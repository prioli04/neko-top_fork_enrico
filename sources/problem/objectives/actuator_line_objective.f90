!> @file actuator_line_objective.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
!> Implements the `actuator_line_objective_t` type.
!
! J = D + (\beta / 2) (L - {L_target})^2
!
module actuator_line_objective
  use vector, only: vector_t
  use matrix, only: matrix_t
  use objective, only: objective_t
  use design, only: design_t
  use actuator_line_design, only: actuator_line_design_t
  use simulation_m, only: simulation_t
  use actuator_line_source_term, only: make_registry_name
  use adjoint_actuator_line_source_term, only: &
       adjoint_actuator_line_source_term_t
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use num_types, only: rp
  use field, only: field_t
  use scratch_registry, only: neko_scratch_registry, scratch_registry_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use mask_ops, only: mask_exterior_const, compute_masked_volume
  use utils, only: neko_error
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use registry, only: neko_registry
  use interpolation, only: interpolator_t
  use registry, only: neko_registry
  use space, only: space_t, GL
  use coefs, only: coef_t
  use math, only: vcross, glsc2, copy, col2, invcol2
  use device_math, only: device_copy, device_glsc2, device_col2, device_invcol2
  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use continuation_scheduler, only: nekotop_continuation
  implicit none
  private

  !> An objective function corresponding to drag plus a quadratic lift deviation penalty
  !! \f$ J = D + (\beta / 2) (L - {L_target})^2 \f$
  type, public, extends(objective_t) :: actuator_line_objective_t
     private

     !> Circulation distribution.
     type(vector_t) :: gamma_vec
     !> weight of the lift deviation penalty term
     real(kind=rp) :: lift_penalty_weight
     !> target lift coefficient
     real(kind=rp) :: L_target

   contains

     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json_sim => &
          actuator_line_init_json_sim
     !> The actual constructor.
     procedure, public, pass(this) :: init_from_attributes => &
          actuator_line_init_attributes
     !> Destructor.
     procedure, public, pass(this) :: free => actuator_line_free
     !> Computes the value of the objective function.
     procedure, public, pass(this) :: update_value => &
          actuator_line_update_value
     !> Computes the sensitivity with respect to the coefficient \f$\chi\f$.
     procedure, public, pass(this) :: update_sensitivity => &
          actuator_line_update_sensitivity

  end type actuator_line_objective_t

contains

  !> The common constructor using a JSON object.
  !! @param this The objective.
  !! @param json the JSON object.
  !! @param design the design.
  !! @param simulation the simulation.
  subroutine actuator_line_init_json_sim(this, json, design, simulation)
    class(actuator_line_objective_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation

    character(len=:), allocatable :: name
    real(kind=rp) :: weight

    call nekotop_continuation%json_get_or_register(json, 'weight', this%weight, weight, 1.0_rp)
    call json_get_or_default(json, "name", name, "Actuator Line")
    call this%init_from_attributes(design, simulation, weight, name)

  end subroutine actuator_line_init_json_sim

  !> The actual constructor.
  !! @param this The objective.
  !! @param design the design.
  !! @param simulation the simulation.
  !! @param weight the weight of the objective function.
  !! @param name the name of the objective.
  subroutine actuator_line_init_attributes(this, design, simulation, &
       weight, name)
    class(actuator_line_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    real(kind=rp), intent(in) :: weight
    character(len=*), intent(in) :: name
    
    type(adjoint_actuator_line_source_term_t) :: actuator_line_adjoint_source

    ! Call the base initializer
    print *, design%size()
    call this%init_base(name, design%size(), weight)

    ! Get the circulation distribution
    call this%gamma_vec%init(design%size())

    select type (design)
    type is (actuator_line_design_t)
       call design%get_values(this%gamma_vec)
       this%lift_penalty_weight = design%lift_penalty_weight
       this%L_target = design%L_target

    class default
       call neko_error('Actuator line objective only works with '// &
            'actuator_line_design')
    end select

  end subroutine actuator_line_init_attributes

  !> Destructor.
  subroutine actuator_line_free(this)
    class(actuator_line_objective_t), intent(inout) :: this
    call this%gamma_vec%free()
    call this%free_base()

  end subroutine actuator_line_free

  !> Compute the objective function.
  !! @param this The objective.
  !! @param design the design.
  subroutine actuator_line_update_value(this, design)
    class(actuator_line_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design

    integer :: alm_id
    real(kind=rp) :: lift, drag, lift_target 
    real(kind=rp), pointer :: force_nondim_factor
    character(len=64) :: force_nondim_factor_name

    ! Get lift and drag coefficients
    select type (design)
    type is (actuator_line_design_t)
       call design%get_resultant_force(lift, drag)
       alm_id = design%alm_id

    class default
       call neko_error('Actuator line objective only works with '// &
            'actuator_line_design')
    end select

    ! Objective: J = drag + (\beta / 2) (lift - lift_target)^2
    this%value = drag + 0.5_rp * this%lift_penalty_weight * (lift - this%L_target)**2

    print *, "actuator_line_update_value"
    print *, this%value

  end subroutine actuator_line_update_value

  !> update_value the sensitivity of the objective function with respect to Gamma
  !! @param this The objective.
  !! @param design the design.
  subroutine actuator_line_update_sensitivity(this, design)
    class(actuator_line_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design

    type(field_t), pointer :: u_adj, v_adj, w_adj, temp_kernel
    type(matrix_t), pointer :: kernel
    type(vector_t), pointer :: u_interp, v_interp, w_interp, resultant_force
    real(kind=rp), allocatable :: velinterp_ex_x(:), velinterp_ex_y(:), velinterp_ex_z(:)
    real(kind=rp), allocatable :: velinterp_ez_x(:), velinterp_ez_y(:), velinterp_ez_z(:)
    real(kind=rp), allocatable :: conv_x(:), conv_y(:), conv_z(:)
    real(kind=rp), allocatable :: velinterp_conv_x(:), velinterp_conv_y(:), velinterp_conv_z(:)
    real(kind=rp), allocatable :: zero_vec(:), one_vec(:)
    real(kind=rp), allocatable :: grad(:, :)
    real(kind=rp) :: delta_L
    integer :: j, n_alm, n_dof, alm_id, temp_index
    character(len=64) :: kernel_name, resultant_force_name

    ! Get interpolated velocities
    select type (design)
    type is (actuator_line_design_t)
       call design%get_adjoint_velocities(u_adj, v_adj, w_adj)
       call design%get_interp_velocities(u_interp, v_interp, w_interp)
       alm_id = design%alm_id

    class default
       call neko_error('Actuator line objective only works with '// &
            'actuator_line_design')
    end select

    n_alm = this%gamma_vec%size()

    ! Allocate arrays
    allocate(velinterp_ex_x(n_alm), velinterp_ex_y(n_alm), velinterp_ex_z(n_alm))
    allocate(velinterp_ez_x(n_alm), velinterp_ez_y(n_alm), velinterp_ez_z(n_alm))
    allocate(conv_x(n_alm), conv_y(n_alm), conv_z(n_alm))
    allocate(velinterp_conv_x(n_alm), velinterp_conv_y(n_alm), velinterp_conv_z(n_alm))
    allocate(zero_vec(n_alm), one_vec(n_alm))
    allocate(grad(n_alm, 3))

    zero_vec = 0.0_rp
    one_vec = 1.0_rp

    ! Get kernel values
    kernel_name = make_registry_name("kernel", alm_id)
    kernel => neko_registry%get_matrix(kernel_name)
    n_dof = size(kernel%x, 1)

    ! Compute lift deviation
    resultant_force_name = make_registry_name("resultant_force", alm_id)
    resultant_force => neko_registry%get_vector(resultant_force_name)
    delta_L = resultant_force%x(1) - this%L_target

    call vcross(velinterp_ex_x, velinterp_ex_y, velinterp_ex_z, &
      u_interp%x, v_interp%x, w_interp%x, one_vec, zero_vec, zero_vec, n_alm) ! Cross product of interpolated velocities and x direction
    call vcross(velinterp_ez_x, velinterp_ez_y, velinterp_ez_z, &
      u_interp%x, v_interp%x, w_interp%x, zero_vec, zero_vec, one_vec, n_alm) ! Cross product of interpolated velocities and z direction

    ! Volumetric convolution of the adjoint velocity with the gaussian kernel
    if (NEKO_BCKND_DEVICE .eq. 1) then
      call neko_scratch_registry%request_field(temp_kernel, temp_index, .false.)

      do j = 1, n_alm
        call device_memcpy(kernel%x(:,j), temp_kernel%x_d, n_dof, HOST_TO_DEVICE, sync=.true.)
        conv_x(j) = device_glsc2(u_adj%x_d, temp_kernel%x_d, n_dof)
        conv_y(j) = device_glsc2(v_adj%x_d, temp_kernel%x_d, n_dof)
        conv_z(j) = device_glsc2(w_adj%x_d, temp_kernel%x_d, n_dof)
      end do

      call neko_scratch_registry%relinquish_field(temp_index)
    else
      do j = 1, n_alm
        conv_x(j) = glsc2(u_adj%x(:,1,1,1), kernel%x(:,j), n_dof)
        conv_y(j) = glsc2(v_adj%x(:,1,1,1), kernel%x(:,j), n_dof)
        conv_z(j) = glsc2(w_adj%x(:,1,1,1), kernel%x(:,j), n_dof)
      end do
    end if

    ! Cross product of interpolated velocities and volumetric convolution
    call vcross(velinterp_conv_x, velinterp_conv_y, velinterp_conv_z, &
      u_interp%x, v_interp%x, w_interp%x, conv_x, conv_y, conv_z, n_alm) 

    ! Compute gradient with respect to gamma_vec
    grad(:, 1) = velinterp_conv_x - velinterp_ex_x - this%lift_penalty_weight * delta_L * velinterp_ez_x
    grad(:, 2) = velinterp_conv_y - velinterp_ex_y - this%lift_penalty_weight * delta_L * velinterp_ez_y
    grad(:, 3) = velinterp_conv_z - velinterp_ex_z - this%lift_penalty_weight * delta_L * velinterp_ez_z

    ! Take only y-component perturbations for now
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(grad(:, 2), this%sensitivity%x_d, n_alm, HOST_TO_DEVICE, sync=.true.)
    else
       call copy(this%sensitivity%x, grad(:, 2), this%sensitivity%size())
    end if

  end subroutine actuator_line_update_sensitivity

end module actuator_line_objective
