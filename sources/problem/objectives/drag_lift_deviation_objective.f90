!> @file drag_lift_deviation_objective.f90
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
!> Implements the `drag_lift_deviation_objective_t` type.
!
! J = D_ind + (\beta / 2) (L - {L_target})^2
!
module drag_lift_deviation_objective
  use global_interpolation, only: global_interpolation_t
  use field_list, only: field_list_t
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use vector, only: vector_t
  use matrix, only: matrix_t
  use objective, only: objective_t
  use design, only: design_t
  use actuator_line_design, only: actuator_line_design_t
  use simulation_m, only: simulation_t
  use adjoint_drag_lift_deviation_source_term, only: adjoint_drag_lift_deviation_source_term_t
  use num_types, only: rp
  use field, only: field_t
  use scratch_registry, only: neko_scratch_registry
  use neko_config, only: NEKO_BCKND_DEVICE
  use utils, only: neko_error
  use json_module, only: json_file
  use json_utils, only: json_get_or_default, json_get
  use registry, only: neko_registry
  use registry, only: neko_registry
  use coefs, only: coef_t
  use math, only: vcross, glsc3, copy, pi, glsc2
  use device_math, only: device_copy, device_glsc3
  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use continuation_scheduler, only: nekotop_continuation
  use comm, only: pe_rank
  implicit none
  private

  !> An objective function corresponding to drag plus a quadratic lift deviation penalty
  !! \f$ J = D_ind + (\beta / 2) (L - {L_target})^2 \f$
  type, public, extends(objective_t) :: drag_lift_deviation_objective_t
     private
     !> Adjoint forcing terms fields
     type(field_list_t) :: fields_adjoint
     !> Circulation distribution.
     type(vector_t) :: gamma_vec
     !> Actuator line elements' length
     type(vector_t), pointer :: dy => null()
     !> Gaussian kernel for the line convolution (Gaussian integrated on segments of constant f(x))
     type(matrix_t), pointer :: kernel_line => null()
     !> Gaussian kernel for the volume convolution (Gaussian evaluated at a point)
     type(matrix_t), pointer :: kernel_volume => null()
     !> Gaussian kernel core size
     real(kind=rp) :: eps
     !> weight of the lift deviation penalty term
     real(kind=rp) :: lift_penalty_weight
     !> target lift
     real(kind=rp) :: L_target
     !> Coefficients for the SEM.
     type(coef_t), pointer :: coef => null()

   contains

     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json_sim => &
          drag_lift_deviation_init_json_sim
     !> The actual constructor.
     procedure, public, pass(this) :: init_from_attributes => &
          drag_lift_deviation_init_attributes
     !> Destructor.
     procedure, public, pass(this) :: free => drag_lift_deviation_free
     !> Computes the value of the objective function.
     procedure, public, pass(this) :: update_value => &
          drag_lift_deviation_update_value
     !> Computes the sensitivity with respect to the coefficient \f$\chi\f$.
     procedure, public, pass(this) :: update_sensitivity => &
          drag_lift_deviation_update_sensitivity

  end type drag_lift_deviation_objective_t

contains

  !> The common constructor using a JSON object.
  !! @param this The objective.
  !! @param json the JSON object.
  !! @param design the design.
  !! @param simulation the simulation.
  subroutine drag_lift_deviation_init_json_sim(this, json, design, simulation)
    class(drag_lift_deviation_objective_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation

    character(len=:), allocatable :: name
    real(kind=rp) :: weight, lift_penalty_weight, L_target

    call nekotop_continuation%json_get_or_register(json, 'weight', this%weight, weight, 1.0_rp)
    call json_get_or_default(json, "name", name, "Drag with lift deviation")
    call json_get(json, "lift_penalty_weight", lift_penalty_weight)
    call json_get(json, "L_target", L_target)
    call this%init_from_attributes(design, simulation, weight, name, lift_penalty_weight, L_target)

  end subroutine drag_lift_deviation_init_json_sim

  !> The actual constructor.
  !! @param this The objective.
  !! @param design the design.
  !! @param simulation the simulation.
  !! @param weight the weight of the objective function.
  !! @param name the name of the objective.
  !! @param lift_penalty_weight weight of the lift deviation penalty term.
  !! @param L_target target lift.
  subroutine drag_lift_deviation_init_attributes(this, design, simulation, &
       weight, name, lift_penalty_weight, L_target)
    class(drag_lift_deviation_objective_t), intent(inout) :: this
    class(design_t), target, intent(in) :: design
    type(simulation_t), target, intent(inout) :: simulation
    real(kind=rp), intent(in) :: weight
    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: lift_penalty_weight
    real(kind=rp), intent(in) :: L_target

    integer :: alm_id
    type(global_interpolation_t), pointer :: interpolator
    type(adjoint_drag_lift_deviation_source_term_t) :: drag_lift_deviation

    ! Call the base initializer
    call this%init_base(name, design%size(), weight)
    this%lift_penalty_weight = lift_penalty_weight
    this%L_target = L_target

    ! Get the circulation distribution
    call this%gamma_vec%init(design%size())

    ! Get adjoint mass matrix coefficients
    this%coef => simulation%adjoint_case%fluid_adj%c_Xh

    select type (design)
    type is (actuator_line_design_t)
       call design%get_values(this%gamma_vec)
       call design%get_dy(this%dy)
       call design%get_kernels(ptr_volume=this%kernel_volume, ptr_line=this%kernel_line)
       this%eps = design%eps
       interpolator => design%interpolator
       alm_id = design%alm_id

    class default
       call neko_error('Drag with lift deviation objective only works with '// &
            'actuator_line_design')
    end select

    ! Init list of source fields used for initializing the adjoint source
    call this%fields_adjoint%init(3)
    call this%fields_adjoint%assign(1, simulation%adjoint_fluid%f_adj_x)
    call this%fields_adjoint%assign(2, simulation%adjoint_fluid%f_adj_y)
    call this%fields_adjoint%assign(3, simulation%adjoint_fluid%f_adj_z)

    call drag_lift_deviation%init_from_components(this%fields_adjoint, simulation%adjoint_fluid%c_Xh, interpolator, &
         this%dy, this%gamma_vec%x, this%kernel_line, this%lift_penalty_weight, this%L_target, alm_id)

    ! append adjoint forcing term based on objective function
    select type (f => simulation%adjoint_fluid)
    type is (adjoint_fluid_pnpn_t)
       call f%source_term%add_source_term(drag_lift_deviation)
    end select

  end subroutine drag_lift_deviation_init_attributes

  !> Destructor.
  subroutine drag_lift_deviation_free(this)
    class(drag_lift_deviation_objective_t), intent(inout) :: this
    call this%gamma_vec%free()
    nullify(this%dy)
    nullify(this%kernel_line)
    nullify(this%kernel_volume)
    nullify(this%coef)
    call this%free_base()

  end subroutine drag_lift_deviation_free

  !> Compute the objective function.
  !! @param this The objective.
  !! @param design the design.
  subroutine drag_lift_deviation_update_value(this, design)
    class(drag_lift_deviation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design

    real(kind=rp) :: lift, drag

    ! Get lift and drag coefficients
    select type (design)
    type is (actuator_line_design_t)
       call design%get_resultant_force(lift, drag)

    class default
       call neko_error('Drag with lift deviation objective only works with '// &
            'actuator_line_design')
    end select

    ! Objective: J = D_ind + (\beta / 2) (L - L_target)^2
    this%value = drag + 0.5_rp * this%lift_penalty_weight * (lift - this%L_target)**2

    print *, "Drag", drag
    print *, "Lift", lift
    print *, "Objective value", this%value
    print *, "Objective weight", this%weight

  end subroutine drag_lift_deviation_update_value

  !> update_value the sensitivity of the objective function with respect to Gamma
  !! @param this The objective.
  !! @param design the design.
  subroutine drag_lift_deviation_update_sensitivity(this, design)
    class(drag_lift_deviation_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design

    type(field_t), pointer :: u_adj, v_adj, w_adj
    type(vector_t), pointer :: u_interp, v_interp, w_interp
    real(kind=rp), allocatable :: velinterp_ex_x(:), velinterp_ex_y(:), velinterp_ex_z(:)
    real(kind=rp), allocatable :: velinterp_ez_x(:), velinterp_ez_y(:), velinterp_ez_z(:)
    real(kind=rp), allocatable :: conv_x(:), conv_y(:), conv_z(:)
    real(kind=rp), allocatable :: velinterp_conv_x(:), velinterp_conv_y(:), velinterp_conv_z(:)
    real(kind=rp), allocatable :: zero_vec(:), one_vec(:)
    real(kind=rp), allocatable :: grad(:, :)
    real(kind=rp) :: delta_L, lift, drag
    integer :: j, n_alm, n_dof

    ! Get interpolated velocities
    select type (design)
    type is (actuator_line_design_t)
       call design%get_adjoint_velocities(u_adj, v_adj, w_adj)
       call design%get_interp_velocities(u_interp, v_interp, w_interp)
       call design%get_resultant_force(lift, drag)

    class default
       call neko_error('Drag with lift deviation objective only works with '// &
            'actuator_line_design')
    end select

    n_alm = this%gamma_vec%size()
    n_dof = this%coef%dof%size()

    ! Allocate arrays
    allocate(velinterp_ex_x(n_alm), velinterp_ex_y(n_alm), velinterp_ex_z(n_alm))
    allocate(velinterp_ez_x(n_alm), velinterp_ez_y(n_alm), velinterp_ez_z(n_alm))
    allocate(conv_x(n_alm), conv_y(n_alm), conv_z(n_alm))
    allocate(velinterp_conv_x(n_alm), velinterp_conv_y(n_alm), velinterp_conv_z(n_alm))
    allocate(zero_vec(n_alm), one_vec(n_alm))
    allocate(grad(n_alm, 3))

    zero_vec = 0.0_rp
    one_vec = 1.0_rp

    ! Compute lift deviation
    delta_L = lift - this%L_target

    call vcross(velinterp_ex_x, velinterp_ex_y, velinterp_ex_z, &
      u_interp%x, v_interp%x, w_interp%x, one_vec, zero_vec, zero_vec, n_alm) ! Cross product of interpolated velocities and x direction
    call vcross(velinterp_ez_x, velinterp_ez_y, velinterp_ez_z, &
      u_interp%x, v_interp%x, w_interp%x, zero_vec, zero_vec, one_vec, n_alm) ! Cross product of interpolated velocities and z direction

    ! Volumetric adjoint convolution of the adjoint velocity with the gaussian kernel
    do j = 1, n_alm
      ! Usage of glsc3 needed in order to take into account the nodes in other ranks (MPI Reduce)
      conv_x(j) = glsc3(u_adj%x(:,1,1,1), this%kernel_volume%x(:,j), this%coef%B, n_dof)
      conv_y(j) = glsc3(v_adj%x(:,1,1,1), this%kernel_volume%x(:,j), this%coef%B, n_dof)
      conv_z(j) = glsc3(w_adj%x(:,1,1,1), this%kernel_volume%x(:,j), this%coef%B, n_dof)
    end do 

    ! Cross product of interpolated velocities and volumetric convolution
    call vcross(velinterp_conv_x, velinterp_conv_y, velinterp_conv_z, &
      u_interp%x, v_interp%x, w_interp%x, conv_x, conv_y, conv_z, n_alm) 

    ! Compute gradient with respect to gamma_vec
    grad(:, 1) = (velinterp_conv_x - velinterp_ex_x - this%lift_penalty_weight*delta_L*velinterp_ez_x) * this%dy%x
    grad(:, 2) = (velinterp_conv_y - velinterp_ex_y - this%lift_penalty_weight*delta_L*velinterp_ez_y) * this%dy%x
    grad(:, 3) = (velinterp_conv_z - velinterp_ex_z - this%lift_penalty_weight*delta_L*velinterp_ez_z) * this%dy%x

    if (pe_rank .eq. 0) then
      print *, "(u_int * ex) * dy:", velinterp_ex_y * this%dy%x
      print *, "(u_dagger * eta)_ z: ", conv_z 
      print *, "(u_int x (u_dagger * eta)_z) * dy: ", velinterp_conv_y * this%dy%x
      print *, "gradient: ", grad(:, 2)
    end if

    ! Take only y-component perturbations for now
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(grad(:, 2), this%sensitivity%x_d, n_alm, HOST_TO_DEVICE, sync=.true.)
    else
       call copy(this%sensitivity%x, grad(:, 2), this%sensitivity%size())
    end if

  end subroutine drag_lift_deviation_update_sensitivity

end module drag_lift_deviation_objective
