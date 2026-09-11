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
! J = CD + (\beta / 2) (CL - {CL_target})^2
!
module actuator_line_objective
  use vector, only: vector_t
  use objective, only: objective_t
  use design, only: design_t
  use actuator_line_design, only: actuator_line_design_t
  use simulation_m, only: simulation_t
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
  use math, only: glsc2, copy, col2, invcol2
  use device_math, only: device_copy, device_glsc2, device_col2, device_invcol2
  use math_ext, only: glsc2_mask
  use field_math, only: field_col3, field_addcol3, field_cmult, field_col2
  use continuation_scheduler, only: nekotop_continuation
  implicit none
  private

  !> An objective function corresponding to drag plus a quadratic lift deviation penalty
  !! \f$ J = CD + (\beta / 2) (CL - {CL_target})^2 \f$
  type, public, extends(objective_t) :: actuator_line_objective_t
     private

     !> Circulation distribution.
     type(vector_t) :: gamma_vec
     !> weight of the lift deviation penalty term
     real(kind=rp) :: lift_penalty_weight
     !> target lift coefficient
     real(kind=rp) :: CL_target

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
    call this%init_base(name, design%size(), weight)

    ! Get the circulation distribution
    call this%gamma_vec%init(design%size())

    select type (design)
    type is (actuator_line_design_t)
       call design%get_values(this%gamma_vec)
       this%lift_penalty_weight = design%lift_penalty_weight
       this%CL_target = design%CL_target

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
    write(force_nondim_factor_name, '("alm_", A, "_", I0)') "force_nondim_factor", alm_id
    force_nondim_factor => neko_registry%get_real_scalar(force_nondim_factor_name)
    lift_target = this%CL_target / force_nondim_factor
    this%value = drag + 0.5_rp * this%lift_penalty_weight * (lift - lift_target)**2

    print *, "actuator_line_update_value"
    print *, this%value

  end subroutine actuator_line_update_value

  !> update_value the sensitivity of the objective function with respect to
  !! \f$chi\f$
  !! @param this The objective.
  !! @param design the design.
  subroutine actuator_line_update_sensitivity(this, design)
    class(actuator_line_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(field_t), pointer :: work
    integer :: temp_indices(1)
    integer :: n_GL, nel
    type(field_t), pointer :: accumulate, fld_GL
    integer :: temp_indices_GL(2)

    ! The Brinkman dissipation adds an extra term in the sensitivity.

   !  call neko_scratch_registry%request_field(work, temp_indices(1), .false.)

   !  if(this%dealias_sensitivity) then
   !     nel = this%c_Xh_GLL%msh%nelv
   !     n_GL = nel * this%Xh_GL%lxyz
   !     call this%scratch_GL%request_field(accumulate, temp_indices_GL(1), &
   !          .false.)
   !     call this%scratch_GL%request_field(fld_GL, temp_indices_GL(2), .false.)

   !     call this%GLL_to_GL%map(fld_GL%x, this%u%x, nel, this%Xh_GL)
   !     call field_col3(accumulate, fld_GL, fld_GL)
   !     call this%GLL_to_GL%map(fld_GL%x, this%v%x, nel, this%Xh_GL)
   !     call field_addcol3(accumulate, fld_GL, fld_GL)
   !     if (this%gdim .eq. 3) then
   !        call this%GLL_to_GL%map(fld_GL%x, this%w%x, nel, this%Xh_GL)
   !        call field_addcol3(accumulate, fld_GL, fld_GL)
   !     end if
   !     ! scale
   !     call field_cmult(accumulate, this%weight * 0.5_rp / this%volume)

   !     ! Evaluate term on GL and preempt the GLL premultiplication
   !     if (NEKO_BCKND_DEVICE .eq. 1) then
   !        call device_col2(accumulate%x_d, this%c_Xh_GL%B_d, n_GL)
   !        call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
   !        call device_invcol2(work%x_d, this%c_Xh_GLL%B_d, work%size())
   !     else
   !        call col2(accumulate%x, this%c_Xh_GL%B, n_GL)
   !        call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
   !        call invcol2(work%x, this%c_Xh_GLL%B, work%size())
   !     end if

   !     call this%scratch_GL%relinquish_field(temp_indices_GL)

   !  else
   !     call field_col3(work, this%u, this%u)
   !     call field_addcol3(work, this%v, this%v)
   !     if (this%gdim .eq. 3) then
   !        call field_addcol3(work, this%w, this%w)
   !     end if
   !     ! scale
   !     call field_cmult(work, this%weight * 0.5_rp / this%volume)
   !  end if

   !  if (NEKO_BCKND_DEVICE .eq. 1) then
   !     call device_copy(this%sensitivity%x_d, work%x_d, this%sensitivity%size())
   !  else
   !     call copy(this%sensitivity%x, work%x, this%sensitivity%size())
   !  end if

   !  call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine actuator_line_update_sensitivity

end module actuator_line_objective
