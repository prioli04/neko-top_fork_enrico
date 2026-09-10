!> @file adjoint_actuator_line_source_term.f90
!! @copyright
!! Copyright (c) 2025-2026, The Neko-TOP Authors
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
!> Implements the `adjoint_actuator_line_source_term` type.
module adjoint_actuator_line_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use field, only: field_t
  use registry, only: neko_registry
  use scratch_registry, only: neko_scratch_registry
  use json_module, only : json_file
  use time_state, only: time_state_t
  use json_utils, only: json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use vector, only: vector_t
  use math, only : rzero, add2s2, vcross
  use field_math, only: field_add2s2, field_copy, field_cadd
  use mask_ops, only: mask_exterior_const, compute_masked_volume
  use point_zone, only: point_zone_t
  implicit none
  private
  public :: adjoint_actuator_line_source_term_allocate

  type, public, extends(source_term_t) :: adjoint_actuator_line_source_term_t

     !> The circulation distribution.
     type(vector_t), pointer :: gamma_vec => null()
     !> The lift deviation.
     real(kind=rp) :: delta_L
     !> The penalty factor for the lift deviation.
     real(kind=rp), pointer :: beta

   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => &
          adjoint_actuator_line_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          adjoint_actuator_line_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => adjoint_actuator_line_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => &
          adjoint_actuator_line_source_term_compute
  end type adjoint_actuator_line_source_term_t

contains
  !> Allocator for the adjoint mixing scalar source term.
  subroutine adjoint_actuator_line_source_term_allocate(obj)
    class(source_term_t), allocatable, intent(inout) :: obj
    allocate(adjoint_actuator_line_source_term_t::obj)
  end subroutine adjoint_actuator_line_source_term_allocate

  !> The common constructor using a JSON object.
  !! @param this The object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable where the source term acts.
  subroutine adjoint_actuator_line_source_term_init_from_json(this, &
       json, fields, coef, variable_name)
    class(adjoint_actuator_line_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name

  end subroutine adjoint_actuator_line_source_term_init_from_json

  !> The constructor from type components.
  !! @param this The source term.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param beta The penalty factor for the lift deviation
  !! @param delta_L The lift deviation
  !! @param gamma_vec The design circulation vector
  subroutine adjoint_actuator_line_source_term_init_from_components(this, fields, coef, gamma_vec)
    class(adjoint_actuator_line_source_term_t), intent(inout) :: this
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in) :: coef
    type(vector_t), pointer, intent(in) :: gamma_vec

    real(kind=rp), pointer :: CL_target
    real(kind=rp) :: start_time, end_time, delta_CL

    ! Mandatory parameters for the general source term
    start_time = 0.0_rp
    end_time = huge(0.0_rp)

    call this%init_base(fields, coef, start_time, end_time)

    ! Point everything in the correct places
    this%gamma_vec => gamma_vec
    this%beta => neko_registry%get_real_scalar("alm_lift_penalty_weight")

    ! Compute lift deviation
    CL_target => neko_registry%get_real_scalar("alm_CL_target")


    delta_CL =  - CL_target

  end subroutine adjoint_actuator_line_source_term_init_from_components

  !> Destructor.
  subroutine adjoint_actuator_line_source_term_free(this)
    class(adjoint_actuator_line_source_term_t), intent(inout) :: this

    nullify(this%gamma_vec)
    call this%free_base()

  end subroutine adjoint_actuator_line_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param this The object.
  !! @param time The time state.
  subroutine adjoint_actuator_line_source_term_compute(this, time)
    class(adjoint_actuator_line_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    
    type(field_t), pointer :: fu, fv, fw
    real(kind=rp), allocatable :: f_dagger(:, :)
    real(kind=rp), allocatable :: gamma_ex_x(:), gamma_ex_y(:), gamma_ex_z(:)
    real(kind=rp), allocatable :: gamma_ez_x(:), gamma_ez_y(:), gamma_ez_z(:)
    real(kind=rp), allocatable :: zero_vec(:), one_vec(:)
    integer :: n, n_alm, j

    n = this%fields%item_size(1)
    n_alm = size(this%gamma_vec%x)
    allocate(f_dagger(n_alm, 3))

    allocate(gamma_ex_x(n_alm), gamma_ex_y(n_alm), gamma_ex_z(n_alm))
    allocate(gamma_ez_x(n_alm), gamma_ez_y(n_alm), gamma_ez_z(n_alm))
    allocate(zero_vec(n_alm), one_vec(n_alm))

    ! Get adjoint RHS fields
    fu => this%fields%get_by_index(1)
    fv => this%fields%get_by_index(2)
    fw => this%fields%get_by_index(3)

    ! Clear old source terms
    call rzero(fu%x, n)
    call rzero(fv%x, n)
    call rzero(fw%x, n)

    ! Update source terms
    call vcross(gamma_ex_x, gamma_ex_y, gamma_ex_z, zero_vec, this%gamma_vec%x(j), zero_vec, one_vec, zero_vec, zero_vec, n_alm) ! Cross product of gamma and x direction
    call vcross(gamma_ez_x, gamma_ez_y, gamma_ez_z, zero_vec, this%gamma_vec%x(j), zero_vec, zero_vec, zero_vec, one_vec, n_alm) ! Cross product of gamma and z direction
    
    ! Volumetric convolution of the adjoint veloctiy with the gaussian kernel
    
    ! Cross product of gamma and volumetric convolution
    
    f_dagger(:, 1) = -gamma_ex_x - this%beta * this%delta_L * gamma_ez_x
    f_dagger(:, 2) = -gamma_ex_y - this%beta * this%delta_L * gamma_ez_y
    f_dagger(:, 3) = -gamma_ex_z - this%beta * this%delta_L * gamma_ez_z

  end subroutine adjoint_actuator_line_source_term_compute

end module adjoint_actuator_line_source_term
