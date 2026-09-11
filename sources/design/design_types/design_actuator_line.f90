!> @file design_actuator_line.f90
!! @copyright
!! Copyright (c) 2024-2026, The Neko-TOP Authors
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

! Implements the `actuator_line_design_t` type.
module actuator_line_design
  use num_types, only: rp, sp, dp
  use field_list, only: field_list_t
  use field, only: field_t
  use global_interpolation, only: global_interpolation_t
  use json_module, only: json_file
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use registry, only : neko_registry
  use design, only: design_t
  use simulation_m, only: simulation_t
  use actuator_line_source_term, only: actuator_line_source_term_t
  use adjoint_actuator_line_source_term, only: adjoint_actuator_line_source_term_t
  use vector, only: vector_t
  use matrix, only: matrix_t
  use math, only: copy
  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use json_utils, only: json_get, json_get_or_default
  use utils, only: neko_error
  implicit none
  private

  !> An actuator line circulation design variable
  type, extends(design_t), public :: actuator_line_design_t
     private

     !> Vector of circulations (the design variables)
     real(kind=rp), allocatable :: gamma_vec(:)
     !> Vector of sensitivities (dJ/dGamma)
     real(kind=rp), allocatable :: sensitivity(:)
     !> Global interpolation object
     type(global_interpolation_t) :: interpolator
     !> Weight of the lift deviation penalty term
     real(kind=rp), public :: lift_penalty_weight
     !> Target lift coefficient
     real(kind=rp), public :: CL_target
     !> Actuator line instance id
     integer, public :: alm_id

   contains

     ! ----------------------------------------------------------------------- !
     ! Initializations

     !> Initialize the design
     generic, public :: init => init_from_json_sim, init_from_components
     !> Initialize the design from a JSON file
     procedure, pass(this) :: init_from_json_sim => &
          actuator_line_design_init_from_json_sim
     !> Initialize the design from components
     procedure, pass(this) :: init_from_components => &
          actuator_line_design_init_from_components
     !> Retrieve lift and drag
     procedure, pass(this) :: get_resultant_force => actuator_line_design_get_resultant_force
     !> Retrieve the design variables
     procedure, pass(this) :: get_values => actuator_line_design_get_design
     !> Retrieve the sensitivity
     procedure, pass(this) :: get_sensitivity => actuator_line_design_get_sensitivity
     !> Retrieve the x location of the design variables
    !  procedure, pass(this) :: design_get_x => actuator_line_design_get_x
    !  !> Retrieve the y location of the design variables
    !  procedure, pass(this) :: design_get_y => actuator_line_design_get_y
    !  !> Retrieve the z location of the design variables
    !  procedure, pass(this) :: design_get_z => actuator_line_design_get_z
     !> Update the design
     procedure, pass(this) :: update_design => actuator_line_design_update_design
     !> Destructor
     procedure, pass(this) :: free => actuator_line_design_free
     procedure, pass(this) :: map_forward => actuator_line_design_map_forward
     procedure, pass(this) :: map_backward => actuator_line_design_map_backward
     procedure, pass(this) :: write => actuator_line_design_write

     end type actuator_line_design_t

contains

  !> Initialize the design from a JSON file
  subroutine actuator_line_design_init_from_json_sim(this, parameters, simulation)
    class(actuator_line_design_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(simulation_t), intent(inout) :: simulation
    type(json_file) :: json_subdict
    character(len=:), allocatable :: name

    integer :: N
    real(kind=rp) :: lift_penalty_weight, CL_target, V_inf, b, AR, eps, x_center, y_center, z_center
    real(kind=rp), allocatable :: gamma_vec(:)

    call json_get_or_default(parameters, 'name', name, 'Actuator Line Design')
    call json_get(parameters, "N", N)
    call json_get(parameters, "lift_penalty_weight", lift_penalty_weight)
    call json_get(parameters, "CL_target", CL_target)
    call json_get(parameters, "V_inf", V_inf)
    call json_get(parameters, 'b', b)
    call json_get(parameters, 'AR', AR)
    call json_get(parameters, 'eps', eps)
    call json_get(parameters, 'xcenter', x_center)
    call json_get(parameters, 'ycenter', y_center)
    call json_get(parameters, 'zcenter', z_center)
    call json_get(parameters, 'gamma', gamma_vec)
    
    ! Initialize and inject into the simulation
    call this%init_from_components(name, simulation, N, lift_penalty_weight, CL_target, V_inf, b, AR, eps, &
    x_center, y_center, z_center, gamma_vec)

  end subroutine actuator_line_design_init_from_json_sim

  !> Free the design
  subroutine actuator_line_design_free(this)
    class(actuator_line_design_t), intent(inout) :: this

    if (allocated(this%gamma_vec)) then
      deallocate(this%gamma_vec)
    end if
    if (allocated(this%sensitivity)) then
      deallocate(this%sensitivity)
    end if
    call this%free_base()

  end subroutine actuator_line_design_free

    subroutine actuator_line_design_init_from_components(this, name, simulation, N, lift_penalty_weight, CL_target, V_inf, b, AR, &
      eps, x_center, y_center, z_center, gamma_vec)
    class(actuator_line_design_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    type(simulation_t), intent(inout) :: simulation
    integer, intent(in) :: N
    real(kind=rp), intent(in) :: lift_penalty_weight
    real(kind=rp), intent(in) :: CL_target
    real(kind=rp), intent(in) :: V_inf
    real(kind=rp), intent(in) :: b
    real(kind=rp), intent(in) :: AR
    real(kind=rp), intent(in) :: eps
    real(kind=rp), intent(in) :: x_center
    real(kind=rp), intent(in) :: y_center
    real(kind=rp), intent(in) :: z_center
    real(kind=rp), intent(in) :: gamma_vec(:)

    type(actuator_line_source_term_t) :: forward_source
    type(adjoint_actuator_line_source_term_t) :: adjoint_source
    type(field_list_t) :: fields_forward, fields_adjoint

    call this%init_base(name, N)
    this%gamma_vec = gamma_vec
    this%lift_penalty_weight = lift_penalty_weight
    this%CL_target = CL_target

    ! Init list of source fields used for initializing Neko's sources
    call fields_forward%init(3)
    call fields_forward%assign(1, simulation%fluid%f_x)
    call fields_forward%assign(2, simulation%fluid%f_y)
    call fields_forward%assign(3, simulation%fluid%f_z)

    ! Init interpolator
    call this%interpolator%init(simulation%fluid%u%dof)

    ! Init the actuator line term for the forward problem
    call forward_source%init_from_components(fields_forward, simulation%fluid%c_Xh, &
         N, CL_target, V_inf, b, AR, eps, x_center, y_center, z_center, gamma_vec, this%interpolator, this%alm_id)
    
    ! Append source term to the forward problem
    call simulation%fluid%source_term%add(forward_source)

    ! Init list of source fields used for initializing the adjoint source
    call fields_adjoint%init(3)
    call fields_adjoint%assign(1, simulation%adjoint_fluid%f_adj_x)
    call fields_adjoint%assign(2, simulation%adjoint_fluid%f_adj_y)
    call fields_adjoint%assign(3, simulation%adjoint_fluid%f_adj_z)

    ! Init the actuator line term for the adjoint
    call adjoint_source%init_from_components(fields_adjoint, simulation%adjoint_fluid%c_Xh, &
         simulation%adjoint_fluid%u_adj, simulation%adjoint_fluid%v_adj, simulation%adjoint_fluid%w_adj, &
         this%gamma_vec, lift_penalty_weight, CL_target, this%alm_id)
         
    ! Append source term to the adjoint problem
    select type (f => simulation%adjoint_fluid)
    type is (adjoint_fluid_pnpn_t)
       call f%source_term%add(adjoint_source)
    class default
    end select

  end subroutine actuator_line_design_init_from_components

  subroutine actuator_line_design_map_forward(this)
    class(actuator_line_design_t), intent(inout) :: this

  end subroutine actuator_line_design_map_forward

  subroutine actuator_line_design_get_design(this, values)
    class(actuator_line_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: values

    if (size(this%gamma_vec) .ne. values%size()) then
       call neko_error('Get design: size mismatch')
    end if

    call copy(values%x, this%gamma_vec, size(this%gamma_vec))

  end subroutine actuator_line_design_get_design

  subroutine actuator_line_design_get_sensitivity(this, values)
    class(actuator_line_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: values

    if (size(this%gamma_vec) .ne. values%size()) then
       call neko_error('Get sensitivity: size mismatch')
    end if

    call copy(values%x, this%sensitivity, size(this%gamma_vec))

  end subroutine actuator_line_design_get_sensitivity

  subroutine actuator_line_design_get_resultant_force(this, lift, drag)
    class(actuator_line_design_t), intent(in) :: this
    real(kind=rp), intent(out) :: lift
    real(kind=rp), intent(out) :: drag

    character(len=64) :: resultant_force_name
    type(vector_t), pointer :: resultant_force

    write(resultant_force_name, '("alm_", A, "_", I0)') "resultant_force", this%alm_id
    resultant_force => neko_registry%get_vector(trim(resultant_force_name))

    lift = resultant_force%x(1)
    drag = resultant_force%x(2)

  end subroutine actuator_line_design_get_resultant_force

  subroutine actuator_line_design_get_x(this, x)
    class(actuator_line_design_t), intent(in) :: this
    real(kind=rp), intent(out), allocatable :: x(:)
    character(len=64) :: name
    type(matrix_t), pointer :: x_vec

    write(name, '("alm_", A, "_", I0)') "x_vec", this%alm_id
    x_vec => neko_registry%get_matrix(trim(name))
    x = x_vec%x(:, 1)

  end subroutine actuator_line_design_get_x

  subroutine actuator_line_design_get_y(this, y)
    class(actuator_line_design_t), intent(in) :: this
    real(kind=rp), intent(out), allocatable :: y(:)
    character(len=64) :: name
    type(matrix_t), pointer :: x_vec

    write(name, '("alm_", A, "_", I0)') "x_vec", this%alm_id
    x_vec => neko_registry%get_matrix(trim(name))
    y = x_vec%x(:, 2)

  end subroutine actuator_line_design_get_y

  subroutine actuator_line_design_get_z(this, z)
    class(actuator_line_design_t), intent(in) :: this
    real(kind=rp), intent(out), allocatable :: z(:)
    character(len=64) :: name
    type(matrix_t), pointer :: x_vec

    write(name, '("alm_", A, "_", I0)') "x_vec", this%alm_id
    x_vec => neko_registry%get_matrix(trim(name))
    z = x_vec%x(:, 3)

  end subroutine actuator_line_design_get_z

  subroutine actuator_line_design_update_design(this, values)
    class(actuator_line_design_t), intent(inout) :: this
    type(vector_t), intent(inout) :: values
    
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(values%x, values%x_d, size(this%gamma_vec), DEVICE_TO_HOST, .true.)
    end if

    call copy(this%gamma_vec, values%x, size(this%gamma_vec))

  end subroutine actuator_line_design_update_design

  subroutine actuator_line_design_map_backward(this, sensitivity)
    class(actuator_line_design_t), intent(inout) :: this
    type(vector_t), intent(in) :: sensitivity

    this%sensitivity = sensitivity%x

  end subroutine actuator_line_design_map_backward

  subroutine actuator_line_design_write(this, idx)
    class(actuator_line_design_t), intent(inout) :: this
    integer, intent(in) :: idx

    print *, "actuator_line_design_write"
    print *, this%gamma_vec

  end subroutine actuator_line_design_write

end module actuator_line_design
