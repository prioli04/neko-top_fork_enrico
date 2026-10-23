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
  use actuator_line_source_term, only: actuator_line_source_term_t, make_registry_name
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
     type(global_interpolation_t), public :: interpolator
     !> Forward forcing terms fields
     type(field_list_t) :: fields_forward
     !> Adjoint forcing terms fields
     type(field_list_t) :: fields_adjoint
     !> u of the adjoint
     type(field_t), pointer :: u_adj => null()
     !> v of the adjoint
     type(field_t), pointer :: v_adj => null()
     !> w of the adjoint
     type(field_t), pointer :: w_adj => null()
     !> Gaussian kernel core size
     real(kind=rp), public :: eps
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
     !> Retrieve the actuator line elements' length
     procedure, pass(this) :: get_dy => actuator_line_design_get_dy
     !> Retrieve the actuator line kernels
     procedure, pass(this) :: get_kernels => actuator_line_design_get_kernels
     !> Retrieve adjoint velocity field
     procedure, pass(this) :: get_adjoint_velocities => actuator_line_design_get_adjoint_velocities
     !> Retrieve interpolated velocities
     procedure, pass(this) :: get_interp_velocities => actuator_line_design_get_interp_velocities
     !> Retrieve the sensitivity
     procedure, pass(this) :: get_sensitivity => actuator_line_design_get_sensitivity
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
    character(len=:), allocatable :: name

    integer :: N
    real(kind=rp) :: b, eps, x_center, y_center, z_center
    real(kind=rp), allocatable :: gamma_vec(:)

    call json_get_or_default(parameters, 'name', name, 'Actuator Line Design')
    call json_get(parameters, "N", N)
    call json_get(parameters, 'b', b)
    call json_get(parameters, 'eps', eps)
    call json_get(parameters, 'xcenter', x_center)
    call json_get(parameters, 'ycenter', y_center)
    call json_get(parameters, 'zcenter', z_center)
    call json_get(parameters, 'gamma', gamma_vec)

    ! Initialize and inject into the simulation
    call this%init_from_components(name, simulation, N, b, eps, x_center, y_center, z_center, gamma_vec)

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
    nullify(this%u_adj)
    nullify(this%v_adj)
    nullify(this%w_adj)
    call this%free_base()

  end subroutine actuator_line_design_free

  subroutine actuator_line_design_init_from_components(this, name, simulation, N, b,&
      eps, x_center, y_center, z_center, gamma_vec)
    class(actuator_line_design_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    type(simulation_t), intent(inout) :: simulation
    integer, intent(in) :: N
    real(kind=rp), intent(in) :: b
    real(kind=rp), intent(in) :: eps
    real(kind=rp), intent(in) :: x_center
    real(kind=rp), intent(in) :: y_center
    real(kind=rp), intent(in) :: z_center
    real(kind=rp), intent(in) :: gamma_vec(:)

    type(actuator_line_source_term_t) :: forward_source
    type(adjoint_actuator_line_source_term_t) :: adjoint_source
    type(matrix_t), pointer :: kernel_line => null()
    type(matrix_t), pointer :: kernel_volume => null()
    type(vector_t), pointer :: dy => null()

    call this%init_base(name, N)
    this%u_adj => simulation%adjoint_fluid%u_adj
    this%v_adj => simulation%adjoint_fluid%v_adj
    this%w_adj => simulation%adjoint_fluid%w_adj
    this%gamma_vec = gamma_vec
    this%eps = eps

    ! Init list of source fields used for initializing Neko's sources
    call this%fields_forward%init(3)
    call this%fields_forward%assign(1, simulation%fluid%f_x)
    call this%fields_forward%assign(2, simulation%fluid%f_y)
    call this%fields_forward%assign(3, simulation%fluid%f_z)

    ! Init interpolator
    call this%interpolator%init(simulation%fluid%u%dof)

    ! Init the actuator line term for the forward problem
    call forward_source%init_from_components(this%fields_forward, simulation%fluid%c_Xh, N, &
      b, this%eps, x_center, y_center, z_center, gamma_vec, this%interpolator, this%alm_id)

    ! Append source term to the forward problem
    call simulation%fluid%source_term%add(forward_source)

    ! Retrieve variables created in the forward source term via Neko's registry
    call this%get_dy(dy)
    call this%get_kernels(ptr_volume=kernel_volume, ptr_line=kernel_line)

    ! Init list of source fields used for initializing the adjoint source
    call this%fields_adjoint%init(3)
    call this%fields_adjoint%assign(1, simulation%adjoint_fluid%f_adj_x)
    call this%fields_adjoint%assign(2, simulation%adjoint_fluid%f_adj_y)
    call this%fields_adjoint%assign(3, simulation%adjoint_fluid%f_adj_z)

    ! Init the actuator line source term for the adjoint
    call adjoint_source%init_from_components(this%fields_adjoint, simulation%adjoint_fluid%c_Xh, this%interpolator, &
         simulation%adjoint_fluid%u_adj, simulation%adjoint_fluid%v_adj, simulation%adjoint_fluid%w_adj, &
         dy, kernel_line, kernel_volume, this%gamma_vec, this%eps)
         
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

    values%x = this%gamma_vec

  end subroutine actuator_line_design_get_design

  subroutine actuator_line_design_get_dy(this, dy)
    class(actuator_line_design_t), intent(in) :: this
    type(vector_t), pointer, intent(out) :: dy

    integer :: i, n_alm
    character(len=64) :: dy_name

    dy_name = make_registry_name("dy", this%alm_id)
    dy => neko_registry%get_vector(trim(dy_name))

  end subroutine actuator_line_design_get_dy

  subroutine actuator_line_design_get_kernels(this, ptr_volume, ptr_line)
    class(actuator_line_design_t), intent(in) :: this
    type(matrix_t), pointer, intent(out), optional :: ptr_volume
    type(matrix_t), pointer, intent(out), optional :: ptr_line

    character(len=64) :: kernel_volume_name, kernel_line_name

    if (present(ptr_volume)) then
      kernel_volume_name = make_registry_name("kernel_volume", this%alm_id)
      ptr_volume => neko_registry%get_matrix(trim(kernel_volume_name))
    end if

    if (present(ptr_line)) then
      kernel_line_name = make_registry_name("kernel_line", this%alm_id)
      ptr_line => neko_registry%get_matrix(trim(kernel_line_name))
    end if

  end subroutine actuator_line_design_get_kernels

  subroutine actuator_line_design_get_adjoint_velocities(this, u_adj, v_adj, w_adj)
    class(actuator_line_design_t), intent(in) :: this
    type(field_t), pointer, intent(out) :: u_adj, v_adj, w_adj

    u_adj => this%u_adj
    v_adj => this%v_adj
    w_adj => this%w_adj

  end subroutine actuator_line_design_get_adjoint_velocities

  subroutine actuator_line_design_get_interp_velocities(this, u_interp, v_interp, w_interp)
    class(actuator_line_design_t), intent(in) :: this
    type(vector_t), pointer, intent(out) :: u_interp, v_interp, w_interp

    character(len=64) :: u_name, v_name, w_name

    u_name = make_registry_name("u_interp", this%alm_id)
    v_name = make_registry_name("v_interp", this%alm_id)
    w_name = make_registry_name("w_interp", this%alm_id)

    u_interp => neko_registry%get_vector(trim(u_name))
    v_interp => neko_registry%get_vector(trim(v_name))
    w_interp => neko_registry%get_vector(trim(w_name))

  end subroutine actuator_line_design_get_interp_velocities

  subroutine actuator_line_design_get_sensitivity(this, values)
    class(actuator_line_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: values

    if (size(this%gamma_vec) .ne. values%size()) then
       call neko_error('Get sensitivity: size mismatch')
    end if

    values%x = this%sensitivity

  end subroutine actuator_line_design_get_sensitivity

  subroutine actuator_line_design_get_resultant_force(this, lift, drag)
    class(actuator_line_design_t), intent(in) :: this
    real(kind=rp), intent(out) :: lift
    real(kind=rp), intent(out) :: drag

    character(len=64) :: resultant_force_name
    type(vector_t), pointer :: resultant_force

    resultant_force_name = make_registry_name("resultant_force", this%alm_id)
    resultant_force => neko_registry%get_vector(trim(resultant_force_name))

    lift = resultant_force%x(1)
    drag = resultant_force%x(2)

  end subroutine actuator_line_design_get_resultant_force

  subroutine actuator_line_design_update_design(this, values)
    class(actuator_line_design_t), intent(inout) :: this
    type(vector_t), intent(inout) :: values
    
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(values%x, values%x_d, size(this%gamma_vec), DEVICE_TO_HOST, .true.)
    end if

    this%gamma_vec = values%x

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
