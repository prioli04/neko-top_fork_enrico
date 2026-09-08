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
  use json_module, only: json_file
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use design, only: design_t
  use simulation_m, only: simulation_t
  use actuator_line_source_term, only: actuator_line_source_term_t
  use adjoint_actuator_line_source_term, only: adjoint_actuator_line_source_term_t
  use vector, only: vector_t
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
     real(kind=rp), allocatable :: gamm(:)
     !> Vector of sensitivities (dJ/dGamma)
     real(kind=rp), allocatable :: sensitivity(:)
     !> Pointer to the actuator line forward source term object
     type(actuator_line_source_term_t), pointer :: forward_source_ptr => null()

   contains

     ! ----------------------------------------------------------------------- !
     ! Initializations

     !> Initialize the design
     generic, public :: init => init_from_json_sim, init_from_components
     !> Initialize the design from a JSON file
     procedure, pass(this), public :: init_from_json_sim => &
          actuator_line_design_init_from_json_sim
     !> Initialize the design from components
     procedure, pass(this), public :: init_from_components => &
          actuator_line_design_init_from_components
     !> Retrieve lift and drag coefficients
     procedure, pass(this), public :: get_coefficients => actuator_line_design_get_coefficients
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
    real(kind=rp) :: CL, b, AR, eps, x_center, y_center, z_center

    call json_get_or_default(parameters, 'name', name, 'Actuator Line Design')
    call json_get(parameters, "N", N)
    call json_get(parameters, 'b', b)
    call json_get(parameters, 'AR', AR)
    call json_get(parameters, 'eps', eps)
    call json_get(parameters, 'xcenter', x_center)
    call json_get(parameters, 'ycenter', y_center)
    call json_get(parameters, 'zcenter', z_center)
    
    allocate(this%gamm(N))
    this%gamm = 1.0_rp

    ! Initialize and inject into the simulation
    call this%init_from_components(name, simulation, N, b, AR, eps, &
        x_center, y_center, z_center)

  end subroutine actuator_line_design_init_from_json_sim

  !> Free the design
  subroutine actuator_line_design_free(this)
    class(actuator_line_design_t), intent(inout) :: this

    if (allocated(this%gamm)) then
      deallocate(this%gamm)
    end if
    if (allocated(this%sensitivity)) then
      deallocate(this%sensitivity)
    end if
    call this%free_base()

  end subroutine actuator_line_design_free

  subroutine actuator_line_design_init_from_components(this, name, simulation, N, b, AR, eps, &
        x_center, y_center, z_center)
    class(actuator_line_design_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(simulation_t), intent(inout) :: simulation
    integer, intent(in) :: N
    real(kind=rp), intent(in) :: b
    real(kind=rp), intent(in) :: AR
    real(kind=rp), intent(in) :: eps
    real(kind=rp), intent(in) :: x_center
    real(kind=rp), intent(in) :: y_center
    real(kind=rp), intent(in) :: z_center

    type(actuator_line_source_term_t), target :: forward_source
    type(adjoint_actuator_line_source_term_t) :: adjoint_source
    type(field_list_t) :: fields

    call this%init_base(name, N)

    ! Init list of source fields used for initializing Neko's sources
    call fields%init(3)
    call fields%assign(1, simulation%fluid%f_x)
    call fields%assign(2, simulation%fluid%f_y)
    call fields%assign(3, simulation%fluid%f_z)

    ! Init the actuator line term for the forward problem
    call forward_source%init_from_compenents(fields, simulation%fluid%c_Xh, &
         N, 1.0_rp, b, AR, eps, x_center, y_center, z_center, this%gamm)
    
    ! Store a pointer to the forward object
    this%forward_source_ptr => forward_source 

    ! Append source term to the forward problem
    call simulation%fluid%source_term%add(forward_source)

    ! Init the actuator line term for the adjoint
    ! call adjoint_source%init_from_components( &
    !      simulation%adjoint_fluid%f_adj_x, &
    !      simulation%adjoint_fluid%f_adj_y, &
    !      simulation%adjoint_fluid%f_adj_z, &
    !      this%brinkman_amplitude, &
    !      simulation%adjoint_fluid%u_adj, &
    !      simulation%adjoint_fluid%v_adj, &
    !      simulation%adjoint_fluid%w_adj, &
    !      simulation%adjoint_fluid%c_Xh, &
    !      simulation%adjoint_fluid%c_Xh_GL, &
    !      simulation%adjoint_fluid%GLL_to_GL, &
    !      dealias, simulation%adjoint_fluid%scratch_GL)
         
    ! ! Append source term to the adjoint problem
    ! select type (f => simulation%adjoint_fluid)
    ! type is (adjoint_fluid_pnpn_t)
    !    call f%source_term%add(adjoint_source)
    ! class default
    ! end select

  end subroutine actuator_line_design_init_from_components

  subroutine actuator_line_design_map_forward(this)
    class(actuator_line_design_t), intent(inout) :: this

  end subroutine actuator_line_design_map_forward

  subroutine actuator_line_design_get_design(this, values)
    class(actuator_line_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: values

    if (size(this%gamm) .ne. values%size()) then
       call neko_error('Get design: size mismatch')
    end if

    call copy(values%x, this%gamm, size(this%gamm))

  end subroutine actuator_line_design_get_design

  subroutine actuator_line_design_get_sensitivity(this, values)
    class(actuator_line_design_t), intent(in) :: this
    type(vector_t), intent(inout) :: values

    if (size(this%gamm) .ne. values%size()) then
       call neko_error('Get sensitivity: size mismatch')
    end if

    call copy(values%x, this%sensitivity, size(this%gamm))

  end subroutine actuator_line_design_get_sensitivity

  subroutine actuator_line_design_get_coefficients(this, CL, CD)
    class(actuator_line_design_t), intent(in) :: this
    real(kind=rp), intent(out) :: CL
    real(kind=rp), intent(out) :: CD

    call this%forward_source_ptr%compute_coefficients(CL, CD)

  end subroutine actuator_line_design_get_coefficients
  
  subroutine actuator_line_design_get_x(this, x)
    class(actuator_line_design_t), intent(in) :: this
    real(kind=rp), intent(out), allocatable :: x(:)
    real(kind=rp), allocatable :: x_vec(:,:)

    call this%forward_source_ptr%get_x_vec(x_vec)
    x = x_vec(:, 1)

  end subroutine actuator_line_design_get_x

  subroutine actuator_line_design_get_y(this, y)
    class(actuator_line_design_t), intent(in) :: this
    real(kind=rp), intent(out), allocatable :: y(:)
    real(kind=rp), allocatable :: x_vec(:,:)

    call this%forward_source_ptr%get_x_vec(x_vec)
    y = x_vec(:, 2)

  end subroutine actuator_line_design_get_y

  subroutine actuator_line_design_get_z(this, z)
    class(actuator_line_design_t), intent(in) :: this
    real(kind=rp), intent(out), allocatable :: z(:)
    real(kind=rp), allocatable :: x_vec(:,:)

    call this%forward_source_ptr%get_x_vec(x_vec)
    z = x_vec(:, 3)

  end subroutine actuator_line_design_get_z

  subroutine actuator_line_design_update_design(this, values)
    class(actuator_line_design_t), intent(inout) :: this
    type(vector_t), intent(inout) :: values
    
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(values%x, values%x_d, size(this%gamm), DEVICE_TO_HOST, .true.)
    end if

    call copy(this%gamm, values%x, size(this%gamm))

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
    print *, this%gamm

  end subroutine actuator_line_design_write

end module actuator_line_design
