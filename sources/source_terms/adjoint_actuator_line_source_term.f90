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
  use neko_config, only: NEKO_BCKND_DEVICE
  use registry, only: neko_registry
  use scratch_registry, only: neko_scratch_registry
  use json_module, only : json_file
  use time_state, only: time_state_t
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use global_interpolation, only : global_interpolation_t
  use vector, only: vector_t
  use matrix, only : matrix_t
  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use math, only : vcross, glsc3, pi, glsc2, invcol2, glmax, add2s2
  use comm, only: pe_rank
  implicit none
  private
  public :: adjoint_actuator_line_source_term_allocate

  type, public, extends(source_term_t) :: adjoint_actuator_line_source_term_t

     !> The circulation distribution.
     real(kind=rp), allocatable :: gamma_vec(:)
     !> Gaussian kernel core size
     real(kind=rp) :: eps
     !> Gaussian kernel for the line convolution (Gaussian integrated on segments of constant f(x))
     type(matrix_t), pointer :: kernel_line => null()
     !> Gaussian kernel for the volume convolution (Gaussian evaluated at a point)
     type(matrix_t), pointer :: kernel_volume => null()
     !> Actuator line elements' length
     type(vector_t), pointer :: dy => null()
     !> u of the adjoint
     type(field_t), pointer :: u_adj => null()
     !> v of the adjoint
     type(field_t), pointer :: v_adj => null()
     !> w of the adjoint
     type(field_t), pointer :: w_adj => null()
     !> Pointer to the interpolator object
     type(global_interpolation_t), pointer :: interp => null()

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
     !> Computes the adjoint of the interpolation operator
     procedure, private, pass(this) :: adjoint_interpolation_compute
     procedure, private, pass(this) :: adjoint_interpolation_compute2

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
  !! @param u_adj X component of the adjoint velocity field.
  !! @param v_adj Y component of the adjoint velocity field.
  !! @param w_adj Z component of the adjoint velocity field.
  !! @param dy Actuator line elements' length.
  !! @param kernel_line Gaussian kernel for the line convolution (Gaussian integrated on segments of constant f(x)).
  !! @param kernel_volume Gaussian kernel for the volume convolution (Gaussian evaluated at a point).
  !! @param gamma_vec The design circulation vector.
  !! @param eps Gaussian kernel core size.
  subroutine adjoint_actuator_line_source_term_init_from_components(this, fields, coef, interp, u_adj, v_adj, w_adj,&
    dy, kernel_line, kernel_volume, gamma_vec, eps)
    class(adjoint_actuator_line_source_term_t), intent(inout) :: this
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in) :: coef
    type(global_interpolation_t), target, intent(in) :: interp
    type(field_t), intent(in), target :: u_adj, v_adj, w_adj
    type(vector_t), pointer, intent(in) :: dy
    type(matrix_t), pointer, intent(in) :: kernel_line
    type(matrix_t), pointer, intent(in) :: kernel_volume
    real(kind=rp), allocatable, intent(in) :: gamma_vec(:)
    real(kind=rp), intent(in) :: eps

    real(kind=rp) :: start_time, end_time

    ! Mandatory parameters for the general source term
    start_time = 0.0_rp
    end_time = huge(0.0_rp)

    call this%init_base(fields, coef, start_time, end_time)

    ! Point everything in the correct places
    this%u_adj => u_adj
    this%v_adj => v_adj
    this%w_adj => w_adj
    this%interp => interp
    this%dy => dy
    this%kernel_line => kernel_line
    this%kernel_volume => kernel_volume
    this%gamma_vec = gamma_vec
    this%eps = eps

  end subroutine adjoint_actuator_line_source_term_init_from_components

  !> Destructor.
  subroutine adjoint_actuator_line_source_term_free(this)
    class(adjoint_actuator_line_source_term_t), intent(inout) :: this
    if (allocated(this%gamma_vec)) then
      deallocate(this%gamma_vec)
    end if
    nullify(this%dy)
    nullify(this%kernel_line)
    nullify(this%kernel_volume)
    nullify(this%u_adj)
    nullify(this%v_adj)
    nullify(this%w_adj)
    nullify(this%interp)
    call this%free_base()

  end subroutine adjoint_actuator_line_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param this The object.
  !! @param time The time state.
  subroutine adjoint_actuator_line_source_term_compute(this, time)
    class(adjoint_actuator_line_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    
    integer :: n_alm, n_dof, j
    real(kind=rp) :: max_udagger_z
    real(kind=rp), allocatable :: conv_x(:), conv_y(:), conv_z(:)
    real(kind=rp), allocatable :: zero_vec(:), one_vec(:)
    real(kind=rp), allocatable :: f_dagger(:, :)

    n_alm = size(this%gamma_vec)
    n_dof = this%coef%dof%size()

    ! Allocate arrays
    allocate(conv_x(n_alm), conv_y(n_alm), conv_z(n_alm))
    allocate(zero_vec(n_alm), one_vec(n_alm))
    allocate(f_dagger(n_alm, 3))

    zero_vec = 0.0_rp
    one_vec = 1.0_rp

    ! Volumetric adjoint convolution of the adjoint velocity with the gaussian kernel
    do j = 1, n_alm
      ! Usage of glsc3 needed in order to take into account the nodes in other ranks (MPI Reduce)
      conv_x(j) = glsc3(this%u_adj%x(:,1,1,1), this%kernel_volume%x(:,j), this%coef%B, n_dof)
      conv_y(j) = glsc3(this%v_adj%x(:,1,1,1), this%kernel_volume%x(:,j), this%coef%B, n_dof)
      conv_z(j) = glsc3(this%w_adj%x(:,1,1,1), this%kernel_volume%x(:,j), this%coef%B, n_dof)
    end do 

    ! Compute adjoint point forcing
    call vcross(f_dagger(:, 1), f_dagger(:, 2), f_dagger(:, 3), zero_vec, this%gamma_vec * this%dy%x, zero_vec,&
    conv_x, conv_y, conv_z, n_alm) 

    ! Perform the interpolation adjoint (scattering operation) only if the calling rank has points to interpolate
    call this%adjoint_interpolation_compute(f_dagger, n_alm)

    max_udagger_z = glmax(abs(this%w_adj%x), size(this%w_adj%x))
    
    if (pe_rank .eq. 0) then
      print *, "(u_dagger * eta)_z", conv_z
      print *, "max(u_dagger_z)", max_udagger_z

      ! print *, "alm f_adj_x: ", f_dagger(:, 1)
      ! print *, "alm f_adj_z: ", f_dagger(:, 3)
    end if


  end subroutine adjoint_actuator_line_source_term_compute

  !> Computes the adjoint of the interpolation operator.
  !! @param this The object.
  !! @param f_dagger The actuator line localized adjoint forcing
  !! @param n_alm Number of actuator line segments
  subroutine adjoint_interpolation_compute2(this, f_dagger, n_alm)
    class(adjoint_actuator_line_source_term_t), intent(inout) :: this
    integer, intent(in) :: n_alm
    real(kind=rp), intent(in) :: f_dagger(n_alm, 3)
    
    type(field_t), pointer :: fu, fv, fw
    integer :: j, n_dof

    n_dof = this%coef%dof%size()

    ! Get adjoint RHS fields
    fu => this%fields%get_by_index(1)
    fv => this%fields%get_by_index(2)
    fw => this%fields%get_by_index(3)

    do j = 1, n_alm
      call add2s2(fu%x, this%kernel_line%x(:,j), -f_dagger(j,1), n_dof)
      call add2s2(fv%x, this%kernel_line%x(:,j), -f_dagger(j,2), n_dof)
      call add2s2(fw%x, this%kernel_line%x(:,j), -f_dagger(j,3), n_dof)
    end do

    ! call invcol2(fu%x, this%coef%B, n_dof)
    ! call invcol2(fv%x, this%coef%B, n_dof)
    ! call invcol2(fw%x, this%coef%B, n_dof)

    ! Copy forcing to the device
    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_memcpy(fu%x, fu%x_d, n_dof, HOST_TO_DEVICE, .true.)
      call device_memcpy(fv%x, fv%x_d, n_dof, HOST_TO_DEVICE, .true.)
      call device_memcpy(fw%x, fw%x_d, n_dof, HOST_TO_DEVICE, .true.)
    end if

  end subroutine adjoint_interpolation_compute2

  !> Computes the adjoint of the interpolation operator.
  !! @param this The object.
  !! @param f_dagger The actuator line localized adjoint forcing
  !! @param n_alm Number of actuator line segments
  subroutine adjoint_interpolation_compute(this, f_dagger, n_alm)
    class(adjoint_actuator_line_source_term_t), intent(inout) :: this
    integer, intent(in) :: n_alm
    real(kind=rp), intent(in) :: f_dagger(n_alm, 3)
    
    type(field_t), pointer :: fu, fv, fw
    real(kind=rp) :: weight
    real(kind=rp), allocatable :: f_dagger_local(:,:)
    integer :: p, p_glb, e, n_dof
    integer :: i, j, k
    ! integer :: rank, n, recv_size
    integer :: recv_size
    ! integer, pointer :: dof_ids(:), ids(:), p_glb_ids(:)
    integer, pointer :: p_glb_ids(:)

    n_dof = this%fields%item_size(1)

    ! Get adjoint RHS fields
    fu => this%fields%get_by_index(1)
    fv => this%fields%get_by_index(2)
    fw => this%fields%get_by_index(3)

    ! Only proceed with the interpolation adjoint if the calling rank has points to interpolate
    if (this%interp%n_points_local > 0) then

      recv_size = this%interp%glb_intrp_comm%recv_dof(pe_rank)%size()
      
      ! Allocate the forcing terms local to this rank
      allocate(f_dagger_local(recv_size, 3))
      p_glb_ids => this%interp%glb_intrp_comm%recv_dof(pe_rank)%array()

      do p = 1, this%interp%glb_intrp_comm%recv_dof(pe_rank)%size()
        p_glb = p_glb_ids(p)
        f_dagger_local(p, :) = f_dagger(p_glb, :)
      end do

      do p = 1, recv_size
        e = this%interp%el_owner0_local(p) + 1

        do k = 1, this%interp%Xh%lz
          do j = 1, this%interp%Xh%ly
            do i = 1, this%interp%Xh%lx
              weight = this%interp%local_interp%weights_r(i, p) * &
                      this%interp%local_interp%weights_s(j, p) * &
                      this%interp%local_interp%weights_t(k, p)

              ! write(*,'(a,i0,a,5(i0,1x),a,es16.8,a,3(es16.8,1x))') &
              !   'rank=', pe_rank, ' (p,e,i,j,k)=', p, e, i, j, k, &
              !   ' weight=', weight, ' xyz=', &
              !   this%coef%dof%x(i,j,k,e), &
              !   this%coef%dof%y(i,j,k,e), &
              !   this%coef%dof%z(i,j,k,e)

              fu%x(i, j, k, e) = fu%x(i, j, k, e) - weight * f_dagger_local(p, 1)
              fv%x(i, j, k, e) = fv%x(i, j, k, e) - weight * f_dagger_local(p, 2)
              fw%x(i, j, k, e) = fw%x(i, j, k, e) - weight * f_dagger_local(p, 3)
            end do
          end do
        end do
      end do

      call invcol2(fu%x, this%coef%B, n_dof)
      call invcol2(fv%x, this%coef%B, n_dof)
      call invcol2(fw%x, this%coef%B, n_dof)

      ! Copy forcing to the device
      if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_memcpy(fu%x, fu%x_d, n_dof, HOST_TO_DEVICE, .true.)
        call device_memcpy(fv%x, fv%x_d, n_dof, HOST_TO_DEVICE, .true.)
        call device_memcpy(fw%x, fw%x_d, n_dof, HOST_TO_DEVICE, .true.)
      end if
    end if

  end subroutine adjoint_interpolation_compute

end module adjoint_actuator_line_source_term
