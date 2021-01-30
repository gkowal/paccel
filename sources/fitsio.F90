!!******************************************************************************
!!
!!  This file is part of the PACCEL source code, a program to integrate
!!  test particle trajectories in fields obtained from Newtonian or
!!  relativistic magnetohydrodynamical simulations.
!!
!!  Copyright (C) 2007-2021 Grzegorz Kowal <grzegorz@amuncode.org>
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!******************************************************************************
!!
!! module: FITSIO
!!
!!  This module provides subroutines to handle FITS files.
!!
!!******************************************************************************
!
module fitsio

  implicit none

! module parameters
!
  character(len=128)              , save :: idir  = './'
  integer                         , save :: ndims = 0
  integer(kind=8)                 , save :: ncells = 0_8
  integer, dimension(3)           , save :: dims  = 1
  real(kind=8), dimension(2,3)    , save :: bnds  = 0.0d+00
  character(len=4)  , dimension(9), save :: vars  = ''
  character(len=128), dimension(9), save :: vpath = ''
  integer           , dimension(9), save :: kinds = 0

  private

  public :: fitsio_init, fitsio_get_ndims, fitsio_get_dims, fitsio_get_bounds
  public :: fitsio_read_var

!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! subroutine FITSIO_INIT:
! ----------------------
!
!   Subroutine initializes the module FITSIO.
!
!===============================================================================
!
  subroutine fitsio_init()

    use parameters, only : get_parameter

    implicit none

! local variables
!
    logical :: info
    integer :: status, io, p, bsize
!
!-------------------------------------------------------------------------------
!
    status = 0

    call get_parameter('idir', idir)
    p = len_trim(idir)
    if (idir(p:p) /= '/') then
      idir = trim(idir) // '/'
    end if
    call get_parameter('xmin_domain', bnds(1,1))
    call get_parameter('xmax_domain', bnds(2,1))
    call get_parameter('ymin_domain', bnds(1,2))
    call get_parameter('ymax_domain', bnds(2,2))
    call get_parameter('zmin_domain', bnds(1,3))
    call get_parameter('zmax_domain', bnds(2,3))

    if (any(bnds(1,:) >= bnds(2,:))) then
      write(*,"('ERROR',5x,': domain limits are not set correctly!')")
      write(*,"('ERROR',5x,': use xmin_domain, xmax_domain, etc.')")
      stop
    end if


! parse filename information
!
    vars(1:6) = [ 'velx', 'vely', 'velz', 'magx', 'magy', 'magz' ]
    call ftgiou(io, status)
    do p = 1, 6

! generate the file path
!
      vpath(p) = trim(idir) // trim(adjustl(vars(p))) // '.fits'

! check if the file exists
!
      inquire(file = vpath(p), exist = info)
      if (.not. info) then
        write(*,"('ERROR',5x,': file ',a,' does not exist!')") trim(vpath(p))
        stop
      end if

! get some info about data stored in FITS format
!
      call ftopen(io, vpath(p), 0, bsize, status)
      call ftgipr(io, 3, kinds(p), ndims, dims, status)
      call ftclos(io, status)

    end do
    call ftfiou(io, status)

! calculate number of elements in array
!
    ncells = product(dims)

  end subroutine fitsio_init
!
!===============================================================================
!
! subroutine FITSIO_GET_NDIMS:
! ---------------------------
!
!   Subroutine returns the number of the domain dimensions.
!
!===============================================================================
!
  subroutine fitsio_get_ndims(dm)

    implicit none

    integer, intent(out) :: dm
!
!-------------------------------------------------------------------------------
!
    dm = ndims

  end subroutine fitsio_get_ndims
!
!===============================================================================
!
! subroutine FITSIO_GET_DIMS:
! --------------------------
!
!   Subroutine returns the domain dimensions.
!
!===============================================================================
!
  subroutine fitsio_get_dims(dm)

    implicit none

    integer, dimension(:), intent(out) :: dm
!
!-------------------------------------------------------------------------------
!
    dm(:) = dims(:)

  end subroutine fitsio_get_dims
!
!===============================================================================
!
! subroutine FITSIO_GET_BOUNDS:
! ----------------------------
!
!   Subroutine returns the domain bounds.
!
!===============================================================================
!
  subroutine fitsio_get_bounds(xmin, xmax, ymin, ymax, zmin, zmax)

    implicit none

    real(kind=8), intent(out) :: xmin, xmax, ymin, ymax, zmin, zmax

!-------------------------------------------------------------------------------
!
    xmin = bnds(1,1)
    xmax = bnds(2,1)
    ymin = bnds(1,2)
    ymax = bnds(2,2)
    zmin = bnds(1,3)
    zmax = bnds(2,3)

  end subroutine fitsio_get_bounds
!
!===============================================================================
!
! subroutine FITSIO_READ_VAR:
! --------------------------
!
!   Subroutine reads dataset for a given variable.
!
!===============================================================================
!
  subroutine fitsio_read_var(var, qty)

    use iso_fortran_env, only : error_unit

    implicit none

    character(len=*)              , intent(in)    :: var
    real(kind=8), dimension(:,:,:), intent(inout) :: qty

    logical :: info
    integer :: status, io, p, bsize

    real(kind=4), dimension(:,:,:), allocatable :: tmp

    character(len=*), parameter :: loc = 'FITSIO::fitsio_read_var()'

!-------------------------------------------------------------------------------
!
    status = 0

    p = 1
    do while(vars(p) /= var .and. p <= size(vars))
      p = p + 1
    end do
    write(*,"('INFO',6x,': reading from ',a)") trim(vpath(p))

    call ftgiou(io, status)
    call ftopen(io, trim(vpath(p)), 0, bsize, status)
    select case(kinds(p))
    case(-32)
      allocate(tmp(dims(1),dims(2),dims(3)))
      call ftgpve(io, 0, 1, ncells, 0.0e+00, tmp, info, status)
      qty = real(tmp, 8)
      deallocate(tmp)
    case(-64)
      call ftgpvd(io, 0, 1, ncells, 0.0d+00, qty, info, status)
    case default
      write(error_unit,"('[',a,']: ',a)") trim(loc),                           &
                       "Unsupported data format in file '" // trim(vpath(p)) //&
                       "!"
      stop
    end select
    call ftclos(io, status)
    call ftfiou(io, status)

!-------------------------------------------------------------------------------
!
  end subroutine fitsio_read_var

!===============================================================================
!
end module fitsio
