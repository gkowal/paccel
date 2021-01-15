!!******************************************************************************
!!
!!  This file is part of the PACCEL source code, a program to integrate
!!  test particle trajectories in fields obtained from Newtonian or
!!  relativistic magnetohydrodynamical simulations.
!!
!!  Copyright (C) 2013-2021 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: COORDINATES
!!
!!  This module provides subroutines to convert between domain coordinates.
!!
!!******************************************************************************
!
module coordinates

  implicit none

! field component dimensions
!
  integer, dimension(3)       , save :: dm     = (/ 1, 1, 1 /)

! domain bounds
!
  real(kind=8), dimension(3,2), save :: bounds = 0.0d+00

! domain sizes
!
  real(kind=8), dimension(3)  , save :: dlen   = 1.0d+00
  real(kind=8), dimension(3)  , save :: plen   = 1.0d+00

! periodic boundary indicator
!
  logical, dimension(3)       , save :: pbnd   = .true.

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_coordinates
  public :: map_position_to_index, map_normalized_to_position
  public :: is_inside
!
!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_COORDINATES:
! ---------------------------------
!
!   Subroutine initializes module parameters and variables.
!
!   Arguments:
!
!     verbose - indicates if it should print any messages;
!     status  - the return value; if it is 0 everything went successfully,
!               otherwise there was a problem;
!
!===============================================================================
!
  subroutine initialize_coordinates(verbose, status)

! import required modules
!
    use fields    , only : get_domain_dimensions, get_domain_bounds
    use parameters, only : get_parameter

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: status

! local variables
!
    integer                :: i
    character(len=255)     :: xbndry = "periodic"
    character(len=255)     :: ybndry = "periodic"
    character(len=255)     :: zbndry = "periodic"
!
!-------------------------------------------------------------------------------
!
    status = 0

! print info
!
    if (verbose) write(*, "('INFO',6x,': initializing coordinates')")

! get the domain dimensions
!
    call get_domain_dimensions(dm(:))

! prepare the index coordinate limits
!
    plen(1)   = 1.0d+00 * dm(1)
    plen(2)   = 1.0d+00 * dm(2)
    plen(3)   = 1.0d+00 * dm(3)

! get domain bounds
!
    call get_domain_bounds(bounds(:,:))

! calculate the domain sizes
!
    dlen(1:3) = bounds(1:3,2) - bounds(1:3,1)

! get the type of boundaries along each direction
!
    call get_parameter("xbndry", xbndry)
    call get_parameter("ybndry", ybndry)
    call get_parameter("zbndry", zbndry)

! set the conversion procedure pointers
!
    pbnd(1) = trim(xbndry) == "periodic"
    pbnd(2) = trim(ybndry) == "periodic"
    pbnd(3) = trim(zbndry) == "periodic"

! print information about the coordinates
!
    if (verbose) then
      write( *, "('INFO',6x,': x-boundary = ', a)" ) trim(xbndry)
      write( *, "('INFO',6x,': y-boundary = ', a)" ) trim(ybndry)
      write( *, "('INFO',6x,': z-boundary = ', a)" ) trim(zbndry)
    end if

!-------------------------------------------------------------------------------
!
  end subroutine initialize_coordinates
!
!===============================================================================
!
! subroutine MAP_POSITION_TO_INDEX:
! --------------------------------
!
!   Subroutine converts a physical position to corresponding array index
!   within the computational domain.
!
!   Arguments:
!
!     x - the position vector in the domain coordinates;
!     p - the index vector in the domain dimensions;
!
!===============================================================================
!
  subroutine map_position_to_index(x, p)

    implicit none

! subroutine arguments
!
    real(kind=8), dimension(3), intent(in)  :: x
    real(kind=8), dimension(3), intent(out) :: p

! local variables
!
    integer                    :: i
    real(kind=8), dimension(3) :: t
!
!-------------------------------------------------------------------------------
!
    t(1:3) = (x(1:3) - bounds(1:3,1)) / dlen(1:3)
    do i = 1, 3
      if (pbnd(i)) then
        t(i) = t(i) - floor(t(i))
      else
        t(i) = min(1.0d+00, max(0.0d+00, t(i)))
      end if
    end do
    p(1:3) = plen(1:3) * t(1:3) + 0.5d+00

!-------------------------------------------------------------------------------
!
  end subroutine map_position_to_index
!
!===============================================================================
!
! subroutine MAP_NORMALIZED_TO_POSITION:
! -------------------------------------
!
!   Subroutine converts a normalized position (0.0 to 1.0) to its physical
!   representation.
!
!   Arguments:
!
!     r - the index vector in the domain dimensions;
!     x - the position vector in the domain coordinates;
!
!===============================================================================
!
  subroutine map_normalized_to_position(r, x)

    implicit none

! output arguments
!
    real(kind=8), dimension(3), intent(in)  :: r
    real(kind=8), dimension(3), intent(out) :: x
!
!-------------------------------------------------------------------------------
!
    x(1:3) = dlen(1:3) * r(1:3) + bounds(1:3,1)

!-------------------------------------------------------------------------------
!
  end subroutine map_normalized_to_position
!
!===============================================================================
!
! function IS_INSIDE:
! ------------------
!
!   Function checks if the point lays inside the selected region.
!
!   Arguments:
!
!     x - the physical position of particle;
!
!===============================================================================
!
  logical function is_inside(x) result(inside)

    implicit none

! subroutine arguments
!
    real(kind=8), dimension(3), intent(in) :: x

! local variables
!
    integer :: i
!
!-------------------------------------------------------------------------------
!
! reset the inside flag
!
    inside = .true.

! check if the position is inside the domain
!
    do i = 1, 3
      if (.not. pbnd(i)) then
        inside = inside .and. (x(i) >= bounds(i,1) .and. x(i) <= bounds(i,2))
      end if
    end do

!-------------------------------------------------------------------------------
!
  end function is_inside

!===============================================================================
!
end module coordinates
