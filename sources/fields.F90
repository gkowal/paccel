!!******************************************************************************
!!
!!  This file is part of the PACCEL source code, a program to integrate
!!  test particle trajectories in fields obtained from Newtonian or
!!  relativistic magnetohydrodynamical simulations.
!!
!!  Copyright (C) 2009-2021 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: FIELDS
!!
!!  This module provides allocates field variables and read them from
!!  the data snapshots.
!!
!!******************************************************************************
!
module fields

! module variables are not implicit by default
!
  implicit none

! file format
!
  character(len=32)     , save :: fformat = "fits"

! field component dimensions
!
  integer, dimension(3) , save :: dm = (/ 1, 1, 1 /)
  integer, dimension(3) , save :: lm = (/ 1, 1, 1 /)
  integer, dimension(3) , save :: um = (/ 1, 1, 1 /)

! the number of ghost layers for interpolation
!
  integer               , save :: nghosts = 4

! domain bounds
!
  real                  , save :: xmin = 0.0e+00, xmax = 1.0e+00
  real                  , save :: ymin = 0.0e+00, ymax = 1.0e+00
  real                  , save :: zmin = 0.0e+00, zmax = 1.0e+00

! resistivity
!
  real(kind=8)          , save :: ufac = 1.0d+00
  real(kind=8)          , save :: bfac = 1.0d+00
  real(kind=8)          , save :: ueta = 0.0d+00

! boundary conditions
!
  character(len=64)     , save :: xbndry = "periodic"
  character(len=64)     , save :: ybndry = "periodic"
  character(len=64)     , save :: zbndry = "periodic"

! the domain size
!
  real   , dimension(3), save :: ln = (/ 1.0, 1.0, 1.0 /)

! the cell size
!
  real   , dimension(3), save :: dh = (/ 1.0, 1.0, 1.0 /)

! arrays for storing the field components
!
  real, dimension(:,:,:), save, allocatable :: ux, uy, uz
  real, dimension(:,:,:), save, allocatable :: bx, by, bz
#ifdef CURRENT
  real, dimension(:,:,:), save, allocatable :: jx, jy, jz
#endif /* CURRENT */

! arrays for acceleration rate
!
  integer     , dimension(:,:,:), save, allocatable :: cn
  real(kind=8), dimension(:,:,:), save, allocatable :: ar, mr

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_fields, finalize_fields, read_fields
  public :: get_domain_dimensions, get_domain_bounds, get_domain_sizes
  public :: get_cell_sizes
  public :: nghosts
  public :: ufac, bfac
  public :: ux, uy, uz
  public :: bx, by, bz
#ifdef CURRENT
  public :: jx, jy, jz
#endif /* CURRENT */
  public :: cn, ar, mr
!
!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_FIELDS:
! ----------------------------
!
!   Subroutine reads variable dimensions and allocates memory to store
!   variables.
!
!   Arguments:
!
!     verbose - indicates if it should print any messages;
!     iret    - the return value; if it is 0 everything went successfully,
!               otherwise there was a problem;
!
!===============================================================================
!
  subroutine initialize_fields(verbose, iret)

! include subroutines and variables from other modules
!
    use dataxml   , only : dataxml_init, dataxml_get_dims, dataxml_get_bounds
    use fitsio    , only : fits_init, fits_get_dims, fits_get_bounds
    use parameters, only : get_parameter

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret
!
!-------------------------------------------------------------------------------
!
! print information
!
    if (verbose) then
      write( *, "('INFO      : ',a)" ) "initializing field components"
    end if

! get module parameters
!
    call get_parameter('fformat', fformat)
    call get_parameter("nghosts", nghosts)
    call get_parameter('ueta'   , ueta   )

! get the type of boundaries along each direction
!
    call get_parameter("xbndry", xbndry)
    call get_parameter("ybndry", ybndry)
    call get_parameter("zbndry", zbndry)

! initialize component dimensions
!
    dm(:) = 1

! obtain array dimensions
!
    select case(fformat)
    case('dataxml')
      call dataxml_init()
      call dataxml_get_dims(dm)
      call dataxml_get_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    case('fits')
      call fits_init('magx')
      call fits_get_dims(dm)
      call fits_get_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    case default
      write( *, "('ERROR     : ',a,1x,a)" ) "unsupported data format:", fformat
      stop
    end select

! prepare lower and upper indices
!
    lm(:) =     1 - nghosts
    um(:) = dm(:) + nghosts

! calculate the domain size
!
    ln(1) = xmax - xmin
    ln(2) = ymax - ymin
    ln(3) = zmax - zmin

! calculate the cell sizes
!
    dh(:) = ln(:) / dm(:)

! allocate space for field components
!
    allocate(ux(lm(1):um(1),lm(2):um(2),lm(3):um(3)))
    allocate(uy(lm(1):um(1),lm(2):um(2),lm(3):um(3)))
    allocate(uz(lm(1):um(1),lm(2):um(2),lm(3):um(3)))
    allocate(bx(lm(1):um(1),lm(2):um(2),lm(3):um(3)))
    allocate(by(lm(1):um(1),lm(2):um(2),lm(3):um(3)))
    allocate(bz(lm(1):um(1),lm(2):um(2),lm(3):um(3)))
#ifdef CURRENT
    allocate(jx(lm(1):um(1),lm(2):um(2),lm(3):um(3)))
    allocate(jy(lm(1):um(1),lm(2):um(2),lm(3):um(3)))
    allocate(jz(lm(1):um(1),lm(2):um(2),lm(3):um(3)))
#endif /* CURRENT */

! allocate space for acceleration rate distribution
!
    allocate(cn(lm(1):um(1),lm(2):um(2),lm(3):um(3)))
    allocate(ar(lm(1):um(1),lm(2):um(2),lm(3):um(3)))
    allocate(mr(lm(1):um(1),lm(2):um(2),lm(3):um(3)))

    cn = 0
    ar = 0.0d+00
    mr = 0.0d+00

!-------------------------------------------------------------------------------
!
  end subroutine initialize_fields
!
!===============================================================================
!
! subroutine READ_FIELDS:
! ----------------------
!
!   Subroutine reads variables.
!
!   Arguments:
!
!     verbose - indicates if it should print any messages;
!     iret    - the return value; if it is 0 everything went successfully,
!               otherwise there was a problem;
!
!===============================================================================
!
  subroutine read_fields(verbose, iret)

! include subroutines and variables from other modules
!
    use dataxml, only : dataxml_read_var
    use fitsio , only : fits_read_var

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret

    character(len=90) :: fmt
    real(kind=8)      :: um, bm
!
!-------------------------------------------------------------------------------
!
! print information
!
    if (verbose) then
      write( *, "('INFO      : ',a)" ) "reading field components"
    end if

! read field components from the file
!
    select case(fformat)
    case('dataxml')
      call dataxml_read_var('velx', ux(1:dm(1),1:dm(2),1:dm(3)))
      call dataxml_read_var('vely', uy(1:dm(1),1:dm(2),1:dm(3)))
      call dataxml_read_var('velz', uz(1:dm(1),1:dm(2),1:dm(3)))
      call dataxml_read_var('magx', bx(1:dm(1),1:dm(2),1:dm(3)))
      call dataxml_read_var('magy', by(1:dm(1),1:dm(2),1:dm(3)))
      call dataxml_read_var('magz', bz(1:dm(1),1:dm(2),1:dm(3)))
#ifdef CURRENT
      call dataxml_read_var('curx', jx(1:dm(1),1:dm(2),1:dm(3)))
      call dataxml_read_var('cury', jy(1:dm(1),1:dm(2),1:dm(3)))
      call dataxml_read_var('curz', jz(1:dm(1),1:dm(2),1:dm(3)))
#endif /* CURRENT */
    case('fits')
      call fits_read_var('velx', ux(1:dm(1),1:dm(2),1:dm(3)))
      call fits_read_var('vely', uy(1:dm(1),1:dm(2),1:dm(3)))
      call fits_read_var('velz', uz(1:dm(1),1:dm(2),1:dm(3)))
      call fits_read_var('magx', bx(1:dm(1),1:dm(2),1:dm(3)))
      call fits_read_var('magy', by(1:dm(1),1:dm(2),1:dm(3)))
      call fits_read_var('magz', bz(1:dm(1),1:dm(2),1:dm(3)))
#ifdef CURRENT
      call fits_read_var('curx', jx(1:dm(1),1:dm(2),1:dm(3)))
      call fits_read_var('cury', jy(1:dm(1),1:dm(2),1:dm(3)))
      call fits_read_var('curz', jz(1:dm(1),1:dm(2),1:dm(3)))
#endif /* CURRENT */
    end select

#ifdef CURRENT
! multiply the current density by resistivity
!
    jx(:,:,:) = ueta * jx(:,:,:)
    jy(:,:,:) = ueta * jy(:,:,:)
    jz(:,:,:) = ueta * jz(:,:,:)
#endif /* CURRENT */

! update the ghosts cells
!
    call expand_array(ux)
    call expand_array(uy)
    call expand_array(uz)
    call expand_array(bx)
    call expand_array(by)
    call expand_array(bz)
#ifdef CURRENT
    call expand_array(jx)
    call expand_array(jy)
    call expand_array(jz)
#endif /* CURRENT */

! renormalize velocity field components to [c]
!
    if (ufac /= 1.0d+00) then
      ux(:,:,:) = ufac * ux(:,:,:)
      uy(:,:,:) = ufac * uy(:,:,:)
      uz(:,:,:) = ufac * uz(:,:,:)
    end if

! renormalize magnetic field components to [Gs]
!
    if (bfac /= 1.0d+00) then
      bx(:,:,:) = bfac * bx(:,:,:)
      by(:,:,:) = bfac * by(:,:,:)
      bz(:,:,:) = bfac * bz(:,:,:)
#ifdef CURRENT
      jx(:,:,:) = bfac * jx(:,:,:)
      jy(:,:,:) = bfac * jy(:,:,:)
      jz(:,:,:) = bfac * jz(:,:,:)
#endif /* CURRENT */
    end if

! get the maximum of velocity and magnetic field
!
    um = sqrt(maxval(ux(:,:,:)**2 + uy(:,:,:)**2 + uz(:,:,:)**2))
    bm = sqrt(maxval(bx(:,:,:)**2 + by(:,:,:)**2 + bz(:,:,:)**2))

! print information about maximum velocity and magnetic field
!
    fmt = "('INFO      : maximum plasma velocity is ', 1es12.6,' [c]')"
    if (verbose) write(*, fmt) um
    fmt = "('INFO      : maximum magnetic field is ', 1es12.6,' [Gs]')"
    if (verbose) write(*, fmt) bm

! check if velocity is physical
!
    if (um >= 1.0d+00) then
      if (verbose) &
        write( *, "('WARNING   : ',a)" ) "non-physical plasma velocities!"
      iret = 111
    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_fields
!
!===============================================================================
!
! subroutine FINALIZE_FIELDS:
! --------------------------
!
!   Subroutine finalizes the fields module.
!
!   Arguments:
!
!     verbose - indicates if it should print any messages;
!
!===============================================================================
!
  subroutine finalize_fields(verbose)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in) :: verbose
!
!-------------------------------------------------------------------------------
!
! print information
!
    if (verbose) then
      write( *, "('INFO      : ',a)" ) "deallocating field components"
    end if

! deallocate field arrays
!
    if (allocated(ux)) deallocate(ux)
    if (allocated(uy)) deallocate(uy)
    if (allocated(uz)) deallocate(uz)
    if (allocated(bx)) deallocate(bx)
    if (allocated(by)) deallocate(by)
    if (allocated(bz)) deallocate(bz)
#ifdef CURRENT
    if (allocated(jx)) deallocate(jx)
    if (allocated(jy)) deallocate(jy)
    if (allocated(jz)) deallocate(jz)
#endif /* CURRENT */
    if (allocated(cn)) deallocate(cn)
    if (allocated(ar)) deallocate(ar)
    if (allocated(mr)) deallocate(mr)

!-------------------------------------------------------------------------------
!
  end subroutine finalize_fields
!
!===============================================================================
!
! subroutine GET_DOMAIN_DIMENSIONS:
! --------------------------------
!
!   Subroutine returns domain dimensions.
!
!   Arguments:
!
!     dims - domain dimensions;
!
!===============================================================================
!
  subroutine get_domain_dimensions(dims)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, dimension(3), intent(out) :: dims
!
!-------------------------------------------------------------------------------
!
    dims(:) = dm(:)
!
!-------------------------------------------------------------------------------
!
  end subroutine get_domain_dimensions
!
!===============================================================================
!
! subroutine GET_DOMAIN_BOUNDS:
! ----------------------------
!
!   Subroutine returns domain bounds.
!
!   Arguments:
!
!     bounds - domain bounds;
!
!===============================================================================
!
  subroutine get_domain_bounds(bounds)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8), dimension(3,2), intent(out) :: bounds
!
!-------------------------------------------------------------------------------
!
    bounds(1,1) = xmin
    bounds(1,2) = xmax
    bounds(2,1) = ymin
    bounds(2,2) = ymax
    bounds(3,1) = zmin
    bounds(3,2) = zmax
!
!-------------------------------------------------------------------------------
!
  end subroutine get_domain_bounds
!
!===============================================================================
!
! subroutine GET_DOMAIN_SIZES:
! ----------------------------
!
!   Subroutine returns the domain sizes.
!
!   Arguments:
!
!     ds(:) - the domain sizes in all directions;
!
!===============================================================================
!
  subroutine get_domain_sizes(ds)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real, dimension(3), intent(out) :: ds
!
!-------------------------------------------------------------------------------
!
    ds(:) = ln(:)

!-------------------------------------------------------------------------------
!
  end subroutine get_domain_sizes
!
!===============================================================================
!
! subroutine GET_CELL_SIZES:
! -------------------------
!
!   Subroutine returns the cell sizes.
!
!   Arguments:
!
!     cs(:) - the domain sizes in all directions;
!
!===============================================================================
!
  subroutine get_cell_sizes(cs)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real, dimension(3), intent(out) :: cs
!
!-------------------------------------------------------------------------------
!
    cs(:) = dh(:)

!-------------------------------------------------------------------------------
!
  end subroutine get_cell_sizes
!
!===============================================================================
!
! subroutine EXPAND_ARRAY:
! -----------------------
!
!   Subroutine expands array by filling its ghost zones depending on the type
!   of the boundary condition.
!
!   Arguments:
!
!     u - the expanded array;
!
!===============================================================================
!
  subroutine expand_array(u)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real, dimension(lm(1):um(1),lm(2):um(2),lm(3):um(3)), intent(inout) :: u

! local variables
!
    integer :: l, ls, lt
!
!-------------------------------------------------------------------------------
!
! update the ghosts cells depending on the boundary type
!
    if (trim(xbndry) == "periodic") then
      do l = 1, nghosts
        ls =         l
        lt = dm(1) + l
        u(lt,1:dm(2),1:dm(3)) = u(ls,1:dm(2),1:dm(3))
        ls = dm(1) + 1 - l
        lt =         1 - l
        u(lt,1:dm(2),1:dm(3)) = u(ls,1:dm(2),1:dm(3))
      end do
    else
      do l = 1, nghosts
        ls = dm(1)
        lt = dm(1) + l
        u(lt,1:dm(2),1:dm(3)) = u(ls,1:dm(2),1:dm(3))
        ls = 1
        lt = 1 - l
        u(lt,1:dm(2),1:dm(3)) = u(ls,1:dm(2),1:dm(3))
      end do
    end if
    if (trim(ybndry) == "periodic") then
      do l = 1, nghosts
        ls =         l
        lt = dm(2) + l
        u(lm(1):um(1),lt,1:dm(3)) = u(lm(1):um(1),ls,1:dm(3))
        ls = dm(2) + 1 - l
        lt =         1 - l
        u(lm(1):um(1),lt,1:dm(3)) = u(lm(1):um(1),ls,1:dm(3))
      end do
    else
      do l = 1, nghosts
        ls = dm(2)
        lt = dm(2) + l
        u(lm(1):um(1),lt,1:dm(3)) = u(lm(1):um(1),ls,1:dm(3))
        ls = 1
        lt = 1 - l
        u(lm(1):um(1),lt,1:dm(3)) = u(lm(1):um(1),ls,1:dm(3))
      end do
    end if
    if (trim(zbndry) == "periodic") then
      do l = 1, nghosts
        ls =         l
        lt = dm(3) + l
        u(lm(1):um(1),lm(2):um(2),lt) = u(lm(1):um(1),lm(2):um(2),ls)
        ls = dm(3) + 1 - l
        lt =         1 - l
        u(lm(1):um(1),lm(2):um(2),lt) = u(lm(1):um(1),lm(2):um(2),ls)
      end do
    else
      do l = 1, nghosts
        ls = dm(3)
        lt = dm(3) + l
        u(lm(1):um(1),lm(2):um(2),lt) = u(lm(1):um(1),lm(2):um(2),ls)
        ls = 1
        lt = 1 - l
        u(lm(1):um(1),lm(2):um(2),lt) = u(lm(1):um(1),lm(2):um(2),ls)
      end do
    end if

!-------------------------------------------------------------------------------
!
  end subroutine expand_array

!===============================================================================
!
end module fields
