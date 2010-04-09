!!******************************************************************************
!!
!! module: fields - subroutines to access velocity and magnetic fields
!!
!! Copyright (C) 2009 Grzegorz Kowal <grzegorz@gkowal.info>
!!
!!******************************************************************************
!!
!!  This file is part of GAccel.
!!
!!  GAccel is free software; you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation; either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  GAccel is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!******************************************************************************
!!
!
module fields

  implicit none

! field component dimensions
!
  integer, dimension(3), save :: dm, qm

! domain bounds
!
  real                 , save :: xmin, xmax, ymin, ymax, zmin, zmax

! arrays for storing the electric and magnetic field components
!
  real, dimension(:,:,:), allocatable, save :: ux, uy, uz, bx, by, bz
!
!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! init_fields: subroutine initializes the field variables
!
!===============================================================================
!
  subroutine init_fields()

    use params  , only : fformat
    use fitsio  , only : fits_init, fits_get_dims, fits_get_bounds, fits_read_var
    use hdf5io  , only : hdf5_init, hdf5_get_dims, hdf5_get_bounds, hdf5_read_var

    implicit none

! local variables
!
    integer            :: p, n, ib, ie, jb, je, kb, ke
    real, dimension(6) :: buf

! parameters
!
#ifdef TRICUB
    integer, parameter :: nguard = 2
#else /* TRICUB */
    integer, parameter :: nguard = 1
#endif /* TRICUB */

! local allocatable arrays
!
    real, dimension(:,:,:), allocatable :: tt
!
!-------------------------------------------------------------------------------
!
! initialize component dimensions
!
    dm(:) = 1

! obtain array dimensions
!
    select case(fformat)
    case('fits')
      call fits_init('magx')
      call fits_get_dims(dm)
      call fits_get_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    case('hdf5')
      call hdf5_init()
      call hdf5_get_dims(dm)
      call hdf5_get_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    case default
      write( *, "('ERROR     : ',a,1x,a)" ) "unsupported data format:", fformat
      stop
    end select

! prepare extended dimension
!
    qm(:) = dm(:) + nguard

    ib =     1 - nguard
    ie = dm(1) + nguard
    jb =     1 - nguard
    je = dm(2) + nguard
    kb =     1 - nguard
    ke = dm(3) + nguard

    if (dm(3) .eq. 1) then
      qm(3) = 1
      kb    = 1
      ke    = 1
    end if

! allocate space for field components
!
    allocate(tt(dm(1),dm(2),dm(3)))
#if DIMS == 2
    allocate(ux(ib:ie,jb:je,dm(3)))
    allocate(uy(ib:ie,jb:je,dm(3)))
    allocate(uz(ib:ie,jb:je,dm(3)))
    allocate(bx(ib:ie,jb:je,dm(3)))
    allocate(by(ib:ie,jb:je,dm(3)))
    allocate(bz(ib:ie,jb:je,dm(3)))
#else /* DIMS == 2 */
    allocate(ux(ib:ie,jb:je,kb:ke))
    allocate(uy(ib:ie,jb:je,kb:ke))
    allocate(uz(ib:ie,jb:je,kb:ke))
    allocate(bx(ib:ie,jb:je,kb:ke))
    allocate(by(ib:ie,jb:je,kb:ke))
    allocate(bz(ib:ie,jb:je,kb:ke))
#endif /* DIMS == 2 */

! read field components from the file
!
    write( *, "('INFO      : ',a)" ) "reading velocity and magnetic field"
    select case(fformat)
    case('fits')
      call fits_read_var('velx', tt)
      call expand_array(tt, ux, ib, jb, kb, nguard)
      call fits_read_var('vely', tt)
      call expand_array(tt, uy, ib, jb, kb, nguard)
      call fits_read_var('velz', tt)
      call expand_array(tt, uz, ib, jb, kb, nguard)
      call fits_read_var('magx', tt)
      call expand_array(tt, bx, ib, jb, kb, nguard)
      call fits_read_var('magy', tt)
      call expand_array(tt, by, ib, jb, kb, nguard)
      call fits_read_var('magz', tt)
      call expand_array(tt, bz, ib, jb, kb, nguard)
    case('hdf5')
      call hdf5_read_var('velx', tt)
      call expand_array(tt, ux, ib, jb, kb, nguard)
      call hdf5_read_var('vely', tt)
      call expand_array(tt, uy, ib, jb, kb, nguard)
      call hdf5_read_var('velz', tt)
      call expand_array(tt, uz, ib, jb, kb, nguard)
      call hdf5_read_var('magx', tt)
      call expand_array(tt, bx, ib, jb, kb, nguard)
      call hdf5_read_var('magy', tt)
      call expand_array(tt, by, ib, jb, kb, nguard)
      call hdf5_read_var('magz', tt)
      call expand_array(tt, bz, ib, jb, kb, nguard)
    end select

! deallocate temporary local array
!
    if (allocated(tt)) deallocate(tt)
!
!-------------------------------------------------------------------------------
!
  end subroutine init_fields
!
!===============================================================================
!
! finit_fields: subroutine deallocates the field variables
!
!===============================================================================
!
  subroutine finit_fields()

    implicit none
!
!-------------------------------------------------------------------------------
!
    if (allocated(ux)) deallocate(ux)
    if (allocated(uy)) deallocate(uy)
    if (allocated(uz)) deallocate(uz)
    if (allocated(bx)) deallocate(bx)
    if (allocated(by)) deallocate(by)
    if (allocated(bz)) deallocate(bz)

!-------------------------------------------------------------------------------
!
  end subroutine finit_fields
!
!===============================================================================
!
! get_dimensions: subroutine returns the variable array dimensions
!
!===============================================================================
!
  subroutine get_dimensions(dims)

    implicit none

! output arguments
!
    integer, dimension(3), intent(out) :: dims
!
!-------------------------------------------------------------------------------
!
    dims(:) = dm(:)
!
!-------------------------------------------------------------------------------
!
  end subroutine get_dimensions
!
!===============================================================================
!
! get_domain_bounds: subroutine returns the domain bounds
!
!===============================================================================
!
  subroutine get_domain_bounds(bounds)

    implicit none

! output arguments
!
    real, dimension(3,2), intent(out) :: bounds
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
! expand_array: subroutine expand array by additional boundary layers
!
!===============================================================================
!
  subroutine expand_array(a, b, ib, jb, kb, ng)

    implicit none

! output arguments
!
    integer                                    , intent(in)  :: ib, jb, kb, ng
    real, dimension( 1:dm(1), 1:dm(2), 1:dm(3)), intent(in)  :: a
#if DIMS == 2
#ifdef TRICUB
    real, dimension(-1:qm(1),-1:qm(2), 1      ), intent(out) :: b
#else /* TRICUB */
    real, dimension( 0:qm(1), 0:qm(2), 1      ), intent(out) :: b
#endif /* TRICUB */
#else /* DIMS = 2 */
#ifdef TRICUB
    real, dimension(-1:qm(1),-1:qm(2),-1:qm(3)), intent(out) :: b
#else /* TRICUB */
    real, dimension( 0:qm(1), 0:qm(2), 0:qm(3)), intent(out) :: b
#endif /* TRICUB */
#endif /* DIMS = 2 */

! local variables
!
    integer :: il, jl, kl, iu, ju, ku
!
!-------------------------------------------------------------------------------
!
! copy the interior
!
    b(1:dm(1),1:dm(2),1:dm(3)) = a(1:dm(1),1:dm(2),1:dm(3))

! calculate indices
!
    il = dm(1) + 1 - ng
    iu = dm(1) + 1
    jl = dm(2) + 1 - ng
    ju = dm(2) + 1

! copy X boundaries
!
    b(ib:0    ,1:dm(2),1:dm(3)) = b(il:dm(1),1:dm(2),1:dm(3))
    b(iu:qm(1),1:dm(2),1:dm(3)) = b( 1:ng   ,1:dm(2),1:dm(3))

! copy Y boundaries
!
    b(ib:qm(1),jb:0    ,1:dm(3)) = b(ib:qm(1),jl:dm(2),1:dm(3))
    b(ib:qm(1),ju:qm(2),1:dm(3)) = b(ib:qm(1), 1:ng   ,1:dm(3))

#if DIMS == 3
! copy Z boundaries
!
    if (dm(3) .gt. 1) then
      kl = dm(3) + 1 - ng
      ku = dm(3) + 1
      b(ib:qm(1),jb:qm(2),kb:0    ) = b(ib:qm(1),jb:qm(2),kl:dm(3))
      b(ib:qm(1),jb:qm(2),ku:qm(3)) = b(ib:qm(1),jb:qm(2), 1:ng   )
    end if
#endif /* DIMS == 3 */
!
!-------------------------------------------------------------------------------
!
  end subroutine expand_array
!
end module fields
