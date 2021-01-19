!!******************************************************************************
!!
!!  This file is part of the PACCEL source code, a program to integrate
!!  test particle trajectories in fields obtained from Newtonian or
!!  relativistic magnetohydrodynamical simulations.
!!
!!  Copyright (C) 2021 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: DataXML
!!
!!  This module provides subroutines to handle DataXML files.
!!
!!******************************************************************************
!
module dataxml

  implicit none

! interfaces to compression algorithms
!
#ifdef ZSTD
  interface

    integer(c_size_t) function zstd_decompress(dst, dstCapacity,               &
                                               src, compressedSize)            &
                               bind(C, name="ZSTD_decompress")
      use iso_c_binding, only: c_size_t, c_int, c_ptr
      implicit none
      integer(kind=c_size_t), value :: dstCapacity, compressedSize
      type(c_ptr)           , value :: src, dst
    end function zstd_decompress

  end interface
#endif /* ZSTD */

! interfaces to local subroutines
!
  interface get_attribute
    module procedure get_attribute_integer
    module procedure get_attribute_double
  end interface get_attribute

! module variables
!
  character(len=128)              , save :: idir  = './'
  integer                         , save :: ndims = 0
  integer, dimension(3)           , save :: dims  = 1
  real(kind=8), dimension(2,3)    , save :: bnds  = 0.0d+00
  character(len=4)  , dimension(9), save :: vars  = ''
  character(len=128), dimension(9), save :: vpath = ''
  integer           , dimension(9), save :: kinds = 0
  integer           , dimension(9), save :: comp  = 0
  integer(kind=8)   , dimension(9), save :: usize = 0
  integer(kind=8)   , dimension(9), save :: csize = 0
  integer(kind=8)   , dimension(9), save :: uhash = 0
  integer(kind=8)   , dimension(9), save :: chash = 0

  private

  public :: dataxml_init, dataxml_get_ndims, dataxml_get_dims
  public :: dataxml_get_bounds
  public :: dataxml_read_var
  public :: dataxml_store_last_state, dataxml_restore_last_state

!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! subroutine DATAXML_INIT:
! -----------------------
!
!   Subroutine initializes the module DataXML.
!
!===============================================================================
!
  subroutine dataxml_init()

    use iso_fortran_env, only : error_unit
    use parameters     , only : get_parameter

    implicit none

    character(len=256) :: fname
    logical            :: status
    integer            :: io, p

    character(len=:), allocatable :: buffer

    character(len=*), parameter :: loc = 'DATAXML::dataxml_init()'

!-------------------------------------------------------------------------------
!
    call get_parameter('idir', idir)
    p = len_trim(idir)
    if (idir(p:p) /= '/') then
      idir = trim(idir) // '/'
    end if

    write(fname,"(a,'datasets.xml')") trim(idir)

! check if the file exists and get its size
!
    inquire(file = fname, exist = status, size = p)
    if (.not. status) then
      write(*,'("ERROR   : file ",a," does not exist!")') trim(fname)
      stop
    end if

! allocate a bufer and store the content of the metadata file in it
!
    allocate(character(len=p) :: buffer)
    open(newunit = io, file = fname, form = 'unformatted', access = 'direct',  &
                       recl = p)
    read(io, rec = 1) buffer(:)
    close(io)

! parse the buffer to get the required information
!
    call get_attribute(buffer, 'ndim'  , ndims)
    call get_attribute(buffer, 'dim(1)', dims(1))
    call get_attribute(buffer, 'dim(2)', dims(2))
    call get_attribute(buffer, 'dim(3)', dims(3))
    call get_attribute(buffer, 'xmin'  , bnds(1,1))
    call get_attribute(buffer, 'xmax'  , bnds(2,1))
    call get_attribute(buffer, 'ymin'  , bnds(1,2))
    call get_attribute(buffer, 'ymax'  , bnds(2,2))
    call get_attribute(buffer, 'zmin'  , bnds(1,3))
    call get_attribute(buffer, 'zmax'  , bnds(2,3))

! parse filename information
!
    vars(1:6) = [ 'velx', 'vely', 'velz', 'magx', 'magy', 'magz' ]
    do p = 1, 6
      call get_filename(buffer, vars(p), vpath(p), kinds(p), comp(p),          &
                                      usize(p), csize(p), uhash(p), chash(p))
      vpath(p) = trim(idir) // trim(adjustl(vpath(p)))
    end do

! deallocate the buffer
!
    deallocate(buffer)

  end subroutine dataxml_init
!
!===============================================================================
!
! subroutine DATAXML_GET_NDIMS:
! ----------------------------
!
!   Subroutine returns the number of the domain dimensions.
!
!===============================================================================
!
  subroutine dataxml_get_ndims(dm)

    implicit none

    integer, intent(out) :: dm
!
!-------------------------------------------------------------------------------
!
    dm = ndims

  end subroutine dataxml_get_ndims
!
!===============================================================================
!
! subroutine DATAXML_GET_DIMS:
! ---------------------------
!
!   Subroutine returns the domain dimensions.
!
!===============================================================================
!
  subroutine dataxml_get_dims(dm)

    implicit none

    integer, dimension(3), intent(out) :: dm
!
!-------------------------------------------------------------------------------
!
    dm(:) = dims(:)

  end subroutine dataxml_get_dims
!
!===============================================================================
!
! subroutine DATAXML_GET_BOUNDS:
! -----------------------------
!
!   Subroutine returns the domain bounds.
!
!===============================================================================
!
  subroutine dataxml_get_bounds(xmin, xmax, ymin, ymax, zmin, zmax)

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

  end subroutine dataxml_get_bounds
!
!===============================================================================
!
! subroutine DATAXML_READ_VAR:
! ---------------------------
!
!   Subroutine reads dataset for a given variable.
!
!===============================================================================
!
  subroutine dataxml_read_var(var, qty)

    use hash           , only : xxh64
    use iso_c_binding  , only : c_int, c_loc
    use iso_fortran_env, only : error_unit

    implicit none

    character(len=*)              , intent(in)    :: var
    real(kind=8), dimension(:,:,:), intent(inout) :: qty

    integer                 :: io, p
    integer(kind=8), target :: rsize

    integer(kind=1), dimension(:), allocatable, target :: buffer, input

    character(len=*), parameter :: loc = 'DATAXML::dataxml_read_var()'

!-------------------------------------------------------------------------------
!
    p = 1
    do while(vars(p) /= var .and. p <= size(vars))
      p = p + 1
    end do
    write(*,"('INFO',6x,': reading from ',a)") trim(vpath(p))

    select case(comp(p))
#ifdef ZSTD
    case(1)
      allocate(input(csize(p)), buffer(usize(p)))

      open(newunit = io, file = trim(vpath(p)), form = 'unformatted',          &
           access = 'direct', recl = csize(p))
      read(io, rec = 1) input(1:csize(p))
      close(io)

      rsize = zstd_decompress(c_loc(buffer), size(buffer, kind=8),             &
                              c_loc(input) , size(input, kind=8))

      if (chash(p) /= xxh64(transfer(input, 1_1, csize(p)))) then
        write(error_unit,"('[',a,']: ',a)") trim(loc),                         &
                         "File '" // trim(vpath(p)) // "' seem to be corrupted!"
      end if

      if (uhash(p) /= xxh64(transfer(buffer, 1_1, usize(p)))) then
        write(error_unit,"('[',a,']: ',a)") trim(loc),                         &
                         "Decompressed data from '" // trim(vpath(p)) //       &
                         "' seem to be corrupted!"
      end if

      if (kinds(p) == 4) then
        qty = real(reshape(transfer(buffer, [ 0.0_4 ]), shape(qty)), 8)
      else if (kinds(p) == 8) then
        qty =      reshape(transfer(buffer, [ 0.0_8 ]), shape(qty))
      end if

      deallocate(input)
      deallocate(buffer)
#endif /* ZSTD */
    case default
      allocate(buffer(usize(p)))

      open(newunit = io, file = trim(vpath(p)), form = 'unformatted',          &
           access = 'direct', recl = usize(p))
      read(io, rec = 1) buffer(1:usize(p))
      close(io)

      if (uhash(p) /= xxh64(transfer(buffer, 1_1, usize(p)))) then
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "'" // trim(vpath(p)) // "' seems to be corrupted!"
      end if

      if (kinds(p) == 4) then
        qty = real(reshape(transfer(buffer, [ 0.0_4 ]), shape(qty)), 8)
      else if (kinds(p) == 8) then
        qty =      reshape(transfer(buffer, [ 0.0_8 ]), shape(qty))
      end if

      deallocate(buffer)
    end select

!-------------------------------------------------------------------------------
!
  end subroutine dataxml_read_var
!
!===============================================================================
!
! subroutine DATAXML_STORE_LAST_STATE:
! -----------------------------------
!
!   Subroutine stores the last particle states.
!
!===============================================================================
!
  subroutine dataxml_store_last_state(state, lsnap, status)

    use hash, only : xxh64

    implicit none

! subroutine arguments
!
    real(kind=8), dimension(:,:), intent(in)  :: state
    integer     , dimension(:)  , intent(in)  :: lsnap
    integer                     , intent(out) :: status

! local variables
!
    character(len=128) :: ss, dd, hh
    integer            :: nf, np, io, reclen
    integer(kind=8)    :: shash, lhash

!-------------------------------------------------------------------------------
!
    status = 0

! get the number of stored fields and particles
!
    nf = size(state, 1)
    np = size(state, 2)

! get the field hashes
!
    shash = xxh64(transfer(state, [ 0_1 ]))
    lhash = xxh64(transfer(lsnap, [ 0_1 ]))

! write the metadata
!
    open(newunit = io, file = "last_state.xml")
    write(io,"(a)") '<?xml version="1.0" ?>'
    write(io,"(a)") '<DataXML  version="1.0" byte_order="LittleEndian">'
    write(io,"(2x,a)") '<Attributes>'
    write(ss,"(i9)") nf
    dd = trim(adjustl(ss))
    write(io,"(4x,a,a,a)") '<Attribute name="number_of_fields" type="int32">', &
                           trim(adjustl(ss)), '</Attribute>'
    write(ss,"(i9)") np
    dd = trim(adjustl(dd)) // " " // trim(adjustl(ss))
    write(io,"(4x,a,a,a)") '<Attribute name="number_of_particles"' //          &
                           ' type="int32">', trim(adjustl(ss)), '</Attribute>'
    write(io,"(2x,a)") '</Attributes>'
    write(io,"(2x,a)") '<DataSets variables="last_state last_snapshot">'
    write(ss,"(i9)") 8 * nf * np
    write(hh,"(1z0.16)") shash
    write(io,"(4x,a)") '<FileName name="last_state" dtype="float64"' //        &
                       ' ndim="3" dims="' // trim(adjustl(dd)) //              &
                       '" size="' // trim(adjustl(ss)) //                      &
                       '" digest_type="xxh64" digest="' // trim(adjustl(hh)) //&
                       '">last_state.bin</FileName>'
    write(dd,"(i9)") np
    write(ss,"(i9)") 4 * np
    write(hh,"(1z0.16)") lhash
    write(io,"(4x,a)") '<FileName name="last_snapshot" dtype="int32"' //       &
                       ' ndim="1" dims="' // trim(adjustl(dd)) //              &
                       '" size="' // trim(adjustl(ss)) //                      &
                       '" digest_type="xxh64" digest="' // trim(adjustl(hh)) //&
                       '">last_snapshot.bin</FileName>'
    write(io,"(2x,a)") '</DataSets>'
    write(io,"(a)") '</DataXML>'
    close(io)

! write the particle data
!
    inquire(iolength = reclen) state(:,:)
    open(newunit = io, file = "last_state.bin", form = 'unformatted',          &
         access = 'direct', recl = reclen)
    write(io, rec = 1) state
    close(io)
    inquire(iolength = reclen) lsnap(:)
    open(newunit = io, file = "last_snapshot.bin", form = 'unformatted',       &
         access = 'direct', recl = reclen)
    write(io, rec = 1) lsnap
    close(io)

!-------------------------------------------------------------------------------
!
  end subroutine dataxml_store_last_state
!
!===============================================================================
!
! subroutine DATAXML_RESTORE_LAST_STATE:
! -------------------------------------
!
!   Subroutine restores the last particle states.
!
!===============================================================================
!
  subroutine dataxml_restore_last_state(state, lsnap, status)

    use hash           , only : xxh64
    use iso_fortran_env, only : error_unit

    implicit none

! subroutine arguments
!
    real(kind=8), dimension(:,:), intent(inout) :: state
    integer     , dimension(:)  , intent(inout) :: lsnap
    integer                     , intent(out)   :: status

! local variables
!
    logical            :: info, spres = .false., lpres = .false.
    character(len=256) :: path
    integer            :: io, nf, np, p, dtype, ctype
    integer(kind=8)    :: sbytes, lbytes, shash, lhash, tmp1, tmp2

    character(len=:), allocatable :: buffer

! local parameters
!
    character(len=*), parameter :: loc = 'DATAXML::dataxml_restore_last_state()'

!-------------------------------------------------------------------------------
!
    status = 0

! check if the file exists and get its size
!
    inquire(file = "last_state.xml", exist = info, size = p)
    if (.not. info) then
      write(*,'("ERROR   : file ",a," does not exist!")') '"last_state.xml"'
      status = 1
      return
    end if

! allocate a bufer and store the content of the metadata file in it
!
    allocate(character(len=p) :: buffer)
    open(newunit = io, file = "last_state.xml", form = 'unformatted',          &
                       access = 'direct', recl = p)
    read(io, rec = 1) buffer(:)
    close(io)

! parse the buffer to get the required information
!
    call get_attribute(buffer, 'number_of_fields'   , nf)
    call get_attribute(buffer, 'number_of_particles', np)

! parse filename information
!
    call get_filename(buffer, 'last_state', path, dtype, ctype,                &
                                                  sbytes, tmp1, shash, tmp2)
    call get_filename(buffer, 'last_snapshot', path, dtype, ctype,             &
                                                  lbytes, tmp1, lhash, tmp2)

! deallocate the buffer
!
    deallocate(buffer)

! check if the dimensions are correct
!
    if (nf /= size(state, 1) .or. np /= size(state, 2) .or.                    &
                                  np /= size(lsnap, 1)) then

      write(error_unit,"('[',a,']: ',a)") trim(loc), "Dimensions do not match!"

      status = 2

      return

    end if

! read binary data
!
    open(newunit = io, file = "last_state.bin", form = 'unformatted',          &
         access = 'direct', recl = sbytes)
    read(io, rec = 1) state
    close(io)

    open(newunit = io, file = "last_snapshot.bin", form = 'unformatted',       &
         access = 'direct', recl = lbytes)
    read(io, rec = 1) lsnap
    close(io)

! check if the hashes of the read data are correct
!
    if (spres) then
      if (shash /= xxh64(transfer(state, [ 0_1 ]))) then
        write(error_unit,"('[',a,']: ',a)") trim(loc),                         &
                         "File 'last_state.bin' seems to be corrupted!"
      end if
    end if
    if (lpres) then
      if (lhash /= xxh64(transfer(lsnap, [ 0_1 ]))) then
        write(error_unit,"('[',a,']: ',a)") trim(loc),                         &
                         "File 'last_snapshot.bin' seems to be corrupted!"
      end if
    end if

    return

!-------------------------------------------------------------------------------
!
  end subroutine dataxml_restore_last_state
!
!===============================================================================
!
! subroutine GET_ATTRIBUTE_INTEGER:
! --------------------------------
!
!   Subroutine returns the value of an integer attribute specified by its name.
!
!===============================================================================
!
  subroutine get_attribute_integer(buffer, aname, avalue)

    implicit none

! subroutine arguments
!
    character(len=*), intent(in)  :: buffer, aname
    integer         , intent(out) :: avalue

! local variables
!
    integer :: ib, ie, lb, le

!-------------------------------------------------------------------------------
!
    ib = index(buffer, '<Attribute name="' // trim(aname) // '"')
    if (ib > 0) then
      ie = index(buffer(ib:),'</Attribute>') + ib + len('</Attribute>') - 2
      lb = index(buffer(ib:ie), '>') + ib
      le = index(buffer(ib:ie), '<', .true.) + ib - 2
      read(buffer(lb:le), fmt=*) avalue
    end if

!-------------------------------------------------------------------------------
!
  end subroutine get_attribute_integer
!
!===============================================================================
!
! subroutine GET_ATTRIBUTE_DOUBLE:
! -------------------------------
!
!   Subroutine returns the value of a double precision attribute specified by
!   its name.
!
!===============================================================================
!
  subroutine get_attribute_double(buffer, aname, avalue)

    implicit none

! subroutine arguments
!
    character(len=*), intent(in)  :: buffer, aname
    real(kind=8)    , intent(out) :: avalue

! local variables
!
    integer :: ib, ie, lb, le

!-------------------------------------------------------------------------------
!
    ib = index(buffer, '<Attribute name="' // trim(aname) // '"')
    if (ib > 0) then
      ie = index(buffer(ib:),'</Attribute>') + ib + len('</Attribute>') - 2
      lb = index(buffer(ib:ie), '>') + ib
      le = index(buffer(ib:ie), '<', .true.) + ib - 2
      read(buffer(lb:le), fmt=*) avalue
    end if

!-------------------------------------------------------------------------------
!
  end subroutine get_attribute_double
!
!===============================================================================
!
! subroutine GET_FILENAME:
! -----------------------
!
!   Subroutine returns the name of the file name with its size, type of data
!   and other attributes.
!
!===============================================================================
!
  subroutine get_filename(buffer, aname, avalue, atype, acomp,                 &
                                                 usize, csize, uhash, chash)

    implicit none

! subroutine arguments
!
    character(len=*), intent(in)    :: buffer, aname
    character(len=*), intent(inout) :: avalue
    integer         , intent(out)   :: atype, acomp
    integer(kind=8) , intent(out)   :: usize, csize, uhash, chash

! local variables
!
    integer :: ib, ie, lb, le

!-------------------------------------------------------------------------------
!
    ib = index(buffer, '<FileName name="' // trim(aname) // '"')
    if (ib > 0) then
      ie = index(buffer(ib:),'</FileName>') + ib + len('</FileName>') - 2
      lb = index(buffer(ib:ie), '>') + ib
      le = index(buffer(ib:ie), '<', .true.) + ib - 2
      read(buffer(lb:le), fmt=*) avalue
      lb = index(buffer(ib:ie), 'dtype="') + len('dtype="') + ib - 1
      le = index(buffer(lb:ie), '"') + lb - 2
      select case(buffer(lb:le))
      case('float64', 'int64')
        atype = 8
      case('float32', 'int32')
        atype = 4
      case default
        atype = 0
      end select
      lb = index(buffer(ib:ie), 'compression="') + len('compression="') + ib - 1
      le = index(buffer(lb:ie), '"') + lb - 2
      select case(buffer(lb:le))
      case('zstd')
        acomp = 1
      case('lz4')
        acomp = 2
      case default
        acomp = 0
      end select
      lb = index(buffer(ib:ie), 'size="') + len('size="') + ib - 1
      le = index(buffer(lb:ie), '"') + lb - 2
      read(buffer(lb:le), fmt=*) usize
      if (index(buffer(ib:ie), 'compressed_size="') > 0) then
        lb = index(buffer(ib:ie), 'compressed_size="') +                       &
                            len('compressed_size="') + ib - 1
        le = index(buffer(lb:ie), '"') + lb - 2
        read(buffer(lb:le), fmt=*) csize
      end if
      if (index(buffer(ib:ie), 'digest="') > 0) then
        lb = index(buffer(ib:ie), 'digest="') + len('digest="') + ib - 1
        le = index(buffer(lb:ie), '"') + lb - 2
        read(buffer(lb:le), fmt='(1z16)') uhash
      end if
      if (index(buffer(ib:ie), 'compressed_digest="') > 0) then
        lb = index(buffer(ib:ie), 'compressed_digest="') +                     &
                              len('compressed_digest="') + ib - 1
        le = index(buffer(lb:ie), '"') + lb - 2
        read(buffer(lb:le), fmt='(1z16)') chash
      end if
    end if

!-------------------------------------------------------------------------------
!
  end subroutine get_filename

!===============================================================================
!
end module dataxml
