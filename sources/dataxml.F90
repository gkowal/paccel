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
! subroutine PARSE_ATTRIBUTE:
! --------------------------
!
!   Subroutine parses the buffer for specific attribute and returns
!   the indices to its value.
!
!===============================================================================
!
  subroutine parse_attribute(buffer, aname, ibeg, iend)

    use iso_fortran_env, only : error_unit

    implicit none

    character(len=*), intent(in)  :: buffer, aname
    integer         , intent(out) :: ibeg, iend

    logical :: found
    integer :: il, iu, ib, ie

    character(len=*), parameter :: loc = 'DATAXML::parse_attribute()'

!-------------------------------------------------------------------------------
!
    ibeg = 0
    iend = 0
    ib   = 0
    ie   = 0

    il = index(buffer,'<Attributes>')
    iu = index(buffer(il:),'</Attributes>')
    if (il > 0 .and. iu > 0) then
      iu = iu - len('</Attributes>') + il + len('</Attributes>') - 2
      il = il + len('<Attributes>')

      found = .false.
      do while(.not. found .and. il < iu)
        ib = index(buffer(il:iu), '<Attribute')
        if (ib > 0) then
          ib = il + ib - 1
          ie = index(buffer(ib:iu), '</Attribute>')                            &
                              + len('</Attribute>') + ib - 2
          found = index(buffer(ib:ie), 'name="' // trim(aname) // '"') > 0
          if (.not. found) il = ie
        else
          il = iu
          write(error_unit,"('[',a,']: ',a)") trim(loc),                       &
                           "Group <Attributes> has no element <Attribute>!"
          stop
        end if
      end do
      if (found) then
        ibeg = index(buffer(ib:ie), '>') + ib
        iend = index(buffer(ib:ie), '<', .true.) + ib - 2
        return
      else
        write(error_unit,"('[',a,']: ',a)") trim(loc),                         &
                         'Attribute "' // trim(aname) // '" could not be found.'
        stop
      end if
    else
      write(error_unit,"('[',a,']: ',a)") trim(loc),                           &
                       "Missing or wrong group <Attributes>!"
      stop
    end if

!-------------------------------------------------------------------------------
!
  end subroutine parse_attribute
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

    character(len=*), intent(in)  :: buffer, aname
    integer         , intent(out) :: avalue

    integer :: ibeg, iend

!-------------------------------------------------------------------------------
!
    call parse_attribute(buffer, aname, ibeg, iend)
    if (ibeg > 0 .and. iend > 0) read(buffer(ibeg:iend), fmt=*) avalue

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

    character(len=*), intent(in)  :: buffer, aname
    real(kind=8)    , intent(out) :: avalue

    integer :: ibeg, iend

!-------------------------------------------------------------------------------
!
    call parse_attribute(buffer, aname, ibeg, iend)
    if (ibeg > 0 .and. iend > 0) read(buffer(ibeg:iend), fmt=*) avalue

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

    use iso_fortran_env, only : error_unit

    implicit none

    character(len=*), intent(in)    :: buffer, aname
    character(len=*), intent(inout) :: avalue
    integer         , intent(out)   :: atype, acomp
    integer(kind=8) , intent(out)   :: usize, csize, uhash, chash

    logical :: found
    integer :: il, iu, ib, ie, lb, le

    character(len=*), parameter :: loc = 'DATAXML::get_filename()'

!-------------------------------------------------------------------------------
!
    il = index(buffer,'<DataSets')
    iu = index(buffer(il:),'</DataSets>')
    if (il > 0 .and. iu > 0) then
      iu = iu + il + len('<DataSets') - len('</DataSets>')
      il = il + index(buffer(il:),'>')

      found = .false.
      do while(.not. found .and. il < iu)
        ib = index(buffer(il:iu), '<FileName')
        if (ib > 0) then
          ib = il + ib - 1
          ie = index(buffer(ib:iu), '</FileName>')                             &
                              + len('</FileName>') + ib - 2
          found = index(buffer(ib:ie), 'name="' // trim(aname) // '"') > 0
          if (.not. found) il = ie
        else
          il = iu
          write(error_unit,"('[',a,']: ',a)") trim(loc),                       &
                           "Group <DataSets> has no element <FileName>!"
          stop
        end if
      end do
      if (found) then

! parse the value
!
        lb = index(buffer(ib:ie), '>') + ib
        le = index(buffer(ib:ie), '<', .true.) + ib - 2
        write(avalue,*) trim(adjustl(buffer(lb:le)))

! parse the datatype
!
        lb = index(buffer(ib:ie), 'dtype="')
        if (lb > 0) then
          lb = lb + len('dtype="') + ib - 1
          le = index(buffer(lb:ie), '"') + lb - 2
          select case(buffer(lb:le))
          case('float64', 'int64')
            atype = 8
          case('float32', 'int32')
            atype = 4
          case default
            atype = 0
          end select
        else
          write(error_unit,"('[',a,']: ',a)") trim(loc),                       &
                           "Element does not have 'dtype' attribute!"
          stop
        end if

! parse the uncompressed size
!
        lb = index(buffer(ib:ie), ' size="')
        if (lb > 0) then
          lb = lb + len(' size="') + ib - 1
          le = index(buffer(lb:ie), '"') + lb - 2
          read(buffer(lb:le), fmt=*) usize
        else
          write(error_unit,"('[',a,']: ',a)") trim(loc),                       &
                           "Element does not have 'size' attribute!"
          stop
        end if

! parse the uncompressed data digest
!
        if (index(buffer(ib:ie), ' digest="') > 0) then
          lb = index(buffer(ib:ie), ' digest="') + len(' digest="') + ib - 1
          le = index(buffer(lb:ie), '"') + lb - 2
          read(buffer(lb:le), fmt='(1z16)') uhash
        end if

! parse the compression
!
        if (index(buffer(ib:ie), 'compression=') > 0) then
          lb = index(buffer(ib:ie), 'compression="')                           &
                              + len('compression="') + ib - 1
          le = index(buffer(lb:ie), '"') + lb - 2
          select case(buffer(lb:le))
          case('zstd')
            acomp = 1
          case('lz4')
            acomp = 2
          case default
            acomp = 0
          end select

! parse the compressed size
!
          if (index(buffer(ib:ie), 'compressed_size="') > 0) then
            lb = index(buffer(ib:ie), 'compressed_size="') +                   &
                                  len('compressed_size="') + ib - 1
            le = index(buffer(lb:ie), '"') + lb - 2
            read(buffer(lb:le), fmt=*) csize
          else
            write(error_unit,"('[',a,']: ',a)") trim(loc),                     &
                  "Compression used but no 'compressed_size' attribute present!"
            stop
          end if

! parse the compressed data digest
!
          if (index(buffer(ib:ie), 'compressed_digest="') > 0) then
            lb = index(buffer(ib:ie), 'compressed_digest="') +                 &
                                  len('compressed_digest="') + ib - 1
            le = index(buffer(lb:ie), '"') + lb - 2
            read(buffer(lb:le), fmt='(1z16)') chash
          end if
        else
          acomp = 0
        end if
      else
        write(error_unit,"('[',a,']: ',a)") trim(loc),                         &
            'Element <FileName name="' // trim(aname) // '" could not be found.'
        stop
      end if
    else
      write(error_unit,"('[',a,']: ',a)") trim(loc),                           &
                       "Missing or wrong group <DataSets>!"
      stop
    end if

!-------------------------------------------------------------------------------
!
  end subroutine get_filename

!===============================================================================
!
end module dataxml
