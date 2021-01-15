!!******************************************************************************
!!
!!  This file is part of the GACCEL source code, a program to integrate
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
  integer           , dimension(9), save :: usize = 0
  integer           , dimension(9), save :: csize = 0

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

    character(len=256) :: fname, line, str
    logical            :: status
    integer            :: p, l, io
    integer(kind=8)    :: s

    character(len=*), parameter :: loc = 'DATAXML::dataxml_init()'

!-------------------------------------------------------------------------------
!
    call get_parameter('idir', idir)
    l = len_trim(idir)
    if (idir(l:l) /= '/') then
      idir = trim(idir) // '/'
    end if

    write(fname,"(a,'datasets.xml')") trim(idir)

    inquire(file = fname, exist = status)
    if (.not. status) then
      write(*,'("ERROR   : file ",a," does not exist!")') trim(fname)
      stop
    end if

    p = 1
    open(newunit = io, file = fname, err = 30)
10  read(io, fmt = "(a)", end = 20) line
    if ((len_trim(line) == 0)                                                  &
                           .or. index(trim(adjustl(line)), '#') == 1) go to 10
    if (index(line, 'ndims') > 0) then
      l = index(line, '</Attribute>') - 1
      read(line(index(line(1:l), '>')+1:l), fmt = *) ndims
    end if
    if (index(line, 'dim(1)') > 0) then
      l = index(line, '</Attribute>') - 1
      read(line(index(line(1:l), '>')+1:l), fmt = *) dims(1)
    end if
    if (index(line, 'dim(2)') > 0) then
      l = index(line, '</Attribute>') - 1
      read(line(index(line(1:l), '>')+1:l), fmt = *) dims(2)
    end if
    if (index(line, 'dim(3)') > 0) then
      l = index(line, '</Attribute>') - 1
      read(line(index(line(1:l), '>')+1:l), fmt = *) dims(3)
    end if
    if (index(line, 'xmin') > 0) then
      l = index(line, '</Attribute>') - 1
      read(line(index(line(1:l), '>')+1:l), fmt = *) bnds(1,1)
    end if
    if (index(line, 'xmax') > 0) then
      l = index(line, '</Attribute>') - 1
      read(line(index(line(1:l), '>')+1:l), fmt = *) bnds(2,1)
    end if
    if (index(line, 'ymin') > 0) then
      l = index(line, '</Attribute>') - 1
      read(line(index(line(1:l), '>')+1:l), fmt = *) bnds(1,2)
    end if
    if (index(line, 'ymax') > 0) then
      l = index(line, '</Attribute>') - 1
      read(line(index(line(1:l), '>')+1:l), fmt = *) bnds(2,2)
    end if
    if (index(line, 'zmin') > 0) then
      l = index(line, '</Attribute>') - 1
      read(line(index(line(1:l), '>')+1:l), fmt = *) bnds(1,3)
    end if
    if (index(line, 'zmax') > 0) then
      l = index(line, '</Attribute>') - 1
      read(line(index(line(1:l), '>')+1:l), fmt = *) bnds(2,3)
    end if
    if (index(line, 'FileName') > 0) then
      l = index(line, '</FileName>') - 1
      read(line(index(line(1:l), '>')+1:l), fmt = *) str
      vars(p)  = str(1:4)
      vpath(p) = trim(idir) // trim(adjustl(str))
      l = index(line, 'dtype="') + 7
      str = trim(adjustl(line(l:)))
      l = index(str, '"') - 1
      if (str(:l) == 'float32') kinds(p) = 4
      if (str(:l) == 'float64') kinds(p) = 8
      l = index(line, 'size="') + 6
      str = trim(adjustl(line(l:)))
      l = index(str, '"') - 1
      read(str(:l), fmt = *) usize(p)
      if (index(line, 'compression') > 0) then
        if (index(line, 'compression="zstd"') > 0) comp(p) = 1
        if (index(line, 'compression="lz4"')  > 0) comp(p) = 2
        l = index(line, 'compressed_size="') + 17
        str = trim(adjustl(line(l:)))
        l = index(str, '"') - 1
        read(str(:l), fmt = *) csize(p)
      end if
      p = p + 1
    end if
    go to 10
20  close(io)

    return

30  write(error_unit,"('[',a,']: ',a)") trim(loc)                              &
                    , "Cannot open the parameter file '" // trim(fname) // "'!"

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

    real(kind=4), intent(out) :: xmin, xmax, ymin, ymax, zmin, zmax

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

    use iso_c_binding, only: c_int, c_loc

    implicit none

    character(len=*)      , intent(in)    :: var
    real, dimension(:,:,:), intent(inout) :: qty

    integer         :: io, p
    integer(kind=8), target :: rsize

    integer(kind=1), dimension(:), allocatable, target :: buffer, input

!-------------------------------------------------------------------------------
!
    p = 1
    do while(vars(p) /= var .and. p <= size(vars))
      p = p + 1
    end do
    write(*,"(a,a)") 'INFO      : reading from ', trim(vpath(p))

    if (comp(p) == 0) then
      allocate(buffer(usize(p)))

      open(newunit = io, file = trim(vpath(p)), form = 'unformatted',          &
           access = 'direct', recl = usize(p))
      read(io, rec = 1) buffer(1:usize(p))
      close(io)

      if (kinds(p) == 4) then
        qty = reshape(transfer(buffer, [ 0.0_4 ]), shape(qty))
      else if (kinds(p) == 8) then
        qty = reshape(transfer(buffer, [ 0.0_8 ]), shape(qty))
      end if

      deallocate(buffer)
#ifdef ZSTD
    else if (comp(p) == 1) then
      allocate(input(csize(p)), buffer(usize(p)))

      open(newunit = io, file = trim(vpath(p)), form = 'unformatted',          &
           access = 'direct', recl = csize(p))
      read(io, rec = 1) input(1:csize(p))
      close(io)

      rsize = zstd_decompress(c_loc(buffer), size(buffer, kind=8),             &
                              c_loc(input) , size(input, kind=8))

      if (kinds(p) == 4) then
        qty = reshape(transfer(buffer, [ 0.0_4 ]), shape(qty))
      else if (kinds(p) == 8) then
        qty = reshape(transfer(buffer, [ 0.0_8 ]), shape(qty))
      end if

      deallocate(input)
      deallocate(buffer)
#endif /* ZSTD */
    end if

  end subroutine dataxml_read_var

!===============================================================================
!
end module dataxml
