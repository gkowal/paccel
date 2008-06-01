!!******************************************************************************
!!
!! module: params - subroutines to read parameters file.
!!
!! Copyright (C) 2007 Grzegorz Kowal <kowal@astro.wisc.edu>
!!
!!******************************************************************************
!!
!!  This file is part of Spectrum.
!!
!!  Godunov-MHD is free software; you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation; either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  Foobar is distributed in the hope that it will be useful,
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
module params

  implicit none

! input parameters
!
  character(len = 128), save :: idir    = "./"         ! input directory
  character(len = 128), save :: odir    = "./"         ! output directory
  character(len =   4), save :: field   = "velo"       ! spectrum for which field
  character(len =   3), save :: frame   = "glo"        ! reference frame 'glo' or 'loc'
  character(len =   1), save :: stype   = "A"          ! kind of spectrum:
                                                       ! 'F' - Fourier
                                                       ! 'A' - Fourier (averaged over integer ks)
                                                       ! 'W' - Wavelet
  character(len =   4), save :: fformat = 'fits'       ! file format:
                                                       ! 'fits' - FITS
                                                       ! 'hdf4' - HDF4
                                                       ! 'hdf5' - HDF5
  character(len =   1), save :: ftype   = 'r'          ! file type: 'r', 'p', 'f'
  character(len =   1), save :: sdir    = 'r'          ! direction along which spectrum is taken
                                                       ! 'R' - radial
                                                       ! 'X', 'Y', 'Z' - x, y, z direction
  integer             , save :: fnumber = 0            ! file number
  integer             , save :: nsteps  = 1            ! number of steps
  real                , save :: xc      = 0.0          ! initial position
  real                , save :: yc      = 0.0
  real                , save :: zc      = 0.0
  real                , save :: dt      = 1.0          ! time steps of integration
  real                , save :: ueta    = 0.0          ! uniform resistivity coeff
  real                , save :: aeta    = 0.0          ! anomalous resistivity coeff
  real                , save :: jcrit   = 1.0e3        ! critical current density
  real                , save :: qom     = 1.0          ! e / m
  real                , save :: c       = 1.0          ! the speed of light in VA

! common variables
!
  character(len =   4), save :: vars(3)
!
!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! read_config: subroutine to read config file and to fill proper parameters
!
!===============================================================================
!
  subroutine read_params
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    character(len =   *), parameter :: configfile = './params.in' ! config file with parameters
    character(len=255) :: line, name, value
    integer            :: l, i, ios
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    open(unit=1, file=configfile, err=100)

10  read(unit=1, fmt="(a)", end=20, err=200 ) line

    call parse_line(line, name, value)

    select case(name)
      case ('idir')
        l = len_trim(value)
        write(idir   , "(a)" ) value(2:l-1)
      case ('odir')
        l = len_trim(value)
        write(odir   , "(a)" ) value(2:l-1)
      case ('stype')
        l = len_trim(value)
        write(stype  , "(a)" ) value(2:l-1)
      case ('frame')
        l = len_trim(value)
        write(frame  , "(a)" ) value(2:l-1)
      case ('field')
        l = len_trim(value)
        write(field  , "(a)" ) value(2:l-1)
      case ('fformat')
        l = len_trim(value)
        write(fformat, "(a)" ) value(2:l-1)
      case ('ftype')
        l = len_trim(value)
        write(ftype  , "(a)" ) value(2:l-1)
      case ('sdir')
        l = len_trim(value)
        write(sdir   , "(a)" ) value(2:l-1)
      case ('fnumber')
        read (value  , "(i6)") fnumber
      case ('nsteps')
        read (value  , "(i9)") nsteps
      case ('xc')
        read (value  , *) xc
      case ('yc')
        read (value  , *) yc
      case ('zc')
        read (value  , *) zc
      case ('dt')
        read (value  , *) dt
      case ('ueta')
        read (value  , *) ueta
      case ('aeta')
        read (value  , *) aeta
      case ('jcrit')
        read (value  , *) jcrit
      case ('qom')
        read (value  , *) qom
      case ('c')
        read (value  , *) c
      case default
    end select

    go to 10

20  close(1)
!
    return
!
100 print *, 'Error opening file ', configfile
    stop
200 print *, 'Error reading file ', configfile
    stop
!
    return
!
  end subroutine read_params
!
!===============================================================================
!
! parse_line: subroutine to parse line into name and value of parameter
!
!===============================================================================
!
  subroutine parse_line(line, name, value)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    character(len=*), intent(in)  :: line
    character(len=*), intent(out) :: name, value

    integer :: l, i
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    l = len_trim(line)

    i = index( line, '=' )

    name  = trim(adjustl(line(1:i-1)))
    value = trim(adjustl(line(i+1:l)))
!
  end subroutine parse_line

end module params
