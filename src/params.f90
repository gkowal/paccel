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
  character(len = 128), save :: idir     = "./"        ! input directory
  character(len = 128), save :: odir     = "./"        ! output directory
  character(len =   4), save :: fformat  = 'fits'      ! file format:
                                                       ! 'fits' - FITS
                                                       ! 'hdf4' - HDF4
                                                       ! 'hdf5' - HDF5
  character(len =   1), save :: ftype    = 'r'         ! file type: 'r', 'p', 'f'
  integer             , save :: fnumber  = 0           ! file number
  character(len =   1), save :: ptype    = 'p'         ! particle type:
                                                       ! 'p' - proton
                                                       ! 'e' - electron
  character(len =   1), save :: tunit    = 's'         ! time units:
                                                       ! 's' - second
                                                       ! 'm' - minute
                                                       ! 'h' - hour
                                                       ! 'd' - day
                                                       ! 'w' - week
                                                       ! 'y' - year
  real                , save :: tmulti   = 1.0         ! time unit count
  real                , save :: c        = 1.0         ! the speed of light in Va
  real                , save :: dens     = 1.0         ! density [1 u/cm^3]
  real                , save :: ueta     = 0.0         ! uniform resistivity coeff
  real                , save :: aeta     = 0.0         ! anomalous resistivity coeff
  real                , save :: jcrit    = 1.0e3       ! critical current density
  integer             , save :: nsteps   = 1           ! number of steps
  real                , save :: xc       = 0.0         ! initial position
  real                , save :: yc       = 0.0
  real                , save :: zc       = 0.0
  real                , save :: vper     = 1.0         ! initial perpendicular speed
  character(len =   1), save :: approx   = 'n'         ! guiding centre approximation
  character(len =   1), save :: periodic = 'y'         ! periodic box or not
  character(len =   1), save :: efield   = 'y'         ! take electric field into account
  character(len =   1), save :: current  = 'y'         ! take current density into account
  real                , save :: cfl      = 0.5         ! cfl condition
  real                , save :: dtout    = 1.0         ! interval between data writing
  real                , save :: tmax     = 1.0         ! maximum time for integration
  real                , save :: ethres   = 1.0         ! energy threshold
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
      case ('fformat')
        l = len_trim(value)
        write(fformat, "(a)" ) value(2:l-1)
      case ('ftype')
        l = len_trim(value)
        write(ftype  , "(a)" ) value(2:l-1)
      case ('fnumber')
        read (value  , "(i6)") fnumber
      case ('ptype')
        l = len_trim(value)
        write(ptype  , "(a)" ) value(2:l-1)
      case ('tunit')
        l = len_trim(value)
        write(tunit  , "(a)" ) value(2:l-1)
      case ('tmulti')
        read (value  , *) tmulti
      case ('c')
        read (value  , *) c
      case ('dens')
        read (value  , *) dens
      case ('ueta')
        read (value  , *) ueta
      case ('aeta')
        read (value  , *) aeta
      case ('jcrit')
        read (value  , *) jcrit
      case ('nsteps')
        read (value  , "(i9)") nsteps
      case ('xc')
        read (value  , *) xc
      case ('yc')
        read (value  , *) yc
      case ('zc')
        read (value  , *) zc
      case ('vper')
        read (value  , *) vper
      case ('approximation')
        l = len_trim(value)
        write(approx  , "(a)" ) value(2:l-1)
      case ('periodic')
        l = len_trim(value)
        write(periodic, "(a)" ) value(2:l-1)
      case ('efield')
        l = len_trim(value)
        write(efield  , "(a)" ) value(2:l-1)
      case ('current')
        l = len_trim(value)
        write(current , "(a)" ) value(2:l-1)
      case ('cfl')
        read (value  , *) cfl
      case ('dtout')
        read (value  , *) dtout
      case ('tmax')
        read (value  , *) tmax
      case ('ethres')
        read (value  , *) ethres
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
