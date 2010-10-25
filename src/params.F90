!!******************************************************************************
!!
!! module: params - subroutines to read parameters file.
!!
!! Copyright (C) 2007-2010 Grzegorz Kowal <grzegorz@gkowal.info>
!!
!!******************************************************************************
!!
!!  This file is part of PAccel.
!!
!!  PAccel is free software; you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation; either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  PAccel is distributed in the hope that it will be useful,
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
  character(len = 128), save :: idir     = "./"          ! input directory
  character(len = 128), save :: odir     = "./"          ! output directory
  character(len =   4), save :: fformat  = 'fits'        ! file format:
                                                         ! 'fits' - FITS
                                                         ! 'hdf4' - HDF4
                                                         ! 'hdf5' - HDF5
  character(len =   1), save :: ftype    = 'r'           ! file type: 'r', 'p', 'f'
  integer             , save :: fnumber  = 0             ! file number
  character(len =   1), save :: ptype    = 'p'           ! particle type:
                                                         ! 'p' - proton
                                                         ! 'e' - electron
  character(len =   1), save :: tunit    = 's'           ! time units:
                                                         ! 's' - second
                                                         ! 'm' - minute
                                                         ! 'h' - hour
                                                         ! 'd' - day
                                                         ! 'w' - week
                                                         ! 'y' - year
  real                , save :: tmulti   = 1.0           ! time unit count
  real                , save :: c        = 1.0           ! the speed of light in Va
  real                , save :: dens     = 1.0           ! density [1 u/cm^3]
  real                , save :: ueta     = 0.0           ! uniform resistivity coeff
  real                , save :: aeta     = 0.0           ! anomalous resistivity coeff
  real                , save :: jcrit    = 1.0e3         ! critical current density
  real                , save :: xc       = 0.0           ! initial position
  real                , save :: yc       = 0.0
  real                , save :: zc       = 0.0
  real                , save :: vpar     = 0.0           ! initial parallel speed [in c]
  real                , save :: vper     = 0.1           ! initial perpendicular speed [in c]
  real                , save :: rho      = 0.5           ! safety coefficient
  real                , save :: tolerance = 1.0e-4       ! integration tolerance
  real                , save :: dtini     = 1.0e-8       ! the initial time step
  real                , save :: dtmax     = 1.0          ! maximum allowed step size

! the parameters controlling the data output
!
  character(len =   1), save :: output    = 'i'          ! the type of output:
                                                         ! 'i' - by the iteration number
                                                         ! 'l' - by the logarithmic time
  integer             , save :: ndumps    = 1000         ! number of steps between subsequent dumps if the output is 'i'
                                                         ! or number of dumps per time decade if the output is 'l'
  real                , save :: tmin      = 1.0e-3       ! minimum time of writing data
  real                , save :: tmax      = 1.0          ! maximum time for integration

#ifdef TEST
! the parameters controlling the test cases
!
  real                , save :: omega    = 1.0           ! omega
  real                , save :: bini     = 1.0           ! mean magnetic field
  real                , save :: bamp     = 0.1           ! amplitude of the magnetic field fluctuations
  real                , save :: vamp     = 0.1           ! amplitude of the velocity field fluctuations
  real                , save :: freq     = 1.0           ! frequency of the field perturbation
  real                , save :: epar     = 0.0           ! constant electric field along parallel direction
#endif /* TEST */
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
    real               :: vv
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
      case ('xc')
        read (value  , *) xc
      case ('yc')
        read (value  , *) yc
      case ('zc')
        read (value  , *) zc
      case ('vpar')
        read (value  , *) vpar
      case ('vper')
        read (value  , *) vper
      case ('rho')
        read (value  , *) rho
      case ('tolerance')
        read (value  , *) tolerance
      case ('dtini')
        read (value  , *) dtini
      case ('dtmax')
        read (value  , *) dtmax

      case ('output')
        l = len_trim(value)
        write(output , "(a)" ) value(2:l-1)
      case ('ndumps')
        read (value  , "(i9)") ndumps
      case ('tmin')
        read (value  , *     ) tmin
      case ('tmax')
        read (value  , *     ) tmax

#ifdef TEST
      case ('omega')
        read (value  , *) omega
      case ('bini')
        read (value  , *) bini
      case ('bamp')
        read (value  , *) bamp
      case ('vamp')
        read (value  , *) vamp
      case ('freq')
        read (value  , *) freq
      case ('epar')
        read (value  , *) epar
#endif /* TEST */
      case default
    end select

    go to 10

20  close(1)

! check input parameters
!
    vv = sqrt(vpar**2 + vper**2)
    if (vv .ge. 1.0) then
      write( *, "('ERROR     : ',a)" ) "absolute speed of the particle is larger than c!"
      write( *, "('ERROR     : ',a,1pe15.8)" ) "|v| = ", vv
      stop
    end if
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
