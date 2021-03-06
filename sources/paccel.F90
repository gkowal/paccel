!!******************************************************************************
!!
!! Program: PAccel - particle acceleration
!!
!! Copyright (C) 2008-2021 Grzegorz Kowal <grzegorz@gkowal.info>
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
program paccel

! import required modules
!
  use coordinates   , only : initialize_coordinates
  use fields        , only : initialize_fields, finalize_fields, read_fields
  use interpolations, only : initialize_interpolations
  use parameters    , only : read_parameters, get_parameter
  use particles     , only : initialize_particles, generate_particle,          &
                             integrate_trajectory_rk4,                         &
                             integrate_trajectory_dop853,                      &
                             integrate_trajectory_si4,                         &
                             integrate_trajectory_si6,                         &
                             integrate_trajectory_si8,                         &
                             integrate_trajectory_si10,                        &
                             integrate_trajectory_si12

  implicit none

  character(len = 32), save :: method = 'rk4' ! the integration method

  integer :: status
  real    :: timer
!
!-------------------------------------------------------------------------------
!
! print header
!
  write (*,'(a)') '------------------------------------------------------------------------------'
  write (*,'(a)') '===        PAccel algorithm started         =================================='
  write (*,'(a)') '===  Copyright (C) 2008-2021 Grzegorz Kowal =================================='
  write (*,*)
  write( *, "('TASK      : ',a)" ) "integrating the trajectory of a charged particle"

  call read_parameters(.true., status)

#ifndef TEST
  call initialize_fields(.true., status)

  call initialize_coordinates(.true., status)

  call initialize_interpolations(.true., status)
#endif /* !TEST */

  call initialize_particles(.true., status)

#ifndef TEST
  call read_fields(.true., status)
#endif /* !TEST */

  call generate_particle(.true., status)

! get the integration method
!
  call get_parameter('method', method)

! take the time of calculations only
!
  timer = secnds(0.0)

! integrate particle trajectories
!
  select case(trim(method))
  case('rk4')
    write( *, "('INFO      : ',a)" ) "integrating the particle trajectory (RK4 method)"
    call integrate_trajectory_rk4()
  case('rk8', 'dop853')
    write( *, "('INFO      : ',a)" ) "integrating the particle trajectory (DOP853 method)"
    call integrate_trajectory_dop853()
  case('si4')
    write( *, "('INFO      : ',a)" ) "integrating the particle trajectory (SI4 method)"
    call integrate_trajectory_si4()
  case('si6')
    write( *, "('INFO      : ',a)" ) "integrating the particle trajectory (SI6 method)"
    call integrate_trajectory_si6()
  case('si8')
    write( *, "('INFO      : ',a)" ) "integrating the particle trajectory (SI8 method)"
    call integrate_trajectory_si8()
  case('si10')
    write( *, "('INFO      : ',a)" ) "integrating the particle trajectory (SI10 method)"
    call integrate_trajectory_si10()
  case('si12')
    write( *, "('INFO      : ',a)" ) "integrating the particle trajectory (SI12 method)"
    call integrate_trajectory_si12()
  end select

! display performance information
!
  timer = secnds(timer)
  write( *, "('COMPUTED  : ',a,1pe12.5,a)" ) 'computing done in ', timer, ' seconds'

! deallocate field variables
!
  call finalize_fields(.true.)

!===============================================================================
!
end program paccel
