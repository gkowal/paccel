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

  use fields    , only : init_fields, finit_fields
  use parameters, only : read_parameters, get_parameter
  use particles , only : init_particle, finit_particle                         &
                       , integrate_trajectory_rk4                              &
                       , integrate_trajectory_si4, integrate_trajectory_si4v   &
                       , integrate_trajectory_si6, integrate_trajectory_si6v   &
                       , integrate_trajectory_si8, integrate_trajectory_si8v

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
  write( *, "('INFO      : ',a)" ) "reading parameters"

! read parameters
!
  call read_parameters(.true., status)

#ifndef TEST
! initialize field variables
!
  write( *, "('INFO      : ',a)" ) "initializing and reading the field components"
  call init_fields()
#endif /* !TEST */

! initiate particle
!
  write( *, "('INFO      : ',a)" ) "initializing the particle positions and velocities"
  call init_particle()

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
  case('si4')
    write( *, "('INFO      : ',a)" ) "integrating the particle trajectory (SI4 method)"
    call integrate_trajectory_si4()
  case('si4v')
    write( *, "('INFO      : ',a)" ) "integrating the particle trajectory (SI4v method)"
    call integrate_trajectory_si4v()
  case('si6')
    write( *, "('INFO      : ',a)" ) "integrating the particle trajectory (SI6 method)"
    call integrate_trajectory_si6()
  case('si6v')
    write( *, "('INFO      : ',a)" ) "integrating the particle trajectory (SI6v method)"
    call integrate_trajectory_si6v()
  case('si8')
    write( *, "('INFO      : ',a)" ) "integrating the particle trajectory (SI8 method)"
    call integrate_trajectory_si8()
  case('si8v')
    write( *, "('INFO      : ',a)" ) "integrating the particle trajectory (SI8v method)"
    call integrate_trajectory_si8v()
  end select

! display performance information
!
  timer = secnds(timer)
  write( *, "('COMPUTED  : ',a,1pe12.5,a)" ) 'computing done in ', timer, ' seconds'

! deallocate particles
!
  write( *, "('INFO      : ',a)" ) "deallocating the particle"
  call finit_particle()

! deallocate field variables
!
  write( *, "('INFO      : ',a)" ) "deallocating the field components"
  call finit_fields()
!
end program paccel
