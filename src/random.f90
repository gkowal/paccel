!!******************************************************************************
!!
!! module: random - handles random number generators by Marsaglia & Tsang
!!
!! references:
!! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
!! random variables', J. Statist. Software, v5(8).
!!
!! Copyright (C) 2007-2009 Grzegorz Kowal <kowal@astro.wisc.edu>
!!
!!******************************************************************************
!!
!!  This file is part of Godunov.
!!
!!  Godunov is free software; you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation; either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  Godunov is distributed in the hope that it will be useful,
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
! N.B. It is assumed that all integers are 32-bit.
! N.B. The value of M2 has been halved to compensate for the lack of
!      unsigned integers in Fortran.

! Latest version - 1 January 2001

! Parallel version - October 2006

! This version has been customised for parallel processing use,
! specifically with OpenMP.  Each thread uses its own pseudo-random
! sequence. (Gib Bogle)
!--------------------------------------------------------------------------

module random

  implicit none

  private

  integer                                   , save :: nseeds
  integer(kind=4), dimension(:), allocatable, save :: seeds

  public :: initgen, nseeds, get_seeds, randomu, randomz, randomn

contains
!
!======================================================================
!
!  Initialize generator
!
subroutine initgen(np, seed)

  integer                           , intent(in) :: np
  integer(kind=4), dimension(0:np-1), intent(in) :: seed

  integer  :: k
!
!----------------------------------------------------------------------
!
  nseeds = np

  if (allocated(seeds)) then
    if (np .ne. nseeds) then
      print *, 'INITGEN: np != nseeds'
      print *, '       : seeds for random numbers are already allocated but their size'
      print *, '       : do not match np!'
    endif
  else
    allocate(seeds(0:np-1))
  endif

  do k = 0, np-1
    seeds(k) = seed(k)
  enddo
  return
end subroutine initgen
!
!======================================================================
!
!  Get seeds
!
subroutine get_seeds(seed)

  integer(kind=4), dimension(0:nseeds-1), intent(out) :: seed

  integer  :: k
!
!----------------------------------------------------------------------
!
  do k = 0, nseeds-1
    seed(k) = seeds(k)
  enddo
  return
end subroutine get_seeds
!
!======================================================================
!
!  Generate uniformly distributed random numbers, sequence kp
!  range 0.0...1.0
!
function randomu(kp) result(val)

  integer(kind=4) :: kp, jz, jsr
  real            :: val
!
!----------------------------------------------------------------------
!
  jsr = seeds(kp)
  jz  = jsr

  jsr = ieor( jsr, ishft( jsr,  13 ) )
  jsr = ieor( jsr, ishft( jsr, -17 ) )
  jsr = ieor( jsr, ishft( jsr,   5 ) )

  seeds(kp) = jsr

  val = 0.5 + 0.23283064365e-9*(jz + jsr)
  return
end function randomu
!
!======================================================================
!
!  Generate uniformly distributed random numbers, sequence kp
!  range -0.5...0.5
!
function randomz(kp) result(val)

  integer(kind=4) :: kp, jz, jsr
  real            :: val
!
!----------------------------------------------------------------------
!
  jsr = seeds(kp)
  jz  = jsr

  jsr = ieor( jsr, ishft( jsr,  13 ) )
  jsr = ieor( jsr, ishft( jsr, -17 ) )
  jsr = ieor( jsr, ishft( jsr,   5 ) )

  seeds(kp) = jsr

  val = 0.23283064365e-9*(jz + jsr)
  return
end function randomz
!
!======================================================================
!
!  Generate uniformly distributed random numbers, sequence kp
!  range -1.0...1.0
!
function randomn(kp) result(val)

  integer(kind=4) :: kp, jz, jsr
  real            :: val
!
!----------------------------------------------------------------------
!
  jsr = seeds(kp)
  jz  = jsr

  jsr = ieor( jsr, ishft( jsr,  13 ) )
  jsr = ieor( jsr, ishft( jsr, -17 ) )
  jsr = ieor( jsr, ishft( jsr,   5 ) )

  seeds(kp) = jsr

  val = 0.46566128730e-9*(jz + jsr)
  return
end function randomn

end module random
