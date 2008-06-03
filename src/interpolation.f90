!!******************************************************************************
!!
!! module: interpolation - subroutines to interpolation 3D data.
!!
!! Copyright (C) 2008 Grzegorz Kowal <kowal@astro.wisc.edu>
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
module interpolation

  implicit none

!
!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! pconst: subroutine for piecewise triconstant interpolation
!
!===============================================================================
!
  real function ptricon(u, x, y, z) result (v)

    implicit none

! input/output arguments
!
    real, dimension(:,:,:)  :: u
    real                    :: x, y, z

! local variables
!
    integer :: i, j, k, d(3)
!
!------------------------------------------------------------------------------
!
    d = size(u)

    i = min(d(1), max(1, int(x)))
    j = min(d(2), max(1, int(y)))
    k = min(d(3), max(1, int(z)))

    v = u(i,j,k)

  end function ptricon
!
!===============================================================================
!
! ptrilin: subroutine for piecewise trilinear interpolation
!
!===============================================================================
!
  real function ptrilin(u, x, y, z) result (v)

    implicit none

! input/output arguments
!
    real, dimension(:,:,:)  :: u
    real                    :: x, y, z

! local variables
!
    integer :: i, j, k, ip1, jp1, kp1, d(3)
    real    :: dx, dy, dz, u1, u2, v1, v2, w1, w2
!
!------------------------------------------------------------------------------
!
    d = size(u)

    x   = min(x, real(d(1)))
    y   = min(y, real(d(2)))
    z   = min(z, real(d(3)))

    i   = int(x)
    j   = int(y)
    k   = int(z)
    ip1 = min(i + 1, d(1))
    jp1 = min(j + 1, d(2))
    kp1 = min(k + 1, d(3))
    dx  = x - i
    dy  = y - j
    dz  = z - k

    u1  = (1.0 - dz) * u(i  ,j  ,k) + dz * u(i  ,j  ,kp1)
    u2  = (1.0 - dz) * u(i  ,jp1,k) + dz * u(i  ,jp1,kp1)
    v1  = (1.0 - dz) * u(ip1,j  ,k) + dz * u(ip1,j  ,kp1)
    v2  = (1.0 - dz) * u(ip1,jp1,k) + dz * u(ip1,jp1,kp1)

    w1  = (1.0 - dy) * u1 + dy * u2
    w2  = (1.0 - dy) * v1 + dy * v2

    v   = (1.0 - dx) * w1 + dx * w2

  end function ptrilin
!
!===============================================================================
!
! ptricub: subroutine for piecewise tricubic interpolation
!
!===============================================================================
!
  real function ptricub(u, x, y, z) result (v)

    implicit none

! input/output arguments
!
    real, dimension(:,:,:)  :: u
    real                    :: x, y, z

! local variables
!
    integer :: i, j, k, d(3)
!
!------------------------------------------------------------------------------
!
    d = size(u)

    i = min(d(1), max(1, int(x)))
    j = min(d(2), max(1, int(y)))
    k = min(d(3), max(1, int(z)))

    v = u(i,j,k)

  end function ptricub

end module
