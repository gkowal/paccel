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

    return
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

    return
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
    integer :: i, j, k, im1, jm1, km1, ip1, jp1, kp1, ip2, jp2, kp2, d(3)
    real    :: dx, dy, dz
    real    :: tm1m1, tm1c, tm1p1, tm1p2, tcm1, tcc, tcp1, tcp2     &
             , tp1m1, tp1c, tp1p1, tp1p2, tp2m1, tp2c, tp2p1, tp2p2 &
             , um1, uc, up1, up2
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
    im1 = max(i - 1, 1)
    jm1 = max(j - 1, 1)
    km1 = max(k - 1, 1)
    ip1 = min(i + 1, d(1))
    jp1 = min(j + 1, d(2))
    kp1 = min(k + 1, d(3))
    ip2 = min(i + 2, d(1))
    jp2 = min(j + 2, d(2))
    kp2 = min(k + 2, d(3))
    dx  = x - i
    dy  = y - j
    dz  = z - k

    tm1m1 = cint(dz, u(im1,jm1,km1), u(im1,jm1,k), u(im1,jm1,kp1), u(im1,jm1,kp2))
    tm1c  = cint(dz, u(im1,j  ,km1), u(im1,j  ,k), u(im1,j  ,kp1), u(im1,j  ,kp2))
    tm1p1 = cint(dz, u(im1,jp1,km1), u(im1,jp1,k), u(im1,jp1,kp1), u(im1,jp1,kp2))
    tm1p2 = cint(dz, u(im1,jp2,km1), u(im1,jp2,k), u(im1,jp2,kp1), u(im1,jp2,kp2))

    tcm1  = cint(dz, u(i  ,jm1,km1), u(i  ,jm1,k), u(i  ,jm1,kp1), u(i  ,jm1,kp2))
    tcc   = cint(dz, u(i  ,j  ,km1), u(i  ,j  ,k), u(i  ,j  ,kp1), u(i  ,j  ,kp2))
    tcp1  = cint(dz, u(i  ,jp1,km1), u(i  ,jp1,k), u(i  ,jp1,kp1), u(i  ,jp1,kp2))
    tcp2  = cint(dz, u(i  ,jp2,km1), u(i  ,jp2,k), u(i  ,jp2,kp1), u(i  ,jp2,kp2))

    tp1m1 = cint(dz, u(ip1,jm1,km1), u(ip1,jm1,k), u(ip1,jm1,kp1), u(ip1,jm1,kp2))
    tp1c  = cint(dz, u(ip1,j  ,km1), u(ip1,j  ,k), u(ip1,j  ,kp1), u(ip1,j  ,kp2))
    tp1p1 = cint(dz, u(ip1,jp1,km1), u(ip1,jp1,k), u(ip1,jp1,kp1), u(ip1,jp1,kp2))
    tp1p2 = cint(dz, u(ip1,jp2,km1), u(ip1,jp2,k), u(ip1,jp2,kp1), u(ip1,jp2,kp2))

    tp2m1 = cint(dz, u(ip2,jm1,km1), u(ip2,jm1,k), u(ip2,jm1,kp1), u(ip2,jm1,kp2))
    tp2c  = cint(dz, u(ip2,j  ,km1), u(ip2,j  ,k), u(ip2,j  ,kp1), u(ip2,j  ,kp2))
    tp2p1 = cint(dz, u(ip2,jp1,km1), u(ip2,jp1,k), u(ip2,jp1,kp1), u(ip2,jp1,kp2))
    tp2p2 = cint(dz, u(ip2,jp2,km1), u(ip2,jp2,k), u(ip2,jp2,kp1), u(ip2,jp2,kp2))

    um1 = cint(dy, tm1m1, tm1c, tm1p1, tm1p2)
    uc  = cint(dy, tcm1 , tcc , tcp1 , tcp2 )
    up1 = cint(dy, tp1m1, tp1c, tp1p1, tp1p2)
    up2 = cint(dy, tp2m1, tp2c, tp2p1, tp2p2)

    v = cint(dx, um1, uc, up1, up2)

    return
  end function ptricub
!
!===============================================================================
!
! cint: function for cubic interpolation
!
!===============================================================================
!
  real function cint(x, um1, uc, up1, up2) result (v)

    implicit none

! input/output arguments
!
    real :: x, um1, uc, up1, up2
!
!------------------------------------------------------------------------------
!
    v =      - x        * (x - 1.0) * (x - 2.0) * um1 &
      + 3.0 * (x + 1.0) * (x - 1.0) * (x - 2.0) * uc  &
      - 3.0 * (x + 1.0) *  x        * (x - 2.0) * up1 &
      +       (x + 1.0) *  x        * (x - 1.0) * up2

    v = v / 6.0

    return
  end function cint

end module
