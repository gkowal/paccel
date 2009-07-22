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
  real function ptricon(d, u, x, y, z) result (v)

    implicit none

! input/output arguments
!
    integer, dimension(3)               :: d
    real   , dimension(d(1),d(2),d(3))  :: u
    real                                :: x, y, z

! local variables
!
    integer :: i, j, k
!
!------------------------------------------------------------------------------
!
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
  real function ptrilin(d, u, x, y, z) result (v)

    implicit none

! input/output arguments
!
    integer, dimension(3)               :: d
    real   , dimension(d(1),d(2),d(3))  :: u
    real(kind=16)                       :: x, y, z

! local variables
!
    integer :: i, j, k, ip1, jp1, kp1
    real    :: dx, dy, dz, u1, u2, v1, v2, w1, w2
    logical :: out = .false.
!
!------------------------------------------------------------------------------
!
      if (x .lt. 1) out = .true.
      if (y .lt. 1) out = .true.
      if (z .lt. 1) out = .true.

      if (x .gt. d(1)) out = .true.
      if (y .gt. d(2)) out = .true.
      if (z .gt. d(3)) out = .true.

    if (out) then
      print *, x, y, z
      stop
    endif

    x   = max(1.0, min(x, real(d(1))))
    y   = max(1.0, min(y, real(d(2))))
    z   = max(1.0, min(z, real(d(3))))

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
  real function ptricub(dm, u, x, y, z) result (v)

    use params  , only : periodic

    implicit none

! input/output arguments
!
    integer, dimension(3)    , intent(in) :: dm
    real   , dimension(:,:,:), intent(in) :: u
    real(kind=16)            , intent(in) :: x, y, z

! local flags
!
    logical, save :: first = .true.
    logical, save :: per   = .false.

! local saved variables
!
    integer(kind=4), dimension(3), save :: pm

! local variables
!
    integer :: i, j, k, im1, jm1, km1, ip1, jp1, kp1, ip2, jp2, kp2
    real    :: dx, dy, dz
    real    :: tm1m1, tm1c, tm1p1, tm1p2, tcm1, tcc, tcp1, tcp2     &
             , tp1m1, tp1c, tp1p1, tp1p2, tp2m1, tp2c, tp2p1, tp2p2 &
             , um1, uc, up1, up2
!
!------------------------------------------------------------------------------
!
    if (first) then
      if (periodic .eq. 'y') per = .true.

      pm = dm - 1

      first = .false.
    endif

! calculate indices
!
    i   = floor(x)
    j   = floor(y)
    k   = floor(z)

    im1 = i - 1
    jm1 = j - 1
    km1 = k - 1

    ip1 = i + 1
    jp1 = j + 1
    kp1 = k + 1

    ip2 = i + 2
    jp2 = j + 2
    kp2 = k + 2

! adjust according to the periodicity of the box
!
    if (per) then
      if (i .eq.     1) im1 = im1 + dm(1)
      if (j .eq.     1) jm1 = jm1 + dm(2)
      if (k .eq.     1) km1 = km1 + dm(3)

      if (i .eq. pm(1)) ip2 = ip2 - dm(1)
      if (j .eq. pm(2)) jp2 = jp2 - dm(2)
      if (k .eq. pm(3)) kp2 = kp2 - dm(3)

      if (i .eq. dm(1)) then
        ip1 = ip1 - dm(1)
        ip2 = ip2 - dm(1)
      endif
      if (j .eq. dm(2)) then
        jp1 = jp1 - dm(2)
        jp2 = jp2 - dm(2)
      endif
      if (k .eq. dm(3)) then
        kp1 = kp1 - dm(3)
        kp2 = kp2 - dm(3)
      endif
    else
      im1 = max(im1,     1)
      jm1 = max(jm1,     1)
      km1 = max(km1,     1)

      ip1 = min(ip1, dm(1))
      jp1 = min(jp1, dm(2))
      kp1 = min(kp1, dm(3))

      ip2 = min(ip2, dm(1))
      jp2 = min(jp2, dm(2))
      kp2 = min(kp2, dm(3))
    endif

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
!
!===============================================================================
!
! pos2index: subroutine converts a give position to the array index
!
!===============================================================================
!
  subroutine pos2index(xp, yp, zp, px, py, pz, out)

    use params  , only : fformat, periodic
    use mod_fits, only : fits_get_dims, fits_get_bounds
    use mod_hdf5, only : hdf5_get_dims, hdf5_get_bounds

    implicit none

! input and output arguments
!
    real(kind=16), intent(in)  :: xp, yp, zp
    real(kind=16), intent(out) :: px, py, pz
    logical      , intent(out) :: out

! local flags
!
    logical, save :: first = .true.
    logical, save :: per   = .false.

! local saved variables
!
    integer(kind=4), dimension(3), save :: dm, pm
    real(kind=4)                 , save :: xmn, xmx, ymn, ymx, zmn, zmx
    real(kind=16)                , save :: xli, yli, zli

! local temporary variables
!
    real(kind=16)                       :: xt, yt, zt
!
!------------------------------------------------------------------------------
!
    if (first) then
      if (periodic .eq. 'y') per = .true.

! obtain the box dimensions and bounds
!
      select case(fformat)
      case('fits')
        call fits_get_dims(dm)
        call fits_get_bounds(xmn, xmx, ymn, ymx, zmn, zmx)
      case('hdf5')
        call hdf5_get_dims(dm)
        call hdf5_get_bounds(xmn, xmx, ymn, ymx, zmn, zmx)
      case default
        stop
      end select

! calculate dimensions needed for convertion
!
      pm(:) = dm(:) - 1

! calculate factors
!
      xli   = 1.0 / (xmx - xmn)
      yli   = 1.0 / (ymx - ymn)
      zli   = 1.0 / (zmx - zmn)

      first = .false.
    endif

! set the out of box flag
!
    out = .false.

! convert position to index
!
    xt = xli * (xp - xmn)
    yt = yli * (yp - ymn)
    zt = zli * (zp - zmn)

! calculate indices
!
    if (per) then
      px = pm(1) * (xt - floor(xt)) + 1
      py = pm(2) * (yt - floor(yt)) + 1
      pz = pm(3) * (zt - floor(zt)) + 1
    else
      px = pm(1) * xt + 1
      py = pm(2) * yt + 1
      pz = pm(3) * zt + 1

! check if the particle left the box
!
      if (px .lt. 1    ) out = .true.
      if (py .lt. 1    ) out = .true.
      if (pz .lt. 1    ) out = .true.

      if (px .gt. dm(1)) out = .true.
      if (py .gt. dm(2)) out = .true.
      if (pz .gt. dm(3)) out = .true.
    endif

  end subroutine pos2index

end module
