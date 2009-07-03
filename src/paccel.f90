!!******************************************************************************
!!
!! Program: PAccel - particle acceleration
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
program paccel

  use params  , only : read_params, idir, odir                                &
                     , ueta, aeta, jcrit, nsteps, xc, yc, zc                  &
                     , ptype, tunit, tmulti                                   &
                     , dens, c, periodic, cfl, dtout, tmax, ethres, current   &
                     , efield, vp
  use mod_hdf5, only : hdf5_init, hdf5_get_dims, hdf5_read_var                &
                     , xmn, ymn, zmn, xmx, ymx, zmx, dxi, dyi, dzi, dx, dy, dz
  use mod_fits, only : fits_put_data_2d, fits_put_data_1d
  use interpolation, only : interpolate => ptricub

  implicit none

! local variables
!
  character(len = 255) :: fl
  integer              :: dm(3), pm(3)
  integer              :: i, j, k, n, nmax
  real                 :: jx, jy, jz, ja, et, xx, dxmin
  real(kind=8)         :: ww, wx, wy, wz, uu, ux, uy, uz, vv
  real(kind=16)        :: ax, ay, az, mu
  real                 :: xli, yli, zli
  real                 :: xt, yt, zt, xi, yi, zi, xr, yr, zr
  real(kind=16)        :: mev, ulen, plen, bavg, qom, va, rg, tg, en
  real(kind=16)        :: xp, yp, zp, vx, vy, vz, xp1, yp1, zp1, vx1, vy1, vz1, vu
  real(kind=16)        :: tm, dt, dtp, dt1, dt2, fc, aa, bb, vc
  logical              :: per = .false., fin = .false.

! allocatable arrays
!
  real, dimension(:,:,:), allocatable :: qt, bx, by, bz, ex, ey, ez
  real, dimension(:,:)  , allocatable :: x, v
  real, dimension(:)    , allocatable :: t, e
!
!-------------------------------------------------------------------------------
!
! read parameters
!
  write (*,'(a)') '------------------------------------------------------------------------------'
  write (*,'(a)') '===        PAccel algorithm started         =================================='
  write (*,'(a)') '===    Copyright (C) 2008 Grzegorz Kowal    =================================='
  write (*,*)
  write( *, "('TASK      : ',a)" ) "integrating trajectories of charged particle"
  write( *, "('INFO      : ',a)" ) "reading parameters"
  call read_params

! print some info
!
  write( *, "('INDIR     : ',a)" ) trim(idir)
  write( *, "('OUDIR     : ',a)" ) trim(odir)

! compute coefficients
!
  mev   = 938.27199893682302445085952058434
  va    = 299792.45799999998416751623153687  / c
  bavg  = 0.13694624848330305688648422801634 * sqrt(dens) / c
  qom   = 9578.8340668294185888953506946564  * bavg * tmulti
  tg    = 0.00065594468630975164436663904510283 / tmulti
  ulen  = va * tmulti
  select case(tunit)
  case('m')
    qom  = 60.0 * qom
    ulen = 60.0 * ulen
    tg   = tg / 60.0
  case('h')
    qom  = 3600.0 * qom
    ulen = 3600.0 * ulen
    tg   = tg / 3600.0
  case('d')
    qom  = 86400.0 * qom
    ulen = 86400.0 * ulen
    tg   = tg / 86400.0
  case('w')
    qom  = 604800.0 * qom
    ulen = 604800.0 * ulen
    tg   = tg / 604800.0
  case('y')
    qom  = 31556925.974678400903940200805664 * qom
    ulen = 31556925.974678400903940200805664 * ulen
    tg   = tg / 31556925.974678400903940200805664
  case default
  end select

  plen = 3.2407792896656063447898874599466e-14 * ulen
  mu   = 0.5d0 * (vp * c)**2 / bavg

  select case(ptype)
  case ('e')
    mev  = 0.00054461702326790686814333986021097 * mev
    qom  = 1836.152667427881851835991255939 * qom
    tg   = tg / 1836.152667427881851835991255939
    write( *, "('INFO      : trajectory for electron')" )
  case default
    write( *, "('INFO      : trajectory for proton')" )
  end select

  write( *, "('INFO      : Va    = ',1pe15.8,' km/s')" ) va
  write( *, "('INFO      : c     = ',1pe15.8,' Va')" ) c
  write( *, "('INFO      : <B>   = ',1pe15.8,' G')" ) bavg
  write( *, "('INFO      : e/m   = ',1pe15.8,' [1 / G ',a1,']')" ) qom, tunit
  write( *, "('INFO      : Om    = ',1pe15.8,' [1 / ',a1,']')" ) qom*bavg, tunit
  write( *, "('INFO      : Tg    = ',1pe15.8,' [',a1,']')" ) tg, tunit
  write( *, "('INFO      : mu    = ',1pe15.8,' [',a1,']')" ) mu, tunit
  write( *, "('INFO      : mu/Om = ',1pe15.8,' [',a1,']')" ) mu/(qom*bavg), tunit
  write( *, "('INFO      : L     = ',1pe15.8,' km')" ) ulen
  write( *, "('INFO      : L     = ',1pe15.8,' pc')" ) plen

! check if periodic box
!
  if (periodic .eq. 'y') per = .true.

! initialize dimensions
!
  dm(:) = 1

! init HDF5 interface and get data dimensions
!
  call hdf5_init()
  call hdf5_get_dims(dm)

! allocate variables
!
  write( *, "('INFO      : ',a)" ) "allocating variables"
  allocate(bx(dm(1),dm(2),dm(3)))
  allocate(by(dm(1),dm(2),dm(3)))
  allocate(bz(dm(1),dm(2),dm(3)))
  if (efield .eq. 'y') then
    allocate(ex(dm(1),dm(2),dm(3)))
    allocate(ey(dm(1),dm(2),dm(3)))
    allocate(ez(dm(1),dm(2),dm(3)))
  endif

! read magnetic field components
!
  write( *, "('INFO      : ',a)" ) "reading magnetic field"
  call hdf5_read_var('magx', bx)
  call hdf5_read_var('magy', by)
  call hdf5_read_var('magz', bz)

  if (efield .eq. 'y') then

! computing current density
!
    write( *, "('INFO      : ',a)" ) "computing electric field"
    ex(:,:,:) = 0.0
    ey(:,:,:) = 0.0
    ez(:,:,:) = 0.0

    if (ueta .gt. 0.0 .and. current .eq. 'y') then
      write( *, "('INFO      : ',a)" ) "computing current density"

      do k = 2, dm(3)-1
        do j = 2, dm(2)-1
          do i = 2, dm(1)-1
            jx = dyi*(bz(i,j+1,k) - bz(i,j-1,k)) - dzi*(by(i,j,k+1) - by(i,j,k-1))
            jy = dzi*(bx(i,j,k+1) - bx(i,j,k-1)) - dxi*(bz(i+1,j,k) - bz(i-1,j,k))
            jz = dxi*(by(i+1,j,k) - by(i-1,j,k)) - dyi*(bx(i,j+1,k) - bx(i,j-1,k))

            et = ueta
            if (aeta .gt. 0.0) then
              ja = sqrt(jx*jx + jy*jy + jz*jz)
              xx = ja / jcrit
              fc = 0.5 * (tanh(20.0*(xx - 1.0)) + 1.0)
              et = ueta + aeta * fc * xx
            endif

            ex(i,j,k) = et * jx
            ey(i,j,k) = et * jy
            ez(i,j,k) = et * jz
          enddo
        enddo
      enddo
    endif

! read velocity field components
!
    allocate(qt(dm(1),dm(2),dm(3)))

    call hdf5_read_var('velx', qt)
    ey = ey + qt * bz
    ez = ez - qt * by

    call hdf5_read_var('vely', qt)
    ex = ex - qt * bz
    ez = ez + qt * bx

    call hdf5_read_var('velz', qt)
    ex = ex + qt * by
    ey = ey - qt * bx

    if (allocated(qt)) deallocate(qt)
  endif

! allocate particle positions & velocities
!
  nmax = int(tmax / dtout) + 2
  allocate(t(nmax))
  allocate(e(nmax))
  allocate(x(3,nmax))
  allocate(v(3,nmax))

! set initial positions and speeds
!
  xp  = xc
  yp  = yc
  zp  = zc

! convert position to index
!
  xt = xli * (xp - xmn)
  yt = yli * (yp - ymn)
  zt = zli * (zp - zmn)
    if (per) then
      xr = pm(1) * (xt - floor(xt)) + 1
      yr = pm(2) * (yt - floor(yt)) + 1
      zr = pm(3) * (zt - floor(zt)) + 1
    else
      xr = pm(1) * xt + 1
      yr = pm(2) * yt + 1
      zr = pm(3) * zt + 1

      if (xr .lt. 1) then
        vx = - vx
        xp = xmn
        xt = xli * (xp - xmn)
        xr = pm(1) * xt + 1
      endif
      if (xr .gt. dm(1)) then
        vx = - vx
        xp = xmx
        xt = xli * (xp - xmn)
        xr = pm(1) * xt + 1
      endif

      if (yr .lt. 1) then
        vy = - vy
        yp = ymn
        yt = yli * (yp - ymn)
        yr = pm(2) * yt + 1
      endif
      if (yr .gt. dm(2)) then
        vy = - vy
        yp = ymx
        yt = yli * (yp - ymn)
        yr = pm(2) * yt + 1
      endif

      if (zr .lt. 1) then
        vz = - vz
        zp = zmn
        zt = zli * (zp - zmn)
        zr = pm(3) * zt + 1
      endif
      if (zr .gt. dm(3)) then
        vz = - vz
        zp = zmx
        zt = zli * (zp - zmn)
        zr = pm(3) * zt + 1
      endif
    endif

! calculate magnetic field components in the initial position
!
  wx = interpolate(dm, bx, xr, yr, zr)
  wy = interpolate(dm, by, xr, yr, zr)
  wz = interpolate(dm, bz, xr, yr, zr)

! calculate the direction of local field
!
  ww = sqrt(wx*wx + wy*wy + wz*wz)
  if (ww .gt. 0.0d0) then
    wx = wx / ww
    wy = wy / ww
    wz = wz / ww
  else
    write( *, "('ERROR     : ',a)" ) "B=0 at the initial position! Choose another one."
    stop
  endif

  ux = 0.0d0
  uy = 0.0d0
  uz = 1.0d0

  vx = uy * wz - uz * wy
  vy = uz * wx - ux * wz
  vz = ux * wy - uy * wx
  vv = sqrt(vx*vx + vy*vy + vz*vz)
  if (vv .gt. 0.0d0) then
    vx = vx / vv
    vy = vy / vv
    vz = vz / vv
  else
    write( *, "('ERROR     : ',a)" ) "V=0 at the initial position! Choose another one."
    stop
  endif

  vx  = vp * c * vx
  vy  = vp * c * vy
  vz  = vp * c * vz

  dtp = 1.0e-16
  tm  = 0.0
  n   = 1

  vu = vx**2 + vy**2 + vz**2
  fc = qom * sqrt(1.0 - min(1.0, vu / c**2))
  vu = sqrt(vu)

  x(1,1) = xp
  x(2,1) = yp
  x(3,1) = zp
  v(1,1) = vx
  v(2,1) = vy
  v(3,1) = vz

! check minimum dx
!
  dxmin = min(dx, dy, dz)
  xli   = 1.0 / (xmx - xmn)
  yli   = 1.0 / (ymx - ymn)
  zli   = 1.0 / (zmx - zmn)

  pm = dm - 1

! print headers
!
  write (*,"('PROGRESS  : ',a8,2x,4(a14))") 'ITER', 'TIME', 'TIMESTEP', 'SPEED (c)', 'ENERGY (MeV)'

! integrate particles
!
  do while (tm .le. tmax .and. n .le. nmax .and. .not. fin)

! update time
!
    tm = tm + dt

!! 1st step of RK integration
!!
! compute relativistic factor
!
    vu = vx**2 + vy**2 + vz**2
    fc = qom * sqrt(1.0 - min(1.0, vu / c**2))
    vu = sqrt(vu)

! convert particle position to arrays index
!
    xt = xli * (xp - xmn)
    yt = yli * (yp - ymn)
    zt = zli * (zp - zmn)

    if (per) then
      xr = pm(1) * (xt - floor(xt)) + 1
      yr = pm(2) * (yt - floor(yt)) + 1
      zr = pm(3) * (zt - floor(zt)) + 1
    else
      xr = pm(1) * xt + 1
      yr = pm(2) * yt + 1
      zr = pm(3) * zt + 1

      if (xr .lt. 1) then
        vx = - vx
        xp = xmn
        xt = xli * (xp - xmn)
        xr = pm(1) * xt + 1
      endif
      if (xr .gt. dm(1)) then
        vx = - vx
        xp = xmx
        xt = xli * (xp - xmn)
        xr = pm(1) * xt + 1
      endif

      if (yr .lt. 1) then
        vy = - vy
        yp = ymn
        yt = yli * (yp - ymn)
        yr = pm(2) * yt + 1
      endif
      if (yr .gt. dm(2)) then
        vy = - vy
        yp = ymx
        yt = yli * (yp - ymn)
        yr = pm(2) * yt + 1
      endif

      if (zr .lt. 1) then
        vz = - vz
        zp = zmn
        zt = zli * (zp - zmn)
        zr = pm(3) * zt + 1
      endif
      if (zr .gt. dm(3)) then
        vz = - vz
        zp = zmx
        zt = zli * (zp - zmn)
        zr = pm(3) * zt + 1
      endif
    endif

! interpolate field at particle position
!
    if (efield .eq. 'y') then
      ux = interpolate(dm, ex, xr, yr, zr)
      uy = interpolate(dm, ey, xr, yr, zr)
      uz = interpolate(dm, ez, xr, yr, zr)
    endif
    wx = interpolate(dm, bx, xr, yr, zr)
    wy = interpolate(dm, by, xr, yr, zr)
    wz = interpolate(dm, bz, xr, yr, zr)

! compute force components
!
  if (efield .eq. 'y') then
    ax = fc * (ux + vy * wz - vz * wy)
    ay = fc * (uy + vz * wx - vx * wz)
    az = fc * (uz + vx * wy - vy * wx)
  else
    ax = fc * (vy * wz - vz * wy)
    ay = fc * (vz * wx - vx * wz)
    az = fc * (vx * wy - vy * wx)
  endif

! integrate velocity
!
    vx1 = vx + dt * ax
    vy1 = vy + dt * ay
    vz1 = vz + dt * az

! integrate position
!
    xp1 = xp + dt * vx
    yp1 = yp + dt * vy
    zp1 = zp + dt * vz

! compute new timestep
!
    aa  = sqrt(ax * ax + ay * ay + az * az)
    bb  = sqrt(wx * wx + wy * wy + wz * wz)
    dt1 = min(tg / max(bb, 1.0e-16), (c - vu) / max(aa, 1.0e-16), sqrt(dxmin / max(aa, 1.0e-16)), dxmin / max(vu, 1.0e-16))


!! 2nd step of RK integration
!!
! compute relativistic factor
!
    vu = vx1**2 + vy1**2 + vz1**2
    fc = qom * sqrt(1.0 - min(1.0, vu / c**2))
    vu = sqrt(vu)

! convert particle position to arrays index
!
    xt = xli * (xp1 - xmn)
    yt = yli * (yp1 - ymn)
    zt = zli * (zp1 - zmn)

    if (per) then
      xr = pm(1) * (xt - floor(xt)) + 1
      yr = pm(2) * (yt - floor(yt)) + 1
      zr = pm(3) * (zt - floor(zt)) + 1
    else
      xr = max(1.0, min(real(dm(1)), pm(1) * xt + 1.0))
      yr = max(1.0, min(real(dm(2)), pm(2) * yt + 1.0))
      zr = max(1.0, min(real(dm(3)), pm(3) * zt + 1.0))
    endif

! interpolate field at particle position
!
    if (efield .eq. 'y') then
      ux = interpolate(dm, ex, xr, yr, zr)
      uy = interpolate(dm, ey, xr, yr, zr)
      uz = interpolate(dm, ez, xr, yr, zr)
    endif
    wx = interpolate(dm, bx, xr, yr, zr)
    wy = interpolate(dm, by, xr, yr, zr)
    wz = interpolate(dm, bz, xr, yr, zr)

! compute force components
!
    if (efield .eq. 'y') then
      ax = fc * (ux + vy * wz - vz * wy)
      ay = fc * (uy + vz * wx - vx * wz)
      az = fc * (uz + vx * wy - vy * wx)
    else
      ax = fc * (vy * wz - vz * wy)
      ay = fc * (vz * wx - vx * wz)
      az = fc * (vx * wy - vy * wx)
    endif

! integrate velocity
!
    vx = 0.5 * (vx + vx1 + dt * ax)
    vy = 0.5 * (vy + vy1 + dt * ay)
    vz = 0.5 * (vz + vz1 + dt * az)

! integrate position
!
    xp = 0.5 * (xp + xp1 + dt * vx)
    yp = 0.5 * (yp + yp1 + dt * vy)
    zp = 0.5 * (zp + zp1 + dt * vz)

! compute new timestep
!
    aa  = sqrt(ax * ax + ay * ay + az * az)
    bb  = sqrt(wx * wx + wy * wy + wz * wz)
    dt2 = min(tg / max(bb, 1.0e-16), (c - vu) / max(aa, 1.0e-16), sqrt(dxmin / max(aa, 1.0e-16)), dxmin / max(vu, 1.0e-16))

! compute new timestep
!
    dt  = min(2.0 * dtp, cfl * min(dt1, dt2))
    dtp = dt

! compute particle energy
!
    vc = (vu/c)**2
    en = mev * vc / sqrt(1.0 - vc)

! copy data to array
!
    if (tm .ge. (n*dtout)) then

      t(n)   = tm

      x(1,n) = xp
      x(2,n) = yp
      x(3,n) = zp

      v(1,n) = vx
      v(2,n) = vy
      v(3,n) = vz

      e(n)   = en

      write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, tm, dt, vu/c, e(n), char(13)
      n = n + 1
    endif

    if (en .ge. ethres) fin = .true.

  enddo

100 continue

  t(n)   = tm

  x(1,n) = xp
  x(2,n) = yp
  x(3,n) = zp

  v(1,n) = vx
  v(2,n) = vy
  v(3,n) = vz

  e(n)   = en

  write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6))") n, tm, dt, vu/c, e(n)

! write positions and speeds to a file
!
  write( *, "('INFO      : ',a)" ) 'writing time, positions and speeds'
  write(fl, '("!",a)') trim(odir) // trim(ptype) // '_tim.fits.gz'
  call fits_put_data_1d(fl, t(1:n))
  write(fl, '("!",a)') trim(odir) // trim(ptype) // '_ene.fits.gz'
  call fits_put_data_1d(fl, e(1:n))
  write(fl, '("!",a)') trim(odir) // trim(ptype) // '_pos.fits.gz'
  call fits_put_data_2d(fl, x(:,1:n))
  write(fl, '("!",a)') trim(odir) // trim(ptype) // '_vel.fits.gz'
  call fits_put_data_2d(fl, v(:,1:n))

! deallocate variables
!
  if (allocated(t )) deallocate(t )
  if (allocated(x )) deallocate(x )
  if (allocated(v )) deallocate(v )
  if (allocated(bx)) deallocate(bx)
  if (allocated(by)) deallocate(by)
  if (allocated(bz)) deallocate(bz)
  if (allocated(ex)) deallocate(ex)
  if (allocated(ey)) deallocate(ey)
  if (allocated(ez)) deallocate(ez)

end program paccel
