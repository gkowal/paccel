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
                     , ptype, tunit, dens, c, periodic, cfl, dtout, tmax
  use mod_hdf5, only : hdf5_init, hdf5_get_dims, hdf5_read_var                &
                     , xmn, ymn, zmn, xmx, ymx, zmx, dxi, dyi, dzi, dx, dy, dz
  use mod_fits, only : fits_put_data_2d, fits_put_data_1d

  implicit none

! local variables
!
  character(len = 255) :: fl
  integer              :: dm(3)
  integer              :: i, j, k, n, nmax
  real                 :: jx, jy, jz, ja, et, xx, fc
  real                 :: xp, yp, zp, vx, vy, vz, vp
  real                 :: ax, ay, az, aa, bb, xi, yi, zi, dxmin
  real                 :: bavg, qom, va, rg, tg, ulen, tm, dt, dtp
  logical              :: per = .false., fin = .false.

! allocatable arrays
!
  real, dimension(:,:,:), allocatable :: vv, bx, by, bz, ex, ey, ez
  real, dimension(:,:)  , allocatable :: x, v
  real, dimension(:)    , allocatable :: t
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
  va    = 299792.45799999998416751623153687  / c
  bavg  = 0.13694624848330305688648422801634 * sqrt(dens) / c
  qom   = 9578.8340668294185888953506946564  * bavg
  tg    = 0.00065594468630975164436663904510283
  ulen  = va
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
  select case(ptype)
  case ('e')
    qom  = 1836.152667427881851835991255939 * qom
    tg   = tg / 1836.152667427881851835991255939
    write( *, "('INFO      : trajectory for electron')" )
  case default
    write( *, "('INFO      : trajectory for proton')" )
  end select
  write( *, "('INFO      : Va  = ',1pe15.8,' km/s')" ) va
  write( *, "('INFO      : c   = ',1pe15.8,' Va')" ) c
  write( *, "('INFO      : <B> = ',1pe15.8,' G')" ) bavg
  write( *, "('INFO      : e/m = ',1pe15.8,' [1 / G ',a1,']')" ) qom, tunit
  write( *, "('INFO      : Tg  = ',1pe15.8,' [',a1,']')" ) tg, tunit
  write( *, "('INFO      : L   = ',1pe15.8,' km')" ) ulen
  write( *, "('INFO      : L   = ',1pe15.8,' pc')" ) ulen * 3.2407793e-14


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

! check minimum dx
!
  dxmin = min(dx, dy, dz)

! allocate variables
!
  write( *, "('INFO      : ',a)" ) "allocating variables"
  allocate(bx(dm(1),dm(2),dm(3)))
  allocate(by(dm(1),dm(2),dm(3)))
  allocate(bz(dm(1),dm(2),dm(3)))
  allocate(ex(dm(1),dm(2),dm(3)))
  allocate(ey(dm(1),dm(2),dm(3)))
  allocate(ez(dm(1),dm(2),dm(3)))

! read magnetic field components
!
  write( *, "('INFO      : ',a)" ) "reading magnetic field"
  call hdf5_read_var('magx', bx)
  call hdf5_read_var('magy', by)
  call hdf5_read_var('magz', bz)

! computing current density
!
  write( *, "('INFO      : ',a)" ) "computing current density"
  ex(:,:,:) = 0.0
  ey(:,:,:) = 0.0
  ez(:,:,:) = 0.0

  if (ueta .gt. 0.0) then
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
  allocate(vv(dm(1),dm(2),dm(3)))

  call hdf5_read_var('velx', vv)
  ey = ey + vv * bz
  ez = ez - vv * by

  call hdf5_read_var('vely', vv)
  ex = ex - vv * bz
  ez = ez + vv * bx

  call hdf5_read_var('velz', vv)
  ex = ex + vv * by
  ey = ey - vv * bx

  if (allocated(vv)) deallocate(vv)

! allocate particle positions & velocities
!
  nmax = int(tmax / dtout) + 1
  allocate(t(nmax))
  allocate(x(3,nmax))
  allocate(v(3,nmax))

! set initial positions and speeds
!
  xp  = xc
  yp  = yc
  zp  = zc
  vx  = 0.0
  vy  = 0.0
  vz  = 0.1
  vp  = 0.1
  dtp = 1.0e-16
  tm  = 0.0
  n   = 1

  x(1,1) = xp
  x(2,1) = yp
  x(3,1) = zp
  v(1,1) = vx
  v(2,1) = vy
  v(3,1) = vz

! integrate particles
!
! F = q/m * (E + VpxB)
  do while (tm .le. tmax .and. n .le. nmax .and. .not. fin)

! interpolate fields at the particle position
!
    xi = (xp - xmn) / (xmx - xmn)
    yi = (yp - ymn) / (ymx - ymn)
    zi = (zp - zmn) / (zmx - zmn)

    if (per) then
      i = (dm(1) - 1) * (xi - floor(xi)) + 1
      j = (dm(2) - 1) * (yi - floor(yi)) + 1
      k = (dm(3) - 1) * (zi - floor(zi)) + 1
    else
      i = (dm(1) - 1) * xi + 1
      j = (dm(2) - 1) * yi + 1
      k = (dm(3) - 1) * zi + 1

      if (i .lt. 1 .or. i .gt. dm(1)) goto 100
      if (j .lt. 1 .or. j .gt. dm(2)) goto 100
      if (k .lt. 1 .or. k .gt. dm(3)) goto 100
    endif

! compute relativistic factor
!
    vp = vx**2 + vy**2 + vz**2
    fc = qom * sqrt(1.0 - min(1.0, vp / c**2))
    vp = sqrt(vp)

! compute force components
!
    ax = fc * (ex(i,j,k) + vy * bz(i,j,k) - vz * by(i,j,k))
    ay = fc * (ey(i,j,k) + vz * bx(i,j,k) - vx * bz(i,j,k))
    az = fc * (ez(i,j,k) + vx * by(i,j,k) - vy * bx(i,j,k))
    aa = sqrt(ax*ax + ay*ay + az*az)
    bb = sqrt(bx(i,j,k)**2 + by(i,j,k)**2 + bz(i,j,k)**2)

! compute new timestep
!
    dt  = cfl * min(tg / bb, dxmin / max(vp, 1.0e-8))
!     print *, dt, bb
!     dt  = cfl * min(sqrt(dxmin / max(aa, 1.0e-8)), dxmin / max(vp, 1.0e-8))
    dt  = min(2.0 * dtp, dt)
    dtp = dt

! update time
!
    tm = tm + dt

! integrate velocity
!
    vx = vx + dt * ax
    vy = vy + dt * ay
    vz = vz + dt * az

! integrate position
!
    xp = xp + dt * vx
    yp = yp + dt * vy
    zp = zp + dt * vz

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

      n = n + 1
      print *, n, tm, dt, vp/c

      if (vp .ge. c) fin = .true.
    endif

  enddo

100 continue

! write positions and speeds to a file
!
  write( *, "('INFO      : ',a)" ) 'writing time, positions and speeds'
  write(fl, '("!",a)') trim(odir) // trim(ptype) // '_tim.fits.gz'
  call fits_put_data_1d(fl, t(1:n-1))
  write(fl, '("!",a)') trim(odir) // trim(ptype) // '_pos.fits.gz'
  call fits_put_data_2d(fl, x(:,1:n-1))
  write(fl, '("!",a)') trim(odir) // trim(ptype) // '_vel.fits.gz'
  call fits_put_data_2d(fl, v(:,1:n-1))

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
