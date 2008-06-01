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
                     , ueta, aeta, jcrit, nsteps, xc, yc, zc, dt, qom, c
  use mod_hdf5, only : hdf5_init, hdf5_get_dims, hdf5_read_var                &
                     , xmn, ymn, zmn, xmx, ymx, zmx, dxi, dyi, dzi
  use mod_fits, only : fits_put_data_2d

  implicit none

! local variables
!
  character(len = 255) :: fl
  integer              :: dm(3)
  integer              :: i, j, k, n
  real                 :: jx, jy, jz, ja, et, xx, fc
  real                 :: xp, yp, zp, vx, vy, vz, vp
  real                 :: xi, yi, zi

! allocatable arrays
!
  real, dimension(:,:,:), allocatable :: vv, bx, by, bz, ex, ey, ez
  real, dimension(:,:)  , allocatable :: x, v
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
  allocate(x(3,nsteps))
  allocate(v(3,nsteps))

! set initial positions and speeds
!
  x(1,1) = xc
  x(2,1) = yc
  x(3,1) = zc
  v(1,1) = 0.0
  v(2,1) = 0.0
  v(3,1) = 0.0

  xp = xc
  yp = yc
  zp = zc
  vx = 0.0
  vy = 0.0
  vz = 0.1

! integrate particles
!
! F = q/m * (E + VpxB)
  do n = 1, nsteps

! interpolate fields at the particle position
!
!     i = int(dxi * (xp - xmn))
!     j = int(dyi * (yp - ymn))
!     k = int(dzi * (zp - zmn))

    xi = (xp - xmn) / (xmx - xmn)
    yi = (yp - ymn) / (ymx - ymn)
    zi = (zp - zmn) / (zmx - zmn)

    i = (dm(1) - 1) * (xi - floor(xi)) + 1
    j = (dm(2) - 1) * (yi - floor(yi)) + 1
    k = (dm(3) - 1) * (zi - floor(zi)) + 1

    if (i .lt. 1 .or. i .gt. dm(1)) goto 100
    if (j .lt. 1 .or. j .gt. dm(2)) goto 100
    if (k .lt. 1 .or. k .gt. dm(3)) goto 100

! compute relativistic factor
!
    vp = vx**2 + vy**2 + vz**2
    fc = dt * qom * sqrt(1.0 - min(1.0, vp / c**2))

! integrate velocity
!
    v(1,n) = vx + fc * (ex(i,j,k) + vy * bz(i,j,k) - vz * by(i,j,k))
    v(2,n) = vy + fc * (ey(i,j,k) + vz * bx(i,j,k) - vx * bz(i,j,k))
    v(3,n) = vz + fc * (ez(i,j,k) + vx * by(i,j,k) - vy * bx(i,j,k))
    vx = v(1,n)
    vy = v(2,n)
    vz = v(3,n)

! integrate position
!
    x(1,n) = xp + dt * vx
    x(2,n) = yp + dt * vy
    x(3,n) = zp + dt * vz

    xp = x(1,n)
    yp = x(2,n)
    zp = x(3,n)
  enddo

100 continue

! write positions and speeds to a file
!
  write( *, "('INFO      : ',a)" ) 'writing positions and speeds'
  write(fl, '("!",a)') trim(odir) // 'ppos.fits.gz'
  call fits_put_data_2d(fl, x(:,1:n-1))
  write(fl, '("!",a)') trim(odir) // 'pvel.fits.gz'
  call fits_put_data_2d(fl, v(:,1:n-1))

! deallocate variables
!
  if (allocated(x )) deallocate(x )
  if (allocated(v )) deallocate(v )
  if (allocated(bx)) deallocate(bx)
  if (allocated(by)) deallocate(by)
  if (allocated(bz)) deallocate(bz)
  if (allocated(ex)) deallocate(ex)
  if (allocated(ey)) deallocate(ey)
  if (allocated(ez)) deallocate(ez)

end program paccel
