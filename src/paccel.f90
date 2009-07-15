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
                     , efield, vpar, vper, approx
  use mod_hdf5, only : hdf5_init, hdf5_get_dims, hdf5_read_var                &
                     , xmn, ymn, zmn, xmx, ymx, zmx, dxi, dyi, dzi, dx, dy, dz, dtc
  use mod_fits, only : fits_put_data_2d, fits_put_data_1d
  use interpolation, only : interpolate => ptricub, pos2index

  implicit none

! local variables
!
  character(len = 255) :: fl
  integer              :: dm(3)
  integer              :: i, j, k, n, nmax
  real(kind=16)        :: va, dn, mu0, bavg, qom, om, tg, vp, pc, rg, gm, mp   &
                        , mu, ln, dr, ts, yy, bt, pu, del, vr, om0
  real(kind=16)        :: rx, ry, rz, xr, yr, zr
  real                 :: jx, jy, jz, ja, et, xx, dxmin
  real(kind=8)         :: ww, wx, wy, wz, uu, ux, uy, uz, vv
  real(kind=16)        :: ax, ay, az
  real(kind=16)        :: r1x, r1y, r1z, r2x, r2y, r2z, r3x, r3y, r3z, r4x, r4y, r4z
  real(kind=16)        :: v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, v4x, v4y, v4z
  real(kind=16)        :: p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, p4x, p4y, p4z
  real(kind=16)        :: k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z, k4x, k4y, k4z
  real                 :: xt, yt, zt, xi, yi, zi
  real(kind=16)        :: mev, ulen, plen, en
  real(kind=16)        :: vx, vy, vz, px, py, pz, xp1, yp1, zp1, vx1, vy1, vz1, vu
  real(kind=16)        :: tm, dt, dtp, dt1, dt2, fc, aa, bb, vc
  logical              :: per = .false., fin = .false., out

! parameters
!
  real(kind=16)        :: cc  = 299792457.99999998416751623153687d0
  real(kind=16)        :: pi  = 3.1415926535897931159979634685442d0
  real(kind=16)        :: pi2 = 6.2831853071795862319959269370884d0

! allocatable arrays
!
  real, dimension(:,:,:), allocatable :: qt, bx, by, bz, ex, ey, ez
!
!-------------------------------------------------------------------------------
!
! read parameters
!
  write (*,'(a)') '------------------------------------------------------------------------------'
  write (*,'(a)') '===        PAccel algorithm started         =================================='
  write (*,'(a)') '===    Copyright (C) 2008-2009 Grzegorz Kowal    ============================='
  write (*,*)
  write( *, "('TASK      : ',a)" ) "integrating trajectories of charged particle"
  write( *, "('INFO      : ',a)" ) "reading parameters"
  call read_params

! check if parameters are correct
!
  if (c .lt. 1.0d0) then
    write( *, "('ERROR     : ',a)" ) "parameter c must be larger than zero!"
    stop
  endif

! initialize dimensions
!
  dm(:) = 1

! init HDF5 interface and get data dimensions
!
  call hdf5_init()
  call hdf5_get_dims(dm)

! print some info
!
  write( *, "('INDIR     : ',a)" ) trim(idir)
  write( *, "('OUDIR     : ',a)" ) trim(odir)

! compute plasma parameters
!                                                       ! c is expressed in Va
  gm    = 1.0d0 / sqrt(1.0d0 - (1.0 / c)**2)            ! Lorentz factor
  va    = gm * cc  / c                                  ! Alfven speed [m/s]
  dn    = 1.6726215850718025379202284485224e-21 * dens  ! density conversion from
                                                        ! protonmass/cm^3 to kg/m^3
  mu0   = 125.66370614359171042906382353976             ! magnetic permeability [Gs^2 m s^2 / kg]
  bavg  = va * sqrt(mu0 * dn)                           ! magnetic field strength [Gs]

! print plasma parametes
!
  write( *, "('INFO      : plasma parameters:')" )
  write( *, "('INFO      : c     =',1pe15.8,' [Va]')"       ) c
  write( *, "('INFO      : Va    =',1pe15.8,' [m / s]')"    ) va
  write( *, "('INFO      : dens  =',1pe15.8,' [u / cm^3] =',1pe15.8,' [kg / m^3]')" ) dens, dn
  write( *, "('INFO      : <B>   =',1pe15.8,' [G]')"        ) bavg

! compute particle parameters
!
  qom   = 9578.8340668294185888953506946564d0           ! e/m [1 / Gs s]
  select case(ptype)
  case ('e')
    qom  = 1836.152667427881851835991255939 * qom
    mp   = 9.1093818871545313708798643833606e-31        ! electron mass [kg]
    mev  = 0.51099890307660134070033564057667           ! rest energy of electron [MeV]
  case ('p')
    mp   = 1.6726215850718025086476640481627e-27        ! proton mass [kg]
    mev  = 938.27199893682302445085952058434            ! rest energy of proton [MeV]
  case default
    mp   = 1.6726215850718025086476640481627e-27        ! proton mass [kg]
    mev  = 938.27199893682302445085952058434            ! rest energy of proton [MeV]
  end select
  vp = cc * vpar                                        ! parallel particle speed
  vr = cc * vper                                        ! perpendicular particle speed
  gm = 1.0d0 / sqrt(1.0d0 - vpar**2 - vper**2)          ! Lorentz factor
  om0 = qom * bavg                                      ! classical gyrofrequency
  om = om0 / gm                                         ! relativistic gyrofrequency
  tg = 1.0d0 / om                                       ! gyroperiod
  tg = 2.0 * pi * tg
  rg = vp / om                                          ! gyroradius (Larmor radius)
  pc = 3.2407792896656065765177783686188e-17            ! 1 meter [pc]
  mu = 0.5d0 * mp * vp**2 / bavg                        ! magnetic moment [kg m^2 / s^2 Gs]
  yy = 3.168876464084018437308447107767e-08             ! 1 second [yr]
  en = gm * mev

! print particle parameters
!
  write( *, "('INFO      : particle parameters:')" )
  select case(ptype)
  case ('e')
    write( *, "('INFO      : trajectory for electron')" )
  case default
    write( *, "('INFO      : trajectory for proton')" )
  end select
  write( *, "('INFO      : e/m   =',1pe15.8,' [1 / G s]')" ) qom
  write( *, "('INFO      : Vpar  =',1pe15.8,' [c] =',1pe15.8,' [m / s]')" ) vpar, vp
  write( *, "('INFO      : Vper  =',1pe15.8,' [c] =',1pe15.8,' [m / s]')" ) vper, vr
  write( *, "('INFO      : gamma =',1pe15.8)"              ) gm
  write( *, "('INFO      : Om    =',1pe15.8,' [1 / s]')"   ) om
  write( *, "('INFO      : Tg    =',1pe15.8,' [s]')"       ) tg
  write( *, "('INFO      : Rg    =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) rg, pc * rg
  write( *, "('INFO      : mu    =',1pe15.8,' [N m / Gs]')") mu
  write( *, "('INFO      : E0    =',1pe15.8,' [MeV]')"     ) en

! calculate geometry parameters
!
  ln = va                                               ! the size of the box
  dr = ln * min(dx, dy, dz)
  ts = dtc

! change time unit
!
  select case(tunit)
  case('u')
    fc = 1.0d-6
  case('s')
    fc = 1.0
  case('m')
    fc = 60.0
  case('h')
    fc = 3600.0
  case('d')
    fc = 86400.0
  case('w')
    fc = 604800.0
  case('y')
    fc = 31556925.974678400903940200805664
  case default
    fc = 1.0
  end select

  fc  = tmulti * fc
  qom = qom    * fc
  ln  = va     * fc
  dr  = dr     * fc
  ts  = ts     * fc

! print geometry parameters
!
  write( *, "('INFO      : geometry parameters:')" )
  write( *, "('INFO      : T     =',1pe15.8,' [s] =',1pe15.8,' [yr]')" ) fc, yy * fc
  write( *, "('INFO      : dt    =',1pe15.8,' [s] =',1pe15.8,' [yr]')" ) ts, yy * ts
  write( *, "('INFO      : L     =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) ln, pc * ln
  write( *, "('INFO      : dx    =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) dr, pc * dr


! calculate conditions
!

! print conditions
!
  write( *, "('INFO      : conditions:')" )
  write( *, "('INFO      : Rg/L  =',1pe15.8)" ) rg / ln
  write( *, "('INFO      : Rg/dx =',1pe15.8)" ) rg / dr
  write( *, "('INFO      : Tg/dt =',1pe15.8)" ) tg / ts

! write parameters to info.txt
!
  open  (10, file = 'info.txt', form = 'formatted', status = 'replace')

! print plasma parametes
!
  write (10, "('INFO      : plasma parameters:')" )
  write (10, "('INFO      : c     =',1pe15.8,' [Va]')"       ) c
  write (10, "('INFO      : Va    =',1pe15.8,' [m / s]')"    ) va
  write (10, "('INFO      : dens  =',1pe15.8,' [u / cm^3] =',1pe15.8,' [kg / m^3]')" ) dens, dn
  write (10, "('INFO      : <B>   =',1pe15.8,' [G]')"        ) bavg

  write (10, "('INFO      : particle parameters:')" )
  select case(ptype)
  case ('e')
    write (10, "('INFO      : trajectory for electron')" )
  case default
    write (10, "('INFO      : trajectory for proton')" )
  end select
  write (10, "('INFO      : e/m   =',1pe15.8,' [1 / G s]')" ) qom
  write (10, "('INFO      : Om    =',1pe15.8,' [1 / s]')"   ) om
  write (10, "('INFO      : Tg    =',1pe15.8,' [s]')"       ) tg
  write (10, "('INFO      : Vpar  =',1pe15.8,' [c] =',1pe15.8,' [m / s]')" ) vpar, vp
  write (10, "('INFO      : Vper  =',1pe15.8,' [c] =',1pe15.8,' [m / s]')" ) vper, vr
  write (10, "('INFO      : Rg    =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) rg, pc * rg
  write (10, "('INFO      : gamma =',1pe15.8)"              ) gm
  write (10, "('INFO      : mu    =',1pe15.8,' [N m / Gs]')") mu
  write (10, "('INFO      : E0    =',1pe15.8,' [MeV]')"     ) en

! print geometry parameters
!
  write (10, "('INFO      : geometry parameters:')" )
  write (10, "('INFO      : T     =',1pe15.8,' [s] =',1pe15.8,' [yr]')" ) fc, yy * fc
  write (10, "('INFO      : dt    =',1pe15.8,' [s] =',1pe15.8,' [yr]')" ) ts, yy * ts
  write (10, "('INFO      : L     =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) ln, pc * ln
  write (10, "('INFO      : dx    =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) dr, pc * dr

! print conditions
!
  write (10, "('INFO      : conditions:')" )
  write (10, "('INFO      : Rg/L  =',1pe15.8)" ) rg / ln
  write (10, "('INFO      : Rg/dx =',1pe15.8)" ) rg / dr
  write (10, "('INFO      : Tg/dt =',1pe15.8)" ) tg / ts

  close (10)

! convert e/m to the units of magnetic field
!
  qom = qom * bavg

! check if periodic box
!
  if (periodic .eq. 'y') per = .true.

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

! relativistic Lorentz force
!
!  dp/dt = q ( E + (u/c x B) )
!
!  p = g m u -> u = p / g m
!
!  g = ( 1 - ( u / c )^2 )^(-1/2)
!
!  dp / dt = q ( E + ( p / (g m c) x B ) )
!
!  1 / c dp / dt = q ( E + ( p / (g m c^2) x B ) )

! allocate particle positions & velocities
!
  nmax = int(tmax / dtout) + 2

! set the initial integration parameters
!
  tm  = 0.0
  dtp = 1.0e-16
  dt  = dtp
  n   = 1

! set the initial particle position
!
  rx  = xc
  ry  = yc
  rz  = zc

! convert position to index
!
  call pos2index(rx, ry, rz, xr, yr, zr, out)

  if (out) then
    write( *, "('ERROR     : ',a)" ) "The initial position is out of the box! Choose another one."
    stop
  endif

! calculate magnetic field components in the initial position
!
  wx = interpolate(dm, bx, xr, yr, zr)
  wy = interpolate(dm, by, xr, yr, zr)
  wz = interpolate(dm, bz, xr, yr, zr)

! calculate the direction of the local magnetic field
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

! calculate the perpendicular velocity component
!
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

! calculate the initial particle velocity in the resting reference frame
!
  vx  = ( vper * vx + vpar * wx ) * c
  vy  = ( vper * vy + vpar * wy ) * c
  vz  = ( vper * vz + vpar * wz ) * c

! calculate the Lorentz factor
!
  vu = sqrt(vx**2 + vy**2 + vz**2)
  bt = vu / c
  gm = 1.0 / sqrt(1.0d0 - min(1.0d0, bt * bt))

! calculate the initial particle momentuum
!
  px = gm * vx
  py = gm * vy
  pz = gm * vz

! calculate the initial particle energy
!
  en = gm * mev

! set the perpendicular speed
!
  vp = vpar * c
  vr = vper * c

! calculate gyroperiod and gyroradius
!
  om = om0 * ww / gm
  tg = pi2 / om
  rg = vp * va / om

! print the progress information
!
  write ( *,"('PROGRESS  : ',a8,2x,4(a14))") 'ITER', 'TIME', 'TIMESTEP', 'SPEED (c)', 'ENERGY (MeV)'
  write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, tm, dt, vu/c, en, char(13)

! print headers and the initial values
!
  open  (10, file = 'output.dat', form = 'formatted', status = 'replace')
  write (10, "('#',1a16,22a18)") 'Time', 'X', 'Y', 'Z', 'Px', 'Py', 'Pz'       &
             , 'Vx', 'Vy', 'Vz', '|V|', '|Vpar|', '|Vper|', '|V|/c'            &
             , '|Vpar|/c', '|Vper|/c', 'gamma', 'En [MeV]', '<B> [Gs]'         &
             , 'Omega [1/s]', 'Tg [s]', 'Rg [m]', 'Rg/L'
  write (10, "(23(1pe18.10))") tm, rx, ry, rz, mp*px, mp*py, mp*pz, vx, vy, vz &
             , vu, vr, vp, vu/c, vr/c, vp/c, gm, en, bavg*ww, om, tg, rg, rg/ln
  close (10)

! iterpolate electric field at the initial position
!
  if (efield .eq. 'y') then
    ux = interpolate(dm, ex, xr, yr, zr)
    uy = interpolate(dm, ey, xr, yr, zr)
    uz = interpolate(dm, ez, xr, yr, zr)
  endif

! check minimum dx
!
  dxmin = min(dx, dy, dz)

! integrate particles
!
  do while (tm .le. tmax .and. n .le. nmax .and. .not. fin)

! update time
!
    tm = tm + dt

!! 1st step of the RK integration
!!
! integrate the momentum
!
    p1x = px
    p1y = py
    p1z = pz

! calculate the Lorentz factor
!
! p2 = gm2 * u2 = u2 / (1 - (u2/c2))
! p2 - p2 u2 / c2 = u2
! p2 = u2 (1 + p2/c2)
! u2 = p2 / (1 + p2/c2)
! u  = p / sqrt(1 + p2/c2)
!
! p = gm u = sqrt(1 + p2/c2) u
! gm = sqrt(1 + p2/c2)
!
    pu = sqrt(p1x**2 + p1y**2 + p1z**2)
    bt = pu / c
    gm = sqrt( 1.0d0 + bt * bt )

! calculate the velocity
!
    v1x = p1x / gm
    v1y = p1y / gm
    v1z = p1z / gm

! integrate the position
!
    r1x = rx
    r1y = ry
    r1z = rz

! convert position to index
!
    call pos2index(r1x, r1y, r1z, xr, yr, zr, out)

    if (out) goto 100

! interpolate field components at the particle position
!
    if (efield .eq. 'y') then
      ux = interpolate(dm, ex, xr, yr, zr)
      uy = interpolate(dm, ey, xr, yr, zr)
      uz = interpolate(dm, ez, xr, yr, zr)
    endif
    wx = interpolate(dm, bx, xr, yr, zr)
    wy = interpolate(dm, by, xr, yr, zr)
    wz = interpolate(dm, bz, xr, yr, zr)

! compute the force components
!
    if (efield .eq. 'y') then
      ax = (ux + v1y * wz - v1z * wy)
      ay = (uy + v1z * wx - v1x * wz)
      az = (uz + v1x * wy - v1y * wx)
    else
      ax = (v1y * wz - v1z * wy)
      ay = (v1z * wx - v1x * wz)
      az = (v1x * wy - v1y * wx)
    endif

! set the first term
!
    k1x = qom * ax
    k1y = qom * ay
    k1z = qom * az

!! 2nd step of the RK integration
!!
! integrate the momentum
!
    p2x = px + 0.5 * dt * k1x
    p2y = py + 0.5 * dt * k1y
    p2z = pz + 0.5 * dt * k1z

! calculate the Lorentz factor
!
    pu = sqrt(p2x**2 + p2y**2 + p2z**2)
    bt = pu / c
    gm = sqrt( 1.0d0 + bt * bt )

! calculate the velocity
!
    v2x = p2x / gm
    v2y = p2y / gm
    v2z = p2z / gm

! integrate the position
!
    r2x = rx + 0.25 * dt * ( vx + v2x )
    r2y = ry + 0.25 * dt * ( vy + v2y )
    r2z = rz + 0.25 * dt * ( vz + v2z )

! convert position to index
!
    call pos2index(r2x, r2y, r2z, xr, yr, zr, out)

    if (out) goto 100

! interpolate field components at the particle position
!
    if (efield .eq. 'y') then
      ux = interpolate(dm, ex, xr, yr, zr)
      uy = interpolate(dm, ey, xr, yr, zr)
      uz = interpolate(dm, ez, xr, yr, zr)
    endif
    wx = interpolate(dm, bx, xr, yr, zr)
    wy = interpolate(dm, by, xr, yr, zr)
    wz = interpolate(dm, bz, xr, yr, zr)

! compute the force components
!
    if (efield .eq. 'y') then
      ax = (ux + v2y * wz - v2z * wy)
      ay = (uy + v2z * wx - v2x * wz)
      az = (uz + v2x * wy - v2y * wx)
    else
      ax = (v2y * wz - v2z * wy)
      ay = (v2z * wx - v2x * wz)
      az = (v2x * wy - v2y * wx)
    endif

! the second term
!
    k2x = qom * ax
    k2y = qom * ay
    k2z = qom * az

!! 3rd step of the RK integration
!!
! integrate the momentum
!
    p3x = px + 0.5 * dt * k2x
    p3y = py + 0.5 * dt * k2y
    p3z = pz + 0.5 * dt * k2z

! calculate the Lorentz factor
!
    pu = sqrt(p3x**2 + p3y**2 + p3z**2)
    bt = pu / c
    gm = sqrt( 1.0d0 + bt * bt )

! calculate the velocity
!
    v3x = p3x / gm
    v3y = p3y / gm
    v3z = p3z / gm

! integrate the position
!
    r3x = rx + 0.25 * dt * ( vx + v2x )
    r3y = ry + 0.25 * dt * ( vy + v2y )
    r3z = rz + 0.25 * dt * ( vz + v2z )

! convert position to index
!
    call pos2index(r3x, r3y, r3z, xr, yr, zr, out)

    if (out) goto 100

! interpolate field components at the particle position
!
    if (efield .eq. 'y') then
      ux = interpolate(dm, ex, xr, yr, zr)
      uy = interpolate(dm, ey, xr, yr, zr)
      uz = interpolate(dm, ez, xr, yr, zr)
    endif
    wx = interpolate(dm, bx, xr, yr, zr)
    wy = interpolate(dm, by, xr, yr, zr)
    wz = interpolate(dm, bz, xr, yr, zr)

!! 3rd step of the RK integration
!!
! compute the force components
!
    if (efield .eq. 'y') then
      ax = (ux + v3y * wz - v3z * wy)
      ay = (uy + v3z * wx - v3x * wz)
      az = (uz + v3x * wy - v3y * wx)
    else
      ax = (v3y * wz - v3z * wy)
      ay = (v3z * wx - v3x * wz)
      az = (v3x * wy - v3y * wx)
    endif

! the third term
!
    k3x = qom * ax
    k3y = qom * ay
    k3z = qom * az

!! 4th step of the RK integration
!!
! integrate the momentum
!
    p4x = px + dt * k3x
    p4y = py + dt * k3y
    p4z = pz + dt * k3z

! calculate the Lorentz factor
!
    pu = sqrt(p4x**2 + p4y**2 + p4z**2)
    bt = pu / c
    gm = sqrt( 1.0d0 + bt * bt )

! calculate the velocity
!
    v4x = p4x / gm
    v4y = p4y / gm
    v4z = p4z / gm

! integrate the position
!
    r4x = rx + dt * v4x
    r4y = ry + dt * v4y
    r4z = rz + dt * v4z

! convert position to index
!
    call pos2index(r4x, r4y, r4z, xr, yr, zr, out)

    if (out) goto 100

! interpolate field components at the particle position
!
    if (efield .eq. 'y') then
      ux = interpolate(dm, ex, xr, yr, zr)
      uy = interpolate(dm, ey, xr, yr, zr)
      uz = interpolate(dm, ez, xr, yr, zr)
    endif
    wx = interpolate(dm, bx, xr, yr, zr)
    wy = interpolate(dm, by, xr, yr, zr)
    wz = interpolate(dm, bz, xr, yr, zr)

! compute the force components
!
    if (efield .eq. 'y') then
      ax = (ux + v4y * wz - v4z * wy)
      ay = (uy + v4z * wx - v4x * wz)
      az = (uz + v4x * wy - v4y * wx)
    else
      ax = (v4y * wz - v4z * wy)
      ay = (v4z * wx - v4x * wz)
      az = (v4x * wy - v4y * wx)
    endif

! the fourth term
!
    k4x = qom * ax
    k4y = qom * ay
    k4z = qom * az

!! the final integration of the particle velocity and position
!!
    px = px + dt * ( k1x + 2.0 * ( k2x + k3x ) + k4x ) / 6.0
    py = py + dt * ( k1y + 2.0 * ( k2y + k3y ) + k4y ) / 6.0
    pz = pz + dt * ( k1z + 2.0 * ( k2z + k3z ) + k4z ) / 6.0

! calculate the Lorentz factor
!
    pu = sqrt(px**2 + py**2 + pz**2)
    bt = pu / c
    gm = sqrt( 1.0d0 + bt * bt )
    vu = pu / gm

! calculate the velocity
!
    vx = px / gm
    vy = py / gm
    vz = pz / gm

! calcuate the particle position
!
    rx = rx + dt * ( vx + 2.0 * ( v1x + v2x ) + v3x ) / 6.0
    ry = ry + dt * ( vy + 2.0 * ( v1y + v2y ) + v3y ) / 6.0
    rz = rz + dt * ( vz + 2.0 * ( v1z + v2z ) + v3z ) / 6.0

! convert position to index
!
    call pos2index(r4x, r4y, r4z, xr, yr, zr, out)

    if (out) goto 100

! interpolate field components at the particle position
!
    if (efield .eq. 'y') then
      ux = interpolate(dm, ex, xr, yr, zr)
      uy = interpolate(dm, ey, xr, yr, zr)
      uz = interpolate(dm, ez, xr, yr, zr)
    endif
    wx = interpolate(dm, bx, xr, yr, zr)
    wy = interpolate(dm, by, xr, yr, zr)
    wz = interpolate(dm, bz, xr, yr, zr)

! compute the force components
!
    if (efield .eq. 'y') then
      ax = (ux + vy * wz - vz * wy)
      ay = (uy + vz * wx - vx * wz)
      az = (uz + vx * wy - vy * wx)
    else
      ax = (vy * wz - vz * wy)
      ay = (vz * wx - vx * wz)
      az = (vx * wy - vy * wx)
    endif

! calculate the particle energy
!
    en = gm * mev

! new time step
!
    k4x = k4x - qom * ax
    k4y = k4y - qom * ay
    k4z = k4z - qom * az

    del = sqrt(k4x*k4x + k4y*k4y + k4z*k4z) / 6.0
    dt1 = dt * sqrt(cfl * 1.0e-4 / del)

    dt  = min(2.0 * dt, dt1)

! copy data to array
!
    if (tm .ge. (n*dtout)) then

! calculate the perpendicular particle speed
!
      ww = sqrt(wx*wx + wy*wy + wz*wz)
      wx = wx / ww
      wy = wy / ww
      wz = wz / ww

      ux = vy * wz - vz * wy
      uy = vz * wx - vx * wz
      uz = vx * wy - vy * wx

      vp = sqrt(ux*ux + uy*uy + uz*uz)
      vr = sqrt(vu*vu - vp*vp)

! calculate gyroperiod and gyroradius
!
      om = om0 * ww / gm
      tg = pi2 / om
      rg = vp * va / om

! write the progress
!
      write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, tm, dt, vu/c, en, char(13)

! write results to the output file
!
      open  (10, file = 'output.dat', form = 'formatted', position = 'append')
      write (10, "(23(1pe18.10))") tm, rx, ry, rz, mp*px, mp*py, mp*pz, vx, vy &
          , vz, vu, vr, vp, vu/c, vr/c, vp/c, gm, en, bavg*ww, om, tg, rg, rg/ln
      close (10)

      n = n + 1
    endif

    if (en .ge. ethres) fin = .true.

  enddo

100 continue

! calculate the perpendicular particle speed
!
  ww = sqrt(wx*wx + wy*wy + wz*wz)
  wx = wx / ww
  wy = wy / ww
  wz = wz / ww

  ux = vy * wz - vz * wy
  uy = vz * wx - vx * wz
  uz = vx * wy - vy * wx

  vp = sqrt(ux*ux + uy*uy + uz*uz)
  vr = sqrt(vu*vu - vp*vp)

! calculate gyroperiod and gyroradius
!
  om = om0 * ww / gm
  tg = pi2 / om
  rg = vp * va / om

! write the progress
!
  write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6))") n, tm, dt, vu/c, en

! write results to the output file
!
  open  (10, file = 'output.dat', form = 'formatted', position = 'append')
  write (10, "(23(1pe18.10))") tm, rx, ry, rz, mp*px, mp*py, mp*pz, vx, vy, vz &
             , vu, vr, vp, vu/c, vr/c, vp/c, gm, en, bavg*ww, om, tg, rg, rg/ln
  close (10)

! deallocate variables
!
  if (allocated(bx)) deallocate(bx)
  if (allocated(by)) deallocate(by)
  if (allocated(bz)) deallocate(bz)
  if (allocated(ex)) deallocate(ex)
  if (allocated(ey)) deallocate(ey)
  if (allocated(ez)) deallocate(ez)

end program paccel
