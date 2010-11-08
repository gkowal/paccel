!!******************************************************************************
!!
!! module: particles - subroutines to prepare and advance particles
!!
!! Copyright (C) 2010 Grzegorz Kowal <grzegorz@gkowal.info>
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
module particles

  implicit none

! domain dimensions
!
  integer, dimension(3)  , save :: dm, qm

! domain bounds
!
  real   , dimension(3,2), save :: bnds

! domain size
!
  real   , dimension(3)  , save :: bsiz

! particle mass and the speed of light
!
  real(kind=8)           , save :: mrest, qom, c2, om0, fc, ln, bavg, bpar

! arrays containing the initial positions and velocities of particle
!
  real(kind=PREC), dimension(3), save :: x0, v0, p0

! array to store the dump times
!
  real(kind=8), dimension(:), save, allocatable :: tt
!
!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! init_particle: subroutine initializes the particle position and momentum
!
!===============================================================================
!
  subroutine init_particle()

    use fields, only : get_dimensions, get_domain_bounds, bx, by, bz
    use params, only : ptype, vpar, vper, c, dens, tunit, tmulti, xc, yc, zc
    use params, only : output, tmin, tmax, ndumps
#ifdef TEST
    use params, only : bini, bamp, vamp, freq
#endif /* TEST */

    implicit none

! local variables
!
    integer      :: p, n
    real(kind=8) :: vp, vr, vv, va
    real(kind=8) :: gm, dn, mu0, om, tg, rg, mu, mp, en, ek, ba
    real(kind=PREC) :: bb, rt

! arrays
!
    real(kind=PREC), dimension(3) :: r0
    real(kind=PREC), dimension(3) :: b, u, w

! position indices
!
    integer        , dimension(4) :: ii, jj, kk
    real(kind=8   ), dimension(4) :: cx, cy, cz
    real(kind=8   ), dimension(3) :: dr

! parameters
!
    real(kind=8) :: pi2 = 6.2831853071795862319959269370884d0
    real(kind=8) :: cc  = 299792457.99999998416751623153687     ! the speed of light [m/s]
    real(kind=8) :: pc  = 3.2407792896656065765177783686188e-17 ! 1 meter [pc]
    real(kind=8) :: sc  = 3.168876464084018437308447107767e-08  ! 1 second [yr]
!
!-------------------------------------------------------------------------------
!
! get dimain dimensions
!
    call get_dimensions(dm)

#ifdef TRICUB
    qm(:) = dm(:) + 2
#else /* TRUCUB */
    qm(:) = dm(:) + 1
#endif /* TRUCUB */
    if (dm(3) .eq. 1) qm(:) = 1

! get domain bounds
!
    call get_domain_bounds(bnds)

! calculate the domain size
!
    do p = 1, 3
      bsiz(p) = bnds(p,2) - bnds(p,1)
    end do

! compute plasma parameters
!                                                        ! c is expressed in Va
    dn   = 1.6726215850718025379202284485224e-21 * dens  ! density conversion from
                                                         ! protonmass/cm^3 to kg/m^3
    gm   = 1.0d0 / sqrt(1.0d0 - (1.0 / c)**2)            ! Lorentz factor
    va   = gm * cc  / c                                  ! Alfven speed [m/s]
    mu0  = 125.66370614359171042906382353976             ! magnetic permeability [Gs^2 m s^2 / kg]
    bavg = va * sqrt(mu0 * dn)                           ! magnetic field strength [Gs]
    c2    = c * c                                        ! square of the speed of light

! initialize particle parameters
!
    select case(ptype)
    case ('e')
      mrest =  0.51099890307660134070033564057667        ! rest energy of electron [MeV]
      qom   = -17588201.72265790030360221862793          ! e/m [1 / Gs s]
      mp    = 9.1093818871545313708798643833606e-31      ! electron mass [kg]
    case default
      mrest =  938.27199893682302445085952058434         ! rest energy of proton   [MeV]
      qom   =  9578.8340668294185888953506946564         ! e/m [1 / Gs s]
      mp    = 1.6726215850718025086476640481627e-27      ! proton mass [kg]
    end select
    vp = cc * vpar                                       ! parallel particle speed
    vr = cc * vper                                       ! perpendicular particle speed
    vv = sqrt(vpar**2 + vper**2)                         ! absolute velocity
    mu = 0.5d0 * mp * vr**2 / bavg                       ! magnetic moment [kg m^2 / s^2 Gs]
    om0   = abs(qom * bavg)                              ! classical gyrofrequency
    om = om0 / gm                                        ! relativistic gyrofrequency
    tg = 1.0d0 / om                                      ! gyroperiod
    tg = pi2 * tg
    rg = vr / om                                         ! gyroradius (Larmor radius)

! print plasma parametes
!
    write( *, "('INFO      : plasma parameters:')" )
    write( *, "('INFO      : c     =',1pe15.8,' [Va]')"       ) c
    write( *, "('INFO      : Va    =',1pe15.8,' [m / s]')"    ) va
    write( *, "('INFO      : dens  =',1pe15.8,' [u / cm^3] =',1pe15.8,' [kg / m^3]')" ) dens, dn
    write( *, "('INFO      : <B>   =',1pe15.8,' [G]')"        ) bavg

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
    write( *, "('INFO      : |V|   =',1pe15.8,' [c] =',1pe15.8,' [m / s]')" ) vv  , vv * cc
    write( *, "('INFO      : gamma =',1pe15.8)"              ) gm
    write( *, "('INFO      : Om    =',1pe15.8,' [1 / s]')"   ) om
    write( *, "('INFO      : Tg    =',1pe15.8,' [s]')"       ) tg
    write( *, "('INFO      : Rg    =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) rg, pc * rg
    write( *, "('INFO      : mu    =',1pe15.8,' [N m / Gs]')") mu
    write( *, "('INFO      : E0    =',1pe15.8,' [MeV]')"     ) mrest

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
    qom = qom * fc

! calculate geometry parameters
!
    ln = va * fc                                         ! the size of the box

! print geometry parameters
!
    write( *, "('INFO      : geometry parameters:')" )
    write( *, "('INFO      : T     =',1pe15.8,' [s] =',1pe15.8,' [yr]')" ) fc, sc * fc
    write( *, "('INFO      : L     =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) ln, pc * ln
    write( *, "('INFO      : Rg/L  =',1pe15.8)" ) rg / ln
    write( *, "('INFO      : Tg/T  =',1pe15.8)" ) tg / fc

    write( *, "('INFO      : code units:')" )
    write( *, "('INFO      : e/m   =',1pe15.8)" ) qom * bavg

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
    write (10, "('INFO      : |V|   =',1pe15.8,' [c] =',1pe15.8,' [m / s]')" ) vv  , vv * cc
    write (10, "('INFO      : Rg    =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) rg, pc * rg
    write (10, "('INFO      : gamma =',1pe15.8)"              ) gm
    write (10, "('INFO      : mu    =',1pe15.8,' [N m / Gs]')") mu
    write (10, "('INFO      : E0    =',1pe15.8,' [MeV]')"     ) mrest

! print geometry parameters
!
    write (10, "('INFO      : geometry parameters:')" )
    write (10, "('INFO      : T     =',1pe15.8,' [s] =',1pe15.8,' [yr]')" ) fc, sc * fc
    write (10, "('INFO      : L     =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) ln, pc * ln
    write (10, "('INFO      : Rg/L  =',1pe15.8)" ) rg / ln
    write (10, "('INFO      : Tg/T  =',1pe15.8)" ) tg / fc

    write (10, "('INFO      : code units:')" )
    write (10, "('INFO      : e/m   =',1pe15.8)" ) qom * bavg

    close (10)

! convert e/m to the units of magnetic field
!
    qom = qom * bavg

! initial position and velocity
!
    x0(:) = (/ xc, yc, zc /)

#ifdef TEST
#ifdef WTEST
    bpar = sqrt(bini**2 - bamp**2)
    b(1) = bpar
    b(2) = bamp * cos(pi2 * freq * xc)
    b(3) = bamp * sin(pi2 * freq * xc)

    u(1) = 0.0
    u(2) = 0.0
    u(3) = 1.0
#endif /* WTEST */

#ifdef ITEST
    rt   = sqrt(xc * xc + yc * yc)

    if (rt .gt. 0.0) then
      b(1) = - bini * yc / rt
      b(2) =   bini * xc / rt
      b(3) = 0.0
    else
      b(:) = 0.0
    end if
#endif /* ITEST */
#else /* TEST */
! convert position to index
!
    call pos2index(x0, r0)

! prepare coefficients for interpolation
!
    call prepare_interpolation(r0, ii, jj, kk, dr, cx, cy, cz)

! interpolate field components at the particle position
!
    b(1) = interpolate(bx, ii, jj, kk, dr, cx, cy, cz)
    b(2) = interpolate(by, ii, jj, kk, dr, cx, cy, cz)
    b(3) = interpolate(bz, ii, jj, kk, dr, cx, cy, cz)
#endif /* TEST */

! calculate the direction of the local magnetic field
!
    bb = sqrt(dot_product(b, b))
    ba = bb
    if (bb .gt. 0.0d0) then
      b(:) = b(:) / bb
    else
      write( *, "('ERROR     : ',a)" ) "B=0 at the initial position! Choose another one."
      stop
    endif

#ifndef TEST
! calculate the perpendicular unit vector
!
    if (dm(3) .eq. 1) then
      w(1) = 0.0
      w(2) = 0.0
      w(3) = 1.0
    else
      call random_number(w)
      w = w - 0.5
    end if

    bb = sqrt(dot_product(w, w))
    if (bb .gt. 0.0d0) then
      w(:) = w(:) / bb
    else
      write( *, "('ERROR     : ',a)" ) "V=0 at the initial position! Choose another one."
      stop
    end if

    u(1) = w(2) * b(3) - w(3) * b(2)
    u(2) = w(3) * b(1) - w(1) * b(3)
    u(3) = w(1) * b(2) - w(2) * b(1)

    bb = sqrt(dot_product(u, u))
    if (bb .gt. 0.0d0) then
      u(:) = u(:) / bb
    else
      write( *, "('ERROR     : ',a)" ) "V=0 at the initial position! Choose another one."
      stop
    end if
#endif /* !TEST */

! calculate the initial velocity
!
    v0(:) = (vpar * b(:) + vper * u(:)) * c

! calculate the Lorentz factor of the initial state
!
#ifdef RELAT
    gm = 1.0 / sqrt(1.0d0 - dot_product(v0, v0) / c2)
#else
    gm = 1.0
#endif

! calculate the initial particle momentuum
!
    p0(:) = gm * v0(:)

! calculate particle energy
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else
    en = 0.5 * (vpar**2 + vper**2) * c2
    ek = en
#endif

! print headers and the initial values
!
    open  (10, file = 'output.dat', form = 'formatted', status = 'replace')
    write (10, "('#',1a16,19a18)") 'Time', 'X', 'Y', 'Z', 'Vx', 'Vy', 'Vz'     &
                                 , '|V| [c]', '|Vpar| [c]', '|Vper| [c]'       &
                                 , 'gamma', 'En [MeV]', 'Ek [MeV]'             &
                                 , '<B> [Gs]', 'Omega [1/s]'                   &
                                 , 'Tg [s]', 'Rg [m]', 'Tg [T]', 'Rg [L]'
    write (10, "(19(1pe18.10))") 0.0, x0(1), x0(2), x0(3), v0(1), v0(2), v0(3) &
                                    , vv, vpar, vper, gm, en, ek               &
                                    , bavg * ba, om, tg, rg, tg / fc, rg / ln
    close (10)

! prepare dump times
!
    if (output .eq. 'l') then
      n = ndumps * (alog10(tmax) - alog10(tmin)) + 1
      allocate(tt(n))
      do p = 1, n
        tt(p) = 10.0d0**((p - 1.0) / ndumps + alog10(tmin))
      enddo
    endif
!
!-------------------------------------------------------------------------------
!
  end subroutine init_particle
!
!===============================================================================
!
! finit_particle: subroutine deallocates the particle variables
!
!===============================================================================
!
  subroutine finit_particle()

    implicit none
!
!-------------------------------------------------------------------------------
!
    if (allocated(tt)) deallocate(tt)

!-------------------------------------------------------------------------------
!
  end subroutine finit_particle
!
!===============================================================================
!
! integrate_trajectory_rk4: subroutine integrates particle trajectory using
!                           the 4th order RK method
!
!===============================================================================
!
  subroutine integrate_trajectory_rk4()

    use params, only : c, tmin, tmax, rho, maxtol, dtini, dtmax, ndumps,    &
                       vpar, vper

    implicit none

! local variables
!
    integer                       :: n, m
    real(kind=PREC)               ::    t1, t2, t3, t4, t5
    real(kind=PREC), dimension(3) :: x, x1, x2, x3, x4, x5
    real(kind=PREC), dimension(3) :: v, v1, v2, v3, v4, v5
    real(kind=PREC), dimension(3) :: p, p1, p2, p3, p4, p5
    real(kind=PREC), dimension(3) ::    k1, k2, k3, k4, k5
    real(kind=PREC), dimension(3) ::    l1, l2, l3, l4, l5
    real(kind=PREC), dimension(3) :: a, u, b
    real(kind=PREC)               :: gamma
    real(kind=8   )               :: delta
    real(kind=8   )               :: ba, vu, vp, vr, en, ek, om, tg, rg
    real(kind=8   )               :: t, dt, dtq, dtn
!
!-------------------------------------------------------------------------------
!
! initialize time
!
    n   = 0
    m   = 0
    t   = 0.0
    dt  = dtini
    dtq = qom * dt

! initial position and velocity
!
    x(:) = x0(:)
    v(:) = v0(:)
    p(:) = p0(:)

! calculate parameters
!
    vu = sqrt(vpar**2 + vper**2)

! calculate the Lorentz factor
!
    gamma = lorentz_factor(p)

! calculate particle energy
!
#ifdef RELAT
    en = gamma * mrest
    ek = en - mrest
#else
    en = 0.5 * vu**2
    ek = en
#endif

! print the progress information
!
    write (*,"('PROGRESS  : ',a8,2x,4(a14))") 'ITER', 'TIME', 'TIMESTEP', 'SPEED (c)', 'ENERGY (MeV)'
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, vu, ek, char(13)

!== INTEGRATION LOOP ==
!
! integrate particles
!
    do while (t .le. tmax)

!! 1st step of the RK integration
!!
! integrate the position and momentum
!
      t1    = t
      x1(:) = x(:)
      p1(:) = p(:)

! calculate the Lorentz factor
!
      gamma = lorentz_factor(p1)

! calculate velocity
!
      v1(:) = p1(:) / gamma

! calculate acceleration for the location x1 and velocity v1
!
      call acceleration(t1, x1, v1, a, u, b)

! calculate the first term
!
      l1(:) = dt  * v1(:)
      k1(:) = dtq * a (:)

!! 2nd step of the RK integration
!!
! integrate the position and momentum
!
      t2    = t    + 0.5 * dt
      x2(:) = x(:) + 0.5 * l1(:)
      p2(:) = p(:) + 0.5 * k1(:)

! calculate the Lorentz factor
!
      gamma = lorentz_factor(p2)

! calculate the velocity
!
      v2(:) = p2(:) / gamma

! calculate acceleration for the location x2 and velocity v2
!
      call acceleration(t2, x2, v2, a, u, b)

! calculate the second term
!
      l2(:) = dt  * v2(:)
      k2(:) = dtq * a (:)

!! 3rd step of the RK integration
!!
! integrate the position and momentum
!
      t3    = t    + 0.5 * dt
      x3(:) = x(:) + 0.5 * l2(:)
      p3(:) = p(:) + 0.5 * k2(:)

! calculate the Lorentz factor
!
      gamma = lorentz_factor(p3(:))

! calculate the velocity
!
      v3(:) = p3(:) / gamma

! calculate acceleration for the location x3 and velocity v3
!
      call acceleration(t3, x3, v3, a, u, b)

! calculate the third term
!
      l3(:) = dt  * v3(:)
      k3(:) = dtq * a (:)

!! 4th step of the RK integration
!!
! integrate the position and momentum
!
      t4    = t    + 0.5 * dt
      x4(:) = x(:) + l3(:)
      p4(:) = p(:) + k3(:)

! calculate the Lorentz factor
!
      gamma = lorentz_factor(p4(:))

! calculate the velocity
!
      v4(:) = p4(:) / gamma

! calculate acceleration for the location x4 and velocity v4
!
      call acceleration(t4, x4, v4, a, u, b)

! calculate the third term
!
      l4(:) = dt  * v4(:)
      k4(:) = dtq * a (:)

!! the final integration of the particle position and momentum
!!
      t5    = t    + dt
      x5(:) = x(:) + ( l1(:) + 2.0 * ( l2(:) + l3(:) ) + l4(:) ) / 6.0
      p5(:) = p(:) + ( k1(:) + 2.0 * ( k2(:) + k3(:) ) + k4(:) ) / 6.0

! calculate the Lorentz factor
!
      gamma = lorentz_factor(p5(:))

! calculate the velocity
!
      v5(:) = p5(:) / gamma

! calculate acceleration at the location x
!
      call acceleration(t5, x5, v5, a, u, b)

! estimate error
!
      l4(:) = l4(:) - dt  * v5(:)
      k4(:) = k4(:) - dtq * a (:)

      delta = sqrt(dot_product(l4, l4) + dot_product(k4, k4)) / 6.0

! estimate new timestep
!
      dtn   = dt * (rho * maxtol / delta)**0.2

      if (delta .gt. maxtol) then

! repeat integration with this timestep
!
        dt  = dtn
        dtq = qom * dt

      else

! update time
!
        t   = t + dt

! update new timestep
!
        dt  = min(2.0 * dt, dtn, dtmax)
        dtq = qom * dt

! update position, velocity and momentum
!
        x(:) = x5(:)
        v(:) = v5(:)
        p(:) = p5(:)

! copy data to array
!
        if (m .eq. ndumps) then

! separate particle velocity into parallel and perpendicular components
!
          call separate_velocity(v, b, ba, vu, vp, vr)

! calculate the particle gyroperiod and gyroradius
!
          call gyro_parameters(gamma, ba, vr, om, tg, rg)

! calculate particle energy
!
#ifdef RELAT
          en = gamma * mrest
          ek = en - mrest
#else
          en = 0.5 * vu**2
          ek = en
#endif

! write the progress
!
          write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, vu / c, ek, char(13)

! write results to the output file
!
          open  (10, file = 'output.dat', form = 'formatted', position = 'append')
          write (10, "(19(1pe18.10))") t, x(1), x(2), x(3), v(1), v(2), v(3)   &
                                     , vu / c, vp / c, vr / c, gamma, en, ek   &
                                     , bavg * ba, om, tg * fc, rg * ln, tg, rg
          close (10)

          n = n + 1
          m = 0

        end if

! increase data write counter
!
        m = m + 1

      end if

    end do

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(v, b, ba, vu, vp, vr)

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gamma, ba, vr, om, tg, rg)

! calculate particle energy
!
#ifdef RELAT
    en = gamma * mrest
    ek = en - mrest
#else
    en = 0.5 * vu**2
    ek = en
#endif

! write the progress
!
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6))") n, t, dt, vu / c, ek

! write results to the output file
!
    open  (10, file = 'output.dat', form = 'formatted', position = 'append')
    write (10, "(19(1pe18.10))") t, x(1), x(2), x(3), v(1), v(2), v(3)         &
                               , vu / c, vp / c, vr / c, gamma, en, ek         &
                               , bavg * ba, om, tg * fc, rg * ln, tg, rg
    close (10)

!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_rk4
!
!===============================================================================
!
! integrate_trajectory_rk4_log: subroutine integrates particle trajectory using
!                               the 4th order method with snapshots equally
!                               distributed in the logarithm of time
!
!===============================================================================
!
  subroutine integrate_trajectory_rk4_log()

    use params, only : c, tmin, tmax, rho, maxtol, dtini, dtmax, ndumps,    &
                       vpar, vper

    implicit none

! local variables
!
    integer                       :: n, m
    real(kind=PREC)               ::    t1, t2, t3, t4, t5
    real(kind=PREC), dimension(3) :: x, x1, x2, x3, x4, x5, xt
    real(kind=PREC), dimension(3) :: v, v1, v2, v3, v4, v5, vt
    real(kind=PREC), dimension(3) :: p, p1, p2, p3, p4, p5, pt
    real(kind=PREC), dimension(3) ::    k1, k2, k3, k4, k5
    real(kind=PREC), dimension(3) ::    l1, l2, l3, l4, l5
    real(kind=PREC), dimension(3) :: a, u, b
    real(kind=PREC)               :: gamma
    real(kind=8   )               :: delta
    real(kind=8   )               :: ba, vu, vp, vr, en, ek, om, tg, rg
    real(kind=8   )               :: t, dt, dtq, dtn, wl, wr, tp
!
!-------------------------------------------------------------------------------
!
! initialize time
!
    n   = 1
    t   = 0.0
    dt  = dtini
    dtq = qom * dt

! initial position and velocity
!
    x(:) = x0(:)
    v(:) = v0(:)
    p(:) = p0(:)

! calculate parameters
!
    vu = sqrt(vpar**2 + vper**2)

! calculate the Lorentz factor
!
    gamma = lorentz_factor(p)

! calculate particle energy
!
#ifdef RELAT
    en = gamma * mrest
    ek = en - mrest
#else
    en = 0.5 * vu**2
    ek = en
#endif

! print the progress information
!
    write (*,"('PROGRESS  : ',a8,2x,4(a14))") 'ITER', 'TIME', 'TIMESTEP', 'SPEED (c)', 'ENERGY (MeV)'
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, vu, ek, char(13)

!== INTEGRATION LOOP ==
!
! integrate particles
!
    do while (t .lt. tmax)

!! 1st step of the RK integration
!!
! integrate the position and momentum
!
      t1    = t
      x1(:) = x(:)
      p1(:) = p(:)

! calculate the Lorentz factor
!
      gamma = lorentz_factor(p1)

! calculate velocity
!
      v1(:) = p1(:) / gamma

! calculate acceleration for the location x1 and velocity v1
!
      call acceleration(t1, x1, v1, a, u, b)

! calculate the first term
!
      l1(:) = dt  * v1(:)
      k1(:) = dtq * a (:)

!! 2nd step of the RK integration
!!
! integrate the position and momentum
!
      t2    = t    + 0.5 * dt
      x2(:) = x(:) + 0.5 * l1(:)
      p2(:) = p(:) + 0.5 * k1(:)

! calculate the Lorentz factor
!
      gamma = lorentz_factor(p2)

! calculate the velocity
!
      v2(:) = p2(:) / gamma

! calculate acceleration for the location x2 and velocity v2
!
      call acceleration(t2, x2, v2, a, u, b)

! calculate the second term
!
      l2(:) = dt  * v2(:)
      k2(:) = dtq * a (:)

!! 3rd step of the RK integration
!!
! integrate the position and momentum
!
      t3    = t    + 0.5 * dt
      x3(:) = x(:) + 0.5 * l2(:)
      p3(:) = p(:) + 0.5 * k2(:)

! calculate the Lorentz factor
!
      gamma = lorentz_factor(p3(:))

! calculate the velocity
!
      v3(:) = p3(:) / gamma

! calculate acceleration for the location x3 and velocity v3
!
      call acceleration(t3, x3, v3, a, u, b)

! calculate the third term
!
      l3(:) = dt  * v3(:)
      k3(:) = dtq * a (:)

!! 4th step of the RK integration
!!
! integrate the position and momentum
!
      t4    = t    + 0.5 * dt
      x4(:) = x(:) + l3(:)
      p4(:) = p(:) + k3(:)

! calculate the Lorentz factor
!
      gamma = lorentz_factor(p4(:))

! calculate the velocity
!
      v4(:) = p4(:) / gamma

! calculate acceleration for the location x4 and velocity v4
!
      call acceleration(t4, x4, v4, a, u, b)

! calculate the third term
!
      l4(:) = dt  * v4(:)
      k4(:) = dtq * a (:)

!! the final integration of the particle position and momentum
!!
      t5    = t    + dt
      x5(:) = x(:) + ( l1(:) + 2.0 * ( l2(:) + l3(:) ) + l4(:) ) / 6.0
      p5(:) = p(:) + ( k1(:) + 2.0 * ( k2(:) + k3(:) ) + k4(:) ) / 6.0

! calculate the Lorentz factor
!
      gamma = lorentz_factor(p5(:))

! calculate the velocity
!
      v5(:) = p5(:) / gamma

! calculate acceleration at the location x
!
      call acceleration(t5, x5, v5, a, u, b)

! estimate error
!
      l4(:) = l4(:) - dt  * v5(:)
      k4(:) = k4(:) - dtq * a (:)

      delta = sqrt(dot_product(l4, l4) + dot_product(k4, k4)) / 6.0

! estimate new timestep
!
      dtn   = dt * (rho * maxtol / delta)**0.2

      if (delta .gt. maxtol) then

! repeat integration with this timestep
!
        dt  = dtn
        dtq = qom * dt

      else

! update time
!
        t   = t + dt

! copy data to array
!
        do while (tt(n) .le. t .and. t .ge. tmin .and. t .lt. tmax)

! calculate the left and right weights
!
          wl = (t - tt(n)) / dt
          wr = 1.0d0 - wl

! interpolate the particle state at the proper time
!
          tp    = wl * (t - dt) + wr * t
          xt(:) = wl * x(:) + wr * x5(:)
          vt(:) = wl * v(:) + wr * v5(:)
          pt(:) = wl * p(:) + wr * p5(:)

! calculate acceleration for the location x4 and velocity v4
!
          call acceleration(tp, xt, vt, a, u, b)

! separate particle velocity into parallel and perpendicular components
!
          call separate_velocity(vt, b, ba, vu, vp, vr)

! calculate the particle gyroperiod and gyroradius
!
          call gyro_parameters(gamma, ba, vr, om, tg, rg)

! calculate particle energy
!
#ifdef RELAT
          en = gamma * mrest
          ek = en - mrest
#else
          en = 0.5 * vu**2
          ek = en
#endif

! write the progress
!
          write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, tp, dt, vu / c, ek, char(13)

! write results to the output file
!
          open  (10, file = 'output.dat', form = 'formatted', position = 'append')
          write (10, "(19(1pe18.10))") tp, xt(1), xt(2), xt(3), vt(1), vt(2), vt(3)  &
                                     , vu / c, vp / c, vr / c, gamma, en, ek   &
                                     , bavg * ba, om, tg * fc, rg * ln, tg, rg
          close (10)

          n = n + 1

        end do

! update position, velocity and momentum
!
        x(:) = x5(:)
        v(:) = v5(:)
        p(:) = p5(:)

! update new timestep
!
        dt  = min(2.0 * dt, dtn, dtmax, max(1.0e-16, tmax - t))
        dtq = qom * dt

      end if

    end do

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(v, b, ba, vu, vp, vr)

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gamma, ba, vr, om, tg, rg)

! calculate particle energy
!
#ifdef RELAT
    en = gamma * mrest
    ek = en - mrest
#else
    en = 0.5 * vu**2
    ek = en
#endif

! write the progress
!
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6))") n, t, dt, vu / c, ek

! write results to the output file
!
    open  (10, file = 'output.dat', form = 'formatted', position = 'append')
    write (10, "(19(1pe18.10))") t, x(1), x(2), x(3), v(1), v(2), v(3)         &
                               , vu / c, vp / c, vr / c, gamma, en, ek         &
                               , bavg * ba, om, tg * fc, rg * ln, tg, rg
    close (10)
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_rk4_log
!
!===============================================================================
!
! integrate_trajectory_si4: subroutine integrates particle trajectory using
!                           the 4th order simplectic method
!
! references: Mackey, Marchand & Kabin, 2006, JGR, 111, A06208
!             Calvo, Laburta & Montijano, 2003, CaMwA, 25, 401
!
!===============================================================================
!
  subroutine integrate_trajectory_si4()

    use params, only : dtini, tmax, c, ndumps

    implicit none

! local variables
!
    integer                         :: n, m
    real(kind=PREC), dimension(2,6) :: y, y0, y1, y2
    real(kind=PREC), dimension(3)   :: x , u , p , a
    real(kind=PREC), dimension(3)   :: x1, u1, p1, a1
    real(kind=PREC), dimension(3)   :: x2, u2, p2, a2
    real(kind=PREC)                 :: gm, g1, g2
    real(kind=8   ), dimension(3)   :: v, b
    real(kind=8   )                 :: t, dt, ds
    real(kind=8   )                 :: en, ek, ua, ba, up, ur, om, tg, rg

! local flags
!
    logical :: flag = .true.

! local parameters
!
    real(kind=8), parameter :: ac1  =    4.82308546376020871664d0              &
                             , ac2  =   67.17691453623979128336d0
    real(kind=8), parameter :: bc1  =    0.21132486540518711775d0              &
                             , bc2  =    0.78867513459481288225d0
    real(kind=8), parameter :: bc11 = -  8.36344801571300286532d0              &
                             , bc12 =    7.48780366869514391996d0              &
                             , bc21 = -115.48780366869514391996d0              &
                             , bc22 =   90.36344801571300286532d0
    real(kind=8), parameter :: gc11 = -  5.79422863405994782084d0              &
                             , gc12 =    2.84678751731759804956d0              &
                             , gc21 = - 50.84678751731759804956d0              &
                             , gc22 =    9.79422863405994782084d0
    real(kind=8), parameter :: ec   =    1.73205080756887729353d0
!
!-------------------------------------------------------------------------------
!
! initialize the iteration number, snapshot number, time, and time steps
!
    n  = 1
    m  = 0
    t  = 0.0d0
    dt = dtini
    ds = qom * dt

! substitute the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = v0(:)
    p(:) = p0(:)

! calculate the acceleration at the starting point
!
    call acceleration(t, x, u, a, v, b)

! calculate the particle speed
!
    ua = sqrt(sum(u(:)**2))

! calculate the Lorentz factor
!
    gm = lorentz_factor(p)

! calculate the particle energy
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else /* RELAT */
    en = 0.5 * ua * ua
    ek = en
#endif /* RELAT */

! print the progress information
!
    write (*,"('PROGRESS  : ',a8,2x,4(a14))") 'ITER', 'TIME', 'TIMESTEP'       &
            , 'SPEED (c)', 'ENERGY (MeV)'
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua, ek        &
            , char(13)

!== INTEGRATION LOOP ==
!
! iterate until the maximum time is reached
!
    do while (t .le. tmax)

! update the previous estimates of Y1 and Y2
!
      y2 = y1
      y1 = y

      y0(1,1:3) = x(:)
      y0(1,4:6) = p(:)
      y0(2,1:3) = x(:)
      y0(2,4:6) = p(:)

! find the initial guess for the states Y1 and Y2
!
      if (flag) then

        y(1,1:3) = x(:) + bc1 * dt * u(:)
        y(1,4:6) = p(:) + bc1 * ds * a(:)
        y(2,1:3) = x(:) + bc2 * dt * u(:)
        y(2,4:6) = p(:) + bc2 * ds * a(:)

      else

        y(1,1:6) = ac1 * y0(1,1:6) + bc11 * y2(1,1:6) + bc12 * y2(2,1:6)       &
                                   + gc11 * y1(1,1:6) + gc12 * y1(2,1:6)
        y(2,1:6) = ac2 * y0(2,1:6) + bc21 * y2(1,1:6) + bc22 * y2(2,1:6)       &
                                   + gc21 * y1(1,1:6) + gc22 * y1(2,1:6)

        flag = .false.
      end if

! estimate the states Y1 and Y2 iteratively
!
!   Y1 = y(n) + h * [ a11 * F(Y1) + a12 * F(Y2) ]
!   Y2 = y(n) + h * [ a21 * F(Y1) + a22 * F(Y2) ]
!
      call estimate(x, p, y, t, dt, ds)

! obtain the positions, momenta, gammas, and velocities for the estimated states
!
      x1(:) = y(1,1:3)
      x2(:) = y(2,1:3)
      p1(:) = y(1,4:6)
      p2(:) = y(2,4:6)
      g1    = lorentz_factor(p1(:))
      g2    = lorentz_factor(p2(:))
      u1(:) = p1(:) / g1
      u2(:) = p2(:) / g2

! calculate the acceleration at the estimated states
!
      call acceleration(t, x1(1:3), u1(1:3), a1(1:3), v, b)
      call acceleration(t, x2(1:3), u2(1:3), a2(1:3), v, b)

! update the solution
!
!   y(n+1) = y(n) + h * [ b1 * F(Y1) + b2 * F(Y2) ]
!
      x(:) = x(:) + 0.5d0 * dt * (u1(1:3) + u2(1:3))
      p(:) = p(:) + 0.5d0 * ds * (a1(1:3) + a2(1:3))

! update the integration time
!
      t = t + dt

! store the particle parameters at a given snapshot
!
      if (m .eq. ndumps) then

! calculate the Lorentz factor and particle velocity
!
        gm   = lorentz_factor(p(:))
        u(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(1:3), u(1:3), a(1:3), v, b)

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(u, b, ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
#ifdef RELAT
        en = gm * mrest
        ek = en - mrest
#else
        en = 0.5 * ua * ua
        ek = en
#endif

! write the progress
!
        write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua / c, ek, char(13)

! write results to the output file
!
        open  (10, file = 'output.dat', form = 'formatted', position = 'append')
        write (10, "(19(1pe18.10))") t, x(1), x(2), x(3), u(1), u(2), u(3)     &
                                   , ua / c, up / c, ur / c, gm, en, ek        &
                                   , bavg * ba, om, tg * fc, rg * ln, tg, rg
        close (10)

        n = n + 1
        m = 0

      end if

! increase data write counter
!
      m = m + 1

! end of iteration
!
    end do

! write the progress
!
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6))") n, t, dt, ua / c, ek
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_si4
!
!===============================================================================
!
! estimate: subroutine estimates the solution for the equation of motion using
!           a simple functional iteration
!
!===============================================================================
!
  subroutine estimate(x, p, y, t, dt, ds)

    use params, only : maxeps, maxtol, maxit, dtmax

    implicit none

! subroutine arguments
!
    real(kind=PREC), dimension(3)  , intent(inout) :: x, p
    real(kind=PREC), dimension(2,6), intent(inout) :: y
    real(kind=PREC)                , intent(inout) :: t, dt, ds

! local variables
!
    integer                         :: it
    real(kind=PREC), dimension(2,6) :: yn
    real(kind=PREC), dimension(6)   :: dh
    real(kind=PREC), dimension(3)   :: x1, p1, u1, a1
    real(kind=PREC), dimension(3)   :: x2, p2, u2, a2
    real(kind=PREC), dimension(3)   :: v, b
    real(kind=PREC)                 :: g1, g2, eps, tol
    real(kind=PREC)                 :: u1a, u2a, a1a, a2a

! local parameter
!
    real(kind=8), parameter :: a11 = 0.25d0, a12 = -0.03867513459481288225d0   &
                             , a22 = 0.25d0, a21 =  0.53867513459481288225d0
    real(kind=8), parameter :: ec  = 1.73205080756887729353d0
!
!-------------------------------------------------------------------------------
!
! initiate the iteration control parameters
!
    it  = 1
    eps = 1.0e+16

! perform the simple functional iteration until the conditions are met
!
    do while ((eps .gt. maxeps .or. abs(maxtol / tol - 1.0d0) .gt. 1.0e-8) .and. it .lt. maxit)

! prepare the initial particle position and momentum
!
      x1(:) = y(1,1:3)
      x2(:) = y(2,1:3)
      p1(:) = y(1,4:6)
      p2(:) = y(2,4:6)

! calculate the Lorentz factors and particle velocity
!
      g1    = lorentz_factor(p1(:))
      g2    = lorentz_factor(p2(:))
      u1(:) = p1(:) / g1
      u2(:) = p2(:) / g2

! calculate the accelerations
!
      call acceleration(t, x1(1:3), u1(1:3), a1(1:3), v, b)
      call acceleration(t, x2(1:3), u2(1:3), a2(1:3), v, b)

! calculate the error and estimate the new time step
!
      u1a    = max(1.0e-16, sqrt(sum(u1(:) * u1(:))))
      u2a    = max(1.0e-16, sqrt(sum(u2(:) * u2(:))))
      a1a    = max(1.0e-16, sqrt(sum(a1(:) * a1(:))))
      a2a    = max(1.0e-16, sqrt(sum(a2(:) * a2(:))))
      dh(1:3) = u2(:) / u2a - u1(:) / u1a
      dh(4:6) = a2(:) / a2a - a1(:) / a1a
      tol     = ec * max(1.0e-16, sqrt(sum(dh(1:3) * dh(1:3)))                 &
                                , sqrt(sum(dh(4:6) * dh(4:6))))
      dt      = min(dt * sqrt(maxtol / tol), dtmax)
      ds      = dt * qom

! update the positions and momenta
!
      yn(1,1:3) = x(1:3) + dt * (a11 * u1(1:3) + a12 * u2(1:3))
      yn(1,4:6) = p(1:3) + ds * (a11 * a1(1:3) + a12 * a2(1:3))
      yn(2,1:3) = x(1:3) + dt * (a21 * u1(1:3) + a22 * u2(1:3))
      yn(2,4:6) = p(1:3) + ds * (a21 * a1(1:3) + a22 * a2(1:3))

! calculate the maximum of residuum
!
      eps = maxval(abs(yn - y))

! substitute the new solution
!
      y = yn

! increase the iteration counter
!
      it = it + 1

    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine estimate
!
!===============================================================================
!
! pos2index: subroutine converts a given position to the array index
!
!===============================================================================
!
  subroutine pos2index(x, r)

    use fields, only : ng

    implicit none

! input and output arguments
!
    real(kind=PREC), dimension(3), intent(in)  :: x
    real(kind=PREC), dimension(3), intent(out) :: r

! local variables
!
    integer         :: i
    real(kind=PREC) :: t
!
!------------------------------------------------------------------------------
!
    do i = 1, DIMS
      t    = (x(i) - bnds(i,1)) / bsiz(i) + 0.5 / dm(i)
      t    = t - floor(t)
      r(i) = dm(i) * t + ng
    end do
!
!------------------------------------------------------------------------------
!
  end subroutine pos2index
!
!===============================================================================
!
! prepare_interpolation: subroutine prepares coeafficients for interpolation
!
!===============================================================================
!
  subroutine prepare_interpolation(x, ii, jj, kk, dr, cx, cy, cz)

    implicit none

! input and output arguments
!
    real(kind=PREC), dimension(3), intent(in)  :: x
    integer        , dimension(4), intent(out) :: ii, jj, kk
    real(kind=8   ), dimension(3), intent(out) :: dr
    real(kind=8   ), dimension(4), intent(out) :: cx, cy, cz
!
!------------------------------------------------------------------------------
!
#if !defined TRILIN && !defined TRICUB
!
!= nearest interpolation =
!
! calculate indices
!
    ii(1) = nint(x(1))
    jj(1) = nint(x(2))
#if DIMS == 3
    kk(1) = nint(x(3))
#endif /* DIMS == 3 */
#endif /* !TRILIN & !TRICUB */
#ifdef TRILIN
!
!= trilinear interpolation =
!
! calculate indices
!
    ii(1) = floor(x(1))
    ii(2) = ii(1) + 1
    jj(1) = floor(x(2))
    jj(2) = jj(1) + 1
#if DIMS == 3
    kk(1) = floor(x(3))
    kk(2) = kk(1) + 1
#endif /* DIMS == 3 */

! calculate intercell position
!
    dr(1) = x(1) - ii(1)
    dr(2) = x(2) - jj(1)
#if DIMS == 3
    dr(3) = x(3) - kk(1)
#endif /* DIMS == 3 */
#endif /* TRILIN */
#ifdef TRICUB
!
!= tricubic interpolation =
!
! calculate indices
!
    ii(2) = floor(x(1))
    ii(1) = ii(2) - 1
    ii(3) = ii(2) + 1
    ii(4) = ii(2) + 2
    jj(2) = floor(x(2))
    jj(1) = jj(2) - 1
    jj(3) = jj(2) + 1
    jj(4) = jj(2) + 2
#if DIMS == 3
    kk(2) = floor(x(3))
    kk(1) = kk(2) - 1
    kk(3) = kk(2) + 1
    kk(4) = kk(2) + 2
#endif /* DIMS == 3 */

! calculate intercell position
!
    dr(1) = x(1) - ii(2)
    dr(2) = x(2) - jj(2)
#if DIMS == 3
    dr(3) = x(3) - kk(2)
#endif /* DIMS == 3 */

! coefficients for dx, dy, and dz
!
    call coefficients_cubic(dr(1), cx(:))
    call coefficients_cubic(dr(2), cy(:))
#if DIMS == 3
    call coefficients_cubic(dr(3), cz(:))
#endif /* DIMS == 3 */
#endif /* TRICUB */
!
!------------------------------------------------------------------------------
!
  end subroutine prepare_interpolation
!
!===============================================================================
!
! interpolate: subroutine interpolates field value for a given position
!
!===============================================================================
!
  real(kind=8) function interpolate(f, ii, jj, kk, dr, cx, cy, cz) result(q)

    implicit none

! input and output arguments
!
    real(kind=4), dimension(:,:,:), intent(in) :: f
    integer     , dimension(4)    , intent(in) :: ii, jj, kk
    real(kind=8), dimension(3)    , intent(in) :: dr
    real(kind=8), dimension(4)    , intent(in) :: cx, cy, cz

! local variables
!
#ifdef TRILIN
    real(kind=4) :: q11, q12, q21, q22, q1, q2
#endif /* TRILIN */
#ifdef TRICUB
    real(kind=4) :: q11, q12, q13, q14, q21, q22, q23, q24                     &
                  , q31, q32, q33, q34, q41, q42, q43, q44, q1, q2, q3, q4
#endif /* TRICUB */
!
!------------------------------------------------------------------------------
!
#if !defined TRILIN && !defined TRICUB
!
!= nearest interpolation =
!
#if DIMS == 2
    q = f(ii(1),jj(1),1)
#else /* DIMS == 2 */
    q = f(ii(1),jj(1),kk(1))
#endif /* DIMS == 2 */
#endif /* !TRILIN & !TRICUB */
#ifdef TRILIN
#if DIMS == 2
!= bilinear interpolation =
!
! interpolate along the Y direction
!
    q1 = plinear(dr(2), f(ii(1),jj(1),1), f(ii(1),jj(2),1))
    q2 = plinear(dr(2), f(ii(2),jj(1),1), f(ii(2),jj(2),1))
#else /* DIMS == 2 */
!= trilinear interpolation =
!
! interpolate along the Z direction
!
    q11 = plinear(dr(3), f(ii(1),jj(1),kk(1)), f(ii(1),jj(1),kk(2)))
    q12 = plinear(dr(3), f(ii(1),jj(2),kk(1)), f(ii(1),jj(2),kk(2)))
    q21 = plinear(dr(3), f(ii(2),jj(1),kk(1)), f(ii(2),jj(1),kk(2)))
    q22 = plinear(dr(3), f(ii(2),jj(2),kk(1)), f(ii(2),jj(2),kk(2)))

! interpolate along the Y direction
!
    q1 = plinear(dr(2), q11, q12)
    q2 = plinear(dr(2), q21, q22)
#endif /* DIMS == 2 */

! interpolate the value at a given position
!
    q  = plinear(dr(1), q1 , q2 )
#endif /* TRILIN */
#ifdef TRICUB
#if DIMS == 2
!= bicubic interpolation =
!
! interpolate along the Y direction
!
    q1  = pcubic(dr(2), cy, f(ii(1),jj(1),1), f(ii(1),jj(2),1)                 &
                          , f(ii(1),jj(3),1), f(ii(1),jj(4),1))
    q2  = pcubic(dr(2), cy, f(ii(2),jj(1),1), f(ii(2),jj(2),1)                 &
                          , f(ii(2),jj(3),1), f(ii(2),jj(4),1))
    q3  = pcubic(dr(2), cy, f(ii(3),jj(1),1), f(ii(3),jj(2),1)                 &
                          , f(ii(3),jj(3),1), f(ii(3),jj(4),1))
    q4  = pcubic(dr(2), cy, f(ii(4),jj(1),1), f(ii(4),jj(2),1)                 &
                          , f(ii(4),jj(3),1), f(ii(4),jj(4),1))
#else /* DIMS == 2 */
!= tricubic interpolation =
!
! interpolate along Z direction
!
    q11 = pcubic(dr(3), cz, f(ii(1),jj(1),kk(1)), f(ii(1),jj(1),kk(2))         &
                          , f(ii(1),jj(1),kk(3)), f(ii(1),jj(1),kk(4)))
    q12 = pcubic(dr(3), cz, f(ii(1),jj(2),kk(1)), f(ii(1),jj(2),kk(2))         &
                          , f(ii(1),jj(2),kk(3)), f(ii(1),jj(2),kk(4)))
    q13 = pcubic(dr(3), cz, f(ii(1),jj(3),kk(1)), f(ii(1),jj(3),kk(2))         &
                          , f(ii(1),jj(3),kk(3)), f(ii(1),jj(3),kk(4)))
    q14 = pcubic(dr(3), cz, f(ii(1),jj(4),kk(1)), f(ii(1),jj(4),kk(2))         &
                          , f(ii(1),jj(4),kk(3)), f(ii(1),jj(4),kk(4)))

    q21 = pcubic(dr(3), cz, f(ii(2),jj(1),kk(1)), f(ii(2),jj(1),kk(2))         &
                          , f(ii(2),jj(1),kk(3)), f(ii(2),jj(1),kk(4)))
    q22 = pcubic(dr(3), cz, f(ii(2),jj(2),kk(1)), f(ii(2),jj(2),kk(2))         &
                          , f(ii(2),jj(2),kk(3)), f(ii(2),jj(2),kk(4)))
    q23 = pcubic(dr(3), cz, f(ii(2),jj(3),kk(1)), f(ii(2),jj(3),kk(2))         &
                          , f(ii(2),jj(3),kk(3)), f(ii(2),jj(3),kk(4)))
    q24 = pcubic(dr(3), cz, f(ii(2),jj(4),kk(1)), f(ii(2),jj(4),kk(2))         &
                          , f(ii(2),jj(4),kk(3)), f(ii(2),jj(4),kk(4)))

    q31 = pcubic(dr(3), cz, f(ii(3),jj(1),kk(1)), f(ii(3),jj(1),kk(2))         &
                          , f(ii(3),jj(1),kk(3)), f(ii(3),jj(1),kk(4)))
    q32 = pcubic(dr(3), cz, f(ii(3),jj(2),kk(1)), f(ii(3),jj(2),kk(2))         &
                          , f(ii(3),jj(2),kk(3)), f(ii(3),jj(2),kk(4)))
    q33 = pcubic(dr(3), cz, f(ii(3),jj(3),kk(1)), f(ii(3),jj(3),kk(2))         &
                          , f(ii(3),jj(3),kk(3)), f(ii(3),jj(3),kk(4)))
    q34 = pcubic(dr(3), cz, f(ii(3),jj(4),kk(1)), f(ii(3),jj(4),kk(2))         &
                          , f(ii(3),jj(4),kk(3)), f(ii(3),jj(4),kk(4)))

    q41 = pcubic(dr(3), cz, f(ii(4),jj(1),kk(1)), f(ii(4),jj(1),kk(2))         &
                          , f(ii(4),jj(1),kk(3)), f(ii(4),jj(1),kk(4)))
    q42 = pcubic(dr(3), cz, f(ii(4),jj(2),kk(1)), f(ii(4),jj(2),kk(2))         &
                          , f(ii(4),jj(2),kk(3)), f(ii(4),jj(2),kk(4)))
    q43 = pcubic(dr(3), cz, f(ii(4),jj(3),kk(1)), f(ii(4),jj(3),kk(2))         &
                          , f(ii(4),jj(3),kk(3)), f(ii(4),jj(3),kk(4)))
    q44 = pcubic(dr(3), cz, f(ii(4),jj(4),kk(1)), f(ii(4),jj(4),kk(2))         &
                          , f(ii(4),jj(4),kk(3)), f(ii(4),jj(4),kk(4)))

! interpolate along the Y direction
!
    q1  = pcubic(dr(2), cy, q11, q12, q13, q14)
    q2  = pcubic(dr(2), cy, q21, q22, q23, q24)
    q3  = pcubic(dr(2), cy, q31, q32, q33, q34)
    q4  = pcubic(dr(2), cy, q41, q42, q43, q44)
#endif /* DIMS == 2 */

! interpolate along the X direction
!
    q   = pcubic(dr(1), cx, q1 , q2 , q3 , q4 )
#endif /* TRICUB */
!
!------------------------------------------------------------------------------
!
  end function interpolate
!
!===============================================================================
!
! plinear: subroutine performs one dimensional linear interpolation
!
!===============================================================================
!
  real(kind=8) function plinear(x, fl, fr) result(q)

    implicit none

! input and output arguments
!
    real(kind=8), intent(in)  :: x
    real(kind=4), intent(in)  :: fl, fr
!
!------------------------------------------------------------------------------
!
    q = fl + x * (fr - fl)
!
!------------------------------------------------------------------------------
!
  end function plinear
!
!===============================================================================
!
! pcubic: subroutine performs one dimensional cubic interpolation
!
!===============================================================================
!
  real(kind=8) function pcubic(x, c, fk, fl, fr, fq) result(q)

    implicit none

! input and output arguments
!
    real(kind=8)              , intent(in)  :: x
    real(kind=8), dimension(4), intent(in)  :: c
    real(kind=4)              , intent(in)  :: fk, fl, fr, fq

#ifdef TVD
! local parameters
!
    real(kind=4) :: dfl, dfr, ds, dl, dr
#endif /* TVD */
!
!------------------------------------------------------------------------------
!
#ifdef TVD
    ds  = (fr - fl)
    dl  = (fl - fk)
    dr  = (fq - fr)
    dfl = sign(1.0, ds) * min(abs(ds), abs(dl))
    dfr = sign(1.0, ds) * min(abs(ds), abs(dr))

    if ((dl * ds) .le. 0.0) then
      dfl = 0.0
    end if

    if ((dr * ds) .le. 0.0) then
      dfr = 0.0
    end if

    q = c(1) * fl + c(2) * fr + c(3) * dfl + c(4) * dfr
#else /* TVD */
    q = c(1) * fk + c(2) * fl + c(3) * fr + c(4) * fq
#endif /* TVD */
!
!------------------------------------------------------------------------------
!
  end function pcubic
!
!===============================================================================
!
! coefficients_cubic: subroutine prepares coeafficients for cubic interpolation
!
!===============================================================================
!
  subroutine coefficients_cubic(x, c)

    implicit none

! input and output arguments
!
    real(kind=8)              , intent(in)  :: x
    real(kind=8), dimension(4), intent(out) :: c

! local variables
!
    real(kind=8) :: x1, x2, x3, xd
!
!------------------------------------------------------------------------------
!
! prepare local variables
!
    x1 = x - 1.0
    x2 = x * x
#ifdef TVD
    x3 = x1 * x1
    xd = 2.0 * x
#else /* TVD */
#endif /* TVD */

! calculate coefficients
!
#ifdef TVD
    c(1) = x3 * (xd + 1.0)
    c(2) = x2 * (3.0 - xd)
    c(3) = x  * x3
    c(4) = x2 * x1
#else /* TVD */
    c(1) = 0.5 * x * ( ( 2.0 - x ) * x - 1.0 )
    c(2) = 0.5 * x2 * ( 3.0 * x - 5.0 ) + 1.0
    c(3) = 0.5 * x * ( ( 4.0 - 3.0 * x ) * x + 1.0 )
    c(4) = 0.5 * x1 * x2
#endif /* TVD */
!
!------------------------------------------------------------------------------
!
  end subroutine coefficients_cubic
!
!===============================================================================
!
! lorentz_factor: subroutine calculates the Lorentz factor
!
!===============================================================================
!
  real(kind=PREC) function lorentz_factor(p) result(gm)

    implicit none

! input and output arguments
!
    real(kind=PREC), dimension(3), intent(in)  :: p
!
!------------------------------------------------------------------------------
!
#ifdef RELAT
    gm = sqrt(1.0d0 + dot_product(p, p) / c2)
#else
    gm = 1.0
#endif
!
!-------------------------------------------------------------------------------
!
  end function lorentz_factor
!
!===============================================================================
!
! acceleration: subroutine calculates acceleration vector at a given location
!
!===============================================================================
!
  subroutine acceleration(t, x, v, a, u, b)

    use fields, only : ux, uy, uz, bx, by, bz
#ifdef TEST
    use params, only : omega, bini, bamp, vamp, freq, epar
#endif /* TEST */

    implicit none

! input and output arguments
!
    real(kind=PREC)              , intent(in)  :: t
    real(kind=PREC), dimension(3), intent(in)  :: x, v
    real(kind=PREC), dimension(3), intent(out) :: a, u, b

#ifdef TEST
! local variables
!
    real(kind=PREC), dimension(3) :: w
#ifdef ITEST
    real(kind=PREC)               :: ra, rb, xt, yt, rt
#endif /* ITEST */

! parameters
!
    real(kind=8   ) :: pi2 = 6.2831853071795862319959269370884d0
#else /* TEST */
! local variables
!
    real(kind=PREC), dimension(3) :: r, w

! position indices
!
    integer        , dimension(4) :: ii, jj, kk
    real(kind=8   ), dimension(4) :: cx, cy, cz
    real(kind=8   ), dimension(3) :: dr
#endif /* TEST */
!
!------------------------------------------------------------------------------
!
#ifndef TEST
! convert position to index
!
      call pos2index(x, r)

#ifdef BNDRY
      if (minval(r) .gt. 4 .or. minval(dm - r) .gt. 4) then
#endif /* BNDRY */

! prepare coefficients for interpolation
!
        call prepare_interpolation(r, ii, jj, kk, dr, cx, cy, cz)
#endif /* !TEST */

! interpolate field components at the particle position
!
#ifdef TEST
#ifdef WTEST
        u(1) = 0.0
        u(2) = - vamp * sin(pi2 * freq * x(1))
        u(3) = 0.0

        b(1) = bpar
        b(2) = bamp * cos(pi2 * freq * x(1))
        b(3) = bamp * sin(pi2 * freq * x(1))
#endif /* WTEST */

#ifdef ITEST
!         ra   = 1.0 + amp * sin(pi2 * omega * t)
!         rb   = amp * pi2 * omega * cos(pi2 * omega * t)
        ra   = 1.0 + bamp

        rt   = sqrt(x(1) * x(1) + x(2) * x(2))

        u(1) = - x(1)
        u(2) =   x(2)
        u(3) = 0.0
!         if (rt .gt. 0.0) then
!           u(1) = - x(1) * rb / ra**2
!           u(2) =   x(2) * rb
!           u(3) = 0.0
!         else
!           u(:) = 0.0
!         end if

        xt   = x(1) / ra
        yt   = x(2) * ra

        rt   = sqrt(xt * xt + yt * yt)

        if (rt .gt. 0.0) then
          b(1) = - yt / rt
          b(2) =   xt / rt
          b(3) = 0.0
        else
          b(:) = 0.0
        end if
#endif /* ITEST */
#else /* TEST */
        u(1) = interpolate(ux, ii, jj, kk, dr, cx, cy, cz)
        u(2) = interpolate(uy, ii, jj, kk, dr, cx, cy, cz)
        u(3) = interpolate(uz, ii, jj, kk, dr, cx, cy, cz)
        b(1) = interpolate(bx, ii, jj, kk, dr, cx, cy, cz)
        b(2) = interpolate(by, ii, jj, kk, dr, cx, cy, cz)
        b(3) = interpolate(bz, ii, jj, kk, dr, cx, cy, cz)
#endif /* TEST */

! subtract the fluid velocity
!
        w(:) = v(:) - u(:)

! compute the acceleration
!
        a(1) = w(2) * b(3) - w(3) * b(2)
        a(2) = w(3) * b(1) - w(1) * b(3)
        a(3) = w(1) * b(2) - w(2) * b(1)
#ifdef TEST
        a(1) = a(1) + epar
#endif /* TEST */
#ifndef TEST
#ifdef BNDRY
      else
        a(:) = 0.0
      endif
#endif /* BNDRY */
#endif /* !TEST */
!
!-------------------------------------------------------------------------------
!
  end subroutine acceleration
!
!===============================================================================
!
! separate_velocity: subroutine separates velocity into two components, parallel
!                    and perpendicular to the local magnetic field
!
!===============================================================================
!
  subroutine separate_velocity(v, b, ba, vu, vp, vr)

    implicit none

! input and output arguments
!
    real(kind=PREC), dimension(3), intent(in)  :: v, b
    real(kind=8   )              , intent(out) :: ba, vu, vp, vr

! local variables
!
    real(kind=8   ), dimension(3) :: p
    real(kind=8   )               :: pp, vv
!
!------------------------------------------------------------------------------
!
! calculate amplitude of magnetic field
!
    ba = sqrt(dot_product(b, b))

! calculate unit vector parallel to B
!
    if (ba .gt. 0.0) then
      p  = b / ba
    else
      p  = 0.0
    end if

! calculate component parallel to B
!
    pp = dot_product(v, p)**2

! calculate amplitude of velocity
!
    vv = dot_product(v, v)

! calculate amplitude of the parallel and perpendicular components of velocity
!
    vu = sqrt(vv)
    vp = sqrt(pp)
    vr = sqrt(vv - pp)
!
!-------------------------------------------------------------------------------
!
  end subroutine separate_velocity
!
!===============================================================================
!
! gyro_parameters: subroutine calculates particle gyrofrequency, gyroperiod and
!                  gyroradius
!
!===============================================================================
!
  subroutine gyro_parameters(gm, ba, vr, om, tg, rg)

    implicit none

! input and output arguments
!
    real(kind=PREC), intent(in)  :: gm
    real(kind=8   ), intent(in)  :: ba, vr
    real(kind=8   ), intent(out) :: om, tg, rg

! parameters
!
    real(kind=8   ) :: pi2 = 6.2831853071795862319959269370884d0
!
!------------------------------------------------------------------------------
!
    om = om0 * ba / gm
    tg = pi2 / om / fc
    rg = vr / om / fc
!
!-------------------------------------------------------------------------------
!
  end subroutine gyro_parameters
!
end module particles
