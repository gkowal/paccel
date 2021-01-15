!!******************************************************************************
!!
!! module: particles - subroutines to prepare and advance particles
!!
!! Copyright (C) 2010-2021 Grzegorz Kowal <grzegorz@gkowal.info>
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
  real(kind=8), dimension(3,2), save :: bnds

! domain size
!
  real(kind=8), dimension(3)  , save :: bsiz

! the units of time, length, velocity and magnetic field
!
  real(kind=8), save :: tunit = 1.0d+00, lunit = 1.0d+00
  real(kind=8), save :: vunit = 1.0d+00, bunit = 1.0d+00
  integer     , save :: nghost = 8 ! the number of ghost pixels near the boundary

! particle mass and the speed of light
!
  real(kind=8), save :: qom   = 9.5788332d+03
  real(kind=8), save :: mrest, bpar

#ifdef TEST
! test problem parameters
!
  real(kind=8), save :: bini = 1.0d+00 ! the mean magnetic field
  real(kind=8), save :: bamp = 0.0d+00 ! the amplitude of the magnetic field fluctuations
  real(kind=8), save :: vamp = 0.0d+00 ! the amplitude of the velocity field fluctuations
  real(kind=8), save :: freq = 1.0d+00 ! the frequency of the field fluctuations
#ifdef ITEST
  real(kind=8), save :: vrat = 1.0d+00 ! the ratio between velocity fluctuations amplitudes in different directions
  real(kind=8), save :: bshr = 0.0d+00 ! the guilde field
  real(kind=8), save :: epar = 0.0d+00 ! the constant electric field along the parallel direction
#endif /* ITEST */
#endif /* TEST */

! output data parameters
!
  integer     , save :: ndumps = 1000    ! number of steps between subsequent dumps
  real(kind=8), save :: tmin   = 1.0d-03 ! minimum time of writing data
  real(kind=8), save :: tmax   = 1.0d+00 ! maximum time for integration

! integration parameters
!
  real(kind=8), save :: safety = 5.0d-01  ! safety coefficient
  real(kind=8), save :: maxtol = 1.0d-04  ! the maximi integration tolerance
  real(kind=8), save :: maxeps = 1.0d-15  ! the maximum iteration error
  real(kind=8), save :: dtini  = 1.0d-08  ! the initial time step
  real(kind=8), save :: dtmax  = 1.0d+00  ! maximum allowed step size
  integer     , save :: maxit  = 1000     ! the limit of iterations


! arrays containing the initial positions and velocities of particle
!
  real(kind=8), dimension(3), save :: x0, u0, p0

! global parameters
!
  character(len=1), parameter :: term = char(13)
  real(kind=8)    , parameter :: pi2 = 6.2831853071795862319959269370884d+00
  real(kind=8)    , parameter :: c   = 2.99792458d+08 ! the speed of light [m/s]
!
!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_PARTICLES:
! -------------------------------
!
!   Subroutine initializes the module.
!
!   Arguments:
!
!     verbose - indicates if it should print any messages;
!     status  - the status value: 0 - success, otherwise there was a problem;
!
!===============================================================================
!
  subroutine initialize_particles(verbose, status)

! import required modules
!
    use fields    , only : ufac, bfac
    use parameters, only : get_parameter

    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: status

! particle parameters
!
    character(len= 1) :: ptype  = 'p'     ! particle type: 'p' - proton,
!                                                          'e' - electron
    character(len=16) :: stunit = 's'     ! time unit:     'u' - microsecond,
!                                                          's' - second,
!                                                          'm' - minute,
!                                                          'h' - hour
!                                                          'd' - dat, etc.
    real(kind=8)      :: tmulti = 1.0d+00 ! the factor multiplying time unit
!
!-------------------------------------------------------------------------------
!
    status = 0

! print info
!
    if (verbose) write(*,"('INFO',6x,': initializing particles')")

! get input parameters
!
    call get_parameter('vunit' , vunit)
    call get_parameter('bunit' , bunit)
    call get_parameter('tunit' , stunit)
    call get_parameter('tmulti', tmulti)
    call get_parameter('ptype' , ptype)
    call get_parameter('ndumps', ndumps)
    call get_parameter('tmin'  , tmin)
    call get_parameter('tmax'  , tmax)
    call get_parameter('safety', safety)
    call get_parameter('maxtol', maxtol)
    call get_parameter('maxeps', maxeps)
    call get_parameter('dtini' , dtini)
    call get_parameter('dtmax' , dtmax)
    call get_parameter('maxit' , maxit)

! get time unit in seconds
!
    select case(trim(stunit))
    case('u')
      tunit = 1.0d-06
    case('s')
      tunit = 1.0d+00
    case('m')
      tunit = 6.0d+01
    case('h')
      tunit = 3.6d+03
    case('d')
      tunit = 8.64d+04
    case('w')
      tunit = 6.048d+05
    case('y')
      tunit = 3.1557d+07
    case default
      tunit = 1.0d+00
    end select
    tunit = tmulti * tunit
    lunit = vunit * tunit
    vunit = vunit / c

! print geometry parameters
!
    write(*,"('INFO',6x,': geometry parameters:')"  )
    write(*,"('INFO',6x,': T     =',1es15.8,' [s]')") tunit
    write(*,"('INFO',6x,': L     =',1es15.8,' [m]')") lunit

! print plasma parameters
!
    write(*,"('INFO',6x,': plasma parameters:')")
    write(*,"('INFO',6x,': V     =',1es15.8,' [m/s] =',1es15.8,' [c]')")       &
                                                      vunit * c, vunit
    write(*,"('INFO',6x,': B     =',1es15.8,' [G]')") bunit

! initialize particle parameters
!
    select case(ptype)
    case ('e')
      mrest =  5.10998949998580864751d-01  ! rest energy of electron [MeV]
      qom   = -1.75882001076384559274d+07  ! e/m [1 / Gs s]
    case default
      mrest =  9.38272088161040869636d+02  ! rest energy of proton   [MeV]
      qom   =  9.57883315593801671639d+03  ! e/m [1 / Gs s]
    end select

! print particle parameters
!
    write(*,"('INFO',6x,': particle parameters:')")
    select case(ptype)
    case ('e')
      write(*,"('INFO',6x,': trajectory for electron')")
    case default
      write(*,"('INFO',6x,': trajectory for proton')")
    end select
    write(*,"('INFO',6x,': q/m   =',1pe15.8,' [1/G·s] =',1pe15.8,' [1/G·T]')" )&
                                                        qom, qom * tunit

! convert q/m to the desired time units
!
    qom   = qom * tunit

! set factors for the plasma velocity and magnetic field
!
    ufac  = vunit
    bfac  = bunit
!
!-------------------------------------------------------------------------------
!
  end subroutine initialize_particles
!
!===============================================================================
!
! subroutine GENERATE_PARTICLE:
! ----------------------------
!
!   Subroutine generate the particle initial state (position and moment).
!
!   Arguments:
!
!     verbose - indicates if it should print any messages;
!     status  - the status value: 0 - success, otherwise there was a problem;
!
!===============================================================================
!
  subroutine generate_particle(verbose, status)

! import required modules
!
    use fields    , only : ux, uy, uz, bx, by, bz
    use fields    , only : get_domain_dimensions, get_domain_bounds
    use parameters, only : get_parameter

    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: status

! local variables
!
    real(kind=8) :: vpar  = 0.0d+00  ! the parallel and
    real(kind=8) :: vper  = 1.0d-02  ! perpendicular component of
                                     ! the initial speed in c
    integer      :: p

! local variables
!
    real(kind=8) :: babs, vabs, lfac, grad, gper, gfrq, ener, ekin
#ifdef ITEST
    real(kind=8) :: xt, yt, rt, dl, ra, rb, ec
#endif /* ITEST */
    real(kind=8) :: bb, ub, uu, ww
    real(kind=8), dimension(3) :: r, b, u, w

! position indices
!
    integer     , dimension(4) :: ii, jj, kk
    real(kind=8), dimension(4) :: cx, cy, cz
    real(kind=8), dimension(3) :: dr
!
!-------------------------------------------------------------------------------
!
    status = 0

! print info
!
    if (verbose) write(*,"('INFO',6x,': initializing particle initial state')")

! get parameters
!
    call get_parameter('xc'    , x0(1))
    call get_parameter('yc'    , x0(2))
    call get_parameter('zc'    , x0(3))
    call get_parameter('vpar'  , vpar)
    call get_parameter('vper'  , vper)
#ifdef TEST
    call get_parameter('bini'  , bini)
    call get_parameter('bamp'  , bamp)
    call get_parameter('vamp'  , vamp)
    call get_parameter('freq'  , freq)
#ifdef ITEST
    call get_parameter('bshr'  , bshr)
    call get_parameter('vrat'  , vrat)
    call get_parameter('epar'  , epar)
#endif /* ITEST */
#endif /* TEST */

! calculate the particle initial parameters
!
    vabs = sqrt(vpar**2 + vper**2)
    lfac = 1.0d+00 / sqrt(1.0d+00 - vabs**2)
    grad = lfac * vper * c / (abs(qom) * bunit)
    gper = pi2 * grad / (vper * c)
    gfrq = abs(qom) * bunit

#ifdef TEST
#ifdef WTEST
    bpar = sqrt(bini**2 - bamp**2)
    b(1) = bpar
    b(2) = bamp * cos(pi2 * freq * x0(1))
    b(3) = bamp * sin(pi2 * freq * x0(1))

    u(1) = 0.0
    u(2) = 0.0
    u(3) = 1.0
#endif /* WTEST */
#ifdef ITEST
! print infor about the island problem
!
    write( *, "('PROBLEM   : motion in the contracting magnetic island:')" )

! check if the initial field is non zero
!
    if (bini .eq. 0.0d0) then
      write (*, "('ERROR     : bini must be not zero!')" )
      stop
    end if

! prepare parameters for the magnetic field topology
!
    dl = bamp / bini

    if (abs(dl) .ge. 1.0) then
      write (*, "('ERROR     : parameter bamp must be smaller than bini!')" )
      stop
    end if

    ra = 1.0d0 + dl
    rb = 1.0d0 - dl

! calculate the eccentricity
!
    ec = dsqrt(1.0d0 - (rb / ra)**2)

! print info about the problem
!
    write (*, "('INFO      : magnetic field strength          =',1pe15.8)") bini
    write (*, "('INFO      : guide field strength             =',1pe15.8)") bshr
    write (*, "('INFO      : magnetic field perturbation      =',1pe15.8)") bamp
    write (*, "('INFO      : magnetic island eccentricity     =',1pe15.8)") ec
    write (*, "('INFO      : horizontal velocity perturbation =',1pe15.8)") vamp
    write (*, "('INFO      : vertical velocity perturbation   =',1pe15.8)") vamp * vrat

! calculate the parallel and perpendicular directions with respect to the
! local magnetic field at the initial particle position
!
    xt   = x0(1) / ra
    yt   = x0(2) / rb

    rt   = dsqrt(xt * xt + yt * yt)

    if (rt .gt. 0.0d0) then

      b(1) =   yt / rb / rt
      b(2) = - xt / ra / rt
      b(3) = bshr

      bb = dsqrt(dot_product(b(:), b(:)))

      if (bb .gt. 0.0d0) then
        b(:) = b(:) / bb
      else
        write (*, "('ERROR     : zero magnetic field at the initial position!')" )
        stop
      end if

      u(1) = 0.0d0
      u(2) = 0.0d0
      u(3) = 1.0d0

      ub   = dot_product(u(:), b(:))

      u(1) = u(1) - b(1) * ub
      u(2) = u(2) - b(2) * ub
      u(3) = u(3) - b(3) * ub

      uu   = dsqrt(dot_product(u(:), u(:)))

      if (uu .gt. 0.0d0) then
        u(:) = u(:) / uu
      else
        write (*, "('ERROR     : cannot determine perpendicular direction!')" )
        stop
      end if
    else
      write (*, "('ERROR     : initial position cannot be located in the origin!')" )
      stop
    end if

#endif /* ITEST */
#else /* TEST */
! get domain dimensions
!
    call get_domain_dimensions(dm)

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

! convert position to index
!
    call pos2index(x0, r)

! prepare coefficients for interpolation
!
    call prepare_interpolation(r, ii, jj, kk, dr, cx, cy, cz)

! interpolate field components at the particle position
!
    b(1) = interpolate(bx, ii, jj, kk, dr, cx, cy, cz)
    b(2) = interpolate(by, ii, jj, kk, dr, cx, cy, cz)
    b(3) = interpolate(bz, ii, jj, kk, dr, cx, cy, cz)
#endif /* TEST */

! calculate the direction of the local magnetic field
!
    babs = sqrt(dot_product(b, b))
    if (babs > 0.0d+00) then
      b(:) = b(:) / babs
    else
      write( *, "('ERROR     : ',a)" ) "B=0 at the initial position! Choose another one."
      stop
    endif

! calculate a unit vector perpendicular to B
!
    if (dm(3) == 1) then
      w(1) = 0.0d+00
      w(2) = 0.0d+00
      w(3) = 1.0d+00
    else
      call random_number(w)
      w = w - 0.5d+00
    end if

    ww = sqrt(dot_product(w, w))
    if (ww > 0.0d+00) then
      w(:) = w(:) / ww
    else
      write( *, "('ERROR     : ',a)" ) "V=0 at the initial position! Choose another one."
      stop
    end if

    u(1) = w(2) * b(3) - w(3) * b(2)
    u(2) = w(3) * b(1) - w(1) * b(3)
    u(3) = w(1) * b(2) - w(2) * b(1)

    uu = sqrt(dot_product(u, u))
    if (uu > 0.0d+00) then
      u(:) = u(:) / uu
    else
      write( *, "('ERROR     : ',a)" ) "V=0 at the initial position! Choose another one."
      stop
    end if

! calculate the initial particle velocity vector
!
    u0(:) = (vpar * b(:) + vper * u(:))

! calculate the initial particle momentuum
!
    lfac = 1.0d+00 / sqrt(1.0d+00 - dot_product(u0, u0))
    p0(:) = lfac * u0(:)

! allow to set the particle moment explicitely
!
    call get_parameter('px', p0(1))
    call get_parameter('py', p0(2))
    call get_parameter('pz', p0(3))

    lfac = lorentz_factor(p0(:))
    u0(:) = p0(:) / lfac

! calculate particle energy
!
    ener = lfac * mrest
    ekin = (lfac - 1.0d+00) * mrest

! print the particle initial parameters
!
    write(*,"('INFO',6x,': Vpar  =',1pe15.8,' [c] =',1pe15.8,' [m / s]')") vpar, vpar * c
    write(*,"('INFO',6x,': Vper  =',1pe15.8,' [c] =',1pe15.8,' [m / s]')") vper, vper * c
    write(*,"('INFO',6x,': |V|   =',1pe15.8,' [c] =',1pe15.8,' [m / s]')") vabs, vabs * c
    write(*,"('INFO',6x,': gamma =',1pe15.8)"           ) lfac
    write(*,"('INFO',6x,': Om    =',1pe15.8,' [1 / s]')") gfrq
    write(*,"('INFO',6x,': Tg    =',1pe15.8,' [s]')"    ) gper
    write(*,"('INFO',6x,': Rg    =',1pe15.8,' [m]')"    ) grad
    write(*,"('INFO',6x,': E0    =',1pe15.8,' [MeV]')"  ) ener
    write(*,"('INFO',6x,': Ek    =',1pe15.8,' [MeV]')"  ) ekin
    write(*,"('INFO',6x,': Rg/L  =',1pe15.8)" ) grad / lunit
    write(*,"('INFO',6x,': Tg/T  =',1pe15.8)" ) gper / tunit
!
!-------------------------------------------------------------------------------
!
  end subroutine generate_particle
!
!===============================================================================
!
! integrate_trajectory_rk4: subroutine integrates particle trajectory using
!                           the 4th order RK method
!
!===============================================================================
!
  subroutine integrate_trajectory_rk4()

    implicit none

! local variables
!
    logical                    :: keepon = .true.
    integer                    :: n, m
    real(kind=8)               ::    t1, t2, t3, t4, t5
    real(kind=8), dimension(3) :: x, x1, x2, x3, x4, x5
    real(kind=8), dimension(3) :: u, u1, u2, u3, u4, u5
    real(kind=8), dimension(3) :: p, p1, p2, p3, p4, p5
    real(kind=8), dimension(3) ::    k1, k2, k3, k4, k5
    real(kind=8), dimension(3) ::    l1, l2, l3, l4, l5
    real(kind=8), dimension(3) :: s, a, v, b
    real(kind=8)               :: gm, t, dt, dtn
    real(kind=8)               :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=8)               :: tol = 0.0d+00
!
!-------------------------------------------------------------------------------
!
! initialize counters, time and timesteps
!
    n  = 0
    m  = 0
    t  = 0.0d0
    dt = dtini

! set the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the initial position
!
    call acceleration(t, x(:), u(:), s(:), a(:), v(:), b(:))

! separate the particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate the particle energies
!
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua, ek, term

! open the output file, print headers and the initial values
!
    open (10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,19a22)") 'Time', 'X', 'Y', 'Z', 'Vx', 'Vy', 'Vz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]',              &
                                 '<B> [Gs]', 'Omega [1/s]',                    &
                                 'Tg [s]', 'Rg [m]', 'Tg [T]', 'Rg [L]',       &
                                 'Tolerance'
    write(10,"(20(1es22.14))") t, x(1), x(2), x(3), u(1), u(2), u(3),          &
                               ua, up, ur, gm, en, ek,                         &
                               bunit * ba, om / tunit, tg * tunit, rg * lunit, &
                               tg, rg, tol

!== INTEGRATION LOOP ==
!
! integrate the trajectory
!
    do while (keepon)

!! 1st step of the RK integration
!!
! integrate the position and momentum
!
      t1    = t
      x1(:) = x(:)
      p1(:) = p(:)

! calculate the Lorentz factor
!
      gm = lorentz_factor(p1(:))

! calculate the velocity
!
      u1(:) = p1(:) / gm

! calculate the acceleration for the location x1 and velocity u1
!
      call acceleration(t1, x1(:), u1(:), s(:), a(:), v(:), b(:))

! calculate the first term
!
      l1(:) = dt * s(:)
      k1(:) = dt * a(:)

!! 2nd step of the RK integration
!!
! integrate the position and momentum
!
      t2    = t    + 0.5d0 * dt
      x2(:) = x(:) + 0.5d0 * l1(:)
      p2(:) = p(:) + 0.5d0 * k1(:)

! calculate the Lorentz factor
!
      gm = lorentz_factor(p2(:))

! calculate the velocity
!
      u2(:) = p2(:) / gm

! calculate the acceleration for the location x2 and velocity u2
!
      call acceleration(t2, x2(:), u2(:), s(:), a(:), v(:), b(:))

! calculate the second term
!
      l2(:) = dt * s(:)
      k2(:) = dt * a(:)

!! 3rd step of the RK integration
!!
! integrate the position and momentum
!
      t3    = t    + 0.5d0 * dt
      x3(:) = x(:) + 0.5d0 * l2(:)
      p3(:) = p(:) + 0.5d0 * k2(:)

! calculate the Lorentz factor
!
      gm = lorentz_factor(p3(:))

! calculate the velocity
!
      u3(:) = p3(:) / gm

! calculate the acceleration for the location x3 and velocity u3
!
      call acceleration(t3, x3(:), u3(:), s(:), a(:), v(:), b(:))

! calculate the third term
!
      l3(:) = dt * s(:)
      k3(:) = dt * a(:)

!! 4th step of the RK integration
!!
! integrate the position and momentum
!
      t4    = t    + dt
      x4(:) = x(:) + l3(:)
      p4(:) = p(:) + k3(:)

! calculate the Lorentz factor
!
      gm = lorentz_factor(p4(:))

! calculate the velocity
!
      u4(:) = p4(:) / gm

! calculate the acceleration for the location x4 and velocity u4
!
      call acceleration(t4, x4(:), u4(:), s(:), a(:), v(:), b(:))

! calculate the third term
!
      l4(:) = dt * s(:)
      k4(:) = dt * a(:)

!! the final integration of the particle position and momentum
!!
      t5    = t    + dt
      x5(:) = x(:) + ( l1(:) + 2.0d0 * ( l2(:) + l3(:) ) + l4(:) ) / 6.0d0
      p5(:) = p(:) + ( k1(:) + 2.0d0 * ( k2(:) + k3(:) ) + k4(:) ) / 6.0d0

! calculate the Lorentz factor
!
      gm = lorentz_factor(p5(:))

! calculate the velocity
!
      u5(:) = p5(:) / gm

! calculate the acceleration at the updated location
!
      call acceleration(t5, x5(:), u5(:), s(:), a(:), v(:), b(:))

! estimate the error for timestep control
!
      l4(:) = l4(:) - dt * s(:)
      k4(:) = k4(:) - dt * a(:)

      tol = sqrt(dot_product(l4(:), l4(:)) + dot_product(k4(:), k4(:))) / 6.0d0

! estimate the new timestep
!
      dtn   = dt * (safety * maxtol / tol)**0.2d0

! check if the error is below desired tolerance
!
      if (tol .gt. maxtol) then

! repeat the integration with a new timestep
!
        dt = dtn

      else

! update the time
!
        t   = t + dt

! check if time exceeded the maximum time
!
        if (t >= tmax) keepon = .false.

! update the new timestep
!
        dt = min(2.0d0 * dt, dtn, dtmax)

! update the position, velocity and momentum
!
        x(:) = x5(:)
        u(:) = u5(:)
        p(:) = p5(:)

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
        if (x(1) < bnds(1,1)) keepon = .false.
        if (x(1) > bnds(1,2)) keepon = .false.
        if (x(2) < bnds(2,1)) keepon = .false.
        if (x(2) > bnds(2,2)) keepon = .false.
        if (x(3) < bnds(3,1)) keepon = .false.
        if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! store the current particle state
!
        if (m .eq. ndumps) then

! separate the particle velocity into the parallel and perpendicular components
!
          call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
          call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate the particle energies
!
          en = gm * mrest
          ek = en - mrest

! print the progress
!
          write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua,  &
                                                            ek, term

! store the particle parameters
!
          write(10,"(20(1es22.14))") t, x(1), x(2), x(3), u(1), u(2), u(3),    &
                                     ua, up, ur, gm, en, ek,                   &
                                     bunit * ba, om / tunit, tg * tunit,       &
                                     rg * lunit, tg, rg, tol

! update the counters
!
          n = n + 1
          m = 0

        end if

! increase the data write counter
!
        m = m + 1

      end if

    end do

! separate the particle velocity into the parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate the particle energies
!
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, ua, ek

! store the particle parameters
!
    write(10,"(20(1es22.14))") t, x(1), x(2), x(3), u(1), u(2), u(3),          &
                               ua, up, ur, gm, en, ek,                         &
                               bunit * ba, om / tunit, tg * tunit,             &
                               rg * lunit, tg, rg, tol
    close (10)

!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_rk4
!
!===============================================================================
!
! integrate_trajectory_si4: subroutine integrates particle trajectory using
!                           the 4th order simplectic method
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!             "High order starting iterates for implicit Runge-Kutta methods:
!              an improvement for variable-step symplectic integrators", 2002,
!              IMA J. of Num. Ana., 22, 153
!
!===============================================================================
!
  subroutine integrate_trajectory_si4()

    implicit none

! local variables
!
    character(len=32)              :: str
    integer                        :: n, m, i = 0, mi, ti, k
    real(kind=8), dimension(2,6)   :: z
    real(kind=8), dimension(5,2,6) :: zp
    real(kind=8), dimension(3)     :: x, u, p, s, a
    real(kind=8), dimension(3)     :: xc, xe, xs
    real(kind=8), dimension(3)     :: pc, pe, ps
    real(kind=8), dimension(3)     :: v, b
    real(kind=8)                   :: gm, t, dt, tc, te, ts
    real(kind=8)                   :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=8)                   :: tol = 0.0d+00

! local flags
!
    logical                        :: keepon = .true.

! local parameters
!
    real(kind=8), parameter :: b1  = 5.0d-01
    real(kind=8), parameter :: ch  = sqrt(3.0d+00) / 6.0d+00
    real(kind=8), parameter :: c1  = b1 - ch, c2  = b1 + ch
    real(kind=8), parameter :: d1   = - sqrt(3.0d+00), d2   = - d1
!
!-------------------------------------------------------------------------------
!
! initialize the iteration number, snapshot number, time, and time steps
!
    n  = 1
    m  = 1
    k  = 5
    mi = 0
    ti = 0
    t  = 0.0d+00
    dt = dtini
    te = 0.0d+00

! reset the initial guess
!
    zp(:,:,:) = 0.0d+00

! reset the vector of the position and momentum errors
!
    xe(:) = 0.0d+00
    pe(:) = 0.0d+00

! substitute the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the starting point
!
    call acceleration(t, x(:), u(:), s(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate the particle energy
!
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua, ek, term

! open the output file, print headers and the initial values
!
    open(10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,20a22)") 'Time', 'X', 'Y', 'Z', 'Vx', 'Vy', 'Vz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]', '<B> [Gs]',  &
                                 'Omega [1/s]', 'Tg [s]', 'Rg [m]', 'Tg [T]',  &
                                 'Rg [L]', 'Tolerance', 'Iterations'
    write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), u(1), u(2), u(3),      &
                                   ua, up, ur, gm, en, ek, bunit * ba,         &
                                   om / tunit, tg * tunit, rg * lunit, tg, rg, &
                                   tol, i

!== INTEGRATION LOOP ==
!
! iterate until the maximum time is reached
!
    do while (keepon)

! find the initial guess for the vector Z using Newton's interpolation formula
!
      zp(5,:,:) = zp(4,:,:)
      zp(4,:,:) = zp(3,:,:)
      zp(3,:,:) = zp(2,:,:)
      zp(2,:,:) = zp(1,:,:)
      zp(1,:,:) = z (  :,:)
      if (k == 0) then
        z (  :,:) = 5.0d+00 * zp(1,:,:) - 1.0d+01 * zp(2,:,:)                  &
                  + 1.0d+01 * zp(3,:,:) - 5.0d+00 * zp(4,:,:) + zp(5,:,:)
      else if (k == 1) then
        z (  :,:) = 4.0d+00 * zp(1,:,:) - 6.0d+00 * zp(2,:,:)                  &
                  + 4.0d+00 * zp(3,:,:) - zp(4,:,:)
        k = k - 1
      else if (k == 2) then
        z (  :,:) = 3.0d+00 * zp(1,:,:) - 3.0d+00 * zp(2,:,:) + zp(3,:,:)
        k = k - 1
      else if (k == 3) then
        z (  :,:) = 2.0d+00 * zp(1,:,:) - zp(2,:,:)
        k = k - 1
      else if (k == 4) then
        z (  :,:) = zp(1,:,:)
        k = k - 1
      else if (k == 5) then
        z (  :,:) = 0.0d+00
        k = k - 1
      end if

! estimate the vector Z (eq. 5.3)
!
!   Z1 = dt * [ a11 * F(y + Z1) + a12 * F(y + Z2) ]
!   Z2 = dt * [ a21 * F(y + Z1) + a22 * F(y + Z2) ]
!
      call estimate_si4(x(:), p(:), z(:,:), t, dt, tol, i)

! update the solution
!
!   y(n+1) = y(n) + [ d1 * Z1 + d2 * Z2 ]
!
      xc(1:3) = (d1 * z(1,1:3) + d2 * z(2,1:3)) - xe(1:3)
      pc(1:3) = (d1 * z(1,4:6) + d2 * z(2,4:6)) - pe(1:3)
      xs(1:3) = x(1:3) + xc(1:3)
      ps(1:3) = p(1:3) + pc(1:3)
      xe(1:3) = (xs(1:3) - x(1:3)) - xc(1:3)
      pe(1:3) = (ps(1:3) - p(1:3)) - pc(1:3)
      x (1:3) = xs(1:3)
      p (1:3) = ps(1:3)

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
      if (x(1) < bnds(1,1)) keepon = .false.
      if (x(1) > bnds(1,2)) keepon = .false.
      if (x(2) < bnds(2,1)) keepon = .false.
      if (x(2) > bnds(2,2)) keepon = .false.
      if (x(3) < bnds(3,1)) keepon = .false.
      if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! update the integration time
!
      tc = dt - te
      ts = t  + tc
      te = (ts - t) - tc
      t  = ts

! check if time exceeded the maximum time
!
      if (t >= tmax) keepon = .false.

! find the maximum number of iteration in the estimator and update the counter
! of the total number of iterations
!
      mi = max(mi, i)
      ti = ti + i

! store the particle parameters at a given snapshot time
!
      if (m == ndumps) then

! calculate the Lorentz factor and particle velocity
!
        gm   = lorentz_factor(p(:))
        u(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(:), u(:), s(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
        en = gm * mrest
        ek = en - mrest

! print the progress
!
        write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua, ek,&
                                                          term

! write results to the output file
!
        write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), u(1), u(2), u(3),  &
                                       ua, up, ur, gm, en, ek, bunit * ba,     &
                                       om / tunit, tg * tunit, rg * lunit, tg, &
                                       rg, tol, i

        n = n + 1
        m = 0

      end if

! increase data write counter
!
      m = m + 1

! end of iteration
!
    end do

    close(10)

! print the progress
!
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, ua, ek

! write info about the estimator
!
    write(str,"(i12)") mi
    write(*,"('INFO      : maximum iterations per step = ',a)"      )          &
          trim(adjustl(str))
    write(*,"('INFO      : average iterations per step = ',1pe12.6)")          &
          real(ti, kind=8) / ((n - 1) * ndumps)
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_si4
!
!===============================================================================
!
! integrate_trajectory_si4: subroutine integrates particle trajectory using
!                           the 4th order simplectic method
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!             "High order starting iterates for implicit Runge-Kutta methods:
!              an improvement for variable-step symplectic integrators", 2002,
!              IMA J. of Num. Ana., 22, 153
!
!===============================================================================
!
  subroutine integrate_trajectory_si4v()

    implicit none

! local variables
!
    character(len=32)              :: str
    integer                        :: n, m, i = 0, mi, ti, k
    real(kind=8), dimension(2,6)   :: z
    real(kind=8), dimension(5,2,6) :: zp
    real(kind=8), dimension(5)     :: hp
    real(kind=8), dimension(3)     :: x, u, p, s, a
    real(kind=8), dimension(3)     :: xc, xe, xs
    real(kind=8), dimension(3)     :: pc, pe, ps
    real(kind=8), dimension(3)     :: v, b
    real(kind=8)                   :: gm, t, dt, dtp, tc, te, ts
    real(kind=8)                   :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=8)                   :: tol = 0.0d+00

! local flags
!
    logical                        :: keepon = .true.

! local parameters
!
    real(kind=8), parameter :: b1  = 5.0d-01
    real(kind=8), parameter :: ch  = sqrt(3.0d+00) / 6.0d+00
    real(kind=8), parameter :: c1  = b1 - ch, c2  = b1 + ch
    real(kind=8), parameter :: d1   = - sqrt(3.0d+00), d2   = - d1
!
!-------------------------------------------------------------------------------
!
! initialize the iteration number, snapshot number, time, and time steps
!
    n  = 1
    m  = 1
    k  = 5
    mi = 0
    ti = 0
    t  = 0.0d+00
    dt = dtini
    te = 0.0d+00

! reset the initial guess
!
    zp(:,:,:) = 0.0d+00

! reset the vector of the position and momentum errors
!
    xe(:) = 0.0d+00
    pe(:) = 0.0d+00

! substitute the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the starting point
!
    call acceleration(t, x(:), u(:), s(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate the particle energy
!
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua, ek, term

! open the output file, print headers and the initial values
!
    open(10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,20a22)") 'Time', 'X', 'Y', 'Z', 'Vx', 'Vy', 'Vz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]', '<B> [Gs]',  &
                                 'Omega [1/s]', 'Tg [s]', 'Rg [m]', 'Tg [T]',  &
                                 'Rg [L]', 'Tolerance', 'Iterations'
    write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), u(1), u(2), u(3),      &
                                   ua, up, ur, gm, en, ek, bunit * ba,         &
                                   om / tunit, tg * tunit, rg * lunit, tg, rg, &
                                   tol, i

!== INTEGRATION LOOP ==
!
! iterate until the maximum time is reached
!
    do while (keepon)

! find the initial guess for the vector Z using Newton's interpolation formula
!
      dtp   = dt
      hp(2) = hp(1)
      hp(1) = dt / dtp
      zp(3,:,:) = zp(2,:,:)
      zp(2,:,:) = zp(1,:,:)
      zp(1,:,:) = z (  :,:)
      if (k == 2) then
        z (  :,:) = zp(1,:,:) * (1.0d+00 + hp(1) + hp(1) * hp(1))              &
                  - zp(2,:,:) * (hp(1) * (1.0d+00 + hp(1) + hp(2)))            &
                  + zp(3,:,:) * (hp(1) * hp(2))
      else if (k == 3) then
        z (  :,:) = zp(1,:,:) + (zp(1,:,:) - zp(2,:,:)) * hp(1)
        k = k - 1
      else if (k == 4) then
        z (  :,:) = zp(1,:,:)
        k = k - 1
      else if (k == 5) then
        z (  :,:) = 0.0d+00
        k = k - 1
      end if

! estimate the vector Z (eq. 5.3)
!
!   Z1 = dt * [ a11 * F(y + Z1) + a12 * F(y + Z2) ]
!   Z2 = dt * [ a21 * F(y + Z1) + a22 * F(y + Z2) ]
!
      call estimate_si4(x(:), p(:), z(:,:), t, dt, tol, i)
      do while(tol >= maxtol)
        dt = dt * min(5.0d+00, max(1.0d-01, 0.8d+00 * (maxtol / tol)**(0.2)))
        call estimate_si4(x(:), p(:), z(:,:), t, dt, tol, i)
      end do

! update the solution
!
!   y(n+1) = y(n) + [ d1 * Z1 + d2 * Z2 ]
!
      xc(1:3) = (d1 * z(1,1:3) + d2 * z(2,1:3)) - xe(1:3)
      pc(1:3) = (d1 * z(1,4:6) + d2 * z(2,4:6)) - pe(1:3)
      xs(1:3) = x(1:3) + xc(1:3)
      ps(1:3) = p(1:3) + pc(1:3)
      xe(1:3) = (xs(1:3) - x(1:3)) - xc(1:3)
      pe(1:3) = (ps(1:3) - p(1:3)) - pc(1:3)
      x (1:3) = xs(1:3)
      p (1:3) = ps(1:3)

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
      if (x(1) < bnds(1,1)) keepon = .false.
      if (x(1) > bnds(1,2)) keepon = .false.
      if (x(2) < bnds(2,1)) keepon = .false.
      if (x(2) > bnds(2,2)) keepon = .false.
      if (x(3) < bnds(3,1)) keepon = .false.
      if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! update the integration time
!
      tc = dt - te
      ts = t  + tc
      te = (ts - t) - tc
      t  = ts

! check if time exceeded the maximum time
!
      if (t >= tmax) keepon = .false.

! find the maximum number of iteration in the estimator and update the counter
! of the total number of iterations
!
      mi = max(mi, i)
      ti = ti + i

! store the particle parameters at a given snapshot time
!
      if (m == ndumps) then

! calculate the Lorentz factor and particle velocity
!
        gm   = lorentz_factor(p(:))
        u(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(:), u(:), s(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
        en = gm * mrest
        ek = en - mrest

! print the progress
!
        write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua, ek,&
                                                          term

! write results to the output file
!
        write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), u(1), u(2), u(3),  &
                                       ua, up, ur, gm, en, ek, bunit * ba,     &
                                       om / tunit, tg * tunit, rg * lunit, tg, &
                                       rg, tol, i

        n = n + 1
        m = 0

      end if

! increase data write counter
!
      m = m + 1

! calculate new timestep
!
      dt = dt * min(2.7d+00 * maxit / (2.0d+00 * maxit + i), 1.0d+00)          &
              * min(6.0d+00, max(1.0d-01, (maxtol / tol)**(2.0d-01)))
      dt = min(dt, dtmax)

! end of iteration
!
    end do

    close(10)

! print the progress
!
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, ua, ek

! write info about the estimator
!
    write(str,"(i12)") mi
    write(*,"('INFO      : maximum iterations per step = ',a)"      )          &
          trim(adjustl(str))
    write(*,"('INFO      : average iterations per step = ',1pe12.6)")          &
          real(ti, kind=8) / ((n - 1) * ndumps)
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_si4v
!
!===============================================================================
!
! estimate_si4: subroutine estimates the solution for the equation of motion
!               using a simple functional iteration (SI4 version)
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!
! description: This subroutines find the solution of the equation 5.3 for
!              the increment Z using the functional iteration
!
!===============================================================================
!
  subroutine estimate_si4(x, p, z, t, dt, tol, n)

    implicit none

! subroutine arguments
!
    real(kind=8), dimension(3)  , intent(in)    :: x, p
    real(kind=8), dimension(2,6), intent(inout) :: z
    real(kind=8)                , intent(in)    :: t
    real(kind=8)                , intent(inout) :: dt, tol
    integer                     , intent(inout) :: n

! local variables
!
    integer                      :: i
    real(kind=8), dimension(2,6) :: zn
    real(kind=8), dimension(2,3) :: ui, si, ai
    real(kind=8), dimension(3)   :: xi, pi, xm, pm, v, b
    real(kind=8), dimension(2)   :: ti
    real(kind=8)                 :: lf

! local parameter
!
    real(kind=8), parameter :: b1  = 5.0d-01, bh = 2.5d-01
    real(kind=8), parameter :: ch  = sqrt(3.0d+00) / 6.0d+00
    real(kind=8), parameter :: c1  = b1 - ch, c2  = b1 + ch
    real(kind=8), parameter :: a11 = bh, a12 = bh - ch, a21 = bh + ch, a22 = bh
!
!-------------------------------------------------------------------------------
!
! initiate the iteration control parameters
!
    n   = 0
    tol = 2.0d+00 * maxtol

! prepare the time moments for intermedia states
!
    ti(1)   = t + c1 * dt
    ti(2)   = t + c2 * dt

! prepare normalized state to calculate tolerance
!
    xm(1:3) = max(1.0d+00, abs(x(1:3)))
    pm(1:3) = max(1.0d+00, abs(p(1:3)))

! perform fixed-point iteration
!
    do while (tol > maxtol .and. n < maxit)

! iterate over intermediate states
!
      do i = 1, 2

! prepare the particle intermediate state (position and momentum)
!
        xi(1:3) = x(:) + z(i,1:3)
        pi(1:3) = p(:) + z(i,4:6)

! convert particle momentum to velocity
!
        lf        = lorentz_factor(pi(1:3))
        ui(i,1:3) = pi(1:3) / lf

! get acceleration for the current state
!
        call acceleration(ti(i), xi(1:3), ui(i,1:3), si(i,1:3), ai(i,1:3), v(:), b(:))

      end do

! get the new increment estimate for the intermediate states
!
      zn(1,1:3) = dt * (a11 * si(1,1:3) + a12 * si(2,1:3))
      zn(2,1:3) = dt * (a21 * si(1,1:3) + a22 * si(2,1:3))
      zn(1,4:6) = dt * (a11 * ai(1,1:3) + a12 * ai(2,1:3))
      zn(2,4:6) = dt * (a21 * ai(1,1:3) + a22 * ai(2,1:3))

! calculate the maximum of residuum of the increment
!
      tol = 0.0d+00
      do i = 1, 2
        tol = max(tol, maxval(abs(zn(i,1:3) - z(i,1:3)) / xm(1:3)))
        tol = max(tol, maxval(abs(zn(i,4:6) - z(i,4:6)) / pm(1:3)))
      end do

! update the intermediate states with the new estimate
!
      z(:,:) = zn(:,:)

! increase the iteration counter
!
      n = n + 1

    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine estimate_si4
!
!===============================================================================
!
! integrate_trajectory_si6: subroutine integrates particle trajectory using
!                           the 6th order simplectic method
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!             "High order starting iterates for implicit Runge-Kutta methods:
!              an improvement for variable-step symplectic integrators", 2002,
!              IMA J. of Num. Ana., 22, 153
!
!===============================================================================
!
  subroutine integrate_trajectory_si6()

    implicit none

! local variables
!
    logical                        :: keepon = .true.
    character(len=32)              :: str
    integer                        :: n, m, i = 0, mi, ti, k
    real(kind=8), dimension(3,6)   :: z
    real(kind=8), dimension(5,3,6) :: zp
    real(kind=8), dimension(3)     :: x, u, p, s, a
    real(kind=8), dimension(3)     :: xc, xe, xs
    real(kind=8), dimension(3)     :: pc, pe, ps
    real(kind=8), dimension(3)     :: v, b
    real(kind=8)                   :: gm, t, dt, tc, te, ts
    real(kind=8)                   :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=8)                   :: tol = 0.0d+00

! local parameters
!
    real(kind=8), parameter :: b1   =   5.0d0 / 3.0d0                       &
                             , b2   = - 4.0d0 / 3.0d0                       &
                             , b3   =   5.0d0 / 3.0d0
    real(kind=8), parameter :: c1   =   0.5d0 - 0.1d0 * dsqrt(15.0d0)       &
                             , c2   =   0.5d0                               &
                             , c3   =   0.5d0 + 0.1d0 * dsqrt(15.0d0)
    real(kind=8), parameter :: d1   =   5.0d0 / 3.0d0                       &
                             , d2   = - 4.0d0 / 3.0d0                       &
                             , d3   =   5.0d0 / 3.0d0
!
!-------------------------------------------------------------------------------
!
! initialize the iteration number, snapshot number, time, and time steps
!
    n  = 1
    m  = 1
    k  = 5
    mi = 0
    ti = 0
    t  = 0.0d+00
    dt = dtini

! reset the initial guess
!
    zp(:,:,:) = 0.0d+00

! substitute the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the starting point
!
    call acceleration(t, x(:), u(:), s(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate the particle energy
!
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua, ek, term

! open the output file, print headers and the initial values
!
    open(10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,20a22)") 'Time', 'X', 'Y', 'Z', 'Vx', 'Vy', 'Vz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]', '<B> [Gs]',  &
                                 'Omega [1/s]', 'Tg [s]', 'Rg [m]', 'Tg [T]',  &
                                 'Rg [L]', 'Tolerance', 'Iterations'
    write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), u(1), u(2), u(3),      &
                                   ua, up, ur, gm, en, ek, bunit * ba,         &
                                   om / tunit, tg * tunit, rg * lunit, tg, rg, &
                                   tol, i

!== INTEGRATION LOOP ==
!
! iterate until the maximum time is reached
!
    do while (keepon)

! find the initial guess for the vector Z using Newton's interpolation formula
!
      zp(5,:,:) = zp(4,:,:)
      zp(4,:,:) = zp(3,:,:)
      zp(3,:,:) = zp(2,:,:)
      zp(2,:,:) = zp(1,:,:)
      zp(1,:,:) = z (  :,:)
      if (k == 0) then
        z (  :,:) = 5.0d+00 * zp(1,:,:) - 1.0d+01 * zp(2,:,:)                  &
                  + 1.0d+01 * zp(3,:,:) - 5.0d+00 * zp(4,:,:) + zp(5,:,:)
      else if (k == 1) then
        z (  :,:) = 4.0d+00 * zp(1,:,:) - 6.0d+00 * zp(2,:,:)                  &
                  + 4.0d+00 * zp(3,:,:) - zp(4,:,:)
        k = k - 1
      else if (k == 2) then
        z (  :,:) = 3.0d+00 * zp(1,:,:) - 3.0d+00 * zp(2,:,:) + zp(3,:,:)
        k = k - 1
      else if (k == 3) then
        z (  :,:) = 2.0d+00 * zp(1,:,:) - zp(2,:,:)
        k = k - 1
      else if (k == 4) then
        z (  :,:) = zp(1,:,:)
        k = k - 1
      else if (k == 5) then
        z (  :,:) = 0.0d+00
        k = k - 1
      end if

! estimate the vector Z (eq. 5.3)
!
!   Z1 = [ a11 * F(y + Z1) + a12 * F(y + Z2) + a13 * F(y + Z3) ]
!   Z2 = [ a21 * F(y + Z1) + a22 * F(y + Z2) + a23 * F(y + Z3) ]
!   Z3 = [ a31 * F(y + Z1) + a32 * F(y + Z2) + a33 * F(y + Z3) ]
!
      call estimate_si6(x(:), p(:), z(:,:), t, dt, tol, i)

! update the solution
!
!   y(n+1) = y(n) + [ d1 * Z1 + d2 * Z2 + d3 * Z3 ]
!
      xc(1:3) = (d1 * z(1,1:3) + d2 * z(2,1:3) + d3 * z(3,1:3)) - xe(1:3)
      pc(1:3) = (d1 * z(1,4:6) + d2 * z(2,4:6) + d3 * z(3,4:6)) - pe(1:3)
      xs(1:3) = x(1:3) + xc(1:3)
      ps(1:3) = p(1:3) + pc(1:3)
      xe(1:3) = (xs(1:3) - x(1:3)) - xc(1:3)
      pe(1:3) = (ps(1:3) - p(1:3)) - pc(1:3)
      x (1:3) = xs(1:3)
      p (1:3) = ps(1:3)

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
      if (x(1) < bnds(1,1)) keepon = .false.
      if (x(1) > bnds(1,2)) keepon = .false.
      if (x(2) < bnds(2,1)) keepon = .false.
      if (x(2) > bnds(2,2)) keepon = .false.
      if (x(3) < bnds(3,1)) keepon = .false.
      if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! update the integration time
!
      tc = dt - te
      ts = t  + tc
      te = (ts - t) - tc
      t  = ts

! check if time exceeded the maximum time
!
      if (t >= tmax) keepon = .false.

! find the maximum number of iteration in the estimator and update the counter
! of the total number of iterations
!
      mi = max(mi, i)
      ti = ti + i

! store the particle parameters at a given snapshot time
!
      if (m == ndumps) then

! calculate the Lorentz factor and particle velocity
!
        gm   = lorentz_factor(p(:))
        u(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(:), u(:), s(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
        en = gm * mrest
        ek = en - mrest

! print the progress
!
        write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua, ek,&
                                                          term

! write results to the output file
!
        write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), u(1), u(2), u(3),  &
                                       ua, up, ur, gm, en, ek, bunit * ba,     &
                                       om / tunit, tg * tunit, rg * lunit, tg, &
                                       rg, tol, i

        n = n + 1
        m = 0

      end if

! increase data write counter
!
      m = m + 1

! end of iteration
!
    end do

    close(10)

! print the progress
!
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, ua, ek

! write info about the estimator
!
    write(str,"(i12)") mi
    write(*,"('INFO      : maximum iterations per step = ',a)"      )          &
          trim(adjustl(str))
    write(*,"('INFO      : average iterations per step = ',1pe12.6)")          &
          real(ti, kind=8) / ((n - 1) * ndumps)
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_si6
!
!===============================================================================
!
! integrate_trajectory_si6v: subroutine integrates particle trajectory using
!                           the 6th order simplectic method
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!             "High order starting iterates for implicit Runge-Kutta methods:
!              an improvement for variable-step symplectic integrators", 2002,
!              IMA J. of Num. Ana., 22, 153
!
!===============================================================================
!
  subroutine integrate_trajectory_si6v()

    implicit none

! local variables
!
    logical                        :: keepon = .true.
    character(len=32)              :: str
    integer                        :: n, m, i = 0, mi, ti, k
    real(kind=8), dimension(3,6)   :: z
    real(kind=8), dimension(5,3,6) :: zp
    real(kind=8), dimension(5)     :: hp
    real(kind=8), dimension(3)     :: x, u, p, s, a
    real(kind=8), dimension(3)     :: xc, xe, xs
    real(kind=8), dimension(3)     :: pc, pe, ps
    real(kind=8), dimension(3)     :: v, b
    real(kind=8)                   :: gm, t, dt, dtp, tc, te, ts
    real(kind=8)                   :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=8)                   :: tol = 0.0d+00

! local parameters
!
    real(kind=8), parameter :: b1   =   5.0d0 / 3.0d0                       &
                             , b2   = - 4.0d0 / 3.0d0                       &
                             , b3   =   5.0d0 / 3.0d0
    real(kind=8), parameter :: c1   =   0.5d0 - 0.1d0 * dsqrt(15.0d0)       &
                             , c2   =   0.5d0                               &
                             , c3   =   0.5d0 + 0.1d0 * dsqrt(15.0d0)
    real(kind=8), parameter :: d1   =   5.0d0 / 3.0d0                       &
                             , d2   = - 4.0d0 / 3.0d0                       &
                             , d3   =   5.0d0 / 3.0d0
!
!-------------------------------------------------------------------------------
!
! initialize the iteration number, snapshot number, time, and time steps
!
    n  = 1
    m  = 1
    k  = 5
    mi = 0
    ti = 0
    t  = 0.0d+00
    dt = dtini

! reset the initial guess
!
    zp(:,:,:) = 0.0d+00

! substitute the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the starting point
!
    call acceleration(t, x(:), u(:), s(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate the particle energy
!
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua, ek, term

! open the output file, print headers and the initial values
!
    open(10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,20a22)") 'Time', 'X', 'Y', 'Z', 'Vx', 'Vy', 'Vz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]', '<B> [Gs]',  &
                                 'Omega [1/s]', 'Tg [s]', 'Rg [m]', 'Tg [T]',  &
                                 'Rg [L]', 'Tolerance', 'Iterations'
    write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), u(1), u(2), u(3),      &
                                   ua, up, ur, gm, en, ek, bunit * ba,         &
                                   om / tunit, tg * tunit, rg * lunit, tg, rg, &
                                   tol, i

!== INTEGRATION LOOP ==
!
! iterate until the maximum time is reached
!
    do while (keepon)

! find the initial guess for the vector Z using Newton's interpolation formula
!
      dtp   = dt
      hp(2) = hp(1)
      hp(1) = dt / dtp
      zp(3,:,:) = zp(2,:,:)
      zp(2,:,:) = zp(1,:,:)
      zp(1,:,:) = z (  :,:)
      if (k == 2) then
        z (  :,:) = zp(1,:,:) * (1.0d+00 + hp(1) + hp(1) * hp(1))              &
                  - zp(2,:,:) * (hp(1) * (1.0d+00 + hp(1) + hp(2)))            &
                  + zp(3,:,:) * (hp(1) * hp(2))
      else if (k == 3) then
        z (  :,:) = zp(1,:,:) + (zp(1,:,:) - zp(2,:,:)) * hp(1)
        k = k - 1
      else if (k == 4) then
        z (  :,:) = zp(1,:,:)
        k = k - 1
      else if (k == 5) then
        z (  :,:) = 0.0d+00
        k = k - 1
      end if

! estimate the vector Z (eq. 5.3)
!
!   Z1 = [ a11 * F(y + Z1) + a12 * F(y + Z2) + a13 * F(y + Z3) ]
!   Z2 = [ a21 * F(y + Z1) + a22 * F(y + Z2) + a23 * F(y + Z3) ]
!   Z3 = [ a31 * F(y + Z1) + a32 * F(y + Z2) + a33 * F(y + Z3) ]
!
      call estimate_si6(x(:), p(:), z(:,:), t, dt, tol, i)
      do while(tol >= maxtol)
        dt = dt * min(5.0d+00, max(1.0d-01, 0.8d+00 * (maxtol / tol)**(0.2)))
        call estimate_si6(x(:), p(:), z(:,:), t, dt, tol, i)
      end do

! update the solution
!
!   y(n+1) = y(n) + [ d1 * Z1 + d2 * Z2 + d3 * Z3 ]
!
      xc(1:3) = (d1 * z(1,1:3) + d2 * z(2,1:3) + d3 * z(3,1:3)) - xe(1:3)
      pc(1:3) = (d1 * z(1,4:6) + d2 * z(2,4:6) + d3 * z(3,4:6)) - pe(1:3)
      xs(1:3) = x(1:3) + xc(1:3)
      ps(1:3) = p(1:3) + pc(1:3)
      xe(1:3) = (xs(1:3) - x(1:3)) - xc(1:3)
      pe(1:3) = (ps(1:3) - p(1:3)) - pc(1:3)
      x (1:3) = xs(1:3)
      p (1:3) = ps(1:3)

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
      if (x(1) < bnds(1,1)) keepon = .false.
      if (x(1) > bnds(1,2)) keepon = .false.
      if (x(2) < bnds(2,1)) keepon = .false.
      if (x(2) > bnds(2,2)) keepon = .false.
      if (x(3) < bnds(3,1)) keepon = .false.
      if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! update the integration time
!
      tc = dt - te
      ts = t  + tc
      te = (ts - t) - tc
      t  = ts

! check if time exceeded the maximum time
!
      if (t >= tmax) keepon = .false.

! find the maximum number of iteration in the estimator and update the counter
! of the total number of iterations
!
      mi = max(mi, i)
      ti = ti + i

! store the particle parameters at a given snapshot time
!
      if (m == ndumps) then

! calculate the Lorentz factor and particle velocity
!
        gm   = lorentz_factor(p(:))
        u(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(:), u(:), s(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
        en = gm * mrest
        ek = en - mrest

! print the progress
!
        write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua, ek,&
                                                          term

! write results to the output file
!
        write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), u(1), u(2), u(3),  &
                                       ua, up, ur, gm, en, ek, bunit * ba,     &
                                       om / tunit, tg * tunit, rg * lunit, tg, &
                                       rg, tol, i

        n = n + 1
        m = 0

      end if

! increase data write counter
!
      m = m + 1

! calculate new timestep
!
      dt = dt * min(2.7d+00 * maxit / (2.0d+00 * maxit + i), 1.0d+00)          &
              * min(6.0d+00, max(1.0d-01, (maxtol / tol)**(2.0d-01)))
      dt = min(dt, dtmax)

! end of iteration
!
    end do

    close(10)

! print the progress
!
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, ua, ek

! write info about the estimator
!
    write(str,"(i12)") mi
    write(*,"('INFO      : maximum iterations per step = ',a)"      )          &
          trim(adjustl(str))
    write(*,"('INFO      : average iterations per step = ',1pe12.6)")          &
          real(ti, kind=8) / ((n - 1) * ndumps)
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_si6v
!
!===============================================================================
!
! estimate_si6: subroutine estimates the solution for the equation of motion
!               using a simple functional iteration (SI6 version)
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!
! description: This subroutines find the solution of the equation 5.3 for
!              the increment Z using the functional iteration
!
!===============================================================================
!
  subroutine estimate_si6(x, p, z, t, dt, tol, n)

    implicit none

! subroutine arguments
!
    real(kind=8), dimension(3)  , intent(in)    :: x, p
    real(kind=8), dimension(3,6), intent(inout) :: z
    real(kind=8)                , intent(in)    :: t
    real(kind=8)                , intent(inout) :: dt, tol
    integer                     , intent(inout) :: n

! local variables
!
    integer                      :: i
    real(kind=8), dimension(3,6) :: zn
    real(kind=8), dimension(3,3) :: ui, si, ai
    real(kind=8), dimension(3)   :: xi, pi, xm, pm, v, b
    real(kind=8), dimension(3)   :: ti
    real(kind=8)                 :: lf

! local parameter
!
    real(kind=8), parameter :: bh = 5.0d-01, ch  = sqrt(1.5d+01)
    real(kind=8), parameter :: a11 = 5.0d+00 / 3.6d+01                         &
                             , a12 = 2.0d+00 / 9.0d+00 - ch / 1.5d+01          &
                             , a13 = 5.0d+00 / 3.6d+01 - ch / 3.0d+01          &
                             , a21 = 5.0d+00 / 3.6d+01 + ch / 2.4d+01          &
                             , a22 = 2.0d+00 / 9.0d+00                         &
                             , a23 = 5.0d+00 / 3.6d+01 - ch / 2.4d+01          &
                             , a31 = 5.0d+00 / 3.6d+01 + ch / 3.0d+01          &
                             , a32 = 2.0d+00 / 9.0d+00 + ch / 1.5d+01          &
                             , a33 = 5.0d+00 / 3.6d+01
    real(kind=8), parameter :: c1  = bh - ch / 1.0d+01                         &
                             , c2  = bh                                        &
                             , c3  = bh + ch / 1.0d+01
!
!-------------------------------------------------------------------------------
!
! initiate the iteration control parameters
!
    n   = 0
    tol = 2.0d+00 * maxtol

! prepare the time moments for intermedia states
!
    ti(1)   = t + c1 * dt
    ti(2)   = t + c2 * dt
    ti(3)   = t + c3 * dt

! prepare normalized state to calculate tolerance
!
    xm(1:3) = max(1.0d+00, abs(x(1:3)))
    pm(1:3) = max(1.0d+00, abs(p(1:3)))

! perform fixed-point iteration
!
    do while (tol > maxtol .and. n < maxit)

! iterate over intermediate states
!
      do i = 1, 3

! prepare the particle intermediate state (position and momentum)
!
        xi(1:3) = x(:) + z(i,1:3)
        pi(1:3) = p(:) + z(i,4:6)

! convert particle momentum to velocity
!
        lf        = lorentz_factor(pi(1:3))
        ui(i,1:3) = pi(1:3) / lf

! get acceleration for the current state
!
        call acceleration(ti(i), xi(1:3), ui(i,1:3), si(i,1:3), ai(i,1:3), v(:), b(:))

      end do

! get the new increment estimate for the intermediate states
!
      zn(1,1:3) = dt * (a11 * si(1,1:3) + a12 * si(2,1:3) + a13 * si(3,1:3))
      zn(2,1:3) = dt * (a21 * si(1,1:3) + a22 * si(2,1:3) + a23 * si(3,1:3))
      zn(3,1:3) = dt * (a31 * si(1,1:3) + a32 * si(2,1:3) + a33 * si(3,1:3))
      zn(1,4:6) = dt * (a11 * ai(1,1:3) + a12 * ai(2,1:3) + a13 * ai(3,1:3))
      zn(2,4:6) = dt * (a21 * ai(1,1:3) + a22 * ai(2,1:3) + a23 * ai(3,1:3))
      zn(3,4:6) = dt * (a31 * ai(1,1:3) + a32 * ai(2,1:3) + a33 * ai(3,1:3))

! calculate the maximum of residuum of the increment
!
      tol = 0.0d+00
      do i = 1, 3
        tol = max(tol, maxval(abs(zn(i,1:3) - z(i,1:3)) / xm(1:3)))
        tol = max(tol, maxval(abs(zn(i,4:6) - z(i,4:6)) / pm(1:3)))
      end do

! update the intermediate states with the new estimate
!
      z(:,:) = zn(:,:)

! increase the iteration counter
!
      n = n + 1

    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine estimate_si6
!
!===============================================================================
!
! integrate_trajectory_si8: subroutine integrates particle trajectory using
!                           the 8th order simplectic method
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!             "High order starting iterates for implicit Runge-Kutta methods:
!              an improvement for variable-step symplectic integrators", 2002,
!              IMA J. of Num. Ana., 22, 153
!
!===============================================================================
!
  subroutine integrate_trajectory_si8()

    implicit none

! local variables
!
    logical                        :: keepon = .true.
    character(len=32)              :: str
    integer                        :: n, m, i = 0, mi, ti, k
    real(kind=8), dimension(4,6)   :: z
    real(kind=8), dimension(5,4,6) :: zp
    real(kind=8), dimension(3)     :: x, u, p, s, a
    real(kind=8), dimension(3)     :: xc, xe, xs
    real(kind=8), dimension(3)     :: pc, pe, ps
    real(kind=8), dimension(3)     :: v, b
    real(kind=8)                   :: gm, t, dt, tc, te, ts
    real(kind=8)                   :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=8)                   :: tol = 0.0d+00

! local parameters, Butcher's coefficients c_i and Sans-Serna & Calvo's
! coefficients d_i
!
    real(kind=8), parameter :: c1 =  6.9431844202973712388026755553596d-02     &
                             , c2 =  3.3000947820757186759866712044838d-01     &
                             , c3 =  6.6999052179242813240133287955162d-01     &
                             , c4 =  9.3056815579702628761197324444641d-01
    real(kind=8), parameter :: d1 = -1.6407053217392567182070402516331d+00     &
                             , d2 =  1.2143939697985776653621798588684d+00     &
                             , d3 = -1.2143939697985776653621798588684d+00     &
                             , d4 =  1.6407053217392567182070402516331d+00
!
!-------------------------------------------------------------------------------
!
! initialize the iteration number, snapshot number, time, and time steps
!
    n  = 1
    m  = 1
    k  = 5
    mi = 0
    ti = 0
    t  = 0.0d+00
    dt = dtini

! reset the initial guess
!
    zp(:,:,:) = 0.0d+00

! substitute the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the starting point
!
    call acceleration(t, x(:), u(:), s(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate the particle energy
!
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua, ek, term

! open the output file, print headers and the initial values
!
    open(10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,20a22)") 'Time', 'X', 'Y', 'Z', 'Vx', 'Vy', 'Vz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]', '<B> [Gs]',  &
                                 'Omega [1/s]', 'Tg [s]', 'Rg [m]', 'Tg [T]',  &
                                 'Rg [L]', 'Tolerance', 'Iterations'
    write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), u(1), u(2), u(3),      &
                                   ua, up, ur, gm, en, ek, bunit * ba,         &
                                   om / tunit, tg * tunit, rg * lunit, tg, rg, &
                                   tol, i

!== INTEGRATION LOOP ==
!
! iterate until the maximum time is reached
!
    do while (keepon)

! find the initial guess for the vector Z using Newton's interpolation formula
!
      zp(5,:,:) = zp(4,:,:)
      zp(4,:,:) = zp(3,:,:)
      zp(3,:,:) = zp(2,:,:)
      zp(2,:,:) = zp(1,:,:)
      zp(1,:,:) = z (  :,:)
      if (k == 0) then
        z (  :,:) = 5.0d+00 * zp(1,:,:) - 1.0d+01 * zp(2,:,:)                  &
                  + 1.0d+01 * zp(3,:,:) - 5.0d+00 * zp(4,:,:) + zp(5,:,:)
      else if (k == 1) then
        z (  :,:) = 4.0d+00 * zp(1,:,:) - 6.0d+00 * zp(2,:,:)                  &
                  + 4.0d+00 * zp(3,:,:) - zp(4,:,:)
        k = k - 1
      else if (k == 2) then
        z (  :,:) = 3.0d+00 * zp(1,:,:) - 3.0d+00 * zp(2,:,:) + zp(3,:,:)
        k = k - 1
      else if (k == 3) then
        z (  :,:) = 2.0d+00 * zp(1,:,:) - zp(2,:,:)
        k = k - 1
      else if (k == 4) then
        z (  :,:) = zp(1,:,:)
        k = k - 1
      else if (k == 5) then
        z (  :,:) = 0.0d+00
        k = k - 1
      end if

! estimate the vector Z (eq. 5.3)
!
!   Z1 = [ a11 * F(y + Z1) + a12 * F(y + Z2) + a13 * F(y + Z3) ]
!   Z2 = [ a21 * F(y + Z1) + a22 * F(y + Z2) + a23 * F(y + Z3) ]
!   Z3 = [ a31 * F(y + Z1) + a32 * F(y + Z2) + a33 * F(y + Z3) ]
!
      call estimate_si8(x(:), p(:), z(:,:), t, dt, tol, i)

! update the solution
!
!   y(n+1) = y(n) + [ d1 * Z1 + d2 * Z2 + d3 * Z3 ]
!
      xc(1:3) = (d1 * z(1,1:3) + d2 * z(2,1:3) + d3 * z(3,1:3)                 &
                                               + d4 * z(4,1:3)) - xe(1:3)
      pc(1:3) = (d1 * z(1,4:6) + d2 * z(2,4:6) + d3 * z(3,4:6)                 &
                                               + d4 * z(4,4:6)) - pe(1:3)
      xs(1:3) = x(1:3) + xc(1:3)
      ps(1:3) = p(1:3) + pc(1:3)
      xe(1:3) = (xs(1:3) - x(1:3)) - xc(1:3)
      pe(1:3) = (ps(1:3) - p(1:3)) - pc(1:3)
      x (1:3) = xs(1:3)
      p (1:3) = ps(1:3)

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
      if (x(1) < bnds(1,1)) keepon = .false.
      if (x(1) > bnds(1,2)) keepon = .false.
      if (x(2) < bnds(2,1)) keepon = .false.
      if (x(2) > bnds(2,2)) keepon = .false.
      if (x(3) < bnds(3,1)) keepon = .false.
      if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! update the integration time
!
      tc = dt - te
      ts = t  + tc
      te = (ts - t) - tc
      t  = ts

! check if time exceeded the maximum time
!
      if (t >= tmax) keepon = .false.

! find the maximum number of iteration in the estimator and update the counter
! of the total number of iterations
!
      mi = max(mi, i)
      ti = ti + i

! store the particle parameters at a given snapshot time
!
      if (m == ndumps) then

! calculate the Lorentz factor and particle velocity
!
        gm   = lorentz_factor(p(:))
        u(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(:), u(:), s(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
        en = gm * mrest
        ek = en - mrest

! print the progress
!
        write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua, ek,&
                                                          term

! write results to the output file
!
        write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), u(1), u(2), u(3),  &
                                       ua, up, ur, gm, en, ek, bunit * ba,     &
                                       om / tunit, tg * tunit, rg * lunit, tg, &
                                       rg, tol, i

        n = n + 1
        m = 0

      end if

! increase data write counter
!
      m = m + 1

! end of iteration
!
    end do

    close(10)

! print the progress
!
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, ua, ek

! write info about the estimator
!
    write(str,"(i12)") mi
    write(*,"('INFO      : maximum iterations per step = ',a)"      )          &
          trim(adjustl(str))
    write(*,"('INFO      : average iterations per step = ',1pe12.6)")          &
          real(ti, kind=8) / ((n - 1) * ndumps)
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_si8
!
!===============================================================================
!
! integrate_trajectory_si8: subroutine integrates particle trajectory using
!                           the 8th order simplectic method
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!             "High order starting iterates for implicit Runge-Kutta methods:
!              an improvement for variable-step symplectic integrators", 2002,
!              IMA J. of Num. Ana., 22, 153
!
!===============================================================================
!
  subroutine integrate_trajectory_si8v()

    implicit none

! local variables
!
    logical                        :: keepon = .true.
    character(len=32)              :: str
    integer                        :: n, m, i = 0, mi, ti, k
    real(kind=8), dimension(4,6)   :: z
    real(kind=8), dimension(5,4,6) :: zp
    real(kind=8), dimension(5)     :: hp
    real(kind=8), dimension(3)     :: x, u, p, s, a
    real(kind=8), dimension(3)     :: xc, xe, xs
    real(kind=8), dimension(3)     :: pc, pe, ps
    real(kind=8), dimension(3)     :: v, b
    real(kind=8)                   :: gm, t, dt, dtp, tc, te, ts
    real(kind=8)                   :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=8)                   :: tol = 0.0d+00

! local parameters, Butcher's coefficients c_i and Sans-Serna & Calvo's
! coefficients d_i
!
    real(kind=8), parameter :: c1 =  6.9431844202973712388026755553596d-02     &
                             , c2 =  3.3000947820757186759866712044838d-01     &
                             , c3 =  6.6999052179242813240133287955162d-01     &
                             , c4 =  9.3056815579702628761197324444641d-01
    real(kind=8), parameter :: d1 = -1.6407053217392567182070402516331d+00     &
                             , d2 =  1.2143939697985776653621798588684d+00     &
                             , d3 = -1.2143939697985776653621798588684d+00     &
                             , d4 =  1.6407053217392567182070402516331d+00
!
!-------------------------------------------------------------------------------
!
! initialize the iteration number, snapshot number, time, and time steps
!
    n  = 1
    m  = 1
    k  = 5
    mi = 0
    ti = 0
    t  = 0.0d+00
    dt = dtini

! reset the initial guess
!
    zp(:,:,:) = 0.0d+00

! substitute the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the starting point
!
    call acceleration(t, x(:), u(:), s(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate the particle energy
!
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua, ek, term

! open the output file, print headers and the initial values
!
    open(10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,20a22)") 'Time', 'X', 'Y', 'Z', 'Vx', 'Vy', 'Vz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]', '<B> [Gs]',  &
                                 'Omega [1/s]', 'Tg [s]', 'Rg [m]', 'Tg [T]',  &
                                 'Rg [L]', 'Tolerance', 'Iterations'
    write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), u(1), u(2), u(3),      &
                                   ua, up, ur, gm, en, ek, bunit * ba,         &
                                   om / tunit, tg * tunit, rg * lunit, tg, rg, &
                                   tol, i

!== INTEGRATION LOOP ==
!
! iterate until the maximum time is reached
!
    do while (keepon)

! find the initial guess for the vector Z using Newton's interpolation formula
!
      dtp   = dt
      hp(2) = hp(1)
      hp(1) = dt / dtp
      zp(3,:,:) = zp(2,:,:)
      zp(2,:,:) = zp(1,:,:)
      zp(1,:,:) = z (  :,:)
      if (k == 2) then
        z (  :,:) = zp(1,:,:) * (1.0d+00 + hp(1) + hp(1) * hp(1))              &
                  - zp(2,:,:) * (hp(1) * (1.0d+00 + hp(1) + hp(2)))            &
                  + zp(3,:,:) * (hp(1) * hp(2))
      else if (k == 3) then
        z (  :,:) = zp(1,:,:) + (zp(1,:,:) - zp(2,:,:)) * hp(1)
        k = k - 1
      else if (k == 4) then
        z (  :,:) = zp(1,:,:)
        k = k - 1
      else if (k == 5) then
        z (  :,:) = 0.0d+00
        k = k - 1
      end if

! estimate the vector Z (eq. 5.3)
!
!   Z1 = [ a11 * F(y + Z1) + a12 * F(y + Z2) + a13 * F(y + Z3) ]
!   Z2 = [ a21 * F(y + Z1) + a22 * F(y + Z2) + a23 * F(y + Z3) ]
!   Z3 = [ a31 * F(y + Z1) + a32 * F(y + Z2) + a33 * F(y + Z3) ]
!
      call estimate_si8(x(:), p(:), z(:,:), t, dt, tol, i)
      do while(tol >= maxtol)
        dt = dt * min(5.0d+00, max(1.0d-01, 0.8d+00 * (maxtol / tol)**(0.2)))
        call estimate_si8(x(:), p(:), z(:,:), t, dt, tol, i)
      end do

! update the solution
!
!   y(n+1) = y(n) + [ d1 * Z1 + d2 * Z2 + d3 * Z3 ]
!
      xc(1:3) = (d1 * z(1,1:3) + d2 * z(2,1:3) + d3 * z(3,1:3)                 &
                                               + d4 * z(4,1:3)) - xe(1:3)
      pc(1:3) = (d1 * z(1,4:6) + d2 * z(2,4:6) + d3 * z(3,4:6)                 &
                                               + d4 * z(4,4:6)) - pe(1:3)
      xs(1:3) = x(1:3) + xc(1:3)
      ps(1:3) = p(1:3) + pc(1:3)
      xe(1:3) = (xs(1:3) - x(1:3)) - xc(1:3)
      pe(1:3) = (ps(1:3) - p(1:3)) - pc(1:3)
      x (1:3) = xs(1:3)
      p (1:3) = ps(1:3)

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
      if (x(1) < bnds(1,1)) keepon = .false.
      if (x(1) > bnds(1,2)) keepon = .false.
      if (x(2) < bnds(2,1)) keepon = .false.
      if (x(2) > bnds(2,2)) keepon = .false.
      if (x(3) < bnds(3,1)) keepon = .false.
      if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! update the integration time
!
      tc = dt - te
      ts = t  + tc
      te = (ts - t) - tc
      t  = ts

! check if time exceeded the maximum time
!
      if (t >= tmax) keepon = .false.

! find the maximum number of iteration in the estimator and update the counter
! of the total number of iterations
!
      mi = max(mi, i)
      ti = ti + i

! store the particle parameters at a given snapshot time
!
      if (m == ndumps) then

! calculate the Lorentz factor and particle velocity
!
        gm   = lorentz_factor(p(:))
        u(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(:), u(:), s(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
        en = gm * mrest
        ek = en - mrest

! print the progress
!
        write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1,$)") n, t, dt, tg, ua, ek,&
                                                          term

! write results to the output file
!
        write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), u(1), u(2), u(3),  &
                                       ua, up, ur, gm, en, ek, bunit * ba,     &
                                       om / tunit, tg * tunit, rg * lunit, tg, &
                                       rg, tol, i

        n = n + 1
        m = 0

      end if

! increase data write counter
!
      m = m + 1

! calculate new timestep
!
      dt = dt * min(2.7d+00 * maxit / (2.0d+00 * maxit + i), 1.0d+00)          &
              * min(6.0d+00, max(1.0d-01, (maxtol / tol)**(2.0d-01)))
      dt = min(dt, dtmax)

! end of iteration
!
    end do

    close(10)

! print the progress
!
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, ua, ek

! write info about the estimator
!
    write(str,"(i12)") mi
    write(*,"('INFO      : maximum iterations per step = ',a)"      )          &
          trim(adjustl(str))
    write(*,"('INFO      : average iterations per step = ',1pe12.6)")          &
          real(ti, kind=8) / ((n - 1) * ndumps)
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_si8v
!
!===============================================================================
!
! estimate_si8: subroutine estimates the solution for the equation of motion
!               using a simple functional iteration (SI8 version)
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!
! description: This subroutines find the solution of the equation 5.3 for
!              the increment Z using the functional iteration
!
!===============================================================================
!
  subroutine estimate_si8(x, p, z, t, dt, tol, n)

    implicit none

! subroutine arguments
!
    real(kind=8), dimension(3)  , intent(in)    :: x, p
    real(kind=8), dimension(4,6), intent(inout) :: z
    real(kind=8)                , intent(in)    :: t
    real(kind=8)                , intent(inout) :: dt, tol
    integer                     , intent(inout) :: n

! local variables
!
    integer                      :: i
    real(kind=8), dimension(4,6) :: zn
    real(kind=8), dimension(4,3) :: ui, si, ai
    real(kind=8), dimension(3)   :: xi, pi, xm, pm, v, b
    real(kind=8), dimension(4)   :: ti
    real(kind=8)                 :: lf

! local parameters, Butcher's coefficients a_ij
!
    real(kind=8), parameter :: a11 =  8.6963711284363464343265987305500d-02    &
                             , a12 = -2.6604180084998793313385130476953d-02    &
                             , a13 =  1.2627462689404724515056880574618d-02    &
                             , a14 = -3.5551496857956831569109818495695d-03    &
                             , a21 =  1.8811811749986807165068554508717d-01    &
                             , a22 =  1.6303628871563653565673401269450d-01    &
                             , a23 = -2.7880428602470895224151106418997d-02    &
                             , a24 =  6.7355005945381555153986690857040d-03    &
                             , a31 =  1.6719192197418877317113330552530d-01    &
                             , a32 =  3.5395300603374396653761913180800d-01    &
                             , a33 =  1.6303628871563653565673401269450d-01    &
                             , a34 = -1.4190694931141142964153570476171d-02    &
                             , a41 =  1.7748257225452261184344295646057d-01    &
                             , a42 =  3.1344511474186834679841114481438d-01    &
                             , a43 =  3.5267675751627186462685315586596d-01    &
                             , a44 =  8.6963711284363464343265987305500d-02
    real(kind=8), parameter :: c1  =  6.9431844202973712388026755553596d-02    &
                             , c2  =  3.3000947820757186759866712044838d-01    &
                             , c3  =  6.6999052179242813240133287955162d-01    &
                             , c4  =  9.3056815579702628761197324444641d-01
!
!-------------------------------------------------------------------------------
!
! initiate the iteration control parameters
!
    n   = 0
    tol = 2.0d+00 * maxtol

! prepare the time moments for intermedia states
!
    ti(1)   = t + c1 * dt
    ti(2)   = t + c2 * dt
    ti(3)   = t + c3 * dt
    ti(4)   = t + c4 * dt

! prepare normalized state to calculate tolerance
!
    xm(1:3) = max(1.0d+00, abs(x(1:3)))
    pm(1:3) = max(1.0d+00, abs(p(1:3)))

! perform fixed-point iteration
!
    do while (tol > maxtol .and. n < maxit)

! iterate over intermediate states
!
      do i = 1, 4

! prepare the particle intermediate state (position and momentum)
!
        xi(1:3) = x(:) + z(i,1:3)
        pi(1:3) = p(:) + z(i,4:6)

! convert particle momentum to velocity
!
        lf        = lorentz_factor(pi(1:3))
        ui(i,1:3) = pi(1:3) / lf

! get acceleration for the current state
!
        call acceleration(ti(i), xi(1:3), ui(i,1:3), si(i,1:3), ai(i,1:3), v(:), b(:))

      end do

! get the new increment estimate for the intermediate states
!
      zn(1,1:3) = dt * (a11 * si(1,1:3) + a12 * si(2,1:3) + a13 * si(3,1:3)    &
                                                          + a14 * si(4,1:3))
      zn(2,1:3) = dt * (a21 * si(1,1:3) + a22 * si(2,1:3) + a23 * si(3,1:3)    &
                                                          + a24 * si(4,1:3))
      zn(3,1:3) = dt * (a31 * si(1,1:3) + a32 * si(2,1:3) + a33 * si(3,1:3)    &
                                                          + a34 * si(4,1:3))
      zn(4,1:3) = dt * (a41 * si(1,1:3) + a42 * si(2,1:3) + a43 * si(3,1:3)    &
                                                          + a44 * si(4,1:3))
      zn(1,4:6) = dt * (a11 * ai(1,1:3) + a12 * ai(2,1:3) + a13 * ai(3,1:3)    &
                                                          + a14 * ai(4,1:3))
      zn(2,4:6) = dt * (a21 * ai(1,1:3) + a22 * ai(2,1:3) + a23 * ai(3,1:3)    &
                                                          + a24 * ai(4,1:3))
      zn(3,4:6) = dt * (a31 * ai(1,1:3) + a32 * ai(2,1:3) + a33 * ai(3,1:3)    &
                                                          + a34 * ai(4,1:3))
      zn(4,4:6) = dt * (a41 * ai(1,1:3) + a42 * ai(2,1:3) + a43 * ai(3,1:3)    &
                                                          + a44 * ai(4,1:3))

! calculate the maximum of residuum of the increment
!
      tol = 0.0d+00
      do i = 1, 4
        tol = max(tol, maxval(abs(zn(i,1:3) - z(i,1:3)) / xm(1:3)))
        tol = max(tol, maxval(abs(zn(i,4:6) - z(i,4:6)) / pm(1:3)))
      end do

! update the intermediate states with the new estimate
!
      z(:,:) = zn(:,:)

! increase the iteration counter
!
      n = n + 1

    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine estimate_si8
!
!===============================================================================
!
! pos2index: subroutine converts a given position to the array index
!
!===============================================================================
!
  subroutine pos2index(x, r)

    use fields, only : nghosts

    implicit none

! input and output arguments
!
    real(kind=8), dimension(3), intent(in)  :: x
    real(kind=8), dimension(3), intent(out) :: r

! local variables
!
    integer      :: i
    real(kind=8) :: t
!
!------------------------------------------------------------------------------
!
    do i = 1, DIMS
      t    = (x(i) - bnds(i,1)) / bsiz(i) + 0.5 / dm(i)
      t    = t - floor(t)
      r(i) = dm(i) * t + nghosts
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
    real(kind=8), dimension(3), intent(in)  :: x
    integer     , dimension(4), intent(out) :: ii, jj, kk
    real(kind=8), dimension(3), intent(out) :: dr
    real(kind=8), dimension(4), intent(out) :: cx, cy, cz
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
    q1  = pcubic(cy, f(ii(1),jj(1),1), f(ii(1),jj(2),1)                        &
                   , f(ii(1),jj(3),1), f(ii(1),jj(4),1))
    q2  = pcubic(cy, f(ii(2),jj(1),1), f(ii(2),jj(2),1)                        &
                   , f(ii(2),jj(3),1), f(ii(2),jj(4),1))
    q3  = pcubic(cy, f(ii(3),jj(1),1), f(ii(3),jj(2),1)                        &
                   , f(ii(3),jj(3),1), f(ii(3),jj(4),1))
    q4  = pcubic(cy, f(ii(4),jj(1),1), f(ii(4),jj(2),1)                        &
                   , f(ii(4),jj(3),1), f(ii(4),jj(4),1))
#else /* DIMS == 2 */
!= tricubic interpolation =
!
! interpolate along Z direction
!
    q11 = pcubic(cz, f(ii(1),jj(1),kk(1)), f(ii(1),jj(1),kk(2))                &
                   , f(ii(1),jj(1),kk(3)), f(ii(1),jj(1),kk(4)))
    q12 = pcubic(cz, f(ii(1),jj(2),kk(1)), f(ii(1),jj(2),kk(2))                &
                   , f(ii(1),jj(2),kk(3)), f(ii(1),jj(2),kk(4)))
    q13 = pcubic(cz, f(ii(1),jj(3),kk(1)), f(ii(1),jj(3),kk(2))                &
                   , f(ii(1),jj(3),kk(3)), f(ii(1),jj(3),kk(4)))
    q14 = pcubic(cz, f(ii(1),jj(4),kk(1)), f(ii(1),jj(4),kk(2))                &
                   , f(ii(1),jj(4),kk(3)), f(ii(1),jj(4),kk(4)))

    q21 = pcubic(cz, f(ii(2),jj(1),kk(1)), f(ii(2),jj(1),kk(2))                &
                   , f(ii(2),jj(1),kk(3)), f(ii(2),jj(1),kk(4)))
    q22 = pcubic(cz, f(ii(2),jj(2),kk(1)), f(ii(2),jj(2),kk(2))                &
                   , f(ii(2),jj(2),kk(3)), f(ii(2),jj(2),kk(4)))
    q23 = pcubic(cz, f(ii(2),jj(3),kk(1)), f(ii(2),jj(3),kk(2))                &
                   , f(ii(2),jj(3),kk(3)), f(ii(2),jj(3),kk(4)))
    q24 = pcubic(cz, f(ii(2),jj(4),kk(1)), f(ii(2),jj(4),kk(2))                &
                   , f(ii(2),jj(4),kk(3)), f(ii(2),jj(4),kk(4)))

    q31 = pcubic(cz, f(ii(3),jj(1),kk(1)), f(ii(3),jj(1),kk(2))                &
                   , f(ii(3),jj(1),kk(3)), f(ii(3),jj(1),kk(4)))
    q32 = pcubic(cz, f(ii(3),jj(2),kk(1)), f(ii(3),jj(2),kk(2))                &
                   , f(ii(3),jj(2),kk(3)), f(ii(3),jj(2),kk(4)))
    q33 = pcubic(cz, f(ii(3),jj(3),kk(1)), f(ii(3),jj(3),kk(2))                &
                   , f(ii(3),jj(3),kk(3)), f(ii(3),jj(3),kk(4)))
    q34 = pcubic(cz, f(ii(3),jj(4),kk(1)), f(ii(3),jj(4),kk(2))                &
                   , f(ii(3),jj(4),kk(3)), f(ii(3),jj(4),kk(4)))

    q41 = pcubic(cz, f(ii(4),jj(1),kk(1)), f(ii(4),jj(1),kk(2))                &
                   , f(ii(4),jj(1),kk(3)), f(ii(4),jj(1),kk(4)))
    q42 = pcubic(cz, f(ii(4),jj(2),kk(1)), f(ii(4),jj(2),kk(2))                &
                   , f(ii(4),jj(2),kk(3)), f(ii(4),jj(2),kk(4)))
    q43 = pcubic(cz, f(ii(4),jj(3),kk(1)), f(ii(4),jj(3),kk(2))                &
                   , f(ii(4),jj(3),kk(3)), f(ii(4),jj(3),kk(4)))
    q44 = pcubic(cz, f(ii(4),jj(4),kk(1)), f(ii(4),jj(4),kk(2))                &
                   , f(ii(4),jj(4),kk(3)), f(ii(4),jj(4),kk(4)))

! interpolate along the Y direction
!
    q1  = pcubic(cy, q11, q12, q13, q14)
    q2  = pcubic(cy, q21, q22, q23, q24)
    q3  = pcubic(cy, q31, q32, q33, q34)
    q4  = pcubic(cy, q41, q42, q43, q44)
#endif /* DIMS == 2 */

! interpolate along the X direction
!
    q   = pcubic(cx, q1 , q2 , q3 , q4 )
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
  real(kind=8) function pcubic(c, fk, fl, fr, fq) result(q)

    implicit none

! input and output arguments
!
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
  real(kind=8) function lorentz_factor(p) result(gm)

    implicit none

    real(kind=8), dimension(3), intent(in) :: p
!
!------------------------------------------------------------------------------
!
    gm = sqrt(1.0d+00 + dot_product(p, p))
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
  subroutine acceleration(t, x, v, s, a, u, b)

    use fields, only : ux, uy, uz, bx, by, bz

    implicit none

! input and output arguments
!
    real(kind=8)              , intent(in)  :: t
    real(kind=8), dimension(3), intent(in)  :: x, v
    real(kind=8), dimension(3), intent(out) :: s, a, u, b

#ifdef TEST
! local variables
!
    real(kind=8), dimension(3) :: w
#ifdef ITEST
    real(kind=8)               :: dl, ra, rb, xt, yt, rt
#endif /* ITEST */
#else /* TEST */
! local variables
!
    real(kind=8), dimension(3) :: r, w

! position indices
!
    integer                    :: dist
    integer     , dimension(4) :: ii, jj, kk
    real(kind=8), dimension(4) :: cx, cy, cz
    real(kind=8), dimension(3) :: dr
#endif /* TEST */
!
!-------------------------------------------------------------------------------
!
#ifdef TEST
#ifdef WTEST
        u(1) =   0.0d+00
        u(2) = - vamp * sin(pi2 * freq * x(1))
        u(3) =   0.0d+00

        b(1) =   bpar
        b(2) =   bamp * cos(pi2 * freq * x(1))
        b(3) =   bamp * sin(pi2 * freq * x(1))
#endif /* WTEST */
#ifdef ITEST
! calculate the local velocity
!
        u(1) =      - vamp * x(1)
        u(2) = vrat * vamp * x(2)
        u(3) = 0.0d0

! calculate the local magnetic field
!
        dl   = bamp / bini
        ra   = 1.0d0 + dl
        rb   = 1.0d0 - dl

        xt   = x(1) / ra
        yt   = x(2) / rb

        rt   = dsqrt(xt * xt + yt * yt)

        if (rt .gt. 0.0d0) then
          b(1) =   yt / rb / rt * bini
          b(2) = - xt / ra / rt * bini
          b(3) = bshr
        else
          b(1) = 0.0d0
          b(2) = 0.0d0
          b(3) = bshr
        end if
#endif /* ITEST */
#else /* TEST */
! convert position to index
!
      call pos2index(x, r)

#ifdef BNDRY
      dist = min(minval(dm(1:DIMS) - r(1:DIMS)), minval(r(1:DIMS)))
      if (dist .gt. nghost) then
#endif /* BNDRY */

! prepare coefficients for interpolation
!
        call prepare_interpolation(r, ii, jj, kk, dr, cx, cy, cz)

! interpolate field components at the particle position
!
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

! normalize the speed
!
        s(:) = v(:) / vunit

! compute the acceleration
!
        a(1) = qom * (w(2) * b(3) - w(3) * b(2))
        a(2) = qom * (w(3) * b(1) - w(1) * b(3))
        a(3) = qom * (w(1) * b(2) - w(2) * b(1))
#if !defined TEST && defined BNDRY
      else
        a(:) = 0.0
      endif
#endif /* !TEST & BNDRY */
!
!-------------------------------------------------------------------------------
!
  end subroutine acceleration
!
!===============================================================================
!
! subroutine FIELDS:
! -----------------
!
!   Subroutine returns plasma fields for a given position.
!
!   Arguments:
!
!     x  - the particle position (in code units), input;
!     u  - the local plasma velocity (in units of c), output;
!     b  - the local magnetic field (renormalized by q/m), output;
!     j  - the local current density (renormalized by q/m), output, optional;
!
!===============================================================================
!
  subroutine fields(x, u, b)

! import required modules
!
    use coordinates   , only : map_position_to_index
    use fields        , only : nghosts
    use fields        , only : ux, uy, uz
    use fields        , only : bx, by, bz
    use interpolations, only : prepare_interpolation, interpolate

    implicit none

! subroutine arguments
!
    real(kind=8), dimension(3), intent(in)            :: x
    real(kind=8), dimension(3), intent(out)           :: u, b

! local variables
!
    integer     , dimension(4,3) :: ip
    real(kind=8), dimension(4,3) :: cc
    real(kind=8), dimension(3)   :: dr
    real(kind=8), dimension(3)   :: r
!
!-------------------------------------------------------------------------------
!
! map position to domain index
!
    call map_position_to_index(x(:), r(:))

! shift the position to domain interior
!
    r(:) = r(:) + nghosts

! prepare interpolation coefficients
!
    call prepare_interpolation(r(:), ip(:,:), dr(:), cc(:,:))

! interpolate plasma field components
!
    u(1) = interpolate(ux(:,:,:), ip(:,:), dr(:), cc(:,:))
    u(2) = interpolate(uy(:,:,:), ip(:,:), dr(:), cc(:,:))
    u(3) = interpolate(uz(:,:,:), ip(:,:), dr(:), cc(:,:))
    b(1) = interpolate(bx(:,:,:), ip(:,:), dr(:), cc(:,:))
    b(2) = interpolate(by(:,:,:), ip(:,:), dr(:), cc(:,:))
    b(3) = interpolate(bz(:,:,:), ip(:,:), dr(:), cc(:,:))

!-------------------------------------------------------------------------------
!
  end subroutine fields
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
    real(kind=8), dimension(3), intent(in)  :: v, b
    real(kind=8)              , intent(out) :: ba, vu, vp, vr

! local variables
!
    real(kind=8), dimension(3) :: p
    real(kind=8)               :: pp, vv
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
    real(kind=8), intent(in)  :: gm
    real(kind=8), intent(in)  :: ba, vr
    real(kind=8), intent(out) :: om, tg, rg
!
!------------------------------------------------------------------------------
!
    om = abs(qom) * ba / gm
    tg = pi2 / om
    rg = vr / (vunit * om)
!
!-------------------------------------------------------------------------------
!
  end subroutine gyro_parameters
!
end module particles
