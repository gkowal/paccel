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

! the units of time, length, velocity and magnetic field
!
  real(kind=8), save :: tunit = 1.0d+00, lunit = 1.0d+00
  real(kind=8), save :: vunit = 1.0d+00, bunit = 1.0d+00

! particle mass and the ratio q/m
!
  real(kind=8), save :: mrest = 9.38272088161040869636d+02
  real(kind=8), save :: qom   = 9.57883315593801671639d+03

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
  real(kind=8), save :: bpar
#endif /* TEST */

! output data parameters
!
  integer     , save :: ndumps = 1000    ! number of steps between subsequent dumps
  real(kind=8), save :: tmin   = 1.0d-03 ! minimum time of writing data
  real(kind=8), save :: tmax   = 1.0d+00 ! maximum time for integration

! integration parameters
!
  real(kind=8), save :: safety = 5.0d-01  ! safety coefficient
  real(kind=8), save :: maxtol = 1.0d-04  ! the maximum integration tolerance
  real(kind=8), save :: dtini  = 1.0d-08  ! the initial time step
  real(kind=8), save :: dtmax  = 1.0d+00  ! maximum allowed step size
  integer     , save :: maxit  = 1000     ! the limit of iterations
  real(kind=8), save :: atol   = 1.0d-06
  real(kind=8), save :: rtol   = 1.0d-06
  real(kind=8), save :: safe   = 9.00d-01
  real(kind=8), save :: beta   = 0.00d+00
  real(kind=8), save :: facmin = 3.33d-01
  real(kind=8), save :: facmax = 6.00d+00

! arrays containing the initial positions and velocities of particle
!
  real(kind=8), dimension(3), save :: x0, p0

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
    call get_parameter('dtini' , dtini)
    call get_parameter('dtmax' , dtmax)
    call get_parameter('maxit' , maxit)
    call get_parameter('atol'  , atol  )
    call get_parameter('rtol'  , rtol  )
    call get_parameter('safe'  , safe  )
    call get_parameter('beta'  , beta  )
    call get_parameter('facmin', facmin)
    call get_parameter('facmax', facmax)

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
    real(kind=8), dimension(3) :: r, b, u, w, v0
!
!-------------------------------------------------------------------------------
!
    status = 0

! print info
!
    if (verbose) write(*,"('INFO',6x,': initializing particle initial state')")

! get parameters
!
    call get_parameter('xp'    , x0(1))
    call get_parameter('yp'    , x0(2))
    call get_parameter('zp'    , x0(3))
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
! get plasma field components
!
    call fields(x0(:), u(:), b(:))
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
#if DIMS == 3
    call random_number(w)
    w = w - 0.5d+00
#else /* DIMS == 3 */
    w(1) = 0.0d+00
    w(2) = 0.0d+00
    w(3) = 1.0d+00
#endif /* DIMS == 3 */

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
    v0(:) = (vpar * b(:) + vper * u(:))

! calculate the initial particle momentuum
!
    lfac = 1.0d+00 / sqrt(1.0d+00 - dot_product(v0, v0))
    p0(:) = lfac * v0(:)

! allow to set the particle moment explicitely
!
    call get_parameter('px', p0(1))
    call get_parameter('py', p0(2))
    call get_parameter('pz', p0(3))

    lfac  = lorentz_factor(p0(:))
    v0(:) = p0(:) / lfac
    vabs  = sqrt(dot_product(v0(:), v0(:)))
    vpar  = abs(dot_product(v0(:), b(:)))
    vper  = sqrt(vabs**2 - vpar**2)

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
! subroutine INTEGRATE_RK4:
! ------------------------
!
!   Subroutine integrates particle trajectory using the classic 4th order
!   explicit Runge-Kutta method.
!
!===============================================================================
!
  subroutine integrate_trajectory_rk4()

! import required modules
!
    use coordinates, only : is_inside

    implicit none

! local variables
!
    logical                        :: keepon = .true.
    integer                        :: n, m
    real(kind=8)                   :: t, dt, dtn, tt, err
    real(kind=8)                   :: gm, ba, va, vp, vr, om, tg, rg, en, ek
    real(kind=8), dimension(3)     :: v, u, b
    real(kind=8), dimension(3,2)   :: ss, si, sr, er
    real(kind=8), dimension(3,2,5) :: ff

! local parameters
!
    real(kind=8), parameter        :: expo = -2.0d-01
!
!-------------------------------------------------------------------------------
!
! initialize counters, time and timesteps
!
    n  = 0
    m  = 0
    t  = 0.0d+00
    dt = dtini

! set the initial position and momentum
!
    si(:,1) = x0(:)
    si(:,2) = p0(:)

! calculate the particle parameters at the initial state
!
    gm   = lorentz_factor(si(:,2))
    v(:) = si(:,2) / gm
    call acceleration(t, si(:,1), si(:,2), ff(:,1,1), ff(:,2,1), u(:), b(:))
    call separate_velocity(v(:), b(:), ba, va, vp, vr)
    call gyro_parameters(gm, ba, vr, om, tg, rg)
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')            &
                                                    n, t, dt, tg, va, ek, term

! open the output file, print headers and the initial values
!
    open (10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,19a22)") 'Time', 'X', 'Y', 'Z', 'Px', 'Py', 'Pz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]',              &
                                 '<B> [Gs]', 'Omega [1/s]',                    &
                                 'Tg [s]', 'Rg [m]', 'Tg [T]', 'Rg [L]',       &
                                 'Tolerance'
    write(10,"(20(1es22.14))") t, si(:,:),                                     &
                               va, vp, vr, gm, en, ek,                         &
                               bunit * ba, om / tunit, tg * tunit, rg * lunit, &
                               tg, rg, err

!== INTEGRATION LOOP ==
!
! integrate the trajectory
!
    do while (keepon)

!! 1st step of the RK integration
!!
      tt      = t
      ss(:,:) = si(:,:)

      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,1), ff(:,2,1), u(:), b(:))

!! 2nd step of the RK integration
!!
      tt      = t       + 5.0d-01 * dt
      ss(:,:) = si(:,:) + 5.0d-01 * dt * ff(:,:,1)

      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,2), ff(:,2,2), u(:), b(:))

!! 3rd step of the RK integration
!!
      tt      = t       + 5.0d-01 * dt
      ss(:,:) = si(:,:) + 5.0d-01 * dt * ff(:,:,2)

      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,3), ff(:,2,3), u(:), b(:))

!! 4th step of the RK integration
!!
      tt      = t       + dt
      ss(:,:) = si(:,:) + dt * ff(:,:,3)

      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,4), ff(:,2,4), u(:), b(:))

!! the final integration of the particle position and momentum
!!
      tt      = t       + dt
      ss(:,:) = si(:,:) + dt * (ff(:,:,1)                                      &
                   + 2.0d+00 * (ff(:,:,2) + ff(:,:,3)) + ff(:,:,4)) / 6.0d+00

! estimate the error for timestep control
!
      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,1), ff(:,2,1), u(:), b(:))
      sr(:,:) = atol + rtol * max(abs(si(:,:)), abs(ss(:,:)))
      er(:,:) = ff(:,:,4) - ff(:,:,1)
      err = abs(dt) * sqrt(sum((er(:,:) / sr(:,:))**2)) / 6.0d+00

! update the solution, if the error is small, otherwise reduce
! the time step and repeat
!
      if (err <= 1.0d+00) then

        t       = tt
        si(:,:) = ss(:,:)

! store the particle state, if desired
!
        if (m >= ndumps) then

          gm   = lorentz_factor(si(:,2))
          v(:) = si(:,2) / gm
          call acceleration(t, si(:,1), si(:,2), ff(:,1,1), ff(:,2,1), u(:), b(:))
          call separate_velocity(v(:), b(:), ba, va, vp, vr)
          call gyro_parameters(gm, ba, vr, om, tg, rg)
          en = gm * mrest
          ek = en - mrest

! print the progress
!
          write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')      &
                                                    n, t, dt, tg, va, ek, term

! store the particle parameters
!
          write(10,"(20(1es22.14))") t, si(:,:),                               &
                                     va, vp, vr, gm, en, ek,                   &
                                     bunit * ba, om / tunit, tg * tunit,       &
                                     rg * lunit, tg, rg, err

! update the counters
!
          n = n + 1
          m = 0

        end if

! increase the data write counter
!
        m = m + 1

! check if the particle time did not exceed the maximum time and
! if the particle is still inside the domain
!
        keepon = (t < tmax) .and. is_inside(si(:,1))

! determine the new timestep
!
        dtn = dt * min(facmax, max(facmin, safe * err**expo))

      else

        dtn = dt *             max(facmin, safe * err**expo)

      end if

! substitute time step
!
      dt = min(dtn, dtmax)

    end do

! calculate the particle parameters at the final state
!
    if (m > 1) then

      gm   = lorentz_factor(si(:,2))
      v(:) = si(:,2) / gm
      call acceleration(t, si(:,1), si(:,2), ff(:,1,1), ff(:,2,1), u(:), b(:))
      call separate_velocity(v(:), b(:), ba, va, vp, vr)
      call gyro_parameters(gm, ba, vr, om, tg, rg)
      en = gm * mrest
      ek = en - mrest

! print the progress
!
      write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, va, ek

! store the particle parameters
!
      write(10,"(20(1es22.14))") t, si(:,:),                                   &
                                 va, vp, vr, gm, en, ek,                       &
                                 bunit * ba, om / tunit, tg * tunit,           &
                                 rg * lunit, tg, rg, err

    end if

    close (10)

!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_rk4
!
!===============================================================================
!
! subroutine INTEGRATE_DOP853:
! ---------------------------
!
!   Subroutine integrates particle trajectory using the 8th order explicit
!   Runge-Kutta method by Dorman & Prince with step control and dense output.
!
!   References:
!
!   [1] "Solving Ordinary Differential Equations I - Nonstiff Problems",
!       E. Hairer, S. P. Nørsett, G. Wanner,
!       Springer-Verlag Berlin Heidelberg, 2008
!
!===============================================================================
!
  subroutine integrate_trajectory_dop853()

! import required modules
!
    use coordinates, only : is_inside

    implicit none

! local variables
!
    logical                         :: keepon = .true., rejected = .false.
    integer                         :: n, m
    real(kind=8)                    :: t, dt, dtn, tt
    real(kind=8)                    :: gm, en, ek, ba, va, vp, vr, om, tg, rg
    real(kind=8)                    :: err, err3, err5, errold, deno, expo
    real(kind=8), dimension(3)      :: v, u, b
    real(kind=8), dimension(3,2)    :: si, ss, sr, er
    real(kind=8), dimension(3,2,10) :: ff

! parameters
!
    real(kind=8), parameter :: a0201 =  5.26001519587677318785587544488d-02    &
                             , a0301 =  1.97250569845378994544595329183d-02    &
                             , a0302 =  5.91751709536136983633785987549d-02    &
                             , a0401 =  2.95875854768068491816892993775d-02    &
                             , a0403 =  8.87627564304205475450678981324d-02    &
                             , a0501 =  2.41365134159266685502369798665d-01    &
                             , a0503 = -8.84549479328286085344864962717d-01    &
                             , a0504 =  9.24834003261792003115737966543d-01    &
                             , a0601 =  3.70370370370370370370370370370d-02    &
                             , a0604 =  1.70828608729473871279604482173d-01    &
                             , a0605 =  1.25467687566822425016691814123d-01    &
                             , a0701 =  3.71093750000000000000000000000d-02    &
                             , a0704 =  1.70252211019544039314978060272d-01    &
                             , a0705 =  6.02165389804559606850219397283d-02    &
                             , a0706 = -1.75781250000000000000000000000d-02    &
                             , a0801 =  3.70920001185047927108779319836d-02    &
                             , a0804 =  1.70383925712239993810214054705d-01    &
                             , a0805 =  1.07262030446373284651809199168d-01    &
                             , a0806 = -1.53194377486244017527936158236d-02    &
                             , a0807 =  8.27378916381402288758473766002d-03    &
                             , a0901 =  6.24110958716075717114429577812d-01    &
                             , a0904 = -3.36089262944694129406857109825d+00    &
                             , a0905 = -8.68219346841726006818189891453d-01    &
                             , a0906 =  2.75920996994467083049415600797d+01    &
                             , a0907 =  2.01540675504778934086186788979d+01    &
                             , a0908 = -4.34898841810699588477366255144d+01    &
                             , a1001 =  4.77662536438264365890433908527d-01    &
                             , a1004 = -2.48811461997166764192642586468d+00    &
                             , a1005 = -5.90290826836842996371446475743d-01    &
                             , a1006 =  2.12300514481811942347288949897d+01    &
                             , a1007 =  1.52792336328824235832596922938d+01    &
                             , a1008 = -3.32882109689848629194453265587d+01    &
                             , a1009 = -2.03312017085086261358222928593d-02    &
                             , a1101 = -9.37142430085987325717040216580d-01    &
                             , a1104 =  5.18637242884406370830023853209d+00    &
                             , a1105 =  1.09143734899672957818500254654d+00    &
                             , a1106 = -8.14978701074692612513997267357d+00    &
                             , a1107 = -1.85200656599969598641566180701d+01    &
                             , a1108 =  2.27394870993505042818970056734d+01    &
                             , a1109 =  2.49360555267965238987089396762d+00    &
                             , a1110 = -3.04676447189821950038236690220d+00    &
                             , a1201 =  2.27331014751653820792359768449d+00    &
                             , a1204 = -1.05344954667372501984066689879d+01    &
                             , a1205 = -2.00087205822486249909675718444d+00    &
                             , a1206 = -1.79589318631187989172765950534d+01    &
                             , a1207 =  2.79488845294199600508499808837d+01    &
                             , a1208 = -2.85899827713502369474065508674d+00    &
                             , a1209 = -8.87285693353062954433549289258d+00    &
                             , a1210 =  1.23605671757943030647266201528d+01    &
                             , a1211 =  6.43392746015763530355970484046d-01

    real(kind=8), parameter :: b01 =  5.42937341165687622380535766363d-02      &
                             , b06 =  4.45031289275240888144113950566d+00      &
                             , b07 =  1.89151789931450038304281599044d+00      &
                             , b08 = -5.80120396001058478146721142270d+00      &
                             , b09 =  3.11164366957819894408916062370d-01      &
                             , b10 = -1.52160949662516078556178806805d-01      &
                             , b11 =  2.01365400804030348374776537501d-01      &
                             , b12 =  4.47106157277725905176885569043d-02

    real(kind=8), parameter :: c2  = 5.260015195876773187855875444880d-02,     &
                               c3  = 7.890022793815159781783813167320d-02,     &
                               c4  = 1.183503419072273967267571975098d-01,     &
                               c5  = 2.816496580927726032732428024902d-01,     &
                               c6  = 3.333333333333333333333333333333d-01,     &
                               c7  = 2.500000000000000000000000000000d-01,     &
                               c8  = 3.076923076923076923076923076923d-01,     &
                               c9  = 6.512820512820512820512820512821d-01,     &
                               c10 = 6.000000000000000000000000000000d-01,     &
                               c11 = 8.571428571428571428571428571429d-01

    real(kind=8), parameter :: bh1 = 0.244094488188976377952755905512d+00      &
                             , bh2 = 0.733846688281611857341361741547d+00      &
                             , bh3 = 0.220588235294117647058823529412d-01

    real(kind=8), parameter :: er01 =  0.1312004499419488073250102996d-01      &
                             , er06 = -0.1225156446376204440720569753d+01      &
                             , er07 = -0.4957589496572501915214079952d+00      &
                             , er08 =  0.1664377182454986536961530415d+01      &
                             , er09 = -0.3503288487499736816886487290d+00      &
                             , er10 =  0.3341791187130174790297318841d+00      &
                             , er11 =  0.8192320648511571246570742613d-01      &
                             , er12 = -0.2235530786388629525884427845d-01
!
!-------------------------------------------------------------------------------
!
! initialize counters, time and timesteps
!
    n  = 0
    m  = 0
    t  = 0.0d+00
    dt = dtini
    err = 0.0d+00

    expo     = 2.0d-01 * beta - 1.25d-01
    errold   = 1.0d-04

    keepon   = .true.
    rejected = .false.

! set the initial position, velocity, and momentum
!
    si(:,1) = x0(:)
    si(:,2) = p0(:)

! determine the initial state of the particle
!
    gm   = lorentz_factor(si(:,2))
    v(:) = si(:,2) / gm
    call acceleration(t, si(:,1), si(:,2), ff(:,1,1), ff(:,2,1), u(:), b(:))
    call separate_velocity(v(:), b(:), ba, va, vp, vr)
    call gyro_parameters(gm, ba, vr, om, tg, rg)
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')            &
                                                    n, t, dt, tg, va, ek, term

! open the output file, print headers and the initial values
!
    open (10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,19a22)") 'Time', 'X', 'Y', 'Z', 'Px', 'Py', 'Pz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]',              &
                                 '<B> [Gs]', 'Omega [1/s]',                    &
                                 'Tg [s]', 'Rg [m]', 'Tg [T]', 'Rg [L]',       &
                                 'Error'
    write(10,"(20(1es22.14))") t, si(:,:),                                     &
                               va, vp, vr, gm, en, ek,                         &
                               bunit * ba, om / tunit, tg * tunit, rg * lunit, &
                               tg, rg, err

!== INTEGRATION LOOP ==
!
! integrate the trajectory
!
    do while (keepon)

! the 12 steps
!
      tt      = t
      ss(:,:) = si(:,:)
      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,1), ff(:,2,1), u(:), b(:))

      tt      = t       + dt * c2
      ss(:,:) = si(:,:) + dt * a0201 * ff(:,:,1)
      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,2), ff(:,2,2), u(:), b(:))

      tt      = t       + dt * c3
      ss(:,:) = si(:,:) + dt * (a0301 * ff(:,:,1) + a0302 * ff(:,:,2))
      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,3), ff(:,2,3), u(:), b(:))

      tt      = t       + dt * c4
      ss(:,:) = si(:,:) + dt * (a0401 * ff(:,:,1) + a0403 * ff(:,:,3))
      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,4), ff(:,2,4), u(:), b(:))

      tt      = t       + dt * c5
      ss(:,:) = si(:,:) + dt * (a0501 * ff(:,:,1) + a0503 * ff(:,:,3)          &
                              + a0504 * ff(:,:,4))
      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,5), ff(:,2,5), u(:), b(:))

      tt      = t       + dt * c6
      ss(:,:)  = si(:,:) + dt * (a0601 * ff(:,:,1) + a0604 * ff(:,:,4)         &
                               + a0605 * ff(:,:,5))
      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,6), ff(:,2,6), u(:), b(:))

      tt      = t       + dt * c7
      ss(:,:)  = si(:,:) + dt * (a0701 * ff(:,:,1) + a0704 * ff(:,:,4)         &
                               + a0705 * ff(:,:,5) + a0706 * ff(:,:,6))
      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,7), ff(:,2,7), u(:), b(:))

      tt      = t       + dt * c8
      ss(:,:)  = si(:,:) + dt * (a0801 * ff(:,:,1) + a0804 * ff(:,:,4)         &
                               + a0805 * ff(:,:,5) + a0806 * ff(:,:,6)         &
                               + a0807 * ff(:,:,7))
      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,8), ff(:,2,8), u(:), b(:))

      tt      = t       + dt * c9
      ss(:,:)  = si(:,:) + dt * (a0901 * ff(:,:,1) + a0904 * ff(:,:,4)         &
                               + a0905 * ff(:,:,5) + a0906 * ff(:,:,6)         &
                               + a0907 * ff(:,:,7) + a0908 * ff(:,:,8))
      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,9), ff(:,2,9), u(:), b(:))

      tt      = t       + dt * c10
      ss(:,:)  = si(:,:) + dt * (a1001 * ff(:,:,1) + a1004 * ff(:,:,4)         &
                               + a1005 * ff(:,:,5) + a1006 * ff(:,:,6)         &
                               + a1007 * ff(:,:,7) + a1008 * ff(:,:,8)         &
                               + a1009 * ff(:,:,9))
      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,10), ff(:,2,10), u(:), b(:))

      tt      = t       + dt * c11
      ss(:,:)  = si(:,:) + dt * (a1101 * ff(:,:,1) + a1104 * ff(:,:, 4)        &
                               + a1105 * ff(:,:,5) + a1106 * ff(:,:, 6)        &
                               + a1107 * ff(:,:,7) + a1108 * ff(:,:, 8)        &
                               + a1109 * ff(:,:,9) + a1110 * ff(:,:,10))
      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,2), ff(:,2,2), u(:), b(:))

      tt      = t       + dt
      ss(:,:)  = si(:,:) + dt * (a1201 * ff(:,:,1) + a1204 * ff(:,:, 4)        &
                               + a1205 * ff(:,:,5) + a1206 * ff(:,:, 6)        &
                               + a1207 * ff(:,:,7) + a1208 * ff(:,:, 8)        &
                               + a1209 * ff(:,:,9) + a1210 * ff(:,:,10)        &
                               + a1211 * ff(:,:,2))
      call acceleration(tt, ss(:,1), ss(:,2), ff(:,1,3), ff(:,2,3), u(:), b(:))

      ff(:,:,4) = b01 * ff(:,:, 1) + b06 * ff(:,:, 6) + b07 * ff(:,:, 7)       &
                + b08 * ff(:,:, 8) + b09 * ff(:,:, 9) + b10 * ff(:,:,10)       &
                + b11 * ff(:,:, 2) + b12 * ff(:,:, 3)
      ss(:,:)  = si(:,:) + dt * ff(:,:,4)

! error estimation
!
      sr(:,:) = atol + rtol * max(abs(si(:,:)), abs(ss(:,:)))
      er(:,:) = ff(:,:,4) - bh1 * ff(:,:,1) - bh2 * ff(:,:,9) - bh3 * ff(:,:,3)
      err3 = sum((er(:,:) / sr(:,:))**2) ! 3rd order error
      er(:,:) = er01 * ff(:,:,1) + er06 * ff(:,:,6) + er07 * ff(:,:, 7)        &
              + er08 * ff(:,:,8) + er09 * ff(:,:,9) + er10 * ff(:,:,10)        &
              + er11 * ff(:,:,2) + er12 * ff(:,:,3)
      err5 = sum((er(:,:) / sr(:,:))**2) ! 5th order error
      deno = err5 + 1.0d-02 * err3
      if (deno <= 0.0d+00) then
        err = abs(dt) * err5 / sqrt(6.0d+00)
      else
        err = abs(dt) * err5 / sqrt(6.0d+00 * deno)
      end if

! update the solution, if the error is small, otherwise reduce
! the time step and repeat
!
      if (err <= 1.0d+00) then

        t       = tt
        si(:,:) = ss(:,:)

! store the particle state, if desired
!
        if (m >= ndumps) then

          gm   = lorentz_factor(si(:,2))
          v(:) = si(:,2) / gm
          call acceleration(t, si(:,1), si(:,2), ff(:,1,1), ff(:,2,1), u(:), b(:))
          call separate_velocity(v(:), b(:), ba, va, vp, vr)
          call gyro_parameters(gm, ba, vr, om, tg, rg)
          en = gm * mrest
          ek = en - mrest

! print the progress
!
          write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')      &
                                                    n, t, dt, tg, va, ek, term

! store the particle parameters
!
          write(10,"(20(1es22.14))") t, si(:,:),                               &
                                     va, vp, vr, gm, en, ek,                   &
                                     bunit * ba, om / tunit, tg * tunit,       &
                                     rg * lunit, tg, rg, err

! update the counters
!
          n = n + 1
          m = 0

        end if

! increase the data write counter
!
        m = m + 1

! check if the particle time did not exceed the maximum time and
! if the particle is still inside the domain
!
        keepon = (t < tmax) .and. is_inside(si(:,1))

! new time step
!
        dtn   = max(facmin, safe * err**expo * errold**beta)
        if (rejected) then
          dtn = dt * min(dtn, 1.0d+00)
        else
          dtn = dt * min(dtn, facmax)
        end if
        errold = max(err, 1.0d-04)

! set rejected flag
!
        rejected = .false.

      else

! new time step
!
        dtn = dt * max(facmin, safe * err**expo)

! set rejected flag
!
        rejected = .true.

      end if

! substitute time step
!
      dt = min(dtn, dtmax)

    end do

! calculate the particle parameters at the final state
!
    if (m > 1) then

      gm   = lorentz_factor(si(:,2))
      v(:) = si(:,2) / gm
      call acceleration(t, si(:,1), si(:,2), ff(:,1,1), ff(:,2,1), u(:), b(:))
      call separate_velocity(v(:), b(:), ba, va, vp, vr)
      call gyro_parameters(gm, ba, vr, om, tg, rg)
      en = gm * mrest
      ek = en - mrest

! print the progress
!
      write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, va, ek

! store the particle parameters
!
      write(10,"(20(1es22.14))") t, si(:,:),                                   &
                                 va, vp, vr, gm, en, ek,                       &
                                 bunit * ba, om / tunit, tg * tunit,           &
                                 rg * lunit, tg, rg, err

    end if

    close (10)

!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_dop853
!
!===============================================================================
!
! subroutine INTEGRATE_SI4:
! ------------------------
!
!   Subroutine integrates the particle trajectory using the 4th order implicit
!   symplectic Gauss-Legendre Runge-Kutta method.
!
!   References:
!
!   [1] "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!       Chapman & Hall, London, New York, 1994
!   [2] "High order starting iterates for implicit Runge-Kutta methods:
!       an improvement for variable-step symplectic integrators", M. P. Calvo,
!       IMA Journal of Numerical Analysis, 2002, vol. 22, pp. 153-166
!
!===============================================================================
!
  subroutine integrate_trajectory_si4()

! import required modules
!
    use coordinates, only : is_inside

    implicit none

! local variables
!
    character(len=32)                :: str
    integer                          :: n, m, i = 0, mi, ti, k
    real(kind=8)                     :: gm, t, dt, tc, te, ts
    real(kind=8)                     :: ba, va, vp, vr, om, tg, rg, en, ek
    real(kind=8)                     :: err
    real(kind=8), dimension(3)       :: v, u, b
    real(kind=8), dimension(3,2)     :: si, sc, se, ss, ff
    real(kind=8), dimension(3,2,2)   :: zi
    real(kind=8), dimension(3,2,2,5) :: zp

! local flags
!
    logical                        :: keepon = .true.

! local parameters
!
    real(kind=8), parameter :: b1 =  5.0000000000000000000d-01,                &
                               b2 =  5.0000000000000000000d-01
    real(kind=8), parameter :: c1 =  2.1132486540518711775d-01,                &
                               c2 =  7.8867513459481288225d-01
    real(kind=8), parameter :: d1 = -1.7320508075688772935d+00,                &
                               d2 =  1.7320508075688772935d+00
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
    te = 0.0d+00
    dt = dtini
    err = 0.0d+00

! reset the initial guess
!
    zp(:,:,:,:) = 0.0d+00
    zi(:,:,:)   = 0.0d+00

! reset the vector of the position and momentum errors
!
    se(:,:) = 0.0d+00

! set the initial position, velocity, and momentum
!
    si(:,1) = x0(:)
    si(:,2) = p0(:)

! determine the initial state of the particle
!
    gm   = lorentz_factor(si(:,2))
    v(:) = si(:,2) / gm
    call acceleration(t, si(:,1), si(:,2), ff(:,1), ff(:,2), u(:), b(:))
    call separate_velocity(v(:), b(:), ba, va, vp, vr)
    call gyro_parameters(gm, ba, vr, om, tg, rg)
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')            &
                                                    n, t, dt, tg, va, ek, term

! open the output file, print headers and the initial values
!
    open(10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,20a22)") 'Time', 'X', 'Y', 'Z', 'Px', 'Py', 'Pz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]', '<B> [Gs]',  &
                                 'Omega [1/s]', 'Tg [s]', 'Rg [m]', 'Tg [T]',  &
                                 'Rg [L]', 'Tolerance', 'Iterations'
    write(10,"(20(1es22.14),i22)") t, si(:,:),                                 &
                                   va, vp, vr, gm, en, ek, bunit * ba,         &
                                   om / tunit, tg * tunit, rg * lunit, tg, rg, &
                                   err, i

!== INTEGRATION LOOP ==
!
! iterate until the maximum time is reached
!
    do while (keepon)

! find the initial guess for the vector Z using Newton's interpolation formula
!
      if (k == 0) then
        zp(:,:,:,5) = zp(:,:,:,4)
        zp(:,:,:,4) = zp(:,:,:,3)
        zp(:,:,:,3) = zp(:,:,:,2)
        zp(:,:,:,2) = zp(:,:,:,1)
        zp(:,:,:,1) = zi(:,:,:  )
        zi(:,:,:)   = 5.0d+00 * zp(:,:,:,1) - 1.0d+01 * zp(:,:,:,2)            &
                    + 1.0d+01 * zp(:,:,:,3) - 5.0d+00 * zp(:,:,:,4)            &
                              + zp(:,:,:,5)
      else if (k == 1) then
        zp(:,:,:,4) = zp(:,:,:,3)
        zp(:,:,:,3) = zp(:,:,:,2)
        zp(:,:,:,2) = zp(:,:,:,1)
        zp(:,:,:,1) = zi(:,:,:  )
        zi(:,:,:)   = 4.0d+00 * zp(:,:,:,1) - 6.0d+00 * zp(:,:,:,2)            &
                    + 4.0d+00 * zp(:,:,:,3) -           zp(:,:,:,4)
        k = k - 1
      else if (k == 2) then
        zp(:,:,:,3) = zp(:,:,:,2)
        zp(:,:,:,2) = zp(:,:,:,1)
        zp(:,:,:,1) = zi(:,:,:  )
        zi(:,:,:)   = 3.0d+00 * zp(:,:,:,1) - 3.0d+00 * zp(:,:,:,2)            &
                              + zp(:,:,:,3)
        k = k - 1
      else if (k == 3) then
        zp(:,:,:,2) = zp(:,:,:,1)
        zp(:,:,:,1) = zi(:,:,:  )
        zi(:,:,:)   = 2.0d+00 * zp(:,:,:,1) - zp(:,:,:,2)
        k = k - 1
      else if (k == 4) then
        zp(:,:,:,1) = zi(:,:,:  )
        zi(:,:,:)   = zp(:,:,:,1)
        k = k - 1
      else if (k == 5) then

! calculate the acceleration at the initial state
!
        call acceleration(t, si(:,1), si(:,2), ff(:,1), ff(:,2), u(:), b(:))

! find the initial guess for the increment Z
!
        zi(:,:,1) = c1 * dt * ff(:,:)
        zi(:,:,2) = c2 * dt * ff(:,:)

        k = k - 1
      end if

! estimate the new increment Z (eq. 5.3)
!
!   Z1 = dt * [ a11 * F(y + Z1) + a12 * F(y + Z2) ]
!   Z2 = dt * [ a21 * F(y + Z1) + a22 * F(y + Z2) ]
!
      call estimate_si4(t, dt, si(:,:), zi(:,:,:), err, i)

! update the solution
!
!   y(n+1) = y(n) + [ b1 * Z1 + b2 * Z2 ]
!
      sc(:,:) = (d1 * zi(:,:,1) + d2 * zi(:,:,2)) - se(:,:)
      ss(:,:) =  si(:,:) + sc(:,:)
      se(:,:) = (ss(:,:) - si(:,:)) - sc(:,:)
      si(:,:) =  ss(:,:)

! update the time
!
      tc =  dt - te
      ts =  t  + tc
      te = (ts - t) - tc
      t  =  ts

! store the particle parameters at a given snapshot time
!
      if (m == ndumps) then

        gm   = lorentz_factor(si(:,2))
        v(:) = si(:,2) / gm
        call acceleration(t, si(:,1), si(:,2), ff(:,1), ff(:,2), u(:), b(:))
        call separate_velocity(v(:), b(:), ba, va, vp, vr)
        call gyro_parameters(gm, ba, vr, om, tg, rg)
        en = gm * mrest
        ek = en - mrest

! print the progress
!
        write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')        &
                                                    n, t, dt, tg, va, ek, term

! write results to the output file
!
        write(10,"(20(1es22.14),i22)") t, si(:,:),                             &
                                       va, vp, vr, gm, en, ek, bunit * ba,     &
                                       om / tunit, tg * tunit, rg * lunit, tg, &
                                       rg, err, i

        n = n + 1
        m = 0

      end if

! increase data write counter
!
      m = m + 1

! check if the particle time did not exceed the maximum time and
! if the particle is still inside the domain
!
      keepon = (t < tmax) .and. is_inside(si(:,1))

! find the maximum number of iteration in the estimator and update the counter
! of the total number of iterations
!
      mi = max(mi, i)
      ti = ti + i

! end of iteration
!
    end do

! calculate the particle parameters at the final state
!
    gm   = lorentz_factor(si(:,2))
    v(:) = si(:,2) / gm
    call acceleration(t, si(:,1), si(:,2), ff(:,1), ff(:,2), u(:), b(:))
    call separate_velocity(v(:), b(:), ba, va, vp, vr)
    call gyro_parameters(gm, ba, vr, om, tg, rg)
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, va, ek

! store the particle parameters
!
    if (m > 1) then
      write(10,"(20(1es22.14))") t, si(:,:),                                   &
                                 va, vp, vr, gm, en, ek,                       &
                                 bunit * ba, om / tunit, tg * tunit,           &
                                 rg * lunit, tg, rg, err
    end if

    close(10)

! write info about the estimator
!
    write(str,"(i12)") mi
    write(*,"('INFO      : maximum iterations per step = ',a)"      )          &
          trim(adjustl(str))
    write(*,"('INFO      : average iterations per step = ',1pe12.6)")          &
          real(ti, kind=8) / ((n - 1) * ndumps)

!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_si4
!
!===============================================================================
!
! subroutine ESTIMATE_SI4:
! -----------------------
!
!   Subroutine estimates the intermediate steps for the SI4 method using
!   a simple functional iteration.
!
!   Arguments:
!
!     s(:,:)   - the initial particle state (position and moment);
!     z(:,:,:) - the intermediate particle states;
!     dt       - the position incremental step;
!
!   References:
!
!   [1] "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!       Chapman & Hall, London, New York, 1994
!
!===============================================================================
!
  subroutine estimate_si4(t, dt, si, zi, err, n)

    implicit none

! subroutine arguments
!
    real(kind=8)                  , intent(in)    :: t, dt
    real(kind=8), dimension(3,2)  , intent(in)    :: si
    real(kind=8), dimension(3,2,2), intent(inout) :: zi
    real(kind=8)                  , intent(out)   :: err
    integer                       , intent(inout) :: n

! local variables
!
    integer                        :: m
    real(kind=8), dimension(3)     :: u, b
    real(kind=8), dimension(2)     :: ti
    real(kind=8), dimension(3,2,2) :: zr, zn
    real(kind=8), dimension(3,2,2) :: fi

! local parameter
!
    real(kind=8), parameter :: b1  =  5.0000000000000000000d-01,               &
                               b2  =  5.0000000000000000000d-01
    real(kind=8), parameter :: c1  =  2.1132486540518711775d-01,               &
                               c2  =  7.8867513459481288225d-01
    real(kind=8), parameter :: d1  = -1.7320508075688772935d+00,               &
                               d2  =  1.7320508075688772935d+00
    real(kind=8), parameter :: a11 =  2.5000000000000000000d-01,               &
                               a12 = -3.8675134594812882255d-02,               &
                               a21 =  5.3867513459481288225d-01,               &
                               a22 =  2.5000000000000000000d-01
!
!-------------------------------------------------------------------------------
!
! initiate the iteration control parameters
!
    n   = 0
    err = huge(err)

! prepare the time moments for intermedia states
!
    ti(1)   = t + c1 * dt
    ti(2)   = t + c2 * dt

! perform fixed-point iteration
!
    do while (err > 1.0d+00 .and. n < maxit)

! iterate over intermediate states
!
      do m = 1, 2
        call acceleration(ti(m), si(:,1) + zi(:,1,m), si(:,2) + zi(:,2,m),     &
                                 fi(:,1,m), fi(:,2,m), u(:), b(:))
      end do

! get the new increment estimate for the intermediate states
!
      zn(:,:,1) = dt * (a11 * fi(:,:,1) + a12 * fi(:,:,2))
      zn(:,:,2) = dt * (a21 * fi(:,:,1) + a22 * fi(:,:,2))

! prepare normalized state to calculate tolerance
!
      zr(:,:,:) = atol + rtol * abs(zi(:,:,:))

! calculate the maximum of residuum of the increment
!
      err = 0.0d+00
      do m = 1, 2
        err = max(err, maxval(abs(zn(:,:,m) - zi(:,:,m)) / zr(:,:,m)))
      end do

! substitute the new solution of the increment
!
      zi(:,:,:) = zn(:,:,:)

! increase the iteration counter
!
      n = n + 1

    end do

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

! import required modules
!
    use coordinates, only : is_inside

    implicit none

! local variables
!
    logical                        :: keepon = .true.
    character(len=32)              :: str
    integer                        :: n, m, i = 0, mi, ti, k
    real(kind=8), dimension(3,6)   :: z
    real(kind=8), dimension(5,3,6) :: zp
    real(kind=8), dimension(3)     :: x, v, p, s, a, u, b
    real(kind=8), dimension(3)     :: xc, xe, xs
    real(kind=8), dimension(3)     :: pc, pe, ps
    real(kind=8)                   :: gm, t, dt, tc, te, ts
    real(kind=8)                   :: ba, va, vp, vr, om, tg, rg, en, ek
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
    p(:) = p0(:)

! calculate the Lorentz factor and particle speed
!
    gm = lorentz_factor(p(:))
    v(:) = p(:) / gm

! calculate the acceleration at the initial position
!
    call acceleration(t, x(:), p(:), s(:), a(:), u(:), b(:))

! separate the particle velocity into parallel and perpendicular components
!
    call separate_velocity(v(:), b(:), ba, va, vp, vr)

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, vr, om, tg, rg)

! calculate the particle energies
!
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')            &
                                                    n, t, dt, tg, va, ek, term

! open the output file, print headers and the initial values
!
    open(10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,20a22)") 'Time', 'X', 'Y', 'Z', 'Px', 'Py', 'Pz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]', '<B> [Gs]',  &
                                 'Omega [1/s]', 'Tg [s]', 'Rg [m]', 'Tg [T]',  &
                                 'Rg [L]', 'Tolerance', 'Iterations'
    write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), p(1), p(2), p(3),      &
                                   va, vp, vr, gm, en, ek, bunit * ba,         &
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

! if the boundaries are not periodic and particle is out of the box, stop
! the integration
!
      keepon = keepon .and. is_inside(x)

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
        v(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(:), p(:), s(:), a(:), u(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(v(:), b(:), ba, va, vp, vr)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, vr, om, tg, rg)

! calculate particle energy
!
        en = gm * mrest
        ek = en - mrest

! print the progress
!
        write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')        &
                                                    n, t, dt, tg, va, ek, term

! write results to the output file
!
        write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), p(1), p(2), p(3),  &
                                       va, vp, vr, gm, en, ek, bunit * ba,     &
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
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, va, ek

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

! import required modules
!
    use coordinates, only : is_inside

    implicit none

! local variables
!
    logical                        :: keepon = .true.
    character(len=32)              :: str
    integer                        :: n, m, i = 0, mi, ti, k
    real(kind=8), dimension(3,6)   :: z
    real(kind=8), dimension(5,3,6) :: zp
    real(kind=8), dimension(5)     :: hp
    real(kind=8), dimension(3)     :: x, v, p, s, a, u, b
    real(kind=8), dimension(3)     :: xc, xe, xs
    real(kind=8), dimension(3)     :: pc, pe, ps
    real(kind=8)                   :: gm, t, dt, dtp, tc, te, ts
    real(kind=8)                   :: ba, va, vp, vr, om, tg, rg, en, ek
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
    p(:) = p0(:)

! calculate the Lorentz factor and particle speed
!
    gm = lorentz_factor(p(:))
    v(:) = p(:) / gm

! calculate the acceleration at the initial position
!
    call acceleration(t, x(:), p(:), s(:), a(:), u(:), b(:))

! separate the particle velocity into parallel and perpendicular components
!
    call separate_velocity(v(:), b(:), ba, va, vp, vr)

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, vr, om, tg, rg)

! calculate the particle energies
!
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')            &
                                                    n, t, dt, tg, va, ek, term

! open the output file, print headers and the initial values
!
    open(10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,20a22)") 'Time', 'X', 'Y', 'Z', 'Px', 'Py', 'Pz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]', '<B> [Gs]',  &
                                 'Omega [1/s]', 'Tg [s]', 'Rg [m]', 'Tg [T]',  &
                                 'Rg [L]', 'Tolerance', 'Iterations'
    write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), p(1), p(2), p(3),      &
                                   va, vp, vr, gm, en, ek, bunit * ba,         &
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

! if the boundaries are not periodic and particle is out of the box, stop
! the integration
!
      keepon = keepon .and. is_inside(x)

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
        v(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(:), p(:), s(:), a(:), u(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(v(:), b(:), ba, va, vp, vr)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, vr, om, tg, rg)

! calculate particle energy
!
        en = gm * mrest
        ek = en - mrest

! print the progress
!
        write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')        &
                                                    n, t, dt, tg, va, ek, term

! write results to the output file
!
        write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), p(1), p(2), p(3),  &
                                       va, vp, vr, gm, en, ek, bunit * ba,     &
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
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, va, ek

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
    real(kind=8), dimension(3,3) :: si, ai
    real(kind=8), dimension(3)   :: xi, pi, xm, pm, u, b
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

! get acceleration for the current state
!
        call acceleration(ti(i), xi(1:3), pi(1:3), si(i,1:3), ai(i,1:3), u(:), b(:))

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

! import required modules
!
    use coordinates, only : is_inside

    implicit none

! local variables
!
    logical                        :: keepon = .true.
    character(len=32)              :: str
    integer                        :: n, m, i = 0, mi, ti, k
    real(kind=8), dimension(4,6)   :: z
    real(kind=8), dimension(5,4,6) :: zp
    real(kind=8), dimension(3)     :: x, v, p, s, a, u, b
    real(kind=8), dimension(3)     :: xc, xe, xs
    real(kind=8), dimension(3)     :: pc, pe, ps
    real(kind=8)                   :: gm, t, dt, tc, te, ts
    real(kind=8)                   :: ba, va, vp, vr, om, tg, rg, en, ek
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
    p(:) = p0(:)

! calculate the Lorentz factor and particle speed
!
    gm = lorentz_factor(p(:))
    v(:) = p(:) / gm

! calculate the acceleration at the initial position
!
    call acceleration(t, x(:), p(:), s(:), a(:), u(:), b(:))

! separate the particle velocity into parallel and perpendicular components
!
    call separate_velocity(v(:), b(:), ba, va, vp, vr)

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, vr, om, tg, rg)

! calculate the particle energies
!
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')            &
                                                    n, t, dt, tg, va, ek, term

! open the output file, print headers and the initial values
!
    open(10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,20a22)") 'Time', 'X', 'Y', 'Z', 'Px', 'Py', 'Pz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]', '<B> [Gs]',  &
                                 'Omega [1/s]', 'Tg [s]', 'Rg [m]', 'Tg [T]',  &
                                 'Rg [L]', 'Tolerance', 'Iterations'
    write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), p(1), p(2), p(3),      &
                                   va, vp, vr, gm, en, ek, bunit * ba,         &
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

! if the boundaries are not periodic and particle is out of the box, stop
! the integration
!
      keepon = keepon .and. is_inside(x)

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
        v(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(:), p(:), s(:), a(:), u(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(v(:), b(:), ba, va, vp, vr)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, vr, om, tg, rg)

! calculate particle energy
!
        en = gm * mrest
        ek = en - mrest

! print the progress
!
        write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')        &
                                                    n, t, dt, tg, va, ek, term

! write results to the output file
!
        write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), p(1), p(2), p(3),  &
                                       va, vp, vr, gm, en, ek, bunit * ba,     &
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
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, va, ek

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

! import required modules
!
    use coordinates, only : is_inside

    implicit none

! local variables
!
    logical                        :: keepon = .true.
    character(len=32)              :: str
    integer                        :: n, m, i = 0, mi, ti, k
    real(kind=8), dimension(4,6)   :: z
    real(kind=8), dimension(5,4,6) :: zp
    real(kind=8), dimension(5)     :: hp
    real(kind=8), dimension(3)     :: x, v, p, s, a, u, b
    real(kind=8), dimension(3)     :: xc, xe, xs
    real(kind=8), dimension(3)     :: pc, pe, ps
    real(kind=8)                   :: gm, t, dt, dtp, tc, te, ts
    real(kind=8)                   :: ba, va, vp, vr, om, tg, rg, en, ek
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
    p(:) = p0(:)

! calculate the Lorentz factor and particle speed
!
    gm = lorentz_factor(p(:))
    v(:) = p(:) / gm

! calculate the acceleration at the initial position
!
    call acceleration(t, x(:), p(:), s(:), a(:), u(:), b(:))

! separate the particle velocity into parallel and perpendicular components
!
    call separate_velocity(v(:), b(:), ba, va, vp, vr)

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, vr, om, tg, rg)

! calculate the particle energies
!
    en = gm * mrest
    ek = en - mrest

! print the progress
!
    write(*,"('PROGRESS  : ',a8,2x,5(a14))") 'ITER', 'TIME', 'TIMESTEP',       &
                                             'GPERIOD', 'SPEED (c)',           &
                                             'ENERGY (MeV)'
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')            &
                                                    n, t, dt, tg, va, ek, term

! open the output file, print headers and the initial values
!
    open(10, file = 'output.dat', form = 'formatted', status = 'replace')
    write(10,"('#',1a20,20a22)") 'Time', 'X', 'Y', 'Z', 'Px', 'Py', 'Pz',      &
                                 '|V| [c]', '|Vpar| [c]', '|Vper| [c]',        &
                                 'gamma', 'En [MeV]', 'Ek [MeV]', '<B> [Gs]',  &
                                 'Omega [1/s]', 'Tg [s]', 'Rg [m]', 'Tg [T]',  &
                                 'Rg [L]', 'Tolerance', 'Iterations'
    write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), p(1), p(2), p(3),      &
                                   va, vp, vr, gm, en, ek, bunit * ba,         &
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

! if the boundaries are not periodic and particle is out of the box, stop
! the integration
!
      keepon = keepon .and. is_inside(x)

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
        v(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(:), p(:), s(:), a(:), u(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(v(:), b(:), ba, va, vp, vr)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, vr, om, tg, rg)

! calculate particle energy
!
        en = gm * mrest
        ek = en - mrest

! print the progress
!
        write(*,"('PROGRESS  : ',i8,2x,5(1es14.6),a1)", advance = 'no')        &
                                                    n, t, dt, tg, va, ek, term

! write results to the output file
!
        write(10,"(20(1es22.14),i22)") t, x(1), x(2), x(3), p(1), p(2), p(3),  &
                                       va, vp, vr, gm, en, ek, bunit * ba,     &
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
    write(*,"('PROGRESS  : ',i8,2x,5(1es14.6))") n, t, dt, tg, va, ek

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
    real(kind=8), dimension(4,3) :: si, ai
    real(kind=8), dimension(3)   :: xi, pi, xm, pm, u, b
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

! get acceleration for the current state
!
        call acceleration(ti(i), xi(1:3), pi(1:3), si(i,1:3), ai(i,1:3), u(:), b(:))

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
! subroutine ACCELERATION:
! -----------------------
!
!   Subroutine calculates the acceleration terms for a given state.
!
!   Arguments:
!
!     t    - the time;
!     x, p - the particle state (position and momentum);
!     s, a - the acceleration terms;
!     u, b - the interpolated values of plasma velocity and magnetic field;
!
!===============================================================================
!
  subroutine acceleration(t, x, p, s, a, u, b)

    implicit none

! subroutine arguments
!
    real(kind=8)              , intent(in)  :: t
    real(kind=8), dimension(3), intent(in)  :: x, p
    real(kind=8), dimension(3), intent(out) :: s, a, u, b

! local variables
!
    real(kind=8), dimension(3) :: v, w
#ifdef ITEST
    real(kind=8)               :: dl, ra, rb, xt, yt, rt
#endif /* ITEST */
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
! get plasma field components
!
    call fields(x(:), u(:), b(:))
#endif /* TEST */

! get the particle speed from its momentum
!
    v(:) = p(:) / lorentz_factor(p(:))

! normalize the speed to get position in the code units
!
    s(:) = v(:) / vunit

! subtract the fluid velocity
!
    w(:) = v(:) - u(:)

! compute the acceleration
!
    a(1) = qom * (w(2) * b(3) - w(3) * b(2))
    a(2) = qom * (w(3) * b(1) - w(1) * b(3))
    a(3) = qom * (w(1) * b(2) - w(2) * b(1))

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
! subroutine LORENTZ_FACTOR:
! -------------------------
!
!   Function calculates the Lorentz factor from the particle moment.
!
!   Arguments:
!
!     p - the particle momentum vector;
!
!===============================================================================
!
  function lorentz_factor(p) result(gm)

    implicit none

! subroutine arguments
!
    real(kind=8), dimension(3), intent(in) :: p
    real(kind=8)                           :: gm
!
!------------------------------------------------------------------------------
!
    gm = sqrt(1.0d+00 + dot_product(p(:), p(:)))

!-------------------------------------------------------------------------------
!
  end function lorentz_factor
!
!===============================================================================
!
! subroutine SEPARATE_VELOCITY:
! ----------------------------
!
!   Subroutine separates velocity into components parallel and perpendicular
!   to the local field.
!
!   Arguments:
!
!     v  - the particle velocity vector;
!     b  - the local magnetic field vector;
!     ba - the local magnetic field strength;
!     vm - the particle velocity magnitude;
!     vp - the particle velocity component parallel to the local field;
!     vr - the particle velocity component perpendicular to the local field;
!
!===============================================================================
!
  subroutine separate_velocity(v, b, ba, vm, vp, vr)

    implicit none

! subroutine arguments
!
    real(kind=8), dimension(3), intent(in)  :: v, b
    real(kind=8)              , intent(out) :: ba, vm, vp, vr

! local variables
!
    real(kind=8), dimension(3) :: u
    real(kind=8)               :: uu, vv
!
!------------------------------------------------------------------------------
!
! calculate amplitude of magnetic field
!
    ba = sqrt(dot_product(b, b))

! calculate component parallel to B
!
    if (ba > 0.0d+00) then
      u(:) = b(:) / ba
    else
      u(:) = 0.0d+00
    end if

! calculate component parallel to B
!
    uu = dot_product(v(:), u(:))**2

! calculate amplitude of velocity
!
    vv = dot_product(v(:), v(:))

! calculate amplitude and the parallel and perpendicular components of velocity
!
    vm = sqrt(vv)
    vp = sqrt(uu)
    vr = sqrt(vv - uu)

!-------------------------------------------------------------------------------
!
  end subroutine separate_velocity
!
!===============================================================================
!
! subroutine GYRO_PARAMETERS:
! --------------------------
!
!   Subroutine calculates particle gyro parameters.
!
!   Arguments:
!
!     gm - the particle Lorentz factor;
!     ba - the local magnetic field strength;
!     vr - the particle velocity component perpendicular to the local field;
!     om - the particle gyrofrequency (in radians per time unit);
!     tg - the particle gyroperiod (in time units);
!     rg - the particle gyroradius (in length units);
!
!===============================================================================
!
  subroutine gyro_parameters(gm, ba, vr, om, tg, rg)

    implicit none

! subroutine arguments
!
    real(kind=8), intent(in)  :: gm, ba, vr
    real(kind=8), intent(out) :: om, tg, rg
!
!------------------------------------------------------------------------------

    om = abs(qom) * ba / gm
    tg = pi2 / om
    rg = vr / (vunit * om)

!-------------------------------------------------------------------------------
!
  end subroutine gyro_parameters

!===============================================================================
!
end module particles
