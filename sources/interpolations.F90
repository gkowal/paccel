!!******************************************************************************
!!
!!  This file is part of the PACCEL source code, a program to integrate
!!  test particle trajectories in fields obtained from Newtonian or
!!  relativistic magnetohydrodynamical simulations.
!!
!!  Copyright (C) 2016-2021 Grzegorz Kowal <grzegorz@amuncode.org>
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!******************************************************************************
!!
!! module: INTERPOLATIONS
!!
!!  This module provides subroutines to interpolate variables at the particle
!!  position.
!!
!!******************************************************************************
!
module interpolations

  implicit none

! interfaces for procedure pointers
!
  abstract interface
    subroutine prepare_iface(pr, pi, dr, cc)
      real(kind=8), dimension(3)  , intent(in)  :: pr
      integer     , dimension(4,3), intent(out) :: pi
      real(kind=8), dimension(3)  , intent(out) :: dr
      real(kind=8), dimension(4,3), intent(out) :: cc
    end subroutine
    function interpolate_iface(f, pi, dr, cc) result(q)
      real(kind=4), dimension(:,:,:), intent(in) :: f
      integer     , dimension(4,3)  , intent(in) :: pi
      real(kind=8), dimension(3)    , intent(in) :: dr
      real(kind=8), dimension(4,3)  , intent(in) :: cc
      real(kind=8)                               :: q
    end function
    subroutine coefficients_cubic_iface(x, c)
      real(kind=8)              , intent(in)  :: x
      real(kind=8), dimension(4), intent(out) :: c
    end subroutine
    function pcubic_iface(c, fk, fl, fr, fq) result(q)
      real(kind=8), dimension(4), intent(in)  :: c
      real(kind=4)              , intent(in)  :: fk, fl, fr, fq
      real(kind=8)                            :: q
    end function
  end interface

! pointers to interpolation methods
!
  procedure(prepare_iface)           , pointer, save :: prepare_interpolation => null()
  procedure(interpolate_iface)       , pointer, save :: interpolate           => null()
  procedure(coefficients_cubic_iface), pointer, save :: coefficients_cubic    => null()
  procedure(pcubic_iface)            , pointer, save :: pcubic                => null()

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_interpolations
  public :: prepare_interpolation, interpolate
!
!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_INTERPOLATIONS:
! ------------------------------------
!
!   Subroutine initializes the interpolation module.
!
!   Arguments:
!
!     verbose - indicates if it should print any messages;
!     status  - the return value; if it is 0 everything went successfully,
!               otherwise there was a problem;
!
!===============================================================================
!
  subroutine initialize_interpolations(verbose, status)

! import required modules
!
    use parameters, only : get_parameter

    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: status

! local variables
!
    character(len=32) :: method        = 'nearest'
    character(len=32) :: method_name   = 'nearest neighbor'
    character(len=32) :: tvd           = 'on'
!
!-------------------------------------------------------------------------------
!
    status = 0

! print info
!
    if (verbose) write(*,"('INFO',6x,': initializing interpolations')")

! get module parameters
!
    call get_parameter('interpolation', method)
    call get_parameter('tvd'          , tvd   )

! select the interpolation method
!
    select case(trim(method))
    case("nearest", "zero")
      method_name           = 'nearest neighbor'
      prepare_interpolation => prepare_nearest
      interpolate           => interpolate_nearest
    case("linear")
      method_name           = 'linear'
      prepare_interpolation => prepare_linear
      interpolate           => interpolate_linear
    case("cubic")
      method_name           = 'cubic'
      prepare_interpolation => prepare_cubic
      interpolate           => interpolate_cubic
    case default
      write(*,"('ERROR',5x,': ',a)") "unsupported interpolation method"
      status = 100
      return
    end select

! check TVD option
!
    select case(trim(tvd))
    case ("on", "ON", "t", "T", "y", "Y", "true", "TRUE", "yes", "YES")
      coefficients_cubic => coefficients_cubic_tvd
      pcubic             => pcubic_tvd
      tvd                = 'on'
    case default
      coefficients_cubic => coefficients_cubic_notvd
      pcubic             => pcubic_notvd
      tvd                = 'off'
    end select

! print information about the interpolation method
!
    if (verbose) then
      write(*,"('INFO',6x,': ', a, ' interpolation')") trim(method_name)
      write(*,"('INFO',6x,': TVD limiting is ' ,a)"  ) trim(tvd)
    end if

!-------------------------------------------------------------------------------
!
  end subroutine initialize_interpolations
!
!===============================================================================
!
! subroutine PREPARE_NEAREST:
! --------------------------
!
!   Subroutine prepares the nearest neighbor interpolation.
!
!   Arguments:
!
!     pr - the index position;
!     pi - the integer indices of cells involved in the interpolation;
!     dr - the position fraction between two inferior cells;
!     cc - the interpolation coefficients (not used);
!
!===============================================================================
!
  subroutine prepare_nearest(pr, pi, dr, cc)

    implicit none

! subroutine arguments
!
    real(kind=8), dimension(3)  , intent(in)  :: pr
    integer     , dimension(4,3), intent(out) :: pi
    real(kind=8), dimension(3)  , intent(out) :: dr
    real(kind=8), dimension(4,3), intent(out) :: cc
!
!-------------------------------------------------------------------------------
!
! calculate indices
!
    pi(1,1) = nint(pr(1))
    pi(1,2) = nint(pr(2))
#if DIMS == 3
    pi(1,3) = nint(pr(3))
#else /* DIMS == 3 */
    pi(1,3) = 1
#endif /* DIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine prepare_nearest
!
!===============================================================================
!
! subroutine PREPARE_LINEAR:
! -------------------------
!
!   Subroutine prepares the linear interpolation.
!
!   Arguments:
!
!     pr - the index position;
!     pi - the integer indices of cells involved in the interpolation;
!     dr - the position fraction between two inferior cells;
!     cc - the interpolation coefficients (not used);
!
!===============================================================================
!
  subroutine prepare_linear(pr, pi, dr, cc)

    implicit none

! subroutine arguments
!
    real(kind=8), dimension(3)  , intent(in)  :: pr
    integer     , dimension(4,3), intent(out) :: pi
    real(kind=8), dimension(3)  , intent(out) :: dr
    real(kind=8), dimension(4,3), intent(out) :: cc
!
!-------------------------------------------------------------------------------
!
! calculate indices
!
    pi(1,1) = floor(pr(1))
    pi(2,1) = pi(1,1) + 1
    pi(1,2) = floor(pr(2))
    pi(2,2) = pi(1,2) + 1
#if DIMS == 3
    pi(1,3) = floor(pr(3))
    pi(2,3) = pi(1,3) + 1
#endif /* DIMS == 3 */

! calculate intercell position
!
    dr(1) = pr(1) - pi(1,1)
    dr(2) = pr(2) - pi(1,2)
#if DIMS == 3
    dr(3) = pr(3) - pi(1,3)
#endif /* DIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine prepare_linear
!
!===============================================================================
!
! subroutine PREPARE_CUBIC:
! -------------------------
!
!   Subroutine prepares the cubic interpolation.
!
!   Arguments:
!
!     pr - the index position;
!     pi - the integer indices of cells involved in the interpolation;
!     dr - the position fraction between two inferior cells;
!     cc - the interpolation coefficients;
!
!===============================================================================
!
  subroutine prepare_cubic(pr, pi, dr, cc)

    implicit none

! subroutine arguments
!
    real(kind=8), dimension(3)  , intent(in)  :: pr
    integer     , dimension(4,3), intent(out) :: pi
    real(kind=8), dimension(3)  , intent(out) :: dr
    real(kind=8), dimension(4,3), intent(out) :: cc
!
!-------------------------------------------------------------------------------
!
! calculate indices
!
    pi(2,1) = floor(pr(1))
    pi(1,1) = pi(2,1) - 1
    pi(3,1) = pi(2,1) + 1
    pi(4,1) = pi(2,1) + 2
    pi(2,2) = floor(pr(2))
    pi(1,2) = pi(2,2) - 1
    pi(3,2) = pi(2,2) + 1
    pi(4,2) = pi(2,2) + 2
#if DIMS == 3
    pi(2,3) = floor(pr(3))
    pi(1,3) = pi(2,3) - 1
    pi(3,3) = pi(2,3) + 1
    pi(4,3) = pi(2,3) + 2
#endif /* DIMS == 3 */

! calculate intercell position
!
    dr(1) = pr(1) - pi(2,1)
    dr(2) = pr(2) - pi(2,2)
#if DIMS == 3
    dr(3) = pr(3) - pi(2,3)
#endif /* DIMS == 3 */

! coefficients for dx, dy, and dz
!
    call coefficients_cubic(dr(1), cc(:,1))
    call coefficients_cubic(dr(2), cc(:,2))
#if DIMS == 3
    call coefficients_cubic(dr(3), cc(:,3))
#endif /* DIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine prepare_cubic
!
!===============================================================================
!
! subroutine INTERPOLATE_NEAREST:
! ------------------------------
!
!   Subroutine interpolates the field using the prepared coefficients for
!   the nearest neighbor method.
!
!   Arguments:
!
!     f  - the interpolated field;
!     pi - the integer indices of cells involved in the interpolation;
!     dr - the position fraction between two inferior cells;
!     cc - the interpolation coefficients (not used);
!
!===============================================================================
!
  function interpolate_nearest(f, pi, dr, cc) result(q)

    implicit none

! subroutine arguments
!
    real(kind=4), dimension(:,:,:), intent(in) :: f
    integer     , dimension(4,3)  , intent(in) :: pi
    real(kind=8), dimension(3)    , intent(in) :: dr
    real(kind=8), dimension(4,3)  , intent(in) :: cc
    real(kind=8)                               :: q
!
!-------------------------------------------------------------------------------
!
    q = f(pi(1,1),pi(1,2),pi(1,3))

!-------------------------------------------------------------------------------
!
  end function interpolate_nearest
!
!===============================================================================
!
! subroutine INTERPOLATE_LINEAR:
! -----------------------------
!
!   Subroutine interpolates the field using the prepared coefficients for
!   the linear method.
!
!   Arguments:
!
!     f  - the interpolated field;
!     pi - the integer indices of cells involved in the interpolation;
!     dr - the position fraction between two inferior cells;
!     cc - the interpolation coefficients (not used);
!
!===============================================================================
!
  function interpolate_linear(f, pi, dr, cc) result(q)

    implicit none

! subroutine arguments
!
    real(kind=4), dimension(:,:,:), intent(in) :: f
    integer     , dimension(4,3)  , intent(in) :: pi
    real(kind=8), dimension(3)    , intent(in) :: dr
    real(kind=8), dimension(4,3)  , intent(in) :: cc
    real(kind=8)                               :: q

! local variables
!
    real(kind=4) :: q11, q12, q21, q22, q1, q2
!
!-------------------------------------------------------------------------------
!
#if DIMS == 2
!= bilinear interpolation =
!
! interpolate along the Y direction
!
    q1  = plinear(dr(2), f(pi(1,1),pi(1,2),1), f(pi(1,1),pi(2,2),1))
    q2  = plinear(dr(2), f(pi(2,1),pi(1,2),1), f(pi(2,1),pi(2,2),1))
#else /* DIMS == 2 */
!= trilinear interpolation =
!
! interpolate along the Z direction
!
    q11 = plinear(dr(3), f(pi(1,1),pi(1,2),pi(1,3)), f(pi(1,1),pi(1,2),pi(2,3)))
    q12 = plinear(dr(3), f(pi(1,1),pi(2,2),pi(1,3)), f(pi(1,1),pi(2,2),pi(2,3)))
    q21 = plinear(dr(3), f(pi(2,1),pi(1,2),pi(1,3)), f(pi(2,1),pi(1,2),pi(2,3)))
    q22 = plinear(dr(3), f(pi(2,1),pi(2,2),pi(1,3)), f(pi(2,1),pi(2,2),pi(2,3)))

! interpolate along the Y direction
!
    q1 = plinear(dr(2), q11, q12)
    q2 = plinear(dr(2), q21, q22)
#endif /* DIMS == 2 */

! interpolate along the X direction
!
    q  = plinear(dr(1), q1 , q2 )

!-------------------------------------------------------------------------------
!
  end function interpolate_linear
!
!===============================================================================
!
! subroutine INTERPOLATE_CUBIC:
! ----------------------------
!
!   Subroutine interpolates the field using the prepared coefficients for
!   the cubic method.
!
!   Arguments:
!
!     f  - the interpolated field;
!     pi - the integer indices of cells involved in the interpolation;
!     dr - the position fraction between two inferior cells;
!     cc - the interpolation coefficients (not used);
!
!===============================================================================
!
  function interpolate_cubic(f, pi, dr, cc) result(q)

    implicit none

! subroutine arguments
!
    real(kind=4), dimension(:,:,:), intent(in) :: f
    integer     , dimension(4,3)  , intent(in) :: pi
    real(kind=8), dimension(3)    , intent(in) :: dr
    real(kind=8), dimension(4,3)  , intent(in) :: cc
    real(kind=8)                               :: q

! local variables
!
    real(kind=4) :: q11, q12, q13, q14, q21, q22, q23, q24                     &
                  , q31, q32, q33, q34, q41, q42, q43, q44, q1, q2, q3, q4
!
!-------------------------------------------------------------------------------
!
#if DIMS == 2
!= bicubic interpolation =
!
! interpolate along the Y direction
!
    q1  = pcubic(cc(:,2), f(pi(1,1),pi(1,2),1), f(pi(1,1),pi(2,2),1)           &
                   , f(pi(1,1),pi(3,2),1), f(pi(1,1),pi(4,2),1))
    q2  = pcubic(cc(:,2), f(pi(2,1),pi(1,2),1), f(pi(2,1),pi(2,2),1)           &
                   , f(pi(2,1),pi(3,2),1), f(pi(2,1),pi(4,2),1))
    q3  = pcubic(cc(:,2), f(pi(3,1),pi(1,2),1), f(pi(3,1),pi(2,2),1)           &
                   , f(pi(3,1),pi(3,2),1), f(pi(3,1),pi(4,2),1))
    q4  = pcubic(cc(:,2), f(pi(4,1),pi(1,2),1), f(pi(4,1),pi(2,2),1)           &
                   , f(pi(4,1),pi(3,2),1), f(pi(4,1),pi(4,2),1))
#else /* DIMS == 2 */
!= tricubic interpolation =
!
! interpolate along the Z direction
!
    q11 = pcubic(cc(:,3), f(pi(1,1),pi(1,2),pi(1,3)), f(pi(1,1),pi(1,2),pi(2,3))&
                        , f(pi(1,1),pi(1,2),pi(3,3)), f(pi(1,1),pi(1,2),pi(4,3)))
    q12 = pcubic(cc(:,3), f(pi(1,1),pi(2,2),pi(1,3)), f(pi(1,1),pi(2,2),pi(2,3))&
                        , f(pi(1,1),pi(2,2),pi(3,3)), f(pi(1,1),pi(2,2),pi(4,3)))
    q13 = pcubic(cc(:,3), f(pi(1,1),pi(3,2),pi(1,3)), f(pi(1,1),pi(3,2),pi(2,3))&
                        , f(pi(1,1),pi(3,2),pi(3,3)), f(pi(1,1),pi(3,2),pi(4,3)))
    q14 = pcubic(cc(:,3), f(pi(1,1),pi(4,2),pi(1,3)), f(pi(1,1),pi(4,2),pi(2,3))&
                        , f(pi(1,1),pi(4,2),pi(3,3)), f(pi(1,1),pi(4,2),pi(4,3)))

    q21 = pcubic(cc(:,3), f(pi(2,1),pi(1,2),pi(1,3)), f(pi(2,1),pi(1,2),pi(2,3))&
                        , f(pi(2,1),pi(1,2),pi(3,3)), f(pi(2,1),pi(1,2),pi(4,3)))
    q22 = pcubic(cc(:,3), f(pi(2,1),pi(2,2),pi(1,3)), f(pi(2,1),pi(2,2),pi(2,3))&
                        , f(pi(2,1),pi(2,2),pi(3,3)), f(pi(2,1),pi(2,2),pi(4,3)))
    q23 = pcubic(cc(:,3), f(pi(2,1),pi(3,2),pi(1,3)), f(pi(2,1),pi(3,2),pi(2,3))&
                        , f(pi(2,1),pi(3,2),pi(3,3)), f(pi(2,1),pi(3,2),pi(4,3)))
    q24 = pcubic(cc(:,3), f(pi(2,1),pi(4,2),pi(1,3)), f(pi(2,1),pi(4,2),pi(2,3))&
                        , f(pi(2,1),pi(4,2),pi(3,3)), f(pi(2,1),pi(4,2),pi(4,3)))

    q31 = pcubic(cc(:,3), f(pi(3,1),pi(1,2),pi(1,3)), f(pi(3,1),pi(1,2),pi(2,3))&
                        , f(pi(3,1),pi(1,2),pi(3,3)), f(pi(3,1),pi(1,2),pi(4,3)))
    q32 = pcubic(cc(:,3), f(pi(3,1),pi(2,2),pi(1,3)), f(pi(3,1),pi(2,2),pi(2,3))&
                        , f(pi(3,1),pi(2,2),pi(3,3)), f(pi(3,1),pi(2,2),pi(4,3)))
    q33 = pcubic(cc(:,3), f(pi(3,1),pi(3,2),pi(1,3)), f(pi(3,1),pi(3,2),pi(2,3))&
                        , f(pi(3,1),pi(3,2),pi(3,3)), f(pi(3,1),pi(3,2),pi(4,3)))
    q34 = pcubic(cc(:,3), f(pi(3,1),pi(4,2),pi(1,3)), f(pi(3,1),pi(4,2),pi(2,3))&
                        , f(pi(3,1),pi(4,2),pi(3,3)), f(pi(3,1),pi(4,2),pi(4,3)))

    q41 = pcubic(cc(:,3), f(pi(4,1),pi(1,2),pi(1,3)), f(pi(4,1),pi(1,2),pi(2,3))&
                        , f(pi(4,1),pi(1,2),pi(3,3)), f(pi(4,1),pi(1,2),pi(4,3)))
    q42 = pcubic(cc(:,3), f(pi(4,1),pi(2,2),pi(1,3)), f(pi(4,1),pi(2,2),pi(2,3))&
                        , f(pi(4,1),pi(2,2),pi(3,3)), f(pi(4,1),pi(2,2),pi(4,3)))
    q43 = pcubic(cc(:,3), f(pi(4,1),pi(3,2),pi(1,3)), f(pi(4,1),pi(3,2),pi(2,3))&
                        , f(pi(4,1),pi(3,2),pi(3,3)), f(pi(4,1),pi(3,2),pi(4,3)))
    q44 = pcubic(cc(:,3), f(pi(4,1),pi(4,2),pi(1,3)), f(pi(4,1),pi(4,2),pi(2,3))&
                        , f(pi(4,1),pi(4,2),pi(3,3)), f(pi(4,1),pi(4,2),pi(4,3)))

! interpolate along the Y direction
!
    q1  = pcubic(cc(:,2), q11, q12, q13, q14)
    q2  = pcubic(cc(:,2), q21, q22, q23, q24)
    q3  = pcubic(cc(:,2), q31, q32, q33, q34)
    q4  = pcubic(cc(:,2), q41, q42, q43, q44)
#endif /* DIMS == 2 */

! interpolate along the X direction
!
    q   = pcubic(cc(:,1), q1 , q2 , q3 , q4 )

!-------------------------------------------------------------------------------
!
  end function interpolate_cubic
!
!===============================================================================
!
! subroutine PLINEAR:
! ------------------
!
!   Subroutine performs 1D linear interpolation.
!
!   Arguments:
!
!     fx     - the position fraction between two inferior cells;
!     fl, fr - the left and right states;
!
!===============================================================================
!
  function plinear(fx, fl, fr) result(q)

    implicit none

! subroutine arguments
!
    real(kind=8), intent(in) :: fx
    real(kind=4), intent(in) :: fl, fr
    real(kind=8)             :: q
!
!-------------------------------------------------------------------------------
!
    q = fl + fx * (fr - fl)

!-------------------------------------------------------------------------------
!
  end function plinear
!
!===============================================================================
!
! subroutine PCUBIC_NOTVD:
! -----------------------
!
!   Subroutine performs 1D not limited cubic interpolation.
!
!   Arguments:
!
!     c              - the cubic interpolation coefficients;
!     fk, fl, fr, fq - the variable states;
!
!===============================================================================
!
  function pcubic_notvd(c, fk, fl, fr, fq) result(q)

    implicit none

! input and output arguments
!
    real(kind=8), dimension(4), intent(in)  :: c
    real(kind=4)              , intent(in)  :: fk, fl, fr, fq
    real(kind=8)                            :: q
!
!-------------------------------------------------------------------------------
!
    q = c(1) * fk + c(2) * fl + c(3) * fr + c(4) * fq

!-------------------------------------------------------------------------------
!
  end function pcubic_notvd
!
!===============================================================================
!
! subroutine PCUBIC_TVD:
! ---------------------
!
!   Subroutine performs 1D TVD limited cubic interpolation.
!
!   Arguments:
!
!     c              - the cubic interpolation coefficients;
!     fk, fl, fr, fq - the variable states;
!
!===============================================================================
!
  function pcubic_tvd(c, fk, fl, fr, fq) result(q)

    implicit none

! input and output arguments
!
    real(kind=8), dimension(4), intent(in)  :: c
    real(kind=4)              , intent(in)  :: fk, fl, fr, fq
    real(kind=8)                            :: q

! local parameters
!
    real(kind=8) :: dfl, dfr, ds, dl, dr
!
!-------------------------------------------------------------------------------
!
! calculate left, middle and right differences
!
    ds  = fr - fl
    dl  = fl - fk
    dr  = fq - fr

! calculate the limited derivatives
!
    dfl = sign(1.0d+00, ds) * min(abs(ds), abs(dl))
    dfr = sign(1.0d+00, ds) * min(abs(ds), abs(dr))

! check monotonicity
!
    if ((dl * ds) <= 0.0d+00) then
      dfl = 0.0d+00
    end if
    if ((dr * ds) <= 0.0d+00) then
      dfr = 0.0d+00
    end if

! perform the final interpolation
!
    q = c(1) * fl + c(2) * fr + c(3) * dfl + c(4) * dfr

!-------------------------------------------------------------------------------
!
  end function pcubic_tvd
!
!===============================================================================
!
! subroutine COEFFICIENTS_CUBIC_NOTVD:
! -----------------------------------
!
!   Subroutine calculates cubic interpolation coefficients without limiting.
!
!   Arguments:
!
!     x - the position fraction between cells;
!     c - the cubic interpolation coefficients;
!
!===============================================================================
!
  subroutine coefficients_cubic_notvd(x, c)

    implicit none

! input and output arguments
!
    real(kind=8)              , intent(in)  :: x
    real(kind=8), dimension(4), intent(out) :: c

! local variables
!
    real(kind=8) :: x1, x2
!
!-------------------------------------------------------------------------------
!
! prepare local variables
!
    x1 = x - 1.0d+00
    x2 = x * x

! calculate coefficients
!
    c(1) = 0.5d+00 * x * ( ( 2.0d+00 - x ) * x - 1.0d+00 )
    c(2) = 0.5d+00 * x2 * ( 3.0d+00 * x - 5.0d+00 ) + 1.0d+00
    c(3) = 0.5d+00 * x * ( ( 4.0d+00 - 3.0d+00 * x ) * x + 1.0d+00 )
    c(4) = 0.5d+00 * x1 * x2

!-------------------------------------------------------------------------------
!
  end subroutine coefficients_cubic_notvd
!
!===============================================================================
!
! subroutine COEFFICIENTS_CUBIC_TVD:
! ---------------------------------
!
!   Subroutine calculates cubic interpolation coefficients with TVD limiting.
!
!   Arguments:
!
!     x - the position fraction between cells;
!     c - the cubic interpolation coefficients;
!
!===============================================================================
!
  subroutine coefficients_cubic_tvd(x, c)

    implicit none

! input and output arguments
!
    real(kind=8)              , intent(in)  :: x
    real(kind=8), dimension(4), intent(out) :: c

! local variables
!
    real(kind=8) :: x1, x2, x3, xd
!
!-------------------------------------------------------------------------------
!
! prepare local variables
!
    x1 = x - 1.0d+00
    x2 = x * x
    x3 = x1 * x1
    xd = 2.0d+00 * x

! calculate coefficients
!
    c(1) = x3 * (xd + 1.0d+00)
    c(2) = x2 * (3.0d+00 - xd)
    c(3) = x  * x3
    c(4) = x2 * x1

!-------------------------------------------------------------------------------
!
  end subroutine coefficients_cubic_tvd

!===============================================================================
!
end module interpolations
