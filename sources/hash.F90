!!******************************************************************************
!!
!!  This file is part of the PACCEL source code, a program to integrate
!!  test particle trajectories in fields obtained from Newtonian or
!!  relativistic magnetohydrodynamical simulations.
!!
!!  Copyright (C) 2012-2020 Yann Collet
!!  Copyright (C) 2020-2021 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: HASH
!!
!!  This module provides 64-bit version of the xxHash64 by Yann Collet.
!!  This is a Fortran implementation based on the XXH64 specification
!!  published at
!!      https://github.com/Cyan4973/xxHash/blob/dev/doc/xxhash_spec.md
!!
!!  For additional info, see
!!      http://www.xxhash.com or https://github.com/Cyan4973/xxHash
!!
!!******************************************************************************
!
module hash

! module variables are not implicit by default
!
  implicit none

! hash parameters
!
  integer(kind=8), parameter :: seed   = 0_8
  integer(kind=8), parameter :: prime1 = -7046029288634856825_8,               &
                                prime2 = -4417276706812531889_8,               &
                                prime3 =  1609587929392839161_8,               &
                                prime4 = -8796714831421723037_8,               &
                                prime5 =  2870177450012600261_8

! by default everything is private
!
  private

! declare public subroutines
!
  public :: xxh64

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
!
!===============================================================================
!
! function XXH64:
! --------------
!
!   Function calculates XXH64 hash for a given sequence of bytes.
!
!   Arguments:
!
!     input - the input sequence of bytes;
!
!===============================================================================
!
  integer(kind=8) function xxh64(input) result(hash)

    implicit none

! subroutine arguments
!
    integer(kind=1), dimension(:), intent(in) :: input

! local variables
!
    integer(kind=8) :: length, remaining, offset

! local arrays
!
    integer(kind=8), dimension(4) :: lane, chunk

!-------------------------------------------------------------------------------
!
    length    = size(input)
    hash      = 0_8
    offset    = 1_8
    remaining = length

    if (remaining >= 32_8) then
      lane(1) = seed + prime1 + prime2
      lane(2) = seed + prime2
      lane(3) = seed + 0_8
      lane(4) = seed - prime1

      do while (remaining >= 32_8)
        chunk(1:4) = transfer(input(offset:offset+31), 1_8, 4)

        lane(1) = xxh64_round(lane(1), chunk(1))
        lane(2) = xxh64_round(lane(2), chunk(2))
        lane(3) = xxh64_round(lane(3), chunk(3))
        lane(4) = xxh64_round(lane(4), chunk(4))

        offset    = offset    + 32_8
        remaining = remaining - 32_8
      end do

      hash = xxh64_rotl(lane(1),  1) + xxh64_rotl(lane(2),  7) +               &
             xxh64_rotl(lane(3), 12) + xxh64_rotl(lane(4), 18)

      hash = xxh64_merge(hash, lane(1))
      hash = xxh64_merge(hash, lane(2))
      hash = xxh64_merge(hash, lane(3))
      hash = xxh64_merge(hash, lane(4))

    else
      hash = seed + prime5
    end if

    hash = hash + length

    do while (remaining >= 8_8)
      chunk(1) = transfer(input(offset:offset+7), 1_8)
      hash = ieor(hash, xxh64_round(0_8, chunk(1)))
      hash = xxh64_rotl(hash, 27)
      hash = hash * prime1 + prime4

      offset    = offset    + 8_8
      remaining = remaining - 8_8
    end do

    if (remaining >= 4_8) then
      chunk(1) = transfer((/ input(offset:offset+3), 0_1, 0_1, 0_1, 0_1 /), 1_8)
      hash = ieor(hash, chunk(1) * prime1)
      hash = xxh64_rotl(hash, 23)
      hash = hash * prime2 + prime3

      offset    = offset    + 4_8
      remaining = remaining - 4_8
    end if

    do while (remaining > 0_8)
      chunk(1) = transfer((/ input(offset), 0_1, 0_1, 0_1,                     &
                                            0_1, 0_1, 0_1, 0_1 /), 1_8)
      hash = ieor(hash, chunk(1) * prime5)
      hash = xxh64_rotl(hash, 11)
      hash = hash * prime1

      offset    = offset    + 1_8
      remaining = remaining - 1_8
    end do

    hash = xxh64_aval(hash)

    return

!-------------------------------------------------------------------------------
!
  end function xxh64
!
!===============================================================================
!
! function XXH64_ROUND:
! --------------------
!
!   Function processes one stripe of the input data updating
!   the correponding lane.
!
!   Arguments:
!
!     lane  - the lane;
!     input - the 8-byte data to process;
!
!===============================================================================
!
  integer(kind=8) function xxh64_round(lane, input)

    implicit none

! subroutine arguments
!
    integer(kind=8), intent(in) :: lane, input

!-------------------------------------------------------------------------------
!
    xxh64_round = lane + (input * prime2)
    xxh64_round = xxh64_rotl(xxh64_round, 31)
    xxh64_round = xxh64_round * prime1
    return

!-------------------------------------------------------------------------------
!
  end function xxh64_round
!
!===============================================================================
!
! function XXH64_MERGE:
! --------------------
!
!   Function performs merging of the given lane in to the hash.
!
!   Arguments:
!
!     hash - the hash to merge to;
!     lane - the lane being merged;
!
!===============================================================================
!
  integer(kind=8) function xxh64_merge(hash, lane)

    implicit none

! subroutine arguments
!
    integer(kind=8), intent(in) :: hash, lane

!-------------------------------------------------------------------------------
!
    xxh64_merge = ieor(hash, xxh64_round(0_8, lane))
    xxh64_merge = xxh64_merge * prime1 + prime4
    return

!-------------------------------------------------------------------------------
!
  end function xxh64_merge
!
!===============================================================================
!
! function XXH64_AVAL:
! -------------------
!
!   Function calculates the final mix of the hash.
!
!   Arguments:
!
!     hash   - the hash to mix;
!
!===============================================================================
!
  integer(kind=8) function xxh64_aval(hash)

    implicit none

! subroutine arguments
!
    integer(kind=8), intent(in) :: hash

!-------------------------------------------------------------------------------
!
    xxh64_aval = hash
    xxh64_aval = ieor(xxh64_aval, ishft(xxh64_aval, -33)) * prime2
    xxh64_aval = ieor(xxh64_aval, ishft(xxh64_aval, -29)) * prime3
    xxh64_aval = ieor(xxh64_aval, ishft(xxh64_aval, -32))
    return

!-------------------------------------------------------------------------------
!
  end function xxh64_aval
!
!===============================================================================
!
! function XXH64_ROTL:
! -------------------
!
!   Function calculates the rotation of the input 8-byte word by a given amount.
!
!   Arguments:
!
!     byte   - the byte to be rotates;
!     amount - the amount by which rotate the input byte;
!
!===============================================================================
!
  integer(kind=8) function xxh64_rotl(byte, amount)

    implicit none

! subroutine arguments
!
    integer(kind=8), intent(in) :: byte
    integer(kind=4), intent(in) :: amount

!-------------------------------------------------------------------------------
!
    xxh64_rotl = ior(ishft(byte, amount), ishft(byte, amount - 64))
    return

!-------------------------------------------------------------------------------
!
  end function xxh64_rotl

!===============================================================================
!
end module hash
