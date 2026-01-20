!> @brief This module defines variable data types
!!
!! This module defines the precision of logical, integer, long integer,
!! and real data types used in MODFLOW 6 and are available to all other
!! modules.
!!
!<
module KindModule
  use, intrinsic :: iso_fortran_env, only: I4B => int32, &
                                         &I8B => int64, &
                                         &LGP => int32, &
                                         &DP => real64

  implicit none

  public :: I4B, I8B, LGP, DP, write_kindinfo

contains

  !> @brief Write variable data types
    !!
    !! This subroutine writes the precision of logical, integer, long integer,
    !! and real data types used in MODFLOW 6.
    !!
  !<
  subroutine write_kindinfo(iout)
    ! -- dummy variables
    integer(I4B), intent(in) :: iout !< file unit to output kind variables
    ! -- local variables
    integer(LGP) :: ldum = 0
    integer(I4B) :: idum = 0
    integer(I8B) :: long_idum = 0
    integer(DP) :: irdum = 0 ! for bit size of real variables
    real(DP) :: rdum = 0._DP
    !
    ! -- write kind information
    write (iout, '(/a)') 'Real Variables'
    write (iout, '(2x,a,i0)') 'KIND: ', DP
    write (iout, '(2x,a,1pg15.6)') 'TINY (smallest non-zero value): ', &
      tiny(rdum)
    write (iout, '(2x,a,1pg15.6)') 'HUGE (largest value): ', huge(rdum)
    write (iout, '(2x,a,i0)') 'PRECISION: ', precision(rdum)
    write (iout, '(2x,a,i0)') 'SIZE IN BITS: ', bit_size(irdum)

    write (iout, '(/a)') 'Integer Variables'
    write (iout, '(2x,a,i0)') 'KIND: ', I4B
    write (iout, '(2x,a,i0)') 'HUGE (largest value): ', huge(idum)
    write (iout, '(2x,a,i0)') 'SIZE IN BITS: ', bit_size(idum)

    write (iout, '(/a)') 'Long Integer Variables'
    write (iout, '(2x,a,i0)') 'KIND: ', I8B
    write (iout, '(2x,a,i0)') 'HUGE (largest value): ', huge(long_idum)
    write (iout, '(2x,a,i0)') 'SIZE IN BITS: ', bit_size(long_idum)

    write (iout, '(/a)') 'Logical Variables'
    write (iout, '(2x,a,i0)') 'KIND: ', LGP
    write (iout, '(2x,a,i0)') 'SIZE IN BITS: ', bit_size(ldum)
  end subroutine write_kindinfo

end module KindModule
