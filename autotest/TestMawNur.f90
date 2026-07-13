module TestMawNur

  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE
  use MathUtilModule, only: is_close
  use testdrive, only: check, error_type, new_unittest, unittest_type
  use MawModule, only: maw_damp_weight

  implicit none
  private
  public :: collect_mawnur

  ! default damping parameters, matching the values used in maw_nur
  real(DP), parameter :: damptheta = 0.7_DP
  real(DP), parameter :: damptol = 0.5_DP
  real(DP), parameter :: weightmin = 0.01_DP
  real(DP), parameter :: recover = 0.2_DP

contains

  subroutine collect_mawnur(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("same_direction_recovers", test_same_direction), &
                new_unittest("decaying_reversal_recovers", &
                             test_decaying_reversal), &
                new_unittest("sustained_reversal_damps", &
                             test_sustained_reversal), &
                new_unittest("weight_floor", test_weight_floor), &
                new_unittest("recover_caps_at_one", test_recover_cap), &
                new_unittest("first_iteration_no_damp", test_first_iteration) &
                ]
  end subroutine collect_mawnur

  !> @brief A head change in the same direction as the last is not damped.
  !<
  subroutine test_same_direction(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: w
    ! same sign change, weight grows back by recover
    w = maw_damp_weight(8.0_DP, 3.0_DP, 0.7_DP, damptheta, damptol, &
                        weightmin, recover)
    call check(error, is_close(w, 0.9_DP))
  end subroutine test_same_direction

  !> @brief A reversal that is already shrinking is not damped.
  !<
  subroutine test_decaying_reversal(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: w
    ! reversal but |dxprop| = 1 < damptol * |dxold| = 0.5 * 8 = 4, so recover
    w = maw_damp_weight(-1.0_DP, 8.0_DP, 1.0_DP, damptheta, damptol, &
                        weightmin, recover)
    call check(error, is_close(w, DONE))
  end subroutine test_decaying_reversal

  !> @brief A reversal that is not shrinking (a real oscillation) is damped.
  !<
  subroutine test_sustained_reversal(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: w
    ! reversal and |dxprop| = 53 > damptol * |dxold| = 0.5 * 56 = 28, so damp
    w = maw_damp_weight(-53.0_DP, 56.0_DP, 1.0_DP, damptheta, damptol, &
                        weightmin, recover)
    call check(error, is_close(w, damptheta))
  end subroutine test_sustained_reversal

  !> @brief Damping never reduces the weight below weightmin.
  !<
  subroutine test_weight_floor(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: w
    ! damptheta * 0.011 = 0.0077, which is below the floor, so weightmin is used
    w = maw_damp_weight(-53.0_DP, 56.0_DP, 0.011_DP, damptheta, damptol, &
                        weightmin, recover)
    call check(error, is_close(w, weightmin))
  end subroutine test_weight_floor

  !> @brief Recovery never increases the weight above one.
  !<
  subroutine test_recover_cap(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: w
    ! 0.95 + recover = 1.15, which is capped at one
    w = maw_damp_weight(5.0_DP, 5.0_DP, 0.95_DP, damptheta, damptol, &
                        weightmin, recover)
    call check(error, is_close(w, DONE))
  end subroutine test_recover_cap

  !> @brief On the first iteration there is no previous change, so nothing is
  !! treated as a reversal and the weight is not damped.
  !<
  subroutine test_first_iteration(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: w
    ! dxold = 0, so dxprop * dxold = 0 (not a reversal) and the weight recovers
    w = maw_damp_weight(-53.0_DP, DZERO, 1.0_DP, damptheta, damptol, &
                        weightmin, recover)
    call check(error, is_close(w, DONE))
  end subroutine test_first_iteration

end module TestMawNur
