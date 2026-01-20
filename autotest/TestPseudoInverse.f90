module TestPseudoInverse
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DSAME
  use testdrive, only: error_type, unittest_type, new_unittest, check, test_failed
  use PseudoInverseModule, only: pinv
  use LinearAlgebraUtilsModule, only: eye
  use TesterUtils, only: check_same_matrix

  implicit none
  private
  public :: collect_pinv

contains

  subroutine collect_pinv(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("identity", test_identity), &
                new_unittest("square", test_square), &
                new_unittest("tall", test_tall), &
                new_unittest("wide", test_wide), &
                new_unittest("rank_def", test_rank_def), &
                new_unittest("last_row_col_zero", test_last_row_col_zero), &
                new_unittest("zero", test_zero) &
                ]
  end subroutine collect_pinv

  !> @brief Test the pseudoinverse of the identity matrix.
  !!
  !! This test verifies that the pseudoinverse of a 3x3 identity matrix is itself.
  !<
  subroutine test_identity(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(3, 3) :: A, B_expected, B_actual

    A = eye(3) ! Create a 3x3 identity matrix
    B_expected = A
    B_actual = pinv(A)

    call check_same_matrix(error, B_actual, B_expected, tolerance=DSAME)
    if (allocated(error)) return
  end subroutine

  !> @brief Test the pseudoinverse of a square invertible matrix.
  !!
  !! This test verifies that the pseudoinverse of a regular 2x2 matrix matches its
  !! analytical inverse.
  !<
  subroutine test_square(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(2, 2) :: A, B_expected, B_actual

    A = reshape([ &
                4.0_DP, 7.0_DP, &
                2.0_DP, 6.0_DP &
                ], [2, 2])
    B_expected = reshape([ &
                         0.6_DP, -0.7_DP, &
                         -0.2_DP, 0.4_DP &
                         ], [2, 2])
    B_actual = pinv(A)

    call check_same_matrix(error, B_actual, B_expected, tolerance=DSAME)
    if (allocated(error)) return
  end subroutine

  !> @brief Test the pseudoinverse of a tall matrix (more rows than columns).
  !!
  !! This test checks that the pseudoinverse of a 3x2 matrix matches the expected result,
  !! as computed analytically.
  !<
  subroutine test_tall(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(3, 2) :: A
    real(DP), dimension(2, 3) :: B_expected, B_actual

    A = reshape([ &
                1.0_DP, 2.0_DP, 3.0_DP, &
                4.0_DP, 5.0_DP, 6.0_DP &
                ], [3, 2])

    B_expected = reshape([ &
                         -17.0_DP / 18.0_DP, 4.0_DP / 9.0_DP, &
                         -1.0_DP / 9.0_DP, 1.0_DP / 9.0_DP, &
                         13.0_DP / 18.0_DP, -2.0_DP / 9.0_DP &
                         ], [2, 3])

    B_actual = pinv(A)

    call check_same_matrix(error, B_actual, B_expected, tolerance=DSAME)
    if (allocated(error)) return
  end subroutine

  !> @brief Test the pseudoinverse of a wide matrix (more columns than rows).
  !!
  !! This test checks that the pseudoinverse of a 2x3 matrix matches the expected result,
  !! as computed analytically.
  !<
  subroutine test_wide(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(2, 3) :: A
    real(DP), dimension(3, 2) :: B_expected, B_actual

    A = reshape([ &
                1.0_DP, 2.0_DP, &
                3.0_DP, 4.0_DP, &
                5.0_DP, 6.0_DP &
                ], [2, 3])
    B_expected = reshape([ &
                         -4.0_DP / 3.0_DP, -1.0_DP / 3.0_DP, 2.0_DP / 3.0_DP, &
                         13.0_DP / 12.0_DP, 1.0_DP / 3.0_DP, -5.0_DP / 12.0_DP &
                         ], [3, 2])
    B_actual = pinv(A)

    call check_same_matrix(error, B_actual, B_expected, tolerance=DSAME)
    if (allocated(error)) return
  end subroutine

  !> @brief Test the pseudoinverse of a rank-deficient matrix.
  !!
  !! This test verifies that the pseudoinverse of a singular (rank-deficient) 2x2 matrix
  !! matches the expected result, ensuring correct handling of singular values.
  !! The test matrix is rank-deficient because its second row is a multiple of the first
  !! (row 2 = 2 × row 1), so its rank is 1 instead of 2. This checks that the pseudoinverse
  !! implementation correctly handles matrices that are not full rank.
  !<
  subroutine test_rank_def(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(2, 2) :: A, B_expected, B_actual

    A = reshape([ &
                1.0_DP, 2.0_DP, &
                2.0_DP, 4.0_DP &
                ], [2, 2])
    B_expected = reshape([ &
                         0.04_DP, 0.08_DP, &
                         0.08_DP, 0.16_DP], [2, 2])
    B_actual = pinv(A)

    call check_same_matrix(error, B_actual, B_expected, tolerance=DSAME)
    if (allocated(error)) return
  end subroutine

  !> @brief Test the pseudoinverse of a 3x3 matrix with the last row and column as zeros.
  !!
  !! This test a common case encountered with the TVD limiter in which the last row and column
  !! of the matrix are zero due to solving a 2D or 1D problem in a 3D context.
  !<
  subroutine test_last_row_col_zero(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(3, 3) :: A, B_expected, B_actual

    A = 0.0_DP
    A(1:2, 1:2) = reshape([4.0_DP, 7.0_DP, 2.0_DP, 6.0_DP], [2, 2])

    B_expected = 0.0_DP
    B_expected(1:2, 1:2) = reshape([0.6_DP, -0.7_DP, -0.2_DP, 0.4_DP], [2, 2])

    B_actual = pinv(A)

    call check_same_matrix(error, B_actual, B_expected, tolerance=DSAME)
    if (allocated(error)) return
  end subroutine test_last_row_col_zero

  !> @brief Test the pseudoinverse of a zero matrix.
  !!
  !! This test checks that the pseudoinverse of a 3x3 zero matrix is also a zero matrix.
  !<
  subroutine test_zero(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(3, 3) :: A, B_expected, B_actual

    A = 0.0_DP
    B_expected = 0.0_DP
    B_actual = pinv(A)

    call check_same_matrix(error, B_actual, B_expected, tolerance=DSAME)
    if (allocated(error)) return
  end subroutine

end module TestPseudoInverse
