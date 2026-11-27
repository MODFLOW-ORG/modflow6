module TestLinearAlgebraUtils
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE, DTWO
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       to_string, unittest_type
  use MathUtilModule, only: is_close
  use LinearAlgebraUtilsModule, only: solve_tridiagonal
  implicit none
  private
  public :: collect_linearalgebrautils

contains

  subroutine collect_linearalgebrautils(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("tridiagonal_laplacian", &
                             test_tridiagonal_laplacian), &
                new_unittest("tridiagonal_single", test_tridiagonal_single), &
                new_unittest("tridiagonal_nonsymmetric", &
                             test_tridiagonal_nonsymmetric) &
                ]
  end subroutine collect_linearalgebrautils

  !> @brief Solve the [-1, 2, -1] system for a known solution
  !<
  subroutine test_tridiagonal_laplacian(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: n = 5
    real(DP), dimension(n) :: tl, td, tu, b, x, w, xtrue
    integer(I4B) :: i

    do i = 1, n
      tl(i) = -DONE
      td(i) = DTWO
      tu(i) = -DONE
      xtrue(i) = real(i, DP)
    end do
    call tridiagonal_product(n, tl, td, tu, xtrue, b)

    call solve_tridiagonal(n, tl, td, tu, b, x, w)
    do i = 1, n
      call check(error, is_close(x(i), xtrue(i), rtol=1.0e-12_DP), &
                 "x("//to_string(i)//") = "//to_string(x(i)))
      if (allocated(error)) return
    end do
  end subroutine test_tridiagonal_laplacian

  !> @brief Solve a one-equation system
  !<
  subroutine test_tridiagonal_single(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(1) :: tl, td, tu, b, x, w

    tl(1) = DZERO
    td(1) = 4.0_DP
    tu(1) = DZERO
    b(1) = DTWO
    call solve_tridiagonal(1, tl, td, tu, b, x, w)
    call check(error, is_close(x(1), 0.5_DP), &
               "x(1) = "//to_string(x(1)))
  end subroutine test_tridiagonal_single

  !> @brief Solve a diagonally dominant nonsymmetric system
  !<
  subroutine test_tridiagonal_nonsymmetric(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: n = 20
    real(DP), dimension(n) :: tl, td, tu, b, x, w, xtrue
    integer(I4B) :: i

    do i = 1, n
      tl(i) = -0.3_DP - 0.01_DP * real(i, DP)
      tu(i) = -1.1_DP + 0.02_DP * real(i, DP)
      td(i) = 3.0_DP + 0.1_DP * real(mod(i, 3), DP)
      xtrue(i) = sin(real(i, DP))
    end do
    call tridiagonal_product(n, tl, td, tu, xtrue, b)

    call solve_tridiagonal(n, tl, td, tu, b, x, w)
    do i = 1, n
      call check(error, is_close(x(i), xtrue(i), atol=1.0e-12_DP), &
                 "x("//to_string(i)//") = "//to_string(x(i)))
      if (allocated(error)) return
    end do
  end subroutine test_tridiagonal_nonsymmetric

  !> @brief Multiply a tridiagonal matrix by a vector
  !<
  subroutine tridiagonal_product(n, tl, td, tu, x, b)
    integer(I4B), intent(in) :: n
    real(DP), dimension(n), intent(in) :: tl, td, tu, x
    real(DP), dimension(n), intent(out) :: b
    integer(I4B) :: i

    do i = 1, n
      b(i) = td(i) * x(i)
      if (i > 1) b(i) = b(i) + tl(i) * x(i - 1)
      if (i < n) b(i) = b(i) + tu(i) * x(i + 1)
    end do
  end subroutine tridiagonal_product

end module TestLinearAlgebraUtils
