module LinearAlgebraUtilsModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE

  implicit none
  private

  public :: eye, outer_product, cross_product, solve_tridiagonal

contains

  pure function eye(n) result(A)
    ! -- dummy
    integer(I4B), intent(in) :: n
    real(DP), dimension(:, :), allocatable :: A
    ! -- locals
    integer(I4B) :: i

    allocate (A(n, n))
    A = 0.0_DP
    do i = 1, n
      A(i, i) = 1.0_DP
    end do
  end function eye

  pure function outer_product(A, B) result(AB)
    ! -- dummy
    real(DP), intent(in) :: A(:), B(:)
    real(DP) :: AB(size(A), size(B))
    ! -- locals
    integer :: nA, nB

    nA = size(A)
    nB = size(B)

    AB = spread(source=A, dim=2, ncopies=nB) * &
         spread(source=B, dim=1, ncopies=nA)
  end function outer_product

  function cross_product(a, b) result(c)
    ! -- return
    real(DP), dimension(3) :: c
    ! -- dummy
    real(DP), dimension(3), intent(in) :: a
    real(DP), dimension(3), intent(in) :: b

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product

  !> @brief Solve a tridiagonal system with the Thomas algorithm
  !<
  subroutine solve_tridiagonal(n, tl, td, tu, b, x, w)
    ! -- dummy variables
    integer(I4B), intent(in) :: n !< number of matrix rows
    real(DP), dimension(n), intent(in) :: tl !< lower matrix terms
    real(DP), dimension(n), intent(in) :: td !< diagonal matrix terms
    real(DP), dimension(n), intent(in) :: tu !< upper matrix terms
    real(DP), dimension(n), intent(in) :: b !< right-hand side vector
    real(DP), dimension(n), intent(inout) :: x !< solution vector
    real(DP), dimension(n), intent(inout) :: w !< work vector
    ! -- local variables
    integer(I4B) :: j
    real(DP) :: bet
    real(DP) :: beti
    !
    ! -- initialize variables
    w(1) = DZERO
    bet = td(1)
    beti = DONE / bet
    x(1) = b(1) * beti
    !
    ! -- decomposition and forward substitution
    do j = 2, n
      w(j) = tu(j - 1) * beti
      bet = td(j) - tl(j) * w(j)
      beti = DONE / bet
      x(j) = (b(j) - tl(j) * x(j - 1)) * beti
    end do
    !
    ! -- backsubstitution
    do j = n - 1, 1, -1
      x(j) = x(j) - w(j + 1) * x(j + 1)
    end do
  end subroutine solve_tridiagonal

end module LinearAlgebraUtilsModule
