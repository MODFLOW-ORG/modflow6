module LinearAlgebraUtilsModule
  use KindModule, only: DP, I4B

  implicit none
  private

  public :: eye, outer_product, cross_product

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

end module LinearAlgebraUtilsModule
