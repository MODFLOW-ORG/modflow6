module PseudoInverseModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DPREC

  use SVDModule, only: SVD

  implicit none
  private

  public :: pinv

contains

  !!> @brief Computes the Moore-Penrose pseudoinverse of a matrix using SVD.
  !!
  !! This function computes the pseudoinverse \(A^+\) of the input matrix \(A\) using
  !! the singular value decomposition (SVD). Small singular values (less than `DPREC`)
  !! are set to zero to improve numerical stability. The pseudoinverse is useful for
  !! solving least-squares problems and for handling rank-deficient or non-square matrices.
  !!
  !! The result satisfies:
  !!     B = V * Sigma^+ * U^T
  !! where \(A = U \Sigma V^T\) is the SVD of \(A\), and \(\Sigma^+\) is the diagonal matrix
  !! with reciprocals of the nonzero singular values.
  !<
  function pinv(A) result(B)
    ! -- dummy
    real(DP), intent(in) :: A(:, :)
    real(DP) :: B(SIZE(A, DIM=2), SIZE(A, DIM=1)) !! The pseudoinverse of A (size n x m, where A is m x n)
    ! -- locals
    integer(I4B) :: m, n
    integer(I4B) :: pos
    real(DP), dimension(:, :), allocatable :: U
    real(DP), dimension(:, :), allocatable :: Vt
    real(DP), dimension(:, :), allocatable :: Sigma

    m = size(A, dim=1)
    n = size(A, dim=2)

    CALL SVD(A, U, Sigma, Vt)

    ! Transform Sigma to its Sigma^+
    do pos = 1, min(m, n)
      if (dabs(Sigma(pos, pos)) < DPREC) then
        Sigma(pos, pos) = 0.0_DP
      else
        Sigma(pos, pos) = 1.0_DP / Sigma(pos, pos)
      end if
    end do

    B = matmul(transpose(Vt), matmul(transpose(Sigma), transpose(U)))

  end function pinv

end module PseudoInverseModule
