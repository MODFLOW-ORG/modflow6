module TestSVD
  use KindModule, only: I4B, DP
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use SVDModule, only: SVD, bidiagonal_decomposition, bidiagonal_qr_decomposition
  use LinearAlgebraUtilsModule, only: eye
  use TesterUtils, only: check_same_matrix, &
                         check_matrix_row_orthogonal, &
                         check_matrix_column_orthogonal, &
                         check_matrix_is_diagonal, &
                         check_matrix_bidiagonal

  implicit none
  private
  public :: collect_svd

contains

  subroutine collect_svd(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("bidiagonal_decomposition", &
                             test_bidiagonal_decomposition), &
                new_unittest("diagonalize_matrix", &
                             test_diagonalize_matrix), &
                new_unittest("svd_input_bidiagonal", &
                             test_svd_input_bidiagonal), &
                new_unittest("svd_input_square_matrix", &
                             test_svd_input_square_matrix), &
                new_unittest("svd_input_tall_matrix", &
                             test_svd_input_tall_matrix), &
                new_unittest("svd_input_wide_matrix", &
                             test_svd_input_wide_matrix), &
                new_unittest("svd_zero_on_last_diagonal", &
                             test_svd_zero_on_last_diagonal), &
                new_unittest("svd_zero_on_inner_diagonal", &
                             test_svd_zero_on_inner_diagonal) &
                ]
  end subroutine collect_svd

  !> @brief Test the bidiagonal decomposition (first phase of SVD)
  !!
  !! This test verifies the correctness of the **bidiagonal decomposition** step,
  !! which is the first phase of the SVD algorithm. It checks that the input matrix
  !! can be reduced to bidiagonal form using orthogonal matrices P and Qt, and that
  !! these matrices are orthogonal. The test does **not** perform full diagonalization.
  !! The focus is on ensuring that the reduction to bidiagonal form is correct and
  !! reversible for various matrix shapes (square, tall, wide).
  !<
  subroutine test_bidiagonal_decomposition(error)
    ! -- dummy
    type(error_type), allocatable, intent(out) :: error
    ! -- locals
    real(DP), dimension(4, 3) :: A1, A1_reconstructed, A1_mod ! A1: more rows than columns (4x3)
    real(DP), dimension(3, 4) :: A2, A2_reconstructed, A2_mod ! A2: more columns than rows (3x4)
    real(DP), dimension(4, 4) :: A3, A3_reconstructed, A3_mod ! A3: square matrix (4x4)

    real(DP), dimension(:, :), allocatable :: P1, Qt1
    real(DP), dimension(:, :), allocatable :: P2, Qt2
    real(DP), dimension(:, :), allocatable :: P3, Qt3

    ! - Arrange.
    A1 = reshape( &
         [1.0_DP, 0.0_DP, 1.0_DP, 1.0_DP, &
          0.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, &
          0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP &
          ], [4, 3])
    A1_mod = A1

    A2 = reshape( &
         [1.0_DP, 0.0_DP, 0.0_DP, &
          0.0_DP, 1.0_DP, 1.0_DP, &
          1.0_DP, 1.0_DP, 1.0_DP, &
          1.0_DP, 1.0_DP, 0.0_DP &
          ], [3, 4])
    A2_mod = A2

    A3 = reshape( &
         [1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, &
          0.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, &
          1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, &
          1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP &
          ], [4, 4])
    A3_mod = A3

    ! - Act.
    call bidiagonal_decomposition(A1_mod, P1, Qt1)
    call bidiagonal_decomposition(A2_mod, P2, Qt2)
    call bidiagonal_decomposition(A3_mod, P3, Qt3)

    ! - Assert.
    ! Test A1
    ! A1_reconstructed = P1 * A1_mod * Qt1
    A1_reconstructed = matmul(P1, matmul(A1_mod, Qt1))
    call check_same_matrix(error, A1_reconstructed, A1)
    if (allocated(error)) return

    call check_matrix_row_orthogonal(error, P1)
    call check_matrix_column_orthogonal(error, P1)
    if (allocated(error)) return

    call check_matrix_row_orthogonal(error, Qt1)
    call check_matrix_column_orthogonal(error, Qt1)
    if (allocated(error)) return

    call check_matrix_bidiagonal(error, A1_mod)
    if (allocated(error)) return

    ! Test A2
    ! A2_reconstructed = P2 * A2_mod * Qt2
    A2_reconstructed = matmul(P2, matmul(A2_mod, Qt2))
    call check_same_matrix(error, A2_reconstructed, A2)
    if (allocated(error)) return

    call check_matrix_row_orthogonal(error, P2)
    call check_matrix_column_orthogonal(error, P2)
    if (allocated(error)) return

    call check_matrix_row_orthogonal(error, Qt2)
    call check_matrix_column_orthogonal(error, Qt2)
    if (allocated(error)) return

    call check_matrix_bidiagonal(error, A2_mod)
    if (allocated(error)) return

    ! Test A3
    ! A3_reconstructed = P3 * A3_mod * Qt3
    A3_reconstructed = matmul(P3, matmul(A3_mod, Qt3))
    call check_same_matrix(error, A3_reconstructed, A3)
    if (allocated(error)) return

    call check_matrix_row_orthogonal(error, P3)
    call check_matrix_column_orthogonal(error, P3)
    if (allocated(error)) return

    call check_matrix_row_orthogonal(error, Qt3)
    call check_matrix_column_orthogonal(error, Qt3)
    if (allocated(error)) return

    call check_matrix_bidiagonal(error, A3_mod)
    if (allocated(error)) return
  end subroutine test_bidiagonal_decomposition

  !> @brief Test the diagonalization of a bidiagonal matrix (second phase of SVD)
  !!
  !! This test verifies the correctness of the **diagonalization** step,
  !! which is the second phase of the SVD algorithm. It assumes the input matrix
  !! is already in bidiagonal form and applies the QR algorithm to diagonalize it.
  !! The test checks that the resulting matrix is diagonal, that the orthogonal
  !! matrices U and Vt remain orthogonal, and that the original matrix can be
  !! reconstructed. This test ensures the full SVD process (from bidiagonal to diagonal)
  !! is correct.
  !<
  subroutine test_diagonalize_matrix(error)
    ! -- dummy
    type(error_type), allocatable, intent(out) :: error
    ! -- locals
    real(DP), dimension(4, 4) :: A1, A1_reconstructed, A1_mod
    real(DP), dimension(3, 4) :: A2, A2_reconstructed, A2_mod
    real(DP), dimension(4, 3) :: A3, A3_reconstructed, A3_mod
    real(DP), dimension(:, :), allocatable :: U1, Vt1
    real(DP), dimension(:, :), allocatable :: U2, Vt2
    real(DP), dimension(:, :), allocatable :: U3, Vt3

    ! - Arrange.
    A1 = reshape( &
         [1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
          1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, &
          0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP, &
          0.0_DP, 0.0_DP, 1.0_DP, 1.0_DP &
          ], [4, 4])
    A1_mod = A1
    U1 = Eye(4)
    Vt1 = Eye(4)

    A2 = reshape( &
         [1.0_DP, 0.0_DP, 0.0_DP, &
          1.0_DP, 1.0_DP, 0.0_DP, &
          0.0_DP, 1.0_DP, 1.0_DP, &
          0.0_DP, 0.0_DP, 1.0_DP &
          ], [3, 4])
    A2_mod = A2
    U2 = Eye(3)
    Vt2 = Eye(4)

    A3 = reshape( &
         [1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
          1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, &
          0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP &
          ], [4, 3])
    A3_mod = A3
    U3 = Eye(4)
    Vt3 = Eye(3)

    ! - Act.
    call bidiagonal_qr_decomposition(A1_mod, U1, Vt1)
    call bidiagonal_qr_decomposition(A2_mod, U2, Vt2)
    call bidiagonal_qr_decomposition(A3_mod, U3, Vt3)

    ! - Assert.
    ! Test A1
    ! A1_reconstructed = U1 * A1_mod * Vt1
    A1_reconstructed = matmul(U1, matmul(A1_mod, Vt1))
    call check_same_matrix(error, A1_reconstructed, A1)
    if (allocated(error)) return

    call check_matrix_row_orthogonal(error, U1)
    call check_matrix_column_orthogonal(error, U1)
    if (allocated(error)) return

    call check_matrix_row_orthogonal(error, Vt1)
    call check_matrix_column_orthogonal(error, Vt1)
    if (allocated(error)) return

    call check_matrix_bidiagonal(error, A1_mod)
    if (allocated(error)) return

    ! Test A2
    ! A2_reconstructed = U2 * A2_mod * Vt2
    A2_reconstructed = matmul(U2, matmul(A2_mod, Vt2))
    call check_same_matrix(error, A2_reconstructed, A2)
    if (allocated(error)) return

    call check_matrix_row_orthogonal(error, U2)
    call check_matrix_column_orthogonal(error, U2)
    if (allocated(error)) return

    call check_matrix_row_orthogonal(error, Vt2)
    call check_matrix_column_orthogonal(error, Vt2)
    if (allocated(error)) return

    call check_matrix_bidiagonal(error, A2_mod)
    if (allocated(error)) return

    ! Test A3
    ! A3_reconstructed = U3 * A3_mod * Vt3
    A3_reconstructed = matmul(U3, matmul(A3_mod, Vt3))
    call check_same_matrix(error, A3_reconstructed, A3)
    if (allocated(error)) return

    call check_matrix_row_orthogonal(error, U3)
    call check_matrix_column_orthogonal(error, U3)
    if (allocated(error)) return

    call check_matrix_row_orthogonal(error, Vt3)
    call check_matrix_column_orthogonal(error, Vt3)
    if (allocated(error)) return

    call check_matrix_bidiagonal(error, A3_mod)
    if (allocated(error)) return

  end subroutine test_diagonalize_matrix

  !> @brief Test the SVD decomposition of a bidiagonal matrix
  !!
  !! This test checks that the SVD decomposition works correctly when the input matrix
  !! is already in bidiagonal form. It verifies that the SVD reconstructs the original
  !! matrix and that the resulting S matrix is diagonal. This ensures that the SVD
  !! implementation can handle bidiagonal matrices as input and produce correct results.
  !<
  subroutine test_svd_input_bidiagonal(error)
    ! -- dummy
    type(error_type), allocatable, intent(out) :: error
    ! -- locals
    real(DP), dimension(4, 4) :: A, A_reconstructed
    real(DP), DIMENSION(:, :), allocatable :: U, S, Vt

    ! - Arrange.
    A = reshape( &
        [1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
         1.0_DP, 2.0_DP, 0.0_DP, 0.0_DP, &
         0.0_DP, 1.0_DP, 3.0_DP, 0.0_DP, &
         0.0_DP, 0.0_DP, 1.0_DP, 4.0_DP &
         ], [4, 4])

    ! - Act.
    call SVD(A, U, S, Vt)

    ! - Assert.
    A_reconstructed = matmul(U, matmul(S, Vt))
    call check_same_matrix(error, A_reconstructed, A)
    if (allocated(error)) return

  end subroutine test_svd_input_bidiagonal

  !> @brief Test the SVD decomposition of a square matrix
  !!
  !! This test checks that the SVD decomposition of a square matrix is correct.
  !! The input matrix is chosen such that the returned S-matrix
  !! has elements on the entire diagonal.
  !<
  subroutine test_svd_input_square_matrix(error)
    ! -- dummy
    type(error_type), allocatable, intent(out) :: error
    ! -- locals
    real(DP), dimension(4, 4) :: A, A_reconstructed
    real(DP), DIMENSION(:, :), allocatable :: U, S, Vt

    ! - Arrange.
    A = reshape( &
        [1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, &
         0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, &
         1.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, &
         1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP &
         ], [4, 4])

    ! - Act.
    call SVD(A, U, S, Vt)

    ! - Assert.
    A_reconstructed = matmul(U, matmul(S, Vt))
    call check_same_matrix(error, A_reconstructed, A)
    if (allocated(error)) return

    call check_matrix_is_diagonal(error, S)
    if (allocated(error)) return

  end subroutine test_svd_input_square_matrix

  !> @brief Test the SVD decomposition of a tall matrix (m > n)
  !!
  !! This test checks that the SVD decomposition works correctly for a **tall matrix**,
  !! i.e., a matrix with more rows than columns (m > n). It verifies that the SVD
  !! reconstructs the original matrix, and that the resulting S matrix is diagonal.
  !! This ensures the SVD implementation handles tall matrices as expected.
  !<
  subroutine test_svd_input_tall_matrix(error)
    ! -- dummy
    type(error_type), allocatable, intent(out) :: error
    ! -- locals
    real(DP), dimension(4, 3) :: A, A_reconstructed
    real(DP), DIMENSION(:, :), allocatable :: U, S, Vt

    ! - Arrange.
    A = reshape( &
        [1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, &
         0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, &
         1.0_DP, 0.0_DP, 1.0_DP, 1.0_DP &
         ], [4, 3])

    ! - Act.
    call SVD(A, U, S, Vt)

    ! - Assert.
    A_reconstructed = matmul(U, matmul(S, Vt))
    call check_same_matrix(error, A_reconstructed, A)
    if (allocated(error)) return

    call check_matrix_is_diagonal(error, S)
    if (allocated(error)) return

  end subroutine test_svd_input_tall_matrix

  !> @brief Test the SVD decomposition of a wide matrix (n > m)
  !!
  !! This test checks that the SVD decomposition works correctly for a **wide matrix**,
  !! i.e., a matrix with more columns than rows (n > m). It verifies that the SVD
  !! reconstructs the original matrix, and that the resulting S matrix is diagonal.
  !! This ensures the SVD implementation handles wide matrices as expected.
  !<
  subroutine test_svd_input_wide_matrix(error)
    ! -- dummy
    type(error_type), allocatable, intent(out) :: error
    ! -- locals
    real(DP), dimension(3, 5) :: A, A_reconstructed
    real(DP), DIMENSION(:, :), allocatable :: U, S, Vt

    ! - Arrange.
    A = reshape( &
        [1.0_DP, 0.0_DP, 1.0_DP, &
         0.0_DP, 1.0_DP, 0.0_DP, &
         1.0_DP, 0.0_DP, 0.0_DP, &
         1.0_DP, 1.0_DP, 1.0_DP, &
         1.0_DP, 1.0_DP, 0.0_DP &
         ], [3, 5])

    ! - Act.
    call SVD(A, U, S, Vt)

    ! - Assert.
    A_reconstructed = matmul(U, matmul(S, Vt))
    call check_same_matrix(error, A_reconstructed, A)
    if (allocated(error)) return

    call check_matrix_is_diagonal(error, S)
    if (allocated(error)) return
  end subroutine test_svd_input_wide_matrix

  !> @brief Test the SVD decomposition of a matrix with a zero on the last diagonal element
  !!
  !! This test checks that the SVD decomposition works correctly when the input matrix
  !! has a zero on its last diagonal element. This situation can cause the matrix to be
  !! split into smaller blocks during the SVD process. The test verifies that the SVD
  !! reconstructs the original matrix and that the resulting S matrix is diagonal.
  !! This ensures the implementation handles matrices with trailing zero diagonals robustly.
  !<
  subroutine test_svd_zero_on_last_diagonal(error)
    ! -- dummy
    type(error_type), allocatable, intent(out) :: error
    ! -- locals
    real(DP), dimension(4, 4) :: A, A_reconstructed
    real(DP), DIMENSION(:, :), allocatable :: U, S, Vt

    A = reshape( &
        [1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, &
         1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, &
         0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP, &
         0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP &
         ], [4, 4])
    ! - Act.
    call SVD(A, U, S, Vt)

    ! - Assert.
    A_reconstructed = matmul(U, matmul(S, Vt))
    call check_same_matrix(error, A_reconstructed, A)
    if (allocated(error)) return

    call check_matrix_is_diagonal(error, S)
    if (allocated(error)) return

  end subroutine test_svd_zero_on_last_diagonal

  !> @brief Test the SVD decomposition of a matrix with a zero on an inner diagonal element
  !!
  !! This test checks that the SVD decomposition works correctly when the input matrix
  !! has a zero on an inner (not last) diagonal element. This case is more challenging,
  !! as it requires the SVD algorithm to correctly split and process the matrix into
  !! independent blocks. The test verifies that the SVD reconstructs the original matrix
  !! and that the resulting S matrix is diagonal. This ensures the implementation is robust
  !! for matrices with inner zero diagonals.
  !<
  subroutine test_svd_zero_on_inner_diagonal(error)
    ! -- dummy
    type(error_type), allocatable, intent(out) :: error
    ! -- locals
    real(DP), dimension(6, 6) :: A, A_reconstructed
    real(DP), DIMENSION(:, :), allocatable :: U, S, Vt

    A = reshape( &
        [1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
         1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
         0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
         0.0_DP, 0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, &
         0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP, &
         0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 1.0_DP &
         ], [6, 6])
    ! - Act.
    call SVD(A, U, S, Vt)

    ! - Assert.
    A_reconstructed = matmul(U, matmul(S, Vt))
    call check_same_matrix(error, A_reconstructed, A)
    if (allocated(error)) return

    call check_matrix_is_diagonal(error, S)
    if (allocated(error)) return

  end subroutine test_svd_zero_on_inner_diagonal

end module TestSVD
