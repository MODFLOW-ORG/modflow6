module TestSVD
    use KindModule, only: I4B, DP
    use ConstantsModule, only: DSAME
    use testdrive, only: error_type, unittest_type, new_unittest, check, test_failed
    use SVDModule, only: SVD2, bidiagonal_decomposition, bidiagonal_qr_decomposition
    use MathUtilModule, only: eye
    use TesterUtils, only: check_same_matrix, check_matrix_row_orthogonal, check_matrix_column_orthogonal, &
      check_matrix_is_diagonal, check_matrix_bidiagonal
    
    implicit none
    private
    public :: collect_svd

    contains

  subroutine collect_svd(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("bidiagonal_decomposition", test_bidiagonal_decomposition),  &
                new_unittest("diagonalize_matrix", test_diagonalize_matrix),  &
                new_unittest("svd_input_bidiagonal", test_svd_input_bidiagonal),  &
                new_unittest("svd_input_square_matrix", test_svd_input_square_matrix),  &
                new_unittest("svd_input_M_GT_N_matrix", test_svd_input_M_GT_N_matrix),  &
                new_unittest("svd_input_N_GT_M_matrix", test_svd_input_N_GT_M_matrix),  &
                new_unittest("svd_zero_on_last_diagonal", test_svd_zero_on_last_diagonal),  &
                new_unittest("svd_zero_on_inner_diagonal", test_svd_zero_on_inner_diagonal)  &
                ]
  end subroutine collect_svd

  !> Test the bidiagonal decomposition
  !!
  !! This test checks that the bidiagonal decomposition is correct.
  !! It does so by checking that the matrices P and Qt are orthogonal,
  !! that the bidiagonal matrix is correct and that the original matrix
  !! can be reconstructed.
  !!
  !! The matrices A1 A2 and A3 are chosen such that all diagonal and non-diagonal
  !! elements are non-zero.
  subroutine test_bidiagonal_decomposition(error)
    type(error_type), allocatable, intent(out) :: error
    ! - locals
    real(DP), dimension(4,3) :: A1, A1_reconstructed, A1_mod
    real(DP), dimension(3,4) :: A2, A2_reconstructed, A2_mod
    real(DP), dimension(4,4) :: A3, A3_reconstructed, A3_mod

    REAL(DP), DIMENSION(:,:), allocatable :: P1, Qt1
    REAL(DP), DIMENSION(:,:), allocatable :: P2, Qt2
    REAL(DP), DIMENSION(:,:), allocatable :: P3, Qt3

    ! - Arrange.
    A1 = reshape( &
      [ 1.0_DP, 0.0_DP, 1.0_DP, 1.0_DP, &
        0.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, &
        0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP &
      ], [4,3])
    A1_mod = A1

    A2 = reshape( &
      [ 1.0_DP, 0.0_DP, 0.0_DP, &
        0.0_DP, 1.0_DP, 1.0_DP, &
        1.0_DP, 1.0_DP, 1.0_DP, &
        1.0_DP, 1.0_DP, 0.0_DP &
      ], [3,4])
      A2_mod = A2

      A3 = reshape( &
      [ 1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, &
        0.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, &
        1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, &
        1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP &
      ], [4,4])
      A3_mod = A3

    ! - Act.
    call bidiagonal_decomposition(A1_mod, P1, Qt1)
    call bidiagonal_decomposition(A2_mod, P2, Qt2)
    call bidiagonal_decomposition(A3_mod, P3, Qt3)

    ! - Assert.
    ! Test A1
    ! A1_reconstructed = P1 * A1_mod * Qt1
    A1_reconstructed = MATMUL(P1, MATMUL(A1_mod, Qt1))
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
    A2_reconstructed = MATMUL(P2, MATMUL(A2_mod, Qt2))
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
    A3_reconstructed = MATMUL(P3, MATMUL(A3_mod, Qt3))
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

  subroutine test_diagonalize_matrix(error)
    type(error_type), allocatable, intent(out) :: error
    ! - locals
    real(DP), dimension(4,4) :: A1, A1_reconstructed, A1_mod
    real(DP), dimension(3,4) :: A2, A2_reconstructed, A2_mod
    real(DP), dimension(4,3) :: A3, A3_reconstructed, A3_mod
    REAL(DP), DIMENSION(:,:), allocatable :: U1, Vt1
    REAL(DP), DIMENSION(:,:), allocatable :: U2, Vt2
    REAL(DP), DIMENSION(:,:), allocatable :: U3, Vt3

    ! - Arrange.
    A1 = reshape( &
      [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
        1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, &
        0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP, &
        0.0_DP, 0.0_DP, 1.0_DP, 1.0_DP &
      ], [4,4])
    A1_mod = A1
    U1 = Eye(4)
    Vt1 = Eye(4)

    A2 = reshape( &
    [ 1.0_DP, 0.0_DP, 0.0_DP, &
      1.0_DP, 1.0_DP, 0.0_DP, &
      0.0_DP, 1.0_DP, 1.0_DP, &
      0.0_DP, 0.0_DP, 1.0_DP &
    ], [3,4])
    A2_mod = A2
    U2 = Eye(3)
    Vt2 = Eye(4)

    A3 = reshape( &
    [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
      1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, &
      0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP &
    ], [4,3])
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
    A1_reconstructed = MATMUL(U1, MATMUL(A1_mod, Vt1))
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
    A2_reconstructed = MATMUL(U2, MATMUL(A2_mod, Vt2))
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
    A3_reconstructed = MATMUL(U3, MATMUL(A3_mod, Vt3))
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

  subroutine test_svd_input_bidiagonal(error)
    type(error_type), allocatable, intent(out) :: error
    ! - locals
    real(DP), dimension(4,4) :: A, A_reconstructed
    REAL(DP), DIMENSION(:,:), allocatable :: U, S, Vt

    ! - Arrange.
    A = reshape( &
    [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
      1.0_DP, 2.0_DP, 0.0_DP, 0.0_DP, &
      0.0_DP, 1.0_DP, 3.0_DP, 0.0_DP, &
      0.0_DP, 0.0_DP, 1.0_DP, 4.0_DP &
    ], [4,4])

    ! - Act.
    call SVD2(A, U, S, Vt)

    ! - Assert.
    A_reconstructed = MATMUL(U, MATMUL(S, Vt))
    call check_same_matrix(error, A_reconstructed, A)
    if (allocated(error)) return

  end subroutine test_svd_input_bidiagonal

  !> Test the SVD decomposition of a square matrix
  !!
  !! This test checks that the SVD decomposition of a square matrix is correct.
  !! The input matrix is choosen such that the returned S-matrix
  !! has elements on the entire diagonal.
  subroutine test_svd_input_square_matrix(error)
    type(error_type), allocatable, intent(out) :: error
    ! - locals
    real(DP), dimension(4,4) :: A, A_reconstructed
    REAL(DP), DIMENSION(:,:), allocatable :: U, S, Vt

    ! - Arrange.
    A = reshape( &
    [ 1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, &
      0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, &
      1.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, &
      1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP &
    ], [4,4])

    ! - Act.
    call SVD2(A, U, S, Vt)

    ! - Assert.
    A_reconstructed = MATMUL(U, MATMUL(S, Vt))
    call check_same_matrix(error, A_reconstructed, A)
    if (allocated(error)) return

    call check_matrix_is_diagonal(error, S)
    if (allocated(error)) return

  end subroutine test_svd_input_square_matrix

  !> Test the SVD decomposition of a square matrix
  !!
  !! This test checks that the SVD decomposition of a square matrix is correct.
  !! The input matrix is choosen such that the returned S-matrix
  !! has elements on the entire diagonal.
  subroutine test_svd_input_M_GT_N_matrix(error)
    type(error_type), allocatable, intent(out) :: error
    ! - locals
    real(DP), dimension(4,3) :: A, A_reconstructed
    REAL(DP), DIMENSION(:,:), allocatable :: U, S, Vt

    ! - Arrange.
    A = reshape( &
    [ 1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, &
      0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, &
      1.0_DP, 0.0_DP, 1.0_DP, 1.0_DP &
    ], [4,3])

    ! - Act.
    call SVD2(A, U, S, Vt)

    ! - Assert.
    A_reconstructed = MATMUL(U, MATMUL(S, Vt))
    call check_same_matrix(error, A_reconstructed, A)
    if (allocated(error)) return

    call check_matrix_is_diagonal(error, S)
    if (allocated(error)) return

  end subroutine test_svd_input_M_GT_N_matrix 

  !> Test the SVD decomposition of a square matrix
  !!
  !! Test zero on last diagonal element
  subroutine test_svd_input_N_GT_M_matrix(error)
    type(error_type), allocatable, intent(out) :: error
    ! - locals
    real(DP), dimension(3,5) :: A, A_reconstructed
    REAL(DP), DIMENSION(:,:), allocatable :: U, S, Vt

    ! - Arrange.
    A = reshape( &
    [ 1.0_DP, 0.0_DP, 1.0_DP, &
      0.0_DP, 1.0_DP, 0.0_DP, &
      1.0_DP, 0.0_DP, 0.0_DP, &
      1.0_DP, 1.0_DP, 1.0_DP, &
      1.0_DP, 1.0_DP, 0.0_DP &
    ], [3,5])

    ! - Act.
    call SVD2(A, U, S, Vt)

    ! - Assert.
    A_reconstructed = MATMUL(U, MATMUL(S, Vt))
    call check_same_matrix(error, A_reconstructed, A)
    if (allocated(error)) return

    call check_matrix_is_diagonal(error, S)
    if (allocated(error)) return
  end subroutine test_svd_input_N_GT_M_matrix


  !> Test the SVD decomposition of a square matrix
  !!
  !! Test zero on last diagonal element
  subroutine test_svd_zero_on_last_diagonal(error)
    type(error_type), allocatable, intent(out) :: error
    ! - locals
    real(DP), dimension(4,4) :: A, A_reconstructed
    REAL(DP), DIMENSION(:,:), allocatable :: U, S, Vt

    A = reshape( &
    [ 1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, &
      1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, &
      0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP, &
      0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP &
    ], [4,4])
    ! - Act.
    call SVD2(A, U, S, Vt)

    ! - Assert.
    A_reconstructed = MATMUL(U, MATMUL(S, Vt))
    call check_same_matrix(error, A_reconstructed, A)
    if (allocated(error)) return

    call check_matrix_is_diagonal(error, S)
    if (allocated(error)) return

  end subroutine test_svd_zero_on_last_diagonal


    !> Test the SVD decomposition of a square matrix
  !!
  !! Test zero on last diagonal element
  subroutine test_svd_zero_on_inner_diagonal(error)
    type(error_type), allocatable, intent(out) :: error
    ! - locals
    real(DP), dimension(6,6) :: A, A_reconstructed
    REAL(DP), DIMENSION(:,:), allocatable :: U, S, Vt

    A = reshape( &
    [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
      1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
      0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
      0.0_DP, 0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, &
      0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP, &
      0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 1.0_DP &
    ], [6,6])
    ! - Act.
    call SVD2(A, U, S, Vt)

    ! - Assert.
    A_reconstructed = MATMUL(U, MATMUL(S, Vt))
    call check_same_matrix(error, A_reconstructed, A)
    if (allocated(error)) return

    call check_matrix_is_diagonal(error, S)
    if (allocated(error)) return

  end subroutine test_svd_zero_on_inner_diagonal

end module TestSVD