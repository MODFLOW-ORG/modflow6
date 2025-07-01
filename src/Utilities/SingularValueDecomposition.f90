module SVDModule
  use KindModule, only: DP, LGP, I4B
  use ConstantsModule, only: DONE, DPREC, DSAME, DEM6
  use LinearAlgebraUtilsModule, only: eye, outer_product

  implicit none
  private

  public :: SVD
  public :: bidiagonal_decomposition
  public :: bidiagonal_qr_decomposition

contains

  !!> @brief Constructs a Householder reflection matrix for a given vector.
  !!
  !<
  pure function house_holder_matrix(x) result(Q)
    ! -- dummy
    real(DP), intent(in) :: x(:)
    real(DP), allocatable, dimension(:, :) :: Q
    ! -- locals
    real(DP) :: x_norm, y_norm
    real(DP), allocatable :: v(:)
    real(DP), allocatable :: w(:)

    x_norm = norm2(x)
    y_norm = norm2(x(2:))

    Q = Eye(size(x))

    if (dabs(y_norm) < DSAME) then
      return
    end if

    v = x
    if (x(1) > 0) then
      ! Parlett (1971) suggested this modification to avoid cancellation
      v(1) = -(y_norm**2) / (X(1) + x_norm)
    else
      v(1) = v(1) - x_norm
    end if

    w = v / norm2(v)

    Q = Q - 2.0_DP * outer_product(w, w)

  end function house_holder_matrix

  !!> @brief Reduces a matrix to bidiagonal form using Householder transformations.
  !!
  !! This subroutine transforms the input matrix `A` into a bidiagonal matrix by applying a sequence of
  !! Householder reflections from the left and right. The orthogonal transformations are accumulated in
  !! the matrices `P` (left transformations) and `Qt` (right transformations).
  !!
  !! The result is:
  !!     P^T * A * Qt = B
  !! where `B` is bidiagonal, `P` and `Qt` are orthogonal matrices.
  !!
  !! This step is the first phase of the Golub-Reinsch SVD algorithm and prepares the matrix for further
  !! diagonalization using the QR algorithm.
  !!
  !<
  pure subroutine bidiagonal_decomposition(A, P, Qt)
    ! -- dummy
    real(DP), intent(inout), dimension(:, :) :: A
    real(DP), intent(out), dimension(:, :), allocatable :: P, Qt
    ! -- locals
    integer(I4B) :: M, N
    integer(I4B) :: I
    real(DP), allocatable, dimension(:, :) :: Qi, Pi
    real(DP), allocatable, dimension(:) :: h

    M = size(A, dim=1) ! Number of rows
    N = size(A, dim=2) ! Number of columns

    Qt = Eye(N)
    P = Eye(M)

    do I = 1, min(M, N)
      ! columns
      h = A(I:M, I)
      Pi = house_holder_matrix(h)
      A(I:M, :) = matmul(Pi, A(I:M, :)) ! Apply householder transformation from left
      P(:, I:M) = matmul(P(:, I:M), Pi)

      ! rows
      if (I < N) then
        h = A(I, I + 1:N)
        Qi = transpose(house_holder_matrix(h))
        A(:, I + 1:N) = matmul(A(:, I + 1:N), Qi) ! Apply householder transformation from right
        Qt(I + 1:N, :) = matmul(Qi, Qt(I + 1:N, :))
      end if

    end do

  end subroutine bidiagonal_decomposition

  !!> @brief Constructs a Givens rotation matrix.
  !!
  !<
  pure function givens_rotation(a, b) result(G)
    ! -- dummy
    real(DP), intent(in) :: a, b
    real(DP), dimension(2, 2) :: G
    ! -- locals
    real(DP) :: c, s, h, d

    if (abs(b) < DPREC) then
      G = Eye(2)
      return
    end if

    h = hypot(a, b)
    d = 1.0 / h
    c = abs(a) * d
    s = sign(d, a) * b

    G(1, 1) = c
    G(1, 2) = s
    G(2, 1) = -s
    G(2, 2) = c

  end function givens_rotation

  !!> @brief Computes the Wilkinson shift for the lower-right 2x2 block of a bidiagonal matrix.
  !!
  !! This function calculates the Wilkinson shift, which is used to accelerate convergence in the implicit QR algorithm
  !! for bidiagonal matrices during SVD. The shift is computed from the lower-right 2x2 block of the input matrix `A`
  !! and is chosen to be the eigenvalue of that block which is closer to the bottom-right element.
  !!
  !! The shift is used to improve numerical stability and speed up the annihilation of off-diagonal elements.
  !!
  !<
  pure function compute_shift(A) result(mu)
    ! -- dummy
    real(DP), intent(in), dimension(:, :) :: A
    real(DP) :: mu
    ! -- locals
    integer(I4B) :: m, n, min_mn
    real(DP) T11, T12, T21, T22
    real(DP) dm, fmmin, fm, dn
    real(DP) :: mean, product, mean_product, mu1, mu2

    m = size(A, dim=1) ! Number of rows
    n = size(A, dim=2) ! Number of columns

    min_mn = min(m, n)

    dn = A(min_mn, min_mn)
    dm = A(min_mn - 1, min_mn - 1)
    fm = A(min_mn - 1, min_mn)
    if (min_mn > 2) then
      fmmin = A(min_mn - 2, min_mn - 1)
    else
      fmmin = 0.0_DP
    end if

    T11 = dm**2 + fmmin**2
    T12 = dm * fm
    T21 = T12
    T22 = fm**2 + dn**2

    mean = (T11 + T22) / 2.0_DP
    product = T11 * T22 - T12 * T21
    mean_product = mean**2 - product
    if (abs(mean_product) < DEM6) then
      mean_product = 0.0_DP
    end if
    mu1 = mean - sqrt(mean_product)
    mu2 = mean + sqrt(mean_product)
    if (abs(T22 - mu1) < abs(T22 - mu2)) then
      mu = mu1
    else
      mu = mu2
    end if

  end function compute_shift

  !!> @brief Performs QR iterations on a bidiagonal matrix as part of the SVD algorithm.
  !!
  !<
  pure subroutine bidiagonal_qr_decomposition(A, U, VT)
    ! -- dummy
    real(DP), intent(inout), dimension(:, :) :: A
    real(DP), intent(inout), dimension(:, :) :: U, Vt
    ! -- locals
    integer(I4B) :: m, n, I, J
    real(DP), dimension(2, 2) :: G
    real(DP) :: t11, t12, mu

    m = size(A, dim=1) ! Number of rows
    n = size(A, dim=2) ! Number of columns

    DO I = 1, n - 1
      J = I + 1
      if (I == 1) then
        ! For the first iteration use the implicit shift
        mu = compute_shift(A)
        ! Apply the shift to the full tri-diagonal matrix.
        t11 = A(1, 1)**2 - mu
        t12 = A(1, 1) * A(1, 2)

        G = givens_rotation(t11, t12)
      else
        G = givens_rotation(A(I - 1, I), A(I - 1, J))
      end if

      A(:, I:J) = matmul(A(:, I:J), transpose(G))
      Vt(I:J, :) = matmul(G, Vt(I:J, :))

      if (j > m) cycle
      G = givens_rotation(A(I, I), A(J, I))
      A(I:J, :) = matmul(G, A(I:J, :))
      U(:, I:J) = matmul(U(:, I:J), transpose(G))
    END DO

  end subroutine bidiagonal_qr_decomposition

  !!> @brief Handles zero (or near-zero) diagonal elements in a bidiagonal matrix during SVD.
  !!
  !! This subroutine detects and processes zero (or near-zero) diagonal elements in the bidiagonal matrix `A`.
  !! If a zero is found on the diagonal, it applies a sequence of Givens rotations to either zero out the
  !! corresponding column (if the zero is at the end of the diagonal) or the corresponding row (otherwise).
  !! The orthogonal transformations are accumulated in the `U` and `VT` matrices as appropriate.
  !! This step is essential for splitting the matrix into smaller blocks and ensuring numerical stability
  !! during the iterative QR step of the SVD algorithm.
  !!
  !<
  pure subroutine handle_zero_diagonal(A, U, VT)
    ! -- dummy
    real(DP), intent(inout), dimension(:, :) :: A
    real(DP), intent(inout), dimension(:, :) :: U, Vt
    ! -- locals
    integer(I4B) :: m, n, I
    real(DP), dimension(2, 2) :: G
    integer(I4B) :: zero_index

    m = size(A, dim=1) ! Number of rows
    n = size(A, dim=2) ! Number of columns

    do I = 1, min(m, n)
      if (abs(A(I, I)) < DPREC) then
        zero_index = I
        exit
      end if
    end do

    if (zero_index == min(n, m)) then
      ! If the zero index is the last element of the diagonal then zero out the column
      do I = zero_index - 1, 1, -1
        G = givens_rotation(A(I, I), A(I, zero_index))
        A(:, [I, zero_index]) = matmul(A(:, [I, zero_index]), transpose(G))
        Vt([I, zero_index], :) = matmul(G, Vt([I, zero_index], :))
      end do
    else
      ! Else zero out the row
      do I = zero_index + 1, n
        G = givens_rotation(A(I, I), A(zero_index, I))
        A([zero_index, I], :) = matmul(transpose(G), A([zero_index, I], :))
        U(:, [zero_index, I]) = matmul(U(:, [zero_index, I]), G)
      end do
    end if
  end subroutine handle_zero_diagonal

  !!> @brief Computes the infinity norm of the superdiagonal elements of a matrix.
  !!
  !<
  pure function superdiagonal_norm(A) result(norm)
    ! -- dummy
    real(DP), intent(in) :: A(:, :)
    real(DP) :: norm
    ! -- locals
    integer(I4B) :: m, n, I

    m = size(A, DIM=1) ! Number of rows
    n = size(A, DIM=2) ! Number of columns

    norm = 0.0_DP
    do I = 1, min(m, n) - 1
      norm = max(norm, abs(A(I, I + 1)))
    end do

  end function superdiagonal_norm

  !!> @brief Finds the range of nonzero superdiagonal elements in a bidiagonal matrix.
  !!
  !! This subroutine scans the superdiagonal of the matrix `A` and determines the indices `p` and `q`
  !! such that all superdiagonal elements outside the range `A(p:q, p:q)` are (near-)zero, while
  !! those within the range may be nonzero. This is useful for identifying active submatrices during
  !! the iterative QR step of the SVD algorithm.
  !!
  !! The search starts from the lower-right corner and moves upward to find the last nonzero
  !! superdiagonal (sets `q`), then moves upward again to find the first zero (sets `p`).
  !<
  pure subroutine find_nonzero_superdiagonal(A, p, q)
    ! -- dummy
    real(DP), intent(in), dimension(:, :) :: A
    integer(I4B), intent(out) :: p ! Index of the first nonzero superdiagonal (start of active block)
    integer(I4B), intent(out) :: q ! Index of the last nonzero superdiagonal (end of active block)
    ! -- locals
    integer(I4B) :: m, n, j, min_mn

    m = size(A, DIM=1) ! Number of rows
    n = size(A, DIM=2) ! Number of columns

    min_mn = min(m, n)
    p = 1
    q = min_mn

    do j = min_mn, 2, -1
      if (abs(A(j - 1, j)) > DPREC) then
        q = j
        exit
      end if
    end do

    do j = q - 1, 2, -1
      if (abs(A(j - 1, j)) < DPREC) then
        p = j
        exit
      end if
    end do
  end subroutine find_nonzero_superdiagonal

  !!> @brief Checks if a matrix has a (near-)zero diagonal element.
  !!
  !<
  pure function has_zero_diagonal(A) result(has_zero)
    ! -- dummy
    real(DP), intent(in) :: A(:, :)
    logical(LGP) :: has_zero
    ! -- locals
    integer(I4B) :: m, n, I

    m = size(A, dim=1) ! Number of rows
    n = size(A, dim=2) ! Number of columns

    has_zero = .FALSE.
    do I = 1, min(m, n)
      if (abs(A(I, I)) < DPREC) then
        has_zero = .TRUE.
        exit
      end if
    end do

  end function has_zero_diagonal

  !!> @brief Reduces a rectangular matrix to square form using Givens rotations.
  !!
  !! This subroutine transforms a rectangular matrix `A` into a square matrix by applying Givens rotations
  !! to eliminate extra columns (when the number of columns exceeds the number of rows).
  !! The corresponding orthogonal transformations are accumulated in the matrix `Qt`.
  !<
  pure subroutine make_matrix_square(A, Qt)
    ! -- dummy
    real(DP), intent(inout), dimension(:, :) :: A
    real(DP), intent(inout), dimension(:, :), allocatable :: Qt
    ! -- locals
    real(DP), dimension(2, 2) :: G
    integer(I4B) :: m, n, I

    m = size(A, dim=1) ! Number of rows
    n = size(A, dim=2) ! Number of columns

    do I = m, 1, -1
      G = givens_rotation(A(I, I), A(I, M + 1))
      A(:, [I, M + 1]) = matmul(A(:, [I, M + 1]), transpose(G))
      Qt([I, M + 1], :) = matmul(G, Qt([I, M + 1], :))
    end do

  end subroutine make_matrix_square

  !!> @brief Sets small superdiagonal elements to zero for numerical stability.
  !!
  !! This subroutine scans the superdiagonal elements of the matrix `A` and sets to zero any element
  !! whose absolute value is less than or equal to a tolerance (relative to the neighboring diagonal elements).
  !! This helps to maintain numerical stability and ensures that the matrix is properly treated as bidiagonal
  !! in subsequent SVD steps.
  !!
  !! The criterion used is:
  !!   |A(j, j+1)| <= DPREC * (|A(j, j)| + |A(j+1, j+1)|)
  !! In such cases, A(j, j+1) is set to zero.
  !<
  pure subroutine clean_superdiagonal(A)
    ! -- dummy
    real(DP), intent(inout), dimension(:, :) :: A
    ! -- locals
    integer(I4B) :: n, m, j

    m = size(A, DIM=1) ! Number of rows
    n = size(A, DIM=2) ! Number of columns

    do j = 1, min(n, m) - 1
      if (abs(A(j, j + 1)) <= DPREC * (abs(A(j, j)) + abs(A(j + 1, j + 1)))) then
        A(j, j + 1) = 0.0_DP
      end if
    end do

  end subroutine clean_superdiagonal

  !> @brief Singular Value Decomposition
  !!
  !! This method decomposes the matrix A into U, S and VT.
  !! It follows the algorithm as described by Golub and Reinsch.
  !!
  !! The first step is to decompose the matrix A into a bidiagonal matrix.
  !! This is done using Householder transformations.
  !! Then second step is to decompose the bidiagonal matrix into U, S and VT
  !! by repetitively applying the QR algorithm.
  !! If there is a zero on the diagonal or superdiagonal the matrix can be split
  !! into two smaller matrices and the QR algorithm can be applied to the smaller
  !! matrices.
  !!
  !! The matrix U is the eigenvectors of A*A^T
  !! The matrix VT is the eigenvectors of A^T*A
  !! The matrix S is the square root of the eigenvalues of A*A^T or A^T*A
  !!
  !<
  pure subroutine SVD(A, U, S, VT)
    ! -- dummy
    real(DP), intent(in), dimension(:, :) :: A
    real(DP), intent(out), dimension(:, :), allocatable :: U
    real(DP), intent(out), dimension(:, :), allocatable :: S
    real(DP), intent(out), dimension(:, :), allocatable :: VT
    ! -- locals
    integer(I4B) :: i ! Loop counter for QR iterations
    integer(I4B) :: max_itr ! Maximum number of QR iterations allowed for convergence
    real(DP) :: error ! Convergence criterion (superdiagonal norm)

    integer(I4B) :: m ! Number of rows in input matrix A
    integer(I4B) :: n ! Number of columns in input matrix A

    integer(I4B) :: r ! Start index of active submatrix (for QR step)
    integer(I4B) :: q ! End index of active submatrix (for QR step)

    max_itr = 100

    m = size(A, dim=1) ! Number of rows
    n = size(A, dim=2) ! Number of columns
    S = A

    call bidiagonal_decomposition(S, U, Vt)

    if (n > m) then
      call make_matrix_square(S, Vt)
    end if

    do i = 1, max_itr
      call find_nonzero_superdiagonal(S, r, q)

      ! find zero diagonal
      if (has_zero_diagonal(S(r:q, r:q))) then
        ! write(*,*) 'Iteration: ', i, ' handle zero diagonal element'
        call handle_zero_diagonal(S(r:q, r:q), U(:, r:q), Vt(r:q, :))
      else
        call bidiagonal_qr_decomposition(S(r:q, r:q), U(:, r:q), Vt(r:q, :))
        call clean_superdiagonal(S(r:q, r:q))

        error = superdiagonal_norm(S)
        ! write(*,*) 'Iteration: ', i, ' Error: ', error
        if (error < DPREC) exit
      end if
    end do

  end subroutine SVD

end module
