module TesterUtils
  use testdrive, only: error_type, test_failed, check
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DSAME

  implicit none
  private

  public :: check_same_matrix
  public :: check_matrix_row_orthogonal
  public :: check_matrix_column_orthogonal
  public :: check_matrix_is_diagonal
  public :: check_matrix_bidiagonal

contains
  !> Check if two matrices A and B are the same within a given tolerance.
  !!
  !! This subroutine compares two matrices element-wise and sets the
  !! error flag if any corresponding elements differ by more than the
  !! specified tolerance.
  subroutine check_same_matrix(error, A, B, tolerance)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(:, :), intent(in) :: A, B
    !- Locals
    integer(I4B) :: i, j
    real(DP), optional :: tolerance
    real(DP) :: tol
    logical :: is_different

    if (.not. present(tolerance)) then
      tol = DSAME
    else
      tol = tolerance
    end if

    if (size(A, 1) /= size(B, 1) .or. size(A, 2) /= size(B, 2)) then
      call test_failed(error, "Matrices have different dimensions")
      return
    end if

    is_different = .false.
    do i = 1, size(A, 1)
      do j = 1, size(A, 2)
        if (abs(A(i, j) - B(i, j)) > tol) then
          call test_failed(error, "Matrices have different elements")
          return
        end if
      end do
    end do

  end subroutine check_same_matrix

  !> @brief Checks if the rows of a matrix are orthogonal within a given tolerance.
  !>
  !> This subroutine verifies whether the rows of the input matrix `A` are orthogonal
  !> to each other within a specified tolerance. If the rows are not orthogonal, the
  !> subroutine sets the `error` flag to true.
  subroutine check_matrix_row_orthogonal(error, A, tolerance)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(:, :), intent(in) :: A
    !- Locals
    integer(I4B) :: i, j
    real(DP), optional :: tolerance
    real(DP) :: tol
    real(DP) :: err

    if (.not. present(tolerance)) then
      tol = DSAME
    else
      tol = tolerance
    end if

    do i = 1, size(A, 1) - 1
      do j = i + 1, size(A, 1)
        err = abs(dot_product(A(i, :), A(j, :)))
        call check(error, err, 0.0_DP, thr=DSAME)
        if (allocated(error)) return
      end do
    end do
  end subroutine check_matrix_row_orthogonal

  !> @brief Checks if the columns of a matrix are orthogonal within a given tolerance.
  !>
  !> This subroutine verifies whether the columns of the input matrix `A` are orthogonal
  !> to each other within a specified tolerance. If the columns are not orthogonal, the
  !> subroutine sets the `error` flag to true.
  subroutine check_matrix_column_orthogonal(error, A, tolerance)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(:, :), intent(in) :: A
    !- Locals
    integer(I4B) :: i, j
    real(DP), optional :: tolerance
    real(DP) :: tol
    real(DP) :: err

    if (.not. present(tolerance)) then
      tol = DSAME
    else
      tol = tolerance
    end if

    do i = 1, size(A, 2) - 1
      do j = i + 1, size(A, 2)
        err = abs(dot_product(A(:, i), A(:, j)))
        call check(error, err, 0.0_DP, thr=DSAME)
        if (allocated(error)) return
      end do
    end do
  end subroutine check_matrix_column_orthogonal

  !>  Checks if the given matrix A is diagonal within a specified tolerance.
  !!
  !! This subroutine verifies whether the input matrix A is diagonal by
  !! comparing its off-diagonal elements to zero within a given tolerance.
  subroutine check_matrix_is_diagonal(error, A, tolerance)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(:, :), intent(in) :: A
    !- Locals
    integer(I4B) :: i, j
    real(DP), optional :: tolerance
    real(DP) :: tol
    real(DP) :: err

    if (.not. present(tolerance)) then
      tol = DSAME
    else
      tol = tolerance
    end if

    do i = 1, size(A, 1)
      do j = 1, size(A, 2)
        if (i /= j) then
          err = abs(A(i, j))
          call check(error, err, 0.0_DP, thr=DSAME)
          if (allocated(error)) return
        end if
      end do
    end do

  end subroutine check_matrix_is_diagonal

  !>  Checks if the given matrix A is bidiagonal within a specified tolerance.
  !!
  !! This subroutine verifies whether the input matrix A is bidiagonal by
  !! comparing its off-diagonal and off-superdiagonal elements to zero within
  !! a given tolerance.
  subroutine check_matrix_bidiagonal(error, A, tolerance)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(:, :), intent(in) :: A
    !- Locals
    integer(I4B) :: i, j
    real(DP), optional :: tolerance
    real(DP) :: tol
    real(DP) :: err

    if (.not. present(tolerance)) then
      tol = DSAME
    else
      tol = tolerance
    end if

    do i = 1, size(A, 1)
      do j = 1, size(A, 2)
        if (i /= j .and. i + 1 /= j) then
          err = abs(A(i, j))
          call check(error, err, 0.0_DP, thr=DSAME)
          if (allocated(error)) return
        end if
      end do
    end do

  end subroutine check_matrix_bidiagonal

end module TesterUtils
