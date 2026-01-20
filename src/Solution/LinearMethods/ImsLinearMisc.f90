MODULE IMSLinearMisc

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE, DEM30

  private
  public :: ims_misc_thomas
  public :: ims_misc_dvscale

CONTAINS

  !> @brief Tridiagonal solve using the Thomas algorithm
    !!
    !! Subroutine to solve tridiagonal linear equations using the
    !! Thomas algorithm.
    !!
  !<
  subroutine ims_misc_thomas(n, tl, td, tu, b, x, w)
    implicit none
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
  end subroutine ims_misc_thomas

  !
  !> @ brief Scale X and RHS
  !!
  !!  Scale X and B to avoid big or small values. Scaling value is the
  !!  maximum ABS(X).
  !<
  SUBROUTINE ims_misc_dvscale(IOPT, NEQ, DSCALE, X, B)
    ! -- dummy variables
    integer(I4B), INTENT(IN) :: IOPT !< flag to scale (0) or unscale the system of equations
    integer(I4B), INTENT(IN) :: NEQ !< number of equations
    real(DP), INTENT(INOUT) :: DSCALE !< scaling value
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: X !< dependent variable
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: B !< right-hand side
    ! -- local variables
    integer(I4B) :: n
    !
    ! -- SCALE SCALE AMAT, X, AND B
    IF (IOPT == 0) THEN
      DSCALE = ABS(DSCALE)
      if (DSCALE < DEM30) then
        DSCALE = DONE
      end if
      DO n = 1, NEQ
        B(n) = B(n) / DSCALE
        X(n) = X(n) / DSCALE
      END DO
      !
      ! -- UNSCALE X, AND B
    ELSE
      DO n = 1, NEQ
        B(n) = B(n) * DSCALE
        X(n) = X(n) * DSCALE
      END DO
    END IF
  end SUBROUTINE ims_misc_dvscale

END MODULE IMSLinearMisc
