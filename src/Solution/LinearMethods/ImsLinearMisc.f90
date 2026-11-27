MODULE IMSLinearMisc

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DEM30

  private
  public :: ims_misc_dvscale

CONTAINS

  !
  !> @brief Scale X and RHS
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
