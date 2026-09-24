MODULE IMSLinearMisc

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE, DEM6, DEM30

  private
  public :: ims_misc_thomas
  public :: ims_misc_dvscale
  public :: ims_base_pccrs
  public :: ims_base_pcilu0
  public :: ims_base_ilu0a
  public :: ims_base_isort

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

  SUBROUTINE ims_base_pcilu0(NJA, NEQ, AMAT, IA, JA, &
                             APC, IAPC, JAPC, IW, W, &
                             RELAX, IPCFLAG, DELTA)
    ! -- dummy variables
    integer(I4B), INTENT(IN) :: NJA !< number of non-zero entries
    integer(I4B), INTENT(IN) :: NEQ !< number of equations
    real(DP), DIMENSION(NJA), INTENT(IN) :: AMAT !< coefficient matrix
    integer(I4B), DIMENSION(NEQ + 1), INTENT(IN) :: IA !< CRS row pointers
    integer(I4B), DIMENSION(NJA), INTENT(IN) :: JA !< CRS column pointers
    real(DP), DIMENSION(NJA), INTENT(INOUT) :: APC !< preconditioned matrix
    integer(I4B), DIMENSION(NEQ + 1), INTENT(INOUT) :: IAPC !< preconditioner CRS row pointers
    integer(I4B), DIMENSION(NJA), INTENT(INOUT) :: JAPC !< preconditioner CRS column pointers
    integer(I4B), DIMENSION(NEQ), INTENT(INOUT) :: IW !< preconditioner integer work vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: W !< preconditioner work vector
    real(DP), INTENT(IN) :: RELAX !< MILU0 preconditioner relaxation factor
    integer(I4B), INTENT(INOUT) :: IPCFLAG !< preconditioner error flag
    real(DP), INTENT(IN) :: DELTA !< factor used to correct non-diagonally dominant matrices
    ! -- local variables
    integer(I4B) :: ic0, ic1
    integer(I4B) :: iic0, iic1
    integer(I4B) :: iu, iiu
    integer(I4B) :: j, n
    integer(I4B) :: jj
    integer(I4B) :: jcol, jw
    integer(I4B) :: jjcol
    real(DP) :: drelax
    real(DP) :: sd1
    real(DP) :: tl
    real(DP) :: rs
    real(DP) :: d
    !
    ! -- initialize local variables
    drelax = RELAX
    DO n = 1, NEQ
      IW(n) = 0
      W(n) = DZERO
    END DO
    MAIN: DO n = 1, NEQ
      ic0 = IA(n)
      ic1 = IA(n + 1) - 1
      DO j = ic0, ic1
        jcol = JA(j)
        IW(jcol) = 1
        W(jcol) = W(jcol) + AMAT(j)
      END DO
      ic0 = IAPC(n)
      ic1 = IAPC(n + 1) - 1
      iu = JAPC(n)
      rs = DZERO
      LOWER: DO j = ic0, iu - 1
        jcol = JAPC(j)
        iic0 = IAPC(jcol)
        iic1 = IAPC(jcol + 1) - 1
        iiu = JAPC(jcol)
        tl = W(jcol) * APC(jcol)
        W(jcol) = tl
        DO jj = iiu, iic1
          jjcol = JAPC(jj)
          jw = IW(jjcol)
          IF (jw .NE. 0) THEN
            W(jjcol) = W(jjcol) - tl * APC(jj)
          ELSE
            rs = rs + tl * APC(jj)
          END IF
        END DO
      END DO LOWER
      !
      ! -- DIAGONAL - CALCULATE INVERSE OF DIAGONAL FOR SOLUTION
      d = W(n)
      tl = (DONE + DELTA) * d - (drelax * rs)
      !
      ! -- ENSURE THAT THE SIGN OF THE DIAGONAL HAS NOT CHANGED AND IS
      sd1 = SIGN(d, tl)
      IF (sd1 .NE. d) THEN
        !
        ! -- USE SMALL VALUE IF DIAGONAL SCALING IS NOT EFFECTIVE FOR
        !    PIVOTS THAT CHANGE THE SIGN OF THE DIAGONAL
        IF (IPCFLAG > 1) THEN
          tl = SIGN(DEM6, d)
          !
          ! -- DIAGONAL SCALING CONTINUES TO BE EFFECTIVE
        ELSE
          IPCFLAG = 1
          EXIT MAIN
        END IF
      END IF
      IF (ABS(tl) == DZERO) THEN
        !
        ! -- USE SMALL VALUE IF DIAGONAL SCALING IS NOT EFFECTIVE FOR
        !    ZERO PIVOTS
        IF (IPCFLAG > 1) THEN
          tl = SIGN(DEM6, d)
          !
          ! -- DIAGONAL SCALING CONTINUES TO BE EFFECTIVE FOR ELIMINATING
        ELSE
          IPCFLAG = 1
          EXIT MAIN
        END IF
      END IF
      APC(n) = DONE / tl
      !
      ! -- RESET POINTER FOR IW TO ZERO
      IW(n) = 0
      W(n) = DZERO
      DO j = ic0, ic1
        jcol = JAPC(j)
        APC(j) = W(jcol)
        IW(jcol) = 0
        W(jcol) = DZERO
      END DO
    END DO MAIN
    !
    ! -- RESET IPCFLAG IF SUCCESSFUL COMPLETION OF MAIN
    IPCFLAG = 0
  end SUBROUTINE ims_base_pcilu0

  SUBROUTINE ims_base_ilu0a(NJA, NEQ, APC, IAPC, JAPC, R, D)
    ! -- dummy variables
    integer(I4B), INTENT(IN) :: NJA !< number of non-zero entries
    integer(I4B), INTENT(IN) :: NEQ !< number of equations
    real(DP), DIMENSION(NJA), INTENT(IN) :: APC !< ILU0/MILU0 preconditioner matrix
    integer(I4B), DIMENSION(NEQ + 1), INTENT(IN) :: IAPC !< ILU0/MILU0 preconditioner CRS row pointers
    integer(I4B), DIMENSION(NJA), INTENT(IN) :: JAPC !< ILU0/MILU0 preconditioner CRS column pointers
    real(DP), DIMENSION(NEQ), INTENT(IN) :: R !< input vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: D !< output vector after applying APC to R
    ! -- local variables
    integer(I4B) :: ic0, ic1
    integer(I4B) :: iu
    integer(I4B) :: jcol
    integer(I4B) :: j, n
    real(DP) :: tv
    !
    ! -- FORWARD SOLVE - APC * D = R
    FORWARD: DO n = 1, NEQ
      tv = R(n)
      ic0 = IAPC(n)
      ic1 = IAPC(n + 1) - 1
      iu = JAPC(n) - 1
      LOWER: DO j = ic0, iu
        jcol = JAPC(j)
        tv = tv - APC(j) * D(jcol)
      END DO LOWER
      D(n) = tv
    END DO FORWARD
    !
    ! -- BACKWARD SOLVE - D = D / U
    BACKWARD: DO n = NEQ, 1, -1
      ic0 = IAPC(n)
      ic1 = IAPC(n + 1) - 1
      iu = JAPC(n)
      tv = D(n)
      UPPER: DO j = iu, ic1
        jcol = JAPC(j)
        tv = tv - APC(j) * D(jcol)
      END DO UPPER
      !
      ! -- COMPUTE D FOR DIAGONAL - D = D / U
      D(n) = tv * APC(n)
    END DO BACKWARD
  end SUBROUTINE ims_base_ilu0a

  SUBROUTINE ims_base_pccrs(NEQ, NJA, IA, JA, &
                            IAPC, JAPC)
    ! -- dummy variables
    integer(I4B), INTENT(IN) :: NEQ !<
    integer(I4B), INTENT(IN) :: NJA !<
    integer(I4B), DIMENSION(NEQ + 1), INTENT(IN) :: IA !<
    integer(I4B), DIMENSION(NJA), INTENT(IN) :: JA !<
    integer(I4B), DIMENSION(NEQ + 1), INTENT(INOUT) :: IAPC !<
    integer(I4B), DIMENSION(NJA), INTENT(INOUT) :: JAPC !<
    ! -- local variables
    integer(I4B) :: n, j
    integer(I4B) :: i0, i1
    integer(I4B) :: nlen
    integer(I4B) :: ic, ip
    integer(I4B) :: jcol
    integer(I4B), DIMENSION(:), ALLOCATABLE :: iarr
    ! -- code
    ip = NEQ + 1
    DO n = 1, NEQ
      i0 = IA(n)
      i1 = IA(n + 1) - 1
      nlen = i1 - i0
      ALLOCATE (iarr(nlen))
      ic = 0
      DO j = i0, i1
        jcol = JA(j)
        IF (jcol == n) CYCLE
        ic = ic + 1
        iarr(ic) = jcol
      END DO
      CALL ims_base_isort(nlen, iarr)
      IAPC(n) = ip
      DO j = 1, nlen
        jcol = iarr(j)
        JAPC(ip) = jcol
        ip = ip + 1
      END DO
      DEALLOCATE (iarr)
    END DO
    IAPC(NEQ + 1) = NJA + 1
    !
    ! -- POSITION OF THE FIRST UPPER ENTRY FOR ROW
    DO n = 1, NEQ
      i0 = IAPC(n)
      i1 = IAPC(n + 1) - 1
      JAPC(n) = IAPC(n + 1)
      DO j = i0, i1
        jcol = JAPC(j)
        IF (jcol > n) THEN
          JAPC(n) = j
          EXIT
        END IF
      END DO
    END DO
  end SUBROUTINE ims_base_pccrs

  SUBROUTINE ims_base_isort(NVAL, IARRAY)
    ! -- dummy variables
    integer(I4B), INTENT(IN) :: NVAL !< length of the integer array
    integer(I4B), DIMENSION(NVAL), INTENT(INOUT) :: IARRAY !< integer array to be sorted
    ! -- local variables
    integer(I4B) :: i, j, itemp
    ! -- code
    DO i = 1, NVAL - 1
      DO j = i + 1, NVAL
        if (IARRAY(i) > IARRAY(j)) then
          itemp = IARRAY(j)
          IARRAY(j) = IARRAY(i)
          IARRAY(i) = itemp
        END IF
      END DO
    END DO
  end SUBROUTINE ims_base_isort

END MODULE IMSLinearMisc
