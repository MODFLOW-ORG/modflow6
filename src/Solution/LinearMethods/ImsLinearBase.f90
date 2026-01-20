
!> @brief This module contains the IMS linear accelerator subroutines
!!
!! This module contains the IMS linear accelerator subroutines used by a
!! MODFLOW 6 solution.
!<
MODULE IMSLinearBaseModule
  ! -- modules
  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH, IZERO, &
                             DZERO, DPREC, DEM6, DEM3, DEM4, DHALF, DONE
  use MathUtilModule, only: is_close
  use BlockParserModule, only: BlockParserType
  use IMSReorderingModule, only: ims_odrv
  use ConvergenceSummaryModule

  IMPLICIT NONE

  type(BlockParserType), private :: parser

contains

  !> @ brief Preconditioned Conjugate Gradient linear accelerator
  !!
  !!  Apply the Preconditioned Conjugate Gradient linear accelerator to
  !!  the current coefficient matrix, right-hand side, using the current
  !!  dependent-variable.
  !!
  !<
  SUBROUTINE ims_base_cg(ICNVG, ITMAX, INNERIT, &
                         NEQ, NJA, NIAPC, NJAPC, &
                         IPC, ICNVGOPT, NORTH, &
                         DVCLOSE, RCLOSE, L2NORM0, EPFACT, &
                         IA0, JA0, A0, IAPC, JAPC, APC, &
                         X, B, D, P, Q, Z, &
                         NJLU, IW, JLU, &
                         NCONV, CONVNMOD, CONVMODSTART, &
                         CACCEL, summary)
    ! -- dummy variables
    integer(I4B), INTENT(INOUT) :: ICNVG !< convergence flag (1) non-convergence (0)
    integer(I4B), INTENT(IN) :: ITMAX !< maximum number of inner iterations
    integer(I4B), INTENT(INOUT) :: INNERIT !< inner iteration count
    integer(I4B), INTENT(IN) :: NEQ !< number of equations
    integer(I4B), INTENT(IN) :: NJA !< number of non-zero entries
    integer(I4B), INTENT(IN) :: NIAPC !< preconditioner number of rows
    integer(I4B), INTENT(IN) :: NJAPC !< preconditioner number of non-zero entries
    integer(I4B), INTENT(IN) :: IPC !< preconditioner option
    integer(I4B), INTENT(IN) :: ICNVGOPT !< flow convergence criteria option
    integer(I4B), INTENT(IN) :: NORTH !< orthogonalization frequency
    real(DP), INTENT(IN) :: DVCLOSE !< dependent-variable closure criteria
    real(DP), INTENT(IN) :: RCLOSE !< flow closure criteria
    real(DP), INTENT(IN) :: L2NORM0 !< initial L-2 norm for system of equations
    real(DP), INTENT(IN) :: EPFACT !< factor for decreasing flow convergence criteria for subsequent Picard iterations
    integer(I4B), DIMENSION(NEQ + 1), INTENT(IN) :: IA0 !< CRS row pointers
    integer(I4B), DIMENSION(NJA), INTENT(IN) :: JA0 !< CRS column pointers
    real(DP), DIMENSION(NJA), INTENT(IN) :: A0 !< coefficient matrix
    integer(I4B), DIMENSION(NIAPC + 1), INTENT(IN) :: IAPC !< preconditioner CRS row pointers
    integer(I4B), DIMENSION(NJAPC), INTENT(IN) :: JAPC !< preconditioner CRS column pointers
    real(DP), DIMENSION(NJAPC), INTENT(IN) :: APC !< preconditioner matrix
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: X !< dependent-variable vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: B !< right-hand side vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: D !< working vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: P !< working vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: Q !< working vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: Z !< working vector
    ! -- ILUT dummy variables
    integer(I4B), INTENT(IN) :: NJLU !< preconditioner length of JLU vector
    integer(I4B), DIMENSION(NIAPC), INTENT(IN) :: IW !< preconditioner integer working vector
    integer(I4B), DIMENSION(NJLU), INTENT(IN) :: JLU !< preconditioner JLU working vector
    ! -- convergence information dummy variables dummy variables
    integer(I4B), INTENT(IN) :: NCONV !< maximum number of inner iterations in a time step (maxiter * maxinner)
    integer(I4B), INTENT(IN) :: CONVNMOD !< number of models in the solution
    integer(I4B), DIMENSION(CONVNMOD + 1), INTENT(INOUT) :: CONVMODSTART !< pointer to the start of each model in the convmod* arrays
    character(len=31), DIMENSION(NCONV), INTENT(INOUT) :: CACCEL !< convergence string
    type(ConvergenceSummaryType), pointer, intent(in) :: summary !< Convergence summary report
    ! -- local variables
    LOGICAL :: lorth
    logical :: lsame
    character(len=31) :: cval
    integer(I4B) :: n
    integer(I4B) :: iiter
    integer(I4B) :: xloc, rloc
    integer(I4B) :: im, im0, im1
    real(DP) :: ddot
    real(DP) :: tv
    real(DP) :: deltax
    real(DP) :: rmax
    real(DP) :: l2norm
    real(DP) :: rcnvg
    real(DP) :: denominator
    real(DP) :: alpha, beta
    real(DP) :: rho, rho0
    !
    ! -- initialize local variables
    rho0 = DZERO
    rho = DZERO
    INNERIT = 0
    !
    ! -- INNER ITERATION
    INNER: DO iiter = 1, itmax
      INNERIT = INNERIT + 1
      summary%iter_cnt = summary%iter_cnt + 1
      !
      ! -- APPLY PRECONDITIONER
      SELECT CASE (IPC)
        !
        ! -- ILU0 AND MILU0
      CASE (1, 2)
        CALL ims_base_ilu0a(NJA, NEQ, APC, IAPC, JAPC, D, Z)
        !
        ! -- ILUT AND MILUT
      CASE (3, 4)
        CALL lusol(NEQ, D, Z, APC, JLU, IW)
      END SELECT
      rho = ddot(NEQ, D, 1, Z, 1)
      !
      ! -- COMPUTE DIRECTIONAL VECTORS
      IF (IITER == 1) THEN
        DO n = 1, NEQ
          P(n) = Z(n)
        END DO
      ELSE
        beta = rho / rho0
        DO n = 1, NEQ
          P(n) = Z(n) + beta * P(n)
        END DO
      END IF
      !
      ! -- COMPUTE ITERATES
      !
      ! -- UPDATE Q
      call amux(NEQ, P, Q, A0, JA0, IA0)
      denominator = ddot(NEQ, P, 1, Q, 1)
      denominator = denominator + SIGN(DPREC, denominator)
      alpha = rho / denominator
      !
      ! -- UPDATE X AND RESIDUAL
      deltax = DZERO
      rmax = DZERO
      l2norm = DZERO
      DO im = 1, CONVNMOD
        summary%locdv(im) = 0
        summary%dvmax(im) = DZERO
        summary%locr(im) = 0
        summary%rmax(im) = DZERO
      END DO
      im = 1
      im0 = CONVMODSTART(1)
      im1 = CONVMODSTART(2)
      DO n = 1, NEQ
        !
        ! -- determine current model index
        if (n == im1) then
          im = im + 1
          im0 = CONVMODSTART(im)
          im1 = CONVMODSTART(im + 1)
        end if
        !
        ! -- identify deltax and rmax
        tv = alpha * P(n)
        X(n) = X(n) + tv
        IF (ABS(tv) > ABS(deltax)) THEN
          deltax = tv
          xloc = n
        END IF
        IF (ABS(tv) > ABS(summary%dvmax(im))) THEN
          summary%dvmax(im) = tv
          summary%locdv(im) = n
        END IF
        tv = D(n)
        tv = tv - alpha * Q(n)
        D(n) = tv
        IF (ABS(tv) > ABS(rmax)) THEN
          rmax = tv
          rloc = n
        END IF
        IF (ABS(tv) > ABS(summary%rmax(im))) THEN
          summary%rmax(im) = tv
          summary%locr(im) = n
        END IF
        l2norm = l2norm + tv * tv
      END DO
      l2norm = SQRT(l2norm)
      !
      ! -- SAVE SOLVER convergence information dummy variables
      IF (NCONV > 1) THEN !<
        n = summary%iter_cnt
        WRITE (cval, '(g15.7)') alpha
        CACCEL(n) = cval
        summary%itinner(n) = iiter
        DO im = 1, CONVNMOD
          summary%convlocdv(im, n) = summary%locdv(im)
          summary%convlocr(im, n) = summary%locr(im)
          summary%convdvmax(im, n) = summary%dvmax(im)
          summary%convrmax(im, n) = summary%rmax(im)
        END DO
      END IF
      !
      ! -- TEST FOR SOLVER CONVERGENCE
      IF (ICNVGOPT == 2 .OR. ICNVGOPT == 3 .OR. ICNVGOPT == 4) THEN
        rcnvg = l2norm
      ELSE
        rcnvg = rmax
      END IF
      CALL ims_base_testcnvg(ICNVGOPT, ICNVG, INNERIT, &
                             deltax, rcnvg, &
                             L2NORM0, EPFACT, DVCLOSE, RCLOSE)
      !
      ! -- CHECK FOR EXACT SOLUTION
      IF (rcnvg == DZERO) ICNVG = 1
      !
      ! -- CHECK FOR STANDARD CONVERGENCE
      IF (ICNVG .NE. 0) EXIT INNER
      !
      ! -- CHECK THAT CURRENT AND PREVIOUS rho ARE DIFFERENT
      lsame = is_close(rho, rho0)
      IF (lsame) THEN
        EXIT INNER
      END IF
      !
      ! -- RECALCULATE THE RESIDUAL
      IF (NORTH > 0) THEN
        lorth = mod(iiter + 1, NORTH) == 0
        IF (lorth) THEN
          call ims_base_residual(NEQ, NJA, X, B, D, A0, IA0, JA0)
        END IF
      END IF
      !
      ! -- exit inner if rho is zero
      if (rho == DZERO) then
        exit inner
      end if
      !
      ! -- SAVE CURRENT INNER ITERATES
      rho0 = rho
    END DO INNER
    !
    ! -- RESET ICNVG
    IF (ICNVG < 0) ICNVG = 0
  end SUBROUTINE ims_base_cg

  !> @ brief Preconditioned BiConjugate Gradient Stabilized linear accelerator
  !!
  !!  Apply the Preconditioned BiConjugate Gradient Stabilized linear
  !!  accelerator to the current coefficient matrix, right-hand side, using
  !!  the currentdependent-variable.
  !!
  !<
  SUBROUTINE ims_base_bcgs(ICNVG, ITMAX, INNERIT, &
                           NEQ, NJA, NIAPC, NJAPC, &
                           IPC, ICNVGOPT, NORTH, ISCL, DSCALE, &
                           DVCLOSE, RCLOSE, L2NORM0, EPFACT, &
                           IA0, JA0, A0, IAPC, JAPC, APC, &
                           X, B, D, P, Q, &
                           T, V, DHAT, PHAT, QHAT, &
                           NJLU, IW, JLU, &
                           NCONV, CONVNMOD, CONVMODSTART, &
                           CACCEL, summary)
    ! -- dummy variables
    integer(I4B), INTENT(INOUT) :: ICNVG !< convergence flag (1) non-convergence (0)
    integer(I4B), INTENT(IN) :: ITMAX !< maximum number of inner iterations
    integer(I4B), INTENT(INOUT) :: INNERIT !< inner iteration count
    integer(I4B), INTENT(IN) :: NEQ !< number of equations
    integer(I4B), INTENT(IN) :: NJA !< number of non-zero entries
    integer(I4B), INTENT(IN) :: NIAPC !< preconditioner number of rows
    integer(I4B), INTENT(IN) :: NJAPC !< preconditioner number of non-zero entries
    integer(I4B), INTENT(IN) :: IPC !< preconditioner option
    integer(I4B), INTENT(IN) :: ICNVGOPT !< flow convergence criteria option
    integer(I4B), INTENT(IN) :: NORTH !< orthogonalization frequency
    integer(I4B), INTENT(IN) :: ISCL !< scaling option
    real(DP), DIMENSION(NEQ), INTENT(IN) :: DSCALE !< scaling vector
    real(DP), INTENT(IN) :: DVCLOSE !< dependent-variable closure criteria
    real(DP), INTENT(IN) :: RCLOSE !< flow closure criteria
    real(DP), INTENT(IN) :: L2NORM0 !< initial L-2 norm for system of equations
    real(DP), INTENT(IN) :: EPFACT !< factor for decreasing flow convergence criteria for subsequent Picard iterations
    integer(I4B), DIMENSION(NEQ + 1), INTENT(IN) :: IA0 !< CRS row pointers
    integer(I4B), DIMENSION(NJA), INTENT(IN) :: JA0 !< CRS column pointers
    real(DP), DIMENSION(NJA), INTENT(IN) :: A0 !< coefficient matrix
    integer(I4B), DIMENSION(NIAPC + 1), INTENT(IN) :: IAPC !< preconditioner CRS row pointers
    integer(I4B), DIMENSION(NJAPC), INTENT(IN) :: JAPC !< preconditioner CRS column pointers
    real(DP), DIMENSION(NJAPC), INTENT(IN) :: APC !< preconditioner matrix
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: X !< dependent-variable vector
    real(DP), DIMENSION(NEQ), INTENT(IN) :: B !< right-hand side vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: D !< preconditioner working vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: P !< preconditioner working vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: Q !< preconditioner working vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: T !< preconditioner working vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: V !< preconditioner working vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: DHAT !< BCGS preconditioner working vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: PHAT !< BCGS preconditioner working vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: QHAT !< BCGS preconditioner working vector
    ! -- ILUT dummy variables
    integer(I4B), INTENT(IN) :: NJLU !< preconditioner length of JLU vector
    integer(I4B), DIMENSION(NIAPC), INTENT(IN) :: IW !< preconditioner integer working vector
    integer(I4B), DIMENSION(NJLU), INTENT(IN) :: JLU !< preconditioner JLU working vector
    ! -- convergence information dummy variables
    integer(I4B), INTENT(IN) :: NCONV !< maximum number of inner iterations in a time step (maxiter * maxinner)
    integer(I4B), INTENT(IN) :: CONVNMOD !< number of models in the solution
    integer(I4B), DIMENSION(CONVNMOD + 1), INTENT(INOUT) :: CONVMODSTART !< pointer to the start of each model in the convmod* arrays
    character(len=31), DIMENSION(NCONV), INTENT(INOUT) :: CACCEL !< convergence string
    type(ConvergenceSummaryType), pointer, intent(in) :: summary !< Convergence summary report
    ! -- local variables
    LOGICAL :: LORTH
    logical :: lsame
    character(len=15) :: cval1, cval2
    integer(I4B) :: n
    integer(I4B) :: iiter
    integer(I4B) :: xloc, rloc
    integer(I4B) :: im, im0, im1
    real(DP) :: ddot
    real(DP) :: tv
    real(DP) :: deltax
    real(DP) :: rmax
    real(DP) :: l2norm
    real(DP) :: rcnvg
    real(DP) :: alpha, alpha0
    real(DP) :: beta
    real(DP) :: rho, rho0
    real(DP) :: omega, omega0
    real(DP) :: numerator, denominator
    !
    ! -- initialize local variables
    INNERIT = 0
    alpha = DZERO
    alpha0 = DZERO
    beta = DZERO
    rho = DZERO
    rho0 = DZERO
    omega = DZERO
    omega0 = DZERO
    !
    ! -- SAVE INITIAL RESIDUAL
    DO n = 1, NEQ
      DHAT(n) = D(n)
    END DO
    !
    ! -- INNER ITERATION
    INNER: DO iiter = 1, itmax
      INNERIT = INNERIT + 1
      summary%iter_cnt = summary%iter_cnt + 1
      !
      ! -- CALCULATE rho
      rho = ddot(NEQ, DHAT, 1, D, 1)
      !
      ! -- COMPUTE DIRECTIONAL VECTORS
      IF (IITER == 1) THEN
        DO n = 1, NEQ
          P(n) = D(n)
        END DO
      ELSE
        beta = (rho / rho0) * (alpha0 / omega0)
        DO n = 1, NEQ
          P(n) = D(n) + beta * (P(n) - omega0 * V(n))
        END DO
      END IF
      !
      ! -- APPLY PRECONDITIONER TO UPDATE PHAT
      SELECT CASE (IPC)
        !
        ! -- ILU0 AND MILU0
      CASE (1, 2)
        CALL ims_base_ilu0a(NJA, NEQ, APC, IAPC, JAPC, P, PHAT)
        !
        ! -- ILUT AND MILUT
      CASE (3, 4)
        CALL lusol(NEQ, P, PHAT, APC, JLU, IW)
      END SELECT
      !
      ! -- COMPUTE ITERATES
      !
      ! -- UPDATE V WITH A AND PHAT
      call amux(NEQ, PHAT, V, A0, JA0, IA0)
      !
      ! -- UPDATE alpha WITH DHAT AND V
      denominator = ddot(NEQ, DHAT, 1, V, 1)
      denominator = denominator + SIGN(DPREC, denominator)
      alpha = rho / denominator
      !
      ! -- UPDATE Q
      DO n = 1, NEQ
        Q(n) = D(n) - alpha * V(n)
      END DO
      !
      ! ! -- CALCULATE INFINITY NORM OF Q - TEST FOR TERMINATION
      ! !    TERMINATE IF rmax IS LESS THAN MACHINE PRECISION (DPREC)
      ! rmax = DZERO
      ! DO n = 1, NEQ
      !     tv = Q(n)
      !     IF (ISCL.NE.0 ) tv = tv / DSCALE(n)
      !     IF (ABS(tv) > ABS(rmax) ) rmax = tv
      ! END DO
      ! IF (ABS(rmax).LE.DPREC) THEN
      !   deltax = DZERO
      !   DO n = 1, NEQ
      !     tv      = alpha * PHAT(n)
      !     IF (ISCL.NE.0) THEN
      !       tv = tv * DSCALE(n)
      !     END IF
      !     X(n)  = X(n) + tv
      !     IF (ABS(tv) > ABS(deltax) ) deltax = tv
      !   END DO
      !   CALL IMSLINEARSUB_TESTCNVG(ICNVGOPT, ICNVG, INNERIT, &
      !                             deltax, rmax, &
      !                             rmax, EPFACT, DVCLOSE, RCLOSE )
      !   IF (ICNVG.NE.0 ) EXIT INNER
      ! END IF
      !
      ! -- APPLY PRECONDITIONER TO UPDATE QHAT
      SELECT CASE (IPC)
        !
        ! -- ILU0 AND MILU0
      CASE (1, 2)
        CALL ims_base_ilu0a(NJA, NEQ, APC, IAPC, JAPC, Q, QHAT)
        !
        ! -- ILUT AND MILUT
      CASE (3, 4)
        CALL lusol(NEQ, Q, QHAT, APC, JLU, IW)
      END SELECT
      !
      ! -- UPDATE T WITH A AND QHAT
      call amux(NEQ, QHAT, T, A0, JA0, IA0)
      !
      ! -- UPDATE omega
      numerator = ddot(NEQ, T, 1, Q, 1)
      denominator = ddot(NEQ, T, 1, T, 1)
      denominator = denominator + SIGN(DPREC, denominator)
      omega = numerator / denominator
      !
      ! -- UPDATE X AND RESIDUAL
      deltax = DZERO
      rmax = DZERO
      l2norm = DZERO
      DO im = 1, CONVNMOD
        summary%dvmax(im) = DZERO
        summary%rmax(im) = DZERO
      END DO
      im = 1
      im0 = CONVMODSTART(1)
      im1 = CONVMODSTART(2)
      DO n = 1, NEQ
        !
        ! -- determine current model index
        if (n == im1) then
          im = im + 1
          im0 = CONVMODSTART(im)
          im1 = CONVMODSTART(im + 1)
        end if
        !
        ! -- X AND DX
        tv = alpha * PHAT(n) + omega * QHAT(n)
        X(n) = X(n) + tv
        IF (ISCL .NE. 0) THEN
          tv = tv * DSCALE(n)
        END IF
        IF (ABS(tv) > ABS(deltax)) THEN
          deltax = tv
          xloc = n
        END IF
        IF (ABS(tv) > ABS(summary%dvmax(im))) THEN
          summary%dvmax(im) = tv
          summary%locdv(im) = n
        END IF
        !
        ! -- RESIDUAL
        tv = Q(n) - omega * T(n)
        D(n) = tv
        IF (ISCL .NE. 0) THEN
          tv = tv / DSCALE(n)
        END IF
        IF (ABS(tv) > ABS(rmax)) THEN
          rmax = tv
          rloc = n
        END IF
        IF (ABS(tv) > ABS(summary%rmax(im))) THEN
          summary%rmax(im) = tv
          summary%locr(im) = n
        END IF
        l2norm = l2norm + tv * tv
      END DO
      l2norm = sqrt(l2norm)
      !
      ! -- SAVE SOLVER convergence information dummy variables
      IF (NCONV > 1) THEN !<
        n = summary%iter_cnt
        WRITE (cval1, '(g15.7)') alpha
        WRITE (cval2, '(g15.7)') omega
        CACCEL(n) = trim(adjustl(cval1))//','//trim(adjustl(cval2))
        summary%itinner(n) = iiter
        DO im = 1, CONVNMOD
          summary%convdvmax(im, n) = summary%dvmax(im)
          summary%convlocdv(im, n) = summary%locdv(im)
          summary%convrmax(im, n) = summary%rmax(im)
          summary%convlocr(im, n) = summary%locr(im)
        END DO
      END IF
      !
      ! -- TEST FOR SOLVER CONVERGENCE
      IF (ICNVGOPT == 2 .OR. ICNVGOPT == 3 .OR. ICNVGOPT == 4) THEN
        rcnvg = l2norm
      ELSE
        rcnvg = rmax
      END IF
      CALL ims_base_testcnvg(ICNVGOPT, ICNVG, INNERIT, &
                             deltax, rcnvg, &
                             L2NORM0, EPFACT, DVCLOSE, RCLOSE)
      !
      ! -- CHECK FOR EXACT SOLUTION
      IF (rcnvg == DZERO) ICNVG = 1
      !
      ! -- CHECK FOR STANDARD CONVERGENCE
      IF (ICNVG .NE. 0) EXIT INNER
      !
      ! -- CHECK THAT CURRENT AND PREVIOUS rho, alpha, AND omega ARE
      !    DIFFERENT
      lsame = is_close(rho, rho0)
      IF (lsame) THEN
        EXIT INNER
      END IF
      lsame = is_close(alpha, alpha0)
      IF (lsame) THEN
        EXIT INNER
      END IF
      lsame = is_close(omega, omega0)
      IF (lsame) THEN
        EXIT INNER
      END IF
      !
      ! -- RECALCULATE THE RESIDUAL
      IF (NORTH > 0) THEN
        LORTH = mod(iiter + 1, NORTH) == 0
        IF (LORTH) THEN
          call ims_base_residual(NEQ, NJA, X, B, D, A0, IA0, JA0)
        END IF
      END IF
      !
      ! -- exit inner if rho or omega are zero
      if (rho * omega == DZERO) then
        exit inner
      end if
      !
      ! -- SAVE CURRENT INNER ITERATES
      rho0 = rho
      alpha0 = alpha
      omega0 = omega
    END DO INNER
    !
    ! -- RESET ICNVG
    IF (ICNVG < 0) ICNVG = 0
  end SUBROUTINE ims_base_bcgs

  !> @ brief Calculate LORDER AND IORDER
  !!
  !!  Calculate LORDER and IORDER for reordering.
  !!
  !<
  SUBROUTINE ims_base_calc_order(IORD, NEQ, NJA, IA, JA, LORDER, IORDER)
    ! -- modules
    use SimModule, only: store_error, count_errors
    ! -- dummy variables
    integer(I4B), INTENT(IN) :: IORD !< reordering option
    integer(I4B), INTENT(IN) :: NEQ !< number of rows
    integer(I4B), INTENT(IN) :: NJA !< number of non-zero entries
    integer(I4B), DIMENSION(NEQ + 1), INTENT(IN) :: IA !< row pointer
    integer(I4B), DIMENSION(NJA), INTENT(IN) :: JA !< column pointer
    integer(I4B), DIMENSION(NEQ), INTENT(INOUT) :: LORDER !< reorder vector
    integer(I4B), DIMENSION(NEQ), INTENT(INOUT) :: IORDER !< inverse of reorder vector
    ! -- local variables
    character(len=LINELENGTH) :: errmsg
    integer(I4B) :: n
    integer(I4B) :: nsp
    integer(I4B), DIMENSION(:), ALLOCATABLE :: iwork0
    integer(I4B), DIMENSION(:), ALLOCATABLE :: iwork1
    integer(I4B) :: iflag
    !
    ! -- initialize lorder and iorder
    DO n = 1, NEQ
      LORDER(n) = IZERO
      IORDER(n) = IZERO
    END DO
    ! ALLOCATE (iwork0(NEQ))
    SELECT CASE (IORD)
    CASE (1)
      CALL genrcm(NEQ, NJA, IA, JA, LORDER)
    CASE (2)
      nsp = 3 * NEQ + 4 * NJA
      allocate (iwork0(NEQ))
      allocate (iwork1(nsp))
      CALL ims_odrv(NEQ, NJA, nsp, IA, JA, LORDER, iwork0, &
                    iwork1, iflag)
      IF (iflag .NE. 0) THEN
        write (errmsg, '(A,1X,A)') &
          'IMSLINEARSUB_CALC_ORDER error creating minimum degree ', &
          'order permutation '
        call store_error(errmsg)
      END IF
      !
      ! -- DEALLOCATE TEMPORARY STORAGE
      deallocate (iwork0, iwork1)
    END SELECT
    !
    ! -- GENERATE INVERSE OF LORDER
    DO n = 1, NEQ
      IORDER(LORDER(n)) = n
    END DO
    !
    ! -- terminate if errors occurred
    if (count_errors() > 0) then
      call parser%StoreErrorUnit()
    end if
  end SUBROUTINE ims_base_calc_order

  !
  !> @ brief Scale the coefficient matrix
  !!
  !!  Scale the coefficient matrix (AMAT), the right-hand side (B),
  !!  and the estimate of the dependent variable (X).
  !!
  !<
  SUBROUTINE ims_base_scale(IOPT, ISCL, NEQ, NJA, IA, JA, AMAT, X, B, &
                            DSCALE, DSCALE2)
    ! -- dummy variables
    integer(I4B), INTENT(IN) :: IOPT !< flag to scale (0) or unscale the system of equations
    integer(I4B), INTENT(IN) :: ISCL !< scaling option (1) symmetric (2) L-2 norm
    integer(I4B), INTENT(IN) :: NEQ !< number of equations
    integer(I4B), INTENT(IN) :: NJA !< number of non-zero entries
    integer(I4B), DIMENSION(NEQ + 1), INTENT(IN) :: IA !< CRS row pointer
    integer(I4B), DIMENSION(NJA), INTENT(IN) :: JA !< CRS column pointer
    real(DP), DIMENSION(NJA), INTENT(INOUT) :: AMAT !< coefficient matrix
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: X !< dependent variable
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: B !< right-hand side
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: DSCALE !< first scaling vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: DSCALE2 !< second scaling vector
    ! -- local variables
    integer(I4B) :: i, n
    integer(I4B) :: id, jc
    integer(I4B) :: i0, i1
    real(DP) :: v, c1, c2
    !
    ! -- SCALE SCALE AMAT, X, AND B
    IF (IOPT == 0) THEN
      !
      ! -- SYMMETRIC SCALING
      SELECT CASE (ISCL)
      CASE (1)
        DO n = 1, NEQ
          id = IA(n)
          v = AMAT(id)
          c1 = DONE / SQRT(ABS(v))
          DSCALE(n) = c1
          DSCALE2(n) = c1
        END DO
        !
        ! -- SCALE AMAT -- AMAT = DSCALE(row) * AMAT(i) * DSCALE2(col)
        DO n = 1, NEQ
          c1 = DSCALE(n)
          i0 = IA(n)
          i1 = IA(n + 1) - 1
          DO i = i0, i1
            jc = JA(i)
            c2 = DSCALE2(jc)
            AMAT(i) = c1 * AMAT(i) * c2
          END DO
        END DO
        !
        ! -- L-2 NORM SCALING
      CASE (2)
        !
        ! -- SCALE EACH ROW SO THAT THE L-2 NORM IS 1
        DO n = 1, NEQ
          c1 = DZERO
          i0 = IA(n)
          i1 = IA(n + 1) - 1
          DO i = i0, i1
            c1 = c1 + AMAT(i) * AMAT(i)
          END DO
          c1 = SQRT(c1)
          IF (c1 == DZERO) THEN
            c1 = DONE
          ELSE
            c1 = DONE / c1
          END IF
          DSCALE(n) = c1
          !
          ! -- INITIAL SCALING OF AMAT -- AMAT = DSCALE(row) * AMAT(i)
          DO i = i0, i1
            AMAT(i) = c1 * AMAT(i)
          END DO
        END DO
        !
        ! -- SCALE EACH COLUMN SO THAT THE L-2 NORM IS 1
        DO n = 1, NEQ
          DSCALE2(n) = DZERO
        END DO
        c2 = DZERO
        DO n = 1, NEQ
          i0 = IA(n)
          i1 = IA(n + 1) - 1
          DO i = i0, i1
            jc = JA(i)
            c2 = AMAT(i)
            DSCALE2(jc) = DSCALE2(jc) + c2 * c2
          END DO
        END DO
        DO n = 1, NEQ
          c2 = DSCALE2(n)
          IF (c2 == DZERO) THEN
            c2 = DONE
          ELSE
            c2 = DONE / SQRT(c2)
          END IF
          DSCALE2(n) = c2
        END DO
        !
        ! -- FINAL SCALING OF AMAT -- AMAT = DSCALE2(col) * AMAT(i)
        DO n = 1, NEQ
          i0 = IA(n)
          i1 = IA(n + 1) - 1
          DO i = i0, i1
            jc = JA(i)
            c2 = DSCALE2(jc)
            AMAT(i) = c2 * AMAT(i)
          END DO
        END DO
      END SELECT
      !
      ! -- SCALE X AND B
      DO n = 1, NEQ
        c1 = DSCALE(n)
        c2 = DSCALE2(n)
        X(n) = X(n) / c2
        B(n) = B(n) * c1
      END DO
      !
      ! -- UNSCALE SCALE AMAT, X, AND B
    ELSE
      DO n = 1, NEQ
        c1 = DSCALE(n)
        i0 = IA(n)
        i1 = IA(n + 1) - 1
        !
        ! -- UNSCALE AMAT
        DO i = i0, i1
          jc = JA(i)
          c2 = DSCALE2(jc)
          AMAT(i) = (DONE / c1) * AMAT(i) * (DONE / c2)
        END DO
        !
        ! -- UNSCALE X AND B
        c2 = DSCALE2(n)
        X(n) = X(n) * c2
        B(n) = B(n) / c1
      END DO
    END IF
  end SUBROUTINE ims_base_scale

  !> @ brief Update the preconditioner
  !!
  !!  Update the preconditioner using the current coefficient matrix.
  !!
  !<
  SUBROUTINE ims_base_pcu(IOUT, NJA, NEQ, NIAPC, NJAPC, IPC, RELAX, &
                          AMAT, IA, JA, APC, IAPC, JAPC, IW, W, &
                          LEVEL, DROPTOL, NJLU, NJW, NWLU, JLU, JW, WLU)
    ! -- modules
    use SimModule, only: store_error, count_errors
    ! -- dummy variables
    integer(I4B), INTENT(IN) :: IOUT !< simulation listing file unit
    integer(I4B), INTENT(IN) :: NJA !< number of non-zero entries
    integer(I4B), INTENT(IN) :: NEQ !< number of equations
    integer(I4B), INTENT(IN) :: NIAPC !< preconditioner number of rows
    integer(I4B), INTENT(IN) :: NJAPC !< preconditioner number of non-zero entries
    integer(I4B), INTENT(IN) :: IPC !< precoditioner (1) ILU0 (2) MILU0 (3) ILUT (4) MILUT
    real(DP), INTENT(IN) :: RELAX !< preconditioner relaxation factor for MILU0 and MILUT
    real(DP), DIMENSION(NJA), INTENT(IN) :: AMAT !< coefficient matrix
    integer(I4B), DIMENSION(NEQ + 1), INTENT(IN) :: IA !< CRS row pointers
    integer(I4B), DIMENSION(NJA), INTENT(IN) :: JA !< CRS column pointers
    real(DP), DIMENSION(NJAPC), INTENT(INOUT) :: APC !< preconditioner matrix
    integer(I4B), DIMENSION(NIAPC + 1), INTENT(INOUT) :: IAPC !< preconditioner CRS row pointers
    integer(I4B), DIMENSION(NJAPC), INTENT(INOUT) :: JAPC !< preconditioner CRS column pointers
    integer(I4B), DIMENSION(NIAPC), INTENT(INOUT) :: IW !< preconditioner integed work vector
    real(DP), DIMENSION(NIAPC), INTENT(INOUT) :: W !< preconditioner work vector
    ! -- ILUT dummy variables
    integer(I4B), INTENT(IN) :: LEVEL !< number of levels of fill for ILUT and MILUT
    real(DP), INTENT(IN) :: DROPTOL !< drop tolerance
    integer(I4B), INTENT(IN) :: NJLU !< length of JLU working vector
    integer(I4B), INTENT(IN) :: NJW !< length of JW working vector
    integer(I4B), INTENT(IN) :: NWLU !< length of WLU working vector
    integer(I4B), DIMENSION(NJLU), INTENT(INOUT) :: JLU !< ILUT/MILUT JLU working vector
    integer(I4B), DIMENSION(NJW), INTENT(INOUT) :: JW !< ILUT/MILUT JW working vector
    real(DP), DIMENSION(NWLU), INTENT(INOUT) :: WLU !< ILUT/MILUT WLU working vector
    ! -- local variables
    character(len=LINELENGTH) :: errmsg
    character(len=100), dimension(5), parameter :: cerr = &
      ["Elimination process has generated a row in L or U whose length is > n.", &
      &"The matrix L overflows the array al.                                  ", &
      &"The matrix U overflows the array alu.                                 ", &
      &"Illegal value for lfil.                                               ", &
      &"Zero row encountered.                                                 "]
    integer(I4B) :: ipcflag
    integer(I4B) :: icount
    integer(I4B) :: ierr
    real(DP) :: delta
    ! -- formats
2000 FORMAT(/, ' MATRIX IS SEVERELY NON-DIAGONALLY DOMINANT.', &
            /, ' ADDED SMALL VALUE TO PIVOT ', i0, ' TIMES IN', &
            ' IMSLINEARSUB_PCU.')
    !
    ! -- initialize local variables
    ipcflag = 0
    icount = 0
    delta = DZERO
    PCSCALE: DO
      SELECT CASE (IPC)
        !
        ! -- ILU0 AND MILU0
      CASE (1, 2)
        CALL ims_base_pcilu0(NJA, NEQ, AMAT, IA, JA, &
                             APC, IAPC, JAPC, IW, W, &
                             RELAX, ipcflag, delta)
        !
        ! -- ILUT AND MILUT
      CASE (3, 4)
        ierr = 0
        CALL ilut(NEQ, AMAT, JA, IA, LEVEL, DROPTOL, &
                  APC, JLU, IW, NJAPC, WLU, JW, ierr, &
                  relax, ipcflag, delta)
        if (ierr /= 0) then
          if (ierr > 0) then
            write (errmsg, '(a,1x,i0,1x,a)') &
              'ILUT: zero pivot encountered at step number', ierr, '.'
          else
            write (errmsg, '(a,1x,a)') 'ILUT:', cerr(-ierr)
          end if
          call store_error(errmsg)
          call parser%StoreErrorUnit()
        end if
        !
        ! -- ADDITIONAL PRECONDITIONERS
      CASE DEFAULT
        ipcflag = 0
      END SELECT
      IF (ipcflag < 1) THEN
        EXIT PCSCALE
      END IF
      delta = 1.5d0 * delta + DEM3
      ipcflag = 0
      IF (delta > DHALF) THEN
        delta = DHALF
        ipcflag = 2
      END IF
      icount = icount + 1
      !
      ! -- terminate pcscale loop if not making progress
      if (icount > 10) then
        exit PCSCALE
      end if

    END DO PCSCALE
    !
    ! -- write error message if small value added to pivot
    if (icount > 0) then
      write (IOUT, 2000) icount
    end if
  end SUBROUTINE ims_base_pcu

  !> @ brief Jacobi preconditioner
  !!
  !!  Calculate the Jacobi preconditioner (inverse of the diagonal) using
  !!  the current coefficient matrix.
  !!
  !<
  SUBROUTINE ims_base_pcjac(NJA, NEQ, AMAT, APC, IA, JA)
    ! -- dummy variables
    integer(I4B), INTENT(IN) :: NJA !< number of non-zero entries
    integer(I4B), INTENT(IN) :: NEQ !< number of equations
    real(DP), DIMENSION(NJA), INTENT(IN) :: AMAT !< coefficient matrix
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: APC !< preconditioner matrix
    integer(I4B), DIMENSION(NEQ + 1), INTENT(IN) :: IA !< CRS row pointers
    integer(I4B), DIMENSION(NJA), INTENT(IN) :: JA !< CRS column pointers
    ! -- local variables
    integer(I4B) :: i, n
    integer(I4B) :: ic0, ic1
    integer(I4B) :: id
    real(DP) :: tv
    ! -- code
    DO n = 1, NEQ
      ic0 = IA(n)
      ic1 = IA(n + 1) - 1
      id = IA(n)
      DO i = ic0, ic1
        IF (JA(i) == n) THEN
          id = i
          EXIT
        END IF
      END DO
      tv = AMAT(id)
      IF (ABS(tv) > DZERO) tv = DONE / tv
      APC(n) = tv
    END DO
  end SUBROUTINE ims_base_pcjac

  !> @ brief Apply the Jacobi preconditioner
  !!
  !!  Apply the Jacobi preconditioner and return the resultant vector.
  !!
  !<
  SUBROUTINE ims_base_jaca(NEQ, A, D1, D2)
    ! -- dummy variables
    integer(I4B), INTENT(IN) :: NEQ !< number of equations
    real(DP), DIMENSION(NEQ), INTENT(IN) :: A !< Jacobi preconditioner
    real(DP), DIMENSION(NEQ), INTENT(IN) :: D1 !< input vector
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: D2 !< resultant vector
    !  -- local variables
    integer(I4B) :: n
    real(DP) :: tv
    ! -- code
    DO n = 1, NEQ
      tv = A(n) * D1(n)
      D2(n) = tv
    END DO
  end SUBROUTINE ims_base_jaca

  !> @ brief Update the ILU0 preconditioner
  !!
  !!  Update the ILU0 preconditioner using the current coefficient matrix.
  !!
  !<
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

  !> @ brief Apply the ILU0 and MILU0 preconditioners
  !!
  !!  Apply the ILU0 and MILU0 preconditioners to the passed vector (R).
  !!
  !<
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

  !> @ brief Test for solver convergence
  !!
  !!  General routine for testing for solver convergence based on the
  !!  user-specified convergence option (Icnvgopt).
  !<
  !
  ! -- TEST FOR SOLVER CONVERGENCE
  SUBROUTINE ims_base_testcnvg(Icnvgopt, Icnvg, Iiter, &
                               Dvmax, Rmax, &
                               Rmax0, Epfact, Dvclose, Rclose)
    ! -- dummy variables
    integer(I4B), INTENT(IN) :: Icnvgopt !< convergence option - see documentation for option
    integer(I4B), INTENT(INOUT) :: Icnvg !< flag indicating if convergence achieved (1) or not (0)
    integer(I4B), INTENT(IN) :: Iiter !< inner iteration number (used for strict convergence option)
    real(DP), INTENT(IN) :: Dvmax !< maximum dependent-variable change
    real(DP), INTENT(IN) :: Rmax !< maximum flow change
    real(DP), INTENT(IN) :: Rmax0 !< initial flow change (initial L2-norm)
    real(DP), INTENT(IN) :: Epfact !< factor for reducing convergence criteria in subsequent Picard iterations
    real(DP), INTENT(IN) :: Dvclose !< Maximum depenendent-variable change allowed
    real(DP), INTENT(IN) :: Rclose !< Maximum flow change allowed
    ! -- code
    IF (Icnvgopt == 0) THEN
      IF (ABS(Dvmax) <= Dvclose .AND. ABS(Rmax) <= Rclose) THEN
        Icnvg = 1
      END IF
    ELSE IF (Icnvgopt == 1) THEN
      IF (ABS(Dvmax) <= Dvclose .AND. ABS(Rmax) <= Rclose) THEN
        IF (iiter == 1) THEN
          Icnvg = 1
        ELSE
          Icnvg = -1
        END IF
      END IF
    ELSE IF (Icnvgopt == 2) THEN
      IF (ABS(Dvmax) <= Dvclose .OR. Rmax <= Rclose) THEN
        Icnvg = 1
      ELSE IF (Rmax <= Rmax0 * Epfact) THEN
        Icnvg = -1
      END IF
    ELSE IF (Icnvgopt == 3) THEN
      IF (ABS(Dvmax) <= Dvclose) THEN
        Icnvg = 1
      ELSE IF (Rmax <= Rmax0 * Rclose) THEN
        Icnvg = -1
      END IF
    ELSE IF (Icnvgopt == 4) THEN
      IF (ABS(Dvmax) <= Dvclose .AND. Rmax <= Rclose) THEN
        Icnvg = 1
      ELSE IF (Rmax <= Rmax0 * Epfact) THEN
        Icnvg = -1
      END IF
    END IF
  end SUBROUTINE ims_base_testcnvg

  subroutine ims_calc_pcdims(neq, nja, ia, level, ipc, &
                             niapc, njapc, njlu, njw, nwlu)
    integer(I4B), intent(in) :: neq !< nr. of rows A
    integer(I4B), intent(in) :: nja !< nr. of nonzeros A
    integer(I4B), dimension(:), intent(in) :: ia !< CSR row pointers A
    integer(I4B), intent(in) :: level !< fill level ILU
    integer(I4B), intent(in) :: ipc !< IMS preconditioner type
    integer(I4B), intent(inout) :: niapc !< work array size
    integer(I4B), intent(inout) :: njapc !< work array size
    integer(I4B), intent(inout) :: njlu !< work array size
    integer(I4B), intent(inout) :: njw !< work array size
    integer(I4B), intent(inout) :: nwlu !< work array size
    ! local
    integer(I4B) :: n, i
    integer(I4B) :: ijlu, ijw, iwlu, iwk

    ijlu = 1
    ijw = 1
    iwlu = 1

    ! ILU0 and MILU0
    niapc = neq
    njapc = nja

    ! ILUT and MILUT
    if (ipc == 3 .or. ipc == 4) then
      niapc = neq
      if (level > 0) then
        iwk = neq * (level * 2 + 1)
      else
        iwk = 0
        do n = 1, neq
          i = ia(n + 1) - ia(n)
          if (i > iwk) then
            iwk = i
          end if
        end do
        iwk = neq * iwk
      end if
      njapc = iwk
      ijlu = iwk
      ijw = 2 * neq
      iwlu = neq + 1
    end if

    njlu = ijlu
    njw = ijw
    nwlu = iwlu

  end subroutine ims_calc_pcdims

  !> @ brief Generate CRS pointers for the preconditioner
  !!
  !!  Generate the CRS row and column pointers for the preconditioner.
  !!  JAPC(1:NEQ) hHas the position of the upper entry for a row,
  !!  JAPC(NEQ+1:NJA) is the column position for entry,
  !!  APC(1:NEQ) is the preconditioned inverse of the diagonal, and
  !!  APC(NEQ+1:NJA) are the preconditioned entries for off diagonals.
  !<
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

  !> @brief In-place sorting for an integer array
  !!
  !! Subroutine sort an integer array in-place.
  !!
  !<
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

  !> @brief Calculate residual
  !!
  !! Subroutine to calculate the residual.
  !!
  !<
  SUBROUTINE ims_base_residual(NEQ, NJA, X, B, D, A, IA, JA)
    ! -- dummy variables
    integer(I4B), INTENT(IN) :: NEQ !< length of vectors
    integer(I4B), INTENT(IN) :: NJA !< length of coefficient matrix
    real(DP), DIMENSION(NEQ), INTENT(IN) :: X !< dependent variable
    real(DP), DIMENSION(NEQ), INTENT(IN) :: B !< right-hand side
    real(DP), DIMENSION(NEQ), INTENT(INOUT) :: D !< residual
    real(DP), DIMENSION(NJA), INTENT(IN) :: A !< coefficient matrix
    integer(I4B), DIMENSION(NEQ + 1), INTENT(IN) :: IA !< CRS row pointers
    integer(I4B), DIMENSION(NJA), INTENT(IN) :: JA !< CRS column pointers
    ! -- local variables
    integer(I4B) :: n
    ! -- code
    !
    ! -- calculate matrix-vector product
    call amux(NEQ, X, D, A, JA, IA)
    !
    ! -- subtract matrix-vector product from right-hand side
    DO n = 1, NEQ
      D(n) = B(n) - D(n)
    END DO
  end SUBROUTINE ims_base_residual

  !> @brief Function returning EPFACT
  !<
  function ims_base_epfact(icnvgopt, kstp) result(epfact)
    integer(I4B) :: icnvgopt !< IMS convergence option
    integer(I4B) :: kstp !< time step number
    real(DP) :: epfact !< factor for decreasing convergence criteria in subsequent Picard iterations

    if (icnvgopt == 2) then
      if (kstp == 1) then
        epfact = 0.01
      else
        epfact = 0.10
      end if
    else if (icnvgopt == 4) then
      epfact = DEM4
    else
      epfact = DONE
    end if

  end function ims_base_epfact

END MODULE IMSLinearBaseModule
