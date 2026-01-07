module PetscConvergenceModule
#include <petsc/finclude/petscksp.h>
  use petscksp
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DPREC, DZERO
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use ListModule
  use ConvergenceSummaryModule
  use ImsLinearSettingsModule
  implicit none
  private

  public :: petsc_cnvg_check, petsc_cnvg_check_internal
  public :: KSPSetConvergenceTest

  ! TODO_MJR: this could be smaller, find a bound
  real(DP), private, parameter :: RNORM_L2_TOL = DPREC

  ! Context for the custom convergence check
  type, public :: PetscCnvgCtxType
    Vec :: x_old !< x vector from the previous iteration
    Vec :: delta_x !< delta in x w.r.t. previous iteration
    Vec :: residual !< the unpreconditoned residual vector (a la IMS)
    integer(I4B) :: icnvg_ims !< IMS convergence number: 1 => converged, -1 => forces next Picard iter
    integer(I4B) :: icnvgopt !< convergence option:
                             !! 0,1,2,3,4,.. for equivalent IMS settings
    real(DP) :: dvclose !< dep. variable closure criterion
    real(DP) :: rclose !< residual closure criterion
    integer(I4B) :: max_its !< maximum number of inner iterations
    real(DP) :: rnorm_L2_init !< the initial L2 norm for (b - Ax)
    type(ConvergenceSummaryType), pointer :: cnvg_summary => null() !< detailed convergence information
    real(DP) :: t_convergence_check !< the time spent convergence checking
  contains
    procedure :: create
    procedure :: destroy
  end type PetscCnvgCtxType

  ! passing our context into PETSc requires an explicit interface
  ! on KSPSetConvergenceTest, it is defined here:
  interface
    subroutine CnvgCheckFunc(ksp, n, rnorm, flag, context, ierr)
      import tKSP, PetscCnvgCtxType
      type(tKSP) :: ksp
      PetscInt :: n
      PetscReal :: rnorm
      KSPConvergedReason :: flag
      class(PetscCnvgCtxType), pointer :: context
      PetscErrorCode :: ierr
    end subroutine

    subroutine CnvgDestroyFunc(context, ierr)
      import PetscCnvgCtxType
      class(PetscCnvgCtxType), pointer :: context
      PetscErrorCode :: ierr
    end subroutine

    subroutine KSPSetConvergenceTest(ksp, check_convergence, context, &
                                     destroy, ierr)
      import tKSP, CnvgCheckFunc, PetscCnvgCtxType, CnvgDestroyFunc
      type(tKSP) :: ksp
      procedure(CnvgCheckFunc) :: check_convergence
      class(PetscCnvgCtxType), pointer :: context
      procedure(CnvgDestroyFunc) :: destroy
      PetscErrorCode :: ierr
    end subroutine
  end interface

contains

  subroutine create(this, mat, settings, summary)
    class(PetscCnvgCtxType) :: this
    Mat, pointer :: mat
    type(ImsLinearSettingsType), pointer :: settings
    type(ConvergenceSummaryType), pointer :: summary
    ! local
    PetscErrorCode :: ierr

    this%icnvg_ims = 0
    this%icnvgopt = settings%icnvgopt
    this%dvclose = settings%dvclose
    this%rclose = settings%rclose
    this%max_its = settings%iter1
    this%cnvg_summary => summary
    call MatCreateVecs(mat, this%x_old, PETSC_NULL_VEC, ierr)
    CHKERRQ(ierr)
    call MatCreateVecs(mat, this%delta_x, PETSC_NULL_VEC, ierr)
    CHKERRQ(ierr)
    call MatCreateVecs(mat, this%residual, PETSC_NULL_VEC, ierr)
    CHKERRQ(ierr)

  end subroutine create

  !> @brief Routine to check the convergence following the configuration
  !< of IMS. (called back from the PETSc solver)
  subroutine petsc_cnvg_check(ksp, n, rnorm_L2, flag, context, ierr)
    KSP :: ksp !< Iterative context
    PetscInt :: n !< Iteration number
    PetscReal :: rnorm_L2 !< 2-norm (preconditioned) residual value
    KSPConvergedReason :: flag !< Converged reason
    class(PetscCnvgCtxType), pointer :: context !< context
    PetscErrorCode :: ierr !< error
    ! local
    PetscReal, parameter :: min_one = -1.0
    PetscReal :: xnorm_inf, rnorm0, rnorm_inf_ims, rnorm_L2_ims
    Vec :: x, res
    type(ConvergenceSummaryType), pointer :: summary

    summary => context%cnvg_summary

    ! NB: KSPBuildResidual needs to have its vector destroyed
    ! to avoid a memory leak, KSPBuildSolution doesn't...
    call KSPBuildSolution(ksp, PETSC_NULL_VEC, x, ierr)
    CHKERRQ(ierr)

    ! for CG the KSPBuildResidual returns the work vector directly,
    ! but BCGS (and possibly others) will do the matrix multiplication
    call KSPBuildResidual(ksp, PETSC_NULL_VEC, PETSC_NULL_VEC, res, ierr)
    CHKERRQ(ierr)

    rnorm0 = huge(rnorm0)
    if (context%icnvgopt == 2 .or. &
        context%icnvgopt == 3 .or. &
        context%icnvgopt == 4) then
      call VecNorm(res, NORM_2, rnorm_L2_ims, ierr)
      rnorm0 = rnorm_L2_ims
      CHKERRQ(ierr)
    end if

    ! n == 0 is before the iteration starts
    if (n == 0) then
      context%rnorm_L2_init = rnorm0
      if (rnorm_L2 < RNORM_L2_TOL) then
        ! exact solution found
        flag = KSP_CONVERGED_HAPPY_BREAKDOWN
      else
        call VecCopy(x, context%x_old, ierr)
        CHKERRQ(ierr)
        flag = KSP_CONVERGED_ITERATING
      end if
      ! early return
      call VecDestroy(res, ierr)
      CHKERRQ(ierr)
      return
    end if

    call VecWAXPY(context%delta_x, min_one, context%x_old, x, ierr)
    CHKERRQ(ierr)

    call VecNorm(context%delta_x, NORM_INFINITY, xnorm_inf, ierr)
    CHKERRQ(ierr)

    rnorm_inf_ims = huge(rnorm_inf_ims)
    if (context%icnvgopt == 0 .or. context%icnvgopt == 1) then
      call VecNorm(res, NORM_INFINITY, rnorm_inf_ims, ierr)
      CHKERRQ(ierr)
    end if

    call VecCopy(x, context%x_old, ierr)
    CHKERRQ(ierr)

    ! fill the summary for reporting
    call fill_cnvg_summary(summary, context%delta_x, res, n)

    if (rnorm_L2 < RNORM_L2_TOL) then
      ! exact solution, set to 'converged'
      flag = KSP_CONVERGED_HAPPY_BREAKDOWN
      context%icnvg_ims = 1
      ! check for strict option
      if (n > 1 .and. context%icnvgopt == 1) then
        context%icnvg_ims = -1
      end if
    else
      ! IMS check on convergence
      flag = apply_check(context, n, xnorm_inf, rnorm_inf_ims, rnorm_L2_ims)
    end if

    if (flag == KSP_CONVERGED_ITERATING) then
      ! not yet converged, max. iters reached? Then stop.
      if (n == context%max_its) then
        flag = KSP_DIVERGED_ITS
      end if
    end if

    call VecDestroy(res, ierr)
    CHKERRQ(ierr)

  end subroutine petsc_cnvg_check

  !> @brief Apply the IMS convergence check
  !<
  function apply_check(ctx, nit, dvmax, rnorm_inf, rnorm_L2) result(flag)
    use TdisModule, only: kstp
    use IMSLinearBaseModule, only: ims_base_testcnvg, ims_base_epfact
    class(PetscCnvgCtxType) :: ctx
    integer(I4B) :: nit !< iteration number
    real(DP) :: dvmax !< infinity norm of dep. var. change
    real(DP) :: rnorm_inf !< infinity norm of residual change
    real(DP) :: rnorm_L2 !< L2-norm of residual change
    KSPConvergedReason :: flag !< the convergence status
    ! local
    real(DP) :: epfact
    real(DP) :: rcnvg

    ! Set to 'not converged'
    flag = KSP_CONVERGED_ITERATING
    ctx%icnvg_ims = 0

    epfact = ims_base_epfact(ctx%icnvgopt, kstp)

    if (ctx%icnvgopt == 2 .or. &
        ctx%icnvgopt == 3 .or. &
        ctx%icnvgopt == 4) then
      rcnvg = rnorm_L2
    else
      rcnvg = rnorm_inf
    end if
    call ims_base_testcnvg(ctx%icnvgopt, ctx%icnvg_ims, nit, &
                           dvmax, rcnvg, ctx%rnorm_L2_init, &
                           epfact, ctx%dvclose, ctx%rclose)

    if (ctx%icnvg_ims /= 0) then
      ! Set to 'converged'
      flag = KSP_CONVERGED_HAPPY_BREAKDOWN
    end if

  end function apply_check

  !> @brief Routine to check the convergence following the configuration
  !< of IMS. (called back from the PETSc solver)
  subroutine petsc_cnvg_check_internal(ksp, n, rnorm_L2, flag, context, ierr)
    KSP :: ksp !< Iterative context
    PetscInt :: n !< Iteration number
    PetscReal :: rnorm_L2 !< 2-norm (preconditioned) residual value
    KSPConvergedReason :: flag !< Converged reason
    class(PetscCnvgCtxType), pointer :: context !< context
    PetscErrorCode :: ierr !< error
    ! local
    PetscReal, parameter :: min_one = -1.0
    PetscReal :: xnorm_inf, rnorm0
    Vec :: x, res
    type(ConvergenceSummaryType), pointer :: summary

    summary => context%cnvg_summary

    ! NB: KSPBuildResidual needs to have its vector destroyed
    ! to avoid a memory leak, KSPBuildSolution doesn't...
    call KSPBuildSolution(ksp, PETSC_NULL_VEC, x, ierr)
    CHKERRQ(ierr)

    rnorm0 = rnorm_L2

    ! n == 0 is before the iteration starts
    if (n == 0) then
      context%rnorm_L2_init = rnorm0
      if (rnorm_L2 < RNORM_L2_TOL) then
        ! exact solution found
        flag = KSP_CONVERGED_HAPPY_BREAKDOWN
        context%icnvg_ims = 1
      else
        call VecCopy(x, context%x_old, ierr)
        CHKERRQ(ierr)
        flag = KSP_CONVERGED_ITERATING
      end if
      ! early return
      return
    end if

    call VecWAXPY(context%delta_x, min_one, context%x_old, x, ierr)
    CHKERRQ(ierr)

    call VecNorm(context%delta_x, NORM_INFINITY, xnorm_inf, ierr)
    CHKERRQ(ierr)

    call VecCopy(x, context%x_old, ierr)
    CHKERRQ(ierr)

    ! fill the summary for reporting
    call fill_cnvg_summary_internal(summary, context%delta_x, n)

    context%icnvg_ims = 0
    if (rnorm_L2 < RNORM_L2_TOL) then
      ! exact solution, set to 'converged'
      flag = KSP_CONVERGED_HAPPY_BREAKDOWN
      context%icnvg_ims = 1
      ! apply 'strict'
      if (n > 1 .and. context%icnvgopt == 1) context%icnvg_ims = -1
    else
      ! L2norm
      if (xnorm_inf <= context%dvclose .and. rnorm_L2 <= context%rclose) then
        flag = KSP_CONVERGED_HAPPY_BREAKDOWN
        context%icnvg_ims = 1
        ! apply 'strict'
        if (n > 1 .and. context%icnvgopt == 1) context%icnvg_ims = -1
      end if
    end if

    if (flag == KSP_CONVERGED_ITERATING) then
      ! not yet converged, max. iters reached? Then stop.
      if (n == context%max_its) then
        flag = KSP_DIVERGED_ITS
      end if
    end if

    call VecDestroy(res, ierr)
    CHKERRQ(ierr)

  end subroutine petsc_cnvg_check_internal

  !> @brief Fill the convergence summary from the context
  !<
  subroutine fill_cnvg_summary(summary, dx, res, n)
    type(ConvergenceSummaryType), pointer :: summary !< the convergence summary
    Vec :: dx !< the vector with changes in x
    Vec :: res !< the residual vector
    PetscInt :: n !< the PETSc iteration number
    ! local
    PetscReal, dimension(:), pointer :: local_dx, local_res
    PetscReal :: dvmax_model, rmax_model
    PetscErrorCode :: ierr
    PetscInt :: idx_dv, idx_r
    PetscInt :: i, j, istart, iend
    PetscInt :: iter_cnt

    ! increment iteration counter
    summary%iter_cnt = summary%iter_cnt + 1
    iter_cnt = summary%iter_cnt

    if (summary%nitermax > 1) then
      summary%itinner(iter_cnt) = n
      do i = 1, summary%convnmod
        summary%convdvmax(i, iter_cnt) = DZERO
        summary%convlocdv(i, iter_cnt) = 0
        summary%convrmax(i, iter_cnt) = DZERO
        summary%convlocr(i, iter_cnt) = 0
      end do
    end if

    ! get dv and dr per local model (readonly!)
    call VecGetArrayReadF90(dx, local_dx, ierr)
    CHKERRQ(ierr)
    call VecGetArrayReadF90(res, local_res, ierr)
    CHKERRQ(ierr)
    do i = 1, summary%convnmod
      ! reset
      dvmax_model = DZERO
      idx_dv = 0
      rmax_model = DZERO
      idx_r = 0
      ! get first and last model index
      istart = summary%model_bounds(i)
      iend = summary%model_bounds(i + 1) - 1
      do j = istart, iend
        if (abs(local_dx(j)) > abs(dvmax_model)) then
          dvmax_model = local_dx(j)
          idx_dv = j
        end if
        if (abs(local_res(j)) > abs(rmax_model)) then
          rmax_model = local_res(j)
          idx_r = j
        end if
      end do
      if (summary%nitermax > 1) then
        summary%convdvmax(i, iter_cnt) = dvmax_model
        summary%convlocdv(i, iter_cnt) = idx_dv
        summary%convrmax(i, iter_cnt) = rmax_model
        summary%convlocr(i, iter_cnt) = idx_r
      end if
    end do
    call VecRestoreArrayF90(dx, local_dx, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(res, local_res, ierr)
    CHKERRQ(ierr)

  end subroutine fill_cnvg_summary

  !> @brief Fill the convergence summary from the context
  !<
  subroutine fill_cnvg_summary_internal(summary, dx, n)
    type(ConvergenceSummaryType), pointer :: summary !< the convergence summary
    Vec :: dx !< the vector with changes in x
    PetscInt :: n !< the PETSc iteration number
    ! local
    PetscReal, dimension(:), pointer :: local_dx
    PetscReal :: dvmax_model
    PetscErrorCode :: ierr
    PetscInt :: idx_dv
    PetscInt :: i, j, istart, iend
    PetscInt :: iter_cnt

    ! increment iteration counter
    summary%iter_cnt = summary%iter_cnt + 1
    iter_cnt = summary%iter_cnt

    if (summary%nitermax > 1) then
      summary%itinner(iter_cnt) = n
      do i = 1, summary%convnmod
        summary%convdvmax(i, iter_cnt) = DZERO
        summary%convlocdv(i, iter_cnt) = 0
        summary%convrmax(i, iter_cnt) = DZERO
        summary%convlocr(i, iter_cnt) = 0
      end do
    end if

    ! get dv and dr per local model (readonly!)
    call VecGetArrayReadF90(dx, local_dx, ierr)
    CHKERRQ(ierr)
    do i = 1, summary%convnmod
      ! reset
      dvmax_model = DZERO
      idx_dv = 0
      ! get first and last model index
      istart = summary%model_bounds(i)
      iend = summary%model_bounds(i + 1) - 1
      do j = istart, iend
        if (abs(local_dx(j)) > abs(dvmax_model)) then
          dvmax_model = local_dx(j)
          idx_dv = j
        end if
      end do
      if (summary%nitermax > 1) then
        summary%convdvmax(i, iter_cnt) = dvmax_model
        summary%convlocdv(i, iter_cnt) = idx_dv
      end if
    end do
    call VecRestoreArrayF90(dx, local_dx, ierr)
    CHKERRQ(ierr)

  end subroutine fill_cnvg_summary_internal

  subroutine destroy(this)
    class(PetscCnvgCtxType) :: this
    ! local
    integer(I4B) :: ierr

    call VecDestroy(this%x_old, ierr)
    CHKERRQ(ierr)
    call VecDestroy(this%delta_x, ierr)
    CHKERRQ(ierr)
    call VecDestroy(this%residual, ierr)
    CHKERRQ(ierr)

  end subroutine destroy

end module PetscConvergenceModule
