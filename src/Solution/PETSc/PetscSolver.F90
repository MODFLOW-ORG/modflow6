module PetscSolverModule
#include <petsc/finclude/petscksp.h>
  use petscksp
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LINELENGTH, LENSOLUTIONNAME, DZERO
  use ListsModule, only: basesolutionlist
  use LinearSolverBaseModule
  use MatrixBaseModule
  use VectorBaseModule
  use PetscMatrixModule
  use PetscVectorModule
  use PetscConvergenceModule
  use PetscImsPreconditionerModule
  use ConvergenceSummaryModule
  use ImsLinearSettingsModule
  use SimVariablesModule, only: iout, simulation_mode, nr_procs, proc_id
  use SimModule, only: store_error, store_warning
  use FeatureFlagsModule, only: developmode

  implicit none
  private

  public :: create_petsc_solver

  !< Linear solver using PETSc KSP
  type, public, extends(LinearSolverBaseType) :: PetscSolverType
    KSP :: ksp_petsc !< the KSP linear solver object
    class(PetscMatrixType), pointer :: matrix => null() !< the system matrix in PETSc compatible format
    Mat, pointer :: mat_petsc => null() !< the PETSc internal matrix format

    type(ImsLinearSettingsType), pointer :: linear_settings => null() !< pointer to linear settings from IMS
    logical(LGP) :: use_ims_pc !< when true, use custom IMS-style preconditioning
    logical(LGP) :: use_ims_cnvgopt !< when true, use IMS convergence check in PETSc solve
    KSPType :: ksp_type !< the KSP solver type (CG, BCGS, ...)
    class(PetscCnvgCtxType), pointer :: petsc_ctx => null() !< context for the PETSc custom convergence check
    type(PcShellCtxType), pointer :: pc_context => null() !< context for the custom (IMS) precondioner
    type(ConvergenceSummaryType), pointer :: convergence_summary => null() !< data structure wrapping the convergence data
    character(len=LENSOLUTIONNAME + 1) :: option_prefix !< prefix for keys in petscrc database in case there are multiple solutions
  contains
    procedure :: initialize => petsc_initialize
    procedure :: solve => petsc_solve
    procedure :: print_summary => petsc_print_summary
    procedure :: destroy => petsc_destroy
    procedure :: create_matrix => petsc_create_matrix

    ! private
    procedure, private :: petsc_check_settings
    procedure, private :: get_options_mf6
    procedure, private :: create_ksp
    procedure, private :: create_convergence_check
    procedure, private :: set_ims_pc
    procedure, private :: print_vec
    procedure, private :: print_petsc_version
  end type PetscSolverType

contains

  !> @brief Create a PETSc solver object
  !<
  function create_petsc_solver(sln_name) result(solver)
    class(LinearSolverBaseType), pointer :: solver !< Uninitialized instance of the PETSc solver
    character(len=LENSOLUTIONNAME) :: sln_name !< the solution name
    ! local
    class(PetscSolverType), pointer :: petsc_solver
    character(len=LINELENGTH) :: errmsg

    allocate (petsc_solver)
    allocate (petsc_solver%petsc_ctx)
    petsc_solver%option_prefix = ""
    if (basesolutionlist%Count() > 1) then
      petsc_solver%option_prefix = trim(sln_name)//"_"
    end if

    solver => petsc_solver
    solver%name = sln_name

    if (simulation_mode /= 'PARALLEL') then
      write (errmsg, '(a,a)') 'PETSc solver not supported for run mode: ', &
        trim(simulation_mode)
      call store_error(errmsg, terminate=.true.)
    end if

  end function create_petsc_solver

  !> @brief Initialize PETSc KSP solver with
  !<  options from the petsc database file
  subroutine petsc_initialize(this, matrix, linear_settings, convergence_summary)
    class(PetscSolverType) :: this !< This solver instance
    class(MatrixBaseType), pointer :: matrix !< The solution matrix as KSP operator
    type(ImsLinearSettingsType), pointer :: linear_settings !< the settings for the linear solver from the .ims file
    type(ConvergenceSummaryType), pointer :: convergence_summary !< a convergence record for diagnostics

    this%linear_settings => linear_settings

    call this%print_petsc_version()

    this%mat_petsc => null()
    select type (pm => matrix)
    class is (PetscMatrixType)
      this%matrix => pm
      this%mat_petsc => pm%mat
    end select

    call this%petsc_check_settings(linear_settings)

    this%use_ims_cnvgopt = .true. ! use IMS convergence check, override with .petscrc
    this%use_ims_pc = .true. ! use IMS preconditioning, override with .petscrc
    allocate (this%pc_context)
    call this%pc_context%create(this%matrix, linear_settings)

    if (linear_settings%ilinmeth == CG_METHOD) then
      this%ksp_type = KSPCG
    else
      this%ksp_type = KSPBCGS
    end if

    ! get MODFLOW options from PETSc database file
    call this%get_options_mf6()

    ! create the solver object
    call this%create_ksp()

    ! Create custom convergence check
    call this%create_convergence_check(convergence_summary)

  end subroutine petsc_initialize

  subroutine petsc_check_settings(this, linear_settings)
    class(PetscSolverType) :: this
    type(ImsLinearSettingsType), pointer :: linear_settings
    ! local
    character(len=LINELENGTH) :: warnmsg, errmsg

    ! errors
    if (linear_settings%ilinmeth /= CG_METHOD .and. &
        linear_settings%ilinmeth /= BCGS_METHOD) then
      write (errmsg, '(a,a)') 'PETSc: unknown linear solver method in ', &
        this%name
      call store_error(errmsg, terminate=.true.)
    end if

    ! warnings
    if (linear_settings%iord > 0) then
      linear_settings%iord = 0
      write (warnmsg, '(a)') 'PETSc: IMS reordering not supported'
      call store_warning(warnmsg)
    end if
    if (linear_settings%iscl > 0) then
      linear_settings%iscl = 0
      write (warnmsg, '(a)') 'PETSc: IMS matrix scaling not supported'
      call store_warning(warnmsg)
    end if
    if (linear_settings%north > 0) then
      linear_settings%north = 0
      write (warnmsg, '(a)') 'PETSc: IMS orthogonalization not supported'
      call store_warning(warnmsg)
    end if

  end subroutine petsc_check_settings

  !> @brief Print PETSc version string from shared lib
  !<
  subroutine print_petsc_version(this)
    class(PetscSolverType) :: this
    ! local
    PetscErrorCode :: ierr
    PetscInt :: major, minor, subminor, release
    character(len=128) :: petsc_version, release_str

    call PetscGetVersionNumber(major, minor, subminor, release, ierr)
    CHKERRQ(ierr)

    if (release == 1) then
      release_str = "(release)"
    else
      release_str = "(unofficial)"
    end if
    write (petsc_version, '(i0,a,i0,a,i0,a,a)') &
      major, ".", minor, ".", subminor, " ", trim(release_str)
    write (iout, '(/,1x,4a,/)') "PETSc Linear Solver will be used for ", &
      trim(this%name), ": version ", petsc_version

  end subroutine print_petsc_version

  !> @brief Get the MODFLOW specific options from the PETSc database
  !<
  subroutine get_options_mf6(this)
    class(PetscSolverType) :: this
    ! local
    PetscErrorCode :: ierr
    logical(LGP) :: found
    logical(LGP) :: use_petsc_pc, use_petsc_cnvg

    use_petsc_pc = .false.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, trim(this%option_prefix), &
                             '-use_petsc_pc', use_petsc_pc, found, ierr)
    CHKERRQ(ierr)
    this%use_ims_pc = .not. use_petsc_pc

    use_petsc_cnvg = .false.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, trim(this%option_prefix), &
                             '-use_petsc_cnvg', use_petsc_cnvg, found, ierr)
    CHKERRQ(ierr)
    this%use_ims_cnvgopt = .not. use_petsc_cnvg

  end subroutine get_options_mf6

  !> @brief Create the PETSc KSP object
  !<
  subroutine create_ksp(this)
    class(PetscSolverType) :: this !< This solver instance
    ! local
    PetscErrorCode :: ierr
    PC :: pc

    call KSPCreate(PETSC_COMM_WORLD, this%ksp_petsc, ierr)
    CHKERRQ(ierr)

    ! set prefix for options database
    call KSPSetOptionsPrefix(this%ksp_petsc, trim(this%option_prefix), ierr)
    CHKERRQ(ierr)

    call KSPSetOperators(this%ksp_petsc, this%mat_petsc, this%mat_petsc, ierr)
    CHKERRQ(ierr)

    call KSPSetInitialGuessNonzero(this%ksp_petsc, .true., ierr)
    CHKERRQ(ierr)

    call KSPSetType(this%ksp_petsc, this%ksp_type, ierr)
    CHKERRQ(ierr)

    if (this%use_ims_pc) then
      call this%set_ims_pc()
    else
      call developmode('PETSc preconditioning is under development, install the &
            &nightly build or compile from source with IDEVELOPMODE = 1.')
      ! The PC options will be set from the .petscrc
      ! file in the call to KSPSetFromOptions below
      call KSPGetPC(this%ksp_petsc, pc, ierr)
      CHKERRQ(ierr)
      call PCSetFromOptions(pc, ierr)
      CHKERRQ(ierr)
    end if

    call KSPSetErrorIfNotConverged(this%ksp_petsc, .false., ierr)
    CHKERRQ(ierr)

    ! finally override these options from the
    ! optional .petscrc file
    call KSPSetFromOptions(this%ksp_petsc, ierr)
    CHKERRQ(ierr)

    call KSPSetup(this%ksp_petsc, ierr)
    CHKERRQ(ierr)

  end subroutine create_ksp

  !> @brief Set up a custom preconditioner following the ones
  !< we have in IMS, i.e. Modified ILU(T)
  subroutine set_ims_pc(this)
    class(PetscSolverType) :: this !< This solver instance
    ! local
    PC :: pc, sub_pc
    KSP, dimension(1) :: sub_ksp
    PetscInt :: n_local, n_first
    PetscErrorCode :: ierr

    call KSPGetPC(this%ksp_petsc, pc, ierr)
    CHKERRQ(ierr)
    call PCSetType(pc, PCBJACOBI, ierr)
    CHKERRQ(ierr)
    call PCSetUp(pc, ierr)
    CHKERRQ(ierr)
    call PCBJacobiGetSubKSP(pc, n_local, n_first, sub_ksp, ierr)
    CHKERRQ(ierr)
    call KSPGetPC(sub_ksp(1), sub_pc, ierr)
    CHKERRQ(ierr)
    call PCSetType(sub_pc, PCSHELL, ierr)
    CHKERRQ(ierr)
    call PCShellSetApply(sub_pc, pcshell_apply, ierr)
    CHKERRQ(ierr)
    call PCShellSetSetUp(sub_pc, pcshell_setup, ierr)
    CHKERRQ(ierr)
    call PCShellSetContext(sub_pc, this%pc_context, ierr)
    CHKERRQ(ierr)

  end subroutine set_ims_pc

  !> @brief Create and assign a custom convergence
  !< check for this solver
  subroutine create_convergence_check(this, convergence_summary)
    use IMSLinearBaseModule, only: ims_base_epfact
    class(PetscSolverType) :: this !< This solver instance
    type(ConvergenceSummaryType), pointer :: convergence_summary
    ! local
    PetscErrorCode :: ierr

    call this%petsc_ctx%create(this%mat_petsc, this%linear_settings, &
                               convergence_summary)

    if (.not. this%use_ims_cnvgopt) then
      ! use PETSc residual L2 norm for convergence
      call developmode('Using PETSc convergence is under development, install &
      &the nightly build or compile from source with IDEVELOPMODE = 1.')
      call KSPSetConvergenceTest(this%ksp_petsc, petsc_cnvg_check_internal, &
                                 this%petsc_ctx, PETSC_NULL_FUNCTION, ierr)
    else
      ! IMS convergence check
      call KSPSetConvergenceTest(this%ksp_petsc, petsc_cnvg_check, &
                                 this%petsc_ctx, PETSC_NULL_FUNCTION, ierr)
    end if
    CHKERRQ(ierr)

  end subroutine create_convergence_check

  subroutine petsc_solve(this, kiter, rhs, x, cnvg_summary)
    class(PetscSolverType) :: this
    integer(I4B) :: kiter
    class(VectorBaseType), pointer :: rhs
    class(VectorBaseType), pointer :: x
    type(ConvergenceSummaryType) :: cnvg_summary
    ! local
    PetscErrorCode :: ierr
    class(PetscVectorType), pointer :: rhs_petsc, x_petsc
    KSPConvergedReason :: cnvg_reason
    integer :: it_number
    character(len=LINELENGTH) :: errmsg

    rhs_petsc => null()
    select type (rhs)
    class is (PetscVectorType)
      rhs_petsc => rhs
    end select

    x_petsc => null()
    select type (x)
    class is (PetscVectorType)
      x_petsc => x
    end select

    this%iteration_number = 0
    this%is_converged = 0
    if (kiter == 1) then
      this%petsc_ctx%cnvg_summary%iter_cnt = 0
    end if

    ! update matrix coefficients
    call this%matrix%update()
    call KSPSolve(this%ksp_petsc, rhs_petsc%vec_impl, x_petsc%vec_impl, ierr)
    CHKERRQ(ierr)

    call KSPGetIterationNumber(this%ksp_petsc, it_number, ierr)
    CHKERRQ(ierr)
    this%iteration_number = it_number

    call KSPGetConvergedReason(this%ksp_petsc, cnvg_reason, ierr)
    CHKERRQ(ierr)
    if (cnvg_reason > 0) then
      if (this%petsc_ctx%icnvg_ims == -1) then
        ! move to next Picard iteration (e.g. with 'STRICT' option)
        this%is_converged = 0
      else
        ! linear convergence reached
        this%is_converged = 1
      end if
    end if

    if (cnvg_reason < 0) then
      if (cnvg_reason == KSP_DIVERGED_BREAKDOWN .or. &
          cnvg_reason == KSP_DIVERGED_ITS) then
        ! out of iterations, or rho/eta became zero,
        ! move to next Picard iteration and try again
        this%is_converged = 0
      else
        ! this is an solver error, terminate
        write (errmsg, '(1x,3a,i0)') "PETSc convergence failure in ", &
          trim(this%name), ": ", cnvg_reason
        call store_error(errmsg, terminate=.true.)
      end if
    end if

  end subroutine petsc_solve

  subroutine petsc_print_summary(this)
    class(PetscSolverType) :: this
    ! local
    character(len=128) :: ksp_str, pc_str, subpc_str, &
                          dvclose_str, rclose_str, relax_str, dtol_str
    character(len=128) :: ksp_logfile
    integer :: ierr
    PC :: pc
    PetscViewer :: ksp_viewer

    if (this%use_ims_pc) then
      call KSPGetType(this%ksp_petsc, ksp_str, ierr)
      CHKERRQ(ierr)
      call KSPGetPC(this%ksp_petsc, pc, ierr)
      CHKERRQ(ierr)
      call PCGetType(pc, pc_str, ierr)
      CHKERRQ(ierr)
      subpc_str = this%pc_context%ims_pc_type

      write (dvclose_str, '(e15.5)') this%linear_settings%dvclose
      write (rclose_str, '(e15.5)') this%linear_settings%rclose
      write (relax_str, '(e15.5)') this%linear_settings%relax
      write (dtol_str, '(e15.5)') this%linear_settings%droptol

      write (iout, '(/,1x,a)') "PETSc linear solver settings: "
      write (iout, '(1x,a)') repeat('-', 66)
      write (iout, '(1x,a,a)') "Linear acceleration method:   ", trim(ksp_str)
      write (iout, '(1x,a,a)') "Preconditioner type:          ", trim(pc_str)
      write (iout, '(1x,a,a)') "Sub-preconditioner type:      ", trim(subpc_str)
      write (iout, '(1x,a,i0)') "Maximum nr. of iterations:    ", &
        this%linear_settings%iter1
      write (iout, '(1x,a,a)') &
        "Dep. var. closure criterion:  ", trim(adjustl(dvclose_str))
      write (iout, '(1x,a,a)') &
        "Residual closure criterion:   ", trim(adjustl(rclose_str))
      if (this%use_ims_cnvgopt) then
        write (iout, '(1x,a,i0)') &
          "Residual convergence option:  ", this%linear_settings%icnvgopt
      else
        write (iout, '(1x,a)') &
          "Residual convergence option:  PETSc L2 norm"
      end if
      write (iout, '(1x,a,a)') &
        "Relaxation factor MILU(T):    ", trim(adjustl(relax_str))
      write (iout, '(1x,a,i0)') &
        "Fill level in factorization:  ", this%linear_settings%level
      write (iout, '(1x,a,a,/)') &
        "Drop tolerance level fill:    ", trim(adjustl(dtol_str))
    else
      ksp_logfile = trim(this%option_prefix)//"ksp_logview.txt"
      write (iout, '(/,1x,a)') "PETSc linear solver settings from .petscrc: "
      write (iout, '(1x,a)') repeat('-', 66)
      write (iout, '(1x,2a)') "see ", trim(ksp_logfile)

      ! collective write
      call PetscViewerASCIIOpen(PETSC_COMM_WORLD, ksp_logfile, ksp_viewer, ierr); 
      CHKERRQ(ierr)
      call PetscViewerPushFormat(ksp_viewer, PETSC_VIEWER_ASCII_INFO_DETAIL, ierr)
      CHKERRQ(ierr)
      call KspView(this%ksp_petsc, ksp_viewer, ierr)
      CHKERRQ(ierr)
      call PetscViewerDestroy(ksp_viewer, ierr)
      CHKERRQ(ierr)
    end if

  end subroutine petsc_print_summary

  subroutine petsc_destroy(this)
    class(PetscSolverType) :: this
    ! local
    PetscErrorCode :: ierr

    call KSPDestroy(this%ksp_petsc, ierr)
    CHKERRQ(ierr)

    ! delete context
    call this%petsc_ctx%destroy()
    deallocate (this%petsc_ctx)

    call this%pc_context%destroy()
    deallocate (this%pc_context)

  end subroutine petsc_destroy

  function petsc_create_matrix(this) result(matrix)
    class(PetscSolverType) :: this
    class(MatrixBaseType), pointer :: matrix
    ! local
    class(PetscMatrixType), pointer :: petsc_matrix

    allocate (petsc_matrix)
    matrix => petsc_matrix

  end function petsc_create_matrix

  subroutine print_vec(this, vec, vec_name, kiter)
    use TdisModule, only: nper, kstp
    class(PetscSolverType) :: this
    class(PetscVectorType) :: vec
    character(len=*) :: vec_name
    integer(I4B) :: kiter
    ! local
    PetscViewer :: viewer
    character(len=24) :: filename
    PetscErrorCode :: ierr

    write (filename, '(2a,i0,a,i0,a,i0,a)') vec_name, '_', nper, &
      '_', kstp, '_', kiter, '.txt'
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, viewer, ierr)
    CHKERRQ(ierr)
    call VecView(vec%vec_impl, viewer, ierr)
    CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr)
    CHKERRQ(ierr)

  end subroutine print_vec

end module PetscSolverModule
