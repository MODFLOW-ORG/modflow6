module LinearSolverBaseModule
  use KindModule, only: I4B, DP
  use ConstantsModule, only: LENSOLUTIONNAME
  use MatrixBaseModule
  use VectorBaseModule
  use ImsLinearSettingsModule
  use ConvergenceSummaryModule

  implicit none
  private

  !> @brief Abstract type for linear solver
  !!
  !! This serves as the base type for our solvers:
  !! sequential, parallel, petsc, block solver, ...
  !<
  type, public, abstract :: LinearSolverBaseType
    character(len=LENSOLUTIONNAME) :: name
    integer(I4B) :: nitermax
    integer(I4B) :: iteration_number
    integer(I4B) :: is_converged !< linear convergence status (INNER_*)
  contains
    procedure(initialize_if), deferred :: initialize
    procedure(print_summary_if), deferred :: print_summary
    procedure(solve_if), deferred :: solve
    procedure(destroy_if), deferred :: destroy

    procedure(create_matrix_if), deferred :: create_matrix
  end type LinearSolverBaseType

  abstract interface
    subroutine initialize_if(this, matrix, linear_settings, convergence_summary)
      import LinearSolverBaseType, MatrixBaseType, &
        ImsLinearSettingsType, ConvergenceSummaryType
      class(LinearSolverBaseType) :: this
      class(MatrixBaseType), pointer :: matrix
      type(ImsLinearSettingsType), pointer :: linear_settings
      type(ConvergenceSummaryType), pointer :: convergence_summary
    end subroutine
    subroutine print_summary_if(this)
      import LinearSolverBaseType
      class(LinearSolverBaseType) :: this
    end subroutine
    subroutine solve_if(this, kiter, rhs, x, cnvg_summary)
      import LinearSolverBaseType, I4B, VectorBaseType, ConvergenceSummaryType
      class(LinearSolverBaseType) :: this
      integer(I4B) :: kiter
      class(VectorBaseType), pointer :: rhs
      class(VectorBaseType), pointer :: x
      type(ConvergenceSummaryType) :: cnvg_summary
    end subroutine
    subroutine destroy_if(this)
      import LinearSolverBaseType
      class(LinearSolverBaseType) :: this
    end subroutine
    function create_matrix_if(this) result(matrix)
      import LinearSolverBaseType, MatrixBaseType
      class(LinearSolverBaseType) :: this
      class(MatrixBaseType), pointer :: matrix
    end function
  end interface

end module LinearSolverBaseModule
