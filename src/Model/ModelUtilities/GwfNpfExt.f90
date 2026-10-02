module GwfNpfFormulationModule
  use KindModule, only: I4B, DP
  use MatrixBaseModule, only: MatrixBaseType
  implicit none
  private

  integer(I4B), public, parameter :: DEFAULT_FLOW = 0
  integer(I4B), public, parameter :: UZR_FLOW = 1
  integer(I4B), public, parameter :: SWI_FLOW = 2
  integer(I4B), public, parameter :: MAX_EXT_FLOW_FORMS = 2

  !> @brief Abstract flow formulation that additively contributes terms
  !!
  !! A formulation owns its own traversal of cells/connections and decides,
  !! per element, whether to add terms. Phases a formulation does not
  !! implement fall back to a no-op so formulations compose additively.
  !<
  type, abstract, public :: GwfNpfFormulationType
  contains
    procedure(fc_if), deferred :: fc
    procedure :: cf => cf_noop
    procedure :: fn => fn_noop
    procedure :: cq => cq_noop
  end type GwfNpfFormulationType

  !> @brief Container to allow arrays of polymorphic extension pointers
  !<
  type, public :: GwfNpfFormContainerType
    class(GwfNpfFormulationType), pointer :: form => null() !< the extended flow calculator
  end type GwfNpfFormContainerType

  abstract interface
    !> @brief Fill coefficients: formulation loops all connections itself
    !<
    subroutine fc_if(this, kiter, matrix_sln, idxglo, rhs, hnew)
      import GwfNpfFormulationType, MatrixBaseType, I4B, DP
      class(GwfNpfFormulationType), intent(inout) :: this
      integer(I4B), intent(in) :: kiter
      class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
      integer(I4B), dimension(:), intent(in) :: idxglo
      real(DP), dimension(:), intent(inout) :: rhs
      real(DP), dimension(:), intent(inout) :: hnew
    end subroutine
  end interface

contains

  !> @brief No-op coefficient calculation; formulation loops cells itself
  !<
  subroutine cf_noop(this, kiter)
    class(GwfNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: kiter
  end subroutine cf_noop

  !> @brief No-op newton terms; formulation loops connections itself
  !<
  subroutine fn_noop(this, kiter, matrix_sln, idxglo, rhs, hnew)
    class(GwfNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(inout) :: rhs
    real(DP), dimension(:), intent(inout) :: hnew
  end subroutine fn_noop

  !> @brief No-op flow calculation; formulation loops connections itself
  !<
  subroutine cq_noop(this, hnew, flowja)
    class(GwfNpfFormulationType), intent(inout) :: this
    real(DP), dimension(:), intent(inout) :: hnew
    real(DP), dimension(:), intent(inout) :: flowja
  end subroutine cq_noop

end module GwfNpfFormulationModule
