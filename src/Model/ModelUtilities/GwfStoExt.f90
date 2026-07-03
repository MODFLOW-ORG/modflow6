module GwfStoExtModule
  use KindModule, only: I4B, LGP, DP
  use MatrixBaseModule, only: MatrixBaseType
  use BudgetModule, only: BudgetType
  implicit none
  private
  
  integer(I4B), public, parameter :: DEFAULT_STORAGE = 0
  integer(I4B), public, parameter :: UZR_STORAGE = 1
  integer(I4B), public, parameter :: SWI_STORAGE = 2
  integer(I4B), public, parameter :: MAX_EXT_STO_FORMS = 2 !< 'cause default is not part of the (external) storage formulations

  type, abstract, public :: GwfStoFormulationType
  contains
    procedure(is_active_if), deferred :: is_active
    procedure(fc_if), deferred :: fc
    procedure(fn_if), deferred :: fn
    procedure(cq_if), deferred :: cq
    procedure(bd_if), deferred :: bd
    procedure(save_flows_if), deferred :: save_flows
  end type GwfStoFormulationType

  !> @brief Container to allow arrays of polymorphic extension pointers
  !<
  type, public :: GwfStoFormContainerType
    logical(LGP) :: is_active = .false.
    class(GwfStoFormulationType), pointer :: form => null() !< the extended flow calculator
  end type GwfStoFormContainerType


  abstract interface
    function is_active_if(this, n) result(is_active)
      import GwfStoFormulationType, I4B, LGP
      class(GwfStoFormulationType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      logical(LGP) :: is_active
    end function
    subroutine fc_if(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
      import GwfStoFormulationType, MatrixBaseType, I4B, DP
      class(GwfStoFormulationType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
      real(DP), dimension(:), intent(inout) :: rhs
      integer(I4B), dimension(:), intent(in) :: idxglo
      real(DP), dimension(:), intent(in) :: h_old
      real(DP), dimension(:), intent(in) :: h_new
    end subroutine
    subroutine fn_if(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
      import GwfStoFormulationType, MatrixBaseType, I4B, DP
      class(GwfStoFormulationType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
      real(DP), dimension(:), intent(inout) :: rhs
      integer(I4B), dimension(:), intent(in) :: idxglo
      real(DP), dimension(:), intent(in) :: h_old
      real(DP), dimension(:), intent(in) :: h_new
    end subroutine
    subroutine cq_if(this, n, flowja, h_new, h_old)
      import GwfStoFormulationType, I4B, DP
      class(GwfStoFormulationType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      real(DP), dimension(:), intent(inout) :: flowja
      real(DP), dimension(:), intent(in) :: h_new
      real(DP), dimension(:), intent(in) :: h_old
    end subroutine
    subroutine bd_if(this, isuppress_output, model_budget)
      import GwfStoFormulationType, BudgetType, I4B
      class(GwfStoFormulationType), intent(inout) :: this
      integer(I4B), intent(in) :: isuppress_output
      type(BudgetType), intent(inout) :: model_budget
    end subroutine bd_if
    subroutine save_flows_if(this, iprint, ibinun)
      import GwfStoFormulationType, I4B
      class(GwfStoFormulationType), intent(inout) :: this
      integer(I4B), intent(in) :: iprint
      integer(I4B), intent(in) :: ibinun
    end subroutine save_flows_if
  end interface

end module GwfStoExtModule
