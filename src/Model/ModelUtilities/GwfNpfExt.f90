module GwfNpfFormulationModule
  use KindModule, only: I4B, LGP, DP
  use MatrixBaseModule, only: MatrixBaseType
  implicit none
  private

  integer(I4B), public, parameter :: DEFAULT_FLOW = 0
  integer(I4B), public, parameter :: UZR_FLOW = 1
  integer(I4B), public, parameter :: SIZE_FORMULATIONS = 2

  type, abstract, public :: GwfNpfFormulationType
  contains
    procedure(is_active_if), deferred :: is_active
    procedure(cf_if), deferred :: cf
    procedure(fc_if), deferred :: fc
    procedure(fn_if), deferred :: fn
    procedure(cq_if), deferred :: cq
  end type GwfNpfFormulationType

  !> @brief Container to allow arrays of polymorphic extension pointers
  !<
  type, public :: GwfNpfFormContainerType
    class(GwfNpfFormulationType), pointer :: form => null() !< the extended flow calculator
  end type GwfNpfFormContainerType

  abstract interface
    function is_active_if(this, n, m) result(is_active)
      import GwfNpfFormulationType, I4B, LGP
      class(GwfNpfFormulationType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      integer(I4B), intent(in) :: m
      logical(LGP) :: is_active
    end function
    subroutine cf_if(this, kiter, n)
      import GwfNpfFormulationType, I4B
      class(GwfNpfFormulationType), intent(inout) :: this
      integer(I4B), intent(in) :: kiter
      integer(I4B), intent(in) :: n
    end subroutine
    subroutine fc_if(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)
      import GwfNpfFormulationType, MatrixBaseType, I4B, DP
      class(GwfNpfFormulationType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      integer(I4B), intent(in) :: m
      integer(I4B), intent(in) :: ipos
      class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
      real(DP), dimension(:), intent(inout) :: rhs
      integer(I4B), dimension(:), intent(in) :: idxglo
      real(DP), dimension(:), intent(in) :: hnew
    end subroutine
    subroutine fn_if(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)
      import GwfNpfFormulationType, MatrixBaseType, I4B, DP
      class(GwfNpfFormulationType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      integer(I4B), intent(in) :: m
      integer(I4B), intent(in) :: ipos
      class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
      real(DP), dimension(:), intent(inout) :: rhs
      integer(I4B), dimension(:), intent(in) :: idxglo
      real(DP), dimension(:), intent(in) :: hnew
    end subroutine
    subroutine cq_if(this, n, m, ipos, flowja, h_new)
      import GwfNpfFormulationType, I4B, DP
      class(GwfNpfFormulationType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      integer(I4B), intent(in) :: m
      integer(I4B), intent(in) :: ipos
      real(DP), dimension(:), intent(inout) :: flowja
      real(DP), dimension(:), intent(in) :: h_new
    end subroutine
  end interface

end module GwfNpfFormulationModule
