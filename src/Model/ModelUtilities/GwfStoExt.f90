module GwfStoExtModule
  use KindModule, only: I4B, LGP, DP
  use MatrixBaseModule, only: MatrixBaseType
  implicit none
  private

  type, abstract, public :: GwfStoExtType
  contains
    procedure(is_active_if), deferred :: is_active
    procedure(fc_if), deferred :: fc
    procedure(fn_if), deferred :: fn
    procedure(cq_if), deferred :: cq
  end type GwfStoExtType

  abstract interface
    function is_active_if(this, n) result(is_active)
      import GwfStoExtType, I4B, LGP
      class(GwfStoExtType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      logical(LGP) :: is_active
    end function
    subroutine fc_if(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
      import GwfStoExtType, MatrixBaseType, I4B, DP
      class(GwfStoExtType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
      real(DP), dimension(:), intent(inout) :: rhs
      integer(I4B), dimension(:), intent(in) :: idxglo
      real(DP), dimension(:), intent(in) :: h_old
      real(DP), dimension(:), intent(in) :: h_new
    end subroutine
    subroutine fn_if(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
      import GwfStoExtType, MatrixBaseType, I4B, DP
      class(GwfStoExtType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
      real(DP), dimension(:), intent(inout) :: rhs
      integer(I4B), dimension(:), intent(in) :: idxglo
      real(DP), dimension(:), intent(in) :: h_old
      real(DP), dimension(:), intent(in) :: h_new
    end subroutine
    subroutine cq_if(this, n, flowja, h_new, h_old)
      import GwfStoExtType, I4B, DP
      class(GwfStoExtType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      real(DP), dimension(:), intent(inout) :: flowja
      real(DP), dimension(:), intent(in) :: h_new
      real(DP), dimension(:), intent(in) :: h_old
    end subroutine
  end interface

end module GwfStoExtModule
