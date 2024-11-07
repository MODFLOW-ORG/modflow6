module GwfNpfExtModule
  use KindModule, only: I4B, LGP, DP
  use MatrixBaseModule, only: MatrixBaseType
  implicit none
  private

  type, abstract, public :: GwfNpfExtType
  contains
    procedure(is_active_if), deferred :: is_active
    procedure(fc_if), deferred :: fc
    procedure(fn_if), deferred :: fn
    procedure(cq_if), deferred :: cq
  end type GwfNpfExtType

  abstract interface
    function is_active_if(this, n, m) result(is_active)
      import GwfNpfExtType, I4B, LGP
      class(GwfNpfExtType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      integer(I4B), intent(in) :: m
      logical(LGP) :: is_active
    end function
    subroutine fc_if(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)
      import GwfNpfExtType, MatrixBaseType, I4B, DP
      class(GwfNpfExtType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      integer(I4B), intent(in) :: m
      integer(I4B), intent(in) :: ipos
      class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
      real(DP), dimension(:), intent(inout) :: rhs
      integer(I4B), dimension(:), intent(in) :: idxglo
      real(DP), dimension(:), intent(in) :: hnew
    end subroutine
    subroutine fn_if(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)
      import GwfNpfExtType, MatrixBaseType, I4B, DP
      class(GwfNpfExtType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      integer(I4B), intent(in) :: m
      integer(I4B), intent(in) :: ipos
      class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
      real(DP), dimension(:), intent(inout) :: rhs
      integer(I4B), dimension(:), intent(in) :: idxglo
      real(DP), dimension(:), intent(in) :: hnew
    end subroutine
    subroutine cq_if(this, n, m, ipos, flowja, h_new)
      import GwfNpfExtType, I4B, DP
      class(GwfNpfExtType), intent(inout) :: this
      integer(I4B), intent(in) :: n
      integer(I4B), intent(in) :: m
      integer(I4B), intent(in) :: ipos
      real(DP), dimension(:), intent(inout) :: flowja
      real(DP), dimension(:), intent(in) :: h_new
    end subroutine
  end interface

end module GwfNpfExtModule
