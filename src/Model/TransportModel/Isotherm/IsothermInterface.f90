module IsothermInterfaceModule
  use KindModule, only: DP, I4B

  implicit none
  private
  public :: IsothermType

  type, abstract :: IsothermType
  contains
    procedure(value_interface), deferred :: value
    procedure(derivative_interface), deferred :: derivative
  end type IsothermType

  abstract interface
    function value_interface(this, c, n) result(val)
      import :: IsothermType, DP, I4B
      class(IsothermType), intent(in) :: this
      real(DP), dimension(:), intent(in) :: c
      integer(I4B), intent(in) :: n
      real(DP) :: val
    end function value_interface

    function derivative_interface(this, c, n) result(derv)
      import :: IsothermType, DP, I4B
      class(IsothermType), intent(in) :: this
      real(DP), dimension(:), intent(in) :: c
      integer(I4B), intent(in) :: n
      real(DP) :: derv
    end function derivative_interface
  end interface

end module IsothermInterfaceModule
