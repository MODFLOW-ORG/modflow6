module IsothermInterfaceModule
  use KindModule, only: DP, I4B

  implicit none
  private
  public :: IsothermType

  !> @brief Interface for isotherms.
  !<
  type, abstract :: IsothermType
  contains
    procedure(value_if), deferred :: value !< Evaluate isotherm value at node n
    procedure(derivative_if), deferred :: derivative !< Evaluate derivative d(value)/dc at node n
  end type IsothermType

  abstract interface
    !> @brief Evaluate the isotherm at a given node
    !!
    !! @return isotherm value
    !<
    function value_if(this, c, n) result(val)
      ! -- import
      import :: IsothermType, DP, I4B
      ! -- return
      real(DP) :: val
      ! -- dummy
      class(IsothermType), intent(in) :: this !< isotherm object
      real(DP), dimension(:), intent(in) :: c !< concentration array
      integer(I4B), intent(in) :: n !< node index
    end function value_if

    !> @brief Evaluate derivative of the isotherm at a given node
    !!
    !! @return derivative d(value)/dc evaluated at c
    !<
    function derivative_if(this, c, n) result(derv)
      ! -- import
      import :: IsothermType, DP, I4B
      ! -- return
      real(DP) :: derv
      ! -- dummy
      class(IsothermType), intent(in) :: this !< isotherm object
      real(DP), dimension(:), intent(in) :: c !< concentration array
      integer(I4B), intent(in) :: n !< node index
    end function derivative_if
  end interface

end module IsothermInterfaceModule
