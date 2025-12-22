module LinearIsothermsModule

  use KindModule, only: DP, I4B
  use IsothermInterfaceModule, only: IsothermType

  implicit none
  private
  public :: LinearIsothermType

  type, extends(IsothermType) :: LinearIsothermType
    real(DP), pointer, dimension(:) :: Kd => null() !< distribution coefficient
  contains
    procedure :: value
    procedure :: derivative
  end type LinearIsothermType

  interface LinearIsothermType
    module procedure constructor
  end interface LinearIsothermType

contains
  function constructor(Kd) Result(isotherm)
    type(LinearIsothermType) :: isotherm
    ! -- dummy
    real(DP), pointer, dimension(:), intent(in) :: Kd
    ! -- local
    isotherm%Kd => Kd

  end function constructor

  function value(this, c, n) result(val)
    class(LinearIsothermType), intent(in) :: this
    real(DP), dimension(:), intent(in) :: c
    integer(I4B), intent(in) :: n
    real(DP) :: val

    val = this%Kd(n) * c(n)
  end function value

  function derivative(this, c, n) result(derv)
    class(LinearIsothermType), intent(in) :: this
    real(DP), dimension(:), intent(in) :: c
    integer(I4B), intent(in) :: n
    real(DP) :: derv

    derv = this%Kd(n)
  end function derivative

end module LinearIsothermsModule
