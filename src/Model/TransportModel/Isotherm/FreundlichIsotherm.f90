Module FreundlichIsothermModule

  use KindModule, only: DP, I4B
  use IsothermInterfaceModule, only: IsothermType

  Implicit None
  Private
  Public :: FreundlichIsothermType

  type, extends(IsothermType) :: FreundlichIsothermType
    real(DP), pointer, dimension(:) :: Kf => null() !< Freundlich constant
    real(DP), pointer, dimension(:) :: a => null() !< Freundlich exponent
  contains
    procedure :: value
    procedure :: derivative
  end type FreundlichIsothermType

  interface FreundlichIsothermType
    module procedure constructor
  end interface FreundlichIsothermType

contains

  function constructor(Kf, a) Result(isotherm)
    type(FreundlichIsothermType) :: isotherm
    ! -- dummy
    real(DP), pointer, dimension(:), intent(in) :: Kf
    real(DP), pointer, dimension(:), intent(in) :: a
    ! -- local
    isotherm%Kf => Kf
    isotherm%a => a

  end function constructor

  function value(this, c, n) result(val)
    class(FreundlichIsothermType), intent(in) :: this
    real(DP), dimension(:), intent(in) :: c
    integer(I4B), intent(in) :: n
    real(DP) :: val

    if (c(n) > 0.0_DP) then
      val = this%Kf(n) * c(n)**this%a(n)
    else
      val = 0.0_DP
    end if
  end function value

  function derivative(this, c, n) result(derv)
    class(FreundlichIsothermType), intent(in) :: this
    real(DP), dimension(:), intent(in) :: c
    integer(I4B), intent(in) :: n
    real(DP) :: derv

    if (c(n) > 0.0_DP) then
      derv = this%a(n) * this%Kf(n) * c(n)**(this%a(n) - 1.0_DP)
    else
      derv = 0.0_DP
    end if
  end function derivative

end module FreundlichIsothermModule
