Module FreundlichIsothermModule

  use KindModule, only: DP, I4B
  use IsothermInterfaceModule, only: IsothermType

  Implicit None
  Private
  Public :: FreundlichIsothermType

  !> @brief Freundlich isotherm implementation of `IsothermType`.
  !>
  !> Sorbed concentration is cs = Kf*c^a.
  !<
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

  !> @brief Constructor for Freundlich isotherm
  !<
  function constructor(Kf, a) Result(isotherm)
    type(FreundlichIsothermType) :: isotherm
    ! -- dummy
    real(DP), pointer, dimension(:), intent(in) :: Kf
    real(DP), pointer, dimension(:), intent(in) :: a
    ! -- local
    isotherm%Kf => Kf
    isotherm%a => a

  end function constructor

  !> @brief Evaluate the isotherm at a given node
  !<
  function value(this, c, n) result(val)
    ! -- return
    real(DP) :: val !< isotherm value
    ! -- dummy
    class(FreundlichIsothermType), intent(in) :: this
    real(DP), dimension(:), intent(in) :: c !< concentration array
    integer(I4B), intent(in) :: n !< node index

    if (c(n) > 0.0_DP) then
      val = this%Kf(n) * c(n)**this%a(n)
    else
      val = 0.0_DP
    end if
  end function value

  !> @brief Evaluate derivative of the isotherm at a given node
  !<
  function derivative(this, c, n) result(derv)
    ! -- return
    real(DP) :: derv !< derivative d(value)/dc evaluated at c
    ! -- dummy
    class(FreundlichIsothermType), intent(in) :: this
    real(DP), dimension(:), intent(in) :: c !< concentration array
    integer(I4B), intent(in) :: n !< node index

    if (c(n) > 0.0_DP) then
      derv = this%a(n) * this%Kf(n) * c(n)**(this%a(n) - 1.0_DP)
    else
      derv = 0.0_DP
    end if
  end function derivative

end module FreundlichIsothermModule
