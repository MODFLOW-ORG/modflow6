Module FreundlichIsothermModule

  use KindModule, only: DP, I4B
  use IsothermInterfaceModule, only: IsothermType
  use ConstantsModule, only: DZERO

  Implicit None
  Private
  Public :: FreundlichIsothermType

  !> @brief Freundlich isotherm implementation of `IsothermType`.
  !>
  !> Sorbed concentration is cs = Kf*c^a.
  !>
  !> However, this expression has a singularity at c = 0 when a < 1,
  !> leading to infinite derivative. To avoid this, the Freundlich isotherm
  !> is modified as follows.
  !>
  !> Modified Sorbed concentration is cs = Kf*(c + eps)^a - Kf*eps^a
  !> where eps = (K/(a*Kf))^(1/(a-1)) and K is a large constant (default 10).
  !> This ensures that the derivative at c = 0 is below K.
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
    real(DP), parameter :: K = 10.0_dp !< constant to limit derivative at c=0
    real(DP) :: eps !< small concentration offset

    eps = (K / (this%a(n) * this%Kf(n)))**(1.0_dp / (this%a(n) - 1.0_dp))

    if (c(n) > DZERO) then
      val = this%Kf(n) * (c(n) + eps)**this%a(n) - this%Kf(n) * eps**this%a(n)
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
    real(DP), parameter :: K = 10.0_dp !< constant to limit derivative at c=0
    real(DP) :: eps !< small concentration offset

    eps = (K / (this%a(n) * this%Kf(n)))**(1.0_dp / (this%a(n) - 1.0_dp))

    if (c(n) > DZERO) then
      derv = this%a(n) * this%Kf(n) * ((c(n) + eps)**(this%a(n) - 1.0_DP))
    else
      derv = 0.0_DP
    end if
  end function derivative

end module FreundlichIsothermModule
