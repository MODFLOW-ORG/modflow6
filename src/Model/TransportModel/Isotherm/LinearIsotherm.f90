module LinearIsothermsModule

  use KindModule, only: DP, I4B
  use IsothermInterfaceModule, only: IsothermType

  implicit none
  private
  public :: LinearIsothermType

  !> @brief Linear (Kd) isotherm implementation of `IsothermType`.
  !>
  !> Sorbed concentration is computed as cs = Kd*c.
  !<
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
  !> @brief Constructor for Linear isotherm
  !<
  function constructor(Kd) Result(isotherm)
    ! -- return
    type(LinearIsothermType) :: isotherm
    ! -- dummy
    real(DP), pointer, dimension(:), intent(in) :: Kd !< distribution coefficient
    ! -- local

    isotherm%Kd => Kd

  end function constructor

  !> @brief Evaluate the isotherm at a given node
  !<
  function value(this, c, n) result(val)
    ! -- return
    real(DP) :: val !< isotherm value
    ! -- dummy
    class(LinearIsothermType), intent(in) :: this
    real(DP), dimension(:), intent(in) :: c !< concentration array
    integer(I4B), intent(in) :: n !< node index

    val = this%Kd(n) * c(n)
  end function value

  !> @brief Evaluate derivative of the isotherm at a given node
  !<
  function derivative(this, c, n) result(derv)
    ! -- return
    real(DP) :: derv !< derivative d(value)/dc evaluated at c
    ! -- dummy
    class(LinearIsothermType), intent(in) :: this
    real(DP), dimension(:), intent(in) :: c !< concentration array
    integer(I4B), intent(in) :: n !< node index
    ! -- local

    derv = this%Kd(n)
  end function derivative

end module LinearIsothermsModule
