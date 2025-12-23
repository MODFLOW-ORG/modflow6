module LangmuirIsothermModule

  use KindModule, only: DP, I4B
  use IsothermInterfaceModule, only: IsothermType

  implicit none
  private
  public :: LangmuirIsothermType

  !> @brief Langmuir isotherm implementation of `IsothermType`.
  !!
  !! Sorbed concentration is cs = (Sbar*Kl*c)/(1 + Kl*c).
  !<
  type, extends(IsothermType) :: LangmuirIsothermType
    real(DP), pointer, dimension(:) :: Kl => null() !< Langmuir constant
    real(DP), pointer, dimension(:) :: Sbar => null() !< Total concentration of sorption sites
  contains
    procedure :: value
    procedure :: derivative
  end type LangmuirIsothermType

  interface LangmuirIsothermType
    module procedure constructor
  end interface LangmuirIsothermType

contains
  !> @brief Constructor for Langmuir isotherm
  !<
  function constructor(Kl, Sbar) Result(isotherm)
    ! -- return
    type(LangmuirIsothermType) :: isotherm
    ! -- dummy
    real(DP), pointer, dimension(:), intent(in) :: Kl !< Langmuir constant
    real(DP), pointer, dimension(:), intent(in) :: Sbar !< Total concentration of sorption sites
    ! -- local

    isotherm%Kl => Kl
    isotherm%Sbar => Sbar

  end function constructor

  !> @brief Evaluate the isotherm at a given node
  !<
  function value(this, c, n) result(val)
    ! -- return
    real(DP) :: val !< isotherm value
    ! -- dummy
    class(LangmuirIsothermType), intent(in) :: this
    real(DP), dimension(:), intent(in) :: c !< concentration array
    integer(I4B), intent(in) :: n !< node index
    ! -- local
    real(DP) :: denom

    if (c(n) > 0.0_DP) then
      denom = 1.0_DP + this%Kl(n) * c(n)
      val = (this%Sbar(n) * this%Kl(n) * c(n)) / denom
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
    class(LangmuirIsothermType), intent(in) :: this
    real(DP), dimension(:), intent(in) :: c !< concentration array
    integer(I4B), intent(in) :: n !< node index
    ! -- local
    real(DP) :: denom

    if (c(n) > 0.0_DP) then
      denom = (1.0_DP + this%Kl(n) * c(n))**2.0_dp
      derv = (this%Sbar(n) * this%Kl(n)) / denom
    else
      derv = 0.0_DP
    end if
  end function derivative

end module LangmuirIsothermModule
