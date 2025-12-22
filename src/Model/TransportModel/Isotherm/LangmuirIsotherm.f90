module LangmuirIsothermModule

  use KindModule, only: DP, I4B
  use IsothermInterfaceModule, only: IsothermType

  implicit none
  private
  public :: LangmuirIsothermType

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
  function constructor(Kl, Sbar) Result(isotherm)
    type(LangmuirIsothermType) :: isotherm
    ! -- dummy
    real(DP), pointer, dimension(:), intent(in) :: Kl
    real(DP), pointer, dimension(:), intent(in) :: Sbar
    ! -- local
    isotherm%Kl => Kl
    isotherm%Sbar => Sbar

  end function constructor

  function value(this, c, n) result(val)
    class(LangmuirIsothermType), intent(in) :: this
    real(DP), dimension(:), intent(in) :: c
    integer(I4B), intent(in) :: n
    real(DP) :: val
    real(DP) :: denom

    if (c(n) > 0.0_DP) then
      denom = 1.0_DP + this%Kl(n) * c(n)
      val = (this%Sbar(n) * this%Kl(n) * c(n)) / denom
    else
      val = 0.0_DP
    end if
  end function value

  function derivative(this, c, n) result(derv)
    class(LangmuirIsothermType), intent(in) :: this
    real(DP), dimension(:), intent(in) :: c
    integer(I4B), intent(in) :: n
    real(DP) :: derv
    real(DP) :: denom

    if (c(n) > 0.0_DP) then
      denom = (1.0_DP + this%Kl(n) * c(n))**2.0_dp
      derv = (this%Sbar(n) * this%Kl(n)) / denom
    else
      derv = 0.0_DP
    end if
  end function derivative

end module LangmuirIsothermModule
