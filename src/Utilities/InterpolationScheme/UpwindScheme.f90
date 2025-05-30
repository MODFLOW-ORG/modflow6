module UpwindSchemeModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO
  use IInterpolationSchemeModule, only: IInterpolationSchemeType, CoefficientsType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use IGradient, only: IGradientType

  implicit none
  private

  public :: UpwindSchemeType

  type, extends(IInterpolationSchemeType) :: UpwindSchemeType
    private
    class(DisBaseType), pointer :: dis
    type(TspFmiType), pointer :: fmi
    class(IGradientType), pointer :: gradient
  contains
    procedure :: compute
  end type UpwindSchemeType

  interface UpwindSchemeType
    module procedure Constructor
  end interface UpwindSchemeType

contains
  function Constructor(dis, fmi, gradient) result(interpolation_scheme)
    ! -- return
    type(UpwindSchemeType) :: interpolation_scheme
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    type(TspFmiType), pointer, intent(in) :: fmi
    class(IGradientType), allocatable, target, intent(in) :: gradient

    interpolation_scheme%dis => dis
    interpolation_scheme%fmi => fmi
    interpolation_scheme%gradient => gradient

  end function Constructor

  function compute(this, n, m, iposnm, phi) result(phi_face)
    !-- return
    type(CoefficientsType) :: phi_face
    ! -- dummy
    class(UpwindSchemeType), target :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: iposnm
    real(DP), intent(in), dimension(:) :: phi
    ! -- local
    real(DP) :: qnm

    ! -- Compute the coefficients for the upwind scheme
    qnm = this%fmi%gwfflowja(iposnm)

    if (qnm < DZERO) then
      phi_face%c_n = 1.0_dp
    else
      phi_face%c_m = 1.0_dp
    end if

  end function compute
end module UpwindSchemeModule