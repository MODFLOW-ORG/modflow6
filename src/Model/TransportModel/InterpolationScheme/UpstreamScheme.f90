module UpstreamSchemeModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO
  use InterpolationSchemeInterfaceModule, only: InterpolationSchemeInterface, &
                                                CoefficientsType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType

  implicit none
  private

  public :: UpstreamSchemeType

  type, extends(InterpolationSchemeInterface) :: UpstreamSchemeType
    private
    class(DisBaseType), pointer :: dis
    type(TspFmiType), pointer :: fmi
  contains
    procedure :: compute
  end type UpstreamSchemeType

  interface UpstreamSchemeType
    module procedure constructor
  end interface UpstreamSchemeType

contains
  function constructor(dis, fmi) result(interpolation_scheme)
    ! -- return
    type(UpstreamSchemeType) :: interpolation_scheme
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    type(TspFmiType), pointer, intent(in) :: fmi

    interpolation_scheme%dis => dis
    interpolation_scheme%fmi => fmi

  end function constructor

  function compute(this, n, m, iposnm, phi) result(phi_face)
    !-- return
    type(CoefficientsType) :: phi_face
    ! -- dummy
    class(UpstreamSchemeType), target :: this
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
end module UpstreamSchemeModule
