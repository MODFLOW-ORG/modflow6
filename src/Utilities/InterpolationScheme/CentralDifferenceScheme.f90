module CentralDifferenceSchemeModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DHALF, DONE
  use IInterpolationSchemeModule, only: IInterpolationSchemeType, CoefficientsType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use IGradient, only: IGradientType

  implicit none
  private

  public :: CentralDifferenceSchemeType

  type, extends(IInterpolationSchemeType) :: CentralDifferenceSchemeType
    private
    class(DisBaseType), pointer :: dis
    type(TspFmiType), pointer :: fmi
    class(IGradientType), pointer :: gradient
  contains
    procedure :: compute
  end type CentralDifferenceSchemeType

  interface CentralDifferenceSchemeType
    module procedure Constructor
  end interface CentralDifferenceSchemeType

contains
  function Constructor(dis, fmi, gradient) result(interpolation_scheme)
    ! -- return
    type(CentralDifferenceSchemeType) :: interpolation_scheme
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
    class(CentralDifferenceSchemeType), target :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: iposnm
    real(DP), intent(in), dimension(:) :: phi
    ! -- local
    real(DP) :: lnm, lmn, omega

    ! -- calculate weight based on distances between nodes and the shared
    !    face of the connection
    if (this%dis%con%ihc(this%dis%con%jas(iposnm)) == 0) then
      ! -- vertical connection; assume cell is fully saturated
      lnm = DHALF * (this%dis%top(n) - this%dis%bot(n))
      lmn = DHALF * (this%dis%top(m) - this%dis%bot(m))
    else
      ! -- horizontal connection
      lnm = this%dis%con%cl1(this%dis%con%jas(iposnm))
      lmn = this%dis%con%cl2(this%dis%con%jas(iposnm))
    end if

    omega = lmn / (lnm + lmn)

    phi_face%c_n = omega
    phi_face%c_m = DONE - omega

  end function compute
end module CentralDifferenceSchemeModule