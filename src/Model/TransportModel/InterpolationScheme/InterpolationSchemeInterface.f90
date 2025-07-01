module InterpolationSchemeInterfaceModule
  use KindModule, only: DP, I4B

  implicit none
  private

  public :: InterpolationSchemeInterface
  public :: CoefficientsType

  type :: CoefficientsType
    real(DP) :: c_n = 0.0_dp
    real(DP) :: c_m = 0.0_dp
    real(DP) :: rhs = 0.0_dp
  end type CoefficientsType

  type, abstract :: InterpolationSchemeInterface
  contains
    procedure(compute_if), deferred :: compute
  end type InterpolationSchemeInterface

  abstract interface

    function compute_if(this, n, m, iposnm, phi) result(phi_face)
      ! -- import
      import DP, I4B
      import InterpolationSchemeInterface
      import CoefficientsType
      ! -- return
      type(CoefficientsType) :: phi_face
      ! -- dummy
      class(InterpolationSchemeInterface), target :: this
      integer(I4B), intent(in) :: n
      integer(I4B), intent(in) :: m
      integer(I4B), intent(in) :: iposnm
      real(DP), intent(in), dimension(:) :: phi
    end function

  end interface

end module InterpolationSchemeInterfaceModule
