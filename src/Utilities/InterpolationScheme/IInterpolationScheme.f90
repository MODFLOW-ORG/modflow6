module IInterpolationSchemeModule
  use KindModule, only: DP, I4B

  implicit none
  private

  public :: IInterpolationSchemeType
  public :: CoefficientsType

  type :: CoefficientsType
    real(DP) :: c_n = 0.0_dp
    real(DP) :: c_m = 0.0_dp
    real(DP) :: rhs = 0.0_dp
  end type CoefficientsType

  type, abstract :: IInterpolationSchemeType
  contains
    procedure(compute_if), deferred :: compute
  end type IInterpolationSchemeType

  abstract interface

    function compute_if(this, n, m, iposnm, phi) result(phi_face)
      ! -- import
      import DP, I4B
      import IInterpolationSchemeType
      import CoefficientsType
      ! -- return
      type(CoefficientsType) :: phi_face
      ! -- dummy
      class(IInterpolationSchemeType), target :: this
      integer(I4B), intent(in) :: n
      integer(I4B), intent(in) :: m
      integer(I4B), intent(in) :: iposnm
      real(DP), intent(in), dimension(:) :: phi
    end function

  end interface

end module IInterpolationSchemeModule
