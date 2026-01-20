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
    procedure(set_field_if), deferred :: set_field
  end type InterpolationSchemeInterface

  abstract interface

    !> @brief Compute interpolation coefficients for a face value
    !!
    !! This method computes the coefficients needed to interpolate a scalar field
    !! value at the face between two connected cells. The method returns coefficients
    !! that define how the face value depends on the cell-centered values.
    !<
    function compute_if(this, n, m, iposnm) result(phi_face)
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
    end function

    !> @brief Set the scalar field for which interpolation will be computed
    !!
    !! This method establishes a pointer to the scalar field data that will be
    !! used for subsequent interpolation computations. Implementations may also
    !! update any dependent cached data to ensure consistent results.
    !<
    subroutine set_field_if(this, phi)
      ! -- import
      import DP
      import InterpolationSchemeInterface
      ! -- dummy
      class(InterpolationSchemeInterface), target :: this
      real(DP), intent(in), dimension(:), pointer :: phi
    end subroutine

  end interface

end module InterpolationSchemeInterfaceModule
