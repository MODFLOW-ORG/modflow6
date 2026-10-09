!> @brief Abstract provider of directional effective hydraulic conductivity.
!!
!! Decouples boundary (and other) packages that need a face/segment
!! conductance from the concrete flow package that owns the conductivity
!! data (e.g. NPF). A concrete provider is injected by the model; callers
!! depend only on this interface.
!<
module ConductanceProviderModule
  use KindModule, only: DP, I4B
  implicit none
  private
  public :: ConductanceProviderType

  type, abstract :: ConductanceProviderType
  contains
    procedure(eff_hy_if), deferred :: eff_hy
    procedure(krel_if), deferred :: krel
    procedure(dkrel_dh_if), deferred :: dkrel_dh
  end type ConductanceProviderType

  abstract interface
    !> @brief Effective hydraulic conductivity of cell n in the direction of
    !! the unit vector vg. ihc selects the vertical (0) or horizontal (1)
    !< anisotropy branch, consistent with MODFLOW connection conventions.
    function eff_hy_if(this, n, ihc, vg) result(hy)
      import :: ConductanceProviderType, DP, I4B
      class(ConductanceProviderType), intent(in) :: this
      integer(I4B), intent(in) :: n !< reduced node number
      integer(I4B), intent(in) :: ihc !< horizontal connection flag
      real(DP), dimension(3), intent(in) :: vg !< unit direction vector
      real(DP) :: hy
    end function eff_hy_if

    !> @brief Current relative permeability of cell n. Returns 1 for a fully
    !! saturated cell; less than 1 where an unsaturated flow formulation (e.g.
    !< UZR) is active. eff_hy returns the saturated K, so callers multiply by this.
    function krel_if(this, n) result(kr)
      import :: ConductanceProviderType, DP, I4B
      class(ConductanceProviderType), intent(in) :: this
      integer(I4B), intent(in) :: n !< reduced node number
      real(DP) :: kr
    end function krel_if

    !> @brief Derivative of the relative permeability of cell n with respect to
    !! head. Returns 0 for saturated cells; nonzero where an unsaturated flow
    !< formulation supplies it. Used to build Newton-Raphson boundary terms.
    function dkrel_dh_if(this, n) result(dkrdh)
      import :: ConductanceProviderType, DP, I4B
      class(ConductanceProviderType), intent(in) :: this
      integer(I4B), intent(in) :: n !< reduced node number
      real(DP) :: dkrdh
    end function dkrel_dh_if
  end interface

end module ConductanceProviderModule
