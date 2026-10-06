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
  end interface

end module ConductanceProviderModule
