!> @brief Abstract provider of cell hydraulic conductivity.
!!
!! Decouples boundary (and other) packages that need a cell's hydraulic
!! conductivity from the concrete flow package that owns the conductivity data
!! (e.g. NPF). A concrete provider is injected by the model; callers depend only
!! on this interface.
!!
!! The provider hands back the full saturated conductivity tensor K (a 3x3,
!! symmetric, model-coordinate tensor), leaving any directional projection or
!! conductance assembly to the consumer. Under an unsaturated (Richards)
!! formulation the current conductivity is krel * K; krel is 1 and dkrel_dh is 0
!! where the cell is fully saturated. A consumer computes the Darcy flux as
!! q = -krel * matmul(K, grad_h) and differentiates through krel/dkrel_dh for
!! Newton-Raphson.
!<
module ConductivityProviderModule
  use KindModule, only: DP, I4B
  implicit none
  private
  public :: ConductivityProviderType

  type, abstract :: ConductivityProviderType
  contains
    procedure(conductivity_if), deferred :: conductivity
    procedure(krel_if), deferred :: krel
    procedure(dkrel_dh_if), deferred :: dkrel_dh
  end type ConductivityProviderType

  abstract interface
    !> @brief Saturated hydraulic conductivity tensor of cell n, expressed in
    !! model coordinates as a 3x3 symmetric tensor (anisotropy and rotation
    !! resolved). Multiply by krel for the current conductivity under
    !< unsaturated conditions.
    function conductivity_if(this, n) result(k)
      import :: ConductivityProviderType, DP, I4B
      class(ConductivityProviderType), intent(in) :: this
      integer(I4B), intent(in) :: n !< reduced node number
      real(DP), dimension(3, 3) :: k !< conductivity tensor in model coordinates
    end function conductivity_if

    !> @brief Current relative permeability of cell n. Returns 1 for a fully
    !! saturated cell; less than 1 where an unsaturated flow formulation (e.g.
    !< UZR) is active. conductivity returns the saturated K, so callers multiply by this.
    function krel_if(this, n) result(kr)
      import :: ConductivityProviderType, DP, I4B
      class(ConductivityProviderType), intent(in) :: this
      integer(I4B), intent(in) :: n !< reduced node number
      real(DP) :: kr
    end function krel_if

    !> @brief Derivative of the relative permeability of cell n with respect to
    !! head. Returns 0 for saturated cells; nonzero where an unsaturated flow
    !< formulation supplies it. Used to build Newton-Raphson boundary terms.
    function dkrel_dh_if(this, n) result(dkrdh)
      import :: ConductivityProviderType, DP, I4B
      class(ConductivityProviderType), intent(in) :: this
      integer(I4B), intent(in) :: n !< reduced node number
      real(DP) :: dkrdh
    end function dkrel_dh_if
  end interface

end module ConductivityProviderModule
