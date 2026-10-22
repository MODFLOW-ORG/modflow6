!> @brief NPF-backed implementation of ConductivityProviderType.
!!
!! Adapter that lets packages obtain a cell's hydraulic conductivity tensor
!! (and relative permeability) from NPF without depending on NPF internals.
!! Because GwfNpfType already extends NumericalPackageType (single inheritance),
!! the capability is exposed through this thin delegating adapter rather than by
!< making NPF extend the abstract provider directly.
module NpfConductivityProviderModule
  use KindModule, only: DP, I4B
  use ConductivityProviderModule, only: ConductivityProviderType
  use GwfNpfModule, only: GwfNpfType
  implicit none
  private
  public :: NpfConductivityProviderType
  public :: create_npf_conductivity_provider

  type, extends(ConductivityProviderType) :: NpfConductivityProviderType
    type(GwfNpfType), pointer :: npf => null() !< flow package supplying the conductivity
  contains
    procedure :: conductivity => npf_conductivity
    procedure :: krel => npf_krel
    procedure :: dkrel_dh => npf_dkrel_dh
  end type NpfConductivityProviderType

contains

  !> @brief Allocate an NPF-backed provider and return it as the abstract type.
  subroutine create_npf_conductivity_provider(provider, npf)
    class(ConductivityProviderType), pointer, intent(out) :: provider
    type(GwfNpfType), pointer, intent(in) :: npf
    ! local
    type(NpfConductivityProviderType), pointer :: npf_provider

    allocate (npf_provider)
    npf_provider%npf => npf
    provider => npf_provider
  end subroutine create_npf_conductivity_provider

  function npf_conductivity(this, n) result(k)
    class(NpfConductivityProviderType), intent(in) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(3, 3) :: k

    k = this%npf%calc_conductivity(n)
  end function npf_conductivity

  function npf_krel(this, n) result(kr)
    class(NpfConductivityProviderType), intent(in) :: this
    integer(I4B), intent(in) :: n
    real(DP) :: kr

    kr = this%npf%krel(n)
  end function npf_krel

  function npf_dkrel_dh(this, n) result(dkrdh)
    class(NpfConductivityProviderType), intent(in) :: this
    integer(I4B), intent(in) :: n
    real(DP) :: dkrdh

    dkrdh = this%npf%dkrdh(n)
  end function npf_dkrel_dh

end module NpfConductivityProviderModule
