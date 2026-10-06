!> @brief NPF-backed implementation of ConductanceProviderType.
!!
!! Adapter that lets packages obtain a directional effective hydraulic
!! conductivity from NPF without depending on NPF internals. Because
!! GwfNpfType already extends NumericalPackageType (single inheritance),
!! the capability is exposed through this thin delegating adapter rather
!< than by making NPF extend the abstract provider directly.
module NpfConductanceProviderModule
  use KindModule, only: DP, I4B
  use ConductanceProviderModule, only: ConductanceProviderType
  use GwfNpfModule, only: GwfNpfType
  implicit none
  private
  public :: NpfConductanceProviderType
  public :: create_npf_conductance_provider

  type, extends(ConductanceProviderType) :: NpfConductanceProviderType
    type(GwfNpfType), pointer :: npf => null() !< flow package supplying the conductivity
  contains
    procedure :: eff_hy => npf_eff_hy
  end type NpfConductanceProviderType

contains

  !> @brief Allocate an NPF-backed provider and return it as the abstract type.
  subroutine create_npf_conductance_provider(provider, npf)
    class(ConductanceProviderType), pointer, intent(out) :: provider
    type(GwfNpfType), pointer, intent(in) :: npf
    ! local
    type(NpfConductanceProviderType), pointer :: npf_provider

    allocate (npf_provider)
    npf_provider%npf => npf
    provider => npf_provider
  end subroutine create_npf_conductance_provider

  function npf_eff_hy(this, n, ihc, vg) result(hy)
    class(NpfConductanceProviderType), intent(in) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: ihc
    real(DP), dimension(3), intent(in) :: vg
    real(DP) :: hy

    hy = this%npf%calc_eff_hy(n, ihc, vg)
  end function npf_eff_hy

end module NpfConductanceProviderModule
