module EtUptakeFeddesModule
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE
  use EtUptakeModelModule, only: EtUptakeModelType
  implicit none
  private

  public :: FeddesUptakeType

  !> @brief Feddes (1978) piecewise-linear pressure-head stress response
  !<
  ! alpha(psi) is zero outside the pressure-head window [h4, h1], unity in the
  ! optimal range [h3, h2], and varies linearly on the wet (h2..h1, oxygen
  ! stress) and dry (h4..h3, water stress) limbs; thresholds satisfy
  ! h1 > h2 > h3 > h4.
  type, extends(EtUptakeModelType) :: FeddesUptakeType
    real(DP) :: h1 = DZERO !< anoxia point (uptake zero at and above)
    real(DP) :: h2 = DZERO !< wet end of optimal range
    real(DP) :: h3 = DZERO !< dry end of optimal range
    real(DP) :: h4 = DZERO !< wilting point (uptake zero at and below)
  contains
    procedure :: stress => feddes_stress
    procedure :: fc_terms => feddes_fc_terms
  end type FeddesUptakeType

contains

  !> @brief Feddes stress response alpha in [0, 1]
  !<
  function feddes_stress(this, psi) result(alpha)
    class(FeddesUptakeType), intent(in) :: this
    real(DP), intent(in) :: psi
    real(DP) :: alpha

    if (psi >= this%h1 .or. psi <= this%h4) then
      alpha = DZERO
    else if (psi >= this%h2) then
      ! wet limb: rises from 0 at h1 to 1 at h2
      alpha = (psi - this%h1) / (this%h2 - this%h1)
    else if (psi >= this%h3) then
      ! optimal range
      alpha = DONE
    else
      ! dry limb: falls from 1 at h3 to 0 at h4
      alpha = (psi - this%h4) / (this%h3 - this%h4)
    end if

  end function feddes_stress

  !> @brief Linearized matrix terms for the active stress limb
  !<
  subroutine feddes_fc_terms(this, h, z, c, hcof, rhs)
    class(FeddesUptakeType), intent(in) :: this
    real(DP), intent(in) :: h
    real(DP), intent(in) :: z
    real(DP), intent(in) :: c
    real(DP), intent(out) :: hcof
    real(DP), intent(out) :: rhs
    ! local
    real(DP) :: psi
    real(DP) :: psi0, psi1

    hcof = DZERO
    rhs = DZERO
    psi = h - z

    if (psi >= this%h1 .or. psi <= this%h4) then
      ! no uptake (anoxia or wilting)
      return
    else if (psi >= this%h2) then
      ! wet limb: alpha = 0 at psi0 = h1, alpha = 1 at psi1 = h2
      psi0 = this%h1
      psi1 = this%h2
    else if (psi >= this%h3) then
      ! optimal range: alpha = 1, constant
      rhs = c
      return
    else
      ! dry limb: alpha = 0 at psi0 = h4, alpha = 1 at psi1 = h3
      psi0 = this%h4
      psi1 = this%h3
    end if

    ! linear limb, exact within the limb: q = -c * (psi - psi0) / (psi1 - psi0)
    hcof = -c / (psi1 - psi0)
    rhs = hcof * (z + psi0)

  end subroutine feddes_fc_terms

end module EtUptakeFeddesModule
