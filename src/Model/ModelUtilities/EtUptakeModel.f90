module EtUptakeModelModule
  use KindModule, only: I4B, DP
  implicit none
  private

  public :: EtUptakeModelType

  !> @brief Abstract base for root-water uptake (transpiration) reduction models
  !<
  type, abstract :: EtUptakeModelType
    integer(I4B) :: id !< the unique uptake model id
  contains
    procedure(stress_if), deferred :: stress
    procedure(fc_terms_if), deferred :: fc_terms
    procedure :: destroy
  end type EtUptakeModelType

  abstract interface
    ! dimensionless stress response in [0, 1] for pressure head psi
    function stress_if(this, psi) result(alpha)
      import EtUptakeModelType, DP
      class(EtUptakeModelType), intent(in) :: this
      real(DP), intent(in) :: psi !< soil pressure head (h - z)
      real(DP) :: alpha !< stress response, 0 (no uptake) to 1 (potential)
    end function

    ! linearized terms so that q = hcof*h - rhs = -c*alpha(h - z), where c is
    ! the potential volumetric transpiration rate (>= 0) and the returned flow
    ! q is the uptake as a flow to the aquifer (negative for extraction)
    subroutine fc_terms_if(this, h, z, c, hcof, rhs)
      import EtUptakeModelType, DP
      class(EtUptakeModelType), intent(in) :: this
      real(DP), intent(in) :: h !< hydraulic head in the cell
      real(DP), intent(in) :: z !< cell-center elevation
      real(DP), intent(in) :: c !< potential volumetric transpiration rate (>= 0)
      real(DP), intent(out) :: hcof !< coefficient on head
      real(DP), intent(out) :: rhs !< right-hand side term
    end subroutine
  end interface

contains

  !> @brief Clean up (base implementation is a no-op)
  !<
  subroutine destroy(this)
    class(EtUptakeModelType), intent(inout) :: this
  end subroutine destroy

end module EtUptakeModelModule
