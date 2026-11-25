!> @brief Relations for compacting porous materials
!!
!! Conversions between porosity and void ratio, the change in thickness and
!! porosity of a compacting layer, and the skeletal storage factor.
!<
module CompactionUtilModule
  use KindModule, only: DP
  use ConstantsModule, only: DZERO, DONE

  implicit none
  private

  public :: porosity_to_void_ratio
  public :: void_ratio_to_porosity
  public :: compact_layer
  public :: adjust_eff_stress
  public :: storage_factor

contains

  !> @brief Calculate the void ratio from the porosity
  !<
  pure function porosity_to_void_ratio(theta) result(void_ratio)
    ! -- dummy variables
    real(DP), intent(in) :: theta !< porosity
    ! -- return
    real(DP) :: void_ratio
    !
    void_ratio = theta / (DONE - theta)
  end function porosity_to_void_ratio

  !> @brief Calculate the porosity from the void ratio
  !<
  pure function void_ratio_to_porosity(void_ratio) result(theta)
    ! -- dummy variables
    real(DP), intent(in) :: void_ratio !< void ratio
    ! -- return
    real(DP) :: theta
    !
    theta = void_ratio / (DONE + void_ratio)
  end function void_ratio_to_porosity

  !> @brief Update the thickness and porosity of a layer after compaction
  !!
  !! The solid thickness, thick / (1 + void ratio), is unchanged.
  !<
  pure subroutine compact_layer(comp, thick, theta)
    ! -- dummy variables
    real(DP), intent(in) :: comp !< compaction
    real(DP), intent(inout) :: thick !< initial and current thickness
    real(DP), intent(inout) :: theta !< initial and current porosity
    ! -- local variables
    real(DP) :: strain
    real(DP) :: void_ratio
    !
    ! -- initialize variables
    strain = DZERO
    void_ratio = porosity_to_void_ratio(theta)
    !
    ! -- calculate strain
    if (thick > DZERO) strain = -comp / thick
    !
    ! -- update void ratio, theta, and thickness
    void_ratio = void_ratio + strain * (DONE + void_ratio)
    theta = void_ratio_to_porosity(void_ratio)
    thick = thick - comp
  end subroutine compact_layer

  !> @brief Calculate the effective stress at elevation z from its value at z0
  !<
  pure function adjust_eff_stress(es0, z0, z, sgs) result(es)
    ! -- dummy variables
    real(DP), intent(in) :: es0 !< effective stress at elevation z0
    real(DP), intent(in) :: z0 !< elevation of es0
    real(DP), intent(in) :: z !< elevation to calculate the effective stress at
    real(DP), intent(in) :: sgs !< saturated specific gravity
    ! -- return
    real(DP) :: es
    !
    es = es0 - (z - z0) * (sgs - DONE)
  end function adjust_eff_stress

  !> @brief Calculate the skeletal storage factor, 1 / ((1 + void ratio) es)
  !!
  !! Returns zero if the effective stress is zero.
  !<
  pure function storage_factor(es, theta) result(fact)
    ! -- dummy variables
    real(DP), intent(in) :: es !< effective stress
    real(DP), intent(in) :: theta !< porosity
    ! -- return
    real(DP) :: fact
    ! -- local variables
    real(DP) :: denom
    !
    fact = DZERO
    denom = es * (DONE + porosity_to_void_ratio(theta))
    if (denom /= DZERO) then
      fact = DONE / denom
    end if
  end function storage_factor

end module CompactionUtilModule
