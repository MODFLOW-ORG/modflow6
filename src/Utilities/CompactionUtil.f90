!> @brief Relations for compacting porous materials
!!
!! Conversions between porosity and void ratio, the change in thickness and
!! porosity of a compacting layer, the skeletal storage factor, and the
!! Newton-Raphson storage terms.
!<
module CompactionUtilModule
  use KindModule, only: DP, LGP
  use ConstantsModule, only: DZERO, DONE

  implicit none
  private

  public :: porosity_to_void_ratio
  public :: void_ratio_to_porosity
  public :: compact_layer
  public :: adjust_eff_stress
  public :: storage_factor
  public :: storage_newton_term
  public :: wcomp_newton_term

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

  !> @brief Calculate the Newton-Raphson term for skeletal storage
  !!
  !! The storage flow is rho sn (es - pcs) + rho0 sn0 (pcs - es0), where
  !! es = gs + bot - hbar; pcs is zero for elastic storage. Returns the
  !! derivative with respect to head, with gs and rho held fixed, less the
  !! Picard term (-rho sn).
  !<
  pure function storage_newton_term(rho, rho0, sn, satderv, hbarderv, es, &
                                    es0, pcs, lagged) result(hcof)
    ! -- dummy variables
    real(DP), intent(in) :: rho !< current storage coefficient without saturation
    real(DP), intent(in) :: rho0 !< previous storage coefficient without saturation
    real(DP), intent(in) :: sn !< current saturation
    real(DP), intent(in) :: satderv !< saturation derivative
    real(DP), intent(in) :: hbarderv !< corrected head derivative
    real(DP), intent(in) :: es !< current effective stress
    real(DP), intent(in) :: es0 !< previous effective stress
    real(DP), intent(in) :: pcs !< preconsolidation stress
    logical(LGP), intent(in) :: lagged !< previous saturation is the current one
    ! -- return
    real(DP) :: hcof
    !
    hcof = rho * sn * (DONE - hbarderv) + rho * (es - pcs) * satderv
    if (lagged) then
      hcof = hcof + rho0 * (pcs - es0) * satderv
    end if
  end function storage_newton_term

  !> @brief Calculate the Newton-Raphson term for water compressibility
  !!
  !! The storage flow is wc0 sn0 h0 - wc sn h. Returns the derivative with
  !! respect to head less the Picard term (-wc sn).
  !<
  pure function wcomp_newton_term(wc, wc0, h, h0, satderv, lagged) &
    result(hcof)
    ! -- dummy variables
    real(DP), intent(in) :: wc !< current water compressibility coefficient
    real(DP), intent(in) :: wc0 !< previous water compressibility coefficient
    real(DP), intent(in) :: h !< current head
    real(DP), intent(in) :: h0 !< previous head
    real(DP), intent(in) :: satderv !< saturation derivative
    logical(LGP), intent(in) :: lagged !< previous saturation is the current one
    ! -- return
    real(DP) :: hcof
    !
    hcof = -wc * h * satderv
    if (lagged) then
      hcof = hcof + wc0 * h0 * satderv
    end if
  end function wcomp_newton_term

end module CompactionUtilModule
