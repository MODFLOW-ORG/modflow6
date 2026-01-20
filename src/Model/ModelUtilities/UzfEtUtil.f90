module UzfETUtilModule
  use KindModule, only: DP
  use ConstantsModule, only: DZERO, DONE, DEM3, DEM4
  use SmoothingModule, only: sCubic

  implicit none
  private
  public :: etfunc_lin
  public :: etfunc_nlin
  public :: calc_lin_scaling_fac

contains

  !> @brief Calculate gwf ET using linear decay ET function from mf-2005
  !<
  function etfunc_lin(efflndsrf, extdp, resid_pet, deriv_et, trhs, thcof, &
                      hgwf, celtop, celbot)
    ! return
    real(DP) :: etfunc_lin
    ! dummy
    real(DP), intent(in) :: efflndsrf !< effective land surface elevation after subtracting off 0.5*surfdep
    real(DP), intent(in) :: extdp !< extinction depth
    real(DP), intent(in) :: resid_pet !< residual pET remaining after subtracting simulated ET in the unsaturated zone
    real(DP), intent(inout) :: deriv_et !< derivative of gw ET for Newton addition to equations in _fn()
    real(DP), intent(inout) :: trhs !< total uzf rhs contribution to GWF model
    real(DP), intent(inout) :: thcof !< total uzf hcof contribution to GWF model
    real(DP), intent(in) :: hgwf !< calculated groundwater head
    real(DP), intent(in) :: celtop !< elevation of the top of the cell
    real(DP), intent(in) :: celbot !< elevation of the bottom of the cell
    ! local
    real(DP) :: etgw
    real(DP) :: range
    real(DP) :: depth, scale, thick, lin_scaling_fac
    !
    ! initialize
    etgw = DZERO
    !
    ! if water table between ET surface and extinction depth
    ! (extdp starts at the bottom of the effective land surface elevation,
    !  where the effective land surface elevation accounts for surfdep)
    if (hgwf > (efflndsrf - extdp) .and. hgwf < efflndsrf) THEN
      lin_scaling_fac = calc_lin_scaling_fac(hgwf, efflndsrf, extdp)
      etgw = resid_pet * lin_scaling_fac
      !
      trhs = resid_pet - resid_pet * efflndsrf / extdp
      thcof = -resid_pet / extdp
      etgw = trhs - (thcof * hgwf)
      !
      ! water table above land surface
    else if (hgwf >= efflndsrf) then
      trhs = resid_pet
      etgw = resid_pet
      !
      ! water table below extinction depth
    else
      etfunc_lin = DZERO
      return
    end if
    !
    ! calculate rate
    depth = hgwf - (efflndsrf - extdp)
    thick = celtop - celbot
    if (depth > thick) depth = thick
    if (depth < DZERO) depth = DZERO
    range = DEM4 * extdp
    call sCubic(depth, range, deriv_et, scale)
    trhs = scale * trhs
    thcof = scale * thcof
    etgw = trhs - (thcof * hgwf)
    deriv_et = -deriv_et * etgw
    etfunc_lin = etgw
    !
  end function etfunc_lin

  !> @brief Calculate gwf ET using a square decay ET function with smoothing
  !! at the specified extinction depth
  !<
  function etfunc_nlin(efflndsrf, extdp, resid_pet, deriv_et, trhs, thcof, hgwf)
    ! -- return
    real(DP) :: etfunc_nlin
    ! -- dummy
    real(DP), intent(in) :: efflndsrf !< effective land surface elevation after subtracting off 0.5*surfdep
    real(DP), intent(in) :: extdp !< extinction depth
    real(DP), intent(in) :: resid_pet !< residual pET remaining after subtracting simulated ET in the unsaturated zone
    real(DP), intent(inout) :: deriv_et !< derivative of gw ET for Newton addition to equations in _fn()
    real(DP), intent(inout) :: trhs !< total uzf rhs contribution to GWF model
    real(DP), intent(inout) :: thcof !< total uzf hcof contribution to GWF model
    real(DP), intent(in) :: hgwf !< calculated groundwater head
    ! -- local
    real(DP) :: etgw
    real(DP) :: range
    real(DP) :: depth, scale
    !
    depth = hgwf - (efflndsrf - extdp)
    if (depth < DZERO) depth = DZERO
    etgw = resid_pet
    range = DEM3 * extdp
    call sCubic(depth, range, deriv_et, scale)
    etgw = etgw * scale
    trhs = etgw
    deriv_et = -deriv_et * etgw
    etfunc_nlin = etgw
  end function etfunc_nlin

  !> @brief Calculate the linear scaling factor
  !<
  pure function calc_lin_scaling_fac(hgwf, lndsrf, extdp) result(sclfac)
    ! dummy
    real(DP), intent(in) :: hgwf !< groundwater head
    real(DP), intent(in) :: lndsrf !< effective land surface (after applying surfdep to land surface)
    real(DP), intent(in) :: extdp !< extinction depth
    ! return
    real(DP) :: sclfac
    !
    sclfac = (hgwf - (lndsrf - extdp)) / extdp
    ! the calculated scaling factor cannot exceed 1.0
    if (sclfac > DONE) sclfac = DONE
  end function calc_lin_scaling_fac

end module UzfETUtilModule
