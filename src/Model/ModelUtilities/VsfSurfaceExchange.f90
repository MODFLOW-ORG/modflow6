!> @brief Land-surface exchange kernel for variably saturated (Richards) flow
!!
!! Evaluates the exchange between a surface pond and the top cell of a Richards
!! cell, given the water that is supplied to the pond (precipitation).  The
!! exchange is driven by the total head difference between the pond and the
!! cell, through a conductance scaled by the relative permeability, so the same
!! expression describes infiltration (positive, into the cell) and
!! exfiltration (negative, out of the cell).
!!
!! Processes and sequencing for one surface cell (flows out of the package to
!! the atmosphere are evaporation; the pond exchanges with the top cell):
!!
!!                 precipitation qp            evaporation to the atmosphere
!!                       |                   (1) E_pond          (4) E_soil
!!     (2) net supply    |                        ^                   ^
!!     qp - E_pond       v                        |                   |
!!            +----------------------+            |                   |
!!            |         POND         |------------+                   |
!!            |    0 <= d <= dmax    |                                |
!!            |   (storage: d_dot)   |----> runoff qro                |
!!            +----------------------+        (to MVR / out of model) |
!!                       |                                            |
!!            (3) qinf   |   infiltration (down) / exfiltration (up)  |
!!                       v                                            |
!!            +----------------------+                                |
!!            |   TOP RICHARDS CELL  |--------------------------------+
!!            |       (aquifer)      |
!!            +----------------------+
!!
!! Only qinf crosses between the pond and the cell; pond storage stays in the
!! pond and runoff leaves the model.  The processes are applied in this order:
!!   (1) Pond evaporation: E_pond = min(E_pot, water available in the pond).
!!       It is head independent, so it is constant during the iterations.
!!   (2) Net supply: the precipitation net of pond evaporation, qp - E_pond,
!!       drives the pond balance below.
!!   (3) Pond balance: the implicit, head-dependent exchange qinf
!!       (infiltration into, or exfiltration out of, the cell), the pond depth
!!       d, and the runoff qro, solved together so the pond balance closes.
!!   (4) Soil evaporation: the residual demand R = E_pot - E_pond is applied to
!!       the soil as bare-soil evaporation E_soil = krel * R, extracted from the
!!       cell; total evaporation E_pond + E_soil never exceeds E_pot.
!!
!! With the conductance frozen at the current iterate, the implicit pond
!! balance is linear in the pond depth and has a closed-form solution that is
!! the minimum of three terms: the water that is available (supply limited),
!! the exchange with a pond that stores water (pond storing), and the exchange
!! with a full pond (pond full).  The minimum is smoothed so that its
!! derivative with respect to head is continuous, and the same smoothed
!! exchange is used to update the pond depth and the runoff, so that the water
!! balance of the pond closes exactly.
!!
!! The kernel is independent of any package and is expressed in terms of the
!! linearization hcof and rhs, with exchange = hcof * h - rhs, and of the
!! additional Newton term that results from the head dependence of the
!! relative permeability.
!<
module VsfSurfaceExchangeModule
  use KindModule, only: DP
  use ConstantsModule, only: DZERO, DHALF, DONE, DQUARTER, DEM30

  implicit none
  private

  public :: SurfaceExchangeResultType
  public :: surface_exchange
  public :: smooth_min

  !> @brief Surface exchange, pond state, and linearization for one surface cell
  !<
  type :: SurfaceExchangeResultType
    real(DP) :: qinf = DZERO !< pond-cell exchange, positive into the cell
    real(DP) :: pond = DZERO !< pond depth at the end of the step
    real(DP) :: qro = DZERO !< runoff (rate)
    real(DP) :: qevap_pond = DZERO !< pond (open-water) evaporation (>= 0, to atmosphere)
    real(DP) :: qevap_soil = DZERO !< bare-soil evaporation (>= 0, extracted from the cell)
    real(DP) :: hcof = DZERO !< coefficient of head, includes pond exchange and soil evaporation
    real(DP) :: rhs = DZERO !< right-hand side, net flow into the cell = hcof * h - rhs
    real(DP) :: deriv = DZERO !< Newton term from the head dependence of krel
  end type SurfaceExchangeResultType

contains

  !> @brief Smooth minimum of two values
  !!
  !! Equal to the minimum when the values differ by at least width and
  !! blended by a quadratic otherwise, so that the result and its first
  !! derivatives are continuous.  The result never exceeds the minimum.  The
  !! weights are the derivatives with respect to a and b.
  !<
  pure subroutine smooth_min(a, b, width, smin, wa, wb)
    real(DP), intent(in) :: a !< first value
    real(DP), intent(in) :: b !< second value
    real(DP), intent(in) :: width !< smoothing width; no smoothing if not positive
    real(DP), intent(out) :: smin !< smoothed minimum
    real(DP), intent(out) :: wa !< derivative of smin with respect to a
    real(DP), intent(out) :: wb !< derivative of smin with respect to b
    ! local
    real(DP) :: x

    x = a - b
    if (width <= DZERO .or. abs(x) >= width) then
      if (x <= DZERO) then
        smin = a
        wa = DONE
        wb = DZERO
      else
        smin = b
        wa = DZERO
        wb = DONE
      end if
    else
      smin = DHALF * (a + b) - (x * x / (4.0_DP * width) + DQUARTER * width)
      wa = DHALF - x / (2.0_DP * width)
      wb = DHALF + x / (2.0_DP * width)
    end if
  end subroutine smooth_min

  !> @brief Calculate the surface exchange, pond state, runoff, and evaporation
  !!
  !! The surface conductance is cs = krel * cbase.  The exchange is
  !! positive when water moves from the pond into the cell.
  !!
  !! Evaporation is removed first from the pond (open-water evaporation, limited
  !! by the water available in the pond), and the remaining demand is applied to
  !! the soil as bare-soil evaporation, reduced by the relative permeability of
  !! the cell.  Pond evaporation is head independent and enters as a constant
  !! reduction of the precipitation supply; soil evaporation is a head-dependent
  !! cell sink and is added to the linearization.
  !<
  pure subroutine surface_exchange(h, zs, dmax, area, cbase, krel, dkrdh, &
                                   qp, epot, d0, delt, hwidth, steady, res)
    real(DP), intent(in) :: h !< head of the cell at the current iterate
    real(DP), intent(in) :: zs !< land-surface elevation
    real(DP), intent(in) :: dmax !< maximum pond depth
    real(DP), intent(in) :: area !< horizontal area of the surface cell
    real(DP), intent(in) :: cbase !< conductance at unit relative permeability
    real(DP), intent(in) :: krel !< relative permeability of the cell
    real(DP), intent(in) :: dkrdh !< head derivative of the relative permeability
    real(DP), intent(in) :: qp !< supply rate to the pond (volume per time)
    real(DP), intent(in) :: epot !< potential evaporation rate (volume per time, >= 0)
    real(DP), intent(in) :: d0 !< pond depth at the start of the time step
    real(DP), intent(in) :: delt !< time step length
    real(DP), intent(in) :: hwidth !< head width for smoothing the minimum
    logical, intent(in) :: steady !< .true. if the pond has no storage
    type(SurfaceExchangeResultType), intent(out) :: res !< results
    ! local
    real(DP) :: cs, dcsdh, ceff, dceffdh, beta, hs0, supply, width
    real(DP) :: q_storing, q_full, q_supply_storing, dunc, hfull
    real(DP) :: wa_supply, wa_storing, w_supply_storing
    real(DP) :: w_supply, w_storing, w_full
    real(DP) :: supply_evap, epond, qp_eff, esoil_pot, esoil

    cs = krel * cbase ! relative permeability times the conductance
    dcsdh = dkrdh * cbase ! its derivative
    hfull = zs + dmax ! head for full pond
    width = cs * hwidth ! smoothing on fluxes, so multiply with conductance

    ! pond (open-water) evaporation, limited by the water available in the pond
    ! this step; it is head independent and so constant during the iterations
    if (steady) then
      supply_evap = qp
    else
      supply_evap = qp + d0 * area / delt
    end if
    epond = min(epot, max(supply_evap, DZERO))
    epond = max(epond, DZERO)
    qp_eff = qp - epond ! precipitation supply net of pond evaporation

    ! residual demand exposed to the drying soil, reduced by krel
    esoil_pot = max(epot - epond, DZERO)
    esoil = krel * esoil_pot

    q_full = cs * (hfull - h) ! rate when ponding at max

    if (steady) then
      ! the pond has no storage: all water that does not exchange is runoff
      call smooth_min(qp_eff, q_full, width, res%qinf, w_supply, w_full)
      res%hcof = -w_full * cs
      res%deriv = w_full * dcsdh * (hfull - h)
      if (cs > DEM30) then
        res%pond = min(max((qp_eff - cs * (zs - h)) / cs, DZERO), dmax)
      else if (qp_eff > DZERO) then
        res%pond = dmax
      else
        res%pond = DZERO
      end if
      res%qro = qp_eff - res%qinf
    else
      beta = delt * cs / area
      ceff = cs / (DONE + beta)
      dceffdh = dcsdh / (DONE + beta)**2
      hs0 = zs + d0 + delt * qp_eff / area
      supply = qp_eff + d0 * area / delt
      q_storing = ceff * (hs0 - h)
      call smooth_min(supply, q_storing, width, q_supply_storing, &
                      wa_supply, wa_storing)
      call smooth_min(q_supply_storing, q_full, width, res%qinf, &
                      w_supply_storing, w_full)
      w_supply = w_supply_storing * wa_supply
      w_storing = w_supply_storing * wa_storing
      res%hcof = -(w_storing * ceff + w_full * cs)
      res%deriv = w_storing * dceffdh * (hs0 - h) + w_full * dcsdh * (hfull - h)
      dunc = d0 + delt * (qp_eff - res%qinf) / area
      res%pond = min(max(dunc, DZERO), dmax)
      res%qro = qp_eff - res%qinf - (res%pond - d0) * area / delt
    end if

    ! guard against round-off, which can leave a tiny negative runoff
    res%qro = max(res%qro, DZERO)

    ! soil evaporation is a head-dependent cell sink through krel; the flow into
    ! the cell is q_soil = -esoil = -krel * esoil_pot.  It contributes a Newton
    ! term and is folded into the net flow into the cell (qinf - esoil).
    res%qevap_pond = epond
    res%qevap_soil = esoil
    res%deriv = res%deriv - esoil_pot * dkrdh
    res%rhs = res%hcof * h - (res%qinf - esoil)
  end subroutine surface_exchange

end module VsfSurfaceExchangeModule
