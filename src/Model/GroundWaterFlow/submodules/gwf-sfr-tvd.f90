submodule(SfrModule) SfrModuleTvd
contains

  !> @brief Kinematic wave routing with TVD flux limiter
  !!
  !! Explicit TVD (Total Variation Diminishing) update for the kinematic
  !! wave continuity equation.  The van Leer flux limiter provides
  !! second-order accuracy in smooth flow regions and reverts to first-order
  !! upwind near flow discontinuities, guaranteeing exact mass conservation.
  !!
  !! The method is activated when ATS_COURANT is specified.  When the
  !! estimated Courant number exceeds 1, the call is delegated to
  !! sfr_calc_transient, whose Cr > 1 branch is already mass-conservative.
  !! Dry reaches (d_old = 0) remain in the TVD path: q_tvd = Q(0) = 0 and
  !! the direct volume balance fills the reach exactly, giving exact closure.
  !!
  !! The old-time depth is taken from stageold(n) - strtop(n), the reach
  !! depth frozen at the start of the time step by sfr_ad.  Using stageold
  !! rather than the current iterate depth(n) makes the TVD update idempotent
  !! across outer Picard iterations, so the outer solver converges in two
  !! passes even for multi-reach networks.
  !!
  !! The pre-computed itvd_upstream(n) array holds the index of the single
  !! upstream reach feeding reach n.  Headwaters and confluences (index = 0)
  !! use first-order upwind for the outflow flux.
  !<
  module procedure sfr_calc_tvd
  use TdisModule, only: delt
  ! -- local
  integer(I4B) :: igwfconn
  integer(I4B) :: i
  integer(I4B) :: j
  integer(I4B) :: iup
  real(DP) :: celerity
  real(DP) :: courant
  real(DP) :: dq
  real(DP) :: d_old
  real(DP) :: d2
  real(DP) :: a_old
  real(DP) :: a_new
  real(DP) :: a2
  real(DP) :: v_new
  real(DP) :: q_out_old
  real(DP) :: q_in_old
  real(DP) :: q_in2_old
  real(DP) :: dq_loc
  real(DP) :: dq_up
  real(DP) :: r
  real(DP) :: phi
  real(DP) :: q_tvd
  real(DP) :: q_in_new
  real(DP) :: q_lat
  real(DP) :: q_avail
  real(DP) :: d1_old
  real(DP) :: tw
  real(DP) :: res
  real(DP) :: delh
  !
  ! -- old depth from start-of-timestep state; frozen by sfr_ad so the TVD
  !    update is idempotent across outer Picard iterations
  d_old = max(this%stageold(n) - this%strtop(n), DZERO)
  a_old = this%calc_area_wet(n, d_old)
  !
  ! -- Manning flow at old depth gives the first-order upwind flux
  call this%sfr_calc_qman(n, d_old, q_out_old)
  q_in_old = this%usinflowold(n)
  !
  ! -- celerity and Courant number from a flow perturbation at old depth
  dq = this%deps
  celerity = DZERO
  if (d_old > DZERO) then
    call this%sfr_calc_reach_depth(n, q_out_old + dq, d2)
    a2 = this%calc_area_wet(n, d2)
    if (a2 > a_old) then
      celerity = dq / (a2 - a_old)
    end if
  end if
  courant = celerity * delt / this%length(n)
  !
  ! -- Cr > 1: delegate to upwind scheme (already mass-conservative)
  if (courant > DONE) then
    call this%sfr_calc_transient(n, d1, hgwf, qu, qi, qfrommvr, &
                                 qr, qe, qro, qgwf, qd)
    return
  end if
  !
  ! -- TVD outflow flux: first-order upwind + van Leer limiter correction
  q_tvd = q_out_old
  iup = this%itvd_upstream(n)
  if (iup > 0) then
    q_in2_old = this%usinflowold(iup)
    dq_loc = q_out_old - q_in_old
    dq_up = q_in_old - q_in2_old
    if (abs(dq_loc) > DEM30) then
      r = dq_up / dq_loc
      ! -- van Leer limiter: phi(r) = (r + |r|) / (1 + |r|)
      phi = (r + abs(r)) / (DONE + abs(r))
      ! -- Lax-Wendroff anti-diffusion correction bounded by limiter
      q_tvd = q_out_old + DHALF * (DONE - courant) * phi * dq_loc
    end if
  end if
  !
  ! -- total new-time volumetric inflows
  q_in_new = qu + qi + qfrommvr
  q_lat = qr + qro - qe
  !
  ! -- GWF exchange (initialized to zero for non-connected reaches)
  igwfconn = this%sfr_gwf_conn(n)
  qgwf = DZERO
  !
  ! -- initial depth estimate from actual stored depth
  d1 = d_old
  d1_old = d1
  !
  ! -- Picard iteration: refine qgwf and depth to convergence
  picard: do i = 1, this%maxsfrpicard
    if (igwfconn == 1) then
      call this%sfr_calc_qgwf(n, d1, hgwf, qgwf)
      qgwf = -qgwf
      q_avail = q_in_new + q_lat
      if (qgwf > q_avail) qgwf = q_avail
    end if
    !
    ! -- conservative volume balance (exact mass conservation)
    v_new = a_old * this%length(n) + delt * (q_in_new + q_lat - qgwf - q_tvd)
    if (v_new < DZERO) then
      ! -- reach goes dry: reduce outflow to available water
      q_tvd = max(q_in_new + q_lat - qgwf + &
                  a_old * this%length(n) / delt, DZERO)
      v_new = DZERO
    end if
    a_new = v_new / this%length(n)
    !
    ! -- Newton solve for depth from area: calc_area_wet(n, d1) = a_new
    d1 = max(d1, DZERO)
    newton_depth: do j = 1, this%maxsfrit
      res = this%calc_area_wet(n, d1) - a_new
      tw = this%calc_top_width_wet(n, d1)
      if (tw > DEM30) then
        d1 = d1 - res / tw
      else
        d1 = DZERO
      end if
      if (d1 < DZERO) d1 = DZERO
      if (abs(res) < this%deps) exit newton_depth
    end do newton_depth
    !
    delh = abs(d1 - d1_old)
    if (i > 1 .and. delh < this%dmaxchg) exit picard
    d1_old = d1
  end do picard
  !
  ! -- set outflow and storage budget term (exact: closes budget by construction)
  qd = max(q_tvd, DZERO)
  this%storage(n) = (a_old - a_new) * this%length(n) / delt
  !
  ! -- accumulate Courant statistics for ATS time step and listing output
  if (courant > DZERO) then
    if (courant < this%crmin(n)) this%crmin(n) = courant
    if (courant > this%crmax(n)) this%crmax(n) = courant
    this%crsum(n) = this%crsum(n) + courant
    this%crcnt(n) = this%crcnt(n) + 1
  end if
  !
  end procedure sfr_calc_tvd

end submodule
