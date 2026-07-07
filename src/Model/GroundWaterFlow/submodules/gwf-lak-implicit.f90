submodule(LakModule) LakModuleImplicit
contains

  !> @brief Sum the non-seepage lake budget terms at a given stage
  !!
  !! Returns B(stage), the sum of every lake-budget term except the groundwater
  !! seepage: rainfall + runoff + specified and external inflow + outlet inflow
  !! and outflow - evaporation - withdrawal + transient storage. The stage-driven
  !! losses are ramped to zero as the stage approaches the lake bottom (over the
  !! surfdep interval) so a drying lake cannot remove more water than it holds.
  !! lak_fc_implicit calls this routine at two stages to finite-difference the
  !! lake-row diagonal and right-hand side.
  !<
  module procedure lak_budget_nogwf
  ! -- modules
  use TdisModule, only: delt
  ! -- local
  real(DP) :: ra, ro, qin, ex, surfin, sa, v0, v1, qout, sf
  !
  call this%lak_calculate_rainfall(n, stage, ra)
  call this%lak_calculate_runoff(n, ro)
  call this%lak_calculate_inflow(n, qin)
  call this%lak_calculate_external(n, ex)
  call this%lak_calculate_outlet_inflow(n, surfin)
  call this%lak_outlet_outflow_rate(n, stage, qout)
  call this%lak_calculate_sarea(n, stage, sa)
  !
  ! -- smooth flow-limiting: ramp the stage-driven losses (evaporation,
  !    withdrawal, outlet outflow) to zero as the lake approaches its bottom
  !    over the surfdep interval, so a drying lake cannot extract more water
  !    than it holds. With surfdep = 0 (or a lake well above its bottom) the
  !    factor is 1 and the losses are unlimited (matching the legacy result
  !    whenever the lake is not water-limited).
  sf = DONE
  if (this%surfdep > DZERO) then
    sf = sQuadraticSaturation(this%lakebot(n) + this%surfdep, &
                              this%lakebot(n), stage)
  end if
  !
  ! -- evaporation and withdrawal are losses (positive demand); qout is the
  !    (negative) outlet outflow; surfin is routed inflow from upstream outlets
  !    (lagged at the previous iterate's stage via simoutrate)
  b = ra + ro + qin + ex + surfin &
      + (qout - this%evaporation(n) * sa - this%withdrawal(n)) * sf
  !
  ! -- transient storage change (v0 - v1)/delt
  if (this%gwfiss /= 1) then
    call this%lak_calculate_vol(n, this%xoldpak(n), v0)
    call this%lak_calculate_vol(n, stage, v1)
    b = b + (v0 - v1) / delt
  end if
  end procedure lak_budget_nogwf

  !> @brief Assemble the lake equations into the groundwater flow matrix
  !!
  !! Implicit formulation. For each active lake the water-balance equation is
  !! assembled directly into the matrix:
  !!   - the lakebed seepage to each connected cell is a conductance coupling
  !!     between the lake-stage row and the cell row. In the seepage, the lake
  !!     stage and the connected-cell head are each kept at or above the lake
  !!     bottom (the same wet/dry cutoff as the default formulation), so a wet
  !!     connection couples both stage and head, while a perched connection keeps
  !!     the stage coupling (a non-zero lake-stage diagonal) but drops the head
  !!     coupling, and
  !!   - the remaining (non-seepage) budget terms -- rainfall, evaporation,
  !!     withdrawal, runoff, specified and external inflow, outlet inflow and
  !!     outflow, and transient storage -- are linearized in stage with a finite
  !!     difference (lak_budget_nogwf) and added to the diagonal and right-hand
  !!     side. This captures the stage-dependent surface area and the storage
  !!     term.
  !! The implicit lake-aquifer coupling is asymmetric (a perched connection
  !! couples the lake to the cell but not the reverse), so the coefficient matrix
  !! is asymmetric and the BICGSTAB linear acceleration is required.
  !! Lake-to-lake outlet routing is lagged one outer iteration through simoutrate.
  !! A constant-stage lake holds its row at the specified stage; an inactive lake
  !! is disconnected from the aquifer. A lake row whose diagonal is still (near)
  !! zero (no connections and no storage) is kept solvable with a small diagonal
  !! term that holds the stage at its current value and becomes zero at
  !! convergence. A lake flagged for the substitution fallback (lak_set_fallback)
  !! is instead solved by the default substitution solver and assembled like a
  !! constant-stage lake at that solved stage.
  !<
  module procedure lak_fc_implicit
  ! -- local
  integer(I4B) :: n, j, ipos, iloc, igwfnode
  logical(LGP) :: lfallback
  real(DP) :: hlak, head, flow, gwfhcof, gwfrhs
  real(DP) :: b0, b1, dbds, avail, sout, csum, deps, adiag
  real(DP) :: dqds, dqdh
  !
  ! -- update simoutrate at the current stages so the routed outlet inflow to
  !    downstream lakes (lak_calculate_outlet_inflow, read in lak_budget_nogwf)
  !    reflects the latest upstream stages (lagged one outer iteration). The
  !    outlet outflow sink itself is linearized per lake in lak_budget_nogwf.
  if (this%noutlets > 0) then
    do n = 1, this%nlakes
      if (this%iboundpak(n) < 1) cycle
      avail = DEP20
      call this%lak_calculate_outlet_outflow(n, this%xnewpak(n), avail, sout)
    end do
    !
    ! -- provide the outlet outflow to the mover (matches the mover-provider
    !    accumulation done at the end of the legacy lak_solve)
    if (this%imover == 1) then
      do n = 1, this%noutlets
        call this%pakmvrobj%accumulate_qformvr(n, -this%simoutrate(n))
      end do
    end if
  end if
  !
  ! -- development option: force every active lake onto the fallback (used to
  !    test the fallback assembly against the legacy formulation)
  if (this%iforcefb /= 0) then
    do n = 1, this%nlakes
      if (this%iboundpak(n) > 0) this%ifallback(n) = 1
    end do
  end if
  !
  ! -- solve the stage of any lake assigned to the substitution fallback
  !    against the current groundwater heads. A fallback lake is then
  !    assembled below like a constant-stage lake at its solved stage, so a
  !    weakly connected or disconnected lake -- which is poorly conditioned as
  !    a matrix unknown -- is handled by the robust 1D substitution instead.
  lfallback = .false.
  do n = 1, this%nlakes
    if (this%ifallback(n) /= 0) then
      lfallback = .true.
      exit
    end if
  end do
  if (lfallback) then
    call this%lak_solve(only_fallback=.true.)
  end if
  !
  ipos = 0
  do n = 1, this%nlakes
    iloc = this%idxlocnode(n)
    hlak = this%xnewpak(n)
    if (this%ifallback(n) /= 0 .and. this%iboundpak(n) > 0) then
      !
      ! -- fallback lake: its stage was solved by substitution above. Hold the
      !    lake row at that stage (unit diagonal) and add the resulting
      !    lake-aquifer exchange (this%hcof/this%rhs, exactly as the legacy
      !    formulation does) to the connected gwf cells.
      call matrix_sln%add_value_pos(this%idxdiag(n), DONE)
      rhs(iloc) = hlak
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        ipos = ipos + 1
        igwfnode = this%cellid(j)
        if (this%ibound(igwfnode) < 1) cycle
        call matrix_sln%add_value_pos(this%idxsymdglo(ipos), this%hcof(j))
        rhs(igwfnode) = rhs(igwfnode) + this%rhs(j)
      end do
    else if (this%iboundpak(n) > 0) then
      !
      ! -- active lake: linearize the non-seepage budget B(stage) about the
      !    current stage with a finite difference. dB/dstage goes on the
      !    lake-row diagonal and (dB/dstage*stage - B) goes on the rhs.
      call this%lak_budget_nogwf(n, hlak, b0)
      call this%lak_budget_nogwf(n, hlak + this%delh, b1)
      dbds = (b1 - b0) / this%delh
      call matrix_sln%add_value_pos(this%idxdiag(n), dbds)
      rhs(iloc) = rhs(iloc) + dbds * hlak - b0
      adiag = dbds
      !
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        ipos = ipos + 1
        igwfnode = this%cellid(j)
        if (this%ibound(igwfnode) < 1) cycle
        head = this%xnew(igwfnode)
        !
        ! -- lakebed seepage and its derivatives. The lake stage and the
        !    connected-cell head are each kept at or above the lake bottom,
        !    the same wet/dry cutoff as the default formulation. A wet
        !    connection (dps = dph = 1) couples both stage and head; a perched
        !    connection (dph = 0) drops the head coupling but keeps the stage
        !    coupling (dqds), so the lake-stage diagonal stays non-zero and the
        !    perched leakage still responds to stage.
        call this%lak_calculate_conn_exchange_deriv(n, j, hlak, head, flow, &
                                                    dqds, dqdh)
        call matrix_sln%add_value_pos(this%idxdiag(n), dqds)
        call matrix_sln%add_value_pos(this%idxoffdglo(ipos), dqdh)
        call matrix_sln%add_value_pos(this%idxsymdglo(ipos), -dqdh)
        call matrix_sln%add_value_pos(this%idxsymoffdglo(ipos), -dqds)
        rhs(iloc) = rhs(iloc) + dqds * hlak + dqdh * head - flow
        rhs(igwfnode) = rhs(igwfnode) - dqdh * head - dqds * hlak + flow
        adiag = adiag + dqds
      end do
      !
      ! -- keep the row solvable when its diagonal is still (near) zero -- a
      !    lake with no wet connections and no storage. Add a small diagonal
      !    term and a matching rhs term so the row holds the stage at its
      !    current value; the contribution is exactly zero at convergence.
      csum = DZERO
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        csum = csum + this%satcond(j)
      end do
      deps = DEM9 * csum + DPREC
      if (abs(adiag) < deps) then
        call matrix_sln%add_value_pos(this%idxdiag(n), -deps)
        rhs(iloc) = rhs(iloc) - deps * hlak
      end if
    else if (this%iboundpak(n) < 0) then
      !
      ! -- constant-stage lake: hold the lake row at the specified stage and
      !    add the seepage to the gwf cells as a constant exchange evaluated at
      !    the fixed stage
      call matrix_sln%add_value_pos(this%idxdiag(n), DONE)
      rhs(iloc) = hlak
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        ipos = ipos + 1
        igwfnode = this%cellid(j)
        if (this%ibound(igwfnode) < 1) cycle
        head = this%xnew(igwfnode)
        call this%lak_calculate_conn_exchange(n, j, hlak, head, flow, &
                                              gwfhcof, gwfrhs)
        call matrix_sln%add_value_pos(this%idxsymdglo(ipos), gwfhcof)
        rhs(igwfnode) = rhs(igwfnode) + gwfrhs
      end do
    else
      !
      ! -- inactive lake (iboundpak == 0): the lake is absent and exchanges no
      !    water with the aquifer. Disconnect the lake row (its global ibound
      !    is 0, so the solver drops the equation) with a unit diagonal and add
      !    NO seepage to the gwf cells. Only advance ipos to keep the
      !    connection index aligned. rhs is left at zero so the large
      !    placeholder value used for an inactive stage is not added to the
      !    linear system.
      call matrix_sln%add_value_pos(this%idxdiag(n), DONE)
      rhs(iloc) = DZERO
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        ipos = ipos + 1
      end do
    end if
  end do
  end procedure lak_fc_implicit

  !> @brief Switch a stalled IMPLICIT lake to the substitution fallback
  !!
  !! Under the IMPLICIT option a weakly connected or disconnected lake (one with
  !! a small lakebed conductance) gives a poorly conditioned lake-stage equation
  !! that the coupled solver may be unable to advance. Such a lake is the
  !! convergence bottleneck: its stage keeps changing by much more than its
  !! connected aquifer heads, which have settled. A well-connected lake, by
  !! contrast, changes its stage in step with its aquifer. When the stage change
  !! exceeds the connected-head change by a large factor for several outer
  !! iterations (see the inline criterion), the lake is flagged (ifallback), and
  !! lak_fc_implicit then solves and assembles it with the robust substitution
  !! solver instead. A lake that changes in step with its aquifer is never
  !! switched. The flag and the per-lake tracking are cleared at the start of
  !! each time step, so a lake switched in one time step is reconsidered for the
  !! implicit formulation in the next, where transient storage or different
  !! stresses may make its stage equation better conditioned.
  !<
  module procedure lak_set_fallback
  ! -- local
  integer(I4B) :: n, j, igwfnode
  real(DP) :: dstage, dhead
  ! -- a lake is treated as the convergence bottleneck when its stage change
  !    exceeds headratio times the largest change in its connected aquifer
  !    heads; more than nstuckmax such consecutive outer iterations switch it to
  !    the fallback. A well-connected lake converges in lockstep with its aquifer
  !    (ratio near one) and is never switched; a weakly connected or disconnected
  !    lake keeps thrashing while its connected heads settle (ratio grows without
  !    bound), so it is switched early -- before the heads commit to the stalled
  !    configuration -- which the fallback then resolves quickly.
  real(DP), parameter :: headratio = DTEN
  integer(I4B), parameter :: nstuckmax = 5
  !
  if (this%iimplicit == 0) return
  !
  ! -- at the start of a time step, snapshot the connected heads and reset the
  !    bottleneck counts. Only snapshot active cells; an inactive cell holds a
  !    placeholder head that would otherwise distort the head-change comparison
  !    if the cell later becomes active.
  if (kiter == 1) then
    do n = 1, this%nlakes
      this%ifallback(n) = 0
      this%nstuck(n) = 0
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        igwfnode = this%cellid(j)
        if (this%ibound(igwfnode) >= 1) then
          this%holdconn(j) = this%xnew(igwfnode)
        end if
      end do
    end do
    return
  end if
  !
  ! -- only evaluate while the model has not converged
  if (icnvgmod /= 0) return
  do n = 1, this%nlakes
    if (this%iboundpak(n) < 1) cycle
    !
    ! -- stage change over this outer iteration and the largest change in the
    !    connected-cell heads, then advance the snapshot
    dstage = abs(this%xnewpak(n) - this%s0(n))
    dhead = DZERO
    do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
      igwfnode = this%cellid(j)
      if (this%ibound(igwfnode) >= 1) then
        dhead = max(dhead, abs(this%xnew(igwfnode) - this%holdconn(j)))
        this%holdconn(j) = this%xnew(igwfnode)
      end if
    end do
    if (this%ifallback(n) /= 0) cycle
    !
    ! -- count outer iterations in which the unconverged stage is changing much
    !    faster than its (settling) aquifer; switch to the fallback after enough
    if (dstage > this%dmaxchg .and. dstage > headratio * dhead) then
      this%nstuck(n) = this%nstuck(n) + 1
    else
      this%nstuck(n) = 0
    end if
    if (this%nstuck(n) > nstuckmax) then
      this%ifallback(n) = 1
    end if
  end do
  end procedure lak_set_fallback

  !> @brief Warn when a fallback lake still prevents IMPLICIT convergence
  !!
  !! Called on the last outer iteration of a solution that did not converge. By
  !! this point any stalled lake has already been switched to the substitution
  !! fallback (lak_set_fallback). If at least one lake is on the fallback and the
  !! solution still did not converge, the implicit formulation needed more outer
  !! iterations than were allowed; a one-time warning reports this and points the
  !! user to the remedies, including the default formulation, which converges
  !! this kind of lake in far fewer outer iterations. No warning is issued when
  !! no lake was switched, because then the lake is not the cause of the failure.
  !<
  module procedure lak_check_disconnected
  ! -- local
  character(len=LINELENGTH) :: warnmsg
  integer(I4B) :: n
  logical(LGP) :: hasfallback
  !
  ! -- only the implicit formulation couples the stage to the solver
  if (this%iimplicit == 0) return
  !
  ! -- only warn if a lake was switched to the fallback (otherwise the lake is
  !    not the cause of the non-convergence)
  hasfallback = .false.
  do n = 1, this%nlakes
    if (this%ifallback(n) /= 0) then
      hasfallback = .true.
      exit
    end if
  end do
  if (.not. hasfallback) return
  !
  ! -- warn with remedies; deduplicated so it is issued at most once
  write (warnmsg, '(a)') &
    "IMPLICIT lake option in package '"//trim(this%packName)//"' did not "// &
    "converge even after switching a stalled lake to the substitution "// &
    "solver. The default solver needs fewer outer iterations for such a "// &
    "lake; consider omitting IMPLICIT, adding an outlet, increasing the "// &
    "leakance, or running transient."
  call store_warning(warnmsg, substring=warnmsg(:40))
  end procedure lak_check_disconnected

end submodule LakModuleImplicit
