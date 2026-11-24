!> @brief CSUB package convergence check and budget
!<
submodule(GwfCsubModule) GwfCsubModuleBudget

  character(len=LENBUDTXT), dimension(4) :: budtxt = & !< text labels for budget terms
                                            [' CSUB-CGELASTIC', &
                                             '   CSUB-ELASTIC', &
                                             ' CSUB-INELASTIC', &
                                             ' CSUB-WATERCOMP']

contains

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module procedure csub_cc
  ! -- modules
  use TdisModule, only: totim, kstp, kper, delt
  ! local variables
  character(len=LENPAKLOC) :: cloc
  integer(I4B) :: icheck
  integer(I4B) :: ipakfail
  integer(I4B) :: ib
  integer(I4B) :: node
  integer(I4B) :: idelay
  integer(I4B) :: locdhmax
  integer(I4B) :: locrmax
  integer(I4B) :: ifirst
  real(DP) :: dhmax
  real(DP) :: rmax
  real(DP) :: dh
  real(DP) :: area
  real(DP) :: hcell
  real(DP) :: hcellold
  real(DP) :: snnew
  real(DP) :: snold
  real(DP) :: stoe
  real(DP) :: stoi
  real(DP) :: dwc
  real(DP) :: tled
  real(DP) :: hcof
  real(DP) :: rhs
  real(DP) :: v1
  real(DP) :: v2
  real(DP) :: df
  !
  ! -- initialize local variables
  icheck = this%iconvchk
  ipakfail = 0
  locdhmax = 0
  locrmax = 0
  ifirst = 1
  dhmax = DZERO
  rmax = DZERO
  !
  ! -- additional checks to see if convergence needs to be checked
  ! -- no convergence check for steady-state stress periods
  if (this%gwfiss /= 0) then
    icheck = 0
  else
    if (icnvgmod == 0) then
      icheck = 0
    end if
  end if
  !
  ! -- perform package convergence check
  if (icheck /= 0) then
    if (DELT > DZERO) then
      tled = DONE / DELT
    else
      tled = DZERO
    end if
    final_check: do ib = 1, this%ninterbeds
      idelay = this%idelay(ib)
      node = this%nodelist(ib)
      !
      ! -- skip nodelay interbeds
      if (idelay == 0) cycle
      !
      ! -- skip inactive cells
      if (this%ibound(node) < 1) cycle
      !
      ! -- evaluate the maximum head change in the interbed
      dh = this%dbdhmax(idelay)
      !
      ! -- evaluate difference between storage changes
      !    in the interbed and exchange between the interbed
      !    and the gwf cell
      area = this%dis%get_area(node)
      hcell = hnew(node)
      hcellold = hold(node)
      !
      ! -- calculate cell saturation
      call this%csub_calc_sat(node, hcell, hcellold, snnew, snold)
      !
      ! -- calculate the change in storage
      call this%csub_delay_calc_dstor(ib, hcell, stoe, stoi)
      v1 = (stoe + stoi) * area * this%rnb(ib) * tled
      !
      ! -- add water compressibility to storage term
      call this%csub_delay_calc_wcomp(ib, dwc)
      v1 = v1 + dwc * area * this%rnb(ib)
      !
      ! -- calculate the flow between the interbed and the cell
      call this%csub_delay_fc(ib, hcof, rhs)
      v2 = (-hcof * hcell - rhs) * area * this%rnb(ib)
      !
      ! -- calculate the difference between the interbed change in
      !    storage and the flow between the interbed and the cell
      df = v2 - v1
      !
      ! -- normalize by cell area and convert to a depth
      df = df * delt / area
      !
      ! -- evaluate magnitude of differences
      if (ifirst == 1) then
        ifirst = 0
        locdhmax = ib
        dhmax = dh
        locrmax = ib
        rmax = df
      else
        if (abs(dh) > abs(dhmax)) then
          locdhmax = ib
          dhmax = dh
        end if
        if (abs(df) > abs(rmax)) then
          locrmax = ib
          rmax = df
        end if
      end if
    end do final_check
    !
    ! -- set dpak and cpak
    ! -- update head error
    if (abs(dhmax) > abs(dpak)) then
      ipak = locdhmax
      dpak = dhmax
      write (cloc, "(a,'-',a)") trim(this%packName), 'head'
      cpak = cloc
    end if
    !
    ! -- update storage error
    if (abs(rmax) > abs(dpak)) then
      ipak = locrmax
      dpak = rmax
      write (cloc, "(a,'-',a)") trim(this%packName), 'storage'
      cpak = cloc
    end if
    !
    ! -- write convergence data to package csv
    if (this%ipakcsv /= 0) then
      !
      ! -- write the data
      call this%pakcsvtab%add_term(innertot)
      call this%pakcsvtab%add_term(totim)
      call this%pakcsvtab%add_term(kper)
      call this%pakcsvtab%add_term(kstp)
      call this%pakcsvtab%add_term(kiter)
      if (this%ndelaybeds > 0) then
        call this%pakcsvtab%add_term(dhmax)
        call this%pakcsvtab%add_term(locdhmax)
        call this%pakcsvtab%add_term(rmax)
        call this%pakcsvtab%add_term(locrmax)
      else
        call this%pakcsvtab%add_term('--')
        call this%pakcsvtab%add_term('--')
        call this%pakcsvtab%add_term('--')
        call this%pakcsvtab%add_term('--')
      end if
      !
      ! -- finalize the package csv
      if (iend == 1) then
        call this%pakcsvtab%finalize_table()
      end if
    end if
  end if
  end procedure csub_cc

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module procedure csub_cq
  ! -- modules
  use TdisModule, only: delt
  use ConstantsModule, only: LENBOUNDNAME, DZERO, DONE
  ! -- local variables
  integer(I4B) :: ib
  integer(I4B) :: idelay
  integer(I4B) :: ielastic
  integer(I4B) :: iconvert
  integer(I4B) :: node
  integer(I4B) :: nn
  integer(I4B) :: n
  integer(I4B) :: idiag
  real(DP) :: es
  real(DP) :: pcs
  real(DP) :: rho1
  real(DP) :: rho2
  real(DP) :: tled
  real(DP) :: tledm
  real(DP) :: es0
  real(DP) :: rrate
  real(DP) :: ratein
  real(DP) :: rateout
  real(DP) :: comp
  real(DP) :: compi
  real(DP) :: compe
  real(DP) :: area
  real(DP) :: h
  real(DP) :: h0
  real(DP) :: snnew
  real(DP) :: snold
  real(DP) :: hcof
  real(DP) :: rhs
  real(DP) :: stoe
  real(DP) :: stoi
  real(DP) :: b
  real(DP) :: q
  real(DP) :: rratewc
  ! -- for observations
  integer(I4B) :: iprobslocal
  ! -- formats
  !
  ! -- Suppress saving of simulated values; they
  !    will be saved at end of this procedure.
  iprobslocal = 0
  ratein = DZERO
  rateout = DZERO
  !
  ! -- coarse-grained coarse-grained storage
  do node = 1, this%dis%nodes
    idiag = this%dis%con%ia(node)
    area = this%dis%get_area(node)
    comp = DZERO
    rrate = DZERO
    rratewc = DZERO
    if (this%gwfiss == 0) then
      if (DELT > DZERO) then
        tled = DONE / DELT
      else
        tled = DZERO
      end if
      if (this%ibound(node) > 0 .and. this%cg_thickini(node) > DZERO) then
        !
        ! -- calculate coarse-grained storage terms
        call this%csub_cg_fc(node, tled, area, hnew(node), hold(node), &
                             hcof, rhs)
        rrate = hcof * hnew(node) - rhs
        !
        ! -- calculate compaction
        call this%csub_cg_calc_comp(node, hnew(node), hold(node), comp)
        !
        ! -- calculate coarse-grained water compressibility storage terms
        call this%csub_cg_wcomp_fc(node, tled, area, hnew(node), hold(node), &
                                   hcof, rhs)
        rratewc = hcof * hnew(node) - rhs
      end if
    end if
    !
    ! -- update coarse-grained storage and water
    !    compression variables
    this%cg_stor(node) = rrate
    this%cell_wcstor(node) = rratewc
    this%cell_thick(node) = this%cg_thick(node)
    !
    ! -- update incremental coarse-grained compaction
    this%cg_comp(node) = comp
    !
    !
    ! -- update states if required
    if (isuppress_output == 0) then
      !
      ! -- calculate strain and change in coarse-grained void ratio and thickness
      ! todo: consider moving error check in csub_cg_update to ot()
      if (this%iupdatematprop /= 0) then
        call this%csub_cg_update(node)
      end if
      !
      ! -- update total compaction
      this%cg_tcomp(node) = this%cg_tcomp(node) + comp
    end if
    !
    ! -- update flowja
    flowja(idiag) = flowja(idiag) + rrate
    flowja(idiag) = flowja(idiag) + rratewc
  end do
  !
  ! -- interbed storage
  !
  ! -- reset delay bed counters for the current time step
  if (this%ndelaybeds > 0) then
    this%idb_nconv_count(1) = 0
  end if
  !
  ! -- initialize tled
  tled = DONE
  !
  ! -- calculate budget terms for each interbed
  do ib = 1, this%ninterbeds
    rratewc = DZERO
    idelay = this%idelay(ib)
    ielastic = this%ielastic(ib)
    !
    ! -- calculate interbed thickness
    ! -- no delay interbeds
    if (idelay == 0) then
      b = this%thick(ib)
      ! -- delay interbeds
    else
      b = this%thick(ib) * this%rnb(ib)
    end if
    !
    ! -- set variables required for no-delay and delay interbeds
    node = this%nodelist(ib)
    idiag = this%dis%con%ia(node)
    area = this%dis%get_area(node)
    !
    ! -- add interbed thickness to cell thickness
    this%cell_thick(node) = this%cell_thick(node) + b
    !
    ! -- update budget terms if transient stress period
    if (this%gwfiss == 0) then
      if (DELT > DZERO) then
        tledm = DONE / DELT
      else
        tledm = DZERO
      end if
      !
      ! -- skip inactive and constant head cells
      if (this%ibound(node) < 1) cycle
      !
      ! -- no delay interbeds
      if (idelay == 0) then
        iconvert = this%iconvert(ib)
        stoi = DZERO
        !
        ! -- calculate compaction
        call this%csub_nodelay_calc_comp(ib, hnew(node), hold(node), comp, &
                                         rho1, rho2)
        !
        ! -- interbed stresses
        es = this%cg_es(node)
        pcs = this%pcs(ib)
        es0 = this%cg_es0(node)
        !
        ! -- calculate inelastic and elastic compaction
        if (ielastic > 0 .or. iconvert == 0) then
          stoe = comp
        else
          stoi = -pcs * rho2 + (rho2 * es)
          stoe = pcs * rho1 - (rho1 * es0)
        end if
        compe = stoe
        compi = stoi
        stoe = stoe * area
        stoi = stoi * area
        this%storagee(ib) = stoe * tledm
        this%storagei(ib) = stoi * tledm
        !
        ! -- update compaction
        this%comp(ib) = comp
        !
        ! -- update states if required
        if (isuppress_output == 0) then
          !
          ! -- calculate strain and change in interbed void ratio and thickness
          if (this%iupdatematprop /= 0) then
            call this%csub_nodelay_update(ib)
          end if
          !
          ! -- update total compaction
          this%tcomp(ib) = this%tcomp(ib) + comp
          this%tcompe(ib) = this%tcompe(ib) + compe
          this%tcompi(ib) = this%tcompi(ib) + compi
        end if
        !
        ! -- delay interbeds
      else
        h = hnew(node)
        h0 = hold(node)
        !
        ! -- calculate cell saturation
        call this%csub_calc_sat(node, h, h0, snnew, snold)
        !
        ! -- calculate inelastic and elastic storage contributions
        call this%csub_delay_calc_dstor(ib, h, stoe, stoi)
        this%storagee(ib) = stoe * area * this%rnb(ib) * tledm
        this%storagei(ib) = stoi * area * this%rnb(ib) * tledm
        !
        ! -- calculate flow across the top and bottom of the delay interbed
        q = this%csub_calc_delay_flow(ib, 1, h) * area * this%rnb(ib)
        this%dbflowtop(idelay) = q
        nn = this%ndelaycells
        q = this%csub_calc_delay_flow(ib, nn, h) * area * this%rnb(ib)
        this%dbflowbot(idelay) = q
        !
        ! -- update states if required
        if (isuppress_output == 0) then
          !
          ! -- calculate sum of compaction in delay interbed
          call this%csub_delay_calc_comp(ib, h, h0, comp, compi, compe)
          !
          ! - calculate strain and change in interbed void ratio and thickness
          ! todo: consider moving error check in csub_delay_update to ot()
          if (this%iupdatematprop /= 0) then
            call this%csub_delay_update(ib)
          end if
          !
          ! -- update total compaction for interbed
          this%tcomp(ib) = this%tcomp(ib) + comp
          this%tcompi(ib) = this%tcompi(ib) + compi
          this%tcompe(ib) = this%tcompe(ib) + compe
          !
          ! -- update total compaction for each delay bed cell
          do n = 1, this%ndelaycells
            this%dbtcomp(n, idelay) = this%dbtcomp(n, idelay) + &
                                      this%dbcomp(n, idelay)
          end do
          !
          ! -- check delay bed heads relative to the top and bottom of each
          !    delay bed cell for convertible and non-convertible gwf cells
          call this%csub_delay_head_check(ib)
        end if
      end if
      !
      ! -- interbed water compressibility
      !
      ! -- no-delay interbed
      if (idelay == 0) then
        call this%csub_nodelay_wcomp_fc(ib, node, tledm, area, &
                                        hnew(node), hold(node), hcof, rhs)
        rratewc = hcof * hnew(node) - rhs
        !
        ! -- delay interbed
      else
        call this%csub_delay_calc_wcomp(ib, q)
        rratewc = q * area * this%rnb(ib)
      end if
      this%cell_wcstor(node) = this%cell_wcstor(node) + rratewc
      !
      ! -- flowja
      flowja(idiag) = flowja(idiag) + rratewc
    else
      this%storagee(ib) = DZERO
      this%storagei(ib) = DZERO
      if (idelay /= 0) then
        this%dbflowtop(idelay) = DZERO
        this%dbflowbot(idelay) = DZERO
      end if
    end if
    !
    ! -- flowja
    flowja(idiag) = flowja(idiag) + this%storagee(ib)
    flowja(idiag) = flowja(idiag) + this%storagei(ib)
  end do
  !
  ! -- terminate if errors encountered when updating material properties
  if (this%iupdatematprop /= 0) then
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end if
  end procedure csub_cq

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module procedure csub_bd
  ! -- modules
  use TdisModule, only: delt
  use ConstantsModule, only: LENBOUNDNAME, DZERO, DONE
  use BudgetModule, only: BudgetType, rate_accumulator
  ! -- local
  real(DP) :: rin
  real(DP) :: rout
  !
  ! -- interbed elastic storage  (this%cg_stor)
  call rate_accumulator(this%cg_stor, rin, rout)
  call model_budget%addentry(rin, rout, delt, budtxt(1), &
                             isuppress_output, '            CSUB')
  if (this%ninterbeds > 0) then
    !
    ! -- interbed elastic storage (this%storagee)
    call rate_accumulator(this%storagee, rin, rout)
    call model_budget%addentry(rin, rout, delt, budtxt(2), &
                               isuppress_output, '            CSUB')
    !
    ! -- interbed elastic storage (this%storagei)
    call rate_accumulator(this%storagei, rin, rout)
    call model_budget%addentry(rin, rout, delt, budtxt(3), &
                               isuppress_output, '            CSUB')
  end if
  call rate_accumulator(this%cell_wcstor, rin, rout)
  call model_budget%addentry(rin, rout, delt, budtxt(4), &
                             isuppress_output, '            CSUB')
  return
  end procedure csub_bd

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module procedure csub_save_model_flows
  ! -- local variables
  character(len=1) :: cdatafmp = ' '
  character(len=1) :: editdesc = ' '
  integer(I4B) :: ibinun
  integer(I4B) :: iprint
  integer(I4B) :: nvaluesp
  integer(I4B) :: nwidthp
  integer(I4B) :: ib
  integer(I4B) :: node
  integer(I4B) :: naux
  real(DP) :: dinact
  real(DP) :: Q
  ! -- formats
  !
  ! -- Set unit number for binary output
  if (this%ipakcb < 0) then
    ibinun = icbcun
  elseif (this%ipakcb == 0) then
    ibinun = 0
  else
    ibinun = this%ipakcb
  end if
  if (icbcfl == 0) ibinun = 0
  !
  ! -- Record the storage rates if requested
  if (ibinun /= 0) then
    iprint = 0
    dinact = DZERO
    !
    ! -- coarse-grained storage (sske)
    call this%dis%record_array(this%cg_stor, this%iout, iprint, -ibinun, &
                               budtxt(1), cdatafmp, nvaluesp, &
                               nwidthp, editdesc, dinact)
    if (this%ninterbeds > 0) then
      naux = 0
      !
      ! -- interbed elastic storage
      call this%dis%record_srcdst_list_header(budtxt(2), &
                                              this%name_model, &
                                              this%name_model, &
                                              this%name_model, &
                                              this%packName, &
                                              naux, &
                                              this%auxname, &
                                              ibinun, &
                                              this%ninterbeds, &
                                              this%iout)
      do ib = 1, this%ninterbeds
        q = this%storagee(ib)
        node = this%nodelist(ib)
        call this%dis%record_mf6_list_entry(ibinun, node, node, q, naux, &
                                            this%auxvar(:, ib))
      end do
      !
      ! -- interbed inelastic storage
      call this%dis%record_srcdst_list_header(budtxt(3), &
                                              this%name_model, &
                                              this%name_model, &
                                              this%name_model, &
                                              this%packName, &
                                              naux, &
                                              this%auxname, &
                                              ibinun, &
                                              this%ninterbeds, &
                                              this%iout)
      do ib = 1, this%ninterbeds
        q = this%storagei(ib)
        node = this%nodelist(ib)
        call this%dis%record_mf6_list_entry(ibinun, node, node, q, naux, &
                                            this%auxvar(:, ib))
      end do
    end if
    !
    ! -- water compressibility
    call this%dis%record_array(this%cell_wcstor, this%iout, iprint, -ibinun, &
                               budtxt(4), cdatafmp, nvaluesp, &
                               nwidthp, editdesc, dinact)
  end if
  end procedure csub_save_model_flows

end submodule GwfCsubModuleBudget
