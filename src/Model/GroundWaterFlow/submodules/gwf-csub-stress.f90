!> @brief CSUB package initial states and stresses
!<
submodule(GwfCsubModule) GwfCsubModuleStress
contains

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module procedure csub_cg_calc_stress
  ! -- local variables
  integer(I4B) :: node
  integer(I4B) :: ii
  integer(I4B) :: nn
  integer(I4B) :: m
  integer(I4B) :: idx_conn
  real(DP) :: gs
  real(DP) :: top
  real(DP) :: bot
  real(DP) :: thick
  real(DP) :: va_scale
  real(DP) :: hcell
  real(DP) :: hbar
  real(DP) :: gs_conn
  real(DP) :: es
  real(DP) :: phead
  real(DP) :: sadd
  !
  ! -- calculate geostatic stress if necessary
  if (this%iupdatestress /= 0) then
    do node = 1, this%dis%nodes
      !
      ! -- calculate geostatic stress for this node
      !    this represents the geostatic stress component
      !    for the cell
      top = this%dis%top(node)
      bot = this%dis%bot(node)
      thick = top - bot
      !
      ! -- calculate cell contribution to geostatic stress
      if (this%ibound(node) /= 0) then
        hcell = hnew(node)
      else
        hcell = bot
      end if
      !
      ! -- calculate corrected head (hbar)
      hbar = sQuadratic0sp(hcell, bot, this%satomega)
      !
      ! -- geostatic stress calculation
      if (hcell < top) then
        gs = (top - hbar) * this%sgm(node) + (hbar - bot) * this%sgs(node)
      else
        gs = thick * this%sgs(node)
      end if
      !
      ! -- cell contribution to geostatic stress
      this%cg_gs(node) = gs
    end do
    !
    ! -- add user specified overlying geostatic stress
    do nn = 1, this%nbound
      node = this%nodelistsig0(nn)
      sadd = this%sig0(nn)
      this%cg_gs(node) = this%cg_gs(node) + sadd
    end do
    !
    ! -- calculate geostatic stress above cell
    do node = 1, this%dis%nodes
      !
      ! -- geostatic stress of cell
      gs = this%cg_gs(node)
      !
      ! -- Add geostatic stress of overlying cells (ihc=0)
      !    m < node = m is vertically above node
      do ii = this%dis%con%ia(node) + 1, this%dis%con%ia(node + 1) - 1
        !
        ! -- Set the m cell number
        m = this%dis%con%ja(ii)
        idx_conn = this%dis%con%jas(ii)
        !
        ! -- vertical connection
        if (this%dis%con%ihc(idx_conn) == 0) then
          !
          ! -- node has an overlying cell
          if (m < node) then
            !
            ! -- dis and disv discretization
            if (this%dis%ndim /= 1) then
              gs = gs + this%cg_gs(m)
              !
              ! -- disu discretization
            else
              va_scale = this%dis%get_area_factor(node, idx_conn)
              gs_conn = this%cg_gs(m)
              gs = gs + (gs_conn * va_scale)
            end if
          end if
        end if
      end do
      !
      ! -- geostatic stress for cell with geostatic stress
      !    of overlying cells
      this%cg_gs(node) = gs
    end do
  end if
  !
  ! -- save effective stress from the last iteration and
  !    calculate the new effective stress for a cell
  do node = 1, this%dis%nodes
    top = this%dis%top(node)
    bot = this%dis%bot(node)
    if (this%ibound(node) /= 0) then
      hcell = hnew(node)
    else
      hcell = bot
    end if
    !
    ! -- calculate corrected head (hbar)
    hbar = sQuadratic0sp(hcell, bot, this%satomega)
    !
    ! -- calculate pressure head
    phead = hbar - bot
    !
    ! -- calculate effective stress
    es = this%cg_gs(node) - phead
    this%cg_es(node) = es
  end do
  end procedure csub_cg_calc_stress

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module procedure csub_set_initial_state
  ! -- local variables
  character(len=LINELENGTH) :: title
  character(len=LINELENGTH) :: tag
  character(len=20) :: cellid
  integer(I4B) :: ib
  integer(I4B) :: node
  integer(I4B) :: n
  integer(I4B) :: idelay
  integer(I4B) :: ntabrows
  integer(I4B) :: ntabcols
  real(DP) :: pcs0
  real(DP) :: pcs
  real(DP) :: fact
  real(DP) :: top
  real(DP) :: bot
  real(DP) :: void_ratio
  real(DP) :: es
  real(DP) :: znode
  real(DP) :: hcell
  real(DP) :: hbar
  real(DP) :: dzhalf
  real(DP) :: zbot
  real(DP) :: dbpcs
  !
  ! -- update geostatic load calculation
  call this%csub_cg_calc_stress(nodes, hnew)
  !
  ! -- initialize coarse-grained material effective stress
  !    for the previous time step and the previous iteration
  do node = 1, nodes
    this%cg_es0(node) = this%cg_es(node)
  end do
  !
  ! -- initialize interbed initial states
  do ib = 1, this%ninterbeds
    idelay = this%idelay(ib)
    node = this%nodelist(ib)
    top = this%dis%top(node)
    bot = this%dis%bot(node)
    hcell = hnew(node)
    pcs = this%pcs(ib)
    pcs0 = pcs
    if (this%ispecified_pcs == 0) then
      ! relative pcs...subtract head (u) from sigma'
      if (this%ipch /= 0) then
        pcs = this%cg_es(node) - pcs0
      else
        pcs = this%cg_es(node) + pcs0
      end if
    else
      ! specified pcs...subtract head (u) from sigma
      if (this%ipch /= 0) then
        pcs = this%cg_gs(node) - (pcs0 - bot)
      end if
      if (pcs < this%cg_es(node)) then
        pcs = this%cg_es(node)
      end if
    end if
    this%pcs(ib) = pcs
    !
    ! -- delay bed initial states
    if (idelay /= 0) then
      dzhalf = DHALF * this%dbdzini(1, idelay)
      !
      ! -- fill delay bed head with aquifer head or offset from aquifer head
      !    heads need to be filled first since used to calculate
      !    the effective stress for each delay bed
      do n = 1, this%ndelaycells
        if (this%ispecified_dbh == 0) then
          this%dbh(n, idelay) = hcell + this%dbh(n, idelay)
        else
          this%dbh(n, idelay) = hcell
        end if
        this%dbh0(n, idelay) = this%dbh(n, idelay)
      end do
      !
      ! -- fill delay bed effective stress
      call this%csub_delay_calc_stress(ib, hcell)
      !
      ! -- fill delay bed pcs
      pcs = this%pcs(ib)
      do n = 1, this%ndelaycells
        zbot = this%dbz(n, idelay) - dzhalf
        ! -- adjust pcs to bottom of each delay bed cell
        !    not using adjust_eff_stress() since smoothing not required
        dbpcs = pcs - (zbot - bot) * (this%sgs(node) - DONE)
        this%dbpcs(n, idelay) = dbpcs
        !
        ! -- initialize effective stress for previous time step
        this%dbes0(n, idelay) = this%dbes(n, idelay)
      end do
    end if
  end do
  !
  ! -- scale coarse-grained materials cr
  do node = 1, nodes
    top = this%dis%top(node)
    bot = this%dis%bot(node)
    !
    ! -- user-specified specific storage
    if (this%istoragec == 1) then
      !
      ! -- retain specific storage values since they are constant
      if (this%lhead_based .EQV. .TRUE.) then
        fact = DONE
        !
        ! -- convert specific storage values since they are simulated to
        !    be a function of the average effective stress
      else
        void_ratio = porosity_to_void_ratio(this%cg_theta(node))
        es = this%cg_es(node)
        hcell = hnew(node)
        !
        ! -- calculate corrected head (hbar)
        hbar = sQuadratic0sp(hcell, bot, this%satomega)
        !
        ! -- calculate znode and factor
        znode = this%csub_calc_znode(top, bot, hbar)
        fact = adjust_eff_stress(es, bot, znode, this%sgs(node))
        fact = fact * (DONE + void_ratio)
      end if
      !
      ! -- user-specified compression indices - multiply by dlog10es
    else
      fact = dlog10es
    end if
    this%cg_ske_cr(node) = this%cg_ske_cr(node) * fact
    !
    ! -- write error message if negative compression indices
    if (fact <= DZERO) then
      call this%dis%noder_to_string(node, cellid)
      write (errmsg, '(a,1x,a,a)') &
        'Negative recompression index calculated for cell', &
        trim(adjustl(cellid)), '.'
      call store_error(errmsg)
    end if
  end do
  !
  ! -- scale interbed cc and cr
  do ib = 1, this%ninterbeds
    idelay = this%idelay(ib)
    node = this%nodelist(ib)
    top = this%dis%top(node)
    bot = this%dis%bot(node)
    !
    ! -- user-specified specific storage
    if (this%istoragec == 1) then
      !
      ! -- retain specific storage values since they are constant
      if (this%lhead_based .EQV. .TRUE.) then
        fact = DONE
        !
        ! -- convert specific storage values since they are simulated to
        !    be a function of the average effective stress
      else
        void_ratio = porosity_to_void_ratio(this%theta(ib))
        es = this%cg_es(node)
        hcell = hnew(node)
        !
        ! -- calculate corrected head (hbar)
        hbar = sQuadratic0sp(hcell, bot, this%satomega)
        !
        ! -- calculate zone and factor
        znode = this%csub_calc_znode(top, bot, hbar)
        fact = adjust_eff_stress(es, bot, znode, this%sgs(node))
        fact = fact * (DONE + void_ratio)
      end if
      !
      ! -- user-specified compression indices - multiply by dlog10es
    else
      fact = dlog10es
    end if
    this%ci(ib) = this%ci(ib) * fact
    this%rci(ib) = this%rci(ib) * fact
    !
    ! -- write error message if negative compression indices
    if (fact <= DZERO) then
      call this%dis%noder_to_string(node, cellid)
      write (errmsg, '(a,1x,i0,2(1x,a),a)') &
        'Negative compression indices calculated for interbed', ib, &
        'in cell', trim(adjustl(cellid)), '.'
      call store_error(errmsg)
    end if
  end do
  !
  ! -- write current stress and initial preconsolidation stress
  if (this%iprpak == 1) then
    ! -- set title
    title = trim(adjustl(this%packName))// &
            ' PACKAGE CALCULATED INITIAL INTERBED STRESSES AT THE CELL BOTTOM'
    !
    ! -- determine the number of columns and rows
    ntabrows = this%ninterbeds
    ntabcols = 5
    if (this%inamedbound /= 0) then
      ntabcols = ntabcols + 1
    end if
    !
    ! -- setup table
    call table_cr(this%inputtab, this%packName, title)
    call this%inputtab%table_df(ntabrows, ntabcols, this%iout)
    !
    ! add columns
    tag = 'INTERBED NUMBER'
    call this%inputtab%initialize_column(tag, 10, alignment=TABLEFT)
    tag = 'CELLID'
    call this%inputtab%initialize_column(tag, 20)
    tag = 'GEOSTATIC STRESS'
    call this%inputtab%initialize_column(tag, 16)
    tag = 'EFFECTIVE STRESS'
    call this%inputtab%initialize_column(tag, 16)
    tag = 'PRECONSOLIDATION STRESS'
    call this%inputtab%initialize_column(tag, 16)
    if (this%inamedbound /= 0) then
      tag = 'BOUNDNAME'
      call this%inputtab%initialize_column(tag, LENBOUNDNAME, &
                                           alignment=TABLEFT)
    end if
    !
    ! -- write the data
    do ib = 1, this%ninterbeds
      node = this%nodelist(ib)
      call this%dis%noder_to_string(node, cellid)
      !
      ! -- write the columns
      call this%inputtab%add_term(ib)
      call this%inputtab%add_term(cellid)
      call this%inputtab%add_term(this%cg_gs(node))
      call this%inputtab%add_term(this%cg_es(node))
      call this%inputtab%add_term(this%pcs(ib))
      if (this%inamedbound /= 0) then
        call this%inputtab%add_term(this%boundname(ib))
      end if
    end do
    !
    ! -- write effective stress and preconsolidation stress
    !    for delay beds
    ! -- set title
    title = trim(adjustl(this%packName))// &
            ' PACKAGE CALCULATED INITIAL DELAY INTERBED STRESSES'
    !
    ! -- determine the number of columns and rows
    ntabrows = 0
    do ib = 1, this%ninterbeds
      idelay = this%idelay(ib)
      if (idelay /= 0) then
        ntabrows = ntabrows + this%ndelaycells
      end if
    end do
    ntabcols = 6
    if (this%inamedbound /= 0) then
      ntabcols = ntabcols + 1
    end if
    !
    ! -- setup table
    call table_cr(this%inputtab, this%packName, title)
    call this%inputtab%table_df(ntabrows, ntabcols, this%iout)
    !
    ! add columns
    tag = 'INTERBED NUMBER'
    call this%inputtab%initialize_column(tag, 10, alignment=TABLEFT)
    tag = 'CELLID'
    call this%inputtab%initialize_column(tag, 20)
    tag = 'DELAY CELL'
    call this%inputtab%initialize_column(tag, 10, alignment=TABLEFT)
    tag = 'GEOSTATIC STRESS'
    call this%inputtab%initialize_column(tag, 16)
    tag = 'EFFECTIVE STRESS'
    call this%inputtab%initialize_column(tag, 16)
    tag = 'PRECONSOLIDATION STRESS'
    call this%inputtab%initialize_column(tag, 16)
    if (this%inamedbound /= 0) then
      tag = 'BOUNDNAME'
      call this%inputtab%initialize_column(tag, LENBOUNDNAME, &
                                           alignment=TABLEFT)
    end if
    !
    ! -- write the data
    do ib = 1, this%ninterbeds
      idelay = this%idelay(ib)
      if (idelay /= 0) then
        node = this%nodelist(ib)
        call this%dis%noder_to_string(node, cellid)
        !
        ! -- write the columns
        do n = 1, this%ndelaycells
          if (n == 1) then
            call this%inputtab%add_term(ib)
            call this%inputtab%add_term(cellid)
          else
            call this%inputtab%add_term(' ')
            call this%inputtab%add_term(' ')
          end if
          call this%inputtab%add_term(n)
          call this%inputtab%add_term(this%dbgeo(n, idelay))
          call this%inputtab%add_term(this%dbes(n, idelay))
          call this%inputtab%add_term(this%dbpcs(n, idelay))
          if (this%inamedbound /= 0) then
            if (n == 1) then
              call this%inputtab%add_term(this%boundname(ib))
            else
              call this%inputtab%add_term(' ')
            end if
          end if
        end do
      end if
    end do
    !
    ! -- write calculated compression indices
    if (this%istoragec == 1) then
      if (this%lhead_based .EQV. .FALSE.) then
        ! -- set title
        title = trim(adjustl(this%packName))// &
                ' PACKAGE COMPRESSION INDICES'
        !
        ! -- determine the number of columns and rows
        ntabrows = this%ninterbeds
        ntabcols = 4
        if (this%inamedbound /= 0) then
          ntabcols = ntabcols + 1
        end if
        !
        ! -- setup table
        call table_cr(this%inputtab, this%packName, title)
        call this%inputtab%table_df(ntabrows, ntabcols, this%iout)
        !
        ! add columns
        tag = 'INTERBED NUMBER'
        call this%inputtab%initialize_column(tag, 10, alignment=TABLEFT)
        tag = 'CELLID'
        call this%inputtab%initialize_column(tag, 20)
        tag = 'CC'
        call this%inputtab%initialize_column(tag, 16)
        tag = 'CR'
        call this%inputtab%initialize_column(tag, 16)
        if (this%inamedbound /= 0) then
          tag = 'BOUNDNAME'
          call this%inputtab%initialize_column(tag, LENBOUNDNAME, &
                                               alignment=TABLEFT)
        end if
        !
        ! -- write the data
        do ib = 1, this%ninterbeds
          fact = DONE / dlog10es
          node = this%nodelist(ib)
          call this%dis%noder_to_string(node, cellid)
          !
          ! -- write the columns
          call this%inputtab%add_term(ib)
          call this%inputtab%add_term(cellid)
          call this%inputtab%add_term(this%ci(ib) * fact)
          call this%inputtab%add_term(this%rci(ib) * fact)
          if (this%inamedbound /= 0) then
            call this%inputtab%add_term(this%boundname(ib))
          end if
        end do
      end if
    end if
  end if
  !
  ! -- terminate if any initialization errors have been detected
  if (count_errors() > 0) then
    call store_error_filename(this%input_fname)
  end if
  !
  ! -- set initialized
  this%initialized = 1
  !
  ! -- set flag to retain initial stresses for entire simulation
  if (this%lhead_based .EQV. .TRUE.) then
    this%iupdatestress = 0
  end if
  end procedure csub_set_initial_state

end submodule GwfCsubModuleStress
