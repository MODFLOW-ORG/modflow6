!> @brief CSUB package output
!<
submodule(GwfCsubModule) GwfCsubModuleOutput

  character(len=LENBUDTXT), dimension(6) :: comptxt = & !< text labels for compaction terms
                                            ['CSUB-COMPACTION', &
                                             ' CSUB-INELASTIC', &
                                             '   CSUB-ELASTIC', &
                                             '  CSUB-INTERBED', &
                                             '    CSUB-COARSE', &
                                             ' CSUB-ZDISPLACE']

contains

  !> @brief Final processing for package
  !!
  !!  Final processing for the  CSUB package. This method generates the final
  !!  strain tables that are output so that the user can evaluate if calculated
  !!  strain rates in coarse-grained sediments and interbeds exceed 1 percent.
  !!
  !<
  module procedure csub_fp
  ! -- local variables
  character(len=LINELENGTH) :: title
  character(len=LINELENGTH) :: tag
  character(len=LINELENGTH) :: msg
  character(len=10) :: ctype
  character(len=20) :: cellid
  character(len=10) :: cflag
  integer(I4B) :: i
  integer(I4B) :: ib
  integer(I4B) :: i0
  integer(I4B) :: i1
  integer(I4B) :: node
  integer(I4B) :: nn
  integer(I4B) :: idelay
  integer(I4B) :: iexceed
  integer(I4B), parameter :: ncells = 20
  integer(I4B) :: nlen
  integer(I4B) :: ntabrows
  integer(I4B) :: ntabcols
  integer(I4B) :: ipos
  real(DP) :: b0
  real(DP) :: b1
  real(DP) :: strain
  real(DP) :: pctcomp
  integer(I4B), dimension(:), allocatable :: imap_sel
  integer(I4B), dimension(:), allocatable :: locs
  real(DP), dimension(:), allocatable :: pctcomp_arr
  !
  ! -- initialize locs
  allocate (locs(this%dis%ndim))
  !
  ! -- calculate and report strain for interbeds
  if (this%ninterbeds > 0) then
    nlen = min(ncells, this%ninterbeds)
    allocate (imap_sel(nlen))
    allocate (pctcomp_arr(this%ninterbeds))
    iexceed = 0
    do ib = 1, this%ninterbeds
      idelay = this%idelay(ib)
      b0 = this%thickini(ib)
      strain = this%tcomp(ib) / b0
      pctcomp = DHUNDRED * strain
      pctcomp_arr(ib) = pctcomp
      if (pctcomp >= DONE) then
        iexceed = iexceed + 1
      end if
    end do
    call selectn(imap_sel, pctcomp_arr, reverse=.TRUE.)
    !
    ! -- summary interbed strain table
    i0 = max(1, this%ninterbeds - ncells + 1)
    i1 = this%ninterbeds
    msg = ''
    if (iexceed /= 0) then
      write (msg, '(1x,a,1x,i0,1x,a,1x,i0,1x,a)') &
        'LARGEST', (i1 - i0 + 1), 'OF', this%ninterbeds, &
        'INTERBED STRAIN VALUES SHOWN'
      call write_message(msg, this%iout, skipbefore=1)
      !
      ! -- interbed strain data
      ! -- set title
      title = trim(adjustl(this%packName))//' PACKAGE INTERBED STRAIN SUMMARY'
      !
      ! -- determine the number of columns and rows
      ntabrows = nlen
      ntabcols = 9
      !
      ! -- setup table
      call table_cr(this%outputtab, this%packName, title)
      call this%outputtab%table_df(ntabrows, ntabcols, this%iout)
      !
      ! add columns
      tag = 'INTERBED NUMBER'
      call this%outputtab%initialize_column(tag, 10, alignment=TABCENTER)
      tag = 'INTERBED TYPE'
      call this%outputtab%initialize_column(tag, 10, alignment=TABCENTER)
      tag = 'CELLID'
      call this%outputtab%initialize_column(tag, 20, alignment=TABLEFT)
      tag = 'INITIAL THICKNESS'
      call this%outputtab%initialize_column(tag, 12, alignment=TABCENTER)
      tag = 'FINAL THICKNESS'
      call this%outputtab%initialize_column(tag, 12, alignment=TABCENTER)
      tag = 'TOTAL COMPACTION'
      call this%outputtab%initialize_column(tag, 12, alignment=TABCENTER)
      tag = 'FINAL STRAIN'
      call this%outputtab%initialize_column(tag, 12, alignment=TABCENTER)
      tag = 'PERCENT COMPACTION'
      call this%outputtab%initialize_column(tag, 12, alignment=TABCENTER)
      tag = 'FLAG'
      call this%outputtab%initialize_column(tag, 10, alignment=TABCENTER)
      !
      ! -- write data
      do i = 1, nlen
        ib = imap_sel(i)
        idelay = this%idelay(ib)
        b0 = this%thickini(ib)
        b1 = this%csub_calc_interbed_thickness(ib)
        if (idelay == 0) then
          ctype = 'no-delay'
        else
          ctype = 'delay'
          b0 = b0 * this%rnb(ib)
        end if
        strain = this%tcomp(ib) / b0
        pctcomp = DHUNDRED * strain
        if (pctcomp >= 5.0_DP) then
          cflag = '**>=5%'
        else if (pctcomp >= DONE) then
          cflag = '*>=1%'
        else
          cflag = ''
        end if
        node = this%nodelist(ib)
        call this%dis%noder_to_string(node, cellid)
        !
        ! -- fill table line
        call this%outputtab%add_term(ib)
        call this%outputtab%add_term(ctype)
        call this%outputtab%add_term(cellid)
        call this%outputtab%add_term(b0)
        call this%outputtab%add_term(b1)
        call this%outputtab%add_term(this%tcomp(ib))
        call this%outputtab%add_term(strain)
        call this%outputtab%add_term(pctcomp)
        call this%outputtab%add_term(cflag)
      end do
      write (this%iout, '(/1X,A,1X,I0,1X,A,1X,I0,1X,A,/1X,A,/1X,A)') &
        'PERCENT COMPACTION IS GREATER THAN OR EQUAL TO 1 PERCENT IN', &
        iexceed, 'OF', this%ninterbeds, 'INTERBED(S).', &
        'USE THE STRAIN_CSV_INTERBED OPTION TO OUTPUT A CSV '// &
        'FILE WITH PERCENT COMPACTION ', 'VALUES FOR ALL INTERBEDS.'
    else
      msg = 'PERCENT COMPACTION WAS LESS THAN 1 PERCENT IN ALL INTERBEDS'
      write (this%iout, '(/1X,A)') trim(adjustl(msg))
    end if
    !
    ! -- write csv file
    if (this%istrainib /= 0) then
      !
      ! -- determine the number of columns and rows
      ntabrows = this%ninterbeds
      ntabcols = 7
      if (this%dis%ndim > 1) then
        ntabcols = ntabcols + 1
      end if
      ntabcols = ntabcols + this%dis%ndim
      !
      ! -- setup table
      call table_cr(this%outputtab, this%packName, '')
      call this%outputtab%table_df(ntabrows, ntabcols, this%istrainib, &
                                   lineseparator=.FALSE., separator=',')
      !
      ! add columns
      tag = 'INTERBED_NUMBER'
      call this%outputtab%initialize_column(tag, 20, alignment=TABRIGHT)
      tag = 'INTERBED_TYPE'
      call this%outputtab%initialize_column(tag, 20, alignment=TABRIGHT)
      tag = 'NODE'
      call this%outputtab%initialize_column(tag, 10, alignment=TABRIGHT)
      if (this%dis%ndim == 2) then
        tag = 'LAYER'
        call this%outputtab%initialize_column(tag, 10, alignment=TABRIGHT)
        tag = 'ICELL2D'
        call this%outputtab%initialize_column(tag, 10, alignment=TABRIGHT)
      else
        tag = 'LAYER'
        call this%outputtab%initialize_column(tag, 10, alignment=TABRIGHT)
        tag = 'ROW'
        call this%outputtab%initialize_column(tag, 10, alignment=TABRIGHT)
        tag = 'COLUMN'
        call this%outputtab%initialize_column(tag, 10, alignment=TABRIGHT)
      end if
      tag = 'INITIAL_THICKNESS'
      call this%outputtab%initialize_column(tag, 20, alignment=TABRIGHT)
      tag = 'FINAL_THICKNESS'
      call this%outputtab%initialize_column(tag, 20, alignment=TABRIGHT)
      tag = 'TOTAL_COMPACTION'
      call this%outputtab%initialize_column(tag, 20, alignment=TABRIGHT)
      tag = 'TOTAL_STRAIN'
      call this%outputtab%initialize_column(tag, 20, alignment=TABRIGHT)
      tag = 'PERCENT_COMPACTION'
      call this%outputtab%initialize_column(tag, 20, alignment=TABRIGHT)
      !
      ! -- write data
      do ib = 1, this%ninterbeds
        idelay = this%idelay(ib)
        b0 = this%thickini(ib)
        b1 = this%csub_calc_interbed_thickness(ib)
        if (idelay == 0) then
          ctype = 'no-delay'
        else
          ctype = 'delay'
          b0 = b0 * this%rnb(ib)
        end if
        strain = this%tcomp(ib) / b0
        pctcomp = DHUNDRED * strain
        node = this%nodelist(ib)
        call this%dis%noder_to_array(node, locs)
        !
        ! -- fill table line
        call this%outputtab%add_term(ib)
        call this%outputtab%add_term(ctype)
        if (this%dis%ndim > 1) then
          call this%outputtab%add_term(this%dis%get_nodeuser(node))
        end if
        do ipos = 1, this%dis%ndim
          call this%outputtab%add_term(locs(ipos))
        end do
        call this%outputtab%add_term(b0)
        call this%outputtab%add_term(b1)
        call this%outputtab%add_term(this%tcomp(ib))
        call this%outputtab%add_term(strain)
        call this%outputtab%add_term(pctcomp)
      end do
    end if
    !
    ! -- deallocate temporary storage
    deallocate (imap_sel)
    deallocate (pctcomp_arr)
  end if
  !
  ! -- calculate and report strain for coarse-grained materials
  nlen = min(ncells, this%dis%nodes)
  allocate (imap_sel(nlen))
  allocate (pctcomp_arr(this%dis%nodes))
  iexceed = 0
  do node = 1, this%dis%nodes
    strain = DZERO
    if (this%cg_thickini(node) > DZERO) then
      strain = this%cg_tcomp(node) / this%cg_thickini(node)
    end if
    pctcomp = DHUNDRED * strain
    pctcomp_arr(node) = pctcomp
    if (pctcomp >= DONE) then
      iexceed = iexceed + 1
    end if
  end do
  call selectn(imap_sel, pctcomp_arr, reverse=.TRUE.)
  !
  ! -- summary coarse-grained strain table
  i0 = max(1, this%dis%nodes - ncells + 1)
  i1 = this%dis%nodes
  msg = ''
  if (iexceed /= 0) then
    write (msg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
      'LARGEST ', (i1 - i0 + 1), 'OF', this%dis%nodes, &
      'CELL COARSE-GRAINED VALUES SHOWN'
    call write_message(msg, this%iout, skipbefore=1)
    !
    ! -- set title
    title = trim(adjustl(this%packName))// &
            ' PACKAGE COARSE-GRAINED STRAIN SUMMARY'
    !
    ! -- determine the number of columns and rows
    ntabrows = nlen
    ntabcols = 7
    !
    ! -- setup table
    call table_cr(this%outputtab, this%packName, title)
    call this%outputtab%table_df(ntabrows, ntabcols, this%iout)
    !
    ! add columns
    tag = 'CELLID'
    call this%outputtab%initialize_column(tag, 20, alignment=TABLEFT)
    tag = 'INITIAL THICKNESS'
    call this%outputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'FINAL THICKNESS'
    call this%outputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'TOTAL COMPACTION'
    call this%outputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'FINAL STRAIN'
    call this%outputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'PERCENT COMPACTION'
    call this%outputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'FLAG'
    call this%outputtab%initialize_column(tag, 10, alignment=TABCENTER)
    ! -- write data
    do nn = 1, nlen
      node = imap_sel(nn)
      if (this%cg_thickini(node) > DZERO) then
        strain = this%cg_tcomp(node) / this%cg_thickini(node)
      else
        strain = DZERO
      end if
      pctcomp = DHUNDRED * strain
      if (pctcomp >= 5.0_DP) then
        cflag = '**>=5%'
      else if (pctcomp >= DONE) then
        cflag = '*>=1%'
      else
        cflag = ''
      end if
      call this%dis%noder_to_string(node, cellid)
      !
      ! -- fill table line
      call this%outputtab%add_term(cellid)
      call this%outputtab%add_term(this%cg_thickini(node))
      call this%outputtab%add_term(this%cg_thick(node))
      call this%outputtab%add_term(this%cg_tcomp(node))
      call this%outputtab%add_term(strain)
      call this%outputtab%add_term(pctcomp)
      call this%outputtab%add_term(cflag)
    end do
    write (this%iout, '(/1X,A,1X,I0,1X,A,1X,I0,1X,A,/1X,A,/1X,A)') &
      'COARSE-GRAINED STORAGE PERCENT COMPACTION IS GREATER THAN OR '// &
      'EQUAL TO 1 PERCENT IN', iexceed, 'OF', this%dis%nodes, 'CELL(S).', &
      'USE THE STRAIN_CSV_COARSE OPTION TO OUTPUT A CSV '// &
      'FILE WITH PERCENT COMPACTION ', 'VALUES FOR ALL CELLS.'
  else
    msg = 'COARSE-GRAINED STORAGE PERCENT COMPACTION WAS LESS THAN '// &
          '1 PERCENT IN ALL CELLS '
    write (this%iout, '(/1X,A)') trim(adjustl(msg))
  end if
  !
  ! -- write csv file
  if (this%istrainsk /= 0) then
    !
    ! -- determine the number of columns and rows
    ntabrows = this%dis%nodes
    ntabcols = 5
    if (this%dis%ndim > 1) then
      ntabcols = ntabcols + 1
    end if
    ntabcols = ntabcols + this%dis%ndim
    !
    ! -- setup table
    call table_cr(this%outputtab, this%packName, '')
    call this%outputtab%table_df(ntabrows, ntabcols, this%istrainsk, &
                                 lineseparator=.FALSE., separator=',')
    !
    ! add columns
    tag = 'NODE'
    call this%outputtab%initialize_column(tag, 10, alignment=TABRIGHT)
    if (this%dis%ndim == 2) then
      tag = 'LAYER'
      call this%outputtab%initialize_column(tag, 10, alignment=TABRIGHT)
      tag = 'ICELL2D'
      call this%outputtab%initialize_column(tag, 10, alignment=TABRIGHT)
    else
      tag = 'LAYER'
      call this%outputtab%initialize_column(tag, 10, alignment=TABRIGHT)
      tag = 'ROW'
      call this%outputtab%initialize_column(tag, 10, alignment=TABRIGHT)
      tag = 'COLUMN'
      call this%outputtab%initialize_column(tag, 10, alignment=TABRIGHT)
    end if
    tag = 'INITIAL_THICKNESS'
    call this%outputtab%initialize_column(tag, 20, alignment=TABRIGHT)
    tag = 'FINAL_THICKNESS'
    call this%outputtab%initialize_column(tag, 20, alignment=TABRIGHT)
    tag = 'TOTAL_COMPACTION'
    call this%outputtab%initialize_column(tag, 20, alignment=TABRIGHT)
    tag = 'TOTAL_STRAIN'
    call this%outputtab%initialize_column(tag, 20, alignment=TABRIGHT)
    tag = 'PERCENT_COMPACTION'
    call this%outputtab%initialize_column(tag, 20, alignment=TABRIGHT)
    !
    ! -- write data
    do node = 1, this%dis%nodes
      if (this%cg_thickini(node) > DZERO) then
        strain = this%cg_tcomp(node) / this%cg_thickini(node)
      else
        strain = DZERO
      end if
      pctcomp = DHUNDRED * strain
      call this%dis%noder_to_array(node, locs)
      !
      ! -- fill table line
      if (this%dis%ndim > 1) then
        call this%outputtab%add_term(this%dis%get_nodeuser(node))
      end if
      do ipos = 1, this%dis%ndim
        call this%outputtab%add_term(locs(ipos))
      end do
      call this%outputtab%add_term(this%cg_thickini(node))
      call this%outputtab%add_term(this%cg_thick(node))
      call this%outputtab%add_term(this%cg_tcomp(node))
      call this%outputtab%add_term(strain)
      call this%outputtab%add_term(pctcomp)
    end do
  end if
  !
  ! -- write a warning message for delay interbeds in non-convertible gwf
  !    cells that violate minimum head assumptions
  if (this%ndelaybeds > 0) then
    if (this%idb_nconv_count(2) > 0) then
      write (warnmsg, '(a,1x,a,1x,i0,1x,a,1x,a)') &
        'Delay interbed cell heads were less than the top of the interbed', &
        'cell in', this%idb_nconv_count(2), 'interbed cells in ', &
        'non-convertible GWF cells for at least one time step during '// &
        'the simulation.'
      call store_warning(warnmsg)
    end if
  end if
  !
  ! -- deallocate temporary storage
  deallocate (imap_sel)
  deallocate (locs)
  deallocate (pctcomp_arr)
  end procedure csub_fp

!> @brief Save and print dependent values for package
!!
!!  Method saves cell-by-cell compaction and z-displacement terms. The method
!!  also calls the method to process observation output.
!!
!<
  module procedure csub_ot_dv
  ! -- local variables
  character(len=1) :: cdatafmp = ' '
  character(len=1) :: editdesc = ' '
  integer(I4B) :: ibinun
  integer(I4B) :: iprint
  integer(I4B) :: nvaluesp
  integer(I4B) :: nwidthp
  integer(I4B) :: ib
  integer(I4B) :: node
  integer(I4B) :: nodem
  integer(I4B) :: nodeu
  integer(I4B) :: i
  integer(I4B) :: ii
  integer(I4B) :: idx_conn
  integer(I4B) :: k
  integer(I4B) :: ncpl
  integer(I4B) :: nlay
  integer(I4B) :: ihc
  real(DP) :: dinact
  real(DP) :: va_scale
  ! -- formats
  character(len=*), parameter :: fmtnconv = &
  "(/4x, 'DELAY INTERBED CELL HEADS IN ', i0, ' INTERBEDS IN', &
  &' NON-CONVERTIBLE GWF CELLS WERE LESS THAN THE TOP OF THE INTERBED CELL')"
  !
  ! -- Save compaction results
  !
  ! -- Set unit number for binary compaction and z-displacement output
  if (this%ioutcomp /= 0 .or. this%ioutzdisp /= 0) then
    ibinun = 1
  else
    ibinun = 0
  end if
  if (idvfl == 0) ibinun = 0
  !
  ! -- save compaction results
  if (ibinun /= 0) then
    iprint = 0
    dinact = DHNOFLO
    !
    ! -- fill buff with total compaction
    do node = 1, this%dis%nodes
      this%buff(node) = this%cg_tcomp(node)
    end do
    do ib = 1, this%ninterbeds
      node = this%nodelist(ib)
      this%buff(node) = this%buff(node) + this%tcomp(ib)
    end do
    !
    ! -- write compaction data to binary file
    if (this%ioutcomp /= 0) then
      ibinun = this%ioutcomp
      call this%dis%record_array(this%buff, this%iout, iprint, ibinun, &
                                 comptxt(1), cdatafmp, nvaluesp, &
                                 nwidthp, editdesc, dinact)
    end if
    !
    ! -- calculate z-displacement (subsidence) and write data to binary file
    if (this%ioutzdisp /= 0) then
      ibinun = this%ioutzdisp
      !
      ! -- initialize buffusr
      do nodeu = 1, this%dis%nodesuser
        this%buffusr(nodeu) = DZERO
      end do
      !
      ! -- fill buffusr with buff
      do node = 1, this%dis%nodes
        nodeu = this%dis%get_nodeuser(node)
        this%buffusr(nodeu) = this%buff(node)
      end do
      !
      ! -- calculate z-displacement
      ncpl = this%dis%get_ncpl()
      !
      ! -- disu
      if (this%dis%ndim == 1) then
        do node = this%dis%nodes, 1, -1
          do ii = this%dis%con%ia(node) + 1, this%dis%con%ia(node + 1) - 1
            !
            ! -- Set the m cell number
            nodem = this%dis%con%ja(ii)
            idx_conn = this%dis%con%jas(ii)
            !
            ! -- vertical connection
            ihc = this%dis%con%ihc(idx_conn)
            if (ihc == 0) then
              !
              ! -- node has an underlying cell
              if (node < nodem) then
                va_scale = this%dis%get_area_factor(node, idx_conn)
                this%buffusr(node) = this%buffusr(node) + &
                                     va_scale * this%buffusr(nodem)
              end if
            end if
          end do
        end do
        ! -- disv or dis
      else
        nlay = this%dis%nodesuser / ncpl
        do k = nlay - 1, 1, -1
          do i = 1, ncpl
            node = (k - 1) * ncpl + i
            nodem = k * ncpl + i
            this%buffusr(node) = this%buffusr(node) + this%buffusr(nodem)
          end do
        end do
      end if
      !
      ! -- fill buff with data from buffusr
      do nodeu = 1, this%dis%nodesuser
        node = this%dis%get_nodenumber_idx1(nodeu, 1)
        if (node > 0) then
          this%buff(node) = this%buffusr(nodeu)
        end if
      end do
      !
      ! -- write z-displacement
      call this%dis%record_array(this%buff, this%iout, iprint, ibinun, &
                                 comptxt(6), cdatafmp, nvaluesp, &
                                 nwidthp, editdesc, dinact)

    end if
  end if
  !
  ! -- Set unit number for binary inelastic interbed compaction
  if (this%ioutcompi /= 0) then
    ibinun = this%ioutcompi
  else
    ibinun = 0
  end if
  if (idvfl == 0) ibinun = 0
  !
  ! -- save inelastic interbed compaction results
  if (ibinun /= 0) then
    iprint = 0
    dinact = DHNOFLO
    !
    ! -- fill buff with inelastic interbed compaction
    do node = 1, this%dis%nodes
      this%buff(node) = DZERO
    end do
    do ib = 1, this%ninterbeds
      node = this%nodelist(ib)
      this%buff(node) = this%buff(node) + this%tcompi(ib)
    end do
    !
    ! -- write inelastic interbed compaction data to binary file
    call this%dis%record_array(this%buff, this%iout, iprint, ibinun, &
                               comptxt(2), cdatafmp, nvaluesp, &
                               nwidthp, editdesc, dinact)
  end if
  !
  ! -- Set unit number for binary elastic interbed compaction
  if (this%ioutcompe /= 0) then
    ibinun = this%ioutcompe
  else
    ibinun = 0
  end if
  if (idvfl == 0) ibinun = 0
  !
  ! -- save elastic interbed compaction results
  if (ibinun /= 0) then
    iprint = 0
    dinact = DHNOFLO
    !
    ! -- fill buff with elastic interbed compaction
    do node = 1, this%dis%nodes
      this%buff(node) = DZERO
    end do
    do ib = 1, this%ninterbeds
      node = this%nodelist(ib)
      this%buff(node) = this%buff(node) + this%tcompe(ib)
    end do
    !
    ! -- write elastic interbed compaction data to binary file
    call this%dis%record_array(this%buff, this%iout, iprint, ibinun, &
                               comptxt(3), cdatafmp, nvaluesp, &
                               nwidthp, editdesc, dinact)
  end if
  !
  ! -- Set unit number for binary interbed compaction
  if (this%ioutcompib /= 0) then
    ibinun = this%ioutcompib
  else
    ibinun = 0
  end if
  if (idvfl == 0) ibinun = 0
  !
  ! -- save interbed compaction results
  if (ibinun /= 0) then
    iprint = 0
    dinact = DHNOFLO
    !
    ! -- fill buff with interbed compaction
    do node = 1, this%dis%nodes
      this%buff(node) = DZERO
    end do
    do ib = 1, this%ninterbeds
      node = this%nodelist(ib)
      this%buff(node) = this%buff(node) + this%tcompe(ib) + this%tcompi(ib)
    end do
    !
    ! -- write interbed compaction data to binary file
    call this%dis%record_array(this%buff, this%iout, iprint, ibinun, &
                               comptxt(4), cdatafmp, nvaluesp, &
                               nwidthp, editdesc, dinact)
  end if
  !
  ! -- Set unit number for binary coarse-grained compaction
  if (this%ioutcomps /= 0) then
    ibinun = this%ioutcomps
  else
    ibinun = 0
  end if
  if (idvfl == 0) ibinun = 0
  !
  ! -- save coarse-grained compaction results
  if (ibinun /= 0) then
    iprint = 0
    dinact = DHNOFLO
    !
    ! -- fill buff with coarse-grained compaction
    do node = 1, this%dis%nodes
      this%buff(node) = this%cg_tcomp(node)
    end do
    !
    ! -- write coarse-grained compaction data to binary file
    call this%dis%record_array(this%buff, this%iout, iprint, ibinun, &
                               comptxt(5), cdatafmp, nvaluesp, &
                               nwidthp, editdesc, dinact)
  end if
  !
  ! -- check that final effective stress values for the time step
  !    are greater than zero
  if (this%gwfiss == 0) then
    call this%csub_cg_chk_stress()
  end if
  !
  ! -- update maximum count of delay interbeds that violate
  !    basic head assumptions for delay beds and write a message
  !    for delay interbeds in non-convertible gwf cells that
  !    violate these head assumptions
  if (this%ndelaybeds > 0) then
    if (this%idb_nconv_count(1) > this%idb_nconv_count(2)) then
      this%idb_nconv_count(2) = this%idb_nconv_count(1)
    end if
    if (this%idb_nconv_count(1) > 0) then
      write (this%iout, fmtnconv) this%idb_nconv_count(1)
    end if
  end if
  end procedure csub_ot_dv

  !> @brief Check effective stress values
  !!
  !!  Method checks calculated effective stress values to ensure that
  !!  effective stress values are positive. An error condition and message are
  !!  issued if calculated effective stress values are less than a small positive
  !!  value (DEM6).
  !!
  !<
  module procedure csub_cg_chk_stress
  ! -- modules
  use TdisModule, only: kper, kstp
  ! -- local variables
  character(len=20) :: cellid
  integer(I4B) :: ierr
  integer(I4B) :: iwarn
  integer(I4B) :: node
  real(DP) :: gs
  real(DP) :: bot
  real(DP) :: hcell
  real(DP) :: es
  real(DP) :: phead
  !
  ! -- initialize variables
  ierr = 0
  iwarn = 0
  !
  ! -- check effective stress in each cell (effective-stress formulation only)
  do node = 1, this%dis%nodes
    if (this%ibound(node) < 1) cycle
    if (this%lhead_based .EQV. .TRUE.) cycle
    bot = this%dis%bot(node)
    gs = this%cg_gs(node)
    es = this%cg_es(node)
    phead = gs - es
    hcell = phead + bot
    if (this%istrict_stress /= 0) then
      ! -- deprecated STRICT_EFFECTIVE_STRESS: terminate on negative stress
      if (es < DEM6) then
        ierr = ierr + 1
        call this%dis%noder_to_string(node, cellid)
        write (errmsg, '(a,g0,a,1x,a,1x,a,4(g0,a))') &
          'Small to negative effective stress (', es, ') in cell', &
          trim(adjustl(cellid)), '. (', es, ' = ', this%cg_gs(node), &
          ' - (', hcell, ' - ', bot, ').'
        call store_error(errmsg)
      end if
    else
      ! -- default: count cells with negative (regularized) effective stress
      if (es < DEM6) then
        iwarn = iwarn + 1
      end if
    end if
  end do
  !
  ! -- STRICT_EFFECTIVE_STRESS: write a summary error message and terminate
  if (ierr > 0) then
    write (errmsg, '(a,1x,i0,3(1x,a))') &
      'Solution: small to negative effective stress values in', ierr, &
      'cells can be eliminated by increasing storage values and/or ', &
      'adding/modifying stress boundaries to prevent water-levels from', &
      'exceeding the top of the model.'
    call store_error(errmsg)
    call store_error_filename(this%input_fname)
  end if
  !
  ! -- default: note the regularized cells and count the time step
  if (iwarn > 0) then
    this%nreg_ts = this%nreg_ts + 1
    write (this%iout, '(1x,a,1x,i0,1x,a,1x,i0,1x,a,1x,i0,a)') &
      'CSUB negative effective stress regularized in', iwarn, &
      'cell(s) in stress period', kper, 'time step', kstp, '.'
  end if
  end procedure csub_cg_chk_stress

end submodule GwfCsubModuleOutput
