!> @brief CSUB package input processing
!<
submodule(GwfCsubModule) GwfCsubModuleInput
contains

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine source_options(this)
    ! -- modules
    use ConstantsModule, only: MAXCHARLEN, DZERO, MNORMAL, DEM3
    use MemoryManagerModule, only: mem_reallocate
    use MemoryManagerExtModule, only: mem_set_value
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: getunit, urdaux, openfile
    use GwfCsubInputModule, only: GwfCsubParamFoundType
    use SourceCommonModule, only: filein_fname
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    ! -- local variables
    integer(I4B), pointer :: ibs
    integer(I4B) :: inobs
    integer(I4B), pointer :: iei_smoothing
    integer(I4B), pointer :: istrict
    character(len=LINELENGTH) :: csv_interbed, csv_coarse
    character(len=LINELENGTH) :: cmp_fn, ecmp_fn, iecmp_fn, ibcmp_fn, cmpcoarse_fn
    character(len=LINELENGTH) :: zdisp_fn, pkg_converge_fn
    type(GwfCsubParamFoundType) :: found
    logical(LGP) :: warn_estress_lag = .false.

    ! -- allocate and initialize variables
    allocate (ibs)
    ibs = 0

    ! -- update defaults from input context
    call mem_set_value(this%inamedbound, 'BOUNDNAMES', this%input_mempath, &
                       found%boundnames)
    call mem_set_value(this%iprpak, 'PRINT_INPUT', this%input_mempath, &
                       found%print_input)
    call mem_set_value(this%ipakcb, 'SAVE_FLOWS', this%input_mempath, &
                       found%save_flows)
    call mem_set_value(this%gammaw, 'GAMMAW', this%input_mempath, found%gammaw)
    call mem_set_value(this%beta, 'BETA', this%input_mempath, found%beta)
    ! -- ELASTIC_INELASTIC_SMOOTHING sets the default smoothing window (DEM3 * pcs)
    allocate (iei_smoothing)
    iei_smoothing = 0
    call mem_set_value(iei_smoothing, 'EI_SMOOTHING', &
                       this%input_mempath, found%ei_smoothing)
    if (found%ei_smoothing) then
      this%pcsomega = DEM3
    end if
    deallocate (iei_smoothing)
    ! -- STRICT_EFFECTIVE_STRESS terminates on negative effective stress
    allocate (istrict)
    istrict = 0
    call mem_set_value(istrict, 'STRICT_STRESS', this%input_mempath, &
                       found%strict_stress)
    if (found%strict_stress) then
      this%istrict_stress = 1
    end if
    deallocate (istrict)
    call mem_set_value(this%ipch, 'HEAD_BASED', this%input_mempath, &
                       found%head_based)
    call mem_set_value(this%ipch, 'PRECON_HEAD', this%input_mempath, &
                       found%precon_head)
    call mem_set_value(this%ndelaycells, 'NDELAYCELLS', this%input_mempath, &
                       found%ndelaycells)
    call mem_set_value(this%istoragec, 'ICOMPRESS', this%input_mempath, &
                       found%icompress)
    call mem_set_value(this%iupdatematprop, 'MATPROP', this%input_mempath, &
                       found%matprop)
    call mem_set_value(this%icellf, 'CELL_FRACTION', this%input_mempath, &
                       found%cell_fraction)
    call mem_set_value(ibs, 'INTERBED_STATE', this%input_mempath, &
                       found%interbed_state)
    call mem_set_value(this%ispecified_pcs, 'PRECON_STRESS', this%input_mempath, &
                       found%precon_stress)
    call mem_set_value(this%ispecified_dbh, 'DELAY_HEAD', this%input_mempath, &
                       found%delay_head)
    call mem_set_value(this%ieslag, 'STRESS_LAG', this%input_mempath, &
                       found%stress_lag)
    call mem_set_value(csv_interbed, 'INTERBEDSTRAINFN', this%input_mempath, &
                       found%interbedstrainfn)
    call mem_set_value(csv_coarse, 'COARSESTRAINFN', this%input_mempath, &
                       found%coarsestrainfn)
    call mem_set_value(cmp_fn, 'CMPFN', this%input_mempath, found%cmpfn)
    call mem_set_value(ecmp_fn, 'ELASTICCMPFN', this%input_mempath, &
                       found%elasticcmpfn)
    call mem_set_value(iecmp_fn, 'INELASTICCMPFN', this%input_mempath, &
                       found%inelasticcmpfn)
    call mem_set_value(ibcmp_fn, 'INTERBEDCMPFN', this%input_mempath, &
                       found%interbedcmpfn)
    call mem_set_value(cmpcoarse_fn, 'CMPCOARSEFN', this%input_mempath, &
                       found%cmpcoarsefn)
    call mem_set_value(zdisp_fn, 'ZDISPFN', this%input_mempath, found%zdispfn)
    call mem_set_value(pkg_converge_fn, 'PKGCONVERGEFN', this%input_mempath, &
                       found%pkgconvergefn)

    ! -- enforce 0 or 1 OBS6_FILENAME entries in option block
    if (filein_fname(this%obs%inputFilename, 'OBS6_FILENAME', &
                     this%input_mempath, this%input_fname)) then
      this%obs%active = .true.
      inobs = GetUnit()
      call openfile(inobs, this%iout, this%obs%inputFilename, 'OBS')
      this%obs%inUnitObs = inobs
      this%inobspkg = inobs
      call this%obs%obs_df(this%iout, this%packName, this%filtyp, this%dis)
      call this%csub_df_obs()
    end if

    ! -- update input dependent internal state
    if (found%save_flows) this%ipakcb = -1
    if (found%head_based) then
      this%lhead_based = .TRUE.
      if (this%ieslag /= 0) then
        this%ieslag = 0
        warn_estress_lag = .true.
      end if
    end if
    if (found%icompress) this%istoragec = 0
    if (found%interbed_state) then
      this%ispecified_pcs = 1
      this%ispecified_dbh = 1
    end if
    if (found%gammaw .or. found%beta) then
      this%brg = this%gammaw * this%beta
    end if

    ! fileout options
    if (found%interbedstrainfn) then
      this%istrainib = getunit()
      call openfile(this%istrainib, this%iout, csv_interbed, 'CSV_OUTPUT', &
                    filstat_opt='REPLACE', mode_opt=MNORMAL)
    end if
    if (found%coarsestrainfn) then
      this%istrainsk = getunit()
      call openfile(this%istrainsk, this%iout, csv_coarse, 'CSV_OUTPUT', &
                    filstat_opt='REPLACE', mode_opt=MNORMAL)
    end if
    if (found%cmpfn) then
      this%ioutcomp = getunit()
      call openfile(this%ioutcomp, this%iout, cmp_fn, 'DATA(BINARY)', &
                    form, access, 'REPLACE', mode_opt=MNORMAL)
    end if
    if (found%elasticcmpfn) then
      this%ioutcompe = getunit()
      call openfile(this%ioutcompe, this%iout, ecmp_fn, &
                    'DATA(BINARY)', form, access, 'REPLACE', &
                    mode_opt=MNORMAL)
    end if
    if (found%inelasticcmpfn) then
      this%ioutcompi = getunit()
      call openfile(this%ioutcompi, this%iout, iecmp_fn, &
                    'DATA(BINARY)', form, access, 'REPLACE', &
                    mode_opt=MNORMAL)
    end if
    if (found%interbedcmpfn) then
      this%ioutcompib = getunit()
      call openfile(this%ioutcompib, this%iout, ibcmp_fn, &
                    'DATA(BINARY)', form, access, 'REPLACE', &
                    mode_opt=MNORMAL)
    end if
    if (found%cmpcoarsefn) then
      this%ioutcomps = getunit()
      call openfile(this%ioutcomps, this%iout, cmpcoarse_fn, &
                    'DATA(BINARY)', form, access, 'REPLACE', &
                    mode_opt=MNORMAL)
    end if
    if (found%zdispfn) then
      this%ioutzdisp = getunit()
      call openfile(this%ioutzdisp, this%iout, zdisp_fn, &
                    'DATA(BINARY)', form, access, 'REPLACE', &
                    mode_opt=MNORMAL)
    end if
    if (found%pkgconvergefn) then
      this%ipakcsv = getunit()
      call openfile(this%ipakcsv, this%iout, pkg_converge_fn, 'CSV', &
                    filstat_opt='REPLACE', mode_opt=MNORMAL)
    end if

    ! -- log user options
    call this%log_options(warn_estress_lag)

    ! -- cleanup
    deallocate (ibs)
  end subroutine source_options

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine log_options(this, warn_estress_lag)
    ! -- modules
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    logical(LGP), intent(in) :: warn_estress_lag
    ! -- local variables
    ! -- formats
    character(len=*), parameter :: fmtts = &
      &"(4x,'TIME-SERIES DATA WILL BE READ FROM FILE: ',a)"
    character(len=*), parameter :: fmtflow = &
      &"(4x,'FLOWS WILL BE SAVED TO FILE: ',a,/4x,'OPENED ON UNIT: ',I7)"
    character(len=*), parameter :: fmtflow2 = &
      &"(4x,'FLOWS WILL BE SAVED TO BUDGET FILE SPECIFIED IN OUTPUT CONTROL')"
    character(len=*), parameter :: fmtssessv = &
      &"(4x,'USING SSE AND SSV INSTEAD OF CR AND CC.')"
    character(len=*), parameter :: fmtoffset = &
      &"(4x,'INITIAL_STRESS TREATED AS AN OFFSET.')"
    character(len=*), parameter :: fmtopt = &
      &"(4x,A)"
    character(len=*), parameter :: fmtopti = &
      &"(4x,A,1X,I0)"
    character(len=*), parameter :: fmtoptr = &
      &"(4x,A,1X,G0)"
    character(len=*), parameter :: fmtfileout = &
      "(4x,'CSUB ',1x,a,1x,' WILL BE SAVED TO FILE: ',a,/4x,&
      &'OPENED ON UNIT: ',I7)"
    !
    ! -- write messages for options
    write (this%iout, '(//2(1X,A))') trim(adjustl(this%packName)), &
      'PACKAGE SETTINGS'
    write (this%iout, fmtopti) 'NUMBER OF DELAY CELLS =', &
      this%ndelaycells
    if (this%lhead_based .EQV. .TRUE.) then
      write (this%iout, '(4x,a)') &
        'HEAD-BASED FORMULATION'
    else
      write (this%iout, '(4x,a)') &
        'EFFECTIVE-STRESS FORMULATION'
    end if
    if (this%istoragec == 0) then
      write (this%iout, '(4x,a,1(/,6x,a))') &
        'COMPRESSION INDICES WILL BE SPECIFIED INSTEAD OF ELASTIC AND', &
        'INELASTIC SPECIFIC STORAGE COEFFICIENTS'
    else
      write (this%iout, '(4x,a,1(/,6x,a))') &
        'ELASTIC AND INELASTIC SPECIFIC STORAGE COEFFICIENTS WILL BE ', &
        'SPECIFIED'
    end if
    if (this%iupdatematprop /= 1) then
      write (this%iout, '(4x,a,1(/,6x,a))') &
        'THICKNESS AND VOID RATIO WILL NOT BE ADJUSTED DURING THE', &
        'SIMULATION'
    else
      write (this%iout, '(4x,a)') &
        'THICKNESS AND VOID RATIO WILL BE ADJUSTED DURING THE SIMULATION'
    end if
    if (this%icellf /= 1) then
      write (this%iout, '(4x,a)') &
        'INTERBED THICKNESS WILL BE SPECIFIED AS A THICKNESS'
    else
      write (this%iout, '(4x,a,1(/,6x,a))') &
        'INTERBED THICKNESS WILL BE SPECIFIED AS A AS A CELL FRACTION'
    end if
    if (this%ispecified_pcs /= 1) then
      if (this%ipch /= 0) then
        write (this%iout, '(4x,a,1(/,6x,a))') &
          'PRECONSOLIDATION HEAD WILL BE SPECIFIED RELATIVE TO INITIAL', &
          'STRESS CONDITIONS'
      else
        write (this%iout, '(4x,a,1(/,6x,a))') &
          'PRECONSOLIDATION STRESS WILL BE SPECIFIED RELATIVE TO INITIAL', &
          'STRESS CONDITIONS'
      end if
    else
      if (this%ipch /= 0) then
        write (this%iout, '(4x,a,1(/,6x,a))') &
          'PRECONSOLIDATION HEAD WILL BE SPECIFIED AS ABSOLUTE VALUES', &
          'INSTEAD OF RELATIVE TO INITIAL HEAD CONDITIONS'
      else
        write (this%iout, '(4x,a,1(/,6x,a))') &
          'PRECONSOLIDATION STRESS WILL BE SPECIFIED AS ABSOLUTE VALUES', &
          'INSTEAD OF RELATIVE TO INITIAL STRESS CONDITIONS'
      end if
    end if
    if (this%ispecified_dbh /= 1) then
      write (this%iout, '(4x,a,1(/,6x,a))') &
        'DELAY INTERBED HEADS WILL BE SPECIFIED RELATIVE TO INITIAL ', &
        'GWF HEADS'
    else
      write (this%iout, '(4x,a,1(/,6x,a))') &
        'DELAY INTERBED HEADS WILL BE SPECIFIED AS ABSOLUTE VALUES INSTEAD', &
        'OF RELATIVE TO INITIAL GWF HEADS'
    end if
    !
    if (this%lhead_based .EQV. .FALSE.) then
      if (this%ieslag /= 0) then
        write (this%iout, '(4x,a,1(/,6x,a))') &
          'SPECIFIC STORAGE VALUES WILL BE CALCULATED USING THE EFFECTIVE', &
          'STRESS FROM THE PREVIOUS TIME STEP'
      else
        write (this%iout, '(4x,a,1(/,6x,a))') &
          'SPECIFIC STORAGE VALUES WILL BE CALCULATED USING THE CURRENT', &
          'EFFECTIVE STRESS'
      end if
      if (this%istrict_stress == 0) then
        write (this%iout, '(4x,a,1(/,6x,a))') &
          'SMALL OR NEGATIVE EFFECTIVE STRESS WILL BE REGULARIZED BY FLOORING', &
          'THE EFFECTIVE STRESS USED TO CALCULATE THE SPECIFIC STORAGE'
      else
        write (this%iout, '(4x,a,1(/,6x,a))') &
          'SMALL OR NEGATIVE EFFECTIVE STRESS WILL TERMINATE THE SIMULATION', &
          '(STRICT_EFFECTIVE_STRESS SPECIFIED)'
      end if
    else if (warn_estress_lag) then
      write (this%iout, '(4x,a,2(/,6x,a))') &
        'EFFECTIVE_STRESS_LAG HAS BEEN SPECIFIED BUT HAS NO EFFECT WHEN', &
        'USING THE HEAD-BASED FORMULATION (HEAD_BASED HAS BEEN SPECIFIED', &
        'IN THE OPTIONS BLOCK)'
    end if
    !
    write (this%iout, fmtoptr) 'GAMMAW        =', this%gammaw
    write (this%iout, fmtoptr) 'BETA          =', this%beta
    write (this%iout, fmtoptr) 'GAMMAW * BETA =', this%brg
    write (this%iout, '((1X,A))') 'END PACKAGE SETTINGS'
  end subroutine log_options

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_source_dimensions(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use GwfCsubInputModule, only: GwfCsubParamFoundType
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    ! -- local variables
    type(GwfCsubParamFoundType) :: found

    ! -- initialize dimensions to -1
    this%ninterbeds = -1

    ! -- update defaults from input context
    call mem_set_value(this%ninterbeds, 'NINTERBEDS', this%input_mempath, &
                       found%ninterbeds)
    call mem_set_value(this%maxsig0, 'MAXBOUND', this%input_mempath, &
                       found%maxbound)

    ! - log dimensions
    write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%packName))// &
      ' DIMENSIONS'
    write (this%iout, '(4x,a,i0)') 'NINTERBEDS = ', this%ninterbeds
    write (this%iout, '(4x,a,i0)') 'MAXSIG0 = ', this%maxsig0
    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%packName))//' DIMENSIONS'

    ! -- verify dimensions were set correctly
    if (.not. found%ninterbeds) then
      write (errmsg, '(a)') &
        'NINTERBEDS is a required dimension.'
      call store_error(errmsg)
      call store_error_filename(this%input_mempath)
    end if

    ! -- Call define_listlabel to construct the list label that is written
    !    when PRINT_INPUT option is used.
    call this%define_listlabel()
  end subroutine csub_source_dimensions

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_source_griddata(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use GwfCsubInputModule, only: GwfCsubParamFoundType
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    ! -- locals
    integer(I4B) :: node
    type(GwfCsubParamFoundType) :: found
    integer(I4B), dimension(:), pointer, contiguous :: map

    ! -- set map to convert user input data into reduced data
    map => null()
    if (this%dis%nodes < this%dis%nodesuser) map => this%dis%nodeuser

    ! -- update defaults from input context
    call mem_set_value(this%cg_ske_cr, 'CG_SKE_CR', this%input_mempath, &
                       map, found%cg_ske_cr)
    call mem_set_value(this%cg_thetaini, 'CG_THETA', this%input_mempath, &
                       map, found%cg_theta)
    call mem_set_value(this%sgm, 'SGM', this%input_mempath, map, found%sgm)
    call mem_set_value(this%sgs, 'SGS', this%input_mempath, map, found%sgs)

    ! -- cg_ske and cg_theta are required input params
    if (.not. found%cg_ske_cr) then
      call store_error('CG_SKE GRIDDATA must be specified.')
      call store_error_filename(this%input_fname)
    end if
    if (.not. found%cg_theta) then
      call store_error('CG_THETA GRIDDATA must be specified.')
      call store_error_filename(this%input_fname)
    end if

    ! -- if sgm and sgs have not been specified assign default values
    if (.not. found%sgm) then
      do node = 1, this%dis%nodes
        this%sgm(node) = 1.7D0
      end do
    end if
    if (.not. found%sgs) then
      do node = 1, this%dis%nodes
        this%sgs(node) = 2.0D0
      end do
    end if
  end subroutine csub_source_griddata

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_source_packagedata(this)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr, mem_allocate
    use MemoryManagerExtModule, only: mem_set_value, memorystore_release
    use CharacterStringModule, only: CharacterStringType
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    integer(I4B), dimension(:), pointer, contiguous :: icsubno
    integer(I4B), dimension(:, :), pointer, contiguous :: cellid_pkgdata
    integer(I4B), dimension(:), pointer :: cellid
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: cdelay
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: boundname
    real(DP), dimension(:), pointer, contiguous :: pcs, thick_frac, rnb
    real(DP), dimension(:), pointer, contiguous :: ssv_cc, sse_cr, theta, kv, h0
    character(len=LINELENGTH) :: cdelaystr
    character(len=LENBOUNDNAME) :: bndname
    character(len=20) :: cellidstr
    real(DP) :: top, botm, baq, q, thick, rval
    integer(I4B) :: idelay, ndelaybeds, csubno
    integer(I4B) :: ib, n, nodeu, noder

    ! -- set input context pointers
    call mem_setptr(icsubno, 'ICSUBNO', this%input_mempath)
    call mem_setptr(cellid_pkgdata, 'CELLID_PKGDATA', this%input_mempath)
    call mem_setptr(cdelay, 'CDELAY', this%input_mempath)
    call mem_setptr(pcs, 'PCS0', this%input_mempath)
    call mem_setptr(thick_frac, 'THICK_FRAC', this%input_mempath)
    call mem_setptr(rnb, 'RNB', this%input_mempath)
    call mem_setptr(ssv_cc, 'SSV_CC', this%input_mempath)
    call mem_setptr(sse_cr, 'SSE_CR', this%input_mempath)
    call mem_setptr(theta, 'THETA', this%input_mempath)
    call mem_setptr(kv, 'KV', this%input_mempath)
    call mem_setptr(h0, 'H0', this%input_mempath)
    call mem_setptr(boundname, 'BOUNDNAME', this%input_mempath)

    ! initialize ndelaybeds
    ndelaybeds = 0

    ! -- update state
    do n = 1, size(icsubno)

      ! -- set cubno
      csubno = icsubno(n)

      ! -- check csubno
      if (csubno < 1 .or. csubno > this%ninterbeds) then
        write (errmsg, '(a,1x,i0,2(1x,a),1x,i0,a)') &
          'Interbed number (', csubno, ') must be greater than 0 and ', &
          'less than or equal to', this%ninterbeds, '.'
        call store_error(errmsg)
        cycle
      end if

      ! -- set cellid
      cellid => cellid_pkgdata(:, n)

      ! -- set node user
      if (this%dis%ndim == 1) then
        nodeu = cellid(1)
      elseif (this%dis%ndim == 2) then
        nodeu = get_node(cellid(1), 1, cellid(2), &
                         this%dis%mshape(1), 1, &
                         this%dis%mshape(2))
      else
        nodeu = get_node(cellid(1), cellid(2), cellid(3), &
                         this%dis%mshape(1), &
                         this%dis%mshape(2), &
                         this%dis%mshape(3))
      end if

      ! -- set node reduced
      noder = this%dis%get_nodenumber(nodeu, 1)
      if (noder <= 0) then
        call this%dis%nodeu_to_string(nodeu, cellidstr)
        write (errmsg, '(a)') &
          'CSUB configured for inactive cell: '// &
          trim(adjustl(cellidstr))//'.'
        call store_error(errmsg)
        cycle
      end if

      ! -- update nodelists
      this%nodelist(csubno) = noder
      this%unodelist(csubno) = nodeu

      ! -- set top, botm, baq
      top = this%dis%top(noder)
      botm = this%dis%bot(noder)
      baq = top - botm

      ! -- set cdelay
      cdelaystr = cdelay(n)
      select case (cdelaystr)
      case ('NODELAY')
        idelay = 0
      case ('DELAY')
        ndelaybeds = ndelaybeds + 1
        idelay = ndelaybeds
      case default
        write (errmsg, '(a,1x,a,1x,i0,1x,a)') &
          'Invalid CDELAY ', trim(adjustl(cdelaystr)), &
          'for packagedata entry', csubno, '.'
        call store_error(errmsg)
        cycle
      end select
      this%idelay(csubno) = idelay

      ! -- set initial preconsolidation stress
      this%pcs(csubno) = pcs(n)

      ! -- set thickness
      if (this%icellf == 0) then
        if (thick_frac(n) < DZERO .or. thick_frac(n) > baq) then
          write (errmsg, '(a,g0,2(a,1x),g0,1x,a,1x,i0,a)') &
            'THICK (', thick_frac(n), ') MUST BE greater than or equal to 0 ', &
            'and less than or equal to than', baq, &
            'for packagedata entry', csubno, '.'
          call store_error(errmsg)
        end if
        thick = thick_frac(n)
      else
        if (thick_frac(n) < DZERO .or. thick_frac(n) > DONE) then
          write (errmsg, '(a,1x,a,1x,i0,a)') &
            'FRAC MUST BE greater than 0 and less than or equal to 1', &
            'for packagedata entry', csubno, '.'
          call store_error(errmsg)
        end if
        thick = thick_frac(n) * baq
      end if
      this%thickini(csubno) = thick
      if (this%iupdatematprop /= 0) then
        this%thick(csubno) = thick
      end if

      ! -- set rnb
      if (idelay > 0) then
        if (rnb(n) < DONE) then
          write (errmsg, '(a,g0,a,1x,a,1x,i0,a)') &
            'RNB (', rnb(n), ') must be greater than or equal to 1', &
            'for packagedata entry', csubno, '.'
          call store_error(errmsg)
        end if
        this%rnb(csubno) = rnb(n)
      else
        this%rnb(csubno) = DONE
      end if

      ! -- set skv or ci
      if (ssv_cc(n) < DZERO) then
        write (errmsg, '(2(a,1x),i0,a)') &
          '(SKV,CI) must be greater than or equal to 0', &
          'for packagedata entry', csubno, '.'
        call store_error(errmsg)
      end if
      this%ci(csubno) = ssv_cc(n)

      ! -- set ske or rci
      if (sse_cr(n) < DZERO) then
        write (errmsg, '(2(a,1x),i0,a)') &
          '(SKE,RCI) must be greater than or equal to 0', &
          'for packagedata entry', csubno, '.'
        call store_error(errmsg)
      end if
      this%rci(csubno) = sse_cr(n)

      ! -- set ielastic
      if (this%ci(csubno) == this%rci(csubno)) then
        this%ielastic(csubno) = 1
      else
        this%ielastic(csubno) = 0
      end if

      ! -- set porosity
      if (theta(n) <= DZERO .or. theta(n) > DONE) then
        write (errmsg, '(a,1x,a,1x,i0,a)') &
          'THETA must be greater than 0 and less than or equal to 1', &
          'for packagedata entry', csubno, '.'
        call store_error(errmsg)
      end if
      this%thetaini(csubno) = theta(n)
      if (this%iupdatematprop /= 0) then
        this%theta(csubno) = theta(n)
      end if

      ! -- set kv
      if (idelay > 0) then
        if (kv(n) <= 0.0) then
          write (errmsg, '(a,1x,i0,a)') &
            'KV must be greater than 0 for packagedata entry', csubno, '.'
          call store_error(errmsg)
        end if
      end if
      this%kv(csubno) = kv(n)

      ! -- set h0
      this%h0(csubno) = h0(n)

      ! -- set bound name
      if (this%inamedbound /= 0) then
        bndname = boundname(n)
        if (len_trim(bndname) < 1) then
          write (errmsg, '(a,1x,i0,a)') &
            'BOUNDNAME must be specified for packagedata entry', csubno, '.'
          call store_error(errmsg)
        end if
        this%boundname(csubno) = bndname
      end if
    end do

    !
    ! -- set the number of delay interbeds
    this%ndelaybeds = ndelaybeds
    !
    ! -- process delay interbeds
    if (ndelaybeds > 0) then
      !
      ! -- reallocate and initialize delay interbed arrays
      call mem_allocate(this%idb_nconv_count, 2, &
                        'IDB_NCONV_COUNT', trim(this%memoryPath))
      call mem_allocate(this%idbconvert, this%ndelaycells, ndelaybeds, &
                        'IDBCONVERT', trim(this%memoryPath))
      call mem_allocate(this%dbdhmax, ndelaybeds, &
                        'DBDHMAX', trim(this%memoryPath))
      call mem_allocate(this%dbz, this%ndelaycells, ndelaybeds, &
                        'DBZ', trim(this%memoryPath))
      call mem_allocate(this%dbrelz, this%ndelaycells, ndelaybeds, &
                        'DBRELZ', trim(this%memoryPath))
      call mem_allocate(this%dbh, this%ndelaycells, ndelaybeds, &
                        'DBH', trim(this%memoryPath))
      call mem_allocate(this%dbh0, this%ndelaycells, ndelaybeds, &
                        'DBH0', trim(this%memoryPath))
      call mem_allocate(this%dbgeo, this%ndelaycells, ndelaybeds, &
                        'DBGEO', trim(this%memoryPath))
      call mem_allocate(this%dbes, this%ndelaycells, ndelaybeds, &
                        'DBES', trim(this%memoryPath))
      call mem_allocate(this%dbes0, this%ndelaycells, ndelaybeds, &
                        'DBES0', trim(this%memoryPath))
      call mem_allocate(this%dbpcs, this%ndelaycells, ndelaybeds, &
                        'DBPCS', trim(this%memoryPath))
      call mem_allocate(this%dbflowtop, ndelaybeds, &
                        'DBFLOWTOP', trim(this%memoryPath))
      call mem_allocate(this%dbflowbot, ndelaybeds, &
                        'DBFLOWBOT', trim(this%memoryPath))
      call mem_allocate(this%dbdzini, this%ndelaycells, ndelaybeds, &
                        'DBDZINI', trim(this%memoryPath))
      call mem_allocate(this%dbthetaini, this%ndelaycells, ndelaybeds, &
                        'DBTHETAINI', trim(this%memoryPath))
      call mem_allocate(this%dbcomp, this%ndelaycells, ndelaybeds, &
                        'DBCOMP', trim(this%memoryPath))
      call mem_allocate(this%dbtcomp, this%ndelaycells, ndelaybeds, &
                        'DBTCOMP', trim(this%memoryPath))
      !
      ! -- allocate delay bed arrays
      if (this%iupdatematprop == 0) then
        call mem_setptr(this%dbdz, 'DBDZINI', trim(this%memoryPath))
        call mem_setptr(this%dbdz0, 'DBDZINI', trim(this%memoryPath))
        call mem_setptr(this%dbtheta, 'DBTHETAINI', trim(this%memoryPath))
        call mem_setptr(this%dbtheta0, 'DBTHETAINI', trim(this%memoryPath))
      else
        call mem_allocate(this%dbdz, this%ndelaycells, ndelaybeds, &
                          'DBDZ', trim(this%memoryPath))
        call mem_allocate(this%dbdz0, this%ndelaycells, ndelaybeds, &
                          'DBDZ0', trim(this%memoryPath))
        call mem_allocate(this%dbtheta, this%ndelaycells, ndelaybeds, &
                          'DBTHETA', trim(this%memoryPath))
        call mem_allocate(this%dbtheta0, this%ndelaycells, ndelaybeds, &
                          'DBTHETA0', trim(this%memoryPath))
      end if
      !
      ! -- allocate delay interbed solution arrays
      call mem_allocate(this%dbal, this%ndelaycells, &
                        'DBAL', trim(this%memoryPath))
      call mem_allocate(this%dbad, this%ndelaycells, &
                        'DBAD', trim(this%memoryPath))
      call mem_allocate(this%dbau, this%ndelaycells, &
                        'DBAU', trim(this%memoryPath))
      call mem_allocate(this%dbrhs, this%ndelaycells, &
                        'DBRHS', trim(this%memoryPath))
      call mem_allocate(this%dbdh, this%ndelaycells, &
                        'DBDH', trim(this%memoryPath))
      call mem_allocate(this%dbaw, this%ndelaycells, &
                        'DBAW', trim(this%memoryPath))
      !
      ! -- initialize delay bed counters
      do n = 1, 2
        this%idb_nconv_count(n) = 0
      end do
      !
      ! -- initialize delay bed storage
      do ib = 1, this%ninterbeds
        idelay = this%idelay(ib)
        if (idelay == 0) then
          cycle
        end if
        !
        ! -- initialize delay interbed variables
        do n = 1, this%ndelaycells
          rval = this%thickini(ib) / real(this%ndelaycells, DP)
          this%dbdzini(n, idelay) = rval
          this%dbh(n, idelay) = this%h0(ib)
          this%dbh0(n, idelay) = this%h0(ib)
          this%dbthetaini(n, idelay) = this%thetaini(ib)
          this%dbgeo(n, idelay) = DZERO
          this%dbes(n, idelay) = DZERO
          this%dbes0(n, idelay) = DZERO
          this%dbpcs(n, idelay) = this%pcs(ib)
          this%dbcomp(n, idelay) = DZERO
          this%dbtcomp(n, idelay) = DZERO
          if (this%iupdatematprop /= 0) then
            this%dbdz(n, idelay) = this%dbdzini(n, idelay)
            this%dbdz0(n, idelay) = this%dbdzini(n, idelay)
            this%dbtheta(n, idelay) = this%theta(ib)
            this%dbtheta0(n, idelay) = this%theta(ib)
          end if
        end do
        !
        ! -- initialize elevation of delay bed cells
        call this%csub_delay_init_zcell(ib)
      end do
      !
      ! -- initialize delay bed solution arrays
      do n = 1, this%ndelaycells
        this%dbal(n) = DZERO
        this%dbad(n) = DZERO
        this%dbau(n) = DZERO
        this%dbrhs(n) = DZERO
        this%dbdh(n) = DZERO
        this%dbaw(n) = DZERO
      end do
    end if
    !
    ! -- check that ndelaycells is odd when using
    !    the effective stress formulation
    if (ndelaybeds > 0) then
      q = MOD(real(this%ndelaycells, DP), DTWO)
      if (q == DZERO) then
        write (errmsg, '(a,i0,a,1x,a)') &
          'NDELAYCELLS (', this%ndelaycells, ') must be an', &
          'odd number when using the effective stress formulation.'
        call store_error(errmsg)
      end if
    end if

    if (this%iprpak /= 0) then
      call this%csub_print_packagedata()
    end if

    ! -- terminate if errors encountered
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if

    call memorystore_release('ICSUBNO', this%input_mempath)
    call memorystore_release('CELLID_PKGDATA', this%input_mempath)
    call memorystore_release('CDELAY', this%input_mempath)
    call memorystore_release('PCS0', this%input_mempath)
    call memorystore_release('THICK_FRAC', this%input_mempath)
    call memorystore_release('RNB', this%input_mempath)
    call memorystore_release('SSV_CC', this%input_mempath)
    call memorystore_release('SSE_CR', this%input_mempath)
    call memorystore_release('THETA', this%input_mempath)
    call memorystore_release('KV', this%input_mempath)
    call memorystore_release('H0', this%input_mempath)
    call memorystore_release('BOUNDNAME', this%input_mempath)
  end subroutine csub_source_packagedata

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_print_packagedata(this)
    ! -- dummy variables
    class(GwfCsubType) :: this
    ! local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: tag
    character(len=10) :: ctype
    character(len=20) :: cellid
    integer(I4B) :: ntabrows
    integer(I4B) :: ntabcols
    integer(I4B) :: ib
    integer(I4b) :: idelay
    integer(I4B) :: node

    ! set title
    title = 'CSUB'//' PACKAGE ('// &
            trim(adjustl(this%packName))//') INTERBED DATA'
    !
    ! determine the number of columns and rows
    ntabrows = this%ninterbeds
    ntabcols = 13
    if (this%inamedbound /= 0) then
      ntabcols = ntabcols + 1
    end if

    ! setup table
    call table_cr(this%inputtab, this%packName, title)
    call this%inputtab%table_df(ntabrows, ntabcols, this%iout)
    !
    ! add columns
    !<icsubno> <cellid> <cdelay> <pcs0> <thick_frac> <rnb> <ssv_cc> <sse_cr> <theta> <kv> <h0> [<boundname>]

    tag = 'INTERBED NUMBER'
    call this%inputtab%initialize_column(tag, 10, alignment=TABCENTER)
    tag = 'CELLID'
    call this%inputtab%initialize_column(tag, 20, alignment=TABLEFT)
    tag = 'INTERBED TYPE'
    call this%inputtab%initialize_column(tag, 10, alignment=TABCENTER)
    tag = 'PCS0'
    call this%inputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'THICK_FRAC'
    call this%inputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'RNB'
    call this%inputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'INTERBED THICKNESS'
    call this%inputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'CELL THICKNESS'
    call this%inputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'SSV_CV'
    call this%inputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'SSE_CR'
    call this%inputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'THETA'
    call this%inputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'KV'
    call this%inputtab%initialize_column(tag, 12, alignment=TABCENTER)
    tag = 'H0'
    call this%inputtab%initialize_column(tag, 12, alignment=TABCENTER)
    if (this%inamedbound /= 0) then
      tag = 'BOUNDNAME'
      call this%inputtab%initialize_column(tag, 40, alignment=TABLEFT)
    end if

    do ib = 1, this%ninterbeds
      idelay = this%idelay(ib)
      node = this%nodelist(ib)
      call this%dis%noder_to_string(node, cellid)
      if (idelay == 0) then
        ctype = 'nodelay'
      else
        ctype = 'delay'
      end if

      ! fill table line
      call this%inputtab%add_term(ib)
      call this%inputtab%add_term(cellid)
      call this%inputtab%add_term(ctype)
      call this%inputtab%add_term(this%pcs(ib))
      call this%inputtab%add_term(this%thickini(ib))
      call this%inputtab%add_term(this%rnb(ib))
      call this%inputtab%add_term(this%thickini(ib) * this%rnb(ib))
      call this%inputtab%add_term(this%dis%top(node) - this%dis%bot(node))
      call this%inputtab%add_term(this%ci(ib))
      call this%inputtab%add_term(this%rci(ib))
      call this%inputtab%add_term(this%theta(ib))
      if (idelay == 0) then
        call this%inputtab%add_term("--")
        call this%inputtab%add_term("--")
      else
        call this%inputtab%add_term(this%kv(ib))
        call this%inputtab%add_term(this%h0(ib))
      end if
      if (this%inamedbound /= 0) then
        call this%inputtab%add_term(this%boundname(ib))
      end if
    end do

  end subroutine csub_print_packagedata

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_initialize_tables(this)
    ! -- dummy variables
    class(GwfCsubType) :: this

    character(len=LINELENGTH) :: tag
    integer(I4B) :: ntabrows
    integer(I4B) :: ntabcols

    if (this%ipakcsv > 0) then
      if (this%ndelaybeds < 1) then
        write (warnmsg, '(a,1x,3a)') &
          'Package convergence data is requested but delay interbeds', &
          'are not included in package (', &
          trim(adjustl(this%packName)), ').'
        call store_warning(warnmsg)
      end if

      ntabrows = 1
      ntabcols = 9

      ! setup table
      call table_cr(this%pakcsvtab, this%packName, '')
      call this%pakcsvtab%table_df(ntabrows, ntabcols, this%ipakcsv, &
                                   lineseparator=.FALSE., separator=',', &
                                   finalize=.FALSE.)

      ! add columns to package csv
      tag = 'total_inner_iterations'
      call this%pakcsvtab%initialize_column(tag, 10, alignment=TABLEFT)
      tag = 'totim'
      call this%pakcsvtab%initialize_column(tag, 10, alignment=TABLEFT)
      tag = 'kper'
      call this%pakcsvtab%initialize_column(tag, 10, alignment=TABLEFT)
      tag = 'kstp'
      call this%pakcsvtab%initialize_column(tag, 10, alignment=TABLEFT)
      tag = 'nouter'
      call this%pakcsvtab%initialize_column(tag, 10, alignment=TABLEFT)
      tag = 'dvmax'
      call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
      tag = 'dvmax_loc'
      call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
      tag = 'dstoragemax'
      call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
      tag = 'dstoragemax_loc'
      call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
    end if

  end subroutine csub_initialize_tables

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine define_listlabel(this)
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    !
    ! -- create the header list label
    this%listlabel = trim(this%filtyp)//' NO.'
    if (this%dis%ndim == 3) then
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'LAYER'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'ROW'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'COL'
    elseif (this%dis%ndim == 2) then
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'LAYER'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'CELL2D'
    else
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'NODE'
    end if
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'SIG0'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
  end subroutine define_listlabel

end submodule GwfCsubModuleInput
