!> @brief CSUB package memory allocation and deallocation
!<
submodule(GwfCsubModule) GwfCsubModuleMemory
contains

  !> @brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the CSUB package. The base model
  !! allocate scalars method is also called.
  !!
  !<
  module procedure csub_allocate_scalars
  ! -- modules
  use MemoryManagerModule, only: mem_allocate
  !
  ! -- call standard NumericalPackageType allocate scalars
  call this%NumericalPackageType%allocate_scalars()
  !
  ! -- allocate character variables
  call mem_allocate(this%listlabel, LENLISTLABEL, 'LISTLABEL', this%memoryPath)
  call mem_allocate(this%stoMemPath, LENMEMPATH, 'STONAME', this%memoryPath)
  !
  ! -- allocate the object and assign values to object variables
  call mem_allocate(this%istounit, 'ISTOUNIT', this%memoryPath)
  call mem_allocate(this%inobspkg, 'INOBSPKG', this%memoryPath)
  call mem_allocate(this%ninterbeds, 'NINTERBEDS', this%memoryPath)
  call mem_allocate(this%maxsig0, 'MAXSIG0', this%memoryPath)
  call mem_allocate(this%nbound, 'NBOUND', this%memoryPath)
  call mem_allocate(this%iscloc, 'ISCLOC', this%memoryPath)
  call mem_allocate(this%iauxmultcol, 'IAUXMULTCOL', this%memoryPath)
  call mem_allocate(this%ndelaycells, 'NDELAYCELLS', this%memoryPath)
  call mem_allocate(this%ndelaybeds, 'NDELAYBEDS', this%memoryPath)
  call mem_allocate(this%initialized, 'INITIALIZED', this%memoryPath)
  call mem_allocate(this%ieslag, 'IESLAG', this%memoryPath)
  call mem_allocate(this%ipch, 'IPCH', this%memoryPath)
  call mem_allocate(this%istrict_stress, 'ISTRICT_STRESS', this%memoryPath)
  call mem_allocate(this%nreg_ts, 'NREG_TS', this%memoryPath)
  call mem_allocate(this%lhead_based, 'LHEAD_BASED', this%memoryPath)
  call mem_allocate(this%iupdatestress, 'IUPDATESTRESS', this%memoryPath)
  call mem_allocate(this%ispecified_pcs, 'ISPECIFIED_PCS', this%memoryPath)
  call mem_allocate(this%ispecified_dbh, 'ISPECIFIED_DBH', this%memoryPath)
  call mem_allocate(this%inamedbound, 'INAMEDBOUND', this%memoryPath)
  call mem_allocate(this%iconvchk, 'ICONVCHK', this%memoryPath)
  call mem_allocate(this%naux, 'NAUX', this%memoryPath)
  call mem_allocate(this%istoragec, 'ISTORAGEC', this%memoryPath)
  call mem_allocate(this%istrainib, 'ISTRAINIB', this%memoryPath)
  call mem_allocate(this%istrainsk, 'ISTRAINSK', this%memoryPath)
  call mem_allocate(this%ioutcomp, 'IOUTCOMP', this%memoryPath)
  call mem_allocate(this%ioutcompi, 'IOUTCOMPI', this%memoryPath)
  call mem_allocate(this%ioutcompe, 'IOUTCOMPE', this%memoryPath)
  call mem_allocate(this%ioutcompib, 'IOUTCOMPIB', this%memoryPath)
  call mem_allocate(this%ioutcomps, 'IOUTCOMPS', this%memoryPath)
  call mem_allocate(this%ioutzdisp, 'IOUTZDISP', this%memoryPath)
  call mem_allocate(this%ipakcsv, 'IPAKCSV', this%memoryPath)
  call mem_allocate(this%iupdatematprop, 'IUPDATEMATPROP', this%memoryPath)
  call mem_allocate(this%epsilon, 'EPSILON', this%memoryPath)
  call mem_allocate(this%cc_crit, 'CC_CRIT', this%memoryPath)
  call mem_allocate(this%gammaw, 'GAMMAW', this%memoryPath)
  call mem_allocate(this%beta, 'BETA', this%memoryPath)
  call mem_allocate(this%brg, 'BRG', this%memoryPath)
  call mem_allocate(this%satomega, 'SATOMEGA', this%memoryPath)
  call mem_allocate(this%pcsomega, 'PCSOMEGA', this%memoryPath)
  call mem_allocate(this%icellf, 'ICELLF', this%memoryPath)
  call mem_allocate(this%gwfiss0, 'GWFISS0', this%memoryPath)
  !
  ! -- allocate text strings
  call mem_allocate(this%auxname, LENAUXNAME, 0, 'AUXNAME', this%memoryPath)
  !
  ! -- initialize values
  this%istounit = 0
  this%inobspkg = 0
  this%ninterbeds = 0
  this%maxsig0 = 0
  this%nbound = 0
  this%iscloc = 0
  this%iauxmultcol = 0
  this%ndelaycells = 19
  this%ndelaybeds = 0
  this%initialized = 0
  this%ieslag = 0
  this%ipch = 0
  this%istrict_stress = 0
  this%nreg_ts = 0
  this%lhead_based = .FALSE.
  this%iupdatestress = 1
  this%ispecified_pcs = 0
  this%ispecified_dbh = 0
  this%inamedbound = 0
  this%iconvchk = 1
  this%naux = 0
  this%istoragec = 1
  this%istrainib = 0
  this%istrainsk = 0
  this%ioutcomp = 0
  this%ioutcompi = 0
  this%ioutcompe = 0
  this%ioutcompib = 0
  this%ioutcomps = 0
  this%ioutzdisp = 0
  this%ipakcsv = 0
  this%iupdatematprop = 0
  this%epsilon = DZERO
  this%cc_crit = DEM7
  this%gammaw = DGRAVITY * 1000._DP
  this%beta = 4.6512e-10_DP
  this%brg = this%gammaw * this%beta
  ! -- fraction of pcs over which ssk blends elastic->inelastic; 0 = hard switch
  this%pcsomega = DZERO
  !
  ! -- set omega value used for saturation calculations
  if (this%inewton /= 0) then
    this%satomega = DEM6
    this%epsilon = DHALF * DEM6
  else
    this%satomega = DZERO
  end if
  this%icellf = 0
  this%ninterbeds = 0
  this%gwfiss0 = 0
  end procedure csub_allocate_scalars

  !> @brief Allocate package arrays
  !!
  !!  Allocate and initialize CSUB package arrays.
  !!
  !<
  module procedure csub_allocate_arrays
  ! -- modules
  use MemoryManagerModule, only: mem_allocate, mem_setptr, mem_checkin
  ! -- local variables
  integer(I4B) :: j
  integer(I4B) :: n
  integer(I4B) :: iblen
  integer(I4B) :: naux
  !
  ! -- grid based data
  if (this%ioutcomp == 0 .and. this%ioutcompi == 0 .and. &
      this%ioutcompe == 0 .and. this%ioutcompib == 0 .and. &
      this%ioutcomps == 0 .and. this%ioutzdisp == 0) then
    call mem_allocate(this%buff, 1, 'BUFF', trim(this%memoryPath))
  else
    call mem_allocate(this%buff, this%dis%nodes, 'BUFF', trim(this%memoryPath))
  end if
  if (this%ioutcomp == 0 .and. this%ioutzdisp == 0) then
    call mem_allocate(this%buffusr, 1, 'BUFFUSR', trim(this%memoryPath))
  else
    call mem_allocate(this%buffusr, this%dis%nodesuser, 'BUFFUSR', &
                      trim(this%memoryPath))
  end if
  call mem_allocate(this%sgm, this%dis%nodes, 'SGM', trim(this%memoryPath))
  call mem_allocate(this%sgs, this%dis%nodes, 'SGS', trim(this%memoryPath))
  call mem_allocate(this%cg_ske_cr, this%dis%nodes, 'CG_SKE_CR', &
                    trim(this%memoryPath))
  call mem_allocate(this%cg_es, this%dis%nodes, 'CG_ES', &
                    trim(this%memoryPath))
  call mem_allocate(this%cg_es0, this%dis%nodes, 'CG_ES0', &
                    trim(this%memoryPath))
  call mem_allocate(this%cg_pcs, this%dis%nodes, 'CG_PCS', &
                    trim(this%memoryPath))
  call mem_allocate(this%cg_comp, this%dis%nodes, 'CG_COMP', &
                    trim(this%memoryPath))
  call mem_allocate(this%cg_tcomp, this%dis%nodes, 'CG_TCOMP', &
                    trim(this%memoryPath))
  call mem_allocate(this%cg_stor, this%dis%nodes, 'CG_STOR', &
                    trim(this%memoryPath))
  call mem_allocate(this%cg_ske, this%dis%nodes, 'CG_SKE', &
                    trim(this%memoryPath))
  call mem_allocate(this%cg_sk, this%dis%nodes, 'CG_SK', &
                    trim(this%memoryPath))
  call mem_allocate(this%cg_thickini, this%dis%nodes, 'CG_THICKINI', &
                    trim(this%memoryPath))
  call mem_allocate(this%cg_thetaini, this%dis%nodes, 'CG_THETAINI', &
                    trim(this%memoryPath))
  if (this%iupdatematprop == 0) then
    call mem_setptr(this%cg_thick, 'CG_THICKINI', trim(this%memoryPath))
    call mem_setptr(this%cg_thick0, 'CG_THICKINI', trim(this%memoryPath))
    call mem_setptr(this%cg_theta, 'CG_THETAINI', trim(this%memoryPath))
    call mem_setptr(this%cg_theta0, 'CG_THETAINI', trim(this%memoryPath))
  else
    call mem_allocate(this%cg_thick, this%dis%nodes, 'CG_THICK', &
                      trim(this%memoryPath))
    call mem_allocate(this%cg_thick0, this%dis%nodes, 'CG_THICK0', &
                      trim(this%memoryPath))
    call mem_allocate(this%cg_theta, this%dis%nodes, 'CG_THETA', &
                      trim(this%memoryPath))
    call mem_allocate(this%cg_theta0, this%dis%nodes, 'CG_THETA0', &
                      trim(this%memoryPath))
  end if
  !
  ! -- cell storage data
  call mem_allocate(this%cell_wcstor, this%dis%nodes, 'CELL_WCSTOR', &
                    trim(this%memoryPath))
  call mem_allocate(this%cell_thick, this%dis%nodes, 'CELL_THICK', &
                    trim(this%memoryPath))
  !
  ! -- interbed data
  iblen = 1
  if (this%ninterbeds > 0) then
    iblen = this%ninterbeds
  end if
  naux = 1
  if (this%naux > 0) then
    naux = this%naux
  end if
  call mem_allocate(this%auxvar, naux, iblen, 'AUXVAR', this%memoryPath)
  do n = 1, iblen
    do j = 1, naux
      this%auxvar(j, n) = DZERO
    end do
  end do
  call mem_allocate(this%unodelist, iblen, 'UNODELIST', trim(this%memoryPath))
  call mem_allocate(this%nodelist, iblen, 'NODELIST', trim(this%memoryPath))
  call mem_allocate(this%cg_gs, this%dis%nodes, 'CG_GS', trim(this%memoryPath))
  call mem_allocate(this%pcs, iblen, 'PCS', trim(this%memoryPath))
  call mem_allocate(this%rnb, iblen, 'RNB', trim(this%memoryPath))
  call mem_allocate(this%kv, iblen, 'KV', trim(this%memoryPath))
  call mem_allocate(this%h0, iblen, 'H0', trim(this%memoryPath))
  call mem_allocate(this%ci, iblen, 'CI', trim(this%memoryPath))
  call mem_allocate(this%rci, iblen, 'RCI', trim(this%memoryPath))
  call mem_allocate(this%idelay, iblen, 'IDELAY', trim(this%memoryPath))
  call mem_allocate(this%ielastic, iblen, 'IELASTIC', trim(this%memoryPath))
  call mem_allocate(this%iconvert, iblen, 'ICONVERT', trim(this%memoryPath))
  call mem_allocate(this%comp, iblen, 'COMP', trim(this%memoryPath))
  call mem_allocate(this%tcomp, iblen, 'TCOMP', trim(this%memoryPath))
  call mem_allocate(this%tcompi, iblen, 'TCOMPI', trim(this%memoryPath))
  call mem_allocate(this%tcompe, iblen, 'TCOMPE', trim(this%memoryPath))
  call mem_allocate(this%storagee, iblen, 'STORAGEE', trim(this%memoryPath))
  call mem_allocate(this%storagei, iblen, 'STORAGEI', trim(this%memoryPath))
  call mem_allocate(this%ske, iblen, 'SKE', trim(this%memoryPath))
  call mem_allocate(this%sk, iblen, 'SK', trim(this%memoryPath))
  call mem_allocate(this%thickini, iblen, 'THICKINI', trim(this%memoryPath))
  call mem_allocate(this%thetaini, iblen, 'THETAINI', trim(this%memoryPath))
  if (this%iupdatematprop == 0) then
    call mem_setptr(this%thick, 'THICKINI', trim(this%memoryPath))
    call mem_setptr(this%thick0, 'THICKINI', trim(this%memoryPath))
    call mem_setptr(this%theta, 'THETAINI', trim(this%memoryPath))
    call mem_setptr(this%theta0, 'THETAINI', trim(this%memoryPath))
  else
    call mem_allocate(this%thick, iblen, 'THICK', trim(this%memoryPath))
    call mem_allocate(this%thick0, iblen, 'THICK0', trim(this%memoryPath))
    call mem_allocate(this%theta, iblen, 'THETA', trim(this%memoryPath))
    call mem_allocate(this%theta0, iblen, 'THETA0', trim(this%memoryPath))
  end if
  !
  ! -- delay bed storage - allocated in csub_source_packagedata
  !    after number of delay beds is defined
  !
  ! -- allocate boundname
  if (this%inamedbound /= 0) then
    call mem_allocate(this%boundname, LENBOUNDNAME, this%ninterbeds, &
                      'BOUNDNAME', trim(this%memoryPath))
  else
    call mem_allocate(this%boundname, LENBOUNDNAME, 1, &
                      'BOUNDNAME', trim(this%memoryPath))

  end if
  !
  ! -- allocate the nodelist and bound arrays
  call mem_allocate(this%nodelistsig0, this%maxsig0, 'NODELISTSIG0', &
                    this%memoryPath)

  ! -- set sig0 input context pointer
  call mem_setptr(this%sig0, 'SIG0', this%input_mempath)
  call mem_checkin(this%sig0, 'SIG0', this%memoryPath, &
                   'SIG0', this%input_mempath)
  !
  ! -- set pointers to gwf variables
  call mem_setptr(this%gwfiss, 'ISS', trim(this%name_model))
  !
  ! -- set pointers to variables in the storage package
  call mem_setptr(this%stoiconv, 'ICONVERT', this%stoMemPath)
  call mem_setptr(this%stoss, 'SS', this%stoMemPath)
  !
  ! -- initialize variables that are not specified by user
  do n = 1, this%dis%nodes
    this%cg_gs(n) = DZERO
    this%cg_es(n) = DZERO
    this%cg_comp(n) = DZERO
    this%cg_tcomp(n) = DZERO
    this%cell_wcstor(n) = DZERO
  end do
  do n = 1, this%ninterbeds
    this%theta(n) = DZERO
    this%tcomp(n) = DZERO
    this%tcompi(n) = DZERO
    this%tcompe(n) = DZERO
  end do
  do n = 1, this%maxsig0
    this%nodelistsig0(n) = 0
  end do
  end procedure csub_allocate_arrays

  !> @brief Deallocate package memory
  !!
  !!  Deallocate CSUB package scalars and arrays.
  !!
  !<
  module procedure csub_da
  ! -- modules
  use MemoryManagerModule, only: mem_deallocate
  !
  ! -- summarize effective-stress regularization for the run
  if (this%nreg_ts > 0) then
    write (warnmsg, '(a,1x,i0,1x,3a)') &
      'CSUB negative effective stress was regularized in', this%nreg_ts, &
      'time step(s); see the model listing file for the number of cells ', &
      'regularized in each time step. This typically occurs in uppermost ', &
      'cells where simulated water levels rise above land surface.'
    call store_warning(warnmsg)
  end if
  !
  ! -- Deallocate arrays if package is active
  if (this%inunit > 0) then
    call mem_deallocate(this%unodelist)
    call mem_deallocate(this%nodelist)
    call mem_deallocate(this%idelay)
    call mem_deallocate(this%ielastic)
    call mem_deallocate(this%iconvert)
    !
    ! -- grid-based storage data
    call mem_deallocate(this%buff)
    call mem_deallocate(this%buffusr)
    call mem_deallocate(this%sgm)
    call mem_deallocate(this%sgs)
    call mem_deallocate(this%cg_ske_cr)
    call mem_deallocate(this%cg_gs)
    call mem_deallocate(this%cg_es)
    call mem_deallocate(this%cg_es0)
    call mem_deallocate(this%cg_pcs)
    call mem_deallocate(this%cg_comp)
    call mem_deallocate(this%cg_tcomp)
    call mem_deallocate(this%cg_stor)
    call mem_deallocate(this%cg_ske)
    call mem_deallocate(this%cg_sk)
    if (this%iupdatematprop == 0) then
      nullify (this%cg_thick)
      nullify (this%cg_thick0)
      nullify (this%cg_theta)
      nullify (this%cg_theta0)
    else
      call mem_deallocate(this%cg_thick)
      call mem_deallocate(this%cg_thick0)
      call mem_deallocate(this%cg_theta)
      call mem_deallocate(this%cg_theta0)
    end if
    call mem_deallocate(this%cg_thickini)
    call mem_deallocate(this%cg_thetaini)
    !
    ! -- cell storage
    call mem_deallocate(this%cell_wcstor)
    call mem_deallocate(this%cell_thick)
    !
    ! -- interbed storage
    call mem_deallocate(this%boundname, 'BOUNDNAME', this%memoryPath)
    call mem_deallocate(this%auxname, 'AUXNAME', this%memoryPath)
    call mem_deallocate(this%auxvar)
    call mem_deallocate(this%ci)
    call mem_deallocate(this%rci)
    call mem_deallocate(this%pcs)
    call mem_deallocate(this%rnb)
    call mem_deallocate(this%kv)
    call mem_deallocate(this%h0)
    call mem_deallocate(this%comp)
    call mem_deallocate(this%tcomp)
    call mem_deallocate(this%tcompi)
    call mem_deallocate(this%tcompe)
    call mem_deallocate(this%storagee)
    call mem_deallocate(this%storagei)
    call mem_deallocate(this%ske)
    call mem_deallocate(this%sk)
    if (this%iupdatematprop == 0) then
      nullify (this%thick)
      nullify (this%thick0)
      nullify (this%theta)
      nullify (this%theta0)
    else
      call mem_deallocate(this%thick)
      call mem_deallocate(this%thick0)
      call mem_deallocate(this%theta)
      call mem_deallocate(this%theta0)
    end if
    call mem_deallocate(this%thickini)
    call mem_deallocate(this%thetaini)
    !
    ! -- delay bed storage
    if (this%ndelaybeds > 0) then
      if (this%iupdatematprop == 0) then
        nullify (this%dbdz)
        nullify (this%dbdz0)
        nullify (this%dbtheta)
        nullify (this%dbtheta0)
      else
        call mem_deallocate(this%dbdz)
        call mem_deallocate(this%dbdz0)
        call mem_deallocate(this%dbtheta)
        call mem_deallocate(this%dbtheta0)
      end if
      call mem_deallocate(this%idb_nconv_count)
      call mem_deallocate(this%idbconvert)
      call mem_deallocate(this%dbdhmax)
      call mem_deallocate(this%dbz)
      call mem_deallocate(this%dbrelz)
      call mem_deallocate(this%dbh)
      call mem_deallocate(this%dbh0)
      call mem_deallocate(this%dbgeo)
      call mem_deallocate(this%dbes)
      call mem_deallocate(this%dbes0)
      call mem_deallocate(this%dbpcs)
      call mem_deallocate(this%dbflowtop)
      call mem_deallocate(this%dbflowbot)
      call mem_deallocate(this%dbdzini)
      call mem_deallocate(this%dbthetaini)
      call mem_deallocate(this%dbcomp)
      call mem_deallocate(this%dbtcomp)
      !
      ! -- delay interbed solution arrays
      call mem_deallocate(this%dbal)
      call mem_deallocate(this%dbad)
      call mem_deallocate(this%dbau)
      call mem_deallocate(this%dbrhs)
      call mem_deallocate(this%dbdh)
      call mem_deallocate(this%dbaw)
    end if
    !
    ! -- period data
    call mem_deallocate(this%nodelistsig0)
    call mem_deallocate(this%sig0, 'SIG0', this%memoryPath)
    !
    ! -- pointers to gwf variables
    nullify (this%gwfiss)
    !
    ! -- pointers to storage variables
    nullify (this%stoiconv)
    nullify (this%stoss)
    !
    ! -- input table
    if (this%iprpak > 0) then
      call this%inputtab%table_da()
      deallocate (this%inputtab)
      nullify (this%inputtab)
    end if
    !
    ! -- output table
    if (associated(this%outputtab)) then
      call this%outputtab%table_da()
      deallocate (this%outputtab)
      nullify (this%outputtab)
    end if
  end if
  !
  ! -- package csv table
  if (this%ipakcsv > 0) then
    call this%pakcsvtab%table_da()
    deallocate (this%pakcsvtab)
    nullify (this%pakcsvtab)
  end if
  !
  ! -- deallocate character variables
  call mem_deallocate(this%listlabel, 'LISTLABEL', this%memoryPath)
  call mem_deallocate(this%stoMemPath, 'STONAME', this%memoryPath)
  !
  ! -- deallocate scalars
  call mem_deallocate(this%istounit)
  call mem_deallocate(this%inobspkg)
  call mem_deallocate(this%ninterbeds)
  call mem_deallocate(this%maxsig0)
  call mem_deallocate(this%nbound)
  call mem_deallocate(this%iscloc)
  call mem_deallocate(this%iauxmultcol)
  call mem_deallocate(this%ndelaycells)
  call mem_deallocate(this%ndelaybeds)
  call mem_deallocate(this%initialized)
  call mem_deallocate(this%ieslag)
  call mem_deallocate(this%ipch)
  call mem_deallocate(this%lhead_based)
  call mem_deallocate(this%iupdatestress)
  call mem_deallocate(this%ispecified_pcs)
  call mem_deallocate(this%ispecified_dbh)
  call mem_deallocate(this%inamedbound)
  call mem_deallocate(this%iconvchk)
  call mem_deallocate(this%naux)
  call mem_deallocate(this%istoragec)
  call mem_deallocate(this%istrainib)
  call mem_deallocate(this%istrainsk)
  call mem_deallocate(this%ioutcomp)
  call mem_deallocate(this%ioutcompi)
  call mem_deallocate(this%ioutcompe)
  call mem_deallocate(this%ioutcompib)
  call mem_deallocate(this%ioutcomps)
  call mem_deallocate(this%ioutzdisp)
  call mem_deallocate(this%ipakcsv)
  call mem_deallocate(this%iupdatematprop)
  call mem_deallocate(this%epsilon)
  call mem_deallocate(this%cc_crit)
  call mem_deallocate(this%gammaw)
  call mem_deallocate(this%beta)
  call mem_deallocate(this%brg)
  call mem_deallocate(this%satomega)
  call mem_deallocate(this%pcsomega)
  call mem_deallocate(this%istrict_stress)
  call mem_deallocate(this%nreg_ts)
  call mem_deallocate(this%icellf)
  call mem_deallocate(this%gwfiss0)
  !
  ! -- deallocate methods on objects
  if (this%inunit > 0) then
    call this%obs%obs_da()
    !
    ! -- deallocate and nullify observations
    deallocate (this%obs)
    nullify (this%obs)
  end if

  !
  ! -- deallocate parent
  call this%NumericalPackageType%da()
  end procedure csub_da

end submodule GwfCsubModuleMemory
