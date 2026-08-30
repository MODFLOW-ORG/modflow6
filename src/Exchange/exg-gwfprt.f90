module GwfPrtExchangeModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENPACKAGENAME
  use ListsModule, only: basemodellist, baseexchangelist
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use BaseModelModule, only: BaseModelType, GetBaseModelFromList
  use GwfModule, only: GwfModelType
  use PrtModule, only: PrtModelType
  use BndModule, only: BndType, GetBndFromList

  implicit none
  public :: GwfPrtExchangeType
  public :: gwfprt_cr

  type, extends(BaseExchangeType) :: GwfPrtExchangeType

    integer(I4B), pointer :: m1id => null()
    integer(I4B), pointer :: m2id => null()

  contains

    procedure :: exg_df
    procedure :: exg_ar
    procedure :: exg_da
    procedure, private :: set_model_pointers
    procedure, private :: allocate_scalars
    procedure, private :: gwfbnd2prtfmi
    ! procedure, private :: gwfconn2prtconn
    ! procedure, private :: link_connections

  end type GwfPrtExchangeType

contains

  !> @brief Create a new GWF to PRT exchange object
  subroutine gwfprt_cr(filename, id, m1id, m2id)
    ! -- modules
    use SimVariablesModule, only: model_loc_idx
    ! -- dummy
    character(len=*), intent(in) :: filename
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: m1id
    integer(I4B), intent(in) :: m2id
    ! -- local
    class(BaseExchangeType), pointer :: baseexchange => null()
    type(GwfPrtExchangeType), pointer :: exchange => null()
    character(len=20) :: cint
    !
    ! -- Create a new exchange and add it to the baseexchangelist container
    allocate (exchange)
    baseexchange => exchange
    call AddBaseExchangeToList(baseexchangelist, baseexchange)
    !
    ! -- Assign id and name
    exchange%id = id
    write (cint, '(i0)') id
    exchange%name = 'GWF-PRT_'//trim(adjustl(cint))
    exchange%memoryPath = exchange%name
    !
    ! -- allocate scalars
    call exchange%allocate_scalars()
    !
    ! -- NB: convert from id to local model index in base model list
    exchange%m1id = model_loc_idx(m1id)
    exchange%m2id = model_loc_idx(m2id)
    !
    ! -- set model pointers
    call exchange%set_model_pointers()
  end subroutine gwfprt_cr

  subroutine set_model_pointers(this)
    ! -- modules
    ! -- dummy
    class(GwfPrtExchangeType) :: this
    ! -- local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(PrtModelType), pointer :: prtmodel => null()
    !
    ! -- set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1id)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set prtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2id)
    select type (mb)
    type is (PrtModelType)
      prtmodel => mb
    end select
    !
    ! -- Verify that GWF model is of the correct type
    if (.not. associated(gwfmodel)) then
      write (errmsg, '(3a)') 'Problem with GWF-PRT exchange ', trim(this%name), &
        '.  Specified GWF Model does not appear to be of the correct type.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Verify that PRT model is of the correct type
    if (.not. associated(prtmodel)) then
      write (errmsg, '(3a)') 'Problem with GWF-PRT exchange ', trim(this%name), &
        '.  Specified PRT Model does not appear to be of the correct type.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Tell particle tracking model fmi flows are not read from file
    prtmodel%fmi%flows_from_file = .false.
    !
    ! -- Set a pointer to the GWF bndlist.  This will allow the transport model
    !    to look through the flow packages and establish a link to GWF flows
    prtmodel%fmi%gwfbndlist => gwfmodel%bndlist
  end subroutine set_model_pointers

  subroutine exg_df(this)
    ! -- modules
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(GwfPrtExchangeType) :: this
    ! -- local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(PrtModelType), pointer :: prtmodel => null()
    integer(I4B) :: ngwfpack, ip
    class(BndType), pointer :: packobj => null()
    !
    !
    ! -- set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1id)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set prtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2id)
    select type (mb)
    type is (PrtModelType)
      prtmodel => mb
    end select
    !
    ! -- Check to make sure that flow is solved before particle tracking and in a
    !    different solution
    if (gwfmodel%idsoln >= prtmodel%idsoln) then
      write (errmsg, '(3a)') 'Problem with GWF-PRT exchange ', trim(this%name), &
        '.  The GWF model must be solved by a different solution than the PRT model. &
        &The IMS specified for GWF must be listed in mfsim.nam &
        &before the EMS for PRT.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Set pointer to flowja
    prtmodel%fmi%gwfflowja => gwfmodel%flowja
    call mem_checkin(prtmodel%fmi%gwfflowja, &
                     'GWFFLOWJA', prtmodel%fmi%memoryPath, &
                     'FLOWJA', gwfmodel%memoryPath)
    !
    ! -- Set the npf flag so that specific discharge is available for
    !    transport calculations if dispersion is active
    if (prtmodel%indsp > 0) then
      gwfmodel%npf%icalcspdis = 1
    end if
    !
    ! -- Set the auxiliary names for gwf flow packages in prt%fmi
    ngwfpack = gwfmodel%bndlist%Count()
    do ip = 1, ngwfpack
      packobj => GetBndFromList(gwfmodel%bndlist, ip)
      call prtmodel%fmi%gwfpackages(ip)%set_auxname(packobj%naux, &
                                                    packobj%auxname)
    end do
  end subroutine exg_df

  subroutine exg_ar(this)
    ! -- modules
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(GwfPrtExchangeType) :: this
    ! -- local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(PrtModelType), pointer :: prtmodel => null()
    ! -- formats
    character(len=*), parameter :: fmtdiserr = &
      "('GWF and PRT Models have incompatible discretizations for exchange&
      & ',a,'.&
      &  GWF Model has ', i0, ' user nodes and ', i0, ' reduced nodes.&
      &  PRT Model has ', i0, ' user nodes and ', i0, ' reduced nodes.&
      &  Ensure grids are identical.')"
    character(len=*), parameter :: fmtidomerr = &
      "('GWF and PRT Models have incompatible discretizations for exchange&
      & ',a,'.&
      &  An inactive cell (IDOMAIN <= 0) in the GWF Model is active &
      &in the PRT Model.&
      &  The PRT active domain must be a subset of the GWF active domain.')"
    !
    ! -- set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1id)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set prtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2id)
    select type (mb)
    type is (PrtModelType)
      prtmodel => mb
    end select
    !
    ! -- Check that the two models share the same underlying grid (same
    !    number of user nodes). The number of reduced nodes may differ:
    !    the PRT domain is allowed to be a subset of the GWF domain.
    if (prtmodel%dis%nodesuser /= gwfmodel%dis%nodesuser) then
      write (errmsg, fmtdiserr) trim(this%name), &
        gwfmodel%dis%nodesuser, &
        gwfmodel%dis%nodes, &
        prtmodel%dis%nodesuser, &
        prtmodel%dis%nodes
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- Check PRT active domain subsets GWF, build node map if not
    call check_and_map_domains(this, gwfmodel, prtmodel, fmtidomerr)
    !
    ! -- setup pointers to gwf variables allocated in gwf_ar
    prtmodel%fmi%gwfhead => gwfmodel%x
    call mem_checkin(prtmodel%fmi%gwfhead, &
                     'GWFHEAD', prtmodel%fmi%memoryPath, &
                     'X', gwfmodel%memoryPath)
    prtmodel%fmi%gwfsat => gwfmodel%npf%sat
    call mem_checkin(prtmodel%fmi%gwfsat, &
                     'GWFSAT', prtmodel%fmi%memoryPath, &
                     'SAT', gwfmodel%npf%memoryPath)
    prtmodel%fmi%gwfspdis => gwfmodel%npf%spdis
    call mem_checkin(prtmodel%fmi%gwfspdis, &
                     'GWFSPDIS', prtmodel%fmi%memoryPath, &
                     'SPDIS', gwfmodel%npf%memoryPath)
    prtmodel%fmi%igwfceltyp = 1
    prtmodel%fmi%gwfceltyp => gwfmodel%npf%icelltype
    call mem_checkin(prtmodel%fmi%gwfceltyp, &
                     'GWFCELTYP', prtmodel%fmi%memoryPath, &
                     'ICELLTYPE', gwfmodel%npf%memoryPath)
    !
    ! -- setup pointers to the flow storage rates. GWF strg arrays are
    !    available after the gwf_ar routine is called.
    if (prtmodel%inmst > 0) then
      if (gwfmodel%insto > 0) then
        prtmodel%fmi%gwfstrgss => gwfmodel%sto%strgss
        prtmodel%fmi%igwfstrgss = 1
        if (gwfmodel%sto%iusesy == 1) then
          prtmodel%fmi%gwfstrgsy => gwfmodel%sto%strgsy
          prtmodel%fmi%igwfstrgsy = 1
        end if
      end if
    end if

    ! -- transfer the boundary package information from gwf to prt
    call this%gwfbnd2prtfmi()

    ! -- if mover package is active, then set a pointer to it's budget object
    if (gwfmodel%inmvr /= 0) &
      prtmodel%fmi%mvrbudobj => gwfmodel%mvr%budobj

    ! -- todo connections
  end subroutine exg_ar

  !> @brief Verify PRT's active domain is a subset of GWF's, and if the two
  !! domains differ, build the node and connection maps needed to
  !! translate between GWF's and PRT's reduced numbering.
  !<
  subroutine check_and_map_domains(this, gwfmodel, prtmodel, fmtidomerr)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwfPrtExchangeType) :: this
    type(GwfModelType), pointer, intent(in) :: gwfmodel
    type(PrtModelType), pointer, intent(inout) :: prtmodel
    character(len=*), intent(in) :: fmtidomerr
    ! -- local
    integer(I4B) :: nu, gn, pn, gm, pm, ipos, jpos
    integer(I4B) :: noder_prt, noder_gwf
    logical :: differs
    !
    differs = .false.
    do nu = 1, prtmodel%dis%nodesuser
      noder_prt = prtmodel%dis%get_nodenumber(nu, 0)
      noder_gwf = gwfmodel%dis%get_nodenumber(nu, 0)
      if (noder_prt /= 0 .and. noder_gwf == 0) then
        write (errmsg, fmtidomerr) trim(this%name)
        call store_error(errmsg, terminate=.TRUE.)
      end if
      if ((noder_prt == 0) .neqv. (noder_gwf == 0)) differs = .true.
    end do
    !
    if (.not. differs) return
    !
    ! -- Build the node maps
    call mem_allocate(prtmodel%fmi%noder_gwf2prt, gwfmodel%dis%nodes, &
                      'NODER_GWF2PRT', prtmodel%fmi%memoryPath)
    call mem_allocate(prtmodel%fmi%noder_prt2gwf, prtmodel%dis%nodes, &
                      'NODER_PRT2GWF', prtmodel%fmi%memoryPath)
    do gn = 1, gwfmodel%dis%nodes
      nu = gwfmodel%dis%get_nodeuser(gn)
      prtmodel%fmi%noder_gwf2prt(gn) = prtmodel%dis%get_nodenumber(nu, 0)
    end do
    do pn = 1, prtmodel%dis%nodes
      nu = prtmodel%dis%get_nodeuser(pn)
      prtmodel%fmi%noder_prt2gwf(pn) = gwfmodel%dis%get_nodenumber(nu, 0)
    end do
    !
    ! -- Build the connection map: for each PRT reduced connection position,
    !    find the corresponding position in GWF's ia/ja. Requires that PRT
    !    active domain is a subset of GWF active domain.
    call mem_allocate(prtmodel%fmi%ipos_prt2gwf, prtmodel%dis%con%nja, &
                      'IPOS_PRT2GWF', prtmodel%fmi%memoryPath)
    do pn = 1, prtmodel%dis%nodes
      gn = prtmodel%fmi%noder_prt2gwf(pn)
      ! -- diagonal position maps directly to GWF's diagonal position
      prtmodel%fmi%ipos_prt2gwf(prtmodel%dis%con%ia(pn)) = gwfmodel%dis%con%ia(gn)
      do ipos = prtmodel%dis%con%ia(pn) + 1, prtmodel%dis%con%ia(pn + 1) - 1
        pm = prtmodel%dis%con%ja(ipos)
        gm = prtmodel%fmi%noder_prt2gwf(pm)
        prtmodel%fmi%ipos_prt2gwf(ipos) = 0
        do jpos = gwfmodel%dis%con%ia(gn) + 1, gwfmodel%dis%con%ia(gn + 1) - 1
          if (gwfmodel%dis%con%ja(jpos) == gm) then
            prtmodel%fmi%ipos_prt2gwf(ipos) = jpos
            exit
          end if
        end do
        if (prtmodel%fmi%ipos_prt2gwf(ipos) == 0) then
          write (errmsg, '(a,a)') 'Programmer error: could not map a PRT &
            &connection onto the corresponding GWF connection for &
            &exchange ', trim(this%name)
          call store_error(errmsg, terminate=.TRUE.)
        end if
      end do
    end do
  end subroutine check_and_map_domains

  ! todo subroutines: gwfconn2prtconn and link_connections

  subroutine exg_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GwfPrtExchangeType) :: this
    ! -- local
    !
    call mem_deallocate(this%m1id)
    call mem_deallocate(this%m2id)
  end subroutine exg_da

  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_checkin
    ! -- dummy
    class(GwfPrtExchangeType) :: this
    ! -- local
    !
    call mem_allocate(this%m1id, 'M1ID', this%memoryPath)
    call mem_allocate(this%m2id, 'M2ID', this%memoryPath)
    this%m1id = 0
    this%m2id = 0
  end subroutine allocate_scalars

  subroutine gwfbnd2prtfmi(this)
    ! -- modules
    ! -- dummy
    class(GwfPrtExchangeType) :: this
    ! -- local
    integer(I4B) :: ngwfpack, ip, iterm, imover
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(PrtModelType), pointer :: prtmodel => null()
    class(BndType), pointer :: packobj => null()
    !
    ! -- set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1id)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set prtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2id)
    select type (mb)
    type is (PrtModelType)
      prtmodel => mb
    end select
    !
    ! -- Call routines in FMI that will set pointers to the necessary flow
    !    data (SIMVALS and SIMTOMVR) stored within each GWF flow package
    ngwfpack = gwfmodel%bndlist%Count()
    iterm = 1
    do ip = 1, ngwfpack
      packobj => GetBndFromList(gwfmodel%bndlist, ip)
      call prtmodel%fmi%gwfpackages(iterm)%set_pointers( &
        'SIMVALS', &
        packobj%memoryPath, &
        packobj%input_mempath)
      iterm = iterm + 1
      !
      ! -- If a mover is active for this package, then establish a separate
      !    pointer link for the mover flows stored in SIMTOMVR
      imover = packobj%imover
      if (packobj%isadvpak /= 0) imover = 0
      if (imover /= 0) then
        call prtmodel%fmi%gwfpackages(iterm)%set_pointers( &
          'SIMTOMVR', &
          packobj%memoryPath, &
          packobj%input_mempath)
        iterm = iterm + 1
      end if
    end do
  end subroutine gwfbnd2prtfmi

end module GwfPrtExchangeModule
