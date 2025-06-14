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
    use DisModule, only: DisType
    use DisvModule, only: DisvType
    use DisuModule, only: DisuType
    ! -- dummy
    class(GwfPrtExchangeType) :: this
    ! -- local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(PrtModelType), pointer :: prtmodel => null()
    ! -- formats
    character(len=*), parameter :: fmtdiserr = &
      "('GWF and PRT Models do not have the same discretization for exchange&
      & ',a,'.&
      &  GWF Model has ', i0, ' user nodes and ', i0, ' reduced nodes.&
      &  PRT Model has ', i0, ' user nodes and ', i0, ' reduced nodes.&
      &  Ensure discretization packages, including IDOMAIN, are identical.')"
    character(len=*), parameter :: fmtidomerr = &
      "('GWF and PRT Models do not have the same discretization for exchange&
      & ',a,'.&
      &  GWF Model and PRT Model have different IDOMAIN arrays.&
      &  Ensure discretization packages, including IDOMAIN, are identical.')"
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
    ! -- Check to make sure sizes are identical
    if (prtmodel%dis%nodes /= gwfmodel%dis%nodes .or. &
        prtmodel%dis%nodesuser /= gwfmodel%dis%nodesuser) then
      write (errmsg, fmtdiserr) trim(this%name), &
        gwfmodel%dis%nodesuser, &
        gwfmodel%dis%nodes, &
        prtmodel%dis%nodesuser, &
        prtmodel%dis%nodes
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- Make sure idomains are identical
    if (.not. all(gwfmodel%dis%idomain == prtmodel%dis%idomain)) then
      write (errmsg, fmtidomerr) trim(this%name)
      call store_error(errmsg, terminate=.TRUE.)
    end if
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
