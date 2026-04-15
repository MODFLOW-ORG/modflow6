module GwfPrtExchangeModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENPACKAGENAME
  use ListsModule, only: basemodellist, baseexchangelist
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use BaseModelModule, only: BaseModelType, GetBaseModelFromList
  use GwfModule, only: GwfModelType
  use PrtModule, only: PrtModelType
  use BndModule, only: BndType, GetBndFromList
  use PrtPrtExchangeModule, only: PrtPrtExchangeType, GetPrtPrtExchangeFromList
  use GwfGwfExchangeModule, only: GwfExchangeType, GetGwfExchangeFromList
  use MemoryManagerModule, only: mem_checkin, mem_allocate, mem_deallocate

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
    procedure, private :: gwfconn2prtconn
    procedure, private :: link_connections

  end type GwfPrtExchangeType

contains

  !> @brief Create a new GWF to PRT exchange object
  subroutine gwfprt_cr(filename, id, m1id, m2id)
    ! modules
    use SimVariablesModule, only: model_loc_idx
    ! dummy
    character(len=*), intent(in) :: filename
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: m1id
    integer(I4B), intent(in) :: m2id
    ! local
    class(BaseExchangeType), pointer :: baseexchange => null()
    type(GwfPrtExchangeType), pointer :: exchange => null()
    character(len=20) :: cint

    ! create a new exchange and add it to the baseexchangelist container
    allocate (exchange)
    baseexchange => exchange
    call AddBaseExchangeToList(baseexchangelist, baseexchange)

    ! assign id and name
    exchange%id = id
    write (cint, '(i0)') id
    exchange%name = 'GWF-PRT_'//trim(adjustl(cint))
    exchange%memoryPath = exchange%name

    ! allocate scalars
    call exchange%allocate_scalars()

    ! NB: convert from id to local model index in base model list
    exchange%m1id = model_loc_idx(m1id)
    exchange%m2id = model_loc_idx(m2id)

    ! set model pointers
    call exchange%set_model_pointers()
  end subroutine gwfprt_cr

  subroutine set_model_pointers(this)
    ! dummy
    class(GwfPrtExchangeType) :: this
    ! local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwf => null()
    type(PrtModelType), pointer :: prt => null()

    ! set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1id)
    select type (mb)
    type is (GwfModelType)
      gwf => mb
    end select

    ! set prtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2id)
    select type (mb)
    type is (PrtModelType)
      prt => mb
    end select

    ! Verify that GWF model is of the correct type
    if (.not. associated(gwf)) then
      write (errmsg, '(3a)') 'Problem with GWF-PRT exchange ', trim(this%name), &
        '.  Specified GWF Model does not appear to be of the correct type.'
      call store_error(errmsg, terminate=.true.)
    end if

    ! Verify that PRT model is of the correct type
    if (.not. associated(prt)) then
      write (errmsg, '(3a)') 'Problem with GWF-PRT exchange ', trim(this%name), &
        '.  Specified PRT Model does not appear to be of the correct type.'
      call store_error(errmsg, terminate=.true.)
    end if

    ! Tell particle tracking model fmi flows are not read from file
    prt%fmi%flows_from_file = .false.

    ! Set a pointer to the GWF bndlist. This will allow the tracking model
    ! to look through the flow packages and establish a link to GWF flows
    prt%fmi%gwfbndlist => gwf%bndlist
  end subroutine set_model_pointers

  subroutine exg_df(this)
    ! dummy
    class(GwfPrtExchangeType) :: this
    ! local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwf => null()
    type(PrtModelType), pointer :: prt => null()
    integer(I4B) :: ngwfpack, ip
    class(BndType), pointer :: packobj => null()

    ! set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1id)
    select type (mb)
    type is (GwfModelType)
      gwf => mb
    end select

    ! set prtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2id)
    select type (mb)
    type is (PrtModelType)
      prt => mb
    end select

    ! Check to make sure that flow is solved before particle tracking and in a
    ! different solution
    if (gwf%idsoln >= prt%idsoln) then
      write (errmsg, '(3a)') 'Problem with GWF-PRT exchange ', trim(this%name), &
        '.  The GWF model must be solved by a different solution than the PRT model. &
        &The IMS specified for GWF must be listed in mfsim.nam &
        &before the EMS for PRT.'
      call store_error(errmsg, terminate=.true.)
    end if

    ! Set pointer to flowja
    prt%fmi%gwfflowja => gwf%flowja
    call mem_checkin(prt%fmi%gwfflowja, &
                     'GWFFLOWJA', prt%fmi%memoryPath, &
                     'FLOWJA', gwf%memoryPath)

    ! Set the npf flag so that specific discharge is available for
    ! transport calculations if dispersion is active
    if (prt%indsp > 0) then
      gwf%npf%icalcspdis = 1
    end if

    ! Set the auxiliary names for gwf flow packages in prt%fmi
    ngwfpack = gwf%bndlist%Count()
    do ip = 1, ngwfpack
      packobj => GetBndFromList(gwf%bndlist, ip)
      call prt%fmi%gwfpackages(ip)%set_auxname(packobj%naux, &
                                               packobj%auxname)
    end do
  end subroutine exg_df

  subroutine exg_ar(this)
    ! modules
    use DisModule, only: DisType
    use DisvModule, only: DisvType
    use DisuModule, only: DisuType
    ! dummy
    class(GwfPrtExchangeType) :: this
    ! local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwf => null()
    type(PrtModelType), pointer :: prt => null()
    ! formats
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

    ! set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1id)
    select type (mb)
    type is (GwfModelType)
      gwf => mb
    end select

    ! set prtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2id)
    select type (mb)
    type is (PrtModelType)
      prt => mb
    end select

    ! Check to make sure sizes are identical
    if (prt%dis%nodes /= gwf%dis%nodes .or. &
        prt%dis%nodesuser /= gwf%dis%nodesuser) then
      write (errmsg, fmtdiserr) trim(this%name), &
        gwf%dis%nodesuser, &
        gwf%dis%nodes, &
        prt%dis%nodesuser, &
        prt%dis%nodes
      call store_error(errmsg, terminate=.TRUE.)
    end if

    ! Make sure idomains are identical
    select type (gwfdis => gwf%dis)
    type is (DisType)
      select type (prtdis => prt%dis)
      type is (DisType)
        if (.not. all(gwfdis%idomain == prtdis%idomain)) then
          write (errmsg, fmtidomerr) trim(this%name)
          call store_error(errmsg, terminate=.TRUE.)
        end if
      end select
    type is (DisvType)
      select type (prtdis => prt%dis)
      type is (DisvType)
        if (.not. all(gwfdis%idomain == prtdis%idomain)) then
          write (errmsg, fmtidomerr) trim(this%name)
          call store_error(errmsg, terminate=.TRUE.)
        end if
      end select
    type is (DisuType)
      select type (prtdis => prt%dis)
      type is (DisuType)
        if (.not. all(gwfdis%idomain == prtdis%idomain)) then
          write (errmsg, fmtidomerr) trim(this%name)
          call store_error(errmsg, terminate=.TRUE.)
        end if
      end select
    end select

    ! setup pointers to gwf variables allocated in gwf_ar
    prt%fmi%gwfhead => gwf%x
    call mem_checkin(prt%fmi%gwfhead, &
                     'GWFHEAD', prt%fmi%memoryPath, &
                     'X', gwf%memoryPath)
    prt%fmi%gwfsat => gwf%npf%sat
    call mem_checkin(prt%fmi%gwfsat, &
                     'GWFSAT', prt%fmi%memoryPath, &
                     'SAT', gwf%npf%memoryPath)
    prt%fmi%gwfspdis => gwf%npf%spdis
    call mem_checkin(prt%fmi%gwfspdis, &
                     'GWFSPDIS', prt%fmi%memoryPath, &
                     'SPDIS', gwf%npf%memoryPath)
    prt%fmi%igwfceltyp = 1
    prt%fmi%gwfceltyp => gwf%npf%icelltype
    call mem_checkin(prt%fmi%gwfceltyp, &
                     'GWFCELTYP', prt%fmi%memoryPath, &
                     'ICELLTYPE', gwf%npf%memoryPath)

    ! setup pointers to the flow storage rates. GWF strg arrays are
    ! available after the gwf_ar routine is called.
    if (prt%inmst > 0) then
      if (gwf%insto > 0) then
        prt%fmi%gwfstrgss => gwf%sto%strgss
        prt%fmi%igwfstrgss = 1
        if (gwf%sto%iusesy == 1) then
          prt%fmi%gwfstrgsy => gwf%sto%strgsy
          prt%fmi%igwfstrgsy = 1
        end if
      end if
    end if

    ! transfer the boundary package information from gwf to prt
    call this%gwfbnd2prtfmi()

    ! if mover package is active, then set a pointer to it's budget object
    if (gwf%inmvr /= 0) &
      prt%fmi%mvrbudobj => gwf%mvr%budobj

    ! link PRT-PRT exchanges to GWF-GWF exchanges
    call this%gwfconn2prtconn(gwf, prt)
  end subroutine exg_ar

  subroutine exg_da(this)
    class(GwfPrtExchangeType) :: this

    call mem_deallocate(this%m1id)
    call mem_deallocate(this%m2id)
  end subroutine exg_da

  subroutine allocate_scalars(this)
    class(GwfPrtExchangeType) :: this

    call mem_allocate(this%m1id, 'M1ID', this%memoryPath)
    call mem_allocate(this%m2id, 'M2ID', this%memoryPath)
    this%m1id = 0
    this%m2id = 0
  end subroutine allocate_scalars

  subroutine gwfbnd2prtfmi(this)
    ! dummy
    class(GwfPrtExchangeType) :: this
    ! local
    integer(I4B) :: ngwfpack, ip, iterm, imover
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwf => null()
    type(PrtModelType), pointer :: prt => null()
    class(BndType), pointer :: packobj => null()

    ! set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1id)
    select type (mb)
    type is (GwfModelType)
      gwf => mb
    end select

    ! set prtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2id)
    select type (mb)
    type is (PrtModelType)
      prt => mb
    end select

    ! call routines in FMI that will set pointers to the necessary flow
    ! data (SIMVALS and SIMTOMVR) stored within each GWF flow package
    ngwfpack = gwf%bndlist%Count()
    iterm = 1
    do ip = 1, ngwfpack
      packobj => GetBndFromList(gwf%bndlist, ip)
      call prt%fmi%gwfpackages(iterm)%set_pointers( &
        'SIMVALS', &
        packobj%memoryPath, &
        packobj%input_mempath)
      iterm = iterm + 1

      ! if a mover is active for this package, then establish a separate
      ! pointer link for the mover flows stored in SIMTOMVR
      imover = packobj%imover
      if (packobj%isadvpak /= 0) imover = 0
      if (imover /= 0) then
        call prt%fmi%gwfpackages(iterm)%set_pointers( &
          'SIMTOMVR', &
          packobj%memoryPath, &
          packobj%input_mempath)
        iterm = iterm + 1
      end if
    end do
  end subroutine gwfbnd2prtfmi

  !> @brief Link PRT-PRT exchange to GWF-GWF exchange
  !!
  !! Connects a PRT-PRT exchange to its corresponding GWF-GWF exchange
  !! to provide flow data for particle tracking across model boundaries.
  !<
  subroutine link_connections(this, prt_exg, gwf_exg)
    use MemoryManagerModule, only: mem_checkin
    ! dummy
    class(GwfPrtExchangeType) :: this
    class(PrtPrtExchangeType), pointer :: prt_exg
    class(GwfExchangeType), pointer :: gwf_exg

    ! Only link if not already linked (prevents duplicate mem_checkin)
    if (.not. associated(prt_exg%gwfsimvals)) then
      prt_exg%gwfsimvals => gwf_exg%simvals
      call mem_checkin(prt_exg%gwfsimvals, &
                       'GWFSIMVALS', prt_exg%memoryPath, &
                       'SIMVALS', gwf_exg%memoryPath)
    end if
  end subroutine link_connections

  !> @brief Find and link PRT-PRT exchanges to GWF-GWF exchanges
  !!
  !! Searches for PRT-PRT exchanges connected to the PRT model and matches
  !! them to corresponding GWF-GWF exchanges based on model names and connection
  !! topology. When a match is found, links the flow data from the GWF-GWF
  !! exchange to the PRT-PRT exchange.
  !<
  subroutine gwfconn2prtconn(this, gwf, prt)
    use SimModule, only: store_error, store_error_filename, count_errors
    use SimVariablesModule, only: iout
    ! dummy
    class(GwfPrtExchangeType) :: this
    type(GwfModelType), pointer :: gwf
    type(PrtModelType), pointer :: prt
    ! local
    class(PrtPrtExchangeType), pointer :: prt_exg => null()
    class(GwfExchangeType), pointer :: gwf_exg => null()
    class(BaseModelType), pointer :: prt_base => null()
    class(BaseModelType), pointer :: gwf_base => null()
    integer(I4B) :: iex_prt, iex_gwf
    logical(LGP) :: eq

    prt_base => prt
    gwf_base => gwf

    ! iterate PRT-PRT exchanges
    do iex_prt = 1, baseexchangelist%Count()
      prt_exg => GetPrtPrtExchangeFromList(baseexchangelist, iex_prt)
      if (.not. associated(prt_exg)) cycle
      if (.not. prt_exg%connects_model(prt_base)) cycle

      ! find matching GWF-GWF exchange
      do iex_gwf = 1, baseexchangelist%Count()
        gwf_exg => GetGwfExchangeFromList(baseexchangelist, iex_gwf)
        if (.not. associated(gwf_exg)) cycle
        if (.not. gwf_exg%connects_model(gwf_base)) cycle

        ! match criteria:
        ! 1. GWF model names match PRT's declared GWF models
        ! 2. exchange arrays (nodem1, nodem2) are identical
        if (gwf_exg%v_model1%name /= prt_exg%gwfmodelname1) cycle
        if (gwf_exg%v_model2%name /= prt_exg%gwfmodelname2) cycle

        eq = (gwf_exg%nexg == prt_exg%nexg)
        if (eq) then
          eq = all(gwf_exg%nodem1 == prt_exg%nodem1)
          eq = eq .and. all(gwf_exg%nodem2 == prt_exg%nodem2)
        end if

        if (eq) then
          write (iout, '(/6a)') 'Linking PRT-PRT exchange ', &
            trim(prt_exg%name), ' to GWF-GWF exchange ', trim(gwf_exg%name), &
            ' for PRT model ', trim(prt%name)
          call this%link_connections(prt_exg, gwf_exg)
          exit
        end if
      end do
    end do
  end subroutine gwfconn2prtconn

end module GwfPrtExchangeModule
