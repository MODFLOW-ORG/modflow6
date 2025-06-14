module GwfGweExchangeModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENPACKAGENAME
  use ListsModule, only: basemodellist, baseexchangelist, &
                         baseconnectionlist
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use SpatialModelConnectionModule, only: SpatialModelConnectionType, &
                                          get_smc_from_list
  use GweGweConnectionModule, only: GweGweConnectionType, CastAsGweGweConnection
  use GwfGwfConnectionModule, only: GwfGwfConnectionType, CastAsGwfGwfConnection
  use GwfGwfExchangeModule, only: GwfExchangeType, &
                                  GetGwfExchangeFromList
  use BaseModelModule, only: BaseModelType, GetBaseModelFromList
  use GwfModule, only: GwfModelType
  use GweModule, only: GweModelType
  use BndModule, only: BndType, GetBndFromList

  implicit none
  public :: GwfGweExchangeType
  public :: gwfgwe_cr

  type, extends(BaseExchangeType) :: GwfGweExchangeType

    integer(I4B), pointer :: m1_idx => null() !< index into the list of base exchanges for model 1
    integer(I4B), pointer :: m2_idx => null() !< index into the list of base exchanges for model 2

  contains

    procedure :: exg_df
    procedure :: exg_ar
    procedure :: exg_da
    procedure, private :: set_model_pointers
    procedure, private :: allocate_scalars
    procedure, private :: gwfbnd2gwefmi
    procedure, private :: gwfconn2gweconn
    procedure, private :: link_connections

  end type GwfGweExchangeType

contains

  !> @brief Create a new GWF to GWE exchange object
  !<
  subroutine gwfgwe_cr(filename, id, m1_id, m2_id)
    ! -- modules
    use SimVariablesModule, only: model_loc_idx
    ! -- dummy
    character(len=*), intent(in) :: filename
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: m1_id
    integer(I4B), intent(in) :: m2_id
    ! -- local
    class(BaseExchangeType), pointer :: baseexchange => null()
    type(GwfGweExchangeType), pointer :: exchange => null()
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
    exchange%name = 'GWF-GWE_'//trim(adjustl(cint))
    exchange%memoryPath = exchange%name
    !
    ! -- allocate scalars
    call exchange%allocate_scalars()
    !
    ! -- NB: convert from id to local model index in base model list
    exchange%m1_idx = model_loc_idx(m1_id)
    exchange%m2_idx = model_loc_idx(m2_id)
    !
    ! -- set model pointers
    call exchange%set_model_pointers()
  end subroutine gwfgwe_cr

  !> @brief Allocate and read
  !<
  subroutine set_model_pointers(this)
    ! -- dummy
    class(GwfGweExchangeType) :: this
    ! -- local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(GweModelType), pointer :: gwemodel => null()
    !
    ! -- set gwfmodel
    gwfmodel => null()
    mb => GetBaseModelFromList(basemodellist, this%m1_idx)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set gwemodel
    gwemodel => null()
    mb => GetBaseModelFromList(basemodellist, this%m2_idx)
    select type (mb)
    type is (GweModelType)
      gwemodel => mb
    end select
    !
    ! -- Verify that gwf model is of the correct type
    if (.not. associated(gwfmodel)) then
      write (errmsg, '(3a)') 'Problem with GWF-GWE exchange ', trim(this%name), &
        '.  Specified GWF Model does not appear to be of the correct type.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Verify that gwe model is of the correct type
    if (.not. associated(gwemodel)) then
      write (errmsg, '(3a)') 'Problem with GWF-GWE exchange ', trim(this%name), &
        '.  Specified GWF Model does not appear to be of the correct type.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Tell transport model fmi flows are not read from file
    gwemodel%fmi%flows_from_file = .false.
    !
    ! -- Set a pointer to the GWF bndlist.  This will allow the transport model
    !    to look through the flow packages and establish a link to GWF flows
    gwemodel%fmi%gwfbndlist => gwfmodel%bndlist
  end subroutine set_model_pointers

  !> @brief Define the GwfGwe Exchange object
  !<
  subroutine exg_df(this)
    ! -- modules
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(GwfGweExchangeType) :: this
    ! -- local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(GweModelType), pointer :: gwemodel => null()
    !
    ! -- set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1_idx)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set gwemodel
    mb => GetBaseModelFromList(basemodellist, this%m2_idx)
    select type (mb)
    type is (GweModelType)
      gwemodel => mb
    end select
    !
    ! -- Check to make sure that flow is solved before transport and in a
    !    different IMS solution
    if (gwfmodel%idsoln >= gwemodel%idsoln) then
      write (errmsg, '(3a)') 'Problem with GWF-GWE exchange ', trim(this%name), &
        '.  The GWF model must be solved by a different IMS than the GWE model. &
        &Furthermore, the IMS specified for GWF must be listed in mfsim.nam &
        &before the IMS for GWE.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Set pointer to flowja
    gwemodel%fmi%gwfflowja => gwfmodel%flowja
    call mem_checkin(gwemodel%fmi%gwfflowja, &
                     'GWFFLOWJA', gwemodel%fmi%memoryPath, &
                     'FLOWJA', gwfmodel%memoryPath)

    !
    ! -- Set the npf flag so that specific discharge is available for
    !    transport calculations if dispersion is active
    if (gwemodel%incnd > 0) then
      gwfmodel%npf%icalcspdis = 1
    end if
  end subroutine exg_df

  !> @brief Allocate and read
  !<
  subroutine exg_ar(this)
    ! -- modules
    use MemoryManagerModule, only: mem_checkin
    use DisModule, only: DisType
    use DisvModule, only: DisvType
    use DisuModule, only: DisuType
    ! -- dummy
    class(GwfGweExchangeType) :: this
    ! -- local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(GweModelType), pointer :: gwemodel => null()
    ! -- formats
    character(len=*), parameter :: fmtdiserr = &
      "('GWF and GWE Models do not have the same discretization for exchange&
      & ',a,'.&
      &  GWF Model has ', i0, ' user nodes and ', i0, ' reduced nodes.&
      &  GWE Model has ', i0, ' user nodes and ', i0, ' reduced nodes.&
      &  Ensure discretization packages, including IDOMAIN, are identical.')"
    character(len=*), parameter :: fmtidomerr = &
      "('GWF and GWE Models do not have the same discretization for exchange&
      & ',a,'.&
      &  GWF Model and GWE Model have different IDOMAIN arrays.&
      &  Ensure discretization packages, including IDOMAIN, are identical.')"
    !
    ! -- set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1_idx)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set gwemodel
    mb => GetBaseModelFromList(basemodellist, this%m2_idx)
    select type (mb)
    type is (GweModelType)
      gwemodel => mb
    end select
    !
    ! -- Check to make sure sizes are identical
    if (gwemodel%dis%nodes /= gwfmodel%dis%nodes .or. &
        gwemodel%dis%nodesuser /= gwfmodel%dis%nodesuser) then
      write (errmsg, fmtdiserr) trim(this%name), &
        gwfmodel%dis%nodesuser, &
        gwfmodel%dis%nodes, &
        gwemodel%dis%nodesuser, &
        gwemodel%dis%nodes
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- Make sure idomains are identical
    if (.not. all(gwfmodel%dis%idomain == gwemodel%dis%idomain)) then
      write (errmsg, fmtidomerr) trim(this%name)
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- setup pointers to gwf variables allocated in gwf_ar
    gwemodel%fmi%gwfhead => gwfmodel%x
    call mem_checkin(gwemodel%fmi%gwfhead, &
                     'GWFHEAD', gwemodel%fmi%memoryPath, &
                     'X', gwfmodel%memoryPath)
    gwemodel%fmi%gwfsat => gwfmodel%npf%sat
    call mem_checkin(gwemodel%fmi%gwfsat, &
                     'GWFSAT', gwemodel%fmi%memoryPath, &
                     'SAT', gwfmodel%npf%memoryPath)
    gwemodel%fmi%gwfspdis => gwfmodel%npf%spdis
    call mem_checkin(gwemodel%fmi%gwfspdis, &
                     'GWFSPDIS', gwemodel%fmi%memoryPath, &
                     'SPDIS', gwfmodel%npf%memoryPath)
    !
    ! -- setup pointers to the flow storage rates. GWF strg arrays are
    !    available after the gwf_ar routine is called.
    if (gwemodel%inest > 0) then
      if (gwfmodel%insto > 0) then
        gwemodel%fmi%gwfstrgss => gwfmodel%sto%strgss
        gwemodel%fmi%igwfstrgss = 1
        if (gwfmodel%sto%iusesy == 1) then
          gwemodel%fmi%gwfstrgsy => gwfmodel%sto%strgsy
          gwemodel%fmi%igwfstrgsy = 1
        end if
      end if
    end if
    !
    ! -- Set a pointer to conc in buy
    if (gwfmodel%inbuy > 0) then
      call gwfmodel%buy%set_concentration_pointer(gwemodel%name, gwemodel%x, &
                                                  gwemodel%ibound)
    end if
    !
    ! -- Set a pointer to conc (which could be a temperature) in vsc
    if (gwfmodel%invsc > 0) then
      call gwfmodel%vsc%set_concentration_pointer(gwemodel%name, gwemodel%x, &
                                                  gwemodel%ibound, 1)
    end if
    !
    ! -- transfer the boundary package information from gwf to gwe
    call this%gwfbnd2gwefmi()
    !
    ! -- if mover package is active, then set a pointer to it's budget object
    if (gwfmodel%inmvr /= 0) then
      gwemodel%fmi%mvrbudobj => gwfmodel%mvr%budobj
    end if
    !
    ! -- connect Connections
    call this%gwfconn2gweconn(gwfmodel, gwemodel)
  end subroutine exg_ar

  !> @brief Link GWE connections to GWF connections or exchanges
  !<
  subroutine gwfconn2gweconn(this, gwfModel, gweModel)
    ! -- modules
    use SimModule, only: store_error
    use SimVariablesModule, only: iout
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(GwfGweExchangeType) :: this !< this exchange
    type(GwfModelType), pointer :: gwfModel !< the flow model
    type(GweModelType), pointer :: gweModel !< the energy transport model
    ! -- local
    class(SpatialModelConnectionType), pointer :: conn => null()
    class(*), pointer :: objPtr => null()
    class(GweGweConnectionType), pointer :: gweConn => null()
    class(GwfGwfConnectionType), pointer :: gwfConn => null()
    class(GwfExchangeType), pointer :: gwfEx => null()
    integer(I4B) :: ic1, ic2, iex
    integer(I4B) :: gwfConnIdx, gwfExIdx
    logical(LGP) :: areEqual
    !
    ! loop over all connections
    gweloop: do ic1 = 1, baseconnectionlist%Count()
      !
      conn => get_smc_from_list(baseconnectionlist, ic1)
      if (.not. associated(conn%owner, gweModel)) cycle gweloop
      !
      ! start with a GWE conn.
      objPtr => conn
      gweConn => CastAsGweGweConnection(objPtr)
      gwfConnIdx = -1
      gwfExIdx = -1
      !
      ! find matching GWF conn. in same list
      gwfloop: do ic2 = 1, baseconnectionlist%Count()
        conn => get_smc_from_list(baseconnectionlist, ic2)
        !
        if (associated(conn%owner, gwfModel)) then
          objPtr => conn
          gwfConn => CastAsGwfGwfConnection(objPtr)
          !
          ! for now, connecting the same nodes nrs will be
          ! sufficient evidence of equality
          areEqual = all(gwfConn%prim_exchange%nodem1 == &
                         gweConn%prim_exchange%nodem1)
          areEqual = areEqual .and. all(gwfConn%prim_exchange%nodem2 == &
                                        gweConn%prim_exchange%nodem2)
          if (areEqual) then
            ! same DIS, same exchange: link and go to next GWE conn.
            write (iout, '(/6a)') 'Linking exchange ', &
              trim(gweConn%prim_exchange%name), &
              ' to ', trim(gwfConn%prim_exchange%name), &
              ' (using interface model) for GWE model ', &
              trim(gweModel%name)
            gwfConnIdx = ic2
            call this%link_connections(gweConn, gwfConn)
            exit gwfloop
          end if
        end if
      end do gwfloop
      !
      ! fallback option: coupling to old gwfgwf exchange,
      ! (this will go obsolete at some point)
      if (gwfConnIdx == -1) then
        gwfloopexg: do iex = 1, baseexchangelist%Count()
          gwfEx => GetGwfExchangeFromList(baseexchangelist, iex)
          !
          ! -- There is no guarantee that iex is a gwfExg, in which case
          !    it will return as null.  cycle if so.
          if (.not. associated(gwfEx)) cycle gwfloopexg
          !
          if (associated(gwfEx%model1, gwfModel) .or. &
              associated(gwfEx%model2, gwfModel)) then

            ! check exchanges have same node counts
            areEqual = size(gwfEx%nodem1) == size(gweConn%prim_exchange%nodem1)
            ! then, connecting the same nodes nrs will be
            ! sufficient evidence of equality
            if (areEqual) &
              areEqual = all(gwfEx%nodem1 == gweConn%prim_exchange%nodem1)
            if (areEqual) &
              areEqual = all(gwfEx%nodem2 == gweConn%prim_exchange%nodem2)
            if (areEqual) then
              ! link exchange to connection
              write (iout, '(/6a)') 'Linking exchange ', &
                trim(gweConn%prim_exchange%name), &
                ' to ', trim(gwfEx%name), ' for GWE model ', &
                trim(gweModel%name)
              gwfExIdx = iex
              if (gweConn%owns_exchange) then
                gweConn%gweExchange%gwfsimvals => gwfEx%simvals
                call mem_checkin(gweConn%gweExchange%gwfsimvals, &
                                 'GWFSIMVALS', gweConn%gweExchange%memoryPath, &
                                 'SIMVALS', gwfEx%memoryPath)
              end if
              !
              !cdl link up mvt to mvr
              if (gwfEx%inmvr > 0) then
                if (gweConn%owns_exchange) then
                  !cdl todo: check and make sure gweEx has mvt active
                  call gweConn%gweExchange%mvt%set_pointer_mvrbudobj( &
                    gwfEx%mvr%budobj)
                end if
              end if
              !
              if (associated(gwfEx%model2, gwfModel)) gweConn%exgflowSign = -1
              gweConn%gweInterfaceModel%fmi%flows_from_file = .false.
              !
              exit gwfloopexg
            end if
          end if
          !
        end do gwfloopexg
      end if
      !
      if (gwfConnIdx == -1 .and. gwfExIdx == -1) then
        ! none found, report
        write (errmsg, '(/6a)') 'Missing GWF-GWF exchange when connecting GWE'// &
          ' model ', trim(gweModel%name), ' with exchange ', &
          trim(gweConn%prim_exchange%name), ' to GWF model ', &
          trim(gwfModel%name)
        call store_error(errmsg, terminate=.true.)
      end if
      !
    end do gweloop
  end subroutine gwfconn2gweconn

  !> @brief Links a GWE connection to its GWF counterpart
  !<
  subroutine link_connections(this, gweConn, gwfConn)
    ! -- modules
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(GwfGweExchangeType) :: this !< this exchange
    class(GweGweConnectionType), pointer :: gweConn !< GWE connection
    class(GwfGwfConnectionType), pointer :: gwfConn !< GWF connection
    !
    !gweConn%exgflowja => gwfConn%exgflowja
    if (gweConn%owns_exchange) then
      gweConn%gweExchange%gwfsimvals => gwfConn%gwfExchange%simvals
      call mem_checkin(gweConn%gweExchange%gwfsimvals, &
                       'GWFSIMVALS', gweConn%gweExchange%memoryPath, &
                       'SIMVALS', gwfConn%gwfExchange%memoryPath)
    end if
    !
    !cdl link up mvt to mvr
    if (gwfConn%gwfExchange%inmvr > 0) then
      if (gweConn%owns_exchange) then
        !cdl todo: check and make sure gweEx has mvt active
        call gweConn%gweExchange%mvt%set_pointer_mvrbudobj( &
          gwfConn%gwfExchange%mvr%budobj)
      end if
    end if
    !
    if (associated(gwfConn%gwfExchange%model2, gwfConn%owner)) then
      gweConn%exgflowSign = -1
    end if
    !
    ! fmi flows are not read from file
    gweConn%gweInterfaceModel%fmi%flows_from_file = .false.
    !
    ! set concentration pointer for buoyancy
    !call gwfConn%gwfInterfaceModel%buy%set_concentration_pointer( &
    !  gweConn%gweModel%name, &
    !  gweConn%conc, &
    !  gweConn%icbound)
  end subroutine link_connections

  !> @brief Deallocate memory
  !<
  subroutine exg_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GwfGweExchangeType) :: this
    !
    call mem_deallocate(this%m1_idx)
    call mem_deallocate(this%m2_idx)
  end subroutine exg_da

  !> @brief Allocate GwfGwe exchange scalars
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwfGweExchangeType) :: this
    !
    call mem_allocate(this%m1_idx, 'M1ID', this%memoryPath)
    call mem_allocate(this%m2_idx, 'M2ID', this%memoryPath)
    this%m1_idx = 0
    this%m2_idx = 0
  end subroutine allocate_scalars

  !> @brief Call routines in FMI that will set pointers to the necessary flow
  !! data (SIMVALS and SIMTOMVR) stored within each GWF flow package
  !<
  subroutine gwfbnd2gwefmi(this)
    ! -- dummy
    class(GwfGweExchangeType) :: this
    ! -- local
    integer(I4B) :: ngwfpack, ip, iterm, imover
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(GweModelType), pointer :: gwemodel => null()
    class(BndType), pointer :: packobj => null()
    !
    ! -- set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1_idx)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set gwemodel
    mb => GetBaseModelFromList(basemodellist, this%m2_idx)
    select type (mb)
    type is (GweModelType)
      gwemodel => mb
    end select
    !
    ! -- Call routines in FMI that will set pointers to the necessary flow
    !    data (SIMVALS and SIMTOMVR) stored within each GWF flow package
    ngwfpack = gwfmodel%bndlist%Count()
    iterm = 1
    do ip = 1, ngwfpack
      packobj => GetBndFromList(gwfmodel%bndlist, ip)
      call gwemodel%fmi%gwfpackages(iterm)%set_pointers( &
        'SIMVALS', &
        packobj%memoryPath, packobj%input_mempath)
      iterm = iterm + 1
      !
      ! -- If a mover is active for this package, then establish a separate
      !    pointer link for the mover flows stored in SIMTOMVR
      imover = packobj%imover
      if (packobj%isadvpak /= 0) imover = 0
      if (imover /= 0) then
        call gwemodel%fmi%gwfpackages(iterm)%set_pointers( &
          'SIMTOMVR', &
          packobj%memoryPath, packobj%input_mempath)
        iterm = iterm + 1
      end if
    end do
  end subroutine gwfbnd2gwefmi

end module GwfGweExchangeModule
