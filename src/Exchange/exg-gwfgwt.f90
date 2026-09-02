module GwfGwtExchangeModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENPACKAGENAME, LINELENGTH
  use ListsModule, only: basemodellist, baseexchangelist, &
                         baseconnectionlist
  use SimModule, only: store_error, store_error_filename, store_warning, &
                       count_errors
  use SimVariablesModule, only: errmsg, warnmsg
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use SpatialModelConnectionModule, only: SpatialModelConnectionType, &
                                          get_smc_from_list
  use GwtGwtConnectionModule, only: GwtGwtConnectionType, CastAsGwtGwtConnection
  use GwfGwfConnectionModule, only: GwfGwfConnectionType, CastAsGwfGwfConnection
  use GwfGwfExchangeModule, only: GwfExchangeType, &
                                  GetGwfExchangeFromList
  use GwtGwtExchangeModule, only: GwtExchangeType
  use DisConnExchangeModule, only: same_exchange_cells
  use BaseModelModule, only: BaseModelType, GetBaseModelFromList
  use GwfModule, only: GwfModelType
  use GwtModule, only: GwtModelType
  use BndModule, only: BndType, GetBndFromList

  implicit none
  public :: GwfGwtExchangeType
  public :: gwfgwt_cr

  type, extends(BaseExchangeType) :: GwfGwtExchangeType

    integer(I4B), pointer :: m1_idx => null() !< index into the list of base exchanges for model 1
    integer(I4B), pointer :: m2_idx => null() !< index into the list of base exchanges for model 2
    character(len=LINELENGTH) :: filename !< the input file for the GWF-GWT exchange

  contains

    procedure :: exg_df
    procedure :: exg_ar
    procedure :: exg_da
    procedure, private :: check_discretization
    procedure, private :: set_model_pointers
    procedure, private :: allocate_scalars
    procedure, private :: gwfbnd2gwtfmi
    procedure, private :: gwfconn2gwtconn
    procedure, private :: link_connections

  end type GwfGwtExchangeType

contains

  !> @brief Create a new GWF to GWT exchange object
  !<
  subroutine gwfgwt_cr(filename, id, m1_id, m2_id)
    ! -- modules
    use SimVariablesModule, only: model_loc_idx
    ! -- dummy
    character(len=*), intent(in) :: filename
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: m1_id
    integer(I4B), intent(in) :: m2_id
    ! -- local
    class(BaseExchangeType), pointer :: baseexchange => null()
    type(GwfGwtExchangeType), pointer :: exchange => null()
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
    exchange%name = 'GWF-GWT_'//trim(adjustl(cint))
    exchange%memoryPath = exchange%name
    exchange%filename = filename
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
  end subroutine gwfgwt_cr

  !> @brief Allocate and read
  !<
  subroutine set_model_pointers(this)
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    ! -- local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(GwtModelType), pointer :: gwtmodel => null()
    !
    ! -- set gwfmodel
    gwfmodel => null()
    mb => GetBaseModelFromList(basemodellist, this%m1_idx)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set gwtmodel
    gwtmodel => null()
    mb => GetBaseModelFromList(basemodellist, this%m2_idx)
    select type (mb)
    type is (GwtModelType)
      gwtmodel => mb
    end select
    !
    ! -- Verify that gwf model is of the correct type
    if (.not. associated(gwfmodel)) then
      write (errmsg, '(3a)') 'Problem with GWF-GWT exchange ', trim(this%name), &
        '.  Specified GWF Model does not appear to be of the correct type.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Verify that gwt model is of the correct type
    if (.not. associated(gwtmodel)) then
      write (errmsg, '(3a)') 'Problem with GWF-GWT exchange ', trim(this%name), &
        '.  Specified GWT Model does not appear to be of the correct type.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Tell transport model fmi flows are not read from file
    gwtmodel%fmi%flows_from_file = .false.
    !
    ! -- Set a pointer to the GWF bndlist.  This will allow the transport model
    !    to look through the flow packages and establish a link to GWF flows
    gwtmodel%fmi%gwfbndlist => gwfmodel%bndlist
  end subroutine set_model_pointers

  !> @brief Define the GwfGwt Exchange object
  !<
  subroutine exg_df(this)
    ! -- modules
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    ! -- local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(GwtModelType), pointer :: gwtmodel => null()
    !
    ! -- set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1_idx)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set gwtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2_idx)
    select type (mb)
    type is (GwtModelType)
      gwtmodel => mb
    end select
    !
    ! -- Check to make sure that flow is solved before transport and in a
    !    different IMS solution
    if (gwfmodel%idsoln >= gwtmodel%idsoln) then
      write (errmsg, '(3a)') 'Problem with GWF-GWT exchange ', trim(this%name), &
        '.  The GWF model must be solved by a different IMS than the GWT model. &
        &Furthermore, the IMS specified for GWF must be listed in mfsim.nam &
        &before the IMS for GWT.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Check the discretization and, when the transport model uses a subset
    !    of the flow model cells, build the maps between the two grids.  This
    !    must be done before the transport model arrays are allocated.
    call this%check_discretization(gwfmodel, gwtmodel)
    !
    ! -- Set pointer to flowja
    if (gwtmodel%fmi%igwfmapped == 0) then
      gwtmodel%fmi%gwfflowja => gwfmodel%flowja
      call mem_checkin(gwtmodel%fmi%gwfflowja, &
                       'GWFFLOWJA', gwtmodel%fmi%memoryPath, &
                       'FLOWJA', gwfmodel%memoryPath)
    end if
    !
    ! -- Set the npf flag so that specific discharge is available for
    !    transport calculations if dispersion is active
    if (gwtmodel%indsp > 0) then
      gwfmodel%npf%icalcspdis = 1
    end if
  end subroutine exg_df

  !> @brief Allocate and read
  !<
  subroutine exg_ar(this)
    ! -- modules
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    ! -- local
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(GwtModelType), pointer :: gwtmodel => null()
    ! -- formats
    character(len=*), parameter :: fmtnocorr = &
      "('GWT Model uses a subset of the GWF Model cells for exchange ',a,'.&
      &  Water that flows across the boundary of the transport domain is not&
      &  accounted for unless FLOW_IMBALANCE_CORRECTION is activated in the&
      &  GWT FMI Package.')"
    !
    ! -- set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1_idx)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set gwtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2_idx)
    select type (mb)
    type is (GwtModelType)
      gwtmodel => mb
    end select
    !
    ! -- setup pointers to gwf variables allocated in gwf_ar.  When the
    !    transport model uses a subset of the flow model cells the arrays are
    !    owned by fmi and refilled through the grid maps each time step.
    if (gwtmodel%fmi%igwfmapped /= 0) then
      call gwtmodel%fmi%set_gwf_sources(gwfmodel%x, gwfmodel%npf%sat, &
                                        gwfmodel%npf%spdis, gwfmodel%flowja)
      if (gwtmodel%fmi%iflowerr == 0) then
        write (warnmsg, fmtnocorr) trim(this%name)
        call store_warning(warnmsg)
      end if
    else
      gwtmodel%fmi%gwfhead => gwfmodel%x
      call mem_checkin(gwtmodel%fmi%gwfhead, &
                       'GWFHEAD', gwtmodel%fmi%memoryPath, &
                       'X', gwfmodel%memoryPath)
      gwtmodel%fmi%gwfsat => gwfmodel%npf%sat
      call mem_checkin(gwtmodel%fmi%gwfsat, &
                       'GWFSAT', gwtmodel%fmi%memoryPath, &
                       'SAT', gwfmodel%npf%memoryPath)
      gwtmodel%fmi%gwfspdis => gwfmodel%npf%spdis
      call mem_checkin(gwtmodel%fmi%gwfspdis, &
                       'GWFSPDIS', gwtmodel%fmi%memoryPath, &
                       'SPDIS', gwfmodel%npf%memoryPath)
    end if
    gwtmodel%fmi%igwfspdis = gwfmodel%npf%icalcspdis
    !
    ! -- setup pointers to the flow storage rates. GWF strg arrays are
    !    available after the gwf_ar routine is called.
    if (gwtmodel%inmst > 0) then
      if (gwfmodel%insto > 0) then
        call gwtmodel%fmi%set_gwf_storage(gwfmodel%sto%strgss, &
                                          gwfmodel%sto%strgsy, &
                                          gwfmodel%sto%iusesy)
      end if
    end if
    !
    ! -- Set a pointer to conc in buy
    if (gwfmodel%inbuy > 0) then
      call gwfmodel%buy%set_concentration_pointer(gwtmodel%name, gwtmodel%x, &
                                                  gwtmodel%ibound)
    end if
    !
    ! -- Set a pointer to conc (which could be a temperature) in vsc
    if (gwfmodel%invsc > 0) then
      call gwfmodel%vsc%set_concentration_pointer(gwtmodel%name, gwtmodel%x, &
                                                  gwtmodel%ibound)
    end if
    !
    ! -- transfer the boundary package information from gwf to gwt
    call this%gwfbnd2gwtfmi()
    !
    ! -- if mover package is active, then set a pointer to it's budget object
    if (gwfmodel%inmvr /= 0) then
      gwtmodel%fmi%mvrbudobj => gwfmodel%mvr%budobj
    end if
    !
    ! -- connect Connections
    call this%gwfconn2gwtconn(gwfmodel, gwtmodel)
  end subroutine exg_ar

  !> @brief Check that the transport grid is a subset of the flow grid
  !!
  !! The two models must use the same user grid.  Every cell that is active in
  !! the transport model must also be active in the flow model, but the
  !! transport model may exclude cells that are active in the flow model.
  !<
  subroutine check_discretization(this, gwfmodel, gwtmodel)
    ! -- modules
    use BaseDisModule, only: DisBaseType
    use SimVariablesModule, only: simulation_mode
    use TspAptModule, only: TspAptType
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    type(GwfModelType), pointer :: gwfmodel !< the flow model
    type(GwtModelType), pointer :: gwtmodel !< the transport model
    ! -- local
    class(BndType), pointer :: packobj => null()
    class(DisBaseType), pointer :: gwfdis => null()
    integer(I4B) :: ip, j
    integer(I4B) :: nu, n, nf, nerr
    character(len=20) :: nodestr
    ! -- parameters
    integer(I4B), parameter :: MAXCELLS = 20
    ! -- formats
    character(len=*), parameter :: fmtdiserr = &
      "('GWF and GWT Models do not have the same discretization for exchange&
      & ',a,'.&
      &  GWF Model has ', i0, ' user nodes and the GWT Model has ', i0, '.&
      &  Ensure the discretization packages define the same grid.')"
    character(len=*), parameter :: fmtidomerr = &
      "('Cell ', a, ' is active in the GWT Model but is not active in the GWF&
      & Model for exchange ',a,'.')"
    character(len=*), parameter :: fmtidomsum = &
      "('IDOMAIN for the GWT Model is not a subset of IDOMAIN for the GWF&
      & Model for exchange ',a,'.  ', i0, ' cells are active in the GWT Model&
      & and inactive in the GWF Model.')"
    character(len=*), parameter :: fmtparerr = &
      "('GWT Model uses a subset of the GWF Model cells for exchange ',a,'.&
      &  This is not supported for a parallel simulation.')"
    character(len=*), parameter :: fmtapterr = &
      "('Advanced transport package ',a,' is connected to cell ',a,', which &
      &is not active in the GWT Model.  Cells connected to an advanced &
      &package must be included in the transport domain.')"
    !
    gwfdis => gwfmodel%dis
    !
    ! -- the user grids must be the same
    if (gwtmodel%dis%nodesuser /= gwfdis%nodesuser) then
      write (errmsg, fmtdiserr) trim(this%name), gwfdis%nodesuser, &
        gwtmodel%dis%nodesuser
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- every active transport cell must be active in the flow model
    nerr = 0
    do nu = 1, gwtmodel%dis%nodesuser
      n = gwtmodel%dis%get_nodenumber(nu, 0)
      if (n <= 0) cycle
      nf = gwfdis%get_nodenumber(nu, 0)
      if (nf > 0) cycle
      nerr = nerr + 1
      if (nerr > MAXCELLS) cycle
      call gwtmodel%dis%nodeu_to_string(nu, nodestr)
      write (errmsg, fmtidomerr) trim(adjustl(nodestr)), trim(this%name)
      call store_error(errmsg)
    end do
    if (nerr > 0) then
      write (errmsg, fmtidomsum) trim(this%name), nerr
      call store_error(errmsg)
      call store_error_filename(this%filename)
    end if
    !
    ! -- nothing more to do if the two models use the same cells
    if (gwtmodel%dis%nodes == gwfdis%nodes) return
    !
    if (simulation_mode == 'PARALLEL') then
      write (errmsg, fmtparerr) trim(this%name)
      call store_error(errmsg)
      call store_error_filename(this%filename)
    end if
    !
    call gwtmodel%fmi%map_gwf_grid(gwfdis)
    !
    ! -- an advanced transport package works from flow model cell numbers
    !    rather than from the mapped nodelist
    do ip = 1, gwtmodel%bndlist%Count()
      packobj => GetBndFromList(gwtmodel%bndlist, ip)
      select type (packobj)
      class is (TspAptType)
        do j = 1, packobj%flowbudptr%budterm(packobj%idxbudgwf)%nlist
          nf = packobj%flowbudptr%budterm(packobj%idxbudgwf)%id2(j)
          if (gwtmodel%fmi%gwfnodeinv(nf) > 0) cycle
          call gwfdis%noder_to_string(nf, nodestr)
          write (errmsg, fmtapterr) trim(packobj%packName), &
            trim(adjustl(nodestr))
          call store_error(errmsg)
        end do
      end select
    end do
    if (count_errors() > 0) call store_error_filename(this%filename)
  end subroutine check_discretization

  !> @brief Link GWT connections to GWF connections or exchanges
  !<
  subroutine gwfconn2gwtconn(this, gwfModel, gwtModel)
    ! -- modules
    use SimModule, only: store_error, store_error_filename, count_errors
    use SimVariablesModule, only: iout
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(GwfGwtExchangeType) :: this !< this exchange
    type(GwfModelType), pointer :: gwfModel !< the flow model
    type(GwtModelType), pointer :: gwtModel !< the transport model
    ! -- local
    class(SpatialModelConnectionType), pointer :: conn => null()
    class(*), pointer :: objPtr => null()
    class(GwtGwtConnectionType), pointer :: gwtConn => null()
    class(GwfGwfConnectionType), pointer :: gwfConn => null()
    class(GwfExchangeType), pointer :: gwfExg => null()
    class(GwtExchangeType), pointer :: gwtExg => null()
    integer(I4B) :: ic1, ic2, iex
    integer(I4B) :: gwfConnIdx, gwfExIdx
    logical(LGP) :: areEqual
    !
    ! loop over all connections
    gwtloop: do ic1 = 1, baseconnectionlist%Count()
      !
      conn => get_smc_from_list(baseconnectionlist, ic1)
      if (.not. associated(conn%owner, gwtModel)) cycle gwtloop
      !
      ! start with a GWT conn.
      objPtr => conn
      gwtConn => CastAsGwtGwtConnection(objPtr)
      gwtExg => gwtConn%gwtExchange
      gwfConnIdx = -1
      gwfExIdx = -1
      !
      ! find matching GWF conn. in same list
      gwfloop: do ic2 = 1, baseconnectionlist%Count()
        conn => get_smc_from_list(baseconnectionlist, ic2)
        !
        if (associated(conn%owner, gwfModel)) then
          !
          objPtr => conn
          gwfConn => CastAsGwfGwfConnection(objPtr)
          gwfExg => gwfConn%gwfExchange
          !
          ! A model can have multiple exchanges, even connecting the same two
          ! models. We have a match if
          !  1. gwtgwt%model1 is connected to gwfgwf%model1
          !  2. gwtgwt%model2 is connected to gwfgwf%model2
          !  3. the list of connected nodes (nodem1, nodem2) is equivalent, such
          !     that it contains the same nodes, appearing in the same order in the
          !     exchange data block
          !
          if (gwfExg%v_model1%name /= gwtExg%gwfmodelname1) cycle
          if (gwfExg%v_model2%name /= gwtExg%gwfmodelname2) cycle
          !
          areEqual = same_exchange_cells(gwfExg, gwtExg)
          if (areEqual) then
            ! same cells, same exchange: link and go to next GWT conn.
            write (iout, '(/6a)') 'Linking exchange ', &
              trim(gwtExg%name), ' to ', trim(gwfExg%name), &
              ' (using interface model) for GWT model ', &
              trim(gwtModel%name)
            gwfConnIdx = ic2
            call this%link_connections(gwtConn, gwfConn)
            exit gwfloop
          end if
        end if
      end do gwfloop
      !
      ! fallback option: coupling to old gwfgwf exchange,
      ! the conditions are equal to what is used above
      ! (this will go obsolete at some point)
      if (gwfConnIdx == -1) then
        gwfloopexg: do iex = 1, baseexchangelist%Count()
          gwfExg => GetGwfExchangeFromList(baseexchangelist, iex)
          !
          if (.not. associated(gwfExg)) cycle gwfloopexg
          !
          if (associated(gwfExg%model1, gwfModel) .or. &
              associated(gwfExg%model2, gwfModel)) then
            !
            if (gwfExg%v_model1%name /= gwtExg%gwfmodelname1) cycle
            if (gwfExg%v_model2%name /= gwtExg%gwfmodelname2) cycle
            !
            areEqual = same_exchange_cells(gwfExg, gwtExg)
            if (areEqual) then
              ! link exchange to connection
              write (iout, '(/6a)') 'Linking exchange ', &
                trim(gwtExg%name), ' to ', trim(gwfExg%name), ' for GWT model ', &
                trim(gwtModel%name)
              gwfExIdx = iex
              if (gwtConn%owns_exchange) then
                gwtExg%gwfsimvals => gwfExg%simvals
                call mem_checkin(gwtExg%gwfsimvals, &
                                 'GWFSIMVALS', gwtExg%memoryPath, &
                                 'SIMVALS', gwfExg%memoryPath)
              end if
              !
              !cdl link up mvt to mvr
              if (gwfExg%inmvr > 0) then
                if (gwtConn%owns_exchange .and. &
                    gwtConn%gwtExchange%inmvt > 0) then
                  call gwtExg%mvt%set_pointer_mvrbudobj(gwfExg%mvr%budobj)
                end if
              end if
              !
              if (associated(gwfExg%model2, gwfModel)) gwtConn%exgflowSign = -1
              gwtConn%gwtInterfaceModel%fmi%flows_from_file = .false.
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
        write (errmsg, *) 'Cannot find GWF-GWF exchange when connecting'// &
          ' GWT model ', trim(gwtModel%name), ' with exchange ', &
          trim(gwtExg%name), ' to GWF model ', trim(gwfModel%name), &
          '. Note: GWF-GWF and GWT-GWT need identical exchange data '// &
          '(both in value and order) for the match to succeed.'
        call store_error(errmsg)
      end if
      !
    end do gwtloop
    !
    ! -- report errors
    if (count_errors() > 0) then
      call store_error_filename(this%filename)
    end if
  end subroutine gwfconn2gwtconn

  !> @brief Links a GWT connection to its GWF counterpart
  !<
  subroutine link_connections(this, gwtConn, gwfConn)
    ! -- modules
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(GwfGwtExchangeType) :: this !< this exchange
    class(GwtGwtConnectionType), pointer :: gwtConn !< GWT connection
    class(GwfGwfConnectionType), pointer :: gwfConn !< GWF connection
    !
    !gwtConn%exgflowja => gwfConn%exgflowja
    if (gwtConn%owns_exchange) then
      gwtConn%gwtExchange%gwfsimvals => gwfConn%gwfExchange%simvals
      call mem_checkin(gwtConn%gwtExchange%gwfsimvals, &
                       'GWFSIMVALS', gwtConn%gwtExchange%memoryPath, &
                       'SIMVALS', gwfConn%gwfExchange%memoryPath)
    end if
    !
    !cdl link up mvt to mvr
    if (gwfConn%gwfExchange%inmvr > 0) then
      if (gwtConn%owns_exchange .and. gwtConn%gwtExchange%inmvt > 0) then
        call gwtConn%gwtExchange%mvt%set_pointer_mvrbudobj( &
          gwfConn%gwfExchange%mvr%budobj)
      end if
    end if
    !
    if (associated(gwfConn%gwfExchange%model2, gwfConn%owner)) then
      gwtConn%exgflowSign = -1
    end if
    !
    ! fmi flows are not read from file
    gwtConn%gwtInterfaceModel%fmi%flows_from_file = .false.
    !
    ! set concentration pointer for buoyancy
    ! call gwfConn%gwfInterfaceModel%buy%set_concentration_pointer( &
    !   gwtConn%gwtModel%name, &
    !   gwtConn%conc, &
    !   gwtConn%icbound)
  end subroutine link_connections

  !> @brief Deallocate memory
  !<
  subroutine exg_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    !
    call mem_deallocate(this%m1_idx)
    call mem_deallocate(this%m2_idx)
  end subroutine exg_da

  !> @brief Allocate package scalars
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    !
    call mem_allocate(this%m1_idx, 'M1ID', this%memoryPath)
    call mem_allocate(this%m2_idx, 'M2ID', this%memoryPath)
    this%m1_idx = 0
    this%m2_idx = 0
  end subroutine allocate_scalars

  !> @brief Call routines in FMI that will set pointers to the necessary flow
  !! data
  !<
  subroutine gwfbnd2gwtfmi(this)
    ! -- dummy
    class(GwfGwtExchangeType) :: this
    ! -- local
    integer(I4B) :: ngwfpack, ip, iterm, imover
    class(BaseModelType), pointer :: mb => null()
    type(GwfModelType), pointer :: gwfmodel => null()
    type(GwtModelType), pointer :: gwtmodel => null()
    class(BndType), pointer :: packobj => null()
    integer(I4B) :: imapped
    !
    ! -- set gwfmodel
    mb => GetBaseModelFromList(basemodellist, this%m1_idx)
    select type (mb)
    type is (GwfModelType)
      gwfmodel => mb
    end select
    !
    ! -- set gwtmodel
    mb => GetBaseModelFromList(basemodellist, this%m2_idx)
    select type (mb)
    type is (GwtModelType)
      gwtmodel => mb
    end select
    !
    ! -- Call routines in FMI that will set pointers to the necessary flow
    !    data (SIMVALS and SIMTOMVR) stored within each GWF flow package
    ngwfpack = gwfmodel%bndlist%Count()
    imapped = gwtmodel%fmi%igwfmapped
    iterm = 1
    do ip = 1, ngwfpack
      packobj => GetBndFromList(gwfmodel%bndlist, ip)
      call gwtmodel%fmi%gwfpackages(iterm)%set_pointers( &
        'SIMVALS', &
        packobj%memoryPath, packobj%input_mempath, imapped)
      iterm = iterm + 1
      !
      ! -- If a mover is active for this package, then establish a separate
      !    pointer link for the mover flows stored in SIMTOMVR
      imover = packobj%imover
      if (packobj%isadvpak /= 0) imover = 0
      if (imover /= 0) then
        call gwtmodel%fmi%gwfpackages(iterm)%set_pointers( &
          'SIMTOMVR', &
          packobj%memoryPath, packobj%input_mempath, imapped)
        iterm = iterm + 1
      end if
    end do
  end subroutine gwfbnd2gwtfmi

end module GwfGwtExchangeModule
