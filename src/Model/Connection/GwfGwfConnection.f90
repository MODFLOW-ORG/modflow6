module GwfGwfConnectionModule
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: DZERO, DONE, DEM6, LENVARNAME, &
                             LENCOMPONENTNAME, LENMEMPATH, LINELENGTH
  use CsrUtilsModule, only: getCSRIndex
  use SparseModule, only: sparsematrix
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use SimModule, only: ustop
  use SpatialModelConnectionModule
  use GwfInterfaceModelModule
  use NumericalModelModule
  use GwfModule, only: GwfModelType, CastAsGwfModel
  use DisConnExchangeModule
  use GwfGwfExchangeModule, only: GwfExchangeType, GetGwfExchangeFromList, &
                                  CastAsGwfExchange
  use GwfNpfModule, only: GwfNpfType
  use GwfBuyModule, only: GwfBuyType
  use BaseDisModule, only: DisBaseType
  use ConnectionsModule, only: ConnectionsType
  use CellWithNbrsModule, only: GlobalCellType
  use DistVariableModule
  use SimStagesModule
  use MatrixBaseModule

  implicit none
  private

  public :: CastAsGwfGwfConnection

  !> Connecting a GWF model to other models in space, implements
  !! NumericalExchangeType so the solution can used this object to determine
  !! the coefficients for the coupling between two adjacent models.
  !!
  !! Two connections are created per exchange between model1 and model2:
  !! one to manage the coefficients in the matrix rows for model1, and
  !! the other to do the same for model2.
  !<
  type, public, extends(SpatialModelConnectionType) :: GwfGwfConnectionType

    class(GwfModelType), pointer :: gwfModel => null() !< the model for which this connection exists
    class(GwfExchangeType), pointer :: gwfExchange => null() !< the primary exchange, cast to its concrete type
    class(GwfInterfaceModelType), pointer :: gwfInterfaceModel => null() !< the interface model
    integer(I4B), pointer :: iXt3dOnExchange => null() !< run XT3D on the interface,
                                                       !! 0 = don't, 1 = matrix, 2 = rhs
    integer(I4B) :: iout = 0 !< the list file for the interface model

  contains
    procedure :: gwfGwfConnection_ctor
    generic, public :: construct => gwfGwfConnection_ctor

    ! overriding NumericalExchangeType
    procedure :: exg_df => gwfgwfcon_df
    procedure :: exg_ar => gwfgwfcon_ar
    procedure :: exg_rp => gwfgwfcon_rp
    procedure :: exg_ad => gwfgwfcon_ad
    procedure :: exg_cf => gwfgwfcon_cf
    procedure :: exg_fc => gwfgwfcon_fc
    procedure :: exg_da => gwfgwfcon_da
    procedure :: exg_cq => gwfgwfcon_cq
    procedure :: exg_bd => gwfgwfcon_bd
    procedure :: exg_ot => gwfgwfcon_ot

    ! overriding 'protected'
    procedure :: validateConnection

    ! local stuff
    procedure, private :: cfg_dist_vars
    procedure, private :: allocateScalars
    procedure, private :: setGridExtent
    procedure, private :: validateGwfExchange
    procedure, private :: setFlowToExchange
    procedure, private :: setFlowToModel
    procedure, private :: setNpfEdgeProps

  end type GwfGwfConnectionType

contains

  !> @brief Basic construction of the connection
  !<
  subroutine gwfGwfConnection_ctor(this, model, gwfEx)
    use NumericalModelModule, only: NumericalModelType
    use InputOutputModule, only: openfile
    class(GwfGwfConnectionType) :: this !< the connection
    class(NumericalModelType), pointer :: model !< the model owning this connection,
                                                !! this must of course be a GwfModelType
    class(DisConnExchangeType), pointer :: gwfEx !< the exchange the interface model is created for
    ! local
    character(len=LINELENGTH) :: fname
    character(len=LENCOMPONENTNAME) :: name
    class(*), pointer :: objPtr
    logical(LGP) :: write_ifmodel_listfile = .false.

    objPtr => model
    this%gwfModel => CastAsGwfModel(objPtr)
    objPtr => gwfEx
    this%gwfExchange => CastAsGwfExchange(objPtr)

    if (gwfEx%v_model1%is_local .and. gwfEx%v_model2%is_local) then
      this%owns_exchange = (gwfEx%v_model1 == model)
    else
      this%owns_exchange = .true.
    end if

    if (gwfEx%v_model1 == model) then
      write (name, '(a,i0)') 'GWFCON1_', gwfEx%id
    else
      write (name, '(a,i0)') 'GWFCON2_', gwfEx%id
    end if

    ! .lst file for interface model
    if (write_ifmodel_listfile) then
      fname = trim(name)//'.im.lst'
      call openfile(this%iout, 0, fname, 'LIST', filstat_opt='REPLACE')
      write (this%iout, '(4a)') 'Creating GWF-GWF connection for model ', &
        trim(this%gwfModel%name), ' from exchange ', &
        trim(gwfEx%name)
    end if

    ! first call base constructor
    call this%SpatialModelConnectionType%spatialConnection_ctor(model, &
                                                                gwfEx, &
                                                                name)

    call this%allocateScalars()

    this%typename = 'GWF-GWF'

    ! determine the required size of the interface grid
    call this%setGridExtent()

    allocate (this%gwfInterfaceModel)
    this%interface_model => this%gwfInterfaceModel

  end subroutine gwfGwfConnection_ctor

  !> @brief Define the connection
  !!
  !! This sets up the GridConnection (for creating the
  !! interface grid), creates and defines the interface
  !< model
  subroutine gwfgwfcon_df(this)
    class(GwfGwfConnectionType) :: this !< this connection
    ! local
    character(len=LENCOMPONENTNAME) :: imName !< the interface model's name
    integer(I4B) :: i

    ! this sets up the GridConnection
    call this%spatialcon_df()

    ! Now grid conn is defined, we create the interface model
    ! here, and the remainder of this routine is define.
    ! we basically follow the logic that is present in sln_df()
    if (this%prim_exchange%v_model1 == this%owner) then
      write (imName, '(a,i0)') 'GWFIM1_', this%gwfExchange%id
    else
      write (imName, '(a,i0)') 'GWFIM2_', this%gwfExchange%id
    end if
    call this%gwfInterfaceModel%gwfifm_cr(imName, this%iout, this%ig_builder)
    call this%gwfInterfaceModel%set_idsoln(this%gwfModel%idsoln)
    this%gwfInterfaceModel%npf%satomega = this%gwfModel%npf%satomega
    this%gwfInterfaceModel%npf%ixt3d = this%iXt3dOnExchange
    this%gwfInterfaceModel%npf%icellavg = this%gwfExchange%icellavg
    call this%gwfInterfaceModel%model_df()

    ! Take these settings from the owning model
    this%gwfInterfaceModel%npf%ik22 = this%gwfModel%npf%ik22
    this%gwfInterfaceModel%npf%ik33 = this%gwfModel%npf%ik33
    this%gwfInterfaceModel%npf%iwetdry = this%gwfModel%npf%iwetdry
    this%gwfInterfaceModel%npf%iangle1 = this%gwfModel%npf%iangle1
    this%gwfInterfaceModel%npf%iangle2 = this%gwfModel%npf%iangle2
    this%gwfInterfaceModel%npf%iangle3 = this%gwfModel%npf%iangle3

    call this%cfg_dist_vars()

    if (this%gwfInterfaceModel%npf%ixt3d > 0) then
      this%gwfInterfaceModel%npf%iangle1 = 1
      this%gwfInterfaceModel%npf%iangle2 = 1
      this%gwfInterfaceModel%npf%iangle3 = 1
    end if

    ! set defaults
    do i = 1, size(this%gwfInterfaceModel%npf%angle1)
      this%gwfInterfaceModel%npf%angle1(i) = 0.0_DP
    end do
    do i = 1, size(this%gwfInterfaceModel%npf%angle2)
      this%gwfInterfaceModel%npf%angle2(i) = 0.0_DP
    end do
    do i = 1, size(this%gwfInterfaceModel%npf%angle3)
      this%gwfInterfaceModel%npf%angle3(i) = 0.0_DP
    end do

    ! point X, RHS, IBOUND to connection
    call this%spatialcon_setmodelptrs()

    ! connect interface model to spatial connection
    call this%spatialcon_connect()

  end subroutine gwfgwfcon_df

  !> @brief Configure distributed variables for this interface model
  !<
  subroutine cfg_dist_vars(this)
    class(GwfGwfConnectionType) :: this !< the connection

    call this%cfg_dv('X', '', SYNC_NDS, &
                     (/STG_BFR_CON_AR, STG_BFR_EXG_AD, STG_BFR_EXG_CF/))
    call this%cfg_dv('IBOUND', '', SYNC_NDS, &
                     (/STG_BFR_CON_AR, STG_BFR_EXG_AD, STG_BFR_EXG_CF/))
    call this%cfg_dv('XOLD', '', SYNC_NDS, (/STG_BFR_EXG_AD, STG_BFR_EXG_CF/))
    call this%cfg_dv('ICELLTYPE', 'NPF', SYNC_NDS, (/STG_BFR_CON_AR/))
    call this%cfg_dv('K11', 'NPF', SYNC_NDS, (/STG_BFR_CON_AR/))
    call this%cfg_dv('K22', 'NPF', SYNC_NDS, (/STG_BFR_CON_AR/))
    call this%cfg_dv('K33', 'NPF', SYNC_NDS, (/STG_BFR_CON_AR/))
    if (this%gwfInterfaceModel%npf%iangle1 == 1) then
      call this%cfg_dv('ANGLE1', 'NPF', SYNC_NDS, (/STG_BFR_CON_AR/))
    end if
    if (this%gwfInterfaceModel%npf%iangle2 == 1) then
      call this%cfg_dv('ANGLE2', 'NPF', SYNC_NDS, (/STG_BFR_CON_AR/))
    end if
    if (this%gwfInterfaceModel%npf%iangle3 == 1) then
      call this%cfg_dv('ANGLE3', 'NPF', SYNC_NDS, (/STG_BFR_CON_AR/))
    end if
    if (this%gwfInterfaceModel%npf%iwetdry == 1) then
      call this%cfg_dv('WETDRY', 'NPF', SYNC_NDS, (/STG_BFR_CON_AR/))
    end if
    call this%cfg_dv('TOP', 'DIS', SYNC_NDS, (/STG_BFR_CON_AR/))
    call this%cfg_dv('BOT', 'DIS', SYNC_NDS, (/STG_BFR_CON_AR/))
    call this%cfg_dv('AREA', 'DIS', SYNC_NDS, (/STG_BFR_CON_AR/))
    if (this%gwfInterfaceModel%inbuy > 0) then
      call this%cfg_dv('DENSE', 'BUY', SYNC_NDS, (/STG_BFR_EXG_CF/))
    end if

  end subroutine cfg_dist_vars

  !> @brief Set the required size of the interface grid from
  !< the configuration
  subroutine setGridExtent(this)
    class(GwfGwfConnectionType) :: this !< the connection
    ! local

    this%iXt3dOnExchange = this%gwfExchange%ixt3d
    if (this%iXt3dOnExchange > 0) then
      this%exg_stencil_depth = 2
      if (this%gwfModel%npf%ixt3d > 0) then
        this%int_stencil_depth = 2
      end if
    end if

  end subroutine setGridExtent

  !> @brief allocation of scalars in the connection
  !<
  subroutine allocateScalars(this)
    use MemoryManagerModule, only: mem_allocate
    class(GwfGwfConnectionType) :: this !< the connection
    ! local

    call mem_allocate(this%iXt3dOnExchange, 'IXT3DEXG', this%memoryPath)

  end subroutine allocateScalars

  !> @brief Allocate and read the connection
  !<
  subroutine gwfgwfcon_ar(this)
    use GridConnectionModule, only: GridConnectionType
    class(GwfGwfConnectionType) :: this !< this connection
    ! local

    ! check if we can construct an interface model
    ! NB: only makes sense after the models' allocate&read have been
    ! called, which is why we do it here
    call this%validateConnection()

    ! allocate and read base
    call this%spatialcon_ar()

    ! ... and now the interface model
    call this%gwfInterfaceModel%model_ar()

    ! AR the movers and obs through the exchange
    if (this%owns_exchange) then
      if (this%gwfExchange%inmvr > 0) then
        call this%gwfExchange%mvr%mvr_ar()
      end if
      if (this%gwfExchange%inobs > 0) then
        call this%gwfExchange%obs%obs_ar()
      end if
    end if

  end subroutine gwfgwfcon_ar

  !> @brief Read time varying data when required
  !<
  subroutine gwfgwfcon_rp(this)
    class(GwfGwfConnectionType) :: this !< this connection

    ! Call exchange rp routines
    if (this%owns_exchange) then
      call this%gwfExchange%exg_rp()
    end if
  end subroutine gwfgwfcon_rp

  !> @brief Advance this connection
  !<
  subroutine gwfgwfcon_ad(this)
    class(GwfGwfConnectionType) :: this !< this connection

    ! this triggers the BUY density calculation
    !if (this%gwfInterfaceModel%inbuy > 0) call this%gwfInterfaceModel%buy%buy_ad()

    if (this%owns_exchange) then
      call this%gwfExchange%exg_ad()
    end if

  end subroutine gwfgwfcon_ad

  subroutine gwfgwfcon_cf(this, kiter)
    class(GwfGwfConnectionType) :: this !< this connection
    integer(I4B), intent(in) :: kiter !< the iteration counter

    call this%SpatialModelConnectionType%spatialcon_cf(kiter)

    ! CF the movers through the exchange
    if (this%owns_exchange) then
      if (this%gwfExchange%inmvr > 0) then
        call this%gwfExchange%mvr%xmvr_cf()
      end if
    end if

  end subroutine gwfgwfcon_cf

  !> @brief Write the calculated coefficients into the global
  !< system matrix and the rhs
  subroutine gwfgwfcon_fc(this, kiter, matrix_sln, rhs_sln, inwtflag)
    class(GwfGwfConnectionType) :: this !< this connection
    integer(I4B), intent(in) :: kiter !< the iteration counter
    class(MatrixBaseType), pointer :: matrix_sln !< global system matrix coefficients
    real(DP), dimension(:), intent(inout) :: rhs_sln !< global right-hand-side
    integer(I4B), optional, intent(in) :: inwtflag !< newton-raphson flag
    ! local

    call this%SpatialModelConnectionType%spatialcon_fc( &
      kiter, matrix_sln, rhs_sln, inwtflag)

    ! FC the movers through the exchange; we cannot call
    ! exg_fc() directly because it calculates matrix terms
    if (this%owns_exchange) then
      if (this%gwfExchange%inmvr > 0) then
        call this%gwfExchange%mvr%mvr_fc()
      end if
    end if

  end subroutine gwfgwfcon_fc

  !> @brief Validate this connection
  !! This is called before proceeding to construct
  !! the interface model
  !<
  subroutine validateConnection(this)
    use SimVariablesModule, only: errmsg
    use SimModule, only: count_errors
    class(GwfGwfConnectionType) :: this !< this connection
    ! local

    ! base validation (geometry/spatial)
    call this%SpatialModelConnectionType%validateConnection()
    call this%validateGwfExchange()

    ! abort on errors
    if (count_errors() > 0) then
      write (errmsg, '(a)') 'Errors occurred while processing exchange(s)'
      call ustop(errmsg)
    end if

  end subroutine validateConnection

  !> @brief Validate the exchange, intercepting those
  !! cases where two models have to be connected with an interface
  !! model, where the individual configurations don't allow this
  !!
  !! Stops with error message on config mismatch
  !<
  subroutine validateGwfExchange(this)
    use SimVariablesModule, only: errmsg, simulation_mode
    use SimModule, only: store_error
    use GwfNpfModule, only: GwfNpfType
    class(GwfGwfConnectionType) :: this !< this connection
    ! local
    class(GwfExchangeType), pointer :: gwfEx
    class(*), pointer :: modelPtr
    class(GwfModelType), pointer :: gwfModel1
    class(GwfModelType), pointer :: gwfModel2
    type(GwfBuyType), pointer :: buy1, buy2
    logical(LGP) :: compatible

    gwfEx => this%gwfExchange

    ! GNC not allowed
    if (gwfEx%ingnc /= 0 .and. gwfEx%ixt3d /= 0) then
      write (errmsg, '(2a)') 'Ghost node correction not supported '// &
        'combined with XT3D for exchange ', trim(gwfEx%name)
      call store_error(errmsg)
    end if
    if (gwfEx%ingnc /= 0 .and. simulation_mode == 'PARALLEL') then
      write (errmsg, '(2a)') 'Ghost node correction not supported '// &
        'in parallel run for exchange ', trim(gwfEx%name)
      call store_error(errmsg)
    end if

    ! we cannot validate the remainder (yet) in parallel mode
    if (.not. gwfEx%v_model1%is_local) return
    if (.not. gwfEx%v_model2%is_local) return

    modelPtr => this%gwfExchange%model1
    gwfModel1 => CastAsGwfModel(modelPtr)
    modelPtr => this%gwfExchange%model2
    gwfModel2 => CastAsGwfModel(modelPtr)

    if ((gwfModel1%inbuy > 0 .and. gwfModel2%inbuy == 0) .or. &
        (gwfModel1%inbuy == 0 .and. gwfModel2%inbuy > 0)) then
      write (errmsg, '(2a)') 'Buoyancy package should be enabled/disabled '// &
        'simultaneously in models connected with the '// &
        'interface model for exchange ', &
        trim(gwfEx%name)
      call store_error(errmsg)

    end if

    if (gwfModel1%inbuy > 0 .and. gwfModel2%inbuy > 0) then
      ! does not work with XT3D
      if (this%iXt3dOnExchange > 0) then
        write (errmsg, '(2a)') 'Connecting models with BUY package not '// &
          'allowed with XT3D enabled on exchange ', &
          trim(gwfEx%name)
        call store_error(errmsg)
      end if

      ! check compatibility of buoyancy
      compatible = .true.
      buy1 => gwfModel1%buy
      buy2 => gwfModel2%buy
      if (buy1%iform /= buy2%iform) compatible = .false.
      if (buy1%denseref /= buy2%denseref) compatible = .false.
      if (buy1%nrhospecies /= buy2%nrhospecies) compatible = .false.
      if (.not. all(buy1%drhodc == buy2%drhodc)) compatible = .false.
      if (.not. all(buy1%crhoref == buy2%crhoref)) compatible = .false.
      if (.not. all(buy1%cauxspeciesname == buy2%cauxspeciesname)) then
        compatible = .false.
      end if

      if (.not. compatible) then
        write (errmsg, '(6a)') 'Buoyancy packages in model ', &
          trim(gwfEx%model1%name), ' and ', &
          trim(gwfEx%model2%name), &
          ' should be equivalent to construct an '// &
          ' interface model for exchange ', &
          trim(gwfEx%name)
        call store_error(errmsg)
      end if

    end if

  end subroutine validateGwfExchange

  !> @brief Deallocate all resources
  !<
  subroutine gwfgwfcon_da(this)
    use KindModule, only: LGP
    class(GwfGwfConnectionType) :: this !< this connection
    ! local
    logical(LGP) :: isOpen

    ! scalars
    call mem_deallocate(this%iXt3dOnExchange)

    call this%gwfInterfaceModel%model_da()
    deallocate (this%gwfInterfaceModel)

    call this%spatialcon_da()

    inquire (this%iout, opened=isOpen)
    if (isOpen) then
      close (this%iout)
    end if

    ! we need to deallocate the baseexchange we own:
    if (this%owns_exchange) then
      call this%gwfExchange%exg_da()
    end if

  end subroutine gwfgwfcon_da

  !> @brief Calculate intra-cell flows
  !! The calculation will be dispatched to the interface
  !! model, and then mapped back to real-world cell ids.
  !<
  subroutine gwfgwfcon_cq(this, icnvg, isuppress_output, isolnid)
    class(GwfGwfConnectionType) :: this !< this connection
    integer(I4B), intent(inout) :: icnvg !< convergence flag
    integer(I4B), intent(in) :: isuppress_output !< suppress output when =1
    integer(I4B), intent(in) :: isolnid !< solution id

    call this%gwfInterfaceModel%model_cq(icnvg, isuppress_output)

    call this%setFlowToExchange()

    call this%setFlowToModel()

    !cdl Could we allow GwfExchange to do this instead, using
    !    simvals?
    ! if needed, we add the edge properties to the model's NPF
    ! package for its spdis calculation:
    if (this%gwfModel%npf%icalcspdis == 1) then
      call this%setNpfEdgeProps()
    end if

    ! Add exchange flows to each model flowja diagonal.  This used
    ! to be done in setNpfEdgeProps, but there was a sign issue
    ! and flowja was only updated if icalcspdis was 1 (it should
    ! always be updated.
    if (this%owns_exchange) then
      call this%gwfExchange%gwf_gwf_add_to_flowja()
    end if

  end subroutine gwfgwfcon_cq

  !> @brief Set the flows (flowja from interface model) to the
  !< simvals in the exchange, leaving the budget calcution in there
  subroutine setFlowToExchange(this)
    use IndexMapModule
    class(GwfGwfConnectionType) :: this !< this connection
    ! local
    integer(I4B) :: i
    class(GwfExchangeType), pointer :: gwfEx
    type(IndexMapSgnType), pointer :: map

    if (this%owns_exchange) then
      gwfEx => this%gwfExchange
      map => this%interface_map%exchange_maps(this%interface_map%prim_exg_idx)

      ! use (half of) the exchange map in reverse:
      do i = 1, size(map%src_idx)
        if (map%sign(i) < 0) cycle ! simvals is defined from exg%m1 => exg%m2
        gwfEx%simvals(map%src_idx(i)) = &
          this%gwfInterfaceModel%flowja(map%tgt_idx(i))
      end do
    end if

  end subroutine setFlowToExchange

  !> @brief Set the flows (flowja from the interface model) to
  !< to the model, update the budget
  subroutine setFlowToModel(this)
    class(GwfGwfConnectionType) :: this !< this connection
    ! local
    integer(I4B) :: n, m, ipos, iposLoc
    integer(I4B) :: nLoc, mLoc
    type(ConnectionsType), pointer :: imCon !< interface model connections
    type(GlobalCellType), dimension(:), pointer :: toGlobal !< map interface index to global cell

    ! for readability
    imCon => this%gwfInterfaceModel%dis%con
    toGlobal => this%ig_builder%idxToGlobal

    do n = 1, this%neq
      if (.not. toGlobal(n)%v_model == this%owner) then
        ! only add flows to own model
        cycle
      end if

      nLoc = toGlobal(n)%index

      do ipos = imCon%ia(n) + 1, imCon%ia(n + 1) - 1
        if (imCon%mask(ipos) < 1) cycle ! skip this connection, it's masked so not determined by us

        m = imCon%ja(ipos)
        mLoc = toGlobal(m)%index
        if (toGlobal(m)%v_model == this%owner) then

          ! internal, need to set flowja for n-m
          iposLoc = getCSRIndex(nLoc, mLoc, this%gwfModel%ia, this%gwfModel%ja)

          ! update flowja with correct value
          this%gwfModel%flowja(iposLoc) = this%gwfInterfaceModel%flowja(ipos)

        end if
      end do
    end do

  end subroutine setFlowToModel

  !> @brief Set flowja as edge properties in the model,
  !< so it can be used for e.g. specific discharge calculation
  subroutine setNpfEdgeProps(this)
    class(GwfGwfConnectionType) :: this !< this connection
    ! local
    integer(I4B) :: n, m, ipos, isym
    integer(I4B) :: nLoc, mLoc
    integer(I4B) :: ihc
    real(DP) :: rrate
    real(DP) :: area
    real(DP) :: satThick
    real(DP) :: nx, ny, nz
    real(DP) :: cx, cy, cz
    real(DP) :: conLen
    real(DP) :: dist
    real(DP) :: cl
    logical :: nozee
    type(ConnectionsType), pointer :: imCon !< interface model connections
    class(GwfNpfType), pointer :: imNpf !< interface model npf package
    class(DisBaseType), pointer :: imDis !< interface model discretization
    type(GlobalCellType), dimension(:), pointer :: toGlobal !< map interface index to global cell

    ! for readability
    imDis => this%gwfInterfaceModel%dis
    imCon => this%gwfInterfaceModel%dis%con
    imNpf => this%gwfInterfaceModel%npf
    toGlobal => this%ig_builder%idxToGlobal

    nozee = .false.
    if (imNpf%ixt3d > 0) then
      nozee = imNpf%xt3d%nozee
    end if

    ! loop over flowja in the interface model and set edge properties
    ! for flows crossing the boundary, and set flowja for internal
    ! flows affected by the connection.
    do n = 1, this%neq
      if (.not. toGlobal(n)%v_model == this%owner) then
        ! only add flows to own model
        cycle
      end if

      nLoc = toGlobal(n)%index

      do ipos = imCon%ia(n) + 1, imCon%ia(n + 1) - 1
        if (imCon%mask(ipos) < 1) then
          ! skip this connection, it's masked so not determined by us
          cycle
        end if

        m = imCon%ja(ipos)
        mLoc = toGlobal(m)%index

        if (.not. toGlobal(m)%v_model == this%owner) then
          ! boundary connection, set edge properties
          isym = imCon%jas(ipos)
          ihc = imCon%ihc(isym)
          area = imCon%hwva(isym)
          satThick = imNpf%calcSatThickness(n, m, ihc)
          rrate = this%gwfInterfaceModel%flowja(ipos)

          call imDis%connection_normal(n, m, ihc, nx, ny, nz, ipos)
          call imDis%connection_vector(n, m, nozee, imNpf%sat(n), imNpf%sat(m), &
                                       ihc, cx, cy, cz, conLen)

          if (ihc == 0) then
            ! check if n is below m
            if (nz > 0) rrate = -rrate
          else
            area = area * satThick
          end if

          cl = imCon%cl1(isym)
          if (m < n) then
            cl = imCon%cl2(isym)
          end if
          dist = conLen * cl / (imCon%cl1(isym) + imCon%cl2(isym))
          call this%gwfModel%npf%set_edge_properties(nLoc, ihc, rrate, area, &
                                                     nx, ny, dist)
        end if
      end do
    end do

  end subroutine setNpfEdgeProps

  !> @brief Calculate the budget terms for this connection, this is
  !! dispatched to the GWF-GWF exchange
  subroutine gwfgwfcon_bd(this, icnvg, isuppress_output, isolnid)
    class(GwfGwfConnectionType) :: this !< this connection
    integer(I4B), intent(inout) :: icnvg !< convergence flag
    integer(I4B), intent(in) :: isuppress_output !< suppress output when =1
    integer(I4B), intent(in) :: isolnid !< solution id
    ! local

    ! call exchange budget routine, also calls bd
    ! for movers.
    if (this%owns_exchange) then
      call this%gwfExchange%exg_bd(icnvg, isuppress_output, isolnid)
    end if

  end subroutine gwfgwfcon_bd

  !> @brief Write output for exchange (and calls
  !< save on the budget)
  subroutine gwfgwfcon_ot(this)
    class(GwfGwfConnectionType) :: this !< this connection
    ! local

    ! Call exg_ot() here as it handles all output processing
    ! based on gwfExchange%simvals(:), which was correctly
    ! filled from gwfgwfcon
    if (this%owns_exchange) then
      call this%gwfExchange%exg_ot()
    end if

  end subroutine gwfgwfcon_ot

  !> @brief Cast to GwfGwfConnectionType
  !<
  function CastAsGwfGwfConnection(obj) result(res)
    implicit none
    class(*), pointer, intent(inout) :: obj !< object to be cast
    class(GwfGwfConnectionType), pointer :: res !< the GwfGwfConnection

    res => null()
    if (.not. associated(obj)) return

    select type (obj)
    class is (GwfGwfConnectionType)
      res => obj
    end select
  end function CastAsGwfGwfConnection

end module GwfGwfConnectionModule
