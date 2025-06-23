module GwtGwtConnectionModule
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LINELENGTH, LENCOMPONENTNAME, DZERO, LENBUDTXT
  use CsrUtilsModule, only: getCSRIndex
  use SimModule, only: ustop
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, mem_checkin
  use SpatialModelConnectionModule
  use NumericalModelModule
  use GwtModule
  use DisConnExchangeModule
  use GwtGwtExchangeModule
  use GwtInterfaceModelModule
  use SparseModule, only: sparsematrix
  use ConnectionsModule, only: ConnectionsType
  use CellWithNbrsModule, only: GlobalCellType
  use DistVariableModule
  use SimStagesModule
  use MatrixBaseModule
  use AdvSchemeEnumModule

  implicit none
  private

  public :: CastAsGwtGwtConnection

  !> Connects a GWT model to other GWT models in space. Derives
  !! from NumericalExchangeType so the solution can use it to
  !! fetch the coefficients for this connection.
  !<
  type, public, extends(SpatialModelConnectionType) :: GwtGwtConnectionType

    class(GwtModelType), pointer :: gwtModel => null() !< the model for which this connection exists
    class(GwtExchangeType), pointer :: gwtExchange => null() !< the primary exchange, cast to GWT-GWT
    class(GwtInterfaceModelType), pointer :: gwtInterfaceModel => null() !< the interface model
    integer(I4B), pointer :: iIfaceAdvScheme => null() !< the advection scheme at the interface:
                                                       !! 0 = upstream, 1 = central, 2 = TVD
    integer(I4B), pointer :: iIfaceXt3d => null() !< XT3D in the interface DSP package: 0 = no, 1 = lhs, 2 = rhs
    integer(I4B), pointer :: exgflowSign => null() !< indicates the flow direction of exgflowja
    real(DP), dimension(:), pointer, contiguous :: exgflowjaGwt => null() !< gwt-flowja at the interface (this is a subset of the GWT
                                                                          !! interface model flowja's)

    real(DP), dimension(:), pointer, contiguous :: gwfflowja => null() !< gwfflowja for the interface model
    real(DP), dimension(:), pointer, contiguous :: gwfsat => null() !< gwfsat for the interface model
    real(DP), dimension(:), pointer, contiguous :: gwfhead => null() !< gwfhead for the interface model
    real(DP), dimension(:, :), pointer, contiguous :: gwfspdis => null() !< gwfspdis for the interface model

    real(DP), dimension(:), pointer, contiguous :: conc => null() !< pointer to concentration array
    integer(I4B), dimension(:), pointer, contiguous :: icbound => null() !< store pointer to gwt ibound array

    integer(I4B) :: iout = 0 !< the list file for the interface model

  contains

    procedure :: gwtGwtConnection_ctor
    generic, public :: construct => gwtGwtConnection_ctor

    procedure :: exg_ar => gwtgwtcon_ar
    procedure :: exg_df => gwtgwtcon_df
    procedure :: exg_rp => gwtgwtcon_rp
    procedure :: exg_ad => gwtgwtcon_ad
    procedure :: exg_fc => gwtgwtcon_fc
    procedure :: exg_da => gwtgwtcon_da
    procedure :: exg_cq => gwtgwtcon_cq
    procedure :: exg_bd => gwtgwtcon_bd
    procedure :: exg_ot => gwtgwtcon_ot

    ! overriding 'protected'
    procedure :: validateConnection

    ! local stuff
    procedure, private :: allocate_scalars
    procedure, private :: allocate_arrays
    procedure, private :: cfg_dist_vars
    procedure, private :: setGridExtent
    procedure, private :: setFlowToExchange

  end type GwtGwtConnectionType

contains

  !> @brief Basic construction of the connection
  !<
  subroutine gwtGwtConnection_ctor(this, model, gwtEx)
    use InputOutputModule, only: openfile
    class(GwtGwtConnectionType) :: this !< the connection
    class(NumericalModelType), pointer :: model !< the model owning this connection,
                                                !! this must be a GwtModelType
    class(DisConnExchangeType), pointer :: gwtEx !< the GWT-GWT exchange the interface model is created for
    ! local
    character(len=LINELENGTH) :: fname
    character(len=LENCOMPONENTNAME) :: name
    class(*), pointer :: objPtr
    logical(LGP) :: write_ifmodel_listfile = .false.

    objPtr => model
    this%gwtModel => CastAsGwtModel(objPtr)
    objPtr => gwtEx
    this%gwtExchange => CastAsGwtExchange(objPtr)

    if (gwtEx%v_model1%is_local .and. gwtEx%v_model2%is_local) then
      this%owns_exchange = associated(model, gwtEx%model1)
    else
      this%owns_exchange = .true.
    end if

    if (gwtEx%v_model1 == model) then
      write (name, '(a,i0)') 'GWTCON1_', gwtEx%id
    else
      write (name, '(a,i0)') 'GWTCON2_', gwtEx%id
    end if

    ! .lst file for interface model
    if (write_ifmodel_listfile) then
      fname = trim(name)//'.im.lst'
      call openfile(this%iout, 0, fname, 'LIST', filstat_opt='REPLACE')
      write (this%iout, '(4a)') 'Creating GWT-GWT connection for model ', &
        trim(this%gwtModel%name), 'from exchange ', &
        trim(gwtEx%name)
    end if

    ! first call base constructor
    call this%SpatialModelConnectionType%spatialConnection_ctor(model, &
                                                                gwtEx, &
                                                                name)

    call this%allocate_scalars()
    this%typename = 'GWT-GWT'
    this%iIfaceAdvScheme = 0
    this%iIfaceXt3d = 0
    this%exgflowSign = 1

    allocate (this%gwtInterfaceModel)
    this%interface_model => this%gwtInterfaceModel

  end subroutine gwtGwtConnection_ctor

  !> @brief Allocate scalar variables for this connection
  !<
  subroutine allocate_scalars(this)
    class(GwtGwtConnectionType) :: this !< the connection

    call mem_allocate(this%iIfaceAdvScheme, 'IADVSCHEME', this%memoryPath)
    call mem_allocate(this%iIfaceXt3d, 'IXT3D', this%memoryPath)
    call mem_allocate(this%exgflowSign, 'EXGFLOWSIGN', this%memoryPath)

  end subroutine allocate_scalars

  !> @brief define the GWT-GWT connection
  !<
  subroutine gwtgwtcon_df(this)
    class(GwtGwtConnectionType) :: this !< the connection
    ! local
    character(len=LENCOMPONENTNAME) :: imName

    ! determine advection scheme (the GWT-GWT exchange
    ! has been read at this point)
    this%iIfaceAdvScheme = this%gwtExchange%iAdvScheme
    !
    ! determine xt3d setting on interface
    this%iIfaceXt3d = this%gwtExchange%ixt3d

    ! turn off when off in the owning model
    if (this%gwtModel%indsp > 0) then
      this%iIfaceXt3d = this%gwtModel%dsp%ixt3d
    end if

    ! determine the required size of the interface model grid
    call this%setGridExtent()

    ! now set up the GridConnection
    call this%spatialcon_df()

    ! we have to 'catch up' and create the interface model
    ! here, then the remainder of this routine will be define
    if (this%prim_exchange%v_model1 == this%owner) then
      write (imName, '(a,i0)') 'GWTIM1_', this%gwtExchange%id
    else
      write (imName, '(a,i0)') 'GWTIM2_', this%gwtExchange%id
    end if
    call this%gwtInterfaceModel%gwtifmod_cr(imName, &
                                            this%iout, &
                                            this%ig_builder)
    call this%gwtInterfaceModel%set_idsoln(this%gwtModel%idsoln)
    this%gwtInterfaceModel%iAdvScheme = this%iIfaceAdvScheme
    this%gwtInterfaceModel%ixt3d = this%iIfaceXt3d
    call this%gwtInterfaceModel%model_df()

    call this%cfg_dist_vars()

    call this%allocate_arrays()
    call this%gwtInterfaceModel%allocate_fmi()

    ! connect X, RHS, IBOUND, and flowja
    call this%spatialcon_setmodelptrs()

    ! connect pointers (used by BUY)
    this%conc => this%gwtInterfaceModel%x
    this%icbound => this%gwtInterfaceModel%ibound

    ! add connections from the interface model to solution matrix
    call this%spatialcon_connect()

  end subroutine gwtgwtcon_df

  !> @brief Configure distributed variables for this interface model
  !<
  subroutine cfg_dist_vars(this)
    class(GwtGwtConnectionType) :: this !< the connection

    call this%cfg_dv('X', '', SYNC_NDS, &
                     (/STG_BFR_CON_AR, STG_BFR_EXG_AD, STG_BFR_EXG_CF/))
    call this%cfg_dv('IBOUND', '', SYNC_NDS, (/STG_BFR_CON_AR/))
    call this%cfg_dv('TOP', 'DIS', SYNC_NDS, (/STG_BFR_CON_AR/))
    call this%cfg_dv('BOT', 'DIS', SYNC_NDS, (/STG_BFR_CON_AR/))
    call this%cfg_dv('AREA', 'DIS', SYNC_NDS, (/STG_BFR_CON_AR/))
    if (this%gwtInterfaceModel%dsp%idiffc > 0) then
      call this%cfg_dv('DIFFC', 'DSP', SYNC_NDS, (/STG_BFR_CON_AR/))
    end if
    if (this%gwtInterfaceModel%dsp%idisp > 0) then
      call this%cfg_dv('ALH', 'DSP', SYNC_NDS, (/STG_BFR_CON_AR/))
      call this%cfg_dv('ALV', 'DSP', SYNC_NDS, (/STG_BFR_CON_AR/))
      call this%cfg_dv('ATH1', 'DSP', SYNC_NDS, (/STG_BFR_CON_AR/))
      call this%cfg_dv('ATH2', 'DSP', SYNC_NDS, (/STG_BFR_CON_AR/))
      call this%cfg_dv('ATV', 'DSP', SYNC_NDS, (/STG_BFR_CON_AR/))
    end if
    call this%cfg_dv('GWFHEAD', 'FMI', SYNC_NDS, (/STG_BFR_EXG_AD/))
    call this%cfg_dv('GWFSAT', 'FMI', SYNC_NDS, (/STG_BFR_EXG_AD/))
    call this%cfg_dv('GWFSPDIS', 'FMI', SYNC_NDS, (/STG_BFR_EXG_AD/))
    call this%cfg_dv('GWFFLOWJA', 'FMI', SYNC_CON, (/STG_BFR_EXG_AD/))
    call this%cfg_dv('GWFFLOWJA', 'FMI', SYNC_EXG, (/STG_BFR_EXG_AD/), &
                     exg_var_name='GWFSIMVALS')
    ! fill thetam from mst packages, needed for dsp
    if (this%gwtModel%indsp > 0 .and. this%gwtModel%inmst > 0) then
      call this%cfg_dv('THETAM', 'MST', SYNC_NDS, (/STG_AFT_CON_AR/))
    end if

  end subroutine cfg_dist_vars

  !> @brief Allocate array variables for this connection
  !<
  subroutine allocate_arrays(this)
    class(GwtGwtConnectionType) :: this !< the connection

    call mem_allocate(this%exgflowjaGwt, this%ig_builder%nrOfBoundaryCells, &
                      'EXGFLOWJAGWT', this%memoryPath)

  end subroutine allocate_arrays

  !> @brief Set required extent of the interface grid from
  !< the configuration
  subroutine setGridExtent(this)
    class(GwtGwtConnectionType) :: this !< the connection
    ! local
    logical(LGP) :: hasAdv, hasDsp

    hasAdv = this%gwtModel%inadv > 0
    hasDsp = this%gwtModel%indsp > 0

    if (hasAdv) then
      if (this%iIfaceAdvScheme == 2) then
        this%exg_stencil_depth = 2
        if (this%gwtModel%adv%iadvwt == ADV_SCHEME_TVD) then
          this%int_stencil_depth = 2
        end if
      end if
    end if

    if (hasDsp) then
      if (this%iIfaceXt3d > 0) then
        this%exg_stencil_depth = 2
        if (this%gwtModel%dsp%ixt3d > 0) then
          this%int_stencil_depth = 2
        end if
      end if
    end if

  end subroutine setGridExtent

  !> @brief allocate and read/set the connection's data structures
  !<
  subroutine gwtgwtcon_ar(this)
    class(GwtGwtConnectionType) :: this !< the connection

    ! check if we can construct an interface model
    ! NB: only makes sense after the models' allocate&read have been
    ! called, which is why we do it here
    call this%validateConnection()

    ! allocate and read base
    call this%spatialcon_ar()

    ! ... and now the interface model
    call this%gwtInterfaceModel%model_ar()

    ! AR the movers and obs through the exchange
    if (this%owns_exchange) then
      !cdl implement this when MVT is ready
      !cdl if (this%gwtExchange%inmvt > 0) then
      !cdl   call this%gwtExchange%mvt%mvt_ar()
      !cdl end if
      if (this%gwtExchange%inobs > 0) then
        call this%gwtExchange%obs%obs_ar()
      end if
    end if

  end subroutine gwtgwtcon_ar

  !> @brief validate this connection prior to constructing
  !< the interface model
  subroutine validateConnection(this)
    use SimVariablesModule, only: errmsg
    use SimModule, only: count_errors, store_error
    class(GwtGwtConnectionType) :: this !< this connection

    ! base validation, the spatial/geometry part
    call this%SpatialModelConnectionType%validateConnection()

    ! we cannot validate this (yet) in parallel mode
    if (.not. this%gwtExchange%v_model1%is_local) return
    if (.not. this%gwtExchange%v_model2%is_local) return

    ! GWT related matters
    if ((this%gwtExchange%gwtmodel1%inadv > 0 .and. &
         this%gwtExchange%gwtmodel2%inadv == 0) .or. &
        (this%gwtExchange%gwtmodel2%inadv > 0 .and. &
         this%gwtExchange%gwtmodel1%inadv == 0)) then
      write (errmsg, '(a,a,a)') 'Cannot connect GWT models in exchange ', &
        trim(this%gwtExchange%name), ' because one model is configured with ADV &
        &and the other one is not'
      call store_error(errmsg)
    end if

    if ((this%gwtExchange%gwtmodel1%indsp > 0 .and. &
         this%gwtExchange%gwtmodel2%indsp == 0) .or. &
        (this%gwtExchange%gwtmodel2%indsp > 0 .and. &
         this%gwtExchange%gwtmodel1%indsp == 0)) then
      write (errmsg, '(a,a,a)') 'Cannot connect GWT models in exchange ', &
        trim(this%gwtExchange%name), ' because one model is configured with DSP &
        &and the other one is not'
      call store_error(errmsg)
    end if

    ! abort on errors
    if (count_errors() > 0) then
      write (errmsg, '(a)') 'Errors occurred while processing exchange(s)'
      call ustop()
    end if

  end subroutine validateConnection

  subroutine gwtgwtcon_rp(this)
    class(GwtGwtConnectionType) :: this !< the connection

    ! Call exchange rp routines
    if (this%owns_exchange) then
      call this%gwtExchange%exg_rp()
    end if

  end subroutine gwtgwtcon_rp

  !> @brief Advance this connection
  !<
  subroutine gwtgwtcon_ad(this)
    class(GwtGwtConnectionType) :: this !< this connection

    ! recalculate dispersion ellipse
    if (this%gwtInterfaceModel%indsp > 0) call this%gwtInterfaceModel%dsp%dsp_ad()

    if (this%owns_exchange) then
      call this%gwtExchange%exg_ad()
    end if

  end subroutine gwtgwtcon_ad

  subroutine gwtgwtcon_fc(this, kiter, matrix_sln, rhs_sln, inwtflag)
    class(GwtGwtConnectionType) :: this !< the connection
    integer(I4B), intent(in) :: kiter !< the iteration counter
    class(MatrixBaseType), pointer :: matrix_sln !< the system matrix
    real(DP), dimension(:), intent(inout) :: rhs_sln !< global right-hand-side
    integer(I4B), optional, intent(in) :: inwtflag !< newton-raphson flag
    !

    call this%SpatialModelConnectionType%spatialcon_fc( &
      kiter, matrix_sln, rhs_sln, inwtflag)
    !
    ! FC the movers through the exchange
    if (this%owns_exchange) then
      if (this%gwtExchange%inmvt > 0) then
        call this%gwtExchange%mvt%mvt_fc(this%gwtExchange%gwtmodel1%x, &
                                         this%gwtExchange%gwtmodel2%x)
      end if
    end if

  end subroutine gwtgwtcon_fc

  subroutine gwtgwtcon_cq(this, icnvg, isuppress_output, isolnid)
    class(GwtGwtConnectionType) :: this !< the connection
    integer(I4B), intent(inout) :: icnvg !< convergence flag
    integer(I4B), intent(in) :: isuppress_output !< suppress output when =1
    integer(I4B), intent(in) :: isolnid !< solution id

    call this%gwtInterfaceModel%model_cq(icnvg, isuppress_output)
    call this%setFlowToExchange()

  end subroutine gwtgwtcon_cq

  !> @brief Set the flows (flowja from interface model) to the
  !< simvals in the exchange, leaving the budget calcution in there
  subroutine setFlowToExchange(this)
    use IndexMapModule
    class(GwtGwtConnectionType) :: this !< this connection
    ! local
    integer(I4B) :: i
    class(GwtExchangeType), pointer :: gwtEx
    type(IndexMapSgnType), pointer :: map

    if (this%owns_exchange) then
      gwtEx => this%gwtExchange
      map => this%interface_map%exchange_maps(this%interface_map%prim_exg_idx)

      ! use (half of) the exchange map in reverse:
      do i = 1, size(map%src_idx)
        if (map%sign(i) < 0) cycle ! simvals is defined from exg%m1 => exg%m2
        gwtEx%simvals(map%src_idx(i)) = &
          this%gwtInterfaceModel%flowja(map%tgt_idx(i))
      end do
    end if

  end subroutine setFlowToExchange

  subroutine gwtgwtcon_bd(this, icnvg, isuppress_output, isolnid)
    use BudgetModule, only: rate_accumulator
    class(GwtGwtConnectionType) :: this !< the connection
    integer(I4B), intent(inout) :: icnvg !< convergence flag
    integer(I4B), intent(in) :: isuppress_output !< suppress output when =1
    integer(I4B), intent(in) :: isolnid !< solution id

    ! call exchange budget routine, also calls bd
    ! for movers.
    if (this%owns_exchange) then
      call this%gwtExchange%exg_bd(icnvg, isuppress_output, isolnid)
    end if

  end subroutine gwtgwtcon_bd

  subroutine gwtgwtcon_ot(this)
    class(GwtGwtConnectionType) :: this !< the connection

    ! Call exg_ot() here as it handles all output processing
    ! based on gwtExchange%simvals(:), which was correctly
    ! filled from gwtgwtcon
    if (this%owns_exchange) then
      call this%gwtExchange%exg_ot()
    end if

  end subroutine gwtgwtcon_ot

  subroutine gwtgwtcon_da(this)
    class(GwtGwtConnectionType) :: this !< the connection
    ! local
    logical(LGP) :: isOpen

    ! scalars
    call mem_deallocate(this%iIfaceAdvScheme)
    call mem_deallocate(this%iIfaceXt3d)
    call mem_deallocate(this%exgflowSign)

    ! arrays
    call mem_deallocate(this%exgflowjaGwt)

    ! interface model
    call this%gwtInterfaceModel%model_da()
    deallocate (this%gwtInterfaceModel)

    ! dealloc base
    call this%spatialcon_da()

    inquire (this%iout, opened=isOpen)
    if (isOpen) then
      close (this%iout)
    end if

    ! we need to deallocate the exchange we own:
    if (this%owns_exchange) then
      call this%gwtExchange%exg_da()
    end if

  end subroutine gwtgwtcon_da

  !> @brief Cast to GwtGwtConnectionType
  !<
  function CastAsGwtGwtConnection(obj) result(res)
    implicit none
    class(*), pointer, intent(inout) :: obj !< object to be cast
    class(GwtGwtConnectionType), pointer :: res !< the GwtGwtConnection

    res => null()
    if (.not. associated(obj)) return

    select type (obj)
    class is (GwtGwtConnectionType)
      res => obj
    end select
  end function CastAsGwtGwtConnection

end module
