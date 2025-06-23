module GweGweConnectionModule
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LINELENGTH, LENCOMPONENTNAME, DZERO, LENBUDTXT
  use CsrUtilsModule, only: getCSRIndex
  use SimModule, only: ustop
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, mem_checkin
  use SpatialModelConnectionModule
  use NumericalModelModule
  use GweModule
  use DisConnExchangeModule
  use GweGweExchangeModule
  use GweInterfaceModelModule
  use SparseModule, only: sparsematrix
  use ConnectionsModule, only: ConnectionsType
  use CellWithNbrsModule, only: GlobalCellType
  use DistVariableModule
  use SimStagesModule
  use MatrixBaseModule
  use AdvSchemeEnumModule

  implicit none
  private

  public :: CastAsGweGweConnection

  !> Connects a GWE model to other GWE models in space. Derives
  !! from NumericalExchangeType so the solution can use it to
  !! fetch the coefficients for this connection.
  !<
  type, public, extends(SpatialModelConnectionType) :: GweGweConnectionType

    class(GweModelType), pointer :: gweModel => null() !< the model for which this connection exists
    class(GweExchangeType), pointer :: gweExchange => null() !< the primary exchange, cast to GWE-GWE
    class(GweInterfaceModelType), pointer :: gweInterfaceModel => null() !< the interface model
    integer(I4B), pointer :: iIfaceAdvScheme => null() !< the advection scheme at the interface:
                                                       !! 0 = upstream, 1 = central, 2 = TVD
    integer(I4B), pointer :: iIfaceXt3d => null() !< XT3D in the interface CND package: 0 = no, 1 = lhs, 2 = rhs
    integer(I4B), pointer :: exgflowSign => null() !< indicates the flow direction of exgflowja
    real(DP), dimension(:), pointer, contiguous :: exgflowjaGwe => null() !< gwe-flowja at the interface (this is a subset of the GWT
                                                                          !! interface model flowja's)

    real(DP), dimension(:), pointer, contiguous :: gwfflowja => null() !< gwfflowja for the interface model
    real(DP), dimension(:), pointer, contiguous :: gwfsat => null() !< gwfsat for the interface model
    real(DP), dimension(:), pointer, contiguous :: gwfhead => null() !< gwfhead for the interface model
    real(DP), dimension(:, :), pointer, contiguous :: gwfspdis => null() !< gwfspdis for the interface model

    real(DP), dimension(:), pointer, contiguous :: conc => null() !< pointer to concentration array
    integer(I4B), dimension(:), pointer, contiguous :: icbound => null() !< store pointer to gwe ibound array

    integer(I4B) :: iout = 0 !< the list file for the interface model

  contains

    procedure, pass(this) :: gweGweConnection_ctor
    generic, public :: construct => gweGweConnection_ctor

    procedure :: exg_ar => gwegwecon_ar
    procedure :: exg_df => gwegwecon_df
    procedure :: exg_rp => gwegwecon_rp
    procedure :: exg_ad => gwegwecon_ad
    procedure :: exg_fc => gwegwecon_fc
    procedure :: exg_da => gwegwecon_da
    procedure :: exg_cq => gwegwecon_cq
    procedure :: exg_bd => gwegwecon_bd
    procedure :: exg_ot => gwegwecon_ot

    ! overriding 'protected'
    procedure :: validateConnection

    ! local stuff
    procedure, private :: allocate_scalars
    procedure, private :: allocate_arrays
    procedure, private :: cfg_dist_vars
    procedure, private :: setGridExtent
    procedure, private :: validateGweExchange
    procedure, private :: setFlowToExchange

  end type GweGweConnectionType

contains

  !> @brief Basic construction of the connection
  !<
  subroutine gweGweConnection_ctor(this, model, gweEx)
    ! modules
    use InputOutputModule, only: openfile
    ! dummy
    class(GweGweConnectionType) :: this !< the connection
    class(NumericalModelType), pointer :: model !< the model owning this connection,
                                                !! this must be a GweModelType
    class(DisConnExchangeType), pointer :: gweEx !< the GWE-GWE exchange the interface model is created for
    ! local
    character(len=LINELENGTH) :: fname
    character(len=LENCOMPONENTNAME) :: name
    class(*), pointer :: objPtr
    logical(LGP) :: write_ifmodel_listfile = .false.

    objPtr => model
    this%gweModel => CastAsGweModel(objPtr)
    objPtr => gweEx
    this%gweExchange => CastAsGweExchange(objPtr)

    if (gweEx%v_model1%is_local .and. gweEx%v_model2%is_local) then
      this%owns_exchange = associated(model, gweEx%model1)
    else
      this%owns_exchange = .true.
    end if

    if (gweEx%v_model1 == model) then
      write (name, '(a,i0)') 'GWECON1_', gweEx%id
    else
      write (name, '(a,i0)') 'GWECON2_', gweEx%id
    end if

    ! .lst file for interface model
    if (write_ifmodel_listfile) then
      fname = trim(name)//'.im.lst'
      call openfile(this%iout, 0, fname, 'LIST', filstat_opt='REPLACE')
      write (this%iout, '(4a)') 'Creating GWE-GWE connection for model ', &
        trim(this%gweModel%name), 'from exchange ', &
        trim(gweEx%name)
    end if

    ! first call base constructor
    call this%SpatialModelConnectionType%spatialConnection_ctor(model, &
                                                                gweEx, &
                                                                name)

    call this%allocate_scalars()
    this%typename = 'GWE-GWE'
    this%iIfaceAdvScheme = 0
    this%iIfaceXt3d = 0
    this%exgflowSign = 1

    allocate (this%gweInterfaceModel)
    this%interface_model => this%gweInterfaceModel

  end subroutine gweGweConnection_ctor

  !> @brief Allocate scalar variables for this connection
  !<
  subroutine allocate_scalars(this)
    ! dummy
    class(GweGweConnectionType) :: this !< the connection

    call mem_allocate(this%iIfaceAdvScheme, 'IADVSCHEME', this%memoryPath)
    call mem_allocate(this%iIfaceXt3d, 'IXT3D', this%memoryPath)
    call mem_allocate(this%exgflowSign, 'EXGFLOWSIGN', this%memoryPath)

  end subroutine allocate_scalars

  !> @brief define the GWE-GWE connection
  !<
  subroutine gwegwecon_df(this)
    ! dummy
    class(GweGweConnectionType) :: this !< the connection
    ! local
    character(len=LENCOMPONENTNAME) :: imName

    ! determine advection scheme (the GWE-GWE exchange
    !    has been read at this point)
    this%iIfaceAdvScheme = this%gweExchange%iAdvScheme

    ! determine xt3d setting on interface
    this%iIfaceXt3d = this%gweExchange%ixt3d

    ! turn off when off in the owning model
    if (this%gweModel%incnd > 0) then
      this%iIfaceXt3d = this%gweModel%cnd%ixt3d
    end if

    ! determine the required size of the interface model grid
    call this%setGridExtent()

    ! now set up the GridConnection
    call this%spatialcon_df()

    ! we have to 'catch up' and create the interface model
    !    here, then the remainder of this routine will be define
    if (this%prim_exchange%v_model1 == this%owner) then
      write (imName, '(a,i0)') 'GWEIM1_', this%gweExchange%id
    else
      write (imName, '(a,i0)') 'GWEIM2_', this%gweExchange%id
    end if
    call this%gweInterfaceModel%gweifmod_cr(imName, &
                                            this%iout, &
                                            this%ig_builder)
    call this%gweInterfaceModel%set_idsoln(this%gweModel%idsoln)
    this%gweInterfaceModel%iAdvScheme = this%iIfaceAdvScheme
    this%gweInterfaceModel%ixt3d = this%iIfaceXt3d
    call this%gweInterfaceModel%model_df()

    call this%cfg_dist_vars()

    call this%allocate_arrays()
    call this%gweInterfaceModel%allocate_fmi()

    ! connect X, RHS, IBOUND, and flowja
    call this%spatialcon_setmodelptrs()

    ! connect pointers (used by BUY)
    this%conc => this%gweInterfaceModel%x
    this%icbound => this%gweInterfaceModel%ibound

    ! add connections from the interface model to solution matrix
    call this%spatialcon_connect()

  end subroutine gwegwecon_df

  !> @brief Configure distributed variables for this interface model
  !<
  subroutine cfg_dist_vars(this)
    ! dummy
    class(GweGweConnectionType) :: this !< the connection

    call this%cfg_dv('X', '', SYNC_NDS, &
                     (/STG_BFR_CON_AR, STG_BFR_EXG_AD, STG_BFR_EXG_CF/))
    call this%cfg_dv('IBOUND', '', SYNC_NDS, (/STG_BFR_CON_AR/))
    call this%cfg_dv('TOP', 'DIS', SYNC_NDS, (/STG_BFR_CON_AR/))
    call this%cfg_dv('BOT', 'DIS', SYNC_NDS, (/STG_BFR_CON_AR/))
    call this%cfg_dv('AREA', 'DIS', SYNC_NDS, (/STG_BFR_CON_AR/))

    if (this%gweInterfaceModel%cnd%idisp > 0) then
      call this%cfg_dv('ALH', 'CND', SYNC_NDS, (/STG_BFR_CON_AR/))
      call this%cfg_dv('ALV', 'CND', SYNC_NDS, (/STG_BFR_CON_AR/))
      call this%cfg_dv('ATH1', 'CND', SYNC_NDS, (/STG_BFR_CON_AR/))
      call this%cfg_dv('ATH2', 'CND', SYNC_NDS, (/STG_BFR_CON_AR/))
      call this%cfg_dv('ATV', 'CND', SYNC_NDS, (/STG_BFR_CON_AR/))
      call this%cfg_dv('KTW', 'CND', SYNC_NDS, (/STG_BFR_CON_AR/))
      call this%cfg_dv('KTS', 'CND', SYNC_NDS, (/STG_BFR_CON_AR/))
    end if
    call this%cfg_dv('GWFHEAD', 'FMI', SYNC_NDS, (/STG_BFR_EXG_AD/))
    call this%cfg_dv('GWFSAT', 'FMI', SYNC_NDS, (/STG_BFR_EXG_AD/))
    call this%cfg_dv('GWFSPDIS', 'FMI', SYNC_NDS, (/STG_BFR_EXG_AD/))
    call this%cfg_dv('GWFFLOWJA', 'FMI', SYNC_CON, (/STG_BFR_EXG_AD/))
    call this%cfg_dv('GWFFLOWJA', 'FMI', SYNC_EXG, (/STG_BFR_EXG_AD/), &
                     exg_var_name='GWFSIMVALS')
    ! fill porosity from est packages, needed for cnd
    if (this%gweModel%incnd > 0 .and. this%gweModel%inest > 0) then
      call this%cfg_dv('POROSITY', 'EST', SYNC_NDS, (/STG_AFT_CON_AR/))
    end if

  end subroutine cfg_dist_vars

  !> @brief Allocate array variables for this connection
  !<
  subroutine allocate_arrays(this)
    class(GweGweConnectionType) :: this !< the connection

    call mem_allocate(this%exgflowjaGwe, this%ig_builder%nrOfBoundaryCells, &
                      'EXGFLOWJAGWT', this%memoryPath)

  end subroutine allocate_arrays

  !> @brief Set required extent of the interface grid from
  !< the configuration
  subroutine setGridExtent(this)
    ! dummy
    class(GweGweConnectionType) :: this !< the connection
    ! local
    logical(LGP) :: hasAdv, hasCnd

    hasAdv = this%gweModel%inadv > 0
    hasCnd = this%gweModel%incnd > 0

    if (hasAdv) then
      if (this%iIfaceAdvScheme == 2) then
        this%exg_stencil_depth = 2
        if (this%gweModel%adv%iadvwt == ADV_SCHEME_TVD) then
          this%int_stencil_depth = 2
        end if
      end if
    end if

    if (hasCnd) then
      if (this%iIfaceXt3d > 0) then
        this%exg_stencil_depth = 2
        if (this%gweModel%cnd%ixt3d > 0) then
          this%int_stencil_depth = 2
        end if
      end if
    end if

  end subroutine setGridExtent

  !> @brief allocate and read/set the connection's data structures
  !<
  subroutine gwegwecon_ar(this)
    class(GweGweConnectionType) :: this !< the connection

    ! check if we can construct an interface model
    ! NB: only makes sense after the models' allocate&read have been
    ! called, which is why we do it here
    call this%validateConnection()

    ! allocate and read base
    call this%spatialcon_ar()

    ! ... and now the interface model
    call this%gweInterfaceModel%model_ar()

    ! set a pointer in the interface model to the gwecommon data
    if (this%gweModel%inest > 0) then
      this%gweInterfaceModel%gwecommon%gwecpw => this%gweModel%gwecommon%gwecpw
      this%gweInterfaceModel%gwecommon%gwerhow => this%gweModel%gwecommon%gwerhow
    end if

    ! set the equation scaling factor in the interface model to that of
    !   underlying GWE model
    if (this%gweModel%incnd > 0) then
      this%gweInterfaceModel%ieqnsclfac = this%gweModel%cnd%eqnsclfac
    end if

    ! AR the movers and obs through the exchange
    if (this%owns_exchange) then
      !cdl implement this when MVT is ready
      !cdl if (this%gweExchange%inmvt > 0) then
      !cdl   call this%gweExchange%mvt%mvt_ar()
      !cdl end if
      if (this%gweExchange%inobs > 0) then
        call this%gweExchange%obs%obs_ar()
      end if
    end if

  end subroutine gwegwecon_ar

  !> @brief validate this connection prior to constructing
  !< the interface model
  subroutine validateConnection(this)
    use SimVariablesModule, only: errmsg
    use SimModule, only: count_errors, store_error
    class(GweGweConnectionType) :: this !< this connection

    ! base validation, the spatial/geometry part
    call this%SpatialModelConnectionType%validateConnection()

    ! we cannot validate this (yet) in parallel mode
    if (.not. this%gweExchange%v_model1%is_local) return
    if (.not. this%gweExchange%v_model2%is_local) return

    ! check specific cross-interface options/values that should be the same
    call this%validateGweExchange()

    ! GWE related matters
    if ((this%gweExchange%gwemodel1%inadv > 0 .and. &
         this%gweExchange%gwemodel2%inadv == 0) .or. &
        (this%gweExchange%gwemodel2%inadv > 0 .and. &
         this%gweExchange%gwemodel1%inadv == 0)) then
      write (errmsg, '(1x,a,a,a)') 'Cannot connect GWE models in exchange ', &
        trim(this%gweExchange%name), ' because one model is configured with ADV &
        &and the other one is not'
      call store_error(errmsg)
    end if

    if ((this%gweExchange%gwemodel1%incnd > 0 .and. &
         this%gweExchange%gwemodel2%incnd == 0) .or. &
        (this%gweExchange%gwemodel2%incnd > 0 .and. &
         this%gweExchange%gwemodel1%incnd == 0)) then
      write (errmsg, '(1x,a,a,a)') 'Cannot connect GWE models in exchange ', &
        trim(this%gweExchange%name), ' because one model is configured with CND &
        &and the other one is not'
      call store_error(errmsg)
    end if

    ! abort on errors
    if (count_errors() > 0) then
      write (errmsg, '(a)') 'Errors occurred while processing exchange(s)'
      call ustop()
    end if

  end subroutine validateConnection

  subroutine gwegwecon_rp(this)
    ! dummy
    class(GweGweConnectionType) :: this !< the connection

    ! call exchange rp routines
    if (this%owns_exchange) then
      call this%gweExchange%exg_rp()
    end if

  end subroutine gwegwecon_rp

  !> @brief Advance this connection
  !<
  subroutine gwegwecon_ad(this)

    class(GweGweConnectionType) :: this !< this connection

    ! recalculate conduction ellipse
    if (this%gweInterfaceModel%incnd > 0) call this%gweInterfaceModel%cnd%cnd_ad()

    if (this%owns_exchange) then
      call this%gweExchange%exg_ad()
    end if

  end subroutine gwegwecon_ad

  subroutine gwegwecon_fc(this, kiter, matrix_sln, rhs_sln, inwtflag)
    ! dummy
    class(GweGweConnectionType) :: this !< the connection
    integer(I4B), intent(in) :: kiter !< the iteration counter
    class(MatrixBaseType), pointer :: matrix_sln !< the system matrix
    real(DP), dimension(:), intent(inout) :: rhs_sln !< global right-hand-side
    integer(I4B), optional, intent(in) :: inwtflag !< newton-raphson flag

    call this%SpatialModelConnectionType%spatialcon_fc( &
      kiter, matrix_sln, rhs_sln, inwtflag)

    ! _fc the movers through the exchange
    if (this%owns_exchange) then
      if (this%gweExchange%inmvt > 0) then
        call this%gweExchange%mvt%mvt_fc(this%gweExchange%gwemodel1%x, &
                                         this%gweExchange%gwemodel2%x)
      end if
    end if

  end subroutine gwegwecon_fc

  subroutine gwegwecon_cq(this, icnvg, isuppress_output, isolnid)
    ! dummy
    class(GweGweConnectionType) :: this !< the connection
    integer(I4B), intent(inout) :: icnvg !< convergence flag
    integer(I4B), intent(in) :: isuppress_output !< suppress output when =1
    integer(I4B), intent(in) :: isolnid !< solution id

    call this%gweInterfaceModel%model_cq(icnvg, isuppress_output)
    call this%setFlowToExchange()

  end subroutine gwegwecon_cq

  !> @brief Set the flows (flowja from interface model) to the
  !< simvals in the exchange, leaving the budget calcution in there
  subroutine setFlowToExchange(this)
    ! modules
    use IndexMapModule
    ! dummy
    class(GweGweConnectionType) :: this !< this connection
    ! local
    integer(I4B) :: i
    class(GweExchangeType), pointer :: gweEx
    type(IndexMapSgnType), pointer :: map

    if (this%owns_exchange) then
      gweEx => this%gweExchange
      map => this%interface_map%exchange_maps(this%interface_map%prim_exg_idx)

      ! use (half of) the exchange map in reverse:
      do i = 1, size(map%src_idx)
        if (map%sign(i) < 0) cycle ! simvals is defined from exg%m1 => exg%m2
        gweEx%simvals(map%src_idx(i)) = &
          this%gweInterfaceModel%flowja(map%tgt_idx(i))
      end do
    end if

  end subroutine setFlowToExchange

  subroutine gwegwecon_bd(this, icnvg, isuppress_output, isolnid)
    ! modules
    use BudgetModule, only: rate_accumulator
    ! dummy
    class(GweGweConnectionType) :: this !< the connection
    integer(I4B), intent(inout) :: icnvg !< convergence flag
    integer(I4B), intent(in) :: isuppress_output !< suppress output when =1
    integer(I4B), intent(in) :: isolnid !< solution id

    ! call exchange budget routine, also calls _bd
    !    for movers.
    if (this%owns_exchange) then
      call this%gweExchange%exg_bd(icnvg, isuppress_output, isolnid)
    end if

  end subroutine gwegwecon_bd

  subroutine gwegwecon_ot(this)
    ! dummy
    class(GweGweConnectionType) :: this !< the connection

    ! call exg_ot() here as it handles all output processing
    !    based on gweExchange%simvals(:), which was correctly
    !    filled from gwegwecon
    if (this%owns_exchange) then
      call this%gweExchange%exg_ot()
    end if

  end subroutine gwegwecon_ot

  !> @brief Validate the exchange, intercepting those
  !! cases where two models have to be connected with an interface
  !! model, where the individual configurations don't allow this
  !!
  !! Stops with error message on config mismatch
  !<
  subroutine validateGweExchange(this)
    ! modules
    use SimVariablesModule, only: errmsg
    use SimModule, only: store_error
    use GweEstModule, only: GweEstType

    ! dummy
    class(GweGweConnectionType) :: this !< this connection

    ! local
    class(GweExchangeType), pointer :: gweEx
    class(*), pointer :: modelPtr
    class(GweModelType), pointer :: gweModel1
    class(GweModelType), pointer :: gweModel2
    type(GweEstType), pointer :: est1, est2
    logical(LGP) :: compatible

    gweEx => this%gweExchange

    ! we cannot validate the remainder (yet) in parallel mode
    if (.not. gweEx%v_model1%is_local) return
    if (.not. gweEx%v_model2%is_local) return

    modelPtr => this%gweExchange%model1
    gweModel1 => CastAsGweModel(modelPtr)
    modelPtr => this%gweExchange%model2
    gweModel2 => CastAsGweModel(modelPtr)

    ! check that EST package usage is the same on both side of the interface
    if ((gweModel1%inest > 0 .and. gweModel2%inest == 0) .or. &
        (gweModel1%inest == 0 .and. gweModel2%inest > 0)) then
      write (errmsg, '(2a)') 'Energy Storage and Transfer package should '// &
        'be enabled/disabled simultaneously in models connected with the '// &
        'interface model for exchange ', &
        trim(gweEx%name)
      call store_error(errmsg)
    end if

    ! conduction options need to be the same in both model
    if ((gweModel1%cnd%ixt3d > 0 .and. gweModel2%cnd%ixt3d == 0) .or. &
        (gweModel1%cnd%ixt3d == 0 .and. gweModel2%cnd%ixt3d > 0)) then
      write (errmsg, '(2a)') 'Use of XT3D to calculate conduction should '// &
        'be the same in both models, either both use XT3D or neither for '// &
        ' exchange '//trim(gweEx%name)
      call store_error(errmsg)
    end if

    ! check compatibility of Energy Storage and Transfer (EST) package
    compatible = .true.
    est1 => gweModel1%est
    est2 => gweModel2%est
    if (est1%rhow /= est2%rhow) compatible = .false.
    if (est1%cpw /= est2%cpw) compatible = .false.
    ! if (est1%nrhospecies /= est2%nrhospecies) compatible = .false.
    ! if (.not. all(buy1%drhodc == buy2%drhodc)) compatible = .false.
    ! if (.not. all(buy1%crhoref == buy2%crhoref)) compatible = .false.
    !if (.not. all(buy1%cauxspeciesname == buy2%cauxspeciesname)) then
    !  compatible = .false.
    !end if
    if (.not. compatible) then
      write (errmsg, '(6a)') 'Energy storage and transfer (EST) packages ', &
        'in model '//trim(gweEx%model1%name), ' and ', &
        trim(gweEx%model2%name), &
        ' should be equivalent to construct an '// &
        ' interface model for exchange ', &
        trim(gweEx%name)
      call store_error(errmsg)
    end if

  end subroutine validateGweExchange

  !> @brief Deallocate all resources
  !<
  subroutine gwegwecon_da(this)
    ! dummy
    class(GweGweConnectionType) :: this !< the connection
    ! local
    logical(LGP) :: isOpen

    ! scalars
    call mem_deallocate(this%iIfaceAdvScheme)
    call mem_deallocate(this%iIfaceXt3d)
    call mem_deallocate(this%exgflowSign)

    ! arrays
    call mem_deallocate(this%exgflowjaGwe)

    ! interface model
    call this%gweInterfaceModel%model_da()
    deallocate (this%gweInterfaceModel)

    ! dealloc base
    call this%spatialcon_da()

    inquire (this%iout, opened=isOpen)
    if (isOpen) then
      close (this%iout)
    end if

    ! we need to deallocate the exchange we own:
    if (this%owns_exchange) then
      call this%gweExchange%exg_da()
    end if

  end subroutine gwegwecon_da

  !> @brief Cast to GweGweConnectionType
  !<
  function CastAsGweGweConnection(obj) result(res)
    implicit none
    ! dummy
    class(*), pointer, intent(inout) :: obj !< object to be cast
    ! return
    class(GweGweConnectionType), pointer :: res !< the GweGweConnection

    res => null()
    if (.not. associated(obj)) return

    select type (obj)
    class is (GweGweConnectionType)
      res => obj
    end select

  end function CastAsGweGweConnection

end module
