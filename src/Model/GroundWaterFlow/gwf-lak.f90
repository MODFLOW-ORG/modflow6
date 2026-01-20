module LakModule
  !
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENBOUNDNAME, LENTIMESERIESNAME, &
                             IWETLAKE, MAXADPIT, &
                             DZERO, DPREC, DEM30, DEM9, DEM6, DEM5, &
                             DEM4, DEM2, DEM1, DHALF, DP7, DP999, DONE, &
                             DTWO, DPI, DTHREE, DEIGHT, DTEN, DHUNDRED, DEP20, &
                             DONETHIRD, DTWOTHIRDS, DFIVETHIRDS, &
                             DGRAVITY, DCD, &
                             NAMEDBOUNDFLAG, LENFTYPE, LENPACKAGENAME, &
                             LENPAKLOC, DNODATA, &
                             TABLEFT, TABCENTER, TABRIGHT, &
                             TABSTRING, TABUCSTRING, TABINTEGER, TABREAL
  use MemoryManagerModule, only: mem_allocate, mem_reallocate, mem_setptr, &
                                 mem_deallocate
  use MemoryHelperModule, only: create_mem_path
  use SmoothingModule, only: sQuadraticSaturation, sQSaturation, &
                             sQuadraticSaturationDerivative, &
                             sQSaturationDerivative
  use BndModule, only: BndType
  use BudgetObjectModule, only: BudgetObjectType, budgetobject_cr
  use TableModule, only: TableType, table_cr
  use ObserveModule, only: ObserveType
  use ObsModule, only: ObsType
  use GeomUtilModule, only: get_node
  use InputOutputModule, only: URWORD, extract_idnum_or_bndname
  use BaseDisModule, only: DisBaseType
  use SimModule, only: count_errors, store_error, store_error_unit, &
                       deprecation_warning
  use MathUtilModule, only: is_close
  use BlockParserModule, only: BlockParserType
  use BaseDisModule, only: DisBaseType
  use SimVariablesModule, only: errmsg, warnmsg
  use MatrixBaseModule
  !
  implicit none
  !
  private
  public :: LakType
  public :: lak_create
  !
  character(len=LENFTYPE) :: ftype = 'LAK'
  character(len=LENPACKAGENAME) :: text = '             LAK'
  !
  type LakTabType
    real(DP), dimension(:), pointer, contiguous :: tabstage => null()
    real(DP), dimension(:), pointer, contiguous :: tabvolume => null()
    real(DP), dimension(:), pointer, contiguous :: tabsarea => null()
    real(DP), dimension(:), pointer, contiguous :: tabwarea => null()
  end type LakTabType
  !
  type, extends(BndType) :: LakType
    ! -- scalars
    ! -- characters
    character(len=16), dimension(:), pointer, contiguous :: clakbudget => NULL()
    character(len=16), dimension(:), pointer, contiguous :: cauxcbc => NULL()
    ! -- control variables
    ! -- integers
    integer(I4B), pointer :: iprhed => null()
    integer(I4B), pointer :: istageout => null()
    integer(I4B), pointer :: ibudgetout => null()
    integer(I4B), pointer :: ibudcsv => null()
    integer(I4B), pointer :: ipakcsv => null()
    integer(I4B), pointer :: cbcauxitems => NULL()
    integer(I4B), pointer :: nlakes => NULL()
    integer(I4B), pointer :: noutlets => NULL()
    integer(I4B), pointer :: ntables => NULL()
    real(DP), pointer :: convlength => NULL()
    real(DP), pointer :: convtime => NULL()
    real(DP), pointer :: outdmax => NULL()
    integer(I4B), pointer :: igwhcopt => NULL()
    integer(I4B), pointer :: iconvchk => NULL()
    integer(I4B), pointer :: iconvresidchk => NULL()
    integer(I4B), pointer :: maxlakit => NULL() !< maximum number of iterations in LAK solve
    real(DP), pointer :: surfdep => NULL()
    real(DP), pointer :: dmaxchg => NULL()
    real(DP), pointer :: delh => NULL()
    real(DP), pointer :: pdmax => NULL()
    integer(I4B), pointer :: check_attr => NULL()
    ! -- for budgets
    integer(I4B), pointer :: bditems => NULL()
    ! -- vectors
    ! -- lake data
    integer(I4B), dimension(:), pointer, contiguous :: nlakeconn => null()
    integer(I4B), dimension(:), pointer, contiguous :: idxlakeconn => null()
    integer(I4B), dimension(:), pointer, contiguous :: ntabrow => null()
    real(DP), dimension(:), pointer, contiguous :: strt => null()
    real(DP), dimension(:), pointer, contiguous :: laketop => null()
    real(DP), dimension(:), pointer, contiguous :: lakebot => null()
    real(DP), dimension(:), pointer, contiguous :: sareamax => null()
    character(len=LENBOUNDNAME), dimension(:), pointer, &
      contiguous :: lakename => null()
    character(len=8), dimension(:), pointer, contiguous :: status => null()
    real(DP), dimension(:), pointer, contiguous :: avail => null()
    real(DP), dimension(:), pointer, contiguous :: lkgwsink => null()
    real(DP), dimension(:), pointer, contiguous :: stage => null()
    real(DP), dimension(:), pointer, contiguous :: rainfall => null()
    real(DP), dimension(:), pointer, contiguous :: evaporation => null()
    real(DP), dimension(:), pointer, contiguous :: runoff => null()
    real(DP), dimension(:), pointer, contiguous :: inflow => null()
    real(DP), dimension(:), pointer, contiguous :: withdrawal => null()
    real(DP), dimension(:, :), pointer, contiguous :: lauxvar => null()
    !
    ! -- table data
    integer(I4B), dimension(:), pointer, contiguous :: ialaktab => null()
    real(DP), dimension(:), pointer, contiguous :: tabstage => null()
    real(DP), dimension(:), pointer, contiguous :: tabvolume => null()
    real(DP), dimension(:), pointer, contiguous :: tabsarea => null()
    real(DP), dimension(:), pointer, contiguous :: tabwarea => null()
    !
    ! -- lake solution data
    integer(I4B), dimension(:), pointer, contiguous :: ncncvr => null()
    real(DP), dimension(:), pointer, contiguous :: surfin => null()
    real(DP), dimension(:), pointer, contiguous :: surfout => null()
    real(DP), dimension(:), pointer, contiguous :: surfout1 => null()
    real(DP), dimension(:), pointer, contiguous :: precip => null()
    real(DP), dimension(:), pointer, contiguous :: precip1 => null()
    real(DP), dimension(:), pointer, contiguous :: evap => null()
    real(DP), dimension(:), pointer, contiguous :: evap1 => null()
    real(DP), dimension(:), pointer, contiguous :: evapo => null()
    real(DP), dimension(:), pointer, contiguous :: withr => null()
    real(DP), dimension(:), pointer, contiguous :: withr1 => null()
    real(DP), dimension(:), pointer, contiguous :: flwin => null()
    real(DP), dimension(:), pointer, contiguous :: flwiter => null()
    real(DP), dimension(:), pointer, contiguous :: flwiter1 => null()
    real(DP), dimension(:), pointer, contiguous :: seep => null()
    real(DP), dimension(:), pointer, contiguous :: seep1 => null()
    real(DP), dimension(:), pointer, contiguous :: seep0 => null()
    real(DP), dimension(:), pointer, contiguous :: stageiter => null()
    real(DP), dimension(:), pointer, contiguous :: chterm => null()
    !
    ! -- lake convergence
    integer(I4B), dimension(:), pointer, contiguous :: iseepc => null()
    integer(I4B), dimension(:), pointer, contiguous :: idhc => null()
    real(DP), dimension(:), pointer, contiguous :: en1 => null()
    real(DP), dimension(:), pointer, contiguous :: en2 => null()
    real(DP), dimension(:), pointer, contiguous :: r1 => null()
    real(DP), dimension(:), pointer, contiguous :: r2 => null()
    real(DP), dimension(:), pointer, contiguous :: dh0 => null()
    real(DP), dimension(:), pointer, contiguous :: s0 => null()
    real(DP), dimension(:), pointer, contiguous :: qgwf0 => null()
    !
    ! -- lake connection data
    integer(I4B), dimension(:), pointer, contiguous :: imap => null()
    integer(I4B), dimension(:), pointer, contiguous :: cellid => null()
    integer(I4B), dimension(:), pointer, contiguous :: nodesontop => null()
    integer(I4B), dimension(:), pointer, contiguous :: ictype => null()
    real(DP), dimension(:), pointer, contiguous :: bedleak => null()
    real(DP), dimension(:), pointer, contiguous :: belev => null()
    real(DP), dimension(:), pointer, contiguous :: telev => null()
    real(DP), dimension(:), pointer, contiguous :: connlength => null()
    real(DP), dimension(:), pointer, contiguous :: connwidth => null()
    real(DP), dimension(:), pointer, contiguous :: sarea => null()
    real(DP), dimension(:), pointer, contiguous :: warea => null()
    real(DP), dimension(:), pointer, contiguous :: satcond => null()
    real(DP), dimension(:), pointer, contiguous :: simcond => null()
    real(DP), dimension(:), pointer, contiguous :: simlakgw => null()
    !
    ! -- lake outlet data
    integer(I4B), dimension(:), pointer, contiguous :: lakein => null()
    integer(I4B), dimension(:), pointer, contiguous :: lakeout => null()
    integer(I4B), dimension(:), pointer, contiguous :: iouttype => null()
    real(DP), dimension(:), pointer, contiguous :: outrate => null()
    real(DP), dimension(:), pointer, contiguous :: outinvert => null()
    real(DP), dimension(:), pointer, contiguous :: outwidth => null()
    real(DP), dimension(:), pointer, contiguous :: outrough => null()
    real(DP), dimension(:), pointer, contiguous :: outslope => null()
    real(DP), dimension(:), pointer, contiguous :: simoutrate => null()
    !
    ! -- lake output data
    real(DP), dimension(:), pointer, contiguous :: qauxcbc => null()
    real(DP), dimension(:), pointer, contiguous :: dbuff => null()
    real(DP), dimension(:), pointer, contiguous :: qleak => null()
    real(DP), dimension(:), pointer, contiguous :: qsto => null()
    !
    ! -- pointer to gwf iss and gwf hk
    integer(I4B), pointer :: gwfiss => NULL()
    real(DP), dimension(:), pointer, contiguous :: gwfk11 => NULL()
    real(DP), dimension(:), pointer, contiguous :: gwfk33 => NULL()
    real(DP), dimension(:), pointer, contiguous :: gwfsat => NULL()
    integer(I4B), pointer :: gwfik33 => NULL()
    !
    ! -- package x, xold, and ibound
    integer(I4B), dimension(:), pointer, contiguous :: iboundpak => null() !package ibound
    real(DP), dimension(:), pointer, contiguous :: xnewpak => null() !package x vector
    real(DP), dimension(:), pointer, contiguous :: xoldpak => null() !package xold vector
    !
    ! -- lake budget object
    type(BudgetObjectType), pointer :: budobj => null()
    !
    ! -- lake table objects
    type(TableType), pointer :: stagetab => null()
    type(TableType), pointer :: pakcsvtab => null()
    !
    ! -- density variables
    integer(I4B), pointer :: idense
    real(DP), dimension(:, :), pointer, contiguous :: denseterms => null()
    !
    ! -- viscosity variables
    real(DP), dimension(:, :), pointer, contiguous :: viscratios => null() !< viscosity ratios (1: lak vsc ratio; 2: gwf vsc ratio)
    !
    ! -- type bound procedures

  contains

    procedure :: lak_allocate_scalars
    procedure :: lak_allocate_arrays
    procedure :: bnd_options => lak_options
    procedure :: read_dimensions => lak_read_dimensions
    procedure :: read_initial_attr => lak_read_initial_attr
    procedure :: set_pointers => lak_set_pointers
    procedure :: bnd_ar => lak_ar
    procedure :: bnd_rp => lak_rp
    procedure :: bnd_ad => lak_ad
    procedure :: bnd_cf => lak_cf
    procedure :: bnd_fc => lak_fc
    procedure :: bnd_fn => lak_fn
    procedure :: bnd_cc => lak_cc
    procedure :: bnd_cq => lak_cq
    procedure :: bnd_ot_model_flows => lak_ot_model_flows
    procedure :: bnd_ot_package_flows => lak_ot_package_flows
    procedure :: bnd_ot_dv => lak_ot_dv
    procedure :: bnd_ot_bdsummary => lak_ot_bdsummary
    procedure :: bnd_da => lak_da
    procedure :: define_listlabel
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => lak_obs_supported
    procedure, public :: bnd_df_obs => lak_df_obs
    procedure, public :: bnd_rp_obs => lak_rp_obs
    procedure, public :: bnd_bd_obs => lak_bd_obs
    ! -- private procedures
    procedure, private :: lak_read_lakes
    procedure, private :: lak_read_lake_connections
    procedure, private :: lak_read_outlets
    procedure, private :: lak_read_tables
    procedure, private :: lak_read_table
    procedure, private :: lak_check_valid
    procedure, private :: lak_set_stressperiod
    procedure, private :: lak_set_attribute_error
    procedure, private :: lak_bound_update
    procedure, private :: lak_calculate_sarea
    procedure, private :: lak_calculate_warea
    procedure, private :: lak_calculate_conn_warea
    procedure, public :: lak_calculate_vol
    procedure, private :: lak_calculate_conductance
    procedure, private :: lak_calculate_cond_head
    procedure, private :: lak_calculate_conn_conductance
    procedure, private :: lak_calculate_exchange
    procedure, private :: lak_calculate_conn_exchange
    procedure, private :: lak_estimate_conn_exchange
    procedure, private :: lak_calculate_storagechange
    procedure, private :: lak_calculate_rainfall
    procedure, private :: lak_calculate_runoff
    procedure, private :: lak_calculate_inflow
    procedure, private :: lak_calculate_external
    procedure, private :: lak_calculate_withdrawal
    procedure, private :: lak_calculate_evaporation
    procedure, private :: lak_calculate_outlet_inflow
    procedure, private :: lak_calculate_outlet_outflow
    procedure, private :: lak_get_internal_inlet
    procedure, private :: lak_get_internal_outlet
    procedure, private :: lak_get_external_outlet
    procedure, private :: lak_get_internal_mover
    procedure, private :: lak_get_external_mover
    procedure, private :: lak_get_outlet_tomover
    procedure, private :: lak_accumulate_chterm
    procedure, private :: lak_vol2stage
    procedure, private :: lak_solve
    procedure, private :: lak_bisection
    procedure, private :: lak_calculate_available
    procedure, private :: lak_calculate_residual
    procedure, private :: lak_linear_interpolation
    procedure, private :: lak_setup_budobj
    procedure, private :: lak_fill_budobj
    procedure, private :: laktables_to_vectors
    ! -- table
    procedure, private :: lak_setup_tableobj
    ! -- density
    procedure :: lak_activate_density
    procedure, private :: lak_calculate_density_exchange
    ! -- viscosity
    procedure :: lak_activate_viscosity
  end type LakType

contains

  !> @brief Create a new LAK Package and point bndobj to the new package
  !<
  subroutine lak_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname)
    ! -- dummy
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    ! -- local
    type(LakType), pointer :: lakobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (lakobj)
    packobj => lakobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype)
    packobj%text = text
    !
    ! -- allocate scalars
    call lakobj%lak_allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()
    !
    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 3
    packobj%iscloc = 0 ! not supported
    packobj%isadvpak = 1
    packobj%ictMemPath = create_mem_path(namemodel, 'NPF')
  end subroutine lak_create

  !> @brief Allocate scalar members
  !<
  subroutine lak_allocate_scalars(this)
    ! -- dummy
    class(LakType), intent(inout) :: this
    !
    ! -- call standard BndType allocate scalars
    call this%BndType%allocate_scalars()
    !
    ! -- allocate the object and assign values to object variables
    call mem_allocate(this%iprhed, 'IPRHED', this%memoryPath)
    call mem_allocate(this%istageout, 'ISTAGEOUT', this%memoryPath)
    call mem_allocate(this%ibudgetout, 'IBUDGETOUT', this%memoryPath)
    call mem_allocate(this%ibudcsv, 'IBUDCSV', this%memoryPath)
    call mem_allocate(this%ipakcsv, 'IPAKCSV', this%memoryPath)
    call mem_allocate(this%nlakes, 'NLAKES', this%memoryPath)
    call mem_allocate(this%noutlets, 'NOUTLETS', this%memoryPath)
    call mem_allocate(this%ntables, 'NTABLES', this%memoryPath)
    call mem_allocate(this%convlength, 'CONVLENGTH', this%memoryPath)
    call mem_allocate(this%convtime, 'CONVTIME', this%memoryPath)
    call mem_allocate(this%outdmax, 'OUTDMAX', this%memoryPath)
    call mem_allocate(this%igwhcopt, 'IGWHCOPT', this%memoryPath)
    call mem_allocate(this%iconvchk, 'ICONVCHK', this%memoryPath)
    call mem_allocate(this%iconvresidchk, 'ICONVRESIDCHK', this%memoryPath)
    call mem_allocate(this%maxlakit, 'MAXLAKIT', this%memoryPath)
    call mem_allocate(this%surfdep, 'SURFDEP', this%memoryPath)
    call mem_allocate(this%dmaxchg, 'DMAXCHG', this%memoryPath)
    call mem_allocate(this%delh, 'DELH', this%memoryPath)
    call mem_allocate(this%pdmax, 'PDMAX', this%memoryPath)
    call mem_allocate(this%check_attr, 'CHECK_ATTR', this%memoryPath)
    call mem_allocate(this%bditems, 'BDITEMS', this%memoryPath)
    call mem_allocate(this%cbcauxitems, 'CBCAUXITEMS', this%memoryPath)
    call mem_allocate(this%idense, 'IDENSE', this%memoryPath)
    !
    ! -- Set values
    this%iprhed = 0
    this%istageout = 0
    this%ibudgetout = 0
    this%ibudcsv = 0
    this%ipakcsv = 0
    this%nlakes = 0
    this%noutlets = 0
    this%ntables = 0
    this%convlength = DONE
    this%convtime = DONE
    this%outdmax = DZERO
    this%igwhcopt = 0
    this%iconvchk = 1
    this%iconvresidchk = 1
    this%maxlakit = MAXADPIT
    this%surfdep = DZERO
    this%dmaxchg = DEM5
    this%delh = DP999 * this%dmaxchg
    this%pdmax = DEM1
    this%bditems = 11
    this%cbcauxitems = 1
    this%idense = 0
    this%ivsc = 0
  end subroutine lak_allocate_scalars

  !> @brief Allocate scalar members
  !<
  subroutine lak_allocate_arrays(this)
    ! -- modules
    ! -- dummy
    class(LakType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i
    !
    ! -- call standard BndType allocate scalars
    call this%BndType%allocate_arrays()
    !
    ! -- allocate character array for budget text
    allocate (this%clakbudget(this%bditems))
    !
    !-- fill clakbudget
    this%clakbudget(1) = '             GWF'
    this%clakbudget(2) = '        RAINFALL'
    this%clakbudget(3) = '     EVAPORATION'
    this%clakbudget(4) = '          RUNOFF'
    this%clakbudget(5) = '      EXT-INFLOW'
    this%clakbudget(6) = '      WITHDRAWAL'
    this%clakbudget(7) = '     EXT-OUTFLOW'
    this%clakbudget(8) = '         STORAGE'
    this%clakbudget(9) = '        CONSTANT'
    this%clakbudget(10) = '        FROM-MVR'
    this%clakbudget(11) = '          TO-MVR'
    !
    ! -- allocate and initialize dbuff
    if (this%istageout > 0) then
      call mem_allocate(this%dbuff, this%nlakes, 'DBUFF', this%memoryPath)
      do i = 1, this%nlakes
        this%dbuff(i) = DZERO
      end do
    else
      call mem_allocate(this%dbuff, 0, 'DBUFF', this%memoryPath)
    end if
    !
    ! -- allocate character array for budget text
    allocate (this%cauxcbc(this%cbcauxitems))
    !
    ! -- allocate and initialize qauxcbc
    call mem_allocate(this%qauxcbc, this%cbcauxitems, 'QAUXCBC', this%memoryPath)
    do i = 1, this%cbcauxitems
      this%qauxcbc(i) = DZERO
    end do
    !
    ! -- allocate qleak and qsto
    call mem_allocate(this%qleak, this%maxbound, 'QLEAK', this%memoryPath)
    do i = 1, this%maxbound
      this%qleak(i) = DZERO
    end do
    call mem_allocate(this%qsto, this%nlakes, 'QSTO', this%memoryPath)
    do i = 1, this%nlakes
      this%qsto(i) = DZERO
    end do
    !
    ! -- allocate denseterms to size 0
    call mem_allocate(this%denseterms, 3, 0, 'DENSETERMS', this%memoryPath)
    !
    ! -- allocate viscratios to size 0
    call mem_allocate(this%viscratios, 2, 0, 'VISCRATIOS', this%memoryPath)
  end subroutine lak_allocate_arrays

  !> @brief Read the dimensions for this package
  !<
  subroutine lak_read_lakes(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error, count_errors, store_error_unit
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(LakType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: text
    character(len=LENBOUNDNAME) :: bndName, bndNameTemp
    character(len=9) :: cno
    character(len=50), dimension(:), allocatable :: caux
    integer(I4B) :: ierr, ival
    logical(LGP) :: isfound, endOfBlock
    integer(I4B) :: n
    integer(I4B) :: ii, jj
    integer(I4B) :: iaux
    integer(I4B) :: itmp
    integer(I4B) :: nlak
    integer(I4B) :: nconn
    integer(I4B), dimension(:), pointer, contiguous :: nboundchk
    real(DP), pointer :: bndElem => null()
    !
    ! -- initialize itmp
    itmp = 0
    !
    ! -- allocate lake data
    call mem_allocate(this%nlakeconn, this%nlakes, 'NLAKECONN', this%memoryPath)
    call mem_allocate(this%idxlakeconn, this%nlakes + 1, 'IDXLAKECONN', &
                      this%memoryPath)
    call mem_allocate(this%ntabrow, this%nlakes, 'NTABROW', this%memoryPath)
    call mem_allocate(this%strt, this%nlakes, 'STRT', this%memoryPath)
    call mem_allocate(this%laketop, this%nlakes, 'LAKETOP', this%memoryPath)
    call mem_allocate(this%lakebot, this%nlakes, 'LAKEBOT', this%memoryPath)
    call mem_allocate(this%sareamax, this%nlakes, 'SAREAMAX', this%memoryPath)
    call mem_allocate(this%stage, this%nlakes, 'STAGE', this%memoryPath)
    call mem_allocate(this%rainfall, this%nlakes, 'RAINFALL', this%memoryPath)
    call mem_allocate(this%evaporation, this%nlakes, 'EVAPORATION', &
                      this%memoryPath)
    call mem_allocate(this%runoff, this%nlakes, 'RUNOFF', this%memoryPath)
    call mem_allocate(this%inflow, this%nlakes, 'INFLOW', this%memoryPath)
    call mem_allocate(this%withdrawal, this%nlakes, 'WITHDRAWAL', this%memoryPath)
    call mem_allocate(this%lauxvar, this%naux, this%nlakes, 'LAUXVAR', &
                      this%memoryPath)
    call mem_allocate(this%avail, this%nlakes, 'AVAIL', this%memoryPath)
    call mem_allocate(this%lkgwsink, this%nlakes, 'LKGWSINK', this%memoryPath)
    call mem_allocate(this%ncncvr, this%nlakes, 'NCNCVR', this%memoryPath)
    call mem_allocate(this%surfin, this%nlakes, 'SURFIN', this%memoryPath)
    call mem_allocate(this%surfout, this%nlakes, 'SURFOUT', this%memoryPath)
    call mem_allocate(this%surfout1, this%nlakes, 'SURFOUT1', this%memoryPath)
    call mem_allocate(this%precip, this%nlakes, 'PRECIP', this%memoryPath)
    call mem_allocate(this%precip1, this%nlakes, 'PRECIP1', this%memoryPath)
    call mem_allocate(this%evap, this%nlakes, 'EVAP', this%memoryPath)
    call mem_allocate(this%evap1, this%nlakes, 'EVAP1', this%memoryPath)
    call mem_allocate(this%evapo, this%nlakes, 'EVAPO', this%memoryPath)
    call mem_allocate(this%withr, this%nlakes, 'WITHR', this%memoryPath)
    call mem_allocate(this%withr1, this%nlakes, 'WITHR1', this%memoryPath)
    call mem_allocate(this%flwin, this%nlakes, 'FLWIN', this%memoryPath)
    call mem_allocate(this%flwiter, this%nlakes, 'FLWITER', this%memoryPath)
    call mem_allocate(this%flwiter1, this%nlakes, 'FLWITER1', this%memoryPath)
    call mem_allocate(this%seep, this%nlakes, 'SEEP', this%memoryPath)
    call mem_allocate(this%seep1, this%nlakes, 'SEEP1', this%memoryPath)
    call mem_allocate(this%seep0, this%nlakes, 'SEEP0', this%memoryPath)
    call mem_allocate(this%stageiter, this%nlakes, 'STAGEITER', this%memoryPath)
    call mem_allocate(this%chterm, this%nlakes, 'CHTERM', this%memoryPath)
    !
    ! -- lake boundary and stages
    call mem_allocate(this%iboundpak, this%nlakes, 'IBOUND', this%memoryPath)
    call mem_allocate(this%xnewpak, this%nlakes, 'XNEWPAK', this%memoryPath)
    call mem_allocate(this%xoldpak, this%nlakes, 'XOLDPAK', this%memoryPath)
    !
    ! -- lake iteration variables
    call mem_allocate(this%iseepc, this%nlakes, 'ISEEPC', this%memoryPath)
    call mem_allocate(this%idhc, this%nlakes, 'IDHC', this%memoryPath)
    call mem_allocate(this%en1, this%nlakes, 'EN1', this%memoryPath)
    call mem_allocate(this%en2, this%nlakes, 'EN2', this%memoryPath)
    call mem_allocate(this%r1, this%nlakes, 'R1', this%memoryPath)
    call mem_allocate(this%r2, this%nlakes, 'R2', this%memoryPath)
    call mem_allocate(this%dh0, this%nlakes, 'DH0', this%memoryPath)
    call mem_allocate(this%s0, this%nlakes, 'S0', this%memoryPath)
    call mem_allocate(this%qgwf0, this%nlakes, 'QGWF0', this%memoryPath)
    !
    ! -- allocate character storage not managed by the memory manager
    allocate (this%lakename(this%nlakes)) ! ditch after boundnames allocated??
    allocate (this%status(this%nlakes))
    !
    do n = 1, this%nlakes
      this%ntabrow(n) = 0
      this%status(n) = 'ACTIVE'
      this%laketop(n) = -DEP20
      this%lakebot(n) = DEP20
      this%sareamax(n) = DZERO
      this%iboundpak(n) = 1
      this%xnewpak(n) = DEP20
      this%xoldpak(n) = DEP20
      !
      ! -- initialize boundary values to zero
      this%rainfall(n) = DZERO
      this%evaporation(n) = DZERO
      this%runoff(n) = DZERO
      this%inflow(n) = DZERO
      this%withdrawal(n) = DZERO
    end do
    !
    ! -- allocate local storage for aux variables
    if (this%naux > 0) then
      allocate (caux(this%naux))
    end if
    !
    ! -- allocate and initialize temporary variables
    allocate (nboundchk(this%nlakes))
    do n = 1, this%nlakes
      nboundchk(n) = 0
    end do
    !
    ! -- read lake well data
    ! -- get lakes block
    call this%parser%GetBlock('PACKAGEDATA', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse locations block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text))// &
        ' PACKAGEDATA'
      nlak = 0
      nconn = 0
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        n = this%parser%GetInteger()
        !
        if (n < 1 .or. n > this%nlakes) then
          write (errmsg, '(a,1x,i0)') 'lakeno MUST BE > 0 and <= ', this%nlakes
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- increment nboundchk
        nboundchk(n) = nboundchk(n) + 1
        !
        ! -- strt
        this%strt(n) = this%parser%GetDouble()
        !
        ! nlakeconn
        ival = this%parser%GetInteger()
        !
        if (ival < 0) then
          write (errmsg, '(a,1x,i0)') 'nlakeconn MUST BE >= 0 for lake ', n
          call store_error(errmsg)
        end if
        !
        nconn = nconn + ival
        this%nlakeconn(n) = ival
        !
        ! -- get aux data
        do iaux = 1, this%naux
          call this%parser%GetString(caux(iaux))
        end do
        !
        ! -- set default bndName
        write (cno, '(i9.9)') n
        bndName = 'Lake'//cno
        !
        ! -- lakename
        if (this%inamedbound /= 0) then
          call this%parser%GetStringCaps(bndNameTemp)
          if (bndNameTemp /= '') then
            bndName = bndNameTemp
          end if
        end if
        this%lakename(n) = bndName
        !
        ! -- fill time series aware data
        ! -- fill aux data
        do jj = 1, this%naux
          text = caux(jj)
          ii = n
          bndElem => this%lauxvar(jj, ii)
          call read_value_or_time_series_adv(text, ii, jj, bndElem, &
                                             this%packName, 'AUX', &
                                             this%tsManager, this%iprpak, &
                                             this%auxname(jj))
        end do
        !
        nlak = nlak + 1
      end do
      !
      ! -- check for duplicate or missing lakes
      do n = 1, this%nlakes
        if (nboundchk(n) == 0) then
          write (errmsg, '(a,1x,i0)') 'NO DATA SPECIFIED FOR LAKE', n
          call store_error(errmsg)
        else if (nboundchk(n) > 1) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
            'DATA FOR LAKE', n, 'SPECIFIED', nboundchk(n), 'TIMES'
          call store_error(errmsg)
        end if
      end do
      !
      write (this%iout, '(1x,a)') 'END OF '//trim(adjustl(this%text))// &
        ' PACKAGEDATA'
    else
      call store_error('REQUIRED PACKAGEDATA BLOCK NOT FOUND.')
    end if
    !
    ! -- terminate if any errors were detected
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- set MAXBOUND
    this%MAXBOUND = nconn
    write (this%iout, '(//4x,a,i7)') 'MAXBOUND = ', this%maxbound
    !
    ! -- set idxlakeconn
    this%idxlakeconn(1) = 1
    do n = 1, this%nlakes
      this%idxlakeconn(n + 1) = this%idxlakeconn(n) + this%nlakeconn(n)
    end do
    !
    ! -- deallocate local storage for aux variables
    if (this%naux > 0) then
      deallocate (caux)
    end if
    !
    ! -- deallocate local storage for nboundchk
    deallocate (nboundchk)
  end subroutine lak_read_lakes

  !> @brief Read the lake connections for this package
  !<
  subroutine lak_read_lake_connections(this)
    use ConstantsModule, only: LINELENGTH, LENVARNAME
    use SimModule, only: store_error, count_errors
    ! -- dummy
    class(LakType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: keyword, cellid
    integer(I4B) :: ierr, ival
    logical(LGP) :: isfound, endOfBlock
    logical(LGP) :: is_lake_bed
    real(DP) :: rval
    integer(I4B) :: j, n
    integer(I4B) :: nn
    integer(I4B) :: ipos, ipos0
    integer(I4B) :: icellid, icellid0
    real(DP) :: top
    real(DP) :: bot
    integer(I4B), dimension(:), pointer, contiguous :: nboundchk
    character(len=LENVARNAME) :: ctypenm
    !
    ! -- allocate local storage
    allocate (nboundchk(this%MAXBOUND))
    do n = 1, this%MAXBOUND
      nboundchk(n) = 0
    end do
    !
    ! -- get connectiondata block
    call this%parser%GetBlock('CONNECTIONDATA', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse connectiondata block if detected
    if (isfound) then
      ! -- allocate connection data using memory manager
      call mem_allocate(this%imap, this%MAXBOUND, 'IMAP', this%memoryPath)
      call mem_allocate(this%cellid, this%MAXBOUND, 'CELLID', this%memoryPath)
      call mem_allocate(this%nodesontop, this%MAXBOUND, 'NODESONTOP', &
                        this%memoryPath)
      call mem_allocate(this%ictype, this%MAXBOUND, 'ICTYPE', this%memoryPath)
      call mem_allocate(this%bedleak, this%MAXBOUND, 'BEDLEAK', this%memoryPath) ! don't need to save this - use a temporary vector
      call mem_allocate(this%belev, this%MAXBOUND, 'BELEV', this%memoryPath)
      call mem_allocate(this%telev, this%MAXBOUND, 'TELEV', this%memoryPath)
      call mem_allocate(this%connlength, this%MAXBOUND, 'CONNLENGTH', &
                        this%memoryPath)
      call mem_allocate(this%connwidth, this%MAXBOUND, 'CONNWIDTH', &
                        this%memoryPath)
      call mem_allocate(this%sarea, this%MAXBOUND, 'SAREA', this%memoryPath)
      call mem_allocate(this%warea, this%MAXBOUND, 'WAREA', this%memoryPath)
      call mem_allocate(this%satcond, this%MAXBOUND, 'SATCOND', this%memoryPath)
      call mem_allocate(this%simcond, this%MAXBOUND, 'SIMCOND', this%memoryPath)
      call mem_allocate(this%simlakgw, this%MAXBOUND, 'SIMLAKGW', this%memoryPath)
      !
      ! -- process the lake connection data
      write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text))// &
        ' LAKE_CONNECTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        n = this%parser%GetInteger()
        !
        if (n < 1 .or. n > this%nlakes) then
          write (errmsg, '(a,1x,i0)') 'lakeno MUST BE > 0 and <= ', this%nlakes
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- read connection number
        ival = this%parser%GetInteger()
        if (ival < 1 .or. ival > this%nlakeconn(n)) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0)') &
            'iconn FOR LAKE ', n, 'MUST BE > 1 and <= ', this%nlakeconn(n)
          call store_error(errmsg)
          cycle
        end if
        !
        j = ival
        ipos = this%idxlakeconn(n) + ival - 1
        !
        ! -- set imap
        this%imap(ipos) = n
        !
        !
        ! -- increment nboundchk
        nboundchk(ipos) = nboundchk(ipos) + 1
        !
        ! -- read gwfnodes from the line
        call this%parser%GetCellid(this%dis%ndim, cellid)
        nn = this%dis%noder_from_cellid(cellid, &
                                        this%parser%iuactive, this%iout)
        !
        ! -- determine if a valid cell location was provided
        if (nn < 1) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0)') &
            'INVALID cellid FOR LAKE ', n, 'connection', j
          call store_error(errmsg)
        end if
        !
        ! -- set gwf cellid for connection
        this%cellid(ipos) = nn
        this%nodesontop(ipos) = nn
        !
        ! -- read ictype
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('VERTICAL')
          this%ictype(ipos) = 0
        case ('HORIZONTAL')
          this%ictype(ipos) = 1
        case ('EMBEDDEDH')
          this%ictype(ipos) = 2
        case ('EMBEDDEDV')
          this%ictype(ipos) = 3
        case default
          write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a,a,a)') &
            'UNKNOWN ctype FOR LAKE ', n, 'connection', j, &
            '(', trim(keyword), ')'
          call store_error(errmsg)
        end select
        write (ctypenm, '(a16)') keyword
        !
        ! -- bed leakance
        !this%bedleak(ipos) = this%parser%GetDouble() !TODO: use this when NONE keyword deprecated
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('NONE')
          is_lake_bed = .FALSE.
          this%bedleak(ipos) = DNODATA
          !
          ! -- create warning message
          write (warnmsg, '(2(a,1x,i0,1x),a,1pe8.1,a)') &
            'BEDLEAK for connection', j, 'in lake', n, 'is specified to '// &
            'be NONE. Lake connections where the lake-GWF connection '// &
            'conductance is solely a function of aquifer properties '// &
            'in the connected GWF cell should be specified with a '// &
            'DNODATA (', DNODATA, ') value.'
          !
          ! -- create deprecation warning
          call deprecation_warning('CONNECTIONDATA', 'bedleak=NONE', '6.4.3', &
                                   warnmsg, this%parser%GetUnit())
        case default
          read (keyword, *) rval
          if (is_close(rval, DNODATA)) then
            is_lake_bed = .FALSE.
          else
            is_lake_bed = .TRUE.
          end if
          this%bedleak(ipos) = rval
        end select
        !
        if (is_lake_bed .and. this%bedleak(ipos) < DZERO) then
          write (errmsg, '(a,1x,i0,1x,a)') 'bedleak FOR LAKE ', n, 'MUST BE >= 0'
          call store_error(errmsg)
        end if
        !
        ! -- belev
        this%belev(ipos) = this%parser%GetDouble()
        !
        ! -- telev
        this%telev(ipos) = this%parser%GetDouble()
        !
        ! -- connection length
        rval = this%parser%GetDouble()
        if (rval <= DZERO) then
          if (this%ictype(ipos) == 1 .or. this%ictype(ipos) == 2 .or. &
              this%ictype(ipos) == 3) then
            write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a,a,1x,a)') &
              'connection length (connlen) FOR LAKE ', n, &
              ', CONNECTION NO.', j, ', MUST BE > 0 FOR SPECIFIED ', &
              'connection type (ctype)', ctypenm
            call store_error(errmsg)
          else
            rval = DZERO
          end if
        end if
        this%connlength(ipos) = rval
        !
        ! -- connection width
        rval = this%parser%GetDouble()
        if (rval < dzero) then
          if (this%ictype(ipos) == 1) then
            write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
              'cell width (connwidth) FOR LAKE ', n, &
              ' HORIZONTAL CONNECTION ', j, 'MUST BE >= 0'
            call store_error(errmsg)
          else
            rval = DZERO
          end if
        end if
        this%connwidth(ipos) = rval
      end do
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' CONNECTIONDATA'
    else
      call store_error('REQUIRED CONNECTIONDATA BLOCK NOT FOUND.')
    end if
    !
    ! -- terminate if any errors were detected
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- check that embedded lakes have only one connection
    do n = 1, this%nlakes
      j = 0
      do ipos = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        if (this%ictype(ipos) /= 2 .and. this%ictype(ipos) /= 3) cycle
        j = j + 1
        if (j > 1) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
            'nlakeconn FOR LAKE', n, 'EMBEDDED CONNECTION', j, ' EXCEEDS 1.'
          call store_error(errmsg)
        end if
      end do
    end do
    ! -- check that an embedded lake is not in the same cell as a lake
    !   with a vertical connection
    do n = 1, this%nlakes
      ipos0 = this%idxlakeconn(n)
      icellid0 = this%cellid(ipos0)
      if (this%ictype(ipos0) /= 2 .and. this%ictype(ipos0) /= 3) cycle
      do nn = 1, this%nlakes
        if (nn == n) cycle
        j = 0
        do ipos = this%idxlakeconn(nn), this%idxlakeconn(nn + 1) - 1
          j = j + 1
          icellid = this%cellid(ipos)
          if (icellid == icellid0) then
            if (this%ictype(ipos) == 0) then
              write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a,1x,i0,1x,a)') &
                'EMBEDDED LAKE', n, &
                'CANNOT COINCIDE WITH VERTICAL CONNECTION', j, &
                'IN LAKE', nn, '.'
              call store_error(errmsg)
            end if
          end if
        end do
      end do
    end do
    !
    ! -- process the data
    do n = 1, this%nlakes
      j = 0
      do ipos = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        j = j + 1
        nn = this%cellid(ipos)
        top = this%dis%top(nn)
        bot = this%dis%bot(nn)
        ! vertical connection
        if (this%ictype(ipos) == 0) then
          this%telev(ipos) = top + this%surfdep
          this%belev(ipos) = top
          this%lakebot(n) = min(this%belev(ipos), this%lakebot(n))
          ! horizontal connection
        else if (this%ictype(ipos) == 1) then
          if (this%belev(ipos) == this%telev(ipos)) then
            this%telev(ipos) = top
            this%belev(ipos) = bot
          else
            if (this%belev(ipos) >= this%telev(ipos)) then
              write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
                'telev FOR LAKE ', n, ' HORIZONTAL CONNECTION ', j, &
                'MUST BE >= belev'
              call store_error(errmsg)
            else if (this%belev(ipos) < bot) then
              write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a,1x,g15.7,1x,a)') &
                'belev FOR LAKE ', n, ' HORIZONTAL CONNECTION ', j, &
                'MUST BE >= cell bottom (', bot, ')'
              call store_error(errmsg)
            else if (this%telev(ipos) > top) then
              write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a,1x,g15.7,1x,a)') &
                'telev FOR LAKE ', n, ' HORIZONTAL CONNECTION ', j, &
                'MUST BE <= cell top (', top, ')'
              call store_error(errmsg)
            end if
          end if
          this%laketop(n) = max(this%telev(ipos), this%laketop(n))
          this%lakebot(n) = min(this%belev(ipos), this%lakebot(n))
          ! embedded connections
        else if (this%ictype(ipos) == 2 .or. this%ictype(ipos) == 3) then
          this%telev(ipos) = top
          this%belev(ipos) = bot
          this%lakebot(n) = bot
        end if
        !
        ! -- check for missing or duplicate lake connections
        if (nboundchk(ipos) == 0) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0)') &
            'NO DATA SPECIFIED FOR LAKE', n, 'CONNECTION', j
          call store_error(errmsg)
        else if (nboundchk(ipos) > 1) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a,1x,i0,1x,a)') &
            'DATA FOR LAKE', n, 'CONNECTION', j, &
            'SPECIFIED', nboundchk(ipos), 'TIMES'
          call store_error(errmsg)
        end if
        !
        ! -- set laketop if it has not been assigned
      end do
      if (this%laketop(n) == -DEP20) then
        this%laketop(n) = this%lakebot(n) + 100.
      end if
    end do
    !
    ! -- deallocate local variable
    deallocate (nboundchk)
    !
    ! -- write summary of lake_connection error messages
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine lak_read_lake_connections

  !> @brief Read the lake tables for this package
  !<
  subroutine lak_read_tables(this)
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error, count_errors
    ! -- dummy
    class(LakType), intent(inout) :: this
    ! -- local
    type(LakTabType), dimension(:), allocatable :: laketables
    character(len=LINELENGTH) :: line
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    logical(LGP) :: isfound, endOfBlock
    integer(I4B) :: n
    integer(I4B) :: iconn
    integer(I4B) :: ntabs
    integer(I4B), dimension(:), pointer, contiguous :: nboundchk
    !
    ! -- skip of no outlets
    if (this%ntables < 1) return
    !
    ! -- allocate and initialize nboundchk
    allocate (nboundchk(this%nlakes))
    do n = 1, this%nlakes
      nboundchk(n) = 0
    end do
    !
    ! -- allocate derived type for table data
    allocate (laketables(this%nlakes))
    !
    ! -- get lake_tables block
    call this%parser%GetBlock('TABLES', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse lake_tables block if detected
    if (isfound) then
      ntabs = 0
      ! -- process the lake table data
      write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text))// &
        ' LAKE_TABLES'
      readtable: do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        n = this%parser%GetInteger()
        !
        if (n < 1 .or. n > this%nlakes) then
          write (errmsg, '(a,1x,i0)') 'lakeno MUST BE > 0 and <= ', this%nlakes
          call store_error(errmsg)
          cycle readtable
        end if
        !
        ! -- increment ntab and nboundchk
        ntabs = ntabs + 1
        nboundchk(n) = nboundchk(n) + 1
        !
        ! -- read FILE keyword
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('TAB6')
          call this%parser%GetStringCaps(keyword)
          if (trim(adjustl(keyword)) /= 'FILEIN') then
            errmsg = 'TAB6 keyword must be followed by "FILEIN" '// &
                     'then by filename.'
            call store_error(errmsg)
            cycle readtable
          end if
          call this%parser%GetString(line)
          call this%lak_read_table(n, line, laketables(n))
        case default
          write (errmsg, '(a,1x,i0,1x,a)') &
            'LAKE TABLE ENTRY for LAKE ', n, 'MUST INCLUDE TAB6 KEYWORD'
          call store_error(errmsg)
          cycle readtable
        end select
      end do readtable
      !
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' LAKE_TABLES'
      !
      ! -- check for missing or duplicate lake connections
      if (ntabs < this%ntables) then
        write (errmsg, '(a,1x,i0,1x,a,1x,i0)') &
          'TABLE DATA ARE SPECIFIED', ntabs, &
          'TIMES BUT NTABLES IS SET TO', this%ntables
        call store_error(errmsg)
      end if
      do n = 1, this%nlakes
        if (this%ntabrow(n) > 0 .and. nboundchk(n) > 1) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
            'TABLE DATA FOR LAKE', n, 'SPECIFIED', nboundchk(n), 'TIMES'
          call store_error(errmsg)
        end if
      end do
    else
      call store_error('REQUIRED TABLES BLOCK NOT FOUND.')
    end if
    !
    ! -- deallocate local storage
    deallocate (nboundchk)
    !
    ! -- write summary of lake_table error messages
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- convert laketables to vectors
    call this%laktables_to_vectors(laketables)
    !
    ! -- destroy laketables
    do n = 1, this%nlakes
      if (this%ntabrow(n) > 0) then
        deallocate (laketables(n)%tabstage)
        deallocate (laketables(n)%tabvolume)
        deallocate (laketables(n)%tabsarea)
        iconn = this%idxlakeconn(n)
        if (this%ictype(iconn) == 2 .or. this%ictype(iconn) == 3) then
          deallocate (laketables(n)%tabwarea)
        end if
      end if
    end do
    deallocate (laketables)
  end subroutine lak_read_tables

  !> @brief Copy the laketables structure data into flattened vectors that are
  !! stored in the memory manager
  !<
  subroutine laktables_to_vectors(this, laketables)
    class(LakType), intent(inout) :: this
    type(LakTabType), intent(in), dimension(:), contiguous :: laketables
    integer(I4B) :: n
    integer(I4B) :: ntabrows
    integer(I4B) :: j
    integer(I4B) :: ipos
    integer(I4B) :: iconn
    !
    ! -- allocate index array for lak tables
    call mem_allocate(this%ialaktab, this%nlakes + 1, 'IALAKTAB', this%memoryPath)
    !
    ! -- Move the laktables structure information into flattened arrays
    this%ialaktab(1) = 1
    do n = 1, this%nlakes
      ! -- ialaktab contains a pointer into the flattened lak table data
      this%ialaktab(n + 1) = this%ialaktab(n) + this%ntabrow(n)
    end do
    !
    ! -- Allocate vectors for storing all lake table data
    ntabrows = this%ialaktab(this%nlakes + 1) - 1
    call mem_allocate(this%tabstage, ntabrows, 'TABSTAGE', this%memoryPath)
    call mem_allocate(this%tabvolume, ntabrows, 'TABVOLUME', this%memoryPath)
    call mem_allocate(this%tabsarea, ntabrows, 'TABSAREA', this%memoryPath)
    call mem_allocate(this%tabwarea, ntabrows, 'TABWAREA', this%memoryPath)
    !
    ! -- Copy data from laketables into vectors
    do n = 1, this%nlakes
      j = 1
      do ipos = this%ialaktab(n), this%ialaktab(n + 1) - 1
        this%tabstage(ipos) = laketables(n)%tabstage(j)
        this%tabvolume(ipos) = laketables(n)%tabvolume(j)
        this%tabsarea(ipos) = laketables(n)%tabsarea(j)
        iconn = this%idxlakeconn(n)
        if (this%ictype(iconn) == 2 .or. this%ictype(iconn) == 3) then
          !
          ! -- tabwarea only filled for ictype 2 and 3
          this%tabwarea(ipos) = laketables(n)%tabwarea(j)
        else
          this%tabwarea(ipos) = DZERO
        end if
        j = j + 1
      end do
    end do
  end subroutine laktables_to_vectors

  !> @brief Read the lake table for this package
  !<
  subroutine lak_read_table(this, ilak, filename, laketable)
    use ConstantsModule, only: LINELENGTH
    use InputOutputModule, only: openfile
    use SimModule, only: store_error, count_errors
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    character(len=*), intent(in) :: filename
    type(LakTabType), intent(inout) :: laketable
    ! -- local
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    logical(LGP) :: isfound, endOfBlock
    integer(I4B) :: iu
    integer(I4B) :: n
    integer(I4B) :: ipos
    integer(I4B) :: j
    integer(I4B) :: jmin
    integer(I4B) :: iconn
    real(DP) :: vol
    real(DP) :: sa
    real(DP) :: wa
    real(DP) :: v
    real(DP) :: v0
    type(BlockParserType) :: parser
    ! -- formats
    character(len=*), parameter :: fmttaberr = &
      &'(a,1x,i0,1x,a,1x,g15.6,1x,a,1x,i0,1x,a,1x,i0,1x,a,1x,g15.6,1x,a)'
    !
    ! -- initialize locals
    n = 0
    j = 0
    !
    ! -- open the table file
    iu = 0
    call openfile(iu, this%iout, filename, 'LAKE TABLE')
    call parser%Initialize(iu, this%iout)
    !
    ! -- get dimensions block
    call parser%GetBlock('DIMENSIONS', isfound, ierr, supportOpenClose=.true.)
    !
    ! -- parse lak table dimensions block if detected
    if (isfound) then
      ! -- process the lake table dimension data
      if (this%iprpak /= 0) then
        write (this%iout, '(/1x,a)') &
          'PROCESSING '//trim(adjustl(this%text))//' DIMENSIONS'
      end if
      readdims: do
        call parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call parser%GetStringCaps(keyword)
        select case (keyword)
        case ('NROW')
          n = parser%GetInteger()

          if (n < 1) then
            write (errmsg, '(a)') 'LAKE TABLE NROW MUST BE > 0'
            call store_error(errmsg)
          end if
        case ('NCOL')
          j = parser%GetInteger()

          if (this%ictype(ilak) == 2 .or. this%ictype(ilak) == 3) then
            jmin = 4
          else
            jmin = 3
          end if
          if (j < jmin) then
            write (errmsg, '(a,1x,i0)') 'LAKE TABLE NCOL MUST BE >= ', jmin
            call store_error(errmsg)
          end if
          !
        case default
          write (errmsg, '(a,a)') &
            'UNKNOWN '//trim(this%text)//' DIMENSIONS KEYWORD: ', trim(keyword)
          call store_error(errmsg)
        end select
      end do readdims
      if (this%iprpak /= 0) then
        write (this%iout, '(1x,a)') &
          'END OF '//trim(adjustl(this%text))//' DIMENSIONS'
      end if
    else
      call store_error('REQUIRED DIMENSIONS BLOCK NOT FOUND.')
    end if
    !
    ! -- check that ncol and nrow have been specified
    if (n < 1) then
      write (errmsg, '(a)') &
        'NROW NOT SPECIFIED IN THE LAKE TABLE DIMENSIONS BLOCK'
      call store_error(errmsg)
    end if
    if (j < 1) then
      write (errmsg, '(a)') &
        'NCOL NOT SPECIFIED IN THE LAKE TABLE DIMENSIONS BLOCK'
      call store_error(errmsg)
    end if
    !
    ! -- only read the lake table data if n and j are specified to be greater
    !    than zero
    if (n * j > 0) then
      !
      ! -- allocate space
      this%ntabrow(ilak) = n
      allocate (laketable%tabstage(n))
      allocate (laketable%tabvolume(n))
      allocate (laketable%tabsarea(n))
      ipos = this%idxlakeconn(ilak)
      if (this%ictype(ipos) == 2 .or. this%ictype(ipos) == 3) then
        allocate (laketable%tabwarea(n))
      end if
      !
      ! -- get table block
      call parser%GetBlock('TABLE', isfound, ierr, supportOpenClose=.true.)
      !
      ! -- parse well_connections block if detected
      if (isfound) then
        !
        ! -- process the table data
        if (this%iprpak /= 0) then
          write (this%iout, '(/1x,a)') &
            'PROCESSING '//trim(adjustl(this%text))//' TABLE'
        end if
        iconn = this%idxlakeconn(ilak)
        ipos = 0
        readtabledata: do
          call parser%GetNextLine(endOfBlock)
          if (endOfBlock) exit
          ipos = ipos + 1
          if (ipos > this%ntabrow(ilak)) then
            cycle readtabledata
          end if
          laketable%tabstage(ipos) = parser%GetDouble()
          laketable%tabvolume(ipos) = parser%GetDouble()
          laketable%tabsarea(ipos) = parser%GetDouble()
          if (this%ictype(iconn) == 2 .or. this%ictype(iconn) == 3) then
            laketable%tabwarea(ipos) = parser%GetDouble()
          end if
        end do readtabledata
        !
        if (this%iprpak /= 0) then
          write (this%iout, '(1x,a)') &
            'END OF '//trim(adjustl(this%text))//' TABLE'
        end if
      else
        call store_error('REQUIRED TABLE BLOCK NOT FOUND.')
      end if
      !
      ! -- error condition if number of rows read are not equal to nrow
      if (ipos /= this%ntabrow(ilak)) then
        write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
          'NROW SET TO', this%ntabrow(ilak), 'BUT', ipos, 'ROWS WERE READ'
        call store_error(errmsg)
      end if
      !
      ! -- set lake bottom based on table if it is an embedded lake
      iconn = this%idxlakeconn(ilak)
      if (this%ictype(iconn) == 2 .or. this%ictype(iconn) == 3) then
        do n = 1, this%ntabrow(ilak)
          vol = laketable%tabvolume(n)
          sa = laketable%tabsarea(n)
          wa = laketable%tabwarea(n)
          vol = vol * sa * wa
          ! -- check if all entries are zero
          if (vol > DZERO) exit
          ! -- set lake bottom
          this%lakebot(ilak) = laketable%tabstage(n)
          this%belev(ilak) = laketable%tabstage(n)
        end do
        ! -- set maximum surface area for rainfall
        n = this%ntabrow(ilak)
        this%sareamax(ilak) = laketable%tabsarea(n)
      end if
      !
      ! -- verify the table data
      do n = 2, this%ntabrow(ilak)
        v = laketable%tabstage(n)
        v0 = laketable%tabstage(n - 1)
        if (v <= v0) then
          write (errmsg, fmttaberr) &
            'TABLE STAGE ENTRY', n, '(', laketable%tabstage(n), ') FOR LAKE ', &
            ilak, 'MUST BE GREATER THAN THE PREVIOUS STAGE ENTRY', &
            n - 1, '(', laketable%tabstage(n - 1), ')'
          call store_error(errmsg)
        end if
        v = laketable%tabvolume(n)
        v0 = laketable%tabvolume(n - 1)
        if (v <= v0) then
          write (errmsg, fmttaberr) &
            'TABLE VOLUME ENTRY', n, '(', laketable%tabvolume(n), &
            ') FOR LAKE ', &
            ilak, 'MUST BE GREATER THAN THE PREVIOUS VOLUME ENTRY', &
            n - 1, '(', laketable%tabvolume(n - 1), ')'
          call store_error(errmsg)
        end if
        v = laketable%tabsarea(n)
        v0 = laketable%tabsarea(n - 1)
        if (v < v0) then
          write (errmsg, fmttaberr) &
            'TABLE SURFACE AREA ENTRY', n, '(', &
            laketable%tabsarea(n), ') FOR LAKE ', ilak, &
            'MUST BE GREATER THAN OR EQUAL TO THE PREVIOUS SURFACE AREA ENTRY', &
            n - 1, '(', laketable%tabsarea(n - 1), ')'
          call store_error(errmsg)
        end if
        iconn = this%idxlakeconn(ilak)
        if (this%ictype(iconn) == 2 .or. this%ictype(iconn) == 3) then
          v = laketable%tabwarea(n)
          v0 = laketable%tabwarea(n - 1)
          if (v < v0) then
            write (errmsg, fmttaberr) &
              'TABLE EXCHANGE AREA ENTRY', n, '(', &
              laketable%tabwarea(n), ') FOR LAKE ', ilak, &
              'MUST BE GREATER THAN OR EQUAL TO THE PREVIOUS EXCHANGE AREA '// &
              'ENTRY', n - 1, '(', laketable%tabwarea(n - 1), ')'
            call store_error(errmsg)
          end if
        end if
      end do
    end if
    !
    ! -- write summary of lake table error messages
    if (count_errors() > 0) then
      call parser%StoreErrorUnit()
    end if
    !
    ! Close the table file and clear other parser members
    call parser%Clear()
  end subroutine lak_read_table

  !> @brief Read the lake outlets for this package
  !<
  subroutine lak_read_outlets(this)
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error, count_errors
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(LakType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: text, keyword
    character(len=LENBOUNDNAME) :: bndName
    character(len=9) :: citem
    integer(I4B) :: ierr, ival
    logical(LGP) :: isfound, endOfBlock
    integer(I4B) :: n
    integer(I4B) :: jj
    integer(I4B), dimension(:), pointer, contiguous :: nboundchk
    real(DP), pointer :: bndElem => null()
    !
    ! -- get well_connections block
    call this%parser%GetBlock('OUTLETS', isfound, ierr, &
                              supportOpenClose=.true., blockRequired=.false.)
    !
    ! -- parse outlets block if detected
    if (isfound) then
      if (this%noutlets > 0) then
        !
        ! -- allocate and initialize local variables
        allocate (nboundchk(this%noutlets))
        do n = 1, this%noutlets
          nboundchk(n) = 0
        end do
        !
        ! -- allocate outlet data using memory manager
        call mem_allocate(this%lakein, this%NOUTLETS, 'LAKEIN', this%memoryPath)
        call mem_allocate(this%lakeout, this%NOUTLETS, 'LAKEOUT', this%memoryPath)
        call mem_allocate(this%iouttype, this%NOUTLETS, 'IOUTTYPE', &
                          this%memoryPath)
        call mem_allocate(this%outrate, this%NOUTLETS, 'OUTRATE', this%memoryPath)
        call mem_allocate(this%outinvert, this%NOUTLETS, 'OUTINVERT', &
                          this%memoryPath)
        call mem_allocate(this%outwidth, this%NOUTLETS, 'OUTWIDTH', &
                          this%memoryPath)
        call mem_allocate(this%outrough, this%NOUTLETS, 'OUTROUGH', &
                          this%memoryPath)
        call mem_allocate(this%outslope, this%NOUTLETS, 'OUTSLOPE', &
                          this%memoryPath)
        call mem_allocate(this%simoutrate, this%NOUTLETS, 'SIMOUTRATE', &
                          this%memoryPath)
        !
        ! -- initialize outlet rate
        do n = 1, this%noutlets
          this%outrate(n) = DZERO
        end do
        !
        ! -- process the lake connection data
        write (this%iout, '(/1x,a)') &
          'PROCESSING '//trim(adjustl(this%text))//' OUTLETS'
        readoutlet: do
          call this%parser%GetNextLine(endOfBlock)
          if (endOfBlock) exit
          n = this%parser%GetInteger()

          if (n < 1 .or. n > this%noutlets) then
            write (errmsg, '(a,1x,i0)') &
              'outletno MUST BE > 0 and <= ', this%noutlets
            call store_error(errmsg)
            cycle readoutlet
          end if
          !
          ! -- increment nboundchk
          nboundchk(n) = nboundchk(n) + 1
          !
          ! -- read outlet lakein
          ival = this%parser%GetInteger()
          if (ival < 1 .or. ival > this%nlakes) then
            write (errmsg, '(a,1x,i0,1x,a,1x,i0)') &
              'lakein FOR OUTLET ', n, 'MUST BE > 0 and <= ', this%nlakes
            call store_error(errmsg)
            cycle readoutlet
          end if
          this%lakein(n) = ival
          !
          ! -- read outlet lakeout
          ival = this%parser%GetInteger()
          if (ival < 0 .or. ival > this%nlakes) then
            write (errmsg, '(a,1x,i0,1x,a,1x,i0)') &
              'lakeout FOR OUTLET ', n, 'MUST BE >= 0 and <= ', this%nlakes
            call store_error(errmsg)
            cycle readoutlet
          end if
          this%lakeout(n) = ival
          !
          ! -- read ictype
          call this%parser%GetStringCaps(keyword)
          select case (keyword)
          case ('SPECIFIED')
            this%iouttype(n) = 0
          case ('MANNING')
            this%iouttype(n) = 1
          case ('WEIR')
            this%iouttype(n) = 2
          case default
            write (errmsg, '(a,1x,i0,1x,a,a,a)') &
              'UNKNOWN couttype FOR OUTLET ', n, '(', trim(keyword), ')'
            call store_error(errmsg)
            cycle readoutlet
          end select
          !
          ! -- build bndname for outlet
          write (citem, '(i9.9)') n
          bndName = 'OUTLET'//citem
          !
          ! -- set a few variables for timeseries aware variables
          jj = 1
          !
          ! -- outlet invert
          call this%parser%GetString(text)
          bndElem => this%outinvert(n)
          call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                             this%packName, 'BND', &
                                             this%tsManager, this%iprpak, &
                                             'INVERT')
          !
          ! -- outlet width
          call this%parser%GetString(text)
          bndElem => this%outwidth(n)
          call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                             this%packName, 'BND', &
                                             this%tsManager, this%iprpak, 'WIDTH')
          !
          ! -- outlet roughness
          call this%parser%GetString(text)
          bndElem => this%outrough(n)
          call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                             this%packName, 'BND', &
                                             this%tsManager, this%iprpak, 'ROUGH')
          !
          ! -- outlet slope
          call this%parser%GetString(text)
          bndElem => this%outslope(n)
          call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                             this%packName, 'BND', &
                                             this%tsManager, this%iprpak, 'SLOPE')
        end do readoutlet
        write (this%iout, '(1x,a)') 'END OF '//trim(adjustl(this%text))// &
          ' OUTLETS'
        !
        ! -- check for duplicate or missing outlets
        do n = 1, this%noutlets
          if (nboundchk(n) == 0) then
            write (errmsg, '(a,1x,i0)') 'NO DATA SPECIFIED FOR OUTLET', n
            call store_error(errmsg)
          else if (nboundchk(n) > 1) then
            write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
              'DATA FOR OUTLET', n, 'SPECIFIED', nboundchk(n), 'TIMES'
            call store_error(errmsg)
          end if
        end do
        !
        ! -- deallocate local storage
        deallocate (nboundchk)
      else
        write (errmsg, '(a,1x,a)') &
          'AN OUTLETS BLOCK SHOULD NOT BE SPECIFIED IF NOUTLETS IS NOT', &
          'SPECIFIED OR IS SPECIFIED TO BE 0.'
        call store_error(errmsg)
      end if
      !
    else
      if (this%noutlets > 0) then
        call store_error('REQUIRED OUTLETS BLOCK NOT FOUND.')
      end if
    end if
    !
    ! -- write summary of lake_connection error messages
    ierr = count_errors()
    if (ierr > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine lak_read_outlets

  !> @brief Read the dimensions for this package
  !<
  subroutine lak_read_dimensions(this)
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error, count_errors
    ! -- dummy
    class(LakType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    logical(LGP) :: isfound, endOfBlock
    !
    ! -- initialize dimensions to -1
    this%nlakes = -1
    this%maxbound = -1
    !
    ! -- get dimensions block
    call this%parser%GetBlock('DIMENSIONS', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse dimensions block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text))// &
        ' DIMENSIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('NLAKES')
          this%nlakes = this%parser%GetInteger()
          write (this%iout, '(4x,a,i7)') 'NLAKES = ', this%nlakes
        case ('NOUTLETS')
          this%noutlets = this%parser%GetInteger()
          write (this%iout, '(4x,a,i7)') 'NOUTLETS = ', this%noutlets
        case ('NTABLES')
          this%ntables = this%parser%GetInteger()
          write (this%iout, '(4x,a,i7)') 'NTABLES = ', this%ntables
        case default
          write (errmsg, '(a,a)') &
            'UNKNOWN '//trim(this%text)//' DIMENSION: ', trim(keyword)
          call store_error(errmsg)
        end select
      end do
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' DIMENSIONS'
    else
      call store_error('REQUIRED DIMENSIONS BLOCK NOT FOUND.')
    end if
    !
    if (this%nlakes < 0) then
      write (errmsg, '(a)') &
        'NLAKES WAS NOT SPECIFIED OR WAS SPECIFIED INCORRECTLY.'
      call store_error(errmsg)
    end if
    !
    ! -- stop if errors were encountered in the DIMENSIONS block
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- read lakes block
    call this%lak_read_lakes()
    !
    ! -- read lake_connections block
    call this%lak_read_lake_connections()
    !
    ! -- read tables block
    call this%lak_read_tables()
    !
    ! -- read outlets block
    call this%lak_read_outlets()
    !
    ! -- Call define_listlabel to construct the list label that is written
    !    when PRINT_INPUT option is used.
    call this%define_listlabel()
    !
    ! -- setup the budget object
    call this%lak_setup_budobj()
    !
    ! -- setup the stage table object
    call this%lak_setup_tableobj()
  end subroutine lak_read_dimensions

  !> @brief Read the initial parameters for this package
  !<
  subroutine lak_read_initial_attr(this)
    use ConstantsModule, only: LINELENGTH
    use MemoryHelperModule, only: create_mem_path
    use SimModule, only: store_error, count_errors
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(LakType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: text
    integer(I4B) :: j, jj, n
    integer(I4B) :: nn
    integer(I4B) :: idx
    real(DP) :: top
    real(DP) :: bot
    real(DP) :: k
    real(DP) :: area
    real(DP) :: length
    real(DP) :: s
    real(DP) :: dx
    real(DP) :: c
    real(DP) :: sa
    real(DP) :: wa
    real(DP) :: v
    real(DP) :: fact
    real(DP) :: c1
    real(DP) :: c2
    real(DP), allocatable, dimension(:) :: clb, caq
    character(len=14) :: cbedleak
    character(len=14) :: cbedcond
    character(len=10), dimension(0:3) :: ctype
    character(len=15) :: nodestr
    real(DP), pointer :: bndElem => null()
    ! -- data
    data ctype(0)/'VERTICAL  '/
    data ctype(1)/'HORIZONTAL'/
    data ctype(2)/'EMBEDDEDH '/
    data ctype(3)/'EMBEDDEDV '/
    !
    ! -- initialize xnewpak and set stage
    do n = 1, this%nlakes
      this%xnewpak(n) = this%strt(n)
      write (text, '(g15.7)') this%strt(n)
      jj = 1 ! For STAGE
      bndElem => this%stage(n)
      call read_value_or_time_series_adv(text, n, jj, bndElem, this%packName, &
                                         'BND', this%tsManager, this%iprpak, &
                                         'STAGE')
    end do
    !
    ! -- initialize status (iboundpak) of lakes to active
    do n = 1, this%nlakes
      if (this%status(n) == 'CONSTANT') then
        this%iboundpak(n) = -1
      else if (this%status(n) == 'INACTIVE') then
        this%iboundpak(n) = 0
      else if (this%status(n) == 'ACTIVE ') then
        this%iboundpak(n) = 1
      end if
    end do
    !
    ! -- set boundname for each connection
    if (this%inamedbound /= 0) then
      do n = 1, this%nlakes
        do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
          this%boundname(j) = this%lakename(n)
        end do
      end do
    end if
    !
    ! -- copy boundname into boundname_cst
    call this%copy_boundname()
    !
    ! -- set pointer to gwf iss and gwf hk
    call mem_setptr(this%gwfiss, 'ISS', create_mem_path(this%name_model))
    call mem_setptr(this%gwfk11, 'K11', create_mem_path(this%name_model, 'NPF'))
    call mem_setptr(this%gwfk33, 'K33', create_mem_path(this%name_model, 'NPF'))
    call mem_setptr(this%gwfik33, 'IK33', create_mem_path(this%name_model, 'NPF'))
    call mem_setptr(this%gwfsat, 'SAT', create_mem_path(this%name_model, 'NPF'))
    !
    ! -- allocate temporary storage
    allocate (clb(this%MAXBOUND))
    allocate (caq(this%MAXBOUND))
    !
    ! -- calculate saturated conductance for each connection
    do n = 1, this%nlakes
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        nn = this%cellid(j)
        top = this%dis%top(nn)
        bot = this%dis%bot(nn)
        ! vertical connection
        if (this%ictype(j) == 0) then
          area = this%dis%area(nn)
          this%sarea(j) = area
          this%warea(j) = area
          this%sareamax(n) = this%sareamax(n) + area
          if (this%gwfik33 == 0) then
            k = this%gwfk11(nn)
          else
            k = this%gwfk33(nn)
          end if
          length = DHALF * (top - bot)
          ! horizontal connection
        else if (this%ictype(j) == 1) then
          area = (this%telev(j) - this%belev(j)) * this%connwidth(j)
          ! -- recalculate area if connected cell is confined and lake
          !    connection top and bot are equal to the cell top and bot
          if (top == this%telev(j) .and. bot == this%belev(j)) then
            if (this%icelltype(nn) == 0) then
              area = this%gwfsat(nn) * (top - bot) * this%connwidth(j)
            end if
          end if
          this%sarea(j) = DZERO
          this%warea(j) = area
          this%sareamax(n) = this%sareamax(n) + DZERO
          k = this%gwfk11(nn)
          length = this%connlength(j)
          ! embedded horizontal connection
        else if (this%ictype(j) == 2) then
          area = DONE
          this%sarea(j) = DZERO
          this%warea(j) = area
          this%sareamax(n) = this%sareamax(n) + DZERO
          k = this%gwfk11(nn)
          length = this%connlength(j)
          ! embedded vertical connection
        else if (this%ictype(j) == 3) then
          area = DONE
          this%sarea(j) = DZERO
          this%warea(j) = area
          this%sareamax(n) = this%sareamax(n) + DZERO
          if (this%gwfik33 == 0) then
            k = this%gwfk11(nn)
          else
            k = this%gwfk33(nn)
          end if
          length = this%connlength(j)
        end if
        if (is_close(this%bedleak(j), DNODATA)) then
          clb(j) = DNODATA
        else if (this%bedleak(j) > DZERO) then
          clb(j) = DONE / this%bedleak(j)
        else
          clb(j) = DZERO
        end if
        if (k > DZERO) then
          caq(j) = length / k
        else
          caq(j) = DZERO
        end if
        if (is_close(this%bedleak(j), DNODATA)) then
          this%satcond(j) = area / caq(j)
        else if (clb(j) * caq(j) > DZERO) then
          this%satcond(j) = area / (clb(j) + caq(j))
        else
          this%satcond(j) = DZERO
        end if
      end do
    end do
    !
    ! -- write a summary of the conductance
    if (this%iprpak > 0) then
      write (this%iout, '(//,29x,a,/)') &
        'INTERFACE CONDUCTANCE BETWEEN LAKE AND AQUIFER CELLS'
      write (this%iout, '(1x,a)') &
     &  '      LAKE CONNECTION                 CONNECTION    LAKEBED'// &
     &  '              C O N D U C T A N C E S        '
      write (this%iout, '(1x,a)') &
     &  '    NUMBER     NUMBER CELLID          DIRECTION    LEAKANCE'// &
     &  '        LAKEBED        AQUIFER       COMBINED'
      write (this%iout, "(1x,108('-'))")
      do n = 1, this%nlakes
        idx = 0
        do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
          idx = idx + 1
          fact = DONE
          if (this%ictype(j) == 1) then
            fact = this%telev(j) - this%belev(j)
            if (ABS(fact) > DZERO) then
              fact = DONE / fact
            end if
          end if
          nn = this%cellid(j)
          area = this%warea(j)
          c1 = DZERO
          if (is_close(clb(j), DNODATA)) then
            cbedleak = '     NONE     '
            cbedcond = '     NONE     '
          else if (clb(j) > DZERO) then
            c1 = area * fact / clb(j)
            write (cbedleak, '(g14.5)') this%bedleak(j)
            write (cbedcond, '(g14.5)') c1
          else
            write (cbedleak, '(g14.5)') c1
            write (cbedcond, '(g14.5)') c1
          end if
          c2 = DZERO
          if (caq(j) > DZERO) then
            c2 = area * fact / caq(j)
          end if
          call this%dis%noder_to_string(nn, nodestr)
          write (this%iout, &
                 '(1x,i10,1x,i10,1x,a15,1x,a10,2(1x,a14),2(1x,g14.5))') &
            n, idx, nodestr, ctype(this%ictype(j)), cbedleak, &
            cbedcond, c2, this%satcond(j) * fact
        end do
      end do
      write (this%iout, "(1x,108('-'))")
      write (this%iout, '(1x,a)') &
        'IF VERTICAL CONNECTION, CONDUCTANCE (L^2/T) IS &
        &BETWEEN AQUIFER CELL AND OVERLYING LAKE CELL.'
      write (this%iout, '(1x,a)') &
        'IF HORIZONTAL CONNECTION, CONDUCTANCES ARE PER &
        &UNIT SATURATED THICKNESS (L/T).'
      write (this%iout, '(1x,a)') &
        'IF EMBEDDED CONNECTION, CONDUCTANCES ARE PER &
        &UNIT EXCHANGE AREA (1/T).'
      !
      !        write(this%iout,*) n, idx, nodestr, this%sarea(j), this%warea(j)
      !
      ! -- calculate stage, surface area, wetted area, volume relation
      do n = 1, this%nlakes
        write (this%iout, '(//1x,a,1x,i10)') 'STAGE/VOLUME RELATION FOR LAKE  ', n
        write (this%iout, '(/1x,5(a14))') '         STAGE', '  SURFACE AREA', &
    &                                    '   WETTED AREA', '   CONDUCTANCE', &
    &                                    '        VOLUME'
        write (this%iout, "(1x,70('-'))")
        dx = (this%laketop(n) - this%lakebot(n)) / 150.
        s = this%lakebot(n)
        do j = 1, 151
          call this%lak_calculate_conductance(n, s, c)
          call this%lak_calculate_sarea(n, s, sa)
          call this%lak_calculate_warea(n, s, wa, s)
          call this%lak_calculate_vol(n, s, v)
          write (this%iout, '(1x,5(E14.5))') s, sa, wa, c, v
          s = s + dx
        end do
        write (this%iout, "(1x,70('-'))")
        !
        write (this%iout, '(//1x,a,1x,i10)') 'STAGE/VOLUME RELATION FOR LAKE  ', n
        write (this%iout, '(/1x,4(a14))') '              ', '              ', &
    &                                    '    CALCULATED', '         STAGE'
        write (this%iout, '(1x,4(a14))') '         STAGE', '        VOLUME', &
    &                                    '         STAGE', '    DIFFERENCE'
        write (this%iout, "(1x,56('-'))")
        s = this%lakebot(n) - dx
        do j = 1, 156
          call this%lak_calculate_vol(n, s, v)
          call this%lak_vol2stage(n, v, c)
          write (this%iout, '(1x,4(E14.5))') s, v, c, s - c
          s = s + dx
        end do
        write (this%iout, "(1x,56('-'))")
      end do
    end if
    !
    ! -- finished with pointer to gwf hydraulic conductivity
    this%gwfk11 => null()
    this%gwfk33 => null()
    this%gwfsat => null()
    this%gwfik33 => null()
    !
    ! -- deallocate temporary storage
    deallocate (clb)
    deallocate (caq)
  end subroutine lak_read_initial_attr

  !> @brief Perform linear interpolation of two vectors.
  !!
  !! Function assumes x data is sorted in ascending order
  !<
  subroutine lak_linear_interpolation(this, n, x, y, z, v)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(n), intent(in) :: x
    real(DP), dimension(n), intent(in) :: y
    real(DP), intent(in) :: z
    real(DP), intent(inout) :: v
    ! -- local
    integer(I4B) :: i
    real(DP) :: dx, dydx
    ! code
    v = DZERO
    ! below bottom of range - set to lowest value
    if (z <= x(1)) then
      v = y(1)
      ! above highest value
      ! slope calculated from interval between n and n-1
    else if (z > x(n)) then
      dx = x(n) - x(n - 1)
      dydx = DZERO
      if (ABS(dx) > DZERO) then
        dydx = (y(n) - y(n - 1)) / dx
      end if
      dx = (z - x(n))
      v = y(n) + dydx * dx
      ! between lowest and highest value in current interval
    else
      do i = 2, n
        dx = x(i) - x(i - 1)
        dydx = DZERO
        if (z >= x(i - 1) .and. z <= x(i)) then
          if (ABS(dx) > DZERO) then
            dydx = (y(i) - y(i - 1)) / dx
          end if
          dx = (z - x(i - 1))
          v = y(i - 1) + dydx * dx
          exit
        end if
      end do
    end if
  end subroutine lak_linear_interpolation

  !> @brief Calculate the surface area of a lake at a given stage
  !<
  subroutine lak_calculate_sarea(this, ilak, stage, sarea)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(in) :: stage
    real(DP), intent(inout) :: sarea
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: ifirst
    integer(I4B) :: ilast
    real(DP) :: topl
    real(DP) :: botl
    real(DP) :: sat
    real(DP) :: sa
    !
    sarea = DZERO
    i = this%ntabrow(ilak)
    if (i > 0) then
      ifirst = this%ialaktab(ilak)
      ilast = this%ialaktab(ilak + 1) - 1
      if (stage <= this%tabstage(ifirst)) then
        sarea = this%tabsarea(ifirst)
      else if (stage >= this%tabstage(ilast)) then
        sarea = this%tabsarea(ilast)
      else
        call this%lak_linear_interpolation(i, this%tabstage(ifirst:ilast), &
                                           this%tabsarea(ifirst:ilast), &
                                           stage, sarea)
      end if
    else
      do i = this%idxlakeconn(ilak), this%idxlakeconn(ilak + 1) - 1
        topl = this%telev(i)
        botl = this%belev(i)
        sat = sQuadraticSaturation(topl, botl, stage)
        sa = sat * this%sarea(i)
        sarea = sarea + sa
      end do
    end if
  end subroutine lak_calculate_sarea

  !> @brief Calculate the wetted area of a lake at a given stage.
  !<
  subroutine lak_calculate_warea(this, ilak, stage, warea, hin)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(in) :: stage
    real(DP), intent(inout) :: warea
    real(DP), optional, intent(inout) :: hin
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: igwfnode
    real(DP) :: head
    real(DP) :: wa
    !
    warea = DZERO
    do i = this%idxlakeconn(ilak), this%idxlakeconn(ilak + 1) - 1
      if (present(hin)) then
        head = hin
      else
        igwfnode = this%cellid(i)
        head = this%xnew(igwfnode)
      end if
      call this%lak_calculate_conn_warea(ilak, i, stage, head, wa)
      warea = warea + wa
    end do
  end subroutine lak_calculate_warea

  !> @brief Calculate the wetted area of a lake connection at a given stage
  !<
  subroutine lak_calculate_conn_warea(this, ilak, iconn, stage, head, wa)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    integer(I4B), intent(in) :: iconn
    real(DP), intent(in) :: stage
    real(DP), intent(in) :: head
    real(DP), intent(inout) :: wa
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: ifirst
    integer(I4B) :: ilast
    integer(I4B) :: node
    real(DP) :: topl
    real(DP) :: botl
    real(DP) :: vv
    real(DP) :: sat
    !
    wa = DZERO
    topl = this%telev(iconn)
    botl = this%belev(iconn)
    call this%lak_calculate_cond_head(iconn, stage, head, vv)
    if (this%ictype(iconn) == 2 .or. this%ictype(iconn) == 3) then
      if (vv > topl) vv = topl
      i = this%ntabrow(ilak)
      ifirst = this%ialaktab(ilak)
      ilast = this%ialaktab(ilak + 1) - 1
      if (vv <= this%tabstage(ifirst)) then
        wa = this%tabwarea(ifirst)
      else if (vv >= this%tabstage(ilast)) then
        wa = this%tabwarea(ilast)
      else
        call this%lak_linear_interpolation(i, this%tabstage(ifirst:ilast), &
                                           this%tabwarea(ifirst:ilast), &
                                           vv, wa)
      end if
    else
      node = this%cellid(iconn)
      ! -- confined cell
      if (this%icelltype(node) == 0) then
        sat = DONE
        ! -- convertible cell
      else
        sat = sQuadraticSaturation(topl, botl, vv)
      end if
      wa = sat * this%warea(iconn)
    end if
  end subroutine lak_calculate_conn_warea

  !> @brief Calculate the volume of a lake at a given stage
  !<
  subroutine lak_calculate_vol(this, ilak, stage, volume)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(in) :: stage
    real(DP), intent(inout) :: volume
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: ifirst
    integer(I4B) :: ilast
    real(DP) :: topl
    real(DP) :: botl
    real(DP) :: ds
    real(DP) :: sa
    real(DP) :: v
    real(DP) :: sat
    !
    volume = DZERO
    i = this%ntabrow(ilak)
    if (i > 0) then
      ifirst = this%ialaktab(ilak)
      ilast = this%ialaktab(ilak + 1) - 1
      if (stage <= this%tabstage(ifirst)) then
        volume = this%tabvolume(ifirst)
      else if (stage >= this%tabstage(ilast)) then
        ds = stage - this%tabstage(ilast)
        sa = this%tabsarea(ilast)
        volume = this%tabvolume(ilast) + ds * sa
      else
        call this%lak_linear_interpolation(i, this%tabstage(ifirst:ilast), &
                                           this%tabvolume(ifirst:ilast), &
                                           stage, volume)
      end if
    else
      do i = this%idxlakeconn(ilak), this%idxlakeconn(ilak + 1) - 1
        topl = this%telev(i)
        botl = this%belev(i)
        sat = sQuadraticSaturation(topl, botl, stage)
        sa = sat * this%sarea(i)
        if (stage < botl) then
          v = DZERO
        else if (stage > botl .and. stage < topl) then
          v = sa * (stage - botl)
        else
          v = sa * (topl - botl) + sa * (stage - topl)
        end if
        volume = volume + v
      end do
    end if
  end subroutine lak_calculate_vol

  !> @brief Calculate the total conductance for a lake at a provided stage
  !<
  subroutine lak_calculate_conductance(this, ilak, stage, conductance)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(in) :: stage
    real(DP), intent(inout) :: conductance
    ! -- local
    integer(I4B) :: i
    real(DP) :: c
    !
    conductance = DZERO
    do i = this%idxlakeconn(ilak), this%idxlakeconn(ilak + 1) - 1
      call this%lak_calculate_conn_conductance(ilak, i, stage, stage, c)
      conductance = conductance + c
    end do
  end subroutine lak_calculate_conductance

  !> @brief Calculate the controlling lake stage or groundwater head used to
  !! calculate the conductance for a lake connection from a provided stage and
  !! groundwater head
  !<
  subroutine lak_calculate_cond_head(this, iconn, stage, head, vv)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: iconn
    real(DP), intent(in) :: stage
    real(DP), intent(in) :: head
    real(DP), intent(inout) :: vv
    ! -- local
    real(DP) :: ss
    real(DP) :: hh
    real(DP) :: topl
    real(DP) :: botl
    !
    topl = this%telev(iconn)
    botl = this%belev(iconn)
    ss = min(stage, topl)
    hh = min(head, topl)
    if (this%igwhcopt > 0) then
      vv = hh
    else if (this%inewton > 0) then
      vv = max(ss, hh)
    else
      vv = DHALF * (ss + hh)
    end if
  end subroutine lak_calculate_cond_head

  !> @brief Calculate the conductance for a lake connection at a provided stage
  !! and groundwater head
  !<
  subroutine lak_calculate_conn_conductance(this, ilak, iconn, stage, head, cond)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    integer(I4B), intent(in) :: iconn
    real(DP), intent(in) :: stage
    real(DP), intent(in) :: head
    real(DP), intent(inout) :: cond
    ! -- local
    integer(I4B) :: node
    !real(DP) :: ss
    !real(DP) :: hh
    real(DP) :: vv
    real(DP) :: topl
    real(DP) :: botl
    real(DP) :: sat
    real(DP) :: wa
    real(DP) :: vscratio
    !
    cond = DZERO
    vscratio = DONE
    topl = this%telev(iconn)
    botl = this%belev(iconn)
    call this%lak_calculate_cond_head(iconn, stage, head, vv)
    sat = sQuadraticSaturation(topl, botl, vv)
    ! vertical connection
    ! use full saturated conductance if top and bottom of the lake connection
    ! are equal
    if (this%ictype(iconn) == 0) then
      if (ABS(topl - botl) < DPREC) then
        sat = DONE
      end if
      ! horizontal connection
      ! use full saturated conductance if the connected cell is not convertible
    else if (this%ictype(iconn) == 1) then
      node = this%cellid(iconn)
      if (this%icelltype(node) == 0) then
        sat = DONE
      end if
      ! embedded connection
    else if (this%ictype(iconn) == 2 .or. this%ictype(iconn) == 3) then
      node = this%cellid(iconn)
      if (this%icelltype(node) == 0) then
        vv = this%telev(iconn)
        call this%lak_calculate_conn_warea(ilak, iconn, vv, vv, wa)
      else
        call this%lak_calculate_conn_warea(ilak, iconn, stage, head, wa)
      end if
      sat = wa
    end if
    !
    ! -- account for viscosity effects (if vsc active)
    if (this%ivsc == 1) then
      ! flow from lake to aquifer
      if (stage > head) then
        vscratio = this%viscratios(1, iconn)
        ! flow from aquifer to lake
      else
        vscratio = this%viscratios(2, iconn)
      end if
    end if
    cond = sat * this%satcond(iconn) * vscratio
  end subroutine lak_calculate_conn_conductance

  !> @brief Calculate the total groundwater-lake flow at a provided stage
  !<
  subroutine lak_calculate_exchange(this, ilak, stage, totflow)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(in) :: stage
    real(DP), intent(inout) :: totflow
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: igwfnode
    real(DP) :: flow
    real(DP) :: hgwf
    !
    totflow = DZERO
    do j = this%idxlakeconn(ilak), this%idxlakeconn(ilak + 1) - 1
      igwfnode = this%cellid(j)
      hgwf = this%xnew(igwfnode)
      call this%lak_calculate_conn_exchange(ilak, j, stage, hgwf, flow)
      totflow = totflow + flow
    end do
  end subroutine lak_calculate_exchange

  !> @brief Calculate the groundwater-lake flow at a provided stage and
  !! groundwater head
  !<
  subroutine lak_calculate_conn_exchange(this, ilak, iconn, stage, head, flow, &
                                         gwfhcof, gwfrhs)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    integer(I4B), intent(in) :: iconn
    real(DP), intent(in) :: stage
    real(DP), intent(in) :: head
    real(DP), intent(inout) :: flow
    real(DP), intent(inout), optional :: gwfhcof
    real(DP), intent(inout), optional :: gwfrhs
    ! -- local
    real(DP) :: botl
    real(DP) :: cond
    real(DP) :: ss
    real(DP) :: hh
    real(DP) :: gwfhcof0
    real(DP) :: gwfrhs0
    !
    flow = DZERO
    call this%lak_calculate_conn_conductance(ilak, iconn, stage, head, cond)
    botl = this%belev(iconn)
    !
    ! -- Set ss to stage or botl
    if (stage >= botl) then
      ss = stage
    else
      ss = botl
    end if
    !
    ! -- set hh to head or botl
    if (head >= botl) then
      hh = head
    else
      hh = botl
    end if
    !
    ! -- calculate flow, positive into lake
    flow = cond * (hh - ss)
    !
    ! -- Calculate gwfhcof and gwfrhs
    if (head >= botl) then
      gwfhcof0 = -cond
      gwfrhs0 = -cond * ss
    else
      gwfhcof0 = DZERO
      gwfrhs0 = flow
    end if
    !
    ! Add density contributions, if active
    if (this%idense /= 0) then
      call this%lak_calculate_density_exchange(iconn, stage, head, cond, botl, &
                                               flow, gwfhcof0, gwfrhs0)
    end if
    !
    ! -- If present update gwfhcof and gwfrhs
    if (present(gwfhcof)) gwfhcof = gwfhcof0
    if (present(gwfrhs)) gwfrhs = gwfrhs0
  end subroutine lak_calculate_conn_exchange

  !> @brief Calculate the groundwater-lake flow at a provided stage and
  !! groundwater head
  !<
  subroutine lak_estimate_conn_exchange(this, iflag, ilak, iconn, idry, stage, &
                                        head, flow, source, gwfhcof, gwfrhs)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: iflag
    integer(I4B), intent(in) :: ilak
    integer(I4B), intent(in) :: iconn
    integer(I4B), intent(inout) :: idry
    real(DP), intent(in) :: stage
    real(DP), intent(in) :: head
    real(DP), intent(inout) :: flow
    real(DP), intent(inout) :: source
    real(DP), intent(inout), optional :: gwfhcof
    real(DP), intent(inout), optional :: gwfrhs
    ! -- local
    real(DP) :: gwfhcof0, gwfrhs0
    !
    flow = DZERO
    idry = 0
    call this%lak_calculate_conn_exchange(ilak, iconn, stage, head, flow, &
                                          gwfhcof0, gwfrhs0)
    if (iflag == 1) then
      if (flow > DZERO) then
        source = source + flow
      end if
    else if (iflag == 2) then
      if (-flow > source) then
        flow = -source
        source = DZERO
        idry = 1
      else if (flow < DZERO) then
        source = source + flow
      end if
    end if
    !
    ! -- Set gwfhcof and gwfrhs if present
    if (present(gwfhcof)) gwfhcof = gwfhcof0
    if (present(gwfrhs)) gwfrhs = gwfrhs0
  end subroutine lak_estimate_conn_exchange

  !> @brief Calculate the storage change in a lake based on provided stages
  !! and a passed delt
  !<
  subroutine lak_calculate_storagechange(this, ilak, stage, stage0, delt, dvr)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(in) :: stage
    real(DP), intent(in) :: stage0
    real(DP), intent(in) :: delt
    real(DP), intent(inout) :: dvr
    ! -- local
    real(DP) :: v
    real(DP) :: v0
    !
    dvr = DZERO
    if (this%gwfiss /= 1) then
      call this%lak_calculate_vol(ilak, stage, v)
      call this%lak_calculate_vol(ilak, stage0, v0)
      dvr = (v0 - v) / delt
    end if
  end subroutine lak_calculate_storagechange

  !> @brief Calculate the rainfall for a lake
  !<
  subroutine lak_calculate_rainfall(this, ilak, stage, ra)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(in) :: stage
    real(DP), intent(inout) :: ra
    ! -- local
    integer(I4B) :: iconn
    real(DP) :: sa
    !
    ! -- rainfall
    iconn = this%idxlakeconn(ilak)
    if (this%ictype(iconn) == 2 .or. this%ictype(iconn) == 3) then
      sa = this%sareamax(ilak)
    else
      call this%lak_calculate_sarea(ilak, stage, sa)
    end if
    ra = this%rainfall(ilak) * sa
  end subroutine lak_calculate_rainfall

  !> @brief Calculate runoff to a lake
  !<
  subroutine lak_calculate_runoff(this, ilak, ro)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(inout) :: ro
    !
    ! -- runoff
    ro = this%runoff(ilak)
  end subroutine lak_calculate_runoff

  !> @brief Calculate specified inflow to a lake
  !<
  subroutine lak_calculate_inflow(this, ilak, qin)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(inout) :: qin
    !
    ! -- inflow to lake
    qin = this%inflow(ilak)
  end subroutine lak_calculate_inflow

  !> @brief Calculate the external flow terms to a lake
  !<
  subroutine lak_calculate_external(this, ilak, ex)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(inout) :: ex
    !
    ! -- If mover is active, add receiver water to rhs and
    !    store available water (as positive value)
    ex = DZERO
    if (this%imover == 1) then
      ex = this%pakmvrobj%get_qfrommvr(ilak)
    end if
  end subroutine lak_calculate_external

  !> @brief Calculate the withdrawal from a lake subject to an available volume
  !<
  subroutine lak_calculate_withdrawal(this, ilak, avail, wr)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(inout) :: avail
    real(DP), intent(inout) :: wr
    !
    ! -- withdrawals - limit to sum of inflows and available volume
    wr = this%withdrawal(ilak)
    if (wr > avail) then
      wr = -avail
    else
      if (wr > DZERO) then
        wr = -wr
      end if
    end if
    avail = avail + wr
  end subroutine lak_calculate_withdrawal

  !> @brief Calculate the evaporation from a lake at a provided stage subject
  !! to an available volume
  !<
  subroutine lak_calculate_evaporation(this, ilak, stage, avail, ev)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(in) :: stage
    real(DP), intent(inout) :: avail
    real(DP), intent(inout) :: ev
    ! -- local
    real(DP) :: sa
    !
    ! -- evaporation - limit to sum of inflows and available volume
    call this%lak_calculate_sarea(ilak, stage, sa)
    ev = sa * this%evaporation(ilak)
    if (ev > avail) then
      if (is_close(avail, DPREC)) then
        ev = DZERO
      else
        ev = -avail
      end if
    else
      ev = -ev
    end if
    avail = avail + ev
  end subroutine lak_calculate_evaporation

  !> @brief Calculate the outlet inflow to a lake
  !<
  subroutine lak_calculate_outlet_inflow(this, ilak, outinf)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(inout) :: outinf
    ! -- local
    integer(I4B) :: n
    !
    outinf = DZERO
    do n = 1, this%noutlets
      if (this%lakeout(n) == ilak) then
        outinf = outinf - this%simoutrate(n)
        if (this%imover == 1) then
          outinf = outinf - this%pakmvrobj%get_qtomvr(n)
        end if
      end if
    end do
  end subroutine lak_calculate_outlet_inflow

  !> @brief Calculate the outlet outflow from a lake
  !<
  subroutine lak_calculate_outlet_outflow(this, ilak, stage, avail, outoutf)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(in) :: stage
    real(DP), intent(inout) :: avail
    real(DP), intent(inout) :: outoutf
    ! -- local
    integer(I4B) :: n
    real(DP) :: g
    real(DP) :: d
    real(DP) :: c
    real(DP) :: gsm
    real(DP) :: rate
    !
    outoutf = DZERO
    do n = 1, this%noutlets
      if (this%lakein(n) == ilak) then
        rate = DZERO
        d = stage - this%outinvert(n)
        if (this%outdmax > DZERO) then
          if (d > this%outdmax) d = this%outdmax
        end if
        g = DGRAVITY * this%convlength * this%convtime * this%convtime
        select case (this%iouttype(n))
          ! specified rate
        case (0)
          rate = this%outrate(n)
          if (-rate > avail) then
            rate = -avail
          end if
          ! manning
        case (1)
          if (d > DZERO) then
            c = (this%convlength**DONETHIRD) * this%convtime
            gsm = DZERO
            if (this%outrough(n) > DZERO) then
              gsm = DONE / this%outrough(n)
            end if
            rate = -c * gsm * this%outwidth(n) * (d**DFIVETHIRDS) * &
                   sqrt(this%outslope(n))
          end if
          ! weir
        case (2)
          if (d > DZERO) then
            rate = -DTWOTHIRDS * DCD * this%outwidth(n) * d * &
                   sqrt(DTWO * g * d)
          end if
        end select
        this%simoutrate(n) = rate
        avail = avail + rate
        outoutf = outoutf + rate
      end if
    end do
  end subroutine lak_calculate_outlet_outflow

  !> @brief Get the outlet inflow to a lake from another lake
  !<
  subroutine lak_get_internal_inlet(this, ilak, outinf)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(inout) :: outinf
    ! -- local
    integer(I4B) :: n
    !
    outinf = DZERO
    do n = 1, this%noutlets
      if (this%lakeout(n) == ilak) then
        outinf = outinf - this%simoutrate(n)
        if (this%imover == 1) then
          outinf = outinf - this%pakmvrobj%get_qtomvr(n)
        end if
      end if
    end do
  end subroutine lak_get_internal_inlet

  !> @brief Get the outlet from a lake to another lake
  !<
  subroutine lak_get_internal_outlet(this, ilak, outoutf)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(inout) :: outoutf
    ! -- local
    integer(I4B) :: n
    !
    outoutf = DZERO
    do n = 1, this%noutlets
      if (this%lakein(n) == ilak) then
        if (this%lakeout(n) < 1) cycle
        outoutf = outoutf + this%simoutrate(n)
      end if
    end do
  end subroutine lak_get_internal_outlet

  !> @brief Get the outlet outflow from a lake to an external boundary
  !<
  subroutine lak_get_external_outlet(this, ilak, outoutf)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(inout) :: outoutf
    ! -- local
    integer(I4B) :: n
    !
    outoutf = DZERO
    do n = 1, this%noutlets
      if (this%lakein(n) == ilak) then
        if (this%lakeout(n) > 0) cycle
        outoutf = outoutf + this%simoutrate(n)
      end if
    end do
  end subroutine lak_get_external_outlet

  !> @brief Get the mover outflow from a lake to an external boundary
  !<
  subroutine lak_get_external_mover(this, ilak, outoutf)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(inout) :: outoutf
    ! -- local
    integer(I4B) :: n
    !
    outoutf = DZERO
    if (this%imover == 1) then
      do n = 1, this%noutlets
        if (this%lakein(n) == ilak) then
          if (this%lakeout(n) > 0) cycle
          outoutf = outoutf + this%pakmvrobj%get_qtomvr(n)
        end if
      end do
    end if
  end subroutine lak_get_external_mover

  !> @brief Get the mover outflow from a lake to another lake
  !<
  subroutine lak_get_internal_mover(this, ilak, outoutf)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(inout) :: outoutf
    ! -- local
    integer(I4B) :: n
    !
    outoutf = DZERO
    if (this%imover == 1) then
      do n = 1, this%noutlets
        if (this%lakein(n) == ilak) then
          if (this%lakeout(n) < 1) cycle
          outoutf = outoutf + this%pakmvrobj%get_qtomvr(n)
        end if
      end do
    end if
  end subroutine lak_get_internal_mover

  !> @brief Get the outlet to mover from a lake
  !<
  subroutine lak_get_outlet_tomover(this, ilak, outoutf)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(inout) :: outoutf
    ! -- local
    integer(I4B) :: n
    !
    outoutf = DZERO
    if (this%imover == 1) then
      do n = 1, this%noutlets
        if (this%lakein(n) == ilak) then
          outoutf = outoutf + this%pakmvrobj%get_qtomvr(n)
        end if
      end do
    end if
  end subroutine lak_get_outlet_tomover

  !> @brief Determine the stage from a provided volume
  !<
  subroutine lak_vol2stage(this, ilak, vol, stage)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(in) :: vol
    real(DP), intent(inout) :: stage
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: ibs
    real(DP) :: s0, s1, sm
    real(DP) :: v0, v1, vm
    real(DP) :: f0, f1, fm
    real(DP) :: sa
    real(DP) :: en0, en1
    real(DP) :: ds, ds0
    real(DP) :: denom
    !
    s0 = this%lakebot(ilak)
    call this%lak_calculate_vol(ilak, s0, v0)
    s1 = this%laketop(ilak)
    call this%lak_calculate_vol(ilak, s1, v1)
    ! -- zero volume
    if (vol <= v0) then
      stage = s0
      ! -- linear relation between stage and volume above top of lake
    else if (vol >= v1) then
      call this%lak_calculate_sarea(ilak, s1, sa)
      stage = s1 + (vol - v1) / sa
      ! -- use combination of secant and bisection
    else
      en0 = s0
      en1 = s1
      ! sm = s1  ! causes divide by zero in 1st line in secantbisection loop
      ! sm = s0  ! causes divide by zero in 1st line in secantbisection loop
      sm = DZERO
      f0 = vol - v0
      f1 = vol - v1
      ibs = 0
      secantbisection: do i = 1, 150
        denom = f1 - f0
        if (denom /= DZERO) then
          ds = f1 * (s1 - s0) / denom
        else
          ibs = 13
        end if
        if (i == 1) then
          ds0 = ds
        end if
        ! -- use bisection if end points are exceeded
        if (sm < en0 .or. sm > en1) ibs = 13
        ! -- use bisection if secant method stagnates or if
        !    ds exceeds previous ds - bisection would occur
        !    after conditions exceeded in 13 iterations
        if (ds * ds0 < DPREC .or. ABS(ds) > ABS(ds0)) ibs = ibs + 1
        if (ibs > 12) then
          ds = DHALF * (s1 - s0)
          ibs = 0
        end if
        sm = s1 - ds
        if (ABS(ds) < DEM6) then
          exit secantbisection
        end if
        call this%lak_calculate_vol(ilak, sm, vm)
        fm = vol - vm
        s0 = s1
        f0 = f1
        s1 = sm
        f1 = fm
        ds0 = ds
      end do secantbisection
      stage = sm
      if (ABS(ds) >= DEM6) then
        write (this%iout, '(1x,a,1x,i0,4(1x,a,1x,g15.6))') &
     &   'LAK_VOL2STAGE failed for lake', ilak, 'volume error =', fm, &
     &   'finding stage (', stage, ') for volume =', vol, &
     &    'final change in stage =', ds
      end if
    end if
  end subroutine lak_vol2stage

  !> @brief Determine if a valid lake or outlet number has been specified
  function lak_check_valid(this, itemno) result(ierr)
    ! -- modules
    use SimModule, only: store_error
    ! -- return
    integer(I4B) :: ierr
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    ! -- local
    integer(I4B) :: ival
    !
    ierr = 0
    ival = abs(itemno)
    if (itemno > 0) then
      if (ival < 1 .or. ival > this%nlakes) then
        write (errmsg, '(a,1x,i0,1x,a,1x,i0,a)') &
          'LAKENO', itemno, 'must be greater than 0 and less than or equal to', &
          this%nlakes, '.'
        call store_error(errmsg)
        ierr = 1
      end if
    else
      if (ival < 1 .or. ival > this%noutlets) then
        write (errmsg, '(a,1x,i0,1x,a,1x,i0,a)') &
          'IOUTLET', itemno, 'must be greater than 0 and less than or equal to', &
          this%noutlets, '.'
        call store_error(errmsg)
        ierr = 1
      end if
    end if
  end function lak_check_valid

  !> @brief Set a stress period attribute for lakweslls(itemno) using keywords
  !<
  subroutine lak_set_stressperiod(this, itemno)
    ! -- modules
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    use SimModule, only: store_error
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    ! -- local
    character(len=LINELENGTH) :: text
    character(len=LINELENGTH) :: caux
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    integer(I4B) :: ii
    integer(I4B) :: jj
    real(DP), pointer :: bndElem => null()
    !
    ! -- read line
    call this%parser%GetStringCaps(keyword)
    select case (keyword)
    case ('STATUS')
      ierr = this%lak_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetStringCaps(text)
      this%status(itemno) = text(1:8)
      if (text == 'CONSTANT') then
        this%iboundpak(itemno) = -1
      else if (text == 'INACTIVE') then
        this%iboundpak(itemno) = 0
      else if (text == 'ACTIVE') then
        this%iboundpak(itemno) = 1
      else
        write (errmsg, '(a,a)') &
          'Unknown '//trim(this%text)//' lak status keyword: ', text//'.'
        call store_error(errmsg)
      end if
    case ('STAGE')
      ierr = this%lak_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1 ! For STAGE
      bndElem => this%stage(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'STAGE')
    case ('RAINFALL')
      ierr = this%lak_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1 ! For RAINFALL
      bndElem => this%rainfall(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'RAINFALL')
      if (this%rainfall(itemno) < DZERO) then
        write (errmsg, '(a,i0,a,G0,a)') &
          'Lake ', itemno, ' was assigned a rainfall value of ', &
          this%rainfall(itemno), '. Rainfall must be positive.'
        call store_error(errmsg)
      end if
    case ('EVAPORATION')
      ierr = this%lak_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1 ! For EVAPORATION
      bndElem => this%evaporation(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'EVAPORATION')
      if (this%evaporation(itemno) < DZERO) then
        write (errmsg, '(a,i0,a,G0,a)') &
          'Lake ', itemno, ' was assigned an evaporation value of ', &
          this%evaporation(itemno), '. Evaporation must be positive.'
        call store_error(errmsg)
      end if
    case ('RUNOFF')
      ierr = this%lak_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1 ! For RUNOFF
      bndElem => this%runoff(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'RUNOFF')
      if (this%runoff(itemno) < DZERO) then
        write (errmsg, '(a,i0,a,G0,a)') &
          'Lake ', itemno, ' was assigned a runoff value of ', &
          this%runoff(itemno), '. Runoff must be positive.'
        call store_error(errmsg)
      end if
    case ('INFLOW')
      ierr = this%lak_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1 ! For specified INFLOW
      bndElem => this%inflow(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'INFLOW')
      if (this%inflow(itemno) < DZERO) then
        write (errmsg, '(a,i0,a,G0,a)') &
          'Lake ', itemno, ' was assigned an inflow value of ', &
          this%inflow(itemno), '. Inflow must be positive.'
        call store_error(errmsg)
      end if
    case ('WITHDRAWAL')
      ierr = this%lak_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1 ! For specified WITHDRAWAL
      bndElem => this%withdrawal(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'WITHDRAWAL')
      if (this%withdrawal(itemno) < DZERO) then
        write (errmsg, '(a,i0,a,G0,a)') &
          'Lake ', itemno, ' was assigned a withdrawal value of ', &
          this%withdrawal(itemno), '. Withdrawal must be positive.'
        call store_error(errmsg)
      end if
    case ('RATE')
      ierr = this%lak_check_valid(-itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1 ! For specified OUTLET RATE
      bndElem => this%outrate(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'RATE')
    case ('INVERT')
      ierr = this%lak_check_valid(-itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1 ! For OUTLET INVERT
      bndElem => this%outinvert(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'INVERT')
    case ('WIDTH')
      ierr = this%lak_check_valid(-itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1 ! For OUTLET WIDTH
      bndElem => this%outwidth(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'WIDTH')
    case ('ROUGH')
      ierr = this%lak_check_valid(-itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1 ! For OUTLET ROUGHNESS
      bndElem => this%outrough(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'ROUGH')
    case ('SLOPE')
      ierr = this%lak_check_valid(-itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1 ! For OUTLET SLOPE
      bndElem => this%outslope(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'SLOPE')
    case ('AUXILIARY')
      ierr = this%lak_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetStringCaps(caux)
      do jj = 1, this%naux
        if (trim(adjustl(caux)) /= trim(adjustl(this%auxname(jj)))) cycle
        call this%parser%GetString(text)
        ii = itemno
        bndElem => this%lauxvar(jj, ii)
        call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                           this%packName, 'AUX', &
                                           this%tsManager, this%iprpak, &
                                           this%auxname(jj))
        exit
      end do
    case default
      write (errmsg, '(2a)') &
        'Unknown '//trim(this%text)//' lak data keyword: ', &
        trim(keyword)//'.'
    end select
    !
    ! -- Return
999 return
  end subroutine lak_set_stressperiod

  !> @brief Issue a parameter error for lakweslls(ilak)
  !!
  !! Read itmp and new boundaries if itmp > 0
  !<
  subroutine lak_set_attribute_error(this, ilak, keyword, msg)
    ! -- modules
    use SimModule, only: store_error
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: ilak
    character(len=*), intent(in) :: keyword
    character(len=*), intent(in) :: msg
    !
    if (len(msg) == 0) then
      write (errmsg, '(a,1x,a,1x,i0,1x,a)') &
        keyword, ' for LAKE', ilak, 'has already been set.'
    else
      write (errmsg, '(a,1x,a,1x,i0,1x,a)') keyword, ' for LAKE', ilak, msg
    end if
    call store_error(errmsg)
  end subroutine lak_set_attribute_error

  !> @brief Set options specific to LakType
  !!
  !! lak_options overrides BndType%bnd_options
  !<
  subroutine lak_options(this, option, found)
    ! -- modules
    use ConstantsModule, only: MAXCHARLEN, DZERO, MNORMAL
    use OpenSpecModule, only: access, form
    use SimModule, only: store_error
    use InputOutputModule, only: urword, getunit, assign_iounit, openfile
    ! -- dummy
    class(LakType), intent(inout) :: this
    character(len=*), intent(inout) :: option
    logical(LGP), intent(inout) :: found
    ! -- local
    character(len=MAXCHARLEN) :: fname, keyword
    real(DP) :: r
    ! -- formats
    character(len=*), parameter :: fmtlengthconv = &
      &"(4x, 'LENGTH CONVERSION VALUE (',g15.7,') SPECIFIED.')"
    character(len=*), parameter :: fmttimeconv = &
      &"(4x, 'TIME CONVERSION VALUE (',g15.7,') SPECIFIED.')"
    character(len=*), parameter :: fmtoutdmax = &
      &"(4x, 'MAXIMUM OUTLET WATER DEPTH (',g15.7,') SPECIFIED.')"
    character(len=*), parameter :: fmtlakeopt = &
      &"(4x, 'LAKE ', a, ' VALUE (',g15.7,') SPECIFIED.')"
    character(len=*), parameter :: fmtlakbin = &
      "(4x, 'LAK ', 1x, a, 1x, ' WILL BE SAVED TO FILE: ', &
      &a, /4x, 'OPENED ON UNIT: ', I0)"
    character(len=*), parameter :: fmtiter = &
      &"(4x, 'MAXIMUM LAK ITERATION VALUE (',i0,') SPECIFIED.')"
    character(len=*), parameter :: fmtdmaxchg = &
      &"(4x, 'MAXIMUM STAGE CHANGE VALUE (',g0,') SPECIFIED.')"
    !
    found = .true.
    select case (option)
    case ('PRINT_STAGE')
      this%iprhed = 1
      write (this%iout, '(4x,a)') trim(adjustl(this%text))// &
        ' STAGES WILL BE PRINTED TO LISTING FILE.'
    case ('STAGE')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        this%istageout = getunit()
        call openfile(this%istageout, this%iout, fname, 'DATA(BINARY)', &
                      form, access, 'REPLACE', mode_opt=MNORMAL)
        write (this%iout, fmtlakbin) 'STAGE', trim(adjustl(fname)), &
          this%istageout
      else
        call store_error('OPTIONAL STAGE KEYWORD MUST BE FOLLOWED BY FILEOUT')
      end if
    case ('BUDGET')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        call assign_iounit(this%ibudgetout, this%inunit, "BUDGET fileout")
        call openfile(this%ibudgetout, this%iout, fname, 'DATA(BINARY)', &
                      form, access, 'REPLACE', mode_opt=MNORMAL)
        write (this%iout, fmtlakbin) 'BUDGET', trim(adjustl(fname)), &
          this%ibudgetout
      else
        call store_error('OPTIONAL BUDGET KEYWORD MUST BE FOLLOWED BY FILEOUT')
      end if
    case ('BUDGETCSV')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        call assign_iounit(this%ibudcsv, this%inunit, "BUDGETCSV fileout")
        call openfile(this%ibudcsv, this%iout, fname, 'CSV', &
                      filstat_opt='REPLACE')
        write (this%iout, fmtlakbin) 'BUDGET CSV', trim(adjustl(fname)), &
          this%ibudcsv
      else
        call store_error('OPTIONAL BUDGETCSV KEYWORD MUST BE FOLLOWED BY &
          &FILEOUT')
      end if
    case ('PACKAGE_CONVERGENCE')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        this%ipakcsv = getunit()
        call openfile(this%ipakcsv, this%iout, fname, 'CSV', &
                      filstat_opt='REPLACE', mode_opt=MNORMAL)
        write (this%iout, fmtlakbin) 'PACKAGE_CONVERGENCE', &
          trim(adjustl(fname)), this%ipakcsv
      else
        call store_error('OPTIONAL PACKAGE_CONVERGENCE KEYWORD MUST BE '// &
                         'FOLLOWED BY FILEOUT')
      end if
    case ('MOVER')
      this%imover = 1
      write (this%iout, '(4x,A)') 'MOVER OPTION ENABLED'
    case ('LENGTH_CONVERSION')
      this%convlength = this%parser%GetDouble()
      write (this%iout, fmtlengthconv) this%convlength
    case ('TIME_CONVERSION')
      this%convtime = this%parser%GetDouble()
      write (this%iout, fmttimeconv) this%convtime
    case ('SURFDEP')
      r = this%parser%GetDouble()
      if (r < DZERO) then
        r = DZERO
      end if
      this%surfdep = r
      write (this%iout, fmtlakeopt) 'SURFDEP', this%surfdep
    case ('MAXIMUM_ITERATIONS')
      this%maxlakit = this%parser%GetInteger()
      write (this%iout, fmtiter) this%maxlakit
    case ('MAXIMUM_STAGE_CHANGE')
      r = this%parser%GetDouble()
      this%dmaxchg = r
      this%delh = DP999 * r
      write (this%iout, fmtdmaxchg) this%dmaxchg
      !
      ! -- right now these are options that are only available in the
      !    development version and are not included in the documentation.
      !    These options are only available when IDEVELOPMODE in
      !    constants module is set to 1
    case ('DEV_GROUNDWATER_HEAD_CONDUCTANCE')
      call this%parser%DevOpt()
      this%igwhcopt = 1
      write (this%iout, '(4x,a)') &
        'CONDUCTANCE FOR HORIZONTAL CONNECTIONS WILL BE CALCULATED &
        &USING THE GROUNDWATER HEAD'
    case ('DEV_MAXIMUM_OUTLET_DEPTH')
      call this%parser%DevOpt()
      this%outdmax = this%parser%GetDouble()
      write (this%iout, fmtoutdmax) this%outdmax
    case ('DEV_NO_FINAL_CHECK')
      call this%parser%DevOpt()
      this%iconvchk = 0
      write (this%iout, '(4x,a)') &
        'A FINAL CONVERGENCE CHECK OF THE CHANGE IN LAKE STAGES &
        &WILL NOT BE MADE'
    case ('DEV_NO_FINAL_RESIDUAL_CHECK')
      call this%parser%DevOpt()
      this%iconvresidchk = 0
      write (this%iout, '(4x,a)') &
        'A FINAL CONVERGENCE CHECK OF THE CHANGE IN LAKE RESIDUALS &
        &WILL NOT BE MADE'
    case ('DEV_MAXIMUM_PERCENT_DIFFERENCE')
      call this%parser%DevOpt()
      r = this%parser%GetDouble()
      if (r < DZERO) then
        r = DEM1
      end if
      this%pdmax = r
      write (this%iout, fmtlakeopt) 'MAXIMUM_PERCENT_DIFFERENCE', this%pdmax
    case default
      !
      ! -- No options found
      found = .false.
    end select
  end subroutine lak_options

  !> @brief Allocate and Read
  !!
  !! Create new LAK package and point bndobj to the new package
  !<
  subroutine lak_ar(this)
    ! -- dummy
    class(LakType), intent(inout) :: this
    !
    call this%obs%obs_ar()
    !
    ! -- Allocate arrays in LAK and in package superclass
    call this%lak_allocate_arrays()
    !
    ! -- read optional initial package parameters
    call this%read_initial_attr()
    !
    ! -- setup pakmvrobj
    if (this%imover /= 0) then
      allocate (this%pakmvrobj)
      call this%pakmvrobj%ar(this%noutlets, this%nlakes, this%memoryPath)
    end if
  end subroutine lak_ar

  !> @brief Read and Prepare
  !!
  !! Read itmp and read new boundaries if itmp > 0
  !<
  subroutine lak_rp(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use TdisModule, only: kper, nper
    use SimModule, only: store_error, count_errors
    ! -- dummy
    class(LakType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: line
    character(len=LINELENGTH) :: text
    logical(LGP) :: isfound
    logical(LGP) :: endOfBlock
    integer(I4B) :: ierr
    integer(I4B) :: node
    integer(I4B) :: n
    integer(I4B) :: itemno
    integer(I4B) :: j
    ! -- formats
    character(len=*), parameter :: fmtblkerr = &
      &"('Looking for BEGIN PERIOD iper.  Found ', a, ' instead.')"
    character(len=*), parameter :: fmtlsp = &
      &"(1X,/1X,'REUSING ',A,'S FROM LAST STRESS PERIOD')"
    !
    ! -- set nbound to maxbound
    this%nbound = this%maxbound
    !
    ! -- Set ionper to the stress period number for which a new block of data
    !    will be read.
    if (this%inunit == 0) return
    !
    ! -- get stress period data
    if (this%ionper < kper) then
      !
      ! -- get period block
      call this%parser%GetBlock('PERIOD', isfound, ierr, &
                                supportOpenClose=.true., &
                                blockRequired=.false.)
      if (isfound) then
        !
        ! -- read ionper and check for increasing period numbers
        call this%read_check_ionper()
      else
        !
        ! -- PERIOD block not found
        if (ierr < 0) then
          ! -- End of file found; data applies for remainder of simulation.
          this%ionper = nper + 1
        else
          ! -- Found invalid block
          call this%parser%GetCurrentLine(line)
          write (errmsg, fmtblkerr) adjustl(trim(line))
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
      end if
    end if
    !
    ! -- Read data if ionper == kper
    if (this%ionper == kper) then
      !
      ! -- setup table for period data
      if (this%iprpak /= 0) then
        !
        ! -- reset the input table object
        title = trim(adjustl(this%text))//' PACKAGE ('// &
                trim(adjustl(this%packName))//') DATA FOR PERIOD'
        write (title, '(a,1x,i6)') trim(adjustl(title)), kper
        call table_cr(this%inputtab, this%packName, title)
        call this%inputtab%table_df(1, 4, this%iout, finalize=.FALSE.)
        text = 'NUMBER'
        call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
        text = 'KEYWORD'
        call this%inputtab%initialize_column(text, 20, alignment=TABLEFT)
        do n = 1, 2
          write (text, '(a,1x,i6)') 'VALUE', n
          call this%inputtab%initialize_column(text, 15, alignment=TABCENTER)
        end do
      end if
      !
      ! -- read the data
      this%check_attr = 1
      stressperiod: do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        !
        ! -- get lake or outlet number
        itemno = this%parser%GetInteger()
        !
        ! -- read data from the rest of the line
        call this%lak_set_stressperiod(itemno)
        !
        ! -- write line to table
        if (this%iprpak /= 0) then
          call this%parser%GetCurrentLine(line)
          call this%inputtab%line_to_columns(line)
        end if
      end do stressperiod
      !
      if (this%iprpak /= 0) then
        call this%inputtab%finalize_table()
      end if
      !
      ! -- using stress period data from the previous stress period
    else
      write (this%iout, fmtlsp) trim(this%filtyp)
    end if
    !
    ! -- write summary of lake stress period error messages
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- fill bound array with lake stage, conductance, and bottom elevation
    do n = 1, this%nlakes
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        node = this%cellid(j)
        this%nodelist(j) = node
        this%bound(1, j) = this%xnewpak(n)
        this%bound(2, j) = this%satcond(j)
        this%bound(3, j) = this%belev(j)
      end do
    end do
    !
    ! -- copy lakein into iprmap so mvr budget contains lake instead of outlet
    if (this%imover == 1) then
      do n = 1, this%noutlets
        this%pakmvrobj%iprmap(n) = this%lakein(n)
      end do
    end if
  end subroutine lak_rp

  !> @brief Add package connection to matrix
  !<
  subroutine lak_ad(this)
    ! -- modules
    use SimVariablesModule, only: iFailedStepRetry
    ! -- dummy
    class(LakType) :: this
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: j
    integer(I4B) :: iaux
    !
    ! -- Advance the time series
    call this%TsManager%ad()
    !
    ! -- update auxiliary variables by copying from the derived-type time
    !    series variable into the bndpackage auxvar variable so that this
    !    information is properly written to the GWF budget file
    if (this%naux > 0) then
      do n = 1, this%nlakes
        do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
          do iaux = 1, this%naux
            if (this%noupdateauxvar(iaux) /= 0) cycle
            this%auxvar(iaux, j) = this%lauxvar(iaux, n)
          end do
        end do
      end do
    end if
    !
    ! -- Update or restore state
    if (iFailedStepRetry == 0) then
      !
      ! -- copy xnew into xold and set xnewpak to stage for
      !    constant stage lakes
      do n = 1, this%nlakes
        this%xoldpak(n) = this%xnewpak(n)
        this%stageiter(n) = this%xnewpak(n)
        if (this%iboundpak(n) < 0) then
          this%xnewpak(n) = this%stage(n)
        end if
        this%seep0(n) = DZERO
      end do
    else
      !
      ! -- copy xold back into xnew as this is a
      !    retry of this time step
      do n = 1, this%nlakes
        this%xnewpak(n) = this%xoldpak(n)
        this%stageiter(n) = this%xnewpak(n)
        if (this%iboundpak(n) < 0) then
          this%xnewpak(n) = this%stage(n)
        end if
        this%seep0(n) = DZERO
      end do
    end if
    !
    ! -- pakmvrobj ad
    if (this%imover == 1) then
      call this%pakmvrobj%ad()
    end if
    !
    ! -- For each observation, push simulated value and corresponding
    !    simulation time from "current" to "preceding" and reset
    !    "current" value.
    call this%obs%obs_ad()
  end subroutine lak_ad

  !> @brief Formulate the HCOF and RHS terms
  !!
  !! Skip if no lakes, otherwise calculate hcof and rhs
  !<
  subroutine lak_cf(this)
    ! -- dummy
    class(LakType) :: this
    ! -- local
    integer(I4B) :: j, n
    integer(I4B) :: igwfnode
    real(DP) :: hlak, bottom_lake
    !
    ! -- save groundwater seepage for lake solution
    do n = 1, this%nlakes
      this%seep0(n) = this%seep(n)
    end do
    !
    ! -- save variables for convergence check
    do n = 1, this%nlakes
      this%s0(n) = this%xnewpak(n)
      call this%lak_calculate_exchange(n, this%s0(n), this%qgwf0(n))
    end do
    !
    ! -- find highest active cell
    do n = 1, this%nlakes
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        ! -- skip horizontal connections
        if (this%ictype(j) /= 0) then
          cycle
        end if
        igwfnode = this%nodesontop(j)
        if (this%ibound(igwfnode) == 0) then
          call this%dis%highest_active(igwfnode, this%ibound)
        end if
        this%nodelist(j) = igwfnode
        this%cellid(j) = igwfnode
      end do
    end do
    !
    ! -- reset ibound for cells where lake stage is above the bottom
    !    of the lake in the cell or the lake is inactive - only applied to
    !    vertical connections
    do n = 1, this%nlakes
      !
      hlak = this%xnewpak(n)
      !
      ! -- Go through lake connections
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        !
        ! -- assign gwf node number
        igwfnode = this%cellid(j)
        !
        ! -- skip inactive or constant head GWF cells
        if (this%ibound(igwfnode) < 1) then
          cycle
        end if
        !
        ! -- skip horizontal connections
        if (this%ictype(j) /= 0) then
          cycle
        end if
        !
        ! -- skip embedded lakes
        if (this%ictype(j) == 2 .or. this%ictype(j) == 3) then
          cycle
        end if
        !
        ! -- Mark ibound for wet lakes or inactive lakes; reset to 1 otherwise
        bottom_lake = this%belev(j)
        if (hlak > bottom_lake .or. this%iboundpak(n) == 0) then
          this%ibound(igwfnode) = IWETLAKE
        else
          this%ibound(igwfnode) = 1
        end if
      end do
      !
    end do
    !
    ! -- Store the lake stage and cond in bound array for other
    !    packages, such as the BUY package
    call this%lak_bound_update()
  end subroutine lak_cf

  !> @brief Copy rhs and hcof into solution rhs and amat
  !<
  subroutine lak_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(LakType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: j, n
    integer(I4B) :: igwfnode
    integer(I4B) :: ipossymd
    !
    ! -- pakmvrobj fc
    if (this%imover == 1) then
      call this%pakmvrobj%fc()
    end if
    !
    ! -- make a stab at a solution
    call this%lak_solve()
    !
    ! -- add terms to the gwf matrix
    do n = 1, this%nlakes
      if (this%iboundpak(n) == 0) cycle
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        igwfnode = this%cellid(j)
        if (this%ibound(igwfnode) < 1) cycle
        ipossymd = idxglo(ia(igwfnode))
        call matrix_sln%add_value_pos(ipossymd, this%hcof(j))
        rhs(igwfnode) = rhs(igwfnode) + this%rhs(j)
      end do
    end do
  end subroutine lak_fc

  !> @brief Fill newton terms
  !<
  subroutine lak_fn(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(LakType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: j, n
    integer(I4B) :: ipos
    integer(I4B) :: igwfnode
    integer(I4B) :: idry
    real(DP) :: hlak
    real(DP) :: avail
    real(DP) :: ra
    real(DP) :: ro
    real(DP) :: qinf
    real(DP) :: ex
    real(DP) :: head
    real(DP) :: q
    real(DP) :: q1
    real(DP) :: rterm
    real(DP) :: drterm
    !
    do n = 1, this%nlakes
      if (this%iboundpak(n) == 0) cycle
      hlak = this%xnewpak(n)
      call this%lak_calculate_available(n, hlak, avail, &
                                        ra, ro, qinf, ex, this%delh)
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        igwfnode = this%cellid(j)
        ipos = ia(igwfnode)
        head = this%xnew(igwfnode)
        if (-this%hcof(j) > DZERO) then
          if (this%ibound(igwfnode) > 0) then
            ! -- estimate lake-aquifer exchange with perturbed groundwater head
            !    exchange is relative to the lake
            !avail = DEP20
            call this%lak_estimate_conn_exchange(2, n, j, idry, hlak, &
                                                 head + this%delh, q1, avail)
            q1 = -q1
            ! -- calculate unperturbed lake-aquifer exchange
            q = this%hcof(j) * head - this%rhs(j)
            ! -- calculate rterm
            rterm = this%hcof(j) * head
            ! -- calculate derivative
            drterm = (q1 - q) / this%delh
            ! -- add terms to convert conductance formulation into
            !    newton-raphson formulation
            call matrix_sln%add_value_pos(idxglo(ipos), drterm - this%hcof(j))
            rhs(igwfnode) = rhs(igwfnode) - rterm + drterm * head
          end if
        end if
      end do
    end do
  end subroutine lak_fn

  !> @brief Final convergence check for package
  !<
  subroutine lak_cc(this, innertot, kiter, iend, icnvgmod, cpak, ipak, dpak)
    ! -- modules
    use TdisModule, only: totim, kstp, kper, delt
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: innertot
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: iend
    integer(I4B), intent(in) :: icnvgmod
    character(len=LENPAKLOC), intent(inout) :: cpak
    integer(I4B), intent(inout) :: ipak
    real(DP), intent(inout) :: dpak
    ! -- local
    character(len=LENPAKLOC) :: cloc
    character(len=LINELENGTH) :: tag
    integer(I4B) :: icheck
    integer(I4B) :: ipakfail
    integer(I4B) :: locdhmax
    integer(I4B) :: locresidmax
    integer(I4B) :: locdgwfmax
    integer(I4B) :: locdqoutmax
    integer(I4B) :: locdqfrommvrmax
    integer(I4B) :: ntabrows
    integer(I4B) :: ntabcols
    integer(I4B) :: n
    real(DP) :: q
    real(DP) :: q0
    real(DP) :: qtolfact
    real(DP) :: area
    real(DP) :: gwf0
    real(DP) :: gwf
    real(DP) :: dh
    real(DP) :: resid
    real(DP) :: dgwf
    real(DP) :: hlak0
    real(DP) :: hlak
    real(DP) :: qout0
    real(DP) :: qout
    real(DP) :: dqout
    real(DP) :: inf
    real(DP) :: ra
    real(DP) :: ro
    real(DP) :: qinf
    real(DP) :: ex
    real(DP) :: dhmax
    real(DP) :: residmax
    real(DP) :: dgwfmax
    real(DP) :: dqoutmax
    real(DP) :: dqfrommvr
    real(DP) :: dqfrommvrmax
    !
    ! -- initialize local variables
    icheck = this%iconvchk
    ipakfail = 0
    locdhmax = 0
    locresidmax = 0
    locdgwfmax = 0
    locdqoutmax = 0
    locdqfrommvrmax = 0
    dhmax = DZERO
    residmax = DZERO
    dgwfmax = DZERO
    dqoutmax = DZERO
    dqfrommvrmax = DZERO
    !
    ! -- if not saving package convergence data on check convergence if
    !    the model is considered converged
    if (this%ipakcsv == 0) then
      if (icnvgmod == 0) then
        icheck = 0
      end if
      !
      ! -- saving package convergence data
    else
      !
      ! -- header for package csv
      if (.not. associated(this%pakcsvtab)) then
        !
        ! -- determine the number of columns and rows
        ntabrows = 1
        ntabcols = 11
        if (this%noutlets > 0) then
          ntabcols = ntabcols + 2
        end if
        if (this%imover == 1) then
          ntabcols = ntabcols + 2
        end if
        !
        ! -- setup table
        call table_cr(this%pakcsvtab, this%packName, '')
        call this%pakcsvtab%table_df(ntabrows, ntabcols, this%ipakcsv, &
                                     lineseparator=.FALSE., separator=',', &
                                     finalize=.FALSE.)
        !
        ! -- add columns to package csv
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
        tag = 'residmax'
        call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
        tag = 'residmax_loc'
        call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
        tag = 'dgwfmax'
        call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
        tag = 'dgwfmax_loc'
        call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
        if (this%noutlets > 0) then
          tag = 'dqoutmax'
          call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
          tag = 'dqoutmax_loc'
          call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
        end if
        if (this%imover == 1) then
          tag = 'dqfrommvrmax'
          call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
          tag = 'dqfrommvrmax_loc'
          call this%pakcsvtab%initialize_column(tag, 16, alignment=TABLEFT)
        end if
      end if
    end if
    !
    ! -- perform package convergence check
    if (icheck /= 0) then
      final_check: do n = 1, this%nlakes
        if (this%iboundpak(n) < 1) cycle
        !
        ! -- set previous and current lake stage
        hlak0 = this%s0(n)
        hlak = this%xnewpak(n)
        !
        ! -- stage difference
        dh = hlak0 - hlak
        !
        ! -- calculate surface area
        call this%lak_calculate_sarea(n, hlak, area)
        !
        ! -- set the Q to length factor
        if (area > DZERO) then
          qtolfact = delt / area
        else
          qtolfact = DZERO
        end if
        !
        ! -- difference in the residual
        call this%lak_calculate_residual(n, hlak, resid)
        resid = resid * qtolfact
        !
        ! -- change in gwf exchange
        dgwf = DZERO
        if (area > DZERO) then
          gwf0 = this%qgwf0(n)
          call this%lak_calculate_exchange(n, hlak, gwf)
          dgwf = (gwf0 - gwf) * qtolfact
        end if
        !
        ! -- change in outflows
        dqout = DZERO
        if (this%noutlets > 0) then
          if (area > DZERO) then
            call this%lak_calculate_available(n, hlak0, inf, ra, ro, qinf, ex)
            call this%lak_calculate_outlet_outflow(n, hlak0, inf, qout0)
            call this%lak_calculate_available(n, hlak, inf, ra, ro, qinf, ex)
            call this%lak_calculate_outlet_outflow(n, hlak, inf, qout)
            dqout = (qout0 - qout) * qtolfact
          end if
        end if
        !
        ! -- q from mvr
        dqfrommvr = DZERO
        if (this%imover == 1) then
          q = this%pakmvrobj%get_qfrommvr(n)
          q0 = this%pakmvrobj%get_qfrommvr0(n)
          dqfrommvr = qtolfact * (q0 - q)
        end if
        !
        ! -- evaluate magnitude of differences
        if (n == 1) then
          locdhmax = n
          dhmax = dh
          locdgwfmax = n
          residmax = resid
          locresidmax = n
          dgwfmax = dgwf
          locdqoutmax = n
          dqoutmax = dqout
          dqfrommvrmax = dqfrommvr
          locdqfrommvrmax = n
        else
          if (abs(dh) > abs(dhmax)) then
            locdhmax = n
            dhmax = dh
          end if
          if (abs(resid) > abs(residmax)) then
            locresidmax = n
            residmax = resid
          end if
          if (abs(dgwf) > abs(dgwfmax)) then
            locdgwfmax = n
            dgwfmax = dgwf
          end if
          if (abs(dqout) > abs(dqoutmax)) then
            locdqoutmax = n
            dqoutmax = dqout
          end if
          if (ABS(dqfrommvr) > abs(dqfrommvrmax)) then
            dqfrommvrmax = dqfrommvr
            locdqfrommvrmax = n
          end if
        end if
      end do final_check
      !
      ! -- set dpak and cpak
      if (ABS(dhmax) > abs(dpak)) then
        ipak = locdhmax
        dpak = dhmax
        write (cloc, "(a,'-',a)") &
          trim(this%packName), 'stage'
        cpak = trim(cloc)
      end if
      if (ABS(residmax) > abs(dpak)) then
        ipak = locresidmax
        dpak = residmax
        write (cloc, "(a,'-',a)") &
          trim(this%packName), 'residual'
        cpak = trim(cloc)
      end if
      if (ABS(dgwfmax) > abs(dpak)) then
        ipak = locdgwfmax
        dpak = dgwfmax
        write (cloc, "(a,'-',a)") &
          trim(this%packName), 'gwf'
        cpak = trim(cloc)
      end if
      if (this%noutlets > 0) then
        if (ABS(dqoutmax) > abs(dpak)) then
          ipak = locdqoutmax
          dpak = dqoutmax
          write (cloc, "(a,'-',a)") &
            trim(this%packName), 'outlet'
          cpak = trim(cloc)
        end if
      end if
      if (this%imover == 1) then
        if (ABS(dqfrommvrmax) > abs(dpak)) then
          ipak = locdqfrommvrmax
          dpak = dqfrommvrmax
          write (cloc, "(a,'-',a)") trim(this%packName), 'qfrommvr'
          cpak = trim(cloc)
        end if
      end if
      !
      ! -- write convergence data to package csv
      if (this%ipakcsv /= 0) then
        !
        ! -- write the data
        call this%pakcsvtab%add_term(innertot)
        call this%pakcsvtab%add_term(totim)
        call this%pakcsvtab%add_term(kper)
        call this%pakcsvtab%add_term(kstp)
        call this%pakcsvtab%add_term(kiter)
        call this%pakcsvtab%add_term(dhmax)
        call this%pakcsvtab%add_term(locdhmax)
        call this%pakcsvtab%add_term(residmax)
        call this%pakcsvtab%add_term(locresidmax)
        call this%pakcsvtab%add_term(dgwfmax)
        call this%pakcsvtab%add_term(locdgwfmax)
        if (this%noutlets > 0) then
          call this%pakcsvtab%add_term(dqoutmax)
          call this%pakcsvtab%add_term(locdqoutmax)
        end if
        if (this%imover == 1) then
          call this%pakcsvtab%add_term(dqfrommvrmax)
          call this%pakcsvtab%add_term(locdqfrommvrmax)
        end if
        !
        ! -- finalize the package csv
        if (iend == 1) then
          call this%pakcsvtab%finalize_table()
        end if
      end if
    end if
  end subroutine lak_cc

  !> @brief Calculate flows
  !<
  subroutine lak_cq(this, x, flowja, iadv)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(LakType), intent(inout) :: this
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    integer(I4B), optional, intent(in) :: iadv
    ! -- local
    real(DP) :: rrate
    real(DP) :: chratin, chratout
    ! -- for budget
    integer(I4B) :: j, n
    real(DP) :: hlak
    real(DP) :: v0, v1
    !
    call this%lak_solve(update=.false.)
    !
    ! -- call base functionality in bnd_cq.  This will calculate lake-gwf flows
    !    and put them into this%simvals
    call this%BndType%bnd_cq(x, flowja, iadv=1)
    !
    ! -- calculate several budget terms
    chratin = DZERO
    chratout = DZERO
    do n = 1, this%nlakes
      this%chterm(n) = DZERO
      if (this%iboundpak(n) == 0) cycle
      hlak = this%xnewpak(n)
      call this%lak_calculate_vol(n, hlak, v1)
      !
      ! -- add budget terms for active lakes
      if (this%iboundpak(n) /= 0) then
        !
        ! -- rainfall
        rrate = this%precip(n)
        call this%lak_accumulate_chterm(n, rrate, chratin, chratout)
        !
        ! -- evaporation
        rrate = this%evap(n)
        call this%lak_accumulate_chterm(n, rrate, chratin, chratout)
        !
        ! -- runoff
        rrate = this%runoff(n)
        call this%lak_accumulate_chterm(n, rrate, chratin, chratout)
        !
        ! -- inflow
        rrate = this%inflow(n)
        call this%lak_accumulate_chterm(n, rrate, chratin, chratout)
        !
        ! -- withdrawals
        rrate = this%withr(n)
        call this%lak_accumulate_chterm(n, rrate, chratin, chratout)
        !
        ! -- add lake storage changes
        rrate = DZERO
        if (this%iboundpak(n) > 0) then
          if (this%gwfiss /= 1) then
            call this%lak_calculate_vol(n, this%xoldpak(n), v0)
            rrate = -(v1 - v0) / delt
            call this%lak_accumulate_chterm(n, rrate, chratin, chratout)
          end if
        end if
        this%qsto(n) = rrate
        !
        ! -- add external outlets
        call this%lak_get_external_outlet(n, rrate)
        call this%lak_accumulate_chterm(n, rrate, chratin, chratout)
        !
        ! -- add mover terms
        if (this%imover == 1) then
          if (this%iboundpak(n) /= 0) then
            rrate = this%pakmvrobj%get_qfrommvr(n)
          else
            rrate = DZERO
          end if
          call this%lak_accumulate_chterm(n, rrate, chratin, chratout)
        end if
      end if
    end do
    !
    ! -- gwf flow and constant flow to lake
    do n = 1, this%nlakes
      rrate = DZERO
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        ! simvals is from aquifer perspective, and so it is positive
        ! for flow into the aquifer.  Need to switch sign for lake
        ! perspective.
        rrate = -this%simvals(j)
        this%qleak(j) = rrate
        if (this%iboundpak(n) /= 0) then
          call this%lak_accumulate_chterm(n, rrate, chratin, chratout)
        end if
      end do
    end do
    !
    ! -- fill the budget object
    call this%lak_fill_budobj()
  end subroutine lak_cq

  !> @brief Output LAK package flow terms
  !<
  subroutine lak_ot_package_flows(this, icbcfl, ibudfl)
    use TdisModule, only: kstp, kper, delt, pertim, totim
    class(LakType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    integer(I4B) :: ibinun
    !
    ! -- write the flows from the budobj
    ibinun = 0
    if (this%ibudgetout /= 0) then
      ibinun = this%ibudgetout
    end if
    if (icbcfl == 0) ibinun = 0
    if (ibinun > 0) then
      call this%budobj%save_flows(this%dis, ibinun, kstp, kper, delt, &
                                  pertim, totim, this%iout)
    end if
    !
    ! -- Print lake flows table
    if (ibudfl /= 0 .and. this%iprflow /= 0) then
      call this%budobj%write_flowtable(this%dis, kstp, kper)
    end if
  end subroutine lak_ot_package_flows

  !> @brief Write flows to binary file and/or print flows to budget
  !<
  subroutine lak_ot_model_flows(this, icbcfl, ibudfl, icbcun, imap)
    class(LakType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    integer(I4B), intent(in) :: icbcun
    integer(I4B), dimension(:), optional, intent(in) :: imap
    !
    ! -- write the flows from the budobj
    call this%BndType%bnd_ot_model_flows(icbcfl, ibudfl, icbcun, this%imap)
  end subroutine lak_ot_model_flows

  !> @brief Save LAK-calculated values to binary file
  !<
  subroutine lak_ot_dv(this, idvsave, idvprint)
    use TdisModule, only: kstp, kper, pertim, totim
    use ConstantsModule, only: DHNOFLO, DHDRY
    use InputOutputModule, only: ulasav
    class(LakType) :: this
    integer(I4B), intent(in) :: idvsave
    integer(I4B), intent(in) :: idvprint
    integer(I4B) :: ibinun
    integer(I4B) :: n
    real(DP) :: v
    real(DP) :: d
    real(DP) :: stage
    real(DP) :: sa
    real(DP) :: wa
    !
    ! -- set unit number for binary dependent variable output
    ibinun = 0
    if (this%istageout /= 0) then
      ibinun = this%istageout
    end if
    if (idvsave == 0) ibinun = 0
    !
    ! -- write lake binary output
    if (ibinun > 0) then
      do n = 1, this%nlakes
        v = this%xnewpak(n)
        d = v - this%lakebot(n)
        if (this%iboundpak(n) == 0) then
          v = DHNOFLO
        else if (d <= DZERO) then
          v = DHDRY
        end if
        this%dbuff(n) = v
      end do
      call ulasav(this%dbuff, '           STAGE', kstp, kper, pertim, totim, &
                  this%nlakes, 1, 1, ibinun)
    end if
    !
    ! -- Print lake stage table
    if (idvprint /= 0 .and. this%iprhed /= 0) then
      !
      ! -- set table kstp and kper
      call this%stagetab%set_kstpkper(kstp, kper)
      !
      ! -- write data
      do n = 1, this%nlakes
        if (this%iboundpak(n) == 0) then
          stage = DHNOFLO
          sa = DHNOFLO
          wa = DHNOFLO
          v = DHNOFLO
        else
          stage = this%xnewpak(n)
          call this%lak_calculate_sarea(n, stage, sa)
          call this%lak_calculate_warea(n, stage, wa)
          call this%lak_calculate_vol(n, stage, v)
        end if
        if (this%inamedbound == 1) then
          call this%stagetab%add_term(this%lakename(n))
        end if
        call this%stagetab%add_term(n)
        call this%stagetab%add_term(stage)
        call this%stagetab%add_term(sa)
        call this%stagetab%add_term(wa)
        call this%stagetab%add_term(v)
      end do
    end if
  end subroutine lak_ot_dv

  !> @brief Write LAK budget to listing file
  !<
  subroutine lak_ot_bdsummary(this, kstp, kper, iout, ibudfl)
    ! -- module
    use TdisModule, only: totim, delt
    ! -- dummy
    class(LakType) :: this !< LakType object
    integer(I4B), intent(in) :: kstp !< time step number
    integer(I4B), intent(in) :: kper !< period number
    integer(I4B), intent(in) :: iout !< flag and unit number for the model listing file
    integer(I4B), intent(in) :: ibudfl !< flag indicating budget should be written
    !
    call this%budobj%write_budtable(kstp, kper, iout, ibudfl, totim, delt)
  end subroutine lak_ot_bdsummary

  !> @brief Deallocate objects
  !<
  subroutine lak_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(LakType) :: this
    !
    ! -- arrays
    deallocate (this%lakename)
    deallocate (this%status)
    deallocate (this%clakbudget)
    call mem_deallocate(this%dbuff)
    deallocate (this%cauxcbc)
    call mem_deallocate(this%qauxcbc)
    call mem_deallocate(this%qleak)
    call mem_deallocate(this%qsto)
    call mem_deallocate(this%denseterms)
    call mem_deallocate(this%viscratios)
    !
    ! -- tables
    if (this%ntables > 0) then
      call mem_deallocate(this%ialaktab)
      call mem_deallocate(this%tabstage)
      call mem_deallocate(this%tabvolume)
      call mem_deallocate(this%tabsarea)
      call mem_deallocate(this%tabwarea)
    end if
    !
    ! -- budobj
    call this%budobj%budgetobject_da()
    deallocate (this%budobj)
    nullify (this%budobj)
    !
    ! -- outlets
    if (this%noutlets > 0) then
      call mem_deallocate(this%lakein)
      call mem_deallocate(this%lakeout)
      call mem_deallocate(this%iouttype)
      call mem_deallocate(this%outrate)
      call mem_deallocate(this%outinvert)
      call mem_deallocate(this%outwidth)
      call mem_deallocate(this%outrough)
      call mem_deallocate(this%outslope)
      call mem_deallocate(this%simoutrate)
    end if
    !
    ! -- stage table
    if (this%iprhed > 0) then
      call this%stagetab%table_da()
      deallocate (this%stagetab)
      nullify (this%stagetab)
    end if
    !
    ! -- package csv table
    if (this%ipakcsv > 0) then
      if (associated(this%pakcsvtab)) then
        call this%pakcsvtab%table_da()
        deallocate (this%pakcsvtab)
        nullify (this%pakcsvtab)
      end if
    end if
    !
    ! -- scalars
    call mem_deallocate(this%iprhed)
    call mem_deallocate(this%istageout)
    call mem_deallocate(this%ibudgetout)
    call mem_deallocate(this%ibudcsv)
    call mem_deallocate(this%ipakcsv)
    call mem_deallocate(this%nlakes)
    call mem_deallocate(this%noutlets)
    call mem_deallocate(this%ntables)
    call mem_deallocate(this%convlength)
    call mem_deallocate(this%convtime)
    call mem_deallocate(this%outdmax)
    call mem_deallocate(this%igwhcopt)
    call mem_deallocate(this%iconvchk)
    call mem_deallocate(this%iconvresidchk)
    call mem_deallocate(this%maxlakit)
    call mem_deallocate(this%surfdep)
    call mem_deallocate(this%dmaxchg)
    call mem_deallocate(this%delh)
    call mem_deallocate(this%pdmax)
    call mem_deallocate(this%check_attr)
    call mem_deallocate(this%bditems)
    call mem_deallocate(this%cbcauxitems)
    call mem_deallocate(this%idense)
    !
    call mem_deallocate(this%nlakeconn)
    call mem_deallocate(this%idxlakeconn)
    call mem_deallocate(this%ntabrow)
    call mem_deallocate(this%strt)
    call mem_deallocate(this%laketop)
    call mem_deallocate(this%lakebot)
    call mem_deallocate(this%sareamax)
    call mem_deallocate(this%stage)
    call mem_deallocate(this%rainfall)
    call mem_deallocate(this%evaporation)
    call mem_deallocate(this%runoff)
    call mem_deallocate(this%inflow)
    call mem_deallocate(this%withdrawal)
    call mem_deallocate(this%lauxvar)
    call mem_deallocate(this%avail)
    call mem_deallocate(this%lkgwsink)
    call mem_deallocate(this%ncncvr)
    call mem_deallocate(this%surfin)
    call mem_deallocate(this%surfout)
    call mem_deallocate(this%surfout1)
    call mem_deallocate(this%precip)
    call mem_deallocate(this%precip1)
    call mem_deallocate(this%evap)
    call mem_deallocate(this%evap1)
    call mem_deallocate(this%evapo)
    call mem_deallocate(this%withr)
    call mem_deallocate(this%withr1)
    call mem_deallocate(this%flwin)
    call mem_deallocate(this%flwiter)
    call mem_deallocate(this%flwiter1)
    call mem_deallocate(this%seep)
    call mem_deallocate(this%seep1)
    call mem_deallocate(this%seep0)
    call mem_deallocate(this%stageiter)
    call mem_deallocate(this%chterm)
    !
    ! -- lake boundary and stages
    call mem_deallocate(this%iboundpak)
    call mem_deallocate(this%xnewpak)
    call mem_deallocate(this%xoldpak)
    !
    ! -- lake iteration variables
    call mem_deallocate(this%iseepc)
    call mem_deallocate(this%idhc)
    call mem_deallocate(this%en1)
    call mem_deallocate(this%en2)
    call mem_deallocate(this%r1)
    call mem_deallocate(this%r2)
    call mem_deallocate(this%dh0)
    call mem_deallocate(this%s0)
    call mem_deallocate(this%qgwf0)
    !
    ! -- lake connection variables
    call mem_deallocate(this%imap)
    call mem_deallocate(this%cellid)
    call mem_deallocate(this%nodesontop)
    call mem_deallocate(this%ictype)
    call mem_deallocate(this%bedleak)
    call mem_deallocate(this%belev)
    call mem_deallocate(this%telev)
    call mem_deallocate(this%connlength)
    call mem_deallocate(this%connwidth)
    call mem_deallocate(this%sarea)
    call mem_deallocate(this%warea)
    call mem_deallocate(this%satcond)
    call mem_deallocate(this%simcond)
    call mem_deallocate(this%simlakgw)
    !
    ! -- pointers to gwf variables
    nullify (this%gwfiss)
    !
    ! -- Parent object
    call this%BndType%bnd_da()
  end subroutine lak_da

  !> @brief Define the list heading that is written to iout when PRINT_INPUT
  !! option is used
  !<
  subroutine define_listlabel(this)
    ! -- modules
    class(LakType), intent(inout) :: this
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
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'STRESS RATE'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
  end subroutine define_listlabel

  !> @brief Set pointers to model arrays and variables so that a package has
  !! access to these things
  !<
  subroutine lak_set_pointers(this, neq, ibound, xnew, xold, flowja)
    ! -- dummy
    class(LakType) :: this
    integer(I4B), pointer :: neq
    integer(I4B), dimension(:), pointer, contiguous :: ibound
    real(DP), dimension(:), pointer, contiguous :: xnew
    real(DP), dimension(:), pointer, contiguous :: xold
    real(DP), dimension(:), pointer, contiguous :: flowja
    !
    ! -- call base BndType set_pointers
    call this%BndType%set_pointers(neq, ibound, xnew, xold, flowja)
    !
    ! -- Set the LAK pointers
    !
    ! -- set package pointers
    !istart = this%dis%nodes + this%ioffset + 1
    !iend = istart + this%nlakes - 1
    !this%iboundpak => this%ibound(istart:iend)
    !this%xnewpak => this%xnew(istart:iend)
    !
    ! -- initialize xnewpak
    !do n = 1, this%nlakes
    !  this%xnewpak(n) = DEP20
    !end do
  end subroutine lak_set_pointers

  !> @brief Procedures related to observations (type-bound)
  !!
  !! Return true because LAK package supports observations. Overrides
  !! BndType%bnd_obs_supported()
  !<
  logical function lak_obs_supported(this)
    ! -- dummy
    class(LakType) :: this
    !
    lak_obs_supported = .true.
  end function lak_obs_supported

  !> @brief Store observation type supported by LAK package. Overrides
  !! BndType%bnd_df_obs
  !<
  subroutine lak_df_obs(this)
    ! -- dummy
    class(LakType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    ! -- Store obs type and assign procedure pointer
    !    for stage observation type.
    call this%obs%StoreObsType('stage', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for ext-inflow observation type.
    call this%obs%StoreObsType('ext-inflow', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for outlet-inflow observation type.
    call this%obs%StoreObsType('outlet-inflow', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for inflow observation type.
    call this%obs%StoreObsType('inflow', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for from-mvr observation type.
    call this%obs%StoreObsType('from-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for rainfall observation type.
    call this%obs%StoreObsType('rainfall', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for runoff observation type.
    call this%obs%StoreObsType('runoff', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for lak observation type.
    call this%obs%StoreObsType('lak', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for evaporation observation type.
    call this%obs%StoreObsType('evaporation', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for withdrawal observation type.
    call this%obs%StoreObsType('withdrawal', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for ext-outflow observation type.
    call this%obs%StoreObsType('ext-outflow', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for to-mvr observation type.
    call this%obs%StoreObsType('to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for storage observation type.
    call this%obs%StoreObsType('storage', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for constant observation type.
    call this%obs%StoreObsType('constant', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for outlet observation type.
    call this%obs%StoreObsType('outlet', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for volume observation type.
    call this%obs%StoreObsType('volume', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for surface-area observation type.
    call this%obs%StoreObsType('surface-area', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for wetted-area observation type.
    call this%obs%StoreObsType('wetted-area', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for conductance observation type.
    call this%obs%StoreObsType('conductance', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => lak_process_obsID
  end subroutine lak_df_obs

  !> @brief Calculate observations this time step and call ObsType%SaveOneSimval
  !! for each LakType observation.
  !<
  subroutine lak_bd_obs(this)
    ! -- dummy
    class(LakType) :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: igwfnode
    integer(I4B) :: j
    integer(I4B) :: jj
    integer(I4B) :: n
    real(DP) :: hgwf
    real(DP) :: hlak
    real(DP) :: v
    real(DP) :: v2
    type(ObserveType), pointer :: obsrv => null()
    !
    ! Write simulated values for all LAK observations
    if (this%obs%npakobs > 0) then
      call this%obs%obs_bd_clear()
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv
        do j = 1, obsrv%indxbnds_count
          v = DNODATA
          jj = obsrv%indxbnds(j)
          select case (obsrv%ObsTypeId)
          case ('STAGE')
            if (this%iboundpak(jj) /= 0) then
              v = this%xnewpak(jj)
            end if
          case ('EXT-INFLOW')
            if (this%iboundpak(jj) /= 0) then
              call this%lak_calculate_inflow(jj, v)
            end if
          case ('OUTLET-INFLOW')
            if (this%iboundpak(jj) /= 0) then
              call this%lak_calculate_outlet_inflow(jj, v)
            end if
          case ('INFLOW')
            if (this%iboundpak(jj) /= 0) then
              call this%lak_calculate_inflow(jj, v)
              call this%lak_calculate_outlet_inflow(jj, v2)
              v = v + v2
            end if
          case ('FROM-MVR')
            if (this%iboundpak(jj) /= 0) then
              if (this%imover == 1) then
                v = this%pakmvrobj%get_qfrommvr(jj)
              end if
            end if
          case ('RAINFALL')
            if (this%iboundpak(jj) /= 0) then
              v = this%precip(jj)
            end if
          case ('RUNOFF')
            if (this%iboundpak(jj) /= 0) then
              v = this%runoff(jj)
            end if
          case ('LAK')
            n = this%imap(jj)
            if (this%iboundpak(n) /= 0) then
              igwfnode = this%cellid(jj)
              hgwf = this%xnew(igwfnode)
              if (this%hcof(jj) /= DZERO) then
                v = -(this%hcof(jj) * (this%xnewpak(n) - hgwf))
              else
                v = -this%rhs(jj)
              end if
            end if
          case ('EVAPORATION')
            if (this%iboundpak(jj) /= 0) then
              v = this%evap(jj)
            end if
          case ('WITHDRAWAL')
            if (this%iboundpak(jj) /= 0) then
              v = this%withr(jj)
            end if
          case ('EXT-OUTFLOW')
            n = this%lakein(jj)
            if (this%iboundpak(n) /= 0) then
              if (this%lakeout(jj) == 0) then
                v = this%simoutrate(jj)
                if (v < DZERO) then
                  if (this%imover == 1) then
                    v = v + this%pakmvrobj%get_qtomvr(jj)
                  end if
                end if
              end if
            end if
          case ('TO-MVR')
            n = this%lakein(jj)
            if (this%iboundpak(n) /= 0) then
              if (this%imover == 1) then
                v = this%pakmvrobj%get_qtomvr(jj)
                if (v > DZERO) then
                  v = -v
                end if
              end if
            end if
          case ('STORAGE')
            if (this%iboundpak(jj) /= 0) then
              v = this%qsto(jj)
            end if
          case ('CONSTANT')
            if (this%iboundpak(jj) /= 0) then
              v = this%chterm(jj)
            end if
          case ('OUTLET')
            n = this%lakein(jj)
            if (this%iboundpak(n) /= 0) then
              v = this%simoutrate(jj)
            end if
          case ('VOLUME')
            if (this%iboundpak(jj) /= 0) then
              call this%lak_calculate_vol(jj, this%xnewpak(jj), v)
            end if
          case ('SURFACE-AREA')
            if (this%iboundpak(jj) /= 0) then
              hlak = this%xnewpak(jj)
              call this%lak_calculate_sarea(jj, hlak, v)
            end if
          case ('WETTED-AREA')
            n = this%imap(jj)
            if (this%iboundpak(n) /= 0) then
              hlak = this%xnewpak(n)
              igwfnode = this%cellid(jj)
              hgwf = this%xnew(igwfnode)
              call this%lak_calculate_conn_warea(n, jj, hlak, hgwf, v)
            end if
          case ('CONDUCTANCE')
            n = this%imap(jj)
            if (this%iboundpak(n) /= 0) then
              hlak = this%xnewpak(n)
              igwfnode = this%cellid(jj)
              hgwf = this%xnew(igwfnode)
              call this%lak_calculate_conn_conductance(n, jj, hlak, hgwf, v)
            end if
          case default
            errmsg = 'Unrecognized observation type: '//trim(obsrv%ObsTypeId)
            call store_error(errmsg)
          end select
          call this%obs%SaveOneSimval(obsrv, v)
        end do
      end do
      !
      ! -- write summary of error messages
      if (count_errors() > 0) then
        call store_error_unit(this%inunit)
      end if
    end if
  end subroutine lak_bd_obs

  !> @brief Process each observation
  !!
  !! Only done the first stress period since boundaries are fixed for the
  !! simulation
  !<
  subroutine lak_rp_obs(this)
    use TdisModule, only: kper
    ! -- dummy
    class(LakType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: nn1
    integer(I4B) :: nn2
    integer(I4B) :: jj
    character(len=LENBOUNDNAME) :: bname
    logical(LGP) :: jfound
    class(ObserveType), pointer :: obsrv => null()
    ! -- formats
10  format('Boundary "', a, '" for observation "', a, &
           '" is invalid in package "', a, '"')
    !
    ! -- process each package observation
    !    only done the first stress period since boundaries are fixed
    !    for the simulation
    if (kper == 1) then
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv
        !
        ! -- get node number 1
        nn1 = obsrv%NodeNumber
        if (nn1 == NAMEDBOUNDFLAG) then
          bname = obsrv%FeatureName
          if (bname /= '') then
            ! -- Observation lake is based on a boundary name.
            !    Iterate through all lakes to identify and store
            !    corresponding index in bound array.
            jfound = .false.
            if (obsrv%ObsTypeId == 'LAK' .or. &
                obsrv%ObsTypeId == 'CONDUCTANCE' .or. &
                obsrv%ObsTypeId == 'WETTED-AREA') then
              do j = 1, this%nlakes
                do jj = this%idxlakeconn(j), this%idxlakeconn(j + 1) - 1
                  if (this%boundname(jj) == bname) then
                    jfound = .true.
                    call obsrv%AddObsIndex(jj)
                  end if
                end do
              end do
            else if (obsrv%ObsTypeId == 'EXT-OUTFLOW' .or. &
                     obsrv%ObsTypeId == 'TO-MVR' .or. &
                     obsrv%ObsTypeId == 'OUTLET') then
              do j = 1, this%noutlets
                jj = this%lakein(j)
                if (this%lakename(jj) == bname) then
                  jfound = .true.
                  call obsrv%AddObsIndex(j)
                end if
              end do
            else
              do j = 1, this%nlakes
                if (this%lakename(j) == bname) then
                  jfound = .true.
                  call obsrv%AddObsIndex(j)
                end if
              end do
            end if
            if (.not. jfound) then
              write (errmsg, 10) &
                trim(bname), trim(obsrv%Name), trim(this%packName)
              call store_error(errmsg)
            end if
          end if
        else
          if (obsrv%indxbnds_count == 0) then
            if (obsrv%ObsTypeId == 'LAK' .or. &
                obsrv%ObsTypeId == 'CONDUCTANCE' .or. &
                obsrv%ObsTypeId == 'WETTED-AREA') then
              nn2 = obsrv%NodeNumber2
              j = this%idxlakeconn(nn1) + nn2 - 1
              call obsrv%AddObsIndex(j)
            else
              call obsrv%AddObsIndex(nn1)
            end if
          else
            errmsg = 'Programming error in lak_rp_obs'
            call store_error(errmsg)
          end if
        end if
        !
        ! -- catch non-cumulative observation assigned to observation defined
        !    by a boundname that is assigned to more than one element
        if (obsrv%ObsTypeId == 'STAGE') then
          if (obsrv%indxbnds_count > 1) then
            write (errmsg, '(a,3(1x,a))') &
              trim(adjustl(obsrv%ObsTypeId)), &
              'for observation', trim(adjustl(obsrv%Name)), &
              ' must be assigned to a lake with a unique boundname.'
            call store_error(errmsg)
          end if
        end if
        !
        ! -- check that index values are valid
        if (obsrv%ObsTypeId == 'TO-MVR' .or. &
            obsrv%ObsTypeId == 'EXT-OUTFLOW' .or. &
            obsrv%ObsTypeId == 'OUTLET') then
          do j = 1, obsrv%indxbnds_count
            nn1 = obsrv%indxbnds(j)
            if (nn1 < 1 .or. nn1 > this%noutlets) then
              write (errmsg, '(a,1x,a,1x,i0,1x,a,1x,i0,a)') &
                trim(adjustl(obsrv%ObsTypeId)), &
                ' outlet must be > 0 and <=', this%noutlets, &
                '(specified value is ', nn1, ')'
              call store_error(errmsg)
            end if
          end do
        else if (obsrv%ObsTypeId == 'LAK' .or. &
                 obsrv%ObsTypeId == 'CONDUCTANCE' .or. &
                 obsrv%ObsTypeId == 'WETTED-AREA') then
          do j = 1, obsrv%indxbnds_count
            nn1 = obsrv%indxbnds(j)
            if (nn1 < 1 .or. nn1 > this%maxbound) then
              write (errmsg, '(a,1x,a,1x,i0,1x,a,1x,i0,a)') &
                trim(adjustl(obsrv%ObsTypeId)), &
                'lake connection number must be > 0 and <=', this%maxbound, &
                '(specified value is ', nn1, ')'
              call store_error(errmsg)
            end if
          end do
        else
          do j = 1, obsrv%indxbnds_count
            nn1 = obsrv%indxbnds(j)
            if (nn1 < 1 .or. nn1 > this%nlakes) then
              write (errmsg, '(a,1x,a,1x,i0,1x,a,1x,i0,a)') &
                trim(adjustl(obsrv%ObsTypeId)), &
                ' lake must be > 0 and <=', this%nlakes, &
                '(specified value is ', nn1, ')'
              call store_error(errmsg)
            end if
          end do
        end if
      end do
      !
      ! -- evaluate if there are any observation errors
      if (count_errors() > 0) then
        call store_error_unit(this%inunit)
      end if
    end if
  end subroutine lak_rp_obs

  !
  ! -- Procedures related to observations (NOT type-bound)

  !> @brief This procedure is pointed to by ObsDataType%ProcesssIdPtr. It
  !! processes the ID string of an observation definition for LAK package
  !! observations.
  !<
  subroutine lak_process_obsID(obsrv, dis, inunitobs, iout)
    ! -- dummy
    type(ObserveType), intent(inout) :: obsrv
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: inunitobs
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B) :: nn1, nn2
    integer(I4B) :: icol, istart, istop
    character(len=LINELENGTH) :: string
    character(len=LENBOUNDNAME) :: bndname
    !
    string = obsrv%IDstring
    ! -- Extract lake number from string and store it.
    !    If 1st item is not an integer(I4B), it should be a
    !    lake name--deal with it.
    icol = 1
    ! -- get lake number or boundary name
    call extract_idnum_or_bndname(string, icol, istart, istop, nn1, bndname)
    if (nn1 == NAMEDBOUNDFLAG) then
      obsrv%FeatureName = bndname
    else
      if (obsrv%ObsTypeId == 'LAK' .or. obsrv%ObsTypeId == 'CONDUCTANCE' .or. &
          obsrv%ObsTypeId == 'WETTED-AREA') then
        call extract_idnum_or_bndname(string, icol, istart, istop, nn2, bndname)
        if (len_trim(bndName) < 1 .and. nn2 < 0) then
          write (errmsg, '(a,1x,a,a,1x,a,1x,a)') &
            'For observation type', trim(adjustl(obsrv%ObsTypeId)), &
            ', ID given as an integer and not as boundname,', &
            'but ID2 (iconn) is missing.  Either change ID to valid', &
            'boundname or supply valid entry for ID2.'
          call store_error(errmsg)
        end if
        if (nn2 == NAMEDBOUNDFLAG) then
          obsrv%FeatureName = bndname
          ! -- reset nn1
          nn1 = nn2
        else
          obsrv%NodeNumber2 = nn2
        end if
      end if
    end if
    ! -- store lake number (NodeNumber)
    obsrv%NodeNumber = nn1
  end subroutine lak_process_obsID

  !
  ! -- private LAK methods
  !

  !> @brief Accumulate constant head terms for budget
  !<
  subroutine lak_accumulate_chterm(this, ilak, rrate, chratin, chratout)
    ! -- dummy
    class(LakType) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(in) :: rrate
    real(DP), intent(inout) :: chratin
    real(DP), intent(inout) :: chratout
    ! -- locals
    real(DP) :: q
    !
    ! code
    if (this%iboundpak(ilak) < 0) then
      q = -rrate
      this%chterm(ilak) = this%chterm(ilak) + q
      !
      ! -- See if flow is into lake or out of lake.
      if (q < DZERO) then
        !
        ! -- Flow is out of lake subtract rate from ratout.
        chratout = chratout - q
      else
        !
        ! -- Flow is into lake; add rate to ratin.
        chratin = chratin + q
      end if
    end if
  end subroutine lak_accumulate_chterm

  !> @brief Store the lake head and connection conductance in the bound array
  !<
  subroutine lak_bound_update(this)
    ! -- dummy
    class(LakType), intent(inout) :: this
    ! -- local
    integer(I4B) :: j, n, node
    real(DP) :: hlak, head, clak
    !
    ! -- Return if no lak lakes
    if (this%nbound == 0) return
    !
    ! -- Calculate hcof and rhs for each lak entry
    do n = 1, this%nlakes
      hlak = this%xnewpak(n)
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        node = this%cellid(j)
        head = this%xnew(node)
        call this%lak_calculate_conn_conductance(n, j, hlak, head, clak)
        this%bound(1, j) = hlak
        this%bound(2, j) = clak
      end do
    end do
  end subroutine lak_bound_update

  !> @brief Solve for lake stage
  !<
  subroutine lak_solve(this, update)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(LakType), intent(inout) :: this
    logical(LGP), intent(in), optional :: update
    ! -- local
    logical(LGP) :: lupdate
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: iicnvg
    integer(I4B) :: iter
    integer(I4B) :: maxiter
    integer(I4B) :: ncnv
    integer(I4B) :: idry
    integer(I4B) :: idry1
    integer(I4B) :: igwfnode
    integer(I4B) :: ibflg
    integer(I4B) :: idhp
    real(DP) :: hlak
    real(DP) :: hlak0
    real(DP) :: v0
    real(DP) :: v1
    real(DP) :: head
    real(DP) :: ra
    real(DP) :: ro
    real(DP) :: qinf
    real(DP) :: ex
    real(DP) :: ev
    real(DP) :: outinf
    real(DP) :: qlakgw
    real(DP) :: qlakgw1
    real(DP) :: gwfhcof
    real(DP) :: gwfrhs
    real(DP) :: avail
    real(DP) :: resid
    real(DP) :: resid1
    real(DP) :: residb
    real(DP) :: wr
    real(DP) :: derv
    real(DP) :: dh
    real(DP) :: adh
    real(DP) :: adh0
    real(DP) :: delh
    real(DP) :: ts
    real(DP) :: area
    real(DP) :: qtolfact
    !
    ! -- set lupdate
    if (present(update)) then
      lupdate = update
    else
      lupdate = .true.
    end if
    !
    ! -- initialize
    avail = DZERO
    delh = this%delh
    !
    ! -- initialize
    do n = 1, this%nlakes
      this%ncncvr(n) = 0
      this%surfin(n) = DZERO
      this%surfout(n) = DZERO
      this%surfout1(n) = DZERO
      if (this%xnewpak(n) < this%lakebot(n)) then
        this%xnewpak(n) = this%lakebot(n)
      end if
      if (this%gwfiss /= 0) then
        this%xoldpak(n) = this%xnewpak(n)
      end if
      ! -- lake iteration items
      this%iseepc(n) = 0
      this%idhc(n) = 0
      this%en1(n) = this%lakebot(n)
      call this%lak_calculate_residual(n, this%en1(n), this%r1(n))
      this%en2(n) = this%laketop(n)
      call this%lak_calculate_residual(n, this%en2(n), this%r2(n))
    end do
    do n = 1, this%noutlets
      this%simoutrate(n) = DZERO
    end do
    !
    ! -- sum up inflows from mover inflows
    do n = 1, this%nlakes
      call this%lak_calculate_outlet_inflow(n, this%surfin(n))
    end do
    !
    ! -- sum up overland runoff, inflows, and external flows into lake
    !    (includes maximum lake volume)
    do n = 1, this%nlakes
      hlak0 = this%xoldpak(n)
      hlak = this%xnewpak(n)
      call this%lak_calculate_runoff(n, ro)
      call this%lak_calculate_inflow(n, qinf)
      call this%lak_calculate_external(n, ex)
      call this%lak_calculate_vol(n, hlak0, v0)
      call this%lak_calculate_vol(n, hlak, v1)
      this%flwin(n) = this%surfin(n) + ro + qinf + ex + &
                      max(v0, v1) / delt
    end do
    !
    ! -- sum up inflows from upstream outlets
    do n = 1, this%nlakes
      call this%lak_calculate_outlet_inflow(n, outinf)
      this%flwin(n) = this%flwin(n) + outinf
    end do
    !
    iicnvg = 0
    maxiter = this%maxlakit
    !
    ! -- outer loop
    converge: do iter = 1, maxiter
      ncnv = 0
      do n = 1, this%nlakes
        if (this%ncncvr(n) == 0) ncnv = 1
      end do
      if (iter == maxiter) ncnv = 0
      if (ncnv == 0) iicnvg = 1
      !
      ! -- initialize variables
      do n = 1, this%nlakes
        this%evap(n) = DZERO
        this%precip(n) = DZERO
        this%precip1(n) = DZERO
        this%seep(n) = DZERO
        this%seep1(n) = DZERO
        this%evap(n) = DZERO
        this%evap1(n) = DZERO
        this%evapo(n) = DZERO
        this%withr(n) = DZERO
        this%withr1(n) = DZERO
        this%flwiter(n) = this%flwin(n)
        this%flwiter1(n) = this%flwin(n)
        if (this%gwfiss /= 0) then
          this%flwiter(n) = DEP20
          this%flwiter1(n) = DEP20
        end if
        do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
          this%hcof(j) = DZERO
          this%rhs(j) = DZERO
        end do
      end do
      !
      estseep: do i = 1, 2
        lakseep: do n = 1, this%nlakes
          ! -- skip inactive lakes
          if (this%iboundpak(n) == 0) then
            cycle lakseep
          end if
          ! - set xoldpak to xnewpak if steady-state
          if (this%gwfiss /= 0) then
            this%xoldpak(n) = this%xnewpak(n)
          end if
          hlak = this%xnewpak(n)
          calcconnseep: do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
            igwfnode = this%cellid(j)
            head = this%xnew(igwfnode)
            if (this%ncncvr(n) /= 2) then
              if (this%ibound(igwfnode) > 0) then
                call this%lak_estimate_conn_exchange(i, n, j, idry, hlak, &
                                                     head, qlakgw, &
                                                     this%flwiter(n), &
                                                     gwfhcof, gwfrhs)
                call this%lak_estimate_conn_exchange(i, n, j, idry1, &
                                                     hlak + delh, head, qlakgw1, &
                                                     this%flwiter1(n))
                !
                ! -- add to gwf matrix
                if (ncnv == 0 .and. i == 2) then
                  if (j == this%maxbound) then
                    this%ncncvr(n) = 2
                  end if
                  if (idry /= 1) then
                    this%hcof(j) = gwfhcof
                    this%rhs(j) = gwfrhs
                  else
                    this%hcof(j) = DZERO
                    this%rhs(j) = qlakgw
                  end if
                end if
                if (i == 2) then
                  this%seep(n) = this%seep(n) + qlakgw
                  this%seep1(n) = this%seep1(n) + qlakgw1
                end if
              end if
            end if
            !
          end do calcconnseep
        end do lakseep
      end do estseep
      !
      laklevel: do n = 1, this%nlakes
        ! -- skip inactive lakes
        if (this%iboundpak(n) == 0) then
          this%ncncvr(n) = 1
          cycle laklevel
        end if
        ibflg = 0
        hlak = this%xnewpak(n)
        if (iter < maxiter) then
          this%stageiter(n) = this%xnewpak(n)
        end if
        call this%lak_calculate_rainfall(n, hlak, ra)
        this%precip(n) = ra
        this%flwiter(n) = this%flwiter(n) + ra
        call this%lak_calculate_rainfall(n, hlak + delh, ra)
        this%precip1(n) = ra
        this%flwiter1(n) = this%flwiter1(n) + ra
        !
        ! -- limit withdrawals to lake inflows and lake storage
        call this%lak_calculate_withdrawal(n, this%flwiter(n), wr)
        this%withr(n) = wr
        call this%lak_calculate_withdrawal(n, this%flwiter1(n), wr)
        this%withr1(n) = wr
        !
        ! -- limit evaporation to lake inflows and lake storage
        call this%lak_calculate_evaporation(n, hlak, this%flwiter(n), ev)
        this%evap(n) = ev
        call this%lak_calculate_evaporation(n, hlak + delh, this%flwiter1(n), ev)
        this%evap1(n) = ev
        !
        ! -- no outlet flow if evaporation consumes all water
        call this%lak_calculate_outlet_outflow(n, hlak + delh, &
                                               this%flwiter1(n), &
                                               this%surfout1(n))
        call this%lak_calculate_outlet_outflow(n, hlak, this%flwiter(n), &
                                               this%surfout(n))
        !
        ! -- update the surface inflow values
        call this%lak_calculate_outlet_inflow(n, this%surfin(n))
        !
        !
        if (ncnv == 1) then
          if (this%iboundpak(n) > 0 .and. lupdate .eqv. .true.) then
            !
            ! -- recalculate flwin
            hlak0 = this%xoldpak(n)
            hlak = this%xnewpak(n)
            call this%lak_calculate_vol(n, hlak0, v0)
            call this%lak_calculate_vol(n, hlak, v1)
            call this%lak_calculate_runoff(n, ro)
            call this%lak_calculate_inflow(n, qinf)
            call this%lak_calculate_external(n, ex)
            this%flwin(n) = this%surfin(n) + ro + qinf + ex + &
                            max(v0, v1) / delt
            !
            ! -- compute new lake stage using Newton's method
            resid = this%precip(n) + this%evap(n) + this%withr(n) + ro + &
                    qinf + ex + this%surfin(n) + &
                    this%surfout(n) + this%seep(n)
            resid1 = this%precip1(n) + this%evap1(n) + this%withr1(n) + ro + &
                     qinf + ex + this%surfin(n) + &
                     this%surfout1(n) + this%seep1(n)
            !
            ! -- add storage changes for transient stress periods
            hlak = this%xnewpak(n)
            if (this%gwfiss /= 1) then
              call this%lak_calculate_vol(n, hlak, v1)
              resid = resid + (v0 - v1) / delt
              call this%lak_calculate_vol(n, hlak + delh, v1)
              resid1 = resid1 + (v0 - v1) / delt
            end if
            !
            ! -- determine the derivative and the stage change
            if (ABS(resid1 - resid) > DZERO) then
              derv = (resid1 - resid) / delh
              dh = DZERO
              if (ABS(derv) > DPREC) then
                dh = resid / derv
              end if
            else
              if (resid < DZERO) then
                resid = DZERO
              end if
              call this%lak_vol2stage(n, resid, dh)
              dh = hlak - dh
              this%ncncvr(n) = 1
            end if
            !
            ! -- determine if the updated stage is outside the endpoints
            ts = hlak - dh
            if (iter == 1) this%dh0(n) = dh
            adh = ABS(dh)
            adh0 = ABS(this%dh0(n))
            if ((ts >= this%en2(n)) .or. (ts < this%en1(n))) then
              ! -- use bisection if dh is increasing or updated stage is below the
              !    bottom of the lake
              if ((adh > adh0) .or. (ts - this%lakebot(n)) < DPREC) then
                residb = resid
                call this%lak_bisection(n, ibflg, hlak, ts, dh, residb)
              end if
            end if
            !
            ! -- set seep0 on the first lake iteration
            if (iter == 1) then
              this%seep0(n) = this%seep(n)
            end if
            !
            ! -- check for slow convergence
            if (this%seep(n) * this%seep0(n) < DPREC) then
              this%iseepc(n) = this%iseepc(n) + 1
            else
              this%iseepc(n) = 0
            end if
            ! -- determine of convergence is slow and oscillating
            idhp = 0
            if (dh * this%dh0(n) < DPREC) idhp = 1
            ! -- determine if stage change is increasing
            adh = ABS(dh)
            if (adh > adh0) idhp = 1
            ! -- increment idhc convergence flag
            if (idhp == 1) then
              this%idhc(n) = this%idhc(n) + 1
            end if
            !
            ! -- switch to bisection when the Newton-Raphson method oscillates
            !    or when convergence is slow
            if (ibflg == 1) then
              if (this%iseepc(n) > 7 .or. this%idhc(n) > 12) then
                call this%lak_bisection(n, ibflg, hlak, ts, dh, residb)
              end if
            end if
          else
            dh = DZERO
          end if
          !
          ! -- update lake stage
          hlak = hlak - dh
          if (hlak < this%lakebot(n)) then
            hlak = this%lakebot(n)
          end if
          !
          ! -- calculate surface area
          call this%lak_calculate_sarea(n, hlak, area)
          !
          ! -- set the Q to length factor
          if (area > DZERO) then
            qtolfact = delt / area
          else
            qtolfact = DZERO
          end if
          !
          ! -- recalculate the residual
          call this%lak_calculate_residual(n, hlak, resid)
          !
          ! -- evaluate convergence
          !if (ABS(dh) < delh) then
          if (ABS(dh) < delh .and. abs(resid) * qtolfact < this%dmaxchg) then
            this%ncncvr(n) = 1
          end if
          this%xnewpak(n) = hlak
          !
          ! -- save iterates for lake
          this%seep0(n) = this%seep(n)
          this%dh0(n) = dh
        end if
      end do laklevel
      !
      if (iicnvg == 1) exit converge
      !
    end do converge
    !
    ! -- Mover terms: store outflow after diversion loss
    !    as qformvr and reduce outflow (qd)
    !    by how much was actually sent to the mover
    if (this%imover == 1) then
      do n = 1, this%noutlets
        call this%pakmvrobj%accumulate_qformvr(n, -this%simoutrate(n))
      end do
    end if
  end subroutine lak_solve

  !> @ brief Lake package bisection method
  !!
  !!  Use bisection method to find lake stage that reduces the residual
  !<
  subroutine lak_bisection(this, n, ibflg, hlak, temporary_stage, dh, residual)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: n !< lake number
    integer(I4B), intent(inout) :: ibflg !< bisection flag
    real(DP), intent(in) :: hlak !< lake stage
    real(DP), intent(inout) :: temporary_stage !< temporary lake stage
    real(DP), intent(inout) :: dh !< lake stage change
    real(DP), intent(inout) :: residual !< lake residual
    ! -- local
    integer(I4B) :: i
    real(DP) :: temporary_stage0
    real(DP) :: residuala
    real(DP) :: endpoint1
    real(DP) :: endpoint2
    ! -- code
    ibflg = 1
    temporary_stage0 = hlak
    endpoint1 = this%en1(n)
    endpoint2 = this%en2(n)
    call this%lak_calculate_residual(n, temporary_stage, residuala)
    if (hlak > endpoint1 .and. hlak < endpoint2) then
      endpoint2 = hlak
    end if
    do i = 1, this%maxlakit
      temporary_stage = DHALF * (endpoint1 + endpoint2)
      call this%lak_calculate_residual(n, temporary_stage, residual)
      if (abs(residual) == DZERO .or. &
          abs(temporary_stage0 - temporary_stage) < this%dmaxchg) then
        exit
      end if
      call this%lak_calculate_residual(n, endpoint1, residuala)
      ! -- change end points
      ! -- root is between temporary_stage and endpoint2
      if (sign(DONE, residuala) == SIGN(DONE, residual)) then
        endpoint1 = temporary_stage
        ! -- root is between endpoint1 and temporary_stage
      else
        endpoint2 = temporary_stage
      end if
      temporary_stage0 = temporary_stage
    end do
    dh = hlak - temporary_stage
  end subroutine lak_bisection

  !> @brief Calculate the available volumetric rate for a lake given a passed
  !! stage
  !<
  subroutine lak_calculate_available(this, n, hlak, avail, &
                                     ra, ro, qinf, ex, headp)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    real(DP), intent(in) :: hlak
    real(DP), intent(inout) :: avail
    real(DP), intent(inout) :: ra
    real(DP), intent(inout) :: ro
    real(DP), intent(inout) :: qinf
    real(DP), intent(inout) :: ex
    real(DP), intent(in), optional :: headp
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: idry
    integer(I4B) :: igwfnode
    real(DP) :: hp
    real(DP) :: head
    real(DP) :: qlakgw
    real(DP) :: v0
    !
    ! -- set hp
    if (present(headp)) then
      hp = headp
    else
      hp = DZERO
    end if
    !
    ! -- initialize
    avail = DZERO
    !
    ! -- calculate the aquifer sources to the lake
    do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
      igwfnode = this%cellid(j)
      if (this%ibound(igwfnode) == 0) cycle
      head = this%xnew(igwfnode) + hp
      call this%lak_estimate_conn_exchange(1, n, j, idry, hlak, head, qlakgw, &
                                           avail)
    end do
    !
    ! -- add rainfall
    call this%lak_calculate_rainfall(n, hlak, ra)
    avail = avail + ra
    !
    ! -- calculate runoff
    call this%lak_calculate_runoff(n, ro)
    avail = avail + ro
    !
    ! -- calculate inflow
    call this%lak_calculate_inflow(n, qinf)
    avail = avail + qinf
    !
    ! -- calculate external flow terms
    call this%lak_calculate_external(n, ex)
    avail = avail + ex
    !
    ! -- calculate volume available in storage
    call this%lak_calculate_vol(n, this%xoldpak(n), v0)
    avail = avail + v0 / delt
  end subroutine lak_calculate_available

  !> @brief Calculate the residual for a lake given a passed stage
  !<
  subroutine lak_calculate_residual(this, n, hlak, resid, headp)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    real(DP), intent(in) :: hlak
    real(DP), intent(inout) :: resid
    real(DP), intent(in), optional :: headp
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: idry
    integer(I4B) :: igwfnode
    real(DP) :: hp
    real(DP) :: avail
    real(DP) :: head
    real(DP) :: ra
    real(DP) :: ro
    real(DP) :: qinf
    real(DP) :: ex
    real(DP) :: ev
    real(DP) :: wr
    real(DP) :: sout
    real(DP) :: sin
    real(DP) :: qlakgw
    real(DP) :: seep
    real(DP) :: hlak0
    real(DP) :: v0
    real(DP) :: v1
    !
    ! -- set hp
    if (present(headp)) then
      hp = headp
    else
      hp = DZERO
    end if
    !
    ! -- initialize
    resid = DZERO
    avail = DZERO
    seep = DZERO
    !
    ! -- calculate the available water
    call this%lak_calculate_available(n, hlak, avail, &
                                      ra, ro, qinf, ex, hp)
    !
    ! -- calculate groundwater seepage
    do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
      igwfnode = this%cellid(j)
      if (this%ibound(igwfnode) == 0) cycle
      head = this%xnew(igwfnode) + hp
      call this%lak_estimate_conn_exchange(2, n, j, idry, hlak, head, qlakgw, &
                                           avail)
      seep = seep + qlakgw
    end do
    !
    ! -- limit withdrawals to lake inflows and lake storage
    call this%lak_calculate_withdrawal(n, avail, wr)
    !
    ! -- limit evaporation to lake inflows and lake storage
    call this%lak_calculate_evaporation(n, hlak, avail, ev)
    !
    ! -- no outlet flow if evaporation consumes all water
    call this%lak_calculate_outlet_outflow(n, hlak, avail, sout)
    !
    ! -- update the surface inflow values
    call this%lak_calculate_outlet_inflow(n, sin)
    !
    ! -- calculate residual
    resid = ra + ev + wr + ro + qinf + ex + sin + sout + seep
    !
    ! -- include storage
    if (this%gwfiss /= 1) then
      hlak0 = this%xoldpak(n)
      call this%lak_calculate_vol(n, hlak0, v0)
      call this%lak_calculate_vol(n, hlak, v1)
      resid = resid + (v0 - v1) / delt
    end if
  end subroutine lak_calculate_residual

  !> @brief Set up the budget object that stores all the lake flows
  !<
  subroutine lak_setup_budobj(this)
    ! -- modules
    use ConstantsModule, only: LENBUDTXT
    ! -- dummy
    class(LakType) :: this
    ! -- local
    integer(I4B) :: nbudterm
    integer(I4B) :: nlen
    integer(I4B) :: j, n, n1, n2
    integer(I4B) :: maxlist, naux
    integer(I4B) :: idx
    real(DP) :: q
    character(len=LENBUDTXT) :: text
    character(len=LENBUDTXT), dimension(1) :: auxtxt
    !
    ! -- Determine the number of lake budget terms. These are fixed for
    !    the simulation and cannot change
    nbudterm = 9
    nlen = 0
    do n = 1, this%noutlets
      if (this%lakein(n) > 0 .and. this%lakeout(n) > 0) then
        nlen = nlen + 1
      end if
    end do
    if (nlen > 0) nbudterm = nbudterm + 1
    if (this%imover == 1) nbudterm = nbudterm + 2
    if (this%naux > 0) nbudterm = nbudterm + 1
    !
    ! -- set up budobj
    call budgetobject_cr(this%budobj, this%packName)
    call this%budobj%budgetobject_df(this%nlakes, nbudterm, 0, 0, &
                                     ibudcsv=this%ibudcsv)
    idx = 0
    !
    ! -- Go through and set up each budget term. nlen is the number
    !    of outlets that discharge into another lake
    if (nlen > 0) then
      text = '    FLOW-JA-FACE'
      idx = idx + 1
      maxlist = 2 * nlen
      naux = 0
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux, ordered_id1=.false.)
      !
      ! -- store connectivity
      call this%budobj%budterm(idx)%reset(2 * nlen)
      q = DZERO
      do n = 1, this%noutlets
        n1 = this%lakein(n)
        n2 = this%lakeout(n)
        if (n1 > 0 .and. n2 > 0) then
          call this%budobj%budterm(idx)%update_term(n1, n2, q)
          call this%budobj%budterm(idx)%update_term(n2, n1, -q)
        end if
      end do
    end if
    !
    ! --
    text = '             GWF'
    idx = idx + 1
    maxlist = this%maxbound
    naux = 1
    auxtxt(1) = '       FLOW-AREA'
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%name_model, &
                                             maxlist, .false., .true., &
                                             naux, auxtxt)
    call this%budobj%budterm(idx)%reset(this%maxbound)
    q = DZERO
    do n = 1, this%nlakes
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        n2 = this%cellid(j)
        call this%budobj%budterm(idx)%update_term(n, n2, q)
      end do
    end do
    !
    ! --
    text = '        RAINFALL'
    idx = idx + 1
    maxlist = this%nlakes
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! --
    text = '     EVAPORATION'
    idx = idx + 1
    maxlist = this%nlakes
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! --
    text = '          RUNOFF'
    idx = idx + 1
    maxlist = this%nlakes
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! --
    text = '      EXT-INFLOW'
    idx = idx + 1
    maxlist = this%nlakes
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! --
    text = '      WITHDRAWAL'
    idx = idx + 1
    maxlist = this%nlakes
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! --
    text = '     EXT-OUTFLOW'
    idx = idx + 1
    maxlist = this%nlakes
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! --
    text = '         STORAGE'
    idx = idx + 1
    maxlist = this%nlakes
    naux = 1
    auxtxt(1) = '          VOLUME'
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux, auxtxt)
    !
    ! --
    text = '        CONSTANT'
    idx = idx + 1
    maxlist = this%nlakes
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! --
    if (this%imover == 1) then
      !
      ! --
      text = '        FROM-MVR'
      idx = idx + 1
      maxlist = this%nlakes
      naux = 0
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux)
      !
      ! --
      text = '          TO-MVR'
      idx = idx + 1
      maxlist = this%noutlets
      naux = 0
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux, ordered_id1=.false.)
      !
      ! -- store to-mvr connection information
      call this%budobj%budterm(idx)%reset(this%noutlets)
      q = DZERO
      do n = 1, this%noutlets
        n1 = this%lakein(n)
        call this%budobj%budterm(idx)%update_term(n1, n1, q)
      end do
    end if
    !
    ! --
    naux = this%naux
    if (naux > 0) then
      !
      ! --
      text = '       AUXILIARY'
      idx = idx + 1
      maxlist = this%nlakes
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux, this%auxname)
    end if
    !
    ! -- if lake flow for each reach are written to the listing file
    if (this%iprflow /= 0) then
      call this%budobj%flowtable_df(this%iout)
    end if
  end subroutine lak_setup_budobj

  !> @brief Copy flow terms into this%budobj
  !<
  subroutine lak_fill_budobj(this)
    ! -- dummy
    class(LakType) :: this
    ! -- local
    integer(I4B) :: naux
    real(DP), dimension(:), allocatable :: auxvartmp
    !integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: n1
    integer(I4B) :: n2
    integer(I4B) :: ii
    integer(I4B) :: jj
    integer(I4B) :: idx
    integer(I4B) :: nlen
    real(DP) :: v, v1
    real(DP) :: q
    real(DP) :: lkstg, gwhead, wa
    !
    ! -- initialize counter
    idx = 0

    ! -- FLOW JA FACE
    nlen = 0
    do n = 1, this%noutlets
      if (this%lakein(n) > 0 .and. this%lakeout(n) > 0) then
        nlen = nlen + 1
      end if
    end do
    if (nlen > 0) then
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(2 * nlen)
      do n = 1, this%noutlets
        n1 = this%lakein(n)
        n2 = this%lakeout(n)
        if (n1 > 0 .and. n2 > 0) then
          q = this%simoutrate(n)
          if (this%imover == 1) then
            q = q + this%pakmvrobj%get_qtomvr(n)
          end if
          call this%budobj%budterm(idx)%update_term(n1, n2, q)
          call this%budobj%budterm(idx)%update_term(n2, n1, -q)
        end if
      end do
    end if
    !
    ! -- GWF (LEAKAGE)
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%maxbound)
    do n = 1, this%nlakes
      do j = this%idxlakeconn(n), this%idxlakeconn(n + 1) - 1
        n2 = this%cellid(j)
        q = this%qleak(j)
        lkstg = this%xnewpak(n)
        ! -- For the case when the lak stage is exactly equal
        !    to the lake bottom, the wetted area is not returned
        !    equal to 0.0
        gwhead = this%xnew(n2)
        call this%lak_calculate_conn_warea(n, j, lkstg, gwhead, wa)
        ! -- For thermal conduction between a lake and a gw cell,
        !    the shared wetted area should be reset to zero when the lake
        !    stage is below the cell bottom
        if (this%belev(j) > lkstg) wa = DZERO
        this%qauxcbc(1) = wa
        call this%budobj%budterm(idx)%update_term(n, n2, q, this%qauxcbc)
      end do
    end do
    !
    ! -- RAIN
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nlakes)
    do n = 1, this%nlakes
      q = this%precip(n)
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- EVAPORATION
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nlakes)
    do n = 1, this%nlakes
      q = this%evap(n)
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- RUNOFF
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nlakes)
    do n = 1, this%nlakes
      q = this%runoff(n)
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- INFLOW
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nlakes)
    do n = 1, this%nlakes
      q = this%inflow(n)
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- WITHDRAWAL
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nlakes)
    do n = 1, this%nlakes
      q = this%withr(n)
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- EXTERNAL OUTFLOW
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nlakes)
    do n = 1, this%nlakes
      call this%lak_get_external_outlet(n, q)
      ! subtract tomover from external outflow
      call this%lak_get_external_mover(n, v)
      q = q + v
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- STORAGE
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nlakes)
    do n = 1, this%nlakes
      call this%lak_calculate_vol(n, this%xnewpak(n), v1)
      q = this%qsto(n)
      this%qauxcbc(1) = v1
      call this%budobj%budterm(idx)%update_term(n, n, q, this%qauxcbc)
    end do
    !
    ! -- CONSTANT FLOW
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nlakes)
    do n = 1, this%nlakes
      q = this%chterm(n)
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- MOVER
    if (this%imover == 1) then
      !
      ! -- FROM MOVER
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%nlakes)
      do n = 1, this%nlakes
        q = this%pakmvrobj%get_qfrommvr(n)
        call this%budobj%budterm(idx)%update_term(n, n, q)
      end do
      !
      ! -- TO MOVER
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%noutlets)
      do n = 1, this%noutlets
        n1 = this%lakein(n)
        q = this%pakmvrobj%get_qtomvr(n)
        if (q > DZERO) then
          q = -q
        end if
        call this%budobj%budterm(idx)%update_term(n1, n1, q)
      end do
      !
    end if
    !
    ! -- AUXILIARY VARIABLES
    naux = this%naux
    if (naux > 0) then
      idx = idx + 1
      allocate (auxvartmp(naux))
      call this%budobj%budterm(idx)%reset(this%nlakes)
      do n = 1, this%nlakes
        q = DZERO
        do jj = 1, naux
          ii = n
          auxvartmp(jj) = this%lauxvar(jj, ii)
        end do
        call this%budobj%budterm(idx)%update_term(n, n, q, auxvartmp)
      end do
      deallocate (auxvartmp)
    end if
    !
    ! --Terms are filled, now accumulate them for this time step
    call this%budobj%accumulate_terms()
  end subroutine lak_fill_budobj

  !> @brief Set up the table object that is used to write the lak stage data
  !!
  !! The terms listed here must correspond in number and order to the ones
  !! written to the stage table in the lak_ot method
  !<
  subroutine lak_setup_tableobj(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH, LENBUDTXT
    ! -- dummy
    class(LakType) :: this
    ! -- local
    integer(I4B) :: nterms
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    !
    ! -- setup stage table
    if (this%iprhed > 0) then
      !
      ! -- Determine the number of lake stage terms. These are fixed for
      !    the simulation and cannot change.  This includes FLOW-JA-FACE
      !    so they can be written to the binary budget files, but these internal
      !    flows are not included as part of the budget table.
      nterms = 5
      if (this%inamedbound == 1) then
        nterms = nterms + 1
      end if
      !
      ! -- set up table title
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(adjustl(this%packName))//') STAGES FOR EACH CONTROL VOLUME'
      !
      ! -- set up stage tableobj
      call table_cr(this%stagetab, this%packName, title)
      call this%stagetab%table_df(this%nlakes, nterms, this%iout, &
                                  transient=.TRUE.)
      !
      ! -- Go through and set up table budget term
      if (this%inamedbound == 1) then
        text = 'NAME'
        call this%stagetab%initialize_column(text, 20, alignment=TABLEFT)
      end if
      !
      ! -- lake number
      text = 'NUMBER'
      call this%stagetab%initialize_column(text, 10, alignment=TABCENTER)
      !
      ! -- lake stage
      text = 'STAGE'
      call this%stagetab%initialize_column(text, 12, alignment=TABCENTER)
      !
      ! -- lake surface area
      text = 'SURFACE AREA'
      call this%stagetab%initialize_column(text, 12, alignment=TABCENTER)
      !
      ! -- lake wetted area
      text = 'WETTED AREA'
      call this%stagetab%initialize_column(text, 12, alignment=TABCENTER)
      !
      ! -- lake volume
      text = 'VOLUME'
      call this%stagetab%initialize_column(text, 12, alignment=TABCENTER)
    end if
  end subroutine lak_setup_tableobj

  !> @brief Activate addition of density terms
  !<
  subroutine lak_activate_density(this)
    ! -- dummy
    class(LakType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i, j
    !
    ! -- Set idense and reallocate denseterms to be of size MAXBOUND
    this%idense = 1
    call mem_reallocate(this%denseterms, 3, this%MAXBOUND, 'DENSETERMS', &
                        this%memoryPath)
    do i = 1, this%maxbound
      do j = 1, 3
        this%denseterms(j, i) = DZERO
      end do
    end do
    write (this%iout, '(/1x,a)') 'DENSITY TERMS HAVE BEEN ACTIVATED FOR LAKE &
      &PACKAGE: '//trim(adjustl(this%packName))
  end subroutine lak_activate_density

  !> @brief Activate viscosity terms
  !!
  !! Method to activate addition of viscosity terms for a LAK package reach.
  !<
  subroutine lak_activate_viscosity(this)
    ! -- modules
    use MemoryManagerModule, only: mem_reallocate
    ! -- dummy variables
    class(LakType), intent(inout) :: this !< LakType object
    ! -- local variables
    integer(I4B) :: i
    integer(I4B) :: j
    !
    ! -- Set ivsc and reallocate viscratios to be of size MAXBOUND
    this%ivsc = 1
    call mem_reallocate(this%viscratios, 2, this%MAXBOUND, 'VISCRATIOS', &
                        this%memoryPath)
    do i = 1, this%maxbound
      do j = 1, 2
        this%viscratios(j, i) = DONE
      end do
    end do
    write (this%iout, '(/1x,a)') 'VISCOSITY HAS BEEN ACTIVATED FOR LAK &
      &PACKAGE: '//trim(adjustl(this%packName))
  end subroutine lak_activate_viscosity

  !> @brief Calculate the groundwater-lake density exchange terms
  !!
  !! Arguments are as follows:
  !!     iconn       : lak-gwf connection number
  !!     stage       : lake stage
  !!     head        : gwf head
  !!     cond        : conductance
  !!     botl        : bottom elevation of this connection
  !!     flow        : calculated flow, updated here with density terms
  !!     gwfhcof     : gwf head coefficient, updated here with density terms
  !!     gwfrhs      : gwf right-hand-side value, updated here with density terms
  !!
  !! Member variable used here
  !!     denseterms  : shape (3, MAXBOUND), filled by buoyancy package
  !!                     col 1 is relative density of lake (denselak / denseref)
  !!                     col 2 is relative density of gwf cell (densegwf / denseref)
  !!                     col 3 is elevation of gwf cell
  !<
  subroutine lak_calculate_density_exchange(this, iconn, stage, head, cond, &
                                            botl, flow, gwfhcof, gwfrhs)
    ! -- dummy
    class(LakType), intent(inout) :: this
    integer(I4B), intent(in) :: iconn
    real(DP), intent(in) :: stage
    real(DP), intent(in) :: head
    real(DP), intent(in) :: cond
    real(DP), intent(in) :: botl
    real(DP), intent(inout) :: flow
    real(DP), intent(inout) :: gwfhcof
    real(DP), intent(inout) :: gwfrhs
    ! -- local
    real(DP) :: ss
    real(DP) :: hh
    real(DP) :: havg
    real(DP) :: rdenselak
    real(DP) :: rdensegwf
    real(DP) :: rdenseavg
    real(DP) :: elevlak
    real(DP) :: elevgwf
    real(DP) :: elevavg
    real(DP) :: d1
    real(DP) :: d2
    logical(LGP) :: stage_below_bot
    logical(LGP) :: head_below_bot
    !
    ! -- Set lak density to lak density or gwf density
    if (stage >= botl) then
      ss = stage
      stage_below_bot = .false.
      rdenselak = this%denseterms(1, iconn) ! lak rel density
    else
      ss = botl
      stage_below_bot = .true.
      rdenselak = this%denseterms(2, iconn) ! gwf rel density
    end if
    !
    ! -- set hh to head or botl
    if (head >= botl) then
      hh = head
      head_below_bot = .false.
      rdensegwf = this%denseterms(2, iconn) ! gwf rel density
    else
      hh = botl
      head_below_bot = .true.
      rdensegwf = this%denseterms(1, iconn) ! lak rel density
    end if
    !
    ! -- todo: hack because denseterms not updated in a cf calculation
    if (rdensegwf == DZERO) return
    !
    ! -- Update flow
    if (stage_below_bot .and. head_below_bot) then
      !
      ! -- flow is zero, so no terms are updated
      !
    else
      !
      ! -- calculate average relative density
      rdenseavg = DHALF * (rdenselak + rdensegwf)
      !
      ! -- Add contribution of first density term:
      !      cond * (denseavg/denseref - 1) * (hgwf - hlak)
      d1 = cond * (rdenseavg - DONE)
      gwfhcof = gwfhcof - d1
      gwfrhs = gwfrhs - d1 * ss
      d1 = d1 * (hh - ss)
      flow = flow + d1
      !
      ! -- Add second density term if stage and head not below bottom
      if (.not. stage_below_bot .and. .not. head_below_bot) then
        !
        ! -- Add contribution of second density term:
        !      cond * (havg - elevavg) * (densegwf - denselak) / denseref
        elevgwf = this%denseterms(3, iconn)
        if (this%ictype(iconn) == 0 .or. this%ictype(iconn) == 3) then
          ! -- vertical or embedded vertical connection
          elevlak = botl
        else
          ! -- horizontal or embedded horizontal connection
          elevlak = elevgwf
        end if
        elevavg = DHALF * (elevlak + elevgwf)
        havg = DHALF * (hh + ss)
        d2 = cond * (havg - elevavg) * (rdensegwf - rdenselak)
        gwfrhs = gwfrhs + d2
        flow = flow + d2
      end if
    end if
  end subroutine lak_calculate_density_exchange

end module LakModule
