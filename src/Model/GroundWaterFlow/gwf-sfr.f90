!> @brief This module contains the SFR package methods
!!
!! This module contains the overridden methods for the streamflow routing (SFR)
!! package. Most of the methods in the base Boundary Package are overridden.
!!
!<
module SfrModule
  !
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENBOUNDNAME, LENTIMESERIESNAME, &
                             MAXADPIT, &
                             DZERO, DPREC, DEM30, DEM6, DEM5, DEM4, DEM2, &
                             DONETHIRD, DHALF, DP6, DTWOTHIRDS, DP7, &
                             DP9, DP99, DP999, &
                             DONE, D1P1, DFIVETHIRDS, DTWO, DPI, DEIGHT, &
                             DHUNDRED, DEP20, &
                             NAMEDBOUNDFLAG, LENBOUNDNAME, LENFTYPE, &
                             LENPACKAGENAME, LENPAKLOC, MAXCHARLEN, &
                             LENBUDTXT, &
                             DHNOFLO, DHDRY, DNODATA, &
                             TABLEFT, TABCENTER, TABRIGHT, &
                             MNORMAL
  use SmoothingModule, only: sQuadraticSaturation, sQSaturation, &
                             sQuadraticSaturationDerivative, &
                             sQSaturationDerivative, &
                             sCubicSaturation, sChSmooth
  use BndModule, only: BndType
  use BudgetObjectModule, only: BudgetObjectType, budgetobject_cr
  use TableModule, only: TableType, table_cr
  use ObserveModule, only: ObserveType
  use InputOutputModule, only: extract_idnum_or_bndname, upcase
  use BaseDisModule, only: DisBaseType
  use SimModule, only: count_errors, store_error, store_error_unit, &
                       store_warning, deprecation_warning
  use SimVariablesModule, only: errmsg, warnmsg
  use GwfSfrCrossSectionUtilsModule, only: get_saturated_topwidth, &
                                           get_wetted_topwidth, &
                                           get_wetted_perimeter, &
                                           get_cross_section_area, &
                                           get_mannings_section
  use dag_module, only: dag
  use MatrixBaseModule
  !
  implicit none
  !
  character(len=LENFTYPE) :: ftype = 'SFR' !< package ftype string
  character(len=LENPACKAGENAME) :: text = '             SFR' !< package budget string
  !
  private
  public :: sfr_create
  public :: SfrType
  !
  type, extends(BndType) :: SfrType
    ! -- scalars
    ! -- for budgets
    ! -- characters
    character(len=16), dimension(:), pointer, contiguous :: csfrbudget => NULL() !< advanced package budget names
    character(len=16), dimension(:), pointer, contiguous :: cauxcbc => NULL() !< aux names
    character(len=LENBOUNDNAME), dimension(:), pointer, &
      contiguous :: sfrname => null() !< internal SFR reach name
    ! -- integers
    integer(I4B), pointer :: istorage => null() !< flag for using kinematic wave approximation
    integer(I4B), pointer :: iprhed => null() !< flag for printing stages to listing file
    integer(I4B), pointer :: istageout => null() !< flag and unit number for binary stage output
    integer(I4B), pointer :: ibudgetout => null() !< flag and unit number for binary sfr budget output
    integer(I4B), pointer :: ibudcsv => null() !< unit number for csv budget output file
    integer(I4B), pointer :: ipakcsv => null() !< flag and unit number for package convergence information
    integer(I4B), pointer :: idiversions => null() !< flag indicating if there are any diversions
    integer(I4B), pointer :: nconn => NULL() !< number of reach connections
    integer(I4B), pointer :: maxsfrpicard => NULL() !< maximum number of Picard iteration calls to SFR solve
    integer(I4B), pointer :: maxsfrit => NULL() !< maximum number of iterations in SFR solve
    integer(I4B), pointer :: bditems => NULL() !< number of SFR budget items
    integer(I4B), pointer :: cbcauxitems => NULL() !< number of aux items in cell-by-cell budget file
    integer(I4B), pointer :: icheck => NULL() !< flag indicating if input should be checked (default is yes)
    integer(I4B), pointer :: iconvchk => NULL() !< flag indicating of final convergence run is executed
    integer(I4B), pointer :: gwfiss => NULL() !< groundwater model steady-state flag
    integer(I4B), pointer :: ianynone => null() !< number of reaches with 'none' connection
    ! -- double precision
    real(DP), pointer :: unitconv => NULL() !< unit conversion factor (SI to model units)
    real(DP), pointer :: lengthconv => NULL() !< length conversion factor (SI to model units)
    real(DP), pointer :: timeconv => NULL() !< time conversion factor (SI to model units)
    real(DP), pointer :: dmaxchg => NULL() !< maximum depth change allowed
    real(DP), pointer :: deps => NULL() !< perturbation value
    real(DP), pointer :: storage_weight => NULL() !< time weighting factor for kinematic wave approximation
    ! -- integer vectors
    integer(I4B), dimension(:), pointer, contiguous :: isfrorder => null() !< sfr reach order determined from DAG of upstream reaches
    integer(I4B), dimension(:), pointer, contiguous :: ia => null() !< CRS row pointer for SFR reaches
    integer(I4B), dimension(:), pointer, contiguous :: ja => null() !< CRS column pointers for SFR reach connections
    ! -- double precision output vectors
    real(DP), dimension(:), pointer, contiguous :: qoutflow => null() !< reach downstream flow
    real(DP), dimension(:), pointer, contiguous :: qextoutflow => null() !< reach discharge to external boundary
    real(DP), dimension(:), pointer, contiguous :: qauxcbc => null() !< aux value
    real(DP), dimension(:), pointer, contiguous :: dbuff => null() !< temporary vector
    !
    ! -- sfr budget object
    type(BudgetObjectType), pointer :: budobj => null() !< SFR budget object
    !
    ! -- sfr table objects
    type(TableType), pointer :: stagetab => null() !< reach stage table written to the listing file
    type(TableType), pointer :: pakcsvtab => null() !< SFR package convergence table
    !
    ! -- sfr reach data
    integer(I4B), dimension(:), pointer, contiguous :: iboundpak => null() !< ibound array for SFR reaches that defines active, inactive, and constant reaches
    integer(I4B), dimension(:), pointer, contiguous :: igwfnode => null() !< groundwater node connected to SFR reaches
    integer(I4B), dimension(:), pointer, contiguous :: igwftopnode => null() !< highest active groundwater node under SFR reaches
    real(DP), dimension(:), pointer, contiguous :: length => null() !< reach length
    real(DP), dimension(:), pointer, contiguous :: width => null() !< reach width
    real(DP), dimension(:), pointer, contiguous :: strtop => null() !< reach bed top elevation
    real(DP), dimension(:), pointer, contiguous :: bthick => null() !< reach bed thickness
    real(DP), dimension(:), pointer, contiguous :: hk => null() !< vertical hydraulic conductivity of reach bed sediments
    real(DP), dimension(:), pointer, contiguous :: slope => null() !< reach slope
    integer(I4B), dimension(:), pointer, contiguous :: nconnreach => null() !< number of connections for each reach
    real(DP), dimension(:), pointer, contiguous :: ustrf => null() !< upstream flow fraction for upstream connections
    real(DP), dimension(:), pointer, contiguous :: ftotnd => null() !< total fraction of connected reaches that are not diversions
    integer(I4B), dimension(:), pointer, contiguous :: ndiv => null() !< number of diversions for each reach
    real(DP), dimension(:), pointer, contiguous :: usflow => null() !< upstream reach flow
    real(DP), dimension(:), pointer, contiguous :: dsflow => null() !< downstream reach flow
    real(DP), dimension(:), pointer, contiguous :: dsflowold => null() !< downstream reach flow for previous time step
    real(DP), dimension(:), pointer, contiguous :: usinflow => null() !< upstream reach flow for previous time step
    real(DP), dimension(:), pointer, contiguous :: usinflowold => null() !< upstream reach flow for previous time step
    real(DP), dimension(:), pointer, contiguous :: depth => null() !< reach depth
    real(DP), dimension(:), pointer, contiguous :: stage => null() !< reach stage
    real(DP), dimension(:), pointer, contiguous :: stageold => null() !< reach stage for last timestep
    real(DP), dimension(:), pointer, contiguous :: gwflow => null() !< flow from groundwater to reach
    real(DP), dimension(:), pointer, contiguous :: simevap => null() !< simulated reach evaporation
    real(DP), dimension(:), pointer, contiguous :: simrunoff => null() !< simulated reach runoff
    real(DP), dimension(:), pointer, contiguous :: stage0 => null() !< previous reach stage iterate
    real(DP), dimension(:), pointer, contiguous :: usflow0 => null() !< previous upstream reach flow iterate
    real(DP), dimension(:), pointer, contiguous :: storage => null() !< previous upstream reach flow iterate
    ! -- cross-section data
    integer(I4B), pointer :: ncrossptstot => null() !< total number of cross-section points
    integer(I4B), dimension(:), pointer, contiguous :: ncrosspts => null() !< number of cross-section points for each reach
    integer(I4B), dimension(:), pointer, contiguous :: iacross => null() !< pointers to cross-section data for each reach
    real(DP), dimension(:), pointer, contiguous :: station => null() !< cross-section station (x-position) data
    real(DP), dimension(:), pointer, contiguous :: xsheight => null() !< cross-section height data
    real(DP), dimension(:), pointer, contiguous :: xsrough => null() !< cross-section roughness data
    ! -- connection data
    integer(I4B), dimension(:), pointer, contiguous :: idir => null() !< reach connection direction
    integer(I4B), dimension(:), pointer, contiguous :: idiv => null() !< reach connection diversion number
    real(DP), dimension(:), pointer, contiguous :: qconn => null() !< reach connection flow
    ! -- boundary data
    real(DP), dimension(:), pointer, contiguous :: rough => null() !< reach Manning's roughness coefficient (SI units)
    real(DP), dimension(:), pointer, contiguous :: rain => null() !< reach rainfall
    real(DP), dimension(:), pointer, contiguous :: evap => null() !< reach potential evaporation
    real(DP), dimension(:), pointer, contiguous :: inflow => null() !< reach upstream inflow
    real(DP), dimension(:), pointer, contiguous :: runoff => null() !< reach maximum runoff
    real(DP), dimension(:), pointer, contiguous :: sstage => null() !< reach specified stage
    ! -- reach aux variables
    real(DP), dimension(:, :), pointer, contiguous :: rauxvar => null() !< reach aux variable
    ! -- diversion data
    integer(I4B), dimension(:), pointer, contiguous :: iadiv => null() !< row pointer for reach diversions
    integer(I4B), dimension(:), pointer, contiguous :: divreach => null() !< diversion reach
    character(len=10), dimension(:), pointer, contiguous :: divcprior => null() !< diversion rule
    real(DP), dimension(:), pointer, contiguous :: divflow => null() !< specified diversion flow value
    real(DP), dimension(:), pointer, contiguous :: divq => null() !< simulated diversion flow
    !
    ! -- density variables
    integer(I4B), pointer :: idense !< flag indicating if density corrections are active
    real(DP), dimension(:, :), pointer, contiguous :: denseterms => null() !< density terms
    !
    ! -- viscosity variables
    real(DP), dimension(:, :), pointer, contiguous :: viscratios => null() !< viscosity ratios (1: sfr vsc ratio; 2: gwf vsc ratio)
    !
    ! -- type bound procedures
  contains
    procedure :: sfr_allocate_scalars
    procedure :: sfr_allocate_arrays
    procedure :: bnd_options => sfr_options
    procedure :: read_dimensions => sfr_read_dimensions
    ! procedure :: set_pointers => sfr_set_pointers
    procedure :: bnd_ar => sfr_ar
    procedure :: bnd_rp => sfr_rp
    procedure :: bnd_ad => sfr_ad
    procedure :: bnd_cf => sfr_cf
    procedure :: bnd_fc => sfr_fc
    procedure :: bnd_fn => sfr_fn
    procedure :: bnd_cc => sfr_cc
    procedure :: bnd_cq => sfr_cq
    procedure :: bnd_ot_package_flows => sfr_ot_package_flows
    procedure :: bnd_ot_dv => sfr_ot_dv
    procedure :: bnd_ot_bdsummary => sfr_ot_bdsummary
    procedure :: bnd_da => sfr_da
    procedure :: define_listlabel
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => sfr_obs_supported
    procedure, public :: bnd_df_obs => sfr_df_obs
    procedure, public :: bnd_rp_obs => sfr_rp_obs
    procedure, public :: bnd_bd_obs => sfr_bd_obs
    ! -- private procedures
    procedure, private :: sfr_set_stressperiod
    procedure, private :: sfr_solve
    procedure, private :: sfr_calc_constant
    procedure, private :: sfr_calc_transient
    procedure, private :: sfr_calc_steady
    procedure, private :: sfr_update_flows
    procedure, private :: sfr_adjust_ro_ev
    procedure, private :: sfr_calc_qgwf
    procedure, private :: sfr_gwf_conn
    procedure, private :: sfr_calc_cond
    procedure, private :: sfr_calc_qman
    procedure, private :: sfr_calc_qd
    procedure, private :: sfr_calc_qsource
    procedure, private :: sfr_calc_div
    ! -- geometry
    procedure, private :: calc_area_wet
    procedure, private :: calc_perimeter_wet
    procedure, private :: calc_surface_area
    procedure, private :: calc_surface_area_wet
    procedure, private :: calc_top_width_wet
    ! -- reading
    procedure, private :: sfr_read_packagedata
    procedure, private :: sfr_read_crossection
    procedure, private :: sfr_read_connectiondata
    procedure, private :: sfr_read_diversions
    procedure, private :: sfr_read_initial_stages
    ! -- calculations
    procedure, private :: sfr_calc_reach_depth
    procedure, private :: sfr_calc_xs_depth
    ! -- error checking
    procedure, private :: sfr_check_conversion
    procedure, private :: sfr_check_storage_weight
    procedure, private :: sfr_check_reaches
    procedure, private :: sfr_check_connections
    procedure, private :: sfr_check_diversions
    procedure, private :: sfr_check_initialstages
    procedure, private :: sfr_check_ustrf
    ! -- budget
    procedure, private :: sfr_setup_budobj
    procedure, private :: sfr_fill_budobj
    ! -- table
    procedure, private :: sfr_setup_tableobj
    ! -- density
    procedure :: sfr_activate_density
    procedure, private :: sfr_calculate_density_exchange
    ! -- viscosity
    procedure :: sfr_activate_viscosity
  end type SfrType

  interface
    module subroutine sfr_calc_steady(this, n, d1, hgwf, &
                                      qu, qi, qfrommvr, qr, qe, qro, &
                                      qgwf, qd)
      class(SfrType) :: this !< SfrType object
      integer(I4B), intent(in) :: n !< reach number
      real(DP), intent(inout) :: d1 !< current reach depth estimate
      real(DP), intent(in) :: hgwf !< head in gw cell
      real(DP), intent(in) :: qu !< reach upstream flow
      real(DP), intent(in) :: qi !< reach specified inflow
      real(DP), intent(in) :: qfrommvr !< reach flow from mover
      real(DP), intent(in) :: qr !< reach rainfall
      real(DP), intent(in) :: qe !< reach evaporation
      real(DP), intent(in) :: qro !< reach runoff flow
      real(DP), intent(inout) :: qgwf !< reach-aquifer exchange
      real(DP), intent(inout) :: qd !< reach outflow
    end subroutine
  end interface

  interface
    module subroutine sfr_calc_transient(this, n, d1, hgwf, &
                                         qu, qi, qfrommvr, qr, qe, qro, &
                                         qgwf, qd)
      class(SfrType) :: this !< SfrType object
      integer(I4B), intent(in) :: n !< reach number
      real(DP), intent(inout) :: d1 !< current reach depth estimate
      real(DP), intent(in) :: hgwf !< head in gw cell
      real(DP), intent(in) :: qu !< reach upstream flow
      real(DP), intent(in) :: qi !< reach specified inflow
      real(DP), intent(in) :: qfrommvr !< reach flow from mover
      real(DP), intent(in) :: qr !< reach rainfall
      real(DP), intent(in) :: qe !< reach evaporation
      real(DP), intent(in) :: qro !< reach runoff flow
      real(DP), intent(inout) :: qgwf !< reach-aquifer exchange
      real(DP), intent(inout) :: qd !< reach outflow
    end subroutine
  end interface

  interface
    module subroutine sfr_calc_constant(this, n, d1, hgwf, qgwf, qd)
      class(SfrType) :: this !< SfrType object
      integer(I4B), intent(in) :: n !< reach number
      real(DP), intent(inout) :: d1 !< current reach depth estimate
      real(DP), intent(in) :: hgwf !< head in gw cell
      real(DP), intent(inout) :: qgwf !< reach-aquifer exchange
      real(DP), intent(inout) :: qd !< reach outflow
    end subroutine
  end interface

contains

  !> @ brief Create a new package object
  !!
  !!  Create a new SFR Package object
  !<
  subroutine sfr_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname)
    ! -- modules
    use MemoryHelperModule, only: create_mem_path
    ! -- dummy
    class(BndType), pointer :: packobj !< pointer to default package type
    integer(I4B), intent(in) :: id !< package id
    integer(I4B), intent(in) :: ibcnum !< boundary condition number
    integer(I4B), intent(in) :: inunit !< unit number of SFR package input file
    integer(I4B), intent(in) :: iout !< unit number of model listing file
    character(len=*), intent(in) :: namemodel !< model name
    character(len=*), intent(in) :: pakname !< package name
    ! -- local
    type(SfrType), pointer :: sfrobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (sfrobj)
    packobj => sfrobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype)
    packobj%text = text
    !
    ! -- allocate scalars
    call sfrobj%sfr_allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()

    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 4
    packobj%iscloc = 0 ! not supported
    packobj%isadvpak = 1
    packobj%ictMemPath = create_mem_path(namemodel, 'NPF')
  end subroutine sfr_create

  !> @ brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the SFR package. The base model
  !! allocate scalars method is also called.
  !<
  subroutine sfr_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    use MemoryHelperModule, only: create_mem_path
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    !
    ! -- call standard BndType allocate scalars
    call this%BndType%allocate_scalars()
    !
    ! -- allocate the object and assign values to object variables
    call mem_allocate(this%istorage, 'ISTORAGE', this%memoryPath)
    call mem_allocate(this%iprhed, 'IPRHED', this%memoryPath)
    call mem_allocate(this%istageout, 'ISTAGEOUT', this%memoryPath)
    call mem_allocate(this%ibudgetout, 'IBUDGETOUT', this%memoryPath)
    call mem_allocate(this%ibudcsv, 'IBUDCSV', this%memoryPath)
    call mem_allocate(this%ipakcsv, 'IPAKCSV', this%memoryPath)
    call mem_allocate(this%idiversions, 'IDIVERSIONS', this%memoryPath)
    call mem_allocate(this%maxsfrpicard, 'MAXSFRPICARD', this%memoryPath)
    call mem_allocate(this%maxsfrit, 'MAXSFRIT', this%memoryPath)
    call mem_allocate(this%bditems, 'BDITEMS', this%memoryPath)
    call mem_allocate(this%cbcauxitems, 'CBCAUXITEMS', this%memoryPath)
    call mem_allocate(this%unitconv, 'UNITCONV', this%memoryPath)
    call mem_allocate(this%lengthconv, 'LENGTHCONV', this%memoryPath)
    call mem_allocate(this%timeconv, 'TIMECONV', this%memoryPath)
    call mem_allocate(this%dmaxchg, 'DMAXCHG', this%memoryPath)
    call mem_allocate(this%deps, 'DEPS', this%memoryPath)
    call mem_allocate(this%storage_weight, 'STORAGE_WEIGHT', this%memoryPath)
    call mem_allocate(this%nconn, 'NCONN', this%memoryPath)
    call mem_allocate(this%icheck, 'ICHECK', this%memoryPath)
    call mem_allocate(this%iconvchk, 'ICONVCHK', this%memoryPath)
    call mem_allocate(this%idense, 'IDENSE', this%memoryPath)
    call mem_allocate(this%ianynone, 'IANYNONE', this%memoryPath)
    call mem_allocate(this%ncrossptstot, 'NCROSSPTSTOT', this%memoryPath)
    !
    ! -- set pointer to gwf iss
    call mem_setptr(this%gwfiss, 'ISS', create_mem_path(this%name_model))
    !
    ! -- Set values
    this%istorage = 0
    this%iprhed = 0
    this%istageout = 0
    this%ibudgetout = 0
    this%ibudcsv = 0
    this%ipakcsv = 0
    this%idiversions = 0
    this%maxsfrpicard = 100
    this%maxsfrit = MAXADPIT
    this%bditems = 8
    this%cbcauxitems = 1
    this%unitconv = DONE
    this%lengthconv = DNODATA
    this%timeconv = DNODATA
    this%dmaxchg = DEM5
    this%deps = DP999 * this%dmaxchg
    this%storage_weight = DNODATA
    this%nconn = 0
    this%icheck = 1
    this%iconvchk = 1
    this%idense = 0
    this%ivsc = 0
    this%ianynone = 0
    this%ncrossptstot = 0
  end subroutine sfr_allocate_scalars

  !> @ brief Allocate arrays
  !!
  !! Allocate and initialize array for the SFR package.
  !<
  subroutine sfr_allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: j
    !
    ! -- allocate character array for budget text
    allocate (this%csfrbudget(this%bditems))
    call mem_allocate(this%sfrname, LENBOUNDNAME, this%maxbound, &
                      'SFRNAME', this%memoryPath)
    !
    ! -- variables originally in SfrDataType
    call mem_allocate(this%iboundpak, this%maxbound, 'IBOUNDPAK', &
                      this%memoryPath)
    call mem_allocate(this%igwfnode, this%maxbound, 'IGWFNODE', this%memoryPath)
    call mem_allocate(this%igwftopnode, this%maxbound, 'IGWFTOPNODE', &
                      this%memoryPath)
    call mem_allocate(this%length, this%maxbound, 'LENGTH', this%memoryPath)
    call mem_allocate(this%width, this%maxbound, 'WIDTH', this%memoryPath)
    call mem_allocate(this%strtop, this%maxbound, 'STRTOP', this%memoryPath)
    call mem_allocate(this%bthick, this%maxbound, 'BTHICK', this%memoryPath)
    call mem_allocate(this%hk, this%maxbound, 'HK', this%memoryPath)
    call mem_allocate(this%slope, this%maxbound, 'SLOPE', this%memoryPath)
    call mem_allocate(this%nconnreach, this%maxbound, 'NCONNREACH', &
                      this%memoryPath)
    call mem_allocate(this%ustrf, this%maxbound, 'USTRF', this%memoryPath)
    call mem_allocate(this%ftotnd, this%maxbound, 'FTOTND', this%memoryPath)
    call mem_allocate(this%ndiv, this%maxbound, 'NDIV', this%memoryPath)
    call mem_allocate(this%usflow, this%maxbound, 'USFLOW', this%memoryPath)
    call mem_allocate(this%dsflow, this%maxbound, 'DSFLOW', this%memoryPath)
    call mem_allocate(this%depth, this%maxbound, 'DEPTH', this%memoryPath)
    call mem_allocate(this%stage, this%maxbound, 'STAGE', this%memoryPath)
    call mem_allocate(this%gwflow, this%maxbound, 'GWFLOW', this%memoryPath)
    call mem_allocate(this%simevap, this%maxbound, 'SIMEVAP', this%memoryPath)
    call mem_allocate(this%simrunoff, this%maxbound, 'SIMRUNOFF', &
                      this%memoryPath)
    call mem_allocate(this%stage0, this%maxbound, 'STAGE0', this%memoryPath)
    call mem_allocate(this%usflow0, this%maxbound, 'USFLOW0', this%memoryPath)
    !
    ! -- stage, usflow, inflow, and dsflow for previous timestep
    if (this%istorage == 1) then
      call mem_allocate(this%stageold, this%maxbound, 'STAGEOLD', &
                        this%memoryPath)
      call mem_allocate(this%usinflow, this%maxbound, 'USINFLOW', &
                        this%memoryPath)
      call mem_allocate(this%usinflowold, this%maxbound, 'USINFLOWOLD', &
                        this%memoryPath)
      call mem_allocate(this%dsflowold, this%maxbound, 'DSFLOWOLD', &
                        this%memoryPath)
      call mem_allocate(this%storage, this%maxbound, 'STORAGE', &
                        this%memoryPath)
    end if
    !
    ! -- reach order and connection data
    call mem_allocate(this%isfrorder, this%maxbound, 'ISFRORDER', &
                      this%memoryPath)
    call mem_allocate(this%ia, this%maxbound + 1, 'IA', this%memoryPath)
    call mem_allocate(this%ja, 0, 'JA', this%memoryPath)
    call mem_allocate(this%idir, 0, 'IDIR', this%memoryPath)
    call mem_allocate(this%idiv, 0, 'IDIV', this%memoryPath)
    call mem_allocate(this%qconn, 0, 'QCONN', this%memoryPath)
    !
    ! -- boundary data
    call mem_allocate(this%rough, this%maxbound, 'ROUGH', this%memoryPath)
    call mem_allocate(this%rain, this%maxbound, 'RAIN', this%memoryPath)
    call mem_allocate(this%evap, this%maxbound, 'EVAP', this%memoryPath)
    call mem_allocate(this%inflow, this%maxbound, 'INFLOW', this%memoryPath)
    call mem_allocate(this%runoff, this%maxbound, 'RUNOFF', this%memoryPath)
    call mem_allocate(this%sstage, this%maxbound, 'SSTAGE', this%memoryPath)
    !
    ! -- aux variables
    call mem_allocate(this%rauxvar, this%naux, this%maxbound, &
                      'RAUXVAR', this%memoryPath)
    !
    ! -- diversion variables
    call mem_allocate(this%iadiv, this%maxbound + 1, 'IADIV', this%memoryPath)
    call mem_allocate(this%divreach, 0, 'DIVREACH', this%memoryPath)
    call mem_allocate(this%divflow, 0, 'DIVFLOW', this%memoryPath)
    call mem_allocate(this%divq, 0, 'DIVQ', this%memoryPath)
    !
    ! -- cross-section data
    call mem_allocate(this%ncrosspts, this%maxbound, 'NCROSSPTS', &
                      this%memoryPath)
    call mem_allocate(this%iacross, this%maxbound + 1, 'IACROSS', &
                      this%memoryPath)
    call mem_allocate(this%station, this%ncrossptstot, 'STATION', &
                      this%memoryPath)
    call mem_allocate(this%xsheight, this%ncrossptstot, 'XSHEIGHT', &
                      this%memoryPath)
    call mem_allocate(this%xsrough, this%ncrossptstot, 'XSROUGH', &
                      this%memoryPath)
    !
    ! -- initialize variables
    this%iacross(1) = 0
    do i = 1, this%maxbound
      this%iboundpak(i) = 1
      this%igwfnode(i) = 0
      this%igwftopnode(i) = 0
      this%length(i) = DZERO
      this%width(i) = DZERO
      this%strtop(i) = DZERO
      this%bthick(i) = DZERO
      this%hk(i) = DZERO
      this%slope(i) = DZERO
      this%nconnreach(i) = 0
      this%ustrf(i) = DZERO
      this%ftotnd(i) = DZERO
      this%ndiv(i) = 0
      this%usflow(i) = DZERO
      this%dsflow(i) = DZERO
      this%depth(i) = DZERO
      this%stage(i) = DZERO
      this%gwflow(i) = DZERO
      this%simevap(i) = DZERO
      this%simrunoff(i) = DZERO
      this%stage0(i) = DZERO
      this%usflow0(i) = DZERO
      !
      ! -- stage
      if (this%istorage == 1) then
        this%stageold(i) = DZERO
        this%usinflow(i) = DZERO
        this%usinflowold(i) = DZERO
        this%dsflowold(i) = DZERO
        this%storage(i) = DZERO
      end if
      !
      ! -- boundary data
      this%rough(i) = DZERO
      this%rain(i) = DZERO
      this%evap(i) = DZERO
      this%inflow(i) = DZERO
      this%runoff(i) = DZERO
      this%sstage(i) = DZERO
      !
      ! -- aux variables
      do j = 1, this%naux
        this%rauxvar(j, i) = DZERO
      end do
      !
      ! -- cross-section data
      this%ncrosspts(i) = 0
      this%iacross(i + 1) = 0
    end do
    !
    ! -- initialize additional cross-section data
    do i = 1, this%ncrossptstot
      this%station(i) = DZERO
      this%xsheight(i) = DZERO
      this%xsrough(i) = DZERO
    end do
    !
    !-- fill csfrbudget
    this%csfrbudget(1) = '        RAINFALL'
    this%csfrbudget(2) = '     EVAPORATION'
    this%csfrbudget(3) = '          RUNOFF'
    this%csfrbudget(4) = '      EXT-INFLOW'
    this%csfrbudget(5) = '             GWF'
    this%csfrbudget(6) = '     EXT-OUTFLOW'
    this%csfrbudget(7) = '        FROM-MVR'
    this%csfrbudget(8) = '          TO-MVR'
    !
    ! -- allocate and initialize budget output data
    call mem_allocate(this%qoutflow, this%maxbound, 'QOUTFLOW', this%memoryPath)
    call mem_allocate(this%qextoutflow, this%maxbound, 'QEXTOUTFLOW', &
                      this%memoryPath)
    do i = 1, this%maxbound
      this%qoutflow(i) = DZERO
      this%qextoutflow(i) = DZERO
    end do
    !
    ! -- allocate and initialize dbuff
    if (this%istageout > 0) then
      call mem_allocate(this%dbuff, this%maxbound, 'DBUFF', this%memoryPath)
      do i = 1, this%maxbound
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
    call mem_allocate(this%qauxcbc, this%cbcauxitems, 'QAUXCBC', &
                      this%memoryPath)
    do i = 1, this%cbcauxitems
      this%qauxcbc(i) = DZERO
    end do
    !
    ! -- fill cauxcbc
    this%cauxcbc(1) = 'FLOW-AREA       '
    !
    ! -- allocate denseterms to size 0
    call mem_allocate(this%denseterms, 3, 0, 'DENSETERMS', this%memoryPath)
    !
    ! -- allocate viscratios to size 0
    call mem_allocate(this%viscratios, 2, 0, 'VISCRATIOS', this%memoryPath)
  end subroutine sfr_allocate_arrays

  !> @ brief Read dimensions for package
  !!
  !!  Read dimensions for the SFR package.
  !<
  subroutine sfr_read_dimensions(this)
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    ! -- local
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    logical(LGP) :: isfound
    logical(LGP) :: endOfBlock
    !
    ! -- initialize dimensions to 0
    this%maxbound = 0
    !
    ! -- get dimensions block
    call this%parser%GetBlock('DIMENSIONS', isFound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse dimensions block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') &
        'PROCESSING '//trim(adjustl(this%text))//' DIMENSIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('NREACHES')
          this%maxbound = this%parser%GetInteger()
          write (this%iout, '(4x,a,i0)') 'NREACHES = ', this%maxbound
        case default
          write (errmsg, '(2a)') &
            'Unknown '//trim(this%text)//' dimension: ', trim(keyword)
          call store_error(errmsg)
        end select
      end do
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' DIMENSIONS'
    else
      call store_error('Required dimensions block not found.')
    end if
    !
    ! -- verify dimensions were set
    if (this%maxbound < 1) then
      write (errmsg, '(a)') &
        'NREACHES was not specified or was specified incorrectly.'
      call store_error(errmsg)
    end if
    !
    ! -- write summary of error messages for block
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- Call define_listlabel to construct the list label that is written
    !    when PRINT_INPUT option is used.
    call this%define_listlabel()
    !
    ! -- Define default cross-section data size
    this%ncrossptstot = this%maxbound
    !
    ! -- Allocate arrays in package superclass
    call this%sfr_allocate_arrays()
    !
    ! -- read package data
    call this%sfr_read_packagedata()
    !
    ! -- read cross-section data
    call this%sfr_read_crossection()
    !
    ! -- read connection data
    call this%sfr_read_connectiondata()
    !
    ! -- read diversion data
    call this%sfr_read_diversions()
    !
    ! -- read initial stage data
    call this%sfr_read_initial_stages()
    !
    ! -- setup the budget object
    call this%sfr_setup_budobj()
    !
    ! -- setup the stage table object
    call this%sfr_setup_tableobj()
  end subroutine sfr_read_dimensions

  !> @ brief Read additional options for package
  !!
  !!  Read additional options for SFR package.
  !<
  subroutine sfr_options(this, option, found)
    ! -- modules
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: getunit, assign_iounit, openfile
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    character(len=*), intent(inout) :: option !< option keyword string
    logical(LGP), intent(inout) :: found !< boolean indicating if option found
    ! -- local
    real(DP) :: r
    character(len=MAXCHARLEN) :: fname
    character(len=MAXCHARLEN) :: keyword
    ! -- formats
    character(len=*), parameter :: fmttimeconv = &
      &"(4x, 'TIME CONVERSION VALUE (',g0,') SPECIFIED.')"
    character(len=*), parameter :: fmtlengthconv = &
      &"(4x, 'LENGTH CONVERSION VALUE (',g0,') SPECIFIED.')"
    character(len=*), parameter :: fmtpicard = &
      &"(4x, 'MAXIMUM SFR PICARD ITERATION VALUE (',i0,') SPECIFIED.')"
    character(len=*), parameter :: fmtiter = &
      &"(4x, 'MAXIMUM SFR ITERATION VALUE (',i0,') SPECIFIED.')"
    character(len=*), parameter :: fmtdmaxchg = &
      &"(4x, 'MAXIMUM DEPTH CHANGE VALUE (',g0,') SPECIFIED.')"
    character(len=*), parameter :: fmtsfrbin = &
      "(4x, 'SFR ', 1x, a, 1x, ' WILL BE SAVED TO FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"
    character(len=*), parameter :: fmtstoweight = &
      &"(4x, 'KINEMATIC STORAGE WEIGHT (',g0,') SPECIFIED.')"
    !
    ! -- Check for SFR options
    found = .true.
    select case (option)
    case ('STORAGE')
      this%istorage = 1
      write (this%iout, '(4x,a)') trim(adjustl(this%text))// &
        ' REACH STORAGE IS ACTIVE.'
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
                      form, access, 'REPLACE', MNORMAL)
        write (this%iout, fmtsfrbin) &
          'STAGE', trim(adjustl(fname)), this%istageout
      else
        call store_error('Optional stage keyword must &
                         &be followed by fileout.')
      end if
    case ('BUDGET')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        call assign_iounit(this%ibudgetout, this%inunit, "BUDGET fileout")
        call openfile(this%ibudgetout, this%iout, fname, 'DATA(BINARY)', &
                      form, access, 'REPLACE', MNORMAL)
        write (this%iout, fmtsfrbin) &
          'BUDGET', trim(adjustl(fname)), this%ibudgetout
      else
        call store_error('Optional budget keyword must be '// &
                         'followed by fileout.')
      end if
    case ('BUDGETCSV')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        call assign_iounit(this%ibudcsv, this%inunit, "BUDGETCSV fileout")
        call openfile(this%ibudcsv, this%iout, fname, 'CSV', &
                      filstat_opt='REPLACE')
        write (this%iout, fmtsfrbin) &
          'BUDGET CSV', trim(adjustl(fname)), this%ibudcsv
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
        write (this%iout, fmtsfrbin) &
          'PACKAGE_CONVERGENCE', trim(adjustl(fname)), this%ipakcsv
      else
        call store_error('Optional package_convergence keyword must be '// &
                         'followed by fileout.')
      end if
    case ('UNIT_CONVERSION')
      this%unitconv = this%parser%GetDouble()
      !
      ! -- create warning message
      write (warnmsg, '(a)') &
        'SETTING UNIT_CONVERSION DIRECTLY'
      !
      ! -- create deprecation warning
      call deprecation_warning('OPTIONS', 'UNIT_CONVERSION', '6.4.2', &
                               warnmsg, this%parser%GetUnit())
    case ('LENGTH_CONVERSION')
      this%lengthconv = this%parser%GetDouble()
      write (this%iout, fmtlengthconv) this%lengthconv
    case ('TIME_CONVERSION')
      this%timeconv = this%parser%GetDouble()
      write (this%iout, fmttimeconv) this%timeconv
    case ('MAXIMUM_PICARD_ITERATIONS')
      this%maxsfrpicard = this%parser%GetInteger()
      write (this%iout, fmtpicard) this%maxsfrpicard
    case ('MAXIMUM_ITERATIONS')
      this%maxsfrit = this%parser%GetInteger()
      write (this%iout, fmtiter) this%maxsfrit
    case ('MAXIMUM_DEPTH_CHANGE')
      r = this%parser%GetDouble()
      this%dmaxchg = r
      this%deps = DP999 * r
      write (this%iout, fmtdmaxchg) this%dmaxchg
    case ('MOVER')
      this%imover = 1
      write (this%iout, '(4x,A)') 'MOVER OPTION ENABLED'
      !
      ! -- right now these are options that are only available in the
      !    development version and are not included in the documentation.
      !    These options are only available when IDEVELOPMODE in
      !    constants module is set to 1
    case ('DEV_NO_CHECK')
      call this%parser%DevOpt()
      this%icheck = 0
      write (this%iout, '(4x,A)') 'SFR CHECKS OF REACH GEOMETRY '// &
        'RELATIVE TO MODEL GRID AND '// &
        'REASONABLE PARAMETERS WILL NOT '// &
        'BE PERFORMED.'
    case ('DEV_NO_FINAL_CHECK')
      call this%parser%DevOpt()
      this%iconvchk = 0
      write (this%iout, '(4x,a)') &
        'A FINAL CONVERGENCE CHECK OF THE CHANGE IN STREAM FLOW ROUTING &
        &STAGES AND FLOWS WILL NOT BE MADE'
    case ('DEV_STORAGE_WEIGHT')
      r = this%parser%GetDouble()
      if (r < DHALF .or. r > DONE) then
        write (errmsg, '(a,g0,a)') &
          "STORAGE_WEIGHT SPECIFIED TO BE '", r, &
          "' BUT CANNOT BE LESS THAN 0.5 OR GREATER THAN 1.0"
        call store_error(errmsg)
      else
        this%storage_weight = r
        write (this%iout, fmtstoweight) this%storage_weight
      end if
      !
      ! -- no valid options found
    case default
      !
      ! -- No options found
      found = .false.
    end select
  end subroutine sfr_options

  !> @ brief Allocate and read method for package
  !!
  !!  Method to read and prepare period data for the SFR package.
  !<
  subroutine sfr_ar(this)
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: ierr
    !
    ! -- allocate and read observations
    call this%obs%obs_ar()
    !
    ! -- call standard BndType allocate scalars
    call this%BndType%allocate_arrays()
    !
    ! -- set boundname for each connection
    if (this%inamedbound /= 0) then
      do n = 1, this%maxbound
        this%boundname(n) = this%sfrname(n)
      end do
    end if
    !
    ! -- copy boundname into boundname_cst
    call this%copy_boundname()
    !
    ! -- copy igwfnode into nodelist
    do n = 1, this%maxbound
      this%nodelist(n) = this%igwfnode(n)
    end do
    !
    ! -- check the sfr unit conversion data
    call this%sfr_check_conversion()
    !
    ! -- check the storage_weight
    call this%sfr_check_storage_weight()
    !
    ! -- check the sfr reach data
    call this%sfr_check_reaches()

    ! -- check the connection data
    call this%sfr_check_connections()

    ! -- check the diversion data
    if (this%idiversions /= 0) then
      call this%sfr_check_diversions()
    end if

    ! -- check the diversion data
    if (this%istorage == 1) then
      call this%sfr_check_initialstages()
    end if
    !
    ! -- terminate if errors were detected in any of the static sfr data
    ierr = count_errors()
    if (ierr > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- setup pakmvrobj
    if (this%imover /= 0) then
      allocate (this%pakmvrobj)
      call this%pakmvrobj%ar(this%maxbound, this%maxbound, this%memoryPath)
    end if
  end subroutine sfr_ar

  !> @ brief Read packagedata for the package
  !!
  !!  Method to read packagedata for each reach for the SFR package.
  !<
  subroutine sfr_read_packagedata(this)
    ! -- modules
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    ! -- local
    character(len=LINELENGTH) :: text
    character(len=LINELENGTH) :: cellid
    character(len=10) :: cnum
    character(len=LENBOUNDNAME) :: bndName
    character(len=LENBOUNDNAME) :: bndNameTemp
    character(len=LENBOUNDNAME) :: hkname
    character(len=LENBOUNDNAME) :: manningname
    character(len=LENBOUNDNAME) :: ustrfname
    character(len=50), dimension(:), allocatable :: caux
    integer(I4B) :: n, ierr, ival
    logical(LGP) :: isfound
    logical(LGP) :: endOfBlock
    integer(I4B) :: i
    integer(I4B) :: ii
    integer(I4B) :: jj
    integer(I4B) :: iaux
    integer(I4B) :: nconzero
    integer(I4B) :: ipos
    integer, allocatable, dimension(:) :: nboundchk
    real(DP), pointer :: bndElem => null()
    !
    ! -- allocate space for checking sfr reach data
    allocate (nboundchk(this%maxbound))
    do i = 1, this%maxbound
      nboundchk(i) = 0
    end do
    nconzero = 0
    !
    ! -- allocate local storage for aux variables
    if (this%naux > 0) then
      allocate (caux(this%naux))
    end if
    !
    ! -- read reach data
    call this%parser%GetBlock('PACKAGEDATA', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse reaches block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text))// &
        ' PACKAGEDATA'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        ! -- read reach number
        n = this%parser%GetInteger()

        if (n < 1 .or. n > this%maxbound) then
          write (errmsg, '(a,1x,a,1x,i0)') &
            'Reach number (rno) must be greater than 0 and less', &
            'than or equal to', this%maxbound
          call store_error(errmsg)
          cycle
        end if

        ! -- increment nboundchk
        nboundchk(n) = nboundchk(n) + 1
        !
        ! -- get model node number
        call this%parser%GetCellid(this%dis%ndim, cellid, flag_string=.true.)
        this%igwfnode(n) = this%dis%noder_from_cellid(cellid, this%inunit, &
                                                      this%iout, &
                                                      flag_string=.true., &
                                                      allow_zero=.true.)
        this%igwftopnode(n) = this%igwfnode(n)
        !
        ! -- read the cellid string and determine if 'none' is specified
        if (this%igwfnode(n) < 1) then
          this%ianynone = this%ianynone + 1
          call upcase(cellid)
          if (cellid == 'NONE') then
            call this%parser%GetStringCaps(cellid)
            !
            ! -- create warning message
            write (cnum, '(i0)') n
            warnmsg = 'CELLID for unconnected reach '//trim(cnum)// &
                      ' specified to be NONE. Unconnected reaches '// &
                      'should be specified with a zero for each grid '// &
                      'dimension. For example, for a DIS grid a CELLID '// &
                      'of 0 0 0 should be specified for unconnected reaches'
            !
            ! -- create deprecation warning
            call deprecation_warning('PACKAGEDATA', 'CELLID=NONE', '6.4.3', &
                                     warnmsg, this%parser%GetUnit())
          else

          end if
        end if
        ! -- get reach length
        this%length(n) = this%parser%GetDouble()
        ! -- get reach width
        this%width(n) = this%parser%GetDouble()
        ! -- get reach slope
        this%slope(n) = this%parser%GetDouble()
        ! -- get reach stream bottom
        this%strtop(n) = this%parser%GetDouble()
        ! -- get reach bed thickness
        this%bthick(n) = this%parser%GetDouble()
        ! -- get reach bed hk
        call this%parser%GetStringCaps(hkname)
        ! -- get reach roughness
        call this%parser%GetStringCaps(manningname)
        ! -- get number of connections for reach
        ival = this%parser%GetInteger()
        this%nconnreach(n) = ival
        this%nconn = this%nconn + ival
        if (ival < 0) then
          write (errmsg, '(a,1x,i0,1x,a,i0,a)') &
            'NCON for reach', n, &
            'must be greater than or equal to 0 (', ival, ').'
          call store_error(errmsg)
        else if (ival == 0) then
          nconzero = nconzero + 1
        end if
        ! -- get upstream fraction for reach
        call this%parser%GetString(ustrfname)
        ! -- get number of diversions for reach
        ival = this%parser%GetInteger()
        this%ndiv(n) = ival
        if (ival > 0) then
          this%idiversions = 1
        else if (ival < 0) then
          ival = 0
        end if

        ! -- get aux data
        do iaux = 1, this%naux
          call this%parser%GetString(caux(iaux))
        end do

        ! -- set default bndName
        write (cnum, '(i10.10)') n
        bndName = 'Reach'//cnum

        ! -- get reach name
        if (this%inamedbound /= 0) then
          call this%parser%GetStringCaps(bndNameTemp)
          if (bndNameTemp /= '') then
            bndName = bndNameTemp
          end if
          !this%boundname(n) = bndName
        end if
        this%sfrname(n) = bndName
        !
        ! -- set reach hydraulic conductivity
        text = hkname
        jj = 1 !for 'BEDK'
        bndElem => this%hk(n)
        call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                           this%packName, 'BND', &
                                           this%tsManager, this%iprpak, &
                                           'BEDK')
        !
        ! -- set Mannings
        text = manningname
        jj = 1 !for 'MANNING'
        bndElem => this%rough(n)
        call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                           this%packName, 'BND', &
                                           this%tsManager, this%iprpak, &
                                           'MANNING')
        !
        ! -- set upstream fraction
        text = ustrfname
        jj = 1 ! For 'USTRF'
        bndElem => this%ustrf(n)
        call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                           this%packName, 'BND', &
                                           this%tsManager, this%iprpak, 'USTRF')
        !
        ! -- get aux data
        do jj = 1, this%naux
          text = caux(jj)
          ii = n
          bndElem => this%rauxvar(jj, ii)
          call read_value_or_time_series_adv(text, ii, jj, bndElem, &
                                             this%packName, 'AUX', &
                                             this%tsManager, this%iprpak, &
                                             this%auxname(jj))
        end do
        !
        ! -- initialize sstage to the top of the reach
        !    this value would be used by simple routing reaches
        !    on kper = 1 and kstp = 1 if a stage is not specified
        !    on the status line for the reach
        this%sstage(n) = this%strtop(n)

      end do
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' PACKAGEDATA'
    else
      call store_error('REQUIRED PACKAGEDATA BLOCK NOT FOUND.')
    end if
    !
    ! -- Check to make sure that every reach is specified and that no reach
    !    is specified more than once.
    do i = 1, this%maxbound
      if (nboundchk(i) == 0) then
        write (errmsg, '(a,i0,1x,a)') &
          'Information for reach ', i, 'not specified in packagedata block.'
        call store_error(errmsg)
      else if (nboundchk(i) > 1) then
        write (errmsg, '(a,1x,i0,1x,a,1x,i0)') &
          'Reach information specified', nboundchk(i), 'times for reach', i
        call store_error(errmsg)
      end if
    end do
    deallocate (nboundchk)
    !
    ! -- Submit warning message if any reach has zero connections
    if (nconzero > 0) then
      write (warnmsg, '(a,1x,a,1x,a,1x,i0,1x, a)') &
        'SFR Package', trim(this%packName), &
        'has', nconzero, 'reach(es) with zero connections.'
      call store_warning(warnmsg)
    end if
    !
    ! -- terminate if errors encountered in reach block
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- initialize the cross-section data
    ipos = 1
    this%iacross(1) = ipos
    do i = 1, this%maxbound
      this%ncrosspts(i) = 1
      this%station(ipos) = this%width(i)
      this%xsheight(ipos) = DZERO
      this%xsrough(ipos) = DONE
      ipos = ipos + 1
      this%iacross(i + 1) = ipos
    end do
    !
    ! -- deallocate local storage for aux variables
    if (this%naux > 0) then
      deallocate (caux)
    end if
  end subroutine sfr_read_packagedata

  !> @ brief Read crosssection block for the package
  !!
  !!  Method to read crosssection data for the SFR package.
  !<
  subroutine sfr_read_crossection(this)
    ! -- modules
    use MemoryManagerModule, only: mem_reallocate
    use sfrCrossSectionManager, only: cross_section_cr, SfrCrossSection
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    ! -- local
    character(len=LINELENGTH) :: keyword
    character(len=LINELENGTH) :: line
    logical(LGP) :: isfound
    logical(LGP) :: endOfBlock
    integer(I4B) :: n
    integer(I4B) :: ierr
    integer(I4B) :: ncrossptstot
    integer, allocatable, dimension(:) :: nboundchk
    type(SfrCrossSection), pointer :: cross_data => null()
    !
    ! -- read cross-section data
    call this%parser%GetBlock('CROSSSECTIONS', isfound, ierr, &
                              supportOpenClose=.true., &
                              blockRequired=.false.)
    !
    ! -- parse reach connectivity block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') &
        'PROCESSING '//trim(adjustl(this%text))//' CROSSSECTIONS'
      !
      ! -- allocate and initialize local variables for reach cross-sections
      allocate (nboundchk(this%maxbound))
      do n = 1, this%maxbound
        nboundchk(n) = 0
      end do
      !
      ! -- create and initialize cross-section data
      call cross_section_cr(cross_data, this%iout, this%iprpak, this%maxbound)
      call cross_data%initialize(this%ncrossptstot, this%ncrosspts, &
                                 this%iacross, &
                                 this%station, this%xsheight, &
                                 this%xsrough)
      !
      ! -- read all of the entries in the block
      readtable: do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        !
        ! -- get reach number
        n = this%parser%GetInteger()
        !
        ! -- check for reach number error
        if (n < 1 .or. n > this%maxbound) then
          write (errmsg, '(a,1x,a,1x,i0)') &
            'SFR reach in crosssections block is less than one or greater', &
            'than NREACHES:', n
          call store_error(errmsg)
          cycle readtable
        end if
        !
        ! -- increment nboundchk
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
          call cross_data%read_table(n, this%width(n), &
                                     trim(adjustl(line)))
        case default
          write (errmsg, '(a,1x,i4,1x,a)') &
            'CROSS-SECTION TABLE ENTRY for REACH ', n, &
            'MUST INCLUDE TAB6 KEYWORD'
          call store_error(errmsg)
          cycle readtable
        end select
      end do readtable

      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' CROSSSECTIONS'

      !
      ! -- check for duplicate sfr crosssections
      do n = 1, this%maxbound
        if (nboundchk(n) > 1) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
            'Cross-section data for reach', n, &
            'specified', nboundchk(n), 'times.'
          call store_error(errmsg)
        end if
      end do
      !
      ! -- terminate if errors encountered in cross-sections block
      if (count_errors() > 0) then
        call this%parser%StoreErrorUnit()
      end if
      !
      ! -- determine the current size of cross-section data
      ncrossptstot = cross_data%get_ncrossptstot()
      !
      ! -- reallocate sfr package cross-section data
      if (ncrossptstot /= this%ncrossptstot) then
        this%ncrossptstot = ncrossptstot
        call mem_reallocate(this%station, this%ncrossptstot, 'STATION', &
                            this%memoryPath)
        call mem_reallocate(this%xsheight, this%ncrossptstot, 'XSHEIGHT', &
                            this%memoryPath)
        call mem_reallocate(this%xsrough, this%ncrossptstot, 'XSROUGH', &
                            this%memoryPath)
      end if
      !
      ! -- write cross-section data to the model listing file
      call cross_data%output(this%width, this%rough)
      !
      ! -- pack cross-section data
      call cross_data%pack(this%ncrossptstot, this%ncrosspts, &
                           this%iacross, &
                           this%station, &
                           this%xsheight, &
                           this%xsrough)
      !
      ! -- deallocate temporary local storage for reach cross-sections
      deallocate (nboundchk)
      call cross_data%destroy()
      deallocate (cross_data)
      nullify (cross_data)
    end if
  end subroutine sfr_read_crossection

  !> @ brief Read connectiondata for the package
  !!
  !!  Method to read connectiondata for each reach for the SFR package.
  !<
  subroutine sfr_read_connectiondata(this)
    ! -- modules
    use MemoryManagerModule, only: mem_reallocate
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    ! -- local
    character(len=LINELENGTH) :: line
    logical(LGP) :: isfound
    logical(LGP) :: endOfBlock
    integer(I4B) :: n
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: jj
    integer(I4B) :: jcol
    integer(I4B) :: jcol2
    integer(I4B) :: nja
    integer(I4B) :: ival
    integer(I4B) :: idir
    integer(I4B) :: ierr
    integer(I4B) :: nconnmax
    integer(I4B) :: nup
    integer(I4B) :: ipos
    integer(I4B) :: istat
    integer(I4B), dimension(:), pointer, contiguous :: rowmaxnnz => null()
    integer, allocatable, dimension(:) :: nboundchk
    integer, allocatable, dimension(:, :) :: iconndata
    type(sparsematrix), pointer :: sparse => null()
    integer(I4B), dimension(:), allocatable :: iup
    integer(I4B), dimension(:), allocatable :: order
    type(dag) :: sfr_dag
    !
    ! -- allocate and initialize local variables for reach connections
    allocate (nboundchk(this%maxbound))
    do n = 1, this%maxbound
      nboundchk(n) = 0
    end do
    !
    ! -- calculate the number of non-zero entries (size of ja maxtrix)
    nja = 0
    nconnmax = 0
    allocate (rowmaxnnz(this%maxbound))
    do n = 1, this%maxbound
      ival = this%nconnreach(n)
      if (ival < 0) ival = 0
      rowmaxnnz(n) = ival + 1
      nja = nja + ival + 1
      if (ival > nconnmax) then
        nconnmax = ival
      end if
    end do
    !
    ! -- reallocate connection data for package
    call mem_reallocate(this%ja, nja, 'JA', this%memoryPath)
    call mem_reallocate(this%idir, nja, 'IDIR', this%memoryPath)
    call mem_reallocate(this%idiv, nja, 'IDIV', this%memoryPath)
    call mem_reallocate(this%qconn, nja, 'QCONN', this%memoryPath)
    !
    ! -- initialize connection data
    do n = 1, nja
      this%idir(n) = 0
      this%idiv(n) = 0
      this%qconn(n) = DZERO
    end do
    !
    ! -- allocate space for iconndata
    allocate (iconndata(nconnmax, this%maxbound))
    !
    ! -- initialize iconndata
    do n = 1, this%maxbound
      do j = 1, nconnmax
        iconndata(j, n) = 0
      end do
    end do
    !
    ! -- allocate space for connectivity
    allocate (sparse)
    !
    ! -- set up sparse
    call sparse%init(this%maxbound, this%maxbound, rowmaxnnz)
    !
    ! -- read connection data
    call this%parser%GetBlock('CONNECTIONDATA', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse reach connectivity block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') &
        'PROCESSING '//trim(adjustl(this%text))//' CONNECTIONDATA'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        !
        ! -- get reach number
        n = this%parser%GetInteger()
        !
        ! -- check for error
        if (n < 1 .or. n > this%maxbound) then
          write (errmsg, '(a,1x,a,1x,i0)') &
            'SFR reach in connectiondata block is less than one or greater', &
            'than NREACHES:', n
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- increment nboundchk
        nboundchk(n) = nboundchk(n) + 1
        !
        ! -- add diagonal connection for reach
        call sparse%addconnection(n, n, 1)
        !
        ! -- fill off diagonals
        do i = 1, this%nconnreach(n)
          !
          ! -- get connected reach
          ival = this%parser%GetInteger()
          !
          ! -- save connection data to temporary iconndata
          iconndata(i, n) = ival
          !
          ! -- determine idir
          if (ival < 0) then
            idir = -1
            ival = abs(ival)
          elseif (ival == 0) then
            call store_error('Missing or zero connection reach in line:')
            call store_error(line)
          else
            idir = 1
          end if
          if (ival > this%maxbound) then
            call store_error('Reach number exceeds NREACHES in line:')
            call store_error(line)
          end if
          !
          ! -- add connection to sparse
          call sparse%addconnection(n, ival, 1)
        end do
      end do

      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' CONNECTIONDATA'

      do n = 1, this%maxbound
        !
        ! -- check for missing or duplicate sfr connections
        if (nboundchk(n) == 0) then
          write (errmsg, '(a,1x,i0)') &
            'No connection data specified for reach', n
          call store_error(errmsg)
        else if (nboundchk(n) > 1) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
            'Connection data for reach', n, &
            'specified', nboundchk(n), 'times.'
          call store_error(errmsg)
        end if
      end do
    else
      call store_error('Required connectiondata block not found.')
    end if
    !
    ! -- terminate if errors encountered in connectiondata block
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- create ia and ja from sparse
    call sparse%filliaja(this%ia, this%ja, ierr, sort=.TRUE.)
    !
    ! -- test for error condition
    if (ierr /= 0) then
      write (errmsg, '(a,3(1x,a))') &
        'Could not fill', trim(this%packName), &
        'package IA and JA connection data.', &
        'Check connectivity data in connectiondata block.'
      call store_error(errmsg)
    end if
    !
    ! -- fill flat connection storage
    do n = 1, this%maxbound
      do j = this%ia(n) + 1, this%ia(n + 1) - 1
        jcol = this%ja(j)
        do jj = 1, this%nconnreach(n)
          jcol2 = iconndata(jj, n)
          if (abs(jcol2) == jcol) then
            idir = 1
            if (jcol2 < 0) then
              idir = -1
            end if
            this%idir(j) = idir
            exit
          end if
        end do
      end do
    end do
    !
    ! -- deallocate temporary local storage for reach connections
    deallocate (rowmaxnnz)
    deallocate (nboundchk)
    deallocate (iconndata)
    !
    ! -- destroy sparse
    call sparse%destroy()
    deallocate (sparse)
    !
    ! -- calculate reach order using DAG
    !
    ! -- initialize the DAG
    call sfr_dag%set_vertices(this%maxbound)
    !
    ! -- fill DAG
    fill_dag: do n = 1, this%maxbound
      !
      ! -- determine the number of upstream reaches
      nup = 0
      do j = this%ia(n) + 1, this%ia(n + 1) - 1
        if (this%idir(j) > 0) then
          nup = nup + 1
        end if
      end do
      !
      ! -- cycle if nu upstream reacches
      if (nup == 0) cycle fill_dag
      !
      ! -- allocate local storage
      allocate (iup(nup))
      !
      ! -- fill local storage
      ipos = 1
      do j = this%ia(n) + 1, this%ia(n + 1) - 1
        if (this%idir(j) > 0) then
          iup(ipos) = this%ja(j)
          ipos = ipos + 1
        end if
      end do
      !
      ! -- add upstream connections to DAG
      call sfr_dag%set_edges(n, iup)
      !
      ! -- clean up local storage
      deallocate (iup)
    end do fill_dag
    !
    ! -- perform toposort on DAG
    call sfr_dag%toposort(order, istat)
    !
    ! -- write warning if circular dependency
    if (istat == -1) then
      write (warnmsg, '(a)') &
        trim(adjustl(this%text))//' PACKAGE ('// &
        trim(adjustl(this%packName))//') cannot calculate a '// &
        'Directed Asyclic Graph for reach connectivity because '// &
        'of circular dependency. Using the reach number for '// &
        'solution ordering.'
      call store_warning(warnmsg)
    end if
    !
    ! -- fill isfrorder
    do n = 1, this%maxbound
      if (istat == 0) then
        this%isfrorder(n) = order(n)
      else
        this%isfrorder(n) = n
      end if
    end do
    !
    ! -- clean up DAG and remaining local storage
    call sfr_dag%destroy()
    if (istat == 0) then
      deallocate (order)
    end if
  end subroutine sfr_read_connectiondata

  !> @ brief Read diversions for the package
  !!
  !!  Method to read diversions for the SFR package.
  !<
  subroutine sfr_read_diversions(this)
    ! -- modules
    use MemoryManagerModule, only: mem_reallocate
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    ! -- local
    character(len=10) :: cnum
    character(len=10) :: cval
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: ierr
    integer(I4B) :: ival
    integer(I4B) :: i0
    integer(I4B) :: ipos
    integer(I4B) :: jpos
    integer(I4B) :: ndiv
    integer(I4B) :: ndiversions
    integer(I4B) :: idivreach
    logical(LGP) :: isfound
    logical(LGP) :: endOfBlock
    integer(I4B) :: idiv
    integer, allocatable, dimension(:) :: iachk
    integer, allocatable, dimension(:) :: nboundchk
    !
    ! -- determine the total number of diversions and fill iadiv
    ndiversions = 0
    i0 = 1
    this%iadiv(1) = i0
    do n = 1, this%maxbound
      ndiversions = ndiversions + this%ndiv(n)
      i0 = i0 + this%ndiv(n)
      this%iadiv(n + 1) = i0
    end do
    !
    ! -- reallocate memory for diversions
    if (ndiversions > 0) then
      call mem_reallocate(this%divreach, ndiversions, 'DIVREACH', &
                          this%memoryPath)
      allocate (this%divcprior(ndiversions))
      call mem_reallocate(this%divflow, ndiversions, 'DIVFLOW', this%memoryPath)
      call mem_reallocate(this%divq, ndiversions, 'DIVQ', this%memoryPath)
    end if
    !
    ! -- initialize diversion flow
    do n = 1, ndiversions
      this%divflow(n) = DZERO
      this%divq(n) = DZERO
    end do
    !
    ! -- read diversions
    call this%parser%GetBlock('DIVERSIONS', isfound, ierr, &
                              supportOpenClose=.true., &
                              blockRequired=.false.)
    !
    ! -- parse reach connectivity block if detected
    if (isfound) then
      if (this%idiversions /= 0) then
        write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text))// &
          ' DIVERSIONS'
        !
        ! -- allocate and initialize local variables for diversions
        ndiv = 0
        do n = 1, this%maxbound
          ndiv = ndiv + this%ndiv(n)
        end do
        allocate (iachk(this%maxbound + 1))
        allocate (nboundchk(ndiv))
        iachk(1) = 1
        do n = 1, this%maxbound
          iachk(n + 1) = iachk(n) + this%ndiv(n)
        end do
        do n = 1, ndiv
          nboundchk(n) = 0
        end do
        !
        ! -- read diversion data
        do
          call this%parser%GetNextLine(endOfBlock)
          if (endOfBlock) exit
          !
          ! -- get reach number
          n = this%parser%GetInteger()
          if (n < 1 .or. n > this%maxbound) then
            write (cnum, '(i0)') n
            errmsg = 'Reach number should be between 1 and '// &
                     trim(cnum)//'.'
            call store_error(errmsg)
            cycle
          end if
          !
          ! -- make sure reach has at least one diversion
          if (this%ndiv(n) < 1) then
            write (cnum, '(i0)') n
            errmsg = 'Diversions cannot be specified '// &
                     'for reach '//trim(cnum)
            call store_error(errmsg)
            cycle
          end if
          !
          ! -- read diversion number
          ival = this%parser%GetInteger()
          if (ival < 1 .or. ival > this%ndiv(n)) then
            write (cnum, '(i0)') n
            errmsg = 'Reach  '//trim(cnum)
            write (cnum, '(i0)') this%ndiv(n)
            errmsg = trim(errmsg)//' diversion number should be between '// &
                     '1 and '//trim(cnum)//'.'
            call store_error(errmsg)
            cycle
          end if

          ! -- increment nboundchk
          ipos = iachk(n) + ival - 1
          nboundchk(ipos) = nboundchk(ipos) + 1

          idiv = ival
          !
          ! -- get target reach for diversion
          ival = this%parser%GetInteger()
          if (ival < 1 .or. ival > this%maxbound) then
            write (cnum, '(i0)') ival
            errmsg = 'Diversion target reach number should be '// &
                     'between 1 and '//trim(cnum)//'.'
            call store_error(errmsg)
            cycle
          end if
          idivreach = ival
          jpos = this%iadiv(n) + idiv - 1
          this%divreach(jpos) = idivreach
          !
          ! -- get cprior
          call this%parser%GetStringCaps(cval)
          ival = -1
          select case (cval)
          case ('UPTO')
            ival = 0
          case ('THRESHOLD')
            ival = -1
          case ('FRACTION')
            ival = -2
          case ('EXCESS')
            ival = -3
          case default
            errmsg = 'Invalid cprior type '//trim(cval)//'.'
            call store_error(errmsg)
          end select
          !
          ! -- set cprior for diversion
          this%divcprior(jpos) = cval
        end do

        write (this%iout, '(1x,a)') 'END OF '//trim(adjustl(this%text))// &
          ' DIVERSIONS'

        do n = 1, this%maxbound
          do j = 1, this%ndiv(n)
            ipos = iachk(n) + j - 1
            !
            ! -- check for missing or duplicate reach diversions
            if (nboundchk(ipos) == 0) then
              write (errmsg, '(a,1x,i0,1x,a,1x,i0)') &
                'No data specified for reach', n, 'diversion', j
              call store_error(errmsg)
            else if (nboundchk(ipos) > 1) then
              write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a,1x,i0,1x,a)') &
                'Data for reach', n, 'diversion', j, &
                'specified', nboundchk(ipos), 'times'
              call store_error(errmsg)
            end if
          end do
        end do
        !
        ! -- deallocate local variables
        deallocate (iachk)
        deallocate (nboundchk)
      else
        !
        ! -- error condition
        write (errmsg, '(a,1x,a)') &
          'A diversions block should not be', &
          'specified if diversions are not specified.'
        call store_error(errmsg)
      end if
    else
      if (this%idiversions /= 0) then
        call store_error('REQUIRED DIVERSIONS BLOCK NOT FOUND.')
      end if
    end if
    !
    ! -- write summary of diversion error messages
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine sfr_read_diversions

  !> @ brief Read initialstages data for the package
  !!
  !!  Method to read initialstages data for each reach for the SFR package.
  !<
  subroutine sfr_read_initial_stages(this)
    ! -- modules
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: ierr
    logical(LGP) :: isfound
    logical(LGP) :: endOfBlock
    integer(I4B) :: i
    real(DP) :: rval
    integer, allocatable, dimension(:) :: nboundchk
    !
    ! -- read data
    call this%parser%GetBlock('INITIALSTAGES', isfound, ierr, &
                              supportOpenClose=.true., &
                              blockRequired=.false.)
    !
    ! -- parse block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') &
        'PROCESSING '//trim(adjustl(this%text))//' INITIALSTAGES'

      allocate (nboundchk(this%maxbound))
      do n = 1, this%maxbound
        nboundchk(n) = 0
      end do

      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit

        ! -- read reach number
        n = this%parser%GetInteger()

        if (n < 1 .or. n > this%maxbound) then
          write (errmsg, '(a,i0,a,1x,i0,a)') &
            'Reach number (', n, ') must be greater than 0 and less &
            &than or equal to', this%maxbound, '.'
          call store_error(errmsg)
          cycle
        end if

        ! -- increment nboundchk
        nboundchk(n) = nboundchk(n) + 1

        rval = this%parser%GetDouble()
        this%stage(n) = rval
        this%depth(n) = rval - this%strtop(n)

        if (rval < this%strtop(n)) then
          write (errmsg, '(a,g0,a,1x,i0,1x,a,g0,a)') &
            'Initial stage (', rval, ') for reach', n, &
            'is less than the reach top (', this%strtop(n), ').'
          call store_error(errmsg)
        end if
      end do

      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' INITIALSTAGES'

      !
      ! -- Check to make sure that every reach is specified and that no reach
      !    is specified more than once.
      do i = 1, this%maxbound
        if (nboundchk(i) == 0) then
          write (errmsg, '(a,i0,1x,a)') &
            'Information for reach ', i, 'not specified in initialstages block.'
          call store_error(errmsg)
        else if (nboundchk(i) > 1) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0)') &
            'Initial stage information specified', &
            nboundchk(i), 'times for reach', i
          call store_error(errmsg)
        end if
      end do
      deallocate (nboundchk)
    else
      ! -- set default initial stage based on a zero depth
      if (this%istorage == 1) then
        do n = 1, this%maxbound
          rval = this%strtop(n)
          this%stage(n) = rval
        end do
      end if
    end if
    !
    ! -- terminate if errors encountered in reach block
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine sfr_read_initial_stages

  !> @ brief Read and prepare period data for package
  !!
  !!  Method to read and prepare period data for the SFR package.
  !<
  subroutine sfr_rp(this)
    ! -- modules
    use TdisModule, only: kper, nper
    use MemoryManagerModule, only: mem_reallocate
    use sfrCrossSectionManager, only: cross_section_cr, SfrCrossSection
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: line
    character(len=LINELENGTH) :: crossfile
    integer(I4B) :: ierr
    integer(I4B) :: n
    integer(I4B) :: ichkustrm
    integer(I4B) :: ichkcross
    integer(I4B) :: ncrossptstot
    logical(LGP) :: isfound
    logical(LGP) :: endOfBlock
    type(SfrCrossSection), pointer :: cross_data => null()
    ! -- formats
    character(len=*), parameter :: fmtblkerr = &
      &"('Looking for BEGIN PERIOD iper.  Found ', a, ' instead.')"
    character(len=*), parameter :: fmtlsp = &
      &"(1X,/1X,'REUSING ',A,'S FROM LAST STRESS PERIOD')"
    character(len=*), parameter :: fmtnbd = &
      "(1X,/1X,'The number of active ',A,'S (',I6, &
      &') is greater than maximum (',I6,')')"
    !
    ! -- initialize flags
    ichkustrm = 0
    ichkcross = 0
    if (kper == 1) then
      ichkustrm = 1
    end if
    !
    ! -- set nbound to maxbound
    this%nbound = this%maxbound
    !
    ! -- Set ionper to the stress period number for which a new block of data
    !    will be read.
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
      ! -- create and initialize cross-section data
      call cross_section_cr(cross_data, this%iout, this%iprpak, this%maxbound)
      call cross_data%initialize(this%ncrossptstot, this%ncrosspts, &
                                 this%iacross, &
                                 this%station, this%xsheight, &
                                 this%xsrough)
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
      ! -- read data
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        n = this%parser%GetInteger()
        if (n < 1 .or. n > this%maxbound) then
          write (errmsg, '(a,1x,a,1x,i0,a)') &
            'Reach number (RNO) must be greater than 0 and', &
            'less than or equal to', this%maxbound, '.'
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- read data from the rest of the line
        call this%sfr_set_stressperiod(n, ichkustrm, crossfile)
        !
        ! -- write line to table
        if (this%iprpak /= 0) then
          call this%parser%GetCurrentLine(line)
          call this%inputtab%line_to_columns(line)
        end if
        !
        ! -- process cross-section file
        if (trim(adjustl(crossfile)) /= 'NONE') then
          call cross_data%read_table(n, this%width(n), &
                                     trim(adjustl(crossfile)))
        end if
      end do
      !
      ! -- write raw period data
      if (this%iprpak /= 0) then
        call this%inputtab%finalize_table()
      end if
      !
      ! -- finalize cross-sections

      !
      ! -- determine the current size of cross-section data
      ncrossptstot = cross_data%get_ncrossptstot()
      !
      ! -- reallocate sfr package cross-section data
      if (ncrossptstot /= this%ncrossptstot) then
        this%ncrossptstot = ncrossptstot
        call mem_reallocate(this%station, this%ncrossptstot, 'STATION', &
                            this%memoryPath)
        call mem_reallocate(this%xsheight, this%ncrossptstot, 'XSHEIGHT', &
                            this%memoryPath)
        call mem_reallocate(this%xsrough, this%ncrossptstot, 'XSROUGH', &
                            this%memoryPath)
      end if
      !
      ! -- write cross-section data to the model listing file
      call cross_data%output(this%width, this%rough, kstp=1, kper=kper)
      !
      ! -- pack cross-section data
      call cross_data%pack(this%ncrossptstot, this%ncrosspts, &
                           this%iacross, &
                           this%station, &
                           this%xsheight, &
                           this%xsrough)
      !
      ! -- deallocate temporary local storage for reach cross-sections
      call cross_data%destroy()
      deallocate (cross_data)
      nullify (cross_data)
      !
      ! -- Reuse data from last stress period
    else
      write (this%iout, fmtlsp) trim(this%filtyp)
    end if
    !
    ! -- check upstream fraction values
    if (ichkustrm /= 0) then
      call this%sfr_check_ustrf()
    end if
    !
    ! -- write summary of package block error messages
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine sfr_rp

  !> @ brief Advance the package
  !!
  !!  Advance data in the SFR package. The method sets advances
  !!  time series, time array series, and observation data.
  !<
  subroutine sfr_ad(this)
    ! -- modules
    use TimeSeriesManagerModule, only: var_timeseries
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: iaux

    ! -- update previous values
    if (this%istorage == 1) then
      do n = 1, this%maxbound
        this%stageold(n) = this%stage(n)
        this%usinflowold(n) = this%usinflow(n)
        this%dsflowold(n) = this%dsflow(n)
      end do
    end if
    !
    ! -- Most advanced package AD routines have to restore state if
    !    the solution failed and the time step is being retried with a smaller
    !    step size.  This is not needed here because there is no old stage
    !    or storage effects in the stream.
    !
    ! -- Advance the time series manager
    call this%TsManager%ad()
    !
    ! -- check upstream fractions if time series are being used to
    !    define this variable
    if (var_timeseries(this%tsManager, this%packName, 'USTRF')) then
      call this%sfr_check_ustrf()
    end if
    !
    ! -- update auxiliary variables by copying from the derived-type time
    !    series variable into the bndpackage auxvar variable so that this
    !    information is properly written to the GWF budget file
    if (this%naux > 0) then
      do n = 1, this%maxbound
        do iaux = 1, this%naux
          if (this%noupdateauxvar(iaux) /= 0) cycle
          this%auxvar(iaux, n) = this%rauxvar(iaux, n)
        end do
      end do
    end if
    !
    ! -- reset upstream flow to zero and set specified stage
    do n = 1, this%maxbound
      this%usflow(n) = DZERO
      if (this%iboundpak(n) < 0) then
        this%stage(n) = this%sstage(n)
      end if
    end do
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
  end subroutine sfr_ad

  !> @ brief Formulate the package hcof and rhs terms.
  !!
  !!  Formulate the hcof and rhs terms for the WEL package that will be
  !!  added to the coefficient matrix and right-hand side vector.
  !<
  subroutine sfr_cf(this)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: igwfnode
    !
    ! -- return if no sfr reaches
    if (this%nbound == 0) return
    !
    ! -- find highest active cell
    do n = 1, this%nbound
      igwfnode = this%igwftopnode(n)
      if (igwfnode > 0) then
        if (this%ibound(igwfnode) == 0) then
          call this%dis%highest_active(igwfnode, this%ibound)
        end if
      end if
      this%igwfnode(n) = igwfnode
      this%nodelist(n) = igwfnode
    end do
  end subroutine sfr_cf

  !> @ brief Copy hcof and rhs terms into solution.
  !!
  !!  Add the hcof and rhs terms for the SFR package to the
  !!  coefficient matrix and right-hand side vector.
  !<
  subroutine sfr_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    real(DP), dimension(:), intent(inout) :: rhs !< right-hand side vector for model
    integer(I4B), dimension(:), intent(in) :: ia !< solution CRS row pointers
    integer(I4B), dimension(:), intent(in) :: idxglo !< mapping vector for model (local) to solution (global)
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: ipos
    integer(I4B) :: node
    real(DP) :: s0
    real(DP) :: ds
    real(DP) :: dsmax
    real(DP) :: hgwf
    real(DP) :: v
    real(DP) :: hhcof
    real(DP) :: rrhs
    !
    ! -- picard iterations for sfr to achieve good solution regardless
    !    of reach order
    sfrpicard: do i = 1, this%maxsfrpicard
      !
      ! -- initialize maximum stage change for iteration to zero
      dsmax = DZERO
      !
      ! -- pakmvrobj fc - reset qformvr to zero
      if (this%imover == 1) then
        call this%pakmvrobj%fc()
      end if
      !
      ! -- solve for each sfr reach
      reachsolve: do j = 1, this%nbound
        n = this%isfrorder(j)
        node = this%igwfnode(n)
        if (node > 0) then
          hgwf = this%xnew(node)
        else
          hgwf = DEP20
        end if
        !
        ! -- save previous stage and upstream flow
        if (i == 1) then
          this%stage0(n) = this%stage(n)
          this%usflow0(n) = this%usflow(n)
        end if
        !
        ! -- set initial stage to calculate stage change
        s0 = this%stage(n)
        !
        ! -- solve for flow in swr
        if (this%iboundpak(n) /= 0) then
          call this%sfr_solve(n, hgwf, hhcof, rrhs)
        else
          this%depth(n) = DZERO
          this%stage(n) = this%strtop(n)
          v = DZERO
          call this%sfr_update_flows(n, v, v)
          hhcof = DZERO
          rrhs = DZERO
        end if
        !
        ! -- set package hcof and rhs
        this%hcof(n) = hhcof
        this%rhs(n) = rrhs
        !
        ! -- calculate stage change
        ds = s0 - this%stage(n)
        !
        ! -- evaluate if stage change exceeds dsmax
        if (abs(ds) > abs(dsmax)) then
          dsmax = ds
        end if

      end do reachsolve
      !
      ! -- evaluate if the sfr picard iterations should be terminated
      if (abs(dsmax) <= this%dmaxchg) then
        exit sfrpicard
      end if

    end do sfrpicard
    !
    ! -- Copy package rhs and hcof into solution rhs and amat
    do n = 1, this%nbound
      node = this%nodelist(n)
      if (node < 1) cycle
      rhs(node) = rhs(node) + this%rhs(n)
      ipos = ia(node)
      call matrix_sln%add_value_pos(idxglo(ipos), this%hcof(n))
    end do
  end subroutine sfr_fc

  !> @ brief Add Newton-Raphson terms for package into solution.
  !!
  !!  Calculate and add the Newton-Raphson terms for the SFR package to the
  !!  coefficient matrix and right-hand side vector.
  !<
  subroutine sfr_fn(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    real(DP), dimension(:), intent(inout) :: rhs !< right-hand side vector for model
    integer(I4B), dimension(:), intent(in) :: ia !< solution CRS row pointers
    integer(I4B), dimension(:), intent(in) :: idxglo !< mapping vector for model (local) to solution (global)
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: ipos
    real(DP) :: rterm
    real(DP) :: drterm
    real(DP) :: rhs1
    real(DP) :: hcof1
    real(DP) :: q1
    real(DP) :: q2
    real(DP) :: hgwf
    !
    ! -- Copy package rhs and hcof into solution rhs and amat
    do j = 1, this%nbound
      i = this%isfrorder(j)
      ! -- skip inactive reaches
      if (this%iboundpak(i) < 1) cycle
      ! -- skip if reach is not connected to gwf
      n = this%nodelist(i)
      if (n < 1) cycle
      ipos = ia(n)
      rterm = this%hcof(i) * this%xnew(n)
      ! -- calculate perturbed head
      hgwf = this%xnew(n) + DEM4
      call this%sfr_solve(i, hgwf, hcof1, rhs1, update=.false.)
      q1 = rhs1 - hcof1 * hgwf
      ! -- calculate unperturbed head
      q2 = this%rhs(i) - this%hcof(i) * this%xnew(n)
      ! -- calculate derivative
      drterm = (q2 - q1) / DEM4
      ! -- add terms to convert conductance formulation into
      !    newton-raphson formulation
      call matrix_sln%add_value_pos(idxglo(ipos), drterm - this%hcof(i))
      rhs(n) = rhs(n) - rterm + drterm * this%xnew(n)
    end do
  end subroutine sfr_fn

  !> @ brief Convergence check for package.
  !!
  !!  Perform additional convergence checks on the flow between the SFR package
  !!  and the model it is attached to.
  !<
  subroutine sfr_cc(this, innertot, kiter, iend, icnvgmod, cpak, ipak, dpak)
    ! -- modules
    use TdisModule, only: totim, kstp, kper, delt
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    integer(I4B), intent(in) :: innertot !< total number of inner iterations
    integer(I4B), intent(in) :: kiter !< Picard iteration number
    integer(I4B), intent(in) :: iend !< flag indicating if this is the last Picard iteration
    integer(I4B), intent(in) :: icnvgmod !< flag inficating if the model has met specific convergence criteria
    character(len=LENPAKLOC), intent(inout) :: cpak !< string for user node
    integer(I4B), intent(inout) :: ipak !< location of the maximum dependent variable change
    real(DP), intent(inout) :: dpak !< maximum dependent variable change
    ! -- local
    character(len=LENPAKLOC) :: cloc
    character(len=LINELENGTH) :: tag
    integer(I4B) :: icheck
    integer(I4B) :: ipakfail
    integer(I4B) :: locdhmax
    integer(I4B) :: locrmax
    integer(I4B) :: locdqfrommvrmax
    integer(I4B) :: ntabrows
    integer(I4B) :: ntabcols
    integer(I4B) :: n
    real(DP) :: q
    real(DP) :: q0
    real(DP) :: qtolfact
    real(DP) :: dh
    real(DP) :: r
    real(DP) :: dhmax
    real(DP) :: rmax
    real(DP) :: dqfrommvr
    real(DP) :: dqfrommvrmax
    !
    ! -- initialize local variables
    icheck = this%iconvchk
    ipakfail = 0
    locdhmax = 0
    locrmax = 0
    r = DZERO
    dhmax = DZERO
    rmax = DZERO
    locdqfrommvrmax = 0
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
        ntabcols = 9
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
        tag = 'dinflowmax'
        call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
        tag = 'dinflowmax_loc'
        call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
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
      final_check: do n = 1, this%maxbound
        if (this%iboundpak(n) == 0) cycle
        !
        ! -- set the Q to length factor
        qtolfact = delt / this%calc_surface_area(n)
        !
        ! -- calculate stage change
        dh = this%stage0(n) - this%stage(n)
        !
        ! -- evaluate flow difference if the time step is transient
        if (this%gwfiss == 0) then
          r = this%usflow0(n) - this%usflow(n)
          !
          ! -- normalize flow difference and convert to a depth
          r = r * qtolfact
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
          locrmax = n
          rmax = r
          dqfrommvrmax = dqfrommvr
          locdqfrommvrmax = n
        else
          if (abs(dh) > abs(dhmax)) then
            locdhmax = n
            dhmax = dh
          end if
          if (abs(r) > abs(rmax)) then
            locrmax = n
            rmax = r
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
        write (cloc, "(a,'-',a)") trim(this%packName), 'stage'
        cpak = trim(cloc)
      end if
      if (ABS(rmax) > abs(dpak)) then
        ipak = locrmax
        dpak = rmax
        write (cloc, "(a,'-',a)") trim(this%packName), 'inflow'
        cpak = trim(cloc)
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
        call this%pakcsvtab%add_term(rmax)
        call this%pakcsvtab%add_term(locrmax)
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
  end subroutine sfr_cc

  !> @ brief Calculate package flows.
  !!
  !!  Calculate the flow between connected SFR package control volumes.
  !<
  subroutine sfr_cq(this, x, flowja, iadv)
    ! -- modules
    use BudgetModule, only: BudgetType
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    real(DP), dimension(:), intent(in) :: x !< current dependent-variable value
    real(DP), dimension(:), contiguous, intent(inout) :: flowja !< flow between two connected control volumes
    integer(I4B), optional, intent(in) :: iadv !< flag that indicates if this is an advance package
    ! -- local
    integer(I4B) :: i
    real(DP) :: qext
    ! -- for budget
    integer(I4B) :: n
    integer(I4B) :: n2
    real(DP) :: qoutflow
    real(DP) :: qfrommvr
    real(DP) :: qtomvr
    !
    ! -- call base functionality in bnd_cq.  This will calculate sfr-gwf flows
    !    and put them into this%simvals
    call this%BndType%bnd_cq(x, flowja, iadv=1)
    !
    ! -- Calculate qextoutflow and qoutflow for subsequent budgets
    do n = 1, this%maxbound
      !
      ! -- mover
      qfrommvr = DZERO
      qtomvr = DZERO
      if (this%imover == 1) then
        qfrommvr = this%pakmvrobj%get_qfrommvr(n)
        qtomvr = this%pakmvrobj%get_qtomvr(n)
        if (qtomvr > DZERO) then
          qtomvr = -qtomvr
        end if
      end if
      !
      ! -- external downstream stream flow
      qext = this%dsflow(n)
      qoutflow = DZERO
      if (qext > DZERO) then
        qext = -qext
      end if
      do i = this%ia(n) + 1, this%ia(n + 1) - 1
        if (this%idir(i) > 0) cycle
        n2 = this%ja(i)
        if (this%iboundpak(n2) == 0) cycle
        qext = DZERO
        exit
      end do
      !
      ! -- adjust external downstream stream flow using qtomvr
      if (qext < DZERO) then
        if (qtomvr < DZERO) then
          qext = qext - qtomvr
        end if
      else
        qoutflow = this%dsflow(n)
        if (qoutflow > DZERO) then
          qoutflow = -qoutflow
        end if
      end if
      !
      ! -- set qextoutflow and qoutflow for cell by cell budget
      !    output and observations
      this%qextoutflow(n) = qext
      this%qoutflow(n) = qoutflow
      !
    end do
    !
    ! -- fill the budget object
    call this%sfr_fill_budobj()
  end subroutine sfr_cq

  !> @ brief Output package flow terms.
  !!
  !!  Output SFR package flow terms.
  !<
  subroutine sfr_ot_package_flows(this, icbcfl, ibudfl)
    ! -- modules
    use TdisModule, only: kstp, kper, delt, pertim, totim
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: icbcfl !< flag and unit number for cell-by-cell output
    integer(I4B), intent(in) :: ibudfl !< flag indication if cell-by-cell data should be saved
    ! -- local
    integer(I4B) :: ibinun
    character(len=20), dimension(:), allocatable :: cellidstr
    integer(I4B) :: n
    integer(I4B) :: node
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
    ! -- Print sfr flows table
    if (ibudfl /= 0 .and. this%iprflow /= 0) then
      !
      ! -- If there are any 'none' gwf connections then need to calculate
      !    a vector of cellids and pass that in to the budget flow table because
      !    the table assumes that there are maxbound gwf entries, which is not
      !    the case if any 'none's are specified.
      if (this%ianynone > 0) then
        allocate (cellidstr(this%maxbound))
        do n = 1, this%maxbound
          node = this%igwfnode(n)
          if (node > 0) then
            call this%dis%noder_to_string(node, cellidstr(n))
          else
            cellidstr(n) = 'NONE'
          end if
        end do
        call this%budobj%write_flowtable(this%dis, kstp, kper, cellidstr)
        deallocate (cellidstr)
      else
        call this%budobj%write_flowtable(this%dis, kstp, kper)
      end if
    end if
  end subroutine sfr_ot_package_flows

  !> @ brief Output package dependent-variable terms.
  !!
  !!  Output SFR boundary package dependent-variable terms.
  !<
  subroutine sfr_ot_dv(this, idvsave, idvprint)
    ! -- modules
    use TdisModule, only: kstp, kper, pertim, totim
    use InputOutputModule, only: ulasav
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: idvsave !< flag and unit number for dependent-variable output
    integer(I4B), intent(in) :: idvprint !< flag indicating if dependent-variable should be written to the model listing file
    ! -- local
    character(len=20) :: cellid
    integer(I4B) :: ibinun
    integer(I4B) :: n
    integer(I4B) :: node
    real(DP) :: d
    real(DP) :: v
    real(DP) :: hgwf
    real(DP) :: sbot
    real(DP) :: depth
    real(DP) :: stage
    real(DP) :: w
    real(DP) :: cond
    real(DP) :: grad
    !
    ! -- set unit number for binary dependent variable output
    ibinun = 0
    if (this%istageout /= 0) then
      ibinun = this%istageout
    end if
    if (idvsave == 0) ibinun = 0
    !
    ! -- write sfr binary output
    if (ibinun > 0) then
      do n = 1, this%maxbound
        d = this%depth(n)
        v = this%stage(n)
        if (this%iboundpak(n) == 0) then
          v = DHNOFLO
        else if (d == DZERO) then
          v = DHDRY
        end if
        this%dbuff(n) = v
      end do
      call ulasav(this%dbuff, '           STAGE', kstp, kper, pertim, totim, &
                  this%maxbound, 1, 1, ibinun)
    end if
    !
    ! -- print sfr stage and depth table
    if (idvprint /= 0 .and. this%iprhed /= 0) then
      !
      ! -- set table kstp and kper
      call this%stagetab%set_kstpkper(kstp, kper)
      !
      ! -- fill stage data
      do n = 1, this%maxbound
        node = this%igwfnode(n)
        if (node > 0) then
          call this%dis%noder_to_string(node, cellid)
          hgwf = this%xnew(node)
        else
          cellid = 'NONE'
        end if
        if (this%inamedbound == 1) then
          call this%stagetab%add_term(this%boundname(n))
        end if
        call this%stagetab%add_term(n)
        call this%stagetab%add_term(cellid)
        if (this%iboundpak(n) /= 0) then
          depth = this%depth(n)
          stage = this%stage(n)
          w = this%calc_top_width_wet(n, depth)
          call this%sfr_calc_cond(n, depth, cond, stage, hgwf)
        else
          depth = DHNOFLO
          stage = DHNOFLO
          w = DHNOFLO
          cond = DHNOFLO
        end if
        if (depth == DZERO) then
          call this%stagetab%add_term(DHDRY)
        else
          call this%stagetab%add_term(stage)
        end if
        call this%stagetab%add_term(depth)
        call this%stagetab%add_term(w)
        if (node > 0) then
          if (this%iboundpak(n) /= 0) then
            sbot = this%strtop(n) - this%bthick(n)
            if (hgwf < sbot) then
              grad = stage - sbot
            else
              grad = stage - hgwf
            end if
            grad = grad / this%bthick(n)
          else
            grad = DHNOFLO
          end if
          call this%stagetab%add_term(hgwf)
          call this%stagetab%add_term(cond)
          call this%stagetab%add_term(grad)
        else
          call this%stagetab%add_term('--')
          call this%stagetab%add_term('--')
          call this%stagetab%add_term('--')
        end if
      end do
    end if
  end subroutine sfr_ot_dv

  !> @ brief Output advanced package budget summary.
  !!
  !!  Output SFR package budget summary.
  !<
  subroutine sfr_ot_bdsummary(this, kstp, kper, iout, ibudfl)
    ! -- module
    use TdisModule, only: totim, delt
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: kstp !< time step number
    integer(I4B), intent(in) :: kper !< period number
    integer(I4B), intent(in) :: iout !< flag and unit number for the model listing file
    integer(I4B), intent(in) :: ibudfl !< flag indicating budget should be written
    !
    call this%budobj%write_budtable(kstp, kper, iout, ibudfl, totim, delt)
  end subroutine sfr_ot_bdsummary

  !> @ brief Deallocate package memory
  !!
  !!  Deallocate SFR package scalars and arrays.
  !<
  subroutine sfr_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    !
    ! -- deallocate arrays
    call mem_deallocate(this%qoutflow)
    call mem_deallocate(this%qextoutflow)
    deallocate (this%csfrbudget)
    call mem_deallocate(this%sfrname, 'SFRNAME', this%memoryPath)
    call mem_deallocate(this%dbuff)
    deallocate (this%cauxcbc)
    call mem_deallocate(this%qauxcbc)
    call mem_deallocate(this%iboundpak)
    call mem_deallocate(this%igwfnode)
    call mem_deallocate(this%igwftopnode)
    call mem_deallocate(this%length)
    call mem_deallocate(this%width)
    call mem_deallocate(this%strtop)
    call mem_deallocate(this%bthick)
    call mem_deallocate(this%hk)
    call mem_deallocate(this%slope)
    call mem_deallocate(this%nconnreach)
    call mem_deallocate(this%ustrf)
    call mem_deallocate(this%ftotnd)
    call mem_deallocate(this%usflow)
    call mem_deallocate(this%dsflow)
    call mem_deallocate(this%depth)
    call mem_deallocate(this%stage)
    call mem_deallocate(this%gwflow)
    call mem_deallocate(this%simevap)
    call mem_deallocate(this%simrunoff)
    call mem_deallocate(this%stage0)
    call mem_deallocate(this%usflow0)
    call mem_deallocate(this%denseterms)
    call mem_deallocate(this%viscratios)
    !
    ! -- stage, usflow, and dsflow for previous timestep
    if (this%istorage == 1) then
      call mem_deallocate(this%stageold)
      call mem_deallocate(this%dsflowold)
      call mem_deallocate(this%storage)
      call mem_deallocate(this%usinflow)
      call mem_deallocate(this%usinflowold)
    end if
    !
    ! -- deallocate reach order and connection data
    call mem_deallocate(this%isfrorder)
    call mem_deallocate(this%ia)
    call mem_deallocate(this%ja)
    call mem_deallocate(this%idir)
    call mem_deallocate(this%idiv)
    call mem_deallocate(this%qconn)
    !
    ! -- deallocate boundary data
    call mem_deallocate(this%rough)
    call mem_deallocate(this%rain)
    call mem_deallocate(this%evap)
    call mem_deallocate(this%inflow)
    call mem_deallocate(this%runoff)
    call mem_deallocate(this%sstage)
    !
    ! -- deallocate aux variables
    call mem_deallocate(this%rauxvar)
    !
    ! -- deallocate diversion variables
    call mem_deallocate(this%iadiv)
    call mem_deallocate(this%divreach)
    if (associated(this%divcprior)) then
      deallocate (this%divcprior)
    end if
    call mem_deallocate(this%divflow)
    call mem_deallocate(this%divq)
    call mem_deallocate(this%ndiv)
    !
    ! -- deallocate cross-section data
    call mem_deallocate(this%ncrosspts)
    call mem_deallocate(this%iacross)
    call mem_deallocate(this%station)
    call mem_deallocate(this%xsheight)
    call mem_deallocate(this%xsrough)
    !
    ! -- deallocate budobj
    call this%budobj%budgetobject_da()
    deallocate (this%budobj)
    nullify (this%budobj)
    !
    ! -- deallocate stage table
    if (this%iprhed > 0) then
      call this%stagetab%table_da()
      deallocate (this%stagetab)
      nullify (this%stagetab)
    end if
    !
    ! -- deallocate package csv table
    if (this%ipakcsv > 0) then
      if (associated(this%pakcsvtab)) then
        call this%pakcsvtab%table_da()
        deallocate (this%pakcsvtab)
        nullify (this%pakcsvtab)
      end if
    end if
    !
    ! -- deallocate scalars
    call mem_deallocate(this%istorage)
    call mem_deallocate(this%storage_weight)
    call mem_deallocate(this%iprhed)
    call mem_deallocate(this%istageout)
    call mem_deallocate(this%ibudgetout)
    call mem_deallocate(this%ibudcsv)
    call mem_deallocate(this%ipakcsv)
    call mem_deallocate(this%idiversions)
    call mem_deallocate(this%maxsfrpicard)
    call mem_deallocate(this%maxsfrit)
    call mem_deallocate(this%bditems)
    call mem_deallocate(this%cbcauxitems)
    call mem_deallocate(this%unitconv)
    call mem_deallocate(this%lengthconv)
    call mem_deallocate(this%timeconv)
    call mem_deallocate(this%dmaxchg)
    call mem_deallocate(this%deps)
    call mem_deallocate(this%nconn)
    call mem_deallocate(this%icheck)
    call mem_deallocate(this%iconvchk)
    call mem_deallocate(this%idense)
    call mem_deallocate(this%ianynone)
    call mem_deallocate(this%ncrossptstot)
    nullify (this%gwfiss)
    !
    ! -- call base BndType deallocate
    call this%BndType%bnd_da()
  end subroutine sfr_da

  !> @ brief Define the list label for the package
  !!
  !!  Method defined the list label for the SFR package. The list label is
  !!  the heading that is written to iout when PRINT_INPUT option is used.
  !<
  subroutine define_listlabel(this)
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
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

  !
  ! -- Procedures related to observations (type-bound)

  !> @brief Determine if observations are supported.
  !!
  !! Function to determine if observations are supported by the SFR package.
  !! Observations are supported by the SFR package.
  !!
  !! @return  sfr_obs_supported       boolean indicating if observations are supported
  !<
  logical function sfr_obs_supported(this)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    !
    ! -- set boolean
    sfr_obs_supported = .true.
  end function sfr_obs_supported

  !> @brief Define the observation types available in the package
  !!
  !! Method to define the observation types available in the SFR package.
  !<
  subroutine sfr_df_obs(this)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    ! -- local
    integer(I4B) :: indx
    !
    ! -- Store obs type and assign procedure pointer
    !    for stage observation type.
    call this%obs%StoreObsType('stage', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for inflow observation type.
    call this%obs%StoreObsType('inflow', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for inflow observation type.
    call this%obs%StoreObsType('ext-inflow', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for rainfall observation type.
    call this%obs%StoreObsType('rainfall', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for runoff observation type.
    call this%obs%StoreObsType('runoff', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for evaporation observation type.
    call this%obs%StoreObsType('evaporation', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for outflow observation type.
    call this%obs%StoreObsType('outflow', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for ext-outflow observation type.
    call this%obs%StoreObsType('ext-outflow', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for to-mvr observation type.
    call this%obs%StoreObsType('to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for sfr-frommvr observation type.
    call this%obs%StoreObsType('from-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for sfr observation type.
    call this%obs%StoreObsType('sfr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for upstream flow observation type.
    call this%obs%StoreObsType('upstream-flow', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for downstream flow observation type.
    call this%obs%StoreObsType('downstream-flow', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for depth observation type.
    call this%obs%StoreObsType('depth', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for wetted-perimeter observation type.
    call this%obs%StoreObsType('wet-perimeter', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for wetted-area observation type.
    call this%obs%StoreObsType('wet-area', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for wetted-width observation type.
    call this%obs%StoreObsType('wet-width', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => sfr_process_obsID
  end subroutine sfr_df_obs

  !> @brief Save observations for the package
  !!
  !! Method to save simulated values for the SFR package.
  !<
  subroutine sfr_bd_obs(this)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: n
    real(DP) :: v
    character(len=100) :: msg
    type(ObserveType), pointer :: obsrv => null()
    !
    ! Write simulated values for all sfr observations
    if (this%obs%npakobs > 0) then
      call this%obs%obs_bd_clear()
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv
        do j = 1, obsrv%indxbnds_count
          n = obsrv%indxbnds(j)
          v = DZERO
          select case (obsrv%ObsTypeId)
          case ('STAGE')
            v = this%stage(n)
          case ('TO-MVR')
            v = DNODATA
            if (this%imover == 1) then
              v = this%pakmvrobj%get_qtomvr(n)
              if (v > DZERO) then
                v = -v
              end if
            end if
          case ('FROM-MVR')
            v = DNODATA
            if (this%imover == 1) then
              v = this%pakmvrobj%get_qfrommvr(n)
            end if
          case ('EXT-INFLOW')
            v = this%inflow(n)
          case ('INFLOW')
            v = this%usflow(n)
          case ('OUTFLOW')
            v = this%qoutflow(n)
          case ('EXT-OUTFLOW')
            v = this%qextoutflow(n)
          case ('RAINFALL')
            if (this%iboundpak(n) /= 0) then
              v = this%rain(n)
            else
              v = DZERO
            end if
          case ('RUNOFF')
            v = this%simrunoff(n)
          case ('EVAPORATION')
            v = this%simevap(n)
          case ('SFR')
            v = this%gwflow(n)
          case ('UPSTREAM-FLOW')
            v = this%usflow(n)
            if (this%imover == 1) then
              v = v + this%pakmvrobj%get_qfrommvr(n)
            end if
          case ('DOWNSTREAM-FLOW')
            v = this%dsflow(n)
            if (v > DZERO) then
              v = -v
            end if
          case ('DEPTH')
            v = this%depth(n)
          case ('WET-PERIMETER')
            v = this%calc_perimeter_wet(n, this%depth(n))
          case ('WET-AREA')
            v = this%calc_area_wet(n, this%depth(n))
          case ('WET-WIDTH')
            v = this%calc_top_width_wet(n, this%depth(n))
          case default
            msg = 'Unrecognized observation type: '//trim(obsrv%ObsTypeId)
            call store_error(msg)
          end select
          call this%obs%SaveOneSimval(obsrv, v)
        end do
      end do
      !
      ! -- write summary of package error messages
      if (count_errors() > 0) then
        call this%parser%StoreErrorUnit()
      end if
    end if
  end subroutine sfr_bd_obs

  !> @brief Read and prepare observations for a package
  !!
  !! Method to read and prepare observations for a SFR package.
  !<
  subroutine sfr_rp_obs(this)
    ! -- modules
    use TdisModule, only: kper
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: nn1
    character(len=LENBOUNDNAME) :: bname
    logical(LGP) :: jfound
    class(ObserveType), pointer :: obsrv => null()
    ! -- formats
10  format('Boundary "', a, '" for observation "', a, &
           '" is invalid in package "', a, '"')
30  format('Boundary name not provided for observation "', a, &
           '" in package "', a, '"')
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
            ! -- Observation location(s) is(are) based on a boundary name.
            !    Iterate through all boundaries to identify and store
            !    corresponding index(indices) in bound array.
            jfound = .false.
            do j = 1, this%maxbound
              if (this%boundname(j) == bname) then
                jfound = .true.
                call obsrv%AddObsIndex(j)
              end if
            end do
            if (.not. jfound) then
              write (errmsg, 10) &
                trim(bname), trim(obsrv%name), trim(this%packName)
              call store_error(errmsg)
            end if
          else
            write (errmsg, 30) trim(obsrv%name), trim(this%packName)
            call store_error(errmsg)
          end if
        else if (nn1 < 1 .or. nn1 > this%maxbound) then
          write (errmsg, '(a,1x,a,1x,i0,1x,a,1x,i0,a)') &
            trim(adjustl(obsrv%ObsTypeId)), &
            'reach must be greater than 0 and less than or equal to', &
            this%maxbound, '(specified value is ', nn1, ')'
          call store_error(errmsg)
        else
          if (obsrv%indxbnds_count == 0) then
            call obsrv%AddObsIndex(nn1)
          else
            errmsg = 'Programming error in sfr_rp_obs'
            call store_error(errmsg)
          end if
        end if
        !
        ! -- catch non-cumulative observation assigned to observation defined
        !    by a boundname that is assigned to more than one element
        if (obsrv%ObsTypeId == 'STAGE' .or. &
            obsrv%ObsTypeId == 'DEPTH' .or. &
            obsrv%ObsTypeId == 'WET-PERIMETER' .or. &
            obsrv%ObsTypeId == 'WET-AREA' .or. &
            obsrv%ObsTypeId == 'WET-WIDTH') then
          nn1 = obsrv%NodeNumber
          if (nn1 == NAMEDBOUNDFLAG) then
            if (obsrv%indxbnds_count > 1) then
              write (errmsg, '(a,3(1x,a))') &
                trim(adjustl(obsrv%ObsTypeId)), &
                'for observation', trim(adjustl(obsrv%Name)), &
                ' must be assigned to a reach with a unique boundname.'
              call store_error(errmsg)
            end if
          end if
        end if
        !
        ! -- check that node number 1 is valid; call store_error if not
        do j = 1, obsrv%indxbnds_count
          nn1 = obsrv%indxbnds(j)
          if (nn1 < 1 .or. nn1 > this%maxbound) then
            write (errmsg, '(a,1x,a,1x,i0,1x,a,1x,i0,a)') &
              trim(adjustl(obsrv%ObsTypeId)), &
              'reach must be greater than 0 and less than or equal to', &
              this%maxbound, '(specified value is ', nn1, ')'
            call store_error(errmsg)
          end if
        end do
      end do
      !
      ! -- evaluate if there are any observation errors
      if (count_errors() > 0) then
        call this%parser%StoreErrorUnit()
      end if
    end if
  end subroutine sfr_rp_obs

  !
  ! -- Procedures related to observations (NOT type-bound)

  !> @brief Process observation IDs for a package
  !!
  !! Method to process observation ID strings for a SFR package.
  !<
  subroutine sfr_process_obsID(obsrv, dis, inunitobs, iout)
    ! -- dummy
    type(ObserveType), intent(inout) :: obsrv !< Observation object
    class(DisBaseType), intent(in) :: dis !< Discretization object
    integer(I4B), intent(in) :: inunitobs !< file unit number for the package observation file
    integer(I4B), intent(in) :: iout !< model listing file unit number
    ! -- local
    integer(I4B) :: nn1
    integer(I4B) :: icol
    integer(I4B) :: istart
    integer(I4B) :: istop
    character(len=LINELENGTH) :: string
    character(len=LENBOUNDNAME) :: bndname
    !
    ! -- initialize local variables
    string = obsrv%IDstring
    !
    ! -- Extract reach number from string and store it.
    !    If 1st item is not an integer(I4B), it should be a
    !    boundary name--deal with it.
    icol = 1
    !
    ! -- get reach number or boundary name
    call extract_idnum_or_bndname(string, icol, istart, istop, nn1, bndname)
    if (nn1 == NAMEDBOUNDFLAG) then
      obsrv%FeatureName = bndname
    end if
    !
    ! -- store reach number (NodeNumber)
    obsrv%NodeNumber = nn1
  end subroutine sfr_process_obsID

  !
  ! -- private sfr methods
  !

  !> @brief Set period data
  !!
  !! Method to read and set period data for a SFR package reach.
  !<
  subroutine sfr_set_stressperiod(this, n, ichkustrm, crossfile)
    ! -- modules
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    integer(I4B), intent(inout) :: ichkustrm !< flag indicating if upstream fraction data specified
    character(len=LINELENGTH), intent(inout) :: crossfile !< cross-section file name
    ! -- local
    character(len=10) :: cnum
    character(len=LINELENGTH) :: text
    character(len=LINELENGTH) :: caux
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ival
    integer(I4B) :: ii
    integer(I4B) :: jj
    integer(I4B) :: idiv
    integer(I4B) :: ixserror
    character(len=10) :: cp
    real(DP) :: divq
    real(DP), pointer :: bndElem => null()
    !
    ! -- initialize variables
    crossfile = 'NONE'
    !
    ! -- read line
    call this%parser%GetStringCaps(keyword)
    select case (keyword)
    case ('STATUS')
      ichkustrm = 1
      call this%parser%GetStringCaps(text)
      if (text == 'INACTIVE') then
        this%iboundpak(n) = 0
      else if (text == 'ACTIVE') then
        this%iboundpak(n) = 1
      else if (text == 'SIMPLE') then
        this%iboundpak(n) = -1
      else
        write (errmsg, '(2a)') &
          'Unknown '//trim(this%text)//' sfr status keyword: ', trim(text)
        call store_error(errmsg)
      end if
    case ('BEDK')
      call this%parser%GetString(text)
      jj = 1 ! For 'BEDK'
      bndElem => this%hk(n)
      call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                         this%packName, 'BND', &
                                         this%tsManager, this%iprpak, &
                                         'BEDK')
    case ('MANNING')
      call this%parser%GetString(text)
      jj = 1 ! For 'MANNING'
      bndElem => this%rough(n)
      call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                         this%packName, 'BND', &
                                         this%tsManager, this%iprpak, &
                                         'MANNING')
    case ('STAGE')
      call this%parser%GetString(text)
      jj = 1 ! For 'STAGE'
      bndElem => this%sstage(n)
      call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                         this%packName, 'BND', &
                                         this%tsManager, this%iprpak, 'STAGE')
    case ('RAINFALL')
      call this%parser%GetString(text)
      jj = 1 ! For 'RAIN'
      bndElem => this%rain(n)
      call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                         this%packName, 'BND', &
                                         this%tsManager, this%iprpak, 'RAIN')
    case ('EVAPORATION')
      call this%parser%GetString(text)
      jj = 1 ! For 'EVAP'
      bndElem => this%evap(n)
      call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                         this%packName, 'BND', &
                                         this%tsManager, this%iprpak, &
                                         'EVAP')
    case ('RUNOFF')
      call this%parser%GetString(text)
      jj = 1 ! For 'RUNOFF'
      bndElem => this%runoff(n)
      call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                         this%packName, 'BND', &
                                         this%tsManager, this%iprpak, &
                                         'RUNOFF')
    case ('INFLOW')
      call this%parser%GetString(text)
      jj = 1 ! For 'INFLOW'
      bndElem => this%inflow(n)
      call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                         this%packName, 'BND', &
                                         this%tsManager, this%iprpak, &
                                         'INFLOW')
    case ('DIVERSION')
      !
      ! -- make sure reach has at least one diversion
      if (this%ndiv(n) < 1) then
        write (cnum, '(i0)') n
        errmsg = 'diversions cannot be specified for reach '//trim(cnum)
        call store_error(errmsg)
      end if
      !
      ! -- read diversion number
      ival = this%parser%GetInteger()
      if (ival < 1 .or. ival > this%ndiv(n)) then
        write (cnum, '(i0)') n
        errmsg = 'Reach  '//trim(cnum)
        write (cnum, '(i0)') this%ndiv(n)
        errmsg = trim(errmsg)//' diversion number should be between 1 '// &
                 'and '//trim(cnum)//'.'
        call store_error(errmsg)
      end if
      idiv = ival
      !
      ! -- read value
      call this%parser%GetString(text)
      ii = this%iadiv(n) + idiv - 1
      jj = 1 ! For 'DIVERSION'
      bndElem => this%divflow(ii)
      call read_value_or_time_series_adv(text, ii, jj, bndElem, &
                                         this%packName, 'BND', &
                                         this%tsManager, this%iprpak, &
                                         'DIVFLOW')
      !
      ! -- if diversion cprior is 'fraction', ensure that 0.0 <= fraction <= 1.0
      cp = this%divcprior(ii)
      divq = this%divflow(ii)
      if (cp == 'FRACTION' .and. (divq < DZERO .or. divq > DONE)) then
        write (errmsg, '(a,1x,i0,a)') &
          'cprior is type FRACTION for diversion no.', ii, &
          ', but divflow not within the range 0.0 to 1.0'
        call store_error(errmsg)
      end if
    case ('UPSTREAM_FRACTION')
      ichkustrm = 1
      call this%parser%GetString(text)
      jj = 1 ! For 'USTRF'
      bndElem => this%ustrf(n)
      call read_value_or_time_series_adv(text, n, jj, bndElem, &
                                         this%packName, 'BND', &
                                         this%tsManager, this%iprpak, 'USTRF')

    case ('CROSS_SECTION')
      ixserror = 0
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
          ixserror = 1
        end if
        if (ixserror == 0) then
          call this%parser%GetString(crossfile)
        end if
      case default
        write (errmsg, '(a,1x,i4,1x,a)') &
          'CROSS-SECTION TABLE ENTRY for REACH ', n, &
          'MUST INCLUDE TAB6 KEYWORD'
        call store_error(errmsg)
      end select

    case ('AUXILIARY')
      call this%parser%GetStringCaps(caux)
      do jj = 1, this%naux
        if (trim(adjustl(caux)) /= trim(adjustl(this%auxname(jj)))) cycle
        call this%parser%GetString(text)
        ii = n
        bndElem => this%rauxvar(jj, ii)
        call read_value_or_time_series_adv(text, ii, jj, bndElem, &
                                           this%packName, 'AUX', &
                                           this%tsManager, this%iprpak, &
                                           this%auxname(jj))
        exit
      end do

    case default
      write (errmsg, '(a,a)') &
        'Unknown '//trim(this%text)//' sfr data keyword: ', &
        trim(keyword)//'.'
      call store_error(errmsg)
    end select
  end subroutine sfr_set_stressperiod

  !> @brief Solve reach continuity equation
  !!
  !! Method to solve the continuity equation for a SFR package reach.
  !<
  subroutine sfr_solve(this, n, h, hcof, rhs, update)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: h !< groundwater head in cell connected to reach
    real(DP), intent(inout) :: hcof !< coefficient term added to the diagonal
    real(DP), intent(inout) :: rhs !< right-hand side term
    logical(LGP), intent(in), optional :: update !< boolean indicating if the reach depth and stage variables should be updated to current iterate
    ! -- local
    logical(LGP) :: lupdate
    integer(I4B) :: i
    integer(I4B) :: ii
    integer(I4B) :: n2
    real(DP) :: hgwf
    real(DP) :: sa
    real(DP) :: sa_wet
    real(DP) :: qu
    real(DP) :: qi
    real(DP) :: qr
    real(DP) :: qe
    real(DP) :: qro
    real(DP) :: qsrc
    real(DP) :: qfrommvr
    real(DP) :: qgwf
    real(DP) :: tp
    real(DP) :: bt
    real(DP) :: hsfr
    real(DP) :: qd
    real(DP) :: d1
    real(DP) :: sumleak
    real(DP) :: sumrch
    real(DP) :: gwfhcof
    real(DP) :: gwfrhs
    !
    ! -- Process optional dummy variables
    if (present(update)) then
      lupdate = update
    else
      lupdate = .true.
    end if

    ! -- initialize variables
    hcof = DZERO
    rhs = DZERO
    !
    if (this%iboundpak(n) == 0) then
      this%depth(n) = DZERO
      this%stage(n) = DHNOFLO
      this%usflow(n) = DZERO
      this%simevap(n) = DZERO
      this%simrunoff(n) = DZERO
      this%dsflow(n) = DZERO
      this%gwflow(n) = DZERO
    else
      hgwf = h
      d1 = DZERO
      qsrc = DZERO
      qgwf = DZERO
      qd = this%dsflow(n)

      ! -- calculate initial depth assuming a wide cross-section and
      !    ignore groundwater leakage
      ! -- calculate upstream flow
      qu = DZERO
      do i = this%ia(n) + 1, this%ia(n + 1) - 1
        if (this%idir(i) < 0) cycle
        n2 = this%ja(i)
        do ii = this%ia(n2) + 1, this%ia(n2 + 1) - 1
          if (this%idir(ii) > 0) cycle
          if (this%ja(ii) /= n) cycle
          qu = qu + this%qconn(ii)
        end do
      end do
      this%usflow(n) = qu

      ! -- calculate remaining terms
      sa = this%calc_surface_area(n)
      sa_wet = this%calc_surface_area_wet(n, this%depth(n))
      qi = this%inflow(n)
      qr = this%rain(n) * sa
      qe = this%evap(n) * sa_wet
      qro = this%runoff(n)

      ! -- Water mover term; assume that it goes in at the upstream end of the reach
      qfrommvr = DZERO
      if (this%imover == 1) then
        qfrommvr = this%pakmvrobj%get_qfrommvr(n)
      end if

      ! -- calculate downstream flow ignoring groundwater leakage
      qsrc = qu + qi + qr - qe + qro + qfrommvr

      ! -- adjust runoff or evaporation if sum of sources is negative
      call this%sfr_adjust_ro_ev(qsrc, qu, qi, qr, qro, qe, qfrommvr)

      ! -- set simulated evaporation and runoff
      this%simevap(n) = qe
      this%simrunoff(n) = qro

      ! -- calculate reach flow using appropriate method
      if (this%iboundpak(n) < 0) then
        call this%sfr_calc_constant(n, d1, hgwf, qgwf, qd)
      else
        if (this%gwfiss == 0 .and. this%istorage == 1) then
          call this%sfr_calc_transient(n, d1, hgwf, qu, qi, &
                                       qfrommvr, qr, qe, qro, &
                                       qgwf, qd)
        else
          call this%sfr_calc_steady(n, d1, hgwf, qu, qi, &
                                    qfrommvr, qr, qe, qro, &
                                    qgwf, qd)
        end if
      end if

      ! -- update sfr stage
      tp = this%strtop(n)
      bt = tp - this%bthick(n)
      hsfr = tp + d1

      ! -- update stored values
      if (lupdate) then
        ! -- save depth and calculate stage
        this%depth(n) = d1
        this%stage(n) = hsfr
        ! -- update flows
        call this%sfr_update_flows(n, qd, qgwf)
      end if

      ! -- calculate sumleak and sumrch
      sumleak = DZERO
      sumrch = DZERO
      if (this%gwfiss == 0) then
        sumleak = qgwf
      else
        sumleak = qgwf
      end if
      if (hgwf < bt) then
        sumrch = qgwf
      end if

      ! -- make final qgwf calculation and obtain
      !    gwfhcof and gwfrhs values
      call this%sfr_calc_qgwf(n, d1, hgwf, qgwf, gwfhcof, gwfrhs)

      ! -- update hcof and rhs terms
      if (abs(sumleak) > DZERO) then
        ! -- stream leakage is not head dependent
        if (hgwf < bt) then
          rhs = rhs - sumrch
          ! -- stream leakage is head dependent
        else if ((sumleak - qsrc) < -DEM30) then
          if (this%gwfiss == 0) then
            rhs = rhs + gwfrhs - sumrch
          else
            rhs = rhs + gwfrhs
          end if
          hcof = gwfhcof
          ! -- place holder for UZF
        else
          if (this%gwfiss == 0) then
            rhs = rhs - sumleak - sumrch
          else
            rhs = rhs - sumleak
          end if
        end if

        ! -- add groundwater leakage
      else if (hgwf < bt) then
        rhs = rhs - sumrch
      end if
    end if
  end subroutine sfr_solve

  !> @brief Update flow terms
  !!
  !! Method to update downstream flow and groundwater leakage terms for
  !! a SFR package reach.
  !<
  subroutine sfr_update_flows(this, n, qd, qgwf)
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(inout) :: qd !< downstream reach flow
    real(DP), intent(in) :: qgwf !< groundwater leakage for reach
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: n2
    integer(I4B) :: idiv
    integer(I4B) :: jpos
    real(DP) :: qdiv
    real(DP) :: f
    !
    ! -- update reach terms
    !
    ! -- save final downstream stream flow
    this%dsflow(n) = qd
    !
    ! -- save groundwater leakage
    this%gwflow(n) = qgwf
    !
    ! -- route downstream flow
    if (qd > DZERO) then
      !
      ! -- route water to diversions
      do i = this%ia(n) + 1, this%ia(n + 1) - 1
        if (this%idir(i) > 0) cycle
        idiv = this%idiv(i)
        if (idiv == 0) cycle
        jpos = this%iadiv(n) + idiv - 1
        call this%sfr_calc_div(n, idiv, qd, qdiv)
        this%qconn(i) = qdiv
        this%divq(jpos) = qdiv
      end do
      !
      ! -- Mover terms: store outflow after diversion loss
      !    as qformvr and reduce outflow (qd)
      !    by how much was actually sent to the mover
      if (this%imover == 1) then
        call this%pakmvrobj%accumulate_qformvr(n, qd)
        qd = MAX(qd - this%pakmvrobj%get_qtomvr(n), DZERO)
      end if
      !
      ! -- route remaining water to downstream reaches
      do i = this%ia(n) + 1, this%ia(n + 1) - 1
        if (this%idir(i) > 0) cycle
        if (this%idiv(i) > 0) cycle
        n2 = this%ja(i)
        if (this%iboundpak(n2) == 0) cycle
        f = this%ustrf(n2) / this%ftotnd(n)
        this%qconn(i) = qd * f
      end do
    else
      do i = this%ia(n) + 1, this%ia(n + 1) - 1
        if (this%idir(i) > 0) cycle
        this%qconn(i) = DZERO
        idiv = this%idiv(i)
        if (idiv == 0) cycle
        jpos = this%iadiv(n) + idiv - 1
        this%divq(jpos) = DZERO
      end do
    end if
  end subroutine sfr_update_flows

  !> @brief Adjust runoff and evaporation
  !!
  !! Method to adjust runoff and evaporation for a SFR package reach
  !! based on the total reach flow.
  !<
  subroutine sfr_adjust_ro_ev(this, qc, qu, qi, qr, qro, qe, qfrommvr)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    real(DP), intent(inout) :: qc !< total reach volumetric flow
    real(DP), intent(in) :: qu !< upstream reach volumetric flow
    real(DP), intent(in) :: qi !< reach volumetric inflow
    real(DP), intent(in) :: qr !< reach volumetric rainfall
    real(DP), intent(inout) :: qro !< reach volumetric runoff
    real(DP), intent(inout) :: qe !< reach volumetric evaporation
    real(DP), intent(in) :: qfrommvr !< reach volumetric flow from mover
    ! -- local
    real(DP) :: qt
    !
    ! -- adjust runoff or evaporation if sum of sources is negative
    if (qc < DZERO) then
      !
      ! -- calculate sources without evaporation
      qt = qu + qi + qr + qro + qfrommvr
      !
      ! -- runoff exceeds sources of water for reach
      if (qt < DZERO) then
        if (qro < DZERO) then
          qro = -(qu + qi + qr + qfrommvr)
          qe = DZERO
        end if
        !
        ! -- evaporation exceeds sources of water for reach
      else
        if (qe > DZERO) then
          qe = qu + qi + qr + qro + qfrommvr
        end if
      end if
      qc = qu + qi + qr - qe + qro + qfrommvr
    end if
  end subroutine sfr_adjust_ro_ev

  !> @brief Calculate downstream flow term
  !!
  !! Method to calculate downstream flow for a SFR package reach.
  !<
  subroutine sfr_calc_qd(this, n, depth, hgwf, qgwf, qd)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: depth !< reach depth
    real(DP), intent(in) :: hgwf !< groundwater head in connected GWF cell
    real(DP), intent(inout) :: qgwf !< groundwater leakage for reach
    real(DP), intent(inout) :: qd !< residual
    ! -- local
    real(DP) :: qsrc
    !
    ! -- initialize residual
    qd = DZERO
    !
    ! -- calculate total water sources excluding groundwater leakage
    call this%sfr_calc_qsource(n, depth, qsrc)
    !
    ! -- estimate groundwater leakage
    call this%sfr_calc_qgwf(n, depth, hgwf, qgwf)
    if (-qgwf > qsrc) qgwf = -qsrc
    !
    ! -- calculate down stream flow
    qd = qsrc + qgwf
    !
    ! -- limit downstream flow to a positive value
    if (qd < DEM30) qd = DZERO
  end subroutine sfr_calc_qd

  !> @brief Calculate sum of sources
  !!
  !! Method to calculate the sum of sources for reach, excluding
  !! reach leakage, for a SFR package reach.
  !<
  subroutine sfr_calc_qsource(this, n, depth, qsrc)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: depth !< reach depth
    real(DP), intent(inout) :: qsrc !< sum of sources for reach
    ! -- local
    real(DP) :: qu
    real(DP) :: qi
    real(DP) :: qr
    real(DP) :: qe
    real(DP) :: qro
    real(DP) :: qfrommvr
    real(DP) :: a
    real(DP) :: ae
    !
    ! -- initialize residual
    qsrc = DZERO
    !
    ! -- calculate flow terms
    qu = this%usflow(n)
    qi = this%inflow(n)
    qro = this%runoff(n)
    !
    ! -- calculate rainfall and evap
    a = this%calc_surface_area(n)
    ae = this%calc_surface_area_wet(n, depth)
    qr = this%rain(n) * a
    qe = this%evap(n) * ae
    !
    ! -- calculate mover term
    qfrommvr = DZERO
    if (this%imover == 1) then
      qfrommvr = this%pakmvrobj%get_qfrommvr(n)
    end if
    !
    ! -- calculate down stream flow
    qsrc = qu + qi + qr - qe + qro + qfrommvr
    !
    ! -- adjust runoff or evaporation if sum of sources is negative
    call this%sfr_adjust_ro_ev(qsrc, qu, qi, qr, qro, qe, qfrommvr)
  end subroutine sfr_calc_qsource

  !> @brief Calculate streamflow
  !!
  !! Method to calculate the streamflow using Manning's equation for a
  !! SFR package reach.
  !<
  subroutine sfr_calc_qman(this, n, depth, qman)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: depth !< reach depth
    real(DP), intent(inout) :: qman !< streamflow
    ! -- local
    integer(I4B) :: npts
    integer(I4B) :: i0
    integer(I4B) :: i1
    real(DP) :: sat
    real(DP) :: derv
    real(DP) :: s
    real(DP) :: r
    real(DP) :: aw
    real(DP) :: wp
    real(DP) :: rh
    !
    ! -- initialize variables
    qman = DZERO
    !
    ! -- calculate Manning's discharge for non-zero depths
    if (depth > DZERO) then
      npts = this%ncrosspts(n)
      !
      ! -- set constant terms for Manning's equation
      call sChSmooth(depth, sat, derv)
      s = this%slope(n)
      !
      ! -- calculate the mannings coefficient that is a
      !    function of depth
      if (npts > 1) then
        !
        ! -- get the location of the cross-section data for the reach
        i0 = this%iacross(n)
        i1 = this%iacross(n + 1) - 1
        !
        ! -- get the Manning's sum of the Manning's discharge
        !    for each section
        qman = get_mannings_section(npts, &
                                    this%station(i0:i1), &
                                    this%xsheight(i0:i1), &
                                    this%xsrough(i0:i1), &
                                    this%rough(n), &
                                    this%unitconv, &
                                    s, &
                                    depth)
      else
        r = this%rough(n)
        aw = this%calc_area_wet(n, depth)
        wp = this%calc_perimeter_wet(n, depth)
        if (wp > DZERO) then
          rh = aw / wp
        else
          rh = DZERO
        end if
        qman = this%unitconv * aw * (rh**DTWOTHIRDS) * sqrt(s) / r
      end if
      !
      ! -- calculate stream flow
      qman = sat * qman
    end if
  end subroutine sfr_calc_qman

  !> @brief Calculate reach-aquifer exchange
  !!
  !! Method to calculate the reach-aquifer exchange for a SFR package reach.
  !! The reach-aquifer exchange is relative to the reach. Calculated flow
  !! is positive if flow is from the aquifer to the reach.
  !<
  subroutine sfr_calc_qgwf(this, n, depth, hgwf, qgwf, gwfhcof, gwfrhs)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: depth !< reach depth
    real(DP), intent(in) :: hgwf !< head in GWF cell connected to reach
    real(DP), intent(inout) :: qgwf !< reach-aquifer exchange
    real(DP), intent(inout), optional :: gwfhcof !< diagonal coefficient term for reach
    real(DP), intent(inout), optional :: gwfrhs !< right-hand side term for reach
    ! -- local
    integer(I4B) :: node
    real(DP) :: tp
    real(DP) :: bt
    real(DP) :: hsfr
    real(DP) :: h_temp
    real(DP) :: cond
    real(DP) :: sat
    real(DP) :: derv
    real(DP) :: gwfhcof0
    real(DP) :: gwfrhs0
    !
    ! -- initialize qgwf
    qgwf = DZERO
    !
    ! -- skip sfr-aquifer exchange in external cells
    node = this%igwfnode(n)
    if (node < 1) return
    !
    ! -- skip sfr-aquifer exchange in inactive cells
    if (this%ibound(node) == 0) return
    !
    ! -- calculate saturation
    call sChSmooth(depth, sat, derv)
    !
    ! -- terms for calculating direction of gradient across streambed
    tp = this%strtop(n)
    bt = tp - this%bthick(n)
    hsfr = tp + depth
    h_temp = hgwf
    if (h_temp < bt) then
      h_temp = bt
    end if
    !
    ! -- calculate conductance
    call this%sfr_calc_cond(n, depth, cond, hsfr, h_temp)
    !
    ! -- calculate groundwater leakage
    qgwf = sat * cond * (h_temp - hsfr)
    gwfrhs0 = -sat * cond * hsfr
    gwfhcof0 = -sat * cond
    !
    ! Add density contributions, if active
    if (this%idense /= 0) then
      call this%sfr_calculate_density_exchange(n, hsfr, hgwf, cond, tp, &
                                               qgwf, gwfhcof0, gwfrhs0)
    end if
    !
    ! -- Set gwfhcof and gwfrhs if present
    if (present(gwfhcof)) gwfhcof = gwfhcof0
    if (present(gwfrhs)) gwfrhs = gwfrhs0
  end subroutine sfr_calc_qgwf

  !> @brief Determine if a reach is connected to a gwf cell
  !!
  !! Function to determine if a reach is connected to a gwf cell. If connected,
  !! the return value is 1. Otherwise, the return value is 0.
  !<
  function sfr_gwf_conn(this, n)
    ! -- return variable
    integer(I4B) :: sfr_gwf_conn !< flag indicating if reach is connected to a gwf cell
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    ! -- local
    integer(I4B) :: node

    sfr_gwf_conn = 0
    node = this%igwfnode(n)
    if (node > 0 .and. this%hk(n) > DZERO) then
      sfr_gwf_conn = 1
    end if
  end function sfr_gwf_conn

  !> @brief Calculate reach-aquifer conductance
  !!
  !! Method to calculate the reach-aquifer conductance for a SFR package reach.
  !<
  subroutine sfr_calc_cond(this, n, depth, cond, hsfr, h_temp)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: depth !< reach depth
    real(DP), intent(inout) :: cond !< reach-aquifer conductance
    real(DP), intent(in), optional :: hsfr !< stream stage
    real(DP), intent(in), optional :: h_temp !< head in gw cell
    ! -- local
    integer(I4B) :: node
    real(DP) :: wp
    real(DP) :: vscratio
    !
    ! -- initialize conductance
    cond = DZERO
    !
    ! -- initial viscosity ratio to 1
    vscratio = DONE
    !
    ! -- calculate conductance if GWF cell is active
    !    rch-gwf flow will not occur if reach connected to an constant head cell
    node = this%igwfnode(n)
    if (node > 0) then
      if (this%ibound(node) > 0) then
        !
        ! -- direction of gradient across streambed determines which vsc ratio
        if (this%ivsc == 1) then
          if (hsfr > h_temp) then
            ! strm stg > gw head
            vscratio = this%viscratios(1, n)
          else
            vscratio = this%viscratios(2, n)
          end if
        end if
        wp = this%calc_perimeter_wet(n, depth)
        cond = this%hk(n) * vscratio * this%length(n) * wp / this%bthick(n)
      end if
    end if
  end subroutine sfr_calc_cond

  !> @brief Calculate diversion flow
  !!
  !! Method to calculate the diversion flow for a diversion connected
  !! to a SFR package reach. The downstream flow for a reach is passed
  !! in and adjusted by the diversion flow amount calculated in this
  !! method.
  !<
  subroutine sfr_calc_div(this, n, i, qd, qdiv)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    integer(I4B), intent(in) :: i !< diversion number in reach
    real(DP), intent(inout) :: qd !< remaining downstream flow for reach
    real(DP), intent(inout) :: qdiv !< diversion flow for diversion i
    ! -- local
    character(len=10) :: cp
    integer(I4B) :: jpos
    integer(I4B) :: n2
    real(DP) :: v
    !
    ! -- set local variables
    jpos = this%iadiv(n) + i - 1
    n2 = this%divreach(jpos)
    cp = this%divcprior(jpos)
    v = this%divflow(jpos)
    !
    ! -- calculate diversion
    select case (cp)
      ! -- flood diversion
    case ('EXCESS')
      if (qd < v) then
        v = DZERO
      else
        v = qd - v
      end if
      ! -- diversion percentage
    case ('FRACTION')
      v = qd * v
      ! -- STR priority algorithm
    case ('THRESHOLD')
      if (qd < v) then
        v = DZERO
      end if
      ! -- specified diversion
    case ('UPTO')
      if (v > qd) then
        v = qd
      end if
    case default
      v = DZERO
    end select
    !
    ! -- update upstream from for downstream reaches
    qd = qd - v
    qdiv = v
  end subroutine sfr_calc_div

  !> @brief Calculate the depth at the midpoint
  !!
  !! Method to calculate the depth at the midpoint of a reach.
  !<
  subroutine sfr_calc_reach_depth(this, n, q1, d1)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: q1 !< streamflow
    real(DP), intent(inout) :: d1 !< stream depth at midpoint of reach
    ! -- local
    real(DP) :: w
    real(DP) :: s
    real(DP) :: r
    real(DP) :: qconst
    !
    ! -- initialize slope and roughness
    s = this%slope(n)
    r = this%rough(n)
    !
    ! -- calculate stream depth at the midpoint
    if (q1 > DZERO) then
      if (this%ncrosspts(n) > 1) then
        call this%sfr_calc_xs_depth(n, q1, d1)
      else
        w = this%station(this%iacross(n))
        qconst = this%unitconv * w * sqrt(s) / r
        d1 = (q1 / qconst)**DP6
      end if
    else
      d1 = DZERO
    end if
  end subroutine sfr_calc_reach_depth

  !> @brief Calculate the depth at the midpoint of a irregular cross-section
  !!
  !! Method to calculate the depth at the midpoint of a reach with a
  !! irregular cross-section using Newton-Raphson.
  !<
  subroutine sfr_calc_xs_depth(this, n, qrch, d)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: qrch !< streamflow
    real(DP), intent(inout) :: d !< stream depth at midpoint of reach
    ! -- local
    integer(I4B) :: iter
    real(DP) :: perturbation
    real(DP) :: q0
    real(DP) :: q1
    real(DP) :: dq
    real(DP) :: derv
    real(DP) :: dd
    real(DP) :: residual
    !
    ! -- initialize variables
    perturbation = this%deps * DTWO
    d = DZERO
    q0 = DZERO
    residual = q0 - qrch
    !
    ! -- Newton-Raphson iteration
    nriter: do iter = 1, this%maxsfrit
      call this%sfr_calc_qman(n, d + perturbation, q1)
      dq = (q1 - q0)
      if (dq /= DZERO) then
        derv = perturbation / (q1 - q0)
      else
        derv = DZERO
      end if
      dd = derv * residual
      d = d - dd
      call this%sfr_calc_qman(n, d, q0)
      residual = q0 - qrch
      !
      ! -- check for convergence
      if (abs(dd) < this%dmaxchg) then
        exit nriter
      end if
    end do nriter
  end subroutine sfr_calc_xs_depth

  !> @brief Check unit conversion data
  !!
  !! Method to check unit conversion data for a SFR package. This method
  !! also calculates unitconv that is used in the Manning's equation.
  !<
  subroutine sfr_check_conversion(this)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    ! -- local
    ! -- formats
    character(len=*), parameter :: fmtunitconv_error = &
      &"('SFR (',a,') UNIT_CONVERSION SPECIFIED VALUE (',g0,') AND', &
      &1x,'LENGTH_CONVERSION OR TIME_CONVERSION SPECIFIED.')"
    character(len=*), parameter :: fmtunitconv = &
      &"(1x,'SFR PACKAGE (',a,') CONVERSION DATA',&
      &/4x,'UNIT CONVERSION VALUE (',g0,').',/)"
    !
    ! -- check the reach data for simple errors
    if (this%lengthconv /= DNODATA .or. this%timeconv /= DNODATA) then
      if (this%unitconv /= DONE) then
        write (errmsg, fmtunitconv_error) &
          trim(adjustl(this%packName)), this%unitconv
        call store_error(errmsg)
      else
        if (this%lengthconv /= DNODATA) then
          this%unitconv = this%unitconv * this%lengthconv**DONETHIRD
        end if
        if (this%timeconv /= DNODATA) then
          this%unitconv = this%unitconv * this%timeconv
        end if
        write (this%iout, fmtunitconv) &
          trim(adjustl(this%packName)), this%unitconv
      end if
    end if
  end subroutine sfr_check_conversion

  !> @brief Check storage weight
  !!
  !! Method to check the kinematic storage weight for a SFR package.
  !! If the kinematic storage weight has not been set it is set to
  !! the default value.
  !<
  subroutine sfr_check_storage_weight(this)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    ! -- formats
    character(len=*), parameter :: fmtweight = &
      &"(1x,'SFR PACKAGE (',a,') SETTING DEFAULT',&
      &/4x,'STORAGE_WEIGHT VALUE (',g0,').',/)"
    !
    ! -- set storage weight if it has not been defined yet
    if (this%istorage == 1) then
      if (this%storage_weight == DNODATA) then
        this%storage_weight = DONE
        write (this%iout, fmtweight) &
          trim(adjustl(this%packName)), this%storage_weight
      end if
    end if
  end subroutine sfr_check_storage_weight

  !> @brief Check reach data
  !!
  !! Method to check specified data for a SFR package. This method
  !! also creates the tables used to print input data, if this
  !! option in enabled in the SFR package.
  !<
  subroutine sfr_check_reaches(this)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    ! -- local
    character(len=5) :: crch
    character(len=10) :: cval
    character(len=30) :: nodestr
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    integer(I4B) :: n
    integer(I4B) :: nn
    real(DP) :: btgwf
    real(DP) :: bt
    !
    ! -- setup inputtab tableobj
    if (this%iprpak /= 0) then
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(adjustl(this%packName))//') STATIC REACH DATA'
      call table_cr(this%inputtab, this%packName, title)
      call this%inputtab%table_df(this%maxbound, 10, this%iout)
      text = 'NUMBER'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'CELLID'
      call this%inputtab%initialize_column(text, 20, alignment=TABLEFT)
      text = 'LENGTH'
      call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
      text = 'WIDTH'
      call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
      text = 'SLOPE'
      call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
      text = 'TOP'
      call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
      text = 'THICKNESS'
      call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
      text = 'HK'
      call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
      text = 'ROUGHNESS'
      call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
      text = 'UPSTREAM FRACTION'
      call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
    end if
    !
    ! -- check the reach data for simple errors
    do n = 1, this%maxbound
      write (crch, '(i5)') n
      nn = this%igwfnode(n)
      if (nn > 0) then
        btgwf = this%dis%bot(nn)
        call this%dis%noder_to_string(nn, nodestr)
      else
        nodestr = 'none'
      end if
      ! -- check reach length
      if (this%length(n) <= DZERO) then
        errmsg = 'Reach '//crch//' length must be greater than 0.0.'
        call store_error(errmsg)
      end if
      ! -- check reach width
      if (this%width(n) <= DZERO) then
        errmsg = 'Reach '//crch//' width must be greater than 0.0.'
        call store_error(errmsg)
      end if
      ! -- check reach slope
      if (this%slope(n) <= DZERO) then
        errmsg = 'Reach '//crch//' slope must be greater than 0.0.'
        call store_error(errmsg)
      end if
      ! -- check bed thickness and bed hk for reaches connected to GWF
      if (nn > 0) then
        bt = this%strtop(n) - this%bthick(n)
        if (bt <= btgwf .and. this%icheck /= 0) then
          write (cval, '(f10.4)') bt
          errmsg = 'Reach '//crch//' bed bottom (rtp-rbth ='// &
                   cval//') must be greater than the bottom of cell ('// &
                   nodestr
          write (cval, '(f10.4)') btgwf
          errmsg = trim(adjustl(errmsg))//'='//cval//').'
          call store_error(errmsg)
        end if
        if (this%hk(n) < DZERO) then
          errmsg = 'Reach '//crch//' hk must be greater than or equal to 0.0.'
          call store_error(errmsg)
        end if
      end if
      ! -- check reach roughness
      if (this%rough(n) <= DZERO) then
        errmsg = 'Reach '//crch//" Manning's roughness "// &
                 'coefficient must be greater than 0.0.'
        call store_error(errmsg)
      end if
      ! -- check reach upstream fraction
      if (this%ustrf(n) < DZERO) then
        errmsg = 'Reach '//crch//' upstream fraction must be greater '// &
                 'than or equal to 0.0.'
        call store_error(errmsg)
      end if
      ! -- write summary of reach information
      if (this%iprpak /= 0) then
        call this%inputtab%add_term(n)
        call this%inputtab%add_term(nodestr)
        call this%inputtab%add_term(this%length(n))
        call this%inputtab%add_term(this%width(n))
        call this%inputtab%add_term(this%slope(n))
        call this%inputtab%add_term(this%strtop(n))
        call this%inputtab%add_term(this%bthick(n))
        call this%inputtab%add_term(this%hk(n))
        call this%inputtab%add_term(this%rough(n))
        call this%inputtab%add_term(this%ustrf(n))
      end if
    end do
  end subroutine sfr_check_reaches

  !> @brief Check connection data
  !!
  !! Method to check connection data for a SFR package. This method
  !! also creates the tables used to print input data, if this
  !! option in enabled in the SFR package.
  !<
  subroutine sfr_check_connections(this)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    ! -- local
    logical(LGP) :: lreorder
    character(len=5) :: crch
    character(len=5) :: crch2
    character(len=LINELENGTH) :: text
    character(len=LINELENGTH) :: title
    integer(I4B) :: n
    integer(I4B) :: nn
    integer(I4B) :: nc
    integer(I4B) :: i
    integer(I4B) :: ii
    integer(I4B) :: j
    integer(I4B) :: ifound
    integer(I4B) :: ierr
    integer(I4B) :: maxconn
    integer(I4B) :: ntabcol
    !
    ! -- determine if the reaches have been reordered
    lreorder = .FALSE.
    do j = 1, this%MAXBOUND
      n = this%isfrorder(j)
      if (n /= j) then
        lreorder = .TRUE.
        exit
      end if
    end do
    !
    ! -- write message that the solution order h
    if (lreorder) then
      write (this%iout, '(/,1x,a)') &
        trim(adjustl(this%text))//' PACKAGE ('// &
        trim(adjustl(this%packName))//') REACH SOLUTION HAS BEEN '// &
        'REORDERED USING A DAG'
      !
      ! -- print table
      if (this%iprpak /= 0) then
        !
        ! -- reset the input table object
        ntabcol = 2
        title = trim(adjustl(this%text))//' PACKAGE ('// &
                trim(adjustl(this%packName))//') REACH SOLUTION ORDER'
        call table_cr(this%inputtab, this%packName, title)
        call this%inputtab%table_df(this%maxbound, ntabcol, this%iout)
        text = 'ORDER'
        call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
        text = 'REACH'
        call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
        !
        ! -- upstream connection data
        do j = 1, this%maxbound
          n = this%isfrorder(j)
          call this%inputtab%add_term(j)
          call this%inputtab%add_term(n)
        end do
      end if
    end if
    !
    ! -- create input table for reach connections data
    if (this%iprpak /= 0) then
      !
      ! -- calculate the maximum number of connections
      maxconn = 0
      do n = 1, this%maxbound
        maxconn = max(maxconn, this%nconnreach(n))
      end do
      ntabcol = 1 + maxconn
      !
      ! -- reset the input table object
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(adjustl(this%packName))//') STATIC REACH CONNECTION DATA'
      call table_cr(this%inputtab, this%packName, title)
      call this%inputtab%table_df(this%maxbound, ntabcol, this%iout)
      text = 'REACH'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      do n = 1, maxconn
        write (text, '(a,1x,i6)') 'CONN', n
        call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      end do
    end if
    !
    ! -- check the reach connections for simple errors
    ! -- connection check
    do n = 1, this%MAXBOUND
      write (crch, '(i5)') n
      eachconn: do i = this%ia(n) + 1, this%ia(n + 1) - 1
        nn = this%ja(i)
        write (crch2, '(i5)') nn
        ifound = 0
        connreach: do ii = this%ia(nn) + 1, this%ia(nn + 1) - 1
          nc = this%ja(ii)
          if (nc == n) then
            ifound = 1
            exit connreach
          end if
        end do connreach
        if (ifound /= 1) then
          errmsg = 'Reach '//crch//' is connected to '// &
                   'reach '//crch2//' but reach '//crch2// &
                   ' is not connected to reach '//crch//'.'
          call store_error(errmsg)
        end if
      end do eachconn
      !
      ! -- write connection data to the table
      if (this%iprpak /= 0) then
        call this%inputtab%add_term(n)
        do i = this%ia(n) + 1, this%ia(n + 1) - 1
          call this%inputtab%add_term(this%ja(i))
        end do
        nn = maxconn - this%nconnreach(n)
        do i = 1, nn
          call this%inputtab%add_term(' ')
        end do
      end if
    end do
    !
    ! -- check for incorrect connections between upstream connections
    !
    ! -- check upstream connections for each reach
    ierr = 0
    do n = 1, this%maxbound
      write (crch, '(i5)') n
      eachconnv: do i = this%ia(n) + 1, this%ia(n + 1) - 1
        !
        ! -- skip downstream connections
        if (this%idir(i) < 0) cycle eachconnv
        nn = this%ja(i)
        write (crch2, '(i5)') nn
        connreachv: do ii = this%ia(nn) + 1, this%ia(nn + 1) - 1
          ! -- skip downstream connections
          if (this%idir(ii) < 0) cycle connreachv
          nc = this%ja(ii)
          !
          ! -- if nc == n then that means reach n is an upstream connection for
          !    reach nn and reach nn is an upstream connection for reach n
          if (nc == n) then
            ierr = ierr + 1
            errmsg = 'Reach '//crch//' is connected to '// &
                     'reach '//crch2//' but streamflow from reach '// &
                     crch//' to reach '//crch2//' is not permitted.'
            call store_error(errmsg)
            exit connreachv
          end if
        end do connreachv
      end do eachconnv
    end do
    !
    ! -- terminate if connectivity errors
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- check that downstream reaches for a reach are
    !    the upstream reaches for the reach
    do n = 1, this%maxbound
      write (crch, '(i5)') n
      eachconnds: do i = this%ia(n) + 1, this%ia(n + 1) - 1
        nn = this%ja(i)
        if (this%idir(i) > 0) cycle eachconnds
        write (crch2, '(i5)') nn
        ifound = 0
        connreachds: do ii = this%ia(nn) + 1, this%ia(nn + 1) - 1
          nc = this%ja(ii)
          if (nc == n) then
            if (this%idir(i) /= this%idir(ii)) then
              ifound = 1
            end if
            exit connreachds
          end if
        end do connreachds
        if (ifound /= 1) then
          errmsg = 'Reach '//crch//' downstream connected reach '// &
                   'is reach '//crch2//' but reach '//crch//' is not'// &
                   ' the upstream connected reach for reach '//crch2//'.'
          call store_error(errmsg)
        end if
      end do eachconnds
    end do
    !
    ! -- create input table for upstream and downstream connections
    if (this%iprpak /= 0) then
      !
      ! -- calculate the maximum number of upstream connections
      maxconn = 0
      do n = 1, this%maxbound
        ii = 0
        do i = this%ia(n) + 1, this%ia(n + 1) - 1
          if (this%idir(i) > 0) then
            ii = ii + 1
          end if
        end do
        maxconn = max(maxconn, ii)
      end do
      ntabcol = 1 + maxconn
      !
      ! -- reset the input table object
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(adjustl(this%packName))//') STATIC UPSTREAM REACH '// &
              'CONNECTION DATA'
      call table_cr(this%inputtab, this%packName, title)
      call this%inputtab%table_df(this%maxbound, ntabcol, this%iout)
      text = 'REACH'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      do n = 1, maxconn
        write (text, '(a,1x,i6)') 'UPSTREAM CONN', n
        call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      end do
      !
      ! -- upstream connection data
      do n = 1, this%maxbound
        call this%inputtab%add_term(n)
        ii = 0
        do i = this%ia(n) + 1, this%ia(n + 1) - 1
          if (this%idir(i) > 0) then
            call this%inputtab%add_term(this%ja(i))
            ii = ii + 1
          end if
        end do
        nn = maxconn - ii
        do i = 1, nn
          call this%inputtab%add_term(' ')
        end do
      end do
      !
      ! -- calculate the maximum number of downstream connections
      maxconn = 0
      do n = 1, this%maxbound
        ii = 0
        do i = this%ia(n) + 1, this%ia(n + 1) - 1
          if (this%idir(i) < 0) then
            ii = ii + 1
          end if
        end do
        maxconn = max(maxconn, ii)
      end do
      ntabcol = 1 + maxconn
      !
      ! -- reset the input table object
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(adjustl(this%packName))//') STATIC DOWNSTREAM '// &
              'REACH CONNECTION DATA'
      call table_cr(this%inputtab, this%packName, title)
      call this%inputtab%table_df(this%maxbound, ntabcol, this%iout)
      text = 'REACH'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      do n = 1, maxconn
        write (text, '(a,1x,i6)') 'DOWNSTREAM CONN', n
        call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      end do
      !
      ! -- downstream connection data
      do n = 1, this%maxbound
        call this%inputtab%add_term(n)
        ii = 0
        do i = this%ia(n) + 1, this%ia(n + 1) - 1
          if (this%idir(i) < 0) then
            call this%inputtab%add_term(this%ja(i))
            ii = ii + 1
          end if
        end do
        nn = maxconn - ii
        do i = 1, nn
          call this%inputtab%add_term(' ')
        end do
      end do
    end if
  end subroutine sfr_check_connections

  !> @brief Check diversions data
  !!
  !! Method to check diversion data for a SFR package. This method
  !! also creates the tables used to print input data, if this
  !! option in enabled in the SFR package.
  !<
  subroutine sfr_check_diversions(this)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    character(len=5) :: crch
    character(len=5) :: cdiv
    character(len=5) :: crch2
    character(len=10) :: cprior
    integer(I4B) :: maxdiv
    integer(I4B) :: n
    integer(I4B) :: nn
    integer(I4B) :: nc
    integer(I4B) :: ii
    integer(I4B) :: idiv
    integer(I4B) :: ifound
    integer(I4B) :: jpos
    ! -- format
10  format('Diversion ', i0, ' of reach ', i0, &
           ' is invalid or has not been defined.')
    !
    ! -- write header
    if (this%iprpak /= 0) then
      !
      ! -- determine the maximum number of diversions
      maxdiv = 0
      do n = 1, this%maxbound
        maxdiv = maxdiv + this%ndiv(n)
      end do
      !
      ! -- reset the input table object
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(adjustl(this%packName))//') REACH DIVERSION DATA'
      call table_cr(this%inputtab, this%packName, title)
      call this%inputtab%table_df(maxdiv, 4, this%iout)
      text = 'REACH'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'DIVERSION'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'REACH 2'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'CPRIOR'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
    end if
    !
    ! -- check that diversion data are correct
    do n = 1, this%maxbound
      if (this%ndiv(n) < 1) cycle
      write (crch, '(i5)') n

      do idiv = 1, this%ndiv(n)
        !
        ! -- determine diversion index
        jpos = this%iadiv(n) + idiv - 1
        !
        ! -- write idiv to cdiv
        write (cdiv, '(i5)') idiv
        !
        !
        nn = this%divreach(jpos)
        write (crch2, '(i5)') nn
        !
        ! -- make sure diversion reach is connected to current reach
        ifound = 0
        if (nn < 1 .or. nn > this%maxbound) then
          write (errmsg, 10) idiv, n
          call store_error(errmsg)
          cycle
        end if
        connreach: do ii = this%ia(nn) + 1, this%ia(nn + 1) - 1
          nc = this%ja(ii)
          if (nc == n) then
            if (this%idir(ii) > 0) then
              ifound = 1
            end if
            exit connreach
          end if
        end do connreach
        if (ifound /= 1) then
          errmsg = 'Reach '//crch//' is not a upstream reach for '// &
                   'reach '//crch2//' as a result diversion '//cdiv// &
                   ' from reach '//crch//' to reach '//crch2// &
                   ' is not possible. Check reach connectivity.'
          call store_error(errmsg)
        end if
        ! -- iprior
        cprior = this%divcprior(jpos)
        !
        ! -- add terms to the table
        if (this%iprpak /= 0) then
          call this%inputtab%add_term(n)
          call this%inputtab%add_term(idiv)
          call this%inputtab%add_term(nn)
          call this%inputtab%add_term(cprior)
        end if
      end do
    end do
  end subroutine sfr_check_diversions

  !> @brief Check initial stage data
  !!
  !! Method to check initial data for a SFR package and calculates
  !! the initial upstream and downstream flows for the reach based
  !! on the initial staalso creates the tables used to print input
  !! data, if this option in enabled in the SFR package.
  !<
  subroutine sfr_check_initialstages(this)
    class(SfrType) :: this !< SfrType object

    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    character(len=5) :: crch
    integer(I4B) :: n
    real(DP) :: qman

    ! skip check if storage is not activated
    if (this%istorage == 0) return

    ! write header
    if (this%iprpak /= 0) then
      !
      ! -- reset the input table object
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(adjustl(this%packName))//') REACH INITIAL STAGE DATA'
      call table_cr(this%inputtab, this%packName, title)
      call this%inputtab%table_df(this%maxbound, 4, this%iout)
      text = 'REACH'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'INITIAL STAGE'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'INITIAL DEPTH'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'INITIAL FLOW'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
    end if
    !
    ! -- check that data are correct
    do n = 1, this%maxbound
      write (crch, '(i5)') n

      ! calculate the initial flows
      call this%sfr_calc_qman(n, this%depth(n), qman)
      this%usinflow(n) = qman
      this%dsflow(n) = qman

      ! add terms to the table
      if (this%iprpak /= 0) then
        call this%inputtab%add_term(n)
        call this%inputtab%add_term(this%stage(n))
        call this%inputtab%add_term(this%depth(n))
        call this%inputtab%add_term(qman)
      end if
    end do
  end subroutine sfr_check_initialstages

  !> @brief Check upstream fraction data
  !!
  !! Method to check upstream fraction data for a SFR package.
  !! This method  also creates the tables used to print input data,
  !! if this option in enabled in the SFR package.
  !<
  subroutine sfr_check_ustrf(this)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    logical(LGP) :: lcycle
    logical(LGP) :: ladd
    character(len=5) :: crch
    character(len=5) :: crch2
    character(len=10) :: cval
    integer(I4B) :: maxcols
    integer(I4B) :: npairs
    integer(I4B) :: ipair
    integer(I4B) :: i
    integer(I4B) :: n
    integer(I4B) :: n2
    integer(I4B) :: idiv
    integer(I4B) :: i0
    integer(I4B) :: i1
    integer(I4B) :: jpos
    integer(I4B) :: ids
    real(DP) :: f
    real(DP) :: rval
    !
    ! -- write table header
    if (this%iprpak /= 0) then
      !
      ! -- determine the maximum number of columns
      npairs = 0
      do n = 1, this%maxbound
        ipair = 0
        ec: do i = this%ia(n) + 1, this%ia(n + 1) - 1
          !
          ! -- skip upstream connections
          if (this%idir(i) > 0) cycle ec
          n2 = this%ja(i)
          !
          ! -- skip inactive downstream reaches
          if (this%iboundpak(n2) == 0) cycle ec
          !
          ! -- increment ipair and see if it exceeds npairs
          ipair = ipair + 1
          npairs = max(npairs, ipair)
        end do ec
      end do
      maxcols = 1 + npairs * 2
      !
      ! -- reset the input table object
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(adjustl(this%packName))//') CONNECTED REACH UPSTREAM '// &
              'FRACTION DATA'
      call table_cr(this%inputtab, this%packName, title)
      call this%inputtab%table_df(this%maxbound, maxcols, this%iout)
      text = 'REACH'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      do i = 1, npairs
        write (cval, '(i10)') i
        text = 'DOWNSTREAM REACH '//trim(adjustl(cval))
        call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
        text = 'FRACTION '//trim(adjustl(cval))
        call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
      end do
    end if
    !
    ! -- fill diversion number for each connection
    do n = 1, this%maxbound
      do idiv = 1, this%ndiv(n)
        i0 = this%iadiv(n)
        i1 = this%iadiv(n + 1) - 1
        do jpos = i0, i1
          do i = this%ia(n) + 1, this%ia(n + 1) - 1
            n2 = this%ja(i)
            if (this%divreach(jpos) == n2) then
              this%idiv(i) = jpos - i0 + 1
              exit
            end if
          end do
        end do
      end do
    end do
    !
    ! -- check that the upstream fraction for reach connected by
    !    a diversion is zero
    do n = 1, this%maxbound
      !
      ! -- determine the number of downstream reaches
      ids = 0
      do i = this%ia(n) + 1, this%ia(n + 1) - 1
        if (this%idir(i) < 0) then
          ids = ids + 1
        end if
      end do
      !
      ! -- evaluate the diversions
      do idiv = 1, this%ndiv(n)
        jpos = this%iadiv(n) + idiv - 1
        n2 = this%divreach(jpos)
        f = this%ustrf(n2)
        if (f /= DZERO) then
          write (errmsg, '(a,2(1x,i0,1x,a),1x,a,g0,a,2(1x,a))') &
            'Reach', n, 'is connected to reach', n2, 'by a diversion', &
            'but the upstream fraction is not equal to zero (', f, '). Check', &
            trim(this%packName), 'package diversion and package data.'
          if (ids > 1) then
            call store_error(errmsg)
          else
            write (warnmsg, '(a,3(1x,a))') &
              trim(warnmsg), &
              'A warning instead of an error is issued because', &
              'the reach is only connected to the diversion reach in the ', &
              'downstream direction.'
            call store_warning(warnmsg)
          end if
        end if
      end do
    end do
    !
    ! -- calculate the total fraction of connected reaches that are
    !    not diversions and check that the sum of upstream fractions
    !    is equal to 1 for each reach
    do n = 1, this%maxbound
      ids = 0
      rval = DZERO
      f = DZERO
      write (crch, '(i5)') n
      if (this%iprpak /= 0) then
        call this%inputtab%add_term(n)
      end if
      ipair = 0
      eachconn: do i = this%ia(n) + 1, this%ia(n + 1) - 1
        lcycle = .FALSE.
        !
        ! -- initialize downstream connection q
        this%qconn(i) = DZERO
        !
        ! -- skip upstream connections
        if (this%idir(i) > 0) then
          lcycle = .TRUE.
        end if
        n2 = this%ja(i)
        !
        ! -- skip inactive downstream reaches
        if (this%iboundpak(n2) == 0) then
          lcycle = .TRUE.
        end if
        if (lcycle) then
          cycle eachconn
        end if
        ipair = ipair + 1
        write (crch2, '(i5)') n2
        ids = ids + 1
        ladd = .true.
        f = f + this%ustrf(n2)
        write (cval, '(f10.4)') this%ustrf(n2)
        !
        ! -- write upstream fractions
        if (this%iprpak /= 0) then
          call this%inputtab%add_term(n2)
          call this%inputtab%add_term(this%ustrf(n2))
        end if
        eachdiv: do idiv = 1, this%ndiv(n)
          jpos = this%iadiv(n) + idiv - 1
          if (this%divreach(jpos) == n2) then
            ladd = .false.
            exit eachdiv
          end if
        end do eachdiv
        if (ladd) then
          rval = rval + this%ustrf(n2)
        end if
      end do eachconn
      this%ftotnd(n) = rval
      !
      ! -- write remaining table columns
      if (this%iprpak /= 0) then
        ipair = ipair + 1
        do i = ipair, npairs
          call this%inputtab%add_term('  ')
          call this%inputtab%add_term('  ')
        end do
      end if
      !
      ! -- evaluate if an error condition has occurred
      !    the sum of fractions is not equal to 1
      if (ids /= 0) then
        if (abs(f - DONE) > DEM6) then
          write (errmsg, '(a,1x,i0,1x,a,g0,a,3(1x,a))') &
            'Upstream fractions for reach ', n, 'is not equal to one (', f, &
            '). Check', trim(this%packName), 'package reach connectivity and', &
            'package data.'
          call store_error(errmsg)
        end if
      end if
    end do
  end subroutine sfr_check_ustrf

  !> @brief Setup budget object for package
  !!
  !! Method to set up the budget object that stores all the sfr flows
  !! The terms listed here must correspond in number and order to the ones
  !! listed in the sfr_fill_budobj method.
  !<
  subroutine sfr_setup_budobj(this)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    ! -- local
    integer(I4B) :: nbudterm
    integer(I4B) :: i
    integer(I4B) :: n
    integer(I4B) :: n1
    integer(I4B) :: n2
    integer(I4B) :: maxlist
    integer(I4B) :: naux
    integer(I4B) :: idx
    real(DP) :: q
    character(len=LENBUDTXT) :: text
    character(len=LENBUDTXT), dimension(1) :: auxtxt
    !
    ! -- Determine the number of sfr budget terms. These are fixed for
    !    the simulation and cannot change.  This includes FLOW-JA-FACE
    !    so they can be written to the binary budget files, but these internal
    !    flows are not included as part of the budget table.
    nbudterm = 8
    if (this%imover == 1) nbudterm = nbudterm + 2
    if (this%naux > 0) nbudterm = nbudterm + 1
    !
    ! -- set up budobj
    call budgetobject_cr(this%budobj, this%packName)
    call this%budobj%budgetobject_df(this%maxbound, nbudterm, 0, 0, &
                                     ibudcsv=this%ibudcsv)
    idx = 0
    !
    ! -- Go through and set up each budget term
    text = '    FLOW-JA-FACE'
    idx = idx + 1
    maxlist = this%nconn
    naux = 1
    auxtxt(1) = '       FLOW-AREA'
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux, auxtxt)
    !
    ! -- store connectivity
    call this%budobj%budterm(idx)%reset(this%nconn)
    q = DZERO
    do n = 1, this%maxbound
      n1 = n
      do i = this%ia(n) + 1, this%ia(n + 1) - 1
        n2 = this%ja(i)
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
      end do
    end do
    !
    ! --
    text = '             GWF'
    idx = idx + 1
    maxlist = this%maxbound - this%ianynone
    naux = 1
    auxtxt(1) = '       FLOW-AREA'
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%name_model, &
                                             maxlist, .false., .true., &
                                             naux, auxtxt)
    call this%budobj%budterm(idx)%reset(maxlist)
    q = DZERO
    do n = 1, this%maxbound
      n2 = this%igwfnode(n)
      if (n2 > 0) then
        call this%budobj%budterm(idx)%update_term(n, n2, q)
      end if
    end do
    !
    ! --
    text = '        RAINFALL'
    idx = idx + 1
    maxlist = this%maxbound
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
    maxlist = this%maxbound
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
    maxlist = this%maxbound
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
    maxlist = this%maxbound
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
    maxlist = this%maxbound
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
    maxlist = this%maxbound
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
    if (this%imover == 1) then
      !
      ! --
      text = '        FROM-MVR'
      idx = idx + 1
      maxlist = this%maxbound
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
      maxlist = this%maxbound
      naux = 0
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux)
    end if
    !
    ! --
    naux = this%naux
    if (naux > 0) then
      !
      ! --
      text = '       AUXILIARY'
      idx = idx + 1
      maxlist = this%maxbound
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux, this%auxname)
    end if
    !
    ! -- if sfr flow for each reach are written to the listing file
    if (this%iprflow /= 0) then
      call this%budobj%flowtable_df(this%iout, cellids='GWF')
    end if
  end subroutine sfr_setup_budobj

  !> @brief Copy flow terms into budget object for package
  !!
  !! Method to copy flows into the budget object that stores all the sfr flows
  !! The terms listed here must correspond in number and order to the ones
  !! added in the sfr_setup_budobj method.
  !<
  subroutine sfr_fill_budobj(this)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    ! -- local
    integer(I4B) :: naux
    integer(I4B) :: i
    integer(I4B) :: n
    integer(I4B) :: n1
    integer(I4B) :: n2
    integer(I4B) :: ii
    integer(I4B) :: idx
    integer(I4B) :: idiv
    integer(I4B) :: jpos
    real(DP) :: q
    real(DP) :: qt
    real(DP) :: d
    real(DP) :: ca
    real(DP) :: a
    real(DP) :: wp
    real(DP) :: l
    !
    ! -- initialize counter
    idx = 0
    !
    ! -- FLOW JA FACE
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nconn)
    do n = 1, this%maxbound
      n1 = n
      q = DZERO
      ca = DZERO
      do i = this%ia(n) + 1, this%ia(n + 1) - 1
        n2 = this%ja(i)
        if (this%iboundpak(n) /= 0) then
          ! flow to downstream reaches
          if (this%idir(i) < 0) then
            qt = this%dsflow(n)
            q = -this%qconn(i)
            ! flow from upstream reaches
          else
            qt = this%usflow(n)
            do ii = this%ia(n2) + 1, this%ia(n2 + 1) - 1
              if (this%idir(ii) > 0) cycle
              if (this%ja(ii) /= n) cycle
              q = this%qconn(ii)
              exit
            end do
          end if
          ! calculate flow area
          call this%sfr_calc_reach_depth(n, qt, d)
          ca = this%calc_area_wet(n, d)
        else
          q = DZERO
          ca = DZERO
        end if
        this%qauxcbc(1) = ca
        call this%budobj%budterm(idx)%update_term(n1, n2, q, this%qauxcbc)
      end do
    end do
    !
    ! -- GWF (LEAKAGE)
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%maxbound - this%ianynone)
    do n = 1, this%maxbound
      n2 = this%igwfnode(n)
      if (n2 > 0) then
        if (this%iboundpak(n) /= 0) then
          ! -- calc_perimeter_wet() does not enforce depth dependence
          if (this%depth(n) > DZERO) then
            wp = this%calc_perimeter_wet(n, this%depth(n))
          else
            wp = DZERO
          end if
          l = this%length(n)
          a = wp * l
          this%qauxcbc(1) = a
          q = -this%gwflow(n)
        else
          this%qauxcbc(1) = DZERO
          q = DZERO
        end if
        call this%budobj%budterm(idx)%update_term(n, n2, q, this%qauxcbc)
      end if
    end do
    !
    ! -- RAIN
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%maxbound)
    do n = 1, this%maxbound
      if (this%iboundpak(n) /= 0) then
        a = this%calc_surface_area(n)
        q = this%rain(n) * a
      else
        q = DZERO
      end if
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- EVAPORATION
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%maxbound)
    do n = 1, this%maxbound
      if (this%iboundpak(n) /= 0) then
        q = -this%simevap(n)
      else
        q = DZERO
      end if
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- RUNOFF
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%maxbound)
    do n = 1, this%maxbound
      if (this%iboundpak(n) /= 0) then
        q = this%simrunoff(n)
      else
        q = DZERO
      end if
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- INFLOW
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%maxbound)
    do n = 1, this%maxbound
      if (this%iboundpak(n) /= 0) then
        q = this%inflow(n)
      else
        q = DZERO
      end if
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- EXTERNAL OUTFLOW
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%maxbound)
    do n = 1, this%maxbound
      q = DZERO
      if (this%iboundpak(n) /= 0) then
        do i = this%ia(n) + 1, this%ia(n + 1) - 1
          if (this%idir(i) > 0) cycle
          idiv = this%idiv(i)
          if (idiv > 0) then
            jpos = this%iadiv(n) + idiv - 1
            q = q + this%divq(jpos)
          else
            q = q + this%qconn(i)
          end if
        end do
        q = q - this%dsflow(n)
        if (this%imover == 1) then
          q = q + this%pakmvrobj%get_qtomvr(n)
        end if
      else
        if (this%imover == 1) then
          q = this%pakmvrobj%get_qfrommvr(n)
        end if
      end if
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- STORAGE
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%maxbound)
    do n = 1, this%maxbound
      q = DZERO
      if (this%iboundpak(n) /= 0) then
        d = this%depth(n)
        a = this%calc_surface_area_wet(n, d)
        this%qauxcbc(1) = a * d
        if (this%gwfiss == 0 .and. this%istorage == 1) then
          q = this%storage(n)
        end if
      else
        q = DZERO
        this%qauxcbc(1) = DZERO
      end if
      call this%budobj%budterm(idx)%update_term(n, n, q, this%qauxcbc)
    end do
    !
    ! -- MOVER
    if (this%imover == 1) then
      !
      ! -- FROM MOVER
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%maxbound)
      do n = 1, this%maxbound
        q = DZERO
        if (this%iboundpak(n) /= 0) then
          q = this%pakmvrobj%get_qfrommvr(n)
        end if
        call this%budobj%budterm(idx)%update_term(n, n, q)
      end do
      !
      ! -- TO MOVER
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%maxbound)
      do n = 1, this%maxbound
        if (this%iboundpak(n) /= 0) then
          q = this%pakmvrobj%get_qtomvr(n)
          if (q > DZERO) then
            q = -q
          end if
        else
          q = DZERO
        end if
        call this%budobj%budterm(idx)%update_term(n, n, q)
      end do
    end if
    !
    ! -- AUXILIARY VARIABLES
    naux = this%naux
    if (naux > 0) then
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%maxbound)
      do n = 1, this%maxbound
        q = DZERO
        call this%budobj%budterm(idx)%update_term(n, n, q, this%auxvar(:, n))
      end do
    end if
    !
    ! --Terms are filled, now accumulate them for this time step
    call this%budobj%accumulate_terms()
  end subroutine sfr_fill_budobj

  !> @brief Setup stage table object for package
  !!
  !! Method to set up the table object that is used to write the sfr
  !! stage data. The terms listed here must correspond in number and
  !! order to the ones written to the stage table in the sfr_ot method.
  !<
  subroutine sfr_setup_tableobj(this)
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    ! -- local
    integer(I4B) :: nterms
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    !
    ! -- setup stage table
    if (this%iprhed > 0) then
      !
      ! -- Determine the number of sfr budget terms. These are fixed for
      !    the simulation and cannot change.  This includes FLOW-JA-FACE
      !    so they can be written to the binary budget files, but these internal
      !    flows are not included as part of the budget table.
      nterms = 8
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
      call this%stagetab%table_df(this%maxbound, nterms, this%iout, &
                                  transient=.TRUE.)
      !
      ! -- Go through and set up table budget term
      if (this%inamedbound == 1) then
        text = 'NAME'
        call this%stagetab%initialize_column(text, LENBOUNDNAME, &
                                             alignment=TABLEFT)
      end if
      !
      ! -- reach number
      text = 'NUMBER'
      call this%stagetab%initialize_column(text, 10, alignment=TABCENTER)
      !
      ! -- cellids
      text = 'CELLID'
      call this%stagetab%initialize_column(text, 20, alignment=TABLEFT)
      !
      ! -- reach stage
      text = 'STAGE'
      call this%stagetab%initialize_column(text, 12, alignment=TABCENTER)
      !
      ! -- reach depth
      text = 'DEPTH'
      call this%stagetab%initialize_column(text, 12, alignment=TABCENTER)
      !
      ! -- reach width
      text = 'WIDTH'
      call this%stagetab%initialize_column(text, 12, alignment=TABCENTER)
      !
      ! -- gwf head
      text = 'GWF HEAD'
      call this%stagetab%initialize_column(text, 12, alignment=TABCENTER)
      !
      ! -- streambed conductance
      text = 'STREAMBED CONDUCTANCE'
      call this%stagetab%initialize_column(text, 12, alignment=TABCENTER)
      !
      ! -- streambed gradient
      text = 'STREAMBED GRADIENT'
      call this%stagetab%initialize_column(text, 12, alignment=TABCENTER)
    end if
  end subroutine sfr_setup_tableobj

  ! -- reach geometry functions

  !> @brief Calculate wetted area
  !!
  !! Function to calculate the wetted area for a SFR package reach.
  !<
  function calc_area_wet(this, n, depth)
    ! -- return variable
    real(DP) :: calc_area_wet !< wetted area
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: depth !< reach depth
    ! -- local
    integer(I4B) :: npts
    integer(I4B) :: i0
    integer(I4B) :: i1
    !
    ! -- Calculate wetted area
    npts = this%ncrosspts(n)
    i0 = this%iacross(n)
    i1 = this%iacross(n + 1) - 1
    if (npts > 1) then
      calc_area_wet = get_cross_section_area(npts, this%station(i0:i1), &
                                             this%xsheight(i0:i1), depth)
    else
      calc_area_wet = this%station(i0) * depth
    end if
  end function calc_area_wet

  !> @brief Calculate wetted perimeter
  !!
  !! Function to calculate the wetted perimeter for a SFR package reach.
  !<
  function calc_perimeter_wet(this, n, depth)
    ! -- return variable
    real(DP) :: calc_perimeter_wet !< wetted perimeter
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: depth !< reach depth
    ! -- local
    integer(I4B) :: npts
    integer(I4B) :: i0
    integer(I4B) :: i1
    !
    ! -- Calculate wetted perimeter
    npts = this%ncrosspts(n)
    i0 = this%iacross(n)
    i1 = this%iacross(n + 1) - 1
    if (npts > 1) then
      calc_perimeter_wet = get_wetted_perimeter(npts, this%station(i0:i1), &
                                                this%xsheight(i0:i1), depth)
    else
      calc_perimeter_wet = this%station(i0) ! no depth dependence in original implementation
    end if
  end function calc_perimeter_wet

  !> @brief Calculate maximum surface area
  !!
  !! Function to calculate the maximum surface area for a SFR package reach.
  !<
  function calc_surface_area(this, n)
    ! -- return variable
    real(DP) :: calc_surface_area !< surface area
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    ! -- local
    integer(I4B) :: npts
    integer(I4B) :: i0
    integer(I4B) :: i1
    real(DP) :: top_width
    !
    ! -- Calculate surface area
    npts = this%ncrosspts(n)
    i0 = this%iacross(n)
    i1 = this%iacross(n + 1) - 1
    if (npts > 1) then
      top_width = get_saturated_topwidth(npts, this%station(i0:i1))
    else
      top_width = this%station(i0)
    end if
    calc_surface_area = top_width * this%length(n)
  end function calc_surface_area

  !> @brief Calculate wetted surface area
  !!
  !! Function to calculate the wetted surface area for a SFR package reach.
  !<
  function calc_surface_area_wet(this, n, depth)
    ! -- return variable
    real(DP) :: calc_surface_area_wet !< wetted surface area
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: depth !< reach depth
    ! -- local
    real(DP) :: top_width
    !
    ! -- Calculate wetted surface area
    top_width = this%calc_top_width_wet(n, depth)
    calc_surface_area_wet = top_width * this%length(n)
  end function calc_surface_area_wet

  !> @brief Calculate wetted top width
  !!
  !! Function to calculate the wetted top width for a SFR package reach.
  !<
  function calc_top_width_wet(this, n, depth)
    ! -- return variable
    real(DP) :: calc_top_width_wet !< wetted top width
    ! -- dummy
    class(SfrType) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: depth !< reach depth
    ! -- local
    integer(I4B) :: npts
    integer(I4B) :: i0
    integer(I4B) :: i1
    real(DP) :: sat
    !
    ! -- Calculate wetted top width
    npts = this%ncrosspts(n)
    i0 = this%iacross(n)
    i1 = this%iacross(n + 1) - 1
    sat = sCubicSaturation(DEM5, DZERO, depth, DEM5)
    if (npts > 1) then
      calc_top_width_wet = sat * get_wetted_topwidth(npts, &
                                                     this%station(i0:i1), &
                                                     this%xsheight(i0:i1), &
                                                     depth)
    else
      calc_top_width_wet = sat * this%station(i0)
    end if
  end function calc_top_width_wet

  !> @brief Activate density terms
  !!
  !! Method to activate addition of density terms for a SFR package reach.
  !<
  subroutine sfr_activate_density(this)
    ! -- modules
    use MemoryManagerModule, only: mem_reallocate
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: j
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
    write (this%iout, '(/1x,a)') 'DENSITY TERMS HAVE BEEN ACTIVATED FOR SFR &
      &PACKAGE: '//trim(adjustl(this%packName))
  end subroutine sfr_activate_density

  !> @brief Activate viscosity terms
  !!
  !! Method to activate addition of viscosity terms for exchange
  !! with groundwater along a SFR package reach.
  !<
  subroutine sfr_activate_viscosity(this)
    ! -- modules
    use MemoryManagerModule, only: mem_reallocate
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    ! -- local
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
    write (this%iout, '(/1x,a)') 'VISCOSITY HAS BEEN ACTIVATED FOR SFR &
      &PACKAGE: '//trim(adjustl(this%packName))
  end subroutine sfr_activate_viscosity

  !> @brief Calculate density terms
  !!
  !! Method to galculate groundwater-reach density exchange terms for a
  !! SFR package reach.
  !!
  !! Member variable used here
  !!   denseterms  : shape (3, MAXBOUND), filled by buoyancy package
  !!                   col 1 is relative density of sfr (densesfr / denseref)
  !!                   col 2 is relative density of gwf cell (densegwf / denseref)
  !!                   col 3 is elevation of gwf cell
  !<
  subroutine sfr_calculate_density_exchange(this, n, stage, head, cond, &
                                            bots, flow, gwfhcof, gwfrhs)
    ! -- dummy
    class(SfrType), intent(inout) :: this !< SfrType object
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: stage !< reach stage
    real(DP), intent(in) :: head !< head in connected GWF cell
    real(DP), intent(in) :: cond !< reach conductance
    real(DP), intent(in) :: bots !< bottom elevation of reach
    real(DP), intent(inout) :: flow !< calculated flow, updated here with density terms
    real(DP), intent(inout) :: gwfhcof !< GWF diagonal coefficient, updated here with density terms
    real(DP), intent(inout) :: gwfrhs !< GWF right-hand-side value, updated here with density terms
    ! -- local
    real(DP) :: ss
    real(DP) :: hh
    real(DP) :: havg
    real(DP) :: rdensesfr
    real(DP) :: rdensegwf
    real(DP) :: rdenseavg
    real(DP) :: elevsfr
    real(DP) :: elevgwf
    real(DP) :: elevavg
    real(DP) :: d1
    real(DP) :: d2
    logical(LGP) :: stage_below_bot
    logical(LGP) :: head_below_bot
    !
    ! -- Set sfr density to sfr density or gwf density
    if (stage >= bots) then
      ss = stage
      stage_below_bot = .false.
      rdensesfr = this%denseterms(1, n) ! sfr rel density
    else
      ss = bots
      stage_below_bot = .true.
      rdensesfr = this%denseterms(2, n) ! gwf rel density
    end if
    !
    ! -- set hh to head or bots
    if (head >= bots) then
      hh = head
      head_below_bot = .false.
      rdensegwf = this%denseterms(2, n) ! gwf rel density
    else
      hh = bots
      head_below_bot = .true.
      rdensegwf = this%denseterms(1, n) ! sfr rel density
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
      rdenseavg = DHALF * (rdensesfr + rdensegwf)
      !
      ! -- Add contribution of first density term:
      !      cond * (denseavg/denseref - 1) * (hgwf - hsfr)
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
        !      cond * (havg - elevavg) * (densegwf - densesfr) / denseref
        elevgwf = this%denseterms(3, n)
        elevsfr = bots
        elevavg = DHALF * (elevsfr + elevgwf)
        havg = DHALF * (hh + ss)
        d2 = cond * (havg - elevavg) * (rdensegwf - rdensesfr)
        gwfrhs = gwfrhs + d2
        flow = flow + d2
      end if
    end if
  end subroutine sfr_calculate_density_exchange

end module SfrModule
