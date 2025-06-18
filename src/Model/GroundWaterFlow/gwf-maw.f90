module MawModule
  !
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENBOUNDNAME, LENTIMESERIESNAME, &
                             LENBUDTXT, DZERO, DEM9, DEM6, DEM4, DEM2, DQUARTER, &
                             DHALF, DP7, DP9, DONE, DTWO, DPI, DTWOPI, DEIGHT, &
                             DHUNDRED, DEP20, NAMEDBOUNDFLAG, LENPACKAGENAME, &
                             LENAUXNAME, LENFTYPE, DHNOFLO, DHDRY, DNODATA, &
                             MAXCHARLEN, TABLEFT, TABCENTER, TABRIGHT, &
                             TABSTRING, TABUCSTRING, TABINTEGER, TABREAL
  use SmoothingModule, only: sQuadraticSaturation, sQSaturation, &
                             sQuadraticSaturationDerivative, &
                             sQSaturationDerivative, &
                             sQuadratic0sp, &
                             sQuadratic0spDerivative
  use BndModule, only: BndType
  use BudgetObjectModule, only: BudgetObjectType, budgetobject_cr
  use TableModule, only: TableType, table_cr
  use ObserveModule, only: ObserveType
  use ObsModule, only: ObsType
  use GeomUtilModule, only: get_node
  use InputOutputModule, only: URWORD, extract_idnum_or_bndname, &
                               GetUnit, openfile
  use BaseDisModule, only: DisBaseType
  use SimModule, only: count_errors, store_error, store_error_unit, &
                       store_warning
  use BlockParserModule, only: BlockParserType
  use SimVariablesModule, only: errmsg, warnmsg
  use MemoryManagerModule, only: mem_allocate, mem_reallocate, mem_setptr, &
                                 mem_deallocate
  use MemoryHelperModule, only: create_mem_path
  use MatrixBaseModule
  !
  implicit none

  public :: MawType

  !
  character(len=LENFTYPE) :: ftype = 'MAW'
  character(len=LENPACKAGENAME) :: text = '             MAW'

  private
  public :: maw_create
  !
  type, extends(BndType) :: MawType
    !
    ! -- scalars
    ! -- characters
    !
    character(len=LENBUDTXT), dimension(:), pointer, &
      contiguous :: cmawbudget => NULL()
    character(len=LENAUXNAME), dimension(:), pointer, &
      contiguous :: cauxcbc => NULL()
    !
    ! -- logical
    logical(LGP), pointer :: correct_flow => null()
    !
    ! -- integers
    integer(I4B), pointer :: iprhed => null()
    integer(I4B), pointer :: iheadout => null()
    integer(I4B), pointer :: ibudgetout => null()
    integer(I4B), pointer :: ibudcsv => null()
    integer(I4B), pointer :: cbcauxitems => NULL()
    integer(I4B), pointer :: iflowingwells => NULL()
    integer(I4B), pointer :: imawiss => NULL()
    integer(I4B), pointer :: imawissopt => NULL()
    integer(I4B), pointer :: nmawwells => NULL()
    integer(I4B), pointer :: check_attr => NULL()
    integer(I4B), pointer :: ishutoffcnt => NULL()
    integer(I4B), pointer :: ieffradopt => NULL()
    integer(I4B), pointer :: ioutredflowcsv => NULL() !< unit number for CSV output file containing MAWs with reduced extraction/injection rates
    real(DP), pointer :: satomega => null()
    !
    ! -- for underrelaxation of estimated well q if using shutoff
    real(DP), pointer :: theta => NULL()
    real(DP), pointer :: kappa => NULL()
    !
    ! -- vector data for each well
    character(len=8), dimension(:), pointer, contiguous :: status => null()
    integer(I4B), dimension(:), pointer, contiguous :: ngwfnodes => null()
    integer(I4B), dimension(:), pointer, contiguous :: ieqn => null()
    integer(I4B), dimension(:), pointer, contiguous :: ishutoff => null()
    integer(I4B), dimension(:), pointer, contiguous :: ifwdischarge => null()
    real(DP), dimension(:), pointer, contiguous :: strt => null()
    real(DP), dimension(:), pointer, contiguous :: radius => null()
    real(DP), dimension(:), pointer, contiguous :: area => null()
    real(DP), dimension(:), pointer, contiguous :: pumpelev => null()
    real(DP), dimension(:), pointer, contiguous :: bot => null()
    real(DP), dimension(:), pointer, contiguous :: ratesim => null()
    real(DP), dimension(:), pointer, contiguous :: reduction_length => null()
    real(DP), dimension(:), pointer, contiguous :: fwelev => null()
    real(DP), dimension(:), pointer, contiguous :: fwcond => null()
    real(DP), dimension(:), pointer, contiguous :: fwrlen => null()
    real(DP), dimension(:), pointer, contiguous :: fwcondsim => null()
    real(DP), dimension(:), pointer, contiguous :: xsto => null()
    real(DP), dimension(:), pointer, contiguous :: xoldsto => null()
    real(DP), dimension(:), pointer, contiguous :: shutoffmin => null()
    real(DP), dimension(:), pointer, contiguous :: shutoffmax => null()
    real(DP), dimension(:), pointer, contiguous :: shutofflevel => null()
    real(DP), dimension(:), pointer, contiguous :: shutoffweight => null()
    real(DP), dimension(:), pointer, contiguous :: shutoffdq => null()
    real(DP), dimension(:), pointer, contiguous :: shutoffqold => null()
    character(len=LENBOUNDNAME), dimension(:), pointer, &
      contiguous :: cmawname => null()
    !
    ! -- time-series aware data
    real(DP), dimension(:), pointer, contiguous :: rate => null()
    real(DP), dimension(:), pointer, contiguous :: well_head => null()
    real(DP), dimension(:, :), pointer, contiguous :: mauxvar => null()
    !
    ! -- ia vector for connections
    integer(I4B), dimension(:), pointer, contiguous :: iaconn => null()
    !
    ! -- vector data for each connections
    integer(I4B), dimension(:), pointer, contiguous :: gwfnodes => NULL()
    real(DP), dimension(:), pointer, contiguous :: sradius => NULL()
    real(DP), dimension(:), pointer, contiguous :: hk => NULL()
    real(DP), dimension(:), pointer, contiguous :: satcond => NULL()
    real(DP), dimension(:), pointer, contiguous :: simcond => NULL()
    real(DP), dimension(:), pointer, contiguous :: topscrn => NULL()
    real(DP), dimension(:), pointer, contiguous :: botscrn => NULL()
    !
    ! -- imap vector
    integer(I4B), dimension(:), pointer, contiguous :: imap => null()
    !
    ! -- maw output data
    real(DP), dimension(:), pointer, contiguous :: qauxcbc => null()
    real(DP), dimension(:), pointer, contiguous :: dbuff => null()
    real(DP), dimension(:), pointer, contiguous :: qleak => null()
    real(DP), dimension(:), pointer, contiguous :: qout => null()
    real(DP), dimension(:), pointer, contiguous :: qfw => null()
    real(DP), dimension(:), pointer, contiguous :: qsto => null()
    real(DP), dimension(:), pointer, contiguous :: qconst => null()
    !
    ! -- for budgets
    integer(I4B), pointer :: bditems => NULL()
    type(BudgetObjectType), pointer :: budobj => null()
    !
    ! -- table objects
    type(TableType), pointer :: headtab => null()
    !
    ! -- pointer to gwf iss, k11, k22.
    integer(I4B), pointer :: gwfiss => NULL()
    real(DP), dimension(:), pointer, contiguous :: gwfk11 => NULL()
    real(DP), dimension(:), pointer, contiguous :: gwfk22 => NULL()
    integer(I4B), pointer :: gwfik22 => NULL()
    real(DP), dimension(:), pointer, contiguous :: gwfsat => NULL()
    !
    ! -- arrays for handling the rows added to the solution matrix
    integer(I4B), dimension(:), pointer, contiguous :: idxlocnode => null() !map position in global rhs and x array of pack entry
    integer(I4B), dimension(:), pointer, contiguous :: idxdglo => null() !map position in global array of package diagonal row entries
    integer(I4B), dimension(:), pointer, contiguous :: idxoffdglo => null() !map position in global array of package off diagonal row entries
    integer(I4B), dimension(:), pointer, contiguous :: idxsymdglo => null() !map position in global array of package diagonal entries to model rows
    integer(I4B), dimension(:), pointer, contiguous :: idxsymoffdglo => null() !map position in global array of package off diagonal entries to model rows
    integer(I4B), dimension(:), pointer, contiguous :: iboundpak => null() !package ibound
    real(DP), dimension(:), pointer, contiguous :: xnewpak => null() !package x vector
    real(DP), dimension(:), pointer, contiguous :: xoldpak => null() !package xold vector
    !
    ! -- density variables
    integer(I4B), pointer :: idense
    real(DP), dimension(:, :), pointer, contiguous :: denseterms => null()
    !
    ! -- viscosity variables
    real(DP), dimension(:, :), pointer, contiguous :: viscratios => null() !< viscosity ratios (1: maw vsc ratio; 2: gwf vsc ratio)
    !
    ! -- type bound procedures

  contains

    procedure :: maw_allocate_scalars
    procedure :: maw_allocate_well_conn_arrays
    procedure :: maw_allocate_arrays
    procedure :: bnd_options => maw_read_options
    procedure :: read_dimensions => maw_read_dimensions
    procedure :: read_initial_attr => maw_read_initial_attr
    procedure :: set_pointers => maw_set_pointers
    procedure :: bnd_ac => maw_ac
    procedure :: bnd_mc => maw_mc
    procedure :: bnd_ar => maw_ar
    procedure :: bnd_rp => maw_rp
    procedure :: bnd_ad => maw_ad
    procedure :: bnd_cf => maw_cf
    procedure :: bnd_fc => maw_fc
    procedure :: bnd_fn => maw_fn
    procedure :: bnd_nur => maw_nur
    procedure :: bnd_cq => maw_cq
    procedure :: bnd_ot_model_flows => maw_ot_model_flows
    procedure :: bnd_ot_package_flows => maw_ot_package_flows
    procedure :: bnd_ot_dv => maw_ot_dv
    procedure :: bnd_ot_bdsummary => maw_ot_bdsummary
    procedure :: bnd_da => maw_da
    procedure :: define_listlabel
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => maw_obs_supported
    procedure, public :: bnd_df_obs => maw_df_obs
    procedure, public :: bnd_rp_obs => maw_rp_obs
    procedure, public :: bnd_bd_obs => maw_bd_obs
    ! -- private procedures
    procedure, private :: maw_read_wells
    procedure, private :: maw_read_well_connections
    procedure, private :: maw_check_attributes
    procedure, private :: maw_set_stressperiod
    procedure, private :: maw_set_attribute_error
    procedure, private :: maw_calculate_saturation
    procedure, private :: maw_calculate_satcond
    procedure, private :: maw_calculate_conn_terms
    procedure, private :: maw_calculate_wellq
    procedure, private :: maw_calculate_qpot
    procedure, private :: maw_cfupdate
    procedure, private :: get_jpos
    procedure, private :: get_gwfnode
    ! -- budget
    procedure, private :: maw_setup_budobj
    procedure, private :: maw_fill_budobj
    ! -- table
    procedure, private :: maw_setup_tableobj
    ! -- density
    procedure :: maw_activate_density
    procedure, private :: maw_calculate_density_exchange
    ! -- MAW reduced flow outputs
    procedure, private :: maw_redflow_csv_init
    procedure, private :: maw_redflow_csv_write
    ! -- viscosity
    procedure :: maw_activate_viscosity
  end type MawType

contains

!> @brief Create a New Multi-Aquifer Well (MAW) Package
!!
!! After creating the package object point bndobj to the new package
!<
  subroutine maw_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname)
    ! -- dummy
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    type(MawType), pointer :: mawobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (mawobj)
    packobj => mawobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype)
    packobj%text = text
    !
    ! -- allocate scalars
    call mawobj%maw_allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()
    !
    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 4
    packobj%iscloc = 0 ! not supported
    packobj%isadvpak = 1
    packobj%ictMemPath = create_mem_path(namemodel, 'NPF')
  end subroutine maw_create

  !> @brief Allocate scalar members
  !<
  subroutine maw_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(MawType), intent(inout) :: this
    !
    ! -- call standard BndType allocate scalars
    call this%BndType%allocate_scalars()
    !
    ! -- allocate the object and assign values to object variables
    call mem_allocate(this%correct_flow, 'CORRECT_FLOW', this%memoryPath)
    call mem_allocate(this%iprhed, 'IPRHED', this%memoryPath)
    call mem_allocate(this%iheadout, 'IHEADOUT', this%memoryPath)
    call mem_allocate(this%ibudgetout, 'IBUDGETOUT', this%memoryPath)
    call mem_allocate(this%ibudcsv, 'IBUDCSV', this%memoryPath)
    call mem_allocate(this%iflowingwells, 'IFLOWINGWELLS', this%memoryPath)
    call mem_allocate(this%imawiss, 'IMAWISS', this%memoryPath)
    call mem_allocate(this%imawissopt, 'IMAWISSOPT', this%memoryPath)
    call mem_allocate(this%nmawwells, 'NMAWWELLS', this%memoryPath)
    call mem_allocate(this%check_attr, 'CHECK_ATTR', this%memoryPath)
    call mem_allocate(this%ishutoffcnt, 'ISHUTOFFCNT', this%memoryPath)
    call mem_allocate(this%ieffradopt, 'IEFFRADOPT', this%memoryPath)
    call mem_allocate(this%ioutredflowcsv, 'IOUTREDFLOWCSV', this%memoryPath) !for writing reduced MAW flows to csv file
    call mem_allocate(this%satomega, 'SATOMEGA', this%memoryPath)
    call mem_allocate(this%bditems, 'BDITEMS', this%memoryPath)
    call mem_allocate(this%theta, 'THETA', this%memoryPath)
    call mem_allocate(this%kappa, 'KAPPA', this%memoryPath)
    call mem_allocate(this%cbcauxitems, 'CBCAUXITEMS', this%memoryPath)
    call mem_allocate(this%idense, 'IDENSE', this%memoryPath)
    !
    ! -- Set values
    this%correct_flow = .FALSE.
    this%nmawwells = 0
    this%iprhed = 0
    this%iheadout = 0
    this%ibudgetout = 0
    this%ibudcsv = 0
    this%iflowingwells = 0
    this%imawiss = 0
    this%imawissopt = 0
    this%ieffradopt = 0
    this%ioutredflowcsv = 0
    this%satomega = DZERO
    this%bditems = 8
    this%theta = DP7
    this%kappa = DEM4
    this%cbcauxitems = 1
    this%idense = 0
    this%ivsc = 0
  end subroutine maw_allocate_scalars

  !> @brief Allocate well arrays
  !<
  subroutine maw_allocate_well_conn_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(MawType), intent(inout) :: this
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: jj
    !
    ! -- allocate character array for budget text
    call mem_allocate(this%cmawbudget, LENBUDTXT, this%bditems, 'CMAWBUDGET', &
                      this%memoryPath)
    !
    !-- fill cmawbudget
    this%cmawbudget(1) = '             GWF'
    this%cmawbudget(2) = '            RATE'
    this%cmawbudget(3) = '         STORAGE'
    this%cmawbudget(4) = '        CONSTANT'
    this%cmawbudget(5) = '         FW-RATE'
    this%cmawbudget(6) = '        FROM-MVR'
    this%cmawbudget(7) = '     RATE-TO-MVR'
    this%cmawbudget(8) = '  FW-RATE-TO-MVR'
    !
    ! -- allocate character arrays
    call mem_allocate(this%cmawname, LENBOUNDNAME, this%nmawwells, 'CMAWNAME', &
                      this%memoryPath)
    call mem_allocate(this%status, 8, this%nmawwells, 'STATUS', this%memoryPath)
    !
    ! -- allocate well data pointers in memory manager
    call mem_allocate(this%ngwfnodes, this%nmawwells, 'NGWFNODES', &
                      this%memoryPath)
    call mem_allocate(this%ieqn, this%nmawwells, 'IEQN', this%memoryPath)
    call mem_allocate(this%ishutoff, this%nmawwells, 'ISHUTOFF', this%memoryPath)
    call mem_allocate(this%ifwdischarge, this%nmawwells, 'IFWDISCHARGE', &
                      this%memoryPath)
    call mem_allocate(this%strt, this%nmawwells, 'STRT', this%memoryPath)
    call mem_allocate(this%radius, this%nmawwells, 'RADIUS', this%memoryPath)
    call mem_allocate(this%area, this%nmawwells, 'AREA', this%memoryPath)
    call mem_allocate(this%pumpelev, this%nmawwells, 'PUMPELEV', this%memoryPath)
    call mem_allocate(this%bot, this%nmawwells, 'BOT', this%memoryPath)
    call mem_allocate(this%ratesim, this%nmawwells, 'RATESIM', this%memoryPath)
    call mem_allocate(this%reduction_length, this%nmawwells, 'REDUCTION_LENGTH', &
                      this%memoryPath)
    call mem_allocate(this%fwelev, this%nmawwells, 'FWELEV', this%memoryPath)
    call mem_allocate(this%fwcond, this%nmawwells, 'FWCONDS', this%memoryPath)
    call mem_allocate(this%fwrlen, this%nmawwells, 'FWRLEN', this%memoryPath)
    call mem_allocate(this%fwcondsim, this%nmawwells, 'FWCONDSIM', &
                      this%memoryPath)
    call mem_allocate(this%xsto, this%nmawwells, 'XSTO', this%memoryPath)
    call mem_allocate(this%xoldsto, this%nmawwells, 'XOLDSTO', this%memoryPath)
    call mem_allocate(this%shutoffmin, this%nmawwells, 'SHUTOFFMIN', &
                      this%memoryPath)
    call mem_allocate(this%shutoffmax, this%nmawwells, 'SHUTOFFMAX', &
                      this%memoryPath)
    call mem_allocate(this%shutofflevel, this%nmawwells, 'SHUTOFFLEVEL', &
                      this%memoryPath)
    call mem_allocate(this%shutoffweight, this%nmawwells, 'SHUTOFFWEIGHT', &
                      this%memoryPath)
    call mem_allocate(this%shutoffdq, this%nmawwells, 'SHUTOFFDQ', &
                      this%memoryPath)
    call mem_allocate(this%shutoffqold, this%nmawwells, 'SHUTOFFQOLD', &
                      this%memoryPath)
    !
    ! -- timeseries aware variables
    call mem_allocate(this%rate, this%nmawwells, 'RATE', this%memoryPath)
    call mem_allocate(this%well_head, this%nmawwells, 'WELL_HEAD', &
                      this%memoryPath)
    if (this%naux > 0) then
      jj = this%naux
    else
      jj = 1
    end if
    call mem_allocate(this%mauxvar, jj, this%nmawwells, 'MAUXVAR', &
                      this%memoryPath)
    !
    ! -- allocate and initialize dbuff
    if (this%iheadout > 0) then
      call mem_allocate(this%dbuff, this%nmawwells, 'DBUFF', this%memoryPath)
    else
      call mem_allocate(this%dbuff, 0, 'DBUFF', this%memoryPath)
    end if
    !
    ! -- allocate iaconn
    call mem_allocate(this%iaconn, this%nmawwells + 1, 'IACONN', this%memoryPath)
    !
    ! -- allocate imap
    call mem_allocate(this%imap, this%MAXBOUND, 'IMAP', this%memoryPath)
    !
    ! -- allocate connection data
    call mem_allocate(this%gwfnodes, this%maxbound, 'GWFNODES', this%memoryPath)
    call mem_allocate(this%sradius, this%maxbound, 'SRADIUS', this%memoryPath)
    call mem_allocate(this%hk, this%maxbound, 'HK', this%memoryPath)
    call mem_allocate(this%satcond, this%maxbound, 'SATCOND', this%memoryPath)
    call mem_allocate(this%simcond, this%maxbound, 'SIMCOND', this%memoryPath)
    call mem_allocate(this%topscrn, this%maxbound, 'TOPSCRN', this%memoryPath)
    call mem_allocate(this%botscrn, this%maxbound, 'BOTSCRN', this%memoryPath)
    !
    ! -- allocate qleak
    call mem_allocate(this%qleak, this%maxbound, 'QLEAK', this%memoryPath)
    !
    ! -- initialize well data
    do n = 1, this%nmawwells
      this%status(n) = 'ACTIVE'
      this%ngwfnodes(n) = 0
      this%ieqn(n) = 0
      this%ishutoff(n) = 0
      this%ifwdischarge(n) = 0
      this%strt(n) = DEP20
      this%radius(n) = DEP20
      this%area(n) = DZERO
      this%pumpelev(n) = DEP20
      this%bot(n) = DEP20
      this%ratesim(n) = DZERO
      this%reduction_length(n) = DEP20
      this%fwelev(n) = DZERO
      this%fwcond(n) = DZERO
      this%fwrlen(n) = DZERO
      this%fwcondsim(n) = DZERO
      this%xsto(n) = DZERO
      this%xoldsto(n) = DZERO
      this%shutoffmin(n) = DZERO
      this%shutoffmax(n) = DZERO
      this%shutofflevel(n) = DEP20
      this%shutoffweight(n) = DONE
      this%shutoffdq(n) = DONE
      this%shutoffqold(n) = DONE
      !
      ! -- timeseries aware variables
      this%rate(n) = DZERO
      this%well_head(n) = DZERO
      do jj = 1, max(1, this%naux)
        this%mauxvar(jj, n) = DZERO
      end do
      !
      ! -- dbuff
      if (this%iheadout > 0) then
        this%dbuff(n) = DZERO
      end if
    end do
    !
    ! -- initialize iaconn
    do n = 1, this%nmawwells + 1
      this%iaconn(n) = 0
    end do
    !
    ! -- allocate character array for budget text
    call mem_allocate(this%cauxcbc, LENAUXNAME, this%cbcauxitems, 'CAUXCBC', &
                      this%memoryPath)
    !
    ! -- allocate and initialize qauxcbc
    call mem_allocate(this%qauxcbc, this%cbcauxitems, 'QAUXCBC', this%memoryPath)
    do j = 1, this%cbcauxitems
      this%qauxcbc(j) = DZERO
    end do
    !
    ! -- allocate flowing well data
    if (this%iflowingwells /= 0) then
      call mem_allocate(this%qfw, this%nmawwells, 'QFW', this%memoryPath)
    else
      call mem_allocate(this%qfw, 1, 'QFW', this%memoryPath)
    end if
    call mem_allocate(this%qout, this%nmawwells, 'QOUT', this%memoryPath)
    call mem_allocate(this%qsto, this%nmawwells, 'QSTO', this%memoryPath)
    call mem_allocate(this%qconst, this%nmawwells, 'QCONST', this%memoryPath)
    !
    ! -- initialize flowing well, storage, and constant flow terms
    do n = 1, this%nmawwells
      if (this%iflowingwells > 0) then
        this%qfw(n) = DZERO
      end if
      this%qsto(n) = DZERO
      this%qconst(n) = DZERO
    end do
    !
    ! -- initialize connection data
    do j = 1, this%maxbound
      this%imap(j) = 0
      this%gwfnodes(j) = 0
      this%sradius(j) = DZERO
      this%hk(j) = DZERO
      this%satcond(j) = DZERO
      this%simcond(j) = DZERO
      this%topscrn(j) = DZERO
      this%botscrn(j) = DZERO
      this%qleak(j) = DZERO
    end do
    !
    ! -- allocate denseterms to size 0
    call mem_allocate(this%denseterms, 3, 0, 'DENSETERMS', this%memoryPath)
    !
    ! -- allocate viscratios to size 0
    call mem_allocate(this%viscratios, 2, 0, 'VISCRATIOS', this%memoryPath)
  end subroutine maw_allocate_well_conn_arrays

  !> @brief Allocate arrays
  !<
  subroutine maw_allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(MawType), intent(inout) :: this
    ! -- local
    !
    ! -- call standard BndType allocate scalars
    call this%BndType%allocate_arrays()
  end subroutine maw_allocate_arrays

  !> @brief Read the packagedata for this package
  !<
  subroutine maw_read_wells(this)
    use ConstantsModule, only: LINELENGTH
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(MawType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: text
    character(len=LINELENGTH) :: keyword
    character(len=LINELENGTH) :: cstr
    character(len=LENBOUNDNAME) :: bndName
    character(len=LENBOUNDNAME) :: bndNameTemp
    character(len=9) :: cno
    logical :: isfound
    logical :: endOfBlock
    integer(I4B) :: ival
    integer(I4B) :: n
    integer(I4B) :: j
    integer(I4B) :: ii
    integer(I4B) :: jj
    integer(I4B) :: ieqn
    integer(I4B) :: itmp
    integer(I4B) :: ierr
    integer(I4B) :: idx
    real(DP) :: rval
    real(DP), pointer :: bndElem => null()
    ! -- local allocatable arrays
    character(len=LINELENGTH), dimension(:), allocatable :: strttext
    character(len=LENBOUNDNAME), dimension(:), allocatable :: nametxt
    character(len=50), dimension(:, :), allocatable :: caux
    integer(I4B), dimension(:), allocatable :: nboundchk
    integer(I4B), dimension(:), allocatable :: wellieqn
    integer(I4B), dimension(:), allocatable :: ngwfnodes
    real(DP), dimension(:), allocatable :: radius
    real(DP), dimension(:), allocatable :: bottom
    ! -- format
    character(len=*), parameter :: fmthdbot = &
      "('well head (', G0, ') must be greater than or equal to the &
      &BOTTOM_ELEVATION (', G0, ').')"
    !
    ! -- allocate and initialize temporary variables
    allocate (strttext(this%nmawwells))
    allocate (nametxt(this%nmawwells))
    if (this%naux > 0) then
      allocate (caux(this%naux, this%nmawwells))
    end if
    allocate (nboundchk(this%nmawwells))
    allocate (wellieqn(this%nmawwells))
    allocate (ngwfnodes(this%nmawwells))
    allocate (radius(this%nmawwells))
    allocate (bottom(this%nmawwells))
    !
    ! -- initialize temporary variables
    do n = 1, this%nmawwells
      nboundchk(n) = 0
    end do
    !
    ! -- initialize itmp
    itmp = 0
    !
    ! -- set npakeq to nmawwells
    this%npakeq = this%nmawwells
    !
    ! -- read maw well data
    ! -- get wells block
    call this%parser%GetBlock('PACKAGEDATA', isfound, ierr, &
                              supportopenclose=.true.)
    !
    ! -- parse locations block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') &
        'PROCESSING '//trim(adjustl(this%text))//' PACKAGEDATA'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        ival = this%parser%GetInteger()
        n = ival

        if (n < 1 .or. n > this%nmawwells) then
          write (errmsg, '(a,1x,i0,a)') &
            'IMAW must be greater than 0 and less than or equal to', &
            this%nmawwells, '.'
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- increment nboundchk
        nboundchk(n) = nboundchk(n) + 1
        !
        ! -- radius
        rval = this%parser%GetDouble()
        if (rval <= DZERO) then
          write (errmsg, '(a,1x,i0,1x,a)') &
            'Radius for well', n, 'must be greater than zero.'
          call store_error(errmsg)
        end if
        radius(n) = rval
        !
        ! -- well bottom
        bottom(n) = this%parser%GetDouble()
        !
        ! -- strt
        call this%parser%GetString(strttext(n))
        !
        ! -- ieqn
        call this%parser%GetStringCaps(keyword)
        if (keyword == 'SPECIFIED') then
          ieqn = 0
        else if (keyword == 'THIEM') then
          ieqn = 1
        else if (keyword == 'THEIM') then ! # codespell:ignore
          ieqn = 1
          write (warnmsg, '(a,a,a,a,a,a)') &
            "CONDEQN in '", trim(this%packName), "' should be ", &
            "corrected from '", trim(keyword), "' to 'THIEM'."
          call store_warning(warnmsg)
        else if (keyword == 'SKIN') then
          ieqn = 2
        else if (keyword == 'CUMULATIVE') then
          ieqn = 3
        else if (keyword == 'MEAN') then
          ieqn = 4
        else
          write (errmsg, '(a,1x,i0,1x,a)') &
            'CONDEQN for well', n, &
            "must be 'CUMULATIVE', 'THIEM', 'MEAN', or 'SKIN'."
        end if
        wellieqn(n) = ieqn
        !
        ! -- ngwnodes
        ival = this%parser%GetInteger()
        if (ival < 1) then
          ival = 0
          write (errmsg, '(a,1x,i0,1x,a)') &
            'NGWFNODES for well', n, 'must be greater than zero.'
          call store_error(errmsg)
        end if

        if (ival > 0) then
          ngwfnodes(n) = ival
        end if
        !
        ! -- increment maxbound
        itmp = itmp + ival
        !
        ! -- get aux data
        do jj = 1, this%naux
          call this%parser%GetString(caux(jj, n))
        end do
        !
        ! -- set default bndName
        write (cno, '(i9.9)') n
        bndName = 'MAWWELL'//cno
        !
        ! -- read well name
        if (this%inamedbound /= 0) then
          call this%parser%GetStringCaps(bndNameTemp)
          if (bndNameTemp /= '') then
            bndName = bndNameTemp
          end if
        end if
        nametxt(n) = bndName
      end do

      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' PACKAGEDATA'
      !
      ! -- check for duplicate or missing wells
      do n = 1, this%nmawwells
        if (nboundchk(n) == 0) then
          write (errmsg, '(a,1x,i0,a)') 'No data specified for maw well', n, '.'
          call store_error(errmsg)
        else if (nboundchk(n) > 1) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
            'Data for maw well', n, 'specified', nboundchk(n), 'times.'
          call store_error(errmsg)
        end if
      end do
    else
      call store_error('Required packagedata block not found.')
    end if
    !
    ! -- terminate if any errors were detected
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- set MAXBOUND
    this%maxbound = itmp
    write (this%iout, '(//4x,a,i7)') 'MAXBOUND = ', this%maxbound
    !
    ! -- allocate well and connection data
    call this%maw_allocate_well_conn_arrays()
    !
    ! -- fill well data with data stored in temporary local arrays
    do n = 1, this%nmawwells
      rval = radius(n)
      this%radius(n) = rval
      this%area(n) = DPI * rval**DTWO
      this%bot(n) = bottom(n)
      this%ieqn(n) = wellieqn(n)
      this%ngwfnodes(n) = ngwfnodes(n)
      this%cmawname(n) = nametxt(n)
      !
      ! fill timeseries aware data
      !
      ! -- well_head and strt
      jj = 1 ! For WELL_HEAD
      bndElem => this%well_head(n)
      call read_value_or_time_series_adv(strttext(n), n, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'WELL_HEAD')
      !
      ! -- set starting head value
      this%strt(n) = this%well_head(n)
      !
      ! -- check for error condition
      if (this%strt(n) < this%bot(n)) then
        write (cstr, fmthdbot) this%strt(n), this%bot(n)
        call this%maw_set_attribute_error(n, 'STRT', trim(cstr))
      end if
      !
      ! -- fill aux data
      do jj = 1, this%naux
        text = caux(jj, n)
        ii = n
        bndElem => this%mauxvar(jj, ii)
        call read_value_or_time_series_adv(text, ii, jj, bndElem, this%packName, &
                                           'AUX', this%tsManager, this%iprpak, &
                                           this%auxname(jj))
      end do
    end do
    !
    ! -- set iaconn and imap for each connection
    idx = 0
    this%iaconn(1) = 1
    do n = 1, this%nmawwells
      do j = 1, this%ngwfnodes(n)
        idx = idx + 1
        this%imap(idx) = n
      end do
      this%iaconn(n + 1) = idx + 1
    end do
    !
    ! -- deallocate local storage
    deallocate (strttext)
    deallocate (nametxt)
    if (this%naux > 0) then
      deallocate (caux)
    end if
    deallocate (nboundchk)
    deallocate (wellieqn)
    deallocate (ngwfnodes)
    deallocate (radius)
    deallocate (bottom)
  end subroutine maw_read_wells

  !> @brief Read the dimensions for this package
  !<
  subroutine maw_read_well_connections(this)
    use ConstantsModule, only: LINELENGTH
    ! -- dummy
    class(MawType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: cellid
    character(len=30) :: nodestr
    logical :: isfound
    logical :: endOfBlock
    integer(I4B) :: ierr
    integer(I4B) :: ival
    integer(I4B) :: j
    integer(I4B) :: jj
    integer(I4B) :: n
    integer(I4B) :: nn
    integer(I4B) :: nn2
    integer(I4B) :: ipos
    integer(I4B) :: jpos
    integer(I4B) :: ireset_scrntop
    integer(I4B) :: ireset_scrnbot
    integer(I4B) :: ireset_wellbot
    real(DP) :: rval
    real(DP) :: topnn
    real(DP) :: botnn
    real(DP) :: botw
    integer(I4B), dimension(:), pointer, contiguous :: nboundchk
    integer(I4B), dimension(:), pointer, contiguous :: iachk
    !
    ! -- initialize counters
    ireset_scrntop = 0
    ireset_scrnbot = 0
    ireset_wellbot = 0
    !
    ! -- allocate and initialize local storage
    allocate (iachk(this%nmawwells + 1))
    iachk(1) = 1
    do n = 1, this%nmawwells
      iachk(n + 1) = iachk(n) + this%ngwfnodes(n)
    end do
    allocate (nboundchk(this%maxbound))
    do n = 1, this%maxbound
      nboundchk(n) = 0
    end do
    !
    ! -- get well_connections block
    call this%parser%GetBlock('CONNECTIONDATA', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse well_connections block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text))// &
        ' CONNECTIONDATA'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        !
        ! -- well number
        ival = this%parser%GetInteger()
        n = ival
        !
        ! -- check for error condition
        if (n < 1 .or. n > this%nmawwells) then
          write (errmsg, '(a,1x,i0,a)') &
            'IMAW must be greater than 0 and less than or equal to ', &
            this%nmawwells, '.'
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- read connection number
        ival = this%parser%GetInteger()
        if (ival < 1 .or. ival > this%ngwfnodes(n)) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0,a)') &
            'JCONN for well ', n, &
            'must be greater than 1 and less than or equal to ', &
            this%ngwfnodes(n), '.'
          call store_error(errmsg)
          cycle
        end if

        ipos = iachk(n) + ival - 1
        nboundchk(ipos) = nboundchk(ipos) + 1

        j = ival
        jpos = this%get_jpos(n, ival)
        !
        ! -- read gwfnodes from the line
        call this%parser%GetCellid(this%dis%ndim, cellid)
        nn = this%dis%noder_from_cellid(cellid, this%inunit, this%iout)
        topnn = this%dis%top(nn)
        botnn = this%dis%bot(nn)
        botw = this%bot(n)
        !
        ! -- set gwf node number for connection
        this%gwfnodes(jpos) = nn
        !
        ! -- top of screen
        rval = this%parser%GetDouble()
        if (this%ieqn(n) /= 4) then
          rval = topnn
        else
          if (rval > topnn) then
            ireset_scrntop = ireset_scrntop + 1
            rval = topnn
          end if
        end if
        this%topscrn(jpos) = rval
        !
        ! -- bottom of screen
        rval = this%parser%GetDouble()
        if (this%ieqn(n) /= 4) then
          rval = botnn
        else
          if (rval < botnn) then
            ireset_scrnbot = ireset_scrnbot + 1
            rval = botnn
          end if
        end if
        this%botscrn(jpos) = rval
        !
        ! -- adjust the bottom of the well for all conductance approaches
        !    except for "mean"
        if (rval < botw) then
          if (this%ieqn(n) /= 4) then
            ireset_wellbot = ireset_wellbot + 1
            botw = rval
            this%bot(n) = rval
          else
            write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a,g0,a,g0,a)') &
              'Screen bottom for maw well', n, 'connection', j, '(', &
              this%botscrn(jpos), ') is less than the well bottom (', &
              this%bot(n), ').'
            call store_error(errmsg)
          end if
        end if
        !
        ! -- hydraulic conductivity or conductance
        rval = this%parser%GetDouble()
        if (this%ieqn(n) == 0) then
          this%satcond(jpos) = rval
        else if (this%ieqn(n) == 2 .OR. this%ieqn(n) == 3 .OR. &
                 this%ieqn(n) == 4) then
          this%hk(jpos) = rval
        end if
        !
        ! -- skin radius
        rval = this%parser%GetDouble()
        if (this%ieqn(n) == 2 .OR. this%ieqn(n) == 3 .OR. &
            this%ieqn(n) == 4) then
          this%sradius(jpos) = rval
          if (this%sradius(jpos) <= this%radius(n)) then
            write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a,g0,a,g0,a)') &
              'Screen radius for maw well', n, 'connection', j, '(', &
              this%sradius(jpos), &
              ') is less than or equal to the well radius (', &
              this%radius(n), ').'
            call store_error(errmsg)
          end if
        end if
      end do
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' CONNECTIONDATA'

      ipos = 0
      do n = 1, this%nmawwells
        do j = 1, this%ngwfnodes(n)
          ipos = ipos + 1
          !
          ! -- check for missing or duplicate maw well connections
          if (nboundchk(ipos) == 0) then
            write (errmsg, '(a,1x,i0,1x,a,1x,i0,a)') &
              'No data specified for maw well', n, 'connection', j, '.'
            call store_error(errmsg)
          else if (nboundchk(ipos) > 1) then
            write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a,1x,i0,1x,a)') &
              'Data for maw well', n, 'connection', j, &
              'specified', nboundchk(n), 'times.'
            call store_error(errmsg)
          end if
        end do
      end do
      !
      ! -- make sure that more than one connection per cell is only specified
      !    wells using the mean conducance type
      do n = 1, this%nmawwells
        if (this%ieqn(n) /= 4) then
          do j = 1, this%ngwfnodes(n)
            nn = this%get_gwfnode(n, j)
            do jj = 1, this%ngwfnodes(n)
              !
              ! -- skip current maw node
              if (jj == j) then
                cycle
              end if
              nn2 = this%get_gwfnode(n, jj)
              if (nn2 == nn) then
                call this%dis%noder_to_string(nn, nodestr)
                write (errmsg, '(a,1x,i0,1x,a,1x,i0,3(1x,a))') &
                  'Only one connection can be specified for maw well', &
                  n, 'connection', j, 'to gwf cell', trim(adjustl(nodestr)), &
                  'unless the mean condeqn is specified.'
                call store_error(errmsg)
              end if
            end do
          end do
        end if
      end do
    else
      call store_error('Required connectiondata block not found.')
    end if
    !
    ! -- deallocate local variable
    deallocate (iachk)
    deallocate (nboundchk)
    !
    ! -- add warning messages
    if (ireset_scrntop > 0) then
      write (warnmsg, '(a,1x,a,1x,a,1x,i0,1x,a)') &
        'The screen tops in multi-aquifer well package', trim(this%packName), &
        'were reset to the top of the connected cell', ireset_scrntop, 'times.'
      call store_warning(warnmsg)
    end if
    if (ireset_scrnbot > 0) then
      write (warnmsg, '(a,1x,a,1x,a,1x,i0,1x,a)') &
        'The screen bottoms in multi-aquifer well package', trim(this%packName), &
        'were reset to the bottom of the connected cell', ireset_scrnbot, &
        'times.'
      call store_warning(warnmsg)
    end if
    if (ireset_wellbot > 0) then
      write (warnmsg, '(a,1x,a,1x,a,1x,i0,1x,a)') &
        'The well bottoms in multi-aquifer well package', trim(this%packName), &
        'were reset to the bottom of the connected cell', ireset_wellbot, &
        'times.'
      call store_warning(warnmsg)
    end if
    !
    ! -- write summary of maw well_connection error messages
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine maw_read_well_connections

  !> @brief Read the dimensions for this package
  !<
  subroutine maw_read_dimensions(this)
    use ConstantsModule, only: LINELENGTH
    ! -- dummy
    class(MawType), intent(inout) :: this
    ! -- local
    character(len=LENBOUNDNAME) :: keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    ! -- format
    !
    ! -- initialize dimensions to -1
    this%nmawwells = -1
    this%maxbound = -1
    !
    ! -- get dimensions block
    call this%parser%GetBlock('DIMENSIONS', isfound, ierr, &
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
        case ('NMAWWELLS')
          this%nmawwells = this%parser%GetInteger()
          write (this%iout, '(4x,a,i0)') 'NMAWWELLS = ', this%nmawwells
        case default
          write (errmsg, '(3a)') &
            'Unknown '//trim(this%text)//' dimension: ', trim(keyword), '.'
          call store_error(errmsg)
        end select
      end do
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' DIMENSIONS'
    else
      call store_error('Required dimensions block not found.', terminate=.TRUE.)
    end if
    !
    ! -- verify dimensions were set correctly
    if (this%nmawwells < 0) then
      write (errmsg, '(a)') &
        'NMAWWELLS was not specified or was specified incorrectly.'
      call store_error(errmsg)
    end if
    !
    ! -- stop if errors were encountered in the DIMENSIONS block
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- read wells block
    call this%maw_read_wells()
    !
    ! -- read well_connections block
    call this%maw_read_well_connections()
    !
    ! -- Call define_listlabel to construct the list label that is written
    !    when PRINT_INPUT option is used.
    call this%define_listlabel()
    !
    ! -- setup the budget object
    call this%maw_setup_budobj()
    !
    ! -- setup the head table object
    call this%maw_setup_tableobj()
  end subroutine maw_read_dimensions

  !> @brief Read the initial parameters for this package
  !<
  subroutine maw_read_initial_attr(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use MemoryManagerModule, only: mem_setptr
    ! -- dummy
    class(MawType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    integer(I4B) :: ntabcols
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: nn
    integer(I4B) :: jpos
    integer(I4B) :: inode
    integer(I4B) :: idx
    real(DP) :: k11
    real(DP) :: k22
    character(len=10), dimension(0:4) :: ccond
    character(len=30) :: nodestr
    ! -- data
    data ccond(0)/'SPECIFIED '/
    data ccond(1)/'THIEM     '/
    data ccond(2)/'SKIN      '/
    data ccond(3)/'CUMULATIVE'/
    data ccond(4)/'MEAN      '/
    ! -- format
    character(len=*), parameter :: fmtwelln = &
      "(1X,//43X,'MULTI-AQUIFER WELL DATA'&
      &/1X,109('-'),&
      &/1X,7(A10,1X),A16)"
    character(len=*), parameter :: fmtwelld = &
      &"(1X,I10,1X,4(G10.3,1X),I10,1X,A10,1X,A16)"
    character(len=*), parameter :: fmtline = &
      &"(1X,119('-'),//)"
    character(len=*), parameter :: fmtwellcn = &
      "(1X,//37X,'MULTI-AQUIFER WELL CONNECTION DATA'&
      &/1X,119('-'),&
      &/1X,2(A10,1X),A20,7(A10,1X))"
    character(len=*), parameter :: fmtwellcd = &
      &"(1X,2(I10,1X),A20,1X,2(G10.3,1X),2(A10,1X),3(G10.3,1X))"
    !
    ! -- initialize xnewpak
    do n = 1, this%nmawwells
      this%xnewpak(n) = this%strt(n)
      this%xsto(n) = this%strt(n)
    end do
    !
    ! -- initialize status (iboundpak) of maw wells to active
    do n = 1, this%nmawwells
      select case (this%status(n))
      case ('CONSTANT')
        this%iboundpak(n) = -1
      case ('INACTIVE')
        this%iboundpak(n) = 0
      case ('ACTIVE')
        this%iboundpak(n) = 1
      end select
    end do
    !
    ! -- set imap and boundname for each connection
    if (this%inamedbound /= 0) then
      idx = 0
      do n = 1, this%nmawwells
        do j = 1, this%ngwfnodes(n)
          idx = idx + 1
          this%boundname(idx) = this%cmawname(n)
          this%imap(idx) = n
        end do
      end do
    else
      do n = 1, this%nmawwells
        this%cmawname(n) = ''
      end do
    end if
    !
    ! -- copy boundname into boundname_cst
    call this%copy_boundname()
    !
    ! -- set pointer to gwf iss and gwf hk
    call mem_setptr(this%gwfiss, 'ISS', create_mem_path(this%name_model))
    call mem_setptr(this%gwfk11, 'K11', create_mem_path(this%name_model, 'NPF'))
    call mem_setptr(this%gwfk22, 'K22', create_mem_path(this%name_model, 'NPF'))
    call mem_setptr(this%gwfik22, 'IK22', create_mem_path(this%name_model, 'NPF'))
    call mem_setptr(this%gwfsat, 'SAT', create_mem_path(this%name_model, 'NPF'))
    !
    ! -- qa data
    call this%maw_check_attributes()
    !
    ! -- Calculate the saturated conductance
    do n = 1, this%nmawwells
      !
      ! -- calculate saturated conductance only if CONDUCTANCE was not
      !    specified for each maw-gwf connection (CONDUCTANCE keyword).
      do j = 1, this%ngwfnodes(n)
        if (this%ieqn(n) /= 0) then
          inode = this%get_gwfnode(n, j)
          call this%maw_calculate_satcond(n, j, inode)
        end if
      end do
    end do
    !
    ! -- write summary of static well data
    ! -- write well data
    if (this%iprpak /= 0) then
      ntabcols = 7
      if (this%inamedbound /= 0) then
        ntabcols = ntabcols + 1
      end if
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(adjustl(this%packName))//') STATIC WELL DATA'
      call table_cr(this%inputtab, this%packName, title)
      call this%inputtab%table_df(this%nmawwells, ntabcols, this%iout)
      text = 'NUMBER'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'RADIUS'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'AREA'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'WELL BOTTOM'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'STARTING HEAD'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'NUMBER OF GWF NODES'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'CONDUCT. EQUATION'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      if (this%inamedbound /= 0) then
        text = 'NAME'
        call this%inputtab%initialize_column(text, 20, alignment=TABLEFT)
      end if
      do n = 1, this%nmawwells
        call this%inputtab%add_term(n)
        call this%inputtab%add_term(this%radius(n))
        call this%inputtab%add_term(this%area(n))
        call this%inputtab%add_term(this%bot(n))
        call this%inputtab%add_term(this%strt(n))
        call this%inputtab%add_term(this%ngwfnodes(n))
        call this%inputtab%add_term(ccond(this%ieqn(n)))
        if (this%inamedbound /= 0) then
          call this%inputtab%add_term(this%cmawname(n))
        end if
      end do
    end if
    !
    ! -- write well connection data
    if (this%iprpak /= 0) then
      ntabcols = 10
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(adjustl(this%packName))//') STATIC WELL CONNECTION DATA'
      call table_cr(this%inputtab, this%packName, title)
      call this%inputtab%table_df(this%maxbound, ntabcols, this%iout)
      text = 'NUMBER'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'WELL CONNECTION'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'CELLID'
      call this%inputtab%initialize_column(text, 20, alignment=TABLEFT)
      text = 'TOP OF SCREEN'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'BOTTOM OF SCREEN'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'SKIN RADIUS'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'SKIN K'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'K11'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'K22'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'SATURATED WELL CONDUCT.'
      call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
      !
      ! -- write the data to the table
      do n = 1, this%nmawwells
        do j = 1, this%ngwfnodes(n)
          call this%inputtab%add_term(n)
          call this%inputtab%add_term(j)
          jpos = this%get_jpos(n, j)
          nn = this%get_gwfnode(n, j)
          call this%dis%noder_to_string(nn, nodestr)
          call this%inputtab%add_term(nodestr)
          call this%inputtab%add_term(this%topscrn(jpos))
          call this%inputtab%add_term(this%botscrn(jpos))
          if (this%ieqn(n) == 2 .or. &
              this%ieqn(n) == 3 .or. &
              this%ieqn(n) == 4) then
            call this%inputtab%add_term(this%sradius(jpos))
            call this%inputtab%add_term(this%hk(jpos))
          else
            call this%inputtab%add_term(' ')
            call this%inputtab%add_term(' ')
          end if
          if (this%ieqn(n) == 1 .or. &
              this%ieqn(n) == 2 .or. &
              this%ieqn(n) == 3) then
            k11 = this%gwfk11(nn)
            if (this%gwfik22 == 0) then
              k22 = this%gwfk11(nn)
            else
              k22 = this%gwfk22(nn)
            end if
            call this%inputtab%add_term(k11)
            call this%inputtab%add_term(k22)
          else
            call this%inputtab%add_term(' ')
            call this%inputtab%add_term(' ')
          end if
          call this%inputtab%add_term(this%satcond(jpos))
        end do
      end do
    end if
    !
    ! -- finished with pointer to gwf hydraulic conductivity
    this%gwfk11 => null()
    this%gwfk22 => null()
    this%gwfik22 => null()
    this%gwfsat => null()
    !
    ! -- check for any error conditions
    if (count_errors() > 0) then
      call store_error_unit(this%inunit)
    end if
  end subroutine maw_read_initial_attr

  !> @brief Set a stress period attribute for mawweslls(imaw) using keywords
  !<
  subroutine maw_set_stressperiod(this, imaw, iheadlimit_warning)
    ! -- modules
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(MawType), intent(inout) :: this
    integer(I4B), intent(in) :: imaw
    integer(I4B), intent(inout) :: iheadlimit_warning
    ! -- local
    character(len=LINELENGTH) :: errmsgr
    character(len=LINELENGTH) :: text
    character(len=LINELENGTH) :: cstr
    character(len=LINELENGTH) :: caux
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ii
    integer(I4B) :: jj
    real(DP) :: rval
    real(DP), pointer :: bndElem => null()
    integer(I4B) :: istat
    ! -- formats
    character(len=*), parameter :: fmthdbot = &
      &"('well head (',G0,') must be >= BOTTOM_ELEVATION (',G0, ').')"
    !
    ! -- read remainder of variables on the line
    call this%parser%GetStringCaps(keyword)
    select case (keyword)
    case ('STATUS')
      call this%parser%GetStringCaps(text)
      this%status(imaw) = text(1:8)
      select case (text)
      case ('CONSTANT')
        this%iboundpak(imaw) = -1
      case ('INACTIVE')
        this%iboundpak(imaw) = 0
      case ('ACTIVE')
        this%iboundpak(imaw) = 1
      case default
        write (errmsg, '(2a)') &
          'Unknown '//trim(this%text)//" maw status keyword: '", &
          trim(text)//"'."
        call store_error(errmsg)
      end select
    case ('RATE')
      call this%parser%GetString(text)
      jj = 1 ! For RATE
      bndElem => this%rate(imaw)
      call read_value_or_time_series_adv(text, imaw, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'RATE')
    case ('WELL_HEAD')
      call this%parser%GetString(text)
      jj = 1 ! For WELL_HEAD
      bndElem => this%well_head(imaw)
      call read_value_or_time_series_adv(text, imaw, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'WELL_HEAD')
      !
      ! -- set xnewpak to well_head
      this%xnewpak(imaw) = this%well_head(imaw)
      !
      ! -- check for error condition
      if (this%well_head(imaw) < this%bot(imaw)) then
        write (cstr, fmthdbot) &
          this%well_head(imaw), this%bot(imaw)
        call this%maw_set_attribute_error(imaw, 'WELL HEAD', trim(cstr))
      end if
    case ('FLOWING_WELL')
      this%fwelev(imaw) = this%parser%GetDouble()
      this%fwcond(imaw) = this%parser%GetDouble()
      this%fwrlen(imaw) = this%parser%GetDouble()
      !
      ! -- test for condition where flowing well data is specified but
      !    flowing_wells is not specified in the options block
      if (this%iflowingwells == 0) then
        this%iflowingwells = -1
        text = 'Flowing well data is specified in the '//trim(this%packName)// &
               ' package but FLOWING_WELL was not specified in the '// &
               'OPTIONS block.'
        call store_warning(text)
      end if
    case ('RATE_SCALING')
      rval = this%parser%GetDouble()
      this%pumpelev(imaw) = rval
      rval = this%parser%GetDouble()
      this%reduction_length(imaw) = rval
      if (rval < DZERO) then
        call this%maw_set_attribute_error(imaw, trim(keyword), &
                                          'must be greater than or equal to 0.')
      end if
    case ('HEAD_LIMIT')
      call this%parser%GetString(text)
      if (trim(text) == 'OFF') then
        this%shutofflevel(imaw) = DEP20
      else
        read (text, *, iostat=istat, iomsg=errmsgr) &
          this%shutofflevel(imaw)
        if (istat /= 0) then
          errmsg = 'Could not read HEAD_LIMIT value. '//trim(errmsgr)
          call store_error(errmsg)
        end if
        if (this%shutofflevel(imaw) <= this%bot(imaw)) then
          iheadlimit_warning = iheadlimit_warning + 1
        end if
      end if
    case ('SHUT_OFF')
      rval = this%parser%GetDouble()
      this%shutoffmin(imaw) = rval
      rval = this%parser%GetDouble()
      this%shutoffmax(imaw) = rval
    case ('AUXILIARY')
      call this%parser%GetStringCaps(caux)
      do jj = 1, this%naux
        if (trim(adjustl(caux)) /= trim(adjustl(this%auxname(jj)))) cycle
        call this%parser%GetString(text)
        ii = imaw
        bndElem => this%mauxvar(jj, ii)
        call read_value_or_time_series_adv(text, ii, jj, bndElem, &
                                           this%packName, 'AUX', &
                                           this%tsManager, this%iprpak, &
                                           this%auxname(jj))
        exit
      end do
    case default
      write (errmsg, '(2a)') &
        'Unknown '//trim(this%text)//" maw data keyword: '", &
        trim(keyword)//"'."
      call store_error(errmsg)
    end select

  end subroutine maw_set_stressperiod

  !> @brief Issue a parameter error for mawweslls(imaw)
  !<
  subroutine maw_set_attribute_error(this, imaw, keyword, msg)
    use SimModule, only: store_error
    ! -- dummy
    class(MawType), intent(inout) :: this
    integer(I4B), intent(in) :: imaw
    character(len=*), intent(in) :: keyword
    character(len=*), intent(in) :: msg
    ! -- local
    ! -- formats
    !
    if (len(msg) == 0) then
      write (errmsg, '(a,1x,a,1x,i0,1x,a)') &
        keyword, ' for MAW well', imaw, 'has already been set.'
    else
      write (errmsg, '(a,1x,a,1x,i0,1x,a)') &
        keyword, ' for MAW well', imaw, msg
    end if
    call store_error(errmsg)
  end subroutine maw_set_attribute_error

  !> @brief Issue parameter errors for mawwells(imaw)
  !<
  subroutine maw_check_attributes(this)
    use SimModule, only: store_error
    ! -- dummy
    class(MawType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: cgwfnode
    integer(I4B) :: idx
    integer(I4B) :: n
    integer(I4B) :: j
    integer(I4B) :: jpos
    ! -- formats
    !
    idx = 1
    do n = 1, this%nmawwells
      if (this%ngwfnodes(n) < 1) then
        call this%maw_set_attribute_error(n, 'NGWFNODES', 'must be greater '// &
                                          'than 0.')
      end if
      if (this%radius(n) == DEP20) then
        call this%maw_set_attribute_error(n, 'RADIUS', 'has not been specified.')
      end if
      if (this%shutoffmin(n) > DZERO) then
        if (this%shutoffmin(n) >= this%shutoffmax(n)) then
          call this%maw_set_attribute_error(n, 'SHUT_OFF', 'shutoffmax must '// &
                                            'be greater than shutoffmin.')
        end if
      end if
      do j = 1, this%ngwfnodes(n)
        !
        ! -- calculate jpos
        jpos = this%get_jpos(n, j)
        !
        ! -- write gwfnode number
        write (cgwfnode, '(a,i0,a)') 'gwfnode(', j, ')'
        !
        ! -- connection screen data
        if (this%botscrn(jpos) >= this%topscrn(jpos)) then
          call this%maw_set_attribute_error(n, 'SCREEN_TOP', 'screen bottom '// &
                                            'must be less than screen top. '// &
                                            trim(cgwfnode))
        end if
        !
        ! -- connection skin hydraulic conductivity
        if (this%ieqn(n) == 2 .OR. this%ieqn(n) == 3 .OR. &
            this%ieqn(n) == 4) then
          if (this%hk(jpos) <= DZERO) then
            call this%maw_set_attribute_error(n, 'HK_SKIN', 'skin hyraulic '// &
                                              'conductivity must be greater '// &
                                              'than zero. '//trim(cgwfnode))
          end if
        else if (this%ieqn(n) == 0) then
          !
          ! -- saturated conductance
          if (this%satcond(jpos) < DZERO) then
            call this%maw_set_attribute_error(n, 'HK_SKIN', &
                                              'skin hyraulic conductivity '// &
                                              'must be greater than or '// &
                                              'equal to zero when using '// &
                                              'SPECIFIED condeqn. '// &
                                              trim(cgwfnode))
          end if
        end if
        idx = idx + 1
      end do
    end do
    ! -- reset check_attr
    this%check_attr = 0
  end subroutine maw_check_attributes

  !> @brief Add package connection to matrix
  !<
  subroutine maw_ac(this, moffset, sparse)
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(MawType), intent(inout) :: this
    integer(I4B), intent(in) :: moffset
    type(sparsematrix), intent(inout) :: sparse
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: jj
    integer(I4B) :: jglo
    integer(I4B) :: nglo
    ! -- format
    !
    ! -- Add package rows to sparse
    do n = 1, this%nmawwells
      nglo = moffset + this%dis%nodes + this%ioffset + n
      call sparse%addconnection(nglo, nglo, 1)
      do j = 1, this%ngwfnodes(n)
        jj = this%get_gwfnode(n, j)
        jglo = jj + moffset
        call sparse%addconnection(nglo, jglo, 1)
        call sparse%addconnection(jglo, nglo, 1)
      end do

    end do
  end subroutine maw_ac

  !> @brief Map package connection to matrix
  !<
  subroutine maw_mc(this, moffset, matrix_sln)
    use SparseModule, only: sparsematrix
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(MawType), intent(inout) :: this
    integer(I4B), intent(in) :: moffset
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: j
    integer(I4B) :: ii
    integer(I4B) :: iglo
    integer(I4B) :: jglo
    integer(I4B) :: ipos
    ! -- format
    !
    ! -- allocate connection mapping vectors
    call mem_allocate(this%idxlocnode, this%nmawwells, 'IDXLOCNODE', &
                      this%memoryPath)
    call mem_allocate(this%idxdglo, this%maxbound, 'IDXDGLO', this%memoryPath)
    call mem_allocate(this%idxoffdglo, this%maxbound, 'IDXOFFDGLO', &
                      this%memoryPath)
    call mem_allocate(this%idxsymdglo, this%maxbound, 'IDXSYMDGLO', &
                      this%memoryPath)
    call mem_allocate(this%idxsymoffdglo, this%maxbound, 'IDXSYMOFFDGLO', &
                      this%memoryPath)
    !
    ! -- Find the position of each connection in the global ia, ja structure
    !    and store them in idxglo.  idxglo allows this model to insert or
    !    retrieve values into or from the global A matrix
    ! -- maw rows
    ipos = 1
    do n = 1, this%nmawwells
      iglo = moffset + this%dis%nodes + this%ioffset + n
      this%idxlocnode(n) = this%dis%nodes + this%ioffset + n
      do ii = 1, this%ngwfnodes(n)
        j = this%get_gwfnode(n, ii)
        jglo = j + moffset
        this%idxdglo(ipos) = matrix_sln%get_position_diag(iglo)
        this%idxoffdglo(ipos) = matrix_sln%get_position(iglo, jglo)
        ipos = ipos + 1
      end do
    end do
    ! -- maw contributions gwf portion of global matrix
    ipos = 1
    do n = 1, this%nmawwells
      do ii = 1, this%ngwfnodes(n)
        iglo = this%get_gwfnode(n, ii) + moffset
        jglo = moffset + this%dis%nodes + this%ioffset + n
        this%idxsymdglo(ipos) = matrix_sln%get_position_diag(iglo)
        this%idxsymoffdglo(ipos) = matrix_sln%get_position(iglo, jglo)
        ipos = ipos + 1
      end do
    end do
  end subroutine maw_mc

  !> @brief Set options specific to MawType.
  !!
  !! Overrides BndType%bnd_options
  !<
  subroutine maw_read_options(this, option, found)
    use ConstantsModule, only: MAXCHARLEN, DZERO, MNORMAL
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: urword, assign_iounit, openfile
    ! -- dummy
    class(MawType), intent(inout) :: this
    character(len=*), intent(inout) :: option
    logical, intent(inout) :: found
    ! -- local
    character(len=MAXCHARLEN) :: fname, keyword
    ! -- formats
    character(len=*), parameter :: fmtflowingwells = &
      &"(4x, 'FLOWING WELLS WILL BE SIMULATED.')"
    character(len=*), parameter :: fmtshutdown = &
      &"(4x, 'SHUTDOWN ', a, ' VALUE (',g15.7,') SPECIFIED.')"
    character(len=*), parameter :: fmtnostoragewells = &
      &"(4x, 'WELL STORAGE WILL NOT BE SIMULATED.')"
    character(len=*), parameter :: fmtmawbin = &
      "(4x, 'MAW ', 1x, a, 1x, ' WILL BE SAVED TO FILE: ', a, /4x, &
     &'OPENED ON UNIT: ', I0)"
    !
    ! -- Check for 'FLOWING_WELLS' and set this%iflowingwells
    found = .true.
    select case (option)
    case ('PRINT_HEAD')
      this%iprhed = 1
      write (this%iout, '(4x,a)') &
        trim(adjustl(this%text))//' heads will be printed to listing file.'
    case ('HEAD')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        call assign_iounit(this%iheadout, this%inunit, "HEAD fileout")
        call openfile(this%iheadout, this%iout, fname, 'DATA(BINARY)', &
                      form, access, 'REPLACE', mode_opt=MNORMAL)
        write (this%iout, fmtmawbin) 'HEAD', trim(adjustl(fname)), &
          this%iheadout
      else
        call store_error('Optional maw stage keyword must be '// &
                         'followed by fileout.')
      end if
    case ('BUDGET')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        call assign_iounit(this%ibudgetout, this%inunit, "BUDGET fileout")
        call openfile(this%ibudgetout, this%iout, fname, 'DATA(BINARY)', &
                      form, access, 'REPLACE', mode_opt=MNORMAL)
        write (this%iout, fmtmawbin) 'BUDGET', trim(adjustl(fname)), &
          this%ibudgetout
      else
        call store_error('Optional maw budget keyword must be '// &
                         'followed by fileout.')
      end if
    case ('BUDGETCSV')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        call assign_iounit(this%ibudcsv, this%inunit, "BUDGETCSV fileout")
        call openfile(this%ibudcsv, this%iout, fname, 'CSV', &
                      filstat_opt='REPLACE')
        write (this%iout, fmtmawbin) 'BUDGET CSV', trim(adjustl(fname)), &
          this%ibudcsv
      else
        call store_error('OPTIONAL BUDGETCSV KEYWORD MUST BE FOLLOWED BY &
          &FILEOUT')
      end if
    case ('FLOWING_WELLS')
      this%iflowingwells = 1
      write (this%iout, fmtflowingwells)
    case ('SHUTDOWN_THETA')
      this%theta = this%parser%GetDouble()
      write (this%iout, fmtshutdown) 'THETA', this%theta
    case ('SHUTDOWN_KAPPA')
      this%kappa = this%parser%GetDouble()
      write (this%iout, fmtshutdown) 'KAPPA', this%kappa
    case ('MOVER')
      this%imover = 1
      write (this%iout, '(4x,A)') 'MOVER OPTION ENABLED'
    case ('NO_WELL_STORAGE')
      this%imawissopt = 1
      write (this%iout, fmtnostoragewells)
    case ('FLOW_CORRECTION')
      this%correct_flow = .TRUE.
      write (this%iout, '(4x,a,/,4x,a)') &
        'MAW-GWF FLOW CORRECTIONS WILL BE APPLIED WHEN MAW HEADS ARE BELOW', &
        'OR GWF HEADS IN CONNECTED CELLS ARE BELOW THE CELL BOTTOM.'
    case ('MAW_FLOW_REDUCE_CSV')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        call this%maw_redflow_csv_init(fname)
      else
        call store_error('OPTIONAL MAW_FLOW_REDUCE_CSV KEYWORD MUST BE &
          &FOLLOWED BY FILEOUT')
      end if
      !
      ! -- right now these are options that are only available in the
      !    development version and are not included in the documentation.
      !    These options are only available when IDEVELOPMODE in
      !    constants module is set to 1
    case ('DEV_PEACEMAN_EFFECTIVE_RADIUS')
      call this%parser%DevOpt()
      this%ieffradopt = 1
      write (this%iout, '(4x,a)') &
        'EFFECTIVE RADIUS FOR STRUCTURED GRIDS WILL BE CALCULATED &
        &USING PEACEMAN 1983'
    case default
      !
      ! -- No options found
      found = .false.
    end select
  end subroutine maw_read_options

  !> @brief Allocate and Read
  !!
  !! Create new MAW package and point bndobj to the new package
  !<
  subroutine maw_ar(this)
    ! -- dummy
    class(MawType), intent(inout) :: this
    ! -- local
    ! -- format
    !
    call this%obs%obs_ar()
    !
    ! -- set omega value used for saturation calculations
    if (this%inewton > 0) then
      this%satomega = DEM6
    end if
    !
    ! -- Allocate connection arrays in MAW and in package superclass
    call this%maw_allocate_arrays()
    !
    ! -- read optional initial package parameters
    call this%read_initial_attr()
    !
    ! -- setup pakmvrobj
    if (this%imover /= 0) then
      allocate (this%pakmvrobj)
      call this%pakmvrobj%ar(this%nmawwells, this%nmawwells, this%memoryPath)
    end if
  end subroutine maw_ar

  !> @brief Read and Prepare
  !!
  !! Read itmp and new boundaries if itmp > 0
  !<
  subroutine maw_rp(this)
    use ConstantsModule, only: LINELENGTH
    use TdisModule, only: kper, nper
    ! -- dummy
    class(MawType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: line
    character(len=LINELENGTH) :: text
    character(len=16) :: csteady
    logical :: isfound
    logical :: endOfBlock
    integer(I4B) :: ierr
    integer(I4B) :: node
    integer(I4B) :: n
    integer(I4B) :: ntabcols
    integer(I4B) :: ntabrows
    integer(I4B) :: imaw
    integer(I4B) :: ibnd
    integer(I4B) :: j
    integer(I4B) :: jpos
    integer(I4B) :: iheadlimit_warning
    ! -- formats
    character(len=*), parameter :: fmtblkerr = &
      &"('Looking for BEGIN PERIOD iper.  Found ', a, ' instead.')"
    character(len=*), parameter :: fmtlsp = &
      &"(1X,/1X,'REUSING ',A,'S FROM LAST STRESS PERIOD')"
    !
    ! -- initialize counters
    iheadlimit_warning = 0
    !
    ! -- set steady-state flag based on gwfiss
    this%imawiss = this%gwfiss
    !
    ! -- reset maw steady flag if 'STEADY-STATE' specified in the OPTIONS block
    if (this%imawissopt == 1) then
      this%imawiss = 1
    end if
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
          call store_error(errmsg, terminate=.TRUE.)
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
        call this%inputtab%table_df(1, 5, this%iout, finalize=.FALSE.)
        text = 'NUMBER'
        call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
        text = 'KEYWORD'
        call this%inputtab%initialize_column(text, 20, alignment=TABLEFT)
        do n = 1, 3
          write (text, '(a,1x,i6)') 'VALUE', n
          call this%inputtab%initialize_column(text, 15, alignment=TABCENTER)
        end do
      end if
      !
      ! -- set flag to check attributes
      this%check_attr = 1
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit

        imaw = this%parser%GetInteger()
        if (imaw < 1 .or. imaw > this%nmawwells) then
          write (errmsg, '(2(a,1x),i0,a)') &
            'IMAW must be greater than 0 and', &
            'less than or equal to ', this%nmawwells, '.'
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- set stress period data
        call this%maw_set_stressperiod(imaw, iheadlimit_warning)
        !
        ! -- write line to table
        if (this%iprpak /= 0) then
          call this%parser%GetCurrentLine(line)
          call this%inputtab%line_to_columns(line)
        end if
      end do
      if (this%iprpak /= 0) then
        call this%inputtab%finalize_table()
      end if
      !
      ! -- using data from the last stress period
    else
      write (this%iout, fmtlsp) trim(this%filtyp)
    end if
    !
    ! -- issue warning messages
    if (iheadlimit_warning > 0) then
      write (warnmsg, '(a,a,a,1x,a,1x,a)') &
        "HEAD_LIMIT in '", trim(this%packName), "' was below the well bottom", &
        "for one or more multi-aquifer well(s). This may result in", &
        "convergence failures for some models."
      call store_warning(warnmsg, substring=warnmsg(:50))
    end if
    !
    ! -- write summary of maw well stress period error messages
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- qa data if necessary
    if (this%check_attr /= 0) then
      call this%maw_check_attributes()

      ! -- write summary of stress period data for MAW
      if (this%iprpak == 1) then
        if (this%imawiss /= 0) then
          csteady = 'STEADY-STATE    '
        else
          csteady = 'TRANSIENT       '
        end if
        !
        ! -- reset the input table object for rate data
        title = trim(adjustl(this%text))//' PACKAGE ('// &
                trim(adjustl(this%packName))//') '//trim(adjustl(csteady))// &
                ' RATE DATA FOR PERIOD'
        write (title, '(a,1x,i6)') trim(adjustl(title)), kper
        ntabcols = 6
        call table_cr(this%inputtab, this%packName, title)
        call this%inputtab%table_df(this%nmawwells, ntabcols, this%iout)
        text = 'NUMBER'
        call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
        text = 'STATUS'
        call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
        text = 'RATE'
        call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
        text = 'SPECIFIED HEAD'
        call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
        text = 'PUMP ELEVATION'
        call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
        text = 'REDUCTION LENGTH'
        call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
        do n = 1, this%nmawwells
          call this%inputtab%add_term(n)
          call this%inputtab%add_term(this%status(n))
          call this%inputtab%add_term(this%rate(n))
          if (this%iboundpak(n) < 0) then
            call this%inputtab%add_term(this%well_head(n))
          else
            call this%inputtab%add_term(' ')
          end if
          call this%inputtab%add_term(this%pumpelev(n))
          if (this%reduction_length(n) /= DEP20) then
            call this%inputtab%add_term(this%reduction_length(n))
          else
            call this%inputtab%add_term(' ')
          end if
        end do
        !
        ! -- flowing wells
        if (this%iflowingwells > 0) then
          !
          ! -- reset the input table object for flowing well data
          title = trim(adjustl(this%text))//' PACKAGE ('// &
                  trim(adjustl(this%packName))//') '//trim(adjustl(csteady))// &
                  ' FLOWING WELL DATA FOR PERIOD'
          write (title, '(a,1x,i6)') trim(adjustl(title)), kper
          ntabcols = 4
          ntabrows = 0
          do n = 1, this%nmawwells
            if (this%fwcond(n) > DZERO) then
              ntabrows = ntabrows + 1
            end if
          end do
          if (ntabrows > 0) then
            call table_cr(this%inputtab, this%packName, title)
            call this%inputtab%table_df(ntabrows, ntabcols, this%iout)
            text = 'NUMBER'
            call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
            text = 'ELEVATION'
            call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
            text = 'CONDUCT.'
            call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
            text = 'REDUCTION LENGTH'
            call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
            do n = 1, this%nmawwells
              if (this%fwcond(n) > DZERO) then
                call this%inputtab%add_term(n)
                call this%inputtab%add_term(this%fwelev(n))
                call this%inputtab%add_term(this%fwcond(n))
                call this%inputtab%add_term(this%fwrlen(n))
              end if
            end do
          end if
        end if
        !
        ! -- reset the input table object for shutoff data
        title = trim(adjustl(this%text))//' PACKAGE ('// &
                trim(adjustl(this%packName))//') '//trim(adjustl(csteady))// &
                ' WELL SHUTOFF DATA FOR PERIOD'
        write (title, '(a,1x,i6)') trim(adjustl(title)), kper
        ntabcols = 4
        ntabrows = 0
        do n = 1, this%nmawwells
          if (this%shutofflevel(n) /= DEP20) then
            ntabrows = ntabrows + 1
          end if
        end do
        if (ntabrows > 0) then
          call table_cr(this%inputtab, this%packName, title)
          call this%inputtab%table_df(ntabrows, ntabcols, this%iout)
          text = 'NUMBER'
          call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
          text = 'ELEVATION'
          call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
          text = 'MINIMUM. Q'
          call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
          text = 'MAXIMUM Q'
          call this%inputtab%initialize_column(text, 12, alignment=TABCENTER)
          do n = 1, this%nmawwells
            if (this%shutofflevel(n) /= DEP20) then
              call this%inputtab%add_term(n)
              call this%inputtab%add_term(this%shutofflevel(n))
              call this%inputtab%add_term(this%shutoffmin(n))
              call this%inputtab%add_term(this%shutoffmax(n))
            end if
          end do
        end if
      end if
    end if
    !
    ! -- fill arrays
    ibnd = 1
    do n = 1, this%nmawwells
      do j = 1, this%ngwfnodes(n)
        jpos = this%get_jpos(n, j)
        node = this%get_gwfnode(n, j)
        this%nodelist(ibnd) = node
        this%bound(1, ibnd) = this%xnewpak(n)
        this%bound(2, ibnd) = this%satcond(jpos)
        this%bound(3, ibnd) = this%botscrn(jpos)
        if (this%iboundpak(n) > 0) then
          this%bound(4, ibnd) = this%rate(n)
        else
          this%bound(4, ibnd) = DZERO
        end if
        ibnd = ibnd + 1
      end do
    end do
  end subroutine maw_rp

  !> @brief Add package connection to matrix
  !<
  subroutine maw_ad(this)
    use TdisModule, only: kper, kstp
    ! -- dummy
    class(MawType) :: this
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: j
    integer(I4B) :: jj
    integer(I4B) :: ibnd
    !
    ! -- Advance the time series
    call this%TsManager%ad()
    !
    ! -- update auxiliary variables by copying from the derived-type time
    !    series variable into the bndpackage auxvar variable so that this
    !    information is properly written to the GWF budget file
    if (this%naux > 0) then
      ibnd = 1
      do n = 1, this%nmawwells
        do j = 1, this%ngwfnodes(n)
          do jj = 1, this%naux
            if (this%noupdateauxvar(jj) /= 0) cycle
            this%auxvar(jj, ibnd) = this%mauxvar(jj, n)
          end do
          ibnd = ibnd + 1
        end do
      end do
    end if
    !
    ! -- copy xnew into xold
    do n = 1, this%nmawwells
      this%xoldpak(n) = this%xnewpak(n)
      this%xoldsto(n) = this%xsto(n)
      if (this%iboundpak(n) < 0) then
        this%xnewpak(n) = this%well_head(n)
      end if
    end do
    !
    !--use the appropriate xoldsto if initial heads are above the
    !  specified flowing well discharge elevation
    if (kper == 1 .and. kstp == 1) then
      do n = 1, this%nmawwells
        if (this%fwcond(n) > DZERO) then
          if (this%xoldsto(n) > this%fwelev(n)) then
            this%xoldsto(n) = this%fwelev(n)
          end if
        end if
      end do
    end if
    !
    ! -- reset ishutoffcnt (equivalent to kiter) to zero
    this%ishutoffcnt = 0
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
  end subroutine maw_ad

  !> @brief Formulate the HCOF and RHS terms
  !!
  !! Skip if no multi-aquifer wells, otherwise, calculate hcof and rhs
  !<
  subroutine maw_cf(this)
    ! -- dummy
    class(MawType) :: this
    ! -- local
    !
    ! -- Calculate maw conductance and update package RHS and HCOF
    call this%maw_cfupdate()
  end subroutine maw_cf

  !> @brief Copy rhs and hcof into solution rhs and amat
  !<
  subroutine maw_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(MawType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: idx
    integer(I4B) :: iloc
    integer(I4B) :: isymloc
    integer(I4B) :: igwfnode
    integer(I4B) :: iposd
    integer(I4B) :: iposoffd
    integer(I4B) :: isymnode
    integer(I4B) :: ipossymd
    integer(I4B) :: ipossymoffd
    integer(I4B) :: jpos
    integer(I4B) :: icflow
    real(DP) :: hmaw
    real(DP) :: hgwf
    real(DP) :: cfw
    real(DP) :: cmaw
    real(DP) :: cterm
    real(DP) :: term
    real(DP) :: scale
    real(DP) :: tp
    real(DP) :: bt
    real(DP) :: rate
    real(DP) :: ratefw
    real(DP) :: flow
    !
    ! -- pakmvrobj fc
    if (this%imover == 1) then
      call this%pakmvrobj%fc()
    end if
    !
    ! -- Copy package rhs and hcof into solution rhs and amat
    idx = 1
    do n = 1, this%nmawwells
      iloc = this%idxlocnode(n)
      !
      ! -- update head value for constant head maw wells
      if (this%iboundpak(n) < 0) then
        this%xnewpak(n) = this%well_head(n)
      end if
      hmaw = this%xnewpak(n)
      !
      ! -- add pumping rate to active or constant maw well
      if (this%iboundpak(n) == 0) then
        this%ratesim(n) = DZERO
      else
        call this%maw_calculate_wellq(n, hmaw, rate)
        this%ratesim(n) = rate
        rhs(iloc) = rhs(iloc) - rate
        !
        ! -- location of diagonal for maw row
        iposd = this%idxdglo(idx)
        !
        ! -- add flowing well
        this%xsto(n) = hmaw
        ratefw = DZERO
        if (this%iflowingwells > 0) then
          if (this%fwcond(n) > DZERO) then
            bt = this%fwelev(n)
            tp = bt + this%fwrlen(n)
            scale = sQSaturation(tp, bt, hmaw)
            cfw = scale * this%fwcond(n)
            this%ifwdischarge(n) = 0
            if (cfw > DZERO) then
              this%ifwdischarge(n) = 1
              this%xsto(n) = bt
            end if
            this%fwcondsim(n) = cfw
            call matrix_sln%add_value_pos(iposd, -cfw)
            rhs(iloc) = rhs(iloc) - cfw * bt
            ratefw = cfw * (bt - hmaw)
          end if
        end if
        !
        ! -- add maw storage changes
        if (this%imawiss /= 1) then
          if (this%ifwdischarge(n) /= 1) then
            call matrix_sln%add_value_pos(iposd, -this%area(n) / delt)
            rhs(iloc) = rhs(iloc) - (this%area(n) * this%xoldsto(n) / delt)
          else
            cterm = this%xoldsto(n) - this%fwelev(n)
            rhs(iloc) = rhs(iloc) - (this%area(n) * cterm / delt)
          end if
        end if
        !
        ! -- If mover is active, add receiver water to rhs and
        !    store available water (as positive value)
        if (this%imover == 1) then
          rhs(iloc) = rhs(iloc) - this%pakmvrobj%get_qfrommvr(n)
          !
          ! -- add pumping rate to mover if not injection
          if (rate < 0) then
            call this%pakmvrobj%accumulate_qformvr(n, -rate) !pumped water
          end if
          !
          ! -- add flowing well flow to mover
          call this%pakmvrobj%accumulate_qformvr(n, -ratefw) !flowing water
        end if
        !
      end if
      !
      ! -- process each maw/gwf connection
      do j = 1, this%ngwfnodes(n)
        if (this%iboundpak(n) /= 0) then
          jpos = this%get_jpos(n, j)
          igwfnode = this%get_gwfnode(n, j)
          hgwf = this%xnew(igwfnode)
          !
          ! -- calculate connection terms
          call this%maw_calculate_conn_terms(n, j, icflow, cmaw, cterm, term, &
                                             flow)
          this%simcond(jpos) = cmaw
          !
          ! -- add to maw row
          iposd = this%idxdglo(idx)
          iposoffd = this%idxoffdglo(idx)
          call matrix_sln%add_value_pos(iposd, -term)
          call matrix_sln%set_value_pos(iposoffd, term)
          !
          ! -- add correction term
          rhs(iloc) = rhs(iloc) - cterm
          !
          ! -- add to gwf row for maw connection
          isymnode = this%get_gwfnode(n, j)
          isymloc = ia(isymnode)
          ipossymd = this%idxsymdglo(idx)
          ipossymoffd = this%idxsymoffdglo(idx)
          call matrix_sln%add_value_pos(ipossymd, -term)
          call matrix_sln%set_value_pos(ipossymoffd, term)
          !
          ! -- add correction term to gwf row
          rhs(isymnode) = rhs(isymnode) + cterm
        end if
        !
        ! -- increment maw connection counter
        idx = idx + 1
      end do
    end do
  end subroutine maw_fc

  !> @brief Fill newton terms
  !<
  subroutine maw_fn(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(MawType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: idx
    integer(I4B) :: iloc
    integer(I4B) :: isymloc
    integer(I4B) :: igwfnode
    integer(I4B) :: iposd
    integer(I4B) :: iposoffd
    integer(I4B) :: isymnode
    integer(I4B) :: ipossymd
    integer(I4B) :: ipossymoffd
    integer(I4B) :: jpos
    integer(I4B) :: icflow
    real(DP) :: hmaw
    real(DP) :: hgwf
    real(DP) :: scale
    real(DP) :: tp
    real(DP) :: bt
    real(DP) :: cfw
    real(DP) :: rate
    real(DP) :: rate2
    real(DP) :: rterm
    real(DP) :: derv
    real(DP) :: drterm
    real(DP) :: cmaw
    real(DP) :: cterm
    real(DP) :: term
    real(DP) :: flow
    real(DP) :: term2
    real(DP) :: rhsterm
    !
    ! -- Calculate Newton-Raphson corrections
    idx = 1
    do n = 1, this%nmawwells
      iloc = this%idxlocnode(n)
      hmaw = this%xnewpak(n)
      !
      ! -- add pumping rate to active or constant maw well
      if (this%iboundpak(n) /= 0) then
        iposd = this%idxdglo(idx)
        scale = DONE
        drterm = DZERO
        rate = this%ratesim(n)
        !
        !-- calculate final derivative for pumping rate
        call this%maw_calculate_wellq(n, hmaw + DEM4, rate2)
        drterm = (rate2 - rate) / DEM4
        !
        !-- fill amat and rhs with newton-raphson terms
        call matrix_sln%add_value_pos(iposd, drterm)
        rhs(iloc) = rhs(iloc) + drterm * hmaw
        !
        ! -- add flowing well
        if (this%iflowingwells > 0) then
          if (this%fwcond(n) > DZERO) then
            bt = this%fwelev(n)
            tp = bt + this%fwrlen(n)
            scale = sQSaturation(tp, bt, hmaw)
            cfw = scale * this%fwcond(n)
            this%ifwdischarge(n) = 0
            if (cfw > DZERO) then
              this%ifwdischarge(n) = 1
            end if
            this%fwcondsim(n) = cfw
            rate = cfw * (bt - hmaw)
            rterm = -cfw * hmaw
            !
            ! --calculate derivative for flowing well
            if (hmaw < tp) then
              derv = sQSaturationDerivative(tp, bt, hmaw)
              drterm = -(cfw + this%fwcond(n) * derv * (hmaw - bt))
              !
              ! -- fill amat and rhs with newton-raphson terms
              call matrix_sln%add_value_pos(iposd, &
                                            -this%fwcond(n) * derv * (hmaw - bt))
              rhs(iloc) = rhs(iloc) - rterm + drterm * hmaw
            end if
          end if
        end if
      end if
      !
      ! -- process each maw/gwf connection
      do j = 1, this%ngwfnodes(n)
        if (this%iboundpak(n) /= 0) then
          jpos = this%get_jpos(n, j)
          igwfnode = this%get_gwfnode(n, j)
          hgwf = this%xnew(igwfnode)
          !
          ! -- add to maw row
          iposd = this%idxdglo(idx)
          iposoffd = this%idxoffdglo(idx)
          !
          ! -- add to gwf row for maw connection
          isymnode = this%get_gwfnode(n, j)
          isymloc = ia(isymnode)
          ipossymd = this%idxsymdglo(idx)
          ipossymoffd = this%idxsymoffdglo(idx)
          !
          ! -- calculate newton terms
          call this%maw_calculate_conn_terms(n, j, icflow, cmaw, cterm, term, &
                                             flow, term2)
          !
          ! -- maw is upstream
          if (hmaw > hgwf) then
            if (icflow /= 0) then
              rhsterm = term2 * hgwf + term * hmaw
              rhs(iloc) = rhs(iloc) + rhsterm
              rhs(isymnode) = rhs(isymnode) - rhsterm
              if (this%iboundpak(n) > 0) then
                call matrix_sln%add_value_pos(iposd, term)
                call matrix_sln%add_value_pos(iposoffd, term2)
              end if
              call matrix_sln%add_value_pos(ipossymd, -term2)
              call matrix_sln%add_value_pos(ipossymoffd, -term)
            else
              rhs(iloc) = rhs(iloc) + term * hmaw
              rhs(isymnode) = rhs(isymnode) - term * hmaw
              call matrix_sln%add_value_pos(iposd, term)
              if (this%ibound(igwfnode) > 0) then
                call matrix_sln%add_value_pos(ipossymoffd, -term)
              end if
            end if
            !
            ! -- gwf is upstream
          else
            if (icflow /= 0) then
              rhsterm = term2 * hmaw + term * hgwf
              rhs(iloc) = rhs(iloc) + rhsterm
              rhs(isymnode) = rhs(isymnode) - rhsterm
              if (this%iboundpak(n) > 0) then
                call matrix_sln%add_value_pos(iposd, term2)
                call matrix_sln%add_value_pos(iposoffd, term)
              end if
              call matrix_sln%add_value_pos(ipossymd, -term)
              call matrix_sln%add_value_pos(ipossymoffd, -term2)
            else
              rhs(iloc) = rhs(iloc) + term * hgwf
              rhs(isymnode) = rhs(isymnode) - term * hgwf
              if (this%iboundpak(n) > 0) then
                call matrix_sln%add_value_pos(iposoffd, term)
              end if
              call matrix_sln%add_value_pos(ipossymd, -term)
            end if
          end if
        end if
        !
        ! -- increment maw connection counter
        idx = idx + 1
      end do
    end do
  end subroutine maw_fn

  !> @brief Calculate under-relaxation of groundwater flow model MAW Package heads
  !! for current outer iteration using the well bottom
  !<
  subroutine maw_nur(this, neqpak, x, xtemp, dx, inewtonur, dxmax, locmax)
    ! -- dummy
    class(MawType), intent(inout) :: this
    integer(I4B), intent(in) :: neqpak
    real(DP), dimension(neqpak), intent(inout) :: x
    real(DP), dimension(neqpak), intent(in) :: xtemp
    real(DP), dimension(neqpak), intent(inout) :: dx
    integer(I4B), intent(inout) :: inewtonur
    real(DP), intent(inout) :: dxmax
    integer(I4B), intent(inout) :: locmax
    ! -- local
    integer(I4B) :: n
    real(DP) :: botw
    real(DP) :: xx
    real(DP) :: dxx
    !
    ! -- Newton-Raphson under-relaxation
    do n = 1, this%nmawwells
      if (this%iboundpak(n) < 1) cycle
      botw = this%bot(n)
      !
      ! -- only apply Newton-Raphson under-relaxation if
      !    solution head is below the bottom of the well
      if (x(n) < botw) then
        inewtonur = 1
        xx = xtemp(n) * (DONE - DP9) + botw * DP9
        dxx = x(n) - xx
        if (abs(dxx) > abs(dxmax)) then
          locmax = n
          dxmax = dxx
        end if
        x(n) = xx
        dx(n) = DZERO
      end if
    end do
  end subroutine maw_nur

  !> @brief Calculate flows
  !<
  subroutine maw_cq(this, x, flowja, iadv)
    ! -- modules
    use TdisModule, only: delt
    use ConstantsModule, only: LENBOUNDNAME
    use BudgetModule, only: BudgetType
    ! -- dummy
    class(MawType), intent(inout) :: this
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    integer(I4B), optional, intent(in) :: iadv
    ! -- local
    real(DP) :: rrate
    ! -- for budget
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: ibnd
    real(DP) :: hmaw
    real(DP) :: cfw
    ! -- for observations
    ! -- formats
    !
    ! -- recalculate package HCOF and RHS terms with latest groundwater and
    !    maw heads prior to calling base budget functionality
    call this%maw_cfupdate()
    !
    ! -- call base functionality in bnd_cq.  This will calculate maw-gwf flows
    !    and put them into this%simvals
    call this%BndType%bnd_cq(x, flowja, iadv=1)
    !
    ! -- calculate maw budget flow and storage terms
    do n = 1, this%nmawwells
      this%qout(n) = DZERO
      this%qsto(n) = DZERO
      if (this%iflowingwells > 0) then
        this%qfw(n) = DZERO
      end if
      if (this%iboundpak(n) == 0) then
        cycle
      end if
      !
      ! -- set hmaw and xsto
      hmaw = this%xnewpak(n)
      this%xsto(n) = hmaw
      !
      ! -- add pumping rate to active maw well
      rrate = this%ratesim(n)
      !
      ! -- If flow is out of maw set qout to rrate.
      if (rrate < DZERO) then
        this%qout(n) = rrate
      end if
      !
      ! -- add flowing well
      if (this%iflowingwells > 0) then
        if (this%fwcond(n) > DZERO) then
          cfw = this%fwcondsim(n)
          this%xsto(n) = this%fwelev(n)
          rrate = cfw * (this%fwelev(n) - hmaw)
          this%qfw(n) = rrate
          !
          ! -- Subtract flowing well rrate from qout.
          this%qout(n) = this%qout(n) + rrate
        end if
      end if
      !
      ! -- Calculate qsto
      if (this%imawiss /= 1) then
        rrate = -this%area(n) * (this%xsto(n) - this%xoldsto(n)) / delt
        this%qsto(n) = rrate
      end if
    end do
    !
    ! -- gwf and constant flow
    ibnd = 1
    do n = 1, this%nmawwells
      hmaw = this%xnewpak(n)
      this%qconst(n) = DZERO
      do j = 1, this%ngwfnodes(n)
        rrate = -this%simvals(ibnd)
        this%qleak(ibnd) = rrate
        if (this%iboundpak(n) < 0) then
          this%qconst(n) = this%qconst(n) - rrate
          !
          ! -- If flow is out increment qout by -rrate.
          if (-rrate < DZERO) then
            this%qout(n) = this%qout(n) - rrate
          end if
        end if
        !
        ! -- increment ibnd counter
        ibnd = ibnd + 1
      end do
      !
      ! -- add additional flow terms to constant head term
      if (this%iboundpak(n) < 0) then
        !
        ! -- add well pumping rate
        this%qconst(n) = this%qconst(n) - this%ratesim(n)
        !
        ! -- add flowing well rate
        if (this%iflowingwells > 0) then
          this%qconst(n) = this%qconst(n) - this%qfw(n)
        end if
        !
        ! -- add storage term
        if (this%imawiss /= 1) then
          this%qconst(n) = this%qconst(n) - this%qsto(n)
        end if
      end if
    end do
    !
    ! -- fill the budget object
    call this%maw_fill_budobj()
  end subroutine maw_cq

  !> @brief Write flows to binary file and/or print flows to budget
  !<
  subroutine maw_ot_model_flows(this, icbcfl, ibudfl, icbcun, imap)
    ! -- dummy
    class(MawType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    integer(I4B), intent(in) :: icbcun
    integer(I4B), dimension(:), optional, intent(in) :: imap
    !
    ! -- write the flows from the budobj
    call this%BndType%bnd_ot_model_flows(icbcfl, ibudfl, icbcun, this%imap)
  end subroutine maw_ot_model_flows

  !> @brief Output MAW package flow terms.
  !<
  subroutine maw_ot_package_flows(this, icbcfl, ibudfl)
    use TdisModule, only: kstp, kper, delt, pertim, totim
    class(MawType) :: this
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
    ! -- Print maw flows table
    if (ibudfl /= 0 .and. this%iprflow /= 0) then
      call this%budobj%write_flowtable(this%dis, kstp, kper)
    end if
  end subroutine maw_ot_package_flows

  !> @brief Save maw-calculated values to binary file
  !<
  subroutine maw_ot_dv(this, idvsave, idvprint)
    use TdisModule, only: kstp, kper, pertim, totim
    use ConstantsModule, only: DHNOFLO, DHDRY
    use InputOutputModule, only: ulasav
    class(MawType) :: this
    integer(I4B), intent(in) :: idvsave
    integer(I4B), intent(in) :: idvprint
    integer(I4B) :: ibinun
    integer(I4B) :: n
    real(DP) :: v
    real(DP) :: d
    !
    ! -- set unit number for binary dependent variable output
    ibinun = 0
    if (this%iheadout /= 0) then
      ibinun = this%iheadout
    end if
    if (idvsave == 0) ibinun = 0
    !
    ! -- write maw binary output
    if (ibinun > 0) then
      do n = 1, this%nmawwells
        v = this%xnewpak(n)
        d = v - this%bot(n)
        if (this%iboundpak(n) == 0) then
          v = DHNOFLO
        else if (d <= DZERO) then
          v = DHDRY
        end if
        this%dbuff(n) = v
      end do
      call ulasav(this%dbuff, '            HEAD', &
                  kstp, kper, pertim, totim, &
                  this%nmawwells, 1, 1, ibinun)
    end if
    !
    ! -- write maw head table
    if (idvprint /= 0 .and. this%iprhed /= 0) then
      !
      ! -- set table kstp and kper
      call this%headtab%set_kstpkper(kstp, kper)
      !
      ! -- fill stage data
      do n = 1, this%nmawwells
        if (this%inamedbound == 1) then
          call this%headtab%add_term(this%cmawname(n))
        end if
        call this%headtab%add_term(n)
        call this%headtab%add_term(this%xnewpak(n))
      end do
    end if
  end subroutine maw_ot_dv

  !> @brief Write MAW budget to listing file
  !<
  subroutine maw_ot_bdsummary(this, kstp, kper, iout, ibudfl)
    ! -- module
    use TdisModule, only: totim, delt
    ! -- dummy
    class(MawType) :: this !< MawType object
    integer(I4B), intent(in) :: kstp !< time step number
    integer(I4B), intent(in) :: kper !< period number
    integer(I4B), intent(in) :: iout !< flag and unit number for the model listing file
    integer(I4B), intent(in) :: ibudfl !< flag indicating budget should be written
    !
    call this%budobj%write_budtable(kstp, kper, iout, ibudfl, totim, delt)
  end subroutine maw_ot_bdsummary

  !> @brief Deallocate memory
  !<
  subroutine maw_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(MawType) :: this
    ! -- local
    !
    ! -- budobj
    call this%budobj%budgetobject_da()
    deallocate (this%budobj)
    nullify (this%budobj)
    !
    ! -- head table
    if (this%iprhed > 0) then
      call this%headtab%table_da()
      deallocate (this%headtab)
      nullify (this%headtab)
    end if
    !
    ! -- character arrays
    call mem_deallocate(this%cmawbudget, 'CMAWBUDGET', this%memoryPath)
    call mem_deallocate(this%cmawname, 'CMAWNAME', this%memoryPath)
    call mem_deallocate(this%status, 'STATUS', this%memoryPath)
    !
    ! -- deallocate well data pointers in memory manager
    call mem_deallocate(this%ngwfnodes)
    call mem_deallocate(this%ieqn)
    call mem_deallocate(this%ishutoff)
    call mem_deallocate(this%ifwdischarge)
    call mem_deallocate(this%strt)
    call mem_deallocate(this%radius)
    call mem_deallocate(this%area)
    call mem_deallocate(this%pumpelev)
    call mem_deallocate(this%bot)
    call mem_deallocate(this%ratesim)
    call mem_deallocate(this%reduction_length)
    call mem_deallocate(this%fwelev)
    call mem_deallocate(this%fwcond)
    call mem_deallocate(this%fwrlen)
    call mem_deallocate(this%fwcondsim)
    call mem_deallocate(this%xsto)
    call mem_deallocate(this%xoldsto)
    call mem_deallocate(this%shutoffmin)
    call mem_deallocate(this%shutoffmax)
    call mem_deallocate(this%shutofflevel)
    call mem_deallocate(this%shutoffweight)
    call mem_deallocate(this%shutoffdq)
    call mem_deallocate(this%shutoffqold)
    !
    ! -- timeseries aware variables
    call mem_deallocate(this%mauxvar)
    call mem_deallocate(this%rate)
    call mem_deallocate(this%well_head)
    !
    ! -- connection data
    call mem_deallocate(this%iaconn)
    call mem_deallocate(this%gwfnodes)
    call mem_deallocate(this%sradius)
    call mem_deallocate(this%hk)
    call mem_deallocate(this%satcond)
    call mem_deallocate(this%simcond)
    call mem_deallocate(this%topscrn)
    call mem_deallocate(this%botscrn)
    !
    ! -- imap vector
    call mem_deallocate(this%imap)
    call mem_deallocate(this%dbuff)
    call mem_deallocate(this%cauxcbc, 'CAUXCBC', this%memoryPath)
    call mem_deallocate(this%qauxcbc)
    call mem_deallocate(this%qleak)
    call mem_deallocate(this%qfw)
    call mem_deallocate(this%qout)
    call mem_deallocate(this%qsto)
    call mem_deallocate(this%qconst)
    call mem_deallocate(this%denseterms)
    call mem_deallocate(this%viscratios)
    call mem_deallocate(this%idxlocnode)
    call mem_deallocate(this%idxdglo)
    call mem_deallocate(this%idxoffdglo)
    call mem_deallocate(this%idxsymdglo)
    call mem_deallocate(this%idxsymoffdglo)
    call mem_deallocate(this%xoldpak)
    !
    ! -- nullify pointers
    call mem_deallocate(this%xnewpak, 'HEAD', this%memoryPath)
    !
    ! -- scalars
    call mem_deallocate(this%correct_flow)
    call mem_deallocate(this%iprhed)
    call mem_deallocate(this%iheadout)
    call mem_deallocate(this%ibudgetout)
    call mem_deallocate(this%ibudcsv)
    call mem_deallocate(this%iflowingwells)
    call mem_deallocate(this%imawiss)
    call mem_deallocate(this%imawissopt)
    call mem_deallocate(this%nmawwells)
    call mem_deallocate(this%check_attr)
    call mem_deallocate(this%ishutoffcnt)
    call mem_deallocate(this%ieffradopt)
    call mem_deallocate(this%ioutredflowcsv)
    call mem_deallocate(this%satomega)
    call mem_deallocate(this%bditems)
    call mem_deallocate(this%theta)
    call mem_deallocate(this%kappa)
    call mem_deallocate(this%cbcauxitems)
    call mem_deallocate(this%idense)
    !
    ! -- pointers to gwf variables
    nullify (this%gwfiss)
    !
    ! -- call standard BndType deallocate
    call this%BndType%bnd_da()
  end subroutine maw_da

  !> @brief Define the list heading that is written to iout when PRINT_INPUT
  !! option is used.
  !<
  subroutine define_listlabel(this)
    class(MawType), intent(inout) :: this
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
  !! has access to these things.
  !<
  subroutine maw_set_pointers(this, neq, ibound, xnew, xold, flowja)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_checkin
    ! -- dummy
    class(MawType) :: this
    integer(I4B), pointer :: neq
    integer(I4B), dimension(:), pointer, contiguous :: ibound
    real(DP), dimension(:), pointer, contiguous :: xnew
    real(DP), dimension(:), pointer, contiguous :: xold
    real(DP), dimension(:), pointer, contiguous :: flowja
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: istart, iend
    !
    ! -- call base BndType set_pointers
    call this%BndType%set_pointers(neq, ibound, xnew, xold, flowja)
    !
    ! -- Set the MAW pointers
    !
    ! -- set package pointers
    istart = this%dis%nodes + this%ioffset + 1
    iend = istart + this%nmawwells - 1
    this%iboundpak => this%ibound(istart:iend)
    this%xnewpak => this%xnew(istart:iend)
    call mem_checkin(this%xnewpak, 'HEAD', this%memoryPath, 'X', &
                     this%memoryPathModel)
    call mem_allocate(this%xoldpak, this%nmawwells, 'XOLDPAK', this%memoryPath)
    !
    ! -- initialize xnewpak
    do n = 1, this%nmawwells
      this%xnewpak(n) = DEP20
    end do
  end subroutine maw_set_pointers

  ! -- Procedures related to observations (type-bound)

  !> @brief Return true because MAW package supports observations
  !!
  !! Overrides BndType%bnd_obs_supported()
  !<
  logical function maw_obs_supported(this)
    class(MawType) :: this
    !
    maw_obs_supported = .true.
  end function maw_obs_supported

  !> @brief Store observation type supported by MAW package
  !!
  !! Overrides BndType%bnd_df_obs
  !<
  subroutine maw_df_obs(this)
    ! -- dummy
    class(MawType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    ! -- Store obs type and assign procedure pointer
    !    for head observation type.
    call this%obs%StoreObsType('head', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => maw_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for frommvr observation type.
    call this%obs%StoreObsType('from-mvr', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => maw_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for conn-rate observation type.
    call this%obs%StoreObsType('maw', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => maw_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for rate observation type.
    call this%obs%StoreObsType('rate', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => maw_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for rate-to-mvr observation type.
    call this%obs%StoreObsType('rate-to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => maw_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for fw-rate observation type.
    call this%obs%StoreObsType('fw-rate', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => maw_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for rate-to-mvr observation type.
    call this%obs%StoreObsType('fw-to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => maw_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for storage observation type.
    call this%obs%StoreObsType('storage', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => maw_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for constant observation type.
    call this%obs%StoreObsType('constant', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => maw_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for cond observation type.
    call this%obs%StoreObsType('conductance', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => maw_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for fw-conductance observation type.
    call this%obs%StoreObsType('fw-conductance', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => maw_process_obsID
  end subroutine maw_df_obs

  !> @brief Calculate observations this time step and call
  !! ObsType%SaveOneSimval for each MawType observation.
  !<
  subroutine maw_bd_obs(this)
    ! -- dummy
    class(MawType) :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: jj
    integer(I4B) :: n
    integer(I4B) :: nn
    integer(I4B) :: jpos
    real(DP) :: cmaw
    real(DP) :: hmaw
    real(DP) :: v
    real(DP) :: qfact
    type(ObserveType), pointer :: obsrv => null()
    !
    ! Calculate, save, and write simulated values for all MAW observations
    if (this%obs%npakobs > 0) then
      call this%obs%obs_bd_clear()
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv
        do j = 1, obsrv%indxbnds_count
          v = DNODATA
          jj = obsrv%indxbnds(j)
          select case (obsrv%ObsTypeId)
          case ('HEAD')
            if (this%iboundpak(jj) /= 0) then
              v = this%xnewpak(jj)
            end if
          case ('FROM-MVR')
            if (this%iboundpak(jj) /= 0) then
              if (this%imover == 1) then
                v = this%pakmvrobj%get_qfrommvr(jj)
              end if
            end if
          case ('MAW')
            n = this%imap(jj)
            if (this%iboundpak(n) /= 0) then
              v = this%qleak(jj)
            end if
          case ('RATE')
            if (this%iboundpak(jj) /= 0) then
              v = this%ratesim(jj)
              if (v < DZERO .and. this%qout(jj) < DZERO) then
                qfact = v / this%qout(jj)
                if (this%imover == 1) then
                  v = v + this%pakmvrobj%get_qtomvr(jj) * qfact
                end if
              end if
            end if
          case ('RATE-TO-MVR')
            if (this%iboundpak(jj) /= 0) then
              if (this%imover == 1) then
                v = this%ratesim(jj)
                qfact = DZERO
                if (v < DZERO .and. this%qout(jj) < DZERO) then
                  qfact = v / this%qout(jj)
                end if
                v = this%pakmvrobj%get_qtomvr(jj) * qfact
                if (v > DZERO) then
                  v = -v
                end if
              end if
            end if
          case ('FW-RATE')
            if (this%iboundpak(jj) /= 0 .and. this%iflowingwells > 0) then
              hmaw = this%xnewpak(jj)
              cmaw = this%fwcondsim(jj)
              v = cmaw * (this%fwelev(jj) - hmaw)
              if (v < DZERO .and. this%qout(jj) < DZERO) then
                qfact = v / this%qout(jj)
                if (this%imover == 1) then
                  v = v + this%pakmvrobj%get_qtomvr(jj) * qfact
                end if
              end if
            end if
          case ('FW-TO-MVR')
            if (this%iboundpak(jj) /= 0 .and. this%iflowingwells > 0) then
              if (this%imover == 1) then
                hmaw = this%xnewpak(jj)
                cmaw = this%fwcondsim(jj)
                v = cmaw * (this%fwelev(jj) - hmaw)
                qfact = DZERO
                if (v < DZERO .and. this%qout(jj) < DZERO) then
                  qfact = v / this%qout(jj)
                end if
                v = this%pakmvrobj%get_qtomvr(jj) * qfact
                if (v > DZERO) then
                  v = -v
                end if
              end if
            end if
          case ('STORAGE')
            if (this%iboundpak(jj) /= 0 .and. this%imawissopt /= 1) then
              v = this%qsto(jj)
            end if
          case ('CONSTANT')
            if (this%iboundpak(jj) /= 0) then
              v = this%qconst(jj)
            end if
          case ('CONDUCTANCE')
            n = this%imap(jj)
            if (this%iboundpak(n) /= 0) then
              nn = jj - this%iaconn(n) + 1
              jpos = this%get_jpos(n, nn)
              v = this%simcond(jpos)
            end if
          case ('FW-CONDUCTANCE')
            if (this%iboundpak(jj) /= 0) then
              v = this%fwcondsim(jj)
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
    !
    ! -- Write the MAW reduced flows to csv file entries for this step
    if (this%ioutredflowcsv > 0) then
      call this%maw_redflow_csv_write()
    end if
  end subroutine maw_bd_obs

  !> @brief Process each observation
  !!
  !! Only done the first stress period since boundaries are fixed for the
  !! simulation
  !<
  subroutine maw_rp_obs(this)
    use TdisModule, only: kper
    ! -- dummy
    class(MawType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: nn1
    integer(I4B) :: nn2
    integer(I4B) :: jj
    character(len=LENBOUNDNAME) :: bname
    logical :: jfound
    class(ObserveType), pointer :: obsrv => null()
    ! -- formats
10  format('Boundary "', a, '" for observation "', a, &
           '" is invalid in package "', a, '"')
    !
    if (kper == 1) then
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv
        !
        ! -- get node number 1
        nn1 = obsrv%NodeNumber
        if (nn1 == NAMEDBOUNDFLAG) then
          bname = obsrv%FeatureName
          if (bname /= '') then
            ! -- Observation maw is based on a boundary name.
            !    Iterate through all multi-aquifer wells to identify and store
            !    corresponding index in bound array.
            jfound = .false.
            if (obsrv%ObsTypeId == 'MAW' .or. &
                obsrv%ObsTypeId == 'CONDUCTANCE') then
              do j = 1, this%nmawwells
                do jj = this%iaconn(j), this%iaconn(j + 1) - 1
                  if (this%boundname(jj) == bname) then
                    jfound = .true.
                    call obsrv%AddObsIndex(jj)
                  end if
                end do
              end do
            else
              do j = 1, this%nmawwells
                if (this%cmawname(j) == bname) then
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
            if (obsrv%ObsTypeId == 'MAW' .or. &
                obsrv%ObsTypeId == 'CONDUCTANCE') then
              nn2 = obsrv%NodeNumber2
              j = this%iaconn(nn1) + nn2 - 1
              call obsrv%AddObsIndex(j)
            else
              call obsrv%AddObsIndex(nn1)
            end if
          else
            errmsg = 'Programming error in maw_rp_obs'
            call store_error(errmsg)
          end if
        end if
        !
        ! -- catch non-cumulative observation assigned to observation defined
        !    by a boundname that is assigned to more than one element
        if (obsrv%ObsTypeId == 'HEAD') then
          if (obsrv%indxbnds_count > 1) then
            write (errmsg, '(a,3(1x,a))') &
              trim(adjustl(obsrv%ObsTypeId)), &
              'for observation', trim(adjustl(obsrv%Name)), &
              'must be assigned to a multi-aquifer well with a unique boundname.'
            call store_error(errmsg)
          end if
        end if
        !
        ! -- check that index values are valid
        if (obsrv%ObsTypeId == 'MAW' .or. &
            obsrv%ObsTypeId == 'CONDUCTANCE') then
          do j = 1, obsrv%indxbnds_count
            nn1 = obsrv%indxbnds(j)
            n = this%imap(nn1)
            nn2 = nn1 - this%iaconn(n) + 1
            jj = this%iaconn(n + 1) - this%iaconn(n)
            if (nn1 < 1 .or. nn1 > this%maxbound) then
              write (errmsg, '(3(a,1x),i0,1x,a,i0,a)') &
                trim(adjustl(obsrv%ObsTypeId)), &
                'multi-aquifer well connection number must be greater than 0', &
                'and less than', jj, '(specified value is ', nn2, ').'
              call store_error(errmsg)
            end if
          end do
        else
          do j = 1, obsrv%indxbnds_count
            nn1 = obsrv%indxbnds(j)
            if (nn1 < 1 .or. nn1 > this%nmawwells) then
              write (errmsg, '(3(a,1x),i0,1x,a,i0,a)') &
                trim(adjustl(obsrv%ObsTypeId)), &
                'multi-aquifer well must be greater than 0 ', &
                'and less than or equal to', this%nmawwells, &
                '(specified value is ', nn1, ').'
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
  end subroutine maw_rp_obs

  !
  ! -- Procedures related to observations (NOT type-bound)

  !> @brief This procedure is pointed to by ObsDataType%ProcesssIdPtr. It
  !! processes the ID string of an observation definition for MAW package
  !! observations.
  !<
  subroutine maw_process_obsID(obsrv, dis, inunitobs, iout)
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
    ! formats
    !
    string = obsrv%IDstring
    ! -- Extract multi-aquifer well number from string and store it.
    !    If 1st item is not an integer(I4B), it should be a
    !    maw name--deal with it.
    icol = 1
    ! -- get multi-aquifer well number or boundary name
    call extract_idnum_or_bndname(string, icol, istart, istop, nn1, bndname)
    if (nn1 == NAMEDBOUNDFLAG) then
      obsrv%FeatureName = bndname
    else
      if (obsrv%ObsTypeId == 'MAW' .or. &
          obsrv%ObsTypeId == 'CONDUCTANCE') then
        call extract_idnum_or_bndname(string, icol, istart, istop, nn2, bndname)
        if (len_trim(bndName) < 1 .and. nn2 < 0) then
          write (errmsg, '(a,1x,a,a,1x,a,1x,a)') &
            'For observation type', trim(adjustl(obsrv%ObsTypeId)), &
            ', ID given as an integer and not as boundname,', &
            'but ID2 (icon) is missing.  Either change ID to valid', &
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
    ! -- store multi-aquifer well number (NodeNumber)
    obsrv%NodeNumber = nn1
  end subroutine maw_process_obsID

  !
  ! -- private MAW methods

  !> @brief Initialize the auto flow reduce csv output file
  !<
  subroutine maw_redflow_csv_init(this, fname)
    ! -- dummy variables
    class(MawType), intent(inout) :: this !< MawType object
    character(len=*), intent(in) :: fname
    ! -- format
    character(len=*), parameter :: fmtredflowcsv = &
      "(4x, 'MAW REDUCED FLOW INFORMATION WILL BE SAVED TO FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"

    this%ioutredflowcsv = getunit()
    call openfile(this%ioutredflowcsv, this%iout, fname, 'CSV', &
                  filstat_opt='REPLACE')
    write (this%iout, fmtredflowcsv) trim(adjustl(fname)), &
      this%ioutredflowcsv
    write (this%ioutredflowcsv, '(a)') &
      'time,period,step,MAWnumber,rate-requested,rate-actual,maw-reduction'
  end subroutine maw_redflow_csv_init

  !> @brief MAW reduced flows only when & where they occur
  !<
  subroutine maw_redflow_csv_write(this)
    ! -- modules
    use TdisModule, only: totim, kstp, kper
    ! -- dummy variables
    class(MawType), intent(inout) :: this !< MawType object
    ! -- local
    integer(I4B) :: n
    !integer(I4B) :: nodereduced
    !integer(I4B) :: nodeuser
    real(DP) :: v
    ! -- format
    do n = 1, this%nmawwells
      !
      ! -- test if node is constant or inactive
      if (this%status(n) .ne. 'ACTIVE') then
        cycle
      end if
      v = this%rate(n) - this%ratesim(n) !reductions in extraction will be negative and reductions in injection will be positive; follows convention of WEL AUTO_FLOW_REDUCE_CSV
      if (abs(v) > DEM9) then !need to check absolute value of difference for both extraction and injection; using 1e-9 as epsilon value but could be tweaked
        write (this%ioutredflowcsv, '(*(G0,:,","))') &
          totim, kper, kstp, n, this%rate(n), this%ratesim(n), v
      end if
    end do
  end subroutine maw_redflow_csv_write

  !> @brief Calculate the appropriate saturated conductance to use based on
  !! aquifer and multi-aquifer well characteristics
  !<
  subroutine maw_calculate_satcond(this, i, j, node)
    ! -- dummy
    class(MawType), intent(inout) :: this
    integer(I4B), intent(in) :: i
    integer(I4B), intent(in) :: j
    integer(I4B), intent(in) :: node
    ! -- local
    integer(I4B) :: iTcontrastErr
    integer(I4B) :: jpos
    real(DP) :: c
    real(DP) :: k11
    real(DP) :: k22
    real(DP) :: sqrtk11k22
    real(DP) :: hks
    real(DP) :: area
    real(DP) :: eradius
    real(DP) :: topw
    real(DP) :: botw
    real(DP) :: tthkw
    real(DP) :: tthka
    real(DP) :: Tcontrast
    real(DP) :: skin
    real(DP) :: ravg
    real(DP) :: slen
    real(DP) :: pavg
    real(DP) :: gwfsat
    real(DP) :: gwftop
    real(DP) :: gwfbot
    real(DP) :: lc1
    real(DP) :: lc2
    real(DP) :: dx
    real(DP) :: dy
    real(DP) :: Txx
    real(DP) :: Tyy
    real(DP) :: T2pi
    real(DP) :: yx4
    real(DP) :: xy4
    ! -- formats
    !
    ! -- initialize conductance variables
    iTcontrastErr = 0
    lc1 = DZERO
    lc2 = DZERO
    !
    ! -- calculate connection position
    jpos = this%get_jpos(i, j)
    !
    ! -- set K11 and K22
    k11 = this%gwfk11(node)
    if (this%gwfik22 == 0) then
      k22 = this%gwfk11(node)
    else
      k22 = this%gwfk22(node)
    end if
    sqrtk11k22 = sqrt(k11 * k22)
    !
    ! -- set gwftop, gwfbot, and gwfsat
    gwftop = this%dis%top(node)
    gwfbot = this%dis%bot(node)
    tthka = gwftop - gwfbot
    gwfsat = this%gwfsat(node)
    !
    ! -- set top and bottom of well screen
    c = DZERO
    topw = this%topscrn(jpos)
    botw = this%botscrn(jpos)
    tthkw = topw - botw
    !
    ! -- scale screen thickness using gwfsat (for NPF Package THICKSTRT)
    if (gwftop == topw .and. gwfbot == botw) then
      if (this%icelltype(node) == 0) then
        tthkw = tthkw * gwfsat
        tthka = tthka * gwfsat
      end if
    end if
    !
    ! -- calculate the aquifer transmissivity (T2pi)
    T2pi = DTWOPI * tthka * sqrtk11k22
    !
    ! -- calculate effective radius
    if (this%dis%ndim == 3 .and. this%ieffradopt /= 0) then
      Txx = k11 * tthka
      Tyy = k22 * tthka
      dx = sqrt(this%dis%area(node))
      dy = dx
      yx4 = (Tyy / Txx)**DQUARTER
      xy4 = (Txx / Tyy)**DQUARTER
      eradius = 0.28_DP * ((yx4 * dx)**DTWO + &
                           (xy4 * dy)**DTWO)**DHALF / (yx4 + xy4)
    else
      area = this%dis%area(node)
      eradius = sqrt(area / (DEIGHT * DPI))
    end if
    !
    ! -- conductance calculations
    ! -- Thiem equation (1) and cumulative Thiem and skin equations (3)
    if (this%ieqn(i) == 1 .or. this%ieqn(i) == 3) then
      lc1 = log(eradius / this%radius(i)) / T2pi
    end if
    !
    ! -- skin equation (2) and cumulative Thiem and skin equations (3)
    if (this%ieqn(i) == 2 .or. this%ieqn(i) == 3) then
      hks = this%hk(jpos)
      if (tthkw * hks > DZERO) then
        Tcontrast = (sqrtk11k22 * tthka) / (hks * tthkw)
        skin = (Tcontrast - DONE) * log(this%sradius(jpos) / this%radius(i))
        !
        ! -- trap invalid transmissvity contrast if using skin equation (2).
        !    Not trapped for cumulative Thiem and skin equations (3)
        !    because the MNW2 package allowed this condition (for
        !    backward compatibility with the MNW2 package for
        !    MODFLOW-2005, MODFLOW-NWT, and MODFLOW-USG).
        if (Tcontrast <= 1 .and. this%ieqn(i) == 2) then
          iTcontrastErr = 1
          write (errmsg, '(a,g0,a,1x,i0,1x,a,1x,i0,a,4(1x,a))') &
            'Invalid calculated transmissivity contrast (', Tcontrast, &
            ') for maw well', i, 'connection', j, '.', 'This happens when the', &
            'skin transmissivity equals or exceeds the aquifer transmissivity.', &
            'Consider decreasing HK_SKIN for the connection or using the', &
            'CUMULATIVE or MEAN conductance equations.'
          call store_error(errmsg)
        else
          lc2 = skin / T2pi
        end if
      end if
    end if
    ! -- conductance using screen elevations, hk, well radius,
    !    and screen radius
    if (this%ieqn(i) == 4) then
      hks = this%hk(jpos)
      ravg = DHALF * (this%radius(i) + this%sradius(jpos))
      slen = this%sradius(jpos) - this%radius(i)
      pavg = DTWOPI * ravg
      c = hks * pavg * tthkw / slen
    end if
    !
    ! -- calculate final conductance for Thiem (1), Skin (2), and
    ! and cumulative Thiem and skin equations (3)
    if (this%ieqn(i) < 4) then
      if (lc1 + lc2 /= DZERO) then
        c = DONE / (lc1 + lc2)
      else
        c = -DNODATA
      end if
    end if
    !
    ! -- ensure that the conductance is not negative. Only write error message
    !    if error condition has not occurred for skin calculations (LC2)
    if (c < DZERO .and. iTcontrastErr == 0) then
      write (errmsg, '(a,g0,a,1x,i0,1x,a,1x,i0,a,4(1x,a))') &
        'Invalid calculated negative conductance (', c, &
        ') for maw well', i, 'connection', j, '.', 'this happens when the', &
        'skin transmissivity equals or exceeds the aquifer transmissivity.', &
        'consider decreasing hk_skin for the connection or using the', &
        'mean conductance equation.'
      call store_error(errmsg)
    end if
    !
    ! -- set saturated conductance
    this%satcond(jpos) = c
  end subroutine maw_calculate_satcond

  !> @brief Calculate the saturation between the aquifer maw well_head
  !<
  subroutine maw_calculate_saturation(this, n, j, node, sat)
    ! -- dummy
    class(MawType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: j
    integer(I4B), intent(in) :: node
    real(DP), intent(inout) :: sat
    ! -- local
    integer(I4B) :: jpos
    real(DP) :: h_temp
    real(DP) :: hwell
    real(DP) :: topw
    real(DP) :: botw
    ! -- formats
    !
    ! -- initialize saturation
    sat = DZERO
    !
    ! -- calculate current saturation for convertible cells
    if (this%icelltype(node) /= 0) then
      !
      ! -- set hwell
      hwell = this%xnewpak(n)
      !
      ! -- set connection position
      jpos = this%get_jpos(n, j)
      !
      ! -- set top and bottom of the well connection
      topw = this%topscrn(jpos)
      botw = this%botscrn(jpos)
      !
      ! -- calculate appropriate saturation
      if (this%inewton /= 1) then
        h_temp = this%xnew(node)
        if (h_temp < botw) then
          h_temp = botw
        end if
        if (hwell < botw) then
          hwell = botw
        end if
        h_temp = DHALF * (h_temp + hwell)
      else
        h_temp = this%xnew(node)
        if (hwell > h_temp) then
          h_temp = hwell
        end if
        if (h_temp < botw) then
          h_temp = botw
        end if
      end if
      ! -- calculate saturation
      sat = sQuadraticSaturation(topw, botw, h_temp, this%satomega)
    else
      sat = DONE
    end if
  end subroutine maw_calculate_saturation

  !> @brief Calculate matrix terms for a multi-aquifer well connection. Terms
  !! for fc and fn methods are calculated based on whether term2 is passed
  !! Arguments are as follows:
  !!     n       : maw well number
  !!     j       : connection number for well n
  !!     icflow  : flag indicating that flow should be corrected
  !!     cmaw    : maw-gwf conducance
  !!     cterm   : correction term for flow to dry cell
  !!     term    : xxx
  !!     flow    : calculated flow for this connection, positive into well
  !!     term2   : xxx
  !<
  subroutine maw_calculate_conn_terms(this, n, j, icflow, cmaw, cterm, term, &
                                      flow, term2)
    ! -- dummy
    class(MawType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: j
    integer(I4B), intent(inout) :: icflow
    real(DP), intent(inout) :: cmaw
    real(DP), intent(inout) :: cterm
    real(DP), intent(inout) :: term
    real(DP), intent(inout) :: flow
    real(DP), intent(inout), optional :: term2
    ! -- local
    logical(LGP) :: correct_flow
    integer(I4B) :: inewton
    integer(I4B) :: jpos
    integer(I4B) :: igwfnode
    real(DP) :: hmaw
    real(DP) :: hgwf
    real(DP) :: hups
    real(DP) :: hdowns
    real(DP) :: sat
    real(DP) :: tmaw
    real(DP) :: bmaw
    real(DP) :: en
    real(DP) :: hbar
    real(DP) :: drterm
    real(DP) :: dhbarterm
    real(DP) :: vscratio
    !
    ! -- initialize terms
    cterm = DZERO
    vscratio = DONE
    icflow = 0
    if (present(term2)) then
      inewton = 1
    else
      inewton = 0
    end if
    !
    ! -- set common terms
    jpos = this%get_jpos(n, j)
    igwfnode = this%get_gwfnode(n, j)
    hgwf = this%xnew(igwfnode)
    hmaw = this%xnewpak(n)
    tmaw = this%topscrn(jpos)
    bmaw = this%botscrn(jpos)
    !
    ! -- if vsc active, select appropriate viscosity ratio
    if (this%ivsc == 1) then
      ! flow out of well (flow is negative)
      if (flow < 0) then
        vscratio = this%viscratios(1, n)
      else
        vscratio = this%viscratios(2, n)
      end if
    end if
    !
    ! -- calculate saturation
    call this%maw_calculate_saturation(n, j, igwfnode, sat)
    cmaw = this%satcond(jpos) * vscratio * sat
    !
    ! -- set upstream head, term, and term2 if returning newton terms
    if (inewton == 1) then
      term = DZERO
      term2 = DZERO
      hups = hmaw
      if (hgwf > hups) then
        hups = hgwf
      end if
      !
      ! -- calculate the derivative of saturation
      drterm = sQuadraticSaturationDerivative(tmaw, bmaw, hups, this%satomega)
    else
      term = cmaw
    end if
    !
    ! -- calculate correction term if flow_correction option specified
    if (this%correct_flow) then
      !
      ! -- set bmaw, determine en, and set correct_flow flag
      en = max(bmaw, this%dis%bot(igwfnode))
      correct_flow = .FALSE.
      if (hmaw < en) then
        correct_flow = .TRUE.
      end if
      if (hgwf < en .and. this%icelltype(igwfnode) /= 0) then
        correct_flow = .TRUE.
      end if
      !
      ! -- if flow should be corrected because hgwf or hmaw is below bottom
      !    then calculate correction term (cterm)
      if (correct_flow) then
        icflow = 1
        hdowns = min(hmaw, hgwf)
        hbar = sQuadratic0sp(hdowns, en, this%satomega)
        if (hgwf > hmaw) then
          cterm = cmaw * (hmaw - hbar)
        else
          cterm = cmaw * (hbar - hgwf)
        end if
      end if
      !
      ! -- if newton formulation then calculate newton terms
      if (inewton /= 0) then
        !
        ! -- maw is upstream
        if (hmaw > hgwf) then
          hbar = sQuadratic0sp(hgwf, en, this%satomega)
          term = drterm * this%satcond(jpos) * vscratio * (hbar - hmaw)
          dhbarterm = sQuadratic0spDerivative(hgwf, en, this%satomega)
          term2 = cmaw * (dhbarterm - DONE)
          !
          ! -- gwf is upstream
        else
          hbar = sQuadratic0sp(hmaw, en, this%satomega)
          term = -drterm * this%satcond(jpos) * vscratio * (hgwf - hbar)
          dhbarterm = sQuadratic0spDerivative(hmaw, en, this%satomega)
          term2 = cmaw * (DONE - dhbarterm)
        end if
      end if
    else
      !
      ! -- flow is not corrected, so calculate term for newton formulation
      if (inewton /= 0) then
        term = drterm * this%satcond(jpos) * vscratio * (hgwf - hmaw)
      end if
    end if
    !
    ! -- calculate flow relative to maw for fc and bd
    flow = DZERO
    if (inewton == 0) then
      flow = term * (hgwf - hmaw) + cterm
    end if
    !
    ! -- add density part here
    if (this%idense /= 0 .and. inewton == 0) then
      call this%maw_calculate_density_exchange(jpos, hmaw, hgwf, cmaw, &
                                               bmaw, flow, term, cterm)
    end if
  end subroutine maw_calculate_conn_terms

  !> @brief Calculate well pumping rate based on constraints
  !<
  subroutine maw_calculate_wellq(this, n, hmaw, q)
    ! -- dummy
    class(MawType) :: this
    integer(I4B), intent(in) :: n
    real(DP), intent(in) :: hmaw
    real(DP), intent(inout) :: q
    ! -- local
    real(DP) :: scale
    real(DP) :: tp
    real(DP) :: bt
    real(DP) :: rate
    real(DP) :: weight
    real(DP) :: dq
    !
    ! -- Initialize q
    q = DZERO
    !
    ! -- Assign rate as the user-provided base pumping rate
    rate = this%rate(n)
    !
    ! -- Assign q differently depending on whether this is an extraction well
    !    (rate < 0) or an injection well (rate > 0).
    if (rate < DZERO) then
      !
      ! -- If well shut off is activated, then turn off well if necessary,
      !    or if shut off is not activated then check to see if rate scaling
      !    is on.
      if (this%shutofflevel(n) /= DEP20) then
        call this%maw_calculate_qpot(n, q)
        if (q < DZERO) q = DZERO
        if (q > -rate) q = -rate

        if (this%ishutoffcnt == 1) then
          this%shutoffweight(n) = DONE
          this%shutoffdq(n) = DZERO
          this%shutoffqold(n) = q
        end if

        dq = q - this%shutoffqold(n)
        weight = this%shutoffweight(n)
        !
        ! -- for flip-flop condition, decrease factor
        if (this%shutoffdq(n) * dq < DZERO) then
          weight = this%theta * this%shutoffweight(n)
          !
          ! -- when change is of same sign, increase factor
        else
          weight = this%shutoffweight(n) + this%kappa
        end if
        if (weight > DONE) weight = DONE

        q = this%shutoffqold(n) + weight * dq

        this%shutoffqold(n) = q
        this%shutoffdq(n) = dq
        this%shutoffweight(n) = weight
        !
        ! -- If shutoffmin and shutoffmax are specified then apply
        !    additional checks for when to shut off the well.
        if (this%shutoffmin(n) > DZERO) then
          if (hmaw < this%shutofflevel(n)) then
            !
            ! -- calculate adjusted well rate subject to constraints
            ! -- well is shutoff
            if (this%ishutoff(n) /= 0) then
              q = DZERO
              !
              ! --- well is not shut off
            else
              ! -- turn off well if q is less than the minimum rate and
              !    reset the ishutoff flag if at least on iteration 3
              if (q < this%shutoffmin(n)) then
                if (this%ishutoffcnt > 2) then
                  this%ishutoff(n) = 1
                end if
                q = DZERO
                !
                ! -- leave well on and use the specified rate
                !    or the potential rate
              end if
            end if
            !
            ! -- try to use the specified rate or the potential rate
          else
            if (q > this%shutoffmax(n)) then
              if (this%ishutoffcnt <= 2) then
                this%ishutoff(n) = 0
              end if
            end if
            if (this%ishutoff(n) /= 0) then
              q = DZERO
            end if
          end if
        end if

        if (q /= DZERO) q = -q

      else
        scale = DONE
        !
        ! -- Apply rate scaling by reducing pumpage when hmaw is less than the
        !    sum of maw pump elevation (pumpelev) and the specified reduction
        !    length.  The rate will go to zero as hmaw drops to the pump
        !    elevation.
        if (this%reduction_length(n) /= DEP20) then
          bt = this%pumpelev(n)
          tp = bt + this%reduction_length(n)
          scale = sQSaturation(tp, bt, hmaw)
        end if
        q = scale * rate
      end if
      !
    else
      !
      ! -- Handle the injection case (rate > 0) differently than extraction.
      q = rate
      if (this%shutofflevel(n) /= DEP20) then
        call this%maw_calculate_qpot(n, q)
        q = -q
        if (q < DZERO) q = DZERO
        if (q > rate) q = rate

        if (this%ishutoffcnt == 1) then
          this%shutoffweight(n) = DONE
          this%shutoffdq(n) = DZERO
          this%shutoffqold(n) = q
        end if

        dq = q - this%shutoffqold(n)
        weight = this%shutoffweight(n)
        !
        ! -- for flip-flop condition, decrease factor
        if (this%shutoffdq(n) * dq < DZERO) then
          weight = this%theta * this%shutoffweight(n)
          !
          ! -- when change is of same sign, increase factor
        else
          weight = this%shutoffweight(n) + this%kappa
        end if
        if (weight > DONE) weight = DONE

        q = this%shutoffqold(n) + weight * dq

        this%shutoffqold(n) = q
        this%shutoffdq(n) = dq
        this%shutoffweight(n) = weight

      else
        scale = DONE
        !
        ! -- Apply rate scaling for an injection well by reducing the
        !    injection rate as hmaw rises above the pump elevation.  The rate
        !    will approach zero as hmaw approaches pumpelev + reduction_length.
        if (this%reduction_length(n) /= DEP20) then
          bt = this%pumpelev(n)
          tp = bt + this%reduction_length(n)
          scale = DONE - sQSaturation(tp, bt, hmaw)
        end if
        q = scale * rate
      end if
    end if
  end subroutine maw_calculate_wellq

  !> @brief Calculate groundwater inflow to a maw well
  !<
  subroutine maw_calculate_qpot(this, n, qnet)
    use TdisModule, only: delt
    ! -- dummy
    class(MawType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    real(DP), intent(inout) :: qnet
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: jpos
    integer(I4B) :: igwfnode
    real(DP) :: bt
    real(DP) :: tp
    real(DP) :: scale
    real(DP) :: cfw
    real(DP) :: hdterm
    real(DP) :: sat
    real(DP) :: cmaw
    real(DP) :: hgwf
    real(DP) :: bmaw
    real(DP) :: h_temp
    real(DP) :: hv
    real(DP) :: vscratio
    ! -- format
    !
    ! -- initialize qnet and h_temp
    qnet = DZERO
    vscratio = DONE
    h_temp = this%shutofflevel(n)
    !
    ! -- if vsc active, select appropriate viscosity ratio
    if (this%ivsc == 1) then
      ! flow out of well (flow is negative)
      if (qnet < 0) then
        vscratio = this%viscratios(1, n)
      else
        vscratio = this%viscratios(2, n)
      end if
    end if
    !
    ! -- calculate discharge to flowing wells
    if (this%iflowingwells > 0) then
      if (this%fwcond(n) > DZERO) then
        bt = this%fwelev(n)
        tp = bt + this%fwrlen(n)
        scale = sQSaturation(tp, bt, h_temp)
        cfw = scale * this%fwcond(n) * this%viscratios(2, n)
        this%ifwdischarge(n) = 0
        if (cfw > DZERO) then
          this%ifwdischarge(n) = 1
          this%xsto(n) = bt
        end if
        qnet = qnet + cfw * (bt - h_temp)
      end if
    end if
    !
    ! -- calculate maw storage changes
    if (this%imawiss /= 1) then
      if (this%ifwdischarge(n) /= 1) then
        hdterm = this%xoldsto(n) - h_temp
      else
        hdterm = this%xoldsto(n) - this%fwelev(n)
      end if
      qnet = qnet - (this%area(n) * hdterm / delt)
    end if
    !
    ! -- calculate inflow from aquifer
    do j = 1, this%ngwfnodes(n)
      jpos = this%get_jpos(n, j)
      igwfnode = this%get_gwfnode(n, j)
      call this%maw_calculate_saturation(n, j, igwfnode, sat)
      cmaw = this%satcond(jpos) * vscratio * sat
      hgwf = this%xnew(igwfnode)
      bmaw = this%botscrn(jpos)
      hv = h_temp
      if (hv < bmaw) then
        hv = bmaw
      end if
      if (hgwf < bmaw) then
        hgwf = bmaw
      end if
      qnet = qnet + cmaw * (hgwf - hv)
    end do
  end subroutine maw_calculate_qpot

  !> @brief Update MAW satcond and package rhs and hcof
  !<
  subroutine maw_cfupdate(this)
    class(MawType) :: this
    ! -- dummy
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: jpos
    integer(I4B) :: icflow
    integer(I4B) :: ibnd
    real(DP) :: flow
    real(DP) :: cmaw
    real(DP) :: hmaw
    real(DP) :: cterm
    real(DP) :: term
    !
    ! -- Return if no maw wells
    if (this%nbound .eq. 0) return
    !
    ! -- Update shutoff count
    this%ishutoffcnt = this%ishutoffcnt + 1
    !
    ! -- Calculate hcof and rhs for each maw entry
    ibnd = 1
    do n = 1, this%nmawwells
      hmaw = this%xnewpak(n)
      do j = 1, this%ngwfnodes(n)
        jpos = this%get_jpos(n, j)
        this%hcof(ibnd) = DZERO
        this%rhs(ibnd) = DZERO
        !
        ! -- set bound, hcof, and rhs components
        !
        ! -- use connection method so the gwf-maw budget flows
        !    are consistent with the maw-gwf budget flows
        if (this%iboundpak(n) == 0) then
          cmaw = DZERO
          term = DZERO
          cterm = DZERO
        else
          call this%maw_calculate_conn_terms(n, j, icflow, cmaw, cterm, &
                                             term, flow)
        end if
        this%simcond(jpos) = cmaw
        this%bound(2, ibnd) = cmaw
        this%hcof(ibnd) = -term
        this%rhs(ibnd) = -term * hmaw + cterm
        !
        ! -- increment boundary number
        ibnd = ibnd + 1
      end do
    end do
  end subroutine maw_cfupdate

  !> @brief Set up the budget object that stores all the maw flows
  !! The terms listed here must correspond in number and order to the ones
  !! listed in the maw_fill_budobj routine.
  !<
  subroutine maw_setup_budobj(this)
    ! -- modules
    use ConstantsModule, only: LENBUDTXT
    ! -- dummy
    class(MawType) :: this
    ! -- local
    integer(I4B) :: nbudterm
    integer(I4B) :: n, j, n2
    real(DP) :: q
    integer(I4B) :: maxlist, naux
    integer(I4B) :: idx
    character(len=LENBUDTXT) :: text
    character(len=LENBUDTXT), dimension(1) :: auxtxt
    !
    ! -- Determine the number of maw budget terms. These are fixed for
    !    the simulation and cannot change.
    ! gwf rate [flowing_well] storage constant_flow [frommvr tomvr tomvrcf [tomvrfw]] [aux]
    nbudterm = 4
    if (this%iflowingwells > 0) then
      nbudterm = nbudterm + 1
    end if
    if (this%imover == 1) then
      nbudterm = nbudterm + 3
      if (this%iflowingwells > 0) then
        nbudterm = nbudterm + 1
      end if
    end if
    if (this%naux > 0) nbudterm = nbudterm + 1
    !
    ! -- set up budobj
    call budgetobject_cr(this%budobj, this%packName)
    call this%budobj%budgetobject_df(this%nmawwells, nbudterm, 0, 0, &
                                     ibudcsv=this%ibudcsv)
    idx = 0
    !
    ! -- Go through and set up each budget term
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
    do n = 1, this%nmawwells
      do j = 1, this%ngwfnodes(n)
        n2 = this%get_gwfnode(n, j)
        call this%budobj%budterm(idx)%update_term(n, n2, q)
      end do
    end do
    !
    ! --
    text = '            RATE'
    idx = idx + 1
    maxlist = this%nmawwells
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
    if (this%iflowingwells > 0) then
      text = '         FW-RATE'
      idx = idx + 1
      maxlist = this%nmawwells
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
    text = '         STORAGE'
    idx = idx + 1
    maxlist = this%nmawwells
    naux = 1
    auxtxt(1) = '          VOLUME'
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%name_model, &
                                             maxlist, .false., .true., &
                                             naux, auxtxt)
    !
    ! --
    text = '        CONSTANT'
    idx = idx + 1
    maxlist = this%nmawwells
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
      maxlist = this%nmawwells
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
      text = '     RATE-TO-MVR'
      idx = idx + 1
      maxlist = this%nmawwells
      naux = 0
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux)
      !
      ! -- constant-head flow to mover
      text = ' CONSTANT-TO-MVR'
      idx = idx + 1
      maxlist = this%nmawwells
      naux = 0
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux)
      !
      ! -- flowing-well flow to mover
      if (this%iflowingwells > 0) then
        !
        ! --
        text = '  FW-RATE-TO-MVR'
        idx = idx + 1
        maxlist = this%nmawwells
        naux = 0
        call this%budobj%budterm(idx)%initialize(text, &
                                                 this%name_model, &
                                                 this%packName, &
                                                 this%name_model, &
                                                 this%packName, &
                                                 maxlist, .false., .false., &
                                                 naux)
      end if
    end if
    !
    ! -- auxiliary variable
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
    ! -- if maw flow for each reach are written to the listing file
    if (this%iprflow /= 0) then
      call this%budobj%flowtable_df(this%iout)
    end if
  end subroutine maw_setup_budobj

  !> @brief Copy flow terms into this%budobj
  !!
  !! terms include a combination of the following:
  !! gwf rate [flowing_well] [storage] constant_flow [frommvr tomvr tomvrcf [tomvrfw]] [aux]
  !<
  subroutine maw_fill_budobj(this)
    ! -- modules
    ! -- dummy
    class(MawType) :: this
    ! -- local
    integer(I4B) :: naux
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: n2
    integer(I4B) :: jpos
    integer(I4B) :: idx
    integer(I4B) :: ibnd
    real(DP) :: q
    real(DP) :: tmaw
    real(DP) :: bmaw
    real(DP) :: sat
    real(DP) :: qfact
    real(DP) :: q2
    real(DP) :: b
    real(DP) :: v
    ! -- formats
    !
    ! -- initialize counter
    idx = 0
    !
    ! -- GWF (LEAKAGE) and connection surface area (aux)
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%maxbound)
    ibnd = 1
    do n = 1, this%nmawwells
      do j = 1, this%ngwfnodes(n)
        jpos = this%get_jpos(n, j)
        n2 = this%get_gwfnode(n, j)
        tmaw = this%topscrn(jpos)
        bmaw = this%botscrn(jpos)
        call this%maw_calculate_saturation(n, j, n2, sat)
        this%qauxcbc(1) = DTWO * DPI * this%radius(n) * sat * (tmaw - bmaw)
        q = this%qleak(ibnd)
        call this%budobj%budterm(idx)%update_term(n, n2, q, this%qauxcbc)
        ibnd = ibnd + 1
      end do
    end do
    !
    ! -- RATE (WITHDRAWAL RATE)
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nmawwells)
    do n = 1, this%nmawwells
      q = this%ratesim(n)
      !
      ! -- adjust if well rate is an outflow
      if (this%imover == 1 .and. q < DZERO) then
        qfact = DONE
        if (this%qout(n) < DZERO) then
          qfact = q / this%qout(n)
        end if
        q = q + qfact * this%pakmvrobj%get_qtomvr(n)
      end if
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- FLOWING WELL
    if (this%iflowingwells > 0) then
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%nmawwells)
      do n = 1, this%nmawwells
        q = this%qfw(n)
        if (this%imover == 1) then
          qfact = DONE
          !
          ! -- adjust if well rate is an outflow
          if (this%qout(n) < DZERO) then
            qfact = q / this%qout(n)
          end if
          q = q + qfact * this%pakmvrobj%get_qtomvr(n)
        end if
        call this%budobj%budterm(idx)%update_term(n, n, q)
      end do
    end if
    !
    ! -- STORAGE (AND VOLUME AS AUX)
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nmawwells)
    do n = 1, this%nmawwells
      b = this%xsto(n) - this%bot(n)
      if (b < DZERO) then
        b = DZERO
      end if
      v = this%area(n) * b
      if (this%imawissopt /= 1) then
        q = this%qsto(n)
      else
        q = DZERO
      end if
      this%qauxcbc(1) = v
      call this%budobj%budterm(idx)%update_term(n, n, q, this%qauxcbc)
    end do
    !
    ! -- CONSTANT FLOW
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nmawwells)
    do n = 1, this%nmawwells
      q = this%qconst(n)
      !
      ! -- adjust if constant-flow rate is an outflow
      if (this%imover == 1 .and. q < DZERO) then
        qfact = DONE
        if (this%qout(n) < DZERO) then
          qfact = q / this%qout(n)
        end if
        q = q + qfact * this%pakmvrobj%get_qtomvr(n)
      end if
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- MOVER
    if (this%imover == 1) then
      !
      ! -- FROM MOVER
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%nmawwells)
      do n = 1, this%nmawwells
        if (this%iboundpak(n) == 0) then
          q = DZERO
        else
          q = this%pakmvrobj%get_qfrommvr(n)
        end if
        call this%budobj%budterm(idx)%update_term(n, n, q)
      end do
      !
      ! -- RATE TO MOVER
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%nmawwells)
      do n = 1, this%nmawwells
        q = this%pakmvrobj%get_qtomvr(n)
        if (q > DZERO) then
          q = -q
          q2 = this%ratesim(n)
          !
          ! -- adjust TO MOVER if well rate is outflow
          if (q2 < DZERO) then
            qfact = q2 / this%qout(n)
            q = q * qfact
          else
            q = DZERO
          end if
        end if
        call this%budobj%budterm(idx)%update_term(n, n, q)
      end do
      !
      ! -- CONSTANT TO MOVER
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%nmawwells)
      do n = 1, this%nmawwells
        q = this%pakmvrobj%get_qtomvr(n)
        if (q > DZERO) then
          q = -q
          q2 = this%qconst(n)
          ! -- adjust TO MOVER if well rate is outflow
          if (q2 < DZERO) then
            qfact = q2 / this%qout(n)
            q = q * qfact
          else
            q = DZERO
          end if
        end if
        call this%budobj%budterm(idx)%update_term(n, n, q)
      end do
      !
      ! -- FLOWING WELL TO MOVER
      if (this%iflowingwells > 0) then
        idx = idx + 1
        call this%budobj%budterm(idx)%reset(this%nmawwells)
        do n = 1, this%nmawwells
          q = this%pakmvrobj%get_qtomvr(n)
          if (q > DZERO) then
            q = -q
            q2 = this%ratesim(n)
            !
            ! -- adjust TO MOVER if well rate is outflow
            qfact = DONE
            if (this%qout(n) < DZERO) then
              qfact = this%qfw(n) / this%qout(n)
            end if
            q = q * qfact
          end if
          call this%budobj%budterm(idx)%update_term(n, n, q)
        end do
      end if

    end if
    !
    ! -- AUXILIARY VARIABLES
    naux = this%naux
    if (naux > 0) then
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%nmawwells)
      do n = 1, this%nmawwells
        q = DZERO
        call this%budobj%budterm(idx)%update_term(n, n, q, this%auxvar(:, n))
      end do
    end if
    !
    ! --Terms are filled, now accumulate them for this time step
    call this%budobj%accumulate_terms()
  end subroutine maw_fill_budobj

  !> @brief Set up the table object that is used to write the maw head data
  !!
  !! The terms listed here must correspond in number and order to the ones
  !! written to the head table in the maw_ot method.
  !<
  subroutine maw_setup_tableobj(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH, LENBUDTXT
    ! -- dummy
    class(MawType) :: this
    ! -- local
    integer(I4B) :: nterms
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    !
    ! -- setup well head table
    if (this%iprhed > 0) then
      !
      ! -- Determine the number of head table columns
      nterms = 2
      if (this%inamedbound == 1) nterms = nterms + 1
      !
      ! -- set up table title
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(adjustl(this%packName))//') HEADS FOR EACH CONTROL VOLUME'
      !
      ! -- set up head tableobj
      call table_cr(this%headtab, this%packName, title)
      call this%headtab%table_df(this%nmawwells, nterms, this%iout, &
                                 transient=.TRUE.)
      !
      ! -- Go through and set up table budget term
      if (this%inamedbound == 1) then
        text = 'NAME'
        call this%headtab%initialize_column(text, 20, alignment=TABLEFT)
      end if
      !
      ! -- reach number
      text = 'NUMBER'
      call this%headtab%initialize_column(text, 10, alignment=TABCENTER)
      !
      ! -- reach stage
      text = 'HEAD'
      call this%headtab%initialize_column(text, 12, alignment=TABCENTER)
    end if
  end subroutine maw_setup_tableobj

  !> @brief Get position of value in connection data
  !<
  function get_jpos(this, n, j) result(jpos)
    ! -- return variable
    integer(I4B) :: jpos
    ! -- dummy
    class(MawType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: j
    ! -- local
    !
    ! -- set jpos
    jpos = this%iaconn(n) + j - 1
  end function get_jpos

  !> @brief Get the gwfnode for connection
  !<
  function get_gwfnode(this, n, j) result(igwfnode)
    ! -- return variable
    integer(I4B) :: igwfnode
    ! -- dummy
    class(MawType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: j
    ! -- local
    integer(I4B) :: jpos
    !
    ! -- set jpos
    jpos = this%get_jpos(n, j)
    igwfnode = this%gwfnodes(jpos)
  end function get_gwfnode

  !> @brief Activate density terms
  !<
  subroutine maw_activate_density(this)
    ! -- dummy
    class(MawType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i, j
    ! -- formats
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
    write (this%iout, '(/1x,a)') 'DENSITY TERMS HAVE BEEN ACTIVATED FOR MAW &
      &PACKAGE: '//trim(adjustl(this%packName))
  end subroutine maw_activate_density

  !> @brief Activate viscosity terms
  !!
  !! Method to activate addition of viscosity terms for a MAW package reach.
  !<
  subroutine maw_activate_viscosity(this)
    ! -- modules
    use MemoryManagerModule, only: mem_reallocate
    ! -- dummy variables
    class(MawType), intent(inout) :: this !< MawType object
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
    write (this%iout, '(/1x,a)') 'VISCOSITY HAS BEEN ACTIVATED FOR MAW &
      &PACKAGE: '//trim(adjustl(this%packName))
  end subroutine maw_activate_viscosity

  !> @brief Calculate the groundwater-maw density exchange terms
  !!
  !! Arguments are as follows:
  !!     iconn       : maw-gwf connection number
  !!     hmaw        : maw head
  !!     hgwf        : gwf head
  !!     cond        : conductance
  !!     bmaw        : bottom elevation of this connection
  !!     flow        : calculated flow, updated here with density terms, + into maw
  !!     hcofterm    : head coefficient term
  !!     rhsterm     : right-hand-side value, updated here with density terms
  !!
  !! Member variable used here
  !!     denseterms  : shape (3, MAXBOUND), filled by buoyancy package
  !!                     col 1 is relative density of maw (densemaw / denseref)
  !!                     col 2 is relative density of gwf cell (densegwf / denseref)
  !!                     col 3 is elevation of gwf cell
  !!
  !! Upon return, amat and rhs for maw row should be updated as:
  !!    amat(idiag) = amat(idiag) - hcofterm
  !!    rhs(n) = rhs(n) + rhsterm
  !<
  subroutine maw_calculate_density_exchange(this, iconn, hmaw, hgwf, cond, &
                                            bmaw, flow, hcofterm, rhsterm)
    ! -- dummy
    class(MawType), intent(inout) :: this
    integer(I4B), intent(in) :: iconn
    real(DP), intent(in) :: hmaw
    real(DP), intent(in) :: hgwf
    real(DP), intent(in) :: cond
    real(DP), intent(in) :: bmaw
    real(DP), intent(inout) :: flow
    real(DP), intent(inout) :: hcofterm
    real(DP), intent(inout) :: rhsterm
    ! -- local
    real(DP) :: t
    real(DP) :: havg
    real(DP) :: rdensemaw
    real(DP) :: rdensegwf
    real(DP) :: rdenseavg
    real(DP) :: elevavg
    ! -- formats
    !
    ! -- assign relative density terms, return if zero which means not avail yet
    rdensemaw = this%denseterms(1, iconn)
    rdensegwf = this%denseterms(2, iconn)
    if (rdensegwf == DZERO) return
    !
    ! -- update rhsterm with density contribution
    if (hmaw > bmaw .and. hgwf > bmaw) then
      !
      ! -- hmaw and hgwf both above bmaw
      rdenseavg = DHALF * (rdensemaw + rdensegwf)
      !
      ! -- update rhsterm with first density term
      t = cond * (rdenseavg - DONE) * (hgwf - hmaw)
      rhsterm = rhsterm + t
      flow = flow + t
      !
      ! -- update rhterm with second density term
      havg = DHALF * (hgwf + hmaw)
      elevavg = this%denseterms(3, iconn)
      t = cond * (havg - elevavg) * (rdensegwf - rdensemaw)
      rhsterm = rhsterm + t
      flow = flow + t
    else if (hmaw > bmaw) then
      !
      ! -- if only hmaw is above bmaw, then increase correction term by density
      t = (rdensemaw - DONE) * rhsterm
      rhsterm = rhsterm + t
      !
    else if (hgwf > bmaw) then
      !
      ! -- if only hgwf is above bmaw, then increase correction term by density
      t = (rdensegwf - DONE) * rhsterm
      rhsterm = rhsterm + t
      !
    else
      !
      ! -- Flow should be zero so do nothing
    end if
  end subroutine maw_calculate_density_exchange

end module MawModule
