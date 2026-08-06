! -- Uzf module
module UzfModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DEM6, DEM4, DEM2, DEM1, DHALF, &
                             DONE, DHUNDRED, &
                             LINELENGTH, LENFTYPE, LENPACKAGENAME, &
                             LENBOUNDNAME, LENBUDTXT, LENPAKLOC, DNODATA, &
                             NAMEDBOUNDFLAG, MAXCHARLEN, &
                             DHNOFLO, DHDRY, &
                             TABLEFT, TABCENTER, TABRIGHT, &
                             TABSTRING, TABUCSTRING, TABINTEGER, TABREAL
  use MemoryManagerModule, only: mem_allocate, mem_reallocate, mem_setptr, &
                                 mem_deallocate
  use MemoryHelperModule, only: create_mem_path
  use SparseModule, only: sparsematrix
  use BndModule, only: BndType
  use UzfCellGroupModule, only: UzfCellGroupType
  use BudgetObjectModule, only: BudgetObjectType, budgetobject_cr
  use BaseDisModule, only: DisBaseType
  use ObserveModule, only: ObserveType
  use ObsModule, only: ObsType
  use InputOutputModule, only: URWORD
  use SimVariablesModule, only: errmsg, warnmsg
  use SimModule, only: count_errors, store_error, store_error_unit, &
                       deprecation_warning
  use BlockParserModule, only: BlockParserType
  use TableModule, only: TableType, table_cr
  use UzfETUtilModule
  use MatrixBaseModule

  implicit none

  character(len=LENFTYPE) :: ftype = 'UZF'
  character(len=LENPACKAGENAME) :: text = '       UZF CELLS'

  private
  public :: uzf_create
  public :: UzfType

  type, extends(BndType) :: UzfType
    ! output integers
    integer(I4B), pointer :: iprwcont => null()
    integer(I4B), pointer :: iwcontout => null()
    integer(I4B), pointer :: ibudgetout => null()
    integer(I4B), pointer :: ibudcsv => null() !< unit number for csv budget output file
    integer(I4B), pointer :: ipakcsv => null()
    !
    type(BudgetObjectType), pointer :: budobj => null()
    integer(I4B), pointer :: bditems => null() !< number of budget items
    integer(I4B), pointer :: nbdtxt => null() !< number of budget text items
    character(len=LENBUDTXT), dimension(:), pointer, &
      contiguous :: bdtxt => null() !< budget items written to cbc file
    character(len=LENBOUNDNAME), dimension(:), pointer, &
      contiguous :: uzfname => null()
    !
    ! -- uzf table objects
    type(TableType), pointer :: pakcsvtab => null()
    !
    ! -- uzf kinematic object
    type(UzfCellGroupType), pointer :: uzfobj => null()
    !
    ! -- pointer to gwf variables
    integer(I4B), pointer :: gwfiss => null()
    real(DP), dimension(:), pointer, contiguous :: gwfhcond => null()
    !
    ! -- uzf data
    integer(I4B), pointer :: nwav_pvar => null()
    integer(I4B), pointer :: ntrail_pvar => null()
    integer(I4B), pointer :: nsets => null()
    integer(I4B), pointer :: nodes => null()
    integer(I4B), pointer :: readflag => null()
    integer(I4B), pointer :: ietflag => null() !< et flag, 0 is off, 1 or 2 are different types
    integer(I4B), pointer :: igwetflag => null()
    integer(I4B), pointer :: iseepflag => null()
    integer(I4B), pointer :: imaxcellcnt => null()
    integer(I4B), pointer :: iuzf2uzf => null()
    ! -- integer vectors
    integer(I4B), dimension(:), pointer, contiguous :: igwfnode => null()
    integer(I4B), dimension(:), pointer, contiguous :: ia => null()
    integer(I4B), dimension(:), pointer, contiguous :: ja => null()
    ! -- double precision output vectors
    real(DP), dimension(:), pointer, contiguous :: appliedinf => null()
    real(DP), dimension(:), pointer, contiguous :: rejinf => null()
    real(DP), dimension(:), pointer, contiguous :: rejinf0 => null()
    real(DP), dimension(:), pointer, contiguous :: rejinftomvr => null()
    real(DP), dimension(:), pointer, contiguous :: infiltration => null()
    real(DP), dimension(:), pointer, contiguous :: gwet_pvar => null()
    real(DP), dimension(:), pointer, contiguous :: uzet => null()
    real(DP), dimension(:), pointer, contiguous :: gwd => null()
    real(DP), dimension(:), pointer, contiguous :: gwd0 => null()
    real(DP), dimension(:), pointer, contiguous :: gwdtomvr => null()
    real(DP), dimension(:), pointer, contiguous :: rch => null()
    real(DP), dimension(:), pointer, contiguous :: rch0 => null()
    real(DP), dimension(:), pointer, contiguous :: qsto => null() !< change in stored mobile water per time for this time step
    real(DP), dimension(:), pointer, contiguous :: wcnew => null() !< water content for this time step
    real(DP), dimension(:), pointer, contiguous :: wcold => null() !< water content for previous time step
    !
    ! -- timeseries aware package variables; these variables with
    !    _pvar have uzfobj counterparts
    real(DP), dimension(:), pointer, contiguous :: sinf_pvar => null()
    real(DP), dimension(:), pointer, contiguous :: pet_pvar => null()
    real(DP), dimension(:), pointer, contiguous :: extdp => null()
    real(DP), dimension(:), pointer, contiguous :: extwc_pvar => null()
    real(DP), dimension(:), pointer, contiguous :: ha_pvar => null()
    real(DP), dimension(:), pointer, contiguous :: hroot_pvar => null()
    real(DP), dimension(:), pointer, contiguous :: rootact_pvar => null()
    !
    ! -- aux variable
    real(DP), dimension(:, :), pointer, contiguous :: uauxvar => null()
    !
    ! -- convergence check
    integer(I4B), pointer :: iconvchk => null()
    !
    ! formulate variables
    real(DP), dimension(:), pointer, contiguous :: deriv => null()
    !
    ! budget variables
    real(DP), pointer :: totfluxtot => null()
    integer(I4B), pointer :: issflag => null()
    integer(I4B), pointer :: issflagold => null()
    integer(I4B), pointer :: istocb => null()
    !
    ! -- uzf cbc budget items
    integer(I4B), pointer :: cbcauxitems => NULL()
    character(len=16), dimension(:), pointer, contiguous :: cauxcbc => NULL()
    real(DP), dimension(:), pointer, contiguous :: qauxcbc => null()

  contains

    procedure :: uzf_allocate_arrays
    procedure :: uzf_allocate_scalars
    procedure :: bnd_options => uzf_options
    procedure :: read_dimensions => uzf_readdimensions
    procedure :: bnd_ar => uzf_ar
    procedure :: bnd_rp => uzf_rp
    procedure :: bnd_ad => uzf_ad
    procedure :: bnd_cf => uzf_cf
    procedure :: bnd_cc => uzf_cc
    procedure :: bnd_cq => uzf_cq
    procedure :: bnd_bd => uzf_bd
    procedure :: bnd_ot_model_flows => uzf_ot_model_flows
    procedure :: bnd_ot_package_flows => uzf_ot_package_flows
    procedure :: bnd_ot_dv => uzf_ot_dv
    procedure :: bnd_ot_bdsummary => uzf_ot_bdsummary
    procedure :: bnd_fc => uzf_fc
    procedure :: bnd_fn => uzf_fn
    procedure :: bnd_da => uzf_da
    procedure :: define_listlabel
    !
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => uzf_obs_supported
    procedure, public :: bnd_df_obs => uzf_df_obs
    procedure, public :: bnd_rp_obs => uzf_rp_obs
    procedure, public :: bnd_bd_obs => uzf_bd_obs
    !
    ! -- methods specific for uzf
    procedure, private :: uzf_solve
    procedure, private :: read_cell_properties
    procedure, private :: print_cell_properties
    procedure, private :: check_cell_area
    !
    ! -- budget
    procedure, private :: uzf_setup_budobj
    procedure, private :: uzf_fill_budobj

  end type UzfType

  interface

    module subroutine read_cell_properties(this)
      class(UzfType), intent(inout) :: this
    end subroutine read_cell_properties

    module subroutine print_cell_properties(this)
      class(UzfType), intent(inout) :: this
    end subroutine print_cell_properties

    module subroutine check_cell_area(this)
      class(UzfType) :: this
    end subroutine check_cell_area

    module logical function uzf_obs_supported(this)
      class(UzfType) :: this
    end function uzf_obs_supported

    module subroutine uzf_df_obs(this)
      class(UzfType) :: this
    end subroutine uzf_df_obs

    module subroutine uzf_rp_obs(this)
      class(UzfType), intent(inout) :: this
    end subroutine uzf_rp_obs

    module subroutine uzf_bd_obs(this)
      class(UzfType) :: this
    end subroutine uzf_bd_obs

    module subroutine uzf_setup_budobj(this)
      class(UzfType) :: this
    end subroutine uzf_setup_budobj

    module subroutine uzf_fill_budobj(this)
      class(UzfType) :: this
    end subroutine uzf_fill_budobj

  end interface

contains

  !> @brief Create a New UZF Package and point packobj to the new package
  !<
  subroutine uzf_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    ! -- local
    type(UzfType), pointer :: uzfobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (uzfobj)
    packobj => uzfobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype)
    packobj%text = text
    !
    ! -- allocate scalars
    call uzfobj%uzf_allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()
    !
    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 1
    packobj%iscloc = 0 ! not supported
    packobj%isadvpak = 1
    packobj%ictMemPath = create_mem_path(namemodel, 'NPF')
  end subroutine uzf_create

  !> @brief Allocate and Read
  !<
  subroutine uzf_ar(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr, mem_reallocate
    ! -- dummy
    class(UzfType), intent(inout) :: this
    ! -- local
    integer(I4B) :: n !< uzf cell number
    integer(I4B) :: node !< gwf node number
    real(DP) :: hgwf
    !
    call this%obs%obs_ar()
    !
    ! -- call standard BndType allocate scalars
    call this%BndType%allocate_arrays()
    !
    ! -- set pointers now that data is available
    call mem_setptr(this%gwfhcond, 'CONDSAT', create_mem_path(this%name_model, &
                                                              'NPF'))
    call mem_setptr(this%gwfiss, 'ISS', create_mem_path(this%name_model))
    !
    ! -- set boundname for each connection
    if (this%inamedbound /= 0) then
      do n = 1, this%nodes
        this%boundname(n) = this%uzfname(n)
      end do
    end if
    !
    ! -- copy boundname into boundname_cst
    call this%copy_boundname()
    !
    ! -- copy igwfnode into nodelist and set water table
    do n = 1, this%nodes
      this%nodelist(n) = this%igwfnode(n)
      node = this%igwfnode(n)
      hgwf = this%xnew(node)
      call this%uzfobj%sethead(n, hgwf)
    end do
    !
    ! -- setup pakmvrobj
    if (this%imover /= 0) then
      allocate (this%pakmvrobj)
      call this%pakmvrobj%ar(this%maxbound, this%maxbound, this%memoryPath)
    end if
  end subroutine uzf_ar

  !> @brief Allocate arrays used for uzf
  !<
  subroutine uzf_allocate_arrays(this)
    ! -- dummy
    class(UzfType), intent(inout) :: this
    ! -- local
    integer(I4B) :: n !< uzf cell number
    integer(I4B) :: node !< gwf node number
    integer(I4B) :: j
    !
    ! -- call standard BndType allocate scalars (now done from AR)
    !call this%BndType%allocate_arrays()
    !
    ! -- allocate uzf specific arrays
    call mem_allocate(this%igwfnode, this%nodes, 'IGWFNODE', this%memoryPath)
    call mem_allocate(this%appliedinf, this%nodes, 'APPLIEDINF', this%memoryPath)
    call mem_allocate(this%rejinf, this%nodes, 'REJINF', this%memoryPath)
    call mem_allocate(this%rejinf0, this%nodes, 'REJINF0', this%memoryPath)
    call mem_allocate(this%rejinftomvr, this%nodes, 'REJINFTOMVR', &
                      this%memoryPath)
    call mem_allocate(this%infiltration, this%nodes, 'INFILTRATION', &
                      this%memoryPath)
    call mem_allocate(this%gwet_pvar, this%nodes, 'GWET_PVAR', this%memoryPath)
    call mem_allocate(this%uzet, this%nodes, 'UZET', this%memoryPath)
    call mem_allocate(this%gwd, this%nodes, 'GWD', this%memoryPath)
    call mem_allocate(this%gwd0, this%nodes, 'GWD0', this%memoryPath)
    call mem_allocate(this%gwdtomvr, this%nodes, 'GWDTOMVR', this%memoryPath)
    call mem_allocate(this%rch, this%nodes, 'RCH', this%memoryPath)
    call mem_allocate(this%rch0, this%nodes, 'RCH0', this%memoryPath)
    call mem_allocate(this%qsto, this%nodes, 'QSTO', this%memoryPath)
    call mem_allocate(this%deriv, this%nodes, 'DERIV', this%memoryPath)
    !
    ! -- integer vectors
    call mem_allocate(this%ia, this%dis%nodes + 1, 'IA', this%memoryPath)
    call mem_allocate(this%ja, this%nodes, 'JA', this%memoryPath)
    !
    ! -- allocate timeseries aware variables
    call mem_allocate(this%sinf_pvar, this%nodes, 'SINF_PVAR', this%memoryPath)
    call mem_allocate(this%pet_pvar, this%nodes, 'PET_PVAR', this%memoryPath)
    call mem_allocate(this%extdp, this%nodes, 'EXDP_PVAR', this%memoryPath)
    call mem_allocate(this%extwc_pvar, this%nodes, 'EXTWC_PVAR', this%memoryPath)
    call mem_allocate(this%ha_pvar, this%nodes, 'HA_PVAR', this%memoryPath)
    call mem_allocate(this%hroot_pvar, this%nodes, 'HROOT_PVAR', this%memoryPath)
    call mem_allocate(this%rootact_pvar, this%nodes, 'ROOTACT_PVAR', &
                      this%memoryPath)
    call mem_allocate(this%uauxvar, this%naux, this%nodes, 'UAUXVAR', &
                      this%memoryPath)
    !
    ! -- initialize
    do n = 1, this%nodes
      this%appliedinf(n) = DZERO
      this%rejinf(n) = DZERO
      this%rejinf0(n) = DZERO
      this%rejinftomvr(n) = DZERO
      this%gwet_pvar(n) = DZERO
      this%uzet(n) = DZERO
      this%gwd(n) = DZERO
      this%gwd0(n) = DZERO
      this%gwdtomvr(n) = DZERO
      this%rch(n) = DZERO
      this%rch0(n) = DZERO
      this%qsto(n) = DZERO
      this%deriv(n) = DZERO
      ! -- integer variables
      this%ja(n) = 0
      ! -- timeseries aware variables
      this%sinf_pvar(n) = DZERO
      this%pet_pvar(n) = DZERO
      this%extdp(n) = DZERO
      this%extwc_pvar(n) = DZERO
      this%ha_pvar(n) = DZERO
      this%hroot_pvar(n) = DZERO
      this%rootact_pvar(n) = DZERO
      do j = 1, this%naux
        if (this%iauxmultcol > 0 .and. j == this%iauxmultcol) then
          this%uauxvar(j, n) = DONE
        else
          this%uauxvar(j, n) = DZERO
        end if
      end do
    end do
    do node = 1, this%dis%nodes + 1
      this%ia(node) = 0
    end do
    !
    ! -- allocate and initialize character array for budget text
    allocate (this%bdtxt(this%nbdtxt))
    this%bdtxt(1) = '         UZF-INF'
    this%bdtxt(2) = '       UZF-GWRCH'
    this%bdtxt(3) = '         UZF-GWD'
    this%bdtxt(4) = '        UZF-GWET'
    this%bdtxt(5) = '  UZF-GWD TO-MVR'
    !
    ! -- allocate and initialize watercontent arrays
    call mem_allocate(this%wcnew, this%nodes, 'WCNEW', this%memoryPath)
    call mem_allocate(this%wcold, this%nodes, 'WCOLD', this%memoryPath)
    do n = 1, this%nodes
      this%wcnew(n) = DZERO
      this%wcold(n) = DZERO
    end do
    !
    ! -- allocate character array for aux budget text
    allocate (this%cauxcbc(this%cbcauxitems))
    allocate (this%uzfname(this%nodes))
    !
    ! -- allocate and initialize qauxcbc
    call mem_allocate(this%qauxcbc, this%cbcauxitems, 'QAUXCBC', this%memoryPath)
    do n = 1, this%cbcauxitems
      this%qauxcbc(n) = DZERO
    end do
  end subroutine uzf_allocate_arrays

  !> @brief Set options specific to UzfType
  !!
  !! Overrides BoundaryPackageType%child_class_options
  !<
  subroutine uzf_options(this, option, found)
    ! -- modules
    use ConstantsModule, only: DZERO, MNORMAL
    use OpenSpecModule, only: access, form
    use SimModule, only: store_error
    use InputOutputModule, only: urword, getunit, assign_iounit, openfile
    implicit none
    ! -- dummy
    class(uzftype), intent(inout) :: this
    character(len=*), intent(inout) :: option
    logical, intent(inout) :: found
    ! -- local
    character(len=MAXCHARLEN) :: fname, keyword
    ! -- formats
    character(len=*), parameter :: fmtnotfound = &
      &"(4x, 'NO UZF OPTIONS WERE FOUND.')"
    character(len=*), parameter :: fmtet = &
      "(4x, 'ET WILL BE SIMULATED WITHIN UZ AND GW ZONES, WITH LINEAR ', &
      &'GWET IF OPTION NOT SPECIFIED OTHERWISE.')"
    character(len=*), parameter :: fmtgwetlin = &
      &"(4x, 'GROUNDWATER ET FUNCTION WILL BE LINEAR.')"
    character(len=*), parameter :: fmtgwetsquare = &
      &"(4x, 'GROUNDWATER ET FUNCTION WILL BE SQUARE WITH SMOOTHING.')"
    character(len=*), parameter :: fmtgwseepout = &
      &"(4x, 'GROUNDWATER DISCHARGE TO LAND SURFACE WILL BE SIMULATED.')"
    character(len=*), parameter :: fmtuzetwc = &
      &"(4x, 'UNSATURATED ET FUNCTION OF WATER CONTENT.')"
    character(len=*), parameter :: fmtuzetae = &
      &"(4x, 'UNSATURATED ET FUNCTION OF AIR ENTRY PRESSURE.')"
    character(len=*), parameter :: fmtuznlay = &
      &"(4x, 'UNSATURATED FLOW WILL BE SIMULATED SEPARATELY IN EACH LAYER.')"
    character(len=*), parameter :: fmtuzfbin = &
      "(4x, 'UZF ', 1x, a, 1x, ' WILL BE SAVED TO FILE: ', &
      &a, /4x, 'OPENED ON UNIT: ', I0)"
    character(len=*), parameter :: fmtuzfopt = &
      &"(4x, 'UZF ', a, ' VALUE (',g15.7,') SPECIFIED.')"
    !
    !
    found = .true.
    select case (option)
      !case ('PRINT_WATER-CONTENT')
      !  this%iprwcont = 1
      !  write(this%iout,'(4x,a)') trim(adjustl(this%text))// &
      !    ' WATERCONTENT WILL BE PRINTED TO LISTING FILE.'
    case ('WATER_CONTENT')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        this%iwcontout = getunit()
        call openfile(this%iwcontout, this%iout, fname, 'DATA(BINARY)', &
                      form, access, 'REPLACE', mode_opt=MNORMAL)
        write (this%iout, fmtuzfbin) 'WATER-CONTENT', trim(adjustl(fname)), &
          this%iwcontout
      else
        call store_error('OPTIONAL WATER_CONTENT KEYWORD &
                         &MUST BE FOLLOWED BY FILEOUT')
      end if
    case ('BUDGET')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        call assign_iounit(this%ibudgetout, this%inunit, "BUDGET fileout")
        call openfile(this%ibudgetout, this%iout, fname, 'DATA(BINARY)', &
                      form, access, 'REPLACE', mode_opt=MNORMAL)
        write (this%iout, fmtuzfbin) 'BUDGET', trim(adjustl(fname)), &
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
        write (this%iout, fmtuzfbin) 'BUDGET CSV', trim(adjustl(fname)), &
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
        write (this%iout, fmtuzfbin) 'PACKAGE_CONVERGENCE', &
          trim(adjustl(fname)), this%ipakcsv
      else
        call store_error('OPTIONAL PACKAGE_CONVERGENCE KEYWORD MUST BE '// &
                         'FOLLOWED BY FILEOUT')
      end if
    case ('SIMULATE_ET')
      this%ietflag = 1 !default
      this%igwetflag = 0
      write (this%iout, fmtet)
    case ('LINEAR_GWET')
      this%igwetflag = 1
      write (this%iout, fmtgwetlin)
    case ('SQUARE_GWET')
      this%igwetflag = 2
      write (this%iout, fmtgwetsquare)
    case ('SIMULATE_GWSEEP')
      this%iseepflag = 1
      write (this%iout, fmtgwseepout)
      !
      ! -- Create warning message
      write (warnmsg, '(a)') &
        'USE DRN PACKAGE TO SIMULATE GROUNDWATER DISCHARGE TO LAND SURFACE '// &
        'INSTEAD'
      !
      ! -- Create deprecation warning
      call deprecation_warning('OPTIONS', 'SIMULATE_GWSEEP', '6.5.0', &
                               warnmsg, this%parser%GetUnit())
    case ('UNSAT_ETWC')
      this%ietflag = 1
      write (this%iout, fmtuzetwc)
    case ('UNSAT_ETAE')
      this%ietflag = 2
      write (this%iout, fmtuzetae)
    case ('MOVER')
      this%imover = 1
      !
      ! -- right now these are options that are available but may not be available in
      !    the release (or in documentation)
    case ('DEV_NO_FINAL_CHECK')
      call this%parser%DevOpt()
      this%iconvchk = 0
      write (this%iout, '(4x,a)') &
        'A FINAL CONVERGENCE CHECK OF THE CHANGE IN UZF RECHARGE &
        &WILL NOT BE MADE'
      !case('DEV_MAXIMUM_PERCENT_DIFFERENCE')
      !  call this%parser%DevOpt()
      !  r = this%parser%GetDouble()
      !  if (r > DZERO) then
      !    this%pdmax = r
      !    write(this%iout, fmtuzfopt) 'MAXIMUM_PERCENT_DIFFERENCE', this%pdmax
      !  else
      !    write(this%iout, fmtuzfopt) 'INVALID MAXIMUM_PERCENT_DIFFERENCE', r
      !    write(this%iout, fmtuzfopt) 'USING DEFAULT MAXIMUM_PERCENT_DIFFERENCE', this%pdmax
      !  end if
    case default
      ! -- No options found
      found = .false.
    end select
  end subroutine uzf_options

  !> @brief Set dimensions specific to UzfType
  !<
  subroutine uzf_readdimensions(this)
    ! -- modules
    use InputOutputModule, only: urword
    use SimModule, only: store_error, count_errors
    class(uzftype), intent(inout) :: this
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    !
    ! -- initialize dimensions to -1
    this%nodes = -1
    this%ntrail_pvar = 0
    this%nsets = 0
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
        case ('NUZFCELLS')
          this%nodes = this%parser%GetInteger()
          write (this%iout, '(4x,a,i0)') 'NUZFCELLS = ', this%nodes
        case ('NTRAILWAVES')
          this%ntrail_pvar = this%parser%GetInteger()
          write (this%iout, '(4x,a,i0)') 'NTRAILWAVES = ', this%ntrail_pvar
        case ('NWAVESETS')
          this%nsets = this%parser%GetInteger()
          write (this%iout, '(4x,a,i0)') 'NTRAILSETS = ', this%nsets
        case default
          write (errmsg, '(a,a)') &
            'Unknown '//trim(this%text)//' dimension: ', trim(keyword)
        end select
      end do
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' DIMENSIONS'
    else
      call store_error('Required dimensions block not found.')
    end if
    !
    ! -- increment maxbound
    this%maxbound = this%maxbound + this%nodes
    this%nbound = this%maxbound
    !
    ! -- verify dimensions were set
    if (this%nodes <= 0) then
      write (errmsg, '(a)') &
        'NUZFCELLS was not specified or was specified incorrectly.'
      call store_error(errmsg)
    end if

    if (this%ntrail_pvar <= 1) then
      write (errmsg, '(a)') &
        'NTRAILWAVES must be greater than 1. A value of 7 is recommended.'
      call store_error(errmsg)
    end if
    !
    if (this%nsets <= 0) then
      write (errmsg, '(a)') &
        'NWAVESETS was not specified or was specified incorrectly.'
      call store_error(errmsg)
    end if
    !
    ! -- terminate if there are dimension errors
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- set the number of waves
    this%nwav_pvar = this%ntrail_pvar * this%nsets
    !
    ! -- Call define_listlabel to construct the list label that is written
    !    when PRINT_INPUT option is used.
    call this%define_listlabel()
    !
    ! -- Allocate arrays in package superclass
    call this%uzf_allocate_arrays()
    !
    ! -- initialize uzf group object
    allocate (this%uzfobj)
    call this%uzfobj%init(this%nodes, this%nwav_pvar, this%memoryPath)
    !
    !--Read uzf cell properties and set values
    call this%read_cell_properties()
    !
    ! -- print cell data
    if (this%iprpak /= 0) then
      call this%print_cell_properties()
    end if
    !
    ! -- setup the budget object
    call this%uzf_setup_budobj()
  end subroutine uzf_readdimensions

  !> @brief Read stress data
  !!
  !!  - check if bc changes
  !!  - read new bc for stress period
  !!  - set kinematic variables to bc values
  !<
  subroutine uzf_rp(this)
    ! -- modules
    use TdisModule, only: kper, nper
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    use InputOutputModule, only: urword
    use SimModule, only: store_error, count_errors
    ! -- dummy
    class(UzfType), intent(inout) :: this
    ! -- local
    character(len=LENBOUNDNAME) :: bndName
    character(len=LINELENGTH) :: text
    character(len=LINELENGTH) :: line
    logical :: isfound
    logical :: endOfBlock
    integer(I4B) :: n
    integer(I4B) :: j
    integer(I4B) :: jj
    integer(I4B) :: ierr
    real(DP), pointer :: bndElem => null()
    ! -- table output
    character(len=20) :: cellid
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: tag
    integer(I4B) :: ntabrows
    integer(I4B) :: ntabcols
    integer(I4B) :: node
    !-- formats
    character(len=*), parameter :: fmtlsp = &
      &"(1X,/1X,'REUSING ',A,'S FROM LAST STRESS PERIOD')"
    character(len=*), parameter :: fmtblkerr = &
      &"('Looking for BEGIN PERIOD iper.  Found ', a, ' instead.')"
    character(len=*), parameter :: fmtisvflow = &
      "(4x,'CELL-BY-CELL FLOW INFORMATION WILL BE SAVED TO BINARY FILE &
      &WHENEVER ICBCFL IS NOT ZERO.')"
    character(len=*), parameter :: fmtflow = &
      &"(4x, 'FLOWS WILL BE SAVED TO FILE: ', a, /4x, 'OPENED ON UNIT: ', I7)"
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
    ! -- set steady-state flag based on gwfiss
    this%issflag = this%gwfiss
    !
    ! -- read data if ionper == kper
    if (this%ionper == kper) then
      !
      ! -- write header
      if (this%iprpak /= 0) then
        !
        ! -- setup inputtab tableobj
        !
        ! -- table dimensions
        ntabrows = 1
        ntabcols = 3
        if (this%ietflag /= 0) then
          ntabcols = ntabcols + 3
          if (this%ietflag == 2) then
            ntabcols = ntabcols + 3
          end if
        end if
        if (this%inamedbound == 1) then
          ntabcols = ntabcols + 1
        end if
        !
        ! -- initialize table and define columns
        title = trim(adjustl(this%text))//' PACKAGE ('// &
                trim(adjustl(this%packName))//') DATA FOR PERIOD'
        write (title, '(a,1x,i6)') trim(adjustl(title)), kper
        call table_cr(this%inputtab, this%packName, title)
        call this%inputtab%table_df(ntabrows, ntabcols, this%iout, &
                                    finalize=.FALSE.)
        tag = 'NUMBER'
        call this%inputtab%initialize_column(tag, 10)
        tag = 'CELLID'
        call this%inputtab%initialize_column(tag, 20, alignment=TABLEFT)
        tag = 'FINF'
        call this%inputtab%initialize_column(tag, 12)
        if (this%ietflag /= 0) then
          tag = 'PET'
          call this%inputtab%initialize_column(tag, 12)
          tag = 'EXTDEP'
          call this%inputtab%initialize_column(tag, 12)
          tag = 'EXTWC'
          call this%inputtab%initialize_column(tag, 12)
          if (this%ietflag == 2) then
            tag = 'HA'
            call this%inputtab%initialize_column(tag, 12)
            tag = 'HROOT'
            call this%inputtab%initialize_column(tag, 12)
            tag = 'ROOTACT'
            call this%inputtab%initialize_column(tag, 12)
          end if
        end if
        if (this%inamedbound == 1) then
          tag = 'BOUNDNAME'
          call this%inputtab%initialize_column(tag, LENBOUNDNAME, &
                                               alignment=TABLEFT)
        end if
      end if
      !
      ! -- read the stress period data
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        !
        ! -- check for valid uzf node
        n = this%parser%GetInteger()
        if (n < 1 .or. n > this%nodes) then
          tag = trim(adjustl(this%text))//' PACKAGE ('// &
                trim(adjustl(this%packName))//') DATA FOR PERIOD'
          write (tag, '(a,1x,i0)') trim(adjustl(tag)), kper

          write (errmsg, '(a,a,i0,1x,a,i0,a)') &
            trim(adjustl(tag)), ': UZFNO ', n, &
            'must be greater than 0 and less than or equal to ', this%nodes, '.'
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- Setup boundname
        if (this%inamedbound > 0) then
          bndName = this%boundname(n)
        else
          bndName = ''
        end if
        !
        ! -- FINF
        call this%parser%GetStringCaps(text)
        jj = 1 ! For SINF
        bndElem => this%sinf_pvar(n)
        call read_value_or_time_series_adv(text, n, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'SINF')
        !
        ! -- PET
        call this%parser%GetStringCaps(text)
        jj = 1 ! For PET
        bndElem => this%pet_pvar(n)
        call read_value_or_time_series_adv(text, n, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'PET')
        !
        ! -- EXTD
        call this%parser%GetStringCaps(text)
        jj = 1 ! For EXTDP
        bndElem => this%extdp(n)
        call read_value_or_time_series_adv(text, n, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'EXTDP')
        !
        ! -- EXTWC
        call this%parser%GetStringCaps(text)
        jj = 1 ! For EXTWC
        bndElem => this%extwc_pvar(n)
        call read_value_or_time_series_adv(text, n, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'EXTWC')
        !
        ! -- HA
        call this%parser%GetStringCaps(text)
        jj = 1 ! For HA
        bndElem => this%ha_pvar(n)
        call read_value_or_time_series_adv(text, n, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'HA')
        !
        ! -- HROOT
        call this%parser%GetStringCaps(text)
        jj = 1 ! For HROOT
        bndElem => this%hroot_pvar(n)
        call read_value_or_time_series_adv(text, n, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'HROOT')
        !
        ! -- ROOTACT
        call this%parser%GetStringCaps(text)
        jj = 1 ! For ROOTACT
        bndElem => this%rootact_pvar(n)
        call read_value_or_time_series_adv(text, n, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'ROOTACT')
        !
        ! -- read auxiliary variables
        do j = 1, this%naux
          call this%parser%GetStringCaps(text)
          bndElem => this%uauxvar(j, n)
          call read_value_or_time_series_adv(text, n, j, bndElem, this%packName, &
                                             'AUX', this%tsManager, this%iprpak, &
                                             this%auxname(j))
        end do
        !
        ! -- write line
        if (this%iprpak /= 0) then
          !
          ! -- get cellid
          node = this%igwfnode(n)
          if (node > 0) then
            call this%dis%noder_to_string(node, cellid)
          else
            cellid = 'none'
          end if
          !
          ! -- write data to the table
          call this%inputtab%add_term(n)
          call this%inputtab%add_term(cellid)
          call this%inputtab%add_term(this%sinf_pvar(n))
          if (this%ietflag /= 0) then
            call this%inputtab%add_term(this%pet_pvar(n))
            call this%inputtab%add_term(this%extdp(n))
            call this%inputtab%add_term(this%extwc_pvar(n))
            if (this%ietflag == 2) then
              call this%inputtab%add_term(this%ha_pvar(n))
              call this%inputtab%add_term(this%hroot_pvar(n))
              call this%inputtab%add_term(this%rootact_pvar(n))
            end if
          end if
          if (this%inamedbound == 1) then
            call this%inputtab%add_term(this%boundname(n))
          end if
        end if

      end do
      !
      ! -- finalize the table
      if (this%iprpak /= 0) then
        call this%inputtab%finalize_table()
      end if
      !
      ! -- using stress period data from the previous stress period
    else
      write (this%iout, fmtlsp) trim(this%filtyp)
    end if
    !
    ! -- write summary of uzf stress period error messages
    ierr = count_errors()
    if (ierr > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- set wave data for first stress period and second that follows SS
    if ((this%issflag == 0 .AND. kper == 1) .or. &
        (kper == 2 .AND. this%issflagold == 1)) then
      do n = 1, this%nodes
        call this%uzfobj%setwaves(n)
      end do
    end if
    !
    ! -- Initialize the water content
    if (kper == 1) then
      do n = 1, this%nodes
        this%wcnew(n) = this%uzfobj%get_wcnew(n)
      end do
    end if
    !
    ! -- Save old ss flag
    this%issflagold = this%issflag
  end subroutine uzf_rp

  !> @brief Advance UZF Package
  !<
  subroutine uzf_ad(this)
    ! -- modules
    use SimVariablesModule, only: iFailedStepRetry
    ! -- dummy
    class(UzfType) :: this
    ! -- locals
    integer(I4B) :: n !< uzf cell number
    integer(I4B) :: ivertflag
    integer(I4B) :: iaux
    real(DP) :: rval1, rval2, rval3
    !
    ! -- Advance the time series
    call this%TsManager%ad()
    !
    ! -- update auxiliary variables by copying from the derived-type time
    !    series variable into the bndpackage auxvar variable so that this
    !    information is properly written to the GWF budget file
    if (this%naux > 0) then
      do n = 1, this%maxbound
        do iaux = 1, this%naux
          if (this%noupdateauxvar(iaux) /= 0) cycle
          this%auxvar(iaux, n) = this%uauxvar(iaux, n)
        end do
      end do
    end if
    !
    ! -- Update or restore state
    if (iFailedStepRetry == 0) then
      !
      ! -- reset old water content to new water content
      do n = 1, this%nodes
        this%wcold(n) = this%wcnew(n)
      end do
    else
      !
      ! -- Copy wcold back into wcnew as this is a retry of this time step.
      !    Note that there is no need to reset the waves as they are not
      !    advanced to their new state until the _ot() method is called,
      !    and that doesn't happen until a successful solution is obtained.
      do n = 1, this%nodes
        this%wcnew(n) = this%wcold(n)
      end do
    end if
    !
    ! -- advance each uzf obj
    do n = 1, this%nodes
      call this%uzfobj%advance(n)
    end do
    !
    ! -- update uzf objects with timeseries aware variables
    do n = 1, this%nodes
      !
      ! -- Set ivertflag
      ivertflag = this%uzfobj%cell_below(n)
      !
      ! -- recalculate uzfarea
      if (this%iauxmultcol > 0) then
        rval1 = this%uauxvar(this%iauxmultcol, n)
        call this%uzfobj%setdatauzfarea(n, rval1)
      end if
      !
      ! -- FINF
      rval1 = this%sinf_pvar(n)
      call this%uzfobj%setdatafinf(n, rval1)
      !
      ! -- PET, EXTDP
      rval1 = this%pet_pvar(n)
      rval2 = this%extdp(n)
      call this%uzfobj%setdataet(n, ivertflag, rval1, rval2)
      !
      ! -- ETWC
      rval1 = this%extwc_pvar(n)
      call this%uzfobj%setdataetwc(n, ivertflag, rval1)
      !
      ! -- HA, HROOT, ROOTACT
      rval1 = this%ha_pvar(n)
      rval2 = this%hroot_pvar(n)
      rval3 = this%rootact_pvar(n)
      call this%uzfobj%setdataetha(n, ivertflag, rval1, rval2, rval3)
    end do
    !
    ! -- check uzfarea
    if (this%iauxmultcol > 0) then
      call this%check_cell_area()
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
  end subroutine uzf_ad

  !> @brief Formulate the HCOF and RHS terms
  !!
  !!   - skip if no UZF cells
  !!   - calculate hcof and rhs
  !<
  subroutine uzf_cf(this)
    ! -- modules
    ! -- dummy
    class(UzfType) :: this
    ! -- locals
    integer(I4B) :: n
    !
    ! -- Return if no UZF cells
    if (this%nodes == 0) return
    !
    ! -- Store values at start of outer iteration to compare with calculated
    !    values for convergence check
    do n = 1, this%maxbound
      this%rejinf0(n) = this%rejinf(n)
      this%rch0(n) = this%rch(n)
      this%gwd0(n) = this%gwd(n)
    end do
  end subroutine uzf_cf

  !> @brief Copy rhs and hcof into solution rhs and amat
  !<
  subroutine uzf_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(UzfType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: n, node, ipos
    !
    ! -- pakmvrobj fc
    if (this%imover == 1) then
      call this%pakmvrobj%fc()
    end if
    !
    ! -- Solve UZF; set reset_state to true so that waves are reset back to
    !    initial position for each outer iteration
    call this%uzf_solve(reset_state=.true.)
    !
    ! -- Copy package rhs and hcof into solution rhs and amat
    do n = 1, this%nodes
      node = this%nodelist(n)
      rhs(node) = rhs(node) + this%rhs(n)
      ipos = ia(node)
      call matrix_sln%add_value_pos(idxglo(ipos), this%hcof(n))
    end do
  end subroutine uzf_fc

  !> @brief Fill newton terms
  !<
  subroutine uzf_fn(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(UzfType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: n, node
    integer(I4B) :: ipos
    !
    ! -- Add derivative terms to rhs and amat
    do n = 1, this%nodes
      node = this%nodelist(n)
      ipos = ia(node)
      call matrix_sln%add_value_pos(idxglo(ipos), this%deriv(n))
      rhs(node) = rhs(node) + this%deriv(n) * this%xnew(node)
    end do
  end subroutine uzf_fn

  !> @brief Final convergence check for package
  !<
  subroutine uzf_cc(this, innertot, kiter, iend, icnvgmod, cpak, ipak, dpak)
    ! -- modules
    use TdisModule, only: totim, kstp, kper, delt
    ! -- dummy
    class(Uzftype), intent(inout) :: this
    integer(I4B), intent(in) :: innertot
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: icnvgmod
    integer(I4B), intent(in) :: iend
    character(len=LENPAKLOC), intent(inout) :: cpak
    integer(I4B), intent(inout) :: ipak
    real(DP), intent(inout) :: dpak
    ! -- local
    character(len=LENPAKLOC) :: cloc
    character(len=LINELENGTH) :: tag
    integer(I4B) :: icheck
    integer(I4B) :: ipakfail
    integer(I4B) :: locdrejinfmax
    integer(I4B) :: locdrchmax
    integer(I4B) :: locdseepmax
    integer(I4B) :: locdqfrommvrmax
    integer(I4B) :: ntabrows
    integer(I4B) :: ntabcols
    integer(I4B) :: n
    real(DP) :: q
    real(DP) :: q0
    real(DP) :: qtolfact
    real(DP) :: drejinf
    real(DP) :: drejinfmax
    real(DP) :: drch
    real(DP) :: drchmax
    real(DP) :: dseep
    real(DP) :: dseepmax
    real(DP) :: dqfrommvr
    real(DP) :: dqfrommvrmax
    !
    ! -- initialize local variables
    icheck = this%iconvchk
    ipakfail = 0
    locdrejinfmax = 0
    locdrchmax = 0
    locdseepmax = 0
    locdqfrommvrmax = 0
    drejinfmax = DZERO
    drchmax = DZERO
    dseepmax = DZERO
    dqfrommvrmax = DZERO
    !
    ! -- if not saving package convergence data on check convergence if
    !    the model is considered converged
    if (this%ipakcsv == 0) then
      if (icnvgmod == 0) then
        icheck = 0
      end if
    else
      !
      ! -- header for package csv
      if (.not. associated(this%pakcsvtab)) then
        !
        ! -- determine the number of columns and rows
        ntabrows = 1
        ntabcols = 9
        if (this%iseepflag == 1) then
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
        tag = 'drejinfmax'
        call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
        tag = 'drejinfmax_loc'
        call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
        tag = 'drchmax'
        call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
        tag = 'drchmax_loc'
        call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
        if (this%iseepflag == 1) then
          tag = 'dseepmax'
          call this%pakcsvtab%initialize_column(tag, 15, alignment=TABLEFT)
          tag = 'dseepmax_loc'
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
      final_check: do n = 1, this%nodes
        !
        ! -- set the Q to length factor
        qtolfact = delt / this%uzfobj%uzfarea(n)
        !
        ! -- rejected infiltration
        drejinf = qtolfact * (this%rejinf0(n) - this%rejinf(n))
        !
        ! -- groundwater recharge
        drch = qtolfact * (this%rch0(n) - this%rch(n))
        !
        ! -- groundwater seepage to the land surface
        dseep = DZERO
        if (this%iseepflag == 1) then
          dseep = qtolfact * (this%gwd0(n) - this%gwd(n))
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
          drejinfmax = drejinf
          locdrejinfmax = n
          drchmax = drch
          locdrchmax = n
          dseepmax = dseep
          locdseepmax = n
          dqfrommvrmax = dqfrommvr
          locdqfrommvrmax = n
        else
          if (ABS(drejinf) > abs(drejinfmax)) then
            drejinfmax = drejinf
            locdrejinfmax = n
          end if
          if (ABS(drch) > abs(drchmax)) then
            drchmax = drch
            locdrchmax = n
          end if
          if (ABS(dseep) > abs(dseepmax)) then
            dseepmax = dseep
            locdseepmax = n
          end if
          if (ABS(dqfrommvr) > abs(dqfrommvrmax)) then
            dqfrommvrmax = dqfrommvr
            locdqfrommvrmax = n
          end if
        end if
      end do final_check
      !
      ! -- set dpak and cpak
      if (ABS(drejinfmax) > abs(dpak)) then
        ipak = locdrejinfmax
        dpak = drejinfmax
        write (cloc, "(a,'-',a)") trim(this%packName), 'rejinf'
        cpak = trim(cloc)
      end if
      if (ABS(drchmax) > abs(dpak)) then
        ipak = locdrchmax
        dpak = drchmax
        write (cloc, "(a,'-',a)") trim(this%packName), 'rech'
        cpak = trim(cloc)
      end if
      if (this%iseepflag == 1) then
        if (ABS(dseepmax) > abs(dpak)) then
          ipak = locdseepmax
          dpak = dseepmax
          write (cloc, "(a,'-',a)") trim(this%packName), 'seep'
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
        call this%pakcsvtab%add_term(drejinfmax)
        call this%pakcsvtab%add_term(locdrejinfmax)
        call this%pakcsvtab%add_term(drchmax)
        call this%pakcsvtab%add_term(locdrchmax)
        if (this%iseepflag == 1) then
          call this%pakcsvtab%add_term(dseepmax)
          call this%pakcsvtab%add_term(locdseepmax)
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
  end subroutine uzf_cc

  !> @brief Calculate flows
  !<
  subroutine uzf_cq(this, x, flowja, iadv)
    ! -- modules
    use TdisModule, only: delt
    use ConstantsModule, only: LENBOUNDNAME, DZERO, DHNOFLO, DHDRY
    use BudgetModule, only: BudgetType
    ! -- dummy
    class(UzfType), intent(inout) :: this
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    integer(I4B), optional, intent(in) :: iadv
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: node
    real(DP) :: qout
    real(DP) :: qfact
    real(DP) :: qtomvr
    real(DP) :: q
    ! -- for observations
    ! -- formats
    character(len=*), parameter :: fmttkk = &
      &"(1X,/1X,A,'   PERIOD ',I0,'   STEP ',I0)"
    !
    ! -- Make uzf solution for budget calculations, and then reset waves.
    !    Final uzf solve will be done as part of ot().
    call this%uzf_solve(reset_state=.true.)
    !
    ! -- call base functionality in bnd_cq.  This will calculate uzf-gwf flows
    !    and put them into this%simvals and this%simvtomvr
    call this%BndType%bnd_cq(x, flowja, iadv=1)
    !
    ! -- Go through and process each UZF cell
    do n = 1, this%nodes
      !
      ! -- Initialize variables
      node = this%nodelist(n)
      !
      ! -- Skip if cell is not active
      if (this%ibound(node) < 1) cycle
      !
      ! -- infiltration terms
      this%appliedinf(n) = this%uzfobj%finf_spec(n) * this%uzfobj%uzfarea(n)
      this%infiltration(n) = this%uzfobj%surf_infil(n) * this%uzfobj%uzfarea(n)
      !
      ! -- qtomvr
      qout = this%rejinf(n) + this%uzfobj%surf_seep(n)
      qtomvr = DZERO
      if (this%imover == 1) then
        qtomvr = this%pakmvrobj%get_qtomvr(n)
      end if
      !
      ! -- rejected infiltration
      qfact = DZERO
      if (qout > DZERO) then
        qfact = this%rejinf(n) / qout
      end if
      q = this%rejinf(n)
      this%rejinftomvr(n) = qfact * qtomvr
      !
      ! -- set rejected infiltration to the remainder
      q = q - this%rejinftomvr(n)
      !
      ! -- values less than zero represent a volumetric error resulting
      !    from qtomvr being greater than water available to the mover
      if (q < DZERO) then
        q = DZERO
      end if
      this%rejinf(n) = q
      !
      ! -- calculate groundwater discharge and what goes to mover
      this%gwd(n) = this%uzfobj%surf_seep(n)
      qfact = DZERO
      if (qout > DZERO) then
        qfact = this%gwd(n) / qout
      end if
      q = this%gwd(n)
      this%gwdtomvr(n) = qfact * qtomvr
      !
      ! -- set groundwater discharge to the remainder
      q = q - this%gwdtomvr(n)
      !
      ! -- values less than zero represent a volumetric error resulting
      !    from qtomvr being greater than water available to the mover
      if (q < DZERO) then
        q = DZERO
      end if
      this%gwd(n) = q
      !
      ! -- calculate and store remaining budget terms
      this%gwet_pvar(n) = this%uzfobj%gwet(n)
      this%uzet(n) = this%uzfobj%et_uz(n) * this%uzfobj%uzfarea(n) / delt
      !
      ! -- End of UZF cell loop
      !
    end do
    !
    ! -- fill the budget object
    call this%uzf_fill_budobj()
  end subroutine uzf_cq

  !> @brief Calculate the change in mobile water stored in the unsaturated zone
  !<
  function get_storage_change(top, bot, carea, hold, hnew, wcold, wcnew, &
                              thtr, delt, iss) result(qsto)
    ! -- dummy
    real(DP), intent(in) :: top
    real(DP), intent(in) :: bot
    real(DP), intent(in) :: hold
    real(DP), intent(in) :: hnew
    real(DP), intent(in) :: wcold
    real(DP), intent(in) :: wcnew
    real(DP), intent(in) :: thtr
    real(DP), intent(in) :: carea
    real(DP), intent(in) :: delt
    integer(I4B) :: iss
    ! -- return
    real(DP) :: qsto
    ! -- local
    real(DP) :: thknew
    real(DP) :: thkold
    if (iss == 0) then
      thknew = top - max(bot, hnew)
      thkold = top - max(bot, hold)
      qsto = DZERO
      if (thknew > DZERO) then
        qsto = qsto + thknew * (wcnew - thtr)
      end if
      if (thkold > DZERO) then
        qsto = qsto - thkold * (wcold - thtr)
      end if
      qsto = qsto * carea / delt
    else
      qsto = DZERO
    end if
  end function get_storage_change

  !> @brief Add package ratin/ratout to model budget
  !<
  subroutine uzf_bd(this, model_budget)
    ! -- add package ratin/ratout to model budget
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType, rate_accumulator
    class(UzfType) :: this
    type(BudgetType), intent(inout) :: model_budget
    real(DP) :: ratin
    real(DP) :: ratout
    integer(I4B) :: isuppress_output
    isuppress_output = 0
    !
    ! -- Calculate flow from uzf to gwf (UZF-GWRCH)
    call rate_accumulator(this%rch, ratin, ratout)
    call model_budget%addentry(ratin, ratout, delt, this%bdtxt(2), &
                               isuppress_output, this%packName)
    !
    ! -- GW discharge and GW discharge to mover
    if (this%iseepflag == 1) then
      call rate_accumulator(-this%gwd, ratin, ratout)
      call model_budget%addentry(ratin, ratout, delt, this%bdtxt(3), &
                                 isuppress_output, this%packName)
      if (this%imover == 1) then
        call rate_accumulator(-this%gwdtomvr, ratin, ratout)
        call model_budget%addentry(ratin, ratout, delt, this%bdtxt(5), &
                                   isuppress_output, this%packName)
      end if
    end if
    !
    ! -- groundwater et (gwet array is positive, so switch ratin/ratout)
    if (this%igwetflag /= 0) then
      call rate_accumulator(-this%gwet_pvar, ratin, ratout)
      call model_budget%addentry(ratin, ratout, delt, this%bdtxt(4), &
                                 isuppress_output, this%packName)
    end if
  end subroutine uzf_bd

  !> @brief Write flows to binary file and/or print flows to budget
  !<
  subroutine uzf_ot_model_flows(this, icbcfl, ibudfl, icbcun, imap)
    ! -- modules
    use ConstantsModule, only: LENBOUNDNAME, DZERO
    use BndModule, only: save_print_model_flows
    ! -- dummy
    class(UzfType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    integer(I4B), intent(in) :: icbcun
    integer(I4B), dimension(:), optional, intent(in) :: imap
    ! -- local
    character(len=LINELENGTH) :: title
    integer(I4B) :: itxt
    !
    ! -- UZF-GWRCH
    itxt = 2
    title = trim(adjustl(this%bdtxt(itxt)))//' PACKAGE ('// &
            trim(this%packName)//') FLOW RATES'
    call save_print_model_flows(icbcfl, ibudfl, icbcun, this%iprflow, &
                                this%outputtab, this%nbound, this%nodelist, &
                                this%rch, this%ibound, title, this%bdtxt(itxt), &
                                this%ipakcb, this%dis, this%naux, &
                                this%name_model, this%name_model, &
                                this%name_model, this%packName, this%auxname, &
                                this%auxvar, this%iout, this%inamedbound, &
                                this%boundname)
    !
    ! -- UZF-GWD
    if (this%iseepflag == 1) then
      itxt = 3
      title = trim(adjustl(this%bdtxt(itxt)))//' PACKAGE ('// &
              trim(this%packName)//') FLOW RATES'
      call save_print_model_flows(icbcfl, ibudfl, icbcun, this%iprflow, &
                                  this%outputtab, this%nbound, this%nodelist, &
                                  -this%gwd, this%ibound, title, &
                                  this%bdtxt(itxt), this%ipakcb, this%dis, &
                                  this%naux, this%name_model, this%name_model, &
                                  this%name_model, this%packName, this%auxname, &
                                  this%auxvar, this%iout, this%inamedbound, &
                                  this%boundname)
      !
      ! -- UZF-GWD TO-MVR
      if (this%imover == 1) then
        itxt = 5
        title = trim(adjustl(this%bdtxt(itxt)))//' PACKAGE ('// &
                trim(this%packName)//') FLOW RATES'
        call save_print_model_flows(icbcfl, ibudfl, icbcun, this%iprflow, &
                                    this%outputtab, this%nbound, this%nodelist, &
                                    -this%gwdtomvr, this%ibound, title, &
                                    this%bdtxt(itxt), this%ipakcb, this%dis, &
                                    this%naux, this%name_model, this%name_model, &
                                    this%name_model, this%packName, &
                                    this%auxname, this%auxvar, this%iout, &
                                    this%inamedbound, this%boundname)
      end if
    end if
    !
    ! -- UZF-GWET
    if (this%igwetflag /= 0) then
      itxt = 4
      title = trim(adjustl(this%bdtxt(itxt)))//' PACKAGE ('// &
              trim(this%packName)//') FLOW RATES'
      call save_print_model_flows(icbcfl, ibudfl, icbcun, this%iprflow, &
                                  this%outputtab, this%nbound, this%nodelist, &
                                  -this%gwet_pvar, this%ibound, title, &
                                  this%bdtxt(itxt), this%ipakcb, this%dis, &
                                  this%naux, this%name_model, this%name_model, &
                                  this%name_model, this%packName, this%auxname, &
                                  this%auxvar, this%iout, this%inamedbound, &
                                  this%boundname)
    end if
  end subroutine uzf_ot_model_flows

  !> @brief Output UZF package flow terms
  !<
  subroutine uzf_ot_package_flows(this, icbcfl, ibudfl)
    ! -- modules
    use TdisModule, only: kstp, kper, delt, pertim, totim
    ! -- dummy
    class(UzfType) :: this
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
  end subroutine uzf_ot_package_flows

  !> @brief Save UZF-calculated values to binary file
  !<
  subroutine uzf_ot_dv(this, idvsave, idvprint)
    ! -- modules
    use TdisModule, only: kstp, kper, pertim, totim
    ! -- dummy
    use InputOutputModule, only: ulasav
    class(UzfType) :: this
    integer(I4B), intent(in) :: idvsave
    integer(I4B), intent(in) :: idvprint
    ! -- local
    integer(I4B) :: ibinun
    !
    ! -- set unit number for binary dependent variable output
    ibinun = 0
    if (this%iwcontout /= 0) then
      ibinun = this%iwcontout
    end if
    if (idvsave == 0) ibinun = 0
    !
    ! -- write uzf binary moisture-content output
    if (ibinun > 0) then
      call ulasav(this%wcnew, '   WATER-CONTENT', kstp, kper, pertim, &
                  totim, this%nodes, 1, 1, ibinun)
    end if
  end subroutine uzf_ot_dv

  !> @brief Write UZF budget to listing file
  !<
  subroutine uzf_ot_bdsummary(this, kstp, kper, iout, ibudfl)
    ! -- module
    use TdisModule, only: totim, delt
    ! -- dummy
    class(UzfType) :: this !< UzfType object
    integer(I4B), intent(in) :: kstp !< time step number
    integer(I4B), intent(in) :: kper !< period number
    integer(I4B), intent(in) :: iout !< flag and unit number for the model listing file
    integer(I4B), intent(in) :: ibudfl !< flag indicating budget should be written
    !
    call this%budobj%write_budtable(kstp, kper, iout, ibudfl, totim, delt)
  end subroutine uzf_ot_bdsummary

  !> @brief Formulate the HCOF and RHS terms
  !!
  !! Called three times per time step with different intent. uzf_fc and uzf_cq
  !! pass reset_state true, so each cell's waves are saved before it is solved
  !! and restored afterwards and the wave state does not advance. Only uzf_bd_obs
  !! passes false, and that call is what commits the waves to their end-of-step
  !! position. It runs after a converged solution, so the committed waves are the
  !! ones consistent with the final head.
  !!
  !! The cell loop cascades water downward: a cell with cell_below sets the
  !! infiltration of the cell beneath it, which is read when that cell is solved
  !! later in the same loop. This is correct only where UZF cells are numbered
  !! from the top down within a vertical stack.
  !<
  subroutine uzf_solve(this, reset_state)
    ! -- modules
    use TdisModule, only: delt
    logical, intent(in) :: reset_state !< flag indicating that waves should be reset after solution
    ! -- dummy
    class(UzfType) :: this
    ! -- locals
    integer(I4B) :: n, ivertflag
    integer(I4B) :: node, m, ierr
    real(DP) :: trhs1, thcof1, trhs2, thcof2
    real(DP) :: hgwf, uzderiv, derivgwet
    real(DP) :: qfrommvr
    real(DP) :: qformvr
    real(DP) :: wc
    real(DP) :: watabold
    !
    ! -- Initialize
    ierr = 0
    do n = 1, this%nodes
      this%uzfobj%pet(n) = this%uzfobj%pet_max(n)
    end do
    !
    ! -- Calculate hcof and rhs for each UZF entry
    do n = 1, this%nodes
      !
      ! -- Initialize hcof/rhs terms
      this%hcof(n) = DZERO
      this%rhs(n) = DZERO
      thcof1 = DZERO
      thcof2 = DZERO
      trhs1 = DZERO
      trhs2 = DZERO
      uzderiv = DZERO
      derivgwet = DZERO
      !
      ! -- Initialize variables
      node = this%nodelist(n)
      ivertflag = this%uzfobj%cell_below(n)
      watabold = this%uzfobj%water_table_old(n)
      !
      if (this%ibound(node) > 0) then
        !
        ! -- Water mover added to infiltration
        qfrommvr = DZERO
        qformvr = DZERO
        if (this%imover == 1) then
          qfrommvr = this%pakmvrobj%get_qfrommvr(n)
        end if
        !
        hgwf = this%xnew(node)
        m = node
        !
        ! -- solve for current uzf cell
        call this%uzfobj%solve(ivertflag, n, &
                               this%totfluxtot, this%ietflag, &
                               this%issflag, this%iseepflag, hgwf, &
                               qfrommvr, ierr, &
                               reset_state=reset_state, &
                               trhs=trhs1, thcof=thcof1, deriv=uzderiv, &
                               watercontent=wc)
        !
        ! -- terminate if an error condition has occurred
        if (ierr > 0) then
          if (ierr == 1) &
            errmsg = 'UZF variable NWAVESETS needs to be increased.'
          call store_error(errmsg, terminate=.TRUE.)
        end if
        !
        ! -- Calculate gwet
        if (this%igwetflag > 0) then
          call this%uzfobj%setgwpet(n)
          call this%uzfobj%simgwet(this%igwetflag, n, hgwf, trhs2, thcof2, &
                                   derivgwet)
        end if
        !
        ! -- distribute PET to deeper cells
        if (this%ietflag > 0) then
          if (this%uzfobj%cell_below(n) > 0) then
            call this%uzfobj%setbelowpet(n, ivertflag)
          end if
        end if
        !
        ! -- store derivative for Newton addition to equations in _fn()
        this%deriv(n) = uzderiv + derivgwet
        !
        ! -- save current rejected infiltration, groundwater recharge, and
        !    groundwater discharge
        this%rejinf(n) = this%uzfobj%finf_rej(n) * this%uzfobj%uzfarea(n)
        this%rch(n) = this%uzfobj%flux_to_wt(n) * this%uzfobj%uzfarea(n) / delt
        this%gwd(n) = this%uzfobj%surf_seep(n)
        !
        ! -- add to hcof and rhs
        this%hcof(n) = thcof1 + thcof2
        this%rhs(n) = -trhs1 + trhs2
        !
        ! -- add spring discharge and rejected infiltration to mover
        if (this%imover == 1) then
          qformvr = this%gwd(n) + this%rejinf(n)
          call this%pakmvrobj%accumulate_qformvr(n, qformvr)
        end if
        !
        ! -- Store water content
        this%wcnew(n) = wc
        !
        ! -- Calculate change in mobile storage
        this%qsto(n) = get_storage_change(this%uzfobj%celtop(n), &
                                          this%uzfobj%celbot(n), &
                                          this%uzfobj%uzfarea(n), &
                                          watabold, &
                                          this%uzfobj%water_table(n), &
                                          this%wcold(n), this%wcnew(n), &
                                     this%uzfobj%theta_res(n), delt, this%issflag)
        !
      end if
    end do
  end subroutine uzf_solve

  !> @brief Define the list heading that is written to iout when PRINT_INPUT
  !!  option is used
  !<
  subroutine define_listlabel(this)
    ! -- dummy
    class(UzfType), intent(inout) :: this
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

  ! -- Procedures related to observations (type-bound)

  ! -- Procedures related to observations (NOT type-bound)

  !> @brief Allocate scalar members
  !<
  subroutine uzf_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(UzfType) :: this
    !
    ! -- call standard BndType allocate scalars
    call this%BndType%allocate_scalars()
    !
    ! -- allocate uzf specific scalars
    call mem_allocate(this%iprwcont, 'IPRWCONT', this%memoryPath)
    call mem_allocate(this%iwcontout, 'IWCONTOUT', this%memoryPath)
    call mem_allocate(this%ibudgetout, 'IBUDGETOUT', this%memoryPath)
    call mem_allocate(this%ibudcsv, 'IBUDCSV', this%memoryPath)
    call mem_allocate(this%ipakcsv, 'IPAKCSV', this%memoryPath)
    call mem_allocate(this%ntrail_pvar, 'NTRAIL_PVAR', this%memoryPath)
    call mem_allocate(this%nsets, 'NSETS', this%memoryPath)
    call mem_allocate(this%nodes, 'NODES', this%memoryPath)
    call mem_allocate(this%istocb, 'ISTOCB', this%memoryPath)
    call mem_allocate(this%nwav_pvar, 'NWAV_PVAR', this%memoryPath)
    call mem_allocate(this%totfluxtot, 'TOTFLUXTOT', this%memoryPath)
    call mem_allocate(this%bditems, 'BDITEMS', this%memoryPath)
    call mem_allocate(this%nbdtxt, 'NBDTXT', this%memoryPath)
    call mem_allocate(this%issflag, 'ISSFLAG', this%memoryPath)
    call mem_allocate(this%issflagold, 'ISSFLAGOLD', this%memoryPath)
    call mem_allocate(this%readflag, 'READFLAG', this%memoryPath)
    call mem_allocate(this%iseepflag, 'ISEEPFLAG', this%memoryPath)
    call mem_allocate(this%imaxcellcnt, 'IMAXCELLCNT', this%memoryPath)
    call mem_allocate(this%ietflag, 'IETFLAG', this%memoryPath)
    call mem_allocate(this%igwetflag, 'IGWETFLAG', this%memoryPath)
    call mem_allocate(this%iuzf2uzf, 'IUZF2UZF', this%memoryPath)
    call mem_allocate(this%cbcauxitems, 'CBCAUXITEMS', this%memoryPath)
    !
    call mem_allocate(this%iconvchk, 'ICONVCHK', this%memoryPath)
    !
    ! -- initialize scalars
    this%iprwcont = 0
    this%iwcontout = 0
    this%ibudgetout = 0
    this%ibudcsv = 0
    this%ipakcsv = 0
    this%istocb = 0
    this%bditems = 7
    this%nbdtxt = 5
    this%issflag = 0
    this%issflagold = 0
    this%ietflag = 0
    this%igwetflag = 0
    this%iseepflag = 0
    this%imaxcellcnt = 0
    this%iuzf2uzf = 0
    this%cbcauxitems = 1
    this%imover = 0
    !
    ! -- convergence check
    this%iconvchk = 1
  end subroutine uzf_allocate_scalars

  !> @brief Deallocate objects
  !<
  subroutine uzf_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(UzfType) :: this
    !
    ! -- deallocate uzf objects
    call this%uzfobj%dealloc()
    deallocate (this%uzfobj)
    nullify (this%uzfobj)
    !
    call this%budobj%budgetobject_da()
    deallocate (this%budobj)
    nullify (this%budobj)
    !
    ! -- character arrays
    deallocate (this%bdtxt)
    deallocate (this%cauxcbc)
    deallocate (this%uzfname)
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
    ! -- deallocate scalars
    call mem_deallocate(this%iprwcont)
    call mem_deallocate(this%iwcontout)
    call mem_deallocate(this%ibudgetout)
    call mem_deallocate(this%ibudcsv)
    call mem_deallocate(this%ipakcsv)
    call mem_deallocate(this%ntrail_pvar)
    call mem_deallocate(this%nsets)
    call mem_deallocate(this%nodes)
    call mem_deallocate(this%istocb)
    call mem_deallocate(this%nwav_pvar)
    call mem_deallocate(this%totfluxtot)
    call mem_deallocate(this%bditems)
    call mem_deallocate(this%nbdtxt)
    call mem_deallocate(this%issflag)
    call mem_deallocate(this%issflagold)
    call mem_deallocate(this%readflag)
    call mem_deallocate(this%iseepflag)
    call mem_deallocate(this%imaxcellcnt)
    call mem_deallocate(this%ietflag)
    call mem_deallocate(this%igwetflag)
    call mem_deallocate(this%iuzf2uzf)
    call mem_deallocate(this%cbcauxitems)
    !
    ! -- convergence check
    call mem_deallocate(this%iconvchk)
    !
    ! -- deallocate arrays
    call mem_deallocate(this%igwfnode)
    call mem_deallocate(this%appliedinf)
    call mem_deallocate(this%rejinf)
    call mem_deallocate(this%rejinf0)
    call mem_deallocate(this%rejinftomvr)
    call mem_deallocate(this%infiltration)
    call mem_deallocate(this%gwet_pvar)
    call mem_deallocate(this%uzet)
    call mem_deallocate(this%gwd)
    call mem_deallocate(this%gwd0)
    call mem_deallocate(this%gwdtomvr)
    call mem_deallocate(this%rch)
    call mem_deallocate(this%rch0)
    call mem_deallocate(this%qsto)
    call mem_deallocate(this%deriv)
    call mem_deallocate(this%qauxcbc)
    call mem_deallocate(this%wcnew)
    call mem_deallocate(this%wcold)
    !
    ! -- deallocate integer arrays
    call mem_deallocate(this%ia)
    call mem_deallocate(this%ja)
    !
    ! -- deallocate timeseries aware variables
    call mem_deallocate(this%sinf_pvar)
    call mem_deallocate(this%pet_pvar)
    call mem_deallocate(this%extdp)
    call mem_deallocate(this%extwc_pvar)
    call mem_deallocate(this%ha_pvar)
    call mem_deallocate(this%hroot_pvar)
    call mem_deallocate(this%rootact_pvar)
    call mem_deallocate(this%uauxvar)
    !
    ! -- Parent object
    call this%BndType%bnd_da()
  end subroutine uzf_da

end module UzfModule
