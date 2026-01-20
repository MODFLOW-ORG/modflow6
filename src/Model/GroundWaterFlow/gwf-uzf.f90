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
    type(UzfCellGroupType) :: uzfobjwork
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
    procedure, private :: findcellabove
    procedure, private :: check_cell_area
    !
    ! -- budget
    procedure, private :: uzf_setup_budobj
    procedure, private :: uzf_fill_budobj

  end type UzfType

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
    integer(I4B) :: n, i
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
    do i = 1, this%nodes
      this%nodelist(i) = this%igwfnode(i)
      n = this%igwfnode(i)
      hgwf = this%xnew(n)
      call this%uzfobj%sethead(i, hgwf)
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
    integer(I4B) :: i
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
    do i = 1, this%nodes
      this%appliedinf(i) = DZERO
      this%rejinf(i) = DZERO
      this%rejinf0(i) = DZERO
      this%rejinftomvr(i) = DZERO
      this%gwet_pvar(i) = DZERO
      this%uzet(i) = DZERO
      this%gwd(i) = DZERO
      this%gwd0(i) = DZERO
      this%gwdtomvr(i) = DZERO
      this%rch(i) = DZERO
      this%rch0(i) = DZERO
      this%qsto(i) = DZERO
      this%deriv(i) = DZERO
      ! -- integer variables
      this%ja(i) = 0
      ! -- timeseries aware variables
      this%sinf_pvar(i) = DZERO
      this%pet_pvar(i) = DZERO
      this%extdp(i) = DZERO
      this%extwc_pvar(i) = DZERO
      this%ha_pvar(i) = DZERO
      this%hroot_pvar(i) = DZERO
      this%rootact_pvar(i) = DZERO
      do j = 1, this%naux
        if (this%iauxmultcol > 0 .and. j == this%iauxmultcol) then
          this%uauxvar(j, i) = DONE
        else
          this%uauxvar(j, i) = DZERO
        end if
      end do
    end do
    do i = 1, this%dis%nodes + 1
      this%ia(i) = 0
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
    do i = 1, this%nodes
      this%wcnew(i) = DZERO
      this%wcold(i) = DZERO
    end do
    !
    ! -- allocate character array for aux budget text
    allocate (this%cauxcbc(this%cbcauxitems))
    allocate (this%uzfname(this%nodes))
    !
    ! -- allocate and initialize qauxcbc
    call mem_allocate(this%qauxcbc, this%cbcauxitems, 'QAUXCBC', this%memoryPath)
    do i = 1, this%cbcauxitems
      this%qauxcbc(i) = DZERO
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

    if (this%ntrail_pvar <= 0) then
      write (errmsg, '(a)') &
        'NTRAILWAVES was not specified or was specified incorrectly.'
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
    call this%uzfobjwork%init(1, this%nwav_pvar)
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
    integer(I4B) :: i
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
        i = this%parser%GetInteger()
        if (i < 1 .or. i > this%nodes) then
          tag = trim(adjustl(this%text))//' PACKAGE ('// &
                trim(adjustl(this%packName))//') DATA FOR PERIOD'
          write (tag, '(a,1x,i0)') trim(adjustl(tag)), kper

          write (errmsg, '(a,a,i0,1x,a,i0,a)') &
            trim(adjustl(tag)), ': UZFNO ', i, &
            'must be greater than 0 and less than or equal to ', this%nodes, '.'
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- Setup boundname
        if (this%inamedbound > 0) then
          bndName = this%boundname(i)
        else
          bndName = ''
        end if
        !
        ! -- FINF
        call this%parser%GetStringCaps(text)
        jj = 1 ! For SINF
        bndElem => this%sinf_pvar(i)
        call read_value_or_time_series_adv(text, i, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'SINF')
        !
        ! -- PET
        call this%parser%GetStringCaps(text)
        jj = 1 ! For PET
        bndElem => this%pet_pvar(i)
        call read_value_or_time_series_adv(text, i, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'PET')
        !
        ! -- EXTD
        call this%parser%GetStringCaps(text)
        jj = 1 ! For EXTDP
        bndElem => this%extdp(i)
        call read_value_or_time_series_adv(text, i, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'EXTDP')
        !
        ! -- EXTWC
        call this%parser%GetStringCaps(text)
        jj = 1 ! For EXTWC
        bndElem => this%extwc_pvar(i)
        call read_value_or_time_series_adv(text, i, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'EXTWC')
        !
        ! -- HA
        call this%parser%GetStringCaps(text)
        jj = 1 ! For HA
        bndElem => this%ha_pvar(i)
        call read_value_or_time_series_adv(text, i, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'HA')
        !
        ! -- HROOT
        call this%parser%GetStringCaps(text)
        jj = 1 ! For HROOT
        bndElem => this%hroot_pvar(i)
        call read_value_or_time_series_adv(text, i, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'HROOT')
        !
        ! -- ROOTACT
        call this%parser%GetStringCaps(text)
        jj = 1 ! For ROOTACT
        bndElem => this%rootact_pvar(i)
        call read_value_or_time_series_adv(text, i, jj, bndElem, this%packName, &
                                           'BND', this%tsManager, this%iprpak, &
                                           'ROOTACT')
        !
        ! -- read auxiliary variables
        do j = 1, this%naux
          call this%parser%GetStringCaps(text)
          bndElem => this%uauxvar(j, i)
          call read_value_or_time_series_adv(text, i, j, bndElem, this%packName, &
                                             'AUX', this%tsManager, this%iprpak, &
                                             this%auxname(j))
        end do
        !
        ! -- write line
        if (this%iprpak /= 0) then
          !
          ! -- get cellid
          node = this%igwfnode(i)
          if (node > 0) then
            call this%dis%noder_to_string(node, cellid)
          else
            cellid = 'none'
          end if
          !
          ! -- write data to the table
          call this%inputtab%add_term(i)
          call this%inputtab%add_term(cellid)
          call this%inputtab%add_term(this%sinf_pvar(i))
          if (this%ietflag /= 0) then
            call this%inputtab%add_term(this%pet_pvar(i))
            call this%inputtab%add_term(this%extdp(i))
            call this%inputtab%add_term(this%extwc_pvar(i))
            if (this%ietflag == 2) then
              call this%inputtab%add_term(this%ha_pvar(i))
              call this%inputtab%add_term(this%hroot_pvar(i))
              call this%inputtab%add_term(this%rootact_pvar(i))
            end if
          end if
          if (this%inamedbound == 1) then
            call this%inputtab%add_term(this%boundname(i))
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
      do i = 1, this%nodes
        call this%uzfobj%setwaves(i)
      end do
    end if
    !
    ! -- Initialize the water content
    if (kper == 1) then
      do i = 1, this%nodes
        this%wcnew(i) = this%uzfobj%get_wcnew(i)
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
    integer(I4B) :: i
    integer(I4B) :: ivertflag
    integer(I4B) :: n, iaux
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
      do i = 1, this%nodes
        this%wcold(i) = this%wcnew(i)
      end do
    else
      !
      ! -- Copy wcold back into wcnew as this is a retry of this time step.
      !    Note that there is no need to reset the waves as they are not
      !    advanced to their new state until the _ot() method is called,
      !    and that doesn't happen until a successful solution is obtained.
      do i = 1, this%nodes
        this%wcnew(i) = this%wcold(i)
      end do
    end if
    !
    ! -- advance each uzf obj
    do i = 1, this%nodes
      call this%uzfobj%advance(i)
    end do
    !
    ! -- update uzf objects with timeseries aware variables
    do i = 1, this%nodes
      !
      ! -- Set ivertflag
      ivertflag = this%uzfobj%ivertcon(i)
      !
      ! -- recalculate uzfarea
      if (this%iauxmultcol > 0) then
        rval1 = this%uauxvar(this%iauxmultcol, i)
        call this%uzfobj%setdatauzfarea(i, rval1)
      end if
      !
      ! -- FINF
      rval1 = this%sinf_pvar(i)
      call this%uzfobj%setdatafinf(i, rval1)
      !
      ! -- PET, EXTDP
      rval1 = this%pet_pvar(i)
      rval2 = this%extdp(i)
      call this%uzfobj%setdataet(i, ivertflag, rval1, rval2)
      !
      ! -- ETWC
      rval1 = this%extwc_pvar(i)
      call this%uzfobj%setdataetwc(i, ivertflag, rval1)
      !
      ! -- HA, HROOT, ROOTACT
      rval1 = this%ha_pvar(i)
      rval2 = this%hroot_pvar(i)
      rval3 = this%rootact_pvar(i)
      call this%uzfobj%setdataetha(i, ivertflag, rval1, rval2, rval3)
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
    integer(I4B) :: i, n, ipos
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
    do i = 1, this%nodes
      n = this%nodelist(i)
      rhs(n) = rhs(n) + this%rhs(i)
      ipos = ia(n)
      call matrix_sln%add_value_pos(idxglo(ipos), this%hcof(i))
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
    integer(I4B) :: i, n
    integer(I4B) :: ipos
    !
    ! -- Add derivative terms to rhs and amat
    do i = 1, this%nodes
      n = this%nodelist(i)
      ipos = ia(n)
      call matrix_sln%add_value_pos(idxglo(ipos), this%deriv(i))
      rhs(n) = rhs(n) + this%deriv(i) * this%xnew(n)
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
    integer(I4B) :: i
    integer(I4B) :: n
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
    do i = 1, this%nodes
      !
      ! -- Initialize variables
      n = this%nodelist(i)
      !
      ! -- Skip if cell is not active
      if (this%ibound(n) < 1) cycle
      !
      ! -- infiltration terms
      this%appliedinf(i) = this%uzfobj%sinf(i) * this%uzfobj%uzfarea(i)
      this%infiltration(i) = this%uzfobj%surflux(i) * this%uzfobj%uzfarea(i)
      !
      ! -- qtomvr
      qout = this%rejinf(i) + this%uzfobj%surfseep(i)
      qtomvr = DZERO
      if (this%imover == 1) then
        qtomvr = this%pakmvrobj%get_qtomvr(i)
      end if
      !
      ! -- rejected infiltration
      qfact = DZERO
      if (qout > DZERO) then
        qfact = this%rejinf(i) / qout
      end if
      q = this%rejinf(i)
      this%rejinftomvr(i) = qfact * qtomvr
      !
      ! -- set rejected infiltration to the remainder
      q = q - this%rejinftomvr(i)
      !
      ! -- values less than zero represent a volumetric error resulting
      !    from qtomvr being greater than water available to the mover
      if (q < DZERO) then
        q = DZERO
      end if
      this%rejinf(i) = q
      !
      ! -- calculate groundwater discharge and what goes to mover
      this%gwd(i) = this%uzfobj%surfseep(i)
      qfact = DZERO
      if (qout > DZERO) then
        qfact = this%gwd(i) / qout
      end if
      q = this%gwd(i)
      this%gwdtomvr(i) = qfact * qtomvr
      !
      ! -- set groundwater discharge to the remainder
      q = q - this%gwdtomvr(i)
      !
      ! -- values less than zero represent a volumetric error resulting
      !    from qtomvr being greater than water available to the mover
      if (q < DZERO) then
        q = DZERO
      end if
      this%gwd(i) = q
      !
      ! -- calculate and store remaining budget terms
      this%gwet_pvar(i) = this%uzfobj%gwet(i)
      this%uzet(i) = this%uzfobj%etact(i) * this%uzfobj%uzfarea(i) / delt
      !
      ! -- End of UZF cell loop
      !
    end do
    !
    ! -- fill the budget object
    call this%uzf_fill_budobj()
  end subroutine uzf_cq

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
  !<
  subroutine uzf_solve(this, reset_state)
    ! -- modules
    use TdisModule, only: delt
    logical, intent(in) :: reset_state !< flag indicating that waves should be reset after solution
    ! -- dummy
    class(UzfType) :: this
    ! -- locals
    integer(I4B) :: i, ivertflag
    integer(I4B) :: n, m, ierr
    real(DP) :: trhs1, thcof1, trhs2, thcof2
    real(DP) :: hgwf, uzderiv, derivgwet
    real(DP) :: qfrommvr
    real(DP) :: qformvr
    real(DP) :: wc
    real(DP) :: watabold
    !
    ! -- Initialize
    ierr = 0
    do i = 1, this%nodes
      this%uzfobj%pet(i) = this%uzfobj%petmax(i)
    end do
    !
    ! -- Calculate hcof and rhs for each UZF entry
    do i = 1, this%nodes
      !
      ! -- Initialize hcof/rhs terms
      this%hcof(i) = DZERO
      this%rhs(i) = DZERO
      thcof1 = DZERO
      thcof2 = DZERO
      trhs1 = DZERO
      trhs2 = DZERO
      uzderiv = DZERO
      derivgwet = DZERO
      !
      ! -- Initialize variables
      n = this%nodelist(i)
      ivertflag = this%uzfobj%ivertcon(i)
      watabold = this%uzfobj%watabold(i)
      !
      if (this%ibound(n) > 0) then
        !
        ! -- Water mover added to infiltration
        qfrommvr = DZERO
        qformvr = DZERO
        if (this%imover == 1) then
          qfrommvr = this%pakmvrobj%get_qfrommvr(i)
        end if
        !
        hgwf = this%xnew(n)
        m = n
        !
        ! -- solve for current uzf cell
        call this%uzfobj%solve(this%uzfobjwork, ivertflag, i, &
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
          call this%uzfobj%setgwpet(i)
          call this%uzfobj%simgwet(this%igwetflag, i, hgwf, trhs2, thcof2, &
                                   derivgwet)
        end if
        !
        ! -- distribute PET to deeper cells
        if (this%ietflag > 0) then
          if (this%uzfobj%ivertcon(i) > 0) then
            call this%uzfobj%setbelowpet(i, ivertflag)
          end if
        end if
        !
        ! -- store derivative for Newton addition to equations in _fn()
        this%deriv(i) = uzderiv + derivgwet
        !
        ! -- save current rejected infiltration, groundwater recharge, and
        !    groundwater discharge
        this%rejinf(i) = this%uzfobj%finf_rej(i) * this%uzfobj%uzfarea(i)
        this%rch(i) = this%uzfobj%totflux(i) * this%uzfobj%uzfarea(i) / delt
        this%gwd(i) = this%uzfobj%surfseep(i)
        !
        ! -- add to hcof and rhs
        this%hcof(i) = thcof1 + thcof2
        this%rhs(i) = -trhs1 + trhs2
        !
        ! -- add spring discharge and rejected infiltration to mover
        if (this%imover == 1) then
          qformvr = this%gwd(i) + this%rejinf(i)
          call this%pakmvrobj%accumulate_qformvr(i, qformvr)
        end if
        !
        ! -- Store water content
        this%wcnew(i) = wc
        !
        ! -- Calculate change in mobile storage
        this%qsto(i) = get_storage_change(this%uzfobj%celtop(i), &
                                          this%uzfobj%celbot(i), &
                                          this%uzfobj%uzfarea(i), &
                                          watabold, &
                                          this%uzfobj%watab(i), &
                                          this%wcold(i), this%wcnew(i), &
                                          this%uzfobj%thtr(i), delt, this%issflag)
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

  !> @brief Identify overlying cell ID based on user-specified mapping
  !<
  subroutine findcellabove(this, n, nml)
    ! -- dummy
    class(UzfType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(inout) :: nml
    ! -- local
    integer(I4B) :: m, ipos
    !
    ! -- Return nml = n if no cell is above it
    nml = n
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)
      if (this%dis%con%ihc(ipos) /= 0) then
        if (n < m) then
          ! -- m is beneath n
        else
          nml = m ! -- m is above n
          exit
        end if
      end if
    end do
  end subroutine findcellabove

  !> @brief Read UZF cell properties and set them for UzfCellGroup type
  !<
  subroutine read_cell_properties(this)
    ! -- modules
    use InputOutputModule, only: urword
    use SimModule, only: store_error, count_errors
    ! -- dummy
    class(UzfType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: cellid
    integer(I4B) :: ierr
    integer(I4B) :: i, n
    integer(I4B) :: j
    integer(I4B) :: ic
    integer(I4B) :: jcol
    logical :: isfound, endOfBlock
    integer(I4B) :: landflag
    integer(I4B) :: ivertcon
    real(DP) :: surfdep, vks, thtr, thts, thti, eps, hgwf
    integer(I4B), dimension(:), allocatable :: rowmaxnnz
    type(sparsematrix) :: sparse
    integer(I4B), dimension(:), allocatable :: nboundchk
    !
    ! -- allocate space for node counter and initialize
    allocate (rowmaxnnz(this%dis%nodes))
    do n = 1, this%dis%nodes
      rowmaxnnz(n) = 0
    end do
    !
    ! -- allocate space for local variables
    allocate (nboundchk(this%nodes))
    do n = 1, this%nodes
      nboundchk(n) = 0
    end do
    !
    ! -- initialize variables
    landflag = 0
    ivertcon = 0
    surfdep = DZERO
    vks = DZERO
    thtr = DZERO
    thts = DZERO
    thti = DZERO
    eps = DZERO
    hgwf = DZERO
    !
    ! -- get uzf properties block
    call this%parser%GetBlock('PACKAGEDATA', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse locations block if detected
    if (isfound) then
      write (this%iout, '(/1x,3a)') 'PROCESSING ', trim(adjustl(this%text)), &
        ' PACKAGEDATA'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        !
        ! -- get uzf cell number
        i = this%parser%GetInteger()

        if (i < 1 .or. i > this%nodes) then
          write (errmsg, '(2(a,1x),i0,a)') &
            'IUZNO must be greater than 0 and less than', &
            'or equal to', this%nodes, '.'
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- increment nboundchk
        nboundchk(i) = nboundchk(i) + 1
        !
        ! -- store the reduced gwf nodenumber in igwfnode
        call this%parser%GetCellid(this%dis%ndim, cellid)
        ic = this%dis%noder_from_cellid(cellid, &
                                        this%parser%iuactive, this%iout)
        this%igwfnode(i) = ic
        rowmaxnnz(ic) = rowmaxnnz(ic) + 1
        !
        ! -- landflag
        landflag = this%parser%GetInteger()
        if (landflag < 0 .OR. landflag > 1) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0,a)') &
            'LANDFLAG for uzf cell', i, &
            'must be 0 or 1 (specified value is', landflag, ').'
          call store_error(errmsg)
        end if
        !
        ! -- ivertcon
        ivertcon = this%parser%GetInteger()
        if (ivertcon < 0 .OR. ivertcon > this%nodes) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0,a)') &
            'IVERTCON for uzf cell', i, &
            'must be 0 or less than NUZFCELLS (specified value is', &
            ivertcon, ').'
          call store_error(errmsg)
        end if
        !
        ! -- surfdep
        surfdep = this%parser%GetDouble()
        if (surfdep <= DZERO .and. landflag > 0) then !need to check for cell thickness
          write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
            'SURFDEP for uzf cell', i, &
            'must be greater than 0 (specified value is', surfdep, ').'
          call store_error(errmsg)
        end if
        if (surfdep >= this%dis%top(ic) - this%dis%bot(ic)) then
          write (errmsg, '(a,1x,i0,1x,a)') &
            'SURFDEP for uzf cell', i, &
            'cannot be greater than the cell thickness.'
          call store_error(errmsg)
        end if
        !
        ! -- vks
        vks = this%parser%GetDouble()
        if (vks <= DZERO) then
          write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
            'VKS for uzf cell', i, &
            'must be greater than 0 (specified value ia', vks, ').'
          call store_error(errmsg)
        end if
        !
        ! -- thtr
        thtr = this%parser%GetDouble()
        if (thtr <= DZERO) then
          write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
            'THTR for uzf cell', i, &
            'must be greater than 0 (specified value is', thtr, ').'
          call store_error(errmsg)
        end if
        !
        ! -- thts
        thts = this%parser%GetDouble()
        if (thts <= thtr) then
          write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
            'THTS for uzf cell', i, &
            'must be greater than THTR (specified value is', thts, ').'
          call store_error(errmsg)
        end if
        !
        ! -- thti
        thti = this%parser%GetDouble()
        if (thti < thtr .OR. thti > thts) then
          write (errmsg, '(a,1x,i0,1x,a,1x,a,1x,g0,a)') &
            'THTI for uzf cell', i, &
            'must be greater than or equal to THTR AND less than THTS', &
            '(specified value is', thti, ').'
          call store_error(errmsg)
        end if
        !
        ! -- eps
        eps = this%parser%GetDouble()
        if (eps < 3.5 .OR. eps > 14) then
          write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
            'EPSILON for uzf cell', i, &
            'must be between 3.5 and 14.0 (specified value is', eps, ').'
          call store_error(errmsg)
        end if
        !
        ! -- boundname
        if (this%inamedbound == 1) then
          call this%parser%GetStringCaps(this%uzfname(i))
        end if
        !
        ! -- set data if there are no data errors
        if (count_errors() == 0) then
          n = this%igwfnode(i)
          call this%uzfobj%setdata(i, this%dis%area(n), this%dis%top(n), &
                                   this%dis%bot(n), surfdep, vks, thtr, thts, &
                                   thti, eps, this%ntrail_pvar, landflag, &
                                   ivertcon)
          if (ivertcon > 0) then
            this%iuzf2uzf = 1
          end if
        end if
        !
      end do
      write (this%iout, '(1x,3a)') &
        'END OF ', trim(adjustl(this%text)), ' PACKAGEDATA'
    else
      call store_error('Required packagedata block not found.')
    end if
    !
    ! -- check for duplicate or missing uzf cells
    do i = 1, this%nodes
      if (nboundchk(i) == 0) then
        write (errmsg, '(a,1x,i0,a)') &
          'No data specified for uzf cell', i, '.'
        call store_error(errmsg)
      else if (nboundchk(i) > 1) then
        write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
          'Data for uzf cell', i, 'specified', nboundchk(i), 'times.'
        call store_error(errmsg)
      end if
    end do
    !
    ! -- write summary of UZF cell property error messages
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- setup sparse for connectivity used to identify multiple uzf cells per
    !    GWF model cell
    call sparse%init(this%dis%nodes, this%dis%nodes, rowmaxnnz)
    ! --
    do i = 1, this%nodes
      ic = this%igwfnode(i)
      call sparse%addconnection(ic, i, 1)
    end do
    !
    ! -- create ia and ja from sparse
    call sparse%filliaja(this%ia, this%ja, ierr)
    !
    ! -- set imaxcellcnt
    do i = 1, this%dis%nodes
      jcol = 0
      do j = this%ia(i), this%ia(i + 1) - 1
        jcol = jcol + 1
      end do
      if (jcol > this%imaxcellcnt) then
        this%imaxcellcnt = jcol
      end if
    end do
    !
    ! -- do an initial evaluation of the sum of uzfarea relative to the
    !    GWF cell area in the case that there is more than one UZF object
    !    in a GWF cell and a auxmult value is not being applied to the
    !    calculate the UZF cell area from the GWF cell area.
    if (this%imaxcellcnt > 1 .and. this%iauxmultcol < 1) then
      call this%check_cell_area()
    end if
    !
    ! -- deallocate local variables
    deallocate (rowmaxnnz)
    deallocate (nboundchk)
  end subroutine read_cell_properties

  !> @brief Read UZF cell properties and set them for UZFCellGroup type
  !<
  subroutine print_cell_properties(this)
    ! -- dummy
    class(UzfType), intent(inout) :: this
    ! -- local
    character(len=20) :: cellid
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: tag
    integer(I4B) :: ntabrows
    integer(I4B) :: ntabcols
    integer(I4B) :: i
    integer(I4B) :: node
    !
    ! -- setup inputtab tableobj
    !
    ! -- table dimensions
    ntabrows = this%nodes
    ntabcols = 10
    if (this%inamedbound == 1) then
      ntabcols = ntabcols + 1
    end if
    !
    ! -- initialize table and define columns
    title = trim(adjustl(this%text))//' PACKAGE ('// &
            trim(adjustl(this%packName))//') STATIC UZF CELL DATA'
    call table_cr(this%inputtab, this%packName, title)
    call this%inputtab%table_df(ntabrows, ntabcols, this%iout)
    tag = 'NUMBER'
    call this%inputtab%initialize_column(tag, 10)
    tag = 'CELLID'
    call this%inputtab%initialize_column(tag, 20, alignment=TABLEFT)
    tag = 'LANDFLAG'
    call this%inputtab%initialize_column(tag, 12)
    tag = 'IVERTCON'
    call this%inputtab%initialize_column(tag, 12)
    tag = 'SURFDEP'
    call this%inputtab%initialize_column(tag, 12)
    tag = 'VKS'
    call this%inputtab%initialize_column(tag, 12)
    tag = 'THTR'
    call this%inputtab%initialize_column(tag, 12)
    tag = 'THTS'
    call this%inputtab%initialize_column(tag, 12)
    tag = 'THTI'
    call this%inputtab%initialize_column(tag, 12)
    tag = 'EPS'
    call this%inputtab%initialize_column(tag, 12)
    if (this%inamedbound == 1) then
      tag = 'BOUNDNAME'
      call this%inputtab%initialize_column(tag, LENBOUNDNAME, alignment=TABLEFT)
    end if
    !
    ! -- write data for each cell
    do i = 1, this%nodes
      !
      ! -- get cellid
      node = this%igwfnode(i)
      if (node > 0) then
        call this%dis%noder_to_string(node, cellid)
      else
        cellid = 'none'
      end if
      !
      ! -- add data
      call this%inputtab%add_term(i)
      call this%inputtab%add_term(cellid)
      call this%inputtab%add_term(this%uzfobj%landflag(i))
      call this%inputtab%add_term(this%uzfobj%ivertcon(i))
      call this%inputtab%add_term(this%uzfobj%surfdep(i))
      call this%inputtab%add_term(this%uzfobj%vks(i))
      call this%inputtab%add_term(this%uzfobj%thtr(i))
      call this%inputtab%add_term(this%uzfobj%thts(i))
      call this%inputtab%add_term(this%uzfobj%thti(i))
      call this%inputtab%add_term(this%uzfobj%eps(i))
      if (this%inamedbound == 1) then
        call this%inputtab%add_term(this%uzfname(i))
      end if
    end do
  end subroutine print_cell_properties

  !> @brief Check UZF cell areas
  !<
  subroutine check_cell_area(this)
    ! -- modules
    use InputOutputModule, only: urword
    use SimModule, only: store_error, count_errors
    ! -- dummy
    class(UzfType) :: this
    ! -- local
    character(len=16) :: cuzf
    character(len=20) :: cellid
    character(len=LINELENGTH) :: cuzfcells
    integer(I4B) :: i
    integer(I4B) :: i2
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: i0
    integer(I4B) :: i1
    real(DP) :: area
    real(DP) :: area2
    real(DP) :: sumarea
    real(DP) :: cellarea
    real(DP) :: d
    !
    ! -- check that the area of vertically connected uzf cells is the equal
    do i = 1, this%nodes
      !
      ! -- Initialize variables
      i2 = this%uzfobj%ivertcon(i)
      area = this%uzfobj%uzfarea(i)
      !
      ! Create pointer to object below
      if (i2 > 0) then
        area2 = this%uzfobj%uzfarea(i2)
        d = abs(area - area2)
        if (d > DEM6) then
          write (errmsg, '(2(a,1x,g0,1x,a,1x,i0,1x),a)') &
            'UZF cell area (', area, ') for cell ', i, &
            'does not equal uzf cell area (', area2, ') for cell ', i2, '.'
          call store_error(errmsg)
        end if
      end if
    end do
    !
    ! -- check that the area of uzf cells in a GWF cell is less than or equal
    !    to the GWF cell area
    do n = 1, this%dis%nodes
      i0 = this%ia(n)
      i1 = this%ia(n + 1)
      ! -- skip gwf cells with no UZF cells
      if ((i1 - i0) < 1) cycle
      sumarea = DZERO
      cellarea = DZERO
      cuzfcells = ''
      do j = i0, i1 - 1
        i = this%ja(j)
        write (cuzf, '(i0)') i
        cuzfcells = trim(adjustl(cuzfcells))//' '//trim(adjustl(cuzf))
        sumarea = sumarea + this%uzfobj%uzfarea(i)
        cellarea = this%uzfobj%cellarea(i)
      end do
      ! -- calculate the difference between the sum of UZF areas and GWF cell area
      d = sumarea - cellarea
      if (d > DEM6) then
        call this%dis%noder_to_string(n, cellid)
        write (errmsg, '(a,1x,g0,1x,a,1x,g0,1x,a,1x,a,1x,a,a,a)') &
          'Total uzf cell area (', sumarea, &
          ') exceeds the gwf cell area (', cellarea, ') of cell', cellid, &
          'which includes uzf cell(s): ', trim(adjustl(cuzfcells)), '.'
        call store_error(errmsg)
      end if
    end do
    !
    ! -- terminate if errors were encountered
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine check_cell_area

  ! -- Procedures related to observations (type-bound)

  !> @brief Return true because uzf package supports observations
  !!
  !! Overrides BndType%bnd_obs_supported
  !<
  logical function uzf_obs_supported(this)
    ! -- dummy
    class(UzfType) :: this
    !
    uzf_obs_supported = .true.
  end function uzf_obs_supported

  !> @brief Implements bnd_df_obs
  !!
  !! Store observation type supported by uzf package.
  !! Overrides BndType%bnd_df_obs
  !<
  subroutine uzf_df_obs(this)
    ! -- dummy
    class(UzfType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    ! -- Store obs type and assign procedure pointer
    !
    !    for recharge observation type.
    call this%obs%StoreObsType('uzf-gwrch', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
    !
    !    for discharge observation type.
    call this%obs%StoreObsType('uzf-gwd', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
    !
    !    for discharge observation type.
    call this%obs%StoreObsType('uzf-gwd-to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
    !
    !    for gwet observation type.
    call this%obs%StoreObsType('uzf-gwet', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
    !
    !    for infiltration observation type.
    call this%obs%StoreObsType('infiltration', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
    !
    !    for from mover observation type.
    call this%obs%StoreObsType('from-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
    !
    !    for rejected infiltration observation type.
    call this%obs%StoreObsType('rej-inf', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
    !
    !    for rejected infiltration to mover observation type.
    call this%obs%StoreObsType('rej-inf-to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
    !
    !    for uzet observation type.
    call this%obs%StoreObsType('uzet', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
    !
    !    for storage observation type.
    call this%obs%StoreObsType('storage', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
    !
    !    for net infiltration observation type.
    call this%obs%StoreObsType('net-infiltration', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
    !
    !    for water-content observation type.
    call this%obs%StoreObsType('water-content', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
  end subroutine uzf_df_obs

  !> @brief Calculate observations this time step and call ObsType%SaveOneSimval
  !! for each UzfType observation
  !<
  subroutine uzf_bd_obs(this)
    ! -- dummy
    class(UzfType) :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: ii
    integer(I4B) :: n
    real(DP) :: v
    type(ObserveType), pointer :: obsrv => null()
    !
    ! -- Make final uzf solution, and do not reset waves.  This will advance
    !    the waves to their new state at the end of the time step.  This should
    !    be the first step of the uzf ot() routines.
    call this%uzf_solve(reset_state=.false.)
    !
    ! Write simulated values for all uzf observations
    if (this%obs%npakobs > 0) then
      call this%obs%obs_bd_clear()
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv
        do ii = 1, obsrv%indxbnds_count
          n = obsrv%indxbnds(ii)
          v = DNODATA
          select case (obsrv%ObsTypeId)
          case ('UZF-GWRCH')
            v = this%rch(n)
          case ('UZF-GWD')
            v = this%gwd(n)
            if (v > DZERO) then
              v = -v
            end if
          case ('UZF-GWD-TO-MVR')
            if (this%imover == 1) then
              v = this%gwdtomvr(n)
              if (v > DZERO) then
                v = -v
              end if
            end if
          case ('UZF-GWET')
            if (this%igwetflag > 0) then
              v = this%gwet_pvar(n)
              if (v > DZERO) then
                v = -v
              end if
            end if
          case ('INFILTRATION')
            v = this%appliedinf(n)
          case ('FROM-MVR')
            if (this%imover == 1) then
              v = this%pakmvrobj%get_qfrommvr(n)
            end if
          case ('REJ-INF')
            v = this%rejinf(n)
            if (v > DZERO) then
              v = -v
            end if
          case ('REJ-INF-TO-MVR')
            if (this%imover == 1) then
              v = this%rejinftomvr(n)
              if (v > DZERO) then
                v = -v
              end if
            end if
          case ('UZET')
            if (this%ietflag /= 0) then
              v = this%uzet(n)
              if (v > DZERO) then
                v = -v
              end if
            end if
          case ('STORAGE')
            v = -this%qsto(n)
          case ('NET-INFILTRATION')
            v = this%infiltration(n)
          case ('WATER-CONTENT')
            v = this%uzfobj%get_water_content_at_depth(n, obsrv%obsDepth)
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
        call this%parser%StoreErrorUnit()
      end if
    end if
  end subroutine uzf_bd_obs

  !> @brief Process each observation
  !!
  !! Only done the first stress period since boundaries are fixed for the
  !! simulation
  !<
  subroutine uzf_rp_obs(this)
    ! -- modules
    use TdisModule, only: kper
    ! -- dummy
    class(UzfType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: nn
    integer(I4B) :: iuzid
    real(DP) :: obsdepth
    real(DP) :: dmax
    character(len=LENBOUNDNAME) :: bname
    class(ObserveType), pointer :: obsrv => null()
    ! -- formats
60  format('Invalid node number in OBS input: ', i0)
    !
    if (kper == 1) then
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv
        !
        ! -- get node number 1
        nn = obsrv%NodeNumber
        if (nn == NAMEDBOUNDFLAG) then
          bname = obsrv%FeatureName
          !
          ! -- Observation location(s) is(are) based on a boundary name.
          !    Iterate through all boundaries to identify and store
          !    corresponding index(indices) in bound array.
          do j = 1, this%nodes
            if (this%boundname(j) == bname) then
              obsrv%BndFound = .true.
              obsrv%CurrentTimeStepEndValue = DZERO
              call obsrv%AddObsIndex(j)
              if (obsrv%indxbnds_count == 1) then
                !
                ! -- Define intPak1 so that obs_theta is stored (for first uzf
                !    cell if multiple cells share the same boundname).
                obsrv%intPak1 = j
              end if
            end if
          end do
        else
          !
          ! -- get node number
          nn = obsrv%NodeNumber
          !
          ! -- put nn (a value meaningful only to UZF) in intPak1
          obsrv%intPak1 = nn
          ! -- check that node number is valid; call store_error if not
          if (nn < 1 .or. nn > this%nodes) then
            write (errmsg, 60) nn
            call store_error(errmsg)
          else
            obsrv%BndFound = .true.
          end if
          obsrv%CurrentTimeStepEndValue = DZERO
          call obsrv%AddObsIndex(nn)
        end if
        !
        ! -- catch non-cumulative observation assigned to observation defined
        !    by a boundname that is assigned to more than one element
        if (obsrv%ObsTypeId == 'WATER-CONTENT') then
          n = obsrv%indxbnds_count
          if (n /= 1) then
            write (errmsg, '(a,3(1x,a))') &
              trim(adjustl(obsrv%ObsTypeId)), 'for observation', &
              trim(adjustl(obsrv%Name)), &
              'must be assigned to a UZF cell with a unique boundname.'
            call store_error(errmsg, terminate=.TRUE.)
          end if
          !
          ! -- check WATER-CONTENT depth
          obsdepth = obsrv%Obsdepth
          !
          ! -- put obsdepth (a value meaningful only to UZF) in dblPak1
          obsrv%dblPak1 = obsdepth
          !
          ! -- determine maximum cell depth
          ! -- This is presently complicated for landflag = 1 cells and surfdep
          !    greater than zero.  In this case, celtop is dis%top - surfdep.
          iuzid = obsrv%intPak1
          dmax = this%uzfobj%celtop(iuzid) - this%uzfobj%celbot(iuzid)
          ! -- check that obs depth is valid; call store_error if not
          ! -- need to think about a way to put bounds on this depth
          ! -- Also, an observation depth of 0.0, whether a landflag == 1 object
          ! -- or a subsurface object, is not legit since this would be at a
          ! -- a layer interface and therefore a discontinuity.
          if (obsdepth <= DZERO .or. obsdepth > dmax) then
            write (errmsg, '(a,3(1x,a),1x,g0,1x,a,1x,g0,a)') &
              trim(adjustl(obsrv%ObsTypeId)), 'for observation', &
              trim(adjustl(obsrv%Name)), 'specified depth (', obsdepth, &
              ') must be greater than 0.0 and less than ', dmax, '.'
            call store_error(errmsg)
          end if
        else
          do j = 1, obsrv%indxbnds_count
            nn = obsrv%indxbnds(j)
            if (nn < 1 .or. nn > this%maxbound) then
              write (errmsg, '(a,2(1x,a),1x,i0,1x,a,1x,i0,a)') &
                trim(adjustl(obsrv%ObsTypeId)), 'uzfno must be greater than 0 ', &
                'and less than or equal to', this%maxbound, &
                '(specified value is ', nn, ').'
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
  end subroutine uzf_rp_obs

  ! -- Procedures related to observations (NOT type-bound)

  !> @brief This procedure is pointed to by ObsDataType%ProcesssIdPtr
  !!
  !! Process the ID string of an observation definition for UZF-package
  !! observations
  !<
  subroutine uzf_process_obsID(obsrv, dis, inunitobs, iout)
    ! -- .
    ! -- dummy
    type(ObserveType), intent(inout) :: obsrv
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: inunitobs
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B) :: n, nn
    real(DP) :: obsdepth
    integer(I4B) :: icol, istart, istop, istat
    real(DP) :: r
    character(len=LINELENGTH) :: string
    ! formats
30  format(i10)
    !
    string = obsrv%IDstring
    ! -- Extract node number from string and store it.
    !    If 1st item is not an integer(I4B), it should be a
    !    feature name--deal with it.
    icol = 1
    ! -- get node number
    call urword(string, icol, istart, istop, 1, n, r, iout, inunitobs)
    read (string(istart:istop), 30, iostat=istat) nn
    if (istat == 0) then
      ! -- store uzf node number (NodeNumber)
      obsrv%NodeNumber = nn
    else
      ! Integer can't be read from string; it's presumed to be a boundary
      ! name (already converted to uppercase)
      obsrv%FeatureName = string(istart:istop)
      !obsrv%FeatureName = trim(adjustl(string))
      ! -- Observation may require summing rates from multiple boundaries,
      !    so assign NodeNumber as a value that indicates observation
      !    is for a named boundary or group of boundaries.
      obsrv%NodeNumber = NAMEDBOUNDFLAG
    end if
    !
    ! -- for soil water observation, store depth
    if (obsrv%ObsTypeId == 'WATER-CONTENT') then
      call urword(string, icol, istart, istop, 3, n, r, iout, inunitobs)
      obsdepth = r
      ! -- store observations depth
      obsrv%Obsdepth = obsdepth
    end if
  end subroutine uzf_process_obsID

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
    call this%uzfobjwork%dealloc()
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

  !> @brief Set up the budget object that stores all the uzf flows
  !!
  !! The terms listed here must correspond in number and order to the ones
  !! listed in the uzf_fill_budobj routine
  !<
  subroutine uzf_setup_budobj(this)
    ! -- modules
    use ConstantsModule, only: LENBUDTXT
    ! -- dummy
    class(UzfType) :: this
    ! -- local
    integer(I4B) :: nbudterm
    integer(I4B) :: maxlist, naux
    integer(I4B) :: idx
    integer(I4B) :: nlen
    integer(I4B) :: n, n1, n2
    integer(I4B) :: ivertflag
    real(DP) :: q
    character(len=LENBUDTXT) :: text
    character(len=LENBUDTXT), dimension(1) :: auxtxt
    !
    ! -- Determine the number of uzf to uzf connections
    nlen = 0
    do n = 1, this%nodes
      ivertflag = this%uzfobj%ivertcon(n)
      if (ivertflag > 0) then
        nlen = nlen + 1
      end if
    end do
    !
    ! -- Determine the number of uzf budget terms. These are fixed for
    !    the simulation and cannot change.  This includes FLOW-JA-FACE
    !    so they can be written to the binary budget files, but these internal
    !    flows are not included as part of the budget table.
    nbudterm = 4
    if (nlen > 0) nbudterm = nbudterm + 1
    if (this%ietflag /= 0) nbudterm = nbudterm + 1
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
    if (nlen > 0) then
      idx = idx + 1
      maxlist = nlen * 2
      naux = 1
      auxtxt(1) = '       FLOW-AREA'
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux, auxtxt, ordered_id1=.false.)
      !
      ! -- store connectivity
      call this%budobj%budterm(idx)%reset(nlen * 2)
      q = DZERO
      do n = 1, this%nodes
        ivertflag = this%uzfobj%ivertcon(n)
        if (ivertflag > 0) then
          n1 = n
          n2 = ivertflag
          call this%budobj%budterm(idx)%update_term(n1, n2, q)
          call this%budobj%budterm(idx)%update_term(n2, n1, -q)
        end if
      end do
    end if
    !
    ! --
    text = '             GWF'
    idx = idx + 1
    maxlist = this%nodes
    naux = 1
    auxtxt(1) = '       FLOW-AREA'
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%name_model, &
                                             maxlist, .false., .true., &
                                             naux, auxtxt)
    call this%budobj%budterm(idx)%reset(this%nodes)
    q = DZERO
    do n = 1, this%nodes
      n2 = this%igwfnode(n)
      this%qauxcbc(1) = this%uzfobj%uzfarea(n)
      call this%budobj%budterm(idx)%update_term(n, n2, q, this%qauxcbc)
    end do
    !
    ! --
    text = '    INFILTRATION'
    idx = idx + 1
    maxlist = this%nodes
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
    text = '         REJ-INF'
    idx = idx + 1
    maxlist = this%nodes
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
    text = '            UZET'
    if (this%ietflag /= 0) then
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
    text = '         STORAGE'
    idx = idx + 1
    maxlist = this%nodes
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
      maxlist = this%nodes
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
      text = '  REJ-INF-TO-MVR'
      idx = idx + 1
      maxlist = this%nodes
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
    ! -- if uzf flow for each reach are written to the listing file
    if (this%iprflow /= 0) then
      call this%budobj%flowtable_df(this%iout, cellids='GWF')
    end if
  end subroutine uzf_setup_budobj

  !> @brief Copy flow terms into this%budobj
  !<
  subroutine uzf_fill_budobj(this)
    ! -- dummy
    class(UzfType) :: this
    ! -- local
    integer(I4B) :: naux
    integer(I4B) :: nlen
    integer(I4B) :: ivertflag
    integer(I4B) :: n, n1, n2
    integer(I4B) :: idx
    real(DP) :: q
    real(DP) :: a
    real(DP) :: top
    real(DP) :: bot
    real(DP) :: thick
    real(DP) :: fm
    real(DP) :: v
    !
    ! -- initialize counter
    idx = 0
    !
    ! -- FLOW JA FACE
    nlen = 0
    do n = 1, this%nodes
      ivertflag = this%uzfobj%ivertcon(n)
      if (ivertflag > 0) then
        nlen = nlen + 1
      end if
    end do
    if (nlen > 0) then
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(nlen * 2)
      do n = 1, this%nodes
        ivertflag = this%uzfobj%ivertcon(n)
        if (ivertflag > 0) then
          a = this%uzfobj%uzfarea(n)
          q = this%uzfobj%surfluxbelow(n) * a
          this%qauxcbc(1) = a
          if (q > DZERO) then
            q = -q
          end if
          n1 = n
          n2 = ivertflag
          call this%budobj%budterm(idx)%update_term(n1, n2, q, this%qauxcbc)
          call this%budobj%budterm(idx)%update_term(n2, n1, -q, this%qauxcbc)
        end if
      end do
    end if
    !
    ! -- GWF (LEAKAGE)
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nodes)
    do n = 1, this%nodes
      this%qauxcbc(1) = this%uzfobj%uzfarea(n)
      n2 = this%igwfnode(n)
      q = -this%rch(n)
      call this%budobj%budterm(idx)%update_term(n, n2, q, this%qauxcbc)
    end do
    !
    ! -- INFILTRATION
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nodes)
    do n = 1, this%nodes
      q = this%appliedinf(n)
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- REJECTED INFILTRATION
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nodes)
    do n = 1, this%nodes
      q = this%rejinf(n)
      if (q > DZERO) then
        q = -q
      end if
      call this%budobj%budterm(idx)%update_term(n, n, q)
    end do
    !
    ! -- UNSATURATED EVT
    if (this%ietflag /= 0) then
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%nodes)
      do n = 1, this%nodes
        q = this%uzet(n)
        if (q > DZERO) then
          q = -q
        end if
        call this%budobj%budterm(idx)%update_term(n, n, q)
      end do
    end if
    !
    ! -- STORAGE
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nodes)
    do n = 1, this%nodes
      q = -this%qsto(n)
      top = this%uzfobj%celtop(n)
      bot = this%uzfobj%watab(n)
      thick = top - bot
      if (thick > DZERO) then
        fm = thick * (this%wcnew(n) - this%uzfobj%thtr(n))
        v = fm * this%uzfobj%uzfarea(n)
      else
        v = DZERO
      end if
      ! -- save mobile water volume into aux variable
      this%qauxcbc(1) = v
      call this%budobj%budterm(idx)%update_term(n, n, q, this%qauxcbc)
    end do
    !
    ! -- MOVER
    if (this%imover == 1) then
      !
      ! -- FROM MOVER
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%nodes)
      do n = 1, this%nodes
        q = this%pakmvrobj%get_qfrommvr(n)
        call this%budobj%budterm(idx)%update_term(n, n, q)
      end do
      !
      ! -- REJ-INF-TO-MVR
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%nodes)
      do n = 1, this%nodes
        q = this%rejinftomvr(n)
        if (q > DZERO) then
          q = -q
        end if
        call this%budobj%budterm(idx)%update_term(n, n, q)
      end do

    end if
    !
    ! -- AUXILIARY VARIABLES
    naux = this%naux
    if (naux > 0) then
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%nodes)
      do n = 1, this%nodes
        q = DZERO
        call this%budobj%budterm(idx)%update_term(n, n, q, this%auxvar(:, n))
      end do
    end if
    !
    ! --Terms are filled, now accumulate them for this time step
    call this%budobj%accumulate_terms()
  end subroutine uzf_fill_budobj

end module UzfModule
