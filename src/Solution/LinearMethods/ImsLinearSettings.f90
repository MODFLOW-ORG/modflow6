module ImsLinearSettingsModule
  use KindModule
  use ConstantsModule
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use BlockParserModule, only: BlockParserType
  use SimModule, only: store_error, store_warning, deprecation_warning
  implicit none
  private

  integer(I4B), public, parameter :: CG_METHOD = 1
  integer(I4B), public, parameter :: BCGS_METHOD = 2

  !> default number of levels in the AMG hierarchy
  integer(I4B), public, parameter :: AMG_LEVELS_DEFAULT = 10

  !> @brief IMS linear preconditioner types
  !<
  enum, bind(C)
    enumerator :: IPC_UNKNOWN = 0 !< preconditioner type not set
    enumerator :: IPC_ILU0 = 1 !< ILU0 (incomplete LU, zero fill)
    enumerator :: IPC_MILU0 = 2 !< modified ILU0
    enumerator :: IPC_ILUT = 3 !< ILUT (incomplete LU with threshold)
    enumerator :: IPC_MILUT = 4 !< modified ILUT
    enumerator :: IPC_AMG = 5 !< algebraic multigrid
  end enum
  public :: IPC_UNKNOWN, IPC_ILU0, IPC_MILU0, IPC_ILUT, IPC_MILUT, IPC_AMG
  public :: resolve_ipc

  !> @brief IMS AMG smoother types
  !<
  enum, bind(C)
    enumerator :: SMOOTHER_ILU0 = 1 !< ILU(0) smoother at finest level (default)
    enumerator :: SMOOTHER_L1GS = 2 !< l1 Gauss-Seidel smoother at all AMG levels
    enumerator :: SMOOTHER_ILU0_ALL = 3 !< ILU(0) smoother at all AMG levels
  end enum
  public :: SMOOTHER_ILU0, SMOOTHER_L1GS, SMOOTHER_ILU0_ALL

  !> @brief IMS matrix scaling methods
  !<
  enum, bind(C)
    enumerator :: SCL_NONE = 0 !< no scaling (default)
    enumerator :: SCL_DIAGONAL = 1 !< diagonal (symmetric) scaling
    enumerator :: SCL_L2NORM = 2 !< L2-norm scaling
  end enum
  public :: SCL_NONE, SCL_DIAGONAL, SCL_L2NORM

  !> @brief IMS matrix reordering methods
  !<
  enum, bind(C)
    enumerator :: ORD_NONE = 0 !< original ordering (default)
    enumerator :: ORD_RCM = 1 !< reverse Cuthill-McKee ordering
    enumerator :: ORD_MINIMUM_DEGREE = 2 !< minimum degree ordering
  end enum
  public :: ORD_NONE, ORD_RCM, ORD_MINIMUM_DEGREE

  !> @brief IMS residual convergence norm options
  !<
  enum, bind(C)
    enumerator :: CNVG_INF_NORM = 0 !< infinity norm (default)
    enumerator :: CNVG_INF_NORM_STRICT = 1 !< strict infinity norm
    enumerator :: CNVG_L2_NORM = 2 !< L2 norm
    enumerator :: CNVG_REL_L2_NORM = 3 !< relative L2 norm (RCLOSE may be 0)
    enumerator :: CNVG_L2_NORM_REL = 4 !< L2 norm with relative component
  end enum
  public :: CNVG_INF_NORM, CNVG_INF_NORM_STRICT, CNVG_L2_NORM, &
            CNVG_REL_L2_NORM, CNVG_L2_NORM_REL

  type, public :: ImsLinearSettingsType
    character(len=LENMEMPATH) :: memory_path
    real(DP), pointer :: dvclose => null() !< dependent variable closure criterion
    real(DP), pointer :: rclose => null() !< residual closure criterion
    integer(I4B), pointer :: icnvgopt => null() !< convergence option
    integer(I4B), pointer :: iter1 => null() !< max. iterations
    integer(I4B), pointer :: ilinmeth => null() !< linear solver method
    integer(I4B), pointer :: iscl => null() !< scaling method
    integer(I4B), pointer :: iord => null() !< reordering method
    integer(I4B), pointer :: north => null() !< number of orthogonalizations
    real(DP), pointer :: relax => null() !< relaxation factor
    integer(I4B), pointer :: level => null() !< nr. of preconditioner levels
    real(DP), pointer :: droptol => null() !< drop tolerance for preconditioner
    integer(I4B), pointer :: ifdparam => null() !< complexity option
    integer(I4B), pointer :: ipc_type => null() !< preconditioner type (0=ILU default, 1=AMG)
    logical(LGP) :: level_specified = .false. !< PRECONDITIONER_LEVELS was read from file
    integer(I4B), pointer :: nsmooth => null() !< AMG smoother iterations per level
    integer(I4B), pointer :: smoother_type => null() !< AMG smoother type
    real(DP), pointer :: strength_threshold => null() !< AMG strength-of-connection threshold
  contains
    procedure :: init
    procedure :: preset_config
    procedure :: read_from_file
    procedure :: check_settings
    procedure :: destroy
  end type

contains

  subroutine init(this, mem_path)
    use MemoryHelperModule, only: create_mem_path
    class(ImsLinearSettingsType) :: this !< linear settings
    character(len=LENMEMPATH) :: mem_path !< solution memory path

    this%memory_path = create_mem_path(mem_path, 'IMSLINEAR')

    call mem_allocate(this%dvclose, 'DVCLOSE', this%memory_path)
    call mem_allocate(this%rclose, 'RCLOSE', this%memory_path)
    call mem_allocate(this%icnvgopt, 'ICNVGOPT', this%memory_path)
    call mem_allocate(this%iter1, 'ITER1', this%memory_path)
    call mem_allocate(this%ilinmeth, 'ILINMETH', this%memory_path)
    call mem_allocate(this%iscl, 'ISCL', this%memory_path)
    call mem_allocate(this%iord, 'IORD', this%memory_path)
    call mem_allocate(this%north, 'NORTH', this%memory_path)
    call mem_allocate(this%relax, 'RELAX', this%memory_path)
    call mem_allocate(this%level, 'LEVEL', this%memory_path)
    call mem_allocate(this%droptol, 'DROPTOL', this%memory_path)
    call mem_allocate(this%ifdparam, 'IDFPARAM', this%memory_path)
    call mem_allocate(this%ipc_type, 'IPCTYPE', this%memory_path)
    call mem_allocate(this%nsmooth, 'NSMOOTH', this%memory_path)
    call mem_allocate(this%smoother_type, 'SMOOTHERTYPE', this%memory_path)
    call mem_allocate(this%strength_threshold, 'STRTHRESH', this%memory_path)

    ! defaults
    this%dvclose = DZERO
    this%rclose = DZERO
    this%icnvgopt = CNVG_INF_NORM
    this%iter1 = 0
    this%ilinmeth = 0
    this%iscl = SCL_NONE
    this%iord = ORD_NONE
    this%north = 0
    this%relax = DZERO
    this%level = 0
    this%droptol = DZERO
    this%ifdparam = 0
    this%ipc_type = 0
    this%nsmooth = 2
    this%smoother_type = SMOOTHER_ILU0
    this%strength_threshold = 0.25_dp

  end subroutine init

  !> @brief Set solver pre-configured settings based on complexity option
  !<
  subroutine preset_config(this, idfparam)
    class(ImsLinearSettingsType) :: this !< linear settings
    integer(I4B) :: idfparam !< complexity option

    this%ifdparam = idfparam

    select case (idfparam)
    case (1) ! Simple option
      this%iter1 = 50
      this%ilinmeth = CG_METHOD
      this%iscl = SCL_NONE
      this%iord = ORD_NONE
      this%dvclose = DEM3
      this%rclose = DEM1
      this%relax = DZERO
      this%level = 0
      this%droptol = DZERO
      this%north = 0
    case (2) ! Moderate
      this%iter1 = 100
      this%ilinmeth = BCGS_METHOD
      this%iscl = SCL_NONE
      this%iord = ORD_NONE
      this%dvclose = DEM2
      this%rclose = DEM1
      this%relax = 0.97d0
      this%level = 0
      this%droptol = DZERO
      this%north = 0
    case (3) ! Complex
      this%iter1 = 500
      this%ilinmeth = BCGS_METHOD
      this%iscl = SCL_NONE
      this%iord = ORD_NONE
      this%dvclose = DEM1
      this%rclose = DEM1
      this%relax = DZERO
      this%level = 5
      this%droptol = DEM4
      this%north = 2
    end select

  end subroutine preset_config

  !> @brief Read the settings for the linear solver from the .ims file,
  !< overriding a possible pre-set configuration with set_complexity
  subroutine read_from_file(this, parser, iout)
    class(ImsLinearSettingsType) :: this !< linear settings
    type(BlockParserType) :: parser !< block parser
    integer(I4B) :: iout !< listing file
    ! local
    logical(LGP) :: block_found, end_of_block
    integer(I4B) :: ierr
    character(len=LINELENGTH) :: errmsg
    character(len=LINELENGTH) :: warnmsg
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: iscaling, iordering

    call parser%GetBlock('LINEAR', block_found, ierr, supportOpenClose=.true., &
                         blockRequired=.false.)

    if (block_found) then
      write (iout, '(/1x,a)') 'PROCESSING LINEAR DATA'
      do
        call parser%GetNextLine(end_of_block)
        if (end_of_block) exit
        call parser%GetStringCaps(keyword)
        ! -- parse keyword
        select case (keyword)
        case ('INNER_DVCLOSE')
          this%dvclose = parser%GetDouble()
        case ('INNER_RCLOSE')
          this%rclose = parser%GetDouble()
          ! -- look for additional key words
          call parser%GetStringCaps(keyword)
          if (keyword == 'STRICT') then
            this%icnvgopt = CNVG_INF_NORM_STRICT
          else if (keyword == 'L2NORM_RCLOSE') then
            this%icnvgopt = CNVG_L2_NORM
          else if (keyword == 'RELATIVE_RCLOSE') then
            this%icnvgopt = CNVG_REL_L2_NORM
          else if (keyword == 'L2NORM_RELATIVE_RCLOSE') then
            this%icnvgopt = CNVG_L2_NORM_REL
          end if
        case ('INNER_MAXIMUM')
          this%iter1 = parser%GetInteger()
        case ('LINEAR_ACCELERATION')
          call parser%GetStringCaps(keyword)
          if (keyword .eq. 'CG') then
            this%ilinmeth = CG_METHOD
          else if (keyword .eq. 'BICGSTAB') then
            this%ilinmeth = BCGS_METHOD
          else
            this%ilinmeth = 0
            write (errmsg, '(3a)') &
              'Unknown IMSLINEAR LINEAR_ACCELERATION method (', &
              trim(keyword), ').'
            call store_error(errmsg)
          end if
        case ('SCALING_METHOD')
          call parser%GetStringCaps(keyword)
          iscaling = SCL_NONE
          if (keyword .eq. 'NONE') then
            iscaling = SCL_NONE
          else if (keyword .eq. 'DIAGONAL') then
            iscaling = SCL_DIAGONAL
          else if (keyword .eq. 'L2NORM') then
            iscaling = SCL_L2NORM
          else
            write (errmsg, '(3a)') &
              'Unknown IMSLINEAR SCALING_METHOD (', trim(keyword), ').'
            call store_error(errmsg)
          end if
          this%iscl = iscaling
        case ('RED_BLACK_ORDERING')
          iordering = ORD_NONE
        case ('REORDERING_METHOD')
          call parser%GetStringCaps(keyword)
          iordering = ORD_NONE
          if (keyword == 'NONE') then
            iordering = ORD_NONE
          else if (keyword == 'RCM') then
            iordering = ORD_RCM
          else if (keyword == 'MD') then
            iordering = ORD_MINIMUM_DEGREE
          else
            write (errmsg, '(3a)') &
              'Unknown IMSLINEAR REORDERING_METHOD (', trim(keyword), ').'
            call store_error(errmsg)
          end if
          this%iord = iordering
        case ('NUMBER_ORTHOGONALIZATIONS')
          this%north = parser%GetInteger()
        case ('RELAXATION_FACTOR')
          this%relax = parser%GetDouble()
        case ('PRECONDITIONER_LEVELS')
          this%level = parser%GetInteger()
          this%level_specified = .true.
          if (this%level < 0) then
            write (errmsg, '(a,1x,a)') &
              'IMSLINEAR PRECONDITIONER_LEVELS must be greater than', &
              'or equal to zero'
            call store_error(errmsg)
          end if
        case ('PRECONDITIONER_DROP_TOLERANCE')
          this%droptol = parser%GetDouble()
          if (this%droptol < DZERO) then
            write (errmsg, '(a,1x,a)') &
              'IMSLINEAR PRECONDITIONER_DROP_TOLERANCE', &
              'must be greater than or equal to zero'
            call store_error(errmsg)
          end if
          !
          ! -- deprecated variables
        case ('INNER_HCLOSE')
          this%dvclose = parser%GetDouble()
          !
          ! -- create warning message
          write (warnmsg, '(a)') &
            'SETTING INNER_DVCLOSE TO INNER_HCLOSE VALUE'
          !
          ! -- create deprecation warning
          call deprecation_warning('LINEAR', 'INNER_HCLOSE', '6.1.1', &
                                   warnmsg, parser%GetUnit())
          !
          ! -- preconditioner type
        case ('PRECONDITIONER_TYPE')
          call parser%GetStringCaps(keyword)
          if (keyword == 'AMG') then
            this%ipc_type = 1
          else if (keyword == 'ILU') then
            this%ipc_type = 0
          else
            write (errmsg, '(3a)') &
              'Unknown IMSLINEAR PRECONDITIONER_TYPE (', &
              trim(keyword), ').'
            call store_error(errmsg)
          end if
          !
          ! -- AMG smoother iterations per level
        case ('NUMBER_SMOOTHING_ITERATIONS')
          this%nsmooth = parser%GetInteger()
          if (this%nsmooth < 1) then
            write (errmsg, '(a)') &
              'IMSLINEAR NUMBER_SMOOTHING_ITERATIONS must be greater than zero'
            call store_error(errmsg)
          end if
          !
          ! -- AMG smoother type
        case ('SMOOTHER_TYPE')
          call parser%GetStringCaps(keyword)
          if (keyword == 'ILU0') then
            this%smoother_type = SMOOTHER_ILU0
          else if (keyword == 'L1GS') then
            this%smoother_type = SMOOTHER_L1GS
          else if (keyword == 'ILU0_ALL') then
            this%smoother_type = SMOOTHER_ILU0_ALL
          else
            write (errmsg, '(3a)') &
              'Unknown IMSLINEAR SMOOTHER_TYPE (', trim(keyword), &
              '). Valid values are ILU0, L1GS, and ILU0_ALL.'
            call store_error(errmsg)
          end if
          !
          ! -- AMG strength-of-connection threshold for aggregation
        case ('PRECONDITIONER_STRENGTH_THRESHOLD')
          this%strength_threshold = parser%GetDouble()
          if (this%strength_threshold < DZERO .or. &
              this%strength_threshold >= DONE) then
            write (errmsg, '(a)') &
              'IMSLINEAR PRECONDITIONER_STRENGTH_THRESHOLD must be '// &
              'greater than or equal to 0.0 and less than 1.0'
            call store_error(errmsg)
          end if
          !
          ! -- default
        case default
          write (errmsg, '(3a)') &
            'Unknown IMSLINEAR keyword (', trim(keyword), ').'
          call store_error(errmsg)
        end select
      end do
      write (iout, '(1x,a)') 'END OF LINEAR DATA'
    else
      if (this%ifdparam == 0) then
        write (errmsg, '(a)') 'NO LINEAR block detected.'
        call store_error(errmsg)
      end if
    end if

  end subroutine read_from_file

  !> @brief Check the settings after reading the configuration from file
  !<
  subroutine check_settings(this)
    class(ImsLinearSettingsType) :: this !< linear settings
    ! local
    character(len=LINELENGTH) :: warnmsg

    if (this%level == 0 .and. this%droptol > 0.0_dp) then
      write (warnmsg, '(a)') "PRECONDITIONER_DROP_TOLERANCE is ignored because &
                             &PRECONDITIONER_LEVELS equals zero."
      call store_warning(warnmsg)
    end if

    ! -- the COMPLEXITY presets set PRECONDITIONER_LEVELS for the ILU fill
    !    level, which is not a useful hierarchy depth for AMG
    if (this%ipc_type == 1 .and. .not. this%level_specified) then
      this%level = AMG_LEVELS_DEFAULT
    end if

  end subroutine check_settings

  subroutine destroy(this)
    class(ImsLinearSettingsType) :: this !< linear settings

    call mem_deallocate(this%dvclose)
    call mem_deallocate(this%rclose)
    call mem_deallocate(this%icnvgopt)
    call mem_deallocate(this%iter1)
    call mem_deallocate(this%ilinmeth)
    call mem_deallocate(this%iscl)
    call mem_deallocate(this%iord)
    call mem_deallocate(this%north)
    call mem_deallocate(this%relax)
    call mem_deallocate(this%level)
    call mem_deallocate(this%droptol)
    call mem_deallocate(this%ifdparam)
    call mem_deallocate(this%ipc_type)
    call mem_deallocate(this%nsmooth)
    call mem_deallocate(this%smoother_type)
    call mem_deallocate(this%strength_threshold)

  end subroutine destroy

  !> @brief Resolve the preconditioner enum from the ILU controls
  !!
  !! Maps the ILU fill level and relaxation factor to a specific IPC_*
  !! preconditioner:
  !!   level > 0   -> IPC_ILUT  (-> IPC_MILUT when relax > 0)
  !!   level == 0  -> IPC_ILU0  (-> IPC_MILU0 when relax > 0)
  !!
  !! Extracted into a pure function so the same resolution can be reused when
  !! preconditioner settings change at runtime (e.g. by stress period).
  !<
  pure function resolve_ipc(ipc_type, level, relax) result(ipc)
    integer(I4B), intent(in) :: ipc_type !< preconditioner-type keyword (0 = ILU, 1 = AMG)
    integer(I4B), intent(in) :: level !< ILU fill level
    real(DP), intent(in) :: relax !< relaxation factor (> 0 selects the modified ILU)
    integer(I4B) :: ipc !< resolved preconditioner enum (IPC_*)
    !
    ! -- AMG is selected explicitly; otherwise map the ILU controls to the
    !    target enum (do not rely on IPC_* being consecutive); relax > 0
    !    selects the modified variant
    if (ipc_type == 1) then
      ipc = IPC_AMG
    else if (level > 0) then
      if (relax > DZERO) then
        ipc = IPC_MILUT
      else
        ipc = IPC_ILUT
      end if
    else
      if (relax > DZERO) then
        ipc = IPC_MILU0
      else
        ipc = IPC_ILU0
      end if
    end if
  end function resolve_ipc

end module
