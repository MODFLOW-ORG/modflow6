module ImsLinearSettingsModule
  use KindModule
  use ConstantsModule
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use BlockParserModule, only: BlockParserType
  use SimModule, only: store_error, store_warning
  implicit none
  private

  !> @brief IMS complexity options
  !<
  enum, bind(C)
    enumerator :: COMPLEXITY_UNSET = 0 !< no complexity option
    enumerator :: COMPLEXITY_SIMPLE = 1 !< SIMPLE
    enumerator :: COMPLEXITY_MODERATE = 2 !< MODERATE
    enumerator :: COMPLEXITY_COMPLEX = 3 !< COMPLEX
  end enum
  public :: COMPLEXITY_UNSET, COMPLEXITY_SIMPLE, COMPLEXITY_MODERATE, &
            COMPLEXITY_COMPLEX

  !> @brief IMS print options
  !<
  enum, bind(C)
    enumerator :: PRINT_NONE = 0 !< no solver output
    enumerator :: PRINT_SUMMARY = 1 !< summary of solver iterations
    enumerator :: PRINT_ALL = 2 !< detailed solver output
  end enum
  public :: PRINT_NONE, PRINT_SUMMARY, PRINT_ALL

  !> @brief IMS linear acceleration methods
  !<
  enum, bind(C)
    enumerator :: LINACC_UNKNOWN = 0 !< acceleration method not set
    enumerator :: LINACC_CG = 1 !< conjugate gradient
    enumerator :: LINACC_BICGSTAB = 2 !< BiCGSTAB
  end enum
  public :: LINACC_UNKNOWN, LINACC_CG, LINACC_BICGSTAB

  !> @brief IMS inner convergence options
  !<
  enum, bind(C)
    enumerator :: CNVG_INF_NORM = 0 !< infinity norm
    enumerator :: CNVG_INF_NORM_STRICT = 1 !< infinity norm, STRICT
    enumerator :: CNVG_L2_NORM = 2 !< L2 norm, L2NORM_RCLOSE
    enumerator :: CNVG_L2_NORM_RELATIVE = 3 !< relative L2 norm, RELATIVE_RCLOSE
    enumerator :: CNVG_L2_NORM_REL_RCLOSE = 4 !< L2NORM_RELATIVE_RCLOSE
  end enum
  public :: CNVG_INF_NORM, CNVG_INF_NORM_STRICT, CNVG_L2_NORM, &
            CNVG_L2_NORM_RELATIVE, CNVG_L2_NORM_REL_RCLOSE

  !> @brief Inner (linear) convergence status
  !<
  enum, bind(C)
    enumerator :: INNER_FORCE_OUTER = -1 !< stop, but force another outer iteration
    enumerator :: INNER_NOT_CONVERGED = 0 !< not converged
    enumerator :: INNER_CONVERGED = 1 !< converged
  end enum
  public :: INNER_FORCE_OUTER, INNER_NOT_CONVERGED, INNER_CONVERGED

  !> @brief IMS matrix scaling methods
  !<
  enum, bind(C)
    enumerator :: SCALING_NONE = 0 !< no scaling
    enumerator :: SCALING_DIAGONAL = 1 !< symmetric diagonal scaling
    enumerator :: SCALING_L2NORM = 2 !< L2 norm scaling
  end enum
  public :: SCALING_NONE, SCALING_DIAGONAL, SCALING_L2NORM

  !> @brief IMS matrix reordering methods
  !<
  enum, bind(C)
    enumerator :: ORDER_NONE = 0 !< original ordering
    enumerator :: ORDER_RCM = 1 !< reverse Cuthill-McKee
    enumerator :: ORDER_MD = 2 !< minimum degree
  end enum
  public :: ORDER_NONE, ORDER_RCM, ORDER_MD

  !> @brief IMS linear preconditioner types
  !<
  enum, bind(C)
    enumerator :: IPC_UNKNOWN = 0 !< preconditioner type not set
    enumerator :: IPC_ILU0 = 1 !< ILU0 (incomplete LU, zero fill)
    enumerator :: IPC_MILU0 = 2 !< modified ILU0
    enumerator :: IPC_ILUT = 3 !< ILUT (incomplete LU with threshold)
    enumerator :: IPC_MILUT = 4 !< modified ILUT
  end enum
  public :: IPC_UNKNOWN, IPC_ILU0, IPC_MILU0, IPC_ILUT, IPC_MILUT
  public :: resolve_ipc
  public :: uses_l2_norm

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

    ! defaults
    this%dvclose = DZERO
    this%rclose = DZERO
    this%icnvgopt = CNVG_INF_NORM
    this%iter1 = 0
    this%ilinmeth = LINACC_UNKNOWN
    this%iscl = SCALING_NONE
    this%iord = ORDER_NONE
    this%north = 0
    this%relax = DZERO
    this%level = 0
    this%droptol = DZERO
    this%ifdparam = COMPLEXITY_UNSET

  end subroutine init

  !> @brief Set solver pre-configured settings based on complexity option
  !<
  subroutine preset_config(this, idfparam)
    class(ImsLinearSettingsType) :: this !< linear settings
    integer(I4B) :: idfparam !< complexity option

    this%ifdparam = idfparam

    select case (idfparam)
    case (COMPLEXITY_SIMPLE)
      this%iter1 = 50
      this%ilinmeth = LINACC_CG
      this%iscl = SCALING_NONE
      this%iord = ORDER_NONE
      this%dvclose = DEM3
      this%rclose = DEM1
      this%relax = DZERO
      this%level = 0
      this%droptol = DZERO
      this%north = 0
    case (COMPLEXITY_MODERATE)
      this%iter1 = 100
      this%ilinmeth = LINACC_BICGSTAB
      this%iscl = SCALING_NONE
      this%iord = ORDER_NONE
      this%dvclose = DEM2
      this%rclose = DEM1
      this%relax = 0.97D0
      this%level = 0
      this%droptol = DZERO
      this%north = 0
    case (COMPLEXITY_COMPLEX)
      this%iter1 = 500
      this%ilinmeth = LINACC_BICGSTAB
      this%iscl = SCALING_NONE
      this%iord = ORDER_NONE
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
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: iscaling, iordering

    call parser%GetBlock('LINEAR', block_found, ierr, supportOpenClose=.true., &
                         blockRequired=.FALSE.)

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
            this%icnvgopt = CNVG_L2_NORM_RELATIVE
          else if (keyword == 'L2NORM_RELATIVE_RCLOSE') then
            this%icnvgopt = CNVG_L2_NORM_REL_RCLOSE
          end if
        case ('INNER_MAXIMUM')
          this%iter1 = parser%GetInteger()
        case ('LINEAR_ACCELERATION')
          call parser%GetStringCaps(keyword)
          if (keyword .eq. 'CG') then
            this%ilinmeth = LINACC_CG
          else if (keyword .eq. 'BICGSTAB') then
            this%ilinmeth = LINACC_BICGSTAB
          else
            this%ilinmeth = LINACC_UNKNOWN
            write (errmsg, '(3a)') &
              'Unknown IMSLINEAR LINEAR_ACCELERATION method (', &
              trim(keyword), ').'
            call store_error(errmsg)
          end if
        case ('SCALING_METHOD')
          call parser%GetStringCaps(keyword)
          iscaling = SCALING_NONE
          if (keyword .eq. 'NONE') then
            iscaling = SCALING_NONE
          else if (keyword .eq. 'DIAGONAL') then
            iscaling = SCALING_DIAGONAL
          else if (keyword .eq. 'L2NORM') then
            iscaling = SCALING_L2NORM
          else
            write (errmsg, '(3a)') &
              'Unknown IMSLINEAR SCALING_METHOD (', trim(keyword), ').'
            call store_error(errmsg)
          end if
          this%iscl = iscaling
        case ('RED_BLACK_ORDERING')
          iordering = ORDER_NONE
        case ('REORDERING_METHOD')
          call parser%GetStringCaps(keyword)
          iordering = ORDER_NONE
          if (keyword == 'NONE') then
            iordering = ORDER_NONE
          else if (keyword == 'RCM') then
            iordering = ORDER_RCM
          else if (keyword == 'MD') then
            iordering = ORDER_MD
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
          ! -- default
        case default
          write (errmsg, '(3a)') &
            'Unknown IMSLINEAR keyword (', trim(keyword), ').'
          call store_error(errmsg)
        end select
      end do
      write (iout, '(1x,a)') 'END OF LINEAR DATA'
    end if

  end subroutine read_from_file

  !> @brief Check the settings after reading the configuration from file
  !<
  subroutine check_settings(this)
    class(ImsLinearSettingsType) :: this !< linear settings
    ! local
    character(len=LINELENGTH) :: warnmsg

    if (this%level == 0 .and. this%droptol > 0.0_DP) then
      write (warnmsg, '(a)') "PRECONDITIONER_DROP_TOLERANCE is ignored because &
                             &PRECONDITIONER_LEVELS equals zero."
      call store_warning(warnmsg)
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
  pure function resolve_ipc(level, relax) result(ipc)
    integer(I4B), intent(in) :: level !< ILU fill level
    real(DP), intent(in) :: relax !< relaxation factor (> 0 selects the modified ILU)
    integer(I4B) :: ipc !< resolved preconditioner enum (IPC_*)
    !
    ! -- map explicitly to the target enum (do not rely on IPC_* being
    !    consecutive); relax > 0 selects the modified variant
    if (level > 0) then
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

  !> @brief Return true if the convergence option tests the L2 norm
  !<
  pure function uses_l2_norm(icnvgopt) result(l2)
    integer(I4B), intent(in) :: icnvgopt !< convergence option (CNVG_*)
    logical(LGP) :: l2 !< true for the L2 norm options

    select case (icnvgopt)
    case (CNVG_L2_NORM, CNVG_L2_NORM_RELATIVE, CNVG_L2_NORM_REL_RCLOSE)
      l2 = .true.
    case default
      l2 = .false.
    end select
  end function uses_l2_norm

end module
