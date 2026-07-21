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
    procedure :: apply_keyword
    procedure :: write_keyword
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
    this%icnvgopt = 0
    this%iter1 = 0
    this%ilinmeth = 0
    this%iscl = 0
    this%iord = 0
    this%north = 0
    this%relax = DZERO
    this%level = 0
    this%droptol = DZERO
    this%ifdparam = 0

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
      this%ilinmeth = 1
      this%iscl = 0
      this%iord = 0
      this%dvclose = DEM3
      this%rclose = DEM1
      this%relax = DZERO
      this%level = 0
      this%droptol = DZERO
      this%north = 0
    case (2) ! Moderate
      this%iter1 = 100
      this%ilinmeth = 2
      this%iscl = 0
      this%iord = 0
      this%dvclose = DEM2
      this%rclose = DEM1
      this%relax = 0.97D0
      this%level = 0
      this%droptol = DZERO
      this%north = 0
    case (3) ! Complex
      this%iter1 = 500
      this%ilinmeth = 2
      this%iscl = 0
      this%iord = 0
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

    call parser%GetBlock('LINEAR', block_found, ierr, supportOpenClose=.true., &
                         blockRequired=.FALSE.)

    if (block_found) then
      write (iout, '(/1x,a)') 'PROCESSING LINEAR DATA'
      do
        call parser%GetNextLine(end_of_block)
        if (end_of_block) exit
        call parser%GetStringCaps(keyword)
        call this%apply_keyword(parser, keyword)
      end do
      write (iout, '(1x,a)') 'END OF LINEAR DATA'
    else
      if (this%ifdparam == 0) THEN
        write (errmsg, '(a)') 'NO LINEAR block detected.'
        call store_error(errmsg)
      end if
    end if

  end subroutine read_from_file

  !> @brief Apply a single LINEAR-block keyword to the settings
  !!
  !! Reads the value(s) for the given keyword from the parser and stores them in
  !! the settings. Shared by the LINEAR-block reader and the period-varying
  !! settings reader so the keyword vocabulary has a single definition.
  !<
  subroutine apply_keyword(this, parser, keyword)
    class(ImsLinearSettingsType) :: this !< linear settings
    type(BlockParserType) :: parser !< block parser
    character(len=*), intent(in) :: keyword !< keyword (already upper-cased)
    ! local
    character(len=LINELENGTH) :: errmsg
    character(len=LINELENGTH) :: kw
    integer(I4B) :: iscaling, iordering

    select case (keyword)
    case ('INNER_DVCLOSE')
      this%dvclose = parser%GetDouble()
    case ('INNER_RCLOSE')
      this%rclose = parser%GetDouble()
      ! -- look for additional key words
      call parser%GetStringCaps(kw)
      if (kw == 'STRICT') then
        this%icnvgopt = 1
      else if (kw == 'L2NORM_RCLOSE') then
        this%icnvgopt = 2
      else if (kw == 'RELATIVE_RCLOSE') then
        this%icnvgopt = 3
      else if (kw == 'L2NORM_RELATIVE_RCLOSE') then
        this%icnvgopt = 4
      end if
    case ('INNER_MAXIMUM')
      this%iter1 = parser%GetInteger()
    case ('LINEAR_ACCELERATION')
      call parser%GetStringCaps(kw)
      if (kw == 'CG') then
        this%ilinmeth = 1
      else if (kw == 'BICGSTAB') then
        this%ilinmeth = 2
      else
        this%ilinmeth = 0
        write (errmsg, '(3a)') &
          'Unknown IMSLINEAR LINEAR_ACCELERATION method (', trim(kw), ').'
        call store_error(errmsg)
      end if
    case ('SCALING_METHOD')
      call parser%GetStringCaps(kw)
      iscaling = 0
      if (kw == 'NONE') then
        iscaling = 0
      else if (kw == 'DIAGONAL') then
        iscaling = 1
      else if (kw == 'L2NORM') then
        iscaling = 2
      else
        write (errmsg, '(3a)') &
          'Unknown IMSLINEAR SCALING_METHOD (', trim(kw), ').'
        call store_error(errmsg)
      end if
      this%iscl = iscaling
    case ('RED_BLACK_ORDERING')
      iordering = 0
    case ('REORDERING_METHOD')
      call parser%GetStringCaps(kw)
      iordering = 0
      if (kw == 'NONE') then
        iordering = 0
      else if (kw == 'RCM') then
        iordering = 1
      else if (kw == 'MD') then
        iordering = 2
      else
        write (errmsg, '(3a)') &
          'Unknown IMSLINEAR REORDERING_METHOD (', trim(kw), ').'
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

  end subroutine apply_keyword

  !> @brief Write the stored value of a single LINEAR-block keyword
  !!
  !! Echoes the value currently held in the settings for the given keyword, so a
  !! setting that apply_keyword can change can also be reported. Shares the
  !! keyword vocabulary with apply_keyword; keywords whose value is not reported
  !! individually are skipped.
  !<
  subroutine write_keyword(this, iout, keyword)
    class(ImsLinearSettingsType) :: this !< linear settings
    integer(I4B), intent(in) :: iout !< listing file unit
    character(len=*), intent(in) :: keyword !< keyword (already upper-cased)
    ! local
    character(len=LINELENGTH) :: cval

    if (iout <= 0) return

    select case (keyword)
    case ('INNER_DVCLOSE', 'INNER_HCLOSE')
      write (cval, '(1pe15.5)') this%dvclose
    case ('INNER_RCLOSE')
      write (cval, '(1pe15.5)') this%rclose
    case ('INNER_MAXIMUM')
      write (cval, '(i0)') this%iter1
    case ('RELAXATION_FACTOR')
      write (cval, '(1pe15.5)') this%relax
    case ('PRECONDITIONER_LEVELS')
      write (cval, '(i0)') this%level
    case ('PRECONDITIONER_DROP_TOLERANCE')
      write (cval, '(1pe15.5)') this%droptol
    case default
      return
    end select

    write (iout, '(4x,3a)') trim(keyword), ' = ', trim(adjustl(cval))
    !
    ! -- INNER_RCLOSE carries an optional trailing keyword that also selects the
    !    residual convergence option, so report that alongside the value
    if (keyword == 'INNER_RCLOSE') then
      write (iout, '(4x,a,i0)') 'RESIDUAL CONVERGENCE OPTION = ', this%icnvgopt
    end if

  end subroutine write_keyword

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

end module
