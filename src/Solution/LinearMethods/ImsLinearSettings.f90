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
    character(len=LINELENGTH) :: warnmsg
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
            this%icnvgopt = 1
          else if (keyword == 'L2NORM_RCLOSE') then
            this%icnvgopt = 2
          else if (keyword == 'RELATIVE_RCLOSE') then
            this%icnvgopt = 3
          else if (keyword == 'L2NORM_RELATIVE_RCLOSE') then
            this%icnvgopt = 4
          end if
        case ('INNER_MAXIMUM')
          this%iter1 = parser%GetInteger()
        case ('LINEAR_ACCELERATION')
          call parser%GetStringCaps(keyword)
          if (keyword .eq. 'CG') then
            this%ilinmeth = 1
          else if (keyword .eq. 'BICGSTAB') then
            this%ilinmeth = 2
          else
            this%ilinmeth = 0
            write (errmsg, '(3a)') &
              'Unknown IMSLINEAR LINEAR_ACCELERATION method (', &
              trim(keyword), ').'
            call store_error(errmsg)
          end if
        case ('SCALING_METHOD')
          call parser%GetStringCaps(keyword)
          iscaling = 0
          if (keyword .eq. 'NONE') then
            iscaling = 0
          else if (keyword .eq. 'DIAGONAL') then
            iscaling = 1
          else if (keyword .eq. 'L2NORM') then
            iscaling = 2
          else
            write (errmsg, '(3a)') &
              'Unknown IMSLINEAR SCALING_METHOD (', trim(keyword), ').'
            call store_error(errmsg)
          end if
          this%iscl = iscaling
        case ('RED_BLACK_ORDERING')
          iordering = 0
        case ('REORDERING_METHOD')
          call parser%GetStringCaps(keyword)
          iordering = 0
          if (keyword == 'NONE') then
            iordering = 0
          else if (keyword == 'RCM') then
            iordering = 1
          else if (keyword == 'MD') then
            iordering = 2
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
          ! -- default
        case default
          write (errmsg, '(3a)') &
            'Unknown IMSLINEAR keyword (', trim(keyword), ').'
          call store_error(errmsg)
        end select
      end do
      write (iout, '(1x,a)') 'END OF LINEAR DATA'
    else
      if (this%ifdparam == 0) THEN
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

end module
