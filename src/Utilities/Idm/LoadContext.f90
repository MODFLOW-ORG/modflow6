!> @brief Load context for IDM generic dynamic loaders.
!!
!! LoadContextType classifies each input), builds the in-scope parameter
!! list for the active block, and manages memory-manager dimension scalars
!! and array pointers.  It is used by all dynamic loaders: ListLoadType,
!! KeystringLoadType, LayerArrayLoadType, GridArrayLoadType, and
!! StructArray-based static loads.
!!
!<
module LoadContextModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, IZERO, LINELENGTH, LENAUXNAME, &
                             LENVARNAME, LENBOUNDNAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use ModflowInputModule, only: ModflowInputType
  use InputDefinitionModule, only: InputParamDefinitionType
  use CharacterStringModule, only: CharacterStringType

  implicit none
  private
  public :: LoadContextType
  public :: ReadStateVarType
  public :: rsv_name
  public :: is_keystring_period

  enum, bind(C)
    enumerator :: LOAD_UNDEF = 0 !< undefined load type
    enumerator :: LIST = 1 !< list load
    enumerator :: LAYERARRAY = 2 !< readasarrays load
    enumerator :: GRIDARRAY = 3 !< readarraygrid load
    enumerator :: KEYSTRING = 4 !< keystring period block load
  end enum

  !> @brief Pointer type for read state variable
  !<
  type ReadStateVarType
    integer(I4B), pointer :: invar
  end type ReadStateVarType

  interface setptr
    module procedure setptr_int, setptr_charstr1d, &
      setptr_auxvar
  end interface setptr

  !> @brief Input load context for generic dynamic loaders and StructArray
  !! based static loads.  Classifies the input, determines in-scope
  !! parameters, and manages memory-manager scalars and array pointers.
  !<
  type :: LoadContextType
    integer(I4B), pointer :: naux => null() !< number of auxiliary variables
    integer(I4B), pointer :: maxbound => null() !< value associated with named_bound
    integer(I4B), pointer :: boundnames => null() !< are bound names optioned
    integer(I4B), pointer :: iprpak => null() !< print input option
    integer(I4B), pointer :: nbound => null() !< number of bounds in period
    integer(I4B), pointer :: ncpl => null() !< ncpl associated with model shape
    integer(I4B), pointer :: nodes => null() !< nodes associated with model shape
    integer(I4B), dimension(:), pointer, contiguous :: mshape => null() !< model shape
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: auxname_cst => null() !< array of auxiliary names
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: boundname_cst => null() !< array of bound names
    real(DP), dimension(:, :), pointer, &
      contiguous :: auxvar => null() !< auxiliary variable array
    integer(I4B) :: loadtype !< load type enum: LIST, LAYERARRAY, GRIDARRAY, KEYSTRING
    logical(LGP) :: set_scalars = .false. !< .true. when dimension scalars must be set
    logical(LGP) :: set_mshape = .false. !< .true. when model shape is load dependency
    logical(LGP) :: is_exchange = .false. !< .true. for exchange contexts
    logical(LGP) :: is_advanced = .false. !< .true. for advanced package KEYSTRING loadtype
    type(InputParamDefinitionType), pointer :: setting_idt => null() !< internal idt for SETTING column
    character(len=LENVARNAME) :: blockname !< load block name
    character(len=LENVARNAME), allocatable :: named_bound(:) !< dimension variable names to sum for maxbound
    integer(I4B) :: nleading = 0 !< count of leading (pre-keystring) columns
    character(len=LINELENGTH), dimension(:), allocatable :: params !< in-scope param tagnames
    integer(I4B), allocatable :: member_nsubs(:) !< nsub per member (1..nmembers)
    type(ModflowInputType) :: mf6_input !< description of input
  contains
    ! --- public interface ---
    procedure :: init
    procedure :: allocate_arrays
    procedure :: allocate_params
    procedure :: rsv_alloc
    procedure :: destroy
    ! --- internal ---
    procedure, private :: resolve_context
    procedure, private :: resolve_loadtype
    procedure, private :: set_params
    procedure, private :: keystring_member_names
    procedure, private :: resolve_dimensions
    procedure, private :: scale_keystring_maxbound
    procedure, private :: allocate_param
    procedure :: check_developmode
    procedure, private :: in_scope
    procedure, private :: option_check
  end type LoadContextType

contains

  !> @brief Initialize the load context.
  !!
  !! Classifies the input, builds the in-scope parameter list,
  !! and resolves memory-manager dimension scalars.
  !<
  subroutine init(this, mf6_input, blockname, named_bound)
    use InputOutputModule, only: upcase
    class(LoadContextType) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), optional, intent(in) :: blockname
    character(len=*), dimension(:), optional, intent(in) :: named_bound
    integer(I4B) :: n

    this%mf6_input = mf6_input

    if (present(blockname)) then
      this%blockname = blockname
      call upcase(this%blockname)
    else
      this%blockname = 'PERIOD'
    end if

    if (present(named_bound)) then
      allocate (this%named_bound(size(named_bound)))
      do n = 1, size(named_bound)
        this%named_bound(n) = named_bound(n)
        call upcase(this%named_bound(n))
      end do
    else
      allocate (this%named_bound(1))
      this%named_bound(1) = 'MAXBOUND'
    end if

    call this%resolve_context()
    call this%resolve_loadtype()
    call this%set_params()
    call this%resolve_dimensions()
  end subroutine init

  !> @brief Set context flags from input load_scope and component metadata.
  !<
  subroutine resolve_context(this)
    class(LoadContextType) :: this

    this%set_scalars = .false.
    this%set_mshape = .false.
    this%is_exchange = .false.

    select case (this%mf6_input%load_scope)
    case ('ROOT')
      ! no memory setup needed for root context
    case ('SIM')
      ! only exchange inputs need scalar setup under SIM scope
      if (this%mf6_input%component_type == 'EXG') then
        this%set_scalars = .true.
        this%is_exchange = .true.
      end if
    case ('MODEL')
      this%set_scalars = .true.
      ! OC and STO are model packages without period block stress data
      if (this%mf6_input%subcomponent_type /= 'OC' .and. &
          this%mf6_input%subcomponent_type /= 'STO') then
        this%set_mshape = .true.
      end if
    case default
      errmsg = 'LoadContext unrecognized load_scope for mempath: '// &
               trim(this%mf6_input%mempath)
      call store_error(errmsg, .true.)
    end select
  end subroutine resolve_context

  !> @brief Determine loadtype from block and param definitions.
  !<
  subroutine resolve_loadtype(this)
    use DefinitionSelectModule, only: idt_default
    class(LoadContextType) :: this
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: n

    this%loadtype = LOAD_UNDEF

    ! detect aggregate (list/keystring) load type
    do n = 1, size(this%mf6_input%block_dfns)
      if (this%mf6_input%block_dfns(n)%blockname == this%blockname) then
        if (this%mf6_input%block_dfns(n)%aggregate) then
          if (this%blockname == 'PERIOD' .and. &
              is_keystring_period(this%mf6_input)) then
            this%loadtype = KEYSTRING
          else
            this%loadtype = LIST
          end if
          exit
        end if
      end if
    end do

    ! check if KEYSTRING package also has a PACKAGEDATA block
    if (this%loadtype == KEYSTRING) then
      do n = 1, size(this%mf6_input%block_dfns)
        if (this%mf6_input%block_dfns(n)%blockname == 'PACKAGEDATA') then
          this%is_advanced = .true.
          this%setting_idt => &
            idt_default(this%mf6_input%component_type, &
                        this%mf6_input%subcomponent_type, &
                        'PERIOD', 'SETTING', 'SETTING', 'STRING')
          exit
        end if
      end do
    end if

    ! detect array-based load
    if (this%loadtype == LOAD_UNDEF) then
      do n = 1, size(this%mf6_input%param_dfns)
        idt => this%mf6_input%param_dfns(n)
        if (idt%blockname == 'OPTIONS') then
          select case (idt%tagname)
          case ('READASARRAYS')
            this%loadtype = LAYERARRAY
          case ('READARRAYGRID')
            this%loadtype = GRIDARRAY
          case default
            ! no-op
          end select
        end if
      end do
    end if
  end subroutine resolve_loadtype

  !> @brief Resolve dimension scalars and scale keystring maxbound.
  !<
  subroutine resolve_dimensions(this)
    use MemoryManagerModule, only: mem_setptr, get_isize
    class(LoadContextType) :: this
    integer(I4B) :: isize
    integer(I4B), pointer :: maxbound_ptr

    if (this%set_scalars) then

      call setptr(this%nbound, 'NBOUND', this%mf6_input%mempath)
      call setval(this%naux, 'NAUX', this%mf6_input%mempath)
      call setval(this%ncpl, 'NCPL', this%mf6_input%mempath)
      call setval(this%nodes, 'NODES', this%mf6_input%mempath)
      call setval(this%boundnames, 'BOUNDNAMES', this%mf6_input%mempath)
      call setval(this%iprpak, 'IPRPAK', this%mf6_input%mempath)

      ! resolve maxbound: sum all named_bound variable values
      allocate (this%maxbound)
      this%maxbound = 0
      call sum_named_bounds(this%named_bound, this%mf6_input%mempath, &
                            this%maxbound)
      ! fallback: try MAXBOUND directly when named_bound tokens yield nothing
      if (this%maxbound == 0) then
        call get_isize('MAXBOUND', this%mf6_input%mempath, isize)
        if (isize > -1) then
          call mem_setptr(maxbound_ptr, 'MAXBOUND', this%mf6_input%mempath)
          this%maxbound = maxbound_ptr
          nullify (maxbound_ptr)
        end if
      end if

      ! reset nbound
      this%nbound = 0
    end if

    if (this%set_mshape .and. &
        this%blockname == 'PERIOD') then
      call mem_setptr(this%mshape, 'MODEL_SHAPE', &
                      this%mf6_input%component_mempath)

      if (this%ncpl == 0) then
        if (size(this%mshape) == 2) then
          this%ncpl = this%mshape(2)
        else if (size(this%mshape) == 3) then
          this%ncpl = this%mshape(2) * this%mshape(3)
        end if
      end if

      if (this%nodes == 0) this%nodes = product(this%mshape)
      if (this%loadtype == KEYSTRING) call this%scale_keystring_maxbound()
    end if
  end subroutine resolve_dimensions

  !> @brief Scale maxbound by the number of KEYSTRING members.
  !!
  !! When a DIMENSIONS block is present (MAW/SFR/LAK/UZF) maxbound is the
  !! feature count; multiply by member count so every feature can use every
  !! setting in one period.  When there is no DIMENSIONS block (TVK/TVS)
  !! fall back to nodes * member count.
  !<
  subroutine scale_keystring_maxbound(this)
    use DefinitionSelectModule, only: get_aggregate_definition_type, &
                                      idt_parse_rectype
    class(LoadContextType) :: this
    type(InputParamDefinitionType), pointer :: aidt, ks_aidt
    character(len=LINELENGTH), allocatable :: cols(:), ks_cols(:)
    integer(I4B) :: nmembers, ncol

    nmembers = 0
    aidt => get_aggregate_definition_type(this%mf6_input%aggregate_dfns, &
                                          this%mf6_input%component_type, &
                                          this%mf6_input%subcomponent_type, &
                                          'PERIOD')
    call idt_parse_rectype(aidt, cols, ncol)
    ks_aidt => find_setting_aggregate(this%mf6_input, cols, ncol)
    if (associated(ks_aidt)) call idt_parse_rectype(ks_aidt, ks_cols, nmembers)
    if (allocated(cols)) deallocate (cols)
    if (allocated(ks_cols)) deallocate (ks_cols)

    if (nmembers > 0) then
      if (this%maxbound == 0) then
        this%maxbound = this%nodes * nmembers
      else
        this%maxbound = this%maxbound * nmembers
      end if
    end if
  end subroutine scale_keystring_maxbound

  !> @brief allocate arrays
  !!
  !! call this routine after input parameters have been allocated,
  !! e.g. after load_params() with create has been called for array
  !! based loaders or after all mem_create_vector() calls have
  !! been made for list based load.
  !!
  !<
  subroutine allocate_arrays(this)
    use MemoryManagerModule, only: mem_allocate, mem_setptr, get_isize
    class(LoadContextType) :: this
    integer(I4B), dimension(:, :), pointer, contiguous :: cellid
    integer(I4B), dimension(:), pointer, contiguous :: nodeulist

    if (this%set_mshape .and. &
        this%blockname == 'PERIOD') then
      ! allocate cellid if this is not list input
      if (this%loadtype == LAYERARRAY .or. &
          this%loadtype == GRIDARRAY) then
        call mem_allocate(cellid, 0, 0, 'CELLID', this%mf6_input%mempath)
      end if

      ! allocate nodeulist for list and layerarray packages only;
      ! keystring and advanced packages do not use a flat nodeulist
      if (this%loadtype /= GRIDARRAY .and. &
          this%loadtype /= KEYSTRING) then
        call mem_allocate(nodeulist, 0, 'NODEULIST', this%mf6_input%mempath)
      end if

      ! set pointers to aux/bound arrays for non-keystring packages;
      ! keystring packages manage aux through struct array columns
      if (this%loadtype /= KEYSTRING) then
        call setptr(this%auxname_cst, 'AUXILIARY', &
                    this%mf6_input%mempath, LENAUXNAME)
        call setptr(this%boundname_cst, 'BOUNDNAME', &
                    this%mf6_input%mempath, LENBOUNDNAME)
        call setptr(this%auxvar, this%mf6_input%mempath)
      end if

    else if (this%is_exchange) then
      ! set pointers to arrays
      call setptr(this%auxname_cst, 'AUXILIARY', &
                  this%mf6_input%mempath, LENAUXNAME)
      call setptr(this%boundname_cst, 'BOUNDNAME', &
                  this%mf6_input%mempath, LENBOUNDNAME)
      call setptr(this%auxvar, this%mf6_input%mempath)
    end if
  end subroutine allocate_arrays

  !> @brief allocate a package dynamic input parameter
  !<
  subroutine allocate_param(this, idt)
    use InputDefinitionModule, only: InputParamDefinitionType
    class(LoadContextType) :: this
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: dimsize

    ! initialize
    dimsize = 0

    if (this%loadtype == LAYERARRAY .or. &
        this%loadtype == GRIDARRAY) then
      select case (idt%shape)
      case ('NCPL', 'NAUX NCPL')
        dimsize = this%ncpl
      case ('NODES', 'NAUX NODES')
        dimsize = this%maxbound
      case default
      end select
    end if

    select case (idt%datatype)
    case ('INTEGER')
      if (this%loadtype == LIST) then
        call allocate_int1d(this%maxbound, idt%mf6varname, &
                            this%mf6_input%mempath)
      end if
    case ('DOUBLE')
      if (this%loadtype == LIST) then
        call allocate_dbl1d(this%maxbound, idt%mf6varname, &
                            this%mf6_input%mempath)
      end if
    case ('STRING')
      if (this%loadtype == LIST) then
        call allocate_charstr1d(LENBOUNDNAME, this%maxbound, idt%mf6varname, &
                                this%mf6_input%mempath)
      end if
    case ('INTEGER1D')
      if (this%loadtype == LIST) then
        if (idt%shape == 'NCELLDIM') then
          call allocate_int2d(size(this%mshape), this%maxbound, &
                              idt%mf6varname, this%mf6_input%mempath)
        end if
      else if (this%loadtype == LAYERARRAY .or. &
               this%loadtype == GRIDARRAY) then
        call allocate_int1d(dimsize, idt%mf6varname, &
                            this%mf6_input%mempath)
      end if
    case ('DOUBLE1D')
      if (idt%shape == 'NAUX') then
        call allocate_dbl2d(this%naux, this%maxbound, &
                            idt%mf6varname, this%mf6_input%mempath)
      else if (this%loadtype == LAYERARRAY .or. &
               this%loadtype == GRIDARRAY) then
        call allocate_dbl1d(dimsize, idt%mf6varname, &
                            this%mf6_input%mempath)
      end if
    case ('DOUBLE2D')
      if (this%loadtype == LAYERARRAY .or. &
          this%loadtype == GRIDARRAY) then
        call allocate_dbl2d(this%naux, dimsize, idt%mf6varname, &
                            this%mf6_input%mempath)
      end if
    case default
    end select
  end subroutine allocate_param

  !> @brief Allocate each in-scope parameter in the memory manager.
  !!
  !! Call after init() for array-based loaders (LAYERARRAY, GRIDARRAY) that
  !! need memory-manager storage allocated for every in-scope parameter.
  !! Loaders access this%params and size(this%params) directly.
  !<
  subroutine allocate_params(this)
    use DefinitionSelectModule, only: get_param_definition_type
    class(LoadContextType) :: this
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: n
    do n = 1, size(this%params)
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       this%blockname, this%params(n), '')
      call this%allocate_param(idt)
    end do
  end subroutine allocate_params

  !> @brief Return .true. if an optional parameter is active for this load.
  !!
  !! Required and structural (KEYSTRING/RECARRAY/RECORD) params are handled
  !! by the caller; this routine only evaluates optional leaf params.
  !!
  !! Generic conditions (AUX, BOUNDNAME, readarray indicator) are checked
  !! first.  Package-specific conditions follow via a select on
  !! subcomponent_type; unrecognized types return .false. (conservative).
  !<
  function in_scope(this, tagname)
    use DefinitionSelectModule, only: get_param_definition_type, idt_datatype
    class(LoadContextType) :: this
    character(len=*), intent(in) :: tagname
    logical(LGP) :: in_scope
    type(InputParamDefinitionType), pointer :: idt
    character(len=LINELENGTH) :: datatype

    idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                     this%mf6_input%component_type, &
                                     this%mf6_input%subcomponent_type, &
                                     this%blockname, tagname, '')
    ! required params are always in scope
    if (idt%required) then
      in_scope = .true.
      return
    else
      in_scope = .false.
    end if

    ! structural container types are never loaded as leaf params
    datatype = idt_datatype(idt)
    if (datatype == 'KEYSTRING' .or. &
        datatype == 'RECARRAY' .or. &
        datatype == 'RECORD') return

    ! --- generic conditions ---
    if (tagname == 'AUXVAR' .or. tagname == 'AUX') then
      in_scope = this%option_check('NAUX', 0)
      return
    end if

    if (tagname == 'BOUNDNAME') then
      in_scope = this%option_check('BOUNDNAMES', 0)
      return
    end if

    ! readarray indicator variable (e.g. IRCH, IRET): in scope for LAYERARRAY only
    if (tagname == 'I'//trim(this%mf6_input%subcomponent_type(1:3))) then
      in_scope = (this%loadtype == LAYERARRAY)
      return
    end if

    ! --- package-specific conditions ---
    select case (this%mf6_input%subcomponent_type)
    case ('EVT')
      if (tagname == 'PXDP' .or. tagname == 'PETM') then
        in_scope = this%option_check('NSEG', 1)
      else if (tagname == 'PETM0') then
        in_scope = this%option_check('SURFRATESPEC', 0)
      end if
    case ('MVR', 'MVT', 'MVE')
      if (tagname == 'MNAME' .or. &
          tagname == 'MNAME1' .or. &
          tagname == 'MNAME2') then
        in_scope = this%option_check('MODELNAMES', 0)
      end if
    case ('NAM')
      in_scope = .true.
    case ('SFR')
      if (tagname == 'IC') in_scope = .true.
    case ('SSM')
      if (tagname == 'MIXED') in_scope = .true.
    case ('SPC', 'SPCA')
      in_scope = .true.
    case default
      if (idt%in_record) then
        in_scope = .true.
      else
        errmsg = 'LoadContext in_scope needs new check for: '// &
                 trim(this%mf6_input%subcomponent_type)//'/'//trim(tagname)
        call store_error(errmsg, .true.)
      end if
    end select
  end function in_scope

  !> @brief Return .true. if a memory-manager integer option variable exceeds a threshold.
  !<
  function option_check(this, varname, threshold)
    use MemoryManagerModule, only: get_isize, mem_setptr
    class(LoadContextType) :: this
    character(len=*), intent(in) :: varname
    integer(I4B), intent(in) :: threshold
    logical(LGP) :: option_check
    integer(I4B) :: isize
    integer(I4B), pointer :: intptr
    option_check = .false.
    call get_isize(varname, this%mf6_input%mempath, isize)
    if (isize > 0) then
      call mem_setptr(intptr, varname, this%mf6_input%mempath)
      if (intptr > threshold) option_check = .true.
    end if
  end function option_check

  !> @brief set set of in scope parameters for package
  !<
  subroutine set_params(this)
    use ArrayHandlersModule, only: expandarray
    use DefinitionSelectModule, only: get_param_definition_type, &
                                      get_aggregate_definition_type, &
                                      idt_parse_rectype
    class(LoadContextType) :: this
    type(InputParamDefinitionType), pointer :: idt, aidt
    character(len=LINELENGTH), dimension(:), allocatable :: param_buf
    character(len=LINELENGTH), dimension(:), allocatable :: cols
    character(len=LINELENGTH), allocatable :: member_names(:)
    integer(I4B), allocatable :: member_nsubs(:)
    integer(I4B) :: keepcnt, iparam, nparam, nmembers, n
    logical(LGP) :: keep, tag_found

    ! initialize
    keepcnt = 0

    if (this%loadtype == LIST .or. &
        this%loadtype == KEYSTRING) then
      ! get aggregate param definition for period block
      aidt => &
        get_aggregate_definition_type(this%mf6_input%aggregate_dfns, &
                                      this%mf6_input%component_type, &
                                      this%mf6_input%subcomponent_type, &
                                      this%blockname)
      ! split recarray definition
      call idt_parse_rectype(aidt, cols, nparam)
    else
      nparam = size(this%mf6_input%param_dfns)
    end if

    ! allocate dfn input params
    do iparam = 1, nparam
      if (this%loadtype == LIST .or. &
          this%loadtype == KEYSTRING) then
        ! use found so keystring placeholders are silently skipped
        idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                         this%mf6_input%component_type, &
                                         this%mf6_input%subcomponent_type, &
                                         this%blockname, cols(iparam), '', &
                                         found=tag_found)
      else
        tag_found = .true.
        idt => this%mf6_input%param_dfns(iparam)
      end if

      if (.not. tag_found) then
        keep = .false.
      else if (idt%blockname /= this%blockname) then
        keep = .false.
      else
        keep = this%in_scope(idt%tagname)
      end if

      if (keep) then
        keepcnt = keepcnt + 1
        call expandarray(param_buf)
        param_buf(keepcnt) = trim(idt%tagname)
      end if
    end do

    ! record leading-column count before member expansion
    if (this%loadtype == LIST .or. &
        this%loadtype == KEYSTRING) this%nleading = keepcnt

    ! for keystring packages: append member names and store associated metadata
    if (this%loadtype == KEYSTRING) then
      call this%keystring_member_names(member_names, member_nsubs, nmembers)
      do n = 1, nmembers
        keepcnt = keepcnt + 1
        call expandarray(param_buf)
        param_buf(keepcnt) = trim(member_names(n))
      end do
      if (allocated(this%member_nsubs)) deallocate (this%member_nsubs)
      if (nmembers > 0) then
        allocate (this%member_nsubs(nmembers))
        this%member_nsubs = member_nsubs
      end if
    end if

    ! update nparam to total (leading + members)
    nparam = keepcnt

    ! allocate and fill params
    allocate (this%params(nparam))
    do iparam = 1, nparam
      this%params(iparam) = trim(param_buf(iparam))
    end do

    ! cleanup
    if (allocated(param_buf)) deallocate (param_buf)
  end subroutine set_params

  !> @brief allocate a read state variable
  !!
  !! Create and set a read state variable, e.g. 'INRECHARGE',
  !! which are updated per iper load as follows:
  !! -1: unset, not in use
  !!  0: not read in most recent period block
  !!  1: numeric input read in most recent period block
  !!  2: time series input read in most recent period block
  !!
  !<
  function rsv_alloc(this, mf6varname) result(varname)
    use ConstantsModule, only: LENVARNAME
    use MemoryManagerModule, only: mem_setptr, mem_allocate
    class(LoadContextType) :: this
    character(len=*), intent(in) :: mf6varname
    character(len=LENVARNAME) :: varname
    integer(I4B), pointer :: intvar
    varname = rsv_name(mf6varname)
    call mem_allocate(intvar, varname, this%mf6_input%mempath)
    intvar = -1
  end function rsv_alloc

  !> @brief destroy input context object
  !<
  subroutine destroy(this)
    class(LoadContextType) :: this

    if (allocated(this%named_bound)) deallocate (this%named_bound)
    if (allocated(this%member_nsubs)) deallocate (this%member_nsubs)
    if (associated(this%setting_idt)) then
      deallocate (this%setting_idt)
      nullify (this%setting_idt)
    end if

    if (this%set_scalars) then
      ! deallocate local
      deallocate (this%naux)
      deallocate (this%ncpl)
      deallocate (this%nodes)
      deallocate (this%maxbound)
      deallocate (this%boundnames)
      deallocate (this%iprpak)
    end if

    ! nullify
    nullify (this%naux)
    nullify (this%nbound)
    nullify (this%ncpl)
    nullify (this%nodes)
    nullify (this%maxbound)
    nullify (this%boundnames)
    nullify (this%iprpak)
    nullify (this%auxname_cst)
    nullify (this%boundname_cst)
    nullify (this%auxvar)
    nullify (this%mshape)
  end subroutine destroy

  !> @brief Return the KEYSTRING aggregate for the SETTING token in rec_cols, or null().
  !<
  function find_setting_aggregate(mf6_input, rec_cols, nrec_col) result(ks_aidt)
    use InputOutputModule, only: upcase
    use DefinitionSelectModule, only: idt_datatype
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=LINELENGTH), intent(in) :: rec_cols(:)
    integer(I4B), intent(in) :: nrec_col
    type(InputParamDefinitionType), pointer :: ks_aidt
    character(len=LINELENGTH) :: token, tagname
    integer(I4B) :: m, n, ilen
    ks_aidt => null()
    do m = 1, nrec_col
      token = trim(rec_cols(m))
      call upcase(token)
      ilen = len_trim(token)
      ! minimum 8 chars: a valid XSETTING token is at least 1 char prefix + 7 for 'SETTING'
      ! if (ilen < 8) cycle
      if (token(ilen - 6:ilen) /= 'SETTING') cycle
      do n = 1, size(mf6_input%aggregate_dfns)
        tagname = mf6_input%aggregate_dfns(n)%tagname
        call upcase(tagname)
        if (trim(tagname) == trim(token)) then
          ks_aidt => mf6_input%aggregate_dfns(n)
          if (idt_datatype(ks_aidt) /= 'KEYSTRING') ks_aidt => null()
          exit
        end if
      end do
      exit
    end do
  end function find_setting_aggregate

  !> @brief Append sub-member column names from a RECORD compound entry to member_names.
  !<
  subroutine expand_record_submembers(mf6_input, rec_idt, member_names, nmembers)
    use InputOutputModule, only: upcase
    use ArrayHandlersModule, only: expandarray
    use DefinitionSelectModule, only: idt_parse_rectype, idt_datatype
    type(ModflowInputType), intent(in) :: mf6_input
    type(InputParamDefinitionType), pointer, intent(in) :: rec_idt
    character(len=LINELENGTH), allocatable, intent(inout) :: member_names(:)
    integer(I4B), intent(inout) :: nmembers
    type(InputParamDefinitionType), pointer :: sub_idt
    character(len=LINELENGTH), allocatable :: sub_cols(:)
    character(len=LINELENGTH) :: token, tagname
    integer(I4B) :: k, j, nsub_col
    call idt_parse_rectype(rec_idt, sub_cols, nsub_col)
    do k = 1, nsub_col
      token = trim(sub_cols(k))
      call upcase(token)
      do j = 1, size(mf6_input%param_dfns)
        sub_idt => mf6_input%param_dfns(j)
        if (sub_idt%blockname /= 'PERIOD') cycle
        tagname = sub_idt%tagname
        call upcase(tagname)
        if (trim(tagname) /= trim(token)) cycle
        if (idt_datatype(sub_idt) == 'RECORD') cycle
        nmembers = nmembers + 1
        call expandarray(member_names)
        member_names(nmembers) = trim(sub_idt%tagname)
        exit
      end do
    end do
    if (allocated(sub_cols)) deallocate (sub_cols)
  end subroutine expand_record_submembers

  !> @brief Return .true. if mf6_input's PERIOD block uses keystring dispatch.
  !<
  function is_keystring_period(mf6_input) result(res)
    use DefinitionSelectModule, only: get_aggregate_definition_type, &
                                      idt_parse_rectype
    type(ModflowInputType), intent(in) :: mf6_input
    logical(LGP) :: res, has_period
    type(InputParamDefinitionType), pointer :: aidt, ks_aidt
    character(len=LINELENGTH), allocatable :: cols(:)
    integer(I4B) :: n, ncol
    res = .false.
    has_period = .false.
    do n = 1, size(mf6_input%block_dfns)
      if (mf6_input%block_dfns(n)%blockname == 'PERIOD') then
        has_period = .true.
      end if
    end do
    if (.not. has_period) return
    aidt => get_aggregate_definition_type(mf6_input%aggregate_dfns, &
                                          mf6_input%component_type, &
                                          mf6_input%subcomponent_type, &
                                          'PERIOD')
    call idt_parse_rectype(aidt, cols, ncol)
    if (ncol >= 2) then
      ks_aidt => find_setting_aggregate(mf6_input, cols, ncol)
      if (associated(ks_aidt)) res = .true.
    end if
    if (allocated(cols)) deallocate (cols)
  end function is_keystring_period

  !> @brief Return keystring member column names and nsub counts.
  !!
  !! Private helper called from set_params. Results are returned via
  !! output parameters; the caller is responsible for storing them.
  !! Column order follows the KEYSTRING aggregate definition token list.
  !! For each token in the aggregate:
  !!   - RECORD compound group: sub-members expanded in RECORD order;
  !!     first entry (KEYWORD dispatch header) gets nsub = sub-member count,
  !!     remaining entries get nsub = 0.
  !!   - direct-dispatch param: appended with nsub = 0.
  !<
  subroutine keystring_member_names(this, member_names, member_nsubs, nmembers)
    use InputOutputModule, only: upcase
    use ArrayHandlersModule, only: expandarray
    use DefinitionSelectModule, only: idt_parse_rectype, idt_datatype, &
                                      get_aggregate_definition_type
    class(LoadContextType) :: this
    character(len=LINELENGTH), allocatable, intent(out) :: member_names(:)
    integer(I4B), allocatable, intent(out) :: member_nsubs(:)
    integer(I4B), intent(out) :: nmembers
    type(InputParamDefinitionType), pointer :: aidt, ks_aidt, idt
    character(len=LINELENGTH), allocatable :: rec_cols(:), ks_cols(:)
    character(len=LINELENGTH) :: rec_token, tagname
    integer(I4B) :: m, n, nrec_col, nks_col, nmembers0, k

    nmembers = 0

    ! get RECARRAY aggregate for period block and parse its column tokens
    aidt => get_aggregate_definition_type(this%mf6_input%aggregate_dfns, &
                                          this%mf6_input%component_type, &
                                          this%mf6_input%subcomponent_type, &
                                          this%blockname)
    call idt_parse_rectype(aidt, rec_cols, nrec_col)

    ! find the KEYSTRING aggregate for the SETTING token
    ks_aidt => find_setting_aggregate(this%mf6_input, rec_cols, nrec_col)
    if (allocated(rec_cols)) deallocate (rec_cols)
    if (.not. associated(ks_aidt)) return

    ! parse the KEYSTRING aggregate to get member token list — canonical order
    call idt_parse_rectype(ks_aidt, ks_cols, nks_col)

    ! walk the keystring token list in aggregate order
    do m = 1, nks_col
      rec_token = trim(ks_cols(m))
      call upcase(rec_token)

      ! locate matching param_dfns entry for this token
      do n = 1, size(this%mf6_input%param_dfns)
        if (this%mf6_input%param_dfns(n)%blockname /= 'PERIOD') cycle
        tagname = this%mf6_input%param_dfns(n)%tagname
        call upcase(tagname)
        if (trim(tagname) /= trim(rec_token)) cycle

        idt => this%mf6_input%param_dfns(n)
        if (idt_datatype(idt) == 'RECORD') then
          ! compound group: expand sub-members in RECORD type order
          nmembers0 = nmembers
          call expand_record_submembers(this%mf6_input, idt, member_names, &
                                        nmembers)
          ! first added entry is the KEYWORD header; remaining are sub-members
          do k = nmembers0 + 1, nmembers
            call expandarray(member_nsubs)
            if (k == nmembers0 + 1) then
              member_nsubs(k) = nmembers - nmembers0 - 1
            else
              member_nsubs(k) = 0
            end if
          end do
        else
          ! direct-dispatch param
          nmembers = nmembers + 1
          call expandarray(member_names)
          call expandarray(member_nsubs)
          member_names(nmembers) = trim(this%mf6_input%param_dfns(n)%tagname)
          member_nsubs(nmembers) = 0
        end if
        exit
      end do
    end do

    if (allocated(ks_cols)) deallocate (ks_cols)
  end subroutine keystring_member_names

  !> @brief Check whether any in-scope parameter is a development-mode feature.
  !<
  subroutine check_developmode(this, input_name)
    use FeatureFlagsModule, only: developmode
    use SimVariablesModule, only: iout
    use DefinitionSelectModule, only: get_param_definition_type
    class(LoadContextType) :: this
    character(len=*), intent(in) :: input_name
    type(InputParamDefinitionType), pointer :: idt
    character(len=LINELENGTH) :: dev_msg
    integer(I4B) :: n

    do n = 1, size(this%params)
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       this%blockname, this%params(n), '')
      if (idt%developmode) then
        dev_msg = 'Input tag "'//trim(idt%tagname)// &
          &'" read from file "'//trim(input_name)// &
          &'" is still under development. Install the &
          &nightly build or compile from source with IDEVELOPMODE = 1.'
        call developmode(dev_msg, iout)
      end if
    end do
  end subroutine check_developmode

  !> @brief create read state variable name
  !<
  function rsv_name(mf6varname) result(varname)
    use ConstantsModule, only: LENVARNAME
    character(len=*), intent(in) :: mf6varname
    character(len=LENVARNAME) :: varname
    integer(I4B) :: ilen
    character(len=2) :: prefix = 'IN'
    ilen = len_trim(mf6varname)
    if (ilen > (LENVARNAME - len(prefix))) then
      varname = prefix//mf6varname(1:(LENVARNAME - len(prefix)))
    else
      varname = prefix//trim(mf6varname)
    end if
  end function rsv_name

  !> @brief allocate character string type array
  !<
  subroutine allocate_charstr1d(strlen, nrow, varname, mempath)
    use MemoryManagerModule, only: mem_allocate
    integer(I4B), intent(in) :: strlen !< string number of characters
    integer(I4B), intent(in) :: nrow !< integer array number of rows
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable mempath
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: charstr1d
    integer(I4B) :: n
    call mem_allocate(charstr1d, strlen, nrow, varname, mempath)
    do n = 1, nrow
      charstr1d(n) = ''
    end do
  end subroutine allocate_charstr1d

  !> @brief allocate int1d
  !<
  subroutine allocate_int1d(nrow, varname, mempath)
    use MemoryManagerModule, only: mem_allocate
    integer(I4B), intent(in) :: nrow !< integer array number of rows
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable mempath
    integer(I4B), dimension(:), pointer, contiguous :: int1d
    integer(I4B) :: n
    call mem_allocate(int1d, nrow, varname, mempath)
    do n = 1, nrow
      int1d(n) = IZERO
    end do
  end subroutine allocate_int1d

  !> @brief allocate int2d
  !<
  subroutine allocate_int2d(ncol, nrow, varname, mempath)
    use MemoryManagerModule, only: mem_allocate
    integer(I4B), intent(in) :: ncol !< integer array number of cols
    integer(I4B), intent(in) :: nrow !< integer array number of rows
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable mempath
    integer(I4B), dimension(:, :), pointer, contiguous :: int2d
    integer(I4B) :: n, m
    call mem_allocate(int2d, ncol, nrow, varname, mempath)
    do m = 1, nrow
      do n = 1, ncol
        int2d(n, m) = IZERO
      end do
    end do
  end subroutine allocate_int2d

  !> @brief allocate dbl1d
  !<
  subroutine allocate_dbl1d(nrow, varname, mempath)
    use MemoryManagerModule, only: mem_allocate
    integer(I4B), intent(in) :: nrow !< integer array number of rows
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable mempath
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    integer(I4B) :: n
    call mem_allocate(dbl1d, nrow, varname, mempath)
    do n = 1, nrow
      dbl1d(n) = DZERO
    end do
  end subroutine allocate_dbl1d

  !> @brief allocate dbl2d
  !<
  subroutine allocate_dbl2d(ncol, nrow, varname, mempath)
    use MemoryManagerModule, only: mem_allocate
    integer(I4B), intent(in) :: ncol !< integer array number of cols
    integer(I4B), intent(in) :: nrow !< integer array number of rows
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mempath !< variable mempath
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B) :: n, m
    call mem_allocate(dbl2d, ncol, nrow, varname, mempath)
    do m = 1, nrow
      do n = 1, ncol
        dbl2d(n, m) = DZERO
      end do
    end do
  end subroutine allocate_dbl2d

  !> @brief sum named dimension variables from mempath
  !!
  !! Loops over each name in named_bound and accumulates its value
  !! from mempath into total.  Variables not present in mempath are
  !! silently skipped.
  !!
  !<
  subroutine sum_named_bounds(named_bound, mempath, total)
    use MemoryManagerModule, only: mem_setptr, get_isize
    character(len=*), dimension(:), intent(in) :: named_bound
    character(len=*), intent(in) :: mempath
    integer(I4B), intent(inout) :: total
    integer(I4B), pointer :: dimptr
    integer(I4B) :: n, isize

    do n = 1, size(named_bound)
      call get_isize(trim(named_bound(n)), mempath, isize)
      if (isize > -1) then
        call mem_setptr(dimptr, trim(named_bound(n)), mempath)
        total = total + dimptr
        nullify (dimptr)
      end if
    end do
  end subroutine sum_named_bounds

  !> @brief allocate intptr and update from input context
  !!
  !<
  subroutine setval(intptr, varname, mempath)
    use MemoryManagerExtModule, only: mem_set_value
    integer(I4B), pointer, intent(inout) :: intptr
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: mempath
    logical(LGP) :: found
    allocate (intptr)
    intptr = 0
    call mem_set_value(intptr, varname, mempath, found, release=.false.)
  end subroutine setval

  !> @brief set intptr to varname
  !!
  !<
  subroutine setptr_int(intptr, varname, mempath)
    use MemoryManagerModule, only: mem_allocate, mem_setptr, get_isize
    integer(I4B), pointer, intent(inout) :: intptr
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: mempath
    integer(I4B) :: isize
    call get_isize(varname, mempath, isize)
    if (isize > -1) then
      call mem_setptr(intptr, varname, mempath)
    else
      call mem_allocate(intptr, varname, mempath)
      intptr = 0
    end if
  end subroutine setptr_int

  !> @brief set charstr1d pointer to varname
  !<
  subroutine setptr_charstr1d(charstr1d, varname, mempath, strlen)
    use MemoryManagerModule, only: mem_allocate, mem_setptr, get_isize
    type(CharacterStringType), dimension(:), pointer, &
      contiguous, intent(inout) :: charstr1d
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: mempath
    integer(I4B), intent(in) :: strlen
    integer(I4B) :: isize
    call get_isize(varname, mempath, isize)
    if (isize > -1) then
      call mem_setptr(charstr1d, varname, mempath)
    else
      call mem_allocate(charstr1d, strlen, 0, varname, mempath)
    end if
  end subroutine setptr_charstr1d

  !> @brief set auxvar pointer
  !!
  !<
  subroutine setptr_auxvar(auxvar, mempath)
    use MemoryManagerModule, only: mem_allocate, mem_setptr, get_isize
    real(DP), dimension(:, :), pointer, &
      contiguous, intent(inout) :: auxvar
    character(len=*), intent(in) :: mempath
    integer(I4B) :: isize
    call get_isize('AUXVAR', mempath, isize)
    if (isize > -1) then
      call mem_setptr(auxvar, 'AUXVAR', mempath)
    else
      call mem_allocate(auxvar, 0, 0, 'AUXVAR', mempath)
    end if
  end subroutine setptr_auxvar

end module LoadContextModule
