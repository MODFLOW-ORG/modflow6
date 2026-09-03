!> @brief Load context for IDM generic dynamic loaders.
!!
!! LoadContextType classifies each input, builds the in-scope parameter
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
  public :: is_advanced_package
  public :: has_dimensions_block
  public :: is_cellid_addressed

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
    logical(LGP) :: is_advanced = .false. !< .true. for advanced package (PACKAGEDATA-paired) KEYSTRING loadtype
    logical(LGP) :: is_dimensions_scoped = .false. !< .true. for DIMENSIONS-block-paired (e.g. SPC) KEYSTRING loadtype
    logical(LGP) :: is_cellid_scoped = .false. !< .true. for CELLID-addressed (e.g. TVK/TVS) KEYSTRING loadtype
    logical(LGP) :: is_feature_indexed = .false. !< .true. when is_advanced .or. is_dimensions_scoped
    logical(LGP) :: has_setting_dispatch = .false. !< .true. when is_feature_indexed .or. is_cellid_scoped
    type(InputParamDefinitionType), pointer :: setting_idt => null() !< internal idt for SETTING column
    character(len=LENVARNAME) :: blockname !< load block name
    character(len=LENVARNAME) :: named_bound !< name of dimension variable for maxbound; defaults to MAXBOUND
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
    character(len=*), optional, intent(in) :: named_bound

    this%mf6_input = mf6_input

    if (present(blockname)) then
      this%blockname = blockname
      call upcase(this%blockname)
    else
      this%blockname = 'PERIOD'
    end if

    if (present(named_bound)) then
      this%named_bound = named_bound
      call upcase(this%named_bound)
    else
      this%named_bound = 'MAXBOUND'
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

    ! mutually exclusive KEYSTRING package subtypes
    if (this%loadtype == KEYSTRING) then
      this%is_advanced = is_advanced_package(this%mf6_input)
      this%is_dimensions_scoped = has_dimensions_block(this%mf6_input)
      this%is_cellid_scoped = is_cellid_addressed(this%mf6_input)
      this%is_feature_indexed = this%is_advanced .or. this%is_dimensions_scoped
      ! is_cellid_scoped has no exclusion guard of its own; without this
      ! check, both flags true would silently double-allocate a permanent
      ! array in Mf6FileKeystring.f90
      if (this%is_feature_indexed .and. this%is_cellid_scoped) then
        errmsg = 'LoadContext: is_feature_indexed and is_cellid_scoped &
                 &cannot both be true for mempath: '// &
                 trim(this%mf6_input%mempath)
        call store_error(errmsg, .true.)
      end if
      this%has_setting_dispatch = &
        this%is_feature_indexed .or. this%is_cellid_scoped
      if (this%has_setting_dispatch) then
        this%setting_idt => &
          idt_default(this%mf6_input%component_type, &
                      this%mf6_input%subcomponent_type, &
                      'PERIOD', 'SETTING', 'SETTING', 'STRING')
      end if
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
    use MemoryManagerModule, only: mem_setptr
    class(LoadContextType) :: this

    if (this%set_scalars) then

      call setptr(this%nbound, 'NBOUND', this%mf6_input%mempath)
      call setval(this%naux, 'NAUX', this%mf6_input%mempath)
      call setval(this%ncpl, 'NCPL', this%mf6_input%mempath)
      call setval(this%nodes, 'NODES', this%mf6_input%mempath)
      call setval(this%boundnames, 'BOUNDNAMES', this%mf6_input%mempath)
      call setval(this%iprpak, 'IPRPAK', this%mf6_input%mempath)
      call setval(this%maxbound, this%named_bound, this%mf6_input%mempath)

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

  !> @brief Scale maxbound (a feature or node count) by the number of
  !! KEYSTRING members, so every feature can use every setting in one
  !! period.
  !<
  subroutine scale_keystring_maxbound(this)
    class(LoadContextType) :: this
    integer(I4B) :: nmembers

    nmembers = 0
    if (allocated(this%member_nsubs)) nmembers = size(this%member_nsubs)

    if (nmembers > 0) then
      if (this%maxbound == 0) then
        if (.not. this%is_feature_indexed) then
          this%maxbound = this%nodes * nmembers
        end if
        ! else: feature-indexed packages with a genuinely zero count stay
        ! at 0 rather than falling back to node count
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
  !!
  !! Called only from allocate_params(), itself called only by the
  !! LAYERARRAY/GRIDARRAY loaders -- this%loadtype is never LIST here.
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
    case ('INTEGER1D')
      if (this%loadtype == LAYERARRAY .or. &
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
    case ('SSM')
      if (tagname == 'MIXED') in_scope = .true.
    case ('SPC', 'SPCA')
      in_scope = .true.
    case ('LAK', 'MAW', 'SFR')
      in_scope = .true.
    case default
      ! Unrecognized subcomponent with an optional param not handled above.
      ! This is a development error — abort with message so developer knows
      ! which package needs a new case.
      errmsg = 'LoadContext in_scope needs new case for: '// &
               trim(this%mf6_input%subcomponent_type)//'/'//trim(tagname)
      call store_error(errmsg, .true.)
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
      if (ilen < 8) cycle
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

  !> @brief Return .true. if mf6_input is an advanced package: a keystring
  !! PERIOD dispatch paired with a PACKAGEDATA block.
  !<
  function is_advanced_package(mf6_input) result(res)
    type(ModflowInputType), intent(in) :: mf6_input
    logical(LGP) :: res
    integer(I4B) :: n
    res = .false.
    if (.not. is_keystring_period(mf6_input)) return
    do n = 1, size(mf6_input%block_dfns)
      if (mf6_input%block_dfns(n)%blockname == 'PACKAGEDATA') then
        res = .true.
        exit
      end if
    end do
  end function is_advanced_package

  !> @brief Return .true. if mf6_input is a keystring PERIOD dispatch paired
  !! with a DIMENSIONS block (e.g. SPC) rather than PACKAGEDATA.
  !!
  !! Only checks that the block exists, not that it declares a usable
  !! feature-count field -- that resolution is KeystringLoadType%ainit's
  !! named_bound (first DIMENSIONS parameter found, whatever its name).
  !<
  function has_dimensions_block(mf6_input) result(res)
    type(ModflowInputType), intent(in) :: mf6_input
    logical(LGP) :: res
    integer(I4B) :: n
    logical(LGP) :: has_dimensions
    res = .false.
    if (.not. is_keystring_period(mf6_input)) return
    has_dimensions = .false.
    do n = 1, size(mf6_input%block_dfns)
      if (mf6_input%block_dfns(n)%blockname == 'DIMENSIONS') then
        has_dimensions = .true.
      end if
      ! PACKAGEDATA-paired (e.g. LAK/MAW/SFR/UZF) takes precedence over a
      ! coincidental DIMENSIONS block (e.g. NOUTLETS, NTABLES)
      if (mf6_input%block_dfns(n)%blockname == 'PACKAGEDATA') return
    end do
    res = has_dimensions
  end function has_dimensions_block

  !> @brief Return .true. if mf6_input is a keystring PERIOD dispatch whose
  !! leading column is CELLID (e.g. TVK/TVS) rather than a stable integer
  !! feature number.
  !<
  function is_cellid_addressed(mf6_input) result(res)
    type(ModflowInputType), intent(in) :: mf6_input
    logical(LGP) :: res
    integer(I4B) :: n
    res = .false.
    if (.not. is_keystring_period(mf6_input)) return
    do n = 1, size(mf6_input%param_dfns)
      if (mf6_input%param_dfns(n)%blockname == 'PERIOD' .and. &
          mf6_input%param_dfns(n)%tagname == 'CELLID') then
        res = .true.
        exit
      end if
    end do
  end function is_cellid_addressed

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
