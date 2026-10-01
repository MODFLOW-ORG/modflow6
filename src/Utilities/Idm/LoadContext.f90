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
                             LENVARNAME, LENBOUNDNAME, DNODATA
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use ModflowInputModule, only: ModflowInputType
  use InputDefinitionModule, only: InputParamDefinitionType
  use CharacterStringModule, only: CharacterStringType
  use IdmDfnSelectorModule, only: idm_is_advanced

  implicit none
  private
  public :: LoadContextType
  public :: ReadStateVarType
  public :: rsv_name
  public :: is_keystring_period
  public :: is_feature_keystring
  public :: is_advanced

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

  ! addressing modes for an applied keystring setting
  integer(I4B), parameter, public :: ADDR_NONE = 0 !< not an applied setting
  integer(I4B), parameter, public :: ADDR_FEATURE = 1 !< addressed by leading id (IFNO/BNDNO)
  integer(I4B), parameter, public :: ADDR_NODE = 2 !< addressed by CELLID
  integer(I4B), parameter, public :: ADDR_SUBINDEX = 3 !< addressed by (id, subindex)

  !> @brief One descriptor per expanded keystring item, so loaders iterate
  !! uniformly instead of re-deriving per-mode state.
  !<
  type, public :: KeystringItemType
    type(InputParamDefinitionType), pointer :: idt => null() !< input definition
    integer(I4B) :: sa_icol = 0 !< struct-array column (SETTING offset resolved)
    logical(LGP) :: is_body = .false. !< RECORD body (non-head)
    integer(I4B) :: head_nbody = 0 !< if a head: number of body members; else 0
    character(len=LENVARNAME) :: head_setting_varname = '' !< body: owning head's SETTING keyword (dispatch key); '' otherwise
    logical(LGP) :: idm_managed = .false. !< IDM allocates+applies a permanent array
    integer(I4B) :: addr_mode = ADDR_NONE !< how a row resolves to an index
    real(DP) :: init_value = DZERO !< DZERO (feature/SPC) | DNODATA (node)
    integer(I4B) :: nfeatures = 0 !< permanent-array length (per-item; LAK lake/outlet)
  end type KeystringItemType

  !> @brief Input load context for generic dynamic loaders and StructArray
  !! based static loads.  Classifies the input, determines in-scope
  !! parameters, and manages memory-manager scalars and array pointers.
  !<
  type :: LoadContextType
    ! memory-manager dimension/option scalars and arrays
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
    ! classification
    integer(I4B) :: loadtype !< load type enum: LIST, LAYERARRAY, GRIDARRAY, KEYSTRING
    logical(LGP) :: set_scalars = .false. !< .true. when dimension scalars must be set
    logical(LGP) :: set_mshape = .false. !< .true. when model shape is load dependency
    logical(LGP) :: is_exchange = .false. !< .true. for exchange contexts
    logical(LGP) :: is_advanced = .false. !< .true. for an advanced package, any loadtype
    ! input handle and in-scope params
    character(len=LENVARNAME) :: blockname !< load block name
    character(len=LENVARNAME) :: named_bound !< name of dimension variable for maxbound; defaults to MAXBOUND
    integer(I4B) :: nleading = 0 !< count of leading (pre-keystring) columns
    character(len=LINELENGTH), dimension(:), allocatable :: params !< in-scope param tagnames
    type(ModflowInputType) :: mf6_input !< description of input
    ! keystring-specific (only set for KEYSTRING loads)
    logical(LGP) :: keystring_by_feature = .false. !< .true. for feature-addressed (IFNO/NUMBER/BNDNO) KEYSTRING loadtype
    logical(LGP) :: keystring_by_node = .false. !< .true. for CELLID-addressed (e.g. TVK/TVS) KEYSTRING loadtype
    logical(LGP) :: has_setting_dispatch = .false. !< .true. when keystring_by_feature .or. keystring_by_node
    type(InputParamDefinitionType), pointer :: setting_idt => null() !< internal idt for SETTING column
    integer(I4B) :: nkeystring_items = 0 !< count of expanded keystring items (0 for non-keystring)
    type(KeystringItemType), allocatable :: keystring_items(:)
    character(len=LENVARNAME) :: feature_id_varname = '' !< leading id column mf6varname (feature-addressed)
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
    procedure, private :: keystring_item_names
    procedure, private :: resolve_dimensions
    procedure, private :: scale_keystring_maxbound
    procedure, private :: allocate_param
    procedure :: check_developmode
    procedure, private :: in_scope
    procedure :: subindex_dependency
    procedure, private :: option_check
    procedure :: build_keystring_items
    procedure, private :: resolve_nfeatures
    procedure :: resolve_item_nfeatures
    procedure, private :: shape_param_is_array
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
    ! build the per-item descriptor table (keystring only)
    if (this%loadtype == KEYSTRING) call this%build_keystring_items()
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
    use DefinitionSelectModule, only: idt_default, idt_parse_rectype, &
                                      get_aggregate_definition_type
    class(LoadContextType) :: this
    type(InputParamDefinitionType), pointer :: idt, aidt
    character(len=LINELENGTH), dimension(:), allocatable :: cols
    integer(I4B) :: n, nparam

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

    this%is_advanced = is_advanced(this%mf6_input)

    ! classify by the KEYSTRING aggregate's own leading column
    if (this%loadtype == KEYSTRING) then
      aidt => &
        get_aggregate_definition_type(this%mf6_input%aggregate_dfns, &
                                      this%mf6_input%component_type, &
                                      this%mf6_input%subcomponent_type, &
                                      this%blockname)
      call idt_parse_rectype(aidt, cols, nparam)
      if (is_feature_tag(cols(1))) then
        this%keystring_by_feature = .true.
      else if (cols(1) == 'CELLID') then
        this%keystring_by_node = .true.
      end if
      if (allocated(cols)) deallocate (cols)
      this%has_setting_dispatch = &
        this%keystring_by_feature .or. this%keystring_by_node
      if (this%has_setting_dispatch) then
        this%setting_idt => &
          idt_default(this%mf6_input%component_type, &
                      this%mf6_input%subcomponent_type, &
                      this%blockname, 'SETTING', 'SETTING', 'STRING')
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
  !! KEYSTRING items, so every feature can use every setting in one
  !! period.
  !<
  subroutine scale_keystring_maxbound(this)
    class(LoadContextType) :: this
    integer(I4B) :: nkeystring_items

    nkeystring_items = this%nkeystring_items

    if (nkeystring_items > 0) then
      if (this%maxbound == 0) then
        if (.not. this%keystring_by_feature) then
          this%maxbound = this%nodes * nkeystring_items
        end if
        ! else: feature-addressed, genuinely zero stays 0
      else
        this%maxbound = this%maxbound * nkeystring_items
      end if
    end if
  end subroutine scale_keystring_maxbound

  !> @brief allocate arrays
  !!
  !! Call after input parameters are allocated: after load_params() with
  !! create for array-based loaders, or after all mem_create_vector()
  !! calls for list-based load.
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
  !! Required/structural params are handled by the caller; generic
  !! conditions are checked first, then package-specific ones by
  !! subcomponent_type.
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

    ! advanced packages: every optional leaf param is in scope
    if (this%is_advanced) then
      in_scope = .true.
      return
    end if

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
    case default
      ! unrecognized subcomponent with an optional param not handled
      ! above is a development error -- abort so the gap is visible
      errmsg = 'LoadContext in_scope needs new case for: '// &
               trim(this%mf6_input%subcomponent_type)//'/'//trim(tagname)
      call store_error(errmsg, .true.)
    end select
  end function in_scope

  !> @brief Hardcoded per-package dependency for a subindex param
  !! with no SHAPE of its own -- same category of special-casing as
  !! in_scope, for a dimension name and subindex field instead.
  !<
  function subindex_dependency(this, tagname, dimname, subindex_tagname) &
    result(found)
    class(LoadContextType) :: this
    character(len=*), intent(in) :: tagname
    character(len=LENVARNAME), intent(out) :: dimname
    character(len=LENVARNAME), intent(out) :: subindex_tagname
    logical(LGP) :: found

    found = .false.
    dimname = ''
    subindex_tagname = ''
    select case (this%mf6_input%subcomponent_type)
    case ('SFR')
      if (tagname == 'DIVFLOW') then
        dimname = 'NDV'
        subindex_tagname = 'IDV'
        found = .true.
      end if
    end select
  end function subindex_dependency

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

  !> @brief Build the per-item descriptor table, the single source of
  !! truth for how each keystring item is allocated and applied.
  !<
  subroutine build_keystring_items(this)
    use DefinitionSelectModule, only: get_param_definition_type
    class(LoadContextType) :: this
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: icol, k, padj, nfeatures, nkeystring_items
    character(len=LENVARNAME) :: dimname, subindex_tagname
    logical(LGP) :: found, node
    character(len=LINELENGTH), allocatable :: item_names(:)
    integer(I4B), allocatable :: head_nbody(:)
    logical(LGP), allocatable :: item_is_body(:)
    character(len=LENVARNAME) :: cur_head !< most recent head's SETTING keyword

    if (.not. this%has_setting_dispatch) return

    ! derive per-item metadata (single computation; not stored as ctx state)
    call this%keystring_item_names(item_names, head_nbody, item_is_body, &
                                   nkeystring_items)
    if (nkeystring_items < 1) return

    ! generic leading-id address varname (feature-addressed packages,
    ! advanced or not, e.g. SPC): the leading column's mf6varname.
    if (this%keystring_by_feature) then
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       this%blockname, this%params(1), '')
      this%feature_id_varname = trim(idt%mf6varname)
    end if

    padj = 1 ! SETTING column present whenever has_setting_dispatch
    node = this%keystring_by_node

    allocate (this%keystring_items(nkeystring_items))
    nfeatures = this%resolve_nfeatures()
    cur_head = ''

    do icol = this%nleading + 1, size(this%params)
      k = icol - this%nleading
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       this%blockname, this%params(icol), '')
      this%keystring_items(k)%idt => idt
      this%keystring_items(k)%sa_icol = icol + padj
      this%keystring_items(k)%is_body = item_is_body(k)
      this%keystring_items(k)%head_nbody = head_nbody(k)

      ! track the current record head's SETTING keyword; a body dispatches
      ! on its owning head's keyword (the SETTING token on the input row),
      ! not on its own mf6varname.
      if (.not. this%keystring_items(k)%is_body) then
        if (this%keystring_items(k)%head_nbody > 0) then
          cur_head = trim(idt%mf6varname) ! opening a RECORD
        else
          cur_head = '' ! standalone item; not a record head
        end if
      else
        this%keystring_items(k)%head_setting_varname = trim(cur_head)
      end if

      ! managed only for DOUBLE + timeseries standalone settings, which
      ! require a permanent array for cross-period timeseries persistence.
      this%keystring_items(k)%idm_managed = &
        (idt%datatype == 'DOUBLE' .and. idt%timeseries .and. &
         .not. this%keystring_items(k)%is_body)

      if (this%keystring_items(k)%idm_managed) then
        if (node) then
          this%keystring_items(k)%addr_mode = ADDR_NODE
          this%keystring_items(k)%init_value = DNODATA
          if (associated(this%nodes)) &
            this%keystring_items(k)%nfeatures = this%nodes
        else
          this%keystring_items(k)%addr_mode = ADDR_FEATURE
          this%keystring_items(k)%init_value = DZERO
          this%keystring_items(k)%nfeatures = &
            this%resolve_item_nfeatures(idt%tagname, idt%shape, nfeatures)
        end if
      end if

      ! subindex (SFR DIVFLOW): classify here; df-time fields
      ! (index_icol/head_icol/offsets/nfeatures) are completed in df().
      if (this%is_advanced .and. this%keystring_items(k)%is_body) then
        found = this%subindex_dependency(idt%tagname, dimname, &
                                         subindex_tagname)
        if (found) then
          this%keystring_items(k)%addr_mode = ADDR_SUBINDEX
          this%keystring_items(k)%idm_managed = .true.
        end if
      end if
    end do
  end subroutine build_keystring_items

  !> @brief Resolve the permanent array's feature count. The explicit
  !! DIMENSIONS dimension (named_bound) is authoritative when present;
  !! PACKAGEDATA row count is the fallback. Both present and disagreeing
  !! is an error.
  !<
  function resolve_nfeatures(this) result(nfeatures)
    use MemoryManagerModule, only: get_isize
    use MemoryManagerExtModule, only: mem_set_value
    use SimModule, only: store_error
    class(LoadContextType) :: this
    integer(I4B) :: nfeatures
    integer(I4B) :: isize, nkeystring_items, nrow, ndim
    integer(I4B), pointer :: dimval
    logical(LGP) :: have_dim, have_nrow
    character(len=LINELENGTH) :: errmsg

    ! read the RAW dimension scalar from the memory manager, not
    ! this%maxbound (a scaled copy used only for read-array sizing)
    ndim = 0
    have_dim = .false.
    allocate (dimval)
    dimval = 0
    call mem_set_value(dimval, this%named_bound, this%mf6_input%mempath, &
                       have_dim, release=.false.)
    if (have_dim) ndim = dimval
    deallocate (dimval)

    call get_isize('PACKAGEDATA_IFNO', this%mf6_input%mempath, isize)
    have_nrow = (isize > 0)
    nrow = isize

    if (have_dim .and. have_nrow .and. ndim /= nrow) then
      write (errmsg, '(a,1x,a,1x,i0,1x,a,1x,i0,a)') &
        'IDM dimension mismatch:', trim(this%named_bound), ndim, &
        'does not match the PACKAGEDATA row count', nrow, '.'
      call store_error(errmsg)
    end if

    if (have_dim) then
      nfeatures = ndim
      return
    else if (have_nrow) then
      nfeatures = nrow
      return
    end if

    ! defensive: only if neither an explicit dimension nor PACKAGEDATA
    ! is available (not expected for current packages)
    nfeatures = 0
    nkeystring_items = this%nkeystring_items
    if (nkeystring_items > 0 .and. associated(this%maxbound)) then
      if (this%maxbound > 0) nfeatures = this%maxbound / nkeystring_items
    end if
  end function resolve_nfeatures

  !> @brief Per-item feature count from a named dimension, falling back to
  !! default_nfeatures if unset. Populated-then-released is an error.
  !<
  function resolve_item_nfeatures(this, item_tagname, dimname, &
                                  default_nfeatures) result(nfeatures)
    use MemoryManagerModule, only: get_isize, mem_setptr
    use SimModule, only: store_error
    class(LoadContextType) :: this
    character(len=*), intent(in) :: item_tagname
    character(len=*), intent(in) :: dimname
    integer(I4B), intent(in) :: default_nfeatures
    integer(I4B) :: nfeatures
    integer(I4B), pointer :: shape_val => null()
    integer(I4B), dimension(:), pointer, contiguous :: shape_arr => null()
    integer(I4B) :: isize
    character(len=LINELENGTH) :: errmsg

    nfeatures = default_nfeatures
    if (dimname == '') return
    call get_isize(trim(dimname), this%mf6_input%mempath, isize)
    if (isize < 0) then
      nfeatures = 0
      return
    else if (isize == 0) then
      write (errmsg, '(a,1x,a,1x,a)') &
        'item', trim(item_tagname)//': DIMENSION', &
        trim(dimname)//' is not defined.'
      call store_error(errmsg)
      nfeatures = 0
      return
    else if (.not. this%shape_param_is_array(trim(dimname))) then
      call mem_setptr(shape_val, trim(dimname), this%mf6_input%mempath)
      nfeatures = shape_val
    else
      call mem_setptr(shape_arr, trim(dimname), this%mf6_input%mempath)
      nfeatures = sum(shape_arr)
    end if
  end function resolve_item_nfeatures

  !> @brief Is shape_varname a per-feature array (PACKAGEDATA) rather than
  !! a package-wide scalar (DIMENSIONS)?
  !<
  function shape_param_is_array(this, shape_varname) result(is_array)
    class(LoadContextType) :: this
    character(len=*), intent(in) :: shape_varname
    logical(LGP) :: is_array
    integer(I4B) :: i

    is_array = .false.
    do i = 1, size(this%mf6_input%param_dfns)
      if (this%mf6_input%param_dfns(i)%component_type == &
          this%mf6_input%component_type .and. &
          this%mf6_input%param_dfns(i)%subcomponent_type == &
          this%mf6_input%subcomponent_type .and. &
          trim(this%mf6_input%param_dfns(i)%mf6varname) == &
          trim(shape_varname)) then
        is_array = &
          (trim(this%mf6_input%param_dfns(i)%blockname) == 'PACKAGEDATA')
        exit
      end if
    end do
  end function shape_param_is_array

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
    character(len=LINELENGTH), allocatable :: item_names(:)
    integer(I4B), allocatable :: head_nbody(:)
    logical(LGP), allocatable :: item_is_body(:)
    integer(I4B) :: keepcnt, iparam, nparam, nkeystring_items, n
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

    ! record leading-column count before item expansion
    if (this%loadtype == LIST .or. &
        this%loadtype == KEYSTRING) this%nleading = keepcnt

    ! for keystring packages: append item names (metadata head_nbody/
    ! item_is_body is derived into the descriptor by build_keystring_items)
    if (this%loadtype == KEYSTRING) then
      call this%keystring_item_names(item_names, head_nbody, &
                                     item_is_body, nkeystring_items)
      this%nkeystring_items = nkeystring_items
      do n = 1, nkeystring_items
        keepcnt = keepcnt + 1
        call expandarray(param_buf)
        param_buf(keepcnt) = trim(item_names(n))
      end do
    end if

    ! update nparam to total (leading + items)
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

  !> @brief Append body column names from a RECORD compound entry to item_names.
  !<
  subroutine expand_record_body(mf6_input, rec_idt, item_names, nkeystring_items)
    use InputOutputModule, only: upcase
    use ArrayHandlersModule, only: expandarray
    use DefinitionSelectModule, only: idt_parse_rectype, idt_datatype
    type(ModflowInputType), intent(in) :: mf6_input
    type(InputParamDefinitionType), pointer, intent(in) :: rec_idt
    character(len=LINELENGTH), allocatable, intent(inout) :: item_names(:)
    integer(I4B), intent(inout) :: nkeystring_items
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
        nkeystring_items = nkeystring_items + 1
        call expandarray(item_names)
        item_names(nkeystring_items) = trim(sub_idt%tagname)
        exit
      end do
    end do
    if (allocated(sub_cols)) deallocate (sub_cols)
  end subroutine expand_record_body

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

  !> @brief .true. if mf6_input's PERIOD block is feature-addressed
  !! (see is_feature_tag). Block-independent, unlike LoadContextType's own
  !! keystring_by_feature, which is only valid for a PERIOD-scoped instance.
  !<
  function is_feature_keystring(mf6_input) result(res)
    use DefinitionSelectModule, only: get_aggregate_definition_type, &
                                      idt_parse_rectype
    type(ModflowInputType), intent(in) :: mf6_input
    logical(LGP) :: res
    type(InputParamDefinitionType), pointer :: aidt
    character(len=LINELENGTH), allocatable :: cols(:)
    integer(I4B) :: ncol
    res = .false.
    if (.not. is_keystring_period(mf6_input)) return
    aidt => get_aggregate_definition_type(mf6_input%aggregate_dfns, &
                                          mf6_input%component_type, &
                                          mf6_input%subcomponent_type, &
                                          'PERIOD')
    call idt_parse_rectype(aidt, cols, ncol)
    res = is_feature_tag(cols(1))
    if (allocated(cols)) deallocate (cols)
  end function is_feature_keystring

  !> @brief .true. if tagname is a record's leading id column (IFNO and
  !! its legacy package-specific aliases). Single source for both
  !! resolve_loadtype and is_feature_keystring.
  !<
  function is_feature_tag(tagname) result(res)
    character(len=*), intent(in) :: tagname
    logical(LGP) :: res
    select case (tagname)
    case ('IFNO', 'NUMBER', 'BNDNO', 'RNO', 'LAKENO', 'MAWNO', 'UZFNO')
      res = .true.
    case default
      res = .false.
    end select
  end function is_feature_tag

  !> @brief Return .true. if mf6_input is an advanced package
  !<
  function is_advanced(mf6_input) result(res)
    type(ModflowInputType), intent(in) :: mf6_input
    logical(LGP) :: res
    res = idm_is_advanced(mf6_input%component_type, mf6_input%subcomponent_type)
  end function is_advanced

  !> @brief Return keystring item column names, per-head body counts, and
  !! body flags. A RECORD group's trailing columns are body members; its head
  !! and any direct-dispatch param are not.
  !<
  subroutine keystring_item_names(this, item_names, head_nbody, &
                                  item_is_body, nkeystring_items)
    use InputOutputModule, only: upcase
    use ArrayHandlersModule, only: expandarray
    use DefinitionSelectModule, only: idt_parse_rectype, idt_datatype, &
                                      get_aggregate_definition_type
    class(LoadContextType) :: this
    character(len=LINELENGTH), allocatable, intent(out) :: item_names(:)
    integer(I4B), allocatable, intent(out) :: head_nbody(:)
    logical(LGP), allocatable, intent(out) :: item_is_body(:)
    integer(I4B), intent(out) :: nkeystring_items
    type(InputParamDefinitionType), pointer :: aidt, ks_aidt, idt
    character(len=LINELENGTH), allocatable :: rec_cols(:), ks_cols(:)
    character(len=LINELENGTH) :: rec_token, tagname
    integer(I4B) :: m, n, nrec_col, nks_col, nitems0, k

    nkeystring_items = 0

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

    ! parse the KEYSTRING aggregate to get item token list — canonical order
    call idt_parse_rectype(ks_aidt, ks_cols, nks_col)

    ! walk the keystring token list in aggregate order
    do m = 1, nks_col
      rec_token = trim(ks_cols(m))
      call upcase(rec_token)

      ! locate matching param_dfns entry for this token
      do n = 1, size(this%mf6_input%param_dfns)
        if (this%mf6_input%param_dfns(n)%blockname /= this%blockname) cycle
        tagname = this%mf6_input%param_dfns(n)%tagname
        call upcase(tagname)
        if (trim(tagname) /= trim(rec_token)) cycle

        idt => this%mf6_input%param_dfns(n)
        if (idt_datatype(idt) == 'RECORD') then
          ! compound group: expand body members in RECORD type order
          nitems0 = nkeystring_items
          call expand_record_body(this%mf6_input, idt, item_names, &
                                  nkeystring_items)
          ! first added entry is the KEYWORD head; remaining are body members
          do k = nitems0 + 1, nkeystring_items
            call expandarray(head_nbody)
            call expandarray(item_is_body)
            if (k == nitems0 + 1) then
              head_nbody(k) = nkeystring_items - nitems0 - 1
              item_is_body(k) = .false.
            else
              head_nbody(k) = 0
              item_is_body(k) = .true.
            end if
          end do
        else
          ! direct-dispatch param
          nkeystring_items = nkeystring_items + 1
          call expandarray(item_names)
          call expandarray(head_nbody)
          call expandarray(item_is_body)
          item_names(nkeystring_items) = &
            trim(this%mf6_input%param_dfns(n)%tagname)
          head_nbody(nkeystring_items) = 0
          item_is_body(nkeystring_items) = .false.
        end if
        exit
      end do
    end do

    if (allocated(ks_cols)) deallocate (ks_cols)
  end subroutine keystring_item_names

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
