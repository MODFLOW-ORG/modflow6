!> @brief This module contains the LoadContextModule
!!
!! This module creates a load context for IDM generic
!! loaders (ListLoadType, LayerArrayLoadType, GridArrayLoadType)
!! that supports consistent package side access. It also
!! determines in scope parameters for the generic dynamic
!! loaders and all structarray based static loads.
!!
!<
module LoadContextModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, IZERO, LINELENGTH, LENAUXNAME, &
                             LENVARNAME, LENBOUNDNAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use ModflowInputModule, only: ModflowInputType
  use InputDefinitionModule, only: InputParamDefinitionType
  use CharacterStringModule, only: CharacterStringType

  implicit none
  private
  public :: LoadContextType
  public :: ReadStateVarType
  public :: rsv_name

  enum, bind(C)
    enumerator :: LOAD_UNDEF = 0 !< undefined load type
    enumerator :: LIST = 1 !< list (structarray) based load
    enumerator :: LAYERARRAY = 2 !< readasarrays load
    enumerator :: GRIDARRAY = 3 !< readarraygrid load
  end enum

  enum, bind(C)
    enumerator :: CONTEXT_UNDEF = 0 !< undefined context type
    enumerator :: ROOT = 1 !< root context type
    enumerator :: SIM = 2 !< sim context type
    enumerator :: MODEL = 3 !< model context type
    enumerator :: MODELPKG = 4 !< model package context type
    enumerator :: EXCHANGE = 5 !< exchange context type
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

  !> @brief derived type for boundary package input context
  !!
  !! Input Load Context for generic dynamic loaders and
  !! StructArray based static loads
  !!
  !<
  type :: LoadContextType
    character(len=LENVARNAME) :: blockname !< load block name
    character(len=LENVARNAME) :: named_bound !< name of dimensions relevant to load
    integer(I4B), pointer :: naux => null() !< number of auxiliary variables
    integer(I4B), pointer :: maxbound => null() !< value associated with named_bound
    integer(I4B), pointer :: boundnames => null() !< are bound names optioned
    integer(I4B), pointer :: iprpak => null() ! print input option
    integer(I4B), pointer :: nbound => null() !< number of bounds in period
    integer(I4B), pointer :: ncpl => null() !< ncpl associated with model shape
    integer(I4B), pointer :: nodes => null() !< nodes associated with model shape
    integer(I4B) :: loadtype !< enum load type
    integer(I4B) :: ctxtype !< enum context type
    logical(LGP) :: readarray !< is this an array based load
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: auxname_cst => null() !< array of auxiliary names
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: boundname_cst => null() !< array of bound names
    real(DP), dimension(:, :), pointer, &
      contiguous :: auxvar => null() !< auxiliary variable array
    integer(I4B), dimension(:), pointer, contiguous :: mshape => null() !< model shape
    character(len=LINELENGTH), dimension(:), allocatable :: params !< in scope param tags
    type(ModflowInputType) :: mf6_input !< description of input
  contains
    procedure :: init
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: allocate_param
    procedure :: tags
    procedure :: in_scope
    procedure :: set_params
    procedure :: rsv_alloc
    procedure :: destroy
  end type LoadContextType

contains

  !> @brief init loader context object
  !<
  subroutine init(this, mf6_input, blockname, named_bound)
    use InputOutputModule, only: upcase
    use ModelPackageInputsModule, only: supported_model
    class(LoadContextType) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), optional, intent(in) :: blockname
    character(len=*), optional, intent(in) :: named_bound
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: n

    this%mf6_input = mf6_input
    this%readarray = .false.
    this%loadtype = LOAD_UNDEF
    this%ctxtype = CONTEXT_UNDEF

    select case (mf6_input%load_scope)
    case ('ROOT')
      this%ctxtype = ROOT
    case ('SIM')
      if (mf6_input%subcomponent_type == 'NAM') then
        this%ctxtype = MODEL
      else if (mf6_input%subcomponent_type == 'TDIS' .or. &
               mf6_input%subcomponent_type == 'HPC') then
        this%ctxtype = SIM
      else if (mf6_input%component_type == 'EXG') then
        this%ctxtype = EXCHANGE
      end if
    case ('MODEL')
      this%ctxtype = MODELPKG
    case default
    end select

    if (this%ctxtype == CONTEXT_UNDEF) then
      errmsg = 'LoadContext unidentified context for mempath: '// &
               trim(mf6_input%mempath)
      call store_error(errmsg, .true.)
    end if

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

    ! determine if list based load
    do n = 1, size(mf6_input%block_dfns)
      if (mf6_input%block_dfns(n)%blockname == this%blockname) then
        if (mf6_input%block_dfns(n)%aggregate) then
          this%loadtype = LIST
        end if
      end if
    end do

    ! determine if array based load
    if (this%loadtype == LOAD_UNDEF) then
      do n = 1, size(mf6_input%param_dfns)
        idt => mf6_input%param_dfns(n)
        if (idt%blockname == 'OPTIONS') then
          select case (idt%tagname)
          case ('READASARRAYS')
            this%loadtype = LAYERARRAY
            this%readarray = .true.
          case ('READARRAYGRID')
            this%loadtype = GRIDARRAY
            this%readarray = .true.
          case default
            ! no-op
          end select
        end if
      end do
    end if

    ! set in scope params for load
    call this%set_params()

    ! allocate load context scalars
    call this%allocate_scalars()
  end subroutine init

  !> @brief allocate scalars
  !<
  subroutine allocate_scalars(this)
    use MemoryManagerModule, only: mem_setptr
    class(LoadContextType) :: this

    if (this%ctxtype == EXCHANGE .or. &
        this%ctxtype == MODELPKG) then

      call setptr(this%nbound, 'NBOUND', this%mf6_input%mempath)
      call setval(this%naux, 'NAUX', this%mf6_input%mempath)
      call setval(this%ncpl, 'NCPL', this%mf6_input%mempath)
      call setval(this%nodes, 'NODES', this%mf6_input%mempath)
      call setval(this%maxbound, this%named_bound, this%mf6_input%mempath)
      call setval(this%boundnames, 'BOUNDNAMES', this%mf6_input%mempath)
      call setval(this%iprpak, 'IPRPAK', this%mf6_input%mempath)

      ! reset nbound
      this%nbound = 0
    end if

    if (this%ctxtype == MODELPKG .and. &
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
    end if
  end subroutine allocate_scalars

  !> @brief allocate arrays
  !!
  !! call this routine after input parameters have been allocated,
  !!   e.g. after load_params() with create has been called for array
  !!   based loaders or after all mem_create_vector() calls have
  !!   been made for list based load.
  !!
  !<
  subroutine allocate_arrays(this)
    use MemoryManagerModule, only: mem_allocate, mem_setptr, get_isize
    class(LoadContextType) :: this
    integer(I4B), dimension(:, :), pointer, contiguous :: cellid
    integer(I4B), dimension(:), pointer, contiguous :: nodeulist

    if (this%ctxtype == MODELPKG .and. &
        this%blockname == 'PERIOD') then
      ! allocate cellid if this is not list input
      if (this%readarray) then
        call mem_allocate(cellid, 0, 0, 'CELLID', this%mf6_input%mempath)
      end if

      ! allocate nodeulist
      if (this%loadtype /= GRIDARRAY) then
        call mem_allocate(nodeulist, 0, 'NODEULIST', this%mf6_input%mempath)
      end if

      ! set pointers to arrays
      call setptr(this%auxname_cst, 'AUXILIARY', &
                  this%mf6_input%mempath, LENAUXNAME)
      call setptr(this%boundname_cst, 'BOUNDNAME', &
                  this%mf6_input%mempath, LENBOUNDNAME)
      call setptr(this%auxvar, this%mf6_input%mempath)

    else if (this%ctxtype == EXCHANGE) then
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

    if (this%readarray) then
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
      else if (this%readarray) then
        call allocate_int1d(dimsize, idt%mf6varname, &
                            this%mf6_input%mempath)
      end if
    case ('DOUBLE1D')
      if (idt%shape == 'NAUX') then
        call allocate_dbl2d(this%naux, this%maxbound, &
                            idt%mf6varname, this%mf6_input%mempath)
      else if (this%readarray) then
        call allocate_dbl1d(dimsize, idt%mf6varname, &
                            this%mf6_input%mempath)
      end if
    case ('DOUBLE2D')
      if (this%readarray) then
        call allocate_dbl2d(this%naux, dimsize, idt%mf6varname, &
                            this%mf6_input%mempath)
      end if
    case default
    end select
  end subroutine allocate_param

  !> @brief get in scope package params
  !!
  !! set input array to tagnames of in scope params, optionally
  !! allocate the parameters based on datatype.
  !!
  !<
  subroutine tags(this, params, nparam, input_name, create)
    use DefinitionSelectModule, only: get_param_definition_type
    class(LoadContextType) :: this
    character(len=LINELENGTH), dimension(:), allocatable, &
      intent(inout) :: params
    integer(I4B), intent(inout) :: nparam
    character(len=*), intent(in) :: input_name
    logical(LGP), optional, intent(in) :: create
    type(InputParamDefinitionType), pointer :: idt
    logical(LGP) :: allocate_params
    integer(I4B) :: n

    ! initialize allocate_params
    allocate_params = .false.

    ! override default if provided
    if (present(create)) then
      allocate_params = create
    end if

    if (allocated(params)) deallocate (params)
    nparam = size(this%params)
    allocate (params(nparam))
    do n = 1, nparam
      params(n) = this%params(n)
    end do

    if (allocate_params) then
      ! allocate dfn input params
      do n = 1, nparam
        idt => &
          get_param_definition_type(this%mf6_input%param_dfns, &
                                    this%mf6_input%component_type, &
                                    this%mf6_input%subcomponent_type, &
                                    this%blockname, params(n), '')
        call this%allocate_param(idt)
      end do
    end if
  end subroutine tags

  !> @brief establish if input parameter is in scope for package load
  !<
  function in_scope(this, mf6_input, blockname, tagname)
    use MemoryManagerModule, only: get_isize, mem_setptr
    use DefinitionSelectModule, only: get_param_definition_type
    class(LoadContextType) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), intent(in) :: blockname
    character(len=*), intent(in) :: tagname
    logical(LGP) :: in_scope
    type(InputParamDefinitionType), pointer :: idt
    character(len=LENVARNAME) :: checkname
    integer(I4B) :: isize, checksize
    integer(I4B), pointer :: intptr

    idt => &
      get_param_definition_type(mf6_input%param_dfns, &
                                mf6_input%component_type, &
                                mf6_input%subcomponent_type, &
                                blockname, tagname, '')
    if (idt%required) then
      in_scope = .true.
      return
    else
      in_scope = .false.
    end if

    ! initialize
    checkname = ''
    checksize = 0

    if (tagname == 'AUXVAR' .or. &
        tagname == 'AUX') then
      checkname = 'NAUX'
    else if (tagname == 'BOUNDNAME') then
      checkname = 'BOUNDNAMES'
    else if (tagname == 'I'//trim(mf6_input%subcomponent_type(1:3))) then
      if (this%loadtype == LAYERARRAY) in_scope = .true.
    else
      select case (mf6_input%subcomponent_type)
      case ('EVT')
        if (tagname == 'PXDP' .or. tagname == 'PETM') then
          checkname = 'NSEG'
          checksize = 1
        else if (tagname == 'PETM0') then
          checkname = 'SURFRATESPEC'
        end if
      case ('MVR', 'MVT', 'MVE')
        if (tagname == 'MNAME' .or. &
            tagname == 'MNAME1' .or. &
            tagname == 'MNAME2') then
          checkname = 'MODELNAMES'
        end if
      case ('NAM')
        in_scope = .true.
      case ('SSM')
        if (tagname == 'MIXED') in_scope = .true.
      case default
        errmsg = 'LoadContext in_scope needs new check for: '// &
                 trim(idt%tagname)
        call store_error(errmsg, .true.)
      end select
    end if

    ! apply checks
    if (.not. in_scope) then
      call get_isize(checkname, mf6_input%mempath, isize)
      if (isize > 0) then
        call mem_setptr(intptr, checkname, mf6_input%mempath)
        if (intptr > checksize) in_scope = .true.
      end if
    end if
  end function in_scope

  !> @brief set set of in scope parameters for package
  !<
  subroutine set_params(this)
    use ArrayHandlersModule, only: expandarray
    use DefinitionSelectModule, only: get_param_definition_type, &
                                      get_aggregate_definition_type, &
                                      idt_parse_rectype
    class(LoadContextType) :: this
    type(InputParamDefinitionType), pointer :: idt, aidt
    character(len=LINELENGTH), dimension(:), allocatable :: tags
    character(len=LINELENGTH), dimension(:), allocatable :: cols
    integer(I4B) :: keepcnt, iparam, nparam
    logical(LGP) :: keep

    ! initialize
    keepcnt = 0

    if (this%loadtype == LIST) then
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
      if (this%loadtype == LIST) then
        idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                         this%mf6_input%component_type, &
                                         this%mf6_input%subcomponent_type, &
                                         this%blockname, cols(iparam), '')
      else
        idt => this%mf6_input%param_dfns(iparam)
      end if

      if (idt%blockname /= this%blockname) then
        keep = .false.
      else
        keep = this%in_scope(this%mf6_input, this%blockname, idt%tagname)
      end if

      if (keep) then
        keepcnt = keepcnt + 1
        call expandarray(tags)
        tags(keepcnt) = trim(idt%tagname)
      end if
    end do

    ! update nparam
    nparam = keepcnt

    ! allocate filtcols
    allocate (this%params(nparam))

    ! set filtcols
    do iparam = 1, nparam
      this%params(iparam) = trim(tags(iparam))
    end do

    ! cleanup
    if (allocated(tags)) deallocate (tags)
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

    if (this%ctxtype == EXCHANGE .or. &
        this%ctxtype == MODELPKG) then
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

  !> @brief allocate intptr and update from input contextset intptr to varname
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
    call mem_set_value(intptr, varname, mempath, found)
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
