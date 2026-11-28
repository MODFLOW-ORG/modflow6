!> @brief Grid array period block loader
!!
!! Loads PERIOD block grid array input for packages with the
!! READARRAYGRID option. The base type reads the block arrays; the
!! storage mode, fixed by the package dfn, stores them and keeps them
!! between periods:
!!  - Compressed mode (GridArrayCompressedLoadModule): the dfn defines
!!    MAXBOUND, e.g. CHDG, WELG, DRNG, RIVG, GHBG. The PERIOD block
!!    defines boundary cells; the package reads them as a list (NODEULIST,
!!    NBOUND) with arrays indexed by boundary.
!!  - Full-grid mode (GridArrayFullLoadModule): the dfn does not define
!!    MAXBOUND. The PERIOD block changes values on a fixed grid; the
!!    package reads arrays indexed by user node.
!!
!! In both modes values are taken only at active cells (IDOMAIN > 0), and
!! an array given as a time-array series is linked to the series and
!! updated each time step until the array is given again or, in
!! compressed mode, an empty block is read.
!<
module GridArrayLoadModule

  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LINELENGTH, LENVARNAME, LENTIMESERIESNAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use InputDefinitionModule, only: InputParamDefinitionType
  use MemoryManagerModule, only: mem_allocate, mem_setptr
  use BlockParserModule, only: BlockParserType
  use ModflowInputModule, only: ModflowInputType
  use CharacterStringModule, only: CharacterStringType
  use LoadContextModule, only: LoadContextType
  use IdmParamMetaModule, only: ReadStateVarType, RS_NONE, RS_VALUES, &
                                RS_SERIES, RS_NETCDF
  use AsciiInputLoadTypeModule, only: AsciiDynamicPkgLoadBaseType
  use TimeArraySeriesManagerModule, only: TimeArraySeriesManagerType, &
                                          tasmanager_cr
  use TimeArraySeriesModule, only: TimeArraySeriesType

  implicit none
  private
  public :: GridArrayLoadType
  public :: GridSlotType
  public :: defines_maxbound

  !> @brief A period block array: a main array or an auxiliary variable
  !<
  type :: GridSlotType
    character(len=LINELENGTH) :: tagname = '' !< array name for messages
    character(len=LENVARNAME) :: mf6varname = '' !< input context name
    integer(I4B) :: iparam = 0 !< dynamic param index (AUX for aux arrays)
    integer(I4B) :: iaux = 0 !< aux index, 0 for a main array
    logical(LGP) :: given = .false. !< given in a period block
    logical(LGP) :: inblock = .false. !< given in the current block
    character(len=LENTIMESERIESNAME) :: tasname = '' !< linked series name, blank if not linked
  end type GridSlotType

  !> @brief Grid array dynamic loader base type
  !!
  !! Extended by each storage mode, which provides mode_init, array_store
  !! and tas_links_update, and block_end and tas_start_array if it needs
  !! them.
  !<
  type, extends(AsciiDynamicPkgLoadBaseType) :: GridArrayLoadType
    type(ReadStateVarType), dimension(:), allocatable :: param_reads !< read states for current load
    type(LoadContextType) :: ctx
    logical(LGP), dimension(:), allocatable :: active !< user node IDOMAIN > 0
    type(GridSlotType), dimension(:), allocatable :: slots !< main arrays, then aux variables
    integer(I4B) :: nmain = 0 !< number of main arrays
    integer(I4B), dimension(:), allocatable :: blk !< block array slots, in input order
    integer(I4B) :: nreads = 0 !< number of block arrays
    type(TimeArraySeriesManagerType), pointer :: tasmanager => null() !< TAS manager
    logical(LGP) :: tas_active = .false. !< TAS6 files provided
    integer(I4B), dimension(:), pointer, &
      contiguous :: aux_reads => null() !< read state of each auxiliary variable
  contains
    procedure :: ainit
    procedure :: df
    procedure :: ts_advance
    procedure :: rp
    procedure :: line_load
    procedure :: destroy
    procedure :: reset
    procedure :: params_alloc
    procedure :: param_dfn
    procedure :: param_load
    procedure :: set_active
    procedure :: slot_index
    procedure :: slot_clear
    procedure :: array_add
    procedure :: tas_load
    procedure :: make_links
    procedure :: mode_init
    procedure :: array_store
    procedure :: block_end
    procedure :: tas_links_update
    procedure :: tas_start_array
  end type GridArrayLoadType

contains

  subroutine ainit(this, mf6_input, component_name, &
                   component_input_name, input_name, &
                   iperblock, parser, iout)
    use MemoryManagerModule, only: get_isize
    use LoadMf6FileModule, only: LoadMf6FileType
    class(GridArrayLoadType), intent(inout) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), intent(in) :: component_name
    character(len=*), intent(in) :: component_input_name
    character(len=*), intent(in) :: input_name
    integer(I4B), intent(in) :: iperblock
    type(BlockParserType), pointer, intent(inout) :: parser
    integer(I4B), intent(in) :: iout
    type(LoadMf6FileType) :: loader
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: tas_fnames
    character(len=LINELENGTH) :: fname
    integer(I4B) :: isize, n
    integer(I4B), pointer :: maxbound

    ! initialize base type
    call this%DynamicPkgLoadType%init(mf6_input, component_name, &
                                      component_input_name, &
                                      input_name, iperblock, iout)

    ! load static input
    call loader%load(parser, mf6_input, this%nc_vars, this%input_name, iout)

    ! maxbound is optional
    call get_isize('MAXBOUND', mf6_input%mempath, isize)
    if (isize < 0) then
      ! set maxbound to grid nodes
      call mem_allocate(maxbound, 'MAXBOUND', mf6_input%mempath)
      maxbound = product(loader%mshape)
    end if

    ! create full grid time-array series manager
    allocate (this%tasmanager)
    call tasmanager_cr(this%tasmanager, &
                       modelname=this%mf6_input%component_name, &
                       iout=this%iout, fullgrid=.true.)

    ! add TAS6 files provided in the OPTIONS block
    call get_isize('TAS6_FILENAME', mf6_input%mempath, isize)
    if (isize > 0) then
      this%tas_active = .true.
      call mem_setptr(tas_fnames, 'TAS6_FILENAME', mf6_input%mempath)
      do n = 1, size(tas_fnames)
        fname = tas_fnames(n)
        call this%tasmanager%add_tasfile(fname)
      end do
    end if

    ! initialize input context memory
    call this%ctx%init(mf6_input)

    ! read state of each auxiliary variable
    call mem_allocate(this%aux_reads, this%ctx%naux, 'INAUXVARS', &
                      mf6_input%mempath)
    this%aux_reads = RS_NONE

    ! set active user nodes
    call this%set_active()

    ! allocate dfn params
    call this%params_alloc()

    ! storage for the mode
    call this%mode_init()
  end subroutine ainit

  subroutine df(this)
    class(GridArrayLoadType), intent(inout) :: this
    call this%tasmanager%tasmanager_df()
  end subroutine df

  subroutine ts_advance(this)
    class(GridArrayLoadType), intent(inout) :: this
    call this%tasmanager%ad()
  end subroutine ts_advance

  subroutine rp(this, parser)
    use IdmLoggerModule, only: idm_log_header, idm_log_close
    class(GridArrayLoadType), intent(inout) :: this
    type(BlockParserType), pointer, intent(inout) :: parser
    logical(LGP) :: endOfBlock

    ! reset for this period
    call this%reset()

    ! log lst file header
    call idm_log_header(this%mf6_input%component_name, &
                        this%mf6_input%subcomponent_name, this%iout)

    ! read array block
    do
      call parser%GetNextLine(endOfBlock)
      if (endOfBlock) exit
      call this%line_load(parser)
    end do

    ! finish the block for the mode
    call this%block_end()

    ! update time-array series links
    call this%tas_links_update(parser%iuactive)

    ! close lst file log block
    call idm_log_close(this%mf6_input%component_name, &
                       this%mf6_input%subcomponent_name, this%iout)
  end subroutine rp

  !> @brief Load one array line of the period block
  !!
  !! The line gives an array name (an auxiliary variable by its name),
  !! then the array values, a NETCDF keyword, or a TIMEARRAYSERIES name.
  !<
  subroutine line_load(this, parser)
    use BlockParserModule, only: BlockParserType
    use InputDefinitionModule, only: InputParamDefinitionType
    use DefinitionSelectModule, only: get_param_definition_type
    use SourceCommonModule, only: ifind_charstr
    class(GridArrayLoadType), intent(inout) :: this
    type(BlockParserType), pointer, intent(inout) :: parser
    logical(LGP) :: netcdf, layered
    character(len=LINELENGTH) :: keyword, param_tag
    character(len=LENTIMESERIESNAME) :: tas_name
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: iaux

    netcdf = .false.
    layered = .false.

    ! read param_tag
    call parser%GetStringCaps(param_tag)

    ! any auxvar corresponds to the definition tag 'AUX'
    iaux = ifind_charstr(this%ctx%auxname_cst, param_tag)
    if (iaux > 0) then
      param_tag = 'AUX'
    else
      iaux = 0
    end if

    ! set input definition
    idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                     this%mf6_input%component_type, &
                                     this%mf6_input%subcomponent_type, &
                                     'PERIOD', param_tag, this%input_name)
    ! auxiliary arrays are specified by auxiliary variable name
    if (iaux == 0 .and. idt%datatype == 'DOUBLE2D') then
      write (errmsg, '(a)') &
        'Input error: '//trim(param_tag)//' is not a valid PERIOD block &
        &array name. Auxiliary variable arrays are specified by &
        &auxiliary variable name.'
      call store_error(errmsg)
      call store_error_filename(this%input_name)
    end if

    ! look for Layered, NetCDF and time-array series keywords
    call parser%GetStringCaps(keyword)
    if (keyword == 'LAYERED' .and. idt%layered) then
      layered = .true.
    else if (keyword == 'NETCDF') then
      netcdf = .true.
    else if (keyword == 'TIMEARRAYSERIES') then
      if (.not. idt%timeseries) then
        write (errmsg, '(a)') &
          'Input error: '//trim(param_tag)//' cannot be specified as a &
          &TIMEARRAYSERIES.'
        call store_error(errmsg)
        call store_error_filename(this%input_name)
      end if
      if (.not. this%tas_active) then
        write (errmsg, '(a)') &
          'Input error: '//trim(param_tag)//' TIMEARRAYSERIES requires a &
          &TAS6 file in the OPTIONS block.'
        call store_error(errmsg)
        call store_error_filename(this%input_name)
      end if
      call parser%GetStringCaps(tas_name)
      call this%tas_load(idt, iaux, tas_name, parser%iuactive)
      return
    end if

    ! read and load the parameter
    call this%param_load(parser, idt, this%mf6_input%mempath, layered, &
                         netcdf, iaux)
  end subroutine line_load

  subroutine destroy(this)
    class(GridArrayLoadType), intent(inout) :: this
    integer(I4B) :: i
    ! series links point at loader arrays and the context nbound
    call this%tasmanager%da()
    deallocate (this%tasmanager)
    nullify (this%tasmanager)
    call this%ctx%destroy()
    if (allocated(this%active)) deallocate (this%active)
    nullify (this%aux_reads)
    if (allocated(this%slots)) then
      do i = 1, size(this%slots)
        call this%slot_clear(i)
      end do
      deallocate (this%slots)
    end if
    if (allocated(this%blk)) deallocate (this%blk)
  end subroutine destroy

  !> @brief Set active user nodes from discretization IDOMAIN input
  !!
  !! Cells with IDOMAIN less than one are inactive. All cells are active
  !! if IDOMAIN is not given.
  !<
  subroutine set_active(this)
    use ConstantsModule, only: LENMEMPATH, DIS, DISV, DISU
    use MemoryManagerModule, only: get_isize
    use MemoryHelperModule, only: create_mem_path
    use SimVariablesModule, only: idm_context
    class(GridArrayLoadType), intent(inout) :: this
    character(len=LENMEMPATH) :: model_mempath, dis_mempath
    integer(I4B), pointer :: disenum
    integer(I4B), dimension(:), pointer, contiguous :: int1d
    integer(I4B), dimension(:, :), pointer, contiguous :: int2d
    integer(I4B), dimension(:, :, :), pointer, contiguous :: int3d
    integer(I4B) :: isize

    allocate (this%active(this%ctx%nodes))
    this%active = .true.

    ! set discretization input mempath
    model_mempath = create_mem_path(component=this%mf6_input%component_name, &
                                    context=idm_context)
    call mem_setptr(disenum, 'DISENUM', model_mempath)
    select case (disenum)
    case (DIS)
      dis_mempath = create_mem_path(this%mf6_input%component_name, 'DIS', &
                                    idm_context)
    case (DISV)
      dis_mempath = create_mem_path(this%mf6_input%component_name, 'DISV', &
                                    idm_context)
    case (DISU)
      dis_mempath = create_mem_path(this%mf6_input%component_name, 'DISU', &
                                    idm_context)
    case default
      return
    end select

    ! IDOMAIN is optional; read before the discretization package
    ! releases its input
    call get_isize('IDOMAIN', dis_mempath, isize)
    if (isize /= this%ctx%nodes) return

    select case (disenum)
    case (DIS)
      call mem_setptr(int3d, 'IDOMAIN', dis_mempath)
      this%active = reshape(int3d, [this%ctx%nodes]) > 0
    case (DISV)
      call mem_setptr(int2d, 'IDOMAIN', dis_mempath)
      this%active = reshape(int2d, [this%ctx%nodes]) > 0
    case (DISU)
      call mem_setptr(int1d, 'IDOMAIN', dis_mempath)
      this%active = int1d > 0
    end select
  end subroutine set_active

  subroutine reset(this)
    class(GridArrayLoadType), intent(inout) :: this
    integer(I4B) :: n

    do n = 1, this%nreads
      call this%slot_clear(this%blk(n))
    end do
    this%nreads = 0

    do n = 1, this%nparam
      ! reset read state
      this%param_reads(n)%invar = RS_NONE
    end do
    this%aux_reads = RS_NONE
  end subroutine reset

  subroutine params_alloc(this)
    use DefinitionSelectModule, only: get_param_definition_type
    class(GridArrayLoadType), intent(inout) :: this
    type(InputParamDefinitionType), pointer :: idt
    character(len=LENVARNAME) :: rs_varname
    character(len=LENVARNAME) :: auxvar
    integer(I4B), pointer :: intvar
    integer(I4B) :: iparam, iaux, iparam_aux, islot

    ! set in-scope param names directly from context
    this%param_names = this%ctx%params
    this%nparam = size(this%ctx%params)
    call this%ctx%allocate_params()
    call this%ctx%check_developmode(this%input_name)
    call this%ctx%allocate_arrays()

    ! allocate and set param_reads pointer array
    allocate (this%param_reads(this%nparam))
    do iparam = 1, this%nparam
      ! allocate and store name of read state variable
      rs_varname = this%ctx%rsv_alloc(this%param_names(iparam))
      call mem_setptr(intvar, rs_varname, this%mf6_input%mempath)
      this%param_reads(iparam)%invar => intvar
      this%param_reads(iparam)%invar = RS_NONE
    end do

    ! count main arrays (1d); the AUX param holds the aux arrays (2d)
    this%nmain = 0
    iparam_aux = 0
    auxvar = ''
    do iparam = 1, this%nparam
      idt => this%param_dfn(iparam)
      if (idt%datatype == 'DOUBLE1D') then
        this%nmain = this%nmain + 1
      else
        iparam_aux = iparam
        auxvar = idt%mf6varname
      end if
    end do

    ! period block arrays: main arrays, then auxiliary variables
    allocate (this%slots(this%nmain + this%ctx%naux))
    allocate (this%blk(size(this%slots)))
    islot = 0
    do iparam = 1, this%nparam
      idt => this%param_dfn(iparam)
      if (idt%datatype /= 'DOUBLE1D') cycle
      islot = islot + 1
      this%slots(islot)%tagname = this%param_names(iparam)
      this%slots(islot)%mf6varname = idt%mf6varname
      this%slots(islot)%iparam = iparam
    end do
    do iaux = 1, this%ctx%naux
      islot = this%nmain + iaux
      this%slots(islot)%tagname = this%ctx%auxname_cst(iaux)
      this%slots(islot)%mf6varname = auxvar
      this%slots(islot)%iparam = iparam_aux
      this%slots(islot)%iaux = iaux
    end do
  end subroutine params_alloc

  !> @brief PERIOD block definition of a dynamic param
  !<
  function param_dfn(this, iparam) result(idt)
    use DefinitionSelectModule, only: get_param_definition_type
    class(GridArrayLoadType), intent(in) :: this
    integer(I4B), intent(in) :: iparam
    type(InputParamDefinitionType), pointer :: idt
    idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                     this%mf6_input%component_type, &
                                     this%mf6_input%subcomponent_type, &
                                     'PERIOD', this%param_names(iparam), &
                                     this%input_name)
  end function param_dfn

  !> @brief Slot of a block array: a main array by tag, an aux array by
  !! aux index
  !<
  function slot_index(this, tagname, iaux) result(islot)
    class(GridArrayLoadType), intent(in) :: this
    character(len=*), intent(in) :: tagname
    integer(I4B), intent(in) :: iaux
    integer(I4B) :: islot
    if (iaux > 0) then
      islot = this%nmain + iaux
      return
    end if
    do islot = 1, this%nmain
      if (this%slots(islot)%tagname == tagname) return
    end do
    errmsg = 'Programming error: '//trim(tagname)//' is not a grid array &
             &period block array.'
    call store_error(errmsg, terminate=.true.)
  end function slot_index

  !> @brief Clear the block state of a slot
  !<
  subroutine slot_clear(this, islot)
    class(GridArrayLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: islot
    this%slots(islot)%inblock = .false.
  end subroutine slot_clear

  subroutine param_load(this, parser, idt, mempath, layered, netcdf, iaux)
    use TdisModule, only: kper
    use InputDefinitionModule, only: InputParamDefinitionType
    use Double1dReaderModule, only: read_dbl1d
    use LayeredArrayReaderModule, only: read_dbl1d_layered
    use LoadNCInputModule, only: netcdf_read_array
    use SourceCommonModule, only: get_layered_shape
    use IdmLoggerModule, only: idm_log_var
    class(GridArrayLoadType), intent(inout) :: this
    type(BlockParserType), intent(in) :: parser
    type(InputParamDefinitionType), intent(in) :: idt
    character(len=*), intent(in) :: mempath
    logical(LGP), intent(in) :: layered
    logical(LGP), intent(in) :: netcdf
    integer(I4B), intent(in) :: iaux
    real(DP), dimension(:), pointer, contiguous :: nodes
    integer(I4B), dimension(:), allocatable :: layer_shape
    integer(I4B) :: islot, nlay
    integer(I4B) :: rs

    islot = this%slot_index(idt%tagname, iaux)

    select case (idt%datatype)
    case ('DOUBLE1D', 'DOUBLE2D')
      allocate (nodes(this%ctx%nodes))
      if (netcdf) then
        if (iaux > 0) then
          call netcdf_read_array(nodes, this%ctx%mshape, idt, &
                                 this%mf6_input, this%nc_vars, &
                                 this%input_name, this%iout, kper, iaux)
        else
          call netcdf_read_array(nodes, this%ctx%mshape, idt, &
                                 this%mf6_input, this%nc_vars, &
                                 this%input_name, this%iout, kper)
        end if
      else if (layered) then
        call get_layered_shape(this%ctx%mshape, nlay, layer_shape)
        call read_dbl1d_layered(parser, nodes, idt%mf6varname, nlay, layer_shape)
      else
        call read_dbl1d(parser, nodes, idt%mf6varname)
      end if

      call idm_log_var(nodes, idt%tagname, mempath, this%iout)

      ! add to the period block arrays
      call this%array_add(islot, nodes)
      deallocate (nodes)
    case default
      errmsg = 'IDM unimplemented. GridArrayLoad::param_load &
               &datatype='//trim(idt%datatype)
      call store_error(errmsg)
      call store_error_filename(this%input_name)
    end select

    ! set read state
    rs = merge(RS_NETCDF, RS_VALUES, netcdf)
    this%param_reads(this%slots(islot)%iparam)%invar = rs
    if (iaux > 0) this%aux_reads(iaux) = rs
  end subroutine param_load

  !> @brief Add an array read in the period block
  !!
  !! nodes is not associated for a series array whose values are set by
  !! the time steps.
  !<
  subroutine array_add(this, islot, nodes)
    class(GridArrayLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: islot
    real(DP), dimension(:), pointer, contiguous, intent(in) :: nodes

    ! an array may be given only once in a block
    if (this%slots(islot)%inblock) then
      write (errmsg, '(a)') &
        'Input error: '//trim(this%slots(islot)%tagname)// &
        ' is specified more than once in the PERIOD block.'
      call store_error(errmsg)
      call store_error_filename(this%input_name)
    end if

    ! add to the block arrays
    this%nreads = this%nreads + 1
    this%blk(this%nreads) = islot
    this%slots(islot)%inblock = .true.
    this%slots(islot)%given = .true.
    this%slots(islot)%tasname = ''

    ! store for the mode
    call this%array_store(islot, nodes)
  end subroutine array_add

  !> @brief Add an array given as a time-array series
  !!
  !! The mode gives the array to add at the period start, if any. The
  !! link is made at the end of the block.
  !<
  subroutine tas_load(this, idt, iaux, tas_name, inunit)
    use IdmLoggerModule, only: idm_log_var
    class(GridArrayLoadType), intent(inout) :: this
    type(InputParamDefinitionType), intent(in) :: idt
    integer(I4B), intent(in) :: iaux
    character(len=*), intent(in) :: tas_name
    integer(I4B), intent(in) :: inunit
    type(TimeArraySeriesType), pointer :: tasptr
    real(DP), dimension(:), pointer, contiguous :: nodes
    integer(I4B) :: islot

    islot = this%slot_index(idt%tagname, iaux)

    tasptr => this%tasmanager%get_series(tas_name, inunit)
    nodes => this%tas_start_array(tasptr)

    ! log variable
    call idm_log_var(this%slots(islot)%tagname, this%mf6_input%mempath, &
                     this%iout, .true., trim(tas_name))

    ! add to the period block arrays
    call this%array_add(islot, nodes)
    this%slots(islot)%tasname = tas_name

    ! set read state: time series input
    if (iaux > 0) then
      this%aux_reads(iaux) = RS_SERIES
    else
      this%param_reads(this%slots(islot)%iparam)%invar = RS_SERIES
    end if
  end subroutine tas_load

  !> @brief Set up storage for the mode, after the params are allocated
  !<
  subroutine mode_init(this)
    class(GridArrayLoadType), intent(inout) :: this
    call mode_required('mode_init')
  end subroutine mode_init

  !> @brief Store an array added to the period block
  !!
  !! nodes is not associated for a series array set by the time steps.
  !<
  subroutine array_store(this, islot, nodes)
    class(GridArrayLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: islot
    real(DP), dimension(:), pointer, contiguous, intent(in) :: nodes
    call mode_required('array_store')
  end subroutine array_store

  !> @brief Update time-array series links after the period block
  !<
  subroutine tas_links_update(this, inunit)
    class(GridArrayLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: inunit
    call mode_required('tas_links_update')
  end subroutine tas_links_update

  !> @brief Finish the period block after all arrays are read
  !!
  !! Nothing to do by default: arrays are stored as they are read.
  !<
  subroutine block_end(this)
    class(GridArrayLoadType), intent(inout) :: this
  end subroutine block_end

  !> @brief Array to add for a series given in the block
  !!
  !! None by default: the time steps set the values.
  !<
  function tas_start_array(this, tasptr) result(nodes)
    class(GridArrayLoadType), intent(inout) :: this
    type(TimeArraySeriesType), pointer, intent(in) :: tasptr
    real(DP), dimension(:), pointer, contiguous :: nodes
    nullify (nodes)
  end function tas_start_array

  !> @brief Stop for a hook a storage mode must provide
  !<
  subroutine mode_required(name)
    character(len=*), intent(in) :: name
    errmsg = 'Programming error: grid array storage mode does not provide '// &
             trim(name)//'.'
    call store_error(errmsg, terminate=.true.)
  end subroutine mode_required

  !> @brief Link every array given as a series
  !!
  !! Element i of a linked array takes series element imap(i), for the
  !! first NBOUND elements. With keep_nodata, a series DNODATA value or
  !! imap(i) < 1 leaves the element unchanged.
  !<
  subroutine make_links(this, inunit, imap, keep_nodata)
    class(GridArrayLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: inunit
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: imap
    logical(LGP), intent(in) :: keep_nodata
    real(DP), dimension(:), pointer :: bndArrayPtr
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    integer(I4B), dimension(:), pointer, contiguous :: nodelist
    character(len=LINELENGTH) :: text
    integer(I4B) :: islot, iaux

    call this%tasmanager%remove_links(this%mf6_input%subcomponent_name)
    nullify (nodelist)
    do islot = 1, size(this%slots)
      if (this%slots(islot)%tasname == '') cycle
      iaux = this%slots(islot)%iaux
      if (iaux > 0) then
        text = this%slots(islot)%tagname
        bndArrayPtr => this%ctx%auxvar(iaux, :)
      else
        text = this%slots(islot)%mf6varname
        call mem_setptr(dbl1d, this%slots(islot)%mf6varname, &
                        this%mf6_input%mempath)
        bndArrayPtr => dbl1d
      end if
      call this%tasmanager%MakeTasLink(this%mf6_input%subcomponent_name, &
                                       bndArrayPtr, this%ctx%iprpak, &
                                       this%slots(islot)%tasname, text, &
                                       .false., nodelist, inunit, &
                                       imap=imap, nmap=this%ctx%nbound, &
                                       keep_nodata=keep_nodata)
    end do
  end subroutine make_links

  !> @brief Return true if the package input defines MAXBOUND
  !<
  function defines_maxbound(mf6_input) result(found)
    type(ModflowInputType), intent(in) :: mf6_input
    logical(LGP) :: found
    integer(I4B) :: n
    found = .false.
    do n = 1, size(mf6_input%param_dfns)
      if (mf6_input%param_dfns(n)%blockname == 'DIMENSIONS' .and. &
          mf6_input%param_dfns(n)%tagname == 'MAXBOUND') then
        found = .true.
        exit
      end if
    end do
  end function defines_maxbound

end module GridArrayLoadModule
