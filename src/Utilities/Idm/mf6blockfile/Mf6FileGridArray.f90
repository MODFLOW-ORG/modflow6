!> @brief This module contains the GridArrayLoadModule
!!
!! This module contains the routines for reading period block
!! array based input associated with the full grid, such as
!! with the GHBA package.
!!
!<
module GridArrayLoadModule

  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: DZERO, IZERO, LINELENGTH, LENVARNAME, &
                             LENTIMESERIESNAME, LENAUXNAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use InputDefinitionModule, only: InputParamDefinitionType
  use MemoryManagerModule, only: mem_allocate, mem_reallocate, mem_setptr
  use CharacterStringModule, only: CharacterStringType
  use BlockParserModule, only: BlockParserType
  use ModflowInputModule, only: ModflowInputType, getModflowInput
  use BoundInputContextModule, only: BoundInputContextType, ReadStateVarType
  use AsciiInputLoadTypeModule, only: AsciiDynamicPkgLoadBaseType

  implicit none
  private
  public :: GridArrayLoadType

  !> @brief Ascii grid based dynamic loader type
  !<
  type, extends(AsciiDynamicPkgLoadBaseType) :: GridArrayLoadType
    type(ReadStateVarType), dimension(:), allocatable :: param_reads !< read states for current load
    type(BoundInputContextType) :: bound_context
  contains
    procedure :: ainit
    procedure :: df
    procedure :: ad
    procedure :: rp
    procedure :: destroy
    procedure :: reset
    procedure :: params_alloc
    procedure :: param_load
  end type GridArrayLoadType

contains

  subroutine ainit(this, mf6_input, component_name, &
                   component_input_name, input_name, &
                   iperblock, parser, iout)
    use MemoryManagerModule, only: get_isize
    use BlockParserModule, only: BlockParserType
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

    ! initialize base type
    call this%DynamicPkgLoadType%init(mf6_input, component_name, &
                                      component_input_name, &
                                      input_name, iperblock, iout)
    ! initialize
    this%iout = iout

    ! load static input
    call loader%load(parser, mf6_input, this%nc_vars, this%input_name, iout)

    ! initialize input context memory
    call this%bound_context%create(mf6_input, this%readarray_grid)

    ! allocate dfn params
    call this%params_alloc()
  end subroutine ainit

  subroutine df(this)
    class(GridArrayLoadType), intent(inout) :: this
  end subroutine df

  subroutine ad(this)
    class(GridArrayLoadType), intent(inout) :: this
  end subroutine ad

  subroutine rp(this, parser)
    use MemoryManagerModule, only: mem_setptr
    use BlockParserModule, only: BlockParserType
    use InputDefinitionModule, only: InputParamDefinitionType
    use DefinitionSelectModule, only: get_param_definition_type
    use ArrayHandlersModule, only: ifind
    use SourceCommonModule, only: ifind_charstr
    use IdmLoggerModule, only: idm_log_header, idm_log_close, idm_log_var
    class(GridArrayLoadType), intent(inout) :: this
    type(BlockParserType), pointer, intent(inout) :: parser
    logical(LGP) :: endOfBlock, netcdf, layered
    character(len=LINELENGTH) :: keyword, param_tag
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: iaux

    ! reset for this period
    call this%reset()

    ! log lst file header
    call idm_log_header(this%mf6_input%component_name, &
                        this%mf6_input%subcomponent_name, this%iout)

    ! read array block
    do
      ! initialize
      iaux = 0
      netcdf = .false.
      layered = .false.

      ! read next line
      call parser%GetNextLine(endOfBlock)
      if (endOfBlock) exit
      ! read param_tag
      call parser%GetStringCaps(param_tag)

      ! is param tag an auxvar?
      iaux = ifind_charstr(this%bound_context%auxname_cst, param_tag)

      ! any auvxar corresponds to the definition tag 'AUX'
      if (iaux > 0) param_tag = 'AUX'

      ! set input definition
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       'PERIOD', param_tag, this%input_name)
      ! look for Layered and NetCDF keywords
      call parser%GetStringCaps(keyword)
      if (keyword == 'LAYERED' .and. idt%layered) then
        layered = .true.
      else if (keyword == 'NETCDF') then
        netcdf = .true.
      end if

      ! read and load the parameter
      call this%param_load(parser, idt, this%mf6_input%mempath, layered, &
                           netcdf, iaux)
    end do

    ! log lst file header
    call idm_log_close(this%mf6_input%component_name, &
                       this%mf6_input%subcomponent_name, this%iout)
  end subroutine rp

  subroutine destroy(this)
    class(GridArrayLoadType), intent(inout) :: this
  end subroutine destroy

  subroutine reset(this)
    class(GridArrayLoadType), intent(inout) :: this
    integer(I4B) :: n, m

    do n = 1, this%nparam
      ! reset read state
      this%param_reads(n)%invar = 0
    end do

    ! explicitly reset auxvar array each period
    do m = 1, this%bound_context%nodes
      do n = 1, this%bound_context%naux
        this%bound_context%auxvar(n, m) = DZERO
      end do
    end do
  end subroutine reset

  subroutine params_alloc(this)
    class(GridArrayLoadType), intent(inout) :: this
    character(len=LENVARNAME) :: rs_varname
    integer(I4B), pointer :: intvar
    integer(I4B) :: iparam

    ! set in scope param names
    call this%bound_context%bound_params(this%param_names, this%nparam, &
                                         this%input_name)
    call this%bound_context%allocate_arrays()

    ! allocate and set param_reads pointer array
    allocate (this%param_reads(this%nparam))

    ! store read state variable pointers
    do iparam = 1, this%nparam
      ! allocate and store name of read state variable
      rs_varname = this%bound_context%rsv_alloc(this%param_names(iparam))
      call mem_setptr(intvar, rs_varname, this%mf6_input%mempath)
      this%param_reads(iparam)%invar => intvar
      this%param_reads(iparam)%invar = 0
    end do
  end subroutine params_alloc

  subroutine param_load(this, parser, idt, mempath, layered, netcdf, iaux)
    use TdisModule, only: kper
    use ConstantsModule, only: DNODATA
    use MemoryManagerModule, only: mem_setptr
    use ArrayHandlersModule, only: ifind
    use InputDefinitionModule, only: InputParamDefinitionType
    use DefinitionSelectModule, only: get_param_definition_type
    use Double1dReaderModule, only: read_dbl1d
    use Double2dReaderModule, only: read_dbl2d
    use Integer1dReaderModule, only: read_int1d
    use LayeredArrayReaderModule, only: read_dbl1d_layered, &
                                        read_int1d_layered
    use LoadNCInputModule, only: netcdf_read_array
    use SourceCommonModule, only: get_shape_from_string, get_layered_shape
    use IdmLoggerModule, only: idm_log_var
    class(GridArrayLoadType), intent(inout) :: this
    type(BlockParserType), intent(in) :: parser
    type(InputParamDefinitionType), intent(in) :: idt
    character(len=*), intent(in) :: mempath
    logical(LGP), intent(in) :: layered
    logical(LGP), intent(in) :: netcdf
    integer(I4B), dimension(:), pointer, contiguous :: int1d
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B), dimension(:), allocatable :: layer_shape
    integer(I4B) :: iaux, iparam, n, nlay

    select case (idt%datatype)
    case ('INTEGER1D')
      call mem_setptr(int1d, idt%mf6varname, mempath)
      if (netcdf) then
        call netcdf_read_array(int1d, this%bound_context%mshape, idt, &
                               this%mf6_input, this%nc_vars, this%input_name, &
                               this%iout, kper)
      else if (layered) then
        call get_layered_shape(this%bound_context%mshape, nlay, layer_shape)
        call read_int1d_layered(parser, int1d, idt%mf6varname, nlay, layer_shape)
      else
        call read_int1d(parser, int1d, idt%mf6varname)
      end if
      call idm_log_var(int1d, idt%tagname, mempath, this%iout)
    case ('DOUBLE1D')
      ! set pointer to managed memory input variable
      call mem_setptr(dbl1d, idt%mf6varname, mempath)

      ! read user input
      if (netcdf) then
        call netcdf_read_array(dbl1d, this%bound_context%mshape, idt, &
                               this%mf6_input, this%nc_vars, this%input_name, &
                               this%iout, kper)
      else if (layered) then
        call get_layered_shape(this%bound_context%mshape, nlay, layer_shape)
        call read_dbl1d_layered(parser, dbl1d, idt%mf6varname, nlay, layer_shape)
      else
        call read_dbl1d(parser, dbl1d, idt%mf6varname)
      end if

      ! log user input
      call idm_log_var(dbl1d, idt%tagname, mempath, this%iout)
    case ('DOUBLE2D')
      ! set pointer to managed memory input variable
      call mem_setptr(dbl2d, idt%mf6varname, mempath)

      ! allocate local array
      allocate (dbl1d(this%bound_context%nodes))

      ! read user input
      if (netcdf) then
        call netcdf_read_array(dbl1d, this%bound_context%mshape, idt, &
                               this%mf6_input, this%nc_vars, this%input_name, &
                               this%iout, kper, iaux)
      else if (layered) then
        call get_layered_shape(this%bound_context%mshape, nlay, layer_shape)
        call read_dbl1d_layered(parser, dbl1d, idt%mf6varname, nlay, layer_shape)
      else
        call read_dbl1d(parser, dbl1d, idt%mf6varname)
      end if

      ! copy into 2d array
      do n = 1, this%bound_context%nodes
        dbl2d(iaux, n) = dbl1d(n)
      end do

      ! log user input
      call idm_log_var(dbl1d, idt%tagname, mempath, this%iout)

      ! cleanup
      deallocate (dbl1d)
    case default
      errmsg = 'IDM unimplemented. GridArrayLoad::param_load &
               &datatype='//trim(idt%datatype)
      call store_error(errmsg)
      call store_error_filename(this%input_name)
    end select

    ! if param is tracked set read state
    iparam = ifind(this%param_names, idt%tagname)
    if (iparam > 0) then
      this%param_reads(iparam)%invar = 1
    end if
  end subroutine param_load

end module GridArrayLoadModule
