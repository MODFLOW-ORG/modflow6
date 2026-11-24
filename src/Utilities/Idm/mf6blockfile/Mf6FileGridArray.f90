!> @brief This module contains the GridArrayLoadModule
!!
!! This module contains the routines for reading period block
!! grid array-based input for stress packages that use the
!! READARRAYGRID option (CHD, WEL, DRN, RIV, GHB).
!!
!<
module GridArrayLoadModule

  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LINELENGTH, LENVARNAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use InputDefinitionModule, only: InputParamDefinitionType
  use MemoryManagerModule, only: mem_allocate, mem_setptr
  use BlockParserModule, only: BlockParserType
  use ModflowInputModule, only: ModflowInputType
  use LoadContextModule, only: LoadContextType, ReadStateVarType, RS_NONE, &
                               RS_VALUES, RS_NETCDF
  use AsciiInputLoadTypeModule, only: AsciiDynamicPkgLoadBaseType

  implicit none
  private
  public :: GridArrayLoadType

  !> @brief Ascii grid based dynamic loader type
  !<
  type, extends(AsciiDynamicPkgLoadBaseType) :: GridArrayLoadType
    type(ReadStateVarType), dimension(:), allocatable :: param_reads !< read states for current load
    type(LoadContextType) :: ctx
    integer(I4B), dimension(:), pointer, contiguous :: nodeulist
    logical(LGP), dimension(:), allocatable :: active !< user node IDOMAIN > 0
    logical(LGP) :: nodeulist_set = .false. !< nodeulist set in this block
    integer(I4B), dimension(:), pointer, contiguous :: aux_reads => null() !< read state of each auxiliary variable
  contains
    procedure :: ainit
    procedure :: df
    procedure :: ts_advance
    procedure :: rp
    procedure :: destroy
    procedure :: reset
    procedure :: params_alloc
    procedure :: param_load
    procedure :: set_active
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
    integer(I4B) :: isize
    integer(I4B), pointer :: maxbound

    ! initialize base type
    call this%DynamicPkgLoadType%init(mf6_input, component_name, &
                                      component_input_name, &
                                      input_name, iperblock, iout)
    ! initialize
    this%iout = iout

    ! load static input
    call loader%load(parser, mf6_input, this%nc_vars, this%input_name, iout)

    ! maxbound is optional
    call get_isize('MAXBOUND', mf6_input%mempath, isize)
    if (isize < 0) then
      ! set maxbound to grid nodes
      call mem_allocate(maxbound, 'MAXBOUND', mf6_input%mempath)
      maxbound = product(loader%mshape)
    end if

    ! initialize input context memory
    call this%ctx%init(mf6_input)

    ! set active user nodes
    call this%set_active()

    ! allocate user nodelist
    call mem_allocate(this%nodeulist, this%ctx%maxbound, &
                      'NODEULIST', mf6_input%mempath)

    ! allocate dfn params
    call this%params_alloc()
  end subroutine ainit

  subroutine df(this)
    class(GridArrayLoadType), intent(inout) :: this
  end subroutine df

  subroutine ts_advance(this)
    class(GridArrayLoadType), intent(inout) :: this
  end subroutine ts_advance

  subroutine rp(this, parser)
    use BlockParserModule, only: BlockParserType
    use InputDefinitionModule, only: InputParamDefinitionType
    use DefinitionSelectModule, only: get_param_definition_type
    use ArrayHandlersModule, only: ifind
    use SourceCommonModule, only: ifind_charstr
    use IdmLoggerModule, only: idm_log_header, idm_log_close
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
      iaux = ifind_charstr(this%ctx%auxname_cst, param_tag)

      ! any auxvar corresponds to the definition tag 'AUX'
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
    use MemoryManagerModule, only: mem_deallocate
    class(GridArrayLoadType), intent(inout) :: this
    call mem_deallocate(this%nodeulist)
    if (allocated(this%active)) deallocate (this%active)
  end subroutine destroy

  !> @brief Set active user nodes from discretization IDOMAIN input
  !!
  !! Cells with IDOMAIN less than one are not boundaries. All cells are
  !! active if IDOMAIN is not given.
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

    this%ctx%nbound = 0
    this%nodeulist_set = .false.

    do n = 1, this%nparam
      ! reset read state
      this%param_reads(n)%invar = RS_NONE
    end do
    this%aux_reads = RS_NONE
  end subroutine reset

  subroutine params_alloc(this)
    class(GridArrayLoadType), intent(inout) :: this
    character(len=LENVARNAME) :: rs_varname
    integer(I4B), pointer :: intvar
    integer(I4B) :: iparam

    ! set in-scope param names directly from context
    this%param_names = this%ctx%params
    this%nparam = size(this%ctx%params)
    call this%ctx%allocate_params()
    call this%ctx%check_developmode(this%input_name)
    call this%ctx%allocate_arrays()

    ! allocate and set param_reads pointer array
    allocate (this%param_reads(this%nparam))

    ! store read state variable pointers
    do iparam = 1, this%nparam
      ! allocate and store name of read state variable
      rs_varname = this%ctx%rsv_alloc(this%param_names(iparam))
      call mem_setptr(intvar, rs_varname, this%mf6_input%mempath)
      this%param_reads(iparam)%invar => intvar
      this%param_reads(iparam)%invar = RS_NONE
    end do

    ! read state of each auxiliary variable
    call mem_allocate(this%aux_reads, this%ctx%naux, 'INAUXVARS', &
                      this%mf6_input%mempath)
    this%aux_reads = RS_NONE
  end subroutine params_alloc

  subroutine param_load(this, parser, idt, mempath, layered, netcdf, iaux)
    use TdisModule, only: kper
    use ConstantsModule, only: DNODATA
    use ArrayHandlersModule, only: ifind
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
    real(DP), dimension(:), pointer, contiguous :: dbl1d, nodes
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B), dimension(:), allocatable :: layer_shape
    integer(I4B) :: iparam, n, nlay, nnode

    nnode = 0

    select case (idt%datatype)
    case ('DOUBLE1D')
      call mem_setptr(dbl1d, idt%mf6varname, mempath)
      allocate (nodes(this%ctx%nodes))
      if (netcdf) then
        call netcdf_read_array(nodes, this%ctx%mshape, idt, &
                               this%mf6_input, this%nc_vars, this%input_name, &
                               this%iout, kper)
      else if (layered) then
        call get_layered_shape(this%ctx%mshape, nlay, layer_shape)
        call read_dbl1d_layered(parser, nodes, idt%mf6varname, nlay, layer_shape)
      else
        call read_dbl1d(parser, nodes, idt%mf6varname)
      end if

      call idm_log_var(nodes, idt%tagname, mempath, this%iout)

      if (this%nodeulist_set) then
        ! nodeulist already established: extract values at known positions
        do n = 1, this%ctx%nbound
          dbl1d(n) = nodes(this%nodeulist(n))
        end do
      else
        ! first array sets nodeulist and nbound (defined, active cells)
        do n = 1, this%ctx%nodes
          if (nodes(n) /= DNODATA .and. this%active(n)) then
            nnode = nnode + 1
            if (nnode > this%ctx%maxbound) then
              write (errmsg, '(a,i0,a)') &
                'Input error: number of defined (non-DNODATA) active cells &
                &exceeds MAXBOUND=', this%ctx%maxbound, '.'
              call store_error(errmsg)
              call store_error_filename(this%input_name)
            end if
            dbl1d(nnode) = nodes(n)
            this%nodeulist(nnode) = n
          end if
        end do
        this%ctx%nbound = nnode
        this%nodeulist_set = .true.
      end if
      deallocate (nodes)
    case ('DOUBLE2D')
      call mem_setptr(dbl2d, idt%mf6varname, mempath)
      allocate (nodes(this%ctx%nodes))

      if (netcdf) then
        call netcdf_read_array(nodes, this%ctx%mshape, idt, &
                               this%mf6_input, this%nc_vars, this%input_name, &
                               this%iout, kper, iaux)
      else if (layered) then
        call get_layered_shape(this%ctx%mshape, nlay, layer_shape)
        call read_dbl1d_layered(parser, nodes, idt%mf6varname, nlay, layer_shape)
      else
        call read_dbl1d(parser, nodes, idt%mf6varname)
      end if

      call idm_log_var(nodes, idt%tagname, mempath, this%iout)

      if (this%nodeulist_set) then
        ! nodeulist already established: extract values at known positions
        do n = 1, this%ctx%nbound
          dbl2d(iaux, n) = nodes(this%nodeulist(n))
        end do
      else
        ! first array sets nodeulist and nbound (defined, active cells)
        do n = 1, this%ctx%nodes
          if (nodes(n) /= DNODATA .and. this%active(n)) then
            nnode = nnode + 1
            if (nnode > this%ctx%maxbound) then
              write (errmsg, '(a,i0,a)') &
                'Input error: number of defined (non-DNODATA) active cells &
                &exceeds MAXBOUND=', this%ctx%maxbound, '.'
              call store_error(errmsg)
              call store_error_filename(this%input_name)
            end if
            dbl2d(iaux, nnode) = nodes(n)
            this%nodeulist(nnode) = n
          end if
        end do
        this%ctx%nbound = nnode
        this%nodeulist_set = .true.
      end if
      deallocate (nodes)
      this%aux_reads(iaux) = merge(RS_NETCDF, RS_VALUES, netcdf)
    case default
      errmsg = 'IDM unimplemented. GridArrayLoad::param_load &
               &datatype='//trim(idt%datatype)
      call store_error(errmsg)
      call store_error_filename(this%input_name)
    end select

    ! if param is tracked set read state
    iparam = ifind(this%param_names, idt%tagname)
    if (iparam > 0) then
      this%param_reads(iparam)%invar = merge(RS_NETCDF, RS_VALUES, netcdf)
    end if
  end subroutine param_load

end module GridArrayLoadModule
