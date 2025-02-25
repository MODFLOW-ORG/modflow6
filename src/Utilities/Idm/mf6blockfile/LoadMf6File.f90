!> @brief This module contains the LoadMf6FileModule
!!
!! This module contains the input data model routines for
!! loading static data from a MODFLOW 6 input file using the
!! block parser.
!!
!<
module LoadMf6FileModule

  use KindModule, only: DP, I4B, LGP
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use ConstantsModule, only: LINELENGTH, LENVARNAME
  use BlockParserModule, only: BlockParserType
  use LayeredArrayReaderModule, only: read_dbl1d_layered, &
                                      read_dbl2d_layered, &
                                      read_dbl3d_layered, &
                                      read_int1d_layered, &
                                      read_int2d_layered, &
                                      read_int3d_layered
  use Double1dReaderModule, only: read_dbl1d
  use Double2dReaderModule, only: read_dbl2d
  use Integer1dReaderModule, only: read_int1d
  use Integer2dReaderModule, only: read_int2d
  use InputOutputModule, only: parseline
  use InputDefinitionModule, only: InputParamDefinitionType
  use DefinitionSelectModule, only: get_param_definition_type, &
                                    get_aggregate_definition_type
  use ModflowInputModule, only: ModflowInputType, getModflowInput
  use MemoryManagerModule, only: mem_allocate, mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use StructArrayModule, only: StructArrayType
  use NCFileVarsModule, only: NCPackageVarsType
  use IdmLoggerModule, only: idm_log_var, idm_log_header, idm_log_close, &
                             idm_export

  implicit none
  private
  public :: LoadMf6FileType
  public :: read_control_record

  !> @brief Static parser based input loader
  !!
  !! This type defines a static input context loader
  !! for traditional mf6 ascii input files.
  !!
  !<
  type :: LoadMf6FileType
    type(BlockParserType), pointer :: parser !< ascii block parser
    integer(I4B), dimension(:), pointer, contiguous :: mshape => null() !< model shape
    type(StructArrayType), pointer :: structarray => null() !< structarray for loading list input
    type(ModflowInputType) :: mf6_input !< description of input
    type(NCPackageVarsType), pointer :: nc_vars => null()
    character(len=LINELENGTH) :: filename !< name of ascii input file
    character(len=LINELENGTH), dimension(:), allocatable :: block_tags !< read block tags
    logical(LGP) :: ts_active !< is timeseries active
    logical(LGP) :: export !< is array export active
    logical(LGP) :: readasarrays
    integer(I4B) :: inamedbound
    integer(I4B) :: iauxiliary
    integer(I4B) :: iout !< inunit for list log
  contains
    procedure :: load
    procedure :: init
    procedure :: load_block
    procedure :: finalize
    procedure :: parse_block
    procedure :: block_post_process
    procedure :: parse_io_tag
    procedure :: parse_keyword_tag
    procedure :: parse_tag
    procedure :: block_index_dfn
    procedure :: parse_structarray_block
  end type LoadMf6FileType

contains

  !> @brief load all static input blocks
  !!
  !! Invoke this routine to load all static input blocks
  !! in single call.
  !!
  !<
  subroutine load(this, parser, mf6_input, nc_vars, filename, iout)
    use MemoryManagerModule, only: get_isize
    class(LoadMf6FileType) :: this
    type(BlockParserType), target, intent(inout) :: parser
    type(ModflowInputType), intent(in) :: mf6_input
    type(NCPackageVarsType), pointer, intent(in) :: nc_vars
    character(len=*), intent(in) :: filename
    integer(I4B), intent(in) :: iout
    integer(I4B) :: iblk

    ! initialize static load
    call this%init(parser, mf6_input, filename, iout)

    ! set netcdf vars
    this%nc_vars => nc_vars

    ! process blocks
    do iblk = 1, size(this%mf6_input%block_dfns)
      ! don't load dynamic input data
      if (this%mf6_input%block_dfns(iblk)%blockname == 'PERIOD') exit
      ! load the block
      call this%load_block(iblk)
    end do

    ! finalize static load
    call this%finalize()
  end subroutine load

  !> @brief init
  !!
  !! init / finalize are only used when load_block() will be called
  !!
  !<
  subroutine init(this, parser, mf6_input, filename, iout)
    use MemoryManagerModule, only: get_isize
    class(LoadMf6FileType) :: this
    type(BlockParserType), target, intent(inout) :: parser
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), intent(in) :: filename
    integer(I4B), intent(in) :: iout
    integer(I4B) :: isize

    this%parser => parser
    this%mf6_input = mf6_input
    this%filename = filename
    this%ts_active = .false.
    this%export = .false.
    this%readasarrays = .false.
    this%inamedbound = 0
    this%iauxiliary = 0
    this%iout = iout

    call get_isize('MODEL_SHAPE', mf6_input%component_mempath, isize)
    if (isize > 0) then
      call mem_setptr(this%mshape, 'MODEL_SHAPE', mf6_input%component_mempath)
    end if

    ! log lst file header
    call idm_log_header(this%mf6_input%component_name, &
                        this%mf6_input%subcomponent_name, this%iout)
  end subroutine init

  !> @brief load a single block
  !!
  !! Assumed in order load of single (next) block. If a
  !! StructArray object is allocated to load this block
  !! it persists until this routine (or finalize) is
  !! called again.
  !!
  !<
  subroutine load_block(this, iblk)
    use StructArrayModule, only: destructStructArray
    class(LoadMf6FileType) :: this
    integer(I4B), intent(in) :: iblk

    ! reset structarray if it was created for previous block
    if (associated(this%structarray)) then
      ! destroy the structured array reader
      call destructStructArray(this%structarray)
    end if

    allocate (this%block_tags(0))
    ! load the block
    call this%parse_block(iblk, .false.)
    ! post process block
    call this%block_post_process(iblk)
    ! cleanup
    deallocate (this%block_tags)
  end subroutine load_block

  !> @brief finalize
  !!
  !! init / finalize are only used when load_block() will be called
  !!
  !<
  subroutine finalize(this)
    use StructArrayModule, only: destructStructArray
    class(LoadMf6FileType) :: this
    ! cleanup
    if (associated(this%structarray)) then
      ! destroy the structured array reader
      call destructStructArray(this%structarray)
    end if
    ! close logging block
    call idm_log_close(this%mf6_input%component_name, &
                       this%mf6_input%subcomponent_name, this%iout)
  end subroutine finalize

  !> @brief Post parse block handling
  !!
  !<
  subroutine block_post_process(this, iblk)
    use ConstantsModule, only: LENBOUNDNAME
    use CharacterStringModule, only: CharacterStringType
    use SourceCommonModule, only: set_model_shape
    class(LoadMf6FileType) :: this
    integer(I4B), intent(in) :: iblk
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: iparam
    integer(I4B), pointer :: intptr

    ! update state based on read tags
    do iparam = 1, size(this%block_tags)
      select case (this%mf6_input%block_dfns(iblk)%blockname)
      case ('OPTIONS')
        if (this%block_tags(iparam) == 'AUXILIARY') then
          this%iauxiliary = 1
        else if (this%block_tags(iparam) == 'BOUNDNAMES') then
          this%inamedbound = 1
        else if (this%block_tags(iparam) == 'READASARRAYS') then
          this%readasarrays = .true.
        else if (this%block_tags(iparam) == 'TS6') then
          this%ts_active = .true.
        else if (this%block_tags(iparam) == 'EXPORT_ARRAY_ASCII') then
          this%export = .true.
        end if
      case default
      end select
    end do

    ! update input context allocations based on dfn set and input
    select case (this%mf6_input%block_dfns(iblk)%blockname)
    case ('OPTIONS')
      ! allocate naux and set to 0 if not allocated
      do iparam = 1, size(this%mf6_input%param_dfns)
        idt => this%mf6_input%param_dfns(iparam)
        if (idt%blockname == 'OPTIONS' .and. &
            idt%tagname == 'AUXILIARY') then
          if (this%iauxiliary == 0) then
            call mem_allocate(intptr, 'NAUX', this%mf6_input%mempath)
            intptr = 0
          end if
          exit
        end if
      end do
    case ('DIMENSIONS')
      ! set model shape if discretization dimensions have been read
      if (this%mf6_input%pkgtype(1:3) == 'DIS') then
        call set_model_shape(this%mf6_input%pkgtype, this%filename, &
                             this%mf6_input%component_mempath, &
                             this%mf6_input%mempath, this%mshape)
      end if
    case default
    end select
  end subroutine block_post_process

  !> @brief parse block
  !!
  !<
  recursive subroutine parse_block(this, iblk, recursive_call)
    use MemoryTypeModule, only: MemoryType
    use MemoryManagerModule, only: get_from_memorystore
    class(LoadMf6FileType) :: this
    integer(I4B), intent(in) :: iblk
    logical(LGP), intent(in) :: recursive_call !< true if recursive call
    logical(LGP) :: isblockfound
    logical(LGP) :: endOfBlock
    logical(LGP) :: supportOpenClose
    integer(I4B) :: ierr
    logical(LGP) :: found, required
    type(MemoryType), pointer :: mt

    ! disu vertices/cell2d blocks are contingent on NVERT dimension
    if (this%mf6_input%pkgtype == 'DISU6' .or. &
        this%mf6_input%pkgtype == 'DISV1D6' .or. &
        this%mf6_input%pkgtype == 'DISV2D6') then
      if (this%mf6_input%block_dfns(iblk)%blockname == 'VERTICES' .or. &
          this%mf6_input%block_dfns(iblk)%blockname == 'CELL2D') then
        call get_from_memorystore('NVERT', this%mf6_input%mempath, mt, found, &
                                  .false.)
        if (.not. found) return
        if (mt%intsclr == 0) return
      end if
    end if

    ! block open/close support
    supportOpenClose = (this%mf6_input%block_dfns(iblk)%blockname /= 'GRIDDATA')

    ! parser search for block
    required = this%mf6_input%block_dfns(iblk)%required .and. .not. recursive_call
    call this%parser%GetBlock(this%mf6_input%block_dfns(iblk)%blockname, &
                              isblockfound, ierr, &
                              supportOpenClose=supportOpenClose, &
                              blockRequired=required)
    ! process block
    if (isblockfound) then
      if (this%mf6_input%block_dfns(iblk)%aggregate) then
        ! process block recarray type, set of variable 1d/2d types
        call this%parse_structarray_block(iblk)
      else
        do
          ! process each line in block
          call this%parser%GetNextLine(endOfBlock)
          if (endOfBlock) exit
          ! process line as tag(s)
          call this%parse_tag(iblk, .false.)
        end do
      end if
    end if

    ! recurse if block is reloadable and was just read
    if (this%mf6_input%block_dfns(iblk)%block_variable) then
      if (isblockfound) then
        call this%parse_block(iblk, .true.)
      end if
    end if
  end subroutine parse_block

  subroutine parse_io_tag(this, iblk, pkgtype, which, tag)
    class(LoadMf6FileType) :: this
    integer(I4B), intent(in) :: iblk
    character(len=*), intent(in) :: pkgtype
    character(len=*), intent(in) :: which
    character(len=*), intent(in) :: tag
    type(InputParamDefinitionType), pointer :: idt !< input data type object describing this record
    ! matches, read and load file name
    idt => &
      get_param_definition_type(this%mf6_input%param_dfns, &
                                this%mf6_input%component_type, &
                                this%mf6_input%subcomponent_type, &
                                this%mf6_input%block_dfns(iblk)%blockname, &
                                tag, this%filename)
    ! load io tag
    call load_io_tag(this%parser, idt, this%mf6_input%mempath, which, this%iout)
  end subroutine parse_io_tag

  subroutine parse_keyword_tag(this, iblk, tag, idt)
    use DefinitionSelectModule, only: split_record_definition
    class(LoadMf6FileType) :: this
    integer(I4B), intent(in) :: iblk
    character(len=LINELENGTH), intent(in) :: tag
    type(InputParamDefinitionType), pointer, intent(in) :: idt
    character(len=40), dimension(:), allocatable :: words
    integer(I4B) :: nwords
    character(len=LINELENGTH) :: io_tag
    logical(LGP) :: found

    ! initialization
    found = .false.

    ! if in record tag check and load if input/output file
    if (idt%in_record) then
      ! get tokens in matching definition
      call split_record_definition(this%mf6_input%param_dfns, &
                                   this%mf6_input%component_type, &
                                   this%mf6_input%subcomponent_type, &
                                   tag, nwords, words)
      ! a filein/fileout record tag definition has 4 tokens
      if (nwords == 4) then
        ! verify third definition token is FILEIN/FILEOUT
        if (words(3) == 'FILEIN' .or. words(3) == 'FILEOUT') then
          ! read 3rd token
          call this%parser%GetStringCaps(io_tag)
          ! check if 3rd token matches definition
          if (io_tag == words(3)) then
            call this%parse_io_tag(iblk, words(2), words(3), words(4))
            found = .true.
          else
            errmsg = 'Expected "'//trim(words(3))//'" following keyword "'// &
                     trim(tag)//'" but instead found "'//trim(io_tag)//'"'
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end if
        end if
      end if

      ! deallocate words
      if (allocated(words)) deallocate (words)
    end if

    if (.not. found) then
      if (idt%shape == ':') then
        call load_string_type(this%parser, idt, this%mf6_input%mempath, this%iout)
      else
        ! load standard keyword tag
        call load_keyword_type(this%parser, idt, this%mf6_input%mempath, &
                               this%iout)
        ! check/set as dev option
        if (idt%tagname(1:4) == 'DEV_' .and. &
            this%mf6_input%block_dfns(iblk)%blockname == 'OPTIONS') then
          call this%parser%DevOpt()
        end if
      end if
    end if
  end subroutine parse_keyword_tag

  !> @brief load an individual input record into memory
  !!
  !! Load an individual input record into the memory
  !! manager.  Allow for recursive calls in the case that multiple
  !! tags are on a single line.
  !!
  !<
  recursive subroutine parse_tag(this, iblk, recursive_call)
    use ArrayHandlersModule, only: expandarray
    class(LoadMf6FileType) :: this
    integer(I4B), intent(in) :: iblk
    logical(LGP), intent(in) :: recursive_call !< true if recursive call
    character(len=LINELENGTH) :: tag
    type(InputParamDefinitionType), pointer :: idt !< input data type object describing this record
    integer(I4B) :: itag

    ! read tag name
    call this%parser%GetStringCaps(tag)
    if (recursive_call) then
      if (tag == '') then
        ! no data on line so return
        return
      end if
    end if

    ! find keyword in input definition
    idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                     this%mf6_input%component_type, &
                                     this%mf6_input%subcomponent_type, &
                                     this%mf6_input%block_dfns(iblk)%blockname, &
                                     tag, this%filename)

    ! allocate and load data type
    select case (idt%datatype)
    case ('KEYWORD')
      call this%parse_keyword_tag(iblk, tag, idt)
      if (idt%shape == ':') return
    case ('STRING')
      if (idt%shape == 'NAUX') then
        call load_auxvar_names(this%parser, idt, this%mf6_input%mempath, &
                               this%iout)
      else
        call load_string_type(this%parser, idt, this%mf6_input%mempath, this%iout)
      end if
    case ('INTEGER')
      call load_integer_type(this%parser, idt, this%mf6_input%mempath, this%iout)
    case ('INTEGER1D')
      call load_integer1d_type(this%parser, idt, this%mf6_input, this%mshape, &
                               this%export, this%nc_vars, this%filename, &
                               this%iout)
    case ('INTEGER2D')
      call load_integer2d_type(this%parser, idt, this%mf6_input, this%mshape, &
                               this%export, this%nc_vars, this%filename, &
                               this%iout)
    case ('INTEGER3D')
      call load_integer3d_type(this%parser, idt, this%mf6_input, this%mshape, &
                               this%export, this%nc_vars, this%filename, &
                               this%iout)
    case ('DOUBLE')
      call load_double_type(this%parser, idt, this%mf6_input%mempath, this%iout)
    case ('DOUBLE1D')
      call load_double1d_type(this%parser, idt, this%mf6_input, this%mshape, &
                              this%export, this%nc_vars, this%filename, this%iout)
    case ('DOUBLE2D')
      call load_double2d_type(this%parser, idt, this%mf6_input, this%mshape, &
                              this%export, this%nc_vars, this%filename, this%iout)
    case ('DOUBLE3D')
      call load_double3d_type(this%parser, idt, this%mf6_input, this%mshape, &
                              this%export, this%nc_vars, this%filename, this%iout)
    case default
      write (errmsg, '(a,a)') 'Failure reading data for tag: ', trim(tag)
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end select

    call expandarray(this%block_tags)
    itag = size(this%block_tags)
    this%block_tags(itag) = trim(idt%tagname)

    ! continue line if in same record
    if (idt%in_record) then
      ! recursively call parse tag again to read rest of line
      call this%parse_tag(iblk, .true.)
    end if
  end subroutine parse_tag

  function block_index_dfn(this, iblk) result(idt)
    class(LoadMf6FileType) :: this
    integer(I4B), intent(in) :: iblk
    type(InputParamDefinitionType) :: idt !< input data type object describing this record
    character(len=LENVARNAME) :: varname
    integer(I4B) :: ilen
    character(len=3) :: block_suffix = 'NUM'

    ! assign first column as the block number
    ilen = len_trim(this%mf6_input%block_dfns(iblk)%blockname)

    if (ilen > (LENVARNAME - len(block_suffix))) then
      varname = &
        this%mf6_input%block_dfns(iblk)% &
        blockname(1:(LENVARNAME - len(block_suffix)))//block_suffix
    else
      varname = trim(this%mf6_input%block_dfns(iblk)%blockname)//block_suffix
    end if

    idt%component_type = trim(this%mf6_input%component_type)
    idt%subcomponent_type = trim(this%mf6_input%subcomponent_type)
    idt%blockname = trim(this%mf6_input%block_dfns(iblk)%blockname)
    idt%tagname = varname
    idt%mf6varname = varname
    idt%datatype = 'INTEGER'
  end function block_index_dfn

  !> @brief parse a structured array record into memory manager
  !!
  !! A structarray is similar to a numpy recarray.  It it used to
  !! load a list of data in which each column in the list may be a
  !! different type.  Each column in the list is stored as a 1d
  !! vector.
  !!
  !<
  subroutine parse_structarray_block(this, iblk)
    use StructArrayModule, only: StructArrayType, constructStructArray
    use DynamicPackageParamsModule, only: DynamicPackageParamsType
    class(LoadMf6FileType) :: this
    integer(I4B), intent(in) :: iblk
    type(DynamicPackageParamsType) :: block_params
    type(InputParamDefinitionType), pointer :: idt !< input data type object describing this record
    type(InputParamDefinitionType), target :: blockvar_idt
    integer(I4B) :: blocknum
    integer(I4B), pointer :: nrow
    integer(I4B) :: nrows, nrowsread
    integer(I4B) :: ibinary, oc_inunit
    integer(I4B) :: icol, iparam
    integer(I4B) :: ncol

    ! initialize package params object
    call block_params%init(this%mf6_input, &
                           this%mf6_input%block_dfns(iblk)%blockname, &
                           this%readasarrays, this%iauxiliary, this%inamedbound)
    ! set input definition for this block
    idt => &
      get_aggregate_definition_type(this%mf6_input%aggregate_dfns, &
                                    this%mf6_input%component_type, &
                                    this%mf6_input%subcomponent_type, &
                                    this%mf6_input%block_dfns(iblk)%blockname)
    ! if block is reloadable read the block number
    if (this%mf6_input%block_dfns(iblk)%block_variable) then
      blocknum = this%parser%GetInteger()
    else
      blocknum = 0
    end if

    ! set ncol
    ncol = block_params%nparam
    ! add col if block is reloadable
    if (blocknum > 0) ncol = ncol + 1
    ! use shape to set the max num of rows
    if (idt%shape /= '') then
      call mem_setptr(nrow, idt%shape, this%mf6_input%mempath)
      nrows = nrow
    else
      nrows = -1
    end if

    ! create a structured array
    this%structarray => constructStructArray(this%mf6_input, ncol, nrows, &
                                             blocknum, this%mf6_input%mempath, &
                                             this%mf6_input%component_mempath)
    ! create structarray vectors for each column
    do icol = 1, ncol
      ! if block is reloadable, block number is first column
      if (blocknum > 0) then
        if (icol == 1) then
          blockvar_idt = this%block_index_dfn(iblk)
          idt => blockvar_idt
          call this%structarray%mem_create_vector(icol, idt)
          ! continue as this column managed by internally SA object
          cycle
        end if
        ! set indexes (where first column is blocknum)
        iparam = icol - 1
      else
        ! set indexes (no blocknum column)
        iparam = icol
      end if
      ! set pointer to input definition for this 1d vector
      idt => &
        get_param_definition_type(this%mf6_input%param_dfns, &
                                  this%mf6_input%component_type, &
                                  this%mf6_input%subcomponent_type, &
                                  this%mf6_input%block_dfns(iblk)%blockname, &
                                  block_params%params(iparam), this%filename)
      ! allocate variable in memory manager
      call this%structarray%mem_create_vector(icol, idt)
    end do

    ! read the block control record
    ibinary = read_control_record(this%parser, oc_inunit, this%iout)

    if (ibinary == 1) then
      ! read from binary
      nrowsread = this%structarray%read_from_binary(oc_inunit, this%iout)
      call this%parser%terminateblock()
      close (oc_inunit)
    else
      ! read from ascii
      nrowsread = this%structarray%read_from_parser(this%parser, this%ts_active, &
                                                    this%iout)
    end if

    ! clean up
    call block_params%destroy()
  end subroutine parse_structarray_block

  !> @brief load type keyword
  !<
  subroutine load_keyword_type(parser, idt, memoryPath, iout)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B), pointer :: intvar
    call mem_allocate(intvar, idt%mf6varname, memoryPath)
    intvar = 1
    call idm_log_var(intvar, idt%tagname, memoryPath, idt%datatype, iout)
  end subroutine load_keyword_type

  !> @brief load type string
  !<
  subroutine load_string_type(parser, idt, memoryPath, iout)
    use ConstantsModule, only: LENBIGLINE
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), intent(in) :: iout !< unit number for output
    character(len=LINELENGTH), pointer :: cstr
    character(len=LENBIGLINE), pointer :: bigcstr
    character(len=:), allocatable :: line
    integer(I4B) :: ilen
    select case (idt%shape)
    case (':')
      ilen = LENBIGLINE
      call parser%GetRemainingLine(line)
      call mem_allocate(bigcstr, ilen, idt%mf6varname, memoryPath)
      bigcstr = line
      call idm_log_var(bigcstr, idt%tagname, memoryPath, iout)
      deallocate (line)
    case ('LENBIGLINE')
      ilen = LENBIGLINE
      call mem_allocate(bigcstr, ilen, idt%mf6varname, memoryPath)
      call parser%GetString(bigcstr, (.not. idt%preserve_case))
      call idm_log_var(bigcstr, idt%tagname, memoryPath, iout)
    case default
      ilen = LINELENGTH
      call mem_allocate(cstr, ilen, idt%mf6varname, memoryPath)
      call parser%GetString(cstr, (.not. idt%preserve_case))
      call idm_log_var(cstr, idt%tagname, memoryPath, iout)
    end select
  end subroutine load_string_type

  !> @brief load io tag
  !<
  subroutine load_io_tag(parser, idt, memoryPath, which, iout)
    use MemoryManagerModule, only: mem_allocate, mem_reallocate, &
                                   mem_setptr, get_isize
    use CharacterStringModule, only: CharacterStringType
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    character(len=*), intent(in) :: which
    integer(I4B), intent(in) :: iout !< unit number for output
    character(len=LINELENGTH) :: cstr
    type(CharacterStringType), dimension(:), pointer, contiguous :: charstr1d
    integer(I4B) :: ilen, isize, idx
    ilen = LINELENGTH
    if (which == 'FILEIN') then
      call get_isize(idt%mf6varname, memoryPath, isize)
      if (isize < 0) then
        call mem_allocate(charstr1d, ilen, 1, idt%mf6varname, memoryPath)
        idx = 1
      else
        call mem_setptr(charstr1d, idt%mf6varname, memoryPath)
        call mem_reallocate(charstr1d, ilen, isize + 1, idt%mf6varname, &
                            memoryPath)
        idx = isize + 1
      end if
      call parser%GetString(cstr, (.not. idt%preserve_case))
      charstr1d(idx) = cstr
    else if (which == 'FILEOUT') then
      call load_string_type(parser, idt, memoryPath, iout)
    end if
  end subroutine load_io_tag

  !> @brief load aux variable names
  !!
  !<
  subroutine load_auxvar_names(parser, idt, memoryPath, iout)
    use ConstantsModule, only: LENAUXNAME, LINELENGTH, LENPACKAGENAME
    use InputOutputModule, only: urdaux
    use CharacterStringModule, only: CharacterStringType
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), intent(in) :: iout !< unit number for output
    character(len=:), allocatable :: line
    character(len=LENAUXNAME), dimension(:), allocatable :: caux
    integer(I4B) :: lloc
    integer(I4B) :: istart
    integer(I4B) :: istop
    integer(I4B) :: i
    character(len=LENPACKAGENAME) :: text = ''
    integer(I4B), pointer :: intvar
    type(CharacterStringType), dimension(:), &
      pointer, contiguous :: acharstr1d !< variable for allocation
    call mem_allocate(intvar, idt%shape, memoryPath)
    intvar = 0
    call parser%GetRemainingLine(line)
    lloc = 1
    call urdaux(intvar, parser%iuactive, iout, lloc, &
                istart, istop, caux, line, text)
    call mem_allocate(acharstr1d, LENAUXNAME, intvar, idt%mf6varname, memoryPath)
    do i = 1, intvar
      acharstr1d(i) = caux(i)
    end do
    deallocate (line)
    deallocate (caux)
  end subroutine load_auxvar_names

  !> @brief load type integer
  !<
  subroutine load_integer_type(parser, idt, memoryPath, iout)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B), pointer :: intvar
    call mem_allocate(intvar, idt%mf6varname, memoryPath)
    intvar = parser%GetInteger()
    call idm_log_var(intvar, idt%tagname, memoryPath, idt%datatype, iout)
  end subroutine load_integer_type

  !> @brief load type 1d integer
  !<
  subroutine load_integer1d_type(parser, idt, mf6_input, mshape, export, &
                                 nc_vars, input_fname, iout)
    use SourceCommonModule, only: get_shape_from_string, get_layered_shape
    use LoadNCInputModule, only: netcdf_read_array
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    type(ModflowInputType), intent(in) :: mf6_input !< description of input
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    logical(LGP), intent(in) :: export !< export to ascii layer files
    type(NCPackageVarsType), pointer, intent(in) :: nc_vars
    character(len=*), intent(in) :: input_fname !< ascii input file name
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B), dimension(:), pointer, contiguous :: int1d
    integer(I4B) :: nlay
    integer(I4B) :: nvals
    integer(I4B), dimension(:), allocatable :: array_shape
    integer(I4B), dimension(:), allocatable :: layer_shape
    character(len=LINELENGTH) :: keyword

    ! Check if it is a full grid sized array (NODES), otherwise use
    ! idt%shape to construct shape from variables in memoryPath
    if (idt%shape == 'NODES') then
      nvals = product(mshape)
    else
      call get_shape_from_string(idt%shape, array_shape, mf6_input%mempath)
      nvals = array_shape(1)
    end if

    ! allocate memory for the array
    call mem_allocate(int1d, nvals, idt%mf6varname, mf6_input%mempath)

    ! read keyword
    keyword = ''
    call parser%GetStringCaps(keyword)

    ! check for "NETCDF" and "LAYERED"
    if (keyword == 'NETCDF') then
      call netcdf_read_array(int1d, mshape, idt, mf6_input, nc_vars, &
                             input_fname, iout)
    else if (keyword == 'LAYERED' .and. idt%layered) then
      call get_layered_shape(mshape, nlay, layer_shape)
      call read_int1d_layered(parser, int1d, idt%mf6varname, nlay, layer_shape)
    else
      call read_int1d(parser, int1d, idt%mf6varname)
    end if

    ! log information on the loaded array to the list file
    call idm_log_var(int1d, idt%tagname, mf6_input%mempath, iout)

    ! create export file for griddata parameters if optioned
    if (export) then
      if (idt%blockname == 'GRIDDATA') then
        call idm_export(int1d, idt%tagname, mf6_input%mempath, idt%shape, iout)
      end if
    end if
  end subroutine load_integer1d_type

  !> @brief load type 2d integer
  !<
  subroutine load_integer2d_type(parser, idt, mf6_input, mshape, export, &
                                 nc_vars, input_fname, iout)
    use SourceCommonModule, only: get_shape_from_string, get_layered_shape
    use LoadNCInputModule, only: netcdf_read_array
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    type(ModflowInputType), intent(in) :: mf6_input !< description of input
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    logical(LGP), intent(in) :: export !< export to ascii layer files
    type(NCPackageVarsType), pointer, intent(in) :: nc_vars
    character(len=*), intent(in) :: input_fname !< ascii input file name
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B), dimension(:, :), pointer, contiguous :: int2d
    integer(I4B) :: nlay
    integer(I4B) :: nsize1, nsize2
    integer(I4B), dimension(:), allocatable :: array_shape
    integer(I4B), dimension(:), allocatable :: layer_shape
    character(len=LINELENGTH) :: keyword

    ! determine the array shape from the input data definition (idt%shape),
    ! which looks like "NCOL, NROW, NLAY"
    call get_shape_from_string(idt%shape, array_shape, mf6_input%mempath)
    nsize1 = array_shape(1)
    nsize2 = array_shape(2)

    ! create a new 3d memory managed variable
    call mem_allocate(int2d, nsize1, nsize2, idt%mf6varname, mf6_input%mempath)

    ! read keyword
    keyword = ''
    call parser%GetStringCaps(keyword)

    ! check for "NETCDF" and "LAYERED"
    if (keyword == 'NETCDF') then
      call netcdf_read_array(int2d, mshape, idt, mf6_input, nc_vars, &
                             input_fname, iout)
    else if (keyword == 'LAYERED' .and. idt%layered) then
      call get_layered_shape(mshape, nlay, layer_shape)
      call read_int2d_layered(parser, int2d, idt%mf6varname, nlay, layer_shape)
    else
      call read_int2d(parser, int2d, idt%mf6varname)
    end if

    ! log information on the loaded array to the list file
    call idm_log_var(int2d, idt%tagname, mf6_input%mempath, iout)

    ! create export file for griddata parameters if optioned
    if (export) then
      if (idt%blockname == 'GRIDDATA') then
        call idm_export(int2d, idt%tagname, mf6_input%mempath, idt%shape, iout)
      end if
    end if
  end subroutine load_integer2d_type

  !> @brief load type 3d integer
  !<
  subroutine load_integer3d_type(parser, idt, mf6_input, mshape, export, &
                                 nc_vars, input_fname, iout)
    use SourceCommonModule, only: get_shape_from_string, get_layered_shape
    use LoadNCInputModule, only: netcdf_read_array
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    type(ModflowInputType), intent(in) :: mf6_input !< description of input
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    logical(LGP), intent(in) :: export !< export to ascii layer files
    type(NCPackageVarsType), pointer, intent(in) :: nc_vars
    character(len=*), intent(in) :: input_fname !< ascii input file name
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B), dimension(:, :, :), pointer, contiguous :: int3d
    integer(I4B) :: nlay
    integer(I4B) :: nsize1, nsize2, nsize3
    integer(I4B), dimension(:), allocatable :: array_shape
    integer(I4B), dimension(:), allocatable :: layer_shape
    integer(I4B), dimension(:), pointer, contiguous :: int1d_ptr
    character(len=LINELENGTH) :: keyword

    ! determine the array shape from the input data definition (idt%shape),
    ! which looks like "NCOL, NROW, NLAY"
    call get_shape_from_string(idt%shape, array_shape, mf6_input%mempath)
    nsize1 = array_shape(1)
    nsize2 = array_shape(2)
    nsize3 = array_shape(3)

    ! create a new 3d memory managed variable
    call mem_allocate(int3d, nsize1, nsize2, nsize3, idt%mf6varname, &
                      mf6_input%mempath)

    ! read keyword
    keyword = ''
    call parser%GetStringCaps(keyword)

    ! check for "NETCDF" and "LAYERED"
    if (keyword == 'NETCDF') then
      call netcdf_read_array(int3d, mshape, idt, mf6_input, nc_vars, &
                             input_fname, iout)
    else if (keyword == 'LAYERED' .and. idt%layered) then
      call get_layered_shape(mshape, nlay, layer_shape)
      call read_int3d_layered(parser, int3d, idt%mf6varname, nlay, &
                              layer_shape)
    else
      int1d_ptr(1:nsize1 * nsize2 * nsize3) => int3d(:, :, :)
      call read_int1d(parser, int1d_ptr, idt%mf6varname)
    end if

    ! log information on the loaded array to the list file
    call idm_log_var(int3d, idt%tagname, mf6_input%mempath, iout)

    ! create export file for griddata parameters if optioned
    if (export) then
      if (idt%blockname == 'GRIDDATA') then
        call idm_export(int3d, idt%tagname, mf6_input%mempath, idt%shape, iout)
      end if
    end if
  end subroutine load_integer3d_type

  !> @brief load type double
  !<
  subroutine load_double_type(parser, idt, memoryPath, iout)
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B), intent(in) :: iout !< unit number for output
    real(DP), pointer :: dblvar
    call mem_allocate(dblvar, idt%mf6varname, memoryPath)
    dblvar = parser%GetDouble()
    call idm_log_var(dblvar, idt%tagname, memoryPath, iout)
  end subroutine load_double_type

  !> @brief load type 1d double
  !<
  subroutine load_double1d_type(parser, idt, mf6_input, mshape, export, &
                                nc_vars, input_fname, iout)
    use SourceCommonModule, only: get_shape_from_string, get_layered_shape
    use LoadNCInputModule, only: netcdf_read_array
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    type(ModflowInputType), intent(in) :: mf6_input !< description of input
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    logical(LGP), intent(in) :: export !< export to ascii layer files
    type(NCPackageVarsType), pointer, intent(in) :: nc_vars
    character(len=*), intent(in) :: input_fname !< ascii input file name
    integer(I4B), intent(in) :: iout !< unit number for output
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    integer(I4B) :: nlay
    integer(I4B) :: nvals
    integer(I4B), dimension(:), allocatable :: array_shape
    integer(I4B), dimension(:), allocatable :: layer_shape
    character(len=LINELENGTH) :: keyword

    ! Check if it is a full grid sized array (NODES)
    if (idt%shape == 'NODES') then
      nvals = product(mshape)
    else
      call get_shape_from_string(idt%shape, array_shape, mf6_input%mempath)
      nvals = array_shape(1)
    end if

    ! allocate memory for the array
    call mem_allocate(dbl1d, nvals, idt%mf6varname, mf6_input%mempath)

    ! read keyword
    keyword = ''
    call parser%GetStringCaps(keyword)

    ! check for "NETCDF" and "LAYERED"
    if (keyword == 'NETCDF') then
      call netcdf_read_array(dbl1d, mshape, idt, mf6_input, nc_vars, &
                             input_fname, iout)
    else if (keyword == 'LAYERED' .and. idt%layered) then
      call get_layered_shape(mshape, nlay, layer_shape)
      call read_dbl1d_layered(parser, dbl1d, idt%mf6varname, nlay, layer_shape)
    else
      call read_dbl1d(parser, dbl1d, idt%mf6varname)
    end if

    ! log information on the loaded array to the list file
    call idm_log_var(dbl1d, idt%tagname, mf6_input%mempath, iout)

    ! create export file for griddata parameters if optioned
    if (export) then
      if (idt%blockname == 'GRIDDATA') then
        call idm_export(dbl1d, idt%tagname, mf6_input%mempath, idt%shape, iout)
      end if
    end if
  end subroutine load_double1d_type

  !> @brief load type 2d double
  !<
  subroutine load_double2d_type(parser, idt, mf6_input, mshape, export, &
                                nc_vars, input_fname, iout)
    use SourceCommonModule, only: get_shape_from_string, get_layered_shape
    use LoadNCInputModule, only: netcdf_read_array
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    type(ModflowInputType), intent(in) :: mf6_input !< description of input
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    logical(LGP), intent(in) :: export !< export to ascii layer files
    type(NCPackageVarsType), pointer, intent(in) :: nc_vars
    character(len=*), intent(in) :: input_fname !< ascii input file name
    integer(I4B), intent(in) :: iout !< unit number for output
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B) :: nlay
    integer(I4B) :: nsize1, nsize2
    integer(I4B), dimension(:), allocatable :: array_shape
    integer(I4B), dimension(:), allocatable :: layer_shape
    character(len=LINELENGTH) :: keyword

    ! determine the array shape from the input data definition (idt%shape),
    ! which looks like "NCOL, NROW, NLAY"
    call get_shape_from_string(idt%shape, array_shape, mf6_input%mempath)
    nsize1 = array_shape(1)
    nsize2 = array_shape(2)

    ! create a new 3d memory managed variable
    call mem_allocate(dbl2d, nsize1, nsize2, idt%mf6varname, mf6_input%mempath)

    ! read keyword
    keyword = ''
    call parser%GetStringCaps(keyword)

    ! check for "NETCDF" and "LAYERED"
    if (keyword == 'NETCDF') then
      call netcdf_read_array(dbl2d, mshape, idt, mf6_input, nc_vars, &
                             input_fname, iout)
    else if (keyword == 'LAYERED' .and. idt%layered) then
      call get_layered_shape(mshape, nlay, layer_shape)
      call read_dbl2d_layered(parser, dbl2d, idt%mf6varname, nlay, layer_shape)
    else
      call read_dbl2d(parser, dbl2d, idt%mf6varname)
    end if

    ! log information on the loaded array to the list file
    call idm_log_var(dbl2d, idt%tagname, mf6_input%mempath, iout)

    ! create export file for griddata parameters if optioned
    if (export) then
      if (idt%blockname == 'GRIDDATA') then
        call idm_export(dbl2d, idt%tagname, mf6_input%mempath, idt%shape, iout)
      end if
    end if
  end subroutine load_double2d_type

  !> @brief load type 3d double
  !<
  subroutine load_double3d_type(parser, idt, mf6_input, mshape, export, &
                                nc_vars, input_fname, iout)
    use SourceCommonModule, only: get_shape_from_string, get_layered_shape
    use LoadNCInputModule, only: netcdf_read_array
    type(BlockParserType), intent(inout) :: parser !< block parser
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    type(ModflowInputType), intent(in) :: mf6_input !< description of input
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    logical(LGP), intent(in) :: export !< export to ascii layer files
    type(NCPackageVarsType), pointer, intent(in) :: nc_vars
    character(len=*), intent(in) :: input_fname !< ascii input file name
    integer(I4B), intent(in) :: iout !< unit number for output
    real(DP), dimension(:, :, :), pointer, contiguous :: dbl3d
    integer(I4B) :: nlay
    integer(I4B) :: nsize1, nsize2, nsize3
    integer(I4B), dimension(:), allocatable :: array_shape
    integer(I4B), dimension(:), allocatable :: layer_shape
    real(DP), dimension(:), pointer, contiguous :: dbl1d_ptr
    character(len=LINELENGTH) :: keyword

    ! determine the array shape from the input data definition (idt%shape),
    ! which looks like "NCOL, NROW, NLAY"
    call get_shape_from_string(idt%shape, array_shape, mf6_input%mempath)
    nsize1 = array_shape(1)
    nsize2 = array_shape(2)
    nsize3 = array_shape(3)

    ! create a new 3d memory managed variable
    call mem_allocate(dbl3d, nsize1, nsize2, nsize3, idt%mf6varname, &
                      mf6_input%mempath)

    ! read keyword
    keyword = ''
    call parser%GetStringCaps(keyword)

    ! check for "NETCDF" and "LAYERED"
    if (keyword == 'NETCDF') then
      call netcdf_read_array(dbl3d, mshape, idt, mf6_input, nc_vars, &
                             input_fname, iout)
    else if (keyword == 'LAYERED' .and. idt%layered) then
      call get_layered_shape(mshape, nlay, layer_shape)
      call read_dbl3d_layered(parser, dbl3d, idt%mf6varname, nlay, &
                              layer_shape)
    else
      dbl1d_ptr(1:nsize1 * nsize2 * nsize3) => dbl3d(:, :, :)
      call read_dbl1d(parser, dbl1d_ptr, idt%mf6varname)
    end if

    ! log information on the loaded array to the list file
    call idm_log_var(dbl3d, idt%tagname, mf6_input%mempath, iout)

    ! create export file for griddata parameters if optioned
    if (export) then
      if (idt%blockname == 'GRIDDATA') then
        call idm_export(dbl3d, idt%tagname, mf6_input%mempath, idt%shape, iout)
      end if
    end if
  end subroutine load_double3d_type

  function read_control_record(parser, oc_inunit, iout) result(ibinary)
    use SimModule, only: store_error_unit
    use InputOutputModule, only: urword
    use InputOutputModule, only: openfile
    use OpenSpecModule, only: form, access
    use ConstantsModule, only: LINELENGTH
    use BlockParserModule, only: BlockParserType
    type(BlockParserType), intent(inout) :: parser
    integer(I4B), intent(inout) :: oc_inunit
    integer(I4B), intent(in) :: iout
    integer(I4B) :: ibinary
    integer(I4B) :: lloc, istart, istop, idum, inunit, itmp, ierr
    integer(I4B) :: nunopn = 99
    character(len=:), allocatable :: line
    character(len=LINELENGTH) :: fname
    logical(LGP) :: exists
    real(DP) :: r
    character(len=*), parameter :: fmtocne = &
      &"('Specified OPEN/CLOSE file ',(A),' does not exist')"
    character(len=*), parameter :: fmtobf = &
      &"(1X,/1X,'OPENING BINARY FILE ON UNIT ',I0,':',/1X,A)"

    ! initialize oc_inunit and ibinary
    oc_inunit = 0
    ibinary = 0
    inunit = parser%getunit()

    ! Read to the first non-commented line
    lloc = 1
    call parser%line_reader%rdcom(inunit, iout, line, ierr)
    call urword(line, lloc, istart, istop, 1, idum, r, iout, inunit)

    if (line(istart:istop) == 'OPEN/CLOSE') then
      ! get filename
      call urword(line, lloc, istart, istop, 0, idum, r, &
                  iout, inunit)
      fname = line(istart:istop)
      ! check to see if file OPEN/CLOSE file exists
      inquire (file=fname, exist=exists)
      if (.not. exists) then
        write (errmsg, fmtocne) line(istart:istop)
        call store_error(errmsg)
        call store_error('Specified OPEN/CLOSE file does not exist')
        call store_error_unit(inunit)
      end if

      ! Check for (BINARY) keyword
      call urword(line, lloc, istart, istop, 1, idum, r, &
                  iout, inunit)

      if (line(istart:istop) == '(BINARY)') ibinary = 1
      ! Open the file depending on ibinary flag
      if (ibinary == 1) then
        oc_inunit = nunopn
        itmp = iout
        if (iout > 0) then
          itmp = 0
          write (iout, fmtobf) oc_inunit, trim(adjustl(fname))
        end if
        call openfile(oc_inunit, itmp, fname, 'OPEN/CLOSE', &
                      fmtarg_opt=form, accarg_opt=access)
      end if
    end if

    if (ibinary == 0) then
      call parser%line_reader%bkspc(parser%getunit())
    end if
  end function read_control_record

end module LoadMf6FileModule
