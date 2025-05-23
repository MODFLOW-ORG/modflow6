module DisuModule

  use ArrayReadersModule, only: ReadArray
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, LENVARNAME, &
                             DZERO, DONE, DHALF
  use ConnectionsModule, only: iac_to_ia
  use InputOutputModule, only: URWORD, ulasav, ulaprufw, ubdsv1, ubdsv06, &
                               getunit, openfile
  use SimModule, only: count_errors, store_error, store_error_unit, &
                       store_error_filename
  use SimVariablesModule, only: errmsg, idm_context
  use BaseDisModule, only: DisBaseType
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, &
                                 mem_reallocate, mem_setptr
  use MemoryManagerExtModule, only: mem_set_value, memorystore_remove
  use TdisModule, only: kstp, kper, pertim, totim, delt
  use DisvGeom, only: line_unit_vector

  implicit none

  private
  public :: DisuType
  public :: disu_cr
  public :: CastAsDisuType

  !> @brief Unstructured grid discretization
  type, extends(DisBaseType) :: DisuType
    integer(I4B), pointer :: njausr => null() ! user-specified nja size
    integer(I4B), pointer :: nvert => null() ! number of x,y vertices
    real(DP), pointer :: voffsettol => null() ! vertical offset tolerance
    real(DP), dimension(:, :), pointer, contiguous :: vertices => null() ! cell vertices stored as 2d array of x and y
    real(DP), dimension(:, :), pointer, contiguous :: cellxy => null() ! cell center stored as 2d array of x and y
    real(DP), dimension(:), pointer, contiguous :: top1d => null() ! (size:nodesuser) cell top elevation
    real(DP), dimension(:), pointer, contiguous :: bot1d => null() ! (size:nodesuser) cell bottom elevation
    real(DP), dimension(:), pointer, contiguous :: area1d => null() ! (size:nodesuser) cell area, in plan view
    integer(I4B), dimension(:), pointer, contiguous :: iainp => null() ! (size:nodesuser+1) user iac converted ia
    integer(I4B), dimension(:), pointer, contiguous :: jainp => null() ! (size:njausr) user-input ja array
    integer(I4B), dimension(:), pointer, contiguous :: ihcinp => null() ! (size:njausr) user-input ihc array
    real(DP), dimension(:), pointer, contiguous :: cl12inp => null() ! (size:njausr) user-input cl12 array
    real(DP), dimension(:), pointer, contiguous :: hwvainp => null() ! (size:njausr) user-input hwva array
    real(DP), dimension(:), pointer, contiguous :: angldegxinp => null() ! (size:njausr) user-input angldegx array
    integer(I4B), pointer :: iangledegx => null() ! =1 when angle information was present in input, 0 otherwise
    integer(I4B), dimension(:), pointer, contiguous :: iavert => null() ! cell vertex pointer ia array
    integer(I4B), dimension(:), pointer, contiguous :: javert => null() ! cell vertex pointer ja array
    logical(LGP) :: readFromFile ! True, when DIS is read from file (almost always)
  contains

    procedure :: dis_df => disu_df
    procedure :: disu_load
    procedure :: dis_da => disu_da
    procedure :: get_dis_type => get_dis_type
    procedure :: get_dis_enum => get_dis_enum
    procedure :: disu_ck
    procedure :: grid_finalize
    procedure :: get_nodenumber_idx1
    procedure :: nodeu_to_string
    procedure :: nodeu_to_array
    procedure :: nodeu_from_string
    procedure :: nodeu_from_cellid
    procedure :: connection_normal
    procedure :: connection_vector
    procedure :: supports_layers
    procedure :: get_ncpl
    procedure, public :: record_array
    procedure, public :: record_srcdst_list_header
    ! -- private
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: allocate_arrays_mem
    procedure :: source_options
    procedure :: source_dimensions
    procedure :: source_griddata
    procedure :: source_connectivity
    procedure :: source_vertices
    procedure :: source_cell2d
    procedure :: log_options
    procedure :: log_dimensions
    procedure :: log_griddata
    procedure :: log_connectivity
    procedure :: define_cellverts
    procedure :: write_grb
    !
    ! -- Read a node-sized model array (reduced or not)
    procedure :: read_int_array
    procedure :: read_dbl_array

  end type DisuType

  type DisuFoundType
    logical :: length_units = .false.
    logical :: nogrb = .false.
    logical :: xorigin = .false.
    logical :: yorigin = .false.
    logical :: angrot = .false.
    logical :: voffsettol = .false.
    logical :: nodes = .false.
    logical :: nja = .false.
    logical :: nvert = .false.
    logical :: top = .false.
    logical :: bot = .false.
    logical :: area = .false.
    logical :: idomain = .false.
    logical :: iac = .false.
    logical :: ja = .false.
    logical :: ihc = .false.
    logical :: cl12 = .false.
    logical :: hwva = .false.
    logical :: angldegx = .false.
    logical :: iv = .false.
    logical :: xv = .false.
    logical :: yv = .false.
    logical :: icell2d = .false.
    logical :: xc = .false.
    logical :: yc = .false.
    logical :: ncvert = .false.
    logical :: icvert = .false.
  end type DisuFoundType

contains

  !> @brief Create a new unstructured discretization object
  !<
  subroutine disu_cr(dis, name_model, input_mempath, inunit, iout)
    ! -- dummy
    class(DisBaseType), pointer :: dis
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    ! -- local
    type(DisuType), pointer :: disnew
    character(len=*), parameter :: fmtheader = &
      "(1X, /1X, 'DISU -- UNSTRUCTURED GRID DISCRETIZATION PACKAGE,', &
      &' VERSION 2 : 3/27/2014 - INPUT READ FROM MEMPATH: ', A, //)"
    !
    ! -- Create a new discretization object
    allocate (disnew)
    dis => disnew
    !
    ! -- Allocate scalars and assign data
    call dis%allocate_scalars(name_model, input_mempath)
    dis%inunit = inunit
    dis%iout = iout
    !
    ! -- If disu is enabled
    if (inunit > 0) then
      !
      ! -- Identify package
      if (iout > 0) then
        write (iout, fmtheader) dis%input_mempath
      end if
      !
      ! -- load disu
      call disnew%disu_load()
    end if
    !
  end subroutine disu_cr

  !> @brief Transfer IDM data into this discretization object
  !<
  subroutine disu_load(this)
    ! -- dummy
    class(DisuType) :: this
    !
    ! -- source input data
    call this%source_options()
    call this%source_dimensions()
    call this%source_griddata()
    call this%source_connectivity()
    !
    ! -- If NVERT specified and greater than 0, then source VERTICES and CELL2D
    if (this%nvert > 0) then
      call this%source_vertices()
      call this%source_cell2d()
    else
      ! -- connection direction information cannot be calculated
      this%icondir = 0
    end if
    !
    ! -- Make some final disu checks on the non-reduced user-provided
    !    input
    call this%disu_ck()
    !
  end subroutine disu_load

  !> @brief Define the discretization
  !<
  subroutine disu_df(this)
    ! -- dummy
    class(DisuType) :: this
    !
    call this%grid_finalize()
    !
  end subroutine disu_df

  !> @brief Finalize the grid
  !<
  subroutine grid_finalize(this)
    ! -- dummy
    class(DisuType) :: this
    ! -- locals
    integer(I4B) :: n
    integer(I4B) :: node
    integer(I4B) :: noder
    integer(I4B) :: nrsize
    ! -- formats
    character(len=*), parameter :: fmtdz = &
      "('CELL (',i0,',',i0,',',i0,') THICKNESS <= 0. ', &
      &'TOP, BOT: ',2(1pg24.15))"
    character(len=*), parameter :: fmtnr = &
      "(/1x, 'The specified IDOMAIN results in a reduced number of cells.',&
      &/1x, 'Number of user nodes: ',I0,&
      &/1X, 'Number of nodes in solution: ', I0, //)"
    !
    ! -- count active cells
    this%nodes = 0
    do n = 1, this%nodesuser
      if (this%idomain(n) > 0) this%nodes = this%nodes + 1
    end do
    !
    ! -- Check to make sure nodes is a valid number
    if (this%nodes == 0) then
      call store_error('Model does not have any active nodes. &
                       &Ensure IDOMAIN array has some values greater &
                       &than zero.')
      call store_error_filename(this%input_fname)
    end if
    !
    ! -- Write message if reduced grid
    if (this%nodes < this%nodesuser) then
      write (this%iout, fmtnr) this%nodesuser, this%nodes
    end if
    !
    ! -- Array size is now known, so allocate
    call this%allocate_arrays()
    !
    ! -- Fill the nodereduced array with the reduced nodenumber, or
    !    a negative number to indicate it is a pass-through cell, or
    !    a zero to indicate that the cell is excluded from the
    !    solution. (negative idomain not supported for disu)
    if (this%nodes < this%nodesuser) then
      noder = 1
      do node = 1, this%nodesuser
        if (this%idomain(node) > 0) then
          this%nodereduced(node) = noder
          noder = noder + 1
        elseif (this%idomain(node) < 0) then
          this%nodereduced(node) = -1
        else
          this%nodereduced(node) = 0
        end if
      end do
    end if
    !
    ! -- Fill nodeuser if a reduced grid
    if (this%nodes < this%nodesuser) then
      noder = 1
      do node = 1, this%nodesuser
        if (this%idomain(node) > 0) then
          this%nodeuser(noder) = node
          noder = noder + 1
        end if
      end do
    end if
    !
    ! -- Move top1d, bot1d, and area1d into top, bot, and area
    do node = 1, this%nodesuser
      noder = node
      if (this%nodes < this%nodesuser) noder = this%nodereduced(node)
      if (noder <= 0) cycle
      this%top(noder) = this%top1d(node)
      this%bot(noder) = this%bot1d(node)
      this%area(noder) = this%area1d(node)
    end do
    !
    ! -- fill cell center coordinates
    if (this%nvert > 0) then
      do node = 1, this%nodesuser
        noder = node
        if (this%nodes < this%nodesuser) noder = this%nodereduced(node)
        if (noder <= 0) cycle
        this%xc(noder) = this%cellxy(1, node)
        this%yc(noder) = this%cellxy(2, node)
      end do
    else
      call mem_reallocate(this%xc, 0, 'XC', this%memoryPath)
      call mem_reallocate(this%yc, 0, 'YC', this%memoryPath)
    end if
    !
    ! -- create and fill the connections object
    nrsize = 0
    if (this%nodes < this%nodesuser) nrsize = this%nodes
    allocate (this%con)
    call this%con%disuconnections(this%name_model, this%nodes, &
                                  this%nodesuser, nrsize, &
                                  this%nodereduced, this%nodeuser, &
                                  this%iainp, this%jainp, &
                                  this%ihcinp, this%cl12inp, &
                                  this%hwvainp, this%angldegxinp, &
                                  this%iangledegx)
    this%nja = this%con%nja
    this%njas = this%con%njas
    !
  end subroutine grid_finalize

  !> @brief Check discretization info
  !<
  subroutine disu_ck(this)
    ! -- dummy
    class(DisuType) :: this
    ! -- local
    integer(I4B) :: n, m
    integer(I4B) :: ipos
    integer(I4B) :: ihc
    real(DP) :: dz
    ! -- formats
    character(len=*), parameter :: fmtidm = &
      &"('Invalid idomain value ', i0, ' specified for node ', i0)"
    character(len=*), parameter :: fmtdz = &
      &"('Cell ', i0, ' with thickness <= 0. Top, bot: ', 2(1pg24.15))"
    character(len=*), parameter :: fmtarea = &
      &"('Cell ', i0, ' with area <= 0. Area: ', 1(1pg24.15))"
    character(len=*), parameter :: fmtjan = &
      &"('Cell ', i0, ' must have its first connection be itself.  Found: ', i0)"
    character(len=*), parameter :: fmtjam = &
      &"('Cell ', i0, ' has invalid connection in JA.  Found: ', i0)"
    character(len=*), parameter :: fmterrmsg = &
      "('Top elevation (', 1pg15.6, ') for cell ', i0, ' is above bottom &
      &elevation (', 1pg15.6, ') for cell ', i0, '. Based on node numbering &
      &rules cell ', i0, ' must be below cell ', i0, '.')"
    !
    ! -- Check connectivity
    do n = 1, this%nodesuser
      !
      ! -- Ensure first connection is to itself, and
      !    that ja(ia(n)) is positive
      ipos = this%iainp(n)
      m = this%jainp(ipos)
      if (m < 0) then
        m = abs(m)
        this%jainp(ipos) = m
      end if
      if (n /= m) then
        write (errmsg, fmtjan) n, m
        call store_error(errmsg)
      end if
      !
      ! -- Check for valid node numbers in connected cells
      do ipos = this%iainp(n) + 1, this%iainp(n + 1) - 1
        m = this%jainp(ipos)
        if (m < 0 .or. m > this%nodesuser) then
          ! -- make sure first connection is to itself
          write (errmsg, fmtjam) n, m
          call store_error(errmsg)
        end if
      end do
    end do
    !
    ! -- terminate if errors found
    if (count_errors() > 0) then
      if (this%inunit > 0) then
        call store_error_filename(this%input_fname)
      end if
    end if
    !
    ! -- Ensure idomain values are valid
    do n = 1, this%nodesuser
      if (this%idomain(n) > 1 .or. this%idomain(n) < 0) then
        write (errmsg, fmtidm) this%idomain(n), n
        call store_error(errmsg)
      end if
    end do
    !
    ! -- Check for zero and negative thickness and zero or negative areas
    !    for cells with idomain == 1
    do n = 1, this%nodesuser
      if (this%idomain(n) == 1) then
        dz = this%top1d(n) - this%bot1d(n)
        if (dz <= DZERO) then
          write (errmsg, fmt=fmtdz) n, this%top1d(n), this%bot1d(n)
          call store_error(errmsg)
        end if
        if (this%area1d(n) <= DZERO) then
          write (errmsg, fmt=fmtarea) n, this%area1d(n)
          call store_error(errmsg)
        end if
      end if
    end do
    !
    ! -- check to make sure voffsettol is >= 0
    if (this%voffsettol < DZERO) then
      write (errmsg, '(a, 1pg15.6)') &
        'Vertical offset tolerance must be greater than zero. Found ', &
        this%voffsettol
      call store_error(errmsg)
      if (this%inunit > 0) then
        call store_error_filename(this%input_fname)
      end if
    end if
    !
    ! -- For cell n, ensure that underlying cells have tops less than
    !    or equal to the bottom of cell n
    do n = 1, this%nodesuser
      do ipos = this%iainp(n) + 1, this%iainp(n + 1) - 1
        m = this%jainp(ipos)
        ihc = this%ihcinp(ipos)
        if (ihc == 0 .and. m > n) then
          dz = this%top1d(m) - this%bot1d(n)
          if (dz > this%voffsettol) then
            write (errmsg, fmterrmsg) this%top1d(m), m, this%bot1d(n), n, m, n
            call store_error(errmsg)
          end if
        end if
      end do
    end do
    !
    ! -- terminate if errors found
    if (count_errors() > 0) then
      if (this%inunit > 0) then
        call store_error_filename(this%input_fname)
      end if
    end if
    !
  end subroutine disu_ck

  !> @brief Deallocate variables
  !<
  subroutine disu_da(this)
    ! -- dummy
    class(DisuType) :: this
    !
    ! -- Deallocate idm memory
    call memorystore_remove(this%name_model, 'DISU', idm_context)
    call memorystore_remove(component=this%name_model, &
                            context=idm_context)
    !
    ! -- scalars
    call mem_deallocate(this%njausr)
    call mem_deallocate(this%nvert)
    call mem_deallocate(this%voffsettol)
    call mem_deallocate(this%iangledegx)
    !
    ! -- arrays
    if (this%readFromFile) then
      call mem_deallocate(this%top1d)
      call mem_deallocate(this%bot1d)
      call mem_deallocate(this%area1d)
      if (associated(this%iavert)) then
        call mem_deallocate(this%iavert)
        call mem_deallocate(this%javert)
      end if
      call mem_deallocate(this%vertices)
      call mem_deallocate(this%iainp)
      call mem_deallocate(this%jainp)
      call mem_deallocate(this%ihcinp)
      call mem_deallocate(this%cl12inp)
      call mem_deallocate(this%hwvainp)
      call mem_deallocate(this%angldegxinp)
    end if
    !
    call mem_deallocate(this%idomain)
    call mem_deallocate(this%cellxy)
    !
    call mem_deallocate(this%nodeuser)
    call mem_deallocate(this%nodereduced)
    !
    ! -- DisBaseType deallocate
    call this%DisBaseType%dis_da()
    !
  end subroutine disu_da

  !> @brief Convert a user nodenumber to a string (nodenumber)
  !<
  subroutine nodeu_to_string(this, nodeu, str)
    ! -- dummy
    class(DisuType) :: this
    integer(I4B), intent(in) :: nodeu
    character(len=*), intent(inout) :: str
    ! -- local
    character(len=10) :: nstr
    !
    write (nstr, '(i0)') nodeu
    str = '('//trim(adjustl(nstr))//')'
    !
  end subroutine nodeu_to_string

  !> @brief Convert a user nodenumber to an array (nodenumber)
  !<
  subroutine nodeu_to_array(this, nodeu, arr)
    class(DisuType) :: this
    integer(I4B), intent(in) :: nodeu
    integer(I4B), dimension(:), intent(inout) :: arr
    ! -- local
    integer(I4B) :: isize
    !
    ! -- check the size of arr
    isize = size(arr)
    if (isize /= this%ndim) then
      write (errmsg, '(a,i0,a,i0,a)') &
        'Program error: nodeu_to_array size of array (', isize, &
        ') is not equal to the discretization dimension (', this%ndim, ')'
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- fill array
    arr(1) = nodeu
    !
  end subroutine nodeu_to_array

  !> @brief Copy options from IDM into package
  !<
  subroutine source_options(this)
    ! -- dummy
    class(DisuType) :: this
    ! -- locals
    character(len=LENVARNAME), dimension(3) :: lenunits = &
      &[character(len=LENVARNAME) :: 'FEET', 'METERS', 'CENTIMETERS']
    type(DisuFoundType) :: found
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%lenuni, 'LENGTH_UNITS', this%input_mempath, &
                       lenunits, found%length_units)
    call mem_set_value(this%nogrb, 'NOGRB', this%input_mempath, found%nogrb)
    call mem_set_value(this%xorigin, 'XORIGIN', this%input_mempath, found%xorigin)
    call mem_set_value(this%yorigin, 'YORIGIN', this%input_mempath, found%yorigin)
    call mem_set_value(this%angrot, 'ANGROT', this%input_mempath, found%angrot)
    call mem_set_value(this%voffsettol, 'VOFFSETTOL', this%input_mempath, &
                       found%voffsettol)
    !
    ! -- log values to list file
    if (this%iout > 0) then
      call this%log_options(found)
    end if
    !
  end subroutine source_options

  !> @brief Write user options to list file
  !<
  subroutine log_options(this, found)
    ! -- dummy
    class(DisuType) :: this
    type(DisuFoundType), intent(in) :: found
    !
    write (this%iout, '(1x,a)') 'Setting Discretization Options'
    !
    if (found%length_units) then
      write (this%iout, '(4x,a,i0)') 'Model length unit [0=UND, 1=FEET, &
      &2=METERS, 3=CENTIMETERS] set as ', this%lenuni
    end if
    !
    if (found%nogrb) then
      write (this%iout, '(4x,a,i0)') 'Binary grid file [0=GRB, 1=NOGRB] &
        &set as ', this%nogrb
    end if
    !
    if (found%xorigin) then
      write (this%iout, '(4x,a,G0)') 'XORIGIN = ', this%xorigin
    end if
    !
    if (found%yorigin) then
      write (this%iout, '(4x,a,G0)') 'YORIGIN = ', this%yorigin
    end if
    !
    if (found%angrot) then
      write (this%iout, '(4x,a,G0)') 'ANGROT = ', this%angrot
    end if
    !
    if (found%voffsettol) then
      write (this%iout, '(4x,a,G0)') 'VERTICAL_OFFSET_TOLERANCE = ', &
        this%voffsettol
    end if
    !
    write (this%iout, '(1x,a,/)') 'End Setting Discretization Options'
    !
  end subroutine log_options

  !> @brief Copy dimensions from IDM into package
  !<
  subroutine source_dimensions(this)
    ! -- dummy
    class(DisuType) :: this
    ! -- locals
    integer(I4B) :: n
    type(DisuFoundType) :: found
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%nodesuser, 'NODES', this%input_mempath, found%nodes)
    call mem_set_value(this%njausr, 'NJA', this%input_mempath, found%nja)
    call mem_set_value(this%nvert, 'NVERT', this%input_mempath, found%nvert)
    !
    ! -- log simulation values
    if (this%iout > 0) then
      call this%log_dimensions(found)
    end if
    !
    ! -- verify dimensions were set
    if (this%nodesuser < 1) then
      call store_error( &
        'NODES was not specified or was specified incorrectly.')
    end if
    if (this%njausr < 1) then
      call store_error( &
        'NJA was not specified or was specified incorrectly.')
    end if
    !
    ! -- terminate if errors were detected
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
    !
    ! -- allocate vectors that are the size of nodesuser
    this%readFromFile = .true.
    call mem_allocate(this%top1d, this%nodesuser, 'TOP1D', this%memoryPath)
    call mem_allocate(this%bot1d, this%nodesuser, 'BOT1D', this%memoryPath)
    call mem_allocate(this%area1d, this%nodesuser, 'AREA1D', this%memoryPath)
    call mem_allocate(this%idomain, this%nodesuser, 'IDOMAIN', this%memoryPath)
    call mem_allocate(this%vertices, 2, this%nvert, 'VERTICES', this%memoryPath)
    call mem_allocate(this%iainp, this%nodesuser + 1, 'IAINP', this%memoryPath)
    call mem_allocate(this%jainp, this%njausr, 'JAINP', this%memoryPath)
    call mem_allocate(this%ihcinp, this%njausr, 'IHCINP', this%memoryPath)
    call mem_allocate(this%cl12inp, this%njausr, 'CL12INP', this%memoryPath)
    call mem_allocate(this%hwvainp, this%njausr, 'HWVAINP', this%memoryPath)
    call mem_allocate(this%angldegxinp, this%njausr, 'ANGLDEGXINP', &
                      this%memoryPath)
    if (this%nvert > 0) then
      call mem_allocate(this%cellxy, 2, this%nodesuser, 'CELLXY', this%memoryPath)
    else
      call mem_allocate(this%cellxy, 2, 0, 'CELLXY', this%memoryPath)
    end if
    !
    ! -- initialize all cells to be active (idomain = 1)
    do n = 1, this%nodesuser
      this%idomain(n) = 1
    end do
    !
  end subroutine source_dimensions

  !> @brief Write dimensions to list file
  !<
  subroutine log_dimensions(this, found)
    class(DisuType) :: this
    type(DisuFoundType), intent(in) :: found
    !
    write (this%iout, '(1x,a)') 'Setting Discretization Dimensions'
    !
    if (found%nodes) then
      write (this%iout, '(4x,a,i0)') 'NODES = ', this%nodesuser
    end if
    !
    if (found%nja) then
      write (this%iout, '(4x,a,i0)') 'NJA = ', this%njausr
    end if
    !
    if (found%nvert) then
      write (this%iout, '(4x,a,i0)') 'NVERT = ', this%nvert
    end if
    !
    write (this%iout, '(1x,a,/)') 'End Setting Discretization Dimensions'
    !
  end subroutine log_dimensions

  !> @brief Copy grid data from IDM into package
  !<
  subroutine source_griddata(this)
    ! -- dummy
    class(DisuType) :: this
    ! -- locals
    type(DisuFoundType) :: found
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%top1d, 'TOP', this%input_mempath, found%top)
    call mem_set_value(this%bot1d, 'BOT', this%input_mempath, found%bot)
    call mem_set_value(this%area1d, 'AREA', this%input_mempath, found%area)
    call mem_set_value(this%idomain, 'IDOMAIN', this%input_mempath, found%idomain)
    !
    ! -- log simulation values
    if (this%iout > 0) then
      call this%log_griddata(found)
    end if
    !
  end subroutine source_griddata

  !> @brief Write griddata found to list file
  !<
  subroutine log_griddata(this, found)
    ! -- dummy
    class(DisuType) :: this
    type(DisuFoundType), intent(in) :: found
    !
    write (this%iout, '(1x,a)') 'Setting Discretization Griddata'
    !
    if (found%top) then
      write (this%iout, '(4x,a)') 'TOP set from input file'
    end if
    !
    if (found%bot) then
      write (this%iout, '(4x,a)') 'BOT set from input file'
    end if
    !
    if (found%area) then
      write (this%iout, '(4x,a)') 'AREA set from input file'
    end if
    !
    if (found%idomain) then
      write (this%iout, '(4x,a)') 'IDOMAIN set from input file'
    end if
    !
    write (this%iout, '(1x,a,/)') 'End Setting Discretization Griddata'
    !
  end subroutine log_griddata

  !> @brief Copy grid connectivity info from IDM into package
  !<
  subroutine source_connectivity(this)
    ! -- dummy
    class(DisuType) :: this
    ! -- locals
    type(DisuFoundType) :: found
    integer(I4B), dimension(:), contiguous, pointer :: iac => null()
    ! -- formats
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%jainp, 'JA', this%input_mempath, found%ja)
    call mem_set_value(this%ihcinp, 'IHC', this%input_mempath, found%ihc)
    call mem_set_value(this%cl12inp, 'CL12', this%input_mempath, found%cl12)
    call mem_set_value(this%hwvainp, 'HWVA', this%input_mempath, found%hwva)
    call mem_set_value(this%angldegxinp, 'ANGLDEGX', this%input_mempath, &
                       found%angldegx)
    !
    ! -- set pointer to iac input array
    call mem_setptr(iac, 'IAC', this%input_mempath)
    !
    ! -- Convert iac to ia
    if (associated(iac)) call iac_to_ia(iac, this%iainp)
    !
    ! -- Set angldegx flag if found
    if (found%angldegx) this%iangledegx = 1
    !
    ! -- log simulation values
    if (this%iout > 0) then
      call this%log_connectivity(found, iac)
    end if
    !
  end subroutine source_connectivity

  !> @brief Write griddata found to list file
  !<
  subroutine log_connectivity(this, found, iac)
    class(DisuType) :: this
    type(DisuFoundType), intent(in) :: found
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: iac
    !
    write (this%iout, '(1x,a)') 'Setting Discretization Connectivity'
    !
    if (associated(iac)) then
      write (this%iout, '(4x,a)') 'IAC set from input file'
    end if
    !
    if (found%ja) then
      write (this%iout, '(4x,a)') 'JA set from input file'
    end if
    !
    if (found%ihc) then
      write (this%iout, '(4x,a)') 'IHC set from input file'
    end if
    !
    if (found%cl12) then
      write (this%iout, '(4x,a)') 'CL12 set from input file'
    end if
    !
    if (found%hwva) then
      write (this%iout, '(4x,a)') 'HWVA set from input file'
    end if
    !
    if (found%angldegx) then
      write (this%iout, '(4x,a)') 'ANGLDEGX set from input file'
    end if
    !
    write (this%iout, '(1x,a,/)') 'End Setting Discretization Connectivity'
    !
  end subroutine log_connectivity

  !> @brief Copy grid vertex data from IDM into package
  !<
  subroutine source_vertices(this)
    ! -- dummy
    class(DisuType) :: this
    ! -- local
    integer(I4B) :: i
    real(DP), dimension(:), contiguous, pointer :: vert_x => null()
    real(DP), dimension(:), contiguous, pointer :: vert_y => null()
    ! -- formats
    !
    ! -- set pointers to memory manager input arrays
    call mem_setptr(vert_x, 'XV', this%input_mempath)
    call mem_setptr(vert_y, 'YV', this%input_mempath)
    !
    ! -- set vertices 2d array
    if (associated(vert_x) .and. associated(vert_y)) then
      do i = 1, this%nvert
        this%vertices(1, i) = vert_x(i)
        this%vertices(2, i) = vert_y(i)
      end do
    else
      call store_error('Required Vertex arrays not found.')
    end if
    !
    ! -- log
    if (this%iout > 0) then
      write (this%iout, '(1x,a)') 'Discretization Vertex data loaded'
    end if
    !
  end subroutine source_vertices

  !> @brief Build data structures to hold cell vertex info
  !<
  subroutine define_cellverts(this, icell2d, ncvert, icvert)
    ! -- modules
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(DisuType) :: this
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: icell2d
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: ncvert
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: icvert
    ! -- locals
    type(sparsematrix) :: vert_spm
    integer(I4B) :: i, j, ierr
    integer(I4B) :: icv_idx, startvert, maxnnz = 5
    !
    ! -- initialize sparse matrix
    call vert_spm%init(this%nodesuser, this%nvert, maxnnz)
    !
    ! -- add sparse matrix connections from input memory paths
    icv_idx = 1
    do i = 1, this%nodesuser
      if (icell2d(i) /= i) call store_error('ICELL2D input sequence violation.')
      do j = 1, ncvert(i)
        call vert_spm%addconnection(i, icvert(icv_idx), 0)
        if (j == 1) then
          startvert = icvert(icv_idx)
        elseif (j == ncvert(i) .and. (icvert(icv_idx) /= startvert)) then
          call vert_spm%addconnection(i, startvert, 0)
        end if
        icv_idx = icv_idx + 1
      end do
    end do
    !
    ! -- allocate and fill iavert and javert
    call mem_allocate(this%iavert, this%nodesuser + 1, 'IAVERT', this%memoryPath)
    call mem_allocate(this%javert, vert_spm%nnz, 'JAVERT', this%memoryPath)
    call vert_spm%filliaja(this%iavert, this%javert, ierr)
    call vert_spm%destroy()
    !
  end subroutine define_cellverts

  !> @brief Copy cell2d data from IDM into package
  !<
  subroutine source_cell2d(this)
    ! -- dummy
    class(DisuType) :: this
    ! -- locals
    integer(I4B), dimension(:), contiguous, pointer :: icell2d => null()
    integer(I4B), dimension(:), contiguous, pointer :: ncvert => null()
    integer(I4B), dimension(:), contiguous, pointer :: icvert => null()
    real(DP), dimension(:), contiguous, pointer :: cell_x => null()
    real(DP), dimension(:), contiguous, pointer :: cell_y => null()
    integer(I4B) :: i
    !
    ! -- set pointers to input path ncvert and icvert
    call mem_setptr(icell2d, 'ICELL2D', this%input_mempath)
    call mem_setptr(ncvert, 'NCVERT', this%input_mempath)
    call mem_setptr(icvert, 'ICVERT', this%input_mempath)
    !
    ! --
    if (associated(icell2d) .and. associated(ncvert) &
        .and. associated(icvert)) then
      call this%define_cellverts(icell2d, ncvert, icvert)
    else
      call store_error('Required cell vertex arrays not found.')
    end if
    !
    ! -- set pointers to cell center arrays
    call mem_setptr(cell_x, 'XC', this%input_mempath)
    call mem_setptr(cell_y, 'YC', this%input_mempath)
    !
    ! -- set cell centers
    if (associated(cell_x) .and. associated(cell_y)) then
      do i = 1, this%nodesuser
        this%cellxy(1, i) = cell_x(i)
        this%cellxy(2, i) = cell_y(i)
      end do
    else
      call store_error('Required cell center arrays not found.')
    end if
    !
    ! -- log
    if (this%iout > 0) then
      write (this%iout, '(1x,a)') 'Discretization Cell2d data loaded'
    end if
    !
  end subroutine source_cell2d

  !> @brief Write a binary grid file
  !<
  subroutine write_grb(this, icelltype)
    ! -- modules
    use OpenSpecModule, only: access, form
    ! -- dummy
    class(DisuType) :: this
    integer(I4B), dimension(:), intent(in) :: icelltype
    ! -- local
    integer(I4B) :: i, iunit, ntxt
    integer(I4B), parameter :: lentxt = 100
    character(len=50) :: txthdr
    character(len=lentxt) :: txt
    character(len=LINELENGTH) :: fname
    ! -- formats
    character(len=*), parameter :: fmtgrdsave = &
      "(4X,'BINARY GRID INFORMATION WILL BE WRITTEN TO:', &
       &/,6X,'UNIT NUMBER: ', I0,/,6X, 'FILE NAME: ', A)"
    !
    ! -- Initialize
    ntxt = 11
    if (this%nvert > 0) ntxt = ntxt + 5
    !
    ! -- Open the file
    fname = trim(this%output_fname)
    iunit = getunit()
    write (this%iout, fmtgrdsave) iunit, trim(adjustl(fname))
    call openfile(iunit, this%iout, trim(adjustl(fname)), 'DATA(BINARY)', &
                  form, access, 'REPLACE')
    !
    ! -- write header information
    write (txthdr, '(a)') 'GRID DISU'
    txthdr(50:50) = new_line('a')
    write (iunit) txthdr
    write (txthdr, '(a)') 'VERSION 1'
    txthdr(50:50) = new_line('a')
    write (iunit) txthdr
    write (txthdr, '(a, i0)') 'NTXT ', ntxt
    txthdr(50:50) = new_line('a')
    write (iunit) txthdr
    write (txthdr, '(a, i0)') 'LENTXT ', lentxt
    txthdr(50:50) = new_line('a')
    write (iunit) txthdr
    !
    ! -- write variable definitions
    write (txt, '(3a, i0)') 'NODES ', 'INTEGER ', 'NDIM 0 # ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'NJA ', 'INTEGER ', 'NDIM 0 # ', this%con%nja
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, 1pg24.15)') 'XORIGIN ', 'DOUBLE ', 'NDIM 0 # ', this%xorigin
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, 1pg24.15)') 'YORIGIN ', 'DOUBLE ', 'NDIM 0 # ', this%yorigin
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, 1pg24.15)') 'ANGROT ', 'DOUBLE ', 'NDIM 0 # ', this%angrot
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'TOP ', 'DOUBLE ', 'NDIM 1 ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'BOT ', 'DOUBLE ', 'NDIM 1 ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'IA ', 'INTEGER ', 'NDIM 1 ', this%nodesuser + 1
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'JA ', 'INTEGER ', 'NDIM 1 ', this%con%nja
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'IDOMAIN ', 'INTEGER ', 'NDIM 1 ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'ICELLTYPE ', 'INTEGER ', 'NDIM 1 ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    !
    ! -- if vertices have been read then write additional header information
    if (this%nvert > 0) then
      write (txt, '(3a, i0)') 'VERTICES ', 'DOUBLE ', 'NDIM 2 2 ', this%nvert
      txt(lentxt:lentxt) = new_line('a')
      write (iunit) txt
      write (txt, '(3a, i0)') 'CELLX ', 'DOUBLE ', 'NDIM 1 ', this%nodesuser
      txt(lentxt:lentxt) = new_line('a')
      write (iunit) txt
      write (txt, '(3a, i0)') 'CELLY ', 'DOUBLE ', 'NDIM 1 ', this%nodesuser
      txt(lentxt:lentxt) = new_line('a')
      write (iunit) txt
      write (txt, '(3a, i0)') 'IAVERT ', 'INTEGER ', 'NDIM 1 ', this%nodesuser + 1
      txt(lentxt:lentxt) = new_line('a')
      write (iunit) txt
      write (txt, '(3a, i0)') 'JAVERT ', 'INTEGER ', 'NDIM 1 ', size(this%javert)
      txt(lentxt:lentxt) = new_line('a')
      write (iunit) txt
    end if
    !
    ! -- write data
    write (iunit) this%nodesuser ! nodes
    write (iunit) this%nja ! nja
    write (iunit) this%xorigin ! xorigin
    write (iunit) this%yorigin ! yorigin
    write (iunit) this%angrot ! angrot
    write (iunit) this%top1d ! top
    write (iunit) this%bot1d ! bot
    write (iunit) this%con%iausr ! ia
    write (iunit) this%con%jausr ! ja
    write (iunit) this%idomain ! idomain
    write (iunit) icelltype ! icelltype
    !
    ! -- if vertices have been read then write additional data
    if (this%nvert > 0) then
      write (iunit) this%vertices ! vertices
      write (iunit) (this%cellxy(1, i), i=1, this%nodesuser) ! cellx
      write (iunit) (this%cellxy(2, i), i=1, this%nodesuser) ! celly
      write (iunit) this%iavert ! iavert
      write (iunit) this%javert ! javert
    end if
    !
    ! -- Close the file
    close (iunit)
    !
  end subroutine write_grb

  !> @brief Get reduced node number from user node number
  !<
  function get_nodenumber_idx1(this, nodeu, icheck) result(nodenumber)
    class(DisuType), intent(in) :: this
    integer(I4B), intent(in) :: nodeu
    integer(I4B), intent(in) :: icheck
    integer(I4B) :: nodenumber
    !
    if (icheck /= 0) then
      if (nodeu < 1 .or. nodeu > this%nodes) then
        write (errmsg, '(a,i0,a,i0,a)') &
          'Node number (', nodeu, ') is less than 1 or greater than nodes (', &
          this%nodes, ').'
        call store_error(errmsg)
      end if
    end if
    !
    ! -- set node number to passed in nodenumber since there is a one to one
    !    mapping for an unstructured grid
    if (this%nodes == this%nodesuser) then
      nodenumber = nodeu
    else
      nodenumber = this%nodereduced(nodeu)
    end if
    !
  end function get_nodenumber_idx1

  !> @brief Get normal vector components between the cell and a given neighbor
  !!
  !! The normal points outward from the shared face between noden and nodem.
  !<
  subroutine connection_normal(this, noden, nodem, ihc, xcomp, ycomp, zcomp, &
                               ipos)
    ! -- dummy
    class(DisuType) :: this
    integer(I4B), intent(in) :: noden !< cell (reduced nn)
    integer(I4B), intent(in) :: nodem !< neighbor (reduced nn)
    integer(I4B), intent(in) :: ihc !< horizontal connection flag
    real(DP), intent(inout) :: xcomp
    real(DP), intent(inout) :: ycomp
    real(DP), intent(inout) :: zcomp
    integer(I4B), intent(in) :: ipos
    ! -- local
    real(DP) :: angle, dmult
    !
    ! -- Set vector components based on ihc
    if (ihc == 0) then
      !
      ! -- connection is vertical
      xcomp = DZERO
      ycomp = DZERO
      if (nodem < noden) then
        !
        ! -- nodem must be above noden, so upward connection
        zcomp = DONE
      else
        !
        ! -- nodem must be below noden, so downward connection
        zcomp = -DONE
      end if
    else
      ! -- find from anglex, since anglex is symmetric, need to flip vector
      !    for lower triangle (nodem < noden)
      angle = this%con%anglex(this%con%jas(ipos))
      dmult = DONE
      if (nodem < noden) dmult = -DONE
      xcomp = cos(angle) * dmult
      ycomp = sin(angle) * dmult
      zcomp = DZERO
    end if
    !
  end subroutine connection_normal

  !> @brief Get unit vector components between the cell and a given neighbor
  !!
  !! Saturation must be provided to compute cell center vertical coordinates.
  !! Also return the straight-line connection length.
  !<
  subroutine connection_vector(this, noden, nodem, nozee, satn, satm, ihc, &
                               xcomp, ycomp, zcomp, conlen)
    ! -- dummy
    class(DisuType) :: this
    integer(I4B), intent(in) :: noden
    integer(I4B), intent(in) :: nodem
    logical, intent(in) :: nozee
    real(DP), intent(in) :: satn
    real(DP), intent(in) :: satm
    integer(I4B), intent(in) :: ihc
    real(DP), intent(inout) :: xcomp
    real(DP), intent(inout) :: ycomp
    real(DP), intent(inout) :: zcomp
    real(DP), intent(inout) :: conlen
    ! -- local
    real(DP) :: xn, xm, yn, ym, zn, zm
    !
    ! -- Terminate with error if requesting unit vector components for problems
    !    without cell data
    if (size(this%cellxy, 2) < 1) then
      write (errmsg, '(a)') &
        'Cannot calculate unit vector components for DISU grid if VERTEX '// &
        'data are not specified'
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- get xy center coords
    xn = this%xc(noden)
    yn = this%yc(noden)
    xm = this%xc(nodem)
    ym = this%yc(nodem)
    !
    ! -- Set vector components based on ihc
    if (ihc == 0) then
      !
      ! -- vertical connection, calculate z as cell center elevation
      zn = this%bot(noden) + DHALF * (this%top(noden) - this%bot(noden))
      zm = this%bot(nodem) + DHALF * (this%top(nodem) - this%bot(nodem))
    else
      !
      ! -- horizontal connection, with possible z component due to cell offsets
      !    and/or water table conditions
      if (nozee) then
        zn = DZERO
        zm = DZERO
      else
        zn = this%bot(noden) + DHALF * satn * (this%top(noden) - this%bot(noden))
        zm = this%bot(nodem) + DHALF * satm * (this%top(nodem) - this%bot(nodem))
      end if
    end if
    !
    ! -- Use coords to find vector components and connection length
    call line_unit_vector(xn, yn, zn, xm, ym, zm, xcomp, ycomp, zcomp, &
                          conlen)
    !
  end subroutine connection_vector

  !> @brief Get the discretization type
  !<
  subroutine get_dis_type(this, dis_type)
    ! -- dummy
    class(DisuType), intent(in) :: this
    character(len=*), intent(out) :: dis_type
    !
    dis_type = "DISU"
    !
  end subroutine get_dis_type

  !> @brief Get the discretization type enumeration
  function get_dis_enum(this) result(dis_enum)
    use ConstantsModule, only: DISU
    class(DisuType), intent(in) :: this
    integer(I4B) :: dis_enum
    dis_enum = DISU
  end function get_dis_enum

  !> @brief Allocate and initialize scalar variables
  !<
  subroutine allocate_scalars(this, name_model, input_mempath)
    ! -- dummy
    class(DisuType) :: this
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    !
    ! -- Allocate parent scalars
    call this%DisBaseType%allocate_scalars(name_model, input_mempath)
    !
    ! -- Allocate variables for DISU
    call mem_allocate(this%njausr, 'NJAUSR', this%memoryPath)
    call mem_allocate(this%nvert, 'NVERT', this%memoryPath)
    call mem_allocate(this%voffsettol, 'VOFFSETTOL', this%memoryPath)
    call mem_allocate(this%iangledegx, 'IANGLEDEGX', this%memoryPath)
    !
    ! -- Set values
    this%ndim = 1
    this%njausr = 0
    this%nvert = 0
    this%voffsettol = DZERO
    this%iangledegx = 0
    this%readFromFile = .false.
    !
  end subroutine allocate_scalars

  !> @brief Allocate and initialize arrays
  !<
  subroutine allocate_arrays(this)
    ! -- dummy
    class(DisuType) :: this
    !
    ! -- Allocate arrays in DisBaseType (mshape, top, bot, area)
    call this%DisBaseType%allocate_arrays()
    !
    ! -- Allocate arrays in DISU
    if (this%nodes < this%nodesuser) then
      call mem_allocate(this%nodeuser, this%nodes, 'NODEUSER', this%memoryPath)
      call mem_allocate(this%nodereduced, this%nodesuser, 'NODEREDUCED', &
                        this%memoryPath)
    else
      call mem_allocate(this%nodeuser, 1, 'NODEUSER', this%memoryPath)
      call mem_allocate(this%nodereduced, 1, 'NODEREDUCED', this%memoryPath)
    end if
    !
    ! -- Initialize
    this%mshape(1) = this%nodesuser
    !
  end subroutine allocate_arrays

  !> @brief Allocate arrays in memory manager
  !<
  subroutine allocate_arrays_mem(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(DisuType) :: this
    !
    call mem_allocate(this%idomain, this%nodes, 'IDOMAIN', this%memoryPath)
    call mem_allocate(this%vertices, 2, this%nvert, 'VERTICES', this%memoryPath)
    call mem_allocate(this%cellxy, 2, this%nodes, 'CELLXY', this%memoryPath)
    !
  end subroutine allocate_arrays_mem

  !> @brief Convert a string to a user nodenumber
  !!
  !! Parse and return user nodenumber.
  !! If flag_string is present and true, the first token may be
  !! non-numeric (e.g. boundary name). In this case, return -2.
  !<
  function nodeu_from_string(this, lloc, istart, istop, in, iout, line, &
                             flag_string, allow_zero) result(nodeu)
    ! -- dummy
    class(DisuType) :: this
    integer(I4B), intent(inout) :: lloc
    integer(I4B), intent(inout) :: istart
    integer(I4B), intent(inout) :: istop
    integer(I4B), intent(in) :: in
    integer(I4B), intent(in) :: iout
    character(len=*), intent(inout) :: line
    logical, optional, intent(in) :: flag_string
    logical, optional, intent(in) :: allow_zero
    integer(I4B) :: nodeu
    ! -- local
    integer(I4B) :: lloclocal, ndum, istat, n
    real(DP) :: r
    !
    if (present(flag_string)) then
      if (flag_string) then
        ! Check to see if first token in line can be read as an integer.
        lloclocal = lloc
        call urword(line, lloclocal, istart, istop, 1, ndum, r, iout, in)
        read (line(istart:istop), *, iostat=istat) n
        if (istat /= 0) then
          ! First token in line is not an integer; return flag to this effect.
          nodeu = -2
          return
        end if
      end if
    end if
    !
    call urword(line, lloc, istart, istop, 2, nodeu, r, iout, in)
    !
    if (nodeu == 0) then
      if (present(allow_zero)) then
        if (allow_zero) then
          return
        end if
      end if
    end if
    !
    if (nodeu < 1 .or. nodeu > this%nodesuser) then
      write (errmsg, '(a,i0,a)') &
        "Node number in list (", nodeu, ") is outside of the grid. "// &
        "Cell number cannot be determined in line '"// &
        trim(adjustl(line))//"'."
      call store_error(errmsg)
      call store_error_unit(in)
    end if
    !
  end function nodeu_from_string

  !> @brief Convert a cellid string to a user nodenumber
  !!
  !! If flag_string is present and true, the first token may be
  !! non-numeric (e.g. boundary name). In this case, return -2.
  !!
  !! If allow_zero is present and true, and all indices are zero, the
  !! result can be zero. If allow_zero is false, a zero in any index is an error.
  !<
  function nodeu_from_cellid(this, cellid, inunit, iout, flag_string, &
                             allow_zero) result(nodeu)
    ! -- return
    integer(I4B) :: nodeu
    ! -- dummy
    class(DisuType) :: this
    character(len=*), intent(inout) :: cellid
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    logical, optional, intent(in) :: flag_string
    logical, optional, intent(in) :: allow_zero
    ! -- local
    integer(I4B) :: lloclocal, istart, istop, ndum, n
    integer(I4B) :: istat
    real(DP) :: r
    !
    if (present(flag_string)) then
      if (flag_string) then
        ! Check to see if first token in cellid can be read as an integer.
        lloclocal = 1
        call urword(cellid, lloclocal, istart, istop, 1, ndum, r, iout, inunit)
        read (cellid(istart:istop), *, iostat=istat) n
        if (istat /= 0) then
          ! First token in cellid is not an integer; return flag to this effect.
          nodeu = -2
          return
        end if
      end if
    end if
    !
    lloclocal = 1
    call urword(cellid, lloclocal, istart, istop, 2, nodeu, r, iout, inunit)
    !
    if (nodeu == 0) then
      if (present(allow_zero)) then
        if (allow_zero) then
          return
        end if
      end if
    end if
    !
    if (nodeu < 1 .or. nodeu > this%nodesuser) then
      write (errmsg, '(a,i0,a)') &
        "Cell number cannot be determined for cellid ("// &
        trim(adjustl(cellid))//") and results in a user "// &
        "node number (", nodeu, ") that is outside of the grid."
      call store_error(errmsg)
      call store_error_unit(inunit)
    end if
    !
  end function nodeu_from_cellid

  !> @brief Indicates whether the grid discretization supports layers
  !<
  logical function supports_layers(this)
    ! -- dummy
    class(DisuType) :: this
    !
    supports_layers = .false.
    !
  end function supports_layers

  !> @brief Get number of cells per layer (total nodes since DISU isn't layered)
  !<
  function get_ncpl(this)
    ! -- return
    integer(I4B) :: get_ncpl
    ! -- dummy
    class(DisuType) :: this
    !
    get_ncpl = this%nodesuser
    !
  end function get_ncpl

  !> @brief Read an integer array
  !<
  subroutine read_int_array(this, line, lloc, istart, istop, iout, in, &
                            iarray, aname)
    ! -- dummy
    class(DisuType), intent(inout) :: this
    character(len=*), intent(inout) :: line
    integer(I4B), intent(inout) :: lloc
    integer(I4B), intent(inout) :: istart
    integer(I4B), intent(inout) :: istop
    integer(I4B), intent(in) :: in
    integer(I4B), intent(in) :: iout
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: iarray
    character(len=*), intent(in) :: aname
    ! -- local
    integer(I4B) :: nval
    integer(I4B), dimension(:), pointer, contiguous :: itemp
    !
    ! -- Point the temporary pointer array, which is passed to the reading
    !    subroutine.  The temporary array will point to ibuff if it is a
    !    reduced structured system, or to iarray if it is an unstructured
    !    model.
    if (this%nodes < this%nodesuser) then
      nval = this%nodesuser
      itemp => this%ibuff
    else
      nval = this%nodes
      itemp => iarray
    end if
    !
    ! -- Read the array
    ! -- Read unstructured input
    call ReadArray(in, itemp, aname, this%ndim, nval, iout, 0)
    !
    ! -- If reduced model, then need to copy from itemp(=>ibuff) to iarray
    if (this%nodes < this%nodesuser) then
      call this%fill_grid_array(itemp, iarray)
    end if
    !
  end subroutine read_int_array

  !> @brief Read a double precision array
  !<
  subroutine read_dbl_array(this, line, lloc, istart, istop, iout, in, &
                            darray, aname)
    ! -- dummy
    class(DisuType), intent(inout) :: this
    character(len=*), intent(inout) :: line
    integer(I4B), intent(inout) :: lloc
    integer(I4B), intent(inout) :: istart
    integer(I4B), intent(inout) :: istop
    integer(I4B), intent(in) :: in
    integer(I4B), intent(in) :: iout
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: darray
    character(len=*), intent(in) :: aname
    ! -- local
    integer(I4B) :: nval
    real(DP), dimension(:), pointer, contiguous :: dtemp
    !
    ! -- Point the temporary pointer array, which is passed to the reading
    !    subroutine.  The temporary array will point to dbuff if it is a
    !    reduced structured system, or to darray if it is an unstructured
    !    model.
    if (this%nodes < this%nodesuser) then
      nval = this%nodesuser
      dtemp => this%dbuff
    else
      nval = this%nodes
      dtemp => darray
    end if
    !
    ! -- Read the array
    call ReadArray(in, dtemp, aname, this%ndim, nval, iout, 0)
    !
    ! -- If reduced model, then need to copy from dtemp(=>dbuff) to darray
    if (this%nodes < this%nodesuser) then
      call this%fill_grid_array(dtemp, darray)
    end if
    !
  end subroutine read_dbl_array

  !> @brief Record a double precision array
  !!
  !! The array is written to a formatted or unformatted external file
  !! depending on the arguments.
  !<
  subroutine record_array(this, darray, iout, iprint, idataun, aname, &
                          cdatafmp, nvaluesp, nwidthp, editdesc, dinact)
    ! -- dummy
    class(DisuType), intent(inout) :: this
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: darray !< double precision array to record
    integer(I4B), intent(in) :: iout !< ascii output unit number
    integer(I4B), intent(in) :: iprint !< whether to print the array
    integer(I4B), intent(in) :: idataun !< binary output unit number
    character(len=*), intent(in) :: aname !< text descriptor
    character(len=*), intent(in) :: cdatafmp ! write format
    integer(I4B), intent(in) :: nvaluesp !< values per line
    integer(I4B), intent(in) :: nwidthp !< number width
    character(len=*), intent(in) :: editdesc !< format type (I, G, F, S, E)
    real(DP), intent(in) :: dinact !< double precision value for cells excluded from model domain
    ! -- local
    integer(I4B) :: k, ifirst
    integer(I4B) :: nlay
    integer(I4B) :: nrow
    integer(I4B) :: ncol
    integer(I4B) :: nval
    integer(I4B) :: nodeu, noder
    integer(I4B) :: istart, istop
    real(DP), dimension(:), pointer, contiguous :: dtemp
    ! -- formats
    character(len=*), parameter :: fmthsv = &
      "(1X,/1X,a,' WILL BE SAVED ON UNIT ',I4, &
       &' AT END OF TIME STEP',I5,', STRESS PERIOD ',I4)"
    !
    ! -- set variables
    nlay = 1
    nrow = 1
    ncol = this%mshape(1)
    !
    ! -- If this is a reduced model, then copy the values from darray into
    !    dtemp.
    if (this%nodes < this%nodesuser) then
      nval = this%nodes
      dtemp => this%dbuff
      do nodeu = 1, this%nodesuser
        noder = this%get_nodenumber(nodeu, 0)
        if (noder <= 0) then
          dtemp(nodeu) = dinact
          cycle
        end if
        dtemp(nodeu) = darray(noder)
      end do
    else
      nval = this%nodes
      dtemp => darray
    end if
    !
    ! -- Print to iout if iprint /= 0
    if (iprint /= 0) then
      istart = 1
      do k = 1, nlay
        istop = istart + nrow * ncol - 1
        call ulaprufw(ncol, nrow, kstp, kper, k, iout, dtemp(istart:istop), &
                      aname, cdatafmp, nvaluesp, nwidthp, editdesc)
        istart = istop + 1
      end do
    end if
    !
    ! -- Save array to an external file.
    if (idataun > 0) then
      ! -- write to binary file by layer
      ifirst = 1
      istart = 1
      do k = 1, nlay
        istop = istart + nrow * ncol - 1
        if (ifirst == 1) write (iout, fmthsv) &
          trim(adjustl(aname)), idataun, &
          kstp, kper
        ifirst = 0
        call ulasav(dtemp(istart:istop), aname, kstp, kper, &
                    pertim, totim, ncol, nrow, k, idataun)
        istart = istop + 1
      end do
    elseif (idataun < 0) then
      !
      ! -- write entire array as one record
      call ubdsv1(kstp, kper, aname, -idataun, dtemp, ncol, nrow, nlay, &
                  iout, delt, pertim, totim)
    end if
    !
  end subroutine record_array

  !> @brief Record list header for imeth=6
  !<
  subroutine record_srcdst_list_header(this, text, textmodel, textpackage, &
                                       dstmodel, dstpackage, naux, auxtxt, &
                                       ibdchn, nlist, iout)
    ! -- dummy
    class(DisuType) :: this
    character(len=16), intent(in) :: text
    character(len=16), intent(in) :: textmodel
    character(len=16), intent(in) :: textpackage
    character(len=16), intent(in) :: dstmodel
    character(len=16), intent(in) :: dstpackage
    integer(I4B), intent(in) :: naux
    character(len=16), dimension(:), intent(in) :: auxtxt
    integer(I4B), intent(in) :: ibdchn
    integer(I4B), intent(in) :: nlist
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B) :: nlay, nrow, ncol
    !
    nlay = 1
    nrow = 1
    ncol = this%mshape(1)
    !
    ! -- Use ubdsv06 to write list header
    call ubdsv06(kstp, kper, text, textmodel, textpackage, dstmodel, dstpackage, &
                 ibdchn, naux, auxtxt, ncol, nrow, nlay, &
                 nlist, iout, delt, pertim, totim)
    !
  end subroutine record_srcdst_list_header

  !> @brief Cast base to DISU
  !<
  function CastAsDisuType(dis) result(disu)
    ! -- dummy
    class(*), pointer :: dis !< base pointer to DISU object
    ! -- return
    class(DisuType), pointer :: disu !< the resulting DISU pointer
    !
    disu => null()
    select type (dis)
    class is (DisuType)
      disu => dis
    end select
    !
  end function CastAsDisuType

end module DisuModule
