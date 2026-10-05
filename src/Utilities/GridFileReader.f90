module GridFileReaderModule

  use KindModule
  use SimModule, only: store_error, store_error_unit
  use SimVariablesModule, only: errmsg
  use ConstantsModule, only: LINELENGTH
  use InputOutputModule, only: urword, openfile
  use HashTableModule, only: HashTableType, hash_table_cr, hash_table_da
  use ArrayHandlersModule, only: ExpandArray

  implicit none

  public :: GridFileReaderType

  integer(I4B), parameter :: TYP_INT = 1 !< integer variable type
  integer(I4B), parameter :: TYP_DBL = 2 !< double precision variable type
  integer(I4B), parameter :: TYP_CHR = 3 !< character variable type

  type :: GridFileReaderType
    private
    integer(I4B), public :: inunit !< file unit
    ! header
    character(len=10), public :: grid_type !< DIS, DISV, DISU, etc
    integer(I4B), public :: version !< binary grid file format version
    integer(I4B) :: ntxt !< number of variables
    integer(I4B) :: lentxt !< header line length per variable
    ! index
    type(HashTableType), pointer :: idx !< map variable name to variable index
    character(len=10), allocatable, public :: keys(:) !< variable names
    integer(I4B), allocatable :: ndims(:) !< variable number of dims
    integer(I4B), allocatable :: typs(:) !< variable type (TYP_INT, TYP_DBL, TYP_CHR)
    integer(I4B), allocatable :: shp_start(:) !< variable shape start in shp
    integer(I8B), allocatable :: pos(:) !< variable position in file
    integer(I4B), allocatable :: shp(:) !< flat array of variable shapes
  contains
    procedure, public :: initialize
    procedure, public :: finalize
    procedure, public :: has_variable
    ! header-reading subroutines
    procedure, private :: read_header
    procedure, private :: read_header_meta
    procedure, private :: read_header_body
    procedure, private :: lookup
    ! scalar read functions
    procedure, public :: read_int
    procedure, public :: read_dbl
    procedure, public :: read_grid_shape
    ! array read functions (allocate and return)
    procedure, public :: read_int_1d
    procedure, public :: read_dbl_1d
    procedure, public :: read_charstr
    ! array read subroutines (populate preallocated)
    procedure, public :: read_int_1d_into
    procedure, public :: read_dbl_1d_into
    procedure, public :: read_charstr_into
  end type GridFileReaderType

contains

  !> @Brief Initialize the grid file reader.
  subroutine initialize(this, iu)
    class(GridFileReaderType) :: this
    integer(I4B), intent(in) :: iu

    this%inunit = iu
    call hash_table_cr(this%idx)
    allocate (this%shp(0))
    call this%read_header()

  end subroutine initialize

  !> @brief Finalize the grid file reader.
  subroutine finalize(this)
    class(GridFileReaderType) :: this

    close (this%inunit)
    call hash_table_da(this%idx)
    if (allocated(this%keys)) deallocate (this%keys)
    if (allocated(this%ndims)) deallocate (this%ndims)
    if (allocated(this%typs)) deallocate (this%typs)
    if (allocated(this%shp_start)) deallocate (this%shp_start)
    if (allocated(this%pos)) deallocate (this%pos)
    if (allocated(this%shp)) deallocate (this%shp)

  end subroutine finalize

  !> @brief Read the file's self-describing header. Internal use only.
  subroutine read_header(this)
    class(GridFileReaderType) :: this
    call this%read_header_meta()
    call this%read_header_body()
  end subroutine read_header

  !> @brief Read self-describing metadata (first four lines). Internal use only.
  subroutine read_header_meta(this)
    ! dummy
    class(GridFileReaderType) :: this
    ! local
    character(len=50) :: line
    integer(I4B) :: lloc, istart, istop
    integer(I4B) :: ival
    real(DP) :: rval

    ! grid type
    read (this%inunit) line
    lloc = 1
    call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
    if (line(istart:istop) /= 'GRID') then
      call store_error('Binary grid file must begin with "GRID". '//&
                       &'Found: '//line(istart:istop))
      call store_error_unit(this%inunit)
    end if
    call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
    this%grid_type = line(istart:istop)

    ! version
    read (this%inunit) line
    lloc = 1
    call urword(line, lloc, istart, istop, 0, ival, rval, 0, 0)
    call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
    this%version = ival

    ! ntxt
    read (this%inunit) line
    lloc = 1
    call urword(line, lloc, istart, istop, 0, ival, rval, 0, 0)
    call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
    this%ntxt = ival

    ! lentxt
    read (this%inunit) line
    lloc = 1
    call urword(line, lloc, istart, istop, 0, ival, rval, 0, 0)
    call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
    this%lentxt = ival

  end subroutine read_header_meta

  !> @brief Read the header body section (text following first
  !< four "meta" lines) and build an index. Internal use only.
  subroutine read_header_body(this)
    ! dummy
    class(GridFileReaderType) :: this
    ! local
    character(len=:), allocatable :: body
    character(len=:), allocatable :: line
    character(len=10) :: key, dtype
    real(DP) :: rval
    integer(I4B) :: i, lloc, istart, istop, ival
    integer(I4B) :: ivar, ndim, dim, ishp, nbytes
    integer(I8B) :: pos

    allocate (this%keys(this%ntxt))
    allocate (this%ndims(this%ntxt))
    allocate (this%typs(this%ntxt))
    allocate (this%shp_start(this%ntxt))
    allocate (this%pos(this%ntxt))
    allocate (character(len=this%lentxt*this%ntxt) :: body)
    allocate (character(len=this%lentxt) :: line)

    read (this%inunit) body
    inquire (this%inunit, pos=pos)
    do ivar = 1, this%ntxt
      i = (ivar - 1) * this%lentxt + 1
      line = body(i:i + this%lentxt - 1)

      ! key
      lloc = 1
      call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
      key = line(istart:istop)
      this%keys(ivar) = key
      call this%idx%add(key, ivar)

      ! type
      call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
      dtype = line(istart:istop)
      select case (dtype)
      case ("INTEGER")
        this%typs(ivar) = TYP_INT
        nbytes = 4
      case ("DOUBLE")
        this%typs(ivar) = TYP_DBL
        nbytes = 8
      case ("CHARACTER")
        this%typs(ivar) = TYP_CHR
        nbytes = 1
      case default
        this%typs(ivar) = 0
        nbytes = 0
      end select

      ! dims
      call urword(line, lloc, istart, istop, 0, ival, rval, 0, 0)
      call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
      ndim = ival
      this%ndims(ivar) = ndim

      ! shape
      this%shp_start(ivar) = 0
      if (ndim > 0) then
        ishp = size(this%shp)
        call ExpandArray(this%shp, increment=ndim)
        do dim = 1, ndim
          call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
          this%shp(ishp + dim) = ival
        end do
        this%shp_start(ivar) = ishp + 1
      end if

      ! position
      this%pos(ivar) = pos
      if (ndim == 0) then
        pos = pos + nbytes
      else
        ishp = this%shp_start(ivar)
        pos = pos + product(int(this%shp(ishp:ishp + ndim - 1), I8B)) * nbytes
      end if
    end do

    rewind (this%inunit)

  end subroutine read_header_body

  !> @brief Look up a variable and check its rank and type.
  !!
  !! Returns the variable's index. Terminates with an error if the
  !! variable does not exist or does not have the expected rank and type.
  !! Internal use only.
  !<
  function lookup(this, key, ndim, typ, desc) result(ivar)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer(I4B), intent(in) :: ndim !< expected number of dims
    integer(I4B), intent(in) :: typ !< expected type
    character(len=*), intent(in) :: desc !< expected kind, for error messages
    integer(I4B) :: ivar

    ivar = this%idx%get(key)
    if (ivar == 0) then
      write (errmsg, '(a)') 'Variable '//trim(key)//' not found'
      call store_error(errmsg, terminate=.TRUE.)
    end if
    if (this%ndims(ivar) /= ndim .or. this%typs(ivar) /= typ) then
      write (errmsg, '(a)') 'Variable '//trim(key)//' is not '//desc
      call store_error(errmsg, terminate=.TRUE.)
    end if
  end function lookup

  !> @brief Read an integer scalar from a grid file.
  function read_int(this, key) result(v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer(I4B) :: v
    ! local
    integer(I4B) :: ivar

    ivar = this%lookup(key, 0, TYP_INT, 'an integer scalar')
    read (this%inunit, pos=this%pos(ivar)) v
    rewind (this%inunit)
  end function read_int

  !> @brief Read a double precision scalar from a grid file.
  function read_dbl(this, key) result(v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(DP) :: v
    ! local
    integer(I4B) :: ivar

    ivar = this%lookup(key, 0, TYP_DBL, 'a double precision scalar')
    read (this%inunit, pos=this%pos(ivar)) v
    rewind (this%inunit)
  end function read_dbl

  !> @brief Read a 1D integer array from a grid file.
  !!
  !! Allocates and returns a new array containing the data.
  !<
  function read_int_1d(this, key) result(v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer(I4B), allocatable :: v(:)
    ! local
    integer(I4B) :: ivar, nvals

    ivar = this%lookup(key, 1, TYP_INT, 'a 1D integer array')
    nvals = this%shp(this%shp_start(ivar))
    allocate (v(nvals))
    read (this%inunit, pos=this%pos(ivar)) v
    rewind (this%inunit)
  end function read_int_1d

  !> @brief Read a 1D integer array into a preallocated array.
  !!
  !! Populates a preallocated array. Array must already be allocated to the
  !! correct size. This version is compatible with both allocatable arrays and
  !! memory-manager-allocated pointer targets.
  !<
  subroutine read_int_1d_into(this, key, v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer(I4B), dimension(:), intent(inout) :: v
    ! local
    integer(I4B) :: ivar, nvals

    ivar = this%lookup(key, 1, TYP_INT, 'a 1D integer array')
    nvals = this%shp(this%shp_start(ivar))
    if (size(v) /= nvals) then
      write (errmsg, '(a,i0,a,i0)') &
        'Array size mismatch for '//trim(key)//': expected ', &
        nvals, ', got ', size(v)
      call store_error(errmsg, terminate=.TRUE.)
    end if
    read (this%inunit, pos=this%pos(ivar)) v
    rewind (this%inunit)
  end subroutine read_int_1d_into

  !> @brief Read a 1D double array from a grid file.
  !!
  !! Allocates and returns a new array containing the data.
  !<
  function read_dbl_1d(this, key) result(v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(DP), allocatable :: v(:)
    ! local
    integer(I4B) :: ivar, nvals

    ivar = this%lookup(key, 1, TYP_DBL, 'a 1D double array')
    nvals = this%shp(this%shp_start(ivar))
    allocate (v(nvals))
    read (this%inunit, pos=this%pos(ivar)) v
    rewind (this%inunit)
  end function read_dbl_1d

  !> @brief Read a 1D double array into a preallocated array.
  !!
  !! Populates a preallocated array. Array must already be allocated to the
  !! correct size. This version is compatible with both allocatable arrays and
  !! memory-manager-allocated pointer targets.
  !<
  subroutine read_dbl_1d_into(this, key, v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(DP), dimension(:), intent(inout) :: v
    ! local
    integer(I4B) :: ivar, nvals

    ivar = this%lookup(key, 1, TYP_DBL, 'a 1D double array')
    nvals = this%shp(this%shp_start(ivar))
    if (size(v) /= nvals) then
      write (errmsg, '(a,i0,a,i0)') &
        'Array size mismatch for '//trim(key)//': expected ', &
        nvals, ', got ', size(v)
      call store_error(errmsg, terminate=.TRUE.)
    end if
    read (this%inunit, pos=this%pos(ivar)) v
    rewind (this%inunit)
  end subroutine read_dbl_1d_into

  !> @brief Read a character string from a grid file.
  !!
  !! Allocates and returns a new character string containing the data.
  !<
  function read_charstr(this, key) result(charstr)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=:), allocatable :: charstr
    ! local
    integer(I4B) :: ivar, nvals

    ivar = this%lookup(key, 1, TYP_CHR, 'a character array')
    nvals = this%shp(this%shp_start(ivar))
    allocate (character(nvals) :: charstr)
    read (this%inunit, pos=this%pos(ivar)) charstr
    rewind (this%inunit)
  end function read_charstr

  !> @brief Read a character string into a preallocated string.
  !!
  !! Populates a preallocated character string. If the string is not allocated
  !! or is the wrong length, it will be (re)allocated to the correct length.
  !<
  subroutine read_charstr_into(this, key, charstr)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=:), allocatable, intent(inout) :: charstr
    ! local
    integer(I4B) :: ivar, nvals

    ivar = this%lookup(key, 1, TYP_CHR, 'a character array')
    nvals = this%shp(this%shp_start(ivar))
    if (allocated(charstr)) then
      if (len(charstr) /= nvals) deallocate (charstr)
    end if
    if (.not. allocated(charstr)) allocate (character(nvals) :: charstr)
    read (this%inunit, pos=this%pos(ivar)) charstr
    rewind (this%inunit)
  end subroutine read_charstr_into

  !> @brief Read the grid shape from a grid file.
  function read_grid_shape(this) result(v)
    ! dummy
    class(GridFileReaderType) :: this
    integer(I4B), allocatable :: v(:)

    select case (this%grid_type)
    case ("DIS")
      allocate (v(3))
      v(1) = this%read_int("NLAY")
      v(2) = this%read_int("NROW")
      v(3) = this%read_int("NCOL")
    case ("DISV")
      allocate (v(2))
      v(1) = this%read_int("NLAY")
      v(2) = this%read_int("NCPL")
    case ("DISU")
      allocate (v(1))
      v(1) = this%read_int("NODES")
    case ("DIS2D")
      allocate (v(2))
      v(1) = this%read_int("NROW")
      v(2) = this%read_int("NCOL")
    case ("DISV2D")
      allocate (v(1))
      v(1) = this%read_int("NODES")
    case ("DISV1D")
      allocate (v(1))
      v(1) = this%read_int("NCELLS")
    end select

  end function read_grid_shape

  !> @brief Check whether the grid file contains a variable.
  function has_variable(this, key) result(has)
    class(GridFileReaderType) :: this
    character(len=*), intent(in) :: key
    logical(LGP) :: has

    has = this%idx%get(key) /= 0
  end function has_variable

end module GridFileReaderModule
