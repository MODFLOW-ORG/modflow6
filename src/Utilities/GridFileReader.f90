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

  type :: GridFileReaderType
    private
    integer(I4B), public :: inunit !< file unit
    ! header
    character(len=10), public :: grid_type !< DIS, DISV, DISU, etc
    integer(I4B), public :: version !< binary grid file format version
    integer(I4B) :: ntxt !< number of variables
    integer(I4B) :: lentxt !< header line length per variable
    ! index
    type(HashTableType), pointer :: dim !< map variable name to number of dims
    type(HashTableType), pointer :: pos !< map variable name to position in file
    type(HashTableType), pointer :: typ !< map variable name to type (1=int, 2=dbl)
    type(HashTableType), pointer :: shp_idx !< map variable name to index in shp
    integer(I4B), allocatable :: shp(:) !< flat array of variable shapes
    character(len=10), allocatable, public :: keys(:) !< variable names
  contains
    procedure, public :: initialize
    procedure, public :: finalize
    procedure, public :: read_int
    procedure, public :: read_dbl
    procedure, public :: read_int_1d
    procedure, public :: read_dbl_1d
    procedure, public :: read_charstr
    procedure, public :: read_grid_shape
    procedure, private :: read_header
    procedure, private :: read_header_meta
    procedure, private :: read_header_body
  end type GridFileReaderType

contains

  !> @Brief Initialize the grid file reader.
  subroutine initialize(this, iu)
    class(GridFileReaderType) :: this
    integer(I4B), intent(in) :: iu

    this%inunit = iu
    call hash_table_cr(this%dim)
    call hash_table_cr(this%pos)
    call hash_table_cr(this%typ)
    call hash_table_cr(this%shp_idx)
    allocate (this%shp(0))
    call this%read_header()

  end subroutine initialize

  !> @brief Finalize the grid file reader.
  subroutine finalize(this)
    class(GridFileReaderType) :: this

    close (this%inunit)
    call hash_table_da(this%dim)
    call hash_table_da(this%pos)
    call hash_table_da(this%typ)
    call hash_table_da(this%shp_idx)
    deallocate (this%shp)

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
    integer(I4B) :: i, lloc, istart, istop, ival, pos
    integer(I4B) :: nvars, ndim, dim, ishp
    integer(I4B), allocatable :: shp(:)

    allocate (this%keys(this%ntxt))
    allocate (character(len=this%lentxt*this%ntxt) :: body)
    allocate (character(len=this%lentxt) :: line)

    nvars = 0
    read (this%inunit) body
    inquire (this%inunit, pos=pos)
    do i = 1, this%lentxt * this%ntxt, this%lentxt
      line = body(i:i + this%lentxt - 1)
      lloc = 1

      ! key
      lloc = 1
      call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
      key = line(istart:istop)
      nvars = nvars + 1
      this%keys(nvars) = key

      ! type
      call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
      dtype = line(istart:istop)
      if (dtype == "INTEGER") then
        call this%typ%add(key, 1)
      else if (dtype == "DOUBLE") then
        call this%typ%add(key, 2)
      else if (dtype == "CHARACTER") then
        call this%typ%add(key, 3)
      end if

      ! dims
      call urword(line, lloc, istart, istop, 0, ival, rval, 0, 0)
      call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
      ndim = ival
      call this%dim%add(key, ndim)

      ! shape
      if (allocated(shp)) deallocate (shp)
      allocate (shp(ndim))
      if (ndim > 0) then
        do dim = 1, ndim
          call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
          shp(dim) = ival
        end do
        ishp = size(this%shp)
        call ExpandArray(this%shp, increment=ndim)
        this%shp(ishp + 1:ishp + ndim) = shp
        call this%shp_idx%add(key, ishp + 1)
      end if

      ! position
      call this%pos%add(key, pos)
      if (ndim == 0) then
        if (dtype == "INTEGER") then
          pos = pos + 4
        else if (dtype == "DOUBLE") then
          pos = pos + 8
        end if
      else
        if (dtype == "INTEGER") then
          pos = pos + (product(shp) * 4)
        else if (dtype == "DOUBLE") then
          pos = pos + (product(shp) * 8)
        else if (dtype == "CHARACTER") then
          pos = pos + (product(shp) * 8)
        end if
      end if
    end do

    rewind (this%inunit)

  end subroutine read_header_body

  !> @brief Read an integer scalar from a grid file.
  function read_int(this, key) result(v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer(I4B) :: v
    ! local
    integer(I4B) :: ndim, pos, typ
    character(len=:), allocatable :: msg

    msg = 'Variable '//trim(key)//' is not an integer scalar'
    ndim = this%dim%get(key)
    if (ndim /= 0) then
      write (errmsg, '(a)') msg
      call store_error(errmsg, terminate=.TRUE.)
    end if
    typ = this%typ%get(key)
    if (typ /= 1) then
      write (errmsg, '(a)') msg
      call store_error(errmsg, terminate=.TRUE.)
    end if
    pos = this%pos%get(key)
    read (this%inunit, pos=pos) v
    rewind (this%inunit)

  end function read_int

  !> @brief Read a double precision scalar from a grid file.
  function read_dbl(this, key) result(v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(DP) :: v
    ! local
    integer(I4B) :: ndim, pos, typ
    character(len=:), allocatable :: msg

    msg = 'Variable '//trim(key)//' is not a double precision scalar'
    ndim = this%dim%get(key)
    if (ndim /= 0) then
      write (errmsg, '(a)') msg
      call store_error(errmsg, terminate=.TRUE.)
    end if
    typ = this%typ%get(key)
    if (typ /= 2) then
      write (errmsg, '(a)') msg
      call store_error(errmsg, terminate=.TRUE.)
    end if
    pos = this%pos%get(key)
    read (this%inunit, pos=pos) v
    rewind (this%inunit)

  end function read_dbl

  !> @brief Read a 1D integer array from a grid file.
  function read_int_1d(this, key) result(v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer(I4B), allocatable :: v(:)
    ! local
    integer(I4B) :: idx, ndim, nvals, pos, typ
    character(len=:), allocatable :: msg

    msg = 'Variable '//trim(key)//' is not a 1D integer array'
    ndim = this%dim%get(key)
    if (ndim /= 1) then
      write (errmsg, '(a)') msg
      call store_error(errmsg, terminate=.TRUE.)
    end if
    typ = this%typ%get(key)
    if (typ /= 1) then
      write (errmsg, '(a)') msg
      call store_error(errmsg, terminate=.TRUE.)
    end if
    idx = this%shp_idx%get(key)
    pos = this%pos%get(key)
    nvals = this%shp(idx)
    allocate (v(nvals))
    read (this%inunit, pos=pos) v
    rewind (this%inunit)

  end function read_int_1d

  !> @brief Read a 1D double array from a grid file.
  function read_dbl_1d(this, key) result(v)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(DP), allocatable :: v(:)
    ! local
    integer(I4B) :: idx, ndim, nvals, pos, typ
    character(len=:), allocatable :: msg

    msg = 'Variable '//trim(key)//' is not a 1D double array'
    ndim = this%dim%get(key)
    if (ndim /= 1) then
      write (errmsg, '(a)') msg
      call store_error(errmsg, terminate=.TRUE.)
    end if
    typ = this%typ%get(key)
    if (typ /= 2) then
      write (errmsg, '(a)') msg
      call store_error(errmsg, terminate=.TRUE.)
    end if
    idx = this%shp_idx%get(key)
    pos = this%pos%get(key)
    nvals = this%shp(idx)
    allocate (v(nvals))
    read (this%inunit, pos=pos) v
    rewind (this%inunit)

  end function read_dbl_1d

  !> @brief Read a character string from a grid file.
  function read_charstr(this, key) result(charstr)
    class(GridFileReaderType), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=:), allocatable :: charstr
    ! local
    integer(I4B) :: idx, ndim, nvals, pos, typ
    character(len=:), allocatable :: msg

    msg = 'Variable '//trim(key)//' is not a character array'
    ndim = this%dim%get(key)
    if (ndim /= 1) then
      write (errmsg, '(a)') msg
      call store_error(errmsg, terminate=.TRUE.)
    end if
    typ = this%typ%get(key)
    if (typ /= 3) then
      write (errmsg, '(a)') msg
      call store_error(errmsg, terminate=.TRUE.)
    end if
    idx = this%shp_idx%get(key)
    pos = this%pos%get(key)
    nvals = this%shp(idx)
    allocate (character(nvals) :: charstr)
    read (this%inunit, pos=pos) charstr
    rewind (this%inunit)
  end function read_charstr

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

end module GridFileReaderModule
