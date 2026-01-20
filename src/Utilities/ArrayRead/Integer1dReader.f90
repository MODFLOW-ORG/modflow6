module Integer1dReaderModule

  use KindModule, only: DP, I4B, LGP
  use BlockParserModule, only: BlockParserType
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_unit
  use ArrayReadersModule, only: read_binary_header, &
                                check_binary_filesize, &
                                BINARY_INT_BYTES, &
                                BINARY_HEADER_BYTES
  use ArrayReaderBaseModule, only: ArrayReaderBaseType

  implicit none
  private
  public :: read_int1d, read_int1d_layered

  type, extends(ArrayReaderBaseType) :: Integer1dReaderType

    integer(I4B) :: constant_array_value = 0
    integer(I4B) :: factor = 1
    integer(I4B), dimension(:), contiguous, pointer :: int1d => null()

  contains

    procedure :: reset_reader
    procedure :: set_constant ! must be overridden
    procedure :: fill_constant ! must be overridden
    procedure :: read_ascii ! must be overridden
    procedure :: read_binary ! must be overridden
    procedure :: set_factor ! must be overridden
    procedure :: apply_factor ! must be overridden

  end type Integer1dReaderType

contains

  subroutine read_int1d(parser, int1d, aname)
    ! -- dummy
    type(BlockParserType), intent(in), target :: parser
    integer(I4B), dimension(:), contiguous, target :: int1d
    character(len=*), intent(in) :: aname
    ! -- local
    type(Integer1dReaderType) :: this

    this%parser => parser
    this%int1d => int1d
    this%array_name = aname

    call this%read_array()

  end subroutine read_int1d

  subroutine read_int1d_layered(parser, int1d, aname, nlay, layer_shape)
    use Integer2dReaderModule, only: read_int2d
    ! -- dummy
    type(BlockParserType), intent(in), target :: parser
    integer(I4B), dimension(:), contiguous, target :: int1d
    character(len=*), intent(in) :: aname
    integer(I4B), intent(in) :: nlay
    integer(I4B), dimension(:), intent(in) :: layer_shape
    ! -- local
    integer(I4B) :: k
    integer(I4B) :: ncpl, nrow, ncol
    integer(I4B) :: index_start, index_stop
    integer(I4B), dimension(:, :), contiguous, pointer :: int2d_ptr

    ncpl = product(layer_shape)
    index_start = 1
    do k = 1, nlay
      index_stop = index_start + ncpl - 1
      if (size(layer_shape) == 2) then
        ncol = layer_shape(1)
        nrow = layer_shape(2)
        int2d_ptr(1:ncol, 1:nrow) => int1d(index_start:index_stop)
        call read_int2d(parser, int2d_ptr, aname)
      else
        call read_int1d(parser, int1d(index_start:index_stop), aname)
      end if
      index_start = index_stop + 1
    end do

  end subroutine read_int1d_layered

  subroutine reset_reader(this)
    class(Integer1dReaderType) :: this
    call this%ArrayReaderBaseType%reset_reader()
    this%constant_array_value = 0
    this%factor = 1
  end subroutine reset_reader

  subroutine set_constant(this)
    class(Integer1dReaderType) :: this
    this%constant_array_value = this%parser%GetInteger()
  end subroutine set_constant

  subroutine fill_constant(this)
    class(Integer1dReaderType) :: this
    integer(I4B) :: i
    do i = 1, size(this%int1d)
      this%int1d(i) = this%constant_array_value
    end do
  end subroutine fill_constant

  subroutine read_ascii(this)
    class(Integer1dReaderType) :: this
    integer(I4B) :: i
    integer(I4B) :: nvals
    integer(I4B) :: istat
    nvals = size(this%int1d)
    read (this%input_unit, *, iostat=istat, iomsg=errmsg) &
      (this%int1d(i), i=1, size(this%int1d))
    if (istat /= 0) then
      errmsg = 'Error reading data for array '//trim(this%array_name)// &
               '.  '//trim(errmsg)
      call store_error(errmsg)
      call store_error_unit(this%input_unit)
    end if
  end subroutine read_ascii

  subroutine read_binary(this)
    class(Integer1dReaderType) :: this
    integer(I4B) :: i
    integer(I4B) :: nvals
    integer(I4B) :: istat
    integer(I4B) :: expected_size
    expected_size = BINARY_HEADER_BYTES + (size(this%int1d) * BINARY_INT_BYTES)
    call read_binary_header(this%input_unit, this%iout, this%array_name, nvals)
    call check_binary_filesize(this%input_unit, expected_size, this%array_name)
    read (this%input_unit, iostat=istat, iomsg=errmsg) &
      (this%int1d(i), i=1, size(this%int1d))
    if (istat /= 0) then
      errmsg = 'Error reading data for array '//trim(this%array_name)// &
               '.  '//trim(errmsg)
      call store_error(errmsg)
      call store_error_unit(this%input_unit)
    end if
  end subroutine read_binary

  subroutine set_factor(this)
    class(Integer1dReaderType) :: this
    this%factor = this%parser%GetInteger()
  end subroutine set_factor

  subroutine apply_factor(this)
    class(Integer1dReaderType) :: this
    integer(I4B) :: i
    if (this%factor /= 0) then
      do i = 1, size(this%int1d)
        this%int1d(i) = this%int1d(i) * this%factor
      end do
    end if
  end subroutine apply_factor

end module Integer1dReaderModule
