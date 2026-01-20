module Integer2dReaderModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE
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
  public :: read_int2d

  type, extends(ArrayReaderBaseType) :: Integer2dReaderType

    integer(I4B) :: constant_array_value = DZERO
    integer(I4B) :: factor = DONE
    integer(I4B), dimension(:, :), contiguous, pointer :: int2d => null()

  contains

    procedure :: reset_reader
    procedure :: set_constant ! must be overridden
    procedure :: fill_constant ! must be overridden
    procedure :: read_ascii ! must be overridden
    procedure :: read_binary ! must be overridden
    procedure :: set_factor ! must be overridden
    procedure :: apply_factor ! must be overridden

  end type Integer2dReaderType

contains

  subroutine read_int2d(parser, int2d, aname)
    ! -- dummy
    type(BlockParserType), intent(in), target :: parser
    integer(I4B), dimension(:, :), contiguous, target :: int2d
    character(len=*), intent(in) :: aname
    ! -- local
    type(Integer2dReaderType) :: this

    this%parser => parser
    this%int2d => int2d
    this%array_name = aname

    call this%read_array()

  end subroutine read_int2d

  subroutine reset_reader(this)
    class(Integer2dReaderType) :: this
    call this%ArrayReaderBaseType%reset_reader()
    this%constant_array_value = 0
    this%factor = 1
  end subroutine reset_reader

  subroutine set_constant(this)
    class(Integer2dReaderType) :: this
    this%constant_array_value = this%parser%GetInteger()
  end subroutine set_constant

  subroutine fill_constant(this)
    class(Integer2dReaderType) :: this
    integer(I4B) :: i, j
    do i = 1, size(this%int2d, dim=2)
      do j = 1, size(this%int2d, dim=1)
        this%int2d(j, i) = this%constant_array_value
      end do
    end do
  end subroutine fill_constant

  subroutine read_ascii(this)
    class(Integer2dReaderType) :: this
    integer(I4B) :: i, j
    integer(I4B) :: istat
    do i = 1, size(this%int2d, dim=2)
      read (this%input_unit, *, iostat=istat, iomsg=errmsg) &
        (this%int2d(j, i), j=1, size(this%int2d, dim=1))
    end do
    if (istat /= 0) then
      errmsg = 'Error reading data for array '//trim(this%array_name)// &
               '.  '//trim(errmsg)
      call store_error(errmsg)
      call store_error_unit(this%input_unit)
    end if
  end subroutine read_ascii

  subroutine read_binary(this)
    class(Integer2dReaderType) :: this
    integer(I4B) :: i, j
    integer(I4B) :: nvals
    integer(I4B) :: istat
    integer(I4B) :: expected_size
    expected_size = BINARY_HEADER_BYTES + (size(this%int2d) * BINARY_INT_BYTES)
    call read_binary_header(this%input_unit, this%iout, this%array_name, nvals)
    call check_binary_filesize(this%input_unit, expected_size, this%array_name)
    read (this%input_unit, iostat=istat, iomsg=errmsg) &
      ((this%int2d(j, i), j=1, size(this%int2d, dim=1)), &
       i=1, size(this%int2d, dim=2))
    if (istat /= 0) then
      errmsg = 'Error reading data for array '//trim(this%array_name)// &
               '.  '//trim(errmsg)
      call store_error(errmsg)
      call store_error_unit(this%input_unit)
    end if
  end subroutine read_binary

  subroutine set_factor(this)
    class(Integer2dReaderType) :: this
    this%factor = this%parser%GetInteger()
  end subroutine set_factor

  subroutine apply_factor(this)
    class(Integer2dReaderType) :: this
    integer(I4B) :: i, j
    if (this%factor /= DZERO) then
      do i = 1, size(this%int2d, dim=2)
        do j = 1, size(this%int2d, dim=1)
          this%int2d(j, i) = this%int2d(j, i) * this%factor
        end do
      end do
    end if
  end subroutine apply_factor

end module Integer2dReaderModule
