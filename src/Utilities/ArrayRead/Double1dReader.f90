module Double1dReaderModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE
  use BlockParserModule, only: BlockParserType
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_unit
  use ArrayReadersModule, only: read_binary_header, &
                                check_binary_filesize, &
                                BINARY_DOUBLE_BYTES, &
                                BINARY_HEADER_BYTES
  use ArrayReaderBaseModule, only: ArrayReaderBaseType

  implicit none
  private
  public :: read_dbl1d

  type, extends(ArrayReaderBaseType) :: Double1dReaderType

    real(DP) :: constant_array_value = DZERO
    real(DP) :: factor = DONE
    real(DP), dimension(:), contiguous, pointer :: dbl1d => null()

  contains

    procedure :: reset_reader
    procedure :: set_constant ! must be overridden
    procedure :: fill_constant ! must be overridden
    procedure :: read_ascii ! must be overridden
    procedure :: read_binary ! must be overridden
    procedure :: set_factor ! must be overridden
    procedure :: apply_factor ! must be overridden

  end type Double1dReaderType

contains

  subroutine read_dbl1d(parser, dbl1d, aname)
    ! -- dummy
    type(BlockParserType), intent(in), target :: parser
    real(DP), dimension(:), contiguous, target :: dbl1d
    character(len=*), intent(in) :: aname
    ! -- local
    type(Double1dReaderType) :: this

    this%parser => parser
    this%dbl1d => dbl1d
    this%array_name = aname

    call this%read_array()

  end subroutine read_dbl1d

  subroutine reset_reader(this)
    class(Double1dReaderType) :: this
    call this%ArrayReaderBaseType%reset_reader()
    this%constant_array_value = DZERO
    this%factor = DONE
  end subroutine reset_reader

  subroutine set_constant(this)
    class(Double1dReaderType) :: this
    this%constant_array_value = this%parser%GetDouble()
  end subroutine set_constant

  subroutine fill_constant(this)
    class(Double1dReaderType) :: this
    integer(I4B) :: i
    do i = 1, size(this%dbl1d)
      this%dbl1d(i) = this%constant_array_value
    end do
  end subroutine fill_constant

  subroutine read_ascii(this)
    class(Double1dReaderType) :: this
    integer(I4B) :: i
    integer(I4B) :: istat
    read (this%input_unit, *, iostat=istat, iomsg=errmsg) &
      (this%dbl1d(i), i=1, size(this%dbl1d))
    if (istat /= 0) then
      errmsg = 'Error reading data for array '//trim(this%array_name)// &
               '.  '//trim(errmsg)
      call store_error(errmsg)
      call store_error_unit(this%input_unit)
    end if
  end subroutine read_ascii

  subroutine read_binary(this)
    class(Double1dReaderType) :: this
    integer(I4B) :: i
    integer(I4B) :: nvals
    integer(I4B) :: istat
    integer(I4B) :: expected_size
    expected_size = BINARY_HEADER_BYTES + (size(this%dbl1d) * BINARY_DOUBLE_BYTES)
    call read_binary_header(this%input_unit, this%iout, this%array_name, nvals)
    call check_binary_filesize(this%input_unit, expected_size, this%array_name)
    read (this%input_unit, iostat=istat, iomsg=errmsg) &
      (this%dbl1d(i), i=1, size(this%dbl1d))
    if (istat /= 0) then
      errmsg = 'Error reading data for array '//trim(this%array_name)// &
               '.  '//trim(errmsg)
      call store_error(errmsg)
      call store_error_unit(this%input_unit)
    end if
  end subroutine read_binary

  subroutine set_factor(this)
    class(Double1dReaderType) :: this
    this%factor = this%parser%GetDouble()
  end subroutine set_factor

  subroutine apply_factor(this)
    class(Double1dReaderType) :: this
    integer(I4B) :: i
    if (this%factor /= DZERO) then
      do i = 1, size(this%dbl1d)
        this%dbl1d(i) = this%dbl1d(i) * this%factor
      end do
    end if
  end subroutine apply_factor

end module Double1dReaderModule
