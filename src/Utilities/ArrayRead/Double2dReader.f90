module Double2dReaderModule

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
  public :: read_dbl2d

  type, extends(ArrayReaderBaseType) :: Double2dReaderType

    real(DP) :: constant_array_value = DZERO
    real(DP) :: factor = DONE
    real(DP), dimension(:, :), contiguous, pointer :: dbl2d => null()

  contains

    procedure :: reset_reader
    procedure :: set_constant ! must be overridden
    procedure :: fill_constant ! must be overridden
    procedure :: read_ascii ! must be overridden
    procedure :: read_binary ! must be overridden
    procedure :: set_factor ! must be overridden
    procedure :: apply_factor ! must be overridden

  end type Double2dReaderType

contains

  subroutine read_dbl2d(parser, dbl2d, aname)
    ! -- dummy
    type(BlockParserType), intent(in), target :: parser
    real(DP), dimension(:, :), contiguous, target :: dbl2d
    character(len=*), intent(in) :: aname
    ! -- local
    type(Double2dReaderType) :: this

    this%parser => parser
    this%dbl2d => dbl2d
    this%array_name = aname

    call this%read_array()

  end subroutine read_dbl2d

  subroutine reset_reader(this)
    class(Double2dReaderType) :: this
    call this%ArrayReaderBaseType%reset_reader()
    this%constant_array_value = DZERO
    this%factor = DONE
  end subroutine reset_reader

  subroutine set_constant(this)
    class(Double2dReaderType) :: this
    this%constant_array_value = this%parser%GetDouble()
  end subroutine set_constant

  subroutine fill_constant(this)
    class(Double2dReaderType) :: this
    integer(I4B) :: i, j
    do i = 1, size(this%dbl2d, dim=2)
      do j = 1, size(this%dbl2d, dim=1)
        this%dbl2d(j, i) = this%constant_array_value
      end do
    end do
  end subroutine fill_constant

  subroutine read_ascii(this)
    class(Double2dReaderType) :: this
    integer(I4B) :: i, j
    integer(I4B) :: istat
    do i = 1, size(this%dbl2d, dim=2)
      read (this%input_unit, *, iostat=istat, iomsg=errmsg) &
        (this%dbl2d(j, i), j=1, size(this%dbl2d, dim=1))
    end do
    if (istat /= 0) then
      errmsg = 'Error reading data for array '//trim(this%array_name)// &
               '.  '//trim(errmsg)
      call store_error(errmsg)
      call store_error_unit(this%input_unit)
    end if
  end subroutine read_ascii

  subroutine read_binary(this)
    class(Double2dReaderType) :: this
    integer(I4B) :: i, j
    integer(I4B) :: nvals
    integer(I4B) :: istat
    integer(I4B) :: expected_size
    expected_size = BINARY_HEADER_BYTES + (size(this%dbl2d) * BINARY_DOUBLE_BYTES)
    call read_binary_header(this%input_unit, this%iout, this%array_name, nvals)
    call check_binary_filesize(this%input_unit, expected_size, this%array_name)
    read (this%input_unit, iostat=istat, iomsg=errmsg) &
      ((this%dbl2d(j, i), j=1, size(this%dbl2d, dim=1)), &
       i=1, size(this%dbl2d, dim=2))
    if (istat /= 0) then
      errmsg = 'Error reading data for array '//trim(this%array_name)// &
               '.  '//trim(errmsg)
      call store_error(errmsg)
      call store_error_unit(this%input_unit)
    end if
  end subroutine read_binary

  subroutine set_factor(this)
    class(Double2dReaderType) :: this
    this%factor = this%parser%GetDouble()
  end subroutine set_factor

  subroutine apply_factor(this)
    class(Double2dReaderType) :: this
    integer(I4B) :: i, j
    if (this%factor /= DZERO) then
      do i = 1, size(this%dbl2d, dim=2)
        do j = 1, size(this%dbl2d, dim=1)
          this%dbl2d(j, i) = this%dbl2d(j, i) * this%factor
        end do
      end do
    end if
  end subroutine apply_factor

end module Double2dReaderModule
