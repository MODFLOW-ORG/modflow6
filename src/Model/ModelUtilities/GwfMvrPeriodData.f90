!> @brief This module contains the GwfMvrPeriodDataModule Module
!!
!! This module contains the code for storing and reading
!! stress period data for the GwfMvr Package.
!!
!<
module GwfMvrPeriodDataModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENMEMPATH, LENMODELNAME, LENPACKAGENAME, &
                             LINELENGTH
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use BlockParserModule, only: BlockParserType

  implicit none
  private
  public GwfMvrPeriodDataType

  !> @brief Derived type for GwfMvrPeriodDataType
  !!
  !! This derived type contains information and methods for
  !! the data read for the GwfMvr Package.
  !!
  !<
  type GwfMvrPeriodDataType
    character(len=LENMODELNAME), &
      dimension(:), pointer, contiguous :: mname1 => null() !< provider model name
    character(len=LENPACKAGENAME), &
      dimension(:), pointer, contiguous :: pname1 => null() !< provider package name
    character(len=LENMODELNAME), &
      dimension(:), pointer, contiguous :: mname2 => null() !< receiver model name
    character(len=LENPACKAGENAME), &
      dimension(:), pointer, contiguous :: pname2 => null() !< receiver package name
    integer(I4B), dimension(:), pointer, contiguous :: id1 => null() !< provider reach number
    integer(I4B), dimension(:), pointer, contiguous :: id2 => null() !< receiver reach number
    integer(I4B), dimension(:), pointer, contiguous :: imvrtype => null() !< mover type (1, 2, 3, 4) corresponds to mvrtypes
    real(DP), dimension(:), pointer, contiguous :: value => null() !< factor or rate depending on mvrtype
  contains
    procedure :: construct
    procedure :: read_from_parser
    procedure :: destroy
  end type GwfMvrPeriodDataType

contains

  !> @ brief Construct arrays
  !!
  !! Allocate maximum space for mover input.
  !!
  !<
  subroutine construct(this, maxsize, memoryPath)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwfMvrPeriodDataType) :: this !< GwfMvrPeriodDataType
    integer(I4B), intent(in) :: maxsize !< size of arrays
    character(len=LENMEMPATH), intent(in) :: memoryPath !< memory manager path

    ! -- character arrays
    allocate (this%mname1(maxsize))
    allocate (this%pname1(maxsize))
    allocate (this%mname2(maxsize))
    allocate (this%pname2(maxsize))

    ! -- integer and real
    call mem_allocate(this%id1, maxsize, 'ID1', memoryPath)
    call mem_allocate(this%id2, maxsize, 'ID2', memoryPath)
    call mem_allocate(this%imvrtype, maxsize, 'IMVRTYPE', memoryPath)
    call mem_allocate(this%value, maxsize, 'VALUE', memoryPath)
  end subroutine construct

  !> @ brief Fill the arrays from parser
  !!
  !! Use the provided block parser to fill the input arrays.
  !!
  !<
  subroutine read_from_parser(this, parser, nmvr, modelname)
    ! -- dummy
    class(GwfMvrPeriodDataType) :: this !< GwfMvrPeriodDataType
    type(BlockParserType), intent(inout) :: parser !< block parser
    integer(I4B), intent(out) :: nmvr !< number of mover entries read
    character(len=LENMODELNAME), intent(in) :: modelname !< name of model or empty string
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: maxmvr
    logical :: endOfBlock
    character(len=LINELENGTH) :: line
    character(len=12) :: mvrtype_char
    !
    ! -- Initialize
    i = 1
    maxmvr = size(this%id1)
    !
    ! -- Read each mover entry
    do
      call parser%GetNextLine(endOfBlock)
      if (endOfBlock) exit
      !
      ! -- Raise error if movers exceeds maxmvr
      if (i > maxmvr) then
        call parser%GetCurrentLine(line)
        write (errmsg, '(a,a)') 'Movers exceed MAXMVR on line: ', &
          trim(adjustl(line))
        call store_error(errmsg)
        call parser%StoreErrorUnit()
      end if
      !
      ! -- modelname, package name, id for provider
      if (modelname == '') then
        call parser%GetStringCaps(this%mname1(i))
      else
        this%mname1(i) = modelname
      end if
      call parser%GetStringCaps(this%pname1(i))
      this%id1(i) = parser%GetInteger()
      !
      ! -- modelname, package name, id for receiver
      if (modelname == '') then
        call parser%GetStringCaps(this%mname2(i))
      else
        this%mname2(i) = modelname
      end if
      call parser%GetStringCaps(this%pname2(i))
      this%id2(i) = parser%GetInteger()
      !
      ! -- Mover type and value
      call parser%GetStringCaps(mvrtype_char)
      select case (mvrtype_char)
      case ('FACTOR')
        this%imvrtype(i) = 1
      case ('EXCESS')
        this%imvrtype(i) = 2
      case ('THRESHOLD')
        this%imvrtype(i) = 3
      case ('UPTO')
        this%imvrtype(i) = 4
      case default
        call store_error('Invalid mover type: '//trim(mvrtype_char))
        call parser%StoreErrorUnit()
      end select
      this%value(i) = parser%GetDouble()
      i = i + 1
    end do
    nmvr = i - 1
  end subroutine read_from_parser

  !> @ brief Destroy memory
  !!
  !! Deallocate memory from the memory manager.
  !!
  !<
  subroutine destroy(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GwfMvrPeriodDataType) :: this !< GwfMvrPeriodDataType

    ! -- character arrays
    deallocate (this%mname1)
    deallocate (this%pname1)
    deallocate (this%mname2)
    deallocate (this%pname2)

    ! -- integer and real
    call mem_deallocate(this%id1)
    call mem_deallocate(this%id2)
    call mem_deallocate(this%imvrtype)
    call mem_deallocate(this%value)
  end subroutine destroy

end module GwfMvrPeriodDataModule

