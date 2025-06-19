module TspSourceInfoProviderModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO

  use TspSsmModule, only: TspSsmType
  use TspFmiModule, only: TspFmiType
  use PtrHashTableModule, only: PtrHashTableType

  implicit none
  private
  public :: TspSourceInfoProviderType

  type TspSourceInfoProviderType
    private
    type(TspSsmType), pointer :: ssm
    type(TspFmiType), pointer :: fmi
    type(PtrHashTableType) :: source_map
  contains
    procedure :: initialize
    procedure :: has_source !< Check if a source term exists for a given node
    procedure :: get_source !< Get the source term for a given node

    procedure, private :: get_source_id
  end type TspSourceInfoProviderType

  type SourceIdType
    integer(I4B) :: ipackage !< package number
    integer(I4B) :: ientry !< bound number
  end type SourceIdType

  interface TspSourceInfoProviderType
    module procedure constructor
  end interface TspSourceInfoProviderType

contains
  function constructor(ssm, fmi) result(source_info_provider)
    ! -- return
    type(TspSourceInfoProviderType) :: source_info_provider
    ! --dummy
    type(TspSsmType), pointer, intent(in) :: ssm
    type(TspFmiType), pointer, intent(in) :: fmi

    source_info_provider%ssm => ssm
    source_info_provider%fmi => fmi

  end function constructor

  subroutine initialize(this)
    ! -- dummy
    class(TspSourceInfoProviderType), intent(inout) :: this
    ! -- local
    integer(I4B) :: nflowpack !< number of flow packages
    integer(I4B) :: ipackage !< package number
    integer(I4B) :: ientry !< bound number
    integer(I4B) :: nbound !< number of bounds in the package
    integer(I4B) :: n !< node number
    character(len=:), allocatable :: key !< key for the source map
    class(*), pointer :: obj !< source id object
    type(SourceIdType), pointer :: source_id !< source id object

    ! -- return if ssm or fmi is not available
    if (.not. associated(this%ssm)) then
      return
    end if

    if (.not. associated(this%fmi)) then
      return
    end if

    ! -- Loop through flow packages and create node to source id mapping
    nflowpack = this%fmi%nflowpack
    do ipackage = 1, nflowpack
      if (this%fmi%iatp(ipackage) /= 0) cycle

      nbound = this%fmi%gwfpackages(ipackage)%nbound
      do ientry = 1, nbound
        n = this%fmi%gwfpackages(ipackage)%nodelist(ientry)
        if (n <= 0) cycle

        allocate (source_id, &
                  source=SourceIdType(ipackage=ipackage, ientry=ientry))
        obj => source_id

        key = int2str(n)
        if (.not. this%source_map%contains(key)) then
          call this%source_map%add(key, obj)
        else
          ! Can a cell contains multiple sources?
          ! If so we should add it as a list and append other sources
        end if
      end do
    end do

  end subroutine initialize

  function has_source(this, node) result(exists)
    ! -- return
    logical(LGP) :: exists !< boolean indicating if a source exists
    ! -- dummy
    class(TspSourceInfoProviderType), intent(in) :: this
    integer(I4B), intent(in) :: node
    ! -- local
    type(SourceIdType), pointer :: source_id
    real(DP) :: qbnd

    exists = .false.
    source_id => this%get_source_id(node)
    if (.not. associated(source_id)) return

    qbnd = this%fmi%gwfpackages(source_id%ipackage)%get_flow(source_id%ientry)
    if (qbnd < DZERO) return ! Not a source, but a sink

    exists = .true.
  end function has_source

  function get_source(this, node) result(res)
    ! -- dummy
    class(TspSourceInfoProviderType), intent(in) :: this
    integer(I4B), intent(in) :: node
    real(DP) :: res
    ! -- local
    type(SourceIdType), pointer :: source_id !< source id object
    logical(LGP) :: lauxmixed
    integer(I4B) :: nbound_flow !< number of bounds in the flow package

    if (.not. this%has_source(node)) then
      res = DZERO
      return
    else
      source_id => this%get_source_id(node)
      nbound_flow = this%fmi%gwfpackages(source_id%ipackage)%nbound

      call this%ssm%get_ssm_conc( &
        source_id%ipackage, source_id%ientry, nbound_flow, res, lauxmixed)
    end if

  end function get_source

  ! Converts an integer to a string
  function int2str(i) result(str)
    integer(I4B), intent(in) :: i
    character(len=:), allocatable :: str
    character(len=32) :: tmp

    write (tmp, '(I0)') i
    str = trim(tmp)
  end function int2str

  function get_source_id(this, node) result(source_id)
    ! -- return
    type(SourceIdType), pointer :: source_id
    ! -- dummy
    class(TspSourceInfoProviderType), intent(in) :: this
    integer(I4B), intent(in) :: node
    ! -- local
    character(len=:), allocatable :: key !< key for the source map
    class(*), pointer :: obj !< source id object

    key = int2str(node)
    if (.not. this%source_map%contains(key)) then
      source_id => null()
      return
    end if

    obj => this%source_map%get(key)
    select type (obj)
    type is (SourceIdType)
      source_id => obj
    class default
      source_id => null()
    end select
  end function get_source_id

end module TspSourceInfoProviderModule
