module TspSourceInfoProviderModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO

  use TspSsmModule, only: TspSsmType
  use TspFmiModule, only: TspFmiType

  implicit none
  private
  public :: TspSourceInfoProviderType

  type TspSourceInfoProviderType
    private
    type(TspSsmType), pointer :: ssm
    type(TspFmiType), pointer :: fmi
    integer(DP), pointer :: source(:) => null() !< Pointer to source term array
  contains
    procedure :: has_source !< Check if a source term exists for a given node
    procedure :: get_source !< Get the source term for a given node

    procedure, private :: get_source_value
  end type TspSourceInfoProviderType

  interface TspSourceInfoProviderType
    module procedure Constructor
  end interface TspSourceInfoProviderType

contains
  function Constructor(ssm, fmi) result(source_info_provider)
    ! -- return
    type(TspSourceInfoProviderType) :: source_info_provider
    ! --dummy
    type(TspSsmType), pointer, intent(in) :: ssm
    type(TspFmiType), pointer, intent(in) :: fmi

    source_info_provider%ssm => ssm
    source_info_provider%fmi => fmi

  end function Constructor

  function has_source(this, node) result(exists)
    ! -- dummy
    class(TspSourceInfoProviderType), intent(in) :: this
    integer(I4B), intent(in) :: node
    logical :: exists

    exists = get_source_value(this, node) > DZERO

  end function has_source

  function get_source(this, node) result(res)
    ! -- dummy
    class(TspSourceInfoProviderType), intent(in) :: this
    integer(I4B), intent(in) :: node
    real(DP) :: res

    res = get_source_value(this, node)

  end function get_source

  function get_source_value(this, node) result(res)
    ! -- return
    real(DP) :: res
    ! -- dummy
    class(TspSourceInfoProviderType) :: this !< TspSsmType object
    integer(I4B), intent(in) :: node !< node number
    ! -- local
    integer(I4B) :: ipackage !< package number
    integer(I4B) :: ientry !< bound number
    logical(LGP) :: lauxmixed
    integer(I4B) :: nbound_flow
    integer(I4B) :: nflowpack
    integer(I4B) :: nbound
    integer(I4B) :: n
    real(DP) :: qbnd

    res = DZERO

    nflowpack = this%fmi%nflowpack
    do ipackage = 1, nflowpack
      if (this%fmi%iatp(ipackage) /= 0) cycle

      nbound = this%fmi%gwfpackages(ipackage)%nbound
      do ientry = 1, nbound
        n = this%fmi%gwfpackages(ipackage)%nodelist(ientry)
        if (n <= 0) cycle
        if (n /= node) cycle
        qbnd = this%fmi%gwfpackages(ipackage)%get_flow(ientry)
        if (qbnd < DZERO) cycle ! Not a source, but a sink

        nbound_flow = this%fmi%gwfpackages(ipackage)%nbound
        call this%ssm%get_ssm_conc(ipackage, ientry, nbound_flow, res, lauxmixed)
        return
      end do
    end do

  end function get_source_value

end module TspSourceInfoProviderModule
