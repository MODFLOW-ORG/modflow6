module TspIcModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENVARNAME
  use GwfIcModule, only: GwfIcType
  use BlockParserModule, only: BlockParserType
  use BaseDisModule, only: DisBaseType

  implicit none
  private
  public :: TspIcType
  public :: ic_cr

  ! -- Most of the TspIcType functionality comes from GwfIcType
  type, extends(GwfIcType) :: TspIcType
    ! -- strings
    character(len=LENVARNAME) :: depvartype = ''
  end type TspIcType

contains

  !> @brief Create a new initial conditions object
  !<
  subroutine ic_cr(ic, name_model, input_mempath, inunit, iout, dis, depvartype)
    ! -- dummy
    type(TspIcType), pointer :: ic
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    class(DisBaseType), pointer, intent(in) :: dis
    character(len=LENVARNAME), intent(in) :: depvartype
    !
    ! -- Create the object
    allocate (ic)
    !
    ! -- create name and memory path
    call ic%set_names(1, name_model, 'IC', 'IC', input_mempath)
    !
    ! -- Allocate scalars
    call ic%allocate_scalars()
    !
    ic%inunit = inunit
    ic%iout = iout
    !
    ! -- set pointers
    ic%dis => dis
    !
    ! -- Give package access to the assigned labelsd based on dependent variable
    ic%depvartype = depvartype
  end subroutine ic_cr

end module TspIcModule
