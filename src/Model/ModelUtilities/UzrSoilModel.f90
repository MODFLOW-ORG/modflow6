module UzrSoilModelModule
  use KindModule, only: I4B, DP
  use ConstantsModule, only: LENMEMPATH
  use SimModule, only: store_error, count_errors, ustop
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use MemoryManagerExtModule, only: mem_set_value
  use GwfUzrInputModule, only: GwfUzrParamFoundType

  implicit none
  private

  ! TODO_UZR: split this file in separate models

  public :: SoilModelType

  type, abstract :: SoilModelType
    integer(I4B) :: id !< the unique soil model id
    character(len=LENMEMPATH) :: input_mem_path !< the path for loading from input context
    character(len=LENMEMPATH) :: mem_path !< memory manager path for storing variables
    integer(I4B) :: nodes !< the number of nodes, dimensioning the model arrays
    integer(I4B), dimension(:), pointer, contiguous :: map => null() !< a map for node reduction
    real(DP), dimension(:), pointer, contiguous :: porosity => null() !< the model porosity per node
    real(DP), dimension(:), pointer, contiguous :: sat_res => null() !< the model residual saturation per node
  contains
    procedure :: create
    procedure(load_if), deferred :: load
    procedure :: destroy
    procedure :: destroy_base
    procedure(saturation_if), deferred :: saturation
    procedure(capacity_if), deferred :: capacity
    procedure(krelative_if), deferred :: krelative
  end type SoilModelType

  abstract interface
    subroutine load_if(this, porosity, sat_res)
      import SoilModelType, DP
      class(SoilModelType), intent(inout) :: this
      real(DP), dimension(:), pointer, contiguous :: porosity
      real(DP), dimension(:), pointer, contiguous :: sat_res
    end subroutine
    function saturation_if(this, psi, i) result(s)
      import SoilModelType, DP, I4B
      class(SoilModelType), intent(inout) :: this
      real(DP) :: psi
      integer(I4B) :: i
      real(DP) :: s
    end function
    function capacity_if(this, psi, i) result(Cm)
      import SoilModelType, DP, I4B
      class(SoilModelType), intent(inout) :: this
      real(DP) :: psi
      integer(I4B) :: i
      real(DP) :: Cm
    end function
    function krelative_if(this, psi, i) result(kr)
      import SoilModelType, DP, I4B
      class(SoilModelType), intent(inout) :: this
      real(DP) :: psi
      integer(I4B) :: i
      real(DP) :: kr
    end function
  end interface

contains

  subroutine create(this, input_mem_path, mem_path, nodes, nodemap)
    class(SoilModelType), intent(inout) :: this
    character(len=*) :: input_mem_path
    character(len=*) :: mem_path
    integer(I4B) :: nodes
    integer(I4B), dimension(:), pointer, contiguous :: nodemap

    this%input_mem_path = input_mem_path
    this%mem_path = mem_path
    this%nodes = nodes
    this%map => nodemap

  end subroutine create

  subroutine destroy(this)
    class(SoilModelType), intent(inout) :: this

    call this%destroy_base()

  end subroutine destroy

  subroutine destroy_base(this)
    class(SoilModelType), intent(inout) :: this

    this%map => null()
    this%porosity => null()
    this%sat_res => null()

  end subroutine destroy_base

end module UzrSoilModelModule
