module UzrSoilModelFactoryModule
  use KindModule, only: I4B, DP
  use ConstantsModule, only: LENMEMPATH, LENVARNAME
  use SimModule, only: ustop
  use DevFeatureModule, only: dev_feature
  use STLVecIntModule
  use UzrSoilModelModule, only: SoilModelType
  use UzrHaverkampModule, only: HaverkampModelType
  use UzrVanGenuchtenModule, only: VanGenuchtenModelType
  use UzrDevSoilModule, only: UzrDevSoilType
  use BaseDisModule, only: DisBaseType
  implicit none
  private

  ! public access
  public :: init_soil_models, get_soil_model, destroy_soil_models
  public :: NR_SOIL_MODELS, soil_model_name
  public :: BROOKS_COREY, HAVERKAMP, VANGENUCHTEN

  ! the available models
  integer(I4B), parameter :: NR_SOIL_MODELS = 4
  character(len=LENVARNAME), parameter, dimension(NR_SOIL_MODELS) &
    :: soil_model_name = &
       [character(len=LENVARNAME) :: 'BROOKS-COREY', 'HAVERKAMP', &
                                      'VANGENUCHTEN', 'DEVSOILMODEL']
  enum, bind(C)
    enumerator :: BROOKS_COREY = 1 !< Brooks-Corey soil model
    enumerator :: HAVERKAMP = 2 !< Haverkamp soil model
    enumerator :: VANGENUCHTEN = 3 !< Van Genuchten model
    enumerator :: DEVSOILMODEL = 4 !< for development
  end enum

  ! Wrapper to make it possible to store array of pointers
  type :: SoilModelWrapperType
    class(SoilModelType), pointer :: instance => null() !< the wrapped instance
  end type SoilModelWrapperType

  type :: SoilModelFactoryType
    character(len=LENMEMPATH) :: mem_path !< memory path
    character(len=LENMEMPATH) :: input_mem_path !< input context memory path
    integer(I4B) :: nodes !< the number of nodes
    integer(I4B), dimension(:), pointer, contiguous :: nodemap => null() !< points to the node map used for loading from the input context
    real(DP), dimension(:), pointer, contiguous :: porosity => null() !< points to the package porosity parameter
    real(DP), dimension(:), pointer, contiguous :: sat_res => null() !< points to the residual/irreducible saturation used in the package
    type(SoilModelWrapperType), dimension(NR_SOIL_MODELS), private :: models !< all possible soil models (not necessarily active)
  contains
    procedure, private :: create_soil_model
  end type SoilModelFactoryType

  type(SoilModelFactoryType), pointer, private :: soil_model_factory => null() !< the private factory instance for creating soil models

contains

  !> @brief Prepare the factory for creation of soil models
  !<
  subroutine init_soil_models(porosity, sat_res, dis, mem_path, input_mem_path)
    real(DP), dimension(:), pointer, contiguous :: porosity
    real(DP), dimension(:), pointer, contiguous :: sat_res
    class(DisBaseType), pointer :: dis
    character(len=*) :: mem_path
    character(len=*) :: input_mem_path

    allocate (soil_model_factory)
    soil_model_factory%mem_path = mem_path
    soil_model_factory%input_mem_path = input_mem_path
    soil_model_factory%nodes = dis%nodes
    if (dis%nodes < dis%nodesuser) then
      soil_model_factory%nodemap => dis%nodeuser
    end if
    soil_model_factory%porosity => porosity
    soil_model_factory%sat_res => sat_res

  end subroutine init_soil_models

  !> @brief Get a soil model instance by id
  !<
  function get_soil_model(soil_model_id) result(soil_model)
    integer(I4B) :: soil_model_id !< type of soil model
    class(SoilModelType), pointer :: soil_model !< the returned soil model
    ! local
    type(SoilModelFactoryType), pointer :: factory !< convenience

    factory => soil_model_factory
    if (.not. associated(factory%models(soil_model_id)%instance)) then
      soil_model => factory%create_soil_model(soil_model_id)
      factory%models(soil_model_id)%instance => soil_model
    else
      soil_model => factory%models(soil_model_id)%instance
    end if

  end function get_soil_model

  !> @brief Factory method: create a soil model per id
  !<
  function create_soil_model(this, id) result(soil_model)
    class(SoilModelFactoryType), intent(inout) :: this !, the soil model factory
    integer(I4B) :: id !< the id of the soil model
    class(SoilModelType), pointer :: soil_model !< the returned soil model

    soil_model => null()
    select case (id)
    case (BROOKS_COREY)
      !allocate (BrooksCoreyModelType :: soil_model)
    case (HAVERKAMP)
      allocate (HaverkampModelType :: soil_model)
      soil_model%id = HAVERKAMP
    case (VANGENUCHTEN)
      allocate (VanGenuchtenModelType :: soil_model)
      soil_model%id = VANGENUCHTEN
    case (DEVSOILMODEL)
      call dev_feature('Setting the soil model to DEVELOPMENT is a development &
        &feature, install the nightly build or compile from source with &
        &IDEVELOPMODE = 1.')
      allocate (UzrDevSoilType :: soil_model)
      soil_model%id = DEVSOILMODEL
    case default
      write (*, *) "Internal error: invalid soil model id"
      call ustop()
    end select

    call soil_model%create(this%input_mem_path, this%mem_path, &
                           this%nodes, this%nodemap)
    call soil_model%load(this%porosity, this%sat_res)

  end function create_soil_model

  subroutine destroy_soil_models()
    integer(I4B) :: i

    do i = 1, NR_SOIL_MODELS
      if (associated(soil_model_factory%models(i)%instance)) then
        call soil_model_factory%models(i)%instance%destroy()
      end if
    end do
    deallocate (soil_model_factory)

  end subroutine destroy_soil_models

end module UzrSoilModelFactoryModule
