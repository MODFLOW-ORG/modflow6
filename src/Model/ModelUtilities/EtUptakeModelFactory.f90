module EtUptakeModelFactoryModule
  use KindModule, only: I4B
  use ConstantsModule, only: LENVARNAME
  use SimModule, only: ustop
  use EtUptakeModelModule, only: EtUptakeModelType
  use EtUptakeFeddesModule, only: FeddesUptakeType
  implicit none
  private

  public :: create_uptake_model, destroy_uptake_model
  public :: NR_UPTAKE_MODELS, uptake_model_name
  public :: FEDDES

  ! the available uptake models
  integer(I4B), parameter :: NR_UPTAKE_MODELS = 1
  character(len=LENVARNAME), parameter, dimension(NR_UPTAKE_MODELS) &
    :: uptake_model_name = &
       [character(len=LENVARNAME) :: 'FEDDES']

  enum, bind(C)
    enumerator :: FEDDES = 1 !< Feddes pressure-head stress response
  end enum

contains

  !> @brief Create an uptake model instance by id
  !<
  function create_uptake_model(id) result(model)
    integer(I4B), intent(in) :: id !< the uptake model id
    class(EtUptakeModelType), pointer :: model !< the returned uptake model

    model => null()
    select case (id)
    case (FEDDES)
      allocate (FeddesUptakeType :: model)
      model%id = FEDDES
    case default
      write (*, *) "Internal error: invalid uptake model id"
      call ustop()
    end select

  end function create_uptake_model

  !> @brief Destroy an uptake model instance
  !<
  subroutine destroy_uptake_model(model)
    class(EtUptakeModelType), pointer, intent(inout) :: model

    if (associated(model)) then
      call model%destroy()
      deallocate (model)
      model => null()
    end if

  end subroutine destroy_uptake_model

end module EtUptakeModelFactoryModule
