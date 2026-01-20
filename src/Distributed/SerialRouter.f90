module SerialRouterModule
  use RouterBaseModule
  use KindModule, only: I4B
  use VirtualSolutionModule
  implicit none
  private

  public :: create_serial_router

  !> @brief Serial router: currently doesn't do anything
  !<
  type, public, extends(RouterBaseType) :: SerialRouterType
  contains
    procedure :: initialize => sr_initialize
    procedure :: route_all => sr_route_all
    procedure :: route_sln => sr_route_sln
    procedure :: finalize => sr_finalize
    procedure :: destroy => sr_destroy
  end type SerialRouterType

contains

  !> Factory method to create serial router
  !<
  function create_serial_router() result(router)
    class(RouterBaseType), pointer :: router
    ! local
    class(SerialRouterType), pointer :: serial_router

    allocate (serial_router)
    router => serial_router

  end function create_serial_router

  subroutine sr_initialize(this)
    class(SerialRouterType) :: this
  end subroutine sr_initialize

  subroutine sr_route_all(this, stage)
    class(SerialRouterType) :: this
    integer(I4B) :: stage

  end subroutine sr_route_all

  subroutine sr_route_sln(this, virtual_sol, stage)
    class(SerialRouterType) :: this
    type(VirtualSolutionType) :: virtual_sol
    integer(I4B) :: stage

  end subroutine sr_route_sln

  subroutine sr_finalize(this)
    class(SerialRouterType) :: this
  end subroutine sr_finalize

  subroutine sr_destroy(this)
    class(SerialRouterType) :: this
  end subroutine sr_destroy

end module SerialRouterModule
