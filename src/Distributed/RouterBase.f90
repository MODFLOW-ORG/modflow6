module RouterBaseModule
  use KindModule, only: I4B, LGP
  use VirtualSolutionModule
  implicit none
  private

  type, abstract, public :: RouterBaseType
    logical(LGP) :: halo_activated !< when true, the halo has been activated
    integer(I4B) :: nr_virt_solutions !< number of virtual solution to be routed
    !< (allowing more efficient routing of virtual data)
  contains
    procedure(initialize_if), deferred :: initialize
    procedure(route_all_if), deferred :: route_all
    procedure(route_sln_if), deferred :: route_sln
    procedure(finalize_if), deferred :: finalize
    procedure(destroy_if), deferred :: destroy
  end type RouterBaseType

  abstract interface
    subroutine initialize_if(this)
      import RouterBaseType
      class(RouterBaseType) :: this
    end subroutine initialize_if
    subroutine route_all_if(this, stage)
      import RouterBaseType, I4B
      class(RouterBaseType) :: this
      integer(I4B) :: stage
    end subroutine route_all_if
    subroutine route_sln_if(this, virtual_sol, stage)
      import RouterBaseType, VirtualSolutionType, I4B
      class(RouterBaseType) :: this
      type(VirtualSolutionType) :: virtual_sol
      integer(I4B) :: stage
    end subroutine route_sln_if
    subroutine finalize_if(this)
      import RouterBaseType
      class(RouterBaseType) :: this
    end subroutine finalize_if
    subroutine destroy_if(this)
      import RouterBaseType
      class(RouterBaseType) :: this
    end subroutine destroy_if
  end interface

end module RouterBaseModule
