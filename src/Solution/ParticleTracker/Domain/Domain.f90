module DomainModule

  implicit none
  private
  public :: DomainType

  !> @brief A tracking domain.
  type, abstract :: DomainType
    character(len=40), pointer :: type !< character string that names the tracking domain type
  contains
    procedure(destroy), deferred :: destroy !< destructor
  end type DomainType

  abstract interface
    subroutine destroy(this)
      import DomainType
      class(DomainType), intent(inout) :: this
    end subroutine
  end interface

end module DomainModule
