module CellModule

  use CellDefnModule, only: CellDefnType
  use KindModule, only: I4B
  implicit none
  private
  public :: CellType

  integer(I4B), public, parameter :: MAX_POLY_CELLS = 10

  !> @brief Base type for grid cells of a concrete type. Contains
  !! a cell-definition which is information shared by cell types.
  type, abstract :: CellType
    character(len=40), pointer :: type ! tracking domain type
    type(CellDefnType), pointer :: defn => null() ! cell defn
  contains
    procedure(destroy), deferred :: destroy !< destroy the cell
  end type CellType

  abstract interface
    subroutine destroy(this)
      import CellType
      class(CellType), intent(inout) :: this
    end subroutine
  end interface

end module CellModule
