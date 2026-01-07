module SubcellModule

  use CellDefnModule, only: CellDefnType
  use DomainModule, only: DomainType
  implicit none
  private
  public :: SubcellType

  !> @brief A subcell of a cell.
  type, abstract, extends(DomainType) :: SubcellType
    integer, public :: isubcell !< index of subcell in the cell
    integer, public :: icell !< index of cell in the source grid
  contains
    procedure(init), deferred :: init !< initializer
  end type SubcellType

  abstract interface
    subroutine init(this)
      import SubcellType
      class(SubcellType), intent(inout) :: this
    end subroutine init
  end interface

end module SubcellModule
