module SubcellModule

  use CellDefnModule, only: CellDefnType
  implicit none
  private
  public :: SubcellType

  !> @brief A subcell of a cell.
  type, abstract :: SubcellType
    private
    character(len=40), pointer, public :: type !< character string that names the tracking domain type
    integer, public :: isubcell !< index of subcell in the cell
    integer, public :: icell !< index of cell in the source grid
  contains
    procedure(destroy), deferred :: destroy !< destructor
    procedure(init), deferred :: init !< initializer
  end type SubcellType

  abstract interface
    subroutine destroy(this)
      import SubcellType
      class(SubcellType), intent(inout) :: this
    end subroutine

    subroutine init(this)
      import SubcellType
      class(SubcellType), intent(inout) :: this
    end subroutine init
  end interface

end module SubcellModule
