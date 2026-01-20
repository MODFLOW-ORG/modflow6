module CellPolyModule

  use CellModule, only: CellType
  use CellDefnModule, only: CellDefnType, create_defn
  implicit none

  private
  public :: CellPolyType
  public :: create_cell_poly

  type, extends(CellType) :: CellPolyType
  contains
    procedure :: destroy => destroy_cell_poly
  end type CellPolyType

contains

  !> @brief Create a new polygonal cell
  subroutine create_cell_poly(cell)
    type(CellPolyType), pointer :: cell
    allocate (cell)
    call create_defn(cell%defn)
    allocate (cell%type)
    cell%type = 'poly'
  end subroutine create_cell_poly

  !> @brief Destroy the polygonal cell
  subroutine destroy_cell_poly(this)
    class(CellPolyType), intent(inout) :: this
    deallocate (this%defn)
    deallocate (this%type)
  end subroutine destroy_cell_poly

end module CellPolyModule
