module CellRectModule

  use CellModule, only: CellType
  use CellDefnModule, only: CellDefnType, create_defn
  implicit none

  private
  public :: CellRectType
  public :: create_cell_rect

  type, extends(CellType) :: CellRectType
    private
    double precision, public :: dx ! dimension of cell in local x direction
    double precision, public :: dy ! dimension of cell in local y direction
    double precision, public :: dz ! dimension of cell in z direction

    double precision, public :: sinrot ! sine of rotation angle for local (x, y)
    double precision, public :: cosrot ! cosine of rotation angle for local (x, y)

    double precision, public :: vx1 ! west-boundary local-x velocity
    double precision, public :: vx2 ! east-boundary local-x velocity
    double precision, public :: vy1 ! south-boundary local-y velocity
    double precision, public :: vy2 ! north-boundary local-y velocity
    double precision, public :: vz1 ! bottom-boundary z velocity
    double precision, public :: vz2 ! top-boundary z velocity

    integer, public :: ipvOrigin ! origin vertex
    double precision, public :: xOrigin ! model x origin for local (x, y)
    double precision, public :: yOrigin ! model y origin for local (x, y)
    double precision, public :: zOrigin ! model z origin for local z

  contains
    procedure :: destroy => destroy_rect ! destructor for the cell
  end type CellRectType

contains

  !> @brief Create a new rectangular cell
  subroutine create_cell_rect(cell)
    type(CellRectType), pointer :: cell
    allocate (cell)
    call create_defn(cell%defn)
    allocate (cell%type)
    cell%type = 'rect'
  end subroutine create_cell_rect

  !> @brief Destroy the rectangular cell
  subroutine destroy_rect(this)
    class(CellRectType), intent(inout) :: this
    deallocate (this%defn)
    deallocate (this%type)
  end subroutine destroy_rect

end module CellRectModule
