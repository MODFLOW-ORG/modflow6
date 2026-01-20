module SubcellTriModule

  use SubcellModule, only: SubcellType
  implicit none

  private
  public :: SubcellTriType
  public :: create_subcell_tri

  type, extends(SubcellType) :: SubcellTriType
    private
    double precision, public :: x0, y0, x1, y1, x2, y2 !< subcell corner coordinates
    double precision, public :: v0x, v0y, v1x, v1y, v2x, v2y !< subcell corner velocities
    double precision, public :: ztop, zbot !< subcell top and bottom elevations
    double precision, public :: dz !< subcell thickness
    double precision, public :: vztop, vzbot !< subcell top and bottom velocities
  contains
    procedure, public :: destroy => destroy_subcell_tri !< destructor for the subcell
    procedure, public :: init => init_subcell_tri !< initializes the triangular subcell
  end type SubcellTriType

contains

  !> @brief Create a new triangular subcell
  subroutine create_subcell_tri(subcell)
    type(SubcellTriType), pointer :: subcell
    allocate (subcell)
    allocate (subcell%type)
    subcell%type = 'subcelltri'
  end subroutine create_subcell_tri

  !> @brief Destructor for a triangular subcell
  subroutine destroy_subcell_tri(this)
    class(SubcellTriType), intent(inout) :: this
    deallocate (this%type)
  end subroutine destroy_subcell_tri

  !> @brief Initialize a triangular subcell
  subroutine init_subcell_tri(this)
    class(SubcellTriType), intent(inout) :: this
  end subroutine init_subcell_tri

end module SubcellTriModule
