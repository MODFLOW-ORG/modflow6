module SubcellRectModule

  use SubcellModule, only: SubcellType
  implicit none

  private
  public :: SubcellRectType
  public :: create_subcell_rect

  type, extends(SubcellType) :: SubcellRectType
    private
    double precision, public :: sinrot !< sine of rotation angle for local (x, y)
    double precision, public :: cosrot !< cosine of rotation angle for local (x, y)
    double precision, public :: xOrigin !< cell x origin for local (x, y)
    double precision, public :: yOrigin !< cell y origin for local (x, y)
    double precision, public :: zOrigin !< cell z origin for local z
    double precision, public :: dx, dy, dz !< subcell dimensions
    double precision, public :: vx1, vx2, vy1, vy2, vz1, vz2 !< subcell face velocities
  contains
    procedure, public :: destroy => destroy_subcell_rect !< destructor for the subcell
    procedure, public :: init => init_subcell_rect !< initializes the rectangular subcell
  end type SubcellRectType

contains

  !> @brief Create a new rectangular subcell
  subroutine create_subcell_rect(subcell)
    type(SubcellRectType), pointer :: subcell
    allocate (subcell)
    allocate (subcell%type)
    subcell%type = 'subcellrect'
  end subroutine create_subcell_rect

  !> @brief Destructor for a rectangular subcell
  subroutine destroy_subcell_rect(this)
    class(SubcellRectType), intent(inout) :: this
    deallocate (this%type)
  end subroutine destroy_subcell_rect

  !> @brief Initialize a rectangular subcell
  subroutine init_subcell_rect(this)
    class(SubcellRectType), intent(inout) :: this
  end subroutine init_subcell_rect

end module SubcellRectModule
