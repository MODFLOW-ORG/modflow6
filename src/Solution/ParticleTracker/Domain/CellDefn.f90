module CellDefnModule
  use KindModule, only: DP, I4B, LGP
  implicit none

  private
  public :: CellDefnType
  public :: create_defn, get_iatop

  !> @brief Base grid cell definition.
  type CellDefnType
    private
    integer(I4B), public :: icell !< index of cell in source grid
    logical(LGP), public :: can_be_rect !< whether cell is representable as a rectangular cell
    logical(LGP), public :: can_be_quad !< whether cell is representable as a rectangular quad cell
    integer(I4B), public :: npolyverts !< number of vertices for cell polygon
    real(DP), public :: porosity !< cell porosity
    real(DP), public :: retfactor !< cell retardation factor
    integer(I4B), public :: ilay !< layer number
    integer(I4B), public :: izone !< cell zone number
    integer(I4B), public :: iweaksink !< weak sink indicator
    integer(I4B), public :: inoexitface !< no exit face indicator
    integer(I4B), public :: iatop !< index of cell top in grid's top/bot arrays (<0 => top array)
    real(DP), public :: top, bot !< top and bottom elevations of cell
    real(DP), public :: sat !< cell saturation
    real(DP), allocatable, public :: polyvert(:, :) !< vertices for cell polygon
    logical(LGP), allocatable, public :: ispv180(:) !< indicator of 180-degree vertices (.true. = 180-degree angle at vertex)
    integer(I4B), allocatable, public :: facenbr(:) !< neighbors that correspond to faces(/vertices)
    real(DP), allocatable, public :: faceflow(:) !< flows that correspond to faces(/vertices)
    real(DP), public :: distflow !< net distributed flow into cell
  contains
    procedure, public :: get_ispv180 !< returns 180-degree indicator for a vertex
    procedure, public :: get_botflow !< returns bottom flow
    procedure, public :: get_topflow !< returns top flow
    procedure, public :: get_distflow !< returns distributed flow
    procedure, public :: get_faceflow !< returns a face flow
  end type CellDefnType

contains

  !> @brief Create a new cell definition object
  subroutine create_defn(cellDefn)
    type(CellDefnType), pointer :: cellDefn
    allocate (cellDefn)
    ! Initially, allocate arrays to size for structured grid tracking method.
    ! They can be (lazily) expanded as necessary for the unstructured method.
    allocate (cellDefn%ispv180(5))
    allocate (cellDefn%facenbr(7))
    allocate (cellDefn%faceflow(7))
  end subroutine create_defn

  !> @brief Get the index corresponding to top elevation of a cell in the grid.
  !! This is -1 if the cell is in the top layer and 1 otherwise.
  function get_iatop(ncpl, icu) result(iatop)
    integer(I4B), intent(in) :: ncpl, icu
    integer(I4B) :: iatop

    if (icu .le. ncpl) then
      iatop = -1
    else
      iatop = 1
    end if
  end function get_iatop

  !> @brief Return 180-degree indicator for a vertex
  function get_ispv180(this, m) result(ispv180)
    class(CellDefnType), intent(inout) :: this
    integer(I4B) :: m
    logical(LGP) :: ispv180
    ispv180 = this%ispv180(m)
  end function get_ispv180

  !> @brief Return the bottom flow
  function get_botflow(this) result(botflow)
    class(CellDefnType), intent(inout) :: this
    real(DP) :: botflow
    botflow = this%faceflow(this%npolyverts + 2)
  end function get_botflow

  !> @brief Return the top flow
  function get_topflow(this) result(topflow)
    class(CellDefnType), intent(inout) :: this
    real(DP) :: topflow
    topflow = this%faceflow(this%npolyverts + 3)
  end function get_topflow

  !> @brief Return the distributed flow
  function get_distflow(this) result(distflow)
    class(CellDefnType), intent(inout) :: this
    real(DP) :: distflow
    distflow = this%distflow
  end function get_distflow

  !> @brief Return a face flow
  function get_faceflow(this, m) result(faceflow)
    class(CellDefnType), intent(inout) :: this
    integer(I4B) :: m
    real(DP) :: faceflow
    faceflow = this%faceflow(m)
  end function get_faceflow

end module CellDefnModule
