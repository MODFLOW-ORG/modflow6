module CellRectQuadModule

  use KindModule, only: DP, I4B, LGP
  use CellModule, only: CellType
  use CellDefnModule, only: CellDefnType, create_defn
  use ConstantsModule, only: DZERO
  implicit none

  private
  public :: CellRectQuadType
  public :: create_cell_rect_quad

  type, extends(CellType) :: CellRectQuadType
    real(DP) :: dx !< dimension of cell in local x direction
    real(DP) :: dy !< dimension of cell in local y direction
    real(DP) :: dz !< dimension of cell in z direction

    real(DP) :: sinrot !< sine of rotation angle for local (x, y)
    real(DP) :: cosrot !< cosine of rotation angle for local (x, y)

    integer(I4B) :: irvOrigin !< origin rectangle vertex
    real(DP) :: xOrigin !< model x origin for local (x, y)
    real(DP) :: yOrigin !< model y origin for local (x, y)
    real(DP) :: zOrigin !< model z origin for local z

    real(DP) :: qextl1(4), qextl2(4), qintl(5) !< external and internal subcell flows for the cell
    integer(I4B), allocatable :: irectvert(:) !< list of indices of the rectangle vertices
    integer(I4B), allocatable :: ipv4irv(:, :) !< list of the polygon vertex indices that correspond to the rectangle vertex indices
    real(DP), allocatable :: rectflow(:, :) !< flow(s) for each rectangle face
  contains
    procedure :: destroy => destroy_cell_rect_quad !< destructor for the cell
    procedure :: init_from !< initializes the cell from an existing cell

    procedure :: load_rect_verts_flows ! loads list of indices of the rectangle vertices and face flows
    procedure :: get_rect_ivert_sw ! gets index of southwest rectangle vertex
    procedure :: get_rect_dim_rot ! gets rectangular dimensions and rotation
    procedure :: get_rect_flow !< returns a rectangle face flow
    procedure :: face_is_refined !< returns whether a rectangle face is refined
  end type CellRectQuadType

contains

  !> @brief Create a new rectangular-quad cell
  subroutine create_cell_rect_quad(cell)
    type(CellRectQuadType), pointer :: cell
    allocate (cell)
    call create_defn(cell%defn)
    allocate (cell%irectvert(5))
    allocate (cell%ipv4irv(2, 4))
    allocate (cell%rectflow(2, 4))
    allocate (cell%type)
    cell%type = 'rectquad'
  end subroutine create_cell_rect_quad

  !> @brief Destroy the rectangular-quad cell
  subroutine destroy_cell_rect_quad(this)
    class(CellRectQuadType), intent(inout) :: this
    deallocate (this%defn)
    deallocate (this%irectvert)
    deallocate (this%type)
  end subroutine destroy_cell_rect_quad

  !> @brief Initialize a rectangular-quad cell from cell definition
  subroutine init_from(this, defn)
    class(CellRectQuadType), intent(inout) :: this
    type(CellDefnType), pointer :: defn
    this%defn => defn
    call this%load_rect_verts_flows()
  end subroutine init_from

  !> @brief Load local polygon vertex indices and rectangular
  !> face flows
  !!
  !! Loads local polygon vertex indices of the four rectangle
  !! vertices and face flows of a rectangular-quad cell.
  !<
  subroutine load_rect_verts_flows(this)
    ! dummy
    class(CellRectQuadType), intent(inout) :: this
    ! local
    integer(I4B) :: n, m

    n = 0
    do m = 1, this%defn%npolyverts
      if (.not. this%defn%get_ispv180(m)) then
        n = n + 1
        this%irectvert(n) = m
        this%ipv4irv(1, n) = m
        this%rectflow(1, n) = this%defn%get_faceflow(m)
        this%ipv4irv(2, n) = 0
        this%rectflow(2, n) = DZERO
      else
        if (n .ne. 0) then
          this%ipv4irv(2, n) = m
          this%rectflow(2, n) = this%defn%get_faceflow(m)
        end if
      end if
    end do

    ! Wrap around for convenience
    this%irectvert(5) = this%irectvert(1)
  end subroutine load_rect_verts_flows

  !> @brief Get index of SW rectangle vertex
  !!
  !! Return the index (1, 2, 3, or 4) of the southwest
  !! rectangle vertex of a rectangular-quad cell
  !<
  function get_rect_ivert_sw(this) result(irv1)
    ! dummy
    class(CellRectQuadType), intent(inout) :: this
    integer(I4B) :: irv1
    ! local
    integer(I4B) :: irv, irv2, irv4, ipv1, ipv2, ipv4
    integer(I4B), dimension(4) :: irvnxt = (/2, 3, 4, 1/)
    real(DP) :: x1, y1, x2, y2, x4, y4

    ! Find the "southwest" rectangle vertex by finding the vertex formed
    ! either by (1) a rectangle edge over which x decreases (going
    ! clockwise) followed by an edge over which x does not increase, or by
    ! (2) a rectangle edge over which y does not decrease (again going
    ! clockwise) followed by a rectangle edge over which y increases. In
    ! the end, ipv1 is the index (1, 2, 3, or 4) of the southwest
    ! rectangle vertex.
    do irv = 1, 4
      irv4 = irv
      irv1 = irvnxt(irv4)
      ipv4 = this%irectvert(irv4)
      ipv1 = this%irectvert(irv1)
      x4 = this%defn%polyvert(1, ipv4)
      y4 = this%defn%polyvert(2, ipv4)
      x1 = this%defn%polyvert(1, ipv1)
      y1 = this%defn%polyvert(2, ipv1)
      if (x1 .lt. x4) then
        irv2 = irvnxt(irv1)
        ipv2 = this%irectvert(irv2)
        x2 = this%defn%polyvert(1, ipv2)
        if (x2 .le. x1) return
      else if (y1 .ge. y4) then
        irv2 = irvnxt(irv1)
        ipv2 = this%irectvert(irv2)
        y2 = this%defn%polyvert(2, ipv2)
        if (y2 .gt. y1) return
      end if
    end do
  end function get_rect_ivert_sw

  !> @brief Get rectangular cell dimensions and rotation
  !!
  !! Compute rectangular dimensions and rotation of
  !! the cell using the specified rectangle vertex
  !! as the origin
  !<
  subroutine get_rect_dim_rot(this)
    ! dummy
    class(CellRectQuadType), intent(inout) :: this
    ! local
    integer(I4B) :: irv2, irv4, ipv1, ipv2, ipv4
    integer(I4B), dimension(4) :: irvnxt = (/2, 3, 4, 1/)
    real(DP) :: x1, y1, x2, y2, x4, y4, dx2, dy2, dx4, dy4

    ! Get rectangle vertex neighbors irv2 and irv4
    irv2 = irvnxt(this%irvOrigin)
    irv4 = irvnxt(irvnxt(irv2))

    ! Get model coordinates at irv1, irv2, and irv4
    ipv1 = this%irectvert(this%irvOrigin)
    x1 = this%defn%polyvert(1, ipv1)
    y1 = this%defn%polyvert(2, ipv1)
    ipv2 = this%irectvert(irv2)
    x2 = this%defn%polyvert(1, ipv2)
    y2 = this%defn%polyvert(2, ipv2)
    ipv4 = this%irectvert(irv4)
    x4 = this%defn%polyvert(1, ipv4)
    y4 = this%defn%polyvert(2, ipv4)

    ! Compute rectangle dimensions
    this%xOrigin = x1
    this%yOrigin = y1
    this%zOrigin = this%defn%bot
    dx2 = x2 - this%xOrigin
    dy2 = y2 - this%yOrigin
    dx4 = x4 - this%xOrigin
    dy4 = y4 - this%yOrigin
    this%dx = dsqrt(dx4 * dx4 + dy4 * dy4)
    this%dy = dsqrt(dx2 * dx2 + dy2 * dy2)
    this%dz = this%defn%top - this%zOrigin

    ! Compute sine and cosine of rotation angle (angle between "southern"
    ! rectangle side irv1-irv4 and the model x axis)
    this%sinrot = dy4 / this%dx
    this%cosrot = dx4 / this%dx
  end subroutine get_rect_dim_rot

  !> @brief Return a rectangle face flow
  function get_rect_flow(this, iq, irv) result(rectflow)
    class(CellRectQuadType), intent(inout) :: this
    integer(I4B) :: iq, irv
    real(DP) :: rectflow
    rectflow = this%rectflow(iq, irv)
  end function get_rect_flow

  !> @brief Return whether a rectangle face is refined
  function face_is_refined(this, i) result(is_refined)
    ! dummy
    class(CellRectQuadType), intent(inout) :: this
    integer(I4B) :: i !< face index
    logical(LGP) :: is_refined

    if (this%ipv4irv(2, i) .ne. 0) then
      is_refined = .true.
    else
      is_refined = .false.
    end if
  end function face_is_refined

end module CellRectQuadModule
