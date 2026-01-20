module DisvGeom

  use KindModule, only: DP, I4B
  use GeomUtilModule, only: get_node, get_jk
  implicit none
  private
  public :: DisvGeomType
  public :: line_unit_vector, shared_edge

  type DisvGeomType

    integer(I4B) :: k
    integer(I4B) :: j
    integer(I4B) :: nodeusr
    integer(I4B) :: nodered
    integer(I4B) :: nlay
    integer(I4B) :: ncpl
    logical :: reduced
    integer(I4B) :: nodes ! number of reduced nodes; nodes = nlay *ncpl when grid is NOT reduced
    real(DP) :: top
    real(DP) :: bot
    real(DP), pointer, dimension(:) :: top_grid => null()
    real(DP), pointer, dimension(:) :: bot_grid => null()
    integer(I4B), pointer, dimension(:) :: iavert => null()
    integer(I4B), pointer, dimension(:) :: javert => null()
    real(DP), pointer, dimension(:, :) :: vertex_grid => null()
    real(DP), pointer, dimension(:, :) :: cellxy_grid => null()
    integer(I4B), pointer, dimension(:, :) :: nodereduced => null() ! nodered = nodereduced(nodeusr)
    integer(I4B), pointer, dimension(:) :: nodeuser => null() ! nodeusr = nodesuser(nodered)

  contains

    procedure :: init
    generic :: set => set_kj, set_nodered
    procedure :: set_kj
    procedure :: set_nodered
    procedure :: cell_setup
    procedure :: cprops
    procedure :: edge_normal
    procedure :: connection_vector
    procedure :: shares_edge
    procedure :: get_area

  end type DisvGeomType

contains

  !> @brief Initialize
  !<
  subroutine init(this, nlay, ncpl, nodes, top_grid, bot_grid, iavert, &
                  javert, vertex_grid, cellxy_grid, nodereduced, nodeuser)
    ! -- dummy
    class(DisvGeomType) :: this
    integer(I4B), intent(in) :: nlay
    integer(I4B), intent(in) :: ncpl
    integer(I4B), intent(in) :: nodes
    real(DP), dimension(nodes), target :: top_grid
    real(DP), dimension(nodes), target :: bot_grid
    integer(I4B), dimension(:), target :: iavert
    integer(I4B), dimension(:), target :: javert
    real(DP), dimension(:, :), target :: vertex_grid
    real(DP), dimension(:, :), target :: cellxy_grid
    integer(I4B), dimension(ncpl, nlay), target :: nodereduced
    integer(I4B), dimension(nodes), target :: nodeuser
    ! -- local
    integer(I4B) :: nodesuser
    !
    this%nlay = nlay
    this%ncpl = ncpl
    this%nodes = nodes
    this%top_grid => top_grid
    this%bot_grid => bot_grid
    this%iavert => iavert
    this%javert => javert
    this%vertex_grid => vertex_grid
    this%cellxy_grid => cellxy_grid
    this%nodereduced => nodereduced
    this%nodeuser => nodeuser
    nodesuser = ncpl * nlay
    !
    if (nodes < nodesuser) then
      this%reduced = .true.
    else
      this%reduced = .false.
    end if
  end subroutine init

  !> @brief Set node IDs
  !<
  subroutine set_kj(this, k, j)
    ! -- dummy
    class(DisvGeomType) :: this
    integer(I4B), intent(in) :: k
    integer(I4B), intent(in) :: j
    !
    this%k = k
    this%j = j
    this%nodeusr = get_node(k, 1, j, this%nlay, 1, this%ncpl)
    if (this%reduced) then
      this%nodered = this%nodereduced(k, j)
    else
      this%nodered = this%nodeusr
    end if
    call this%cell_setup()
  end subroutine set_kj

  !> @brief Set reduced node number
  !<
  subroutine set_nodered(this, nodered)
    ! -- dummy
    class(DisvGeomType) :: this
    integer(I4B), intent(in) :: nodered
    !
    this%nodered = nodered
    !
    if (this%reduced) then
      this%nodeusr = this%nodeuser(nodered)
    else
      this%nodeusr = nodered
    end if
    !
    call get_jk(this%nodeusr, this%ncpl, this%nlay, this%j, this%k)
    call this%cell_setup()
  end subroutine set_nodered

  !> @brief Set top and bottom elevations of grid cell
  !<
  subroutine cell_setup(this)
    ! -- dummy
    class(DisvGeomType) :: this
    !
    this%top = this%top_grid(this%nodered)
    this%bot = this%bot_grid(this%nodered)
  end subroutine cell_setup

  subroutine cprops(this, cell2, hwva, cl1, cl2, ax, ihc)
    ! -- module
    use ConstantsModule, only: DZERO, DHALF, DONE
    ! -- dummy
    class(DisvGeomType) :: this
    type(DisvGeomType) :: cell2
    real(DP), intent(out) :: hwva
    real(DP), intent(out) :: cl1
    real(DP), intent(out) :: cl2
    real(DP), intent(out) :: ax
    integer(I4B), intent(out) :: ihc
    ! -- local
    integer(I4B) :: ivert1, ivert2
    integer(I4B) :: istart1, istart2, istop1, istop2
    real(DP) :: x0, y0, x1, y1, x2, y2
    !
    if (this%j == cell2%j) then
      !
      ! -- Cells share same j index, so must be a vertical connection
      ihc = 0
      hwva = this%get_area()
      cl1 = DHALF * (this%top - this%bot)
      cl2 = DHALF * (cell2%top - cell2%bot)
      ax = DZERO
    else
      !
      ! -- Must be horizontal connection
      ihc = 1
      istart1 = this%iavert(this%j)
      istop1 = this%iavert(this%j + 1) - 1
      istart2 = cell2%iavert(cell2%j)
      istop2 = this%iavert(cell2%j + 1) - 1
      call shared_edge(this%javert(istart1:istop1), &
                       this%javert(istart2:istop2), &
                       ivert1, ivert2)
      if (ivert1 == 0 .or. ivert2 == 0) then
        !
        ! -- Cells do not share an edge
        hwva = DZERO
        cl1 = DONE
        cl2 = DONE
      else
        x1 = this%vertex_grid(1, ivert1)
        y1 = this%vertex_grid(2, ivert1)
        x2 = this%vertex_grid(1, ivert2)
        y2 = this%vertex_grid(2, ivert2)
        hwva = distance(x1, y1, x2, y2)
        !
        ! -- cl1
        x0 = this%cellxy_grid(1, this%j)
        y0 = this%cellxy_grid(2, this%j)
        cl1 = distance_normal(x0, y0, x1, y1, x2, y2)
        !
        ! -- cl2
        x0 = this%cellxy_grid(1, cell2%j)
        y0 = this%cellxy_grid(2, cell2%j)
        cl2 = distance_normal(x0, y0, x1, y1, x2, y2)
        !
        ! -- anglex
        x1 = this%vertex_grid(1, ivert1)
        y1 = this%vertex_grid(2, ivert1)
        x2 = this%vertex_grid(1, ivert2)
        y2 = this%vertex_grid(2, ivert2)
        ax = anglex(x1, y1, x2, y2)
      end if
    end if
  end subroutine cprops

  !> @brief Return the x and y components of an outward normal facing vector
  !<
  subroutine edge_normal(this, cell2, xcomp, ycomp)
    ! -- module
    use ConstantsModule, only: DZERO, DHALF, DONE
    ! -- dummy
    class(DisvGeomType) :: this
    type(DisvGeomType) :: cell2
    real(DP), intent(out) :: xcomp
    real(DP), intent(out) :: ycomp
    ! -- local
    integer(I4B) :: ivert1, ivert2
    integer(I4B) :: istart1, istart2, istop1, istop2
    real(DP) :: x1, y1, x2, y2
    !
    istart1 = this%iavert(this%j)
    istop1 = this%iavert(this%j + 1) - 1
    istart2 = cell2%iavert(cell2%j)
    istop2 = this%iavert(cell2%j + 1) - 1
    call shared_edge(this%javert(istart1:istop1), &
                     this%javert(istart2:istop2), &
                     ivert1, ivert2)
    x1 = this%vertex_grid(1, ivert1)
    y1 = this%vertex_grid(2, ivert1)
    x2 = this%vertex_grid(1, ivert2)
    y2 = this%vertex_grid(2, ivert2)
    !
    call line_unit_normal(x1, y1, x2, y2, xcomp, ycomp)
  end subroutine edge_normal

  !> @brief Return the x y and z components of a unit vector that points from
  !! from the center of this to the center of cell2, and the straight-line
  !! connection length
  !<
  subroutine connection_vector(this, cell2, nozee, satn, satm, xcomp, &
                               ycomp, zcomp, conlen)
    ! -- module
    use ConstantsModule, only: DZERO, DHALF, DONE
    ! -- dummy
    class(DisvGeomType) :: this
    type(DisvGeomType) :: cell2
    logical, intent(in) :: nozee
    real(DP), intent(in) :: satn
    real(DP), intent(in) :: satm
    real(DP), intent(out) :: xcomp
    real(DP), intent(out) :: ycomp
    real(DP), intent(out) :: zcomp
    real(DP), intent(out) :: conlen
    ! -- local
    real(DP) :: x1, y1, z1, x2, y2, z2
    !
    x1 = this%cellxy_grid(1, this%j)
    y1 = this%cellxy_grid(2, this%j)
    x2 = this%cellxy_grid(1, cell2%j)
    y2 = this%cellxy_grid(2, cell2%j)
    if (nozee) then
      z1 = DZERO
      z2 = DZERO
    else
      z1 = this%bot + DHALF * satn * (this%top - this%bot)
      z2 = cell2%bot + DHALF * satm * (cell2%top - cell2%bot)
    end if
    !
    call line_unit_vector(x1, y1, z1, x2, y2, z2, xcomp, ycomp, zcomp, &
                          conlen)
  end subroutine connection_vector

  !> @brief Return true if this shares a horizontal edge with cell2
  !<
  function shares_edge(this, cell2) result(l)
    ! -- dummy
    class(DisvGeomType) :: this
    type(DisvGeomType) :: cell2
    ! -- return
    logical l
    ! -- local
    integer(I4B) :: istart1, istop1, istart2, istop2
    integer(I4B) :: ivert1, ivert2
    !
    istart1 = this%iavert(this%j)
    istop1 = this%iavert(this%j + 1) - 1
    istart2 = cell2%iavert(cell2%j)
    istop2 = this%iavert(cell2%j + 1) - 1
    call shared_edge(this%javert(istart1:istop1), &
                     this%javert(istart2:istop2), &
                     ivert1, ivert2)
    l = .true.
    if (ivert1 == 0 .or. ivert2 == 0) then
      l = .false.
    end if
  end function shares_edge

  !> @brief Find two common vertices shared by cell1 and cell2.
  !!
  !! Return 0 if there are no shared edges.  Proceed forward through ivlist1
  !! and backward through ivlist2 as a clockwise face in cell1 must correspond
  !! to a counter clockwise face in cell2.
  !<
  subroutine shared_edge(ivlist1, ivlist2, ivert1, ivert2)
    ! -- dummy
    integer(I4B), dimension(:) :: ivlist1
    integer(I4B), dimension(:) :: ivlist2
    integer(I4B), intent(out) :: ivert1
    integer(I4B), intent(out) :: ivert2
    ! -- local
    integer(I4B) :: nv1
    integer(I4B) :: nv2
    integer(I4B) :: il1
    integer(I4B) :: il2
    logical :: found
    !
    found = .false.
    nv1 = size(ivlist1)
    nv2 = size(ivlist2)
    ivert1 = 0
    ivert2 = 0
    outerloop: do il1 = 1, nv1 - 1
      do il2 = nv2, 2, -1
        if (ivlist1(il1) == ivlist2(il2) .and. &
            ivlist1(il1 + 1) == ivlist2(il2 - 1)) then
          found = .true.
          ivert1 = ivlist1(il1)
          ivert2 = ivlist1(il1 + 1)
          exit outerloop
        end if
      end do
      if (found) exit
    end do outerloop
  end subroutine shared_edge

  !> @brief Calculate and return the area of the cell
  !!
  !! a = 1/2 *[(x1*y2 + x2*y3 + x3*y4 + ... + xn*y1) -
  !!             (x2*y1 + x3*y2 + x4*y3 + ... + x1*yn)]
  !<
  function get_area(this) result(area)
    ! -- module
    use ConstantsModule, only: DZERO, DHALF
    ! -- dummy
    class(DisvGeomType) :: this
    ! -- return
    real(DP) :: area
    ! -- local
    integer(I4B) :: ivert
    integer(I4B) :: nvert
    integer(I4B) :: icount
    integer(I4B) :: iv1
    real(DP) :: x
    real(DP) :: y
    real(DP) :: x1
    real(DP) :: y1
    !
    area = DZERO
    nvert = this%iavert(this%j + 1) - this%iavert(this%j)
    icount = 1
    iv1 = this%javert(this%iavert(this%j))
    x1 = this%vertex_grid(1, iv1)
    y1 = this%vertex_grid(2, iv1)
    do ivert = this%iavert(this%j), this%iavert(this%j + 1) - 1
      x = this%vertex_grid(1, this%javert(ivert))
      if (icount < nvert) then
        y = this%vertex_grid(2, this%javert(ivert + 1))
      else
        y = this%vertex_grid(2, this%javert(this%iavert(this%j)))
      end if
      area = area + (x - x1) * (y - y1)
      icount = icount + 1
    end do
    !
    icount = 1
    do ivert = this%iavert(this%j), this%iavert(this%j + 1) - 1
      y = this%vertex_grid(2, this%javert(ivert))
      if (icount < nvert) then
        x = this%vertex_grid(1, this%javert(ivert + 1))
      else
        x = this%vertex_grid(1, this%javert(this%iavert(this%j)))
      end if
      area = area - (x - x1) * (y - y1)
      icount = icount + 1
    end do
    !
    area = abs(area) * DHALF
  end function get_area

  !> @brief Calculate the angle that the x-axis makes with a line that is
  !! normal to the two points.
  !!
  !! This assumes that vertices are numbered clockwise so that the angle is for
  !! the normal outward of cell n.
  !<
  function anglex(x1, y1, x2, y2) result(ax)
    ! -- modules
    use ConstantsModule, only: DZERO, DTWO, DPI
    ! -- dummy
    real(DP), intent(in) :: x1
    real(DP), intent(in) :: x2
    real(DP), intent(in) :: y1
    real(DP), intent(in) :: y2
    ! -- return
    real(DP) :: ax
    ! -- local
    real(DP) :: dx
    real(DP) :: dy
    !
    dx = x2 - x1
    dy = y2 - y1
    ax = atan2(dx, -dy)
    if (ax < DZERO) ax = DTWO * DPI + ax
  end function anglex

  !> @brief Calculate distance between two points
  !<
  function distance(x1, y1, x2, y2) result(d)
    ! -- dummy
    real(DP), intent(in) :: x1
    real(DP), intent(in) :: x2
    real(DP), intent(in) :: y1
    real(DP), intent(in) :: y2
    ! -- return
    real(DP) :: d
    !
    d = (x1 - x2)**2 + (y1 - y2)**2
    d = sqrt(d)
  end function distance

  !> @brief Calculate normal distance from point (x0, y0) to line defined by
  !! two points, (x1, y1), (x2, y2).
  !<
  function distance_normal(x0, y0, x1, y1, x2, y2) result(d)
    ! -- dummy
    real(DP), intent(in) :: x0
    real(DP), intent(in) :: y0
    real(DP), intent(in) :: x1
    real(DP), intent(in) :: y1
    real(DP), intent(in) :: x2
    real(DP), intent(in) :: y2
    ! -- return
    real(DP) :: d
    !
    d = abs((x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1))
    d = d / distance(x1, y1, x2, y2)
  end function distance_normal

  !> @brief Calculate the normal vector components (xcomp and ycomp) for a line
  !! defined by two points, (x0, y0), (x1, y1).
  !<
  subroutine line_unit_normal(x0, y0, x1, y1, xcomp, ycomp)
    ! -- dummy
    real(DP), intent(in) :: x0
    real(DP), intent(in) :: y0
    real(DP), intent(in) :: x1
    real(DP), intent(in) :: y1
    real(DP), intent(out) :: xcomp
    real(DP), intent(out) :: ycomp
    ! -- local
    real(DP) :: dx, dy, vmag
    !
    dx = x1 - x0
    dy = y1 - y0
    vmag = sqrt(dx**2 + dy**2)
    xcomp = -dy / vmag
    ycomp = dx / vmag
  end subroutine line_unit_normal

  !> @brief Calculate the vector components (xcomp, ycomp, and zcomp) for a
  !! line defined by two points, (x0, y0, z0), (x1, y1, z1).
  !!
  !! Also return the magnitude of the original vector, vmag.
  !<
  subroutine line_unit_vector(x0, y0, z0, x1, y1, z1, &
                              xcomp, ycomp, zcomp, vmag)
    ! -- dummy
    real(DP), intent(in) :: x0
    real(DP), intent(in) :: y0
    real(DP), intent(in) :: z0
    real(DP), intent(in) :: x1
    real(DP), intent(in) :: y1
    real(DP), intent(in) :: z1
    real(DP), intent(out) :: xcomp
    real(DP), intent(out) :: ycomp
    real(DP), intent(out) :: zcomp
    real(DP) :: vmag
    ! -- local
    real(DP) :: dx, dy, dz
    !
    dx = x1 - x0
    dy = y1 - y0
    dz = z1 - z0
    vmag = sqrt(dx**2 + dy**2 + dz**2)
    xcomp = dx / vmag
    ycomp = dy / vmag
    zcomp = dz / vmag
  end subroutine line_unit_vector

end module DisvGeom
