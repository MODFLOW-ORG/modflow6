module BoundaryFacesModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, DSAME
  use SimModule, only: store_warning

  use BaseDisModule, only: DisBaseType
  use DisuModule, only: DisuType
  use DisInfoModule, only: number_connected_faces, number_boundary_faces, &
                           number_faces, cell_centroid
  use MathUtilModule, only: cross_product

  implicit none
  private
  public :: BoundaryFacesType

  type BoundaryFacesType
    ! public
    integer(I4B), public, dimension(:), allocatable :: ia
    ! private
    type(BoundaryFaceType), private, allocatable, dimension(:) :: faces
    real(DP), private, allocatable, dimension(:, :) :: normal ! Normal vector of the boundary
  contains
    ! public
    procedure, public :: get_normal
    ! private
    procedure, private :: create_boundary_cells
  end type BoundaryFacesType

  type, private :: BoundaryFaceType
    ! public
    real(DP), dimension(3) :: normal
  end type BoundaryFaceType

  interface BoundaryFacesType
    module procedure constructor
  end interface BoundaryFacesType

contains

  function constructor(dis) result(boundary_cells)
    ! -- return
    type(BoundaryFacesType) :: boundary_cells
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    logical :: has_warnings

    has_warnings = .false.
    select type (dis)
    class is (DisuType)
      ! Check if the discretization package is compatible with the TVD scheme
      if (dis%icondir == 0) then
        has_warnings = .true.
        call store_warning('Vertices not specified for discretization ' // &
        'package, but TVD is active. Vertices must be specified in ' // &
        'discretization package in order to use TVD. Failing to do so may ' // &
        'result in inaccurate results at the boundary.')
      end if
      if (dis%iangledegx == 0) then
        has_warnings = .true.
        call store_warning('ANGLDEGX not specified for discretization ' // &
        'package, but TVD is active. ANGLDEGX must be specified in ' // &
        'discretization package in order to use TVD. Failing to do so may ' // &
        'result in inaccurate results at the boundary.')
      end if

      ! If the discretization doesn't has all the info needed we don't create the boundary faces.
      ! This may impact the results because the cell gradient at the boundary is not computed
      ! correctly.
      if (has_warnings) then  
        call create_empty_boundary_cells(boundary_cells, dis)
        return
      end if
    end select

    call create_boundary_cells(boundary_cells, dis)

  end function constructor

  function get_normal(this, ipos) result(normal)
    ! -- return
    real(DP), dimension(3) :: normal
    ! -- dummy
    class(BoundaryFacesType) :: this
    integer(I4B), intent(in) :: ipos

    normal = this%faces(ipos)%normal
  end function get_normal

  subroutine create_empty_boundary_cells(this, dis)
    ! -- dummy
    class(BoundaryFacesType) :: this
    class(DisBaseType), intent(in) :: dis
    ! -- local
    integer(I4B) :: nodes

    nodes = dis%nodes
    allocate (this%ia(nodes + 1))
    this%ia = 1
  end subroutine create_empty_boundary_cells

  subroutine create_boundary_cells(this, dis)
    ! -- dummy
    class(BoundaryFacesType) :: this
    class(DisBaseType), intent(in) :: dis
    ! -- local
    integer(I4B) :: nodes
    integer(I4B) :: n, ipos
    integer(I4B) :: boundary_face_count, local_boundary_face_count
    real(DP), dimension(:, :, :), allocatable :: boundary_edges
    integer(I4B) :: boundary_cell_id, size_boundary_edges
    real(DP), dimension(:), allocatable :: v1, v2
    real(DP), dimension(3) :: centroid

    ! Allocate memory
    nodes = dis%nodes
    boundary_face_count = number_boundary_faces(dis)
    allocate (this%ia(nodes + 1))
    allocate (this%faces(boundary_face_count))

    ! Create boundary cells
    this%ia(1) = 1
    boundary_cell_id = 1
    do n = 1, nodes
      ! Increment ia for the new boundary cells
      local_boundary_face_count = number_boundary_faces(dis, n)
      this%ia(n + 1) = this%ia(n) + local_boundary_face_count

      if (local_boundary_face_count == 0) cycle

      ! Get the edges belonging to the ghost cells
      boundary_edges = find_boundary_edges(dis, n)

      ! Compute the normal of the boundary
      centroid = cell_centroid(dis, n)
      size_boundary_edges = size(boundary_edges, dim=1)

      do ipos = 1, size_boundary_edges
        v1 = boundary_edges(ipos, 1, :)
        v2 = boundary_edges(ipos, 2, :)

        this%faces(boundary_cell_id)%normal = compute_normal(v1, v2, centroid)

        boundary_cell_id = boundary_cell_id + 1
      end do
    end do
  end subroutine create_boundary_cells

  function compute_normal(v1, v2, centroid) result(normal)
    ! -- return
    real(DP), dimension(3) :: normal
    ! -- dummy
    real(DP), dimension(3), intent(in) :: v1, v2, centroid
    ! -- local

    normal = compute_point_to_line_vector(v1, v2, centroid)

    if (norm2(normal) < DSAME) then
      ! If the normal is zero then the cell centroid is on the boundary.
      ! This occurs when using Voronoi cells.
      ! In this case we compute the normal using the cross product.
      normal = cross_product(v2 - v1, [0.0_dp, 0.0_dp, -1.0_dp])
    end if

    normal = normal / norm2(normal)
  end function compute_normal

  function compute_point_to_line_vector(v1, v2, p) result(vec)
    ! -- return
    real(DP), dimension(3) :: vec
    ! -- dummy
    real(DP), dimension(3) :: v1, v2, p
    ! -- local
    real(DP), dimension(3) :: PV1, V21

    V21 = v2 - v1
    PV1 = p - v1
    vec = dot_product(V21, PV1) * V21 - PV1
  end function compute_point_to_line_vector

  function find_boundary_edges(dis, n) result(boundary_edges)
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    real(DP), allocatable :: boundary_edges(:, :, :)
    ! -- local
    real(DP), dimension(:, :), allocatable :: polyverts
    integer(I4B) :: m, num_sides, iedge, ipos
    logical, dimension(:), allocatable :: is_ghost_boundary
    real(DP) :: x_dir, y_dir, z_dir, prod
    real(DP), dimension(3) :: edge_dir, cross
    integer(I4B) :: num_ghostcells
    integer(I4B) :: isympos, ihc
    real(DP) :: node_z

    call dis%get_polyverts(n, polyverts)
    num_sides = size(polyverts, 2)

    ! We start of by flagging all boundaries as ghost cells.
    ! Later we will check if they are connected to a real cell
    allocate (is_ghost_boundary(num_sides + 2))
    is_ghost_boundary = .true.
    node_z = (dis%top(n) + dis%bot(n)) / 2.0_dp

    ! Check to which edge the connection belongs.
    ! When we find the edge we unflag it as a ghost cell
    do ipos = dis%con%ia(n) + 1, dis%con%ia(n + 1) - 1
      m = dis%con%ja(ipos)
      call dis%connection_normal(n, m, 1, x_dir, y_dir, z_dir, ipos)

      ! Determine connection direction.
      ! Horizontal and vertical connections are treated differently
      isympos = dis%con%jas(ipos)
      ihc = dis%con%ihc(isympos)

      if (ihc == 1) then
        ! Check if the horizontal connection is a ghost cell
        ! We first find the edge to which the connection belongs
        ! Then we take the cross product of the edge direction and the z unit vector
        ! This will give a vector that is perpendicular to the edge
        ! If the dot product of this vector and the connection direction is 1
        ! the connection is not a ghost cell
        do iedge = 1, num_sides
          edge_dir(1:2) = polyverts(:, mod(iedge, num_sides) + 1) &
                          - polyverts(:, iedge)
          edge_dir(3) = DZERO
          edge_dir = edge_dir / norm2(edge_dir)
          cross = cross_product(edge_dir, [0.0_dp, 0.0_dp, -1.0_dp])

          prod = dot_product(cross, [x_dir, y_dir, z_dir])
          if (abs(prod - DONE) < DSAME) then
            is_ghost_boundary(iedge) = .false.
            exit
          end if
        end do
      else
        ! Check if the vertical connection is a ghost cell
        ! By definition the last two edges are the bottom and top
        if (abs(dis%bot(n) - dis%top(m)) < DSAME) then
          is_ghost_boundary(num_sides + 1) = .false.
        end if

        if (abs(dis%top(n) - dis%bot(m)) < DSAME) then
          is_ghost_boundary(num_sides + 2) = .false.
        end if

      end if
    end do

    num_ghostcells = count(is_ghost_boundary)
    allocate (boundary_edges(num_ghostcells, 2, 3))

    ! Collect all the ghost cell edges
    ipos = 1
    boundary_edges = 0
    do iedge = 1, num_sides
      if (is_ghost_boundary(iedge)) then
        boundary_edges(ipos, 1, 1:2) = polyverts(:, iedge)
        boundary_edges(ipos, 2, 1:2) = polyverts(:, mod(iedge, num_sides) + 1)
        boundary_edges(ipos, 1:2, 3) = node_z
        ipos = ipos + 1
      end if
    end do

    ! The top and the bottom are not really edges but we treat them as such
    ! We only need a line on the top/bottom to reflect the centroid over
    if (is_ghost_boundary(num_sides + 1)) then
      boundary_edges(ipos, 1, :) = [dis%xc(n), dis%yc(n), dis%bot(n)]
      boundary_edges(ipos, 2, :) = [dis%xc(n), dis%yc(n), dis%bot(n)]
      ipos = ipos + 1
    end if

    if (is_ghost_boundary(num_sides + 2)) then
      boundary_edges(ipos, 1, :) = [dis%xc(n), dis%yc(n), dis%top(n)]
      boundary_edges(ipos, 2, :) = [dis%xc(n), dis%yc(n), dis%top(n)]
      ipos = ipos + 1
    end if

  end function find_boundary_edges

end module BoundaryFacesModule
