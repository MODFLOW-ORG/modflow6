module DisInfoModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DSAME, DONE
  use BaseDisModule, only: DisBaseType

  implicit none
  private

  !!> @brief Returns the number of boundary faces for a cell or for the entire discretization.
  !!
  !! This generic interface provides two functionalities:
  !! - When called with a cell index, it returns the number of boundary faces for that specific cell.
  !! - When called without a cell index, it returns the total number of boundary faces in the entire discretization.
  !! A boundary face is defined as a face of a cell that is not connected to any neighboring cell.
  !<
  interface number_boundary_faces
    module procedure number_local_boundary_faces, number_global_boundary_faces
  end interface

  public :: number_connected_faces
  public :: cell_centroid
  public :: node_distance
  public :: number_boundary_faces

  private :: number_unique_connected_faces
  private :: number_faces
  private :: number_local_boundary_faces
  private :: number_global_boundary_faces

contains

  !!> @brief Returns the total number of faces for a given cell.
  !!
  !! This function computes the total number of faces for cell `n` in the discretization.
  !! The number of faces is determined by the number of polygon sides for the cell,
  !! plus two additional faces representing the top and bottom of the cell.
  !<
  function number_faces(dis, n) result(faces_count)
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    integer(I4B) :: faces_count
    ! -- local
    real(DP), dimension(:, :), allocatable :: polyverts

    call dis%get_polyverts(n, polyverts)
    faces_count = size(polyverts, 2) + 2 ! We add 2 for the top and bottom
  end function number_faces

  !!> @brief Returns the number of connected faces for a given cell.
  !!
  !! This function computes the number of faces of cell `n` that are connected to neighboring cells
  !! in the discretization. The value is determined from the connectivity information in the
  !! connection arrays, and does not include boundary faces (faces not connected to another cell).
  !<
  function number_connected_faces(dis, n) result(connected_faces_count)
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    integer(I4B) :: connected_faces_count

    connected_faces_count = dis%con%ia(n + 1) - dis%con%ia(n) - 1
  end function number_connected_faces

  !!> @brief Returns the number of unique connected faces for a given cell.
  !!
  !! This function computes the number of unique faces of cell `n` that are connected to neighboring cells
  !! in the discretization. In some grid types, a single face can be connected to multiple other cells
  !! (i.e., there can be multiple connections per face). This function returns the number of distinct faces
  !! that are connected, counting each face only once regardless of how many cells it is connected to.
  !! This is useful for identifying the true number of connected faces without counting duplicate connections.
  !<
  function number_unique_connected_faces(dis, n) result(connected_faces_count)
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    integer(I4B) :: connected_faces_count
    ! -- local
    integer(I4B) :: ipos, count, idx
    integer(I4B) :: m, isympos, ihc
    integer(I4B) :: number_connections
    real(DP), dimension(:, :), allocatable :: face_normals
    real(DP) :: x_dir, y_dir, z_dir
    logical :: normal_exists

    number_connections = number_connected_faces(dis, n)
    allocate (face_normals(3, number_connections))

    count = 0
    do ipos = dis%con%ia(n) + 1, dis%con%ia(n + 1) - 1
      m = dis%con%ja(ipos)
      isympos = dis%con%jas(ipos)
      ihc = dis%con%ihc(isympos)

      call dis%connection_normal(n, m, ihc, x_dir, y_dir, z_dir, ipos)

      normal_exists = .false.
      do idx = 1, count
        if (abs(x_dir - face_normals(1, idx)) < DSAME .and. &
            abs(y_dir - face_normals(2, idx)) < DSAME .and. &
            abs(z_dir - face_normals(3, idx)) < DSAME) then

          normal_exists = .true.
          exit
        end if
      end do

      if (.not. normal_exists) then
        count = count + 1
        face_normals(1, count) = x_dir
        face_normals(2, count) = y_dir
        face_normals(3, count) = z_dir
      end if
    end do

    connected_faces_count = count

  end function number_unique_connected_faces

  function number_local_boundary_faces(dis, n) result(boundary_faces_count)
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    integer(I4B) :: boundary_faces_count
    ! -- local
    integer(I4B) :: number_cell_faces, number_connections

    number_cell_faces = number_faces(dis, n)
    number_connections = number_unique_connected_faces(dis, n)

    boundary_faces_count = number_cell_faces - number_connections

  end function number_local_boundary_faces

  function number_global_boundary_faces(dis) result(count)
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B) :: count
    ! -- local
    integer(I4B) :: n

    count = 0
    do n = 1, dis%nodes
      count = count + number_boundary_faces(dis, n)
    end do

  end function number_global_boundary_faces

  !!> @brief Returns the centroid coordinates of a given cell.
  !!
  !! This function computes the centroid (geometric center) of cell `n` in the discretization.
  !! The centroid is returned as a 3-element vector containing the x, y, and z coordinates.
  !! The x and y coordinates are taken directly from the cell center arrays, while the z coordinate
  !! is computed as the average of the cell's top and bottom elevations.
  !<
  function cell_centroid(dis, n) result(centroid)
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    real(DP), dimension(3) :: centroid

    centroid(1) = dis%xc(n)
    centroid(2) = dis%yc(n)
    centroid(3) = (dis%top(n) + dis%bot(n)) / 2.0_dp
  end function cell_centroid

  !!> @brief Returns the vector distance between two connected nodes.
  !!
  !! This function computes the vector from the center of cell `n` to the center of cell `m`
  !! in the discretization. If the cells are directly connected, the vector is computed along
  !! the connection direction, taking into account cell geometry and, if available, cell saturation.
  !! If the cells are not directly connected, the vector is simply the difference between their centroids.
  !<
  function node_distance(dis, fmi, n, m) result(d)
    !-- modules
    use TspFmiModule, only: TspFmiType
    ! -- return
    real(DP), dimension(3) :: d
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    class(TspFmiType), pointer, intent(in) :: fmi
    integer(I4B), intent(in) :: n, m
    ! -- local
    real(DP) :: x_dir, y_dir, z_dir, length
    real(DP) :: satn, satm
    integer(I4B) :: ipos, isympos, ihc
    real(DP), dimension(3) :: xn, xm

    ! -- Find the connection position (isympos) between cell n and cell m
    isympos = -1
    do ipos = dis%con%ia(n) + 1, dis%con%ia(n + 1) - 1
      if (dis%con%ja(ipos) == m) then
        isympos = dis%con%jas(ipos)
        exit
      end if
    end do

    ! -- if the connection is not found, then return the distance between the two nodes
    ! -- Thi can happen when using an extended stencil (neighbours-of-neigbhours) to compute the gradients
    if (isympos == -1) then
      xn = cell_centroid(dis, n)
      xm = cell_centroid(dis, m)

      d = xm - xn
      return
    end if

    ! -- Account for the saturation levels if available
    ihc = dis%con%ihc(isympos)
    if (associated(fmi%gwfsat)) then
      satn = fmi%gwfsat(n)
      satm = fmi%gwfsat(m)
    else
      satn = DONE
      satm = DONE
    end if

    ! -- Get the connection direction and length
    call dis%connection_vector(n, m, .true., satn, satm, ihc, x_dir, &
                               y_dir, z_dir, length)

    ! -- Compute the distance vector
    d(1) = x_dir * length
    d(2) = y_dir * length
    d(3) = z_dir * length

  end function node_distance

end module DisInfoModule
