module DisUtilsModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DSAME, DONE
  use BaseDisModule, only: DisBaseType

  implicit none
  private

  public :: number_connected_faces
  public :: cell_center
  public :: node_distance

contains

  !> @brief Returns the number of connected faces for a given cell.
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

  !> @brief Returns the vector distance from cell n to cell m.
  !!
  !! This function computes the vector from the center of cell `n` to the center of cell `m`
  !! in the discretization. If the cells are directly connected, the vector is computed along
  !! the connection direction, taking into account cell geometry and, if available, cell saturation.
  !! If the cells are not directly connected (e.g., when using an extended stencil such as neighbors-of-neighbors),
  !! the vector is simply the difference between their centroids: `d = centroid(m) - centroid(n)`.
  !! The returned vector always points from cell `n` to cell `m`.
  !<
  function node_distance(dis, n, m) result(d)
    !-- modules
    use TspFmiModule, only: TspFmiType
    ! -- return
    real(DP), dimension(3) :: d
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n, m
    ! -- local
    real(DP) :: x_dir, y_dir, z_dir, length
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
    ! -- This can happen when using an extended stencil (neighbours-of-neigbhours) to compute the gradients
    if (isympos == -1) then
      xn = cell_center(dis, n)
      xm = cell_center(dis, m)

      d = xm - xn
      return
    end if

    ! -- Get the connection direction and length
    ihc = dis%con%ihc(isympos)
    call dis%connection_vector(n, m, .false., DONE, DONE, ihc, x_dir, &
                               y_dir, z_dir, length)

    ! -- Compute the distance vector
    d(1) = x_dir * length
    d(2) = y_dir * length
    d(3) = z_dir * length

  end function node_distance

  !> @brief Returns the center coordinates of a given cell.
  !!
  !! This function computes the center of cell `n` in the discretization.
  !! The center is returned as a 3-element vector containing the x, y, and z coordinates.
  !! The x and y coordinates are taken directly from the cell center arrays, while the z coordinate
  !! is computed as the average of the cell's top and bottom elevations.
  !<
  function cell_center(dis, n) result(center)
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    real(DP), dimension(3) :: center

    center(1) = dis%xc(n)
    center(2) = dis%yc(n)
    center(3) = (dis%top(n) + dis%bot(n)) / 2.0_dp
  end function cell_center
end module DisUtilsModule
