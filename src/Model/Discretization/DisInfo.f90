module DisInfoModule
  use KindModule, only: DP, I4B
  use BaseDisModule, only: DisBaseType

  implicit none
  private

  interface number_boundary_faces
    module procedure number_local_boundary_faces, number_global_boundary_faces
  end interface

  public :: number_connected_faces
  public :: number_faces
  public :: number_boundary_faces
  public :: cell_centroid
  private :: number_local_boundary_faces
  private :: number_global_boundary_faces

contains

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

  function number_connected_faces(dis, n) result(connected_faces_count)
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    integer(I4B) :: connected_faces_count

    connected_faces_count = dis%con%ia(n + 1) - dis%con%ia(n) - 1
  end function number_connected_faces

  function number_local_boundary_faces(dis, n) result(boundary_faces_count)
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    integer(I4B) :: boundary_faces_count
    ! -- local
    integer(I4B) :: number_sides, number_connections

    number_sides = number_faces(dis, n)
    number_connections = number_connected_faces(dis, n)

    boundary_faces_count = number_sides - number_connections

  end function number_local_boundary_faces

  function number_global_boundary_faces(dis) result(count)
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B) :: count
    ! -- local
    integer(I4B) :: nodes, n

    count = 0
    nodes = dis%nodes
    do n = 1, nodes
      count = count + number_boundary_faces(dis, n)
    end do

  end function number_global_boundary_faces

  function cell_centroid(dis, n) result(centroid)
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    real(DP), dimension(3) :: centroid

    centroid(1) = dis%xc(n)
    centroid(2) = dis%yc(n)
    centroid(3) = (dis%top(n) + dis%bot(n)) / 2.0_dp
  end function cell_centroid

end module DisInfoModule
