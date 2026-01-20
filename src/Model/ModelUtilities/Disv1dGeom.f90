module Disv1dGeom

  use KindModule, only: DP, I4B
  implicit none
  private
  public :: calcdist, line_unit_vector

contains

  !> @brief Calculate distance between two vertices
  !<
  function calcdist(vertices, ivert1, ivert2) result(dist)
    ! -- dummy
    real(DP), dimension(:, :), intent(in) :: vertices
    integer(I4B), intent(in) :: ivert1
    integer(I4B), intent(in) :: ivert2
    real(DP) :: dist, xdist, ydist
    ! -- local
    !
    ! -- calc distance
    xdist = abs(vertices(1, ivert1) - vertices(1, ivert2))
    ydist = abs(vertices(2, ivert1) - vertices(2, ivert2))
    dist = sqrt(xdist * xdist + ydist * ydist)
  end function calcdist

  !> @brief Calculate distance between two vertices
  !!
  !! Calculate the vector components (xcomp, ycomp, and zcomp)
  !! for a line defined by two points, (x0, y0, z0), (x1, y1, z1). Also
  !! return the magnitude of the original vector, vmag.
  !!
  !<
  subroutine line_unit_vector(x0, y0, x1, y1, &
                              xcomp, ycomp, vmag)
    real(DP), intent(in) :: x0
    real(DP), intent(in) :: y0
    real(DP), intent(in) :: x1
    real(DP), intent(in) :: y1
    real(DP), intent(out) :: xcomp
    real(DP), intent(out) :: ycomp
    real(DP) :: dx, dy, vmag
    !
    dx = x1 - x0
    dy = y1 - y0
    vmag = sqrt(dx**2 + dy**2)
    xcomp = dx / vmag
    ycomp = dy / vmag
  end subroutine line_unit_vector

end module Disv1dGeom
