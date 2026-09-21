module TestMawNonVertical

  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DTWO, DPIO180
  use MathUtilModule, only: is_close
  use testdrive, only: check, error_type, new_unittest, unittest_type
  use GeomUtilModule, only: polygon_extent
  use MawModule, only: maw_screen_length

  implicit none
  private
  public :: collect_mawnonvertical

  ! a cell 20 length units long and 10 length units thick, with a screen that
  ! spans the cell thickness, matching the multi-aquifer well non-vertical
  ! connection chapter of the supplemental technical information
  real(DP), parameter :: dx = 20.0_DP
  real(DP), parameter :: dz = 10.0_DP
  real(DP), parameter :: radius = 0.25_DP

contains

  subroutine collect_mawnonvertical(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("screen_length_vertical", test_screen_vertical), &
                new_unittest("screen_length_tilted", test_screen_tilted), &
                new_unittest("screen_length_grows_near_horizontal", &
                             test_screen_near_horizontal), &
                new_unittest("extent_dis_cell", test_extent_dis), &
                new_unittest("extent_disv_cell", test_extent_disv), &
                new_unittest("extent_disu_cell", test_extent_disu), &
                new_unittest("connection_within_cell", test_within_cell), &
                new_unittest("connection_beyond_cell", test_beyond_cell) &
                ]
  end subroutine collect_mawnonvertical

  !> @brief Vertices of a dis cell, which is always a rectangle.
  !<
  subroutine dis_cell(xv, yv, delr, delc)
    real(DP), dimension(4), intent(out) :: xv
    real(DP), dimension(4), intent(out) :: yv
    real(DP), intent(in) :: delr
    real(DP), intent(in) :: delc
    xv = [DZERO, DZERO, delr, delr]
    yv = [DZERO, delc, delc, DZERO]
  end subroutine dis_cell

  !> @brief Vertices of a disv cell, here a regular hexagon.
  !<
  subroutine disv_cell(xv, yv, r)
    real(DP), dimension(6), intent(out) :: xv
    real(DP), dimension(6), intent(out) :: yv
    real(DP), intent(in) :: r
    integer(I4B) :: i
    real(DP) :: a
    do i = 1, 6
      a = real(i - 1, DP) * 6.0e1_DP * DPIO180
      xv(i) = r * cos(a)
      yv(i) = r * sin(a)
    end do
  end subroutine disv_cell

  !> @brief The horizontal distance spanned by a connection.
  !<
  function horizontal_distance(angle) result(hlen)
    real(DP), intent(in) :: angle !< tilt angle from vertical, in degrees
    real(DP) :: hlen
    real(DP) :: omega
    omega = angle * DPIO180
    hlen = maw_screen_length(dz, radius, omega) * sin(omega)
  end function horizontal_distance

  !> @brief A vertical connection is as long as the screen is thick.
  !<
  subroutine test_screen_vertical(error)
    type(error_type), allocatable, intent(out) :: error
    call check(error, is_close(maw_screen_length(dz, radius, DZERO), dz))
  end subroutine test_screen_vertical

  !> @brief A tilted connection is longer than the screen is thick.
  !<
  subroutine test_screen_tilted(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: lw
    ! (10 - 2 * 0.25 * sin(60)) / cos(60)
    lw = maw_screen_length(dz, radius, 6.0e1_DP * DPIO180)
    call check(error, is_close(lw, 19.13397459621556_DP))
  end subroutine test_screen_tilted

  !> @brief The connection length grows without bound near horizontal.
  !<
  subroutine test_screen_near_horizontal(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: lw
    ! (10 - 2 * 0.25 * sin(85)) / cos(85), more than five times the cell length
    lw = maw_screen_length(dz, radius, 8.5e1_DP * DPIO180)
    call check(error, is_close(lw, 109.02210630531792_DP))
  end subroutine test_screen_near_horizontal

  !> @brief The extent of a dis cell is the diagonal of the rectangle.
  !<
  subroutine test_extent_dis(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(4) :: xv, yv
    call dis_cell(xv, yv, dx, 0.5_DP * dx)
    ! sqrt(20^2 + 10^2)
    call check(error, is_close(polygon_extent(xv, yv), 22.360679774997898_DP))
  end subroutine test_extent_dis

  !> @brief The extent of a disv cell is the longest diagonal of the polygon.
  !<
  subroutine test_extent_disv(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(6) :: xv, yv
    ! the longest diagonal of a regular hexagon is twice the circumradius
    call disv_cell(xv, yv, 0.5_DP * dx)
    call check(error, is_close(polygon_extent(xv, yv), dx))
  end subroutine test_extent_disv

  !> @brief A disu cell has no vertices, so its extent comes from the area.
  !<
  subroutine test_extent_disu(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(4) :: xv, yv
    real(DP) :: from_area
    ! the extent taken from the area of a square cell must be the extent
    ! calculated from the vertices of the same cell
    call dis_cell(xv, yv, dx, dx)
    from_area = sqrt(DTWO * dx * dx)
    call check(error, is_close(polygon_extent(xv, yv), from_area))
  end subroutine test_extent_disu

  !> @brief A 60 degree connection stays within a square cell.
  !<
  subroutine test_within_cell(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(4) :: xv, yv
    call dis_cell(xv, yv, dx, dx)
    call check(error, horizontal_distance(6.0e1_DP) < polygon_extent(xv, yv))
  end subroutine test_within_cell

  !> @brief An 85 degree connection extends well beyond a square cell.
  !<
  subroutine test_beyond_cell(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(4) :: xv, yv
    real(DP) :: hlen
    call dis_cell(xv, yv, dx, dx)
    hlen = horizontal_distance(8.5e1_DP)
    call check(error, hlen > polygon_extent(xv, yv))
    ! the connection spans more than three cell diagonals
    call check(error, hlen > 3.0_DP * polygon_extent(xv, yv))
  end subroutine test_beyond_cell

end module TestMawNonVertical
