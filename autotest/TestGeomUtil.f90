module TestGeomUtil
  use KindModule, only: I4B, DP
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       to_string, unittest_type
  use GeomUtilModule, only: get_node, get_ijk, get_jk, point_in_polygon, &
                            skew, area, polygon_extent, shared_face
  use ConstantsModule, only: LINELENGTH, DEM7, DPIO180
  use MathUtilModule, only: is_close
  use DisvGeom, only: shared_edge
  implicit none
  private
  public :: collect_geomutil

contains

  subroutine collect_geomutil(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("get_node_get_ijk", &
                             test_get_node_get_ijk), &
                new_unittest("point_in_polygon_sq", &
                             test_point_in_polygon_sq), &
                new_unittest("point_in_polygon_tri", &
                             test_point_in_polygon_tri), &
                new_unittest("point_in_polygon_irr", &
                             test_point_in_polygon_irr), &
                new_unittest("point_in_polygon_tol", &
                             test_point_in_polygon_tol), &
                new_unittest("skew", test_skew), &
                new_unittest("area", test_area), &
                new_unittest("polygon_extent", test_polygon_extent), &
                new_unittest("polygon_extent_degenerate", &
                             test_polygon_extent_degenerate), &
                new_unittest("shared_face", test_shared_face), &
                new_unittest("shared_face_large", &
                             test_shared_face_large) &
                ]
  end subroutine collect_geomutil

  ! 2D arrays for polygons and check points use column-major indexing

  subroutine test_get_node_get_ijk(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: ilay
    integer :: irow
    integer :: icol
    integer :: nlay
    integer :: nrow
    integer :: ncol
    integer :: nnum
    integer :: ncls
    integer :: k, i, j

    ! trivial grid with 1 cell
    nnum = get_node(1, 1, 1, 1, 1, 1)
    call get_ijk(nnum, 1, 1, 1, ilay, irow, icol)
    call check(error, nnum == 1)
    call check(error, ilay == 1)
    call check(error, irow == 1)
    call check(error, icol == 1)
    if (allocated(error)) return

    ! small grid, 3x4x5
    nlay = 3
    nrow = 4
    ncol = 5
    ncls = nlay * nrow * ncol
    do k = 1, nlay
      do i = 1, nrow
        do j = 1, ncol
          ! node number from ijk
          nnum = get_node(k, i, j, nlay, nrow, ncol)
          call check(error, nnum == (k - 1) * nrow * ncol + (i - 1) * ncol + j)
          if (allocated(error)) return

          ! ijk from node number
          call get_ijk(nnum, nrow, ncol, nlay, irow, icol, ilay)
          call check(error, ilay == k)
          call check(error, irow == i)
          call check(error, icol == j)
          if (allocated(error)) return
        end do
      end do
    end do
  end subroutine test_get_node_get_ijk

  subroutine test_point_in_polygon(error, shape, &
                                   poly, in_pts, out_pts, vert_pts, face_pts)
    type(error_type), allocatable, intent(inout) :: error
    character(len=*), intent(in) :: shape
    real(DP), allocatable, intent(in) :: poly(:, :)
    real(DP), allocatable, intent(in) :: in_pts(:, :)
    real(DP), allocatable, intent(in) :: out_pts(:, :)
    real(DP), allocatable, intent(in) :: vert_pts(:, :)
    real(DP), allocatable, intent(in) :: face_pts(:, :)
    integer(I4B) :: i
    real(DP) :: x, y

    ! test inside points
    do i = 1, size(in_pts, 2)
      x = in_pts(1, i)
      y = in_pts(2, i)
      call check(error, point_in_polygon(x, y, poly), &
                 "point inside "//shape//" failed: " &
                 //to_string(x)//", "//to_string(y))
      if (allocated(error)) return
    end do

    ! test outside points
    do i = 1, size(out_pts, 2)
      x = out_pts(1, i)
      y = out_pts(2, i)
      call check(error, (.not. point_in_polygon(x, y, poly)), &
                 "point outside "//shape//" failed: " &
                 //to_string(x)//", "//to_string(y))
      if (allocated(error)) return
    end do

    ! test vertex points
    do i = 1, size(vert_pts, 2)
      x = vert_pts(1, i)
      y = vert_pts(2, i)
      call check(error, point_in_polygon(x, y, poly), &
                 "point on "//shape//" vertex failed: " &
                 //to_string(x)//", "//to_string(y))
      if (allocated(error)) return
    end do

    ! test face points
    do i = 1, size(face_pts, 2)
      x = face_pts(1, i)
      y = face_pts(2, i)
      call check(error, point_in_polygon(x, y, poly), &
                 "point on "//shape//" face failed: " &
                 //to_string(x)//", "//to_string(y))
      if (allocated(error)) return
    end do
  end subroutine test_point_in_polygon

  !> @brief Test a unit square
  subroutine test_point_in_polygon_sq(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), allocatable :: poly(:, :)
    real(DP), allocatable :: in_pts(:, :)
    real(DP), allocatable :: out_pts(:, :)
    real(DP), allocatable :: vert_pts(:, :)
    real(DP), allocatable :: face_pts(:, :)

    allocate (poly(2, 4))

    allocate (in_pts(2, 3))
    in_pts(:, 1) = (/0.99_DP, 0.01_DP/)
    in_pts(:, 2) = (/0.5_DP, 0.5_DP/)
    in_pts(:, 3) = (/0.0001_DP, 0.9999_DP/)

    allocate (out_pts(2, 2))
    out_pts(:, 1) = (/0.5_DP, 1.00001_DP/)
    out_pts(:, 2) = (/-0.5_DP, 34.0_DP/)

    allocate (vert_pts(2, 4))
    vert_pts(:, 1) = (/0.0_DP, 0.0_DP/)
    vert_pts(:, 2) = (/1.0_DP, 0.0_DP/)
    vert_pts(:, 3) = (/0.0_DP, 1.0_DP/)
    vert_pts(:, 4) = (/1.0_DP, 1.0_DP/)

    allocate (face_pts(2, 4))
    face_pts(:, 1) = (/0.0_DP, 0.5_DP/)
    face_pts(:, 2) = (/0.5_DP, 0.0_DP/)
    face_pts(:, 3) = (/1.0_DP, 0.5_DP/)
    face_pts(:, 4) = (/0.5_DP, 1.0_DP/)

    poly(:, 1) = (/0.0_DP, 0.0_DP/)
    poly(:, 2) = (/0.0_DP, 1.0_DP/)
    poly(:, 3) = (/1.0_DP, 1.0_DP/)
    poly(:, 4) = (/1.0_DP, 0.0_DP/)
    call test_point_in_polygon(error, "clockwise square", &
                               poly, in_pts, out_pts, vert_pts, face_pts)
    if (allocated(error)) return

    poly(:, 1) = (/0.0_DP, 0.0_DP/)
    poly(:, 2) = (/1.0_DP, 0.0_DP/)
    poly(:, 3) = (/1.0_DP, 1.0_DP/)
    poly(:, 4) = (/0.0_DP, 1.0_DP/)
    call test_point_in_polygon(error, "counter-clockwise square", &
                               poly, in_pts, out_pts, vert_pts, face_pts)
    if (allocated(error)) return

    deallocate (poly)
    deallocate (in_pts)
    deallocate (out_pts)
    deallocate (vert_pts)
    deallocate (face_pts)
  end subroutine test_point_in_polygon_sq

  !> @brief Test a right triangle
  subroutine test_point_in_polygon_tri(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), allocatable :: poly(:, :)
    real(DP), allocatable :: in_pts(:, :)
    real(DP), allocatable :: out_pts(:, :)
    real(DP), allocatable :: vert_pts(:, :)
    real(DP), allocatable :: face_pts(:, :)

    allocate (poly(2, 3))

    allocate (in_pts(2, 3))
    in_pts(:, 1) = (/0.8_DP, 0.0001_DP/)
    in_pts(:, 2) = (/0.5_DP, 0.49999_DP/)
    in_pts(:, 3) = (/0.0001_DP, 0.8_DP/)

    allocate (out_pts(2, 2))
    out_pts(:, 1) = (/0.5_DP, 0.50001_DP/)
    out_pts(:, 2) = (/-0.5_DP, 34.0_DP/)

    allocate (vert_pts(2, 3))
    vert_pts(:, 1) = (/0.0_DP, 0.0_DP/)
    vert_pts(:, 2) = (/1.0_DP, 0.0_DP/)
    vert_pts(:, 3) = (/0.0_DP, 1.0_DP/)

    allocate (face_pts(2, 3))
    face_pts(:, 1) = (/0.0_DP, 0.5_DP/)
    face_pts(:, 2) = (/0.5_DP, 0.0_DP/)
    face_pts(:, 3) = (/0.5_DP, 0.5_DP/)

    poly(:, 1) = (/0.0_DP, 0.0_DP/)
    poly(:, 2) = (/0.0_DP, 1.0_DP/)
    poly(:, 3) = (/1.0_DP, 0.0_DP/)
    call test_point_in_polygon(error, "clockwise triangle", &
                               poly, in_pts, out_pts, vert_pts, face_pts)
    if (allocated(error)) return

    poly(:, 1) = (/0.0_DP, 0.0_DP/)
    poly(:, 2) = (/1.0_DP, 0.0_DP/)
    poly(:, 3) = (/0.0_DP, 1.0_DP/)
    call test_point_in_polygon(error, "counter-clockwise triangle", &
                               poly, in_pts, out_pts, vert_pts, face_pts)
    if (allocated(error)) return

    deallocate (poly)
    deallocate (in_pts)
    deallocate (out_pts)
    deallocate (vert_pts)
    deallocate (face_pts)
  end subroutine test_point_in_polygon_tri

  !> @brief Test an irregular polygon
  subroutine test_point_in_polygon_irr(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), allocatable :: poly(:, :)
    real(DP), allocatable :: in_pts(:, :)
    real(DP), allocatable :: out_pts(:, :)
    real(DP), allocatable :: vert_pts(:, :)
    real(DP), allocatable :: face_pts(:, :)

    allocate (poly(2, 5))

    allocate (in_pts(2, 3))
    in_pts(:, 1) = (/0.5_DP, 0.1_DP/)
    in_pts(:, 2) = (/0.5_DP, 0.49_DP/)
    in_pts(:, 3) = (/1.999_DP, 1.999_DP/)

    allocate (out_pts(2, 3))
    out_pts(:, 1) = (/0.5_DP, -0.1_DP/)
    out_pts(:, 2) = (/0.5_DP, 0.51_DP/)
    out_pts(:, 3) = (/-0.5_DP, 34.0_DP/)

    allocate (vert_pts(2, 5))
    vert_pts(:, 1) = (/0.0_DP, 0.0_DP/)
    vert_pts(:, 2) = (/1.0_DP, 1.0_DP/)
    vert_pts(:, 3) = (/1.0_DP, 2.0_DP/)
    vert_pts(:, 4) = (/2.0_DP, 2.0_DP/)
    vert_pts(:, 5) = (/2.0_DP, 0.0_DP/)

    allocate (face_pts(2, 3))
    face_pts(:, 1) = (/0.5_DP, 0.5_DP/)
    face_pts(:, 2) = (/2.0_DP, 1.0_DP/)
    face_pts(:, 3) = (/1.5_DP, 2.0_DP/)

    poly(:, 1) = (/0.0_DP, 0.0_DP/)
    poly(:, 2) = (/1.0_DP, 1.0_DP/)
    poly(:, 3) = (/1.0_DP, 2.0_DP/)
    poly(:, 4) = (/2.0_DP, 2.0_DP/)
    poly(:, 5) = (/2.0_DP, 0.0_DP/)
    call test_point_in_polygon(error, &
                               "clockwise irregular polygon", &
                               poly, in_pts, out_pts, vert_pts, face_pts)
    if (allocated(error)) return

    poly(:, 1) = (/0.0_DP, 0.0_DP/)
    poly(:, 2) = (/2.0_DP, 0.0_DP/)
    poly(:, 3) = (/2.0_DP, 2.0_DP/)
    poly(:, 4) = (/1.0_DP, 2.0_DP/)
    poly(:, 5) = (/1.0_DP, 1.0_DP/)
    call test_point_in_polygon(error, &
                               "counter-clockwise irregular polygon", &
                               poly, in_pts, out_pts, vert_pts, face_pts)
    if (allocated(error)) return

    deallocate (poly)
    deallocate (in_pts)
    deallocate (out_pts)
    deallocate (vert_pts)
    deallocate (face_pts)
  end subroutine test_point_in_polygon_irr

  !> @brief Test the optional edge tolerance.
  subroutine test_point_in_polygon_tol(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), allocatable :: poly(:, :)
    real(DP) :: xll, yll, dx, dy
    real(DP) :: xoff, yoff, cellsize, tol

    xll = 512345.678901_DP
    yll = 4567890.123456_DP
    dx = 10.0_DP
    dy = 10.0_DP

    allocate (poly(2, 4))
    poly(:, 1) = (/xll, yll/)
    poly(:, 2) = (/xll, yll + dy/)
    poly(:, 3) = (/xll + dx, yll + dy/)
    poly(:, 4) = (/xll + dx, yll/)

    ! points exactly on a vertex or edge are always in the polygon,
    ! tolerance or not
    call check(error, point_in_polygon(xll, yll, poly), &
               "point on vertex failed")
    if (allocated(error)) return
    call check(error, point_in_polygon(xll + dx, yll + 5.0_DP, poly), &
               "point on edge failed")
    if (allocated(error)) return

    ! a point slightly off the edge
    xoff = xll + dx + 4.0_DP * spacing(xll + dx)
    yoff = yll + 5.0_DP

    ! without a tolerance, the near-edge point is (correctly) rejected
    call check(error,.not. point_in_polygon(xoff, yoff, poly), &
               "near-edge point wrongly accepted without tolerance")
    if (allocated(error)) return

    ! a tolerance scaled to the cell accepts the same near-edge point
    cellsize = max(maxval(poly(1, :)) - minval(poly(1, :)), &
                   maxval(poly(2, :)) - minval(poly(2, :)))
    tol = cellsize * cellsize * DEM7
    call check(error, point_in_polygon(xoff, yoff, poly, tol), &
               "near-edge point wrongly rejected with tolerance")
    if (allocated(error)) return

    ! a point clearly outside the cell is still rejected, even with
    ! the tolerance
    call check(error, &
               .not. point_in_polygon(xll + dx + 0.001_DP, yoff, &
                                      poly, tol), &
               "clearly-outside point wrongly accepted with tolerance")
    if (allocated(error)) return

    ! the same checks, repeated for a point near a horizontal edge
    ! (top edge, ya == yb) rather than a vertical one, since the
    ! on-edge check is derived separately for horizontal edges
    xoff = xll + 5.0_DP
    yoff = yll + dy + 4.0_DP * spacing(yll + dy)

    call check(error,.not. point_in_polygon(xoff, yoff, poly), &
               "near-horizontal-edge point wrongly accepted " &
               //"without tolerance")
    if (allocated(error)) return

    call check(error, point_in_polygon(xoff, yoff, poly, tol), &
               "near-horizontal-edge point wrongly rejected " &
               //"with tolerance")
    if (allocated(error)) return

    call check(error, &
               .not. point_in_polygon(xoff, yll + dy + 0.001_DP, &
                                      poly, tol), &
               "clearly-outside point wrongly accepted with tolerance " &
               //"(horizontal edge)")
    if (allocated(error)) return

    deallocate (poly)
  end subroutine test_point_in_polygon_tol

  subroutine test_skew(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: v(2)

    ! shear to right
    v = (/1.0_DP, 1.0_DP/)
    v = skew(v, (/1.0_DP, 1.0_DP, 1.0_DP/))
    call check(error, v(1) == 2.0_DP .and. v(2) == 1.0_DP)
    if (allocated(error)) return

    v = (/2.0_DP, 2.0_DP/)
    v = skew(v, (/1.0_DP, 0.5_DP, 1.0_DP/))
    call check(error, v(1) == 3.0_DP .and. v(2) == 2.0_DP)
    if (allocated(error)) return

    ! collapse x dim
    v = (/2.0_DP, 2.0_DP/)
    v = skew(v, (/0.0_DP, 0.5_DP, 1.0_DP/))
    call check(error, v(1) == 1.0_DP .and. v(2) == 2.0_DP, to_string(v(1)))
    if (allocated(error)) return

    ! mirror over x axis
    v = (/2.0_DP, 2.0_DP/)
    v = skew(v, (/-1.0_DP, 0.0_DP, 1.0_DP/))
    call check(error, v(1) == -2.0_DP .and. v(2) == 2.0_DP, to_string(v(1)))
    if (allocated(error)) return

    ! mirror over x and y axis
    v = (/2.0_DP, 2.0_DP/)
    v = skew(v, (/-1.0_DP, 0.0_DP, -1.0_DP/))
    call check(error, v(1) == -2.0_DP .and. v(2) == -2.0_DP, to_string(v(1)))
    if (allocated(error)) return
  end subroutine test_skew

  subroutine test_area(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), allocatable :: poly(:, :)
    real(DP) :: a

    allocate (poly(2, 5))

    poly(:, 1) = (/0.0_DP, 0.0_DP/)
    poly(:, 2) = (/1.0_DP, 1.0_DP/)
    poly(:, 3) = (/1.0_DP, 2.0_DP/)
    poly(:, 4) = (/2.0_DP, 2.0_DP/)
    poly(:, 5) = (/2.0_DP, 0.0_DP/)

    a = area(poly(1, :), poly(2, :))
    call check(error, a == 2.5_DP, to_string(a))
    if (allocated(error)) return

    poly(:, 5) = (/0.0_DP, 0.0_DP/)
    poly(:, 4) = (/1.0_DP, 1.0_DP/)
    poly(:, 3) = (/1.0_DP, 2.0_DP/)
    poly(:, 2) = (/2.0_DP, 2.0_DP/)
    poly(:, 1) = (/2.0_DP, 0.0_DP/)

    a = area(poly(1, :), poly(2, :), cw=.false.)
    call check(error, a == 2.5_DP, to_string(a))
    if (allocated(error)) return

    deallocate (poly)

  end subroutine

  subroutine test_polygon_extent(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), allocatable :: poly(:, :)
    real(DP) :: e
    integer(I4B) :: i

    ! a rectangle, as a dis cell always is: the extent is the diagonal
    allocate (poly(2, 4))
    poly(:, 1) = (/0.0_DP, 0.0_DP/)
    poly(:, 2) = (/0.0_DP, 10.0_DP/)
    poly(:, 3) = (/20.0_DP, 10.0_DP/)
    poly(:, 4) = (/20.0_DP, 0.0_DP/)

    e = polygon_extent(poly(1, :), poly(2, :))
    call check(error, is_close(e, 22.360679774997898_DP), to_string(e))
    if (allocated(error)) return

    ! the vertices may be given in any order
    poly = cshift(poly, shift=2, dim=2)
    e = polygon_extent(poly(1, :), poly(2, :))
    call check(error, is_close(e, 22.360679774997898_DP), to_string(e))
    if (allocated(error)) return
    deallocate (poly)

    ! a regular hexagon, one shape a disv cell takes: the extent is the
    ! longest diagonal, which is twice the circumradius
    allocate (poly(2, 6))
    do i = 1, 6
      poly(1, i) = 10.0_DP * cos(real(i - 1, DP) * DPIO180 * 60.0_DP)
      poly(2, i) = 10.0_DP * sin(real(i - 1, DP) * DPIO180 * 60.0_DP)
    end do

    e = polygon_extent(poly(1, :), poly(2, :))
    call check(error, is_close(e, 20.0_DP), to_string(e))
    if (allocated(error)) return
    deallocate (poly)

    ! a triangle: the extent is the longest side
    allocate (poly(2, 3))
    poly(:, 1) = (/0.0_DP, 0.0_DP/)
    poly(:, 2) = (/3.0_DP, 0.0_DP/)
    poly(:, 3) = (/0.0_DP, 4.0_DP/)

    e = polygon_extent(poly(1, :), poly(2, :))
    call check(error, is_close(e, 5.0_DP), to_string(e))
    if (allocated(error)) return
    deallocate (poly)

  end subroutine

  subroutine test_polygon_extent_degenerate(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: e

    ! a single vertex has no extent
    e = polygon_extent((/1.0_DP/), (/2.0_DP/))
    call check(error, is_close(e, 0.0_DP), to_string(e))
    if (allocated(error)) return

    ! two vertices are a distance apart
    e = polygon_extent((/0.0_DP, 3.0_DP/), (/0.0_DP, 4.0_DP/))
    call check(error, is_close(e, 5.0_DP), to_string(e))
    if (allocated(error)) return

  end subroutine

  subroutine test_shared_face(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: iverts1(5), iverts2(5)
    integer(I4B) :: iface

    iface = 0
    iverts1 = (/1, 2, 5, 4, 1/)
    iverts2 = (/2, 3, 6, 5, 2/)

    call shared_face(iverts1, iverts2, iface)
    call check(error, iface == 2)
    if (allocated(error)) return

  end subroutine

  function rand_int(a, b) result(i)
    integer(I4B) :: a, b, i
    real(DP) :: u

    call random_number(u)
    i = a + floor((b + 1 - a) * u)

  end function

  subroutine test_shared_face_large(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: iverts1(33), iverts2(33)
    integer(I4B) :: iface, i, iv1, iv2
    real(DP) :: tstart, tstop

    iface = 0
    iverts1 = (/1, 2, 3, 4, 10, 14, 18, &
                22, 33, 36, 40, 44, 48, &
                52, 58, 57, 56, 55, 51, &
                47, 43, 39, 35, 32, 30, &
                28, 26, 24, 21, 17, 13, &
                9, 1/)
    iverts2 = (/5, 6, 7, 8, 12, 16, 20, &
                23, 25, 27, 29, 31, 34, &
                38, 42, 46, 50, 54, 62, &
                61, 60, 59, 53, 49, 45, &
                41, 37, 33, 22, 19, 15, &
                11, 5/)

    call shared_face(iverts1, iverts2, iface)
    call check(error, iface == 8, to_string(iface))
    if (allocated(error)) return

    call cpu_time(tstart)
    do i = 1, 1000
      iverts1 = cshift(iverts1, shift=rand_int(0, 32))
      call shared_face(iverts1, iverts2, iface)
    end do
    call cpu_time(tstop)
    print *, 'shared_face took: ', tstop - tstart

    call cpu_time(tstart)
    do i = 1, 1000
      iverts1 = cshift(iverts1, shift=rand_int(0, 32))
      call shared_edge(iverts1, iverts2, iv1, iv2)
    end do
    call cpu_time(tstop)
    print *, 'shared_edge took: ', tstop - tstart

  end subroutine

end module TestGeomUtil
