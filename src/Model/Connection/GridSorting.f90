module GridSorting
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: DHALF, DZERO
  use CellWithNbrsModule, only: GlobalCellType
  use MathUtilModule, only: is_close
  use BaseDisModule, only: dis_transform_xy
  implicit none
  private

  public :: quickSortGrid

contains
  ! Sort an array of integers
  subroutine quickSortGrid(array, arraySize, idxToGlobal, z_only)
    integer, intent(inout), dimension(:) :: array
    integer, intent(in) :: arraySize
    type(GlobalCellType), dimension(:), pointer :: idxToGlobal
    logical(LGP) :: z_only !< only sort by z coordinate
    ! local
    integer :: QSORT_THRESHOLD = 8
    include "qsort_inline.inc"

  contains
    subroutine init()
    end subroutine init

    ! Compare two grid cells, this doesn't work as
    ! smooth for staggered discretizations though...
    function lessThan(n, m) result(isLess)
      integer(I4B), intent(in) :: n
      integer(I4B), intent(in) :: m
      logical(LGP) :: isLess
      ! local
      real(DP), dimension(3) :: xyz_n, xyz_m

      ! get coordinates as 3-vectors
      xyz_n = get_global_xyz(idxToGlobal(array(n)), z_only)
      xyz_m = get_global_xyz(idxToGlobal(array(m)), z_only)

      ! compare
      if (.not. is_close(xyz_n(3), xyz_m(3), 10 * epsilon(xyz_n(3)))) then
        isLess = xyz_n(3) > xyz_m(3)
      else if (.not. is_close(xyz_n(2), xyz_m(2), 10 * epsilon(xyz_n(2)))) then
        isLess = xyz_n(2) > xyz_m(2)
      else if (.not. is_close(xyz_n(1), xyz_m(1), 10 * epsilon(xyz_n(1)))) then
        isLess = xyz_n(1) < xyz_m(1)
      else
        isLess = .false.
      end if

    end function lessThan

    !> @brief Utility function to convert global cell
    !< id to global x,y,z coordinates
    function get_global_xyz(gc, use_only_z) result(global_xyz)
      type(GlobalCellType) :: gc !< the global cell id
      logical(LGP) :: use_only_z !< only z coordinate is needed or available, skip transform
      real(DP), dimension(3) :: global_xyz !< return xyz
      ! local
      real(DP) :: x, y, z
      real(DP) :: xc, yc, xo, yo, angrot

      z = DHALF * (gc%v_model%dis_top%get(gc%index) + &
                   gc%v_model%dis_bot%get(gc%index))

      x = DZERO
      y = DZERO
      if (.not. use_only_z) then
        xc = gc%v_model%dis_xc%get(gc%index)
        yc = gc%v_model%dis_yc%get(gc%index)
        xo = gc%v_model%dis_xorigin%get()
        yo = gc%v_model%dis_yorigin%get()
        angrot = gc%v_model%dis_angrot%get()
        call dis_transform_xy(xc, yc, xo, yo, angrot, x, y)
      end if

      global_xyz = [x, y, z]

    end function get_global_xyz

    ! swap indices
    subroutine swap(a, b)
      integer, intent(in) :: a, b
      integer :: hold

      hold = array(a)
      array(a) = array(b)
      array(b) = hold

    end subroutine swap

    ! circular shift-right by one
    subroutine rshift(left, right)
      integer, intent(in) :: left, right
      integer :: hold

      hold = array(right)
      array(left + 1:right) = array(left:right - 1)
      array(left) = hold

    end subroutine rshift
  end subroutine quickSortGrid
end module GridSorting
