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
      type(GlobalCellType), pointer :: gcn, gcm
      real(DP) :: xn, yn, zn, xm, ym, zm
      real(DP), dimension(:), pointer, contiguous :: dis_top_n, dis_bot_n, &
                                                     dis_top_m, dis_bot_m
      real(DP), dimension(:), pointer, contiguous :: dis_xc_n, dis_yc_n, &
                                                     dis_xc_m, dis_yc_m
      real(DP) :: xorigin_n, yorigin_n, angrot_n, &
                  xorigin_m, yorigin_m, angrot_m

      ! get coordinates
      gcn => idxToGlobal(array(n))
      gcm => idxToGlobal(array(m))

      ! for n:
      if (.not. z_only) then
        dis_xc_n => gcn%v_model%dis_xc%get_array()
        dis_yc_n => gcn%v_model%dis_yc%get_array()
        xorigin_n = gcn%v_model%dis_xorigin%get()
        yorigin_n = gcn%v_model%dis_yorigin%get()
        angrot_n = gcn%v_model%dis_angrot%get()
        ! convert coordinates
        call dis_transform_xy(dis_xc_n(gcn%index), dis_yc_n(gcn%index), &
                              xorigin_n, yorigin_n, angrot_n, &
                              xn, yn)
      else
        xn = DZERO
        yn = DZERO
      end if
      dis_top_n => gcn%v_model%dis_top%get_array()
      dis_bot_n => gcn%v_model%dis_bot%get_array()
      zn = DHALF * (dis_top_n(gcn%index) + dis_bot_n(gcn%index))

      ! for m:
      if (.not. z_only) then
        dis_xc_m => gcm%v_model%dis_xc%get_array()
        dis_yc_m => gcm%v_model%dis_yc%get_array()
        xorigin_m = gcm%v_model%dis_xorigin%get()
        yorigin_m = gcm%v_model%dis_yorigin%get()
        angrot_m = gcm%v_model%dis_angrot%get()
        call dis_transform_xy(dis_xc_m(gcm%index), dis_yc_m(gcm%index), &
                              xorigin_m, yorigin_m, angrot_m, &
                              xm, ym)
      else
        xm = DZERO
        ym = DZERO
      end if
      dis_top_m => gcm%v_model%dis_top%get_array()
      dis_bot_m => gcm%v_model%dis_bot%get_array()
      zm = DHALF * (dis_top_m(gcm%index) + dis_bot_m(gcm%index))

      ! compare
      if (.not. is_close(zn, zm, 10 * epsilon(zn))) then
        isLess = zn > zm
      else if (.not. is_close(yn, ym, 10 * epsilon(yn))) then
        isLess = yn > ym
      else if (.not. is_close(xn, xm, 10 * epsilon(xn))) then
        isLess = xn < xm
      else
        isLess = .false.
      end if

    end function lessThan

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
