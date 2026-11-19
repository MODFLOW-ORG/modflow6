module TestSort
  use KindModule, only: I4B
  use testdrive, only: check, error_type, new_unittest, unittest_type
  use SortModule, only: qsort, search_sorted
  implicit none
  private
  public :: collect_sort

contains

  subroutine collect_sort(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("search_sorted", test_search_sorted), &
                new_unittest("search_sorted_empty", test_search_sorted_empty) &
                ]
  end subroutine collect_sort

  subroutine test_search_sorted(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), dimension(5) :: vals, svals, sidx
    integer(I4B) :: i

    vals = [40, 10, 50, 20, 30]
    do i = 1, 5
      svals(i) = vals(i)
      sidx(i) = i
    end do
    call qsort(sidx, svals)
    ! -- every value is found at its original position
    do i = 1, 5
      call check(error, search_sorted(svals, sidx, vals(i)) == i, &
                 "value should be found at its original position")
      if (allocated(error)) return
    end do
    ! -- values not in the array, below, between, and above
    call check(error, search_sorted(svals, sidx, 5) == 0, &
               "value below the range should not be found")
    if (allocated(error)) return
    call check(error, search_sorted(svals, sidx, 25) == 0, &
               "value between entries should not be found")
    if (allocated(error)) return
    call check(error, search_sorted(svals, sidx, 60) == 0, &
               "value above the range should not be found")
  end subroutine test_search_sorted

  subroutine test_search_sorted_empty(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), dimension(0) :: svals, sidx

    call check(error, search_sorted(svals, sidx, 1) == 0, &
               "nothing should be found in an empty array")
  end subroutine test_search_sorted_empty

end module TestSort
