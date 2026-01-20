module TestArrayHandlers
  use KindModule, only: I4B, DP, LGP
  use testdrive, only: error_type, unittest_type, new_unittest, check, &
                       test_failed, to_string
  use ArrayHandlersModule, only: ExpandArray, ExpandArray2D, ExtendPtrArray, &
                                 remove_character
  use ConstantsModule, only: LINELENGTH
  implicit none
  private
  public :: collect_arrayhandlers

contains

  subroutine collect_arrayhandlers(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("ExpandArray_int", &
                             test_ExpandArray_int), &
                new_unittest("ExpandArray_dbl", &
                             test_ExpandArray_dbl), &
                new_unittest("ExpandArray_lgp", &
                             test_ExpandArray_lgp), &
                new_unittest("ExpandArray2D_int", &
                             test_ExpandArray2D_int), &
                new_unittest("ExpandArray2D_dbl", &
                             test_ExpandArray2D_dbl), &
                new_unittest("remove_character", &
                             test_remove_character) &
                ]
  end subroutine collect_arrayhandlers

  !> @brief Test 1D int array expansion
  subroutine test_ExpandArray_int(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), allocatable :: a(:)
    integer(I4B) :: i, lb, n1, n2

    n1 = 2 ! starting size
    n2 = 5 ! expanded size
    do lb = -1, 1 ! test default lower bound (1) as well as 0 and -1
      ! allocate/populate array
      allocate (a(lb:(lb + n1 - 1)))
      a(lb) = lb
      a(lb + 1) = lb + 1

      ! resize array and check new size and bounds
      call ExpandArray(a, n2 - n1)
      call check(error, size(a, 1) == n2, &
                 "unexpected size: "//to_string(size(a, 1)))
      call check(error, lbound(a, 1) == lb, &
                 "unexpected lower bound: "//to_string(lbound(a, 1)))
      call check(error, ubound(a, 1) == lb + n2 - 1, &
                 "unexpected upper bound: "//to_string(ubound(a, 1)))
      if (allocated(error)) return

      ! set new array elements and check new/old contents
      do i = lb + n1 - 1, lb + n2 - 1
        a(i) = i
      end do
      do i = lb, lb + n2 - 1
        call check(error, a(i) == i, &
                   "unexpected value "//to_string(a(i)) &
                   //" at i="//to_string(i))
        if (allocated(error)) return
      end do
      deallocate (a)
    end do
  end subroutine test_ExpandArray_int

  !> @brief Test 1D dbl array expansion
  subroutine test_ExpandArray_dbl(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), allocatable :: a(:)
    integer(I4B) :: i, lb, n1, n2

    n1 = 2 ! starting size
    n2 = 5 ! expanded size
    do lb = -1, 1 ! test with default lower bound (1) as well as 0 and -1
      ! allocate/populate array
      allocate (a(lb:(lb + n1 - 1)))
      a(lb) = real(lb)
      a(lb + 1) = real(lb + 1)

      ! resize array and check new size and bounds
      call ExpandArray(a, n2 - n1)
      call check(error, size(a, 1) == n2, &
                 "unexpected size: "//to_string(size(a, 1)))
      call check(error, lbound(a, 1) == lb, &
                 "unexpected lower bound: "//to_string(lbound(a, 1)))
      call check(error, ubound(a, 1) == lb + n2 - 1, &
                 "unexpected upper bound: "//to_string(ubound(a, 1)))
      if (allocated(error)) return

      ! set new array elements and check new/old contents
      do i = lb + n1 - 1, lb + n2 - 1
        a(i) = real(i)
      end do
      do i = lb, lb + n2 - 1
        call check(error, a(i) == real(i), &
                   "unexpected value "//to_string(a(i)) &
                   //" at i="//to_string(i))
        if (allocated(error)) return
      end do
      deallocate (a)
    end do
  end subroutine test_ExpandArray_dbl

  !> @brief Test 1D logical array expansion
  subroutine test_ExpandArray_lgp(error)
    type(error_type), allocatable, intent(out) :: error
    logical(LGP), allocatable :: a(:)
    integer(I4B) :: i, lb, n1, n2

    n1 = 2 ! starting size
    n2 = 5 ! expanded size
    do lb = -1, 1 ! test with default lower bound (1) as well as 0 and -1
      ! allocate/populate array (alternate T/F starting with false)
      allocate (a(lb:(lb + n1 - 1)))
      a(lb) = mod(lb, 2) == 0
      a(lb + 1) = mod(lb + 1, 2) == 0

      ! resize array and check new size and bounds
      call ExpandArray(a, n2 - n1)
      call check(error, size(a, 1) == n2, &
                 "unexpected size: "//to_string(size(a, 1)))
      call check(error, lbound(a, 1) == lb, &
                 "unexpected lower bound: "//to_string(lbound(a, 1)))
      call check(error, ubound(a, 1) == lb + n2 - 1, &
                 "unexpected upper bound: "//to_string(ubound(a, 1)))
      if (allocated(error)) return

      ! set new array elements and check new/old contents
      do i = lb + n1 - 1, lb + n2 - 1
        a(i) = mod(i, 2) == 0
      end do
      do i = lb, lb + n2 - 1
        call check(error, a(i) .eqv. (mod(i, 2) == 0), &
                   "unexpected value "// &
                   merge('t', 'f', a(i)) &
                   //" at i="//to_string(i))
        if (allocated(error)) return
      end do
      deallocate (a)
    end do
  end subroutine test_ExpandArray_lgp

  !> @brief Test 2D int array expansion
  subroutine test_ExpandArray2D_int(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), allocatable :: a(:, :)
    integer(I4B) :: i, lb, n1, n2

    n1 = 2 ! starting size
    n2 = 5 ! expanded size
    do lb = -1, 1 ! test with default lower bound (1) as well as 0 and -1
      ! allocate/populate array and check initial size, with
      ! same lower bound and starting/new size for both dims
      allocate (a(lb:(lb + n1 - 1), lb:(lb + n1 - 1)))
      a(lb, :) = lb
      a(lb + 1, :) = lb + 1
      call check(error, size(a, 1) == n1 .and. size(a, 2) == n1)
      if (allocated(error)) return

      ! resize array and check new size and bounds
      call ExpandArray2D(a, n2 - n1, n2 - n1)
      call check(error, size(a, 1) == n2, &
                 "unexpected dim1 size: "//to_string(size(a, 1)))
      call check(error, size(a, 1) == n2, &
                 "unexpected dim2 size: "//to_string(size(a, 1)))
      call check(error, lbound(a, 1) == lb, &
                 "unexpected dim1 lower bound:"//to_string(lbound(a, 1)))
      call check(error, ubound(a, 1) == lb + n2 - 1, &
                 "unexpected dim1 upper bound:"//to_string(ubound(a, 1)))
      call check(error, lbound(a, 2) == lb, &
                 "unexpected dim2 lower bound:"//to_string(lbound(a, 2)))
      call check(error, ubound(a, 2) == lb + n2 - 1, &
                 "unexpected dim2 upper bound:"//to_string(ubound(a, 2)))
      if (allocated(error)) return

      ! set new elements starting from the new region, check new/old contents
      do i = lb + n1 - 1, lb + n2 - 1
        a(i, :) = i
      end do
      do i = lb, lb + n2 - 1
        if (i < (lb + n1 - 1)) then
          ! old contents, expect uninitialized values in new slots
          call check(error, all(a(i, lb:(lb + n1 - 1)) == i), &
                     "unexpected value "//to_string(a(i, i)) &
                     //" at i="//to_string(i))
        else
          ! new contents, expect all values as set in prior loop
          call check(error, all(a(i, :) == i), &
                     "unexpected value "//to_string(a(i, i)) &
                     //" at i="//to_string(i))
        end if
        if (allocated(error)) return
      end do
      deallocate (a)
    end do
  end subroutine test_ExpandArray2D_int

  !> @brief Test 2D dbl array expansion
  subroutine test_ExpandArray2D_dbl(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), allocatable :: a(:, :)
    integer(I4B) :: i, lb, n1, n2

    n1 = 2 ! starting size
    n2 = 5 ! expanded size
    do lb = -1, 1 ! test with default lower bound (1) as well as 0 and -1
      ! allocate/populate array and check initial size, with
      ! same lower bound and starting/new size for both dims
      allocate (a(lb:(lb + n1 - 1), lb:(lb + n1 - 1)))
      a(lb, :) = real(lb)
      a(lb + 1, :) = real(lb + 1)
      call check(error, size(a, 1) == n1 .and. size(a, 2) == n1)
      if (allocated(error)) return

      ! resize array and check new size and bounds
      call ExpandArray2D(a, n2 - n1, n2 - n1)
      call check(error, size(a, 1) == n2, &
                 "unexpected dim1 size: "//to_string(size(a, 1)))
      call check(error, size(a, 1) == n2, &
                 "unexpected dim2 size: "//to_string(size(a, 1)))
      call check(error, lbound(a, 1) == lb, &
                 "unexpected dim1 lower bound:"//to_string(lbound(a, 1)))
      call check(error, ubound(a, 1) == lb + n2 - 1, &
                 "unexpected dim1 upper bound:"//to_string(ubound(a, 1)))
      call check(error, lbound(a, 2) == lb, &
                 "unexpected dim2 lower bound:"//to_string(lbound(a, 2)))
      call check(error, ubound(a, 2) == lb + n2 - 1, &
                 "unexpected dim2 upper bound:"//to_string(ubound(a, 2)))
      if (allocated(error)) return

      ! set new elements starting from the new region, check new/old contents
      do i = lb + n1 - 1, lb + n2 - 1
        a(i, :) = real(i)
      end do
      do i = lb, lb + n2 - 1
        if (i < (lb + n1 - 1)) then
          ! old contents, expect uninitialized values in new slots
          call check(error, all(a(i, lb:(lb + n1 - 1)) == real(i)), &
                     "unexpected value "//to_string(a(i, i)) &
                     //" at i="//to_string(i))
        else
          ! new contents, expect all values as set in prior loop
          call check(error, all(a(i, :) == real(i)), &
                     "unexpected value "//to_string(a(i, i)) &
                     //" at i="//to_string(i))
        end if
        if (allocated(error)) return
      end do
      deallocate (a)
    end do
  end subroutine test_ExpandArray2D_dbl

  subroutine test_remove_character(error)
    type(error_type), allocatable, intent(out) :: error
    character(len=11), allocatable :: s(:)
    allocate (s(2))
    s(1) = "hello world"
    s(2) = "hello earth"
    call remove_character(s, 1)
    call check(error, s(1) == "hello earth")
  end subroutine test_remove_character

end module TestArrayHandlers
