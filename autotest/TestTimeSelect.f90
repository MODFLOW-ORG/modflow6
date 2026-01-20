module TestTimeSelect
  use KindModule, only: I4B, DP, LGP
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       to_string, unittest_type
  use TimeSelectModule, only: TimeSelectType
  use ConstantsModule, only: LINELENGTH
  implicit none
  private
  public :: collect_timeselect

contains
  subroutine collect_timeselect(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("is_increasing", test_is_increasing), &
                new_unittest("select", test_select), &
                new_unittest("extend_and_sort", &
                             test_extend_and_sort) &
                ]
  end subroutine collect_timeselect

  subroutine test_is_increasing(error)
    type(error_type), allocatable, intent(out) :: error
    type(TimeSelectType) :: ts

    call ts%init()
    call ts%expand(3)

    ! increasing
    ts%times = (/0.0_DP, 1.0_DP, 2.0_DP/)
    call check(error, ts%increasing())

    ! not decreasing (duplicates)
    ts%times = (/0.0_DP, 0.0_DP, 2.0_DP/)
    call check(error,.not. ts%increasing())

    ! decreasing
    ts%times = (/2.0_DP, 1.0_DP, 0.0_DP/)
    call check(error,.not. ts%increasing())
  end subroutine

  subroutine test_select(error)
    type(error_type), allocatable, intent(out) :: error
    type(TimeSelectType) :: ts
    logical(LGP) :: changed

    call ts%init()
    call ts%expand(3)
    ts%times = (/0.0_DP, 1.0_DP, 2.0_DP/)
    call check( &
      error, &
      size(ts%times) == 3, &
      "expected size 3, got"//to_string(size(ts%times)))

    ! empty slice
    call ts%select(1.1_DP, 1.9_DP)
    call check( &
      error, &
      ts%selection(1) == -1 .and. ts%selection(2) == -1, &
      "empty slice failed, got ["// &
      to_string(ts%selection(1))//","//to_string(ts%selection(2))//"]")

    ! single-item slice
    call ts%select(0.5_DP, 1.5_DP)
    call check( &
      error, &
      ts%selection(1) == 2 .and. ts%selection(2) == 2, &
      "1-item slice failed, got ["// &
      to_string(ts%selection(1))//","//to_string(ts%selection(2))//"]")

    ! multi-item slice
    changed = .false.
    call ts%select(0.5_DP, 2.5_DP, changed=changed)
    call check(error, changed)
    call check( &
      error, &
      ts%selection(1) == 2 .and. ts%selection(2) == 3, &
      "2-item slice failed, got ["// &
      to_string(ts%selection(1))//","//to_string(ts%selection(2))//"]")

    ! no-change
    call ts%select(0.1_DP, 2.5_DP, changed=changed)
    call check(error,.not. changed)
    call check( &
      error, &
      ts%selection(1) == 2 .and. ts%selection(2) == 3, &
      "2-item slice failed, got ["// &
      to_string(ts%selection(1))//","//to_string(ts%selection(2))//"]")

    ! lower bound equal to a time value (exclusive lower bound)
    call ts%select(0.0_DP, 2.5_DP)
    call check( &
      error, &
      ts%selection(1) == 2 .and. ts%selection(2) == 3, &
      "lb eq slice failed, got [" &
      //to_string(ts%selection(1))//","//to_string(ts%selection(2))//"]")

    ! upper bound equal to a time value
    call ts%select(-0.5_DP, 2.0_DP)
    call check( &
      error, &
      ts%selection(1) == 1 .and. ts%selection(2) == 3, &
      "ub eq slice failed, got [" &
      //to_string(ts%selection(1))//","//to_string(ts%selection(2))//"]")

    ! both bounds equal to a time value (lower exclusive, upper inclusive)
    call ts%select(0.0_DP, 2.0_DP)
    call check( &
      error, &
      ts%selection(1) == 2 .and. ts%selection(2) == 3, &
      "lb ub eq slice failed, got [" &
      //to_string(ts%selection(1))//","//to_string(ts%selection(2))//"]")

    ! Test explicit boundary behavior: interval (0.0, 1.0]
    ! should include 1.0 but not 0.0
    call ts%select(0.0_DP, 1.0_DP)
    call check( &
      error, &
      ts%selection(1) == 2 .and. ts%selection(2) == 2, &
      "boundary (0,1] slice failed, got [" &
      //to_string(ts%selection(1))//","//to_string(ts%selection(2))//"]")

    ! Test explicit boundary behavior: interval (1.0, 2.0]
    ! should include 2.0 but not 1.0
    call ts%select(1.0_DP, 2.0_DP)
    call check( &
      error, &
      ts%selection(1) == 3 .and. ts%selection(2) == 3, &
      "boundary (1,2] slice failed, got [" &
      //to_string(ts%selection(1))//","//to_string(ts%selection(2))//"]")

  end subroutine test_select

  subroutine test_extend_and_sort(error)
    type(error_type), allocatable, intent(out) :: error
    type(TimeSelectType) :: ts
    real(DP) :: a(3)

    a = (/0.0_DP, 2.0_DP, 1.0_DP/)

    call ts%init()
    call ts%extend(a)
    call check(error, ts%increasing())

  end subroutine test_extend_and_sort

end module TestTimeSelect
