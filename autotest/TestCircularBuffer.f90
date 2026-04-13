module TestCircularBuffer
  use KindModule, only: I4B, DP
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use CircularBufferModule
  use SimModule, only: store_error, count_errors, initial_message

  implicit none
  private
  public :: collect_circular_buffer

contains

  subroutine collect_circular_buffer(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("push_back", test_push_back), &
      new_unittest("pop_back", test_pop_back), &
      new_unittest("front_and_back", test_front_and_back), &
      new_unittest("front_at", test_front_at), &
      new_unittest("back_at", test_back_at), &
      new_unittest("is_full_and_is_empty", test_is_full_and_is_empty), &
      new_unittest("wrap_around", test_wrap_around), &
      new_unittest("pop_back_underflow", test_pop_back_underflow), &
      new_unittest("index_out_of_bounds", test_index_out_of_bounds) &
    ]
  end subroutine collect_circular_buffer

  subroutine test_push_back(error)
    type(error_type), allocatable, intent(out) :: error
    type(CircularBufferType) :: buffer
    real(DP) :: element(2)
    real(DP), pointer :: retrieved(:)

    buffer = CircularBufferType(3, 2, "cbuf_pba", "TestPath")

    element = [1.0_DP, 2.0_DP]
    call buffer%push_back(element)
    element = [3.0_DP, 4.0_DP]
    call buffer%push_back(element)

    retrieved => buffer%front(1)
    call check(error, all(retrieved == [1.0_DP, 2.0_DP]), &
               "front(1) should return the oldest element")
    if (allocated(error)) return

    retrieved => buffer%front(2)
    call check(error, all(retrieved == [3.0_DP, 4.0_DP]), &
               "front(2) should return the newest element")
  end subroutine test_push_back

  subroutine test_front_at(error)
    type(error_type), allocatable, intent(out) :: error
    type(CircularBufferType) :: buffer
    real(DP) :: element(2)
    real(DP), pointer :: retrieved(:)

    buffer = CircularBufferType(3, 2, "cbuf_fwa", "TestPath")

    element = [1.0_DP, 2.0_DP]; call buffer%push_back(element)
    element = [3.0_DP, 4.0_DP]; call buffer%push_back(element)
    element = [5.0_DP, 6.0_DP]; call buffer%push_back(element)

    retrieved => buffer%front(1)
    call check(error, all(retrieved == [1.0_DP, 2.0_DP]), &
               "front(1) should return the oldest element")
    if (allocated(error)) return

    retrieved => buffer%front(2)
    call check(error, all(retrieved == [3.0_DP, 4.0_DP]), &
               "front(2) should return the second oldest element")
    if (allocated(error)) return

    retrieved => buffer%front(3)
    call check(error, all(retrieved == [5.0_DP, 6.0_DP]), &
               "front(3) should return the newest element")
    if (allocated(error)) return

    ! no-arg front() must still work via the generic
    retrieved => buffer%front()
    call check(error, all(retrieved == [1.0_DP, 2.0_DP]), &
               "front() should still return the oldest element")
  end subroutine test_front_at

  subroutine test_back_at(error)
    type(error_type), allocatable, intent(out) :: error
    type(CircularBufferType) :: buffer
    real(DP) :: element(2)
    real(DP), pointer :: retrieved(:)

    buffer = CircularBufferType(3, 2, "cbuf_ba", "TestPath")

    element = [1.0_DP, 2.0_DP]; call buffer%push_back(element)
    element = [3.0_DP, 4.0_DP]; call buffer%push_back(element)
    element = [5.0_DP, 6.0_DP]; call buffer%push_back(element)

    retrieved => buffer%back(1)
    call check(error, all(retrieved == [5.0_DP, 6.0_DP]), &
               "back(1) should return the newest element")
    if (allocated(error)) return

    retrieved => buffer%back(2)
    call check(error, all(retrieved == [3.0_DP, 4.0_DP]), &
               "back(2) should return the second newest element")
    if (allocated(error)) return

    retrieved => buffer%back(3)
    call check(error, all(retrieved == [1.0_DP, 2.0_DP]), &
               "back(3) should return the oldest element")
    if (allocated(error)) return

    ! no-arg back() must still work via the generic
    retrieved => buffer%back()
    call check(error, all(retrieved == [5.0_DP, 6.0_DP]), &
               "back() should still return the newest element")
  end subroutine test_back_at

  subroutine test_pop_back(error)
    type(error_type), allocatable, intent(out) :: error
    type(CircularBufferType) :: buffer
    real(DP) :: element(2)
    real(DP), pointer :: retrieved(:)

    buffer = CircularBufferType(3, 2, "cbuf_pb", "TestPath")

    element = [1.0_DP, 2.0_DP]
    call buffer%push_back(element)
    element = [3.0_DP, 4.0_DP]
    call buffer%push_back(element)

    call buffer%pop_back()

    call check(error, buffer%size() == 1, "size should be 1 after pop_back")
    if (allocated(error)) return

    retrieved => buffer%back()
    call check(error, all(retrieved == [1.0_DP, 2.0_DP]), &
               "back after pop_back should be the remaining element")
  end subroutine test_pop_back

  subroutine test_front_and_back(error)
    type(error_type), allocatable, intent(out) :: error
    type(CircularBufferType) :: buffer
    real(DP) :: element(2)
    real(DP), pointer :: retrieved(:)

    buffer = CircularBufferType(3, 2, "cbuf_fb", "TestPath")

    element = [1.0_DP, 2.0_DP]
    call buffer%push_back(element)
    element = [3.0_DP, 4.0_DP]
    call buffer%push_back(element)
    element = [5.0_DP, 6.0_DP]
    call buffer%push_back(element)

    retrieved => buffer%front()
    call check(error, all(retrieved == [1.0_DP, 2.0_DP]), &
               "front should return the oldest element")
    if (allocated(error)) return

    retrieved => buffer%back()
    call check(error, all(retrieved == [5.0_DP, 6.0_DP]), &
               "back should return the newest element")
  end subroutine test_front_and_back

  subroutine test_is_full_and_is_empty(error)
    type(error_type), allocatable, intent(out) :: error
    type(CircularBufferType) :: buffer
    real(DP) :: element(2)

    buffer = CircularBufferType(2, 2, "cbuf_fie", "TestPath")

    call check(error, buffer%is_empty(), "buffer should be empty initially")
    if (allocated(error)) return
    call check(error, .not. buffer%is_full(), "buffer should not be full initially")
    if (allocated(error)) return

    element = [1.0_DP, 2.0_DP]
    call buffer%push_back(element)
    call buffer%push_back(element)

    call check(error, buffer%is_full(), "buffer should be full after push to capacity")
    if (allocated(error)) return
    call check(error, .not. buffer%is_empty(), "buffer should not be empty when full")
  end subroutine test_is_full_and_is_empty

  subroutine test_wrap_around(error)
    type(error_type), allocatable, intent(out) :: error
    type(CircularBufferType) :: buffer
    real(DP) :: element(1)
    real(DP), pointer :: retrieved(:)

    ! Fill a capacity-3 buffer then push a 4th element: oldest is overwritten.
    buffer = CircularBufferType(3, 1, "cbuf_wrap", "TestPath")

    element = [1.0_DP]; call buffer%push_back(element)
    element = [2.0_DP]; call buffer%push_back(element)
    element = [3.0_DP]; call buffer%push_back(element)
    element = [4.0_DP]; call buffer%push_back(element)

    call check(error, buffer%size() == 3, &
               "size should remain at capacity after wrap-around push_back")
    if (allocated(error)) return

    retrieved => buffer%front()
    call check(error, retrieved(1) == 2.0_DP, &
               "front after wrap should be the second inserted element")
    if (allocated(error)) return

    retrieved => buffer%back()
    call check(error, retrieved(1) == 4.0_DP, &
               "back after wrap should be the last inserted element")
  end subroutine test_wrap_around

  subroutine test_pop_back_underflow(error)
    type(error_type), allocatable, intent(out) :: error
    type(CircularBufferType) :: buffer
    real(DP) :: element(1)

    buffer = CircularBufferType(3, 1, "cbuf_ufl", "TestPath")
    call initial_message()

    element = [1.0_DP]; call buffer%push_back(element)
    element = [2.0_DP]; call buffer%push_back(element)

    call buffer%pop_back()
    call buffer%pop_back()

    ! Buffer is now empty; one more pop should store an error
    call buffer%pop_back()
    call check(error, count_errors() == 1, &
               "pop_back on empty buffer should store one error")
  end subroutine test_pop_back_underflow

  subroutine test_index_out_of_bounds(error)
    type(error_type), allocatable, intent(out) :: error
    type(CircularBufferType) :: buffer
    real(DP) :: element(1)
    real(DP), pointer :: retrieved(:)

    buffer = CircularBufferType(3, 1, "cbuf_oob", "TestPath")

    element = [1.0_DP]; call buffer%push_back(element)
    element = [2.0_DP]; call buffer%push_back(element)

    ! front(n) with n > size()
    call initial_message()
    retrieved => buffer%front(5)
    call check(error, count_errors() == 1, &
               "front(5) on size-2 buffer should store one error")
    if (allocated(error)) return

    ! back(n) with n > size()
    call initial_message()
    retrieved => buffer%back(5)
    call check(error, count_errors() == 1, &
               "back(5) on size-2 buffer should store one error")
    if (allocated(error)) return

    ! front(0) — zero index is invalid (1-based)
    call initial_message()
    retrieved => buffer%front(0)
    call check(error, count_errors() == 1, &
               "front(0) should store one error (1-based indexing)")
    if (allocated(error)) return

    ! back(0) — zero index is invalid (1-based)
    call initial_message()
    retrieved => buffer%back(0)
    call check(error, count_errors() == 1, &
               "back(0) should store one error (1-based indexing)")
  end subroutine test_index_out_of_bounds

end module TestCircularBuffer