module TestList
  use KindModule, only: I4B
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use ConstantsModule, only: LINELENGTH
  use ListModule, only: ListType
  implicit none
  private
  public :: collect_list

  type :: IntNodeType
    integer :: value
  end type IntNodeType

contains

  subroutine collect_list(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("add_count_get_item", &
                             test_add_count_get_item), &
                new_unittest("get_get_index_contains", &
                             test_get_index_contains), &
                new_unittest("get_next_previous_item_reset", &
                             test_get_next_previous_item_reset), &
                new_unittest("insert_after", &
                             test_insert_after), &
                new_unittest("remove_node", &
                             test_remove_node) &
                ]
  end subroutine collect_list

  subroutine test_add_count_get_item(error)
    type(error_type), allocatable, intent(out) :: error
    type(ListType), pointer :: list
    type(IntNodeType), pointer :: n
    class(*), pointer :: p

    allocate (list)
    allocate (n)

    ! empty
    call check(error, list%Count() == 0, "count should be 0")
    if (allocated(error)) return

    ! add one node
    n%value = 1
    p => n
    call list%Add(p)

    ! check count
    call check(error, list%Count() == 1, "count should be 1")
    if (allocated(error)) return

    ! retrieve item
    p => list%GetItem(1)
    call check(error, associated(p, n))
    select type (item => p)
    type is (IntNodeType)
      call check(error, item%value == 1, "wrong value")
    class default
      call check(error, .false., "wrong node type")
    end select
    if (allocated(error)) return

    call list%Clear()
    deallocate (list)
    deallocate (n)
  end subroutine test_add_count_get_item

  subroutine test_get_index_contains(error)
    type(error_type), allocatable, intent(out) :: error
    type(ListType), pointer :: list
    type(IntNodeType), pointer :: n1, n2
    class(*), pointer :: p
    integer(I4B) :: i

    allocate (list)
    allocate (n1)
    allocate (n2)

    ! add nodes
    n1%value = 1
    n2%value = 2
    p => n1
    call list%Add(p)
    p => n2
    call list%Add(p)

    ! check count
    call check(error, list%Count() == 2, "count should be 1")
    if (allocated(error)) return

    ! check get index
    i = list%GetIndex(p)
    call check(error, i == 2, "wrong index")
    if (allocated(error)) return

    ! check contains
    p => n1
    call check(error, list%ContainsObject(p), "should contain n1")
    if (allocated(error)) return
    p => n2
    call check(error, list%ContainsObject(p), "should contain n2")
    if (allocated(error)) return

    call list%Clear()
    deallocate (list)
    deallocate (n1)
    deallocate (n2)
  end subroutine test_get_index_contains

  subroutine test_get_next_previous_item_reset(error)
    type(error_type), allocatable, intent(out) :: error
    type(ListType), pointer :: list
    type(IntNodeType), pointer :: n1, n2, n3
    class(*), pointer :: p

    allocate (list)
    allocate (n1)
    allocate (n2)
    allocate (n3)

    ! add nodes
    n1%value = 1
    n2%value = 2
    n3%value = 3
    p => n1
    call list%Add(p)
    p => n2
    call list%Add(p)
    p => n3
    call list%Add(p)

    ! check count
    call check(error, list%Count() == 3, "count should be 3")
    if (allocated(error)) return

    ! check get next/previous item
    p => list%GetNextItem()
    call check(error, associated(p, n1))
    p => list%GetNextItem()
    call check(error, associated(p, n2))
    p => list%GetPreviousItem()
    call check(error, associated(p, n1))
    p => list%GetNextItem()
    call check(error, associated(p, n2))
    p => list%GetNextItem()
    call check(error, associated(p, n3))
    p => list%GetNextItem()
    call check(error, (.not. associated(p)))
    call list%Reset()
    p => list%GetPreviousItem()
    call check(error, (.not. associated(p)))

    call list%Clear()
    deallocate (list)
    deallocate (n1)
    deallocate (n2)
    deallocate (n3)
  end subroutine test_get_next_previous_item_reset

  subroutine test_insert_after(error)
    type(error_type), allocatable, intent(out) :: error
    type(ListType), pointer :: list
    type(IntNodeType), pointer :: n1, n2, n3
    class(*), pointer :: p

    allocate (list)
    allocate (n1)
    allocate (n2)
    allocate (n3)

    ! add nodes 1 and 3
    n1%value = 1
    n2%value = 2
    n3%value = 3
    p => n1
    call list%Add(p)
    p => n3
    call list%Add(p)

    ! check count
    call check(error, list%Count() == 2, "count should be 2")
    if (allocated(error)) return

    ! insert item after first item
    p => n2
    call list%InsertAfter(p, 1)

    ! check count
    call check(error, list%Count() == 3, "count should be 3")
    if (allocated(error)) return

    ! check get next/previous item
    call list%Reset()
    p => list%GetNextItem()
    call check(error, associated(p, n1))
    p => list%GetNextItem()
    call check(error, associated(p, n2))
    p => list%GetNextItem()
    call check(error, associated(p, n3))
    if (allocated(error)) return

    deallocate (list)
    deallocate (n1)
    deallocate (n2)
    deallocate (n3)
  end subroutine test_insert_after

  subroutine test_remove_node(error)
    type(error_type), allocatable, intent(out) :: error
    type(ListType), pointer :: list
    type(IntNodeType), pointer :: n1, n2, n3
    class(*), pointer :: p

    allocate (list)
    allocate (n1)
    allocate (n2)
    allocate (n3)

    ! add nodes
    n1%value = 1
    n2%value = 2
    n3%value = 3
    p => n1
    call list%Add(p)
    p => n2
    call list%Add(p)
    p => n3
    call list%Add(p)

    ! check count
    call check(error, list%Count() == 3, "count should be 3")
    if (allocated(error)) return

    ! remove first node
    call list%RemoveNode(1, .false.)
    call check(error, list%Count() == 2, "count should be 2")
    p => list%GetItem(1)
    call check(error, associated(p, n2))
    p => list%GetItem(2)
    call check(error, associated(p, n3))

    ! remove last node
    call list%RemoveNode(2, .false.)
    call check(error, list%Count() == 1, "count should be 1")
    p => list%GetItem(1)
    call check(error, associated(p, n2))

    call list%Clear()
    deallocate (list)
    deallocate (n1)
    deallocate (n2)
    deallocate (n3)
  end subroutine test_remove_node

end module TestList
