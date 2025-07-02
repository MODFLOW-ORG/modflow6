module TestKeyValueList
  use KindModule, only: I4B
  use KeyValueListModule, only: KeyValueListType
  use testdrive, only: error_type, unittest_type, new_unittest, check

contains
  subroutine collect_keyvaluelist(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("add_get_values", test_add_get_values), &
                new_unittest("get_nonexisting_value", &
                             test_get_nonexisting_value), &
                new_unittest("count_items", test_count_items), &
                new_unittest("clear_list", test_clear_list) &
                ]
  end subroutine collect_keyvaluelist

  !> @brief Test adding and getting values from the list
  !!
  !! Two values are added using a key to identify them
  !! These values are then retrieved by using their keys
  !! and asserted against the expected values
  !<
  subroutine test_add_get_values(error)
    type(error_type), allocatable, intent(out) :: error
    !- Locals
    type(KeyValueListType) :: list
    class(*), pointer :: val_ptr => null()

    character(len=5) :: name1 = "item1"
    character(len=5) :: name2 = "item2"
    integer(I4B), target :: expected_value1 = 9
    integer(I4B), target :: expected_value2 = 4

    integer(I4B), pointer :: value1 => null()
    integer(I4B), pointer :: value2 => null()

    !- Arrange
    val_ptr => expected_value1
    call list%add(name1, val_ptr)
    val_ptr => expected_value2
    call list%add(name2, val_ptr)

    !- Act
    select type (val => list%get(name1))
    type is (integer(I4B))
      value1 => val
    end select

    select type (val => list%get(name2))
    type is (integer(I4B))
      value2 => val
    end select

    !- Assert
    call check(error, value1 == expected_value1)
    if (allocated(error)) return

    call check(error, value2 == expected_value2)
    if (allocated(error)) return

    ! - Clean up
    call list%clear()

  end subroutine test_add_get_values

  !> @brief Test retrieving a value using a non-existing key
  !!
  !! When the key can't be found the list should return a null pointer
  !<
  subroutine test_get_nonexisting_value(error)
    type(error_type), allocatable, intent(out) :: error
    !- Locals
    type(KeyValueListType) :: list
    class(*), pointer :: val_ptr => null()

    character(len=17) :: name = "non-existing-item"

    !- Arrange
    !- Act
    val_ptr => list%get(name)

    !- Assert
    call check(error,.not. associated(val_ptr))
    if (allocated(error)) return

    ! - Clean up
    call list%clear()

  end subroutine test_get_nonexisting_value

  !> @brief Test retrieving the number of items in the list
  !!
  !! Three items are added to the list before calling the count method
  !<
  subroutine test_count_items(error)
    type(error_type), allocatable, intent(out) :: error
    !- Locals
    type(KeyValueListType) :: list
    class(*), pointer :: val_ptr => null()

    integer(I4B), target :: value1 = 1
    integer(I4B), target :: value2 = 2
    integer(I4B), target :: value3 = 3

    integer(I4B) :: cnt = 0

    !- Arrange
    val_ptr => value1
    call list%add("Item1", val_ptr)

    val_ptr => value2
    call list%add("Item2", val_ptr)

    val_ptr => value3
    call list%add("Item3", val_ptr)

    !- Act
    cnt = list%count()

    !- Assert
    call check(error, cnt == 3)
    if (allocated(error)) return

    ! - Clean up
    call list%clear()

  end subroutine test_count_items

  !> @brief Test clearing the list
  !!
  !! When the list is cleared non of the previous added items are gettable anymore.
  !! Furthermore the count of the list is reset to 0.
  !<
  subroutine test_clear_list(error)
    type(error_type), allocatable, intent(out) :: error
    !- Locals
    type(KeyValueListType) :: list
    class(*), pointer :: val_ptr => null()

    integer(I4B), target :: value1 = 1
    integer(I4B), target :: value2 = 2
    integer(I4B), target :: value3 = 3

    integer(I4B) :: cnt = 0

    !- Arrange
    val_ptr => value1
    call list%add("Item1", val_ptr)

    val_ptr => value2
    call list%add("Item2", val_ptr)

    val_ptr => value3
    call list%add("Item3", val_ptr)

    !- Act
    call list%clear()

    !- Assert
    val_ptr => list%get("Item1")
    call check(error,.not. associated(val_ptr))
    if (allocated(error)) return

    val_ptr => list%get("Item2")
    call check(error,.not. associated(val_ptr))
    if (allocated(error)) return

    val_ptr => list%get("Item3")
    call check(error,.not. associated(val_ptr))
    if (allocated(error)) return

    cnt = list%count()
    call check(error, cnt == 0)
    if (allocated(error)) return

  end subroutine test_clear_list

end module TestKeyValueList
