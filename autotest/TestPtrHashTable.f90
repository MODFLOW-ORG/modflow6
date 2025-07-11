module TestPtrHashTable
  use KindModule, only: I4B
  use PtrHashTableModule, only: PtrHashTableType
  use SimModule, only: count_warnings
  use testdrive, only: error_type, unittest_type, new_unittest, check

contains
  subroutine collect_ptrhashtable(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("add_get_values", test_add_get_values), &
                new_unittest("get_nonexisting_value", &
                             test_get_nonexisting_value), &
                new_unittest("count_items", test_count_items), &
                new_unittest("clear_hashtable", test_clear_hashtable), &
                new_unittest("contains_item", test_contains_items), &
                new_unittest("add_duplicate_key", test_add_duplicate_key) &
                ]
  end subroutine collect_ptrhashtable

  !> @brief Test adding and getting values from the hashtable
  !!
  !! Two values are added using a key to identify them
  !! These values are then retrieved by using their keys
  !! and asserted against the expected values
  !<
  subroutine test_add_get_values(error)
    type(error_type), allocatable, intent(out) :: error
    !- Locals
    type(PtrHashTableType) :: hashtable
    class(*), pointer :: val_ptr => null()

    character(len=5) :: name1 = "item1"
    character(len=5) :: name2 = "item2"
    integer(I4B), target :: expected_value1 = 9
    integer(I4B), target :: expected_value2 = 4

    integer(I4B), pointer :: value1 => null()
    integer(I4B), pointer :: value2 => null()

    !- Arrange
    val_ptr => expected_value1
    call hashtable%add(name1, val_ptr)
    val_ptr => expected_value2
    call hashtable%add(name2, val_ptr)

    !- Act
    select type (val => hashtable%get(name1))
    type is (integer(I4B))
      value1 => val
    end select

    select type (val => hashtable%get(name2))
    type is (integer(I4B))
      value2 => val
    end select

    !- Assert
    call check(error, value1 == expected_value1)
    if (allocated(error)) return

    call check(error, value2 == expected_value2)
    if (allocated(error)) return

    ! - Clean up
    call hashtable%clear()

  end subroutine test_add_get_values

  !> @brief Test retrieving a value using a non-existing key
  !!
  !! When the key can't be found the hashtable should return a null pointer
  !<
  subroutine test_get_nonexisting_value(error)
    type(error_type), allocatable, intent(out) :: error
    !- Locals
    type(PtrHashTableType) :: hashtable
    class(*), pointer :: val_ptr => null()

    character(len=17) :: name = "non-existing-item"

    !- Arrange
    !- Act
    val_ptr => hashtable%get(name)

    !- Assert
    call check(error,.not. associated(val_ptr))
    if (allocated(error)) return

    ! - Clean up
    call hashtable%clear()

  end subroutine test_get_nonexisting_value

  !> @brief Test retrieving the number of items in the hashtable
  !!
  !! Three items are added to the hashtable before calling the count method
  !<
  subroutine test_count_items(error)
    type(error_type), allocatable, intent(out) :: error
    !- Locals
    type(PtrHashTableType) :: hashtable
    class(*), pointer :: val_ptr => null()

    integer(I4B), target :: value1 = 1
    integer(I4B), target :: value2 = 2
    integer(I4B), target :: value3 = 3

    integer(I4B) :: cnt = 0

    !- Arrange
    val_ptr => value1
    call hashtable%add("Item1", val_ptr)

    val_ptr => value2
    call hashtable%add("Item2", val_ptr)

    val_ptr => value3
    call hashtable%add("Item3", val_ptr)

    !- Act
    cnt = hashtable%count()

    !- Assert
    call check(error, cnt == 3)
    if (allocated(error)) return

    ! - Clean up
    call hashtable%clear()

  end subroutine test_count_items

  !> @brief Test clearing the hashtable
  !!
  !! When the hashtable is cleared non of the previous added items are gettable anymore.
  !! Furthermore the count of the hashtable is reset to 0.
  !<
  subroutine test_clear_hashtable(error)
    type(error_type), allocatable, intent(out) :: error
    !- Locals
    type(PtrHashTableType) :: hashtable
    class(*), pointer :: val_ptr => null()

    integer(I4B), target :: value1 = 1
    integer(I4B), target :: value2 = 2
    integer(I4B), target :: value3 = 3

    integer(I4B) :: cnt = 0

    !- Arrange
    val_ptr => value1
    call hashtable%add("Item1", val_ptr)

    val_ptr => value2
    call hashtable%add("Item2", val_ptr)

    val_ptr => value3
    call hashtable%add("Item3", val_ptr)

    !- Act
    call hashtable%clear()

    !- Assert
    val_ptr => hashtable%get("Item1")
    call check(error,.not. associated(val_ptr))
    if (allocated(error)) return

    val_ptr => hashtable%get("Item2")
    call check(error,.not. associated(val_ptr))
    if (allocated(error)) return

    val_ptr => hashtable%get("Item3")
    call check(error,.not. associated(val_ptr))
    if (allocated(error)) return

    cnt = hashtable%count()
    call check(error, cnt == 0)
    if (allocated(error)) return

  end subroutine test_clear_hashtable

  !> @brief Test the existence of an item by its key
  !!
  !! This test retrieves the status of 2 items. One of them exists
  !! and the other one doesn't
  !<
  subroutine test_contains_items(error)
    type(error_type), allocatable, intent(out) :: error
    !- Locals
    type(PtrHashTableType) :: hashtable
    class(*), pointer :: val_ptr => null()

    integer(I4B), target :: value = 1
    logical :: found_item1 = .false.
    logical :: found_item2 = .false.

    !- Arrange
    val_ptr => value
    call hashtable%add("Item", val_ptr)

    !- Act
    found_item1 = hashtable%contains("Item")
    found_item2 = hashtable%contains("NonexistingItem")

    !- Assert
    call check(error, found_item1)
    if (allocated(error)) return

    call check(error,.not. found_item2)
    if (allocated(error)) return

    ! - Clean up
    call hashtable%clear()

  end subroutine test_contains_items

  !> @brief Test adding duplicate keys
  !!
  !! When an item is being inserted with an already existing key
  !! a warning is raised. The item is still added but won't be able
  !! to be retrieved using the get method.
  !<
  subroutine test_add_duplicate_key(error)
    type(error_type), allocatable, intent(out) :: error
    !- Locals
    type(PtrHashTableType) :: hashtable
    class(*), pointer :: val_ptr => null()

    integer(I4B), target :: value1 = 1
    integer(I4B), target :: value2 = 2

    !- Arrange
    val_ptr => value1
    call hashtable%add("duplicate_key", val_ptr)

    !- Act
    val_ptr => value2
    call hashtable%add("duplicate_key", val_ptr)

    !- Assert
    call check(error, count_warnings() == 1)
    if (allocated(error)) return

    ! - Clean up
    call hashtable%clear()

  end subroutine test_add_duplicate_key

end module TestPtrHashTable
