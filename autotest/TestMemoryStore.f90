module TestMemoryStore
  use MemoryStoreModule, only: MemoryStoreType
  use MemoryTypeModule, only: MemoryType
  use testdrive, only: error_type, unittest_type, new_unittest, check

contains
  subroutine collect_memorystore(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("add_get_values", test_add_get_values), &
                new_unittest("get_nonexisting_value", &
                             test_get_nonexisting_value) &
                ]
  end subroutine collect_memorystore

  !> @brief Test adding and getting values from the container
  !!
  !! One MemoryType is added. A key is constructed from its
  !! name and path. The MemoryType is then stored and
  !! retrieved using that key from the internal container
  !<
  subroutine test_add_get_values(error)
    type(error_type), allocatable, intent(out) :: error
    !- Locals
    type(MemoryStoreType) :: container

    type(MemoryType), target :: mt
    type(MemoryType), pointer :: mt_ptr

    character(*), parameter :: name = "TestName"
    character(*), parameter :: path = "TestPath"

    !- Arrange
    mt%name = name
    mt%path = path
    mt_ptr => mt

    !- Act
    call container%add(mt_ptr)

    mt_ptr => null()
    mt_ptr => container%get(name, path)

    !- Assert
    call check(error, associated(mt_ptr))
    if (allocated(error)) return

    call check(error, mt_ptr%name == name)
    if (allocated(error)) return

    call check(error, mt_ptr%path == path)
    if (allocated(error)) return

    !- Cleanup
    call container%clear()
  end subroutine test_add_get_values

  !> @brief Test retrieving a MemoryType using a non-existing name and path
  !!
  !! When the key can't be found the MemoryStore should return a null pointer
  !<
  subroutine test_get_nonexisting_value(error)
    type(error_type), allocatable, intent(out) :: error
    !- Locals
    type(MemoryStoreType) :: container
    type(MemoryType), pointer :: mt_ptr
    character(*), parameter :: name = "FakeName"
    character(*), parameter :: path = "FakePath"

    !- Arrange
    !- Act
    mt_ptr => container%get(name, path)

    !- Assert
    call check(error,.not. associated(mt_ptr))
    if (allocated(error)) return

    !- Cleanup
    call container%clear()
  end subroutine test_get_nonexisting_value

end module TestMemoryStore
