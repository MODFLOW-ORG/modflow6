module TestMemoryContainerIterator
  use KindModule, only: I4B, LGP
  use MemoryContainerIteratorModule, only: MemoryContainerIteratorType
  use MemoryStoreModule, only: MemoryStoreType
  use MemoryTypeModule, only: MemoryType
  use testdrive, only: error_type, unittest_type, new_unittest, check

contains
  subroutine collect_memorycontaineriterator(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("iterate_through_container", &
                             test_iterate_through_container) &
                ]
  end subroutine collect_memorycontaineriterator

  !> @brief Iterate through a MemoryContainer
  !!
  !! This test creates an iterator for a container containing 3 MemoryTypes.
  !! It iterates though the container and validates that each type is reached.
  !!
  !! Because the order of the iterator doesn't have to match the order in which
  !! the MemoryTypes have been added an 'iterated' array is used.  A flag in the
  !! array is set to true and at the end it is validated that the entire array
  !! is set to to true (indicating that all memory types have been reached)
  !<
  subroutine test_iterate_through_container(error)
    type(error_type), allocatable, intent(out) :: error
    !- Locals
    type(MemoryStoreType) :: memory_container
    type(MemoryContainerIteratorType), allocatable :: itr

    type(MemoryType), target :: mt1
    type(MemoryType), target :: mt2
    type(MemoryType), target :: mt3

    type(MemoryType), pointer :: current_mt
    integer(I4B) :: mt_index = 0

    logical(LGP) :: iterated(3) = .false.

    !- Arrange.
    mt1%name = "TestName1"
    mt1%path = "TestPath1"
    mt2%name = "TestName2"
    mt2%path = "TestPath2"
    mt3%name = "TestName3"
    mt3%path = "TestPath3"

    current_mt => mt1
    call memory_container%add(current_mt)

    current_mt => mt2
    call memory_container%add(current_mt)

    current_mt => mt3
    call memory_container%add(current_mt)

    itr = memory_container%iterator()

    !- Act.
    current_mt => null()
    do while (itr%has_next())
      call itr%next()
      current_mt => itr%value()

      read (current_mt%name(len_trim(current_mt%name):), '(i1)') mt_index
      iterated(mt_index) = .true.
    end do

    !- Assert.
    call check(error, all(iterated .eqv. .true.))
    if (allocated(error)) return

    !- Clean up.
    call memory_container%clear()

  end subroutine test_iterate_through_container

end module TestMemoryContainerIterator
