module TestHashTable
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DNODATA, DZERO
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       to_string, unittest_type
  use HashTableModule, only: HashTableType, hash_table_cr, hash_table_da
  implicit none
  private
  public :: collect_hashtable

contains

  subroutine collect_hashtable(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("add_and_get_value", &
                             test_add_and_get_value) &
                ]
  end subroutine collect_hashtable

  subroutine test_add_and_get_value(error)
    type(error_type), allocatable, intent(out) :: error
    type(HashTableType), pointer :: map
    integer(I4B) :: i, n

    call hash_table_cr(map)

    n = 3
    do i = 1, n
      call map%add(to_string(i), i)
    end do

    do i = 1, n
      call check(error, map%get(to_string(i)) == i, &
                 'wrong value for '//to_string(i))
    end do

    call hash_table_da(map)

  end subroutine test_add_and_get_value

end module TestHashTable
