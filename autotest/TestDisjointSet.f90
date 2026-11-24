module TestDisjointSet
  use KindModule, only: I4B
  use testdrive, only: check, error_type, new_unittest, unittest_type
  use DisjointSetModule, only: DisjointSetType
  implicit none
  private
  public :: collect_disjointset

contains

  subroutine collect_disjointset(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("singletons", test_singletons), &
                new_unittest("smallest_root", test_smallest_root), &
                new_unittest("separate_sets", test_separate_sets), &
                new_unittest("reinit", test_reinit) &
                ]
  end subroutine collect_disjointset

  subroutine test_singletons(error)
    type(error_type), allocatable, intent(out) :: error
    type(DisjointSetType) :: ds
    integer(I4B) :: i

    call ds%init(4)
    do i = 1, 4
      call check(error, ds%find(i) == i, "each element should be its own root")
      if (allocated(error)) return
    end do
    call ds%destroy()
  end subroutine test_singletons

  subroutine test_smallest_root(error)
    type(error_type), allocatable, intent(out) :: error
    type(DisjointSetType) :: ds
    integer(I4B) :: i

    ! -- join from the large end so the smallest member is joined last
    call ds%init(6)
    call ds%join(6, 5)
    call ds%join(5, 4)
    call ds%join(4, 2)
    call ds%join(2, 2)
    do i = 2, 6
      if (i == 3) cycle
      call check(error, ds%find(i) == 2, "root should be the smallest member")
      if (allocated(error)) return
    end do
    call check(error, ds%find(1) == 1 .and. ds%find(3) == 3, &
               "elements 1 and 3 should be on their own")
    call ds%destroy()
  end subroutine test_smallest_root

  subroutine test_separate_sets(error)
    type(error_type), allocatable, intent(out) :: error
    type(DisjointSetType) :: ds

    call ds%init(6)
    call ds%join(1, 3)
    call ds%join(6, 4)
    call check(error, ds%find(3) == 1, "set {1, 3} should have root 1")
    if (allocated(error)) return
    call check(error, ds%find(6) == 4, "set {4, 6} should have root 4")
    if (allocated(error)) return
    ! -- joining two sets keeps the smaller root
    call ds%join(6, 3)
    call check(error, ds%find(4) == 1 .and. ds%find(6) == 1, &
               "merged set should have root 1")
    if (allocated(error)) return
    call check(error, ds%find(2) == 2 .and. ds%find(5) == 5, &
               "elements 2 and 5 should be on their own")
    call ds%destroy()
  end subroutine test_separate_sets

  subroutine test_reinit(error)
    type(error_type), allocatable, intent(out) :: error
    type(DisjointSetType) :: ds

    call ds%init(3)
    call ds%join(1, 3)
    call ds%init(3)
    call check(error, ds%find(3) == 3, "init should reset the sets")
    if (allocated(error)) return
    call ds%destroy()
    call check(error,.not. allocated(ds%parent), &
               "destroy should deallocate the parents")
  end subroutine test_reinit

end module TestDisjointSet
