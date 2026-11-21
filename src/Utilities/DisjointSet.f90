!> @brief Disjoint set (union-find) of the integers 1 to n
!!
!! Each set is identified by its smallest member, so the root of a set does
!! not depend on the order in which sets are joined.
!<
module DisjointSetModule

  use KindModule, only: I4B

  implicit none
  private
  public :: DisjointSetType

  !> @brief Disjoint set of the integers 1 to n
  !<
  type :: DisjointSetType
    integer(I4B), dimension(:), allocatable :: parent !< parent of each element
  contains
    procedure :: init
    procedure :: find
    procedure :: join
    procedure :: destroy
  end type DisjointSetType

contains

  !> @brief Make n sets of one element each
  !<
  subroutine init(this, n)
    class(DisjointSetType), intent(inout) :: this !< disjoint set
    integer(I4B), intent(in) :: n !< number of elements
    ! -- local
    integer(I4B) :: i

    if (allocated(this%parent)) deallocate (this%parent)
    allocate (this%parent(n))
    do i = 1, n
      this%parent(i) = i
    end do
  end subroutine init

  !> @brief Return the smallest member of the set holding element i
  !<
  function find(this, i) result(r)
    class(DisjointSetType), intent(inout) :: this !< disjoint set
    integer(I4B), intent(in) :: i !< element
    integer(I4B) :: r

    r = i
    do while (this%parent(r) /= r)
      ! -- path halving
      this%parent(r) = this%parent(this%parent(r))
      r = this%parent(r)
    end do
  end function find

  !> @brief Join the sets holding elements i and j
  !<
  subroutine join(this, i, j)
    class(DisjointSetType), intent(inout) :: this !< disjoint set
    integer(I4B), intent(in) :: i, j !< elements
    ! -- local
    integer(I4B) :: ri, rj

    ri = this%find(i)
    rj = this%find(j)
    ! -- the smaller root stays the root
    if (ri < rj) then
      this%parent(rj) = ri
    else if (rj < ri) then
      this%parent(ri) = rj
    end if
  end subroutine join

  !> @brief Deallocate the disjoint set
  !<
  subroutine destroy(this)
    class(DisjointSetType), intent(inout) :: this !< disjoint set

    if (allocated(this%parent)) deallocate (this%parent)
  end subroutine destroy

end module DisjointSetModule
