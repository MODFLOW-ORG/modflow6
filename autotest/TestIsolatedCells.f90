module TestIsolatedCells
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE
  use testdrive, only: check, error_type, new_unittest, unittest_type
  use IsolatedCellsModule, only: row_is_controlled, find_groups, &
                                 boundary_links, merge_groups, &
                                 IsolatedTrackerType
  implicit none
  private
  public :: collect_isolatedcells

  integer(I4B), parameter :: NCHAIN = 6 !< cells in the test chain

contains

  subroutine collect_isolatedcells(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("row_controlled", test_row_controlled), &
                new_unittest("one_group_determined", test_one_group_determined), &
                new_unittest("zero_term_splits", test_zero_term_splits), &
                new_unittest("one_sided_term_joins", &
                             test_one_sided_term_joins), &
                new_unittest("row_offset", test_row_offset), &
                new_unittest("boundary_links", test_boundary_links), &
                new_unittest("merge_determined", test_merge_determined), &
                new_unittest("merge_undetermined", test_merge_undetermined), &
                new_unittest("merge_no_link", test_merge_no_link), &
                new_unittest("tracker_status", test_tracker_status), &
                new_unittest("tracker_cause", test_tracker_cause) &
                ]
  end subroutine collect_isolatedcells

  !> @brief Build rows r1 to r2 of a chain of cells
  !!
  !! Row i has conductance cond(i-1) to cell i-1 and cond(i) to cell i+1, and
  !! storage sto(i), with the diagonal the negative sum, as MODFLOW assembles
  !! it. Columns are global cell numbers.
  !<
  subroutine build_chain(r1, r2, cond, sto, ia, ja, a)
    integer(I4B), intent(in) :: r1, r2 !< first and last row
    real(DP), dimension(NCHAIN - 1), intent(in) :: cond !< conductances
    real(DP), dimension(NCHAIN), intent(in) :: sto !< storage terms
    integer(I4B), dimension(:), allocatable, intent(out) :: ia, ja
    real(DP), dimension(:), allocatable, intent(out) :: a
    integer(I4B) :: i, k, n

    allocate (ia(r2 - r1 + 2), ja(3 * (r2 - r1 + 1)), a(3 * (r2 - r1 + 1)))
    k = 0
    ia(1) = 1
    do i = r1, r2
      n = i - r1 + 1
      k = k + 1
      ja(k) = i
      a(k) = -sto(i)
      if (i > 1) then
        a(k) = a(k) - cond(i - 1)
        k = k + 1
        ja(k) = i - 1
        a(k) = cond(i - 1)
      end if
      if (i < NCHAIN) then
        a(ia(n)) = a(ia(n)) - cond(i)
        k = k + 1
        ja(k) = i + 1
        a(k) = cond(i)
      end if
      ia(n + 1) = k + 1
    end do
  end subroutine build_chain

  !> @brief Set the control flag of each row
  !<
  subroutine row_controls(nrow, ia, a, row_control)
    integer(I4B), intent(in) :: nrow
    integer(I4B), dimension(nrow + 1), intent(in) :: ia
    real(DP), dimension(*), intent(in) :: a
    integer(I4B), dimension(nrow), intent(out) :: row_control
    integer(I4B) :: n

    do n = 1, nrow
      row_control(n) = row_is_controlled(ia(n + 1) - ia(n), a(ia(n)))
    end do
  end subroutine row_controls

  subroutine test_row_controlled(error)
    type(error_type), allocatable, intent(out) :: error

    ! -- conductance only: the row sums to zero
    call check(error, row_is_controlled(3, [-3.0_DP, 1.0_DP, 2.0_DP]) == 0, &
               "conductance-only row should not be controlled")
    if (allocated(error)) return
    ! -- storage makes the row sum nonzero
    call check(error, row_is_controlled(3, [-3.5_DP, 1.0_DP, 2.0_DP]) == 1, &
               "row with storage should be determined")
    if (allocated(error)) return
    ! -- a constant head or inactive cell is an identity row
    call check(error, row_is_controlled(1, [DONE]) == 1, &
               "identity row should be controlled")
    if (allocated(error)) return
    ! -- an all-zero row has nothing to control it
    call check(error, row_is_controlled(2, [DZERO, DZERO]) == 0, &
               "zero row should not be controlled")
    if (allocated(error)) return
    ! -- a sum at rounding level is not storage
    call check(error, &
               row_is_controlled(3, [-0.3_DP, 0.1_DP, 0.2_DP]) == 0, &
               "rounding-level sum should not be controlled")
  end subroutine test_row_controlled

  subroutine test_one_group_determined(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), dimension(:), allocatable :: ia, ja
    real(DP), dimension(:), allocatable :: a
    integer(I4B) :: row_control(NCHAIN), igid(NCHAIN), group_determined(NCHAIN), n
    real(DP) :: cond(NCHAIN - 1), sto(NCHAIN)

    do n = 1, NCHAIN - 1
      cond(n) = DONE
    end do
    do n = 1, NCHAIN
      sto(n) = DZERO
    end do
    sto(1) = 0.1_DP
    call build_chain(1, NCHAIN, cond, sto, ia, ja, a)
    call row_controls(NCHAIN, ia, a, row_control)
    call find_groups(NCHAIN, 0, ia, ja, a, row_control, igid, group_determined)
    do n = 1, NCHAIN
      call check(error, igid(n) == 1, "all cells should be in group 1")
      if (allocated(error)) return
    end do
    call check(error, group_determined(1) == 1, "group 1 should be determined")
  end subroutine test_one_group_determined

  subroutine test_zero_term_splits(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), dimension(:), allocatable :: ia, ja
    real(DP), dimension(:), allocatable :: a
    integer(I4B) :: row_control(NCHAIN), igid(NCHAIN), group_determined(NCHAIN), n
    real(DP) :: cond(NCHAIN - 1), sto(NCHAIN)

    ! -- no conductance between cells 3 and 4, as between two dry cells
    do n = 1, NCHAIN - 1
      cond(n) = DONE
    end do
    cond(3) = DZERO
    do n = 1, NCHAIN
      sto(n) = DZERO
    end do
    sto(1) = 0.1_DP
    call build_chain(1, NCHAIN, cond, sto, ia, ja, a)
    call row_controls(NCHAIN, ia, a, row_control)
    call find_groups(NCHAIN, 0, ia, ja, a, row_control, igid, group_determined)
    do n = 1, 3
      call check(error, igid(n) == 1, "cells 1-3 should be in group 1")
      if (allocated(error)) return
    end do
    do n = 4, NCHAIN
      call check(error, igid(n) == 4, "cells 4-6 should be in group 4")
      if (allocated(error)) return
    end do
    call check(error, group_determined(1) == 1, "group 1 should be determined")
    if (allocated(error)) return
    call check(error, group_determined(4) == 0, &
               "group 4 should not be determined")
  end subroutine test_zero_term_splits

  subroutine test_one_sided_term_joins(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), dimension(:), allocatable :: ia, ja
    real(DP), dimension(:), allocatable :: a
    integer(I4B) :: row_control(NCHAIN), igid(NCHAIN), group_determined(NCHAIN)
    integer(I4B) :: n, k
    real(DP) :: cond(NCHAIN - 1), sto(NCHAIN)

    do n = 1, NCHAIN - 1
      cond(n) = DONE
    end do
    do n = 1, NCHAIN
      sto(n) = DZERO
    end do
    sto(1) = 0.1_DP
    call build_chain(1, NCHAIN, cond, sto, ia, ja, a)
    ! -- remove the term for cell 4 in row 3 only, as a nonsymmetric matrix
    !    can have; row 4 still links the two cells
    do k = ia(3), ia(4) - 1
      if (ja(k) == 4) then
        a(ia(3)) = a(ia(3)) + a(k)
        a(k) = DZERO
      end if
    end do
    call row_controls(NCHAIN, ia, a, row_control)
    call find_groups(NCHAIN, 0, ia, ja, a, row_control, igid, group_determined)
    do n = 1, NCHAIN
      call check(error, igid(n) == 1, "all cells should be in group 1")
      if (allocated(error)) return
    end do
    call check(error, group_determined(1) == 1, "group 1 should be determined")
  end subroutine test_one_sided_term_joins

  subroutine test_row_offset(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), dimension(:), allocatable :: ia, ja
    real(DP), dimension(:), allocatable :: a
    integer(I4B) :: row_control(NCHAIN), igid(NCHAIN), group_determined(NCHAIN)
    integer(I4B) :: n, roff
    real(DP) :: cond(NCHAIN - 1), sto(NCHAIN)

    ! -- the same split chain numbered from 11, as rows of a later process
    roff = 10
    do n = 1, NCHAIN - 1
      cond(n) = DONE
    end do
    cond(3) = DZERO
    do n = 1, NCHAIN
      sto(n) = DZERO
    end do
    call build_chain(1, NCHAIN, cond, sto, ia, ja, a)
    do n = 1, size(ja)
      ja(n) = ja(n) + roff
    end do
    call row_controls(NCHAIN, ia, a, row_control)
    call find_groups(NCHAIN, roff, ia, ja, a, row_control, igid, group_determined)
    call check(error, igid(1) == 11 .and. igid(3) == 11, &
               "cells 1-3 should be in group 11")
    if (allocated(error)) return
    call check(error, igid(4) == 14 .and. igid(NCHAIN) == 14, &
               "cells 4-6 should be in group 14")
    if (allocated(error)) return
    call check(error, group_determined(1) == 0 .and. group_determined(4) == 0, &
               "neither group should be determined")
  end subroutine test_row_offset

  !> @brief Split the chain between two processes at cells 3 and 4
  !<
  subroutine split_chain(cond, sto, determined_a, determined_b, nb, brow, bgid, &
                         ng, gid, determined, nlink, lgid, lrow)
    real(DP), dimension(NCHAIN - 1), intent(in) :: cond
    real(DP), dimension(NCHAIN), intent(in) :: sto
    integer(I4B), dimension(3), intent(out) :: determined_a, determined_b
    integer(I4B), intent(out) :: nb, ng, nlink
    integer(I4B), dimension(:), allocatable, intent(out) :: brow, bgid
    integer(I4B), dimension(:), allocatable, intent(out) :: gid, determined
    integer(I4B), dimension(:), allocatable, intent(out) :: lgid, lrow
    integer(I4B), dimension(:), allocatable :: ia, ja, br1, bg1, lg1, lr1
    integer(I4B), dimension(:), allocatable :: br2, bg2, lg2, lr2
    real(DP), dimension(:), allocatable :: a
    integer(I4B) :: row_control(3), igid(3), nb1, nb2, nl1, nl2

    call build_chain(1, 3, cond, sto, ia, ja, a)
    call row_controls(3, ia, a, row_control)
    call find_groups(3, 0, ia, ja, a, row_control, igid, determined_a)
    call boundary_links(3, 0, ia, ja, a, igid, nb1, br1, bg1, nl1, lg1, lr1)
    call build_chain(4, NCHAIN, cond, sto, ia, ja, a)
    call row_controls(3, ia, a, row_control)
    call find_groups(3, 3, ia, ja, a, row_control, igid, determined_b)
    call boundary_links(3, 3, ia, ja, a, igid, nb2, br2, bg2, nl2, lg2, lr2)
    ! -- the gathered lists, with one group on each side of the split
    nb = nb1 + nb2
    brow = [br1, br2]
    bgid = [bg1, bg2]
    ng = 2
    gid = [1, 4]
    determined = [determined_a(1), determined_b(1)]
    nlink = nl1 + nl2
    lgid = [lg1, lg2]
    lrow = [lr1, lr2]
  end subroutine split_chain

  subroutine test_boundary_links(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), dimension(:), allocatable :: brow, bgid, gid, determined
    integer(I4B), dimension(:), allocatable :: lgid, lrow
    integer(I4B) :: determined_a(3), determined_b(3), nb, ng, nlink, n
    real(DP) :: cond(NCHAIN - 1), sto(NCHAIN)

    do n = 1, NCHAIN - 1
      cond(n) = DONE
    end do
    do n = 1, NCHAIN
      sto(n) = DZERO
    end do
    call split_chain(cond, sto, determined_a, determined_b, nb, brow, bgid, &
                     ng, gid, determined, nlink, lgid, lrow)
    call check(error, nb == 2, "rows 3 and 4 should be boundary rows")
    if (allocated(error)) return
    call check(error, brow(1) == 3 .and. bgid(1) == 1, &
               "row 3 should be listed in group 1")
    if (allocated(error)) return
    call check(error, brow(2) == 4 .and. bgid(2) == 4, &
               "row 4 should be listed in group 4")
    if (allocated(error)) return
    call check(error, nlink == 2, "each side should have one link")
    if (allocated(error)) return
    call check(error, lgid(1) == 1 .and. lrow(1) == 4, &
               "group 1 should link to row 4")
    if (allocated(error)) return
    call check(error, lgid(2) == 4 .and. lrow(2) == 3, &
               "group 4 should link to row 3")
  end subroutine test_boundary_links

  subroutine test_merge_determined(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), dimension(:), allocatable :: brow, bgid, gid, determined
    integer(I4B), dimension(:), allocatable :: lgid, lrow
    integer(I4B) :: determined_a(3), determined_b(3), nb, ng, nlink, n
    real(DP) :: cond(NCHAIN - 1), sto(NCHAIN)

    ! -- storage on one side determines the cells on the other
    do n = 1, NCHAIN - 1
      cond(n) = DONE
    end do
    do n = 1, NCHAIN
      sto(n) = DZERO
    end do
    sto(1) = 0.1_DP
    call split_chain(cond, sto, determined_a, determined_b, nb, brow, bgid, &
                     ng, gid, determined, nlink, lgid, lrow)
    call check(error, determined(1) == 1 .and. determined(2) == 0, &
               "only the first side should be determined before merging")
    if (allocated(error)) return
    call merge_groups(nb, brow, bgid, ng, gid, determined, nlink, lgid, lrow)
    call check(error, determined(1) == 1 .and. determined(2) == 1, &
               "both sides should be determined after merging")
  end subroutine test_merge_determined

  subroutine test_merge_undetermined(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), dimension(:), allocatable :: brow, bgid, gid, determined
    integer(I4B), dimension(:), allocatable :: lgid, lrow
    integer(I4B) :: determined_a(3), determined_b(3), nb, ng, nlink, n
    real(DP) :: cond(NCHAIN - 1), sto(NCHAIN)

    do n = 1, NCHAIN - 1
      cond(n) = DONE
    end do
    do n = 1, NCHAIN
      sto(n) = DZERO
    end do
    call split_chain(cond, sto, determined_a, determined_b, nb, brow, bgid, &
                     ng, gid, determined, nlink, lgid, lrow)
    call merge_groups(nb, brow, bgid, ng, gid, determined, nlink, lgid, lrow)
    call check(error, determined(1) == 0 .and. determined(2) == 0, &
               "neither side should be determined")
  end subroutine test_merge_undetermined

  subroutine test_merge_no_link(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), dimension(:), allocatable :: brow, bgid, gid, determined
    integer(I4B), dimension(:), allocatable :: lgid, lrow
    integer(I4B) :: determined_a(3), determined_b(3), nb, ng, nlink, n
    real(DP) :: cond(NCHAIN - 1), sto(NCHAIN)

    ! -- no conductance across the split, so storage on one side does not
    !    determine the other
    do n = 1, NCHAIN - 1
      cond(n) = DONE
    end do
    cond(3) = DZERO
    do n = 1, NCHAIN
      sto(n) = DZERO
    end do
    sto(1) = 0.1_DP
    call split_chain(cond, sto, determined_a, determined_b, nb, brow, bgid, &
                     ng, gid, determined, nlink, lgid, lrow)
    call check(error, nb == 2 .and. nlink == 0, &
               "rows 3 and 4 should be boundary rows with no links")
    if (allocated(error)) return
    call merge_groups(nb, brow, bgid, ng, gid, determined, nlink, lgid, lrow)
    call check(error, determined(1) == 1 .and. determined(2) == 0, &
               "only the first side should be determined")
  end subroutine test_merge_no_link

  subroutine test_tracker_status(error)
    type(error_type), allocatable, intent(out) :: error
    type(IsolatedTrackerType) :: track
    integer(I4B) :: nnew

    ! -- row 1 found in every check, row 2 only in the first, row 3 only in
    !    the last
    call track%record(3, [1, 1, 0], [0, 0, 0], 1, 1, nnew)
    call check(error, nnew == 2, "two rows should be new in check 1")
    if (allocated(error)) return
    call track%record(3, [1, 0, 0], [0, 0, 0], 1, 2, nnew)
    call check(error, nnew == 0, "no rows should be new in check 2")
    if (allocated(error)) return
    call track%record(3, [1, 0, 1], [0, 0, 0], 2, 1, nnew)
    call check(error, nnew == 1, "one row should be new in check 3")
    if (allocated(error)) return
    call check(error, track%status(1) == 'PERSISTENT', &
               "row 1 should be persistent")
    if (allocated(error)) return
    call check(error, track%status(2) == 'RECOVERED', &
               "row 2 should have recovered")
    if (allocated(error)) return
    call check(error, track%status(3) == 'NOT RECHECKED', &
               "row 3 should not have been rechecked")
    if (allocated(error)) return
    call check(error, track%nfound(1) == 3 .and. track%nfound(2) == 1, &
               "rows 1 and 2 should be found in 3 checks and 1 check")
    if (allocated(error)) return
    call check(error, track%kper_first(3) == 2 .and. &
               track%kstp_last(1) == 1 .and. track%kper_last(1) == 2, &
               "first and last stress periods and time steps")
  end subroutine test_tracker_status

  subroutine test_tracker_cause(error)
    type(error_type), allocatable, intent(out) :: error
    type(IsolatedTrackerType) :: track
    integer(I4B) :: nnew

    ! -- nothing found: no arrays
    call track%record(2, [0, 0], [0, 0], 1, 1, nnew)
    call check(error,.not. allocated(track%nfound), &
               "nothing should be stored until a row is found")
    if (allocated(error)) return
    ! -- the cause is taken when a row is first found
    call track%record(2, [1, 1], [1, 0], 1, 2, nnew)
    call track%record(2, [1, 1], [0, 1], 1, 3, nnew)
    call check(error, track%cause(1) == 'DRY', "row 1 should be dry")
    if (allocated(error)) return
    call check(error, track%cause(2) == 'INPUT', &
               "row 2 should be cut off by the input")
  end subroutine test_tracker_cause

end module TestIsolatedCells
