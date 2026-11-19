!> @brief Find and track cells with an undetermined dependent variable
!!
!! A row is controlled when its sum is nonzero, from storage, a
!! dependent-variable boundary, or a constant value. A group of connected rows
!! with no controlled row has no unique solution.
!<
module IsolatedCellsModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO
  use SortModule, only: qsort, search_sorted
  use DisjointSetModule, only: DisjointSetType

  implicit none
  private
  public :: row_is_controlled
  public :: find_groups
  public :: boundary_links
  public :: merge_groups
  public :: IsolatedTrackerType

  !> a row is controlled when |row sum| exceeds this fraction of its largest
  !! term; conductance terms alone cancel to a few rounding errors
  real(DP), parameter :: CONTROL_TOL = 1.0e-12_DP

  !> @brief Record of the cells found with an undetermined dependent variable
  !<
  type :: IsolatedTrackerType
    integer(I4B) :: ncheck = 0 !< checks made
    integer(I4B), dimension(:), allocatable :: nfound !< checks that found each row
    integer(I4B), dimension(:), allocatable :: icheck_first !< check that first found each row
    integer(I4B), dimension(:), allocatable :: icheck_last !< check that last found each row
    integer(I4B), dimension(:), allocatable :: kper_first !< stress period first found
    integer(I4B), dimension(:), allocatable :: kstp_first !< time step first found
    integer(I4B), dimension(:), allocatable :: kper_last !< stress period last found
    integer(I4B), dimension(:), allocatable :: kstp_last !< time step last found
    integer(I4B), dimension(:), allocatable :: idry !< dry when first found (1)
  contains
    procedure :: record
    procedure :: cause
    procedure :: status
  end type IsolatedTrackerType

contains

  !> @brief Return 1 if the row sum is nonzero relative to the row's terms
  !<
  function row_is_controlled(nval, vals) result(row_control)
    integer(I4B), intent(in) :: nval !< number of terms in the row
    real(DP), dimension(nval), intent(in) :: vals !< row terms
    integer(I4B) :: row_control
    ! -- local
    integer(I4B) :: k
    real(DP) :: rsum, rmax

    rsum = DZERO
    rmax = DZERO
    do k = 1, nval
      rsum = rsum + vals(k)
      if (abs(vals(k)) > rmax) rmax = abs(vals(k))
    end do
    row_control = 0
    if (rmax > DZERO) then
      if (abs(rsum) > CONTROL_TOL * rmax) row_control = 1
    end if
  end function row_is_controlled

  !> @brief Group the local rows connected by nonzero terms
  !!
  !! Each row gets the global number of the first row in its group, which
  !! identifies the group across processes. group_determined is set for each
  !! group from its rows and is indexed by that number less the row offset.
  !<
  subroutine find_groups(nrow, roff, ia, ja, a, row_control, igid, &
                         group_determined)
    integer(I4B), intent(in) :: nrow !< number of local rows
    integer(I4B), intent(in) :: roff !< global number of the row before the first
    integer(I4B), dimension(nrow + 1), intent(in) :: ia !< CSR row pointers
    integer(I4B), dimension(*), intent(in) :: ja !< global column numbers
    real(DP), dimension(*), intent(in) :: a !< matrix terms
    integer(I4B), dimension(nrow), intent(in) :: row_control !< controlled (1)
    integer(I4B), dimension(nrow), intent(inout) :: igid !< group of each row
    integer(I4B), dimension(nrow), intent(inout) :: &
      group_determined !< group has a controlled row (1)
    ! -- local
    type(DisjointSetType) :: groups
    integer(I4B) :: n, k, m

    call groups%init(nrow)
    do n = 1, nrow
      do k = ia(n), ia(n + 1) - 1
        m = ja(k) - roff
        if (m < 1 .or. m > nrow .or. m == n) cycle
        if (a(k) == DZERO) cycle
        call groups%join(n, m)
      end do
    end do
    do n = 1, nrow
      group_determined(n) = 0
    end do
    do n = 1, nrow
      m = groups%find(n)
      igid(n) = m + roff
      if (row_control(n) /= 0) group_determined(m) = 1
    end do
    call groups%destroy()
  end subroutine find_groups

  !> @brief List the links from local groups to rows on other processes
  !!
  !! brow lists the local rows with any off-process column, whether or not the
  !! term is zero, so every linked row can be looked up; lgid and lrow list
  !! the nonzero terms as (local group, remote row).
  !<
  subroutine boundary_links(nrow, roff, ia, ja, a, igid, nbrow, brow, bgid, &
                            nlink, lgid, lrow)
    integer(I4B), intent(in) :: nrow !< number of local rows
    integer(I4B), intent(in) :: roff !< global number of the row before the first
    integer(I4B), dimension(nrow + 1), intent(in) :: ia !< CSR row pointers
    integer(I4B), dimension(*), intent(in) :: ja !< global column numbers
    real(DP), dimension(*), intent(in) :: a !< matrix terms
    integer(I4B), dimension(nrow), intent(in) :: igid !< group of each row
    integer(I4B), intent(out) :: nbrow !< number of boundary rows
    integer(I4B), dimension(:), allocatable, intent(out) :: brow !< global row
    integer(I4B), dimension(:), allocatable, intent(out) :: bgid !< its group
    integer(I4B), intent(out) :: nlink !< number of links
    integer(I4B), dimension(:), allocatable, intent(out) :: lgid !< local group
    integer(I4B), dimension(:), allocatable, intent(out) :: lrow !< remote row
    ! -- local
    integer(I4B) :: n, k, m, ipass
    logical :: isbnd

    ! -- count on the first pass, fill on the second
    do ipass = 1, 2
      nbrow = 0
      nlink = 0
      do n = 1, nrow
        isbnd = .false.
        do k = ia(n), ia(n + 1) - 1
          m = ja(k) - roff
          if (m >= 1 .and. m <= nrow) cycle
          isbnd = .true.
          if (a(k) == DZERO) cycle
          nlink = nlink + 1
          if (ipass == 2) then
            lgid(nlink) = igid(n)
            lrow(nlink) = ja(k)
          end if
        end do
        if (isbnd) then
          nbrow = nbrow + 1
          if (ipass == 2) then
            brow(nbrow) = n + roff
            bgid(nbrow) = igid(n)
          end if
        end if
      end do
      if (ipass == 1) then
        allocate (brow(nbrow), bgid(nbrow), lgid(nlink), lrow(nlink))
      end if
    end do
  end subroutine boundary_links

  !> @brief Merge groups linked across processes
  !!
  !! The arrays hold the boundary data from every process. On return,
  !! group_determined is 1 for each listed group that belongs to a merged group
  !! with a controlled row.
  !<
  subroutine merge_groups(nb, brow, bgid, ng, gid, group_determined, nlink, &
                          lgid, lrow)
    integer(I4B), intent(in) :: nb !< number of boundary rows
    integer(I4B), dimension(nb), intent(in) :: brow !< global boundary rows
    integer(I4B), dimension(nb), intent(in) :: bgid !< group of each row
    integer(I4B), intent(in) :: ng !< number of groups
    integer(I4B), dimension(ng), intent(in) :: gid !< groups, no duplicates
    integer(I4B), dimension(ng), intent(inout) :: &
      group_determined !< group has a controlled row (1)
    integer(I4B), intent(in) :: nlink !< number of links
    integer(I4B), dimension(nlink), intent(in) :: lgid !< group on one side
    integer(I4B), dimension(nlink), intent(in) :: lrow !< row on the other side
    ! -- local
    type(DisjointSetType) :: groups
    integer(I4B), dimension(:), allocatable :: gsort, gidx, rsort, ridx
    integer(I4B) :: i, ga, gb, ir

    call groups%init(ng)
    allocate (gsort(ng), gidx(ng), rsort(nb), ridx(nb))
    do i = 1, ng
      gsort(i) = gid(i)
      gidx(i) = i
    end do
    do i = 1, nb
      rsort(i) = brow(i)
      ridx(i) = i
    end do
    call qsort(gidx, gsort)
    call qsort(ridx, rsort)
    do i = 1, nlink
      ir = search_sorted(rsort, ridx, lrow(i))
      if (ir == 0) cycle
      ga = search_sorted(gsort, gidx, lgid(i))
      gb = search_sorted(gsort, gidx, bgid(ir))
      if (ga == 0 .or. gb == 0) cycle
      call groups%join(ga, gb)
    end do
    ! -- a merged group is determined if any member is
    do i = 1, ng
      if (group_determined(i) /= 0) group_determined(groups%find(i)) = 1
    end do
    do i = 1, ng
      group_determined(i) = group_determined(groups%find(i))
    end do
    call groups%destroy()
    deallocate (gsort, gidx, rsort, ridx)
  end subroutine merge_groups

  !> @brief Record the rows found in one check
  !!
  !! nnew returns the number of rows found for the first time.
  !<
  subroutine record(this, nrow, ifound, idry, kper, kstp, nnew)
    class(IsolatedTrackerType) :: this !< tracker
    integer(I4B), intent(in) :: nrow !< number of local rows
    integer(I4B), dimension(nrow), intent(in) :: ifound !< rows found (1)
    integer(I4B), dimension(nrow), intent(in) :: idry !< dry rows (1)
    integer(I4B), intent(in) :: kper !< stress period
    integer(I4B), intent(in) :: kstp !< time step
    integer(I4B), intent(out) :: nnew !< rows found for the first time
    ! -- local
    integer(I4B) :: n

    this%ncheck = this%ncheck + 1
    nnew = 0
    do n = 1, nrow
      if (ifound(n) == 0) cycle
      if (.not. allocated(this%nfound)) call allocate_tracker(this, nrow)
      if (this%nfound(n) == 0) then
        nnew = nnew + 1
        this%icheck_first(n) = this%ncheck
        this%kper_first(n) = kper
        this%kstp_first(n) = kstp
        this%idry(n) = idry(n)
      end if
      this%nfound(n) = this%nfound(n) + 1
      this%icheck_last(n) = this%ncheck
      this%kper_last(n) = kper
      this%kstp_last(n) = kstp
    end do
  end subroutine record

  !> @brief Return why a row was found: dry, or cut off by the input
  !<
  function cause(this, n) result(text)
    class(IsolatedTrackerType) :: this !< tracker
    integer(I4B), intent(in) :: n !< row
    character(len=5) :: text

    if (this%idry(n) /= 0) then
      text = 'DRY'
    else
      text = 'INPUT'
    end if
  end function cause

  !> @brief Return whether a row was still found by the last check
  !<
  function status(this, n) result(text)
    class(IsolatedTrackerType) :: this !< tracker
    integer(I4B), intent(in) :: n !< row
    character(len=13) :: text

    if (this%icheck_first(n) == this%ncheck) then
      text = 'NOT RECHECKED'
    else if (this%icheck_last(n) == this%ncheck) then
      text = 'PERSISTENT'
    else
      text = 'RECOVERED'
    end if
  end function status

  !> @brief Allocate and zero the tracker arrays
  !<
  subroutine allocate_tracker(this, nrow)
    type(IsolatedTrackerType), intent(inout) :: this !< tracker
    integer(I4B), intent(in) :: nrow !< number of local rows
    ! -- local
    integer(I4B) :: n

    allocate (this%nfound(nrow), this%icheck_first(nrow), &
              this%icheck_last(nrow), this%kper_first(nrow), &
              this%kstp_first(nrow), this%kper_last(nrow), &
              this%kstp_last(nrow), this%idry(nrow))
    do n = 1, nrow
      this%nfound(n) = 0
      this%icheck_first(n) = 0
      this%icheck_last(n) = 0
      this%kper_first(n) = 0
      this%kstp_first(n) = 0
      this%kper_last(n) = 0
      this%kstp_last(n) = 0
      this%idry(n) = 0
    end do
  end subroutine allocate_tracker

end module IsolatedCellsModule
