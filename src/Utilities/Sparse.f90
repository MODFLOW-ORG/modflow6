module SparseModule
!cdl (8/12/2013) sparse matrix module for managing the structure
!of a matrix.  Module uses FORTRAN 2003 extensions to manage
!the data structures in an object oriented fashion.

  use KindModule, only: DP, I4B, LGP
  implicit none

  type rowtype
    integer(I4B) :: nnz ! number of nonzero entries in the row
    integer(I4B), allocatable, dimension(:) :: icolarray ! array of column numbers
  end type rowtype

  type, public :: sparsematrix
    integer(I4B) :: offset !< global offset for first row in this matrix (default = 0)
    integer(I4B) :: nrow !< number of rows in the matrix
    integer(I4B) :: ncol !< number of columns in the matrix
    integer(I4B) :: nnz !< number of nonzero matrix entries
    integer(I4B) :: nnz_od !< number of off-diagonal nonzero matrix entries
    type(rowtype), allocatable, dimension(:) :: row !< one rowtype for each matrix row
  contains
    generic :: init => initialize, initializefixed
    procedure :: addconnection
    procedure :: filliaja
    procedure :: sort
    procedure :: destroy

    procedure, private :: initializefixed
    procedure, private :: initialize
  end type sparsematrix

contains

  subroutine destroy(this)
    class(sparsematrix), intent(inout) :: this
    deallocate (this%row)
  end subroutine destroy

  subroutine initialize(this, nrow, ncol, rowmaxnnz)
    !initial the sparse matrix.  This subroutine
    !acts a method for a sparse matrix by initializing
    !the row data.  It presently requires one maximum
    !value for all rows, however, this can be changed
    !to a vector of maximum values with one value for
    !each row.
    ! -- dummy
    class(sparsematrix), intent(inout) :: this
    integer(I4B), intent(in) :: nrow, ncol
    integer(I4B), intent(in), dimension(nrow) :: rowmaxnnz
    ! -- local
    integer(I4B) :: i
    ! -- code
    this%offset = 0
    this%nrow = nrow
    this%ncol = ncol
    this%nnz = 0
    this%nnz_od = 0
    allocate (this%row(nrow))
    do i = 1, nrow
      allocate (this%row(i)%icolarray(rowmaxnnz(i)))
      this%row(i)%icolarray = 0
      this%row(i)%nnz = 0
    end do
  end subroutine initialize

  ! overload
  subroutine initializefixed(this, nrow, ncol, maxnnz)
    implicit none
    class(sparsematrix), intent(inout) :: this
    integer(I4B), intent(in) :: nrow, ncol
    integer(I4B), intent(in) :: maxnnz
    ! local
    integer(I4B), dimension(:), allocatable :: rowmaxnnz
    integer(I4B) :: i

    allocate (rowmaxnnz(nrow))

    do i = 1, nrow
      rowmaxnnz(i) = maxnnz
    end do

    call this%initialize(nrow, ncol, rowmaxnnz)
    deallocate (rowmaxnnz)

  end subroutine initializefixed

  subroutine filliaja(this, ia, ja, ierror, sort)
    !allocate and fill the ia and ja arrays using information
    !from the sparsematrix.
    !ierror is returned as:
    !  0 if no error
    !  1 if ia is not the correct size
    !  2 if ja is not the correct size
    !  3 if both ia and ja are not correct size
    ! -- dummy
    class(sparsematrix), intent(inout) :: this
    integer(I4B), dimension(:), intent(inout) :: ia, ja
    integer(I4B), intent(inout) :: ierror
    logical, intent(in), optional :: sort
    ! -- local
    logical :: sortja
    integer(I4B) :: i, j, ipos
    ! -- code
    !
    ! -- process optional dummy variables
    if (present(sort)) then
      sortja = sort
    else
      sortja = .FALSE.
    end if
    !
    ! -- initialize error variable
    ierror = 0
    !
    ! -- check for error conditions
    if (ubound(ia, dim=1) /= this%nrow + 1) then
      ierror = 1
    end if
    if (ubound(ja, dim=1) /= this%nnz) then
      ierror = ierror + 2
    end if
    if (ierror /= 0) then
      return
    end if
    !
    ! -- sort this
    if (sortja) then
      call this%sort()
    end if
    !
    ! -- fill ia and ja
    ipos = 1
    ia(1) = ipos
    do i = 1, this%nrow
      do j = 1, this%row(i)%nnz
        ja(ipos) = this%row(i)%icolarray(j)
        ipos = ipos + 1
      end do
      ia(i + 1) = ipos
    end do
  end subroutine filliaja

  subroutine addconnection(this, i, j, inodup, iaddop)
    !add a connection to the sparsematrix.  if inodup
    !(for no duplicates) is 1, then j is added only
    !if it is unique.
    ! -- dummy
    class(sparsematrix), intent(inout) :: this
    integer(I4B), intent(in) :: i, j, inodup
    integer(I4B), optional, intent(inout) :: iaddop
    ! -- local
    integer(I4B) :: irow_local
    integer(I4B) :: iadded
    ! -- code
    !
    ! -- when distributed system, reduce row numbers to local range
    irow_local = i - this%offset
    !
    call insert(j, this%row(irow_local), inodup, iadded)
    this%nnz = this%nnz + iadded
    if (j < this%offset + 1 .or. j > this%offset + this%nrow) then
      ! count the off-diagonal entries separately
      this%nnz_od = this%nnz_od + iadded
    end if
    if (present(iaddop)) iaddop = iadded
  end subroutine addconnection

  subroutine insert(j, thisrow, inodup, iadded)
    !insert j into thisrow (for row i)
    !inodup=1 means do not include duplicate connections
    !iadded is 1 if a new entry was added (meaning that nnz for the row was increased)
    !iadded is 0 if duplicate and not added.  Used to track total number of connections
    ! -- dummy
    integer(I4B), intent(in) :: j, inodup
    type(rowtype), intent(inout) :: thisrow
    integer(I4B), allocatable, dimension(:) :: iwk
    integer(I4B), intent(inout) :: iadded
    ! -- local
    integer(I4B) :: jj, maxnnz
    ! -- code
    iadded = 0
    maxnnz = ubound(thisrow%icolarray, dim=1)
    if (thisrow%icolarray(1) == 0) then
      thisrow%icolarray(1) = j
      thisrow%nnz = thisrow%nnz + 1
      iadded = 1
      return
    end if
    if (thisrow%nnz == maxnnz) then
      ! -- increase size of the row
      allocate (iwk(thisrow%nnz))
      iwk = thisrow%icolarray
      deallocate (thisrow%icolarray)
      ! -- Specify how to increase the size of the icolarray.  Adding 1
      !    will be most memory conservative, but may be a little slower
      !    due to frequent allocate/deallocate.  Another option would be
      !    to double the size: maxnnz=maxnnz*2
      maxnnz = maxnnz + 1
      allocate (thisrow%icolarray(maxnnz))
      thisrow%icolarray(1:thisrow%nnz) = iwk(1:thisrow%nnz)
      thisrow%icolarray(thisrow%nnz + 1:maxnnz) = 0
    end if
    if (inodup == 1) then
      do jj = 1, thisrow%nnz
        if (thisrow%icolarray(jj) == j) then
          return
        end if
      end do
    end if
    !
    ! -- add the connection to end
    thisrow%nnz = thisrow%nnz + 1
    thisrow%icolarray(thisrow%nnz) = j
    iadded = 1
  end subroutine insert

  subroutine sort(this, with_csr)
    !sort the icolarray for each row, but do not include
    !the diagonal position in the sort so that it stays in front
    ! -- dummy
    class(sparsematrix), intent(inout) :: this
    logical(LGP), optional :: with_csr
    ! -- local
    integer(I4B) :: i, nval, start_idx
    ! -- code
    start_idx = 2
    if (present(with_csr)) then
      if (with_csr) then
        ! CSR: don't put diagonal up front
        start_idx = 1
      end if
    end if

    do i = 1, this%nrow
      nval = this%row(i)%nnz
      call sortintarray(nval - start_idx + 1, &
                        this%row(i)%icolarray(start_idx:nval))
    end do
  end subroutine sort

  subroutine sortintarray(nval, iarray)
    !simple subroutine for sorting an array
    !in place.  It is not the fastest sort function
    !but should suffice for relatively short nodelists.
    ! -- dummy
    integer(I4B), intent(in) :: nval
    integer(I4B), intent(inout), dimension(nval) :: iarray
    ! -- local
    integer(I4B) :: i, j, itemp
    ! -- code
    do i = 1, nval - 1
      do j = i + 1, nval
        if (iarray(i) > iarray(j)) then
          itemp = iarray(j)
          iarray(j) = iarray(i)
          iarray(i) = itemp
        end if
      end do
    end do
  end subroutine sortintarray

  subroutine csr_diagsum(ia, flowja)
    !Add up the off diagonal terms and put the sum in the
    !diagonal position
    ! -- dummy
    integer(I4B), dimension(:), contiguous :: ia
    real(DP), dimension(:), contiguous :: flowja
    ! -- local
    integer(I4B) :: nodes
    integer(I4B) :: n
    integer(I4B) :: iposdiag
    integer(I4B) :: ipos
    ! -- code
    nodes = size(ia) - 1
    do n = 1, nodes
      iposdiag = ia(n)
      do ipos = ia(n) + 1, ia(n + 1) - 1
        flowja(iposdiag) = flowja(iposdiag) + flowja(ipos)
      end do
    end do
  end subroutine csr_diagsum

end module SparseModule
