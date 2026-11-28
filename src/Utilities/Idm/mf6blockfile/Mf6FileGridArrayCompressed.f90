!> @brief Grid array period block loader, compressed mode
!!
!! For packages whose dfn defines MAXBOUND (e.g. CHDG, WELG, DRNG, RIVG,
!! GHBG). The PERIOD block defines boundary cells; the package reads them
!! as a list (NODEULIST, NBOUND) with arrays indexed by boundary.
!!
!! Period block rules:
!!  - Boundary cells are the active cells (IDOMAIN > 0) where every main
!!    (non-aux) array given in the block has a defined (non-DNODATA)
!!    value. Values at other cells are ignored. Main arrays with defined
!!    values but no common active cell are an input error.
!!  - An array not given in the block carries forward per cell; a new
!!    boundary cell is then an input error. A main array must be given
!!    before it can carry forward; an aux array never given is zero.
!!  - An aux array given in the block must be defined at every boundary
!!    cell. A block with only aux arrays keeps the boundary cells.
!!  - An empty block removes all boundaries and time-array series links.
!!  - An array given as a time-array series is linked to the series and
!!    updated each time step until the array is given again or an empty
!!    block is read. A linked main array sets boundary cells from the
!!    series array at the period start.
!<
module GridArrayCompressedLoadModule

  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: DNODATA, DZERO
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename, count_errors
  use MemoryManagerModule, only: mem_allocate, mem_setptr
  use TimeArraySeriesModule, only: TimeArraySeriesType
  use GridArrayLoadModule, only: GridArrayLoadType

  implicit none
  private
  public :: GridArrayCompressedLoadType

  !> @brief storage state of an array read in the current period block
  !<
  enum, bind(C)
    enumerator :: STORED = 1 !< values stored at block boundary cells
    enumerator :: DEFERRED = 2 !< defined values kept until cells are known
  end enum

  !> @brief Block storage of a period block array
  !<
  type :: BlockArrayType
    integer(I4B) :: state = 0 !< block storage state
    integer(I4B) :: ndef = 0 !< block number of defined active cells
    integer(I4B), dimension(:), allocatable :: cells !< DEFERRED user nodes (ascending)
    real(DP), dimension(:), allocatable :: vals !< DEFERRED values
  end type BlockArrayType

  !> @brief Grid array dynamic loader, compressed mode
  !<
  type, extends(GridArrayLoadType) :: GridArrayCompressedLoadType
    integer(I4B), dimension(:), pointer, contiguous :: nodeulist => null() !< boundary cells (user nodes)
    type(BlockArrayType), dimension(:), allocatable :: barr !< block storage of each slot
    integer(I4B), dimension(:), allocatable :: bcells !< block boundary cells
    integer(I4B) :: nbcells = 0 !< number of block boundary cells
    logical(LGP) :: bcells_set = .false. !< block boundary cells established
  contains
    procedure :: mode_init
    procedure :: array_store
    procedure :: block_end
    procedure :: tas_links_update
    procedure :: tas_start_array
    procedure :: reset
    procedure :: slot_clear
    procedure :: destroy
    procedure :: defer
    procedure :: bcells_narrow
    procedure :: bcells_compact
    procedure :: bcells_intersect
    procedure :: store_values
    procedure :: store_deferred
    procedure :: apply_deferred
    procedure :: carry_forward
  end type GridArrayCompressedLoadType

contains

  !> @brief Allocate the boundary cell list and block storage
  !<
  subroutine mode_init(this)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    call mem_allocate(this%nodeulist, this%ctx%maxbound, 'NODEULIST', &
                      this%mf6_input%mempath)
    allocate (this%bcells(this%ctx%maxbound))
    allocate (this%barr(size(this%slots)))
  end subroutine mode_init

  subroutine reset(this)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    ! nbound and nodeulist are kept until the block is finalized so
    ! arrays not given in the block can carry forward
    this%nbcells = 0
    this%bcells_set = .false.
    call this%GridArrayLoadType%reset()
  end subroutine reset

  !> @brief Clear the block state of a slot
  !<
  subroutine slot_clear(this, islot)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: islot
    call this%GridArrayLoadType%slot_clear(islot)
    if (.not. allocated(this%barr)) return
    if (allocated(this%barr(islot)%cells)) deallocate (this%barr(islot)%cells)
    if (allocated(this%barr(islot)%vals)) deallocate (this%barr(islot)%vals)
    this%barr(islot)%state = 0
    this%barr(islot)%ndef = 0
  end subroutine slot_clear

  subroutine destroy(this)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    call this%GridArrayLoadType%destroy()
    nullify (this%nodeulist)
    if (allocated(this%barr)) deallocate (this%barr)
    if (allocated(this%bcells)) deallocate (this%bcells)
  end subroutine destroy

  !> @brief Store an array read in the period block
  !!
  !! A main array establishes or narrows the block boundary cells.
  !! Values are stored at the boundary cells once they are known; until
  !! then the defined values are deferred.
  !<
  subroutine array_store(this, islot, nodes)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: islot
    real(DP), dimension(:), pointer, contiguous, intent(in) :: nodes
    integer(I4B) :: iaux, n, nd

    iaux = this%slots(islot)%iaux

    ! count defined active cells
    nd = 0
    do n = 1, size(nodes)
      if (nodes(n) /= DNODATA .and. this%active(n)) nd = nd + 1
    end do
    this%barr(islot)%ndef = nd

    if (this%bcells_set) then
      ! narrow boundary cells to the overlap (main) and store values
      if (iaux == 0) call this%bcells_narrow(nodes)
      call this%store_values(islot, nodes)
    else if (iaux == 0 .and. nd <= this%ctx%maxbound) then
      ! establish boundary cells, then apply deferred arrays
      this%nbcells = 0
      do n = 1, size(nodes)
        if (nodes(n) == DNODATA .or. .not. this%active(n)) cycle
        this%nbcells = this%nbcells + 1
        this%bcells(this%nbcells) = n
      end do
      this%bcells_set = .true.
      call this%store_values(islot, nodes)
      call this%apply_deferred()
    else
      ! an aux array, or a main array with more defined cells than
      ! maxbound: defer until boundary cells are known
      call this%defer(islot, nodes, nd)
    end if
  end subroutine array_store

  !> @brief Series array at the period start, added like a read array
  !!
  !! A main array given as a series sets or narrows the boundary cells.
  !<
  function tas_start_array(this, tasptr) result(nodes)
    use TdisModule, only: totim
    class(GridArrayCompressedLoadType), intent(inout) :: this
    type(TimeArraySeriesType), pointer, intent(in) :: tasptr
    real(DP), dimension(:), pointer, contiguous :: nodes
    ! totim before the first step of the period
    nodes => tasptr%GetArrayAtTime(totim)
  end function tas_start_array

  !> @brief Update time-array series links after a period block
  !!
  !! An array given in the block with values ends its link; an array
  !! given as a series is linked. Links for arrays not given are kept. An
  !! empty block removes all links. Links update the boundary cells.
  !<
  subroutine tas_links_update(this, inunit)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: inunit
    integer(I4B) :: islot

    if (.not. this%tas_active) return
    if (this%nreads == 0) then
      do islot = 1, size(this%slots)
        this%slots(islot)%tasname = ''
      end do
    end if
    call this%make_links(inunit, this%nodeulist, .false.)
  end subroutine tas_links_update

  !> @brief Keep the defined active values of an array until the block
  !! boundary cells are known
  !<
  subroutine defer(this, islot, nodes, nd)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: islot
    real(DP), dimension(:), contiguous, intent(in) :: nodes
    integer(I4B), intent(in) :: nd !< number of defined active cells
    integer(I4B) :: n, k

    allocate (this%barr(islot)%cells(nd))
    allocate (this%barr(islot)%vals(nd))
    k = 0
    do n = 1, size(nodes)
      if (nodes(n) == DNODATA .or. .not. this%active(n)) cycle
      k = k + 1
      this%barr(islot)%cells(k) = n
      this%barr(islot)%vals(k) = nodes(n)
    end do
    this%barr(islot)%state = DEFERRED
  end subroutine defer

  !> @brief Narrow block boundary cells to cells defined in nodes
  !<
  subroutine bcells_narrow(this, nodes)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    real(DP), dimension(:), contiguous, intent(in) :: nodes
    logical(LGP), dimension(:), allocatable :: keep
    integer(I4B) :: k

    allocate (keep(this%nbcells))
    do k = 1, this%nbcells
      keep(k) = (nodes(this%bcells(k)) /= DNODATA)
    end do
    call this%bcells_compact(keep)
  end subroutine bcells_narrow

  !> @brief Keep the block boundary cells where keep is true
  !!
  !! Values already stored for the block are compacted to the
  !! remaining cells.
  !<
  subroutine bcells_compact(this, keep)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    logical(LGP), dimension(:), intent(in) :: keep !< one per boundary cell
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B) :: i, islot, iaux, nkeep

    nkeep = count(keep)
    if (nkeep == this%nbcells) return

    ! compact values stored for the block
    do i = 1, this%nreads
      islot = this%blk(i)
      if (this%barr(islot)%state /= STORED) cycle
      iaux = this%slots(islot)%iaux
      if (iaux == 0) then
        call mem_setptr(dbl1d, this%slots(islot)%mf6varname, &
                        this%mf6_input%mempath)
        dbl1d(1:nkeep) = pack(dbl1d(1:this%nbcells), keep)
      else
        call mem_setptr(dbl2d, this%slots(islot)%mf6varname, &
                        this%mf6_input%mempath)
        dbl2d(iaux, 1:nkeep) = pack(dbl2d(iaux, 1:this%nbcells), keep)
      end if
    end do

    ! compact boundary cells
    this%bcells(1:nkeep) = pack(this%bcells(1:this%nbcells), keep)
    this%nbcells = nkeep
  end subroutine bcells_compact

  !> @brief Set block boundary cells to the cells defined in every
  !! deferred main array
  !<
  subroutine bcells_intersect(this)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    integer(I4B), dimension(:), allocatable :: cells, idx
    integer(I4B) :: islot, imin

    ! start from the main array with the fewest defined cells
    imin = 0
    do islot = 1, this%nmain
      if (.not. this%slots(islot)%inblock) cycle
      if (imin == 0) then
        imin = islot
      else if (this%barr(islot)%ndef < this%barr(imin)%ndef) then
        imin = islot
      end if
    end do
    cells = this%barr(imin)%cells

    ! keep the cells defined in every other main array
    do islot = 1, this%nmain
      if (.not. this%slots(islot)%inblock .or. islot == imin) cycle
      allocate (idx(size(cells)))
      call merge_index(this%barr(islot)%cells, cells, idx)
      cells = pack(cells, idx > 0)
      deallocate (idx)
    end do

    if (size(cells) > this%ctx%maxbound) then
      write (errmsg, '(a,i0,a)') &
        'Input error: number of defined (non-DNODATA) active cells &
        &exceeds MAXBOUND=', this%ctx%maxbound, '.'
      call store_error(errmsg)
      call store_error_filename(this%input_name)
      return
    end if
    this%nbcells = size(cells)
    this%bcells(1:this%nbcells) = cells
  end subroutine bcells_intersect

  !> @brief Store block array values at block boundary cells
  !<
  subroutine store_values(this, islot, nodes)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: islot
    real(DP), dimension(:), contiguous, intent(in) :: nodes
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B) :: k, iaux

    iaux = this%slots(islot)%iaux
    if (iaux == 0) then
      call mem_setptr(dbl1d, this%slots(islot)%mf6varname, &
                      this%mf6_input%mempath)
      do k = 1, this%nbcells
        dbl1d(k) = nodes(this%bcells(k))
      end do
    else
      call mem_setptr(dbl2d, this%slots(islot)%mf6varname, &
                      this%mf6_input%mempath)
      do k = 1, this%nbcells
        dbl2d(iaux, k) = nodes(this%bcells(k))
      end do
    end if
    this%barr(islot)%state = STORED
  end subroutine store_values

  !> @brief Store deferred values at block boundary cells
  !!
  !! Boundary cells without a value are set to DNODATA.
  !<
  subroutine store_deferred(this, islot)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: islot
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B), dimension(:), allocatable :: idx
    real(DP) :: v
    integer(I4B) :: k, iaux

    iaux = this%slots(islot)%iaux
    if (iaux == 0) then
      call mem_setptr(dbl1d, this%slots(islot)%mf6varname, &
                      this%mf6_input%mempath)
    else
      call mem_setptr(dbl2d, this%slots(islot)%mf6varname, &
                      this%mf6_input%mempath)
    end if
    allocate (idx(this%nbcells))
    call merge_index(this%barr(islot)%cells, this%bcells(1:this%nbcells), idx)
    do k = 1, this%nbcells
      v = DNODATA
      if (idx(k) > 0) v = this%barr(islot)%vals(idx(k))
      if (iaux == 0) then
        dbl1d(k) = v
      else
        dbl2d(iaux, k) = v
      end if
    end do
    deallocate (this%barr(islot)%cells)
    deallocate (this%barr(islot)%vals)
    this%barr(islot)%state = STORED
  end subroutine store_deferred

  !> @brief Apply deferred block arrays once boundary cells are known
  !!
  !! Deferred main arrays narrow the boundary cells first (unless narrow
  !! is false), then all deferred arrays are stored at the remaining
  !! cells.
  !<
  subroutine apply_deferred(this, narrow)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    logical(LGP), optional, intent(in) :: narrow !< narrow by deferred main arrays
    integer(I4B), dimension(:), allocatable :: idx
    integer(I4B) :: i, islot
    logical(LGP) :: lnarrow

    lnarrow = .true.
    if (present(narrow)) lnarrow = narrow

    if (lnarrow) then
      do i = 1, this%nreads
        islot = this%blk(i)
        if (this%barr(islot)%state /= DEFERRED .or. &
            this%slots(islot)%iaux /= 0) cycle
        allocate (idx(this%nbcells))
        call merge_index(this%barr(islot)%cells, &
                         this%bcells(1:this%nbcells), idx)
        call this%bcells_compact(idx > 0)
        deallocate (idx)
      end do
    end if

    do i = 1, this%nreads
      islot = this%blk(i)
      if (this%barr(islot)%state == DEFERRED) call this%store_deferred(islot)
    end do
  end subroutine apply_deferred

  !> @brief End a compressed mode period block
  !!
  !! Sets the boundary cells, checks main and aux arrays at them, and
  !! carries forward arrays not given in the block.
  !<
  subroutine block_end(this)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B) :: i, k, islot, nmain_defined, nmain_empty, iaux, nerr
    logical(LGP) :: main_read, lnarrow

    ! errors stored before this block
    nerr = count_errors()

    ! empty block removes all boundaries
    if (this%nreads == 0) then
      this%ctx%nbound = 0
      return
    end if

    ! establish boundary cells if not set by a main array
    if (.not. this%bcells_set) then
      main_read = any(this%slots(1:this%nmain)%inblock)
      this%nbcells = 0
      lnarrow = .true.
      if (main_read) then
        ! all main arrays deferred: cells defined in every main array
        lnarrow = .false.
        call this%bcells_intersect()
      else
        ! aux only block keeps current boundary cells
        this%nbcells = this%ctx%nbound
        this%bcells(1:this%nbcells) = this%nodeulist(1:this%nbcells)
      end if
      this%bcells_set = .true.
      call this%apply_deferred(lnarrow)
    end if

    ! main arrays with defined values must share an active cell; name
    ! any main array with no defined values
    nmain_defined = 0
    nmain_empty = 0
    do islot = 1, this%nmain
      if (.not. this%slots(islot)%inblock) cycle
      if (this%barr(islot)%ndef > 0) then
        nmain_defined = nmain_defined + 1
      else
        nmain_empty = nmain_empty + 1
      end if
    end do
    if (this%nbcells == 0 .and. nmain_defined > 0) then
      if (nmain_empty > 0) then
        do i = 1, this%nreads
          islot = this%blk(i)
          if (this%slots(islot)%iaux /= 0) cycle
          if (this%barr(islot)%ndef == 0) then
            write (errmsg, '(a)') &
              'Input error: '//trim(this%slots(islot)%tagname)//' has no &
              &defined (non-DNODATA) values in active cells.'
            call store_error(errmsg)
          end if
        end do
      else
        write (errmsg, '(a)') &
          'Input error: main arrays have no defined (non-DNODATA) values in &
          &the same active cells.'
        call store_error(errmsg)
      end if
    end if

    ! aux arrays given in the block must be defined at boundary cells
    do i = 1, this%nreads
      islot = this%blk(i)
      iaux = this%slots(islot)%iaux
      if (iaux == 0) cycle
      call mem_setptr(dbl2d, this%slots(islot)%mf6varname, &
                      this%mf6_input%mempath)
      do k = 1, this%nbcells
        if (dbl2d(iaux, k) == DNODATA) then
          write (errmsg, '(a,i0,a)') &
            'Input error: auxiliary variable '// &
            trim(this%slots(islot)%tagname)// &
            ' is not defined (DNODATA) at boundary cell (user node ', &
            this%bcells(k), ').'
          call store_error(errmsg)
          exit
        end if
      end do
    end do

    ! carry forward arrays not given in the block
    do islot = 1, size(this%slots)
      if (this%slots(islot)%inblock) cycle
      iaux = this%slots(islot)%iaux
      if (this%slots(islot)%given) then
        call this%carry_forward(islot)
      else if (iaux == 0) then
        ! a main array must be given before it can carry forward
        if (this%nbcells > 0) then
          write (errmsg, '(a)') &
            'Input error: required array '// &
            trim(this%slots(islot)%tagname)//' has not been specified.'
          call store_error(errmsg)
        end if
      else
        ! aux never given is zero
        call mem_setptr(dbl2d, this%slots(islot)%mf6varname, &
                        this%mf6_input%mempath)
        dbl2d(iaux, 1:this%nbcells) = DZERO
      end if
    end do

    if (count_errors() > nerr) then
      call store_error_filename(this%input_name)
    end if

    ! set period boundary cells
    this%nodeulist(1:this%nbcells) = this%bcells(1:this%nbcells)
    this%ctx%nbound = this%nbcells
  end subroutine block_end

  !> @brief Carry forward values for an array not given in the block
  !!
  !! Values move from the previous period boundary cells to the block
  !! boundary cells. A new boundary cell without a previous value is an
  !! input error.
  !<
  subroutine carry_forward(this, islot)
    class(GridArrayCompressedLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: islot
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    real(DP), dimension(:), allocatable :: prev
    integer(I4B), dimension(:), allocatable :: idx
    integer(I4B) :: k, nprev, iaux

    iaux = this%slots(islot)%iaux
    nprev = this%ctx%nbound
    if (iaux == 0) then
      call mem_setptr(dbl1d, this%slots(islot)%mf6varname, &
                      this%mf6_input%mempath)
      prev = dbl1d(1:nprev)
    else
      call mem_setptr(dbl2d, this%slots(islot)%mf6varname, &
                      this%mf6_input%mempath)
      prev = dbl2d(iaux, 1:nprev)
    end if

    allocate (idx(this%nbcells))
    call merge_index(this%nodeulist(1:nprev), this%bcells(1:this%nbcells), idx)
    do k = 1, this%nbcells
      if (idx(k) == 0) then
        write (errmsg, '(a,i0,a)') &
          'Input error: '//trim(this%slots(islot)%tagname)// &
          ' must be specified for new boundary cell (user node ', &
          this%bcells(k), ').'
        call store_error(errmsg)
        exit
      end if
      if (iaux == 0) then
        dbl1d(k) = prev(idx(k))
      else
        dbl2d(iaux, k) = prev(idx(k))
      end if
    end do
  end subroutine carry_forward

  !> @brief Positions of cells in a cell list
  !!
  !! Both lists are in ascending order. idx(k) is the position of dst(k)
  !! in src, or 0 if dst(k) is not in src.
  !<
  subroutine merge_index(src, dst, idx)
    integer(I4B), dimension(:), intent(in) :: src
    integer(I4B), dimension(:), intent(in) :: dst
    integer(I4B), dimension(:), intent(out) :: idx
    integer(I4B) :: j, k

    j = 1
    do k = 1, size(dst)
      idx(k) = 0
      do while (j <= size(src))
        if (src(j) >= dst(k)) exit
        j = j + 1
      end do
      if (j <= size(src)) then
        if (src(j) == dst(k)) idx(k) = j
      end if
    end do
  end subroutine merge_index

end module GridArrayCompressedLoadModule
