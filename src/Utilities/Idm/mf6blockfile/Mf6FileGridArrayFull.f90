!> @brief Grid array period block loader, full-grid mode
!!
!! For packages whose dfn does not define MAXBOUND (e.g. TVKG). The
!! PERIOD block changes values on a fixed grid; the package reads arrays
!! indexed by user node.
!!
!! Period block rules:
!!  - Arrays are stored by user node, DNODATA where never set.
!!  - A defined (non-DNODATA) value at an active cell replaces the stored
!!    value; DNODATA and inactive cells leave it unchanged.
!!  - An array not given in the block, and an empty block, change nothing.
!!  - An array given as a time-array series is linked and updated each
!!    time step at active cells until the array is given again; a series
!!    DNODATA value leaves the cell unchanged.
!!
!! Array values are stamped when written (see IdmParamMetaModule):
!! RS_VALUES when given in a block, RS_SERIES when written by a series.
!! Each array has a stamp column, named by its input name for a main
!! array and by its auxiliary variable name for an aux array.
!<
module GridArrayFullLoadModule

  use KindModule, only: I4B, DP
  use ConstantsModule, only: LENVARNAME, DNODATA
  use MemoryManagerModule, only: mem_allocate, mem_setptr
  use IdmParamMetaModule, only: ParamStampsType, RS_VALUES, RS_SERIES
  use GridArrayLoadModule, only: GridArrayLoadType

  implicit none
  private
  public :: GridArrayFullLoadType

  !> @brief Grid array dynamic loader, full-grid mode
  !<
  type, extends(GridArrayLoadType) :: GridArrayFullLoadType
    integer(I4B), dimension(:), pointer, contiguous :: activemap => null() !< user node if active, else 0
    type(ParamStampsType) :: written !< write stamps, one column per array
  contains
    procedure :: mode_init
    procedure :: df
    procedure :: ts_advance
    procedure :: array_store
    procedure :: tas_links_update
    procedure :: destroy
    procedure :: stamp_links
    procedure, private :: stamp_name
  end type GridArrayFullLoadType

contains

  !> @brief Set arrays to unset (DNODATA) and map active cells
  !<
  subroutine mode_init(this)
    class(GridArrayFullLoadType), intent(inout) :: this
    integer(I4B), dimension(:), pointer, contiguous :: nodeulist
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    integer(I4B) :: islot, n

    ! no boundary cell list; every user node holds a value
    call mem_allocate(nodeulist, 0, 'NODEULIST', this%mf6_input%mempath)
    this%ctx%nbound = this%ctx%nodes

    do islot = 1, this%nmain
      call mem_setptr(dbl1d, this%slots(islot)%mf6varname, &
                      this%mf6_input%mempath)
      dbl1d = DNODATA
    end do
    if (this%ctx%naux > 0) this%ctx%auxvar = DNODATA

    ! series links update active cells only
    allocate (this%activemap(this%ctx%nodes))
    do n = 1, this%ctx%nodes
      if (this%active(n)) then
        this%activemap(n) = n
      else
        this%activemap(n) = 0
      end if
    end do
  end subroutine mode_init

  !> @brief Allocate write stamps, one column per array (slot)
  !<
  subroutine df(this)
    class(GridArrayFullLoadType), intent(inout) :: this
    character(len=LENVARNAME), dimension(:), allocatable :: varnames
    integer(I4B) :: islot
    call this%GridArrayLoadType%df()
    if (size(this%slots) > 0) then
      allocate (varnames(size(this%slots)))
      do islot = 1, size(this%slots)
        varnames(islot) = this%stamp_name(islot)
      end do
      call this%written%allocate(varnames, this%ctx%nodes, &
                                 this%mf6_input%mempath)
    end if
  end subroutine df

  !> @brief Advance time-array series and stamp the values they write
  !<
  subroutine ts_advance(this)
    class(GridArrayFullLoadType), intent(inout) :: this
    call this%GridArrayLoadType%ts_advance()
    if (associated(this%written%istep)) call this%stamp_links()
  end subroutine ts_advance

  subroutine destroy(this)
    class(GridArrayFullLoadType), intent(inout) :: this
    call this%GridArrayLoadType%destroy()
    if (associated(this%activemap)) deallocate (this%activemap)
    call this%written%destroy()
  end subroutine destroy

  !> @brief Store an array read in the period block
  !!
  !! A defined value at an active cell replaces the stored value; other
  !! cells are unchanged. A series array (nodes not associated) is set by
  !! the time steps.
  !<
  subroutine array_store(this, islot, nodes)
    class(GridArrayFullLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: islot
    real(DP), dimension(:), pointer, contiguous, intent(in) :: nodes
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    integer(I4B) :: n, iaux

    if (.not. associated(nodes)) return
    iaux = this%slots(islot)%iaux
    if (iaux == 0) then
      call mem_setptr(dbl1d, this%slots(islot)%mf6varname, &
                      this%mf6_input%mempath)
    end if

    do n = 1, size(nodes)
      if (nodes(n) == DNODATA .or. .not. this%active(n)) cycle
      if (iaux == 0) then
        dbl1d(n) = nodes(n)
      else
        this%ctx%auxvar(iaux, n) = nodes(n)
      end if
      call this%written%set(islot, n, RS_VALUES)
    end do
  end subroutine array_store

  !> @brief Update time-array series links after a period block
  !!
  !! An array given in the block with values ends its link; an array
  !! given as a series is linked. Links for arrays not given are kept, and
  !! an empty block changes nothing. Links update active cells by user
  !! node and keep the current value where the series is DNODATA.
  !<
  subroutine tas_links_update(this, inunit)
    class(GridArrayFullLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: inunit
    if (.not. this%tas_active) return
    if (this%nreads == 0) return
    call this%make_links(inunit, this%activemap, .true.)
  end subroutine tas_links_update

  !> @brief Stamp the array values written by series links this step
  !<
  subroutine stamp_links(this)
    use TimeArraySeriesLinkModule, only: TimeArraySeriesLinkType
    use InputOutputModule, only: same_word
    class(GridArrayFullLoadType), intent(inout) :: this
    type(TimeArraySeriesLinkType), pointer :: tasLink
    integer(I4B) :: i, islot, n

    do i = 1, this%tasmanager%CountLinks()
      tasLink => this%tasmanager%GetLink(i)
      if (.not. associated(tasLink)) cycle
      if (.not. allocated(tasLink%written)) cycle
      do islot = 1, size(this%slots)
        if (.not. same_word(tasLink%Text, this%stamp_name(islot))) cycle
        do n = 1, size(tasLink%written)
          if (tasLink%written(n)) call this%written%set(islot, n, RS_SERIES)
        end do
      end do
    end do
  end subroutine stamp_links

  !> @brief Stamp column name of an array: its input name, or its
  !! auxiliary variable name for an aux array
  !<
  function stamp_name(this, islot) result(name)
    class(GridArrayFullLoadType), intent(in) :: this
    integer(I4B), intent(in) :: islot
    character(len=LENVARNAME) :: name
    if (this%slots(islot)%iaux > 0) then
      name = this%slots(islot)%tagname(1:LENVARNAME)
    else
      name = this%slots(islot)%mf6varname
    end if
  end function stamp_name

end module GridArrayFullLoadModule
