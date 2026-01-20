! A table term is the information needed to describe flow.
! The table object contains an array of table terms.
! For an advanced package.  The table object describes all of
! the flows.
module TableTermModule

  use KindModule, only: I4B, DP
  use ConstantsModule, only: LINELENGTH, LENBUDTXT, DZERO, &
                             TABLEFT, TABCENTER, TABRIGHT, &
                             TABSTRING, TABUCSTRING, TABINTEGER, TABREAL
  use InputOutputModule, only: UPCASE, parseline

  implicit none

  public :: TableTermType

  type :: TableTermType
    character(len=LINELENGTH), pointer :: tag => null()
    integer(I4B), pointer :: width => null()
    integer(I4B), pointer :: alignment => null()
    integer(I4B), pointer :: nheader_lines => null()

    character(len=LINELENGTH), dimension(:), pointer :: initial_lines => null()
    character(len=LINELENGTH), dimension(:), pointer :: header_lines => null()

  contains

    procedure :: initialize
    procedure, private :: allocate_scalars
    procedure :: get_width
    procedure :: get_alignment
    procedure :: get_header_lines
    procedure :: set_header
    procedure :: get_header
    procedure :: da

  end type TableTermType

contains

  !< @brief initialize the table term
  !<
  subroutine initialize(this, tag, width, alignment)
    ! -- modules
    ! -- dummy
    class(TableTermType) :: this
    character(len=*), intent(in) :: tag
    integer(I4B), intent(in) :: width
    integer(I4B), intent(in), optional :: alignment
    ! -- local
    character(len=LINELENGTH) :: string
    character(len=LINELENGTH) :: tstring
    character(len=LINELENGTH), allocatable, dimension(:) :: words
    integer(I4B) :: nwords
    integer(I4B) :: ilen
    integer(I4B) :: i
    integer(I4B) :: j

    !
    ! -- allocate scalars
    call this%allocate_scalars()

    ! -- process dummy variables
    this%tag = tag

    if (present(alignment)) then
      this%alignment = alignment
    else
      this%alignment = TABCENTER
    end if

    this%width = width
    !
    ! -- parse tag into words
    call parseline(tag, nwords, words, 0)
    !
    ! -- abbreviate any words that exceed the specified width
    !    and trim trailing characters
    do i = 1, nwords
      ilen = len(trim(words(i)))
      if (ilen > width) then
        words(i) (width:width) = '.'
        do j = width + 1, ilen
          words(i) (j:j) = ' '
        end do
      end if
    end do
    !
    ! -- combine words that fit into width
    i = 0
    do
      i = i + 1
      if (i > nwords) then
        exit
      end if
      string = trim(adjustl(words(i)))
      tstring = string
      do j = i + 1, nwords
        if (len(trim(adjustl(string))) > 0) then
          tstring = trim(adjustl(tstring))//' '//trim(adjustl(words(j)))
        else
          tstring = trim(adjustl(words(j)))
        end if
        ilen = len(trim(adjustl(tstring)))
        if (ilen == 0) then
          continue
        else if (ilen <= width) then
          words(j) = ' '
          string = tstring
        else
          exit
        end if
      end do
      words(i) = string
    end do
    !
    ! -- calculate the number of header lines
    do i = 1, nwords
      ilen = len(trim(adjustl(words(i))))
      if (ilen > 0) then
        this%nheader_lines = this%nheader_lines + 1
      end if
    end do
    !
    ! allocate initial_lines and fill with words
    allocate (this%initial_lines(this%nheader_lines))
    do i = 1, this%nheader_lines
      this%initial_lines(i) = words(i) (1:width)
    end do
    !
    ! -- deallocate words
    deallocate (words)
    !
    ! -- return
    return

  end subroutine initialize

  !< @brief get column width
  !<
  function get_width(this)
    ! -- return variable
    integer(I4B) :: get_width
    ! -- modules
    ! -- dummy
    class(TableTermType) :: this
    ! -- local
    get_width = this%width
  end function get_width

  !< @brief get column alignment
  !<
  function get_alignment(this)
    ! -- return variable
    integer(I4B) :: get_alignment
    ! -- modules
    ! -- dummy
    class(TableTermType) :: this
    ! -- local
    get_alignment = this%alignment
  end function get_alignment

  !< @brief get the number of lines in initial_lines
  !<
  function get_header_lines(this)
    ! -- return variable
    integer(I4B) :: get_header_lines
    ! -- modules
    ! -- dummy
    class(TableTermType) :: this
    ! -- local
    get_header_lines = this%nheader_lines
  end function get_header_lines

  !< @brief allocate table term scalars
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    ! -- dummy
    class(TableTermType) :: this
    !
    ! -- allocate scalars
    allocate (this%tag)
    allocate (this%alignment)
    allocate (this%width)
    allocate (this%nheader_lines)
    !
    ! -- initialize scalars
    this%nheader_lines = 0
  end subroutine allocate_scalars

  !< @brief deallocate table terms
  !<
  subroutine da(this)
    ! -- modules
    ! -- dummy
    class(TableTermType) :: this
    ! -- local
    !integer(I4B) :: n
    !
    ! -- deallocate scalars
    deallocate (this%tag)
    deallocate (this%alignment)
    deallocate (this%width)
    deallocate (this%nheader_lines)
    deallocate (this%header_lines)
    !
    ! -- return
  end subroutine da

  !< @brief set final header lines for table term
  !<
  subroutine set_header(this, nlines)
    ! -- modules
    ! -- dummy
    class(TableTermType) :: this
    integer(I4B), intent(in) :: nlines
    ! -- local
    character(len=this%width) :: string
    integer(I4B) :: idiff
    integer(I4B) :: i0
    integer(I4B) :: i
    integer(I4B) :: j
    !
    ! -- initialize variables
    string = ' '
    !
    ! allocate header_lines
    allocate (this%header_lines(nlines))
    !
    ! -- initialize header lines
    do i = 1, nlines
      this%header_lines(i) = string
    end do
    !
    ! -- fill header_lines with initial_lines from
    !    bottom to top
    idiff = nlines - this%nheader_lines
    i0 = 1 - idiff
    do i = this%nheader_lines, 1, -1
      j = i + idiff
      this%header_lines(j) = this%initial_lines(i)
    end do
    !
    ! -- deallocate temporary header lines
    deallocate (this%initial_lines)
    !
    ! -- reinitialize nheader_lines
    this%nheader_lines = nlines
    !
    ! -- return
  end subroutine set_header

  !< @brief get header entry for table term iline
  !<
  subroutine get_header(this, iline, cval)
    ! -- modules
    ! -- dummy
    class(TableTermType) :: this
    integer(I4B), intent(in) :: iline
    character(len=*), intent(inout) :: cval
    ! -- return variable
    ! -- local
    !
    ! -- set return value
    cval = this%header_lines(iline) (1:this%width)
    !
    ! -- return
  end subroutine get_header

end module TableTermModule
