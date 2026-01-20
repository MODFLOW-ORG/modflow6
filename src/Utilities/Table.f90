! Comprehensive table object that stores all of the
! intercell flows, and the inflows and the outflows for
! an advanced package.
module TableModule

  use KindModule, only: I4B, I8B, DP
  use ConstantsModule, only: LINELENGTH, LENBUDTXT, &
                             TABSTRING, TABUCSTRING, TABINTEGER, TABREAL, &
                             TABCENTER, &
                             DHNOFLO, DHDRY
  use TableTermModule, only: TableTermType
  use InputOutputModule, only: UWWORD, parseline
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg

  implicit none

  public :: TableType
  public :: table_cr

  type :: TableType
    !
    ! -- name, number of control volumes, and number of table terms
    character(len=LENBUDTXT) :: name
    character(len=LINELENGTH) :: title
    character(len=1), pointer :: sep => null()
    logical, pointer :: write_csv => null()
    logical, pointer :: first_entry => null()
    logical, pointer :: transient => null()
    logical, pointer :: add_linesep => null()
    logical, pointer :: allow_finalization => null()
    integer(I4B), pointer :: iout => null()
    integer(I4B), pointer :: maxbound => null()
    integer(I4B), pointer :: nheaderlines => null()
    integer(I4B), pointer :: nlinewidth => null()
    integer(I4B), pointer :: ntableterm => null()
    integer(I4B), pointer :: ientry => null()
    integer(I4B), pointer :: iloc => null()
    integer(I4B), pointer :: icount => null()
    integer(I4B), pointer :: kstp => null()
    integer(I4B), pointer :: kper => null()
    !
    ! -- array of table terms, with one separate entry for each term
    !    such as rainfall, et, leakage, etc.
    type(TableTermType), dimension(:), pointer :: tableterm => null()
    !
    ! -- table table object, for writing the typical MODFLOW table
    type(TableType), pointer :: table => null()

    character(len=LINELENGTH), pointer :: linesep => null()
    character(len=LINELENGTH), pointer :: dataline => null()
    character(len=LINELENGTH), dimension(:), pointer :: header => null()

  contains

    procedure :: table_df
    procedure :: table_da
    procedure :: initialize_column
    procedure :: line_to_columns
    procedure :: finalize_table
    procedure :: set_maxbound
    procedure :: set_kstpkper
    procedure :: set_title
    procedure :: set_iout
    procedure :: print_list_entry
    procedure :: print_separator

    procedure, private :: allocate_strings
    procedure, private :: set_header
    procedure, private :: write_header
    procedure, private :: write_line
    procedure, private :: finalize
    procedure, private :: add_error
    procedure, private :: reset

    generic, public :: add_term => add_integer, add_long_integer, &
      add_real, add_string
    procedure, private :: add_integer, add_long_integer, add_real, add_string

  end type TableType

contains

  !< @brief Create a new table object
  !<
  subroutine table_cr(this, name, title)
    ! -- modules
    ! -- dummy
    type(TableType), pointer :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: title
    ! -- local
    !
    ! -- check if table already associated and reset if necessary
    if (associated(this)) then
      call this%table_da()
      deallocate (this)
      nullify (this)
    end if
    !
    ! -- Create the object
    allocate (this)
    !
    ! -- initialize variables
    this%name = name
    this%title = title
  end subroutine table_cr

  !< @brief Define the new table object
  subroutine table_df(this, maxbound, ntableterm, iout, transient, &
                      lineseparator, separator, finalize)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    integer(I4B), intent(in) :: maxbound
    integer(I4B), intent(in) :: ntableterm
    integer(I4B), intent(in) :: iout
    logical, intent(in), optional :: transient
    logical, intent(in), optional :: lineseparator
    character(len=1), intent(in), optional :: separator
    logical, intent(in), optional :: finalize
    !
    ! -- allocate scalars
    allocate (this%sep)
    allocate (this%write_csv)
    allocate (this%first_entry)
    allocate (this%transient)
    allocate (this%add_linesep)
    allocate (this%allow_finalization)
    allocate (this%iout)
    allocate (this%maxbound)
    allocate (this%nheaderlines)
    allocate (this%nlinewidth)
    allocate (this%ntableterm)
    allocate (this%ientry)
    allocate (this%iloc)
    allocate (this%icount)
    !
    ! -- allocate space for tableterm
    allocate (this%tableterm(ntableterm))
    !
    ! -- initialize values based on optional dummy variables
    if (present(transient)) then
      this%transient = transient
      allocate (this%kstp)
      allocate (this%kper)
    else
      this%transient = .FALSE.
    end if
    if (present(separator)) then
      this%sep = separator
      if (separator == ',') then
        this%write_csv = .TRUE.
      else
        this%write_csv = .FALSE.
      end if
    else
      this%sep = ' '
      this%write_csv = .FALSE.
    end if
    if (present(lineseparator)) then
      this%add_linesep = lineseparator
    else
      this%add_linesep = .TRUE.
    end if
    if (present(finalize)) then
      this%allow_finalization = finalize
    else
      this%allow_finalization = .TRUE.
    end if
    !
    ! -- initialize variables
    this%first_entry = .TRUE.
    this%iout = iout
    this%maxbound = maxbound
    this%ntableterm = ntableterm
    this%ientry = 0
    this%icount = 0
  end subroutine table_df

  !< @brief Initialize data for a column
  !<
  subroutine initialize_column(this, text, width, alignment)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    character(len=*), intent(in) :: text
    integer(I4B), intent(in) :: width
    integer(I4B), intent(in), optional :: alignment
    ! -- local
    integer(I4B) :: idx
    integer(I4B) :: ialign
    !
    ! -- process optional dummy variables
    if (present(alignment)) then
      ialign = alignment
    else
      ialign = TABCENTER
    end if
    !
    ! -- update index for tableterm
    this%ientry = this%ientry + 1
    idx = this%ientry
    !
    ! -- check that ientry is in bounds
    if (this%ientry > this%ntableterm) then
      write (errmsg, '(a,a,a,i0,a,1x,a,1x,a,a,a,1x,i0,1x,a)') &
        'Trying to add column "', trim(adjustl(text)), '" (', &
        this%ientry, ') in the', trim(adjustl(this%name)), 'table ("', &
        trim(adjustl(this%title)), '") that only has', this%ntableterm, &
        'columns.'
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- initialize table term
    call this%tableterm(idx)%initialize(text, width, alignment=ialign)
    !
    ! -- create header when all terms have been specified
    if (this%ientry == this%ntableterm) then
      call this%set_header()
      !
      ! -- reset ientry
      this%ientry = 0
    end if
  end subroutine initialize_column

  !< @brief Set the table object header
  !<
  subroutine set_header(this)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    ! -- local
    character(len=LINELENGTH) :: cval
    integer(I4B) :: width
    integer(I4B) :: alignment
    integer(I4B) :: nlines
    integer(I4B) :: iloc
    integer(I4B) :: ival
    real(DP) :: rval
    integer(I4B) :: j
    integer(I4B) :: n
    integer(I4B) :: nn
    !
    ! -- initialize variables
    width = 0
    nlines = 0
    !
    ! -- determine total width and maximum number of lines
    do n = 1, this%ntableterm
      width = width + this%tableterm(n)%get_width()
      nlines = max(nlines, this%tableterm(n)%get_header_lines())
    end do
    !
    ! -- add length of separators
    width = width + this%ntableterm - 1
    !
    ! -- allocate the header and line separator
    call this%allocate_strings(width, nlines)
    !
    ! -- build final header lines
    do n = 1, this%ntableterm
      call this%tableterm(n)%set_header(nlines)
    end do
    !
    ! -- build header
    do n = 1, nlines
      iloc = 1
      this%iloc = 1
      if (this%add_linesep) then
        nn = n + 1
      else
        nn = n
      end if
      do j = 1, this%ntableterm
        width = this%tableterm(j)%get_width()
        alignment = this%tableterm(j)%get_alignment()
        call this%tableterm(j)%get_header(n, cval)
        if (this%write_csv) then
          if (j == 1) then
            write (this%header(nn), '(a)') trim(adjustl(cval))
          else
            write (this%header(nn), '(a,",",G0)') &
              trim(this%header(nn)), trim(adjustl(cval))
          end if
        else
          if (j == this%ntableterm) then
            call UWWORD(this%header(nn), iloc, width, TABUCSTRING, &
                        cval(1:width), ival, rval, ALIGNMENT=alignment)
          else
            call UWWORD(this%header(nn), iloc, width, TABUCSTRING, &
                        cval(1:width), ival, rval, ALIGNMENT=alignment, &
                        SEP=this%sep)
          end if
        end if
      end do
    end do
  end subroutine set_header

  !< @brief Allocate allocatable character arrays
  !<
  subroutine allocate_strings(this, width, nlines)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    integer(I4B), intent(in) :: width
    integer(I4B), intent(in) :: nlines
    ! -- local
    character(len=width) :: string
    character(len=width) :: linesep
    integer(I4B) :: n
    !
    ! -- initialize local variables
    string = ''
    linesep = repeat('-', width)
    !
    ! -- initialize variables
    this%nheaderlines = nlines
    if (this%add_linesep) then
      this%nheaderlines = this%nheaderlines + 2
    end if
    this%nlinewidth = width
    !
    ! -- allocate deferred length strings
    allocate (this%header(this%nheaderlines))
    allocate (this%linesep)
    allocate (this%dataline)
    !
    ! -- initialize lines
    this%linesep = linesep(1:width)
    this%dataline = string(1:width)
    do n = 1, this%nheaderlines
      this%header(n) = string(1:width)
    end do
    !
    ! -- fill first and last header line with
    !    linesep
    if (this%add_linesep) then
      this%header(1) = linesep(1:width)
      this%header(nlines + 2) = linesep(1:width)
    end if
  end subroutine allocate_strings

  !< @brief Write the table header
  !<
  subroutine write_header(this)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    ! -- local
    character(len=LINELENGTH) :: title
    integer(I4B) :: width
    integer(I4B) :: n
    !
    ! -- initialize local variables
    width = this%nlinewidth
    !
    ! -- write the table header
    if (this%first_entry) then
      ! -- write title
      title = this%title
      if (this%transient) then
        write (title, '(a,a,i6)') trim(adjustl(title)), '   PERIOD ', this%kper
        write (title, '(a,a,i8)') trim(adjustl(title)), '   STEP ', this%kstp
      end if
      if (len_trim(title) > 0) then
        write (this%iout, '(/,1x,a)') trim(adjustl(title))
      end if
      !
      ! -- write header
      do n = 1, this%nheaderlines
        write (this%iout, '(1x,a)') this%header(n) (1:width)
      end do
    end if
    !
    ! -- reinitialize variables
    this%first_entry = .FALSE.
    this%ientry = 0
    this%icount = 0
  end subroutine write_header

  !< @brief Write the data line
  !<
  subroutine write_line(this)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    ! -- local
    integer(I4B) :: width
    !
    ! -- initialize local variables
    width = this%nlinewidth
    !
    ! -- write the dataline
    write (this%iout, '(1x,a)') this%dataline(1:width)
    !
    ! -- update column and line counters
    this%ientry = 0
    this%iloc = 1
    this%icount = this%icount + 1
  end subroutine write_line

  !< @brief Private method that test for last line. If last line the
  !!        public finalize_table method is called
  !<
  subroutine finalize(this)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    ! -- local
    !
    ! -- finalize table if last entry
    if (this%icount == this%maxbound) then
      call this%finalize_table()
    end if
  end subroutine finalize

  !< @brief Public method to finalize the table
  !<
  subroutine finalize_table(this)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    ! -- local
    !
    ! -- write the final table separator
    call this%print_separator(iextralines=1)
    !
    ! -- flush file
    flush (this%iout)
    !
    ! -- reinitialize variables
    call this%reset()
  end subroutine finalize_table

  !< @brief deallocate
  !<
  subroutine table_da(this)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    ! -- dummy
    integer(I4B) :: i
    !
    ! -- deallocate each table term
    do i = 1, this%ntableterm
      call this%tableterm(i)%da()
    end do
    !
    ! -- deallocate space for tableterm
    deallocate (this%tableterm)
    !
    ! -- deallocate character scalars and arrays
    deallocate (this%linesep)
    deallocate (this%dataline)
    deallocate (this%header)
    !
    ! -- deallocate scalars
    if (this%transient) then
      deallocate (this%kstp)
      deallocate (this%kper)
    end if
    deallocate (this%sep)
    deallocate (this%write_csv)
    deallocate (this%first_entry)
    deallocate (this%transient)
    deallocate (this%add_linesep)
    deallocate (this%allow_finalization)
    deallocate (this%iout)
    deallocate (this%maxbound)
    deallocate (this%nheaderlines)
    deallocate (this%nlinewidth)
    deallocate (this%ntableterm)
    deallocate (this%ientry)
    deallocate (this%iloc)
    deallocate (this%icount)
  end subroutine table_da

  !< @brief convert a line to the correct number of columns
  !<
  subroutine line_to_columns(this, line)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    character(len=LINELENGTH), intent(in) :: line
    ! -- local
    character(len=LINELENGTH), allocatable, dimension(:) :: words
    integer(I4B) :: nwords
    integer(I4B) :: icols
    integer(I4B) :: i
    !
    ! -- write header
    if (this%icount == 0 .and. this%ientry == 0) then
      call this%write_header()
    end if
    !
    ! -- parse line into words
    call parseline(line, nwords, words, 0)
    !
    ! -- calculate the number of entries in line but
    !    limit it to the maximum number of columns if
    !    the number of words exceeds ntableterm
    icols = this%ntableterm
    icols = min(nwords, icols)
    !
    ! -- add data (as strings) to line
    do i = 1, icols
      call this%add_term(words(i))
    end do
    !
    ! -- add empty strings to complete the line
    do i = icols + 1, this%ntableterm
      call this%add_term(' ')
    end do
    !
    ! -- clean up local allocatable array
    deallocate (words)
  end subroutine line_to_columns

  !< @brief evaluate if error condition occurs when adding data to dataline
  !<
  subroutine add_error(this)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    ! -- local
    !
    ! -- check that ientry is within bounds
    if (this%ientry > this%ntableterm) then
      write (errmsg, '(a,1x,i0,5(1x,a),1x,i0,1x,a)') &
        'Trying to add data to column ', this%ientry, 'in the', &
        trim(adjustl(this%name)), 'table (', trim(adjustl(this%title)), &
        ') that only has', this%ntableterm, 'columns.'
      call store_error(errmsg, terminate=.TRUE.)
    end if
  end subroutine add_error

  !< @brief add integer value to the dataline
  !<
  subroutine add_integer(this, ival)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    integer(I4B), intent(in) :: ival
    ! -- local
    logical :: line_end
    character(len=LINELENGTH) :: cval
    real(DP) :: rval
    integer(I4B) :: width
    integer(I4B) :: alignment
    integer(I4B) :: j
    !
    ! -- write header
    if (this%icount == 0 .and. this%ientry == 0) then
      call this%write_header()
    end if
    !
    ! -- update index for tableterm
    this%ientry = this%ientry + 1
    !
    ! -- check that ientry is within bounds
    call this%add_error()
    !
    ! -- initialize local variables
    j = this%ientry
    width = this%tableterm(j)%get_width()
    alignment = this%tableterm(j)%get_alignment()
    line_end = .FALSE.
    if (j == this%ntableterm) then
      line_end = .TRUE.
    end if
    !
    ! -- add data to line
    if (this%write_csv) then
      if (j == 1) then
        write (this%dataline, '(G0)') ival
      else
        write (this%dataline, '(a,",",G0)') trim(this%dataline), ival
      end if
    else
      if (j == this%ntableterm) then
        call UWWORD(this%dataline, this%iloc, width, TABINTEGER, &
                    cval, ival, rval, ALIGNMENT=alignment)
      else
        call UWWORD(this%dataline, this%iloc, width, TABINTEGER, &
                    cval, ival, rval, ALIGNMENT=alignment, SEP=this%sep)
      end if
    end if
    !
    ! -- write the data line, if necessary
    if (line_end) then
      call this%write_line()
    end if
    !
    ! -- finalize the table, if necessary
    if (this%allow_finalization) then
      call this%finalize()
    end if
  end subroutine add_integer

  !< @brief add long integer value to the dataline
  !<
  subroutine add_long_integer(this, long_ival)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    integer(I8B), intent(in) :: long_ival
    ! -- local
    logical :: line_end
    character(len=LINELENGTH) :: cval
    real(DP) :: rval
    integer(I4B) :: ival
    integer(I4B) :: width
    integer(I4B) :: alignment
    integer(I4B) :: j
    !
    ! -- write header
    if (this%icount == 0 .and. this%ientry == 0) then
      call this%write_header()
    end if
    !
    ! -- update index for tableterm
    this%ientry = this%ientry + 1
    !
    ! -- check that ientry is within bounds
    call this%add_error()
    !
    ! -- initialize local variables
    j = this%ientry
    width = this%tableterm(j)%get_width()
    alignment = this%tableterm(j)%get_alignment()
    line_end = .FALSE.
    if (j == this%ntableterm) then
      line_end = .TRUE.
    end if
    !
    ! -- add data to line
    if (this%write_csv) then
      if (j == 1) then
        write (this%dataline, '(G0)') long_ival
      else
        write (this%dataline, '(a,",",G0)') trim(this%dataline), long_ival
      end if
    else
      write (cval, '(i0)') long_ival
      if (j == this%ntableterm) then
        call UWWORD(this%dataline, this%iloc, width, TABSTRING, &
                    trim(cval), ival, rval, ALIGNMENT=alignment)
      else
        call UWWORD(this%dataline, this%iloc, width, TABSTRING, &
                    trim(cval), ival, rval, ALIGNMENT=alignment, SEP=this%sep)
      end if
    end if
    !
    ! -- write the data line, if necessary
    if (line_end) then
      call this%write_line()
    end if
    !
    ! -- finalize the table, if necessary
    if (this%allow_finalization) then
      call this%finalize()
    end if
  end subroutine add_long_integer

  !< @brief add real value to the dataline
  !<
  subroutine add_real(this, rval)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    real(DP), intent(in) :: rval
    ! -- local
    logical :: line_end
    character(len=LINELENGTH) :: cval
    integer(I4B) :: ival
    integer(I4B) :: j
    integer(I4B) :: width
    integer(I4B) :: alignment

    if (rval == DHNOFLO) then
      call this%add_string("INACTIVE")
    else if (rval == DHDRY) then
      call this%add_string("DRY")
    else
      !
      ! -- write header
      if (this%icount == 0 .and. this%ientry == 0) then
        call this%write_header()
      end if
      !
      ! -- update index for tableterm
      this%ientry = this%ientry + 1
      !
      ! -- check that ientry is within bounds
      call this%add_error()
      !
      ! -- initialize local variables
      j = this%ientry
      width = this%tableterm(j)%get_width()
      alignment = this%tableterm(j)%get_alignment()
      line_end = .FALSE.
      if (j == this%ntableterm) then
        line_end = .TRUE.
      end if
      !
      ! -- add data to line
      if (this%write_csv) then
        if (j == 1) then
          write (this%dataline, '(G0)') rval
        else
          write (this%dataline, '(a,",",G0)') trim(this%dataline), rval
        end if
      else
        if (j == this%ntableterm) then
          call UWWORD(this%dataline, this%iloc, width, TABREAL, &
                      cval, ival, rval, ALIGNMENT=alignment)
        else
          call UWWORD(this%dataline, this%iloc, width, TABREAL, &
                      cval, ival, rval, ALIGNMENT=alignment, SEP=this%sep)
        end if
      end if
      !
      ! -- write the data line, if necessary
      if (line_end) then
        call this%write_line()
      end if
      !
      ! -- finalize the table, if necessary
      if (this%allow_finalization) then
        call this%finalize()
      end if
    end if
  end subroutine add_real

  !< @brief add string value to the dataline
  !<
  subroutine add_string(this, cval)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    character(len=*) :: cval
    ! -- local
    logical :: line_end
    integer(I4B) :: j
    integer(I4B) :: ival
    real(DP) :: rval
    integer(I4B) :: width
    integer(I4B) :: alignment
    !
    ! -- write header
    if (this%icount == 0 .and. this%ientry == 0) then
      call this%write_header()
    end if
    !
    ! -- update index for tableterm
    this%ientry = this%ientry + 1
    !
    ! -- check that ientry is within bounds
    call this%add_error()
    !
    ! -- initialize local variables
    j = this%ientry
    width = this%tableterm(j)%get_width()
    alignment = this%tableterm(j)%get_alignment()
    line_end = .FALSE.
    if (j == this%ntableterm) then
      line_end = .TRUE.
    end if
    !
    ! -- add data to line
    if (this%write_csv) then
      if (j == 1) then
        write (this%dataline, '(a)') trim(adjustl(cval))
      else
        write (this%dataline, '(a,",",a)') &
          trim(this%dataline), trim(adjustl(cval))
      end if
    else
      if (j == this%ntableterm) then
        call UWWORD(this%dataline, this%iloc, width, TABSTRING, &
                    cval, ival, rval, ALIGNMENT=alignment)
      else
        call UWWORD(this%dataline, this%iloc, width, TABSTRING, &
                    cval, ival, rval, ALIGNMENT=alignment, SEP=this%sep)
      end if
    end if
    !
    ! -- write the data line, if necessary
    if (line_end) then
      call this%write_line()
    end if
    !
    ! -- finalize the table, if necessary
    if (this%allow_finalization) then
      call this%finalize()
    end if
  end subroutine add_string

  !< @brief reset maxbound
  !<
  subroutine set_maxbound(this, maxbound)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    integer(I4B), intent(in) :: maxbound
    ! -- local
    !
    ! -- set maxbound
    this%maxbound = maxbound
    !
    ! -- reset counters
    call this%reset()
  end subroutine set_maxbound

  !< @brief reset kstp and kper
  !<
  subroutine set_kstpkper(this, kstp, kper)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    integer(I4B), intent(in) :: kstp
    integer(I4B), intent(in) :: kper
    ! -- local
    !
    ! -- set maxbound
    this%kstp = kstp
    this%kper = kper
  end subroutine set_kstpkper

  !< @brief reset title
  !<
  subroutine set_title(this, title)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    character(len=*), intent(in) :: title
    ! -- local
    !
    ! -- set maxbound
    this%title = title
  end subroutine set_title

  !< @brief reset iout
  !<
  subroutine set_iout(this, iout)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    integer(I4B), intent(in) :: iout
    ! -- local
    !
    ! -- set iout
    this%iout = iout
  end subroutine set_iout

  !< @brief print list entry
  !<
  subroutine print_list_entry(this, i, nodestr, q, bname)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    integer(I4B), intent(in) :: i
    character(len=*), intent(in) :: nodestr
    real(DP), intent(in) :: q
    character(len=*), intent(in) :: bname
    ! -- local
    !
    ! -- fill table terms
    call this%add_term(i)
    call this%add_term(nodestr)
    call this%add_term(q)
    if (this%ntableterm > 3) then
      call this%add_term(bname)
    end if
  end subroutine print_list_entry

  !< @brief print separator
  !<
  subroutine print_separator(this, iextralines)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    integer(I4B), optional :: iextralines
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: iextra
    integer(I4B) :: width
    !
    ! -- process optional variables
    if (present(iextralines)) then
      iextra = iextralines
    else
      iextra = 0
    end if
    !
    ! -- initialize local variables
    width = this%nlinewidth
    !
    ! -- print line separator
    if (this%add_linesep) then
      write (this%iout, '(1x,a)') this%linesep(1:width)
      do i = 1, iextra
        write (this%iout, '(/)')
      end do
    end if
  end subroutine print_separator

  !< @brief Private method to reset table counters
  !<
  subroutine reset(this)
    ! -- modules
    ! -- dummy
    class(TableType) :: this
    ! -- local
    !
    ! -- reset counters
    this%ientry = 0
    this%icount = 0
    this%first_entry = .TRUE.
  end subroutine reset

end module TableModule
