!> @brief Generic List Reader Module
!<
module ListReaderModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH, LENBOUNDNAME, LENTIMESERIESNAME, &
                             LENAUXNAME, LENLISTLABEL, DONE
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, count_errors, store_error_unit
  use LongLineReaderModule, only: LongLineReaderType
  use GeomUtilModule, only: get_ijk, get_jk, get_node

  implicit none
  private
  public ListReaderType

  type :: ListReaderType
    integer(I4B) :: in = 0 !< unit number of file containing control record
    integer(I4B) :: inlist = 0 !< unit number of file from which list will be read
    integer(I4B) :: iout = 0 !< unit number to output messages
    integer(I4B) :: inamedbound = 0 !< flag indicating boundary names are to be read
    integer(I4B) :: ierr = 0 !< error flag
    integer(I4B) :: nlist = 0 !< number of entries in list.  -1 indicates number will be automatically determined
    integer(I4B) :: ibinary = 0 !< flag indicating to read binary list
    integer(I4B) :: istart = 0 !< string starting location
    integer(I4B) :: istop = 0 !< string ending location
    integer(I4B) :: lloc = 0 !< entry number in line
    integer(I4B) :: iclose = 0 !< flag indicating whether or not to close file
    integer(I4B) :: ndim = 0 !< number of dimensions in model
    integer(I4B) :: ntxtrlist = 0 !< number of text entries found in rlist
    integer(I4B) :: ntxtauxvar = 0 !< number of text entries found in auxvar
    character(len=LENLISTLABEL) :: label = '' !< label for printing list
    character(len=:), allocatable, private :: line !< current line
    integer(I4B), dimension(:), pointer, contiguous :: mshape => null() !< pointer to model shape
    integer(I4B), dimension(:), pointer, contiguous :: nodelist => null() !< pointer to nodelist
    real(DP), dimension(:, :), pointer, contiguous :: rlist => null() !< pointer to rlist
    real(DP), dimension(:, :), pointer, contiguous :: auxvar => null() !< pointer to auxvar
    character(len=16), dimension(:), pointer :: auxname => null() !< pointer to aux names
    character(len=LENBOUNDNAME), dimension(:), pointer, &
      contiguous :: boundname => null() !< pointer to boundname
    integer(I4B), dimension(:), allocatable :: idxtxtrow !< row locations of text in rlist
    integer(I4B), dimension(:), allocatable :: idxtxtcol !< col locations of text in rlist
    integer(I4B), dimension(:), allocatable :: idxtxtauxrow !< row locations of text in auxvar
    integer(I4B), dimension(:), allocatable :: idxtxtauxcol !< col locations of text in auxvar
    character(len=LENTIMESERIESNAME), dimension(:), allocatable :: txtrlist !< text found in rlist
    character(len=LENTIMESERIESNAME), dimension(:), allocatable :: txtauxvar !< text found in auxvar
    type(LongLineReaderType), pointer :: line_reader => null()

  contains

    procedure :: read_list
    procedure :: write_list
    procedure, private :: read_control_record
    procedure, private :: read_data
    procedure, private :: set_openclose
    procedure, private :: read_ascii
    procedure, private :: read_binary
  end type ListReaderType

contains

  !> @brief Initialize the reader
  !<
  subroutine read_list(this, line_reader, in, iout, nlist, inamedbound, &
                       mshape, nodelist, rlist, auxvar, auxname, boundname, &
                       label)
    ! -- modules
    use ConstantsModule, only: LENBOUNDNAME
    ! -- dummy
    class(ListReaderType) :: this
    type(LongLineReaderType), intent(inout), target :: line_reader
    integer(I4B), intent(in) :: in
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(inout) :: nlist
    integer(I4B), intent(in) :: inamedbound
    integer(I4B), dimension(:), intent(in), contiguous, pointer :: mshape
    integer(I4B), dimension(:), intent(inout), contiguous, pointer :: nodelist
    real(DP), dimension(:, :), intent(inout), contiguous, pointer :: rlist
    real(DP), dimension(:, :), intent(inout), contiguous, pointer :: auxvar
    character(len=LENAUXNAME), dimension(:), intent(inout), target :: auxname
    character(len=LENBOUNDNAME), &
      dimension(:), pointer, contiguous, intent(inout) :: boundname
    character(len=LENLISTLABEL), intent(in) :: label
    !
    ! -- Copy variables
    this%in = in
    this%iout = iout
    this%nlist = nlist
    this%inamedbound = inamedbound
    this%ndim = size(mshape)
    this%label = label
    !
    ! -- Set pointers
    this%mshape => mshape
    this%nodelist => nodelist
    this%rlist => rlist
    this%auxvar => auxvar
    this%auxname => auxname
    this%boundname => boundname
    this%line_reader => line_reader
    !
    ! -- Allocate arrays for storing text and text locations
    if (.not. allocated(this%idxtxtrow)) allocate (this%idxtxtrow(0))
    if (.not. allocated(this%idxtxtcol)) allocate (this%idxtxtcol(0))
    if (.not. allocated(this%idxtxtauxrow)) allocate (this%idxtxtauxrow(0))
    if (.not. allocated(this%idxtxtauxcol)) allocate (this%idxtxtauxcol(0))
    if (.not. allocated(this%txtrlist)) allocate (this%txtrlist(0))
    if (.not. allocated(this%txtauxvar)) allocate (this%txtauxvar(0))
    !
    ! -- Read control record
    call this%read_control_record()
    !
    ! -- Read data
    call this%read_data()
    !
    ! -- Set nlist for return
    nlist = this%nlist
  end subroutine read_list

  !> @brief Check for a control record, and parse if found
  !<
  subroutine read_control_record(this)
    ! -- modules
    use InputOutputModule, only: urword
    ! -- dummy
    class(ListReaderType) :: this
    ! -- local
    integer(I4B) :: idum
    real(DP) :: r
    ! -- formats
    character(len=*), parameter :: fmtlsf = &
                                   "(1X,'LIST SCALING FACTOR=',1PG12.5)"
    !
    ! -- Set default values, which may be changed by control record
    this%inlist = this%in
    this%iclose = 0
    this%ibinary = 0
    !
    ! -- Read to the first non-commented line
    call this%line_reader%rdcom(this%in, this%iout, this%line, this%ierr)
    this%lloc = 1
    call urword(this%line, this%lloc, this%istart, this%istop, 1, idum, r, &
                this%iout, this%in)
    !
    ! -- Parse record
    select case (this%line(this%istart:this%istop))
    case ('OPEN/CLOSE')
      call this%set_openclose()
    end select
  end subroutine read_control_record

  !> @brief Set up for open/close file
  !!
  !! OPEN/CLOSE fname [(BINARY)]
  !<
  subroutine set_openclose(this)
    ! -- modules
    use InputOutputModule, only: urword, openfile
    use OpenSpecModule, only: form, access
    use ConstantsModule, only: LINELENGTH
    ! -- dummy
    class(ListReaderType) :: this
    ! -- local
    integer(I4B) :: idum, itmp
    real(DP) :: r
    logical :: exists
    integer(I4B) :: nunopn = 99
    character(len=LINELENGTH) :: fname
    ! -- formats
    character(len=*), parameter :: fmtocne = &
      &"('Specified OPEN/CLOSE file ',(A),' does not exist')"
    character(len=*), parameter :: fmtobf = &
      &"(1X,/1X,'OPENING BINARY FILE ON UNIT ',I0,':',/1X,A)"
    character(len=*), parameter :: fmtobfnlist = &
      &"(1X, 'TO READ ', I0, ' RECORDS.')"
    character(len=*), parameter :: fmtofnlist = &
      &"(1x,'TO READ ', I0, ' RECORDS.')"
    character(len=*), parameter :: fmtof = &
      &"(1X,/1X,'OPENING FILE ON UNIT ',I0,':',/1X,A)"
    !
    ! -- Get filename
    call urword(this%line, this%lloc, this%istart, this%istop, 0, idum, r, &
                this%iout, this%in)
    fname = this%line(this%istart:this%istop)
    !
    ! -- Check to see if file OPEN/CLOSE file exists
    inquire (file=fname, exist=exists)
    if (.not. exists) then
      write (errmsg, fmtocne) this%line(this%istart:this%istop)
      call store_error(errmsg)
      call store_error('Specified OPEN/CLOSE file does not exist')
      call store_error_unit(this%in)
    end if
    !
    ! -- Check for (BINARY) keyword
    call urword(this%line, this%lloc, this%istart, this%istop, 1, idum, r, &
                this%iout, this%in)
    if (this%line(this%istart:this%istop) == '(BINARY)') this%ibinary = 1
    !
    ! -- Open the file depending on ibinary flag
    this%inlist = nunopn
    if (this%ibinary == 1) then
      itmp = this%iout
      if (this%iout > 0) then
        itmp = 0
        write (this%iout, fmtobf) this%inlist, trim(adjustl(fname))
        if (this%nlist > 0) write (this%iout, fmtobfnlist) this%nlist
      end if
      call openfile(this%inlist, itmp, fname, 'OPEN/CLOSE', fmtarg_opt=form, &
                    accarg_opt=access)
    else
      itmp = this%iout
      if (this%iout > 0) then
        itmp = 0
        write (this%iout, fmtof) this%inlist, trim(adjustl(fname))
        if (this%nlist > 0) write (this%iout, fmtofnlist) this%nlist
      end if
      call openfile(this%inlist, itmp, fname, 'OPEN/CLOSE')
    end if
    !
    ! -- Set iclose to 1 because it is open/close, to indicate that the
    !    file needs to be closed after the list is read
    this%iclose = 1
    !
    ! -- Read the first line from inlist to be consistent with how the list is
    !    read when it is included in the package input file
    if (this%ibinary /= 1) &
      call this%line_reader%rdcom(this%inlist, this%iout, this%line, &
                                  this%ierr)
  end subroutine set_openclose

  !> @brief Read the data
  !<
  subroutine read_data(this)
    ! -- dummy
    class(ListReaderType) :: this
    !
    ! -- Read the list
    if (this%ibinary == 1) then
      call this%read_binary()
    else
      call this%read_ascii()
    end if
    !
    ! -- If open/close, then close file
    if (this%iclose == 1) then
      close (this%inlist)
    end if
  end subroutine read_data

  !> @brief Read the data from a binary file
  !<
  subroutine read_binary(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH, LENBIGLINE
    ! -- dummy
    class(ListReaderType) :: this
    ! -- local
    integer(I4B) :: mxlist, ldim, naux, nod, ii, jj
    character(len=LINELENGTH) :: fname
    integer(I4B), dimension(:), allocatable :: cellid
    ! -- formats
    character(len=*), parameter :: fmtmxlsterronly = &
      "('ERROR READING LIST FROM FILE: ',&
      &a,' ON UNIT: ',I0,&
      &' THE NUMBER OF RECORDS ENCOUNTERED EXCEEDS THE MAXIMUM NUMBER &
      &OF RECORDS.  TRY INCREASING MAXBOUND FOR THIS LIST.&
      &  NUMBER OF RECORDS: ',I0,' MAXBOUND: ',I0)"
    character(len=*), parameter :: fmtlsterronly = &
      "('ERROR READING LIST FROM FILE: ',&
      &1x,a,1x,' ON UNIT: ',I0)"
    !
    ! -- Determine array sizes
    mxlist = size(this%rlist, 2)
    ldim = size(this%rlist, 1)
    naux = size(this%auxvar, 1)
    !
    ! -- Allocate arrays
    allocate (cellid(this%ndim))
    !
    ii = 1
    readloop: do
      !
      ! -- Read layer, row, col, or cell number
      read (this%inlist, iostat=this%ierr) cellid
      !
      ! -- If not end of record, then store nodenumber, else
      !    calculate lstend and nlist, and exit readloop
      select case (this%ierr)
      case (0)
        !
        ! -- Ensure cellid is valid, store an error otherwise
        call check_cellid(ii, cellid, this%mshape, this%ndim)
        !
        ! -- Check range
        if (ii > mxlist) then
          inquire (unit=this%inlist, name=fname)
          write (errmsg, fmtmxlsterronly) fname, this%inlist, ii, mxlist
          call store_error(errmsg, terminate=.TRUE.)
        end if
        !
        ! -- Calculate and store user node number
        if (this%ndim == 1) then
          nod = cellid(1)
        elseif (this%ndim == 2) then
          nod = get_node(cellid(1), 1, cellid(2), &
                         this%mshape(1), 1, this%mshape(2))
        else
          nod = get_node(cellid(1), cellid(2), cellid(3), &
                         this%mshape(1), this%mshape(2), this%mshape(3))
        end if
        this%nodelist(ii) = nod
        !
        ! -- Read remainder of record
        read (this%inlist, iostat=this%ierr) (this%rlist(jj, ii), jj=1, ldim), &
          (this%auxvar(jj, ii), jj=1, naux)
        if (this%ierr /= 0) then
          inquire (unit=this%inlist, name=fname)
          write (errmsg, fmtlsterronly) trim(adjustl(fname)), this%inlist
          call store_error(errmsg, terminate=.TRUE.)
        end if
        !
      case (:-1)
        !
        ! -- End of record was encountered
        this%nlist = ii - 1
        exit readloop
        !
      case (1:)
        !
        ! -- Error
        inquire (unit=this%inlist, name=fname)
        write (errmsg, fmtlsterronly) trim(adjustl(fname)), this%inlist
        call store_error(errmsg, terminate=.TRUE.)
        !
      end select
      !
      ! -- If nlist is known, then exit when nlist values have been read
      if (this%nlist > 0) then
        if (ii == this%nlist) exit readloop
      end if
      !
      ! -- Increment ii
      ii = ii + 1
      !
    end do readloop
    !
    ! -- Stop if errors were detected
    if (count_errors() > 0) then
      call store_error_unit(this%inlist)
    end if
  end subroutine read_binary

  !> @brief Read the data from an ascii file
  !<
  subroutine read_ascii(this)
    ! -- modules
    use ConstantsModule, only: LENBOUNDNAME, LINELENGTH, DZERO
    use InputOutputModule, only: urword
    use ArrayHandlersModule, only: ExpandArray
    use TdisModule, only: kper
    ! -- dummy
    class(ListReaderType) :: this
    ! -- local
    integer(I4B) :: mxlist, ldim, naux
    integer(I4B) :: ii, jj, idum, nod, istat, increment
    real(DP) :: r
    integer(I4B), dimension(:), allocatable :: cellid
    character(len=LINELENGTH) :: fname
    ! -- formats
    character(len=*), parameter :: fmtmxlsterronly = &
      "('Error reading list. The number of records encountered exceeds &
       &the maximum number of records. Number of records found is ',I0,&
       &' but MAXBOUND is ', I0, '. Try increasing MAXBOUND for this list. &
       &Error occurred reading the following line: ', a, 5x, '>>> ', a)"
    !
    ! -- Determine array sizes
    mxlist = size(this%rlist, 2)
    ldim = size(this%rlist, 1)
    naux = size(this%auxvar, 1)
    this%ntxtrlist = 0
    this%ntxtauxvar = 0
    !
    ! -- Allocate arrays
    allocate (cellid(this%ndim))
    !
    ii = 1
    readloop: do
      !
      ! -- First line was already read, so don't read again
      if (ii /= 1) &
        call this%line_reader%rdcom(this%inlist, 0, this%line, this%ierr)
      !
      ! -- If this is an unknown-length list, then check for END.
      !    If found, then backspace, set nlist, and exit readloop.
      if (this%nlist < 0) then
        this%lloc = 1
        call urword(this%line, this%lloc, this%istart, this%istop, 1, idum, r, &
                    this%iout, this%inlist)
        if (this%line(this%istart:this%istop) == 'END' .or. this%ierr < 0) then
          ! -- If END was found then call line_reader backspace
          !    emulator so that caller can proceed with reading END.
          if (this%ierr == 0) then
            call this%line_reader%bkspc(this%inlist)
          end if
          this%nlist = ii - 1
          exit readloop
        end if
      end if
      !
      ! -- Check range
      if (ii > mxlist) then
        inquire (unit=this%inlist, name=fname)
        write (errmsg, fmtmxlsterronly) &
          ii, mxlist, new_line("A"), trim(this%line)
        call store_error(errmsg)
        call store_error_unit(this%inlist)
      end if
      !
      ! -- Initialize locator
      this%lloc = 1
      !
      ! -- Read cellid
      call urword(this%line, this%lloc, this%istart, this%istop, 2, &
                  cellid(1), r, this%iout, this%inlist)
      if (this%ndim > 1) then
        call urword(this%line, this%lloc, this%istart, this%istop, 2, &
                    cellid(2), r, this%iout, this%inlist)
      end if
      if (this%ndim > 2) then
        call urword(this%line, this%lloc, this%istart, this%istop, 2, &
                    cellid(3), r, this%iout, this%inlist)
      end if
      !
      ! -- Ensure cellid is valid, store an error otherwise
      call check_cellid(ii, cellid, this%mshape, this%ndim)
      !
      ! -- Calculate user node number
      if (this%ndim == 3) then
        nod = get_node(cellid(1), cellid(2), cellid(3), &
                       this%mshape(1), this%mshape(2), this%mshape(3))
      elseif (this%ndim == 2) then
        nod = get_node(cellid(1), 1, cellid(2), &
                       this%mshape(1), 1, this%mshape(2))
      else
        nod = cellid(1)
      end if
      !
      ! -- Assign nod to nodelist
      this%nodelist(ii) = nod
      !
      ! -- Read rlist
      do jj = 1, ldim
        call urword(this%line, this%lloc, this%istart, this%istop, 0, idum, &
                    r, this%iout, this%inlist)
        read (this%line(this%istart:this%istop), *, iostat=istat) r
        !
        ! -- If a double precision value, then store in rlist, otherwise store
        !    the text name and location
        if (istat == 0) then
          this%rlist(jj, ii) = r
        else
          this%rlist(jj, ii) = DZERO
          this%ntxtrlist = this%ntxtrlist + 1
          if (this%ntxtrlist > size(this%txtrlist)) then
            increment = int(size(this%txtrlist) * 0.2)
            increment = max(100, increment)
            call ExpandArray(this%txtrlist, increment)
            call ExpandArray(this%idxtxtrow, increment)
            call ExpandArray(this%idxtxtcol, increment)
          end if
          this%txtrlist(this%ntxtrlist) = this%line(this%istart:this%istop)
          this%idxtxtrow(this%ntxtrlist) = ii
          this%idxtxtcol(this%ntxtrlist) = jj
        end if
        !
      end do
      !
      ! -- Read auxvar
      do jj = 1, naux
        call urword(this%line, this%lloc, this%istart, this%istop, 0, idum, &
                    r, this%iout, this%inlist)
        read (this%line(this%istart:this%istop), *, iostat=istat) r
        !
        ! -- If a double precision value, then store in auxvar, otherwise store
        !    the text name and location
        if (istat == 0) then
          this%auxvar(jj, ii) = r
        else
          this%auxvar(jj, ii) = DZERO
          this%ntxtauxvar = this%ntxtauxvar + 1
          if (this%ntxtauxvar > size(this%txtauxvar)) then
            increment = int(size(this%txtauxvar) * 0.2)
            increment = max(100, increment)
            call ExpandArray(this%txtauxvar, increment)
            call ExpandArray(this%idxtxtauxrow, increment)
            call ExpandArray(this%idxtxtauxcol, increment)
          end if
          this%txtauxvar(this%ntxtauxvar) = this%line(this%istart:this%istop)
          this%idxtxtauxrow(this%ntxtauxvar) = ii
          this%idxtxtauxcol(this%ntxtauxvar) = jj
          if (len_trim(this%txtauxvar(this%ntxtauxvar)) == 0) then
            write (errmsg, '(a,i0,a)') 'Auxiliary data or time series name &
                                        &expected but not found in period &
                                        &block "', kper, '".'
            call store_error(errmsg)
            call store_error_unit(this%inlist)
          end if
        end if
        !
      end do
      !
      ! -- Read the boundary names (only supported for ascii input)
      if (this%inamedbound > 0) then
        call urword(this%line, this%lloc, this%istart, this%istop, 1, idum, r, &
                    this%iout, this%inlist)
        this%boundname(ii) = this%line(this%istart:this%istop)
      end if
      !
      ! -- If nlist is known, then exit when nlist values have been read
      if (this%nlist > 0) then
        if (ii == this%nlist) exit readloop
      end if
      !
      ! -- Increment ii row counter
      ii = ii + 1
      !
    end do readloop
    !
    ! -- Stop if errors were detected
    if (count_errors() > 0) then
      call store_error_unit(this%inlist)
    end if
  end subroutine read_ascii

  !> @brief Check for valid cellid
  !<
  subroutine check_cellid(ii, cellid, mshape, ndim)
    ! -- dummy
    integer(I4B), intent(in) :: ii
    integer(I4B), dimension(:), intent(in) :: cellid !< cellid
    integer(I4B), dimension(:), intent(in) :: mshape !< model shape
    integer(I4B), intent(in) :: ndim !< size of mshape
    ! -- local
    character(len=20) :: cellstr, mshstr
    ! -- formats
    character(len=*), parameter :: fmterr = &
      "('List entry ',i0,' contains cellid ',a,' but this cellid is invalid &
      &for model with shape ', a)"
    character(len=*), parameter :: fmtndim1 = &
                                   "('(',i0,')')"
    character(len=*), parameter :: fmtndim2 = &
                                   "('(',i0,',',i0,')')"
    character(len=*), parameter :: fmtndim3 = &
                                   "('(',i0,',',i0,',',i0,')')"
    !
    if (ndim == 1) then
      if (cellid(1) < 1 .or. cellid(1) > mshape(1)) then
        write (cellstr, fmtndim1) cellid(1)
        write (mshstr, fmtndim1) mshape(1)
        write (errmsg, fmterr) ii, trim(adjustl(cellstr)), trim(adjustl(mshstr))
        call store_error(errmsg)
      end if
    else if (ndim == 2) then
      if (cellid(1) < 1 .or. cellid(1) > mshape(1) .or. &
          cellid(2) < 1 .or. cellid(2) > mshape(2)) then
        write (cellstr, fmtndim2) cellid(1), cellid(2)
        write (mshstr, fmtndim2) mshape(1), mshape(2)
        write (errmsg, fmterr) ii, trim(adjustl(cellstr)), trim(adjustl(mshstr))
        call store_error(errmsg)
      end if
    else if (ndim == 3) then
      if (cellid(1) < 1 .or. cellid(1) > mshape(1) .or. &
          cellid(2) < 1 .or. cellid(2) > mshape(2) .or. &
          cellid(3) < 1 .or. cellid(3) > mshape(3)) then
        write (cellstr, fmtndim3) cellid(1), cellid(2), cellid(3)
        write (mshstr, fmtndim3) mshape(1), mshape(2), mshape(3)
        write (errmsg, fmterr) ii, trim(adjustl(cellstr)), trim(adjustl(mshstr))
        call store_error(errmsg)
      end if
    end if
  end subroutine check_cellid

  !> @brief Write input data to a list
  !<
  subroutine write_list(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH, LENBOUNDNAME, &
                               TABLEFT, TABCENTER
    use InputOutputModule, only: ulstlb
    use TableModule, only: TableType, table_cr
    ! -- dummy
    class(ListReaderType) :: this
    ! -- local
    character(len=10) :: cpos
    character(len=LINELENGTH) :: tag
    character(len=LINELENGTH), allocatable, dimension(:) :: words
    integer(I4B) :: ntabrows
    integer(I4B) :: ntabcols
    integer(I4B) :: ipos
    integer(I4B) :: ii, jj, i, j, k, nod
    integer(I4B) :: ldim
    integer(I4B) :: naux
    type(TableType), pointer :: inputtab => null()
    ! -- formats
    character(len=LINELENGTH) :: fmtlstbn
    !
    ! -- Determine sizes
    ldim = size(this%rlist, 1)
    naux = size(this%auxvar, 1)
    !
    ! -- Dimension table
    ntabrows = this%nlist
    !
    ! -- Start building format statement to parse this%label, which
    !    contains the column headers (except for boundname and auxnames)
    ipos = index(this%label, 'NO.')
    if (ipos /= 0) then
      write (cpos, '(i10)') ipos + 3
      fmtlstbn = '(a'//trim(adjustl(cpos))
    else
      fmtlstbn = '(a7'
    end if
    ! -- Sequence number, layer, row, and column.
    if (size(this%mshape) == 3) then
      ntabcols = 4
      fmtlstbn = trim(fmtlstbn)//',a7,a7,a7'
      !
      ! -- Sequence number, layer, and cell2d.
    else if (size(this%mshape) == 2) then
      ntabcols = 3
      fmtlstbn = trim(fmtlstbn)//',a7,a7'
      !
      ! -- Sequence number and node.
    else
      ntabcols = 2
      fmtlstbn = trim(fmtlstbn)//',a7'
    end if
    !
    ! -- Add fields for non-optional real values
    ntabcols = ntabcols + ldim
    do i = 1, ldim
      fmtlstbn = trim(fmtlstbn)//',a16'
    end do
    !
    ! -- Add field for boundary name
    if (this%inamedbound == 1) then
      ntabcols = ntabcols + 1
      fmtlstbn = trim(fmtlstbn)//',a16'
    end if
    !
    ! -- Add fields for auxiliary variables
    ntabcols = ntabcols + naux
    do i = 1, naux
      fmtlstbn = trim(fmtlstbn)//',a16'
    end do
    fmtlstbn = trim(fmtlstbn)//')'
    !
    ! -- Allocate words
    allocate (words(ntabcols))
    !
    ! -- Parse this%label into words
    read (this%label, fmtlstbn) (words(i), i=1, ntabcols)
    !
    ! -- Initialize the input table object
    call table_cr(inputtab, ' ', ' ')
    call inputtab%table_df(ntabrows, ntabcols, this%iout)
    !
    ! -- Add the columns
    ipos = 1
    call inputtab%initialize_column(words(ipos), 10, alignment=TABCENTER)
    !
    ! -- Discretization
    do i = 1, size(this%mshape)
      ipos = ipos + 1
      call inputtab%initialize_column(words(ipos), 7, alignment=TABCENTER)
    end do
    !
    ! -- Non-optional variables
    do i = 1, ldim
      ipos = ipos + 1
      call inputtab%initialize_column(words(ipos), 16, alignment=TABCENTER)
    end do
    !
    ! -- Boundname
    if (this%inamedbound == 1) then
      ipos = ipos + 1
      tag = 'BOUNDNAME'
      call inputtab%initialize_column(tag, LENBOUNDNAME, alignment=TABLEFT)
    end if
    !
    ! -- Aux variables
    do i = 1, naux
      call inputtab%initialize_column(this%auxname(i), 16, alignment=TABCENTER)
    end do
    !
    ! -- Write the table
    do ii = 1, this%nlist
      call inputtab%add_term(ii)
      !
      ! -- Discretization
      if (size(this%mshape) == 3) then
        nod = this%nodelist(ii)
        call get_ijk(nod, this%mshape(2), this%mshape(3), this%mshape(1), &
                     i, j, k)
        call inputtab%add_term(k)
        call inputtab%add_term(i)
        call inputtab%add_term(j)
      else if (size(this%mshape) == 2) then
        nod = this%nodelist(ii)
        call get_ijk(nod, 1, this%mshape(2), this%mshape(1), i, j, k)
        call inputtab%add_term(k)
        call inputtab%add_term(j)
      else
        nod = this%nodelist(ii)
        call inputtab%add_term(nod)
      end if
      !
      ! -- Non-optional variables
      do jj = 1, ldim
        call inputtab%add_term(this%rlist(jj, ii))
      end do
      !
      ! -- Boundname
      if (this%inamedbound == 1) then
        call inputtab%add_term(this%boundname(ii))
      end if
      !
      ! -- Aux variables
      do jj = 1, naux
        call inputtab%add_term(this%auxvar(jj, ii))
      end do
    end do
    !
    ! -- Deallocate the local variables
    call inputtab%table_da()
    deallocate (inputtab)
    nullify (inputtab)
    deallocate (words)
  end subroutine write_list

end module ListReaderModule
