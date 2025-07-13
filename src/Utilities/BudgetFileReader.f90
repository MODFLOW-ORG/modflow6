module BudgetFileReaderModule

  use KindModule
  use SimModule, only: store_error, store_error_unit
  use ConstantsModule, only: LINELENGTH, LENHUGELINE
  use BinaryFileReaderModule, only: BinaryFileReaderType, BinaryFileHeaderType

  implicit none

  private
  public :: BudgetFileReaderType, BudgetFileHeaderType

  type, extends(BinaryFileHeaderType) :: BudgetFileHeaderType
    character(len=16) :: budtxt
    integer(I4B) :: nval, idum1, idum2, imeth
    real(DP) :: delt
    character(len=16) :: srcmodelname, srcpackagename
    character(len=16) :: dstmodelname, dstpackagename
    integer(I4B) :: ndat, naux, nlist
    character(len=16), dimension(:), allocatable :: auxtxt
  contains
    procedure :: get_str
  end type BudgetFileHeaderType

  type, extends(BinaryFileReaderType) :: BudgetFileReaderType
    logical :: hasimeth1flowja = .false.
    integer(I4B) :: nbudterms
    character(len=16), dimension(:), allocatable :: budtxtarray
    integer(I4B), dimension(:), allocatable :: imetharray
    integer(I4B), dimension(:), allocatable :: nauxarray
    character(len=16), dimension(:, :), allocatable :: auxtxtarray
    real(DP), dimension(:), allocatable :: flowja
    integer(I4B), dimension(:), allocatable :: nodesrc
    integer(I4B), dimension(:), allocatable :: nodedst
    real(DP), dimension(:), allocatable :: flow
    real(DP), dimension(:, :), allocatable :: auxvar
    character(len=16), dimension(:), allocatable :: dstpackagenamearray
  contains
    procedure :: initialize
    procedure :: read_header
    procedure :: read_record
    procedure :: rewind
    procedure :: finalize
  end type BudgetFileReaderType

contains

  !< @brief initialize
  !<
  subroutine initialize(this, iu, iout, ncrbud)
    ! -- dummy
    class(BudgetFileReaderType) :: this
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(out) :: ncrbud
    ! -- local
    integer(I4B) :: ibudterm
    integer(I4B) :: kstp_last, kper_last
    integer(I4B) :: maxaux
    logical :: success
    !
    this%inunit = iu
    this%nbudterms = 0
    ncrbud = 0
    maxaux = 0
    call this%rewind()
    !
    ! -- Determine number of budget terms within a time step
    if (iout > 0) &
      write (iout, '(a)') &
      'Reading budget file to determine number of terms per time step.'
    !
    ! -- Read through the first set of data for time step 1 and stress period 1
    do
      call this%read_record(success)
      if (.not. success) exit
      this%nbudterms = this%nbudterms + 1
      select type (h => this%header)
      type is (BudgetFileHeaderType)
        if (h%naux > maxaux) maxaux = h%naux
      end select

      if (this%endoffile) exit
      if (this%header%kstp /= this%headernext%kstp .or. &
          this%header%kper /= this%headernext%kper) &
        exit
    end do
    kstp_last = this%header%kstp
    kper_last = this%header%kper
    allocate (this%budtxtarray(this%nbudterms))
    allocate (this%imetharray(this%nbudterms))
    allocate (this%dstpackagenamearray(this%nbudterms))
    allocate (this%nauxarray(this%nbudterms))
    allocate (this%auxtxtarray(maxaux, this%nbudterms))
    this%auxtxtarray(:, :) = ''
    call this%rewind()
    !
    ! -- Now read through again and store budget text names
    do ibudterm = 1, this%nbudterms
      call this%read_record(success, iout)
      if (.not. success) exit
      select type (h => this%header)
      type is (BudgetFileHeaderType)
        this%budtxtarray(ibudterm) = h%budtxt
        this%imetharray(ibudterm) = h%imeth
        this%dstpackagenamearray(ibudterm) = h%dstpackagename
        this%nauxarray(ibudterm) = h%naux
        if (h%naux > 0) then
          this%auxtxtarray(1:h%naux, ibudterm) = h%auxtxt(:)
        end if
        if (h%srcmodelname == h%dstmodelname) then
          if (allocated(this%nodesrc)) ncrbud = max(ncrbud, maxval(this%nodesrc))
        end if
      end select
    end do
    call this%rewind()
    if (iout > 0) &
      write (iout, '(a, i0, a)') 'Detected ', this%nbudterms, &
      ' unique flow terms in budget file.'
  end subroutine initialize

  !< @brief read header only
  !<
  subroutine read_header(this, success, iout)
    ! -- dummy
    class(BudgetFileReaderType), intent(inout) :: this
    logical, intent(out) :: success
    integer(I4B), intent(in), optional :: iout
    ! -- local
    integer(I4B) :: iostat
    character(len=LINELENGTH) :: errmsg
    !
    success = .true.
    select type (h => this%header)
    type is (BudgetFileHeaderType)
      h%kstp = 0
      h%kper = 0
      h%budtxt = ''
      h%nval = 0
      h%idum1 = 0
      h%idum2 = 0
      h%imeth = 0
      h%srcmodelname = ''
      h%srcpackagename = ''
      h%dstmodelname = ''
      h%dstpackagename = ''
      h%ndat = 0
      h%naux = 0
      h%nlist = 0
      if (allocated(h%auxtxt)) deallocate (h%auxtxt)

      inquire (unit=this%inunit, pos=h%pos)
      read (this%inunit, iostat=iostat) h%kstp, h%kper, &
        h%budtxt, h%nval, h%idum1, h%idum2
      if (iostat /= 0) then
        success = .false.
        if (iostat < 0) this%endoffile = .true.
        return
      end if
      read (this%inunit) h%imeth, h%delt, h%pertim, h%totim
      if (h%imeth == 6) then
        read (this%inunit) h%srcmodelname
        read (this%inunit) h%srcpackagename
        read (this%inunit) h%dstmodelname
        read (this%inunit) h%dstpackagename
        read (this%inunit) h%ndat
        h%naux = h%ndat - 1
        if (allocated(h%auxtxt)) deallocate (h%auxtxt)
        allocate (h%auxtxt(h%naux))
        read (this%inunit) h%auxtxt
        read (this%inunit) h%nlist
      elseif (h%imeth /= 1) then
        write (errmsg, '(a, a)') 'ERROR READING: ', trim(h%budtxt)
        call store_error(errmsg)
        write (errmsg, '(a, i0)') 'INVALID METHOD CODE DETECTED: ', h%imeth
        call store_error(errmsg)
        call store_error_unit(this%inunit)
      end if
    end select
  end subroutine read_header

  !< @brief read record
  !<
  subroutine read_record(this, success, iout)
    ! -- modules
    use InputOutputModule, only: fseek_stream
    ! -- dummy
    class(BudgetFileReaderType), intent(inout) :: this
    logical, intent(out) :: success
    integer(I4B), intent(in), optional :: iout
    ! -- local
    integer(I4B) :: i, n, iout_opt
    !
    if (present(iout)) then
      iout_opt = iout
    else
      iout_opt = 0
    end if
    !
    call this%read_header(success, iout_opt)
    if (.not. success) return
    !
    select type (h => this%header)
    type is (BudgetFileHeaderType)
      if (h%imeth == 1) then
        if (trim(adjustl(h%budtxt)) == 'FLOW-JA-FACE') then
          if (allocated(this%flowja)) deallocate (this%flowja)
          allocate (this%flowja(h%nval))
          read (this%inunit) this%flowja
          this%hasimeth1flowja = .true.
        else
          h%nval = h%nval * h%idum1 * abs(h%idum2)
          if (allocated(this%flow)) deallocate (this%flow)
          allocate (this%flow(h%nval))
          if (allocated(this%nodesrc)) deallocate (this%nodesrc)
          allocate (this%nodesrc(h%nval))
          read (this%inunit) this%flow
          do i = 1, h%nval
            this%nodesrc(i) = i
          end do
        end if
      elseif (h%imeth == 6) then
        if (allocated(this%nodesrc)) deallocate (this%nodesrc)
        allocate (this%nodesrc(h%nlist))
        if (allocated(this%nodedst)) deallocate (this%nodedst)
        allocate (this%nodedst(h%nlist))
        if (allocated(this%flow)) deallocate (this%flow)
        allocate (this%flow(h%nlist))
        if (allocated(this%auxvar)) deallocate (this%auxvar)
        allocate (this%auxvar(h%naux, h%nlist))
        read (this%inunit) (this%nodesrc(n), this%nodedst(n), this%flow(n), &
                            (this%auxvar(i, n), i=1, h%naux), n=1, h%nlist)
      end if

      if (iout_opt > 0) then
        write (iout_opt, '(1pg15.6, a, 1x, a)') h%totim, h%budtxt, &
          h%dstpackagename
      end if
    end select
    !
    call this%peek_record()
  end subroutine read_record

  !< @brief finalize
  !<
  subroutine finalize(this)
    class(BudgetFileReaderType) :: this
    close (this%inunit)
    if (allocated(this%flowja)) deallocate (this%flowja)
    if (allocated(this%nodesrc)) deallocate (this%nodesrc)
    if (allocated(this%nodedst)) deallocate (this%nodedst)
    if (allocated(this%flow)) deallocate (this%flow)
    if (allocated(this%auxvar)) deallocate (this%auxvar)
    if (allocated(this%header)) deallocate (this%header)
    if (allocated(this%headernext)) deallocate (this%headernext)
  end subroutine finalize

  !> @brief Get a string representation of the budget file header.
  function get_str(this) result(str)
    class(BudgetFileHeaderType), intent(in) :: this
    character(len=:), allocatable :: str
    character(len=LENHUGELINE) :: temp

    write (temp, '(*(G0))') &
      'Budget file header (pos: ', this%pos, &
      ', kper: ', this%kper, &
      ', kstp: ', this%kstp, &
      ', delt: ', this%delt, &
      ', pertim: ', this%pertim, &
      ', totim: ', this%totim, &
      ', budtxt: ', trim(this%budtxt), &
      ', nval: ', this%nval, &
      ', idum1: ', this%idum1, &
      ', idum2: ', this%idum2, &
      ', imeth: ', this%imeth, &
      ', srcmodel: ', trim(this%srcmodelname), &
      ', srcpackage: ', trim(this%srcpackagename), &
      ', dstmodel: ', trim(this%dstmodelname), &
      ', dstpackage: ', trim(this%dstpackagename), &
      ', ndat: ', this%ndat, &
      ', naux: ', this%naux, &
      ', nlist: ', this%nlist, &
      ')'
    str = trim(temp)
  end function get_str

  subroutine rewind (this)
    class(BudgetFileReaderType), intent(inout) :: this

    rewind (this%inunit)
    this%endoffile = .false.
    if (allocated(this%header)) deallocate (this%header)
    if (allocated(this%headernext)) deallocate (this%headernext)
    allocate (BudgetFileHeaderType :: this%header)
    allocate (BudgetFileHeaderType :: this%headernext)
    this%header%pos = 1
    this%headernext%pos = 1
  end subroutine rewind

end module BudgetFileReaderModule
