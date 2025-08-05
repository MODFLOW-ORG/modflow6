module PrtOcModule

  use BaseDisModule, only: DisBaseType
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENMODELNAME, MNORMAL
  use OutputControlModule, only: OutputControlType
  use OutputControlDataModule, only: OutputControlDataType, ocd_cr
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg, warnmsg
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, mem_reallocate
  use MemoryHelperModule, only: create_mem_path
  use BlockParserModule, only: BlockParserType
  use InputOutputModule, only: urword, openfile
  use TimeSelectModule, only: TimeSelectType
  use LongLineReaderModule, only: LongLineReaderType

  implicit none
  private
  public PrtOcType, oc_cr

  !> @ brief Output control for particle tracking models
  type, extends(OutputControlType) :: PrtOcType
    integer(I4B), pointer :: itrkout => null() !< binary output file
    integer(I4B), pointer :: itrkhdr => null() !< output header file
    integer(I4B), pointer :: itrkcsv => null() !< CSV output file
    integer(I4B), pointer :: itrktls => null() !< track time list input file
    logical(LGP), pointer :: trackrelease => null() !< whether to track release events
    logical(LGP), pointer :: trackcellexit => null() !< whether to track cell exit events
    logical(LGP), pointer :: tracktimestep => null() !< whether to track timestep events
    logical(LGP), pointer :: trackterminate => null() !< whether to track termination events
    logical(LGP), pointer :: trackweaksink => null() !< whether to track weak sink exit events
    logical(LGP), pointer :: trackusertime => null() !< whether to track user-specified times
    integer(I4B), pointer :: ntracktimes => null() !< number of user-specified tracking times
    logical(LGP), pointer :: dump_event_trace => null() !< whether to dump event trace for debugging
    type(TimeSelectType), pointer :: tracktimes !< user-specified tracking times

  contains
    procedure :: oc_ar
    procedure :: oc_da => prt_oc_da
    procedure :: allocate_scalars => prt_oc_allocate_scalars
    procedure :: read_options => prt_oc_read_options
    procedure, private :: prt_oc_read_dimensions
    procedure, private :: prt_oc_read_tracktimes

  end type PrtOcType

contains

  !> @ brief Create an output control object
  subroutine oc_cr(ocobj, name_model, inunit, iout)
    type(PrtOcType), pointer :: ocobj !< PrtOcType object
    character(len=*), intent(in) :: name_model !< name of the model
    integer(I4B), intent(in) :: inunit !< unit number for input
    integer(I4B), intent(in) :: iout !< unit number for output

    ! Create the object
    allocate (ocobj)

    ! Allocate scalars
    call ocobj%allocate_scalars(name_model)

    ! Save unit numbers
    ocobj%inunit = inunit
    ocobj%iout = iout

    ! Initialize block parser
    call ocobj%parser%Initialize(inunit, iout)
  end subroutine oc_cr

  subroutine prt_oc_allocate_scalars(this, name_model)
    class(PrtOcType) :: this
    character(len=*), intent(in) :: name_model !< name of model

    this%memoryPath = create_mem_path(name_model, 'OC')

    allocate (this%name_model)
    call mem_allocate(this%dump_event_trace, 'DUMP_EVENT_TRACE', this%memoryPath)
    call mem_allocate(this%inunit, 'INUNIT', this%memoryPath)
    call mem_allocate(this%iout, 'IOUT', this%memoryPath)
    call mem_allocate(this%ibudcsv, 'IBUDCSV', this%memoryPath)
    call mem_allocate(this%iperoc, 'IPEROC', this%memoryPath)
    call mem_allocate(this%iocrep, 'IOCREP', this%memoryPath)
    call mem_allocate(this%itrkout, 'ITRKOUT', this%memoryPath)
    call mem_allocate(this%itrkhdr, 'ITRKHDR', this%memoryPath)
    call mem_allocate(this%itrkcsv, 'ITRKCSV', this%memoryPath)
    call mem_allocate(this%itrktls, 'ITRKTLS', this%memoryPath)
    call mem_allocate(this%trackrelease, 'ITRACKRLS', this%memoryPath)
    call mem_allocate(this%trackcellexit, 'ITRACKTRS', this%memoryPath)
    call mem_allocate(this%tracktimestep, 'ITRACKTST', this%memoryPath)
    call mem_allocate(this%trackterminate, 'ITRACKTER', this%memoryPath)
    call mem_allocate(this%trackweaksink, 'ITRACKWSK', this%memoryPath)
    call mem_allocate(this%trackusertime, 'ITRACKTLS', this%memoryPath)
    call mem_allocate(this%ntracktimes, 'NTRACKTIMES', this%memoryPath)

    this%name_model = name_model
    this%dump_event_trace = .false.
    this%inunit = 0
    this%iout = 0
    this%ibudcsv = 0
    this%iperoc = 0
    this%iocrep = 0
    this%itrkout = 0
    this%itrkhdr = 0
    this%itrkcsv = 0
    this%itrktls = 0
    this%trackrelease = .false.
    this%trackcellexit = .false.
    this%tracktimestep = .false.
    this%trackterminate = .false.
    this%trackweaksink = .false.
    this%trackusertime = .false.
    this%ntracktimes = 0

  end subroutine prt_oc_allocate_scalars

  !> @ brief Setup output control variables.
  subroutine oc_ar(this, dis, dnodata)
    ! dummy
    class(PrtOcType) :: this !< PrtOcType object
    class(DisBaseType), pointer, intent(in) :: dis !< model discretization package
    real(DP), intent(in) :: dnodata !< no data value
    ! local
    integer(I4B) :: i, nocdobj, inodata
    type(OutputControlDataType), pointer :: ocdobjptr
    real(DP), dimension(:), pointer, contiguous :: nullvec => null()

    ! Allocate and initialize variables
    allocate (this%tracktimes)
    call this%tracktimes%init()
    inodata = 0
    nocdobj = 1
    allocate (this%ocds(nocdobj))
    do i = 1, nocdobj
      call ocd_cr(ocdobjptr)
      select case (i)
      case (1)
        call ocdobjptr%init_dbl('BUDGET', nullvec, dis, 'PRINT LAST ', &
                                'COLUMNS 10 WIDTH 11 DIGITS 4 GENERAL ', &
                                this%iout, dnodata)
      end select
      this%ocds(i) = ocdobjptr
      deallocate (ocdobjptr)
    end do

    ! Read options, dimensions, and tracktimes
    ! blocks if this package is enabled
    if (this%inunit <= 0) return
    call this%read_options()
    call this%prt_oc_read_dimensions()
    call this%prt_oc_read_tracktimes()

  end subroutine oc_ar

  subroutine prt_oc_da(this)
    ! dummy
    class(PrtOcType) :: this
    ! local
    integer(I4B) :: i

    call this%tracktimes%deallocate()

    do i = 1, size(this%ocds)
      call this%ocds(i)%ocd_da()
    end do
    deallocate (this%ocds)

    deallocate (this%name_model)
    call mem_deallocate(this%dump_event_trace)
    call mem_deallocate(this%inunit)
    call mem_deallocate(this%iout)
    call mem_deallocate(this%ibudcsv)
    call mem_deallocate(this%iperoc)
    call mem_deallocate(this%iocrep)
    call mem_deallocate(this%itrkout)
    call mem_deallocate(this%itrkhdr)
    call mem_deallocate(this%itrkcsv)
    call mem_deallocate(this%itrktls)
    call mem_deallocate(this%trackrelease)
    call mem_deallocate(this%trackcellexit)
    call mem_deallocate(this%tracktimestep)
    call mem_deallocate(this%trackterminate)
    call mem_deallocate(this%trackweaksink)
    call mem_deallocate(this%trackusertime)
    call mem_deallocate(this%ntracktimes)

  end subroutine prt_oc_da

  subroutine prt_oc_read_options(this)
    ! modules
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: getunit, openfile, lowcase
    use ConstantsModule, only: LINELENGTH
    use ParticleTracksModule, only: TRACKHEADER, TRACKDTYPES
    use SimModule, only: store_error, store_error_unit
    use InputOutputModule, only: openfile, getunit
    ! dummy
    class(PrtOcType) :: this
    ! local
    character(len=LINELENGTH) :: keyword
    character(len=LINELENGTH) :: keyword2
    character(len=LINELENGTH) :: fname
    character(len=:), allocatable :: line
    integer(I4B) :: ierr, ipos
    logical(LGP) :: block_found, param_found, event_found, eob
    type(OutputControlDataType), pointer :: ocdobjptr
    ! formats
    character(len=*), parameter :: fmttrkbin = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO BINARY FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"
    character(len=*), parameter :: fmttrkcsv = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO CSV FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"

    ! get options block
    call this%parser%GetBlock('OPTIONS', block_found, ierr, &
                              supportOpenClose=.true., blockRequired=.false.)

    ! parse options block if detected
    if (block_found) then
      write (this%iout, '(/,1x,a,/)') 'PROCESSING OC OPTIONS'
      event_found = .false.
      do
        call this%parser%GetNextLine(eob)
        if (eob) exit
        call this%parser%GetStringCaps(keyword)
        param_found = .false.
        select case (keyword)
        case ('BUDGETCSV')
          call this%parser%GetStringCaps(keyword2)
          if (keyword2 /= 'FILEOUT') then
            errmsg = "BUDGETCSV must be followed by FILEOUT and then budget &
              &csv file name.  Found '"//trim(keyword2)//"'."
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetString(fname)
          this%ibudcsv = GetUnit()
          call openfile(this%ibudcsv, this%iout, fname, 'CSV', &
                        filstat_opt='REPLACE')
          param_found = .true.
        case ('TRACK')
          call this%parser%GetStringCaps(keyword)
          if (keyword == 'FILEOUT') then
            ! parse filename
            call this%parser%GetString(fname)
            ! open binary track output file
            this%itrkout = getunit()
            call openfile(this%itrkout, this%iout, fname, 'DATA(BINARY)', &
                          form, access, filstat_opt='REPLACE', &
                          mode_opt=MNORMAL)
            write (this%iout, fmttrkbin) trim(adjustl(fname)), this%itrkout
            ! open and write ascii track header file
            this%itrkhdr = getunit()
            fname = trim(fname)//'.hdr'
            call openfile(this%itrkhdr, this%iout, fname, 'CSV', &
                          filstat_opt='REPLACE', mode_opt=MNORMAL)
            write (this%itrkhdr, '(a,/,a)') TRACKHEADER, TRACKDTYPES
          else
            call store_error('OPTIONAL TRACK KEYWORD MUST BE '// &
                             'FOLLOWED BY FILEOUT')
          end if
          param_found = .true.
        case ('TRACKCSV')
          call this%parser%GetStringCaps(keyword)
          if (keyword == 'FILEOUT') then
            ! parse filename
            call this%parser%GetString(fname)
            ! open CSV track output file and write headers
            this%itrkcsv = getunit()
            call openfile(this%itrkcsv, this%iout, fname, 'CSV', &
                          filstat_opt='REPLACE')
            write (this%iout, fmttrkcsv) trim(adjustl(fname)), this%itrkcsv
            write (this%itrkcsv, '(a)') TRACKHEADER
          else
            call store_error('OPTIONAL TRACKCSV KEYWORD MUST BE &
              &FOLLOWED BY FILEOUT')
          end if
          param_found = .true.
        case ('TRACK_RELEASE')
          this%trackrelease = .true.
          event_found = .true.
          param_found = .true.
        case ('TRACK_EXIT')
          this%trackcellexit = .true.
          event_found = .true.
          param_found = .true.
        case ('TRACK_TIMESTEP')
          this%tracktimestep = .true.
          event_found = .true.
          param_found = .true.
        case ('TRACK_TERMINATE')
          this%trackterminate = .true.
          event_found = .true.
          param_found = .true.
        case ('TRACK_WEAKSINK')
          this%trackweaksink = .true.
          event_found = .true.
          param_found = .true.
        case ('TRACK_USERTIME')
          this%trackusertime = .true.
          event_found = .true.
          param_found = .true.
        case ('DEV_DUMP_EVENT_TRACE')
          this%dump_event_trace = .true.
          param_found = .true.
        case default
          param_found = .false.
        end select

        ! check if we're done
        if (.not. param_found) then
          do ipos = 1, size(this%ocds)
            ocdobjptr => this%ocds(ipos)
            if (keyword == trim(ocdobjptr%cname)) then
              param_found = .true.
              exit
            end if
          end do
          if (.not. param_found) then
            errmsg = "UNKNOWN OC OPTION '"//trim(keyword)//"'."
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetRemainingLine(line)
          call ocdobjptr%set_option(line, this%parser%iuactive, this%iout)
        end if
      end do

      ! default to all events
      if (.not. event_found) then
        this%trackrelease = .true.
        this%trackcellexit = .true.
        this%tracktimestep = .true.
        this%trackterminate = .true.
        this%trackweaksink = .true.
        this%trackusertime = .true.
      end if

      ! logging
      write (this%iout, '(1x,a)') 'END OF OC OPTIONS'
    end if
  end subroutine prt_oc_read_options

  !> @brief Read the dimensions block.
  subroutine prt_oc_read_dimensions(this)
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error, count_errors
    ! dummy
    class(PrtOcType), intent(inout) :: this
    ! local
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    logical(LGP) :: isfound, endOfBlock

    ! initialize dimensions to -1
    this%ntracktimes = -1

    ! get dimensions block
    call this%parser%GetBlock('DIMENSIONS', isfound, ierr, &
                              supportOpenClose=.true., &
                              blockRequired=.false.)

    ! parse dimensions block if detected
    if (.not. isfound) return
    write (this%iout, '(/1x,a)') &
      'PROCESSING OUTPUT CONTROL DIMENSIONS'
    do
      call this%parser%GetNextLine(endOfBlock)
      if (endOfBlock) exit
      call this%parser%GetStringCaps(keyword)
      select case (keyword)
      case ('NTRACKTIMES')
        this%ntracktimes = this%parser%GetInteger()
        write (this%iout, '(4x,a,i7)') 'NTRACKTIMES = ', this%ntracktimes
      case default
        write (errmsg, '(a,a)') &
          'UNKNOWN OUTPUT CONTROL DIMENSION: ', trim(keyword)
        call store_error(errmsg)
      end select
    end do
    write (this%iout, '(1x,a)') &
      'END OF OUTPUT CONTROL DIMENSIONS'

    if (this%ntracktimes < 0) then
      write (errmsg, '(a)') &
        'NTRACKTIMES WAS NOT SPECIFIED OR WAS SPECIFIED INCORRECTLY.'
      call store_error(errmsg)
    end if

    ! stop if errors were encountered in the block
    if (count_errors() > 0) &
      call this%parser%StoreErrorUnit()

  end subroutine prt_oc_read_dimensions

  !> @brief Read the tracking times block.
  subroutine prt_oc_read_tracktimes(this)
    ! dummy
    class(PrtOcType), intent(inout) :: this
    ! local
    integer(I4B) :: i, ierr
    logical(LGP) :: eob, found, success
    real(DP) :: t

    ! get tracktimes block
    call this%parser%GetBlock('TRACKTIMES', found, ierr, &
                              supportOpenClose=.true., &
                              blockRequired=.false.)

    ! raise an error if tracktimes has a dimension
    ! but no block was found, otherwise return early
    if (.not. found) then
      if (this%ntracktimes <= 0) return
      write (errmsg, '(a, i0)') &
        "Expected TRACKTIMES with length ", this%ntracktimes
      call store_error(errmsg)
      call this%parser%StoreErrorUnit(terminate=.true.)
    end if

    ! allocate time selection
    call this%tracktimes%expand(this%ntracktimes)

    ! read the block
    write (this%iout, '(/1x,a)') &
      'PROCESSING OUTPUT CONTROL TRACKTIMES'
    do i = 1, this%ntracktimes
      call this%parser%GetNextLine(eob)
      if (eob) exit
      call this%parser%TryGetDouble(t, success)
      if (.not. success) then
        errmsg = "Failed to read double precision value"
        call store_error(errmsg)
        call this%parser%StoreErrorUnit(terminate=.true.)
      end if
      this%tracktimes%times(i) = t
    end do

    ! make sure times strictly increase
    if (.not. this%tracktimes%increasing()) then
      errmsg = "TRACKTIMES must strictly increase"
      call store_error(errmsg)
      call this%parser%StoreErrorUnit(terminate=.true.)
    end if

  end subroutine prt_oc_read_tracktimes

end module PrtOcModule
