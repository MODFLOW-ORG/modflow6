module PrtOcModule

  use BaseDisModule, only: DisBaseType
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENMODELNAME, MNORMAL
  use OutputControlModule, only: OutputControlType
  use OutputControlDataModule, only: OutputControlDataType, ocd_cr
  use SimModule, only: store_error, store_error_filename
  use SimVariablesModule, only: errmsg, warnmsg
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, mem_reallocate
  use MemoryHelperModule, only: create_mem_path
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
    logical(LGP), pointer :: trackfeatexit => null() !< whether to track grid-scale feature exit events
    logical(LGP), pointer :: tracktimestep => null() !< whether to track timestep events
    logical(LGP), pointer :: trackterminate => null() !< whether to track termination events
    logical(LGP), pointer :: trackweaksink => null() !< whether to track weak sink exit events
    logical(LGP), pointer :: trackusertime => null() !< whether to track user-specified times
    logical(LGP), pointer :: tracksubfexit => null() !< whether to track sub-grid-scale feature exit events
    logical(LGP), pointer :: trackdropped => null() !< whether to track drops to water table
    integer(I4B), pointer :: ntracktimes => null() !< number of user-specified tracking times
    logical(LGP), pointer :: dump_event_trace => null() !< whether to dump event trace for debugging
    type(TimeSelectType), pointer :: tracktimes !< user-specified tracking times

  contains
    procedure :: oc_ar
    procedure :: oc_da => prt_oc_da
    procedure :: allocate_scalars => prt_oc_allocate_scalars
    procedure :: source_options => prt_oc_source_options
    procedure, private :: prt_oc_source_dimensions
    procedure, private :: prt_oc_source_tracktimes

  end type PrtOcType

contains

  !> @ brief Create an output control object
  subroutine oc_cr(ocobj, name_model, input_mempath, inunit, iout)
    type(PrtOcType), pointer :: ocobj !< PrtOcType object
    character(len=*), intent(in) :: name_model !< name of the model
    character(len=*), intent(in) :: input_mempath !< input mempath of the package
    integer(I4B), intent(in) :: inunit !< unit number for input
    integer(I4B), intent(in) :: iout !< unit number for output

    ! Create the object
    allocate (ocobj)

    ! Allocate scalars
    call ocobj%allocate_scalars(name_model, input_mempath)

    ! Save unit numbers
    ocobj%inunit = inunit
    ocobj%iout = iout
  end subroutine oc_cr

  subroutine prt_oc_allocate_scalars(this, name_model, input_mempath)
    use MemoryManagerExtModule, only: mem_set_value
    class(PrtOcType) :: this
    character(len=*), intent(in) :: name_model !< name of model
    character(len=*), intent(in) :: input_mempath !< input mempath of the package
    logical(LGP) :: found

    this%memoryPath = create_mem_path(name_model, 'OC')

    allocate (this%name_model)
    allocate (this%input_fname)
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
    call mem_allocate(this%trackrelease, 'ITRACKRELEASE', this%memoryPath)
    call mem_allocate(this%trackfeatexit, 'ITRACKFEATEXIT', this%memoryPath)
    call mem_allocate(this%tracktimestep, 'ITRACKTIMESTEP', this%memoryPath)
    call mem_allocate(this%trackterminate, 'ITRACKTERMINATE', this%memoryPath)
    call mem_allocate(this%trackweaksink, 'ITRACKWEAKSINK', this%memoryPath)
    call mem_allocate(this%trackusertime, 'ITRACKUSERTIME', this%memoryPath)
    call mem_allocate(this%tracksubfexit, 'ITRACKSUBFEXIT', this%memoryPath)
    call mem_allocate(this%trackdropped, 'ITRACKDROPPED', this%memoryPath)
    call mem_allocate(this%ntracktimes, 'NTRACKTIMES', this%memoryPath)

    this%name_model = name_model
    this%input_mempath = input_mempath
    this%input_fname = ''
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
    this%trackfeatexit = .false.
    this%tracktimestep = .false.
    this%trackterminate = .false.
    this%trackweaksink = .false.
    this%trackusertime = .false.
    this%tracksubfexit = .false.
    this%trackdropped = .false.
    this%ntracktimes = 0

    if (this%input_mempath /= '') then
      call mem_set_value(this%input_fname, 'INPUT_FNAME', &
                         this%input_mempath, found)
    end if
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
    if (this%input_mempath == '') return
    call this%source_options()
    call this%prt_oc_source_dimensions()
    call this%prt_oc_source_tracktimes()
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
    call mem_deallocate(this%trackfeatexit)
    call mem_deallocate(this%tracktimestep)
    call mem_deallocate(this%trackterminate)
    call mem_deallocate(this%trackweaksink)
    call mem_deallocate(this%trackusertime)
    call mem_deallocate(this%tracksubfexit)
    call mem_deallocate(this%trackdropped)
    call mem_deallocate(this%ntracktimes)

  end subroutine prt_oc_da

  subroutine prt_oc_source_options(this)
    ! -- modules
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: getunit, openfile
    use ConstantsModule, only: LINELENGTH
    use ParticleTracksModule, only: TRACKHEADER, TRACKDTYPES
    use MemoryManagerExtModule, only: mem_set_value
    use PrtOcInputModule, only: PrtOcParamFoundType
    ! -- dummy
    class(PrtOcType) :: this
    ! -- local
    character(len=LINELENGTH) :: trackfile, trackcsv
    type(PrtOcParamFoundType) :: found
    integer(I4B), pointer :: evinput
    ! formats
    character(len=*), parameter :: fmttrkbin = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO BINARY FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"
    character(len=*), parameter :: fmttrkcsv = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO CSV FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"

    allocate (evinput)

    write (this%iout, '(/,1x,a,/)') 'PROCESSING OC OPTIONS'
    !
    ! -- source base class options
    call this%OutPutControlType%source_options()
    !
    ! -- source options
    call mem_set_value(trackfile, 'TRACKFILE', this%input_mempath, &
                       found%trackfile)
    call mem_set_value(trackcsv, 'TRACKCSVFILE', this%input_mempath, &
                       found%trackcsvfile)
    call mem_set_value(evinput, 'TRACK_RELEASE', this%input_mempath, &
                       found%track_release)
    call mem_set_value(evinput, 'TRACK_EXIT', this%input_mempath, &
                       found%track_exit)
    call mem_set_value(evinput, 'TRACK_SUBF_EXIT', this%input_mempath, &
                       found%track_subf_exit)
    call mem_set_value(evinput, 'TRACK_DROPPED', this%input_mempath, &
                       found%track_dropped)
    call mem_set_value(evinput, 'TRACK_TIMESTEP', this%input_mempath, &
                       found%track_timestep)
    call mem_set_value(evinput, 'TRACK_TERMINATE', this%input_mempath, &
                       found%track_terminate)
    call mem_set_value(evinput, 'TRACK_WEAKSINK', this%input_mempath, &
                       found%track_weaksink)
    call mem_set_value(evinput, 'TRACK_USERTIME', this%input_mempath, &
                       found%track_usertime)
    call mem_set_value(evinput, 'DEV_DUMP_EVTRACE', this%input_mempath, &
                       found%dev_dump_evtrace)

    if (found%track_release) this%trackrelease = .true.
    if (found%track_exit) this%trackfeatexit = .true.
    if (found%track_subf_exit) this%tracksubfexit = .true.
    if (found%track_dropped) this%trackdropped = .true.
    if (found%track_timestep) this%tracktimestep = .true.
    if (found%track_terminate) this%trackterminate = .true.
    if (found%track_weaksink) this%trackweaksink = .true.
    if (found%track_usertime) this%trackusertime = .true.
    if (found%dev_dump_evtrace) this%dump_event_trace = .true.

    ! default to all events
    if (.not. (found%track_release .or. &
               found%track_exit .or. &
               found%track_timestep .or. &
               found%track_terminate .or. &
               found%track_weaksink .or. &
               found%track_usertime .or. &
               found%track_dropped)) then
      this%trackrelease = .true.
      this%trackfeatexit = .true.
      this%tracktimestep = .true.
      this%trackterminate = .true.
      this%trackweaksink = .true.
      this%trackusertime = .true.
      this%trackdropped = .true.
    end if

    if (found%trackfile) then
      ! open binary track output file
      this%itrkout = getunit()
      call openfile(this%itrkout, this%iout, trackfile, 'DATA(BINARY)', &
                    form, access, filstat_opt='REPLACE', &
                    mode_opt=MNORMAL)
      write (this%iout, fmttrkbin) trim(adjustl(trackfile)), this%itrkout
      ! open and write ascii track header file
      this%itrkhdr = getunit()
      trackfile = trim(trackfile)//'.hdr'
      call openfile(this%itrkhdr, this%iout, trackfile, 'CSV', &
                    filstat_opt='REPLACE', mode_opt=MNORMAL)
      write (this%itrkhdr, '(a,/,a)') TRACKHEADER, TRACKDTYPES
    end if

    if (found%trackcsvfile) then
      this%itrkcsv = getunit()
      call openfile(this%itrkcsv, this%iout, trackcsv, 'CSV', &
                    filstat_opt='REPLACE')
      write (this%iout, fmttrkcsv) trim(adjustl(trackcsv)), this%itrkcsv
      write (this%itrkcsv, '(a)') TRACKHEADER
    end if

    write (this%iout, '(1x,a)') 'END OF OC OPTIONS'
    deallocate (evinput)
  end subroutine prt_oc_source_options

  !> @brief source the dimensions block.
  subroutine prt_oc_source_dimensions(this)
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error, count_errors
    use MemoryManagerExtModule, only: mem_set_value
    use PrtOcInputModule, only: PrtOcParamFoundType
    ! dummy
    class(PrtOcType), intent(inout) :: this
    ! local
    type(PrtOcParamFoundType) :: found
    write (this%iout, '(/1x,a)') &
      'PROCESSING OUTPUT CONTROL DIMENSIONS'
    call mem_set_value(this%ntracktimes, 'NTRACKTIMES', this%input_mempath, &
                       found%ntracktimes)
    if (found%ntracktimes) then
      write (this%iout, '(4x,a,i7)') 'NTRACKTIMES = ', this%ntracktimes
    end if
    write (this%iout, '(1x,a)') &
      'END OF OUTPUT CONTROL DIMENSIONS'

    if (this%ntracktimes < 0) then
      write (errmsg, '(a)') &
        'NTRACKTIMES WAS NOT SPECIFIED OR WAS SPECIFIED INCORRECTLY.'
      call store_error(errmsg)
    end if
  end subroutine prt_oc_source_dimensions

  !> @brief source the tracking times block.
  subroutine prt_oc_source_tracktimes(this)
    use ParticleTracksModule, only: TRACKHEADER, TRACKDTYPES
    use MemoryManagerExtModule, only: mem_set_value
    use MemoryManagerModule, only: mem_setptr, get_isize
    ! dummy
    class(PrtOcType), intent(inout) :: this
    ! local
    real(DP), dimension(:), pointer, contiguous :: tracktimes
    integer(I4B) :: n, asize

    if (this%ntracktimes <= 0) return

    call get_isize('TIME', this%input_mempath, asize)

    if (asize /= this%ntracktimes) then
      write (errmsg, '(a, i0)') &
        "Expected TRACKTIMES with length ", this%ntracktimes
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    else
      call mem_setptr(tracktimes, 'TIME', this%input_mempath)

      ! allocate time selection
      call this%tracktimes%expand(this%ntracktimes)

      do n = 1, this%ntracktimes
        this%tracktimes%times(n) = tracktimes(n)
      end do
    end if

    ! make sure times strictly increase
    if (.not. this%tracktimes%increasing()) then
      errmsg = "TRACKTIMES must strictly increase"
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if
  end subroutine prt_oc_source_tracktimes

end module PrtOcModule
