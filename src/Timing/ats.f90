! Outstanding issues for future work:
!   CSUB state advance/restore
!   Add courant time step constraint and other stability controls for GWT model
module AdaptiveTimeStepModule

  use KindModule, only: DP, I4B, LGP
  use SimVariablesModule, only: iout, errmsg, warnmsg
  use SimModule, only: store_error, count_errors, store_warning
  use BlockParserModule, only: BlockParserType
  use ConstantsModule, only: DZERO, DONE, LINELENGTH, DNODATA, TABLEFT, &
                             TABCENTER

  implicit none
  private

  public :: AtsType

  type :: AtsType
    integer(I4B), pointer :: nper => null() !< set equal to nper
    integer(I4B), pointer :: maxats => null() !< number of ats entries
    real(DP), public, pointer :: dtstable => null() !< delt value required for stability
    integer(I4B), dimension(:), pointer, contiguous :: kperats => null() !< array of stress period numbers to apply ats (size NPER)
    integer(I4B), dimension(:), pointer, contiguous :: iperats => null() !< array of stress period numbers to apply ats (size MAXATS)
    real(DP), dimension(:), pointer, contiguous :: dt0 => null() !< input array of initial time step sizes
    real(DP), dimension(:), pointer, contiguous :: dtmin => null() !< input array of minimum time step sizes
    real(DP), dimension(:), pointer, contiguous :: dtmax => null() !< input array of maximum time step sizes
    real(DP), dimension(:), pointer, contiguous :: dtadj => null() !< input array of time step factors for shortening or increasing
    real(DP), dimension(:), pointer, contiguous :: dtfailadj => null() !< input array of time step factors for shortening due to nonconvergence
    type(BlockParserType) :: parser !< block parser for reading input file
  contains
    procedure, public :: ats_init
    procedure, public :: isAdaptivePeriod
    procedure, public :: ats_submit_delt
    procedure, public :: ats_set_delt
    procedure, public :: ats_reset_delt
    procedure, public :: ats_period_message
    procedure, public :: ats_set_endofperiod
    procedure, public :: ats_da
    ! private
    procedure, private :: ats_allocate_scalars
    procedure, private :: ats_allocate_arrays
    procedure, private :: ats_read_options
    procedure, private :: ats_read_dimensions
    procedure, private :: ats_read_timing
    procedure, private :: ats_process_input
    procedure, private :: ats_input_table
    procedure, private :: ats_check_timing
  end type AtsType

contains

  !> @ brief Determine if period is adaptive
  !!
  !!  Check settings and determine if kper is an adaptive
  !!  stress period.
  !!
  !<
  function isAdaptivePeriod(this, kper) result(lv)
    class(AtsType) :: this
    integer(I4B), intent(in) :: kper
    logical(LGP) :: lv
    lv = .false.
    if (associated(this%kperats)) then
      if (this%kperats(kper) > 0) then
        lv = .true.
      end if
    end if
  end function isAdaptivePeriod

  !> @ brief Create ATS object
  !!
  !!  Create a new ATS object, and read and check input.
  !!
  !<
  subroutine ats_init(this, inunit, nper_tdis)
    ! -- modules
    ! -- dummy
    class(AtsType) :: this
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: nper_tdis
    ! -- local
    ! -- formats
    character(len=*), parameter :: fmtheader = &
        "(1X,/1X,'ATS -- ADAPTIVE TIME STEP PACKAGE,',   /                  &
        &' VERSION 1 : 03/18/2021 - INPUT READ FROM UNIT ',I0)"
    !
    ! -- Allocate the scalar variables
    call this%ats_allocate_scalars()
    !
    ! -- Identify package
    write (iout, fmtheader) inunit
    !
    ! -- Initialize block parser
    call this%parser%initialize(inunit, iout)
    !
    ! -- Read options
    call this%ats_read_options()
    !
    ! -- store tdis nper in nper
    this%nper = nper_tdis
    !
    ! -- Read dimensions and then allocate arrays
    call this%ats_read_dimensions()
    call this%ats_allocate_arrays()
    !
    ! -- Read timing
    call this%ats_read_timing()
    !
    ! -- Echo input data to table
    call this%ats_input_table()
    !
    ! -- Check timing
    call this%ats_check_timing()
    !
    ! -- Process input
    call this%ats_process_input()
    !
    ! -- Close the file
    call this%parser%Clear()
  end subroutine ats_init

  !> @ brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the ATS package.
  !!
  !<
  subroutine ats_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(AtsType) :: this
    ! -- memory manager variables
    call mem_allocate(this%nper, 'NPER', 'ATS')
    call mem_allocate(this%maxats, 'MAXATS', 'ATS')
    call mem_allocate(this%dtstable, 'DTSTABLE', 'ATS')
    !
    ! -- Initialize variables
    this%nper = 0
    this%maxats = 0
    this%dtstable = DNODATA
  end subroutine ats_allocate_scalars

  !> @ brief Allocate arrays
  !!
  !! Allocate and initialize arrays for the ATS package.
  !!
  !<
  subroutine ats_allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(AtsType) :: this
    ! -- local
    integer(I4B) :: n
    !
    call mem_allocate(this%kperats, this%nper, 'KPERATS', 'ATS')
    call mem_allocate(this%iperats, this%maxats, 'IPERATS', 'ATS')
    call mem_allocate(this%dt0, this%maxats, 'DT0', 'ATS')
    call mem_allocate(this%dtmin, this%maxats, 'DTMIN', 'ATS')
    call mem_allocate(this%dtmax, this%maxats, 'DTMAX', 'ATS')
    call mem_allocate(this%dtadj, this%maxats, 'DTADJ', 'ATS')
    call mem_allocate(this%dtfailadj, this%maxats, 'DTFAILADJ', 'ATS')
    !
    ! -- initialize kperats
    do n = 1, this%nper
      this%kperats(n) = 0
    end do
    !
    ! -- initialize
    do n = 1, this%maxats
      this%iperats(n) = 0
      this%dt0(n) = DZERO
      this%dtmin(n) = DZERO
      this%dtmax(n) = DZERO
      this%dtadj(n) = DZERO
      this%dtfailadj(n) = DZERO
    end do
  end subroutine ats_allocate_arrays

  !> @ brief Deallocate variables
  !!
  !! Deallocate all ATS variables.
  !!
  !<
  subroutine ats_da(this)
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(AtsType) :: this
    ! -- Scalars
    call mem_deallocate(this%nper)
    call mem_deallocate(this%maxats)
    call mem_deallocate(this%dtstable)
    !
    ! -- Arrays
    call mem_deallocate(this%kperats)
    call mem_deallocate(this%iperats)
    call mem_deallocate(this%dt0)
    call mem_deallocate(this%dtmin)
    call mem_deallocate(this%dtmax)
    call mem_deallocate(this%dtadj)
    call mem_deallocate(this%dtfailadj)
  end subroutine ats_da

  !> @ brief Read options
  !!
  !! Read options from ATS input file.
  !!
  !<
  subroutine ats_read_options(this)
    ! -- dummy
    class(AtsType) :: this
    ! -- local
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    ! -- formats
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, &
                              supportOpenClose=.true., blockRequired=.false.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write (iout, '(1x,a)') 'PROCESSING ATS OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case default
          write (errmsg, '(a,a)') 'Unknown ATS option: ', &
            trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (iout, '(1x,a)') 'END OF ATS OPTIONS'
    end if
  end subroutine ats_read_options

  !> @ brief Read dimensions
  !!
  !! Read dimensions from ATS input file.
  !!
  !<
  subroutine ats_read_dimensions(this)
    ! -- dummy
    class(AtsType) :: this
    ! -- local
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    ! -- formats
    character(len=*), parameter :: fmtmaxats = &
      &"(1X,I0,' ADAPTIVE TIME STEP RECORDS(S) WILL FOLLOW IN PERIODDATA')"
    !
    ! -- get DIMENSIONS block
    call this%parser%GetBlock('DIMENSIONS', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse block if detected
    if (isfound) then
      write (iout, '(1x,a)') 'PROCESSING ATS DIMENSIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('MAXATS')
          this%maxats = this%parser%GetInteger()
          write (iout, fmtmaxats) this%maxats
        case default
          write (errmsg, '(a,a)') 'Unknown ATS dimension: ', &
            trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (iout, '(1x,a)') 'END OF ATS DIMENSIONS'
    else
      write (errmsg, '(a)') 'Required DIMENSIONS block not found.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
  end subroutine ats_read_dimensions

  !> @ brief Read timing
  !!
  !! Read timing information from ATS input file.
  !!
  !<
  subroutine ats_read_timing(this)
    ! -- modules
    ! -- dummy
    class(AtsType) :: this
    ! -- local
    integer(I4B) :: ierr
    integer(I4B) :: n
    logical :: isfound, endOfBlock
    ! -- formats
    !
    ! -- get PERIODDATA block
    call this%parser%GetBlock('PERIODDATA', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse block if detected
    if (isfound) then
      write (iout, '(1x,a)') 'READING ATS PERIODDATA'
      do n = 1, this%maxats
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        !
        ! -- fill the ats data arrays
        this%iperats(n) = this%parser%GetInteger()
        this%dt0(n) = this%parser%GetDouble()
        this%dtmin(n) = this%parser%GetDouble()
        this%dtmax(n) = this%parser%GetDouble()
        this%dtadj(n) = this%parser%GetDouble()
        this%dtfailadj(n) = this%parser%GetDouble()
      end do
      !
      ! -- Close the block
      call this%parser%terminateblock()
      !
      ! -- Check for errors
      if (count_errors() > 0) then
        call this%parser%StoreErrorUnit()
      end if
      write (iout, '(1x,a)') 'END READING ATS PERIODDATA'
    else
      write (errmsg, '(a)') 'Required PERIODDATA block not found.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
  end subroutine ats_read_timing

  !> @ brief Process input
  !!
  !! Process ATS input by filling the kperats array.
  !!
  !<
  subroutine ats_process_input(this)
    ! -- dummy
    class(AtsType) :: this

    integer(I4B) :: kkper
    integer(I4B) :: n
    !
    ! -- fill kperats for valid iperats values
    do n = 1, this%maxats
      kkper = this%iperats(n)
      if (kkper > 0 .and. kkper <= this%nper) then
        this%kperats(kkper) = n
      end if
    end do
  end subroutine ats_process_input

  !> @ brief Write input table
  !!
  !! Write a table showing the ATS input read from the perioddata block.
  !!
  !<
  subroutine ats_input_table(this)
    use TableModule, only: TableType, table_cr
    ! -- dummy
    class(AtsType) :: this

    integer(I4B) :: n
    character(len=LINELENGTH) :: tag
    type(TableType), pointer :: inputtab => null()
    !
    ! -- setup table
    call table_cr(inputtab, 'ATS', 'ATS PERIOD DATA')
    call inputtab%table_df(this%maxats, 7, iout)
    !
    ! add columns
    tag = 'RECORD'
    call inputtab%initialize_column(tag, 10, alignment=TABLEFT)
    tag = 'IPERATS'
    call inputtab%initialize_column(tag, 10, alignment=TABLEFT)
    tag = 'DT0'
    call inputtab%initialize_column(tag, 10, alignment=TABCENTER)
    tag = 'DTMIN'
    call inputtab%initialize_column(tag, 10, alignment=TABCENTER)
    tag = 'DTMAX'
    call inputtab%initialize_column(tag, 10, alignment=TABCENTER)
    tag = 'DTADJ'
    call inputtab%initialize_column(tag, 10, alignment=TABCENTER)
    tag = 'DTFAILADJ'
    call inputtab%initialize_column(tag, 10, alignment=TABCENTER)
    !
    ! -- write the data
    do n = 1, this%maxats
      call inputtab%add_term(n)
      call inputtab%add_term(this%iperats(n))
      call inputtab%add_term(this%dt0(n))
      call inputtab%add_term(this%dtmin(n))
      call inputtab%add_term(this%dtmax(n))
      call inputtab%add_term(this%dtadj(n))
      call inputtab%add_term(this%dtfailadj(n))
    end do
    !
    ! -- deallocate the table
    call inputtab%table_da()
    deallocate (inputtab)
    nullify (inputtab)
  end subroutine ats_input_table

  !> @ brief Check timing
  !!
  !! Perform a check on the input data to make sure values are within
  !! required ranges.
  !!
  !<
  subroutine ats_check_timing(this)
    ! -- dummy
    class(AtsType) :: this

    integer(I4B) :: n
    write (iout, '(1x,a)') 'PROCESSING ATS INPUT'
    do n = 1, this%maxats
      !
      ! -- check iperats
      if (this%iperats(n) < 1) then
        write (errmsg, '(a, i0, a, i0)') &
          'IPERATS must be greater than zero.  Found ', this%iperats(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if
      if (this%iperats(n) > this%nper) then
        write (warnmsg, '(a, i0, a, i0)') &
          'IPERATS greater than NPER.  Found ', this%iperats(n), &
          ' for ATS PERIODDATA record ', n
        call store_warning(warnmsg)
      end if
      !
      ! -- check dt0
      if (this%dt0(n) < DZERO) then
        write (errmsg, '(a, g15.7, a, i0)') &
          'DT0 must be >= zero.  Found ', this%dt0(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if
      !
      ! -- check dtmin
      if (this%dtmin(n) <= DZERO) then
        write (errmsg, '(a, g15.7, a, i0)') &
          'DTMIN must be > zero.  Found ', this%dtmin(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if
      !
      ! -- check dtmax
      if (this%dtmax(n) <= DZERO) then
        write (errmsg, '(a, g15.7, a, i0)') &
          'DTMAX must be > zero.  Found ', this%dtmax(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if
      !
      ! -- check dtmin <= dtmax
      if (this%dtmin(n) > this%dtmax(n)) then
        write (errmsg, '(a, 2g15.7, a, i0)') &
          'DTMIN must be < dtmax.  Found ', this%dtmin(n), this%dtmax(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if
      !
      ! -- check dtadj
      if (this%dtadj(n) .ne. DZERO .and. this%dtadj(n) < DONE) then
        write (errmsg, '(a, g15.7, a, i0)') &
          'DTADJ must be 0 or >= 1.0.  Found ', this%dtadj(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if
      !
      ! -- check dtfailadj
      if (this%dtfailadj(n) .ne. DZERO .and. this%dtfailadj(n) < DONE) then
        write (errmsg, '(a, g15.7, a, i0)') &
          'DTFAILADJ must be 0 or >= 1.0.  Found ', this%dtfailadj(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if

    end do
    !
    ! -- Check for errors
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    write (iout, '(1x,a)') 'DONE PROCESSING ATS INPUT'
  end subroutine ats_check_timing

  !> @ brief Write period message
  !!
  !! Write message to mfsim.lst file with information on ATS settings
  !! for this period.
  !!
  !<
  subroutine ats_period_message(this, kper)
    ! -- dummy
    class(AtsType) :: this
    integer(I4B), intent(in) :: kper
    ! -- local
    integer(I4B) :: n
    character(len=*), parameter :: fmtspts = &
    "(28X,'ATS IS OVERRIDING TIME STEPPING FOR THIS PERIOD',/                  &
      &28X,'INITIAL TIME STEP SIZE                 (DT0) = ',G15.7,/           &
      &28X,'MINIMUM TIME STEP SIZE               (DTMIN) = ',G15.7,/           &
      &28X,'MAXIMUM TIME STEP SIZE               (DTMAX) = ',G15.7,/           &
      &28X,'MULTIPLIER/DIVIDER FOR TIME STEP     (DTADJ) = ',G15.7,/           &
      &28X,'DIVIDER FOR FAILED TIME STEP     (DTFAILADJ) = ',G15.7,/           &
      &)"
    n = this%kperats(kper)
    write (iout, fmtspts) this%dt0(n), this%dtmin(n), this%dtmax(n), &
      this%dtadj(n), this%dtfailadj(n)
  end subroutine ats_period_message

  !> @ brief Allow and external caller to submit preferred time step
  !!
  !!  Submit a preferred time step length.  Alternatively, if idir is
  !!  is passed, then either increase or decrease the submitted time
  !!  step by the dtadj input variable.
  !!
  !<
  subroutine ats_submit_delt(this, kstp, kper, dt, sloc, idir)
    ! -- dummy
    class(AtsType) :: this
    integer(I4B), intent(in) :: kstp
    integer(I4B), intent(in) :: kper
    real(DP), intent(in) :: dt
    character(len=*), intent(in) :: sloc
    integer(I4B), intent(in), optional :: idir
    ! -- local
    integer(I4B) :: n
    real(DP) :: tsfact
    real(DP) :: dt_temp
    character(len=*), parameter :: fmtdtsubmit = &
      &"(1x, 'ATS: ', A,' submitted a preferred time step size of ', G15.7)"

    if (this%isAdaptivePeriod(kper)) then
      n = this%kperats(kper)
      tsfact = this%dtadj(n)
      if (tsfact > DONE) then
        !
        ! -- if idir is present, then dt is a length that should be adjusted
        !    (divided by or multiplied by) by dtadj.  If idir is not present
        !    then dt is the submitted time step.
        if (present(idir)) then
          dt_temp = DZERO
          if (idir == -1) then
            dt_temp = dt / tsfact
          else if (idir == 1) then
            dt_temp = dt * tsfact
          end if
        else
          dt_temp = dt
        end if
        if (kstp > 1 .and. dt_temp > DZERO) then
          write (iout, fmtdtsubmit) trim(adjustl(sloc)), dt_temp
        end if
        if (dt_temp > DZERO .and. dt_temp < this%dtstable) then
          ! -- Reset dtstable to a smaller value
          this%dtstable = dt_temp
        end if
      end if
    end if
  end subroutine ats_submit_delt

  !> @ brief Set time step
  !!
  !! Set the time step length (delt) for this time step using the ATS
  !! controls.
  !!
  !<
  subroutine ats_set_delt(this, kstp, kper, pertim, perlencurrent, delt)
    ! -- modules
    ! -- dummy
    class(AtsType) :: this
    integer(I4B), intent(in) :: kstp
    integer(I4B), intent(in) :: kper
    real(DP), intent(inout) :: pertim
    real(DP), intent(in) :: perlencurrent
    real(DP), intent(inout) :: delt
    ! -- local
    integer(I4B) :: n
    real(DP) :: tstart
    ! -- formats
    character(len=*), parameter :: fmtdt = &
      "(1x, 'ATS: time step set to ', G15.7, ' for step ', i0, &
      &' and period ', i0)"
    !
    ! -- initialize the record position (n) for this stress period
    n = this%kperats(kper)
    !
    ! -- set tstart to the end of the last time step.
    tstart = pertim
    !
    ! -- Calculate delt
    !
    ! -- Setup new stress period if kstp is 1
    if (kstp == 1) then
      !
      ! -- Assign first value of delt for this stress period
      if (this%dt0(n) /= DZERO) then
        delt = this%dt0(n)
      else
        ! leave delt the way it was
      end if
    else
      !
      ! -- Assign delt based on stability
      if (this%dtstable /= DNODATA) then
        delt = this%dtstable
        this%dtstable = DNODATA
      end if
    end if
    !
    ! -- Ensure  tsmin < delt < tsmax
    if (delt < this%dtmin(n)) then
      delt = this%dtmin(n)
    end if
    if (delt > this%dtmax(n)) then
      delt = this%dtmax(n)
    end if
    !
    ! -- Cut timestep down to meet end of period
    if (tstart + delt > perlencurrent - this%dtmin(n)) then
      delt = perlencurrent - tstart
    end if
    !
    ! -- Write time step size information
    write (iout, fmtdt) delt, kstp, kper
  end subroutine ats_set_delt

  !> @ brief Reset time step because failure has occurred
  !!
  !!  Reset the time step using dtfailadj because the time step
  !!  did not converge.
  !!
  !<
  subroutine ats_reset_delt(this, kstp, kper, lastStepFailed, &
                            delt, finishedTrying)
    ! -- modules
    ! -- dummy
    class(AtsType) :: this
    integer(I4B), intent(in) :: kstp
    integer(I4B), intent(in) :: kper
    integer(I4B), intent(in) :: lastStepFailed
    real(DP), intent(inout) :: delt
    logical, intent(inout) :: finishedTrying
    ! -- local
    integer(I4B) :: n
    real(DP) :: delt_temp
    real(DP) :: tsfact
    ! -- formats
    character(len=*), parameter :: fmttsi = &
      "(1X, 'Failed solution for step ', i0, ' and period ', i0, &
      &' will be retried using time step of ', G15.7)"
    if (this%isAdaptivePeriod(kper)) then
      if (lastStepFailed /= 0) then
        delt_temp = delt
        n = this%kperats(kper)
        tsfact = this%dtfailadj(n)
        if (tsfact > DONE) then
          delt_temp = delt / tsfact
          if (delt_temp >= this%dtmin(n)) then
            finishedTrying = .false.
            delt = delt_temp
            write (iout, fmttsi) kstp, kper, delt
          end if
        end if

      end if
    end if
  end subroutine ats_reset_delt

  !> @ brief Set end of period indicator
  !!
  !! Determine if it is the end of the stress period and set the endofperiod
  !! logical variable if so.
  !!
  !<
  subroutine ats_set_endofperiod(this, kper, pertim, perlencurrent, endofperiod)
    ! -- dummy
    class(AtsType) :: this
    integer(I4B), intent(in) :: kper
    real(DP), intent(inout) :: pertim
    real(DP), intent(in) :: perlencurrent
    logical(LGP), intent(inout) :: endofperiod
    ! -- local
    integer(I4B) :: n
    !
    ! -- End of stress period and/or simulation?
    n = this%kperats(kper)
    if (abs(pertim - perlencurrent) < this%dtmin(n)) then
      endofperiod = .true.
    end if
  end subroutine ats_set_endofperiod

end module AdaptiveTimeStepModule
