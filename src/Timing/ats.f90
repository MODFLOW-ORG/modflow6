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
  public :: isAdaptivePeriod
  public :: ats_set_delt
  public :: ats_period_message
  public :: ats_set_endofperiod
  public :: ats_submit_delt
  public :: ats_reset_delt
  public :: ats_cr
  public :: ats_da

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

  !> @ brief Determine if period is adaptive
  !!
  !!  Check settings and determine if kper is an adaptive
  !!  stress period.
  !!
  !<
  function isAdaptivePeriod(kper) result(lv)
    integer(I4B), intent(in) :: kper
    logical(LGP) :: lv
    lv = .false.
    if (associated(kperats)) then
      if (kperats(kper) > 0) then
        lv = .true.
      end if
    end if
  end function isAdaptivePeriod

  !> @ brief Create ATS object
  !!
  !!  Create a new ATS object, and read and check input.
  !!
  !<
  subroutine ats_cr(inunit, nper_tdis)
    ! -- modules
    ! -- dummy
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: nper_tdis
    ! -- local
    ! -- formats
    character(len=*), parameter :: fmtheader = &
        "(1X,/1X,'ATS -- ADAPTIVE TIME STEP PACKAGE,',   /                  &
        &' VERSION 1 : 03/18/2021 - INPUT READ FROM UNIT ',I0)"
    !
    ! -- Allocate the scalar variables
    call ats_allocate_scalars()
    !
    ! -- Identify package
    write (iout, fmtheader) inunit
    !
    ! -- Initialize block parser
    call parser%initialize(inunit, iout)
    !
    ! -- Read options
    call ats_read_options()
    !
    ! -- store tdis nper in nper
    nper = nper_tdis
    !
    ! -- Read dimensions and then allocate arrays
    call ats_read_dimensions()
    call ats_allocate_arrays()
    !
    ! -- Read timing
    call ats_read_timing()
    !
    ! -- Echo input data to table
    call ats_input_table()
    !
    ! -- Check timing
    call ats_check_timing()
    !
    ! -- Process input
    call ats_process_input()
    !
    ! -- Close the file
    call parser%Clear()
  end subroutine ats_cr

  !> @ brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the ATS package.
  !!
  !<
  subroutine ats_allocate_scalars()
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    !
    ! -- memory manager variables
    call mem_allocate(nper, 'NPER', 'ATS')
    call mem_allocate(maxats, 'MAXATS', 'ATS')
    call mem_allocate(dtstable, 'DTSTABLE', 'ATS')
    !
    ! -- Initialize variables
    nper = 0
    maxats = 0
    dtstable = DNODATA
  end subroutine ats_allocate_scalars

  !> @ brief Allocate arrays
  !!
  !! Allocate and initialize arrays for the ATS package.
  !!
  !<
  subroutine ats_allocate_arrays()
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- local
    integer(I4B) :: n
    !
    call mem_allocate(kperats, nper, 'KPERATS', 'ATS')
    call mem_allocate(iperats, maxats, 'IPERATS', 'ATS')
    call mem_allocate(dt0, maxats, 'DT0', 'ATS')
    call mem_allocate(dtmin, maxats, 'DTMIN', 'ATS')
    call mem_allocate(dtmax, maxats, 'DTMAX', 'ATS')
    call mem_allocate(dtadj, maxats, 'DTADJ', 'ATS')
    call mem_allocate(dtfailadj, maxats, 'DTFAILADJ', 'ATS')
    !
    ! -- initialize kperats
    do n = 1, nper
      kperats(n) = 0
    end do
    !
    ! -- initialize
    do n = 1, maxats
      iperats(n) = 0
      dt0(n) = DZERO
      dtmin(n) = DZERO
      dtmax(n) = DZERO
      dtadj(n) = DZERO
      dtfailadj(n) = DZERO
    end do
  end subroutine ats_allocate_arrays

  !> @ brief Deallocate variables
  !!
  !! Deallocate all ATS variables.
  !!
  !<
  subroutine ats_da()
    use MemoryManagerModule, only: mem_deallocate
    !
    ! -- Scalars
    call mem_deallocate(nper)
    call mem_deallocate(maxats)
    call mem_deallocate(dtstable)
    !
    ! -- Arrays
    call mem_deallocate(kperats)
    call mem_deallocate(iperats)
    call mem_deallocate(dt0)
    call mem_deallocate(dtmin)
    call mem_deallocate(dtmax)
    call mem_deallocate(dtadj)
    call mem_deallocate(dtfailadj)
  end subroutine ats_da

  !> @ brief Read options
  !!
  !! Read options from ATS input file.
  !!
  !<
  subroutine ats_read_options()
    ! -- dummy
    ! -- local
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    ! -- formats
    !
    ! -- get options block
    call parser%GetBlock('OPTIONS', isfound, ierr, &
                         supportOpenClose=.true., blockRequired=.false.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write (iout, '(1x,a)') 'PROCESSING ATS OPTIONS'
      do
        call parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call parser%GetStringCaps(keyword)
        select case (keyword)
        case default
          write (errmsg, '(a,a)') 'Unknown ATS option: ', &
            trim(keyword)
          call store_error(errmsg)
          call parser%StoreErrorUnit()
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
  subroutine ats_read_dimensions()
    ! -- dummy
    ! -- local
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    ! -- formats
    character(len=*), parameter :: fmtmaxats = &
      &"(1X,I0,' ADAPTIVE TIME STEP RECORDS(S) WILL FOLLOW IN PERIODDATA')"
    !
    ! -- get DIMENSIONS block
    call parser%GetBlock('DIMENSIONS', isfound, ierr, &
                         supportOpenClose=.true.)
    !
    ! -- parse block if detected
    if (isfound) then
      write (iout, '(1x,a)') 'PROCESSING ATS DIMENSIONS'
      do
        call parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call parser%GetStringCaps(keyword)
        select case (keyword)
        case ('MAXATS')
          maxats = parser%GetInteger()
          write (iout, fmtmaxats) maxats
        case default
          write (errmsg, '(a,a)') 'Unknown ATS dimension: ', &
            trim(keyword)
          call store_error(errmsg)
          call parser%StoreErrorUnit()
        end select
      end do
      write (iout, '(1x,a)') 'END OF ATS DIMENSIONS'
    else
      write (errmsg, '(a)') 'Required DIMENSIONS block not found.'
      call store_error(errmsg)
      call parser%StoreErrorUnit()
    end if
  end subroutine ats_read_dimensions

  !> @ brief Read timing
  !!
  !! Read timing information from ATS input file.
  !!
  !<
  subroutine ats_read_timing()
    ! -- modules
    ! -- dummy
    ! -- local
    integer(I4B) :: ierr
    integer(I4B) :: n
    logical :: isfound, endOfBlock
    ! -- formats
    !
    ! -- get PERIODDATA block
    call parser%GetBlock('PERIODDATA', isfound, ierr, &
                         supportOpenClose=.true.)
    !
    ! -- parse block if detected
    if (isfound) then
      write (iout, '(1x,a)') 'READING ATS PERIODDATA'
      do n = 1, maxats
        call parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        !
        ! -- fill the ats data arrays
        iperats(n) = parser%GetInteger()
        dt0(n) = parser%GetDouble()
        dtmin(n) = parser%GetDouble()
        dtmax(n) = parser%GetDouble()
        dtadj(n) = parser%GetDouble()
        dtfailadj(n) = parser%GetDouble()
      end do
      !
      ! -- Close the block
      call parser%terminateblock()
      !
      ! -- Check for errors
      if (count_errors() > 0) then
        call parser%StoreErrorUnit()
      end if
      write (iout, '(1x,a)') 'END READING ATS PERIODDATA'
    else
      write (errmsg, '(a)') 'Required PERIODDATA block not found.'
      call store_error(errmsg)
      call parser%StoreErrorUnit()
    end if
  end subroutine ats_read_timing

  !> @ brief Process input
  !!
  !! Process ATS input by filling the kperats array.
  !!
  !<
  subroutine ats_process_input()
    integer(I4B) :: kkper
    integer(I4B) :: n
    !
    ! -- fill kperats for valid iperats values
    do n = 1, maxats
      kkper = iperats(n)
      if (kkper > 0 .and. kkper <= nper) then
        kperats(kkper) = n
      end if
    end do
  end subroutine ats_process_input

  !> @ brief Write input table
  !!
  !! Write a table showing the ATS input read from the perioddata block.
  !!
  !<
  subroutine ats_input_table()
    use TableModule, only: TableType, table_cr
    integer(I4B) :: n
    character(len=LINELENGTH) :: tag
    type(TableType), pointer :: inputtab => null()
    !
    ! -- setup table
    call table_cr(inputtab, 'ATS', 'ATS PERIOD DATA')
    call inputtab%table_df(maxats, 7, iout)
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
    do n = 1, maxats
      call inputtab%add_term(n)
      call inputtab%add_term(iperats(n))
      call inputtab%add_term(dt0(n))
      call inputtab%add_term(dtmin(n))
      call inputtab%add_term(dtmax(n))
      call inputtab%add_term(dtadj(n))
      call inputtab%add_term(dtfailadj(n))
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
  subroutine ats_check_timing()
    integer(I4B) :: n
    write (iout, '(1x,a)') 'PROCESSING ATS INPUT'
    do n = 1, maxats
      !
      ! -- check iperats
      if (iperats(n) < 1) then
        write (errmsg, '(a, i0, a, i0)') &
          'IPERATS must be greater than zero.  Found ', iperats(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if
      if (iperats(n) > nper) then
        write (warnmsg, '(a, i0, a, i0)') &
          'IPERATS greater than NPER.  Found ', iperats(n), &
          ' for ATS PERIODDATA record ', n
        call store_warning(warnmsg)
      end if
      !
      ! -- check dt0
      if (dt0(n) < DZERO) then
        write (errmsg, '(a, g15.7, a, i0)') &
          'DT0 must be >= zero.  Found ', dt0(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if
      !
      ! -- check dtmin
      if (dtmin(n) <= DZERO) then
        write (errmsg, '(a, g15.7, a, i0)') &
          'DTMIN must be > zero.  Found ', dtmin(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if
      !
      ! -- check dtmax
      if (dtmax(n) <= DZERO) then
        write (errmsg, '(a, g15.7, a, i0)') &
          'DTMAX must be > zero.  Found ', dtmax(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if
      !
      ! -- check dtmin <= dtmax
      if (dtmin(n) > dtmax(n)) then
        write (errmsg, '(a, 2g15.7, a, i0)') &
          'DTMIN must be < dtmax.  Found ', dtmin(n), dtmax(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if
      !
      ! -- check dtadj
      if (dtadj(n) .ne. DZERO .and. dtadj(n) < DONE) then
        write (errmsg, '(a, g15.7, a, i0)') &
          'DTADJ must be 0 or >= 1.0.  Found ', dtadj(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if
      !
      ! -- check dtfailadj
      if (dtfailadj(n) .ne. DZERO .and. dtfailadj(n) < DONE) then
        write (errmsg, '(a, g15.7, a, i0)') &
          'DTFAILADJ must be 0 or >= 1.0.  Found ', dtfailadj(n), &
          ' for ATS PERIODDATA record ', n
        call store_error(errmsg)
      end if

    end do
    !
    ! -- Check for errors
    if (count_errors() > 0) then
      call parser%StoreErrorUnit()
    end if
    write (iout, '(1x,a)') 'DONE PROCESSING ATS INPUT'
  end subroutine ats_check_timing

  !> @ brief Write period message
  !!
  !! Write message to mfsim.lst file with information on ATS settings
  !! for this period.
  !!
  !<
  subroutine ats_period_message(kper)
    ! -- dummy
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
    n = kperats(kper)
    write (iout, fmtspts) dt0(n), dtmin(n), dtmax(n), dtadj(n), dtfailadj(n)
  end subroutine ats_period_message

  !> @ brief Allow and external caller to submit preferred time step
  !!
  !!  Submit a preferred time step length.  Alternatively, if idir is
  !!  is passed, then either increase or decrease the submitted time
  !!  step by the dtadj input variable.
  !!
  !<
  subroutine ats_submit_delt(kstp, kper, dt, sloc, idir)
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

    if (isAdaptivePeriod(kper)) then
      n = kperats(kper)
      tsfact = dtadj(n)
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
        if (dt_temp > DZERO .and. dt_temp < dtstable) then
          ! -- Reset dtstable to a smaller value
          dtstable = dt_temp
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
  subroutine ats_set_delt(kstp, kper, pertim, perlencurrent, delt)
    ! -- modules
    ! -- dummy
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
    n = kperats(kper)
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
      if (dt0(n) /= DZERO) then
        delt = dt0(n)
      else
        ! leave delt the way it was
      end if
    else
      !
      ! -- Assign delt based on stability
      if (dtstable /= DNODATA) then
        delt = dtstable
        dtstable = DNODATA
      end if
    end if
    !
    ! -- Ensure  tsmin < delt < tsmax
    if (delt < dtmin(n)) then
      delt = dtmin(n)
    end if
    if (delt > dtmax(n)) then
      delt = dtmax(n)
    end if
    !
    ! -- Cut timestep down to meet end of period
    if (tstart + delt > perlencurrent - dtmin(n)) then
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
  subroutine ats_reset_delt(kstp, kper, lastStepFailed, delt, finishedTrying)
    ! -- modules
    ! -- dummy
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
    if (isAdaptivePeriod(kper)) then
      if (lastStepFailed /= 0) then
        delt_temp = delt
        n = kperats(kper)
        tsfact = dtfailadj(n)
        if (tsfact > DONE) then
          delt_temp = delt / tsfact
          if (delt_temp >= dtmin(n)) then
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
  subroutine ats_set_endofperiod(kper, pertim, perlencurrent, endofperiod)
    integer(I4B), intent(in) :: kper
    real(DP), intent(inout) :: pertim
    real(DP), intent(in) :: perlencurrent
    logical(LGP), intent(inout) :: endofperiod
    ! -- local
    integer(I4B) :: n
    !
    ! -- End of stress period and/or simulation?
    n = kperats(kper)
    if (abs(pertim - perlencurrent) < dtmin(n)) then
      endofperiod = .true.
    end if
  end subroutine ats_set_endofperiod

end module AdaptiveTimeStepModule
