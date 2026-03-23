!stress periods and time stepping is handled by these routines
!convert this to a derived type?  May not be necessary since only
!one of them is needed.

module TdisModule

  use KindModule, only: DP, I4B, LGP
  use SimVariablesModule, only: iout, isim_level
  use ConstantsModule, only: LINELENGTH, LENDATETIME, LENMEMPATH, VALL
  !
  implicit none
  !
  private
  public :: tdis_cr
  public :: tdis_set_counters
  public :: tdis_set_timestep
  public :: tdis_delt_reset
  public :: tdis_ot
  public :: tdis_da
  !
  integer(I4B), public, pointer :: nper => null() !< number of stress period
  integer(I4B), public, pointer :: itmuni => null() !< flag indicating time units
  integer(I4B), public, pointer :: kper => null() !< current stress period number
  integer(I4B), public, pointer :: kstp => null() !< current time step number
  integer(I4B), public, pointer :: inats => null() !< flag indicating ats active for simulation
  logical(LGP), public, pointer :: readnewdata => null() !< flag indicating time to read new data
  logical(LGP), public, pointer :: endofperiod => null() !< flag indicating end of stress period
  logical(LGP), public, pointer :: endofsimulation => null() !< flag indicating end of simulation
  real(DP), public, pointer :: delt => null() !< length of the current time step
  real(DP), public, pointer :: pertim => null() !< time relative to start of stress period
  real(DP), public, pointer :: topertim => null() !< simulation time at start of stress period
  real(DP), public, pointer :: totim => null() !< time relative to start of simulation
  real(DP), public, pointer :: totimc => null() !< simulation time at start of time step
  real(DP), public, pointer :: deltsav => null() !< saved value for delt, used for subtiming
  real(DP), public, pointer :: totimsav => null() !< saved value for totim, used for subtiming
  real(DP), public, pointer :: pertimsav => null() !< saved value for pertim, used for subtiming
  real(DP), public, pointer :: totalsimtime => null() !< time at end of simulation
  real(DP), public, dimension(:), pointer, contiguous :: perlen => null() !< length of each stress period
  integer(I4B), public, dimension(:), pointer, contiguous :: nstp => null() !< number of time steps in each stress period
  real(DP), public, dimension(:), pointer, contiguous :: tsmult => null() !< time step multiplier for each stress period
  character(len=LENDATETIME), public, pointer :: datetime0 => null() !< starting date and time for the simulation
  character(len=LENMEMPATH), pointer :: input_mempath => null() !< input context mempath for tdis
  character(len=LINELENGTH), pointer :: input_fname => null() !< input filename for tdis
  !
contains

  !> @brief Create temporal discretization
  !<
  subroutine tdis_cr(fname, inmempath)
    ! -- modules
    use InputOutputModule, only: getunit, openfile
    use ConstantsModule, only: LINELENGTH, DZERO
    use AdaptiveTimeStepModule, only: ats_cr
    ! -- dummy
    character(len=*), intent(in) :: fname
    character(len=*), intent(in) :: inmempath
    ! -- local
    ! -- formats
    character(len=*), parameter :: fmtheader = &
     "(1X,/1X,'TDIS -- TEMPORAL DISCRETIZATION PACKAGE,',   /                  &
      &' VERSION 1 : 11/13/2014 - INPUT READ FROM MEMPATH: ', A)"
    !
    ! -- Allocate the scalar variables
    call tdis_allocate_scalars()
    !
    ! -- set input context and fname
    input_fname = fname
    input_mempath = inmempath
    !
    ! -- Identify package
    write (iout, fmtheader) input_mempath
    !
    ! -- Source options
    call tdis_source_options()
    !
    ! -- Source dimensions and then allocate arrays
    call tdis_source_dimensions()
    call tdis_allocate_arrays()
    !
    ! -- Source timing
    call tdis_source_timing()
    !
    if (inats > 0) then
      call ats_cr(inats, nper)
    end if
  end subroutine tdis_cr

  !> @brief Set kstp and kper
  !<
  subroutine tdis_set_counters()
    ! -- modules
    use ConstantsModule, only: DONE, DZERO, MNORMAL, MVALIDATE, DNODATA
    use SimVariablesModule, only: isim_mode
    use MessageModule, only: write_message
    use AdaptiveTimeStepModule, only: isAdaptivePeriod, dtstable, &
                                      ats_period_message
    ! -- local
    character(len=LINELENGTH) :: line
    character(len=4) :: cpref
    character(len=10) :: cend
    ! -- formats
    character(len=*), parameter :: fmtspts = &
      &"(a, 'Solving:  Stress period: ',i5,4x, 'Time step: ',i5,4x, a)"
    character(len=*), parameter :: fmtvspts = &
      &"(' Validating:  Stress period: ',i5,4x,'Time step: ',i5,4x)"
    character(len=*), parameter :: fmtspi = &
      "('1',/1X,'STRESS PERIOD NO. ',I0,', LENGTH =',G15.7,/ &
      &1X,42('-'))"
    character(len=*), parameter :: fmtspits = &
      "(1X,'NUMBER OF TIME STEPS = ',I0,/ &
      &1X,'MULTIPLIER FOR DELT =',F10.3)"
    !
    ! -- Initialize variables for this step
    if (inats > 0) dtstable = DNODATA
    readnewdata = .false.
    cpref = '    '
    cend = ''
    !
    ! -- Increment kstp and kper
    if (endofperiod) then
      kstp = 1
      kper = kper + 1
      readnewdata = .true.
    else
      kstp = kstp + 1
    end if
    !
    ! -- Print stress period and time step to console
    select case (isim_mode)
    case (MVALIDATE)
      write (line, fmtvspts) kper, kstp
    case (MNORMAL)
      write (line, fmtspts) cpref, kper, kstp, trim(cend)
    end select
    if (isim_level >= VALL) &
      call write_message(line)
    call write_message(line, iunit=iout, skipbefore=1, skipafter=1)
    !
    ! -- Write message if first time step
    if (kstp == 1) then
      write (iout, fmtspi) kper, perlen(kper)
      if (isAdaptivePeriod(kper)) then
        call ats_period_message(kper)
      else
        write (iout, fmtspits) nstp(kper), tsmult(kper)
      end if
    end if
  end subroutine tdis_set_counters

  !> @brief Set time step length
  !<
  subroutine tdis_set_timestep()
    ! -- modules
    use ConstantsModule, only: DONE, DZERO
    use AdaptiveTimeStepModule, only: isAdaptivePeriod, &
                                      ats_set_delt, &
                                      ats_set_endofperiod
    ! -- local
    logical(LGP) :: adaptivePeriod
    ! -- format
    character(len=*), parameter :: fmttsi = &
                                   "(1X,'INITIAL TIME STEP SIZE =',G15.7)"
    !
    ! -- Initialize
    adaptivePeriod = isAdaptivePeriod(kper)
    if (kstp == 1) then
      pertim = DZERO
      topertim = DZERO
    end if
    !
    ! -- Set delt
    if (adaptivePeriod) then
      call ats_set_delt(kstp, kper, pertim, perlen(kper), delt)
    else
      call tdis_set_delt()
      if (kstp == 1) then
        write (iout, fmttsi) delt
      end if
    end if
    !
    ! -- Advance timers and update totim and pertim based on delt
    totimsav = totim
    pertimsav = pertim
    totimc = totimsav
    totim = totimsav + delt
    pertim = pertimsav + delt
    !
    ! -- Set end of period indicator
    endofperiod = .false.
    if (adaptivePeriod) then
      call ats_set_endofperiod(kper, pertim, perlen(kper), endofperiod)
    else
      if (kstp == nstp(kper)) then
        endofperiod = .true.
      end if
    end if
    if (endofperiod) then
      pertim = perlen(kper)
    end if
    !
    ! -- Set end of simulation indicator
    if (endofperiod .and. kper == nper) then
      endofsimulation = .true.
    end if
  end subroutine tdis_set_timestep

  !> @brief Reset delt and update timing variables and indicators
  !!
  !! This routine is called when a timestep fails to converge, and so it is
  !! retried using a smaller time step (deltnew).
  !<
  subroutine tdis_delt_reset(deltnew)
    ! -- modules
    use ConstantsModule, only: DONE, DZERO
    use AdaptiveTimeStepModule, only: isAdaptivePeriod, &
                                      ats_set_delt, &
                                      ats_set_endofperiod
    ! -- dummy
    real(DP), intent(in) :: deltnew
    ! -- local
    logical(LGP) :: adaptivePeriod
    !
    ! -- Set values
    adaptivePeriod = isAdaptivePeriod(kper)
    delt = deltnew
    totim = totimsav + delt
    pertim = pertimsav + delt
    !
    ! -- Set end of period indicator
    endofperiod = .false.
    if (adaptivePeriod) then
      call ats_set_endofperiod(kper, pertim, perlen(kper), endofperiod)
    else
      if (kstp == nstp(kper)) then
        endofperiod = .true.
      end if
    end if
    !
    ! -- Set end of simulation indicator
    if (endofperiod .and. kper == nper) then
      endofsimulation = .true.
      totim = totalsimtime
    end if
  end subroutine tdis_delt_reset

  !> @brief Set time step length
  !<
  subroutine tdis_set_delt()
    ! -- modules
    use ConstantsModule, only: DONE
    !
    if (kstp == 1) then
      ! -- Calculate the first value of delt for this stress period
      topertim = totim
      if (tsmult(kper) /= DONE) then
        ! -- Timestep length has a geometric progression
        delt = perlen(kper) * (DONE - tsmult(kper)) / &
               (DONE - tsmult(kper)**nstp(kper))
      else
        ! -- Timestep length is constant
        delt = perlen(kper) / float(nstp(kper))
      end if
    elseif (kstp == nstp(kper)) then
      ! -- Calculate exact last delt to avoid accumulation errors
      delt = topertim + perlen(kper) - totim
    else
      delt = tsmult(kper) * delt
    end if
  end subroutine tdis_set_delt

  !> @brief Print simulation time
  !<
  subroutine tdis_ot(iout)
    ! -- modules
    use ConstantsModule, only: DZERO, DONE, DSIXTY, DSECPERHR, DHRPERDAY, &
                               DDYPERYR, DSECPERDY, DSECPERYR
    ! -- dummy
    integer(I4B), intent(in) :: iout
    ! -- local
    real(DP) :: cnv, delsec, totsec, persec, delmn, delhr, totmn, tothr, &
                totdy, totyr, permn, perhr, perdy, peryr, deldy, delyr
    ! -- format
    character(len=*), parameter :: fmttmsmry = "(1X, ///9X, &
      &'TIME SUMMARY AT END OF TIME STEP', I5,' IN STRESS PERIOD ', I4)"
    character(len=*), parameter :: fmttmstpmsg = &
      &"(21X, '     TIME STEP LENGTH =', G15.6 / &
      & 21X, '   STRESS PERIOD TIME =', G15.6 / &
      & 21X, 'TOTAL SIMULATION TIME =', G15.6)"
    character(len=*), parameter :: fmttottmmsg = &
      &"(19X, ' SECONDS     MINUTES      HOURS', 7X, &
      &'DAYS        YEARS'/20X, 59('-'))"
    character(len=*), parameter :: fmtdelttm = &
      &"(1X, '  TIME STEP LENGTH', 1P, 5G12.5)"
    character(len=*), parameter :: fmtpertm = &
      &"(1X, 'STRESS PERIOD TIME', 1P, 5G12.5)"
    character(len=*), parameter :: fmttottm = &
      &"(1X, '        TOTAL TIME', 1P, 5G12.5,/)"
    !
    ! -- Write header message for the information that follows
    write (iout, fmttmsmry) kstp, kper
    !
    ! -- Use time unit indicator to get factor to convert to seconds
    cnv = DZERO
    if (itmuni == 1) cnv = DONE
    if (itmuni == 2) cnv = DSIXTY
    if (itmuni == 3) cnv = DSECPERHR
    if (itmuni == 4) cnv = DSECPERDY
    if (itmuni == 5) cnv = DSECPERYR
    !
    ! -- If FACTOR=0 then time units are non-standard
    if (cnv == DZERO) then
      ! -- Print times in non-standard time units
      write (iout, fmttmstpmsg) delt, pertim, totim
    else
      ! -- Calculate length of time step & elapsed time in seconds
      delsec = cnv * delt
      totsec = cnv * totim
      persec = cnv * pertim
      !
      ! -- Calculate times in minutes, hours, days, and years
      delmn = delsec / DSIXTY
      delhr = delmn / DSIXTY
      deldy = delhr / DHRPERDAY
      delyr = deldy / DDYPERYR
      totmn = totsec / DSIXTY
      tothr = totmn / DSIXTY
      totdy = tothr / DHRPERDAY
      totyr = totdy / DDYPERYR
      permn = persec / DSIXTY
      perhr = permn / DSIXTY
      perdy = perhr / DHRPERDAY
      peryr = perdy / DDYPERYR
      !
      ! -- Print time step length and elapsed times in all time units
      write (iout, fmttottmmsg)
      write (iout, fmtdelttm) delsec, delmn, delhr, deldy, delyr
      write (iout, fmtpertm) persec, permn, perhr, perdy, peryr
      write (iout, fmttottm) totsec, totmn, tothr, totdy, totyr
    end if
  end subroutine tdis_ot

  !> @brief Deallocate memory
  !<
  subroutine tdis_da()
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    use AdaptiveTimeStepModule, only: ats_da
    !
    ! -- ats
    if (inats > 0) call ats_da()
    !
    ! -- Scalars
    call mem_deallocate(nper)
    call mem_deallocate(itmuni)
    call mem_deallocate(kper)
    call mem_deallocate(kstp)
    call mem_deallocate(inats)
    call mem_deallocate(readnewdata)
    call mem_deallocate(endofperiod)
    call mem_deallocate(endofsimulation)
    call mem_deallocate(delt)
    call mem_deallocate(pertim)
    call mem_deallocate(topertim)
    call mem_deallocate(totim)
    call mem_deallocate(totimc)
    call mem_deallocate(deltsav)
    call mem_deallocate(totimsav)
    call mem_deallocate(pertimsav)
    call mem_deallocate(totalsimtime)
    !
    ! -- strings
    deallocate (datetime0)
    deallocate (input_mempath)
    deallocate (input_fname)
    !
    ! -- Arrays
    call mem_deallocate(perlen)
    call mem_deallocate(nstp)
    call mem_deallocate(tsmult)
  end subroutine tdis_da

  !> @brief Source the timing discretization options
  !<
  subroutine tdis_source_options()
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use InputOutputModule, only: GetUnit, openfile
    use SimModule, only: store_error, store_error_filename
    use MemoryManagerExtModule, only: mem_set_value
    use SourceCommonModule, only: filein_fname
    use SimTdisInputModule, only: SimTdisParamFoundType
    ! -- local
    type(SimTdisParamFoundType) :: found
    character(len=LINELENGTH), dimension(6) :: time_units = &
      &[character(len=LINELENGTH) :: 'UNDEFINED', 'SECONDS', 'MINUTES', 'HOURS', &
                                     'DAYS', 'YEARS']
    character(len=LINELENGTH) :: fname
    ! -- formats
    character(len=*), parameter :: fmtitmuni = &
      &"(4x,'SIMULATION TIME UNIT IS ',A)"
    character(len=*), parameter :: fmtdatetime0 = &
      &"(4x,'SIMULATION STARTING DATE AND TIME IS ',A)"
    !
    ! -- initialize time unit to undefined
    itmuni = 0
    !
    ! -- source options from input context
    call mem_set_value(itmuni, 'TIME_UNITS', input_mempath, time_units, &
                       found%time_units)
    call mem_set_value(datetime0, 'START_DATE_TIME', input_mempath, &
                       found%start_date_time)
    !
    if (found%time_units) then
      if (itmuni == 0) then
        call store_error('Unrecognized input value for TIME_UNITS option.')
        call store_error_filename(input_fname)
      else
        !
        ! -- adjust to 0-based indexing for itmuni
        itmuni = itmuni - 1
      end if
    end if
    !
    ! -- enforce 0 or 1 ATS6_FILENAME entries in option block
    if (filein_fname(fname, 'ATS6_FILENAME', input_mempath, &
                     input_fname)) then
      inats = GetUnit()
      call openfile(inats, iout, fname, 'ATS')
    end if
    !
    ! -- log values to list file
    write (iout, '(1x,a)') 'PROCESSING TDIS OPTIONS'
    !
    if (found%time_units) then
      select case (itmuni)
      case (0)
        write (iout, fmtitmuni) 'UNDEFINED'
      case (1)
        write (iout, fmtitmuni) 'SECONDS'
      case (2)
        write (iout, fmtitmuni) 'MINUTES'
      case (3)
        write (iout, fmtitmuni) 'HOURS'
      case (4)
        write (iout, fmtitmuni) 'DAYS'
      case (5)
        write (iout, fmtitmuni) 'YEARS'
      case default
      end select
    else
      write (iout, fmtitmuni) 'UNDEFINED'
    end if
    !
    if (found%start_date_time) then
      write (iout, fmtdatetime0) datetime0
    end if
    !
    write (iout, '(1x,a)') 'END OF TDIS OPTIONS'
  end subroutine tdis_source_options

  !> @brief Allocate tdis scalars
  !<
  subroutine tdis_allocate_scalars()
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use ConstantsModule, only: DZERO
    !
    ! -- memory manager variables
    call mem_allocate(nper, 'NPER', 'TDIS')
    call mem_allocate(itmuni, 'ITMUNI', 'TDIS')
    call mem_allocate(kper, 'KPER', 'TDIS')
    call mem_allocate(kstp, 'KSTP', 'TDIS')
    call mem_allocate(inats, 'INATS', 'TDIS')
    call mem_allocate(readnewdata, 'READNEWDATA', 'TDIS')
    call mem_allocate(endofperiod, 'ENDOFPERIOD', 'TDIS')
    call mem_allocate(endofsimulation, 'ENDOFSIMULATION', 'TDIS')
    call mem_allocate(delt, 'DELT', 'TDIS')
    call mem_allocate(pertim, 'PERTIM', 'TDIS')
    call mem_allocate(topertim, 'TOPERTIM', 'TDIS')
    call mem_allocate(totim, 'TOTIM', 'TDIS')
    call mem_allocate(totimc, 'TOTIMC', 'TDIS')
    call mem_allocate(deltsav, 'DELTSAV', 'TDIS')
    call mem_allocate(totimsav, 'TOTIMSAV', 'TDIS')
    call mem_allocate(pertimsav, 'PERTIMSAV', 'TDIS')
    call mem_allocate(totalsimtime, 'TOTALSIMTIME', 'TDIS')
    !
    ! -- strings
    allocate (datetime0)
    allocate (input_mempath)
    allocate (input_fname)
    !
    ! -- Initialize variables
    nper = 0
    itmuni = 0
    kper = 0
    kstp = 0
    inats = 0
    readnewdata = .true.
    endofperiod = .true.
    endofsimulation = .false.
    delt = DZERO
    pertim = DZERO
    topertim = DZERO
    totim = DZERO
    totimc = DZERO
    deltsav = DZERO
    totimsav = DZERO
    pertimsav = DZERO
    totalsimtime = DZERO
    datetime0 = ''
  end subroutine tdis_allocate_scalars

  !> @brief Allocate tdis arrays
  !<
  subroutine tdis_allocate_arrays()
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    !
    call mem_allocate(perlen, nper, 'PERLEN', 'TDIS')
    call mem_allocate(nstp, nper, 'NSTP', 'TDIS')
    call mem_allocate(tsmult, nper, 'TSMULT', 'TDIS')
  end subroutine tdis_allocate_arrays

  !> @brief Source dimension NPER
  !<
  subroutine tdis_source_dimensions()
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use MemoryManagerExtModule, only: mem_set_value
    use SourceCommonModule, only: filein_fname
    use SimTdisInputModule, only: SimTdisParamFoundType
    ! -- local
    type(SimTdisParamFoundType) :: found
    ! -- formats
    character(len=*), parameter :: fmtnper = &
                                   "(1X,I4,' STRESS PERIOD(S) IN SIMULATION')"
    !
    ! -- source dimensions from input context
    call mem_set_value(nper, 'NPER', input_mempath, found%nper)
    !
    ! -- log values to list file
    write (iout, '(1x,a)') 'PROCESSING TDIS DIMENSIONS'
    !
    if (found%nper) then
      write (iout, fmtnper) nper
    end if
    !
    write (iout, '(1x,a)') 'END OF TDIS DIMENSIONS'
  end subroutine tdis_source_dimensions

  !> @brief Source timing information
  !<
  subroutine tdis_source_timing()
    ! -- modules
    use ConstantsModule, only: LINELENGTH, DZERO
    use SimModule, only: store_error_filename, count_errors
    use MemoryManagerExtModule, only: mem_set_value
    use SourceCommonModule, only: filein_fname
    use SimTdisInputModule, only: SimTdisParamFoundType
    ! -- local
    type(SimTdisParamFoundType) :: found
    integer(I4B) :: n
    ! -- formats
    character(len=*), parameter :: fmtheader = &
      "(1X,//1X,'STRESS PERIOD     LENGTH       TIME STEPS', &
       &'     MULTIPLIER FOR DELT',/1X,76('-'))"
    character(len=*), parameter :: fmtrow = &
                                   "(1X,I8,1PG21.7,I7,0PF25.3)"
    !
    ! -- source perioddata from input context
    call mem_set_value(perlen, 'PERLEN', input_mempath, found%perlen)
    call mem_set_value(nstp, 'NSTP', input_mempath, found%nstp)
    call mem_set_value(tsmult, 'TSMULT', input_mempath, found%tsmult)
    !
    ! -- Check timing information
    call check_tdis_timing(nper, perlen, nstp, tsmult)
    !
    ! -- Check for errors
    if (count_errors() > 0) then
      call store_error_filename(input_fname)
    end if
    !
    ! -- log timing
    write (iout, '(1x,a)') 'PROCESSING TDIS PERIODDATA'
    write (iout, fmtheader)
    !
    do n = 1, size(perlen)
      write (iout, fmtrow) n, perlen(n), nstp(n), tsmult(n)
      totalsimtime = totalsimtime + perlen(n)
    end do
    !
    write (iout, '(1x,a)') 'END OF TDIS PERIODDATA'
  end subroutine tdis_source_timing

  !> @brief Check the tdis timing information
  !!
  !! Return back to tdis_read_timing if an error condition is found and let the
  !! ustop routine be called there instead so the StoreErrorUnit routine can be
  !! called to assign the correct file name.
  !<
  subroutine check_tdis_timing(nper, perlen, nstp, tsmult)
    ! -- modules
    use ConstantsModule, only: LINELENGTH, DZERO, DONE
    use SimModule, only: store_error
    ! -- dummy
    integer(I4B), intent(in) :: nper
    real(DP), dimension(:), contiguous, intent(in) :: perlen
    integer(I4B), dimension(:), contiguous, intent(in) :: nstp
    real(DP), dimension(:), contiguous, intent(in) :: tsmult
    ! -- local
    integer(I4B) :: kper, kstp
    real(DP) :: tstart, tend, dt
    character(len=LINELENGTH) :: errmsg
    ! -- formats
    character(len=*), parameter :: fmtpwarn = &
      "(1X,/1X,'PERLEN is zero for stress period ', I0, &
      &'. PERLEN must not be zero for transient periods.')"
    character(len=*), parameter :: fmtsperror = &
      &"(A,' for stress period ', I0)"
    character(len=*), parameter :: fmtdterror = &
      "('Time step length of ', G0, ' is too small in period ', I0, &
      &' and time step ', I0)"
    !
    ! -- Initialize
    tstart = DZERO
    !
    ! -- Go through and check each stress period
    do kper = 1, nper
      !
      ! -- Error if nstp less than or equal to zero
      if (nstp(kper) <= 0) then
        write (errmsg, fmtsperror) 'Number of time steps less than one ', kper
        call store_error(errmsg)
        return
      end if
      !
      ! -- Warn if perlen is zero
      if (perlen(kper) == DZERO) then
        write (iout, fmtpwarn) kper
        return
      end if
      !
      ! -- Error if tsmult is less than zero
      if (tsmult(kper) <= DZERO) then
        write (errmsg, fmtsperror) 'TSMULT must be greater than 0.0 ', kper
        call store_error(errmsg)
        return
      end if
      !
      ! -- Error if negative period length
      if (perlen(kper) < DZERO) then
        write (errmsg, fmtsperror) 'PERLEN cannot be less than 0.0 ', kper
        call store_error(errmsg)
        return
      end if
      !
      ! -- Go through all time step lengths and make sure they are valid
      do kstp = 1, nstp(kper)
        if (kstp == 1) then
          dt = perlen(kper) / float(nstp(kper))
          if (tsmult(kper) /= DONE) &
            dt = perlen(kper) * (DONE - tsmult(kper)) / &
                 (DONE - tsmult(kper)**nstp(kper))
        else
          dt = dt * tsmult(kper)
        end if
        tend = tstart + dt
        !
        ! -- Error condition if tstart == tend
        if (tstart == tend) then
          write (errmsg, fmtdterror) dt, kper, kstp
          call store_error(errmsg)
          return
        end if
      end do
      !
      ! -- reset tstart = tend
      tstart = tend
      !
    end do
  end subroutine check_tdis_timing

end module TdisModule

