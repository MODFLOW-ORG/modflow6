module TimerModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH, DZERO, DONE, DSECPERDY, DSECPERHR, &
                             DSIXTY
  use MessageModule, only: write_message
  implicit none
  private
  public :: print_start_time
  public :: elapsed_time
  public :: code_timer
  integer(I4B), dimension(8) :: ibdt

contains

  !> @brief Start simulation timer
  !<
  subroutine print_start_time()
    ! -- local
    character(len=LINELENGTH) :: line
    integer(I4B) :: i
    ! -- formats
    character(len=*), parameter :: fmtdt = &
      "(1X,'Run start date and time (yyyy/mm/dd hh:mm:ss): ', &
      &I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2)"
    !
    ! -- Get current date and time, assign to IBDT, and write to screen
    call date_and_time(values=ibdt)
    write (line, fmtdt) (ibdt(i), i=1, 3), (ibdt(i), i=5, 7)
    call write_message(line, skipafter=1)
  end subroutine print_start_time

  !> @brief Get end time and calculate elapsed time
  !<
  subroutine elapsed_time(iout, iprtim)
    ! -- dummy
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(in) :: iprtim
    ! -- local
    character(len=LINELENGTH) :: line
    integer(I4B) :: IEDT(8), IDPM(12)
    integer(I4B) :: NSPD
    integer(I4B) :: i
    integer(I4B) :: ndays, leap, ibd, ied, mb, me, nm, mc, m
    integer(I4B) :: nhours, nmins, nsecs, msecs, nrsecs
    real(DP) :: elsec, rsecs
    data IDPM/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/ !< Days per month
    data NSPD/86400/ !< Seconds per day
    ! -- formats
    character(len=*), parameter :: fmtdt = &
      "(1x,'Run end date and time (yyyy/mm/dd hh:mm:ss): ', &
      &I4,'/',I2.2,'/',I2.2,1x,I2,':',I2.2,':',I2.2)"
    character(len=*), parameter :: fmttma = &
      "(1x,'Elapsed run time: ',I3,' Days, ',I2,' Hours, ',I2, &
      &' Minutes, ', I2, ' Seconds')"
    character(len=*), parameter :: fmttmb = &
      &"(1x,'Elapsed run time: ',I2,' Hours, ',I2,' Minutes, ',I2,' Seconds')"
    character(len=*), parameter :: fmttmc = &
      &"(1x,'Elapsed run time: ',I2,' Minutes, ',I2,'.',I3.3,' Seconds')"
    character(len=*), parameter :: fmttmd = &
      &"(1x,'Elapsed run time: ',I2,'.',I3.3,' Seconds')"
    !
    ! -- Get current date and time, assign to IEDT, and write.
    call date_and_time(values=IEDT)
    !
    ! -- Write elapsed time to stdout
    write (line, fmtdt) (IEDT(I), I=1, 3), (IEDT(I), I=5, 7)
    call write_message(line, skipbefore=1)
    !
    ! -- Write elapsted time to iout
    if (iprtim > 0) then
      call write_message(line, iunit=iout, skipbefore=1)
    end if
    !
    ! -- Calculate elapsed time in days and seconds
    ndays = 0
    leap = 0
    if (mod(IEDT(1), 4) == 0) leap = 1
    ibd = ibdt(3) !< Begin Day
    ied = IEDT(3) !< End Day
    !
    ! -- Find days
    if (ibdt(2) /= IEDT(2)) then
      ! -- Months differ
      mb = ibdt(2) !< Begin month
      me = IEDT(2) !< End month
      nm = me - mb + 1 !< Number of months to look at
      if (mb > me) nm = nm + 12
      mc = mb - 1
      do m = 1, nm
        mc = mc + 1 !< mc is current month
        if (mc == 13) mc = 1
        if (mc == mb) then
          ndays = ndays + IDPM(mc) - ibd
          if (mc == 2) ndays = ndays + leap
        elseif (mc == me) then
          ndays = ndays + ied
        else
          ndays = ndays + IDPM(mc)
          if (mc == 2) ndays = ndays + leap
        end if
      end do
    elseif (ibd < ied) then
      ! -- Start and end in same month, only account for days
      ndays = ied - ibd
    end if
    elsec = ndays * NSPD
    !
    ! -- Add or subtract seconds
    elsec = elsec + (IEDT(5) - ibdt(5)) * 3600.0
    elsec = elsec + (IEDT(6) - ibdt(6)) * 60.0
    elsec = elsec + (IEDT(7) - ibdt(7))
    elsec = elsec + (IEDT(8) - ibdt(8)) * 0.001
    !
    ! -- Convert seconds to days, hours, minutes, and seconds
    ndays = int(elsec / nspd)
    rsecs = mod(elsec, DSECPERDY)
    nhours = int(rsecs / 3600.0)
    rsecs = mod(rsecs, DSECPERHR)
    nmins = int(rsecs / 60.0)
    rsecs = mod(rsecs, DSIXTY)
    nsecs = int(rsecs)
    rsecs = mod(rsecs, DONE)
    msecs = nint(rsecs * 1000.0)
    nrsecs = nsecs
    if (rsecs > 0.5) nrsecs = nrsecs + 1
    !
    ! -- Write elapsed time to screen
    if (ndays > 0) then
      write (line, fmttma) ndays, nhours, nmins, nrsecs
    elseif (nhours > 0) then
      write (line, fmttmb) nhours, nmins, nrsecs
    elseif (nmins > 0) then
      write (line, fmttmc) nmins, nsecs, msecs
    else
      write (line, fmttmd) nsecs, msecs
    end if
    call write_message(line, skipafter=1)
    !
    ! -- Write times to file if requested
    if (iprtim > 0) then
      IF (NDAYS > 0) THEN
        WRITE (IOUT, fmttma) NDAYS, NHOURS, NMINS, NRSECS
      ELSEIF (NHOURS > 0) THEN
        WRITE (IOUT, fmttmb) NHOURS, NMINS, NRSECS
      ELSEIF (NMINS > 0) THEN
        WRITE (IOUT, fmttmc) NMINS, NSECS, MSECS
      ELSE
        WRITE (IOUT, fmttmd) NSECS, MSECS
      END IF
    END IF
  end subroutine elapsed_time

  !> @brief Get end time and calculate elapsed time
  !!
  !! Timer for subroutines
  !<
  subroutine code_timer(it, t1, ts)
    ! -- dummy
    integer(I4B), intent(in) :: it
    real(DP), intent(inout) :: t1
    real(DP), intent(inout) :: ts
    ! -- local
    real(DP) :: dt
    !
    if (it == 0) then
      call CPU_TIME(t1)
    else
      call CPU_TIME(dt)
      ts = ts + dt - t1
    end if
  end subroutine code_timer

end module TimerModule
