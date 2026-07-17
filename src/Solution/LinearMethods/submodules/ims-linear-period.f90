!> @brief Implementation of stress-period-varying IMS linear settings
!!
!! Submodule of ImsLinearPeriodModule holding the period-data reader. PERIOD
!! blocks are read on demand as the simulation advances, and their keyword
!! overrides are applied to the running linear settings in place, so a setting
!! persists to later periods until a subsequent PERIOD block changes it.
!<
submodule(ImsLinearPeriodModule) ImsLinearPeriodImpl
  use ConstantsModule, only: LINELENGTH
  use SimModule, only: store_error
  use TdisModule, only: nper
contains

  !> @brief Initialize the period-settings utility from an open file unit
  !<
  module procedure imslinperiod_init
  this%inunit = inunit
  this%iout = iout
  this%active = inunit > 0
  this%ionper = 0
  this%lastonper = 0
  if (this%active) then
    call this%parser%Initialize(this%inunit, this%iout)
  end if
  end procedure imslinperiod_init

  !> @brief Apply this stress period's linear-setting overrides (sticky)
  !<
  module procedure imslinperiod_read_period
  ! -- local variables
  integer(I4B) :: ierr
  logical(LGP) :: isfound, end_of_block
  character(len=LINELENGTH) :: keyword, line, errmsg
  !
  changed = .false.
  if (.not. this%active) return
  !
  ! -- advance to the PERIOD block for the current (or a later) period
  if (this%ionper < kper) then
    call this%parser%GetBlock('PERIOD', isfound, ierr, &
                              supportOpenClose=.true., blockRequired=.false.)
    if (isfound) then
      this%lastonper = this%ionper
      this%ionper = this%parser%GetInteger()
      if (this%ionper <= this%lastonper) then
        write (errmsg, '(a,i0,a,i0,a,i0,a)') &
          'IMSLINEAR PERIOD data period numbers must be increasing. In &
          &stress period ', kper, ' found period ', this%ionper, &
          ' but the previous block was period ', this%lastonper, '.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      end if
    else if (ierr < 0) then
      ! -- end of file: no more period blocks, keep current settings
      this%ionper = nper + 1
    else
      call this%parser%GetCurrentLine(line)
      write (errmsg, '(3a)') &
        'Looking for BEGIN PERIOD.  Found "', trim(adjustl(line)), &
        '" instead.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
  end if
  !
  ! -- apply this period's overrides when the pending block matches kper;
  !    settings are mutated in place so they persist to later periods
  if (this%ionper == kper) then
    write (this%iout, '(/1x,a,i0)') &
      'PROCESSING LINEAR PERIOD DATA FOR STRESS PERIOD ', kper
    do
      call this%parser%GetNextLine(end_of_block)
      if (end_of_block) exit
      call this%parser%GetStringCaps(keyword)
      select case (keyword)
      case ('INNER_DVCLOSE', 'INNER_RCLOSE', 'INNER_MAXIMUM', &
            'PRECONDITIONER_LEVELS', 'PRECONDITIONER_DROP_TOLERANCE', &
            'RELAXATION_FACTOR')
        call settings%apply_keyword(this%parser, keyword)
        changed = .true.
      case default
        write (errmsg, '(3a)') &
          'Keyword "', trim(keyword), &
          '" is not allowed in an IMSLINEAR PERIOD block.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      end select
    end do
    write (this%iout, '(1x,a,i0)') &
      'END PROCESSING LINEAR PERIOD DATA FOR STRESS PERIOD ', kper
  end if
  end procedure imslinperiod_read_period

  !> @brief Deallocate the period-settings utility
  !<
  module procedure imslinperiod_da
  if (this%active) then
    call this%parser%Clear()
  end if
  this%active = .false.
  this%inunit = 0
  this%ionper = 0
  this%lastonper = 0
  end procedure imslinperiod_da

end submodule ImsLinearPeriodImpl
