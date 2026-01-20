!> @brief This module contains the BudgetModule
!!
!! New entries can be added for each time step, however, the same number of
!! entries must be provided, and they must be provided in the same order.  If not,
!! the module will terminate with an error.
!!
!! Maxsize is required as part of the df method and the arrays will be allocated
!! to maxsize.  If additional entries beyond maxsize are added, the arrays
!! will dynamically increase in size, however, to avoid allocation and copying,
!! it is best to set maxsize large enough up front.
!!
!! vbvl(1, :) contains cumulative rate in
!! vbvl(2, :) contains cumulative rate out
!! vbvl(3, :) contains rate in
!! vbvl(4, :) contains rate out
!! vbnm(:)    contains a LENBUDTXT character text string for each entry
!! rowlabel(:) contains a LENBUDROWLABEL character text string to write as a label for each entry
!!
!<
module BudgetModule

  use KindModule, only: DP, I4B
  use SimModule, only: store_error, count_errors
  use ConstantsModule, only: LINELENGTH, LENBUDTXT, LENBUDROWLABEL, DZERO, &
                             DTWO, DHUNDRED

  implicit none
  private
  public :: BudgetType
  public :: budget_cr
  public :: rate_accumulator

  !> @brief Derived type for the Budget object
  !!
  !! This derived type stores and prints information about a
  !! model budget.
  !!
  !<
  type BudgetType
    integer(I4B), pointer :: msum => null()
    integer(I4B), pointer :: maxsize => null()
    real(DP), pointer :: budperc => null()
    logical, pointer :: written_once => null()
    real(DP), dimension(:, :), pointer :: vbvl => null()
    character(len=LENBUDTXT), dimension(:), pointer, contiguous :: vbnm => null()
    character(len=20), pointer :: bdtype => null()
    character(len=5), pointer :: bddim => null()
    character(len=LENBUDROWLABEL), &
      dimension(:), pointer, contiguous :: rowlabel => null()
    character(len=16), pointer :: labeltitle => null()
    character(len=20), pointer :: bdzone => null()
    logical, pointer :: labeled => null()
    !
    ! -- csv output
    integer(I4B), pointer :: ibudcsv => null()
    integer(I4B), pointer :: icsvheader => null()

  contains
    procedure :: budget_df
    procedure :: budget_ot
    procedure :: budget_da
    procedure :: set_ibudcsv
    procedure :: reset
    procedure :: add_single_entry
    procedure :: add_multi_entry
    generic :: addentry => add_single_entry, add_multi_entry
    procedure :: finalize_step
    procedure :: writecsv
    ! -- private
    procedure :: allocate_scalars
    procedure, private :: allocate_arrays
    procedure, private :: resize
    procedure, private :: write_csv_header
  end type BudgetType

contains

  !> @ brief Create a new budget object
  !!
  !!  Create a new budget object.
  !!
  !<
  subroutine budget_cr(this, name_model)
    ! -- modules
    ! -- dummy
    type(BudgetType), pointer :: this !< BudgetType object
    character(len=*), intent(in) :: name_model !< name of the model
    !
    ! -- Create the object
    allocate (this)
    !
    ! -- Allocate scalars
    call this%allocate_scalars(name_model)
  end subroutine budget_cr

  !> @ brief Define information for this object
  !!
  !!  Allocate arrays and set member variables
  !!
  !<
  subroutine budget_df(this, maxsize, bdtype, bddim, labeltitle, bdzone)
    class(BudgetType) :: this !< BudgetType object
    integer(I4B), intent(in) :: maxsize !< maximum size of budget arrays
    character(len=*), optional :: bdtype !< type of budget, default is VOLUME
    character(len=*), optional :: bddim !< dimensions of terms, default is L**3
    character(len=*), optional :: labeltitle !< budget label, default is PACKAGE NAME
    character(len=*), optional :: bdzone !< corresponding zone, default is ENTIRE MODEL
    !
    ! -- Set values
    this%maxsize = maxsize
    !
    ! -- Allocate arrays
    call this%allocate_arrays()
    !
    ! -- Set the budget type
    if (present(bdtype)) then
      this%bdtype = bdtype
    else
      this%bdtype = 'VOLUME'
    end if
    !
    ! -- Set the budget dimension
    if (present(bddim)) then
      this%bddim = bddim
    else
      this%bddim = 'L**3'
    end if
    !
    ! -- Set the budget zone
    if (present(bdzone)) then
      this%bdzone = bdzone
    else
      this%bdzone = 'ENTIRE MODEL'
    end if
    !
    ! -- Set the label title
    if (present(labeltitle)) then
      this%labeltitle = labeltitle
    else
      this%labeltitle = 'PACKAGE NAME'
    end if
  end subroutine budget_df

  !> @ brief Convert a number to a string
  !!
  !!  This is sometimes needed to avoid numbers that do not fit
  !!  correctly into a text string
  !!
  !<
  subroutine value_to_string(val, string, big, small)
    real(DP), intent(in) :: val !< value to convert
    character(len=*), intent(out) :: string !< string to fill
    real(DP), intent(in) :: big !< big value
    real(DP), intent(in) :: small !< small value
    real(DP) :: absval
    !
    absval = abs(val)
    if (val /= DZERO .and. (absval >= big .or. absval < small)) then
      if (absval >= 1.D100 .or. absval <= 1.D-100) then
        ! -- if exponent has 3 digits, then need to explicitly use the ES
        !    format to force writing the E character
        write (string, '(es17.4E3)') val
      else
        write (string, '(1pe17.4)') val
      end if
    else
      ! -- value is within range where number looks good with F format
      write (string, '(f17.4)') val
    end if
  end subroutine value_to_string

  !> @ brief Output the budget table
  !!
  !!  Write the budget table for the current set of budget
  !!  information.
  !!
  !<
  subroutine budget_ot(this, kstp, kper, iout)
    class(BudgetType) :: this !< BudgetType object
    integer(I4B), intent(in) :: kstp !< time step
    integer(I4B), intent(in) :: kper !< stress period
    integer(I4B), intent(in) :: iout !< output unit number
    character(len=17) :: val1, val2
    integer(I4B) :: msum1, l
    real(DP) :: two, hund, bigvl1, bigvl2, small, &
                totrin, totrot, totvin, totvot, diffr, adiffr, &
                pdiffr, pdiffv, avgrat, diffv, adiffv, avgvol
    !
    ! -- Set constants
    two = 2.d0
    hund = 100.d0
    bigvl1 = 9.99999d11
    bigvl2 = 9.99999d10
    small = 0.1d0
    !
    ! -- Determine number of individual budget entries.
    this%budperc = DZERO
    msum1 = this%msum - 1
    if (msum1 <= 0) return
    !
    ! -- Clear rate and volume accumulators.
    totrin = DZERO
    totrot = DZERO
    totvin = DZERO
    totvot = DZERO
    !
    ! -- Add rates and volumes (in and out) to accumulators.
    do l = 1, msum1
      totrin = totrin + this%vbvl(3, l)
      totrot = totrot + this%vbvl(4, l)
      totvin = totvin + this%vbvl(1, l)
      totvot = totvot + this%vbvl(2, l)
    end do
    !
    ! -- Print time step number and stress period number.
    if (this%labeled) then
      write (iout, 261) trim(adjustl(this%bdtype)), trim(adjustl(this%bdzone)), &
        kstp, kper
      write (iout, 266) trim(adjustl(this%bdtype)), trim(adjustl(this%bddim)), &
        trim(adjustl(this%bddim)), this%labeltitle
    else
      write (iout, 260) trim(adjustl(this%bdtype)), trim(adjustl(this%bdzone)), &
        kstp, kper
      write (iout, 265) trim(adjustl(this%bdtype)), trim(adjustl(this%bddim)), &
        trim(adjustl(this%bddim))
    end if
    !
    ! -- Print individual inflow rates and volumes and their totals.
    do l = 1, msum1
      call value_to_string(this%vbvl(1, l), val1, bigvl1, small)
      call value_to_string(this%vbvl(3, l), val2, bigvl1, small)
      if (this%labeled) then
        write (iout, 276) this%vbnm(l), val1, this%vbnm(l), val2, this%rowlabel(l)
      else
        write (iout, 275) this%vbnm(l), val1, this%vbnm(l), val2
      end if
    end do
    call value_to_string(totvin, val1, bigvl1, small)
    call value_to_string(totrin, val2, bigvl1, small)
    write (iout, 286) val1, val2
    !
    ! -- Print individual outflow rates and volumes and their totals.
    write (iout, 287)
    do l = 1, msum1
      call value_to_string(this%vbvl(2, l), val1, bigvl1, small)
      call value_to_string(this%vbvl(4, l), val2, bigvl1, small)
      if (this%labeled) then
        write (iout, 276) this%vbnm(l), val1, this%vbnm(l), val2, this%rowlabel(l)
      else
        write (iout, 275) this%vbnm(l), val1, this%vbnm(l), val2
      end if
    end do
    call value_to_string(totvot, val1, bigvl1, small)
    call value_to_string(totrot, val2, bigvl1, small)
    write (iout, 298) val1, val2
    !
    ! -- Calculate the difference between inflow and outflow.
    !
    ! -- Calculate difference between rate in and rate out.
    diffr = totrin - totrot
    adiffr = abs(diffr)
    !
    ! -- Calculate percent difference between rate in and rate out.
    pdiffr = DZERO
    avgrat = (totrin + totrot) / two
    if (avgrat /= DZERO) pdiffr = hund * diffr / avgrat
    this%budperc = pdiffr
    !
    ! -- Calculate difference between volume in and volume out.
    diffv = totvin - totvot
    adiffv = abs(diffv)
    !
    ! -- Get percent difference between volume in and volume out.
    pdiffv = DZERO
    avgvol = (totvin + totvot) / two
    if (avgvol /= DZERO) pdiffv = hund * diffv / avgvol
    !
    ! -- Print differences and percent differences between input
    ! -- and output rates and volumes.
    call value_to_string(diffv, val1, bigvl2, small)
    call value_to_string(diffr, val2, bigvl2, small)
    write (iout, 299) val1, val2
    write (iout, 300) pdiffv, pdiffr
    !
    ! -- flush the file
    flush (iout)
    !
    ! -- set written_once to .true.
    this%written_once = .true.
    !
    ! -- formats
260 FORMAT(//2X, a, ' BUDGET FOR ', a, ' AT END OF' &
            , ' TIME STEP', I5, ', STRESS PERIOD', I4 / 2X, 78('-'))
261 FORMAT(//2X, a, ' BUDGET FOR ', a, ' AT END OF' &
            , ' TIME STEP', I5, ', STRESS PERIOD', I4 / 2X, 99('-'))
265 FORMAT(1X, /5X, 'CUMULATIVE ', a, 6X, a, 7X &
           , 'RATES FOR THIS TIME STEP', 6X, a, '/T'/5X, 18('-'), 17X, 24('-') &
           //11X, 'IN:', 38X, 'IN:'/11X, '---', 38X, '---')
266 FORMAT(1X, /5X, 'CUMULATIVE ', a, 6X, a, 7X &
           , 'RATES FOR THIS TIME STEP', 6X, a, '/T', 10X, A16, &
           /5X, 18('-'), 17X, 24('-'), 21X, 16('-') &
           //11X, 'IN:', 38X, 'IN:'/11X, '---', 38X, '---')
275 FORMAT(1X, 3X, A16, ' =', A17, 6X, A16, ' =', A17)
276 FORMAT(1X, 3X, A16, ' =', A17, 6X, A16, ' =', A17, 5X, A)
286 FORMAT(1X, /12X, 'TOTAL IN =', A, 14X, 'TOTAL IN =', A)
287 FORMAT(1X, /10X, 'OUT:', 37X, 'OUT:'/10X, 4('-'), 37X, 4('-'))
298 FORMAT(1X, /11X, 'TOTAL OUT =', A, 13X, 'TOTAL OUT =', A)
299 FORMAT(1X, /12X, 'IN - OUT =', A, 14X, 'IN - OUT =', A)
300 FORMAT(1X, /1X, 'PERCENT DISCREPANCY =', F15.2 &
           , 5X, 'PERCENT DISCREPANCY =', F15.2/)
  end subroutine budget_ot

  !> @ brief Deallocate memory
  !!
  !!  Deallocate budget memory
  !!
  !<
  subroutine budget_da(this)
    class(BudgetType) :: this !< BudgetType object
    !
    ! -- Scalars
    deallocate (this%msum)
    deallocate (this%maxsize)
    deallocate (this%budperc)
    deallocate (this%written_once)
    deallocate (this%labeled)
    deallocate (this%bdtype)
    deallocate (this%bddim)
    deallocate (this%labeltitle)
    deallocate (this%bdzone)
    deallocate (this%ibudcsv)
    deallocate (this%icsvheader)
    !
    ! -- Arrays
    deallocate (this%vbvl)
    deallocate (this%vbnm)
    deallocate (this%rowlabel)
  end subroutine budget_da

  !> @ brief Reset the budget object
  !!
  !!  Reset the budget object in preparation for next set of entries
  !!
  !<
  subroutine reset(this)
    ! -- modules
    ! -- dummy
    class(BudgetType) :: this !< BudgetType object
    ! -- local
    integer(I4B) :: i
    !
    this%msum = 1
    do i = 1, this%maxsize
      this%vbvl(3, i) = DZERO
      this%vbvl(4, i) = DZERO
    end do
  end subroutine reset

  !> @ brief Add a single row of information
  !!
  !!  Add information corresponding to one row in the budget table
  !!    rin the inflow rate
  !!    rout is the outflow rate
  !!    delt is the time step length
  !!    text is the name of the entry
  !!    isupress_accumulate is an optional flag.  If specified as 1, then
  !!      the volume is NOT added to the accumulators on vbvl(1, :) and vbvl(2, :).
  !!    rowlabel is a LENBUDROWLABEL character text entry that is written to the
  !!      right of the table.  It can be used for adding package names to budget
  !!      entries.
  !!
  !<
  subroutine add_single_entry(this, rin, rout, delt, text, &
                              isupress_accumulate, rowlabel)
    ! -- dummy
    class(BudgetType) :: this !< BudgetType object
    real(DP), intent(in) :: rin !< inflow rate
    real(DP), intent(in) :: rout !< outflow rate
    real(DP), intent(in) :: delt !< time step length
    character(len=LENBUDTXT), intent(in) :: text !< name of the entry
    integer(I4B), optional, intent(in) :: isupress_accumulate !< accumulate flag
    character(len=*), optional, intent(in) :: rowlabel !< row label
    ! -- local
    character(len=LINELENGTH) :: errmsg
    character(len=*), parameter :: fmtbuderr = &
      &"('Error in MODFLOW 6.', 'Entries do not match: ', (a), (a) )"
    integer(I4B) :: iscv
    integer(I4B) :: maxsize
    !
    iscv = 0
    if (present(isupress_accumulate)) then
      iscv = isupress_accumulate
    end if
    !
    ! -- ensure budget arrays are large enough
    maxsize = this%msum
    if (maxsize > this%maxsize) then
      call this%resize(maxsize)
    end if
    !
    ! -- If budget has been written at least once, then make sure that the present
    !    text entry matches the last text entry
    if (this%written_once) then
      if (trim(adjustl(this%vbnm(this%msum))) /= trim(adjustl(text))) then
        write (errmsg, fmtbuderr) trim(adjustl(this%vbnm(this%msum))), &
          trim(adjustl(text))
        call store_error(errmsg, terminate=.TRUE.)
      end if
    end if
    !
    this%vbvl(3, this%msum) = rin
    this%vbvl(4, this%msum) = rout
    this%vbnm(this%msum) = adjustr(text)
    if (present(rowlabel)) then
      this%rowlabel(this%msum) = adjustl(rowlabel)
      this%labeled = .true.
    end if
    this%msum = this%msum + 1
  end subroutine add_single_entry

  !> @ brief Add multiple rows of information
  !!
  !!  Add information corresponding to one multiple rows in the budget table
  !!    budterm is an array with inflow in column 1 and outflow in column 2
  !!    delt is the time step length
  !!    budtxt is the name of the entries.  It should have one entry for each
  !!      row in budterm
  !!    isupress_accumulate is an optional flag.  If specified as 1, then
  !!      the volume is NOT added to the accumulators on vbvl(1, :) and vbvl(2, :).
  !!    rowlabel is a LENBUDROWLABEL character text entry that is written to the
  !!      right of the table.  It can be used for adding package names to budget
  !!      entries. For multiple entries, the same rowlabel is used for each entry.
  !!
  !<
  subroutine add_multi_entry(this, budterm, delt, budtxt, &
                             isupress_accumulate, rowlabel)
    ! -- dummy
    class(BudgetType) :: this !< BudgetType object
    real(DP), dimension(:, :), intent(in) :: budterm !< array of budget terms
    real(DP), intent(in) :: delt !< time step length
    character(len=LENBUDTXT), dimension(:), intent(in) :: budtxt !< name of the entries
    integer(I4B), optional, intent(in) :: isupress_accumulate !< suppress accumulate
    character(len=*), optional, intent(in) :: rowlabel !< row label
    ! -- local
    character(len=LINELENGTH) :: errmsg
    character(len=*), parameter :: fmtbuderr = &
      &"('Error in MODFLOW 6.', 'Entries do not match: ', (a), (a) )"
    integer(I4B) :: iscv, i
    integer(I4B) :: nbudterms, maxsize
    !
    iscv = 0
    if (present(isupress_accumulate)) then
      iscv = isupress_accumulate
    end if
    !
    ! -- ensure budget arrays are large enough
    nbudterms = size(budtxt)
    maxsize = this%msum - 1 + nbudterms
    if (maxsize > this%maxsize) then
      call this%resize(maxsize)
    end if
    !
    ! -- Process each of the multi-entry budget terms
    do i = 1, size(budtxt)
      !
      ! -- If budget has been written at least once, then make sure that the present
      !    text entry matches the last text entry
      if (this%written_once) then
        if (trim(adjustl(this%vbnm(this%msum))) /= &
            trim(adjustl(budtxt(i)))) then
          write (errmsg, fmtbuderr) trim(adjustl(this%vbnm(this%msum))), &
            trim(adjustl(budtxt(i)))
          call store_error(errmsg)
        end if
      end if
      !
      this%vbvl(3, this%msum) = budterm(1, i)
      this%vbvl(4, this%msum) = budterm(2, i)
      this%vbnm(this%msum) = adjustr(budtxt(i))
      if (present(rowlabel)) then
        this%rowlabel(this%msum) = adjustl(rowlabel)
        this%labeled = .true.
      end if
      this%msum = this%msum + 1
      !
    end do
    !
    ! -- Check for errors
    if (count_errors() > 0) then
      call store_error('Could not add multi-entry', terminate=.TRUE.)
    end if
  end subroutine add_multi_entry

  !> @ brief Update accumulators
  !!
  !! This must be called before any output is written
  !! in order to update the accumulators in vbvl(1,:)
  !! and vbl(2,:).
  !<
  subroutine finalize_step(this, delt)
    ! -- modules
    ! -- dummy
    class(BudgetType) :: this !< BudgetType object
    real(DP), intent(in) :: delt
    ! -- local
    integer(I4B) :: i
    !
    do i = 1, this%msum - 1
      this%vbvl(1, i) = this%vbvl(1, i) + this%vbvl(3, i) * delt
      this%vbvl(2, i) = this%vbvl(2, i) + this%vbvl(4, i) * delt
    end do
  end subroutine finalize_step

  !> @ brief allocate scalar variables
  !!
  !!  Allocate scalar variables of this budget object
  !!
  !<
  subroutine allocate_scalars(this, name_model)
    ! -- modules
    ! -- dummy
    class(BudgetType) :: this !< BudgetType object
    character(len=*), intent(in) :: name_model !< name of the model
    !
    allocate (this%msum)
    allocate (this%maxsize)
    allocate (this%budperc)
    allocate (this%written_once)
    allocate (this%labeled)
    allocate (this%bdtype)
    allocate (this%bddim)
    allocate (this%labeltitle)
    allocate (this%bdzone)
    allocate (this%ibudcsv)
    allocate (this%icsvheader)
    !
    ! -- Initialize values
    this%msum = 0
    this%maxsize = 0
    this%written_once = .false.
    this%labeled = .false.
    this%bdtype = ''
    this%bddim = ''
    this%labeltitle = ''
    this%bdzone = ''
    this%ibudcsv = 0
    this%icsvheader = 0
  end subroutine allocate_scalars

  !> @ brief allocate array variables
  !!
  !!  Allocate array variables of this budget object
  !!
  !<
  subroutine allocate_arrays(this)
    ! -- modules
    ! -- dummy
    class(BudgetType) :: this !< BudgetType object
    !
    ! -- If redefining, then need to deallocate/reallocate
    if (associated(this%vbvl)) then
      deallocate (this%vbvl)
      nullify (this%vbvl)
    end if
    if (associated(this%vbnm)) then
      deallocate (this%vbnm)
      nullify (this%vbnm)
    end if
    if (associated(this%rowlabel)) then
      deallocate (this%rowlabel)
      nullify (this%rowlabel)
    end if
    !
    ! -- Allocate
    allocate (this%vbvl(4, this%maxsize))
    allocate (this%vbnm(this%maxsize))
    allocate (this%rowlabel(this%maxsize))
    !
    ! -- Initialize values
    this%vbvl(:, :) = DZERO
    this%vbnm(:) = ''
    this%rowlabel(:) = ''
  end subroutine allocate_arrays

  !> @ brief Resize the budget object
  !!
  !!  If the size wasn't allocated to be large enough, then the budget object
  !!  we reallocate itself to a larger size.
  !!
  !<
  subroutine resize(this, maxsize)
    ! -- modules
    ! -- dummy
    class(BudgetType) :: this !< BudgetType object
    integer(I4B), intent(in) :: maxsize !< maximum size
    ! -- local
    real(DP), dimension(:, :), allocatable :: vbvl
    character(len=LENBUDTXT), dimension(:), allocatable :: vbnm
    character(len=LENBUDROWLABEL), dimension(:), allocatable :: rowlabel
    integer(I4B) :: maxsizeold
    !
    ! -- allocate and copy into local storage
    maxsizeold = this%maxsize
    allocate (vbvl(4, maxsizeold))
    allocate (vbnm(maxsizeold))
    allocate (rowlabel(maxsizeold))
    vbvl(:, :) = this%vbvl(:, :)
    vbnm(:) = this%vbnm(:)
    rowlabel(:) = this%rowlabel(:)
    !
    ! -- Set new size and reallocate
    this%maxsize = maxsize
    call this%allocate_arrays()
    !
    ! -- Copy from local back into member variables
    this%vbvl(:, 1:maxsizeold) = vbvl(:, 1:maxsizeold)
    this%vbnm(1:maxsizeold) = vbnm(1:maxsizeold)
    this%rowlabel(1:maxsizeold) = rowlabel(1:maxsizeold)
    !
    ! - deallocate local copies
    deallocate (vbvl)
    deallocate (vbnm)
    deallocate (rowlabel)
  end subroutine resize

  !> @ brief Rate accumulator subroutine
  !!
  !!  Routing for tallying inflows and outflows of an array
  !!
  !<
  subroutine rate_accumulator(flow, rin, rout)
    ! -- modules
    ! -- dummy
    real(DP), dimension(:), contiguous, intent(in) :: flow !< array of flows
    real(DP), intent(out) :: rin !< calculated sum of inflows
    real(DP), intent(out) :: rout !< calculated sum of outflows
    integer(I4B) :: n
    !
    rin = DZERO
    rout = DZERO
    do n = 1, size(flow)
      if (flow(n) < DZERO) then
        rout = rout - flow(n)
      else
        rin = rin + flow(n)
      end if
    end do
  end subroutine rate_accumulator

  !> @ brief Set unit number for csv output file
  !!
  !!  This routine can be used to activate csv output
  !!  by passing in a valid unit number opened for output
  !!
  !<
  subroutine set_ibudcsv(this, ibudcsv)
    ! -- modules
    ! -- dummy
    class(BudgetType) :: this !< BudgetType object
    integer(I4B), intent(in) :: ibudcsv !< unit number for csv budget output
    this%ibudcsv = ibudcsv
  end subroutine set_ibudcsv

  !> @ brief Write csv output
  !!
  !!  This routine will write a row of output to the
  !!  csv file, if it is available for output.  Upon first
  !!  call, it will write the csv header.
  !!
  !<
  subroutine writecsv(this, totim)
    ! -- modules
    ! -- dummy
    class(BudgetType) :: this !< BudgetType object
    real(DP), intent(in) :: totim !< time corresponding to this data
    ! -- local
    integer(I4B) :: i
    real(DP) :: totrin
    real(DP) :: totrout
    real(DP) :: diffr
    real(DP) :: pdiffr
    real(DP) :: avgrat
    !
    if (this%ibudcsv > 0) then
      !
      ! -- write header
      if (this%icsvheader == 0) then
        call this%write_csv_header()
        this%icsvheader = 1
      end if
      !
      ! -- Calculate in and out
      totrin = DZERO
      totrout = DZERO
      do i = 1, this%msum - 1
        totrin = totrin + this%vbvl(3, i)
        totrout = totrout + this%vbvl(4, i)
      end do
      !
      ! -- calculate percent difference
      diffr = totrin - totrout
      pdiffr = DZERO
      avgrat = (totrin + totrout) / DTWO
      if (avgrat /= DZERO) then
        pdiffr = DHUNDRED * diffr / avgrat
      end if
      !
      ! -- write data
      write (this%ibudcsv, '(*(G0,:,","))') &
        totim, &
        (this%vbvl(3, i), i=1, this%msum - 1), &
        (this%vbvl(4, i), i=1, this%msum - 1), &
        totrin, totrout, pdiffr
      !
      ! -- flush the file
      flush (this%ibudcsv)
    end if
  end subroutine writecsv

  !> @ brief Write csv header
  !!
  !!  This routine will write the csv header based on the
  !!  names in vbnm
  !!
  !<
  subroutine write_csv_header(this)
    ! -- modules
    ! -- dummy
    class(BudgetType) :: this !< BudgetType object
    ! -- local
    integer(I4B) :: l
    character(len=LINELENGTH) :: txt, txtl
    write (this%ibudcsv, '(a)', advance='NO') 'time,'
    !
    ! -- first write IN
    do l = 1, this%msum - 1
      txt = this%vbnm(l)
      txtl = ''
      if (this%labeled) then
        txtl = '('//trim(adjustl(this%rowlabel(l)))//')'
      end if
      txt = trim(adjustl(txt))//trim(adjustl(txtl))//'_IN,'
      write (this%ibudcsv, '(a)', advance='NO') trim(adjustl(txt))
    end do
    !
    ! -- then write OUT
    do l = 1, this%msum - 1
      txt = this%vbnm(l)
      txtl = ''
      if (this%labeled) then
        txtl = '('//trim(adjustl(this%rowlabel(l)))//')'
      end if
      txt = trim(adjustl(txt))//trim(adjustl(txtl))//'_OUT,'
      write (this%ibudcsv, '(a)', advance='NO') trim(adjustl(txt))
    end do
    write (this%ibudcsv, '(a)') 'TOTAL_IN,TOTAL_OUT,PERCENT_DIFFERENCE'
  end subroutine write_csv_header

end module BudgetModule
