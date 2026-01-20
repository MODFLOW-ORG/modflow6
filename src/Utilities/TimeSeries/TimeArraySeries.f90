module TimeArraySeriesModule

  use ArrayReadersModule, only: ReadArray
  use BlockParserModule, only: BlockParserType
  use ConstantsModule, only: LINELENGTH, UNDEFINED, STEPWISE, LINEAR, &
                             LENTIMESERIESNAME, LENMODELNAME, DZERO, DONE
  use MathUtilModule, only: is_close
  use InputOutputModule, only: GetUnit, openfile
  use KindModule, only: DP, I4B
  use ListModule, only: ListType
  use ListNodeModule, only: ListNodeType
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, store_error_unit
  use TimeArrayModule, only: TimeArrayType, ConstructTimeArray, &
                             AddTimeArrayToList, CastAsTimeArrayType, &
                             GetTimeArrayFromList
  use, intrinsic :: iso_fortran_env, only: IOSTAT_END

  implicit none
  private
  public :: TimeArraySeriesType, ConstructTimeArraySeries, &
            CastAsTimeArraySeriesType, GetTimeArraySeriesFromList

  type TimeArraySeriesType
    ! -- Public members
    character(len=LENTIMESERIESNAME), public :: Name = ''
    ! -- Private members
    integer(I4B), private :: inunit = 0
    integer(I4B), private :: iout = 0
    integer(I4B), private :: iMethod = UNDEFINED
    real(DP), private :: sfac = DONE
    character(len=LINELENGTH), private :: dataFile = ''
    logical, private :: autoDeallocate = .true.
    type(ListType), pointer, private :: list => null()
    character(len=LENMODELNAME) :: modelname
    type(BlockParserType), private :: parser

  contains

    ! -- Public procedures
    procedure, public :: tas_init
    procedure, public :: GetAverageValues
    procedure, public :: GetInunit
    procedure, public :: da => tas_da
    ! -- Private procedures
    procedure, private :: get_integrated_values
    procedure, private :: get_latest_preceding_node
    procedure, private :: get_values_at_time
    procedure, private :: get_surrounding_records
    procedure, private :: read_next_array
    procedure, private :: DeallocateBackward
  end type TimeArraySeriesType

contains

  ! -- Constructor for TimeArraySeriesType

  !> @brief Allocate a new TimeArraySeriesType object.
  !<
  subroutine ConstructTimeArraySeries(newTas, filename)
    ! -- dummy
    type(TimeArraySeriesType), pointer, intent(out) :: newTas
    character(len=*), intent(in) :: filename
    ! -- local
    logical :: lex
    ! -- formats
10  format('Error: Time-array-series file "', a, '" does not exist.')
    !
    ! -- Allocate a new object of type TimeArraySeriesType
    allocate (newTas)
    allocate (newTas%list)
    !
    ! -- Ensure that input file exists
    inquire (file=filename, exist=lex)
    if (.not. lex) then
      write (errmsg, 10) trim(filename)
      call store_error(errmsg, terminate=.TRUE.)
    end if
    newTas%datafile = filename
  end subroutine ConstructTimeArraySeries

  ! -- Public procedures

  !> @brief Initialize the time array series
  !<
  subroutine tas_init(this, fname, modelname, iout, tasname, autoDeallocate)
    ! -- dummy
    class(TimeArraySeriesType), intent(inout) :: this
    character(len=*), intent(in) :: fname
    character(len=*), intent(in) :: modelname
    integer(I4B), intent(in) :: iout
    character(len=*), intent(inout) :: tasname
    logical, optional, intent(in) :: autoDeallocate
    ! -- local
    integer(I4B) :: istatus
    integer(I4B) :: ierr
    integer(I4B) :: inunit
    character(len=40) :: keyword, keyvalue
    logical :: found, continueread, endOfBlock
    !
    ! -- initialize some variables
    if (present(autoDeallocate)) this%autoDeallocate = autoDeallocate
    this%dataFile = fname
    allocate (this%list)
    !
    ! -- assign members
    this%modelname = modelname
    this%iout = iout
    !
    ! -- open time-array series input file
    inunit = GetUnit()
    this%inunit = inunit
    call openfile(inunit, 0, fname, 'TAS6')
    !
    ! -- initialize block parser
    call this%parser%Initialize(this%inunit, this%iout)
    !
    ! -- read ATTRIBUTES block
    continueread = .false.
    ierr = 0
    !
    ! -- get BEGIN line of ATTRIBUTES block
    call this%parser%GetBlock('ATTRIBUTES', found, ierr, &
                              supportOpenClose=.true.)
    if (.not. found) then
      errmsg = 'Error: Attributes block not found in file: '// &
               trim(fname)
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- parse ATTRIBUTES entries
    do
      ! -- read line from input
      call this%parser%GetNextLine(endOfBlock)
      if (endOfBlock) exit
      !
      ! -- get the keyword
      call this%parser%GetStringCaps(keyword)
      !
      ! -- get the word following the keyword (the key value)
      call this%parser%GetStringCaps(keyvalue)
      select case (keyword)
      case ('NAME')
        this%Name = keyvalue
        tasname = keyvalue
      case ('METHOD')
        select case (keyvalue)
        case ('STEPWISE')
          this%iMethod = STEPWISE
        case ('LINEAR')
          this%iMethod = LINEAR
        case default
          errmsg = 'Unknown interpolation method: "'//trim(keyvalue)//'"'
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      case ('AUTODEALLOCATE')
        this%autoDeallocate = (keyvalue == 'TRUE')
      case ('SFAC')
        read (keyvalue, *, iostat=istatus) this%sfac
        if (istatus /= 0) then
          errmsg = 'Error reading numeric SFAC value from "'//trim(keyvalue) &
                   //'"'
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
      case default
        errmsg = 'Unknown option found in ATTRIBUTES block: "'// &
                 trim(keyword)//'"'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      end select
    end do
    !
    ! -- ensure that NAME and METHOD have been specified
    if (this%Name == '') then
      errmsg = 'Name not specified for time array series in file: '// &
               trim(this%dataFile)
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    if (this%iMethod == UNDEFINED) then
      errmsg = 'Interpolation method not specified for time'// &
               ' array series in file: '//trim(this%dataFile)
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- handle any errors encountered so far
    if (count_errors() > 0) then
      errmsg = 'Error(s) encountered initializing time array series from file: ' &
               //trim(this%dataFile)
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- try to read first time array into linked list
    if (.not. this%read_next_array()) then
      errmsg = 'Error encountered reading time-array data from file: '// &
               trim(this%dataFile)
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
  end subroutine tas_init

  !> @brief Populate an array time-weighted average value for a specified time
  !! span
  !<
  subroutine GetAverageValues(this, nvals, values, time0, time1)
    ! -- dummy
    class(TimeArraySeriesType), intent(inout) :: this
    integer(I4B), intent(in) :: nvals
    real(DP), dimension(nvals), intent(inout) :: values
    real(DP), intent(in) :: time0
    real(DP), intent(in) :: time1
    ! -- local
    integer(I4B) :: i
    real(DP) :: timediff
    !
    timediff = time1 - time0
    if (timediff > 0) then
      call this%get_integrated_values(nvals, values, time0, time1)
      do i = 1, nvals
        values(i) = values(i) / timediff
      end do
    else
      ! -- time0 and time1 are the same, so skip the integration step.
      call this%get_values_at_time(nvals, values, time0)
    end if
  end subroutine GetAverageValues

  !> @brief Return unit number
  !<
  function GetInunit(this)
    ! -- return
    integer(I4B) :: GetInunit
    ! -- dummy
    class(TimeArraySeriesType) :: this
    !
    GetInunit = this%inunit
  end function GetInunit

  ! -- Private procedures

  !> @brief Get surrounding records
  !<
  subroutine get_surrounding_records(this, time, taEarlier, taLater)
    ! -- dummy
    class(TimeArraySeriesType), intent(inout) :: this
    real(DP), intent(in) :: time
    type(TimeArrayType), pointer, intent(inout) :: taEarlier
    type(TimeArrayType), pointer, intent(inout) :: taLater
    ! -- local
    real(DP) :: time0, time1
    type(ListNodeType), pointer :: currNode => null()
    type(ListNodeType), pointer :: node0 => null()
    type(ListNodeType), pointer :: node1 => null()
    type(TimeArrayType), pointer :: ta => null(), ta0 => null(), ta1 => null()
    class(*), pointer :: obj
    !
    taEarlier => null()
    taLater => null()
    !
    if (associated(this%list%firstNode)) then
      currNode => this%list%firstNode
    end if
    !
    ! -- If the next node is earlier than time of interest, advance along
    !    linked list until the next node is later than time of interest.
    do
      if (associated(currNode)) then
        if (associated(currNode%nextNode)) then
          obj => currNode%nextNode%GetItem()
          ta => CastAsTimeArrayType(obj)
          if (ta%taTime <= time) then
            currNode => currNode%nextNode
          else
            exit
          end if
        else
          ! -- read another array
          if (.not. this%read_next_array()) exit
        end if
      else
        exit
      end if
    end do
    !
    if (associated(currNode)) then
      !
      ! -- find earlier record
      node0 => currNode
      obj => node0%GetItem()
      ta0 => CastAsTimeArrayType(obj)
      time0 = ta0%taTime
      do while (time0 > time)
        if (associated(node0%prevNode)) then
          node0 => node0%prevNode
          obj => node0%GetItem()
          ta0 => CastAsTimeArrayType(obj)
          time0 = ta0%taTime
        else
          exit
        end if
      end do
      !
      ! -- find later record
      node1 => currNode
      obj => node1%GetItem()
      ta1 => CastAsTimeArrayType(obj)
      time1 = ta1%taTime
      do while (time1 < time)
        if (associated(node1%nextNode)) then
          node1 => node1%nextNode
          obj => node1%GetItem()
          ta1 => CastAsTimeArrayType(obj)
          time1 = ta1%taTime
        else
          ! -- get next array
          if (.not. this%read_next_array()) then
            ! -- end of file reached, so exit loop
            exit
          end if
        end if
      end do
      !
    end if
    !
    if (time0 <= time) taEarlier => ta0
    if (time1 >= time) taLater => ta1
  end subroutine get_surrounding_records

  !> @brief Read next time array from input file and append to list
  !<
  logical function read_next_array(this)
    ! -- modules
    use ConstantsModule, only: LENMEMPATH
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path
    ! -- dummy
    class(TimeArraySeriesType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i, ierr, istart, istat, istop, lloc, nrow, ncol, &
                    nodesperlayer
    logical :: lopen, isFound
    type(TimeArrayType), pointer :: ta => null()
    character(len=LENMEMPATH) :: mempath
    integer(I4B), dimension(:), contiguous, pointer :: mshape
    !
    ! -- initialize
    istart = 1
    istat = 0
    istop = 1
    lloc = 1
    nullify (mshape)
    !
    ! -- create mempath
    mempath = create_mem_path(component=this%modelname, subcomponent='DIS')
    !
    ! -- set mshape pointer
    call mem_setptr(mshape, 'MSHAPE', mempath)
    !
    ! Get dimensions for supported discretization type
    if (size(mshape) == 2) then
      nodesperlayer = mshape(2)
      nrow = 1
      ncol = mshape(2)
    else if (size(mshape) == 3) then
      nodesperlayer = mshape(2) * mshape(3)
      nrow = mshape(2)
      ncol = mshape(3)
    else
      errmsg = 'Time array series is not supported for selected &
               &discretization type.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    !
    read_next_array = .false.
    inquire (unit=this%inunit, opened=lopen)
    if (lopen) then
      call ConstructTimeArray(ta, this%modelname)
      ! -- read a time and an array from the input file
      ! -- Get a TIME block and read the time
      call this%parser%GetBlock('TIME', isFound, ierr, &
                                supportOpenClose=.false.)
      if (isFound) then
        ta%taTime = this%parser%GetDouble()
        ! -- Read the array
        call ReadArray(this%parser%iuactive, ta%taArray, this%Name, &
                       size(mshape), ncol, nrow, 1, nodesperlayer, &
                       this%iout, 0, 0)
        !
        ! -- multiply values by sfac
        do i = 1, nodesperlayer
          ta%taArray(i) = ta%taArray(i) * this%sfac
        end do
        !
        ! -- append the new time array to the list
        call AddTimeArrayToList(this%list, ta)
        read_next_array = .true.
        !
        ! -- make sure block is closed
        call this%parser%terminateblock()
      end if
    end if
  end function read_next_array

  !> @brief Return an array of values for a specified time, same units as
  !! time-series values
  !<
  subroutine get_values_at_time(this, nvals, values, time)
    ! -- dummy
    class(TimeArraySeriesType), intent(inout) :: this
    integer(I4B), intent(in) :: nvals
    real(DP), dimension(nvals), intent(inout) :: values
    real(DP), intent(in) :: time ! time of interest
    ! -- local
    integer(I4B) :: i, ierr
    real(DP) :: ratio, time0, time1, timediff, timediffi, val0, val1, &
                valdiff
    type(TimeArrayType), pointer :: taEarlier => null()
    type(TimeArrayType), pointer :: taLater => null()
    ! -- formats
10  format('Error getting array at time ', g10.3, &
           ' for time-array series "', a, '"')
    !
    ierr = 0
    call this%get_surrounding_records(time, taEarlier, taLater)
    if (associated(taEarlier)) then
      if (associated(taLater)) then
        ! -- values are available for both earlier and later times
        if (this%iMethod == STEPWISE) then
          ! -- Just populate values from elements of earlier time array
          values = taEarlier%taArray
        elseif (this%iMethod == LINEAR) then
          ! -- perform linear interpolation
          time0 = taEarlier%taTime
          time1 = taLater%tatime
          timediff = time1 - time0
          timediffi = time - time0
          if (timediff > 0) then
            ratio = timediffi / timediff
          else
            ! -- should not happen if TS does not contain duplicate times
            ratio = 0.5d0
          end if
          ! -- Iterate through all elements and perform interpolation.
          do i = 1, nvals
            val0 = taEarlier%taArray(i)
            val1 = taLater%taArray(i)
            valdiff = val1 - val0
            values(i) = val0 + (ratio * valdiff)
          end do
        else
          ierr = 1
        end if
      else
        if (is_close(taEarlier%taTime, time)) then
          values = taEarlier%taArray
        else
          ! -- Only earlier time is available, and it is not time of interest;
          !    however, if method is STEPWISE, use value for earlier time.
          if (this%iMethod == STEPWISE) then
            values = taEarlier%taArray
          else
            ierr = 1
          end if
        end if
      end if
    else
      if (associated(taLater)) then
        if (is_close(taLater%taTime, time)) then
          values = taLater%taArray
        else
          ! -- only later time is available, and it is not time of interest
          ierr = 1
        end if
      else
        ! -- Neither earlier nor later time is available.
        !    This should never happen!
        ierr = 1
      end if
    end if
    !
    if (ierr > 0) then
      write (errmsg, 10) time, trim(this%Name)
      call store_error(errmsg)
      call store_error_unit(this%inunit)
    end if
  end subroutine get_values_at_time

  !> @brief Populates an array with integrated values for a specified time span
  !!
  !! Units: (ts-value-unit)*time
  !<
  subroutine get_integrated_values(this, nvals, values, time0, time1)
    ! -- dummy
    class(TimeArraySeriesType), intent(inout) :: this
    integer(I4B), intent(in) :: nvals
    real(DP), dimension(nvals), intent(inout) :: values
    real(DP), intent(in) :: time0
    real(DP), intent(in) :: time1
    ! -- local
    integer(I4B) :: i
    real(DP) :: area, currTime, nextTime, ratio0, ratio1, t0, &
                t01, t1, timediff, value, value0, value1, valuediff
    logical :: ldone
    type(ListNodeType), pointer :: precNode => null()
    type(ListNodeType), pointer :: currNode => null(), nextNode => null()
    type(TimeArrayType), pointer :: currRecord => null(), nextRecord => null()
    class(*), pointer :: currObj => null(), nextObj => null()
    ! -- formats
10  format('Error encountered while performing integration', &
           ' for time-array series "', a, '" for time interval: ', &
           g12.5, ' to ', g12.5)
    !
    values = DZERO
    value = DZERO
    ldone = .false.
    t1 = -DONE
    call this%get_latest_preceding_node(time0, precNode)
    if (associated(precNode)) then
      currNode => precNode
      do while (.not. ldone)
        currObj => currNode%GetItem()
        currRecord => CastAsTimeArrayType(currObj)
        currTime = currRecord%taTime
        if (currTime < time1) then
          if (.not. associated(currNode%nextNode)) then
            ! -- try to read the next array
            if (.not. this%read_next_array()) then
              write (errmsg, 10) trim(this%Name), time0, time1
              call store_error(errmsg)
              call store_error_unit(this%inunit)
            end if
          end if
          if (associated(currNode%nextNode)) then
            nextNode => currNode%nextNode
            nextObj => nextNode%GetItem()
            nextRecord => CastAsTimeArrayType(nextObj)
            nextTime = nextRecord%taTime
            ! -- determine lower and upper limits of time span of interest
            !    within current interval
            if (currTime >= time0) then
              t0 = currTime
            else
              t0 = time0
            end if
            if (nextTime <= time1) then
              t1 = nextTime
            else
              t1 = time1
            end if
            ! -- For each element, find area of rectangle
            !    or trapezoid delimited by t0 and t1.
            t01 = t1 - t0
            select case (this%iMethod)
            case (STEPWISE)
              do i = 1, nvals
                ! -- compute area of a rectangle
                value0 = currRecord%taArray(i)
                area = value0 * t01
                ! -- add area to integrated value
                values(i) = values(i) + area
              end do
            case (LINEAR)
              do i = 1, nvals
                ! -- compute area of a trapezoid
                timediff = nextTime - currTime
                ratio0 = (t0 - currTime) / timediff
                ratio1 = (t1 - currTime) / timediff
                valuediff = nextRecord%taArray(i) - currRecord%taArray(i)
                value0 = currRecord%taArray(i) + ratio0 * valuediff
                value1 = currRecord%taArray(i) + ratio1 * valuediff
                area = 0.5d0 * t01 * (value0 + value1)
                ! -- add area to integrated value
                values(i) = values(i) + area
              end do
            end select
          else
            write (errmsg, 10) trim(this%Name), time0, time1
            call store_error(errmsg)
            call store_error('(Probable programming error)', terminate=.TRUE.)
          end if
        else
          ! Current node time = time1 so should be done
          ldone = .true.
        end if
        !
        ! -- Are we done yet?
        if (t1 >= time1) then
          ldone = .true.
        else
          if (.not. associated(currNode%nextNode)) then
            ! -- try to read the next array
            if (.not. this%read_next_array()) then
              write (errmsg, 10) trim(this%Name), time0, time1
              call store_error(errmsg)
              call this%parser%StoreErrorUnit()
            end if
          end if
          if (associated(currNode%nextNode)) then
            currNode => currNode%nextNode
          else
            write (errmsg, 10) trim(this%Name), time0, time1
            call store_error(errmsg)
            call store_error('(Probable programming error)', terminate=.TRUE.)
          end if
        end if
      end do
    end if
    !
    if (this%autoDeallocate) then
      if (associated(precNode)) then
        if (associated(precNode%prevNode)) then
          call this%DeallocateBackward(precNode%prevNode)
        end if
      end if
    end if
  end subroutine get_integrated_values

  !> @brief Deallocate fromNode and all previous nodes in list;
  !! reassign firstNode
  !<
  subroutine DeallocateBackward(this, fromNode)
    ! -- dummy
    class(TimeArraySeriesType), intent(inout) :: this
    type(ListNodeType), pointer, intent(inout) :: fromNode
    !
    ! -- local
    type(ListNodeType), pointer :: current => null()
    type(ListNodeType), pointer :: prev => null()
    type(TimeArrayType), pointer :: ta => null()
    class(*), pointer :: obj => null()
    !
    if (associated(fromNode)) then
      ! -- reassign firstNode
      if (associated(fromNode%nextNode)) then
        this%list%firstNode => fromNode%nextNode
      else
        this%list%firstNode => null()
      end if
      ! -- deallocate fromNode and all previous nodes
      current => fromNode
      do while (associated(current))
        prev => current%prevNode
        obj => current%GetItem()
        ta => CastAsTimeArrayType(obj)
        ! -- Deallocate the contents of this time array,
        !    then remove it from the list
        call ta%da()
        call this%list%RemoveNode(current, .true.)
        current => prev
      end do
      fromNode => null()
    end if
  end subroutine DeallocateBackward

  !> @brief Return pointer to ListNodeType object for the node representing
  !! the latest preceding time in the time series
  !<
  subroutine get_latest_preceding_node(this, time, tslNode)
    ! -- dummy
    class(TimeArraySeriesType), intent(inout) :: this
    real(DP), intent(in) :: time
    type(ListNodeType), pointer, intent(inout) :: tslNode
    ! -- local
    real(DP) :: time0
    type(ListNodeType), pointer :: currNode => null()
    type(ListNodeType), pointer :: node0 => null()
    type(TimeArrayType), pointer :: ta => null()
    type(TimeArrayType), pointer :: ta0 => null()
    class(*), pointer :: obj => null()
    !
    tslNode => null()
    if (associated(this%list%firstNode)) then
      currNode => this%list%firstNode
    else
      call store_error('probable programming error in &
                       &get_latest_preceding_node', &
                       terminate=.TRUE.)
    end if
    !
    continue
    ! -- If the next node is earlier than time of interest, advance along
    !    linked list until the next node is later than time of interest.
    do
      if (associated(currNode)) then
        if (associated(currNode%nextNode)) then
          obj => currNode%nextNode%GetItem()
          ta => CastAsTimeArrayType(obj)
          if (ta%taTime < time .or. is_close(ta%taTime, time)) then
            currNode => currNode%nextNode
          else
            exit
          end if
        else
          ! -- read another record
          if (.not. this%read_next_array()) exit
        end if
      else
        exit
      end if
    end do
    !
    if (associated(currNode)) then
      !
      ! -- find earlier record
      node0 => currNode
      obj => node0%GetItem()
      ta0 => CastAsTimeArrayType(obj)
      time0 = ta0%taTime
      do while (time0 > time)
        if (associated(node0%prevNode)) then
          node0 => node0%prevNode
          obj => node0%GetItem()
          ta0 => CastAsTimeArrayType(obj)
          time0 = ta0%taTime
        else
          exit
        end if
      end do
    end if
    !
    if (time0 <= time) tslNode => node0
  end subroutine get_latest_preceding_node

  !> @brief Deallocate memory
  !<
  subroutine tas_da(this)
    ! -- dummy
    class(TimeArraySeriesType), intent(inout) :: this
    ! -- local
    integer :: i, n
    type(TimeArrayType), pointer :: ta => null()
    !
    ! -- Deallocate contents of each time array in list
    n = this%list%Count()
    do i = 1, n
      ta => GetTimeArrayFromList(this%list, i)
      call ta%da()
    end do
    !
    ! -- Deallocate the list of time arrays
    call this%list%Clear(.true.)
    deallocate (this%list)
  end subroutine tas_da

  ! -- Procedures not type-bound

  !> @brief Cast an unlimited polymorphic object as class(TimeArraySeriesType)
  !<
  function CastAsTimeArraySeriesType(obj) result(res)
    ! -- dummy
    class(*), pointer, intent(inout) :: obj
    ! -- return
    type(TimeArraySeriesType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    type is (TimeArraySeriesType)
      res => obj
    end select
  end function CastAsTimeArraySeriesType

  !> @brief Get time array from list
  !<
  function GetTimeArraySeriesFromList(list, indx) result(res)
    ! -- dummy
    type(ListType), intent(inout) :: list
    integer, intent(in) :: indx
    ! -- return
    type(TimeArraySeriesType), pointer :: res
    ! -- local
    class(*), pointer :: obj
    !
    obj => list%GetItem(indx)
    res => CastAsTimeArraySeriesType(obj)
  end function GetTimeArraySeriesFromList

end module TimeArraySeriesModule
