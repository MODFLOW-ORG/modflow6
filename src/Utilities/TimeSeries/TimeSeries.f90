module TimeSeriesModule

  use KindModule, only: DP, I4B
  use BlockParserModule, only: BlockParserType
  use ConstantsModule, only: LINELENGTH, UNDEFINED, STEPWISE, LINEAR, &
                             LINEAREND, LENTIMESERIESNAME, LENHUGELINE, &
                             DZERO, DONE, DNODATA
  use MathUtilModule, only: is_close
  use InputOutputModule, only: GetUnit, openfile, ParseLine, upcase
  use ListModule, only: ListType
  use ListNodeModule, only: ListNodeType
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, &
                       store_error_unit
  use TimeSeriesRecordModule, only: TimeSeriesRecordType, &
                                    ConstructTimeSeriesRecord, &
                                    CastAsTimeSeriesRecordType, &
                                    AddTimeSeriesRecordToList

  private
  public :: TimeSeriesType, TimeSeriesFileType, ConstructTimeSeriesFile, &
            TimeSeriesContainerType, AddTimeSeriesFileToList, &
            GetTimeSeriesFileFromList, CastAsTimeSeriesFileClass, &
            SameTimeSeries

  type TimeSeriesType
    ! -- Public members
    integer(I4B), public :: iMethod = UNDEFINED
    character(len=LENTIMESERIESNAME), public :: Name = ''
    ! -- Private members
    real(DP), private :: sfac = DONE
    logical, public :: autoDeallocate = .true.
    type(ListType), pointer, private :: list => null()
    class(TimeSeriesFileType), pointer, private :: tsfile => null()

  contains

    ! -- Public procedures
    procedure, public :: AddTimeSeriesRecord
    procedure, public :: Clear
    procedure, public :: FindLatestTime
    procedure, public :: get_surrounding_records
    procedure, public :: get_surrounding_nodes
    procedure, public :: GetCurrentTimeSeriesRecord
    procedure, public :: GetNextTimeSeriesRecord
    procedure, public :: GetPreviousTimeSeriesRecord
    procedure, public :: GetTimeSeriesRecord
    procedure, public :: GetValue
    procedure, public :: InitializeTimeSeries => initialize_time_series
    procedure, public :: InsertTsr
    procedure, public :: Reset
    ! -- Private procedures
    procedure, private :: da => ts_da
    procedure, private :: get_average_value
    procedure, private :: get_integrated_value
    procedure, private :: get_latest_preceding_node
    procedure, private :: get_value_at_time
    procedure, private :: initialize_time_series
    procedure, private :: read_next_record
  end type TimeSeriesType

  type TimeSeriesFileType
    ! -- Private members
    integer(I4B), public :: inunit = 0
    integer(I4B), public :: iout = 0
    integer(I4B), public :: nTimeSeries = 0
    logical, public :: finishedReading = .false.
    character(len=LINELENGTH), public :: datafile = ''
    type(TimeSeriesType), dimension(:), &
      pointer, contiguous, public :: timeSeries => null()
    type(BlockParserType), pointer, public :: parser

  contains

    ! -- Public procedures
    procedure, public :: Count
    procedure, public :: Initializetsfile
    procedure, public :: GetTimeSeries
    procedure, public :: da => tsf_da
    ! -- Private procedures
    procedure, private :: read_tsfile_line
  end type TimeSeriesFileType

  type TimeSeriesContainerType
    ! -- Public members
    type(TimeSeriesType), pointer, public :: timeSeries => null()
  end type TimeSeriesContainerType

contains

  ! -- non-type-bound procedures

  !> @brief Construct time series file
  !<
  subroutine ConstructTimeSeriesFile(newTimeSeriesFile)
    ! -- dummy
    type(TimeSeriesFileType), pointer, intent(inout) :: newTimeSeriesFile
    !
    allocate (newTimeSeriesFile)
    allocate (newTimeSeriesFile%parser)
  end subroutine ConstructTimeSeriesFile

  !> @brief Cast an unlimited polymorphic object as class(TimeSeriesFileType)
  !<
  function CastAsTimeSeriesFileType(obj) result(res)
    ! -- dummy
    class(*), pointer, intent(inout) :: obj
    ! -- return
    type(TimeSeriesFileType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    type is (TimeSeriesFileType)
      res => obj
    end select
  end function CastAsTimeSeriesFileType

  !> @brief Cast an unlimited polymorphic object as class(TimeSeriesFileType)
  !<
  function CastAsTimeSeriesFileClass(obj) result(res)
    ! -- dummy
    class(*), pointer, intent(inout) :: obj
    ! -- return
    type(TimeSeriesFileType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    class is (TimeSeriesFileType)
      res => obj
    end select
  end function CastAsTimeSeriesFileClass

  !> @brief Add time series file to list
  !<
  subroutine AddTimeSeriesFileToList(list, tsfile)
    ! -- dummy
    type(ListType), intent(inout) :: list
    class(TimeSeriesFileType), pointer, intent(inout) :: tsfile
    ! -- local
    class(*), pointer :: obj => null()
    !
    obj => tsfile
    call list%Add(obj)
  end subroutine AddTimeSeriesFileToList

  !> @brief Get time series from list
  !<
  function GetTimeSeriesFileFromList(list, idx) result(res)
    ! -- dummy
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: idx
    ! -- return
    type(TimeSeriesFileType), pointer :: res
    ! -- local
    class(*), pointer :: obj => null()
    !
    obj => list%GetItem(idx)
    res => CastAsTimeSeriesFileType(obj)
    !
    if (.not. associated(res)) then
      res => CastAsTimeSeriesFileClass(obj)
    end if
  end function GetTimeSeriesFileFromList

  !> @brief Compare two time series; if they are identical, return true
  !<
  function SameTimeSeries(ts1, ts2) result(same)
    ! -- dummy
    type(TimeSeriesType), intent(in) :: ts1
    type(TimeSeriesType), intent(in) :: ts2
    ! -- return
    logical :: same
    ! -- local
    integer :: i, n1, n2
    type(TimeSeriesRecordType), pointer :: tsr1, tsr2
    !
    same = .false.
    n1 = ts1%list%Count()
    n2 = ts2%list%Count()
    if (n1 /= n2) return
    !
    call ts1%Reset()
    call ts2%Reset()
    !
    do i = 1, n1
      tsr1 => ts1%GetNextTimeSeriesRecord()
      tsr2 => ts2%GetNextTimeSeriesRecord()
      if (tsr1%tsrTime /= tsr2%tsrTime) return
      if (tsr1%tsrValue /= tsr2%tsrValue) return
    end do
    !
    same = .true.
  end function SameTimeSeries

  ! Type-bound procedures of TimeSeriesType

  !> @brief Get time series value
  !!
  !! If iMethod is STEPWISE or LINEAR:
  !!     Return a time-weighted average value for a specified time span.
  !! If iMethod is LINEAREND:
  !!     Return value at time1. Time0 argument is ignored.
  !! Units: (ts-value-unit)
  !<
  function GetValue(this, time0, time1, extendToEndOfSimulation)
    ! -- return
    real(DP) :: GetValue
    ! -- dummy
    class(TimeSeriesType), intent(inout) :: this
    real(DP), intent(in) :: time0
    real(DP), intent(in) :: time1
    logical, intent(in), optional :: extendToEndOfSimulation
    ! -- local
    logical :: extend
    !
    if (present(extendToEndOfSimulation)) then
      extend = extendToEndOfSimulation
    else
      extend = .false.
    end if
    !
    select case (this%iMethod)
    case (STEPWISE, LINEAR)
      GetValue = this%get_average_value(time0, time1, extend)
    case (LINEAREND)
      GetValue = this%get_value_at_time(time1, extend)
    end select
  end function GetValue

  !> @brief Initialize time series
  !!
  !! Open time-series file and read options and first time-series record.
  !<
  subroutine initialize_time_series(this, tsfile, name, autoDeallocate)
    ! -- dummy
    class(TimeSeriesType), intent(inout) :: this
    class(TimeSeriesFileType), target :: tsfile
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: autoDeallocate
    ! -- local
    character(len=LENTIMESERIESNAME) :: tsNameTemp
    !
    ! -- Assign the time-series tsfile, name, and autoDeallocate
    this%tsfile => tsfile
    ! Store time-series name as all caps
    tsNameTemp = name
    call UPCASE(tsNameTemp)
    this%Name = tsNameTemp
    !
    this%iMethod = UNDEFINED
    !
    if (present(autoDeallocate)) this%autoDeallocate = autoDeallocate
    !
    ! -- allocate the list
    allocate (this%list)
    !
    ! -- ensure that NAME has been specified
    if (this%Name == '') then
      errmsg = 'Name not specified for time series.'
      call store_error(errmsg, terminate=.TRUE.)
    end if
  end subroutine initialize_time_series

  !> @brief Get surrounding records
  !<
  subroutine get_surrounding_records(this, time, tsrecEarlier, tsrecLater)
    ! -- dummy
    class(TimeSeriesType), intent(inout) :: this
    real(DP), intent(in) :: time
    type(TimeSeriesRecordType), pointer, intent(inout) :: tsrecEarlier
    type(TimeSeriesRecordType), pointer, intent(inout) :: tsrecLater
    ! -- local
    real(DP) :: time0, time1
    type(ListNodeType), pointer :: currNode => null()
    type(ListNodeType), pointer :: tsNode0 => null()
    type(ListNodeType), pointer :: tsNode1 => null()
    type(TimeSeriesRecordType), pointer :: tsr => null(), tsrec0 => null()
    type(TimeSeriesRecordType), pointer :: tsrec1 => null()
    class(*), pointer :: obj => null()
    !
    tsrecEarlier => null()
    tsrecLater => null()
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
          tsr => CastAsTimeSeriesRecordType(obj)
          if (tsr%tsrTime < time .and. .not. is_close(tsr%tsrTime, time)) then
            currNode => currNode%nextNode
          else
            exit
          end if
        else
          ! -- read another record
          if (.not. this%read_next_record()) exit
        end if
      else
        exit
      end if
    end do
    !
    if (associated(currNode)) then
      !
      ! -- find earlier record
      tsNode0 => currNode
      obj => tsNode0%GetItem()
      tsrec0 => CastAsTimeSeriesRecordType(obj)
      time0 = tsrec0%tsrTime
      do while (time0 > time)
        if (associated(tsNode0%prevNode)) then
          tsNode0 => tsNode0%prevNode
          obj => tsNode0%GetItem()
          tsrec0 => CastAsTimeSeriesRecordType(obj)
          time0 = tsrec0%tsrTime
        else
          exit
        end if
      end do
      !
      ! -- find later record
      tsNode1 => currNode
      obj => tsNode1%GetItem()
      tsrec1 => CastAsTimeSeriesRecordType(obj)
      time1 = tsrec1%tsrTime
      do while (time1 < time .and. .not. is_close(time1, time))
        if (associated(tsNode1%nextNode)) then
          tsNode1 => tsNode1%nextNode
          obj => tsNode1%GetItem()
          tsrec1 => CastAsTimeSeriesRecordType(obj)
          time1 = tsrec1%tsrTime
        else
          ! -- get next record
          if (.not. this%read_next_record()) then
            ! -- end of file reached, so exit loop
            exit
          end if
        end if
      end do
      !
    end if
    !
    if (time0 < time .or. is_close(time0, time)) tsrecEarlier => tsrec0
    if (time1 > time .or. is_close(time1, time)) tsrecLater => tsrec1
  end subroutine get_surrounding_records

  !> @brief Get surrounding nodes
  !!
  !! This subroutine is for working with time series already entirely stored
  !! in memory -- it does not read data from a file.
  !<
  subroutine get_surrounding_nodes(this, time, nodeEarlier, nodeLater)
    ! -- dummy
    class(TimeSeriesType), intent(inout) :: this
    real(DP), intent(in) :: time
    type(ListNodeType), pointer, intent(inout) :: nodeEarlier
    type(ListNodeType), pointer, intent(inout) :: nodeLater
    ! -- local
    real(DP) :: time0, time1
    type(ListNodeType), pointer :: currNode => null()
    type(ListNodeType), pointer :: tsNode0 => null()
    type(ListNodeType), pointer :: tsNode1 => null()
    type(TimeSeriesRecordType), pointer :: tsr => null(), tsrec0 => null()
    type(TimeSeriesRecordType), pointer :: tsrec1 => null()
    type(TimeSeriesRecordType), pointer :: tsrecEarlier
    type(TimeSeriesRecordType), pointer :: tsrecLater
    class(*), pointer :: obj => null()
    !
    tsrecEarlier => null()
    tsrecLater => null()
    nodeEarlier => null()
    nodeLater => null()
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
          tsr => CastAsTimeSeriesRecordType(obj)
          if (tsr%tsrTime < time .and. .not. is_close(tsr%tsrTime, time)) then
            currNode => currNode%nextNode
          else
            exit
          end if
        else
          exit
        end if
      else
        exit
      end if
    end do
    !
    if (associated(currNode)) then
      !
      ! -- find earlier record
      tsNode0 => currNode
      obj => tsNode0%GetItem()
      tsrec0 => CastAsTimeSeriesRecordType(obj)
      time0 = tsrec0%tsrTime
      do while (time0 > time)
        if (associated(tsNode0%prevNode)) then
          tsNode0 => tsNode0%prevNode
          obj => tsNode0%GetItem()
          tsrec0 => CastAsTimeSeriesRecordType(obj)
          time0 = tsrec0%tsrTime
        else
          exit
        end if
      end do
      !
      ! -- find later record
      tsNode1 => currNode
      obj => tsNode1%GetItem()
      tsrec1 => CastAsTimeSeriesRecordType(obj)
      time1 = tsrec1%tsrTime
      do while (time1 < time .and. .not. is_close(time1, time))
        if (associated(tsNode1%nextNode)) then
          tsNode1 => tsNode1%nextNode
          obj => tsNode1%GetItem()
          tsrec1 => CastAsTimeSeriesRecordType(obj)
          time1 = tsrec1%tsrTime
        else
          exit
        end if
      end do
      !
    end if
    !
    if (time0 < time .or. is_close(time0, time)) then
      tsrecEarlier => tsrec0
      nodeEarlier => tsNode0
    end if
    if (time1 > time .or. is_close(time1, time)) then
      tsrecLater => tsrec1
      nodeLater => tsNode1
    end if
  end subroutine get_surrounding_nodes

  !> @brief Read next record
  !!
  !! Read next time-series record from input file
  !<
  logical function read_next_record(this)
    ! -- dummy
    class(TimeSeriesType), intent(inout) :: this
    !
    ! -- If we have already encountered the end of the TIMESERIES block, do not try to read any further
    if (this%tsfile%finishedReading) then
      read_next_record = .false.
      return
    end if
    !
    read_next_record = this%tsfile%read_tsfile_line()
    if (.not. read_next_record) then
      this%tsfile%finishedReading = .true.
    end if
  end function read_next_record

  !> @brief Get value for a time
  !!
  !! Return a value for a specified time, same units as time-series values
  !<
  function get_value_at_time(this, time, extendToEndOfSimulation)
    ! -- return
    real(DP) :: get_value_at_time
    ! -- dummy
    class(TimeSeriesType), intent(inout) :: this
    real(DP), intent(in) :: time ! time of interest
    logical, intent(in) :: extendToEndOfSimulation
    ! -- local
    integer(I4B) :: ierr
    real(DP) :: ratio, time0, time1, timediff, timediffi, val0, val1, &
                valdiff
    type(TimeSeriesRecordType), pointer :: tsrEarlier => null()
    type(TimeSeriesRecordType), pointer :: tsrLater => null()
    ! -- formats
10  format('Error getting value at time ', g10.3, ' for time series "', a, '"')
    !
    ierr = 0
    call this%get_surrounding_records(time, tsrEarlier, tsrLater)
    if (associated(tsrEarlier)) then
      if (associated(tsrLater)) then
        ! -- values are available for both earlier and later times
        if (this%iMethod == STEPWISE) then
          get_value_at_time = tsrEarlier%tsrValue
        elseif (this%iMethod == LINEAR .or. this%iMethod == LINEAREND) then
          ! -- For get_value_at_time, result is the same for either
          !    linear method.
          ! -- Perform linear interpolation.
          time0 = tsrEarlier%tsrTime
          time1 = tsrLater%tsrtime
          timediff = time1 - time0
          timediffi = time - time0
          if (timediff > 0) then
            ratio = timediffi / timediff
          else
            ! -- should not happen if TS does not contain duplicate times
            ratio = 0.5d0
          end if
          val0 = tsrEarlier%tsrValue
          val1 = tsrLater%tsrValue
          valdiff = val1 - val0
          get_value_at_time = val0 + (ratio * valdiff)
        else
          ierr = 1
        end if
      else
        if (extendToEndOfSimulation .or. is_close(tsrEarlier%tsrTime, time)) then
          get_value_at_time = tsrEarlier%tsrValue
        else
          ! -- Only earlier time is available, and it is not time of interest;
          !    however, if method is STEPWISE, use value for earlier time.
          if (this%iMethod == STEPWISE) then
            get_value_at_time = tsrEarlier%tsrValue
          else
            ierr = 1
          end if
        end if
      end if
    else
      if (associated(tsrLater)) then
        if (is_close(tsrLater%tsrTime, time)) then
          get_value_at_time = tsrLater%tsrValue
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
      call store_error(errmsg, terminate=.TRUE.)
    end if
  end function get_value_at_time

  !> @brief Get integrated value
  !!
  !! Return an integrated value for a specified time span.
  !! Units: (ts-value-unit)*time
  !<
  function get_integrated_value(this, time0, time1, extendToEndOfSimulation)
    ! -- return
    real(DP) :: get_integrated_value
    ! -- dummy
    class(TimeSeriesType), intent(inout) :: this
    real(DP), intent(in) :: time0
    real(DP), intent(in) :: time1
    logical, intent(in) :: extendToEndOfSimulation
    ! -- local
    real(DP) :: area, currTime, nextTime, ratio0, ratio1, t0, t01, t1, &
                timediff, value, value0, value1, valuediff, currVal, nextVal
    logical :: ldone, lprocess
    type(ListNodeType), pointer :: tslNodePreceding => null()
    type(ListNodeType), pointer :: currNode => null(), nextNode => null()
    type(TimeSeriesRecordType), pointer :: currRecord => null()
    type(TimeSeriesRecordType), pointer :: nextRecord => null()
    class(*), pointer :: currObj => null(), nextObj => null()
    ! -- formats
10  format('Error encountered while performing integration', &
           ' for time series "', a, '" for time interval: ', g12.5, ' to ', g12.5)
    !
    value = DZERO
    ldone = .false.
    t1 = -DONE
    call this%get_latest_preceding_node(time0, tslNodePreceding)
    if (associated(tslNodePreceding)) then
      currNode => tslNodePreceding
      do while (.not. ldone)
        currObj => currNode%GetItem()
        currRecord => CastAsTimeSeriesRecordType(currObj)
        currTime = currRecord%tsrTime
        if (is_close(currTime, time1)) then
          ! Current node time = time1 so should be ldone
          ldone = .true.
        elseif (currTime < time1) then
          if (.not. associated(currNode%nextNode)) then
            ! -- try to read the next record
            if (.not. this%read_next_record()) then
              if (.not. extendToEndOfSimulation) then
                write (errmsg, 10) trim(this%Name), time0, time1
                call store_error(errmsg, terminate=.TRUE.)
              end if
            end if
          end if
          !
          currVal = currRecord%tsrValue
          lprocess = .false.
          if (associated(currNode%nextNode)) then
            nextNode => currNode%nextNode
            nextObj => nextNode%GetItem()
            nextRecord => CastAsTimeSeriesRecordType(nextObj)
            nextTime = nextRecord%tsrTime
            nextVal = nextRecord%tsrValue
            lprocess = .true.
          elseif (extendToEndOfSimulation) then
            ! -- Last time series value extends forever, so integrate the final value over all simulation time after the end of the series
            nextTime = time1
            nextVal = currVal
            lprocess = .true.
          end if
          !
          if (lprocess) then
            ! -- determine lower and upper limits of time span of interest
            !    within current interval
            if (currTime > time0 .or. is_close(currTime, time0)) then
              t0 = currTime
            else
              t0 = time0
            end if
            if (nextTime < time1 .or. is_close(nextTime, time1)) then
              t1 = nextTime
            else
              t1 = time1
            end if
            ! -- find area of rectangle or trapezoid delimited by t0 and t1
            t01 = t1 - t0
            select case (this%iMethod)
            case (STEPWISE)
              ! -- compute area of a rectangle
              value0 = currVal
              area = value0 * t01
            case (LINEAR, LINEAREND)
              ! -- compute area of a trapezoid
              timediff = nextTime - currTime
              ratio0 = (t0 - currTime) / timediff
              ratio1 = (t1 - currTime) / timediff
              valuediff = nextVal - currVal
              value0 = currVal + ratio0 * valuediff
              value1 = currVal + ratio1 * valuediff
              if (this%iMethod == LINEAR) then
                area = 0.5d0 * t01 * (value0 + value1)
              elseif (this%iMethod == LINEAREND) then
                area = DZERO
                value = value1
              end if
            end select
            ! -- add area to integrated value
            value = value + area
          end if
        end if
        !
        ! -- Are we done yet?
        if (t1 > time1) then
          ldone = .true.
        elseif (is_close(t1, time1)) then
          ldone = .true.
        else
          ! -- We are not done yet
          if (.not. associated(currNode%nextNode)) then
            ! -- Not done and no more data, so try to read the next record
            if (.not. this%read_next_record()) then
              write (errmsg, 10) trim(this%Name), time0, time1
              call store_error(errmsg, terminate=.TRUE.)
            end if
          elseif (associated(currNode%nextNode)) then
            currNode => currNode%nextNode
          end if
        end if
      end do
    end if
    !
    get_integrated_value = value
    if (this%autoDeallocate) then
      if (associated(tslNodePreceding)) then
        if (associated(tslNodePreceding%prevNode)) then
          call this%list%DeallocateBackward(tslNodePreceding%prevNode)
        end if
      end if
    end if
  end function get_integrated_value

  !> @brief Get average value
  !!
  !! Return a time-weighted average value for a specified time span.
  !! Units: (ts-value-unit)
  !<
  function get_average_value(this, time0, time1, extendToEndOfSimulation)
    ! -- return
    real(DP) :: get_average_value
    ! -- dummy
    class(TimeSeriesType), intent(inout) :: this
    real(DP), intent(in) :: time0
    real(DP), intent(in) :: time1
    logical, intent(in) :: extendToEndOfSimulation
    ! -- local
    real(DP) :: timediff, value, valueIntegrated
    !
    timediff = time1 - time0
    if (timediff > 0) then
      valueIntegrated = this%get_integrated_value(time0, time1, &
                                                  extendToEndOfSimulation)
      if (this%iMethod == LINEAREND) then
        value = valueIntegrated
      else
        value = valueIntegrated / timediff
      end if
    else
      ! -- time0 and time1 are the same
      value = this%get_value_at_time(time0, extendToEndOfSimulation)
    end if
    get_average_value = value
  end function get_average_value

  !> @brief Get latest preceding node
  !!
  !! Return pointer to ListNodeType object for the node representing the
  !! latest preceding time in the time series
  !<
  subroutine get_latest_preceding_node(this, time, tslNode)
    ! -- dummy
    class(TimeSeriesType), intent(inout) :: this
    real(DP), intent(in) :: time
    type(ListNodeType), pointer, intent(inout) :: tslNode
    ! -- local
    real(DP) :: time0
    type(ListNodeType), pointer :: currNode => null()
    type(ListNodeType), pointer :: tsNode0 => null()
    type(TimeSeriesRecordType), pointer :: tsr => null()
    type(TimeSeriesRecordType), pointer :: tsrec0 => null()
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
    ! -- If the next node is earlier than time of interest, advance along
    !    linked list until the next node is later than time of interest.
    do
      if (associated(currNode)) then
        if (associated(currNode%nextNode)) then
          obj => currNode%nextNode%GetItem()
          tsr => CastAsTimeSeriesRecordType(obj)
          if (tsr%tsrTime < time .or. is_close(tsr%tsrTime, time)) then
            currNode => currNode%nextNode
          else
            exit
          end if
        else
          ! -- read another record
          if (.not. this%read_next_record()) exit
        end if
      else
        exit
      end if
    end do
    !
    if (associated(currNode)) then
      !
      ! -- find earlier record
      tsNode0 => currNode
      obj => tsNode0%GetItem()
      tsrec0 => CastAsTimeSeriesRecordType(obj)
      time0 = tsrec0%tsrTime
      do while (time0 > time)
        if (associated(tsNode0%prevNode)) then
          tsNode0 => tsNode0%prevNode
          obj => tsNode0%GetItem()
          tsrec0 => CastAsTimeSeriesRecordType(obj)
          time0 = tsrec0%tsrTime
        else
          exit
        end if
      end do
    end if
    !
    if (time0 < time .or. is_close(time0, time)) tslNode => tsNode0
  end subroutine get_latest_preceding_node

  !> @brief Deallocate
  !<
  subroutine ts_da(this)
    ! -- dummy
    class(TimeSeriesType), intent(inout) :: this
    !
    if (associated(this%list)) then
      call this%list%Clear(.true.)
      deallocate (this%list)
    end if
  end subroutine ts_da

  !> @brief Add ts record
  !<
  subroutine AddTimeSeriesRecord(this, tsr)
    ! -- dummy
    class(TimeSeriesType) :: this
    type(TimeSeriesRecordType), pointer, intent(inout) :: tsr
    ! -- local
    class(*), pointer :: obj => null()
    !
    obj => tsr
    call this%list%Add(obj)
  end subroutine AddTimeSeriesRecord

  !> @brief Get current ts record
  !<
  function GetCurrentTimeSeriesRecord(this) result(res)
    ! -- dummy
    class(TimeSeriesType) :: this
    ! -- result
    type(TimeSeriesRecordType), pointer :: res
    ! -- local
    class(*), pointer :: obj => null()
    !
    obj => null()
    res => null()
    obj => this%list%GetItem()
    if (associated(obj)) then
      res => CastAsTimeSeriesRecordType(obj)
    end if
  end function GetCurrentTimeSeriesRecord

  !> @brief Get previous ts record
  !<
  function GetPreviousTimeSeriesRecord(this) result(res)
    ! -- dummy
    class(TimeSeriesType) :: this
    ! -- result
    type(TimeSeriesRecordType), pointer :: res
    ! -- local
    class(*), pointer :: obj => null()
    !
    obj => null()
    res => null()
    obj => this%list%GetPreviousItem()
    if (associated(obj)) then
      res => CastAsTimeSeriesRecordType(obj)
    end if
  end function GetPreviousTimeSeriesRecord

  !> @brief Get next ts record
  !<
  function GetNextTimeSeriesRecord(this) result(res)
    ! -- dummy
    class(TimeSeriesType) :: this
    ! -- result
    type(TimeSeriesRecordType), pointer :: res
    ! -- local
    class(*), pointer :: obj => null()
    !
    obj => null()
    res => null()
    obj => this%list%GetNextItem()
    if (associated(obj)) then
      res => CastAsTimeSeriesRecordType(obj)
    end if
  end function GetNextTimeSeriesRecord

  !> @brief Get ts record
  !<
  function GetTimeSeriesRecord(this, time, epsi) result(res)
    ! -- dummy
    class(TimeSeriesType) :: this
    double precision, intent(in) :: time
    double precision, intent(in) :: epsi
    ! -- result
    type(TimeSeriesRecordType), pointer :: res
    ! -- local
    type(TimeSeriesRecordType), pointer :: tsr
    !
    call this%list%Reset()
    res => null()
    do
      tsr => this%GetNextTimeSeriesRecord()
      if (associated(tsr)) then
        if (is_close(tsr%tsrTime, time)) then
          res => tsr
          exit
        end if
        if (tsr%tsrTime > time) exit
      else
        exit
      end if
    end do
  end function GetTimeSeriesRecord

  !> @brief Reset
  !<
  subroutine Reset(this)
    ! -- dummy
    class(TimeSeriesType) :: this
    !
    call this%list%Reset()
  end subroutine Reset

  !> @brief Insert a time series record
  !<
  subroutine InsertTsr(this, tsr)
    ! -- dummy
    class(TimeSeriesType), intent(inout) :: this
    type(TimeSeriesRecordType), pointer, intent(inout) :: tsr
    ! -- local
    double precision :: badtime, time, time0, time1
    type(TimeSeriesRecordType), pointer :: tsrEarlier, tsrLater
    type(ListNodeType), pointer :: nodeEarlier, nodeLater
    class(*), pointer :: obj => null()
    !
    badtime = -9.0d30
    time0 = badtime
    time1 = badtime
    time = tsr%tsrTime
    call this%get_surrounding_nodes(time, nodeEarlier, nodeLater)
    !
    if (associated(nodeEarlier)) then
      obj => nodeEarlier%GetItem()
      tsrEarlier => CastAsTimeSeriesRecordType(obj)
      if (associated(tsrEarlier)) then
        time0 = tsrEarlier%tsrTime
      end if
    end if
    !
    if (associated(nodeLater)) then
      obj => nodeLater%GetItem()
      tsrLater => CastAsTimeSeriesRecordType(obj)
      if (associated(tsrLater)) then
        time1 = tsrLater%tsrTime
      end if
    end if
    !
    if (time0 > badtime) then
      ! Time0 is valid
      if (time1 > badtime) then
        ! Both time0 and time1 are valid
        if (time > time0 .and. time < time1) then
          ! Insert record between two list nodes
          obj => tsr
          call this%list%InsertBefore(obj, nodeLater)
        else
          ! No need to insert a time series record, but if existing record
          ! for time of interest has NODATA as tsrValue, replace tsrValue
          if (time == time0 .and. tsrEarlier%tsrValue == DNODATA .and. &
              tsr%tsrValue /= DNODATA) then
            tsrEarlier%tsrValue = tsr%tsrValue
          elseif (time == time1 .and. tsrLater%tsrValue == DNODATA .and. &
                  tsr%tsrValue /= DNODATA) then
            tsrLater%tsrValue = tsr%tsrValue
          end if
        end if
      else
        ! Time0 is valid and time1 is invalid. Just add tsr to the list.
        call this%AddTimeSeriesRecord(tsr)
      end if
    else
      ! Time0 is invalid, so time1 must be for first node in list
      if (time1 > badtime) then
        ! Time 1 is valid
        if (time < time1) then
          ! Insert tsr at beginning of list
          obj => tsr
          call this%list%InsertBefore(obj, nodeLater)
        elseif (time == time1) then
          ! No need to insert a time series record, but if existing record
          ! for time of interest has NODATA as tsrValue, replace tsrValue
          if (tsrLater%tsrValue == DNODATA .and. tsr%tsrValue /= DNODATA) then
            tsrLater%tsrValue = tsr%tsrValue
          end if
        end if
      else
        ! Both time0 and time1 are invalid. Just add tsr to the list.
        call this%AddTimeSeriesRecord(tsr)
      end if
    end if
  end subroutine InsertTsr

  !> @brief Find latest time
  !<
  function FindLatestTime(this, readToEnd) result(endtime)
    ! -- dummy
    class(TimeSeriesType), intent(inout) :: this
    logical, intent(in), optional :: readToEnd
    ! -- local
    integer :: nrecords
    type(TimeSeriesRecordType), pointer :: tsr
    class(*), pointer :: obj => null()
    ! -- return
    double precision :: endtime
    !
    ! -- If the caller requested the very last time in the series (readToEnd is true), check that we have first read all records
    if (present(readToEnd)) then
      if (readToEnd) then
        do while (this%read_next_record())
        end do
      end if
    end if
    !
    nrecords = this%list%Count()
    obj => this%list%GetItem(nrecords)
    tsr => CastAsTimeSeriesRecordType(obj)
    endtime = tsr%tsrTime
  end function FindLatestTime

  !> @brief Clear the list of time series records
  !<
  subroutine Clear(this, destroy)
    ! -- dummy
    class(TimeSeriesType), intent(inout) :: this
    logical, optional, intent(in) :: destroy
    !
    call this%list%Clear(destroy)
  end subroutine Clear

! Type-bound procedures of TimeSeriesFileType

  !> @brief Count number of time series
  !<
  function Count(this)
    ! -- return
    integer(I4B) :: Count
    ! -- dummy
    class(TimeSeriesFileType) :: this
    !
    if (associated(this%timeSeries)) then
      Count = size(this%timeSeries)
    else
      Count = 0
    end if
  end function Count

  !> @brief Get time series
  !<
  function GetTimeSeries(this, indx) result(res)
    ! -- dummy
    class(TimeSeriesFileType) :: this
    integer(I4B), intent(in) :: indx
    ! -- return
    type(TimeSeriesType), pointer :: res
    !
    res => null()
    if (indx > 0 .and. indx <= this%nTimeSeries) then
      res => this%timeSeries(indx)
    end if
  end function GetTimeSeries

  !> @brief Open time-series tsfile file and read options and first record,
  !! which may contain data to define multiple time series.
  !<
  subroutine Initializetsfile(this, filename, iout, autoDeallocate)
    ! -- dummy
    class(TimeSeriesFileType), target, intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer(I4B), intent(in) :: iout
    logical, optional, intent(in) :: autoDeallocate
    ! -- local
    integer(I4B) :: iMethod, istatus, j, nwords
    integer(I4B) :: ierr, inunit
    logical :: autoDeallocateLocal = .true.
    logical :: continueread, found, endOfBlock
    logical :: methodWasSet
    real(DP) :: sfaclocal
    character(len=40) :: keyword, keyvalue
    character(len=:), allocatable :: line
    character(len=LENTIMESERIESNAME), allocatable, dimension(:) :: words
    !
    ! -- Initialize some variables
    if (present(autoDeallocate)) autoDeallocateLocal = autoDeallocate
    iMethod = UNDEFINED
    methodWasSet = .false.
    !
    ! -- Assign members
    this%iout = iout
    this%datafile = filename
    !
    ! -- Open the time-series tsfile input file
    this%inunit = GetUnit()
    inunit = this%inunit
    call openfile(inunit, 0, filename, 'TS6')
    !
    ! -- Initialize block parser
    call this%parser%Initialize(this%inunit, this%iout)
    !
    ! -- Read the ATTRIBUTES block and count time series
    continueread = .false.
    ierr = 0
    !
    ! -- get BEGIN line of ATTRIBUTES block
    call this%parser%GetBlock('ATTRIBUTES', found, ierr, &
                              supportOpenClose=.true.)
    if (ierr /= 0) then
      ! end of file
      errmsg = 'End-of-file encountered while searching for'// &
               ' ATTRIBUTES in time-series '// &
               'input file "'//trim(this%datafile)//'"'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    elseif (.not. found) then
      errmsg = 'ATTRIBUTES block not found in time-series '// &
               'tsfile input file "'//trim(this%datafile)//'"'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- parse ATTRIBUTES entries
    do
      ! -- read a line from input
      call this%parser%GetNextLine(endOfBlock)
      if (endOfBlock) exit
      !
      ! -- get the keyword
      call this%parser%GetStringCaps(keyword)
      !
      ! support either NAME or NAMES as equivalent keywords
      if (keyword == 'NAMES') keyword = 'NAME'
      !
      if (keyword /= 'NAME' .and. keyword /= 'METHODS' .and. &
          keyword /= 'SFACS') then
        ! -- get the word following the keyword (the key value)
        call this%parser%GetStringCaps(keyvalue)
      end if
      !
      select case (keyword)
      case ('NAME')
!        line = line(istart:linelen)
        call this%parser%GetRemainingLine(line)
        call ParseLine(line, nwords, words, this%parser%iuactive)
        this%nTimeSeries = nwords
        ! -- Allocate the timeSeries array and initialize each
        !    time series.
        allocate (this%timeSeries(this%nTimeSeries))
        do j = 1, this%nTimeSeries
          call this%timeSeries(j)%initialize_time_series(this, words(j), &
                                                         autoDeallocateLocal)
        end do
      case ('METHOD')
        methodWasSet = .true.
        if (this%nTimeSeries == 0) then
          errmsg = 'Error: NAME attribute not provided before METHOD in file: ' &
                   //trim(filename)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
        select case (keyvalue)
        case ('STEPWISE')
          iMethod = STEPWISE
        case ('LINEAR')
          iMethod = LINEAR
        case ('LINEAREND')
          iMethod = LINEAREND
        case default
          errmsg = 'Unknown interpolation method: "'//trim(keyvalue)//'"'
          call store_error(errmsg)
        end select
        do j = 1, this%nTimeSeries
          this%timeSeries(j)%iMethod = iMethod
        end do
      case ('METHODS')
        methodWasSet = .true.
        if (this%nTimeSeries == 0) then
          errmsg = 'Error: NAME attribute not provided before METHODS in file: ' &
                   //trim(filename)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
        call this%parser%GetRemainingLine(line)
        call ParseLine(line, nwords, words, this%parser%iuactive)
        if (nwords < this%nTimeSeries) then
          errmsg = 'METHODS attribute does not list a method for'// &
                   ' all time series.'
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
        do j = 1, this%nTimeSeries
          call upcase(words(j))
          select case (words(j))
          case ('STEPWISE')
            iMethod = STEPWISE
          case ('LINEAR')
            iMethod = LINEAR
          case ('LINEAREND')
            iMethod = LINEAREND
          case default
            errmsg = 'Unknown interpolation method: "'//trim(words(j))//'"'
            call store_error(errmsg)
          end select
          this%timeSeries(j)%iMethod = iMethod
        end do
      case ('SFAC')
        if (this%nTimeSeries == 0) then
          errmsg = 'NAME attribute not provided before SFAC in file: ' &
                   //trim(filename)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
        read (keyvalue, *, iostat=istatus) sfaclocal
        if (istatus /= 0) then
          errmsg = 'Error reading numeric value from: "'//trim(keyvalue)//'"'
          call store_error(errmsg)
        end if
        do j = 1, this%nTimeSeries
          this%timeSeries(j)%sfac = sfaclocal
        end do
      case ('SFACS')
        if (this%nTimeSeries == 0) then
          errmsg = 'NAME attribute not provided before SFACS in file: ' &
                   //trim(filename)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
        do j = 1, this%nTimeSeries
          sfaclocal = this%parser%GetDouble()
          this%timeSeries(j)%sfac = sfaclocal
        end do
      case ('AUTODEALLOCATE')
        do j = 1, this%nTimeSeries
          this%timeSeries(j)%autoDeallocate = (keyvalue == 'TRUE')
        end do
      case default
        errmsg = 'Unknown option found in ATTRIBUTES block: "'// &
                 trim(keyword)//'"'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      end select
    end do
    !
    ! -- Get TIMESERIES block
    call this%parser%GetBlock('TIMESERIES', found, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- Read the first line of time-series data
    if (.not. this%read_tsfile_line()) then
      errmsg = 'Error: No time-series data contained in file: '// &
               trim(this%datafile)
      call store_error(errmsg)
    end if
    !
    ! -- Ensure method was set
    if (.not. methodWasSet) then
      errmsg = 'Interpolation method was not set.  METHOD or METHODS &
      &must be specified in the ATTRIBUTES block for this time series file.'
      call store_error(errmsg)
    end if
    !
    ! -- Clean up and return
    if (allocated(words)) deallocate (words)
    !
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine Initializetsfile

  !> @brief Read time series file line
  logical function read_tsfile_line(this)
    ! -- dummy
    class(TimeSeriesFileType), intent(inout) :: this
    ! -- local
    real(DP) :: tsrTime, tsrValue
    integer(I4B) :: i
    logical :: endOfBlock
    type(TimeSeriesRecordType), pointer :: tsRecord => null()
    !
    read_tsfile_line = .false.
    !
    ! -- Get an arbitrary length, non-comment, non-blank line
    !    from the input file.
    call this%parser%GetNextLine(endOfBlock)
    !
    ! -- Check if we've reached the end of the TIMESERIES block
    if (endOfBlock) then
      return
    end if
    !
    ! -- Get the time
    tsrTime = this%parser%GetDouble()
    !
    ! -- Construct a new record and append a new node to each time series
    tsloop: do i = 1, this%nTimeSeries
      tsrValue = this%parser%GetDouble()
      if (tsrValue == DNODATA) cycle tsloop
      ! -- multiply value by sfac
      tsrValue = tsrValue * this%timeSeries(i)%sfac
      call ConstructTimeSeriesRecord(tsRecord, tsrTime, tsrValue)
      call AddTimeSeriesRecordToList(this%timeSeries(i)%list, tsRecord)
    end do tsloop
    read_tsfile_line = .true.
  end function read_tsfile_line

  !> @brief Deallocate memory
  !<
  subroutine tsf_da(this)
    ! -- dummy
    class(TimeSeriesFileType), intent(inout) :: this
    ! -- local
    integer :: i, n
    type(TimeSeriesType), pointer :: ts => null()
    !
    n = this%Count()
    do i = 1, n
      ts => this%GetTimeSeries(i)
      if (associated(ts)) then
        call ts%da()
!        deallocate(ts)
      end if
    end do
    !
    deallocate (this%timeSeries)
    deallocate (this%parser)
  end subroutine tsf_da

end module TimeSeriesModule
