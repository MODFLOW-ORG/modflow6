module TimeSeriesRecordModule

  use KindModule, only: DP, I4B
  use ListModule, only: ListType

  private
  public :: TimeSeriesRecordType, ConstructTimeSeriesRecord, &
            CastAsTimeSeriesRecordType, AddTimeSeriesRecordToList

  type :: TimeSeriesRecordType
    ! -- Public members
    real(DP), public :: tsrTime
    real(DP), public :: tsrValue
  end type TimeSeriesRecordType

contains

  !> @brief Allocate and assign members of a new TimeSeriesRecordType object
  !<
  subroutine ConstructTimeSeriesRecord(newTsRecord, time, value)
    implicit none
    ! -- dummy
    type(TimeSeriesRecordType), pointer, intent(out) :: newTsRecord
    real(DP), intent(in) :: time, value
    !
    allocate (newTsRecord)
    newTsRecord%tsrTime = time
    newTsRecord%tsrValue = value
  end subroutine ConstructTimeSeriesRecord

  !> @brief Cast an unlimited polymorphic object as TimeSeriesRecordType
  !<
  function CastAsTimeSeriesRecordType(obj) result(res)
    implicit none
    ! -- dummy
    class(*), pointer, intent(inout) :: obj
    ! -- return
    type(TimeSeriesRecordType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    type is (TimeSeriesRecordType)
      res => obj
    end select
  end function CastAsTimeSeriesRecordType

  !> @brief Add time series record to list
  !<
  subroutine AddTimeSeriesRecordToList(list, tsrecord)
    implicit none
    ! -- dummy
    type(ListType), intent(inout) :: list
    type(TimeSeriesRecordType), pointer, intent(inout) :: tsrecord
    ! -- local
    class(*), pointer :: obj => null()
    !
    obj => tsrecord
    call list%Add(obj)
  end subroutine AddTimeSeriesRecordToList

end module TimeSeriesRecordModule
