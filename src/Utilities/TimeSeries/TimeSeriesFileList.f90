module TimeSeriesFileListModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH
  use ListModule, only: ListType
  use TimeSeriesModule, only: TimeSeriesFileType, &
                              ConstructTimeSeriesFile, &
                              GetTimeSeriesFileFromList, &
                              AddTimeSeriesFileToList

  implicit none

  private
  public :: TimeSeriesFileListType

  type :: TimeSeriesFileListType
    ! -- Public members
    integer(I4B), public :: numtsfiles = 0
    type(ListType), public :: tsfileList

  contains

    ! -- Public procedures
    procedure, public :: Add
    procedure, public :: Counttsfiles
    procedure, public :: CountTimeSeries
    procedure, public :: Gettsfile
    procedure, public :: Clear
    procedure, public :: da => tsfl_da
    procedure, public :: add_time_series_tsfile
  end type TimeSeriesFileListType

contains

  ! -- Public procedures

  subroutine Add(this, filename, iout, tsfile)
    implicit none
    ! -- dummy
    class(TimeSeriesFileListType), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer(I4B), intent(in) :: iout
    class(TimeSeriesFileType), pointer, intent(inout) :: tsfile
    ! -- local
    type(TimeSeriesFileType), pointer :: tsf
    !
    ! -- Construct and initialize a new time-series tsfile
    call ConstructTimeSeriesFile(tsf)
    tsfile => tsf
    call tsfile%Initializetsfile(filename, iout, .true.)
    !
    ! -- Add the time-series tsfile to the list
    call this%add_time_series_tsfile(tsfile)
  end subroutine Add

  subroutine Clear(this)
    implicit none
    ! -- dummy
    class(TimeSeriesFileListType), intent(inout) :: this
    !
    call this%tsfileList%Clear()
  end subroutine Clear

  function Counttsfiles(this)
    implicit none
    ! -- return
    integer(I4B) :: Counttsfiles
    ! -- dummy
    class(TimeSeriesFileListType) :: this
    !
    Counttsfiles = this%tsfileList%Count()
    !
  end function Counttsfiles

  function CountTimeSeries(this)
    implicit none
    ! -- return
    integer(I4B) :: CountTimeSeries
    ! -- dummy
    class(TimeSeriesFileListType) :: this
    ! -- local
    integer(I4B) :: i, numtsfiles
    type(TimeSeriesFileType), pointer :: tsfile
    !
    numtsfiles = this%Counttsfiles()
    CountTimeSeries = 0
    do i = 1, numtsfiles
      tsfile => this%Gettsfile(i)
      if (associated(tsfile)) then
        CountTimeSeries = CountTimeSeries + tsfile%Count()
      end if
    end do
  end function CountTimeSeries

  function Gettsfile(this, indx) result(res)
    implicit none
    ! -- dummy
    class(TimeSeriesFileListType) :: this
    integer(I4B), intent(in) :: indx
    ! -- return
    type(TimeSeriesFileType), pointer :: res
    !
    res => GetTimeSeriesFileFromList(this%tsfileList, indx)
  end function Gettsfile

  ! -- Private procedures

  subroutine add_time_series_tsfile(this, tsfile)
    implicit none
    ! -- dummy
    class(TimeSeriesFileListType), intent(inout) :: this
    class(TimeSeriesFileType), pointer, intent(inout) :: tsfile
    !
    call AddTimeSeriesFileToList(this%tsfileList, tsfile)
    this%numtsfiles = this%numtsfiles + 1
  end subroutine add_time_series_tsfile

  subroutine tsfl_da(this)
    ! -- dummy
    class(TimeSeriesFileListType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i, n
    type(TimeSeriesFileType), pointer :: tsf => null()
    !
    n = this%Counttsfiles()
    do i = 1, n
      tsf => this%Gettsfile(i)
      call tsf%da()
    end do
    !
    call this%tsfileList%Clear(.true.)
    !
  end subroutine tsfl_da

end module TimeSeriesFileListModule
