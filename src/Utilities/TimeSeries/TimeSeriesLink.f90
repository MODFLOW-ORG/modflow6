module TimeSeriesLinkModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, LENBOUNDNAME, LENPACKAGENAME, &
                             LENTIMESERIESTEXT
  use InputOutputModule, only: UPCASE
  use ListModule, only: ListType
  use TimeSeriesModule, only: TimeSeriesType

  implicit none

  private
  public :: TimeSeriesLinkType, ConstructTimeSeriesLink, &
            GetTimeSeriesLinkFromList, AddTimeSeriesLinkToList
  private :: CastAsTimeSeriesLinkType

  type :: TimeSeriesLinkType
    ! -- Public members
    integer(I4B), public :: IRow = 0 ! row index (2nd dim) in bound or auxval array
    integer(I4B), public :: JCol = 0 ! column index (1st dim) in bound or auxval array
    integer(I4B), public :: Iprpak = 1
    ! BndElement can point to an element in either the bound or auxval
    ! array of BndType, or any other double precision variable or array
    ! element that contains a value that could be controlled by a time series.
    real(DP), pointer, public :: BndElement => null()
    real(DP), pointer, public :: RMultiplier => null()
    real(DP), public :: CellArea = DZERO
    character(len=LENPACKAGENAME), public :: PackageName = ''
    character(len=3), public :: AuxOrBnd = ''
    character(len=LENTIMESERIESTEXT), public :: Text = ''
    character(len=LENBOUNDNAME), public :: BndName = ''
    logical, public :: Active = .true.
    logical, public :: UseDefaultProc = .true.
    logical, public :: ConvertFlux = .false.
    type(TimeSeriesType), pointer, public :: TimeSeries => null()
  end type TimeSeriesLinkType

contains

  !> @brief Construct time series link
  !<
  subroutine ConstructTimeSeriesLink(newTsLink, timeSeries, pkgName, &
                                     auxOrBnd, bndElem, iRow, jCol, iprpak, &
                                     text)
    implicit none
    ! -- dummy
    type(TimeSeriesLinkType), pointer, intent(out) :: newTsLink
    type(TimeSeriesType), pointer, intent(in) :: timeSeries
    integer(I4B), intent(in) :: iRow, jCol
    character(len=*), intent(in) :: pkgName
    character(len=3), intent(in) :: auxOrBnd
    real(DP), pointer, intent(in) :: bndElem
    integer(I4B), intent(in) :: iprpak
    character(len=*), optional :: text
    ! -- local
    character(len=LENPACKAGENAME) :: pkgNameTemp
    !
    allocate (newTsLink)
    !
    ! Store package name as all caps
    pkgNameTemp = pkgName
    call upcase(pkgNameTemp)
    newTsLink%PackageName = pkgNameTemp
    newTsLink%AuxOrBnd = auxOrBnd
    newTsLink%timeSeries => timeSeries
    newTsLink%iRow = iRow
    newTsLink%jCol = jCol
    newTsLink%BndElement => bndElem
    newTsLink%Iprpak = iprpak
    !
    if (present(text)) then
      newTsLink%Text = text
    end if
  end subroutine ConstructTimeSeriesLink

  !> @brief Cast an unlimited polymorphic object as TimeSeriesLinkType
  !<
  function CastAsTimeSeriesLinkType(obj) result(res)
    implicit none
    ! -- dummy
    class(*), pointer, intent(inout) :: obj
    ! -- return
    type(TimeSeriesLinkType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    type is (TimeSeriesLinkType)
      res => obj
    class default
      continue
    end select
  end function CastAsTimeSeriesLinkType

  !> @brief Get time series link from a list
  !<
  function GetTimeSeriesLinkFromList(list, indx) result(tsLink)
    implicit none
    ! -- dummy
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: indx
    ! -- return
    type(TimeSeriesLinkType), pointer :: tsLink
    ! -- local
    class(*), pointer :: obj
    !
    tsLink => null()
    obj => list%GetItem(indx)
    tsLink => CastAsTimeSeriesLinkType(obj)
  end function GetTimeSeriesLinkFromList

  !> @brief Add time series link to a list
  !<
  subroutine AddTimeSeriesLinkToList(list, tslink)
    implicit none
    ! -- dummy
    type(ListType), intent(inout) :: list
    type(TimeSeriesLinkType), pointer, intent(inout) :: tslink
    ! -- local
    class(*), pointer :: obj
    !
    obj => tslink
    call list%Add(obj)
  end subroutine AddTimeSeriesLinkToList

end module TimeSeriesLinkModule
