module TimeArraySeriesLinkModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENPACKAGENAME, LENTIMESERIESTEXT
  use InputOutputModule, only: UPCASE
  use ListModule, only: ListType
  use TimeArraySeriesModule, only: TimeArraySeriesType

  implicit none

  private
  public :: TimeArraySeriesLinkType, ConstructTimeArraySeriesLink, &
            AddTimeArraySeriesLinkToList, GetTimeArraySeriesLinkFromList
  private :: CastAsTimeArraySeriesLinkType

  type :: TimeArraySeriesLinkType
    ! -- Public members
    character(len=LENPACKAGENAME), public :: PackageName = ''
    character(len=LENTIMESERIESTEXT), public :: Text = ''
    integer(I4B), public :: Iprpak = 1
    logical, public :: UseDefaultProc = .true.
    logical, public :: ConvertFlux = .false.
    integer(I4B), dimension(:), pointer, contiguous, public :: nodelist => null()
    ! BndArray can point to an array in either the bound or auxval
    ! array of BndType, or any other double precision variable or array
    ! element that contains a value that could be controlled by a time series.
    real(DP), dimension(:), pointer, public :: BndArray => null()
    real(DP), dimension(:), pointer, public :: RMultArray => null()
    type(TimeArraySeriesType), pointer, public :: TimeArraySeries => null()

  contains

    procedure, public :: da => tasl_da
  end type TimeArraySeriesLinkType

contains

  !> @brief Deallocate
  !<
  subroutine tasl_da(this)
    ! -- dummy
    class(TimeArraySeriesLinkType), intent(inout) :: this
    !
    this%nodelist => null()
    this%bndarray => null()
    this%rmultarray => null()
    this%TimeArraySeries => null()
  end subroutine tasl_da

  !> @brief Construct a time series of arrays that are linked
  !<
  subroutine ConstructTimeArraySeriesLink(newTasLink, timeArraySeries, &
                                          pkgName, bndArray, iprpak, text)
    ! -- dummy
    type(TimeArraySeriesLinkType), pointer, intent(out) :: newTasLink
    type(TimeArraySeriesType), pointer, intent(in) :: timeArraySeries
    character(len=*), intent(in) :: pkgName
    real(DP), dimension(:), pointer, intent(in) :: bndArray
    integer(I4B), intent(in) :: iprpak
    character(len=*), intent(in) :: text
    ! -- local
    character(len=LENPACKAGENAME) :: pkgNameTemp
    !
    allocate (newTasLink)
    ! Store package name as all caps
    pkgNameTemp = pkgName
    call UPCASE(pkgNameTemp)
    newTasLink%PackageName = pkgNameTemp
    newTasLink%timeArraySeries => timeArraySeries
    newTasLink%BndArray => bndArray
    newTasLink%Iprpak = iprpak
    newTasLink%Text = text
  end subroutine ConstructTimeArraySeriesLink

  !> @brief Cast an unlimited polymorphic object as TimeArraySeriesLinkType
  !<
  function CastAsTimeArraySeriesLinkType(obj) result(res)
    ! -- dummy
    class(*), pointer, intent(inout) :: obj
    ! -- return
    type(TimeArraySeriesLinkType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    type is (TimeArraySeriesLinkType)
      res => obj
    end select
  end function CastAsTimeArraySeriesLinkType

  !> @brief Add time-array series to list
  !<
  subroutine AddTimeArraySeriesLinkToList(list, tasLink)
    ! -- dummy
    type(ListType), intent(inout) :: list
    ! -- return
    type(TimeArraySeriesLinkType), pointer, intent(inout) :: tasLink
    ! -- local
    class(*), pointer :: obj
    !
    obj => tasLink
    call list%Add(obj)
  end subroutine AddTimeArraySeriesLinkToList

  !> @brief Get time-array series from a list and return
  !<
  function GetTimeArraySeriesLinkFromList(list, idx) result(res)
    ! -- dummy
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: idx
    ! -- return
    type(TimeArraySeriesLinkType), pointer :: res
    ! -- local
    class(*), pointer :: obj
    !
    obj => list%GetItem(idx)
    res => CastAsTimeArraySeriesLinkType(obj)
  end function GetTimeArraySeriesLinkFromList

end module TimeArraySeriesLinkModule

