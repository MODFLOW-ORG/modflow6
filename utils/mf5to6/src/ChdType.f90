module ChdModule

  use ConstantsModule, only: LENTIMESERIESNAME, MAXCHARLEN
  use ListModule, only: ListType
  use SimPHMFModule, only: store_error, ustop
  use TimeSeriesModule, only: TimeSeriesType

  implicit none

  private
  public :: ChdType, CastAsChdType, ConstructChdType, AddChdToList, &
            GetChdFromList, CellInChdList

  type :: ChdType
    integer                                     :: jcol
    integer                                     :: irow
    integer                                     :: klay
    double precision                            :: head
    double precision, allocatable, dimension(:) :: auxvars
    logical,          allocatable, dimension(:) :: ActiveByStressPeriod
    logical,          allocatable, dimension(:) :: ActiveInMf2005
    character(len=MAXCHARLEN)                   :: tsfilename = ''
    character(len=LENTIMESERIESNAME)            :: AlternateTsName = ''
    character(len=LENTIMESERIESNAME)            :: OriginalTsName = ''
    type(TimeSeriesType), pointer               :: timeSeries => null()
  contains
    procedure :: CopyTo
  end type ChdType

  interface ConstructChdType
    module procedure ConstructChdType1, ConstructChdType2
  end interface

contains

  subroutine CopyTo(this, newchd)
    ! dummy
    class(ChdType) :: this
    type(ChdType), pointer :: newchd
    ! local
    integer :: i, ndim
    !
    if (.not. associated(newchd)) then
      call store_error('programmer error in call to ChdType%CopyTo.')
      call ustop()
    endif
    !
    newchd%head = this%head
    newchd%irow = this%irow
    newchd%jcol = this%jcol
    newchd%klay = this%klay
    !
    if (allocated(this%auxvars)) then
      ndim = size(this%auxvars)
      allocate(newchd%auxvars(ndim))
      do i=1,ndim
        newchd%auxvars(i) = this%auxvars(i)
      enddo
    endif
    !
    return
  end subroutine CopyTo

  subroutine ConstructChdType1(newChd)
    implicit none
    type(ChdType), pointer, intent(inout) :: newChd
    !
    allocate(newChd)
    allocate(newChd%timeSeries)
    return
  end subroutine ConstructChdType1

  subroutine ConstructChdType2(newChd, nVars)
    implicit none
    type(ChdType), pointer, intent(inout) :: newChd
    integer, intent(in) :: nVars
    !
    allocate(newChd)
    allocate(newChd%timeSeries)
    allocate(newChd%auxvars(nVars))
    return
  end subroutine ConstructChdType2

  function CastAsChdType(obj) result(res)
    ! Cast an unlimited polymorphic object as class(ChdType)
    implicit none
    ! dummy
    class(*), pointer :: obj
    class(ChdType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    class is (ChdType)
      res => obj
    end select
    !
    return
  end function CastAsChdType

  subroutine AddChdToList(chdList, chd)
    ! dummy
    type(ListType) :: chdList
    type(ChdType), pointer :: chd
    ! local
    class(*), pointer :: obj => null()
    !
    obj => chd
    call chdList%Add(obj)
    !
    return
  end subroutine AddChdToList

  function GetChdFromList(chdList, idx) result (res)
    ! dummy
    type(ListType), intent(inout) :: chdList
    integer, intent(in) :: idx
    type(ChdType), pointer :: res
    ! local
    class(*), pointer :: obj => null()
    !
    obj => chdList%GetItem(idx)
    res => CastAsChdType(obj)
    !
    return
  end function GetChdFromList

  logical function CellInChdList(chdList, klay, irow, jcol)
    ! dummy
    type(ListType), intent(inout) :: chdList
    integer, intent(in) :: klay, irow, jcol
    ! local
    integer :: i, n
    type(ChdType), pointer :: chd
    !
    CellInChdList = .false.
    n = chdList%Count()
    do i=1,n
      chd => GetChdFromList(chdList, i)
      if (chd%klay == klay .and. chd%irow == irow .and. chd%jcol == jcol) then
        CellInChdList = .true.
        exit
      endif
    enddo
    !
    return
  end function CellInChdList

end module ChdModule
