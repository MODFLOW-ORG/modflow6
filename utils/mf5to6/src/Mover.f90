module MoverModule

  use ConstantsModule, only: DONE, DZERO, LENMODELNAME, LENPACKAGENAME, &
                             MAXCHARLEN
  use ListModule, only: ListType
  use SimPHMFModule, only: store_error, ustop
  use utl7module, only: UPCASE

  private
  public :: MoverType, ConstructWaterMover, CastAsMoverType, &
            AddMoverToList, GetMoverFromList

  type :: MoverType
    character(len=LENPACKAGENAME) :: ProvPkgName = ''
    character(len=LENPACKAGENAME) :: RecPkgName = ''
    character(len=LENMODELNAME) :: ProvModelName = ''
    character(len=LENMODELNAME) :: RecModelName = ''
    character(len=9) :: MvrType = ''
    character(len=4) :: ProvPkgType = ''
    character(len=4) :: RecPkgType = ''
    integer :: IdProvider = 0
    integer :: IdReceiver = 0
    integer :: IgridProvider = 0
    integer :: IgridReceiver = 0
    double precision :: Valu = DZERO
  end type MoverType

contains

  subroutine ConstructWaterMover(newMover, mvrType, provModelName, &
                                 recModelName, provPkgName, &
                                 recPkgName, idProvider, idReceiver, &
                                 igridProvider, igridReceiver, provPkgType, &
                                 recPkgType, valu)
    implicit none
    ! dummy
    type(MoverType), pointer, intent(inout)   :: newMover
    character(len=*), intent(in)              :: mvrType
    character(len=*), intent(in) :: provPkgName, recPkgName
    character(len=*), intent(in) :: provModelName, recModelName
    integer, intent(in) :: idProvider
    integer, intent(in) :: idReceiver
    integer, intent(in) :: igridProvider
    integer, intent(in) :: igridReceiver
    character(len=*), intent(in) :: provPkgType, recPkgType
    double precision, optional, intent(in):: valu
    ! local
    character(len=9) :: mvrTypeCaps
    character(len=MAXCHARLEN) :: ermsg
    !
    allocate(newMover)
    ! Provider data
    newMover%ProvModelName = provModelName
    newMover%IgridProvider = igridProvider
    newMover%ProvPkgName = provPkgName
    newMover%ProvPkgType = provPkgType
    newMover%IdProvider = idProvider
    ! Receiver data
    newMover%RecModelName = recModelName
    newMover%IgridReceiver = igridReceiver
    newMover%RecPkgName = recPkgName
    newMover%RecPkgType = recPkgType
    newMover%IdReceiver = idReceiver
    !
    mvrTypeCaps = mvrType
    call UPCASE(mvrTypeCaps)
    select case (mvrTypeCaps)
    case ('FACTOR', 'EXCESS', 'THRESHOLD', 'UPTO')
      newMover%MvrType = mvrTypeCaps
    case default
      ermsg = 'Invalid mvrType in ConstructWaterMover: ' // trim(mvrType)
      call store_error(ermsg)
      call ustop()
    end select
    !
    if (present(valu)) then
      newMover%Valu = valu
    endif
    !
    return
  end subroutine ConstructWaterMover

  function CastAsMoverType(obj) result (res)
    implicit none
    class(*), pointer, intent(inout) :: obj
    type(MoverType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    type is (MoverType)
      res => obj
    end select
    return
  end function CastAsMoverType

  subroutine AddMoverToList(moverList, mover)
    implicit none
    ! dummy
    type(ListType) :: moverList
    type(MoverType), pointer :: mover
    ! local
    class(*), pointer :: obj => null()
    !
    obj => mover
    call moverList%Add(obj)
    !
    return
  end subroutine AddMoverToList

  function GetMoverFromList(moverList, idx) result (res)
    implicit none
    ! dummy
    type(ListType) :: moverList
    integer, intent(in) :: idx
    type(MoverType), pointer :: res
    ! local
    class(*), pointer :: obj => null()
    !
    obj => moverList%GetItem(idx)
    res => CastAsMoverType(obj)
    !
    return
  end function GetMoverFromList

end module MoverModule
