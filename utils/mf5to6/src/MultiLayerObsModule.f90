module MultiLayerObs

  use ConstantsModule, only: DONE, MAXCHARLEN
  use ConstantsPHMFModule, only: LENOBSNAMENEW
  use MathUtilModule, only: is_close
  use ListModule, only: ListType
  use SimPHMFModule, only: store_error, ustop

  implicit none
  private
  public :: LayerObsType, MLObsType, ConstructLayerObs, AddLayerObsToList, &
            GetLayerObsFromList, ConstructMLObs, AddMLObsToList, &
            GetMLObsFromList

  type :: LayerObsType
    character(len=LENOBSNAMENEW) :: lobsname = ''
    integer :: layer = 0
    double precision :: weight = 0.0d0
  end type

  type :: MLObsType
    character(len=LENOBSNAMENEW) :: mlobsname = ''
    type(ListType) :: LayerObsList
    logical :: summ = .false.
  contains
    procedure, public :: CheckWeightSum
  end type

contains

  ! Non-type-bound procedures

  function CastAsLayerObsType(obj) result(res)
    class(*), pointer, intent(inout) :: obj
    type(LayerObsType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    type is (LayerObsType)
      res => obj
    end select
    return
  end function CastAsLayerObsType

  subroutine ConstructLayerObs(newLayerObs, layobname, layer, weight)
    ! dummy
    type(LayerObsType), pointer :: newLayerObs
    character(len=*), intent(in) :: layobname
    integer, intent(in) :: layer
    double precision, intent(in) :: weight
    !
    allocate (newLayerObs)
    newLayerObs%lobsname = layobname
    newLayerObs%layer = layer
    newLayerObs%weight = weight
    !
    return
  end subroutine ConstructLayerObs

  subroutine AddLayerObsToList(list, layerobs)
    ! dummy
    type(ListType), intent(inout) :: list
    type(LayerObsType), pointer, intent(inout) :: layerobs
    ! local
    class(*), pointer :: obj => null()
    !
    obj => layerobs
    call list%Add(obj)
    !
    return
  end subroutine AddLayerObsToList

  function GetLayerObsFromList(list, indx) result(res)
    ! dummy
    type(ListType), intent(inout) :: list
    integer, intent(in) :: indx
    type(LayerObsType), pointer :: res
    ! local
    class(*), pointer :: obj => null()
    !
    obj => list%GetItem(indx)
    select type (obj)
    type is (LayerObsType)
      res => obj
    end select
    !
    return
  end function GetLayerObsFromList

  subroutine ConstructMLObs(newMLObs, obsname)
    ! dummy
    type(MLObsType), pointer, intent(inout) :: newMLObs
    character(len=LENOBSNAMENEW), intent(in) :: obsname
    !
    allocate (newMLObs)
    newMLObs%mlobsname = obsname
    !
    return
  end subroutine ConstructMLObs

  subroutine AddMLObsToList(list, mlobs)
    ! dummy
    type(ListType), intent(inout) :: list
    type(MLObsType), pointer, intent(inout) :: mlobs
    ! local
    class(*), pointer :: obj => null()
    !
    obj => mlobs
    call list%Add(obj)
    !
    return
  end subroutine AddMLObsToList

  function GetMLObsFromList(list, indx) result(res)
    ! dummy
    type(ListType), intent(inout) :: list
    integer, intent(in) :: indx
    type(MLObsType), pointer :: res
    ! local
    class(*), pointer :: obj => null()
    !
    obj => list%GetItem(indx)
    select type (obj)
    type is (MLObsType)
      res => obj
    end select
    !
    return
  end function GetMLObsFromList

  ! Type-bound procedures of MLObsType

  subroutine CheckWeightSum(this)
    ! dummy
    class(MLObsType) :: this
    ! local
    double precision :: weightsum
    integer :: i, nlayers
    type(LayerObsType), pointer :: layobs => null()
    character(len=MAXCHARLEN) :: ermsg
    ! formats
10  format('Weights of layer observations do not sum to 1.0 for', &
           ' multilayer observation: ', a)
    !
    if (this%summ) return
    !
    weightsum = 0.0d0
    nlayers = this%LayerObsList%Count()
    do i = 1, nlayers
      layobs => GetLayerObsFromList(this%LayerObsList, i)
      weightsum = weightsum + layobs%weight
    end do
    !
    if (.not. is_close(weightsum, DONE)) then
      write (ermsg, 10) trim(this%mlobsname)
      call store_error(ermsg)
      call ustop()
    end if
    !
    return
  end subroutine

end module MultiLayerObs
