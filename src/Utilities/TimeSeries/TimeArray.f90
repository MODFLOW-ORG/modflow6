module TimeArrayModule

  use KindModule, only: DP, I4B
  use ListModule, only: ListType
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error

  implicit none
  private
  public :: TimeArrayType, ConstructTimeArray, &
            AddTimeArrayToList, CastAsTimeArrayType, &
            GetTimeArrayFromList

  type :: TimeArrayType
    ! -- Public members
    real(DP), public :: taTime
    real(DP), dimension(:), pointer, contiguous, public :: taArray => null()

  contains

    ! -- Public procedures
    ! -- When gfortran adds support for finalization, the
    !    following declaration could be: final :: finalize
    procedure, public :: da => ta_da
  end type TimeArrayType

contains

  !> @brief Construct time array
  !!
  !! Allocate and assign members of a new TimeArrayType object. Allocate space
  !! for the array so that this subroutine can be called repeatedly with the
  !! same array (but with different contents).
  !<
  subroutine ConstructTimeArray(newTa, modelname)
    ! -- modules
    use ConstantsModule, only: LENMEMPATH
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path
    ! -- dummy
    type(TimeArrayType), pointer, intent(out) :: newTa
    character(len=*), intent(in) :: modelname
    ! -- local
    integer(I4B), dimension(:), contiguous, &
      pointer :: mshape
    character(len=LENMEMPATH) :: mempath
    integer(I4B) :: isize
    !
    ! -- initialize
    nullify (mshape)
    !
    ! -- create mempath
    mempath = create_mem_path(component=modelname, subcomponent='DIS')
    !
    ! -- set mshape pointer
    call mem_setptr(mshape, 'MSHAPE', mempath)
    !
    ! Get dimensions for supported discretization type
    if (size(mshape) == 2) then
      isize = mshape(2)
    else if (size(mshape) == 3) then
      isize = mshape(2) * mshape(3)
    else
      errmsg = 'Time array series is not supported for discretization type'
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    allocate (newTa)
    allocate (newTa%taArray(isize))
  end subroutine ConstructTimeArray

  !> @brief Cast an unlimited polymorphic object as TimeArrayType
  !<
  function CastAsTimeArrayType(obj) result(res)
    ! -- dummy
    class(*), pointer, intent(inout) :: obj
    ! -- return
    type(TimeArrayType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    type is (TimeArrayType)
      res => obj
    end select
  end function CastAsTimeArrayType

  !> @brief Add a time array to a to list
  !<
  subroutine AddTimeArrayToList(list, timearray)
    ! -- dummy
    type(ListType), intent(inout) :: list
    type(TimeArrayType), pointer, intent(inout) :: timearray
    ! -- local
    class(*), pointer :: obj
    !
    obj => timearray
    call list%Add(obj)
  end subroutine AddTimeArrayToList

  !> @brief Retrieve a time array from a list
  !<
  function GetTimeArrayFromList(list, indx) result(res)
    ! -- dummy
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: indx
    ! -- return
    type(TimeArrayType), pointer :: res
    ! -- local
    class(*), pointer :: obj
    !
    obj => list%GetItem(indx)
    res => CastAsTimeArrayType(obj)
  end function GetTimeArrayFromList

  !> @brief Deallocate memory
  !<
  subroutine ta_da(this)
    ! -- dummy
    class(TimeArrayType) :: this
    !
    deallocate (this%taArray)
    this%taArray => null()
  end subroutine ta_da

end module TimeArrayModule
