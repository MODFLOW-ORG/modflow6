module BaseExchangeModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENEXCHANGENAME, LENMEMPATH
  use ListModule, only: ListType
  use BaseModelModule, only: BaseModelType

  implicit none

  private
  public :: BaseExchangeType, AddBaseExchangeToList, GetBaseExchangeFromList
  private :: CastAsBaseExchangeClass

  type, abstract :: BaseExchangeType
    character(len=LENEXCHANGENAME) :: name !< the name of this exchange
    character(len=LENMEMPATH) :: memoryPath !< the location in the memory manager where the variables are stored
    character(len=LENMEMPATH) :: input_mempath
    integer(I4B) :: id

  contains

    procedure(exg_df), deferred :: exg_df
    procedure(exg_ar), deferred :: exg_ar
    procedure :: exg_rp
    procedure :: exg_dt
    procedure :: exg_ot
    procedure :: exg_fp
    procedure :: exg_da
    procedure :: connects_model
  end type BaseExchangeType

  abstract interface

    subroutine exg_df(this)
      import BaseExchangeType
      class(BaseExchangeType) :: this
    end subroutine

    subroutine exg_ar(this)
      import BaseExchangeType
      class(BaseExchangeType) :: this
    end subroutine

  end interface

contains

  !> @brief Read and prepare
  !<
  subroutine exg_rp(this)
    ! -- modules
    use TdisModule, only: readnewdata
    ! -- dummy
    class(BaseExchangeType) :: this
    !
    ! -- Check with TDIS on whether or not it is time to RP
    if (.not. readnewdata) return
    !
    ! -- Nothing to do for RP
  end subroutine exg_rp

  !> @brief Calculate time step length
  !<
  subroutine exg_dt(this)
    ! -- dummy
    class(BaseExchangeType) :: this
    !
    ! -- Nothing to do for TU
  end subroutine exg_dt

  !> @brief Run output routines
  !<
  subroutine exg_ot(this)
    ! -- dummy
    class(BaseExchangeType) :: this
  end subroutine exg_ot

  !> @brief Final processing
  !<
  subroutine exg_fp(this)
    ! -- dummy
    class(BaseExchangeType) :: this
  end subroutine exg_fp

  !> @brief Deallocate memory
  !<
  subroutine exg_da(this)
    ! -- dummy
    class(BaseExchangeType) :: this
  end subroutine exg_da

  !> @brief Should return true when the exchange should be added to the
  !! solution where the model resides
  !<
  function connects_model(this, model) result(is_connected)
    ! -- dummy
    class(BaseExchangeType) :: this !< the instance of the exchange
    class(BaseModelType), pointer, intent(in) :: model !< the model to which the exchange might hold a connection
    ! -- return
    logical(LGP) :: is_connected !< true, when connected
    !
    is_connected = .false.
  end function

  !> @brief Cast the object passed in as BaseExchangeType and return it
  !<
  function CastAsBaseExchangeClass(obj) result(res)
    ! -- dummy
    class(*), pointer, intent(inout) :: obj
    ! -- return
    class(BaseExchangeType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    class is (BaseExchangeType)
      res => obj
    end select
  end function CastAsBaseExchangeClass

  !> @brief Add the exchange object (BaseExchangeType) to a list
  !<
  subroutine AddBaseExchangeToList(list, exchange)
    ! -- dummy
    type(ListType), intent(inout) :: list
    class(BaseExchangeType), pointer, intent(inout) :: exchange
    ! -- local
    class(*), pointer :: obj
    !
    obj => exchange
    call list%Add(obj)
  end subroutine AddBaseExchangeToList

  !> @brief Retrieve a specific BaseExchangeType object from a list
  !<
  function GetBaseExchangeFromList(list, idx) result(res)
    ! -- dummy
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: idx
    class(BaseExchangeType), pointer :: res
    ! -- local
    class(*), pointer :: obj
    !
    obj => list%GetItem(idx)
    res => CastAsBaseExchangeClass(obj)
  end function GetBaseExchangeFromList

end module BaseExchangeModule
