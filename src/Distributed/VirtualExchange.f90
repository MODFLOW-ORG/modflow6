module VirtualExchangeModule
  use VirtualBaseModule
  use VirtualDataContainerModule
  use VirtualModelModule, only: VirtualModelType, get_virtual_model
  use KindModule, only: I4B, LGP
  use ListModule, only: ListType
  use STLVecIntModule
  use ConstantsModule, only: LENEXCHANGENAME
  use SimStagesModule
  implicit none
  private

  public :: get_virtual_exchange
  public :: get_virtual_exchange_from_list
  private :: cast_as_virtual_exchange

  !> The Virtual Exchange is based on two Virtual Models
  !! and is therefore not always strictly local or remote.
  !! We have to consider three different cases:
  !!
  !! 1) both virtual models are local
  !!
  !!   RECV: In this case this virtual data container will have
  !!         no data items to receive from other processes.
  !!   SEND: Whenever it is called to send its virtual data items
  !!         to other processes, it simply sends everything.
  !!
  !! 2) one model is local, one model is remote
  !!
  !!   Consequently, there is another exchange which
  !!   has the reverse, we call this our _dual_ exchange.
  !!
  !!   RECV: The sender is our dual exchange, and we have all data
  !!         except its list of reduced model node numbers, either
  !!         this%nodem1 or this%nodem2. We receive the missing
  !!         array. Receiving from a sender that is not the dual
  !!         exchange cannot occur.
  !!
  !!   SEND: here we have to consider two cases
  !!      a) The receiver is our dual exchange, we return the favor
  !!         and send the list of model node numbers that is present
  !!         on this process, this
  !!         would be either this%nodem1 or this%nodem2
  !!      b) The receiver is not the dual exchange. And here we will
  !!         send everything.
  !!
  !! 3) both models are remote
  !!
  !!   RECV: we will receive everything. In case the source
  !!         exchange is fully local, i.e. type 1) above, we get
  !!         all the data at the first attempt. Otherwise, it will
  !!         take a second attempt before all the data is in.
  !!         (To allow for two attempts, the nodem1 and nodem2
  !!         arrays are registered to be synchronized at two
  !!         consecutive stages)
  !!
  !!   SEND: nothing to be sent.
  !!
  !! Exchange mover data follows the pattern described above for nodem1/m2
  !! except that when both models are remote, none of the mover data is
  !! will be synchronized.
  !!
  !! This behavior is different from the general VirtualDataContainer,
  !! so the get_send_items and get_recv_items subroutines are
  !! overridden accordingly.
  !! Additionally, for case 2) the container will have a mix of
  !< local and remote virtual data items.
  type, public, extends(VirtualDataContainerType) :: VirtualExchangeType
    class(VirtualModelType), pointer :: v_model1 => null()
    class(VirtualModelType), pointer :: v_model2 => null()
    ! scalars
    type(VirtualIntType), pointer :: nexg => null()
    type(VirtualIntType), pointer :: naux => null()
    type(VirtualIntType), pointer :: ianglex => null()
    ! arrays
    type(VirtualInt1dType), pointer :: nodem1 => null()
    type(VirtualInt1dType), pointer :: nodem2 => null()
    type(VirtualInt1dType), pointer :: ihc => null()
    type(VirtualDbl1dType), pointer :: cl1 => null()
    type(VirtualDbl1dType), pointer :: cl2 => null()
    type(VirtualDbl1dType), pointer :: hwva => null()
    type(VirtualDbl2dType), pointer :: auxvar => null()
  contains
    procedure :: create => vx_create
    procedure :: prepare_stage => vx_prepare_stage
    procedure :: get_send_items => vx_get_send_items
    procedure :: get_recv_items => vx_get_recv_items
    procedure :: has_mover => vx_has_mover
    procedure :: destroy => vx_destroy
    ! protected
    procedure :: add_vdi_for_stage
    ! private
    procedure, private :: init_virtual_data
    procedure, private :: allocate_data
    procedure, private :: deallocate_data
  end type VirtualExchangeType

contains

  !> @brief Create the virtual exchange base
  !<
  subroutine vx_create(this, name, exg_id, m1_id, m2_id)
    class(VirtualExchangeType) :: this
    character(len=*) :: name
    integer(I4B) :: exg_id
    integer(I4B) :: m1_id
    integer(I4B) :: m2_id
    ! local
    logical(LGP) :: is_local

    this%v_model1 => get_virtual_model(m1_id)
    this%v_model2 => get_virtual_model(m2_id)

    ! 1) both models local:  is_local = true
    ! 2) only one of them:   is_local = true
    ! 3) both models remote: is_local = false
    is_local = this%v_model1%is_local .or. this%v_model2%is_local
    call this%VirtualDataContainerType%vdc_create(name, exg_id, is_local)

    call this%allocate_data()
    call this%init_virtual_data()

  end subroutine vx_create

  subroutine init_virtual_data(this)
    class(VirtualExchangeType) :: this
    ! local
    logical(LGP) :: is_nodem1_local
    logical(LGP) :: is_nodem2_local

    ! exchanges can be hybrid with both local and remote
    ! fields, nodem1/2 array only local when corresponding
    ! model sits on the same process
    is_nodem1_local = this%v_model1%is_local
    is_nodem2_local = this%v_model2%is_local
    call this%set(this%nexg%base(), 'NEXG', '', MAP_ALL_TYPE)
    call this%set(this%naux%base(), 'NAUX', '', MAP_ALL_TYPE)
    call this%set(this%ianglex%base(), 'IANGLEX', '', MAP_ALL_TYPE)
    call this%set(this%nodem1%base(), 'NODEM1', '', &
                  MAP_ALL_TYPE, is_nodem1_local)
    call this%set(this%nodem2%base(), 'NODEM2', '', &
                  MAP_ALL_TYPE, is_nodem2_local)
    call this%set(this%ihc%base(), 'IHC', '', MAP_ALL_TYPE)
    call this%set(this%cl1%base(), 'CL1', '', MAP_ALL_TYPE)
    call this%set(this%cl2%base(), 'CL2', '', MAP_ALL_TYPE)
    call this%set(this%hwva%base(), 'HWVA', '', MAP_ALL_TYPE)
    call this%set(this%auxvar%base(), 'AUXVAR', '', MAP_ALL_TYPE)

  end subroutine init_virtual_data

  subroutine vx_prepare_stage(this, stage)
    class(VirtualExchangeType) :: this
    integer(I4B) :: stage
    ! local
    integer(I4B) :: nexg, naux

    if (stage == STG_AFT_EXG_DF) then

      call this%map(this%nexg%base(), (/STG_AFT_EXG_DF/))
      call this%map(this%naux%base(), (/STG_AFT_EXG_DF/))
      call this%map(this%ianglex%base(), (/STG_AFT_EXG_DF/))

    else if (stage == STG_AFT_CON_CR) then

      nexg = this%nexg%get()
      naux = this%naux%get()
      call this%map(this%nodem1%base(), nexg, (/STG_AFT_CON_CR, &
                                                STG_BFR_CON_DF/))
      call this%map(this%nodem2%base(), nexg, (/STG_AFT_CON_CR, &
                                                STG_BFR_CON_DF/))
      call this%map(this%ihc%base(), nexg, (/STG_AFT_CON_CR/))
      call this%map(this%cl1%base(), nexg, (/STG_AFT_CON_CR/))
      call this%map(this%cl2%base(), nexg, (/STG_AFT_CON_CR/))
      call this%map(this%hwva%base(), nexg, (/STG_AFT_CON_CR/))
      call this%map(this%auxvar%base(), naux, nexg, (/STG_AFT_CON_CR/))

    end if

  end subroutine vx_prepare_stage

  subroutine vx_get_recv_items(this, stg, rank, vi)
    class(VirtualExchangeType) :: this
    integer(I4B) :: stg
    integer(I4B) :: rank
    type(STLVecInt) :: vi

    if (this%is_local .and. rank == this%orig_rank) then
      ! treat the primary exchange case independently, we
      ! have all data available except for nodem1 or nodem2
      if (stg < STG_BFR_CON_DF) then
        if (this%nodem1%is_remote) then
          call this%add_vdi_for_stage(this%nodem1%base(), stg, vi)
        end if
        if (this%nodem2%is_remote) then
          call this%add_vdi_for_stage(this%nodem2%base(), stg, vi)
        end if
      end if
    else
      ! send/receive all
      call this%add_vdi_for_stage(this%nexg%base(), stg, vi)
      call this%add_vdi_for_stage(this%naux%base(), stg, vi)
      call this%add_vdi_for_stage(this%ianglex%base(), stg, vi)
      call this%add_vdi_for_stage(this%nodem1%base(), stg, vi)
      call this%add_vdi_for_stage(this%nodem2%base(), stg, vi)
      call this%add_vdi_for_stage(this%ihc%base(), stg, vi)
      call this%add_vdi_for_stage(this%cl1%base(), stg, vi)
      call this%add_vdi_for_stage(this%cl2%base(), stg, vi)
      call this%add_vdi_for_stage(this%hwva%base(), stg, vi)
      call this%add_vdi_for_stage(this%auxvar%base(), stg, vi)
    end if

  end subroutine vx_get_recv_items

  subroutine vx_get_send_items(this, stg, rank, vi)
    class(VirtualExchangeType) :: this
    integer(I4B) :: stg
    integer(I4B) :: rank
    type(STLVecInt) :: vi

    if (this%is_local .and. rank == this%orig_rank) then
      ! this is a primary exchange, all we need to send are
      ! the node numbers nodem1 or nodem2
      if (stg < STG_BFR_CON_DF) then
        if (.not. this%nodem1%is_remote) then
          call this%add_vdi_for_stage(this%nodem1%base(), stg, vi)
        end if
        if (.not. this%nodem2%is_remote) then
          call this%add_vdi_for_stage(this%nodem2%base(), stg, vi)
        end if
      end if
    else
      ! send/receive all
      call this%add_vdi_for_stage(this%nexg%base(), stg, vi)
      call this%add_vdi_for_stage(this%naux%base(), stg, vi)
      call this%add_vdi_for_stage(this%ianglex%base(), stg, vi)
      call this%add_vdi_for_stage(this%nodem1%base(), stg, vi)
      call this%add_vdi_for_stage(this%nodem2%base(), stg, vi)
      call this%add_vdi_for_stage(this%ihc%base(), stg, vi)
      call this%add_vdi_for_stage(this%cl1%base(), stg, vi)
      call this%add_vdi_for_stage(this%cl2%base(), stg, vi)
      call this%add_vdi_for_stage(this%hwva%base(), stg, vi)
      call this%add_vdi_for_stage(this%auxvar%base(), stg, vi)
    end if

  end subroutine vx_get_send_items

  !> @brief Convenience routine to add virtual data item to a list
  !< after checking the stage
  subroutine add_vdi_for_stage(this, vdata_item, stage, virtual_items)
    class(VirtualExchangeType) :: this
    class(VirtualDataType), pointer :: vdata_item
    integer(I4B) :: stage
    type(STLVecInt) :: virtual_items
    ! local
    class(*), pointer :: vdi
    integer(I4B) :: idx

    vdi => vdata_item
    idx = this%virtual_data_list%GetIndex(vdi)
    if (vdata_item%check_stage(stage)) then
      call virtual_items%push_back(idx)
    end if

  end subroutine add_vdi_for_stage

  !> @brief Checks if there is an active mover in the exchange
  !<
  function vx_has_mover(this) result(has_mover)
    class(VirtualExchangeType) :: this
    logical(LGP) :: has_mover

    has_mover = .false.

  end function vx_has_mover

  subroutine vx_destroy(this)
    class(VirtualExchangeType) :: this

    call this%VirtualDataContainerType%destroy()
    call this%deallocate_data()

  end subroutine vx_destroy

  subroutine allocate_data(this)
    class(VirtualExchangeType) :: this

    allocate (this%nexg)
    allocate (this%naux)
    allocate (this%ianglex)
    allocate (this%nodem1)
    allocate (this%nodem2)
    allocate (this%ihc)
    allocate (this%cl1)
    allocate (this%cl2)
    allocate (this%hwva)
    allocate (this%auxvar)

  end subroutine allocate_data

  subroutine deallocate_data(this)
    class(VirtualExchangeType) :: this

    deallocate (this%nexg)
    deallocate (this%naux)
    deallocate (this%ianglex)
    deallocate (this%nodem1)
    deallocate (this%nodem2)
    deallocate (this%ihc)
    deallocate (this%cl1)
    deallocate (this%cl2)
    deallocate (this%hwva)
    deallocate (this%auxvar)

  end subroutine deallocate_data

  !> @brief Returns a virtual exchange with the specified id
  !< from the global list
  function get_virtual_exchange(exg_id) result(virtual_exg)
    use VirtualDataListsModule, only: virtual_exchange_list
    integer(I4B) :: exg_id
    class(VirtualExchangeType), pointer :: virtual_exg
    ! local
    integer(I4B) :: i
    class(*), pointer :: ve

    virtual_exg => null()
    do i = 1, virtual_exchange_list%Count()
      ve => virtual_exchange_list%GetItem(i)
      select type (ve)
      class is (VirtualExchangeType)
        if (ve%id == exg_id) then
          virtual_exg => ve
          return
        end if
      end select
    end do

  end function get_virtual_exchange

  function get_virtual_exchange_from_list(list, idx) result(virtual_exg)
    type(ListType) :: list
    integer(I4B) :: idx
    class(VirtualExchangeType), pointer :: virtual_exg
    ! local
    class(*), pointer :: obj_ptr

    obj_ptr => list%GetItem(idx)
    virtual_exg => cast_as_virtual_exchange(obj_ptr)

  end function get_virtual_exchange_from_list

  function cast_as_virtual_exchange(obj_ptr) result(virtual_exg)
    class(*), pointer :: obj_ptr
    class(VirtualExchangeType), pointer :: virtual_exg

    virtual_exg => null()
    select type (obj_ptr)
    class is (VirtualExchangeType)
      virtual_exg => obj_ptr
    end select

  end function cast_as_virtual_exchange

end module VirtualExchangeModule
