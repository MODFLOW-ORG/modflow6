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
    procedure :: destroy => vx_destroy
    ! private
    procedure, private :: create_virtual_fields
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

    ! - if both models are local, is_local = true
    ! - if both are remote, is_local = false
    ! - if only one of them is remote, is_local = true and only the
    ! remote nodem1/2 array will get its property is_local = false
    is_local = this%v_model1%is_local .or. this%v_model2%is_local
    call this%VirtualDataContainerType%vdc_create(name, exg_id, is_local)

    ! allocate fields
    call this%create_virtual_fields()

  end subroutine vx_create

  subroutine create_virtual_fields(this)
    class(VirtualExchangeType) :: this
    ! local
    logical(LGP) :: is_nodem1_local
    logical(LGP) :: is_nodem2_local

    ! exchanges can be hybrid with both local and remote fields
    is_nodem1_local = this%is_local .and. this%v_model1%is_local
    is_nodem2_local = this%is_local .and. this%v_model2%is_local

    allocate (this%nexg)
    call this%create_field(this%nexg%to_base(), 'NEXG', '')
    allocate (this%naux)
    call this%create_field(this%naux%to_base(), 'NAUX', '')
    allocate (this%ianglex)
    call this%create_field(this%ianglex%to_base(), 'IANGLEX', '')
    allocate (this%nodem1)
    call this%create_field(this%nodem1%to_base(), 'NODEM1', '', is_nodem1_local)
    allocate (this%nodem2)
    call this%create_field(this%nodem2%to_base(), 'NODEM2', '', is_nodem2_local)
    allocate (this%ihc)
    call this%create_field(this%ihc%to_base(), 'IHC', '')
    allocate (this%cl1)
    call this%create_field(this%cl1%to_base(), 'CL1', '')
    allocate (this%cl2)
    call this%create_field(this%cl2%to_base(), 'CL2', '')
    allocate (this%hwva)
    call this%create_field(this%hwva%to_base(), 'HWVA', '')
    allocate (this%auxvar)
    call this%create_field(this%auxvar%to_base(), 'AUXVAR', '')

  end subroutine create_virtual_fields

  subroutine vx_prepare_stage(this, stage)
    class(VirtualExchangeType) :: this
    integer(I4B) :: stage
    ! local
    integer(I4B) :: nexg, naux

    if (stage == STG_AFTER_EXG_DF) then

      call this%map(this%nexg%to_base(), (/STG_AFTER_EXG_DF/), MAP_ALL_TYPE)
      call this%map(this%naux%to_base(), (/STG_AFTER_EXG_DF/), MAP_ALL_TYPE)
      call this%map(this%ianglex%to_base(), (/STG_AFTER_EXG_DF/), MAP_ALL_TYPE)

    else if (stage == STG_BEFORE_DF) then

      nexg = this%nexg%get()
      naux = this%naux%get()
      call this%map(this%nodem1%to_base(), nexg, (/STG_BEFORE_DF/), MAP_ALL_TYPE)
      call this%map(this%nodem2%to_base(), nexg, (/STG_BEFORE_DF/), MAP_ALL_TYPE)
      call this%map(this%ihc%to_base(), nexg, (/STG_BEFORE_DF/), MAP_ALL_TYPE)
      call this%map(this%cl1%to_base(), nexg, (/STG_BEFORE_DF/), MAP_ALL_TYPE)
      call this%map(this%cl2%to_base(), nexg, (/STG_BEFORE_DF/), MAP_ALL_TYPE)
      call this%map(this%hwva%to_base(), nexg, (/STG_BEFORE_DF/), MAP_ALL_TYPE)
      call this%map(this%auxvar%to_base(), naux, nexg, &
                    (/STG_BEFORE_DF/), MAP_ALL_TYPE)

    end if

  end subroutine vx_prepare_stage

  subroutine vx_get_recv_items(this, stage, rank, virtual_items)
    class(VirtualExchangeType) :: this
    integer(I4B) :: stage
    integer(I4B) :: rank
    type(STLVecInt) :: virtual_items
    ! local
    integer(I4B) :: i, nodem1_idx, nodem2_idx
    class(*), pointer :: virtual_data_item

    nodem1_idx = -1
    nodem2_idx = -1
    do i = 1, this%virtual_data_list%Count()
      virtual_data_item => this%virtual_data_list%GetItem(i)
      if (associated(virtual_data_item, this%nodem1)) then
        nodem1_idx = i
      end if
      if (associated(virtual_data_item, this%nodem2)) then
        nodem2_idx = i
      end if
    end do

    if (.not. this%v_model1%is_local .and. .not. this%v_model2%is_local) then
      ! receive all using base
      call this%VirtualDataContainerType%get_recv_items(stage, rank, &
                                                        virtual_items)
    else if (.not. this%v_model2%is_local) then
      ! receive for model2
      if (this%nodem2%check_stage(stage)) then
        call virtual_items%push_back(nodem2_idx)
      end if
    else if (.not. this%v_model1%is_local) then
      ! receive for model1
      if (this%nodem1%check_stage(stage)) then
        call virtual_items%push_back(nodem1_idx)
      end if
    end if

  end subroutine vx_get_recv_items

  subroutine vx_get_send_items(this, stage, rank, virtual_items)
    class(VirtualExchangeType) :: this
    integer(I4B) :: stage
    integer(I4B) :: rank
    type(STLVecInt) :: virtual_items
    ! local
    integer(I4B) :: i, nodem1_idx, nodem2_idx
    class(*), pointer :: virtual_data_item

    nodem1_idx = -1
    nodem2_idx = -1
    do i = 1, this%virtual_data_list%Count()
      virtual_data_item => this%virtual_data_list%GetItem(i)
      if (associated(virtual_data_item, this%nodem1)) then
        nodem1_idx = i
      end if
      if (associated(virtual_data_item, this%nodem2)) then
        nodem2_idx = i
      end if
    end do

    if (this%v_model1%is_local .and. this%v_model2%is_local) then
      ! send all using base
      call this%VirtualDataContainerType%get_send_items(stage, rank, &
                                                        virtual_items)
    else if (this%v_model1%is_local) then
      ! send for model1
      if (this%nodem1%check_stage(stage)) then
        call virtual_items%push_back(nodem1_idx)
      end if
    else if (this%v_model2%is_local) then
      ! send for model2
      if (this%nodem2%check_stage(stage)) then
        call virtual_items%push_back(nodem2_idx)
      end if
    end if

  end subroutine vx_get_send_items

  subroutine vx_destroy(this)
    class(VirtualExchangeType) :: this

    call this%VirtualDataContainerType%destroy()
    call this%deallocate_data()

  end subroutine vx_destroy

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

  !> @brief Returs a virtual exchange with the specified id
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
