module VirtualTspExchangeModule
  use KindModule, only: I4B, LGP
  use SimStagesModule
  use VirtualBaseModule
  use VirtualDataListsModule, only: virtual_exchange_list
  use VirtualDataContainerModule, only: VDC_GWTEXG_TYPE, VDC_GWEEXG_TYPE
  use VirtualExchangeModule
  use STLVecIntModule
  implicit none
  private

  public :: add_virtual_tsp_exchange

  !> GWE and GWT work fully analogously, so we can do
  !! with only one virtual exchange for both, at least for now
  !!
  !< Performance-TODO: why synchronize movers when the exchange is not primary?!
  type, public, extends(VirtualExchangeType) :: VirtualTspExchangeType
    type(VirtualIntType), pointer :: inmvt => null()
    type(VirtualDbl1dType), pointer :: gwfsimvals => null()
    type(VirtualIntType), pointer :: mvt_maxmvt => null()
    type(VirtualDbl1dType), pointer :: mvt_quantity_m1 => null()
    type(VirtualDbl1dType), pointer :: mvt_quantity_m2 => null()
    logical(LGP) :: quantities_mapped
    ! private
    logical(LGP), private :: has_mvt !< backing field for function
  contains
    procedure :: create => vtx_create
    procedure :: destroy => vtx_destroy
    procedure :: prepare_stage => vtx_prepare_stage
    procedure :: get_send_items => vtx_get_send_items
    procedure :: get_recv_items => vtx_get_recv_items
    procedure :: has_mover => vtx_has_mover
    ! private
    procedure, private :: init_virtual_data
    procedure, private :: allocate_data
    procedure, private :: deallocate_data
  end type VirtualTspExchangeType

contains

!> @brief Add a virtual GWT-GWT or GWE-GWE exchange to the simulation
!<
  subroutine add_virtual_tsp_exchange(name, exchange_id, m1_id, m2_id, qtype)
    character(len=*) :: name
    integer(I4B) :: exchange_id
    integer(I4B) :: m1_id !< id model 1
    integer(I4B) :: m2_id !< id model 2
    character(len=*) :: qtype !< quantity type (for GWE and GWT)
    ! local
    class(VirtualTspExchangeType), pointer :: v_exg
    class(*), pointer :: obj_ptr

    allocate (v_exg)

    call v_exg%create(name, exchange_id, m1_id, m2_id)
    if (qtype == "concentration") then
      v_exg%container_type = VDC_GWTEXG_TYPE
    else if (qtype == "temperature") then
      v_exg%container_type = VDC_GWEEXG_TYPE
    end if

    obj_ptr => v_exg
    call virtual_exchange_list%Add(obj_ptr)
  end subroutine add_virtual_tsp_exchange

!> @brief Create a virtual GWT-GWT exchange
!<
  subroutine vtx_create(this, name, exg_id, m1_id, m2_id)
    class(VirtualTspExchangeType) :: this
    character(len=*) :: name
    integer(I4B) :: exg_id
    integer(I4B) :: m1_id
    integer(I4B) :: m2_id

    ! create base
    call this%VirtualExchangeType%create(name, exg_id, m1_id, m2_id)

    call this%allocate_data()
    call this%init_virtual_data()

    this%quantities_mapped = .false.
    this%has_mvt = .false.

  end subroutine vtx_create

  subroutine init_virtual_data(this)
    class(VirtualTspExchangeType) :: this

    call this%set(this%inmvt%base(), 'INMVT', '', MAP_ALL_TYPE)
    call this%set(this%gwfsimvals%base(), 'GWFSIMVALS', '', MAP_ALL_TYPE)
    call this%set(this%mvt_maxmvt%base(), 'MAXMVT', 'MVT', MAP_ALL_TYPE)
    call this%set(this%mvt_quantity_m1%base(), 'QUANTITY_M1', 'MVT', &
                  MAP_ALL_TYPE, this%v_model1%is_local)
    call this%set(this%mvt_quantity_m2%base(), 'QUANTITY_M2', 'MVT', &
                  MAP_ALL_TYPE, this%v_model2%is_local)

  end subroutine init_virtual_data

  subroutine vtx_prepare_stage(this, stage)
    class(VirtualTspExchangeType) :: this
    integer(I4B) :: stage
    ! local
    integer(I4B) :: nexg, nmax

    ! prepare base exchange data items
    call this%VirtualExchangeType%prepare_stage(stage)

    if (stage == STG_AFT_EXG_DF) then

      ! always synchronize mover flag
      call this%map(this%inmvt%base(), (/STG_AFT_EXG_DF/))

    else if (stage == STG_AFT_CON_CR) then

      ! at this point we know:
      if (this%inmvt%get() > 0) then
        this%has_mvt = .true.
      end if

    else if (stage == STG_BFR_CON_AR) then

      nexg = this%nexg%get()
      call this%map(this%gwfsimvals%base(), nexg, (/STG_BFR_EXG_AD/))

      ! only when MVT is locally active (i.e. primary exchange)
      if (this%has_mvt .and. this%is_local) then
        call this%map(this%mvt_maxmvt%base(), (/STG_BFR_CON_AR/))
      end if

    else if (stage == STG_BFR_EXG_AD) then

      ! only when MVT is locally active
      nmax = 0
      if (this%has_mvt .and. this%is_local) nmax = this%mvt_maxmvt%get()

      ! only map the arrays once, after the first read prepare is done
      if (.not. this%quantities_mapped) then
        if (nmax > 0) then
          call this%map(this%mvt_quantity_m1%base(), nmax, (/STG_BFR_EXG_FC/))
          call this%map(this%mvt_quantity_m2%base(), nmax, (/STG_BFR_EXG_FC/))
        else
          call this%map(this%mvt_quantity_m1%base(), 0, (/STG_NEVER/))
          call this%map(this%mvt_quantity_m2%base(), 0, (/STG_NEVER/))
        end if
        this%quantities_mapped = .true.
      end if
    end if

  end subroutine vtx_prepare_stage

  subroutine vtx_get_recv_items(this, stg, rank, vi)
    class(VirtualTspExchangeType) :: this
    integer(I4B) :: stg !< stage
    integer(I4B) :: rank !< rank of remote process
    type(STLVecInt) :: vi !< virtual data items

    ! get base items to receive
    call this%VirtualExchangeType%get_recv_items(stg, rank, vi)

    ! add more MVT items that follow nodem1/nodem2 pattern,
    ! see comments in VirtualExchange for more details
    if (this%is_local .and. rank == this%orig_rank) then
      if (this%mvt_quantity_m1%is_remote) then
        ! only receive for model1
        call this%add_vdi_for_stage(this%mvt_quantity_m1%base(), stg, vi)
      end if
      if (this%mvt_quantity_m2%is_remote) then
        ! only receive for model2
        call this%add_vdi_for_stage(this%mvt_quantity_m2%base(), stg, vi)
      end if
    end if

  end subroutine vtx_get_recv_items

  subroutine vtx_get_send_items(this, stg, rank, vi)
    class(VirtualTspExchangeType) :: this
    integer(I4B) :: stg !< stage
    integer(I4B) :: rank !< rank of remote process
    type(STLVecInt) :: vi !< virtual data items

    ! get base items to send
    call this%VirtualExchangeType%get_send_items(stg, rank, vi)

    ! add more MVT items that follow nodem1/nodem2 pattern,
    ! see comments in VirtualExchange for more details
    if (this%is_local .and. rank == this%orig_rank) then
      if (.not. this%mvt_quantity_m1%is_remote) then
        ! only send for model1
        call this%add_vdi_for_stage(this%mvt_quantity_m1%base(), stg, vi)
      end if
      if (.not. this%mvt_quantity_m2%is_remote) then
        ! only send for model2
        call this%add_vdi_for_stage(this%mvt_quantity_m2%base(), stg, vi)
      end if
    end if

  end subroutine vtx_get_send_items

  !> @brief Override
  !<
  function vtx_has_mover(this) result(has_mover)
    class(VirtualTspExchangeType) :: this
    logical(LGP) :: has_mover

    has_mover = this%has_mvt

  end function vtx_has_mover

  subroutine vtx_destroy(this)
    class(VirtualTspExchangeType) :: this

    call this%VirtualExchangeType%destroy()
    call this%deallocate_data()

  end subroutine vtx_destroy

  subroutine allocate_data(this)
    class(VirtualTspExchangeType) :: this

    allocate (this%inmvt)
    allocate (this%gwfsimvals)
    allocate (this%mvt_maxmvt)
    allocate (this%mvt_quantity_m1)
    allocate (this%mvt_quantity_m2)

  end subroutine allocate_data

  subroutine deallocate_data(this)
    class(VirtualTspExchangeType) :: this

    deallocate (this%inmvt)
    deallocate (this%gwfsimvals)
    deallocate (this%mvt_maxmvt)
    deallocate (this%mvt_quantity_m1)
    deallocate (this%mvt_quantity_m2)

  end subroutine deallocate_data

end module VirtualTspExchangeModule
