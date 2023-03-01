module MapperModule
  use KindModule, only: I4B, LGP
  use ConstantsModule, only: LENVARNAME, LENMEMPATH
  use MemoryHelperModule, only: create_mem_path
  use IndexMapModule
  use VirtualModelModule, only: VirtualModelType, get_virtual_model
  use VirtualExchangeModule, only: VirtualExchangeType, get_virtual_exchange
  use InterfaceMapModule
  use DistVariableModule
  use MappedMemoryModule
  use ListModule
  implicit none
  private

  public :: MapperType

  type :: MapperType
    type(ListType) :: mapped_data_list
  contains
    procedure :: init
    procedure :: add_exchange_vars
    procedure :: add_interface_vars
    procedure :: scatter
    procedure :: destroy

    procedure, private :: add_dist_vars
    procedure, private :: map_model_data
    procedure, private :: map_exg_data
    procedure, private :: map_data
    procedure, private :: map_data_full
  end type MapperType

contains

  subroutine init(this)
    class(MapperType) :: this

  end subroutine init

  !> @brief Add virtual exchange variables
  !<
  subroutine add_exchange_vars(this)
    use SimStagesModule
    use ListsModule, only: baseconnectionlist
    use SpatialModelConnectionModule, only: SpatialModelConnectionType, &
                                            get_smc_from_list
    use VirtualExchangeModule, only: VirtualExchangeType, get_virtual_exchange
    class(MapperType) :: this
    ! local
    integer(I4B) :: iconn
    class(SpatialModelConnectionType), pointer :: conn
    class(VirtualExchangeType), pointer :: virt_exg
    character(len=LENMEMPATH) :: virt_mem_path

    do iconn = 1, baseconnectionlist%Count()
      conn => get_smc_from_list(baseconnectionlist, iconn)
      virt_exg => get_virtual_exchange(conn%prim_exchange%id)
      if (.not. virt_exg%v_model1%is_local) then
        virt_mem_path = virt_exg%get_vrt_mem_path('NODEM1', '')
        call this%map_data_full(0, 'NODEM1', conn%prim_exchange%memoryPath, &
                                'NODEM1', virt_mem_path, (/STG_BEFORE_DF/))
      end if
      if (.not. virt_exg%v_model2%is_local) then
        virt_mem_path = virt_exg%get_vrt_mem_path('NODEM2', '')
        call this%map_data_full(0, 'NODEM2', conn%prim_exchange%memoryPath, &
                                'NODEM2', virt_mem_path, (/STG_BEFORE_DF/))
      end if
    end do

  end subroutine add_exchange_vars

  !> @brief Add distributed interface variables as memory mapped items
  !<
  subroutine add_interface_vars(this)
    use ListsModule, only: baseconnectionlist
    use SpatialModelConnectionModule, only: SpatialModelConnectionType, &
                                            get_smc_from_list
    class(MapperType) :: this
    ! local
    integer(I4B) :: iconn
    class(SpatialModelConnectionType), pointer :: conn

    do iconn = 1, baseconnectionlist%Count()
      conn => get_smc_from_list(baseconnectionlist, iconn)
      ! add the variables for this interface model to our mapper
      call this%add_dist_vars(conn%owner%idsoln, &
                              conn%iface_dist_vars, &
                              conn%interface_map)
    end do

  end subroutine add_interface_vars

  subroutine add_dist_vars(this, sol_id, var_list, interface_map)
    class(MapperType) :: this
    integer(I4B) :: sol_id
    type(ListType) :: var_list
    type(InterfaceMapType) :: interface_map
    ! local
    integer(I4B) :: i, m, e
    type(DistVarType), pointer :: distvar
    type(IndexMapType), pointer :: idx_map

    ! loop over variables
    do i = 1, var_list%Count()
      distvar => GetDistVarFromList(var_list, i)
      if (distvar%map_type == SYNC_NODES .or. &
          distvar%map_type == SYNC_CONNECTIONS) then
        ! map data for all models in this interface
        do m = 1, interface_map%nr_models

          ! pick the right index map: connection based or node based
          if (distvar%map_type == SYNC_NODES) then
            idx_map => interface_map%node_map(m)
          else if (distvar%map_type == SYNC_CONNECTIONS) then
            idx_map => interface_map%connection_map(m)
          end if

          ! and map ...
          call this%map_model_data(sol_id, &
                                   distvar%comp_name, &
                                   distvar%subcomp_name, &
                                   distvar%var_name, &
                                   interface_map%model_ids(m), &
                                   idx_map, &
                                   distvar%sync_stages)
        end do
      else if (distvar%map_type == SYNC_EXCHANGES) then
        ! map data from the exchanges to the interface
        do e = 1, interface_map%nr_exchanges
          call this%map_exg_data(sol_id, &
                                 distvar%comp_name, &
                                 distvar%subcomp_name, &
                                 distvar%var_name, &
                                 interface_map%exchange_ids(e), &
                                 distvar%exg_var_name, &
                                 interface_map%exchange_map(e), &
                                 distvar%sync_stages)
        end do
      end if
    end do

  end subroutine add_dist_vars

  !> @brief Map data from model memory to a target memory entry,
  !! with the specified map. The source and target items have
  !< the same name and (optionally) subcomponent name.
  subroutine map_model_data(this, controller_id, tgt_model_name, &
                            tgt_subcomp_name, tgt_var_name, src_model_id, &
                            index_map, stages)
    use SimModule, only: ustop
    use MemoryManagerModule, only: get_from_memorylist
    class(MapperType) :: this
    integer(I4B) :: controller_id !< e.g. the numerical solution where synchr. is controlled
    character(len=*), intent(in) :: tgt_model_name
    character(len=*), intent(in) :: tgt_subcomp_name
    character(len=*), intent(in) :: tgt_var_name
    integer(I4B), intent(in) :: src_model_id
    type(IndexMapType), intent(in) :: index_map
    integer(I4B), dimension(:), intent(in) :: stages !< array with 1 or multiple stages for synchronization
    ! local
    character(len=LENVARNAME) :: src_var_name
    character(len=LENMEMPATH) :: src_mem_path, tgt_mem_path
    class(VirtualModelType), pointer :: v_model

    v_model => get_virtual_model(src_model_id)

    if (len_trim(tgt_subcomp_name) > 0) then
      tgt_mem_path = create_mem_path(tgt_model_name, tgt_subcomp_name)
    else
      tgt_mem_path = create_mem_path(tgt_model_name)
    end if

    src_var_name = tgt_var_name
    src_mem_path = v_model%get_vrt_mem_path(src_var_name, tgt_subcomp_name)
    call this%map_data(controller_id, &
                       tgt_var_name, tgt_mem_path, index_map%tgt_idx, &
                       src_var_name, src_mem_path, index_map%src_idx, &
                       null(), stages)

  end subroutine map_model_data

  !> @brief Map memory from a Exchange to the specified memory entry,
  !< using the index map
  subroutine map_exg_data(this, controller_id, tgt_model_name, &
                          tgt_subcomp_name, tgt_var_name, src_exg_id, &
                          src_var_name, index_map_sgn, stages)
    use SimModule, only: ustop
    use MemoryManagerModule, only: get_from_memorylist
    class(MapperType) :: this
    integer(I4B) :: controller_id !< e.g. the numerical solution where synchr. is controlled
    character(len=*), intent(in) :: tgt_model_name
    character(len=*), intent(in) :: tgt_subcomp_name
    character(len=*), intent(in) :: tgt_var_name
    integer(I4B), intent(in) :: src_exg_id
    character(len=*), intent(in) :: src_var_name
    type(IndexMapSgnType), intent(in) :: index_map_sgn
    integer(I4B), dimension(:), intent(in) :: stages !< array with 1 or multiple stages for synchronization
    ! local
    character(len=LENMEMPATH) :: src_mem_path, tgt_mem_path
    class(VirtualExchangeType), pointer :: v_exchange

    v_exchange => get_virtual_exchange(src_exg_id)

    if (len_trim(tgt_subcomp_name) > 0) then
      tgt_mem_path = create_mem_path(tgt_model_name, tgt_subcomp_name)
    else
      tgt_mem_path = create_mem_path(tgt_model_name)
    end if

    src_mem_path = v_exchange%get_vrt_mem_path(src_var_name, '')
    call this%map_data(controller_id, &
                       tgt_var_name, tgt_mem_path, index_map_sgn%tgt_idx, &
                       src_var_name, src_mem_path, index_map_sgn%src_idx, &
                       index_map_sgn%sign, stages)

  end subroutine map_exg_data

  !> @brief Full copy between two variables in memory
  subroutine map_data_full(this, controller_id, tgt_name, tgt_path, &
                           src_name, src_path, stages)
    class(MapperType) :: this
    integer(I4B) :: controller_id
    character(len=*), intent(in) :: tgt_name
    character(len=*), intent(in) :: tgt_path
    character(len=*), intent(in) :: src_name
    character(len=*), intent(in) :: src_path
    integer(I4B), dimension(:), intent(in) :: stages

    call this%map_data(controller_id, tgt_name, tgt_path, null(), &
                       src_name, src_path, null(), &
                       null(), stages)

  end subroutine map_data_full

  !> @brief Generic mapping between two variables in memory, using
  !< an optional sign conversion
  subroutine map_data(this, controller_id, tgt_name, tgt_path, tgt_idx, &
                      src_name, src_path, src_idx, sign_array, stages)
    class(MapperType) :: this
    integer(I4B) :: controller_id
    character(len=*), intent(in) :: tgt_name
    character(len=*), intent(in) :: tgt_path
    integer(I4B), dimension(:), pointer :: tgt_idx
    character(len=*), intent(in) :: src_name
    character(len=*), intent(in) :: src_path
    integer(I4B), dimension(:), pointer :: src_idx
    integer(I4B), dimension(:), pointer :: sign_array
    integer(I4B), dimension(:), intent(in) :: stages
    ! local
    integer(I4B) :: istage, i
    type(MappedMemoryType), pointer :: mapped_data
    class(*), pointer :: obj

    ! loop and set stage bits
    istage = 0
    do i = 1, size(stages)
      istage = ibset(istage, stages(i))
    end do

    ! create MappedVariable and add to list
    allocate (mapped_data)
    mapped_data%controller_id = controller_id
    mapped_data%sync_stage = istage
    mapped_data%src_name = src_name
    mapped_data%src_path = src_path
    mapped_data%src => null()
    mapped_data%tgt_name = tgt_name
    mapped_data%tgt_path = tgt_path
    mapped_data%tgt => null()
    mapped_data%copy_all = .not. associated(src_idx)
    mapped_data%src_idx => src_idx
    mapped_data%tgt_idx => tgt_idx
    mapped_data%sign => sign_array
    obj => mapped_data
    call this%mapped_data_list%Add(obj)

  end subroutine map_data

  !> @brief Scatter the mapped memory, typically into
  !< the memory space of the interface models
  subroutine scatter(this, controller_id, stage)
    class(MapperType) :: this
    integer(I4B) :: controller_id
    integer(I4B), intent(in) :: stage
    ! local
    integer(I4B) :: i
    class(*), pointer :: obj
    class(MappedMemoryType), pointer :: mapped_data

    ! sync all variables (src => tgt) for a given stage
    do i = 1, this%mapped_data_list%Count()
      obj => this%mapped_data_list%GetItem(i)
      mapped_data => CastAsMappedData(obj)
      if (controller_id > 0 .and. &
          mapped_data%controller_id /= controller_id) cycle
      if (.not. check_stage(mapped_data%sync_stage, stage)) cycle

      ! copy data
      call mapped_data%sync()
    end do

  end subroutine scatter

  function check_stage(var_stage, current_stage) result(is_sync)
    integer(I4B) :: var_stage
    integer(I4B) :: current_stage
    logical(LGP) :: is_sync

    is_sync = iand(var_stage, ibset(0, current_stage)) == ibset(0, current_stage)

  end function check_stage

  subroutine destroy(this)
    class(MapperType) :: this

    call this%mapped_data_list%Clear(destroy=.true.)

  end subroutine destroy

end module MapperModule
