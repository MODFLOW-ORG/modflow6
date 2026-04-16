module MapperModule
  use KindModule, only: I4B, LGP
  use ConstantsModule, only: LENVARNAME, LENMEMPATH
  use MemoryHelperModule, only: create_mem_path
  use IndexMapModule
  use VirtualBaseModule, only: VirtualDataType, MAP_NODE_TYPE, MAP_CONN_TYPE
  use VirtualDataContainerModule, only: VDC_GWFEXG_TYPE, VDC_GWTEXG_TYPE
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
    class(VirtualExchangeType), pointer :: vx
    character(len=LENMEMPATH) :: virt_mem_path, local_mem_path

    do iconn = 1, baseconnectionlist%Count()
      conn => get_smc_from_list(baseconnectionlist, iconn)
      vx => get_virtual_exchange(conn%prim_exchange%id)
      if (.not. vx%v_model1%is_local) then
        virt_mem_path = vx%get_vrt_mem_path('NODEM1', '')
        call this%map_data_full(0, 'NODEM1', conn%prim_exchange%memoryPath, &
                                'NODEM1', virt_mem_path, (/STG_BFR_CON_DF/))

        ! these are only present when there is a MVR:
        if (vx%has_mover() .and. vx%container_type == VDC_GWFEXG_TYPE) then
          local_mem_path = create_mem_path(vx%name, 'MVR')
          virt_mem_path = vx%get_vrt_mem_path('QPACTUAL_M1', 'MVR')
          call this%map_data_full(conn%owner%idsoln, 'QPACTUAL_M1', &
                                  local_mem_path, 'QPACTUAL_M1', &
                                  virt_mem_path, (/STG_BFR_EXG_FC/))
          virt_mem_path = vx%get_vrt_mem_path('QAVAILABLE_M1', 'MVR')
          call this%map_data_full(conn%owner%idsoln, 'QAVAILABLE_M1', &
                                  local_mem_path, 'QAVAILABLE_M1', &
                                  virt_mem_path, (/STG_BFR_EXG_FC/))
          virt_mem_path = vx%get_vrt_mem_path('ID_MAPPED_M1', 'MVR')
          call this%map_data_full(conn%owner%idsoln, 'ID_MAPPED_M1', &
                                  local_mem_path, 'ID_MAPPED_M1', &
                                  virt_mem_path, (/STG_AFT_CON_RP/))
        end if
        ! the same for MVT
        if (vx%has_mover() .and. vx%container_type == VDC_GWTEXG_TYPE) then
          local_mem_path = create_mem_path(vx%name, 'MVT')
          virt_mem_path = vx%get_vrt_mem_path('QUANTITY_M1', 'MVT')
          call this%map_data_full(conn%owner%idsoln, 'QUANTITY_M1', &
                                  local_mem_path, 'QUANTITY_M1', &
                                  virt_mem_path, (/STG_BFR_EXG_FC/))
        end if
      end if
      if (.not. vx%v_model2%is_local) then
        virt_mem_path = vx%get_vrt_mem_path('NODEM2', '')
        call this%map_data_full(0, 'NODEM2', conn%prim_exchange%memoryPath, &
                                'NODEM2', virt_mem_path, (/STG_BFR_CON_DF/))

        ! these are only present when there is a mover:
        if (vx%has_mover() .and. vx%container_type == VDC_GWFEXG_TYPE) then
          local_mem_path = create_mem_path(vx%name, 'MVR')
          virt_mem_path = vx%get_vrt_mem_path('QPACTUAL_M2', 'MVR')
          call this%map_data_full(conn%owner%idsoln, 'QPACTUAL_M2', &
                                  local_mem_path, 'QPACTUAL_M2', &
                                  virt_mem_path, (/STG_BFR_EXG_FC/))
          virt_mem_path = vx%get_vrt_mem_path('QAVAILABLE_M2', 'MVR')
          call this%map_data_full(conn%owner%idsoln, 'QAVAILABLE_M2', &
                                  local_mem_path, 'QAVAILABLE_M2', &
                                  virt_mem_path, (/STG_BFR_EXG_FC/))
          virt_mem_path = vx%get_vrt_mem_path('ID_MAPPED_M2', 'MVR')
          call this%map_data_full(conn%owner%idsoln, 'ID_MAPPED_M2', &
                                  local_mem_path, 'ID_MAPPED_M2', &
                                  virt_mem_path, (/STG_AFT_CON_RP/))
        end if
        ! the same for MVT
        if (vx%has_mover() .and. vx%container_type == VDC_GWTEXG_TYPE) then
          local_mem_path = create_mem_path(vx%name, 'MVT')
          virt_mem_path = vx%get_vrt_mem_path('QUANTITY_M2', 'MVT')
          call this%map_data_full(conn%owner%idsoln, 'QUANTITY_M2', &
                                  local_mem_path, 'QUANTITY_M2', &
                                  virt_mem_path, (/STG_BFR_EXG_FC/))
        end if
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

  subroutine add_dist_vars(this, sol_id, var_list, iface_map)
    class(MapperType) :: this
    integer(I4B) :: sol_id
    type(ListType) :: var_list
    type(InterfaceMapType), pointer :: iface_map
    ! local
    integer(I4B) :: i, m, e
    type(DistVarType), pointer :: dist_var

    ! loop over variables
    do i = 1, var_list%Count()
      dist_var => GetDistVarFromList(var_list, i)
      if (dist_var%map_type == SYNC_NDS .or. & ! models
          dist_var%map_type == SYNC_CON) then
        do m = 1, iface_map%nr_models
          call this%map_model_data(sol_id, iface_map, m, dist_var)
        end do
      else if (dist_var%map_type == SYNC_EXG) then ! exchanges
        do e = 1, iface_map%nr_exchanges
          call this%map_exg_data(sol_id, iface_map, e, dist_var)
        end do
      end if
    end do

  end subroutine add_dist_vars

  !> @brief Map data from model memory to a target memory entry,
  !! with the specified map. The source and target items have
  !< the same name and (optionally) subcomponent name.
  !call this%map_model_data(sol_id, interface_map%model_ids(m), &
  !dist_var, idx_map)
  subroutine map_model_data(this, sol_id, iface_map, model_idx, dist_var)
    use SimModule, only: ustop
    class(MapperType) :: this !< this mapper instance
    integer(I4B) :: sol_id !< the numerical solution where synchr. is controlled
    type(InterfaceMapType), pointer :: iface_map !< the full interface map
    integer(I4B) :: model_idx !< the model index (not id) in the interface map
    type(DistVarType), pointer :: dist_var !< the distributed variable to map
    ! local
    character(len=LENVARNAME) :: src_var_name
    character(len=LENMEMPATH) :: src_mem_path, tgt_mem_path
    class(VirtualModelType), pointer :: v_model
    type(IndexMapType), pointer :: idx_map
    integer(I4B), dimension(:), pointer, contiguous :: lookup_table
    class(VirtualDataType), pointer :: vd

    v_model => get_virtual_model(iface_map%model_ids(model_idx))
    vd => v_model%get_virtual_data(dist_var%var_name, dist_var%subcomp_name)

    ! pick the right index map: connection based or node based,
    ! and reduced data items require a lookup table
    lookup_table => null()
    if (dist_var%map_type == SYNC_NDS) then
      idx_map => iface_map%node_maps(model_idx)
      if (vd%is_reduced) then
        lookup_table => v_model%element_luts(MAP_NODE_TYPE)%remote_to_virtual
      end if
    else if (dist_var%map_type == SYNC_CON) then
      idx_map => iface_map%conn_maps(model_idx)
      if (vd%is_reduced) then
        lookup_table => v_model%element_luts(MAP_CONN_TYPE)%remote_to_virtual
      end if
    else
      write (*, *) "Unknown map type for distributed variable ", dist_var%var_name
      call ustop()
    end if

    if (len_trim(dist_var%subcomp_name) > 0) then
      tgt_mem_path = create_mem_path(dist_var%comp_name, dist_var%subcomp_name)
    else
      tgt_mem_path = create_mem_path(dist_var%comp_name)
    end if

    src_var_name = dist_var%var_name
    src_mem_path = v_model%get_vrt_mem_path(src_var_name, dist_var%subcomp_name)
    call this%map_data(sol_id, &
                       src_var_name, tgt_mem_path, idx_map%tgt_idx, &
                       src_var_name, src_mem_path, idx_map%src_idx, &
                       null(), lookup_table, dist_var%sync_stages)

  end subroutine map_model_data

  !> @brief Map memory from a Exchange to the specified memory entry,
  !< using the index map
  subroutine map_exg_data(this, sol_id, iface_map, exg_idx, dist_var)
    class(MapperType) :: this
    integer(I4B) :: sol_id !< the numerical solution where synchr. is controlled
    type(InterfaceMapType), pointer :: iface_map !< the full interface map
    integer(I4B), intent(in) :: exg_idx !< the index (not id) for the exchange
    type(DistVarType), pointer :: dist_var !< the distributed variable to map
    ! local
    character(len=LENMEMPATH) :: src_mem_path, tgt_mem_path
    class(VirtualExchangeType), pointer :: v_exchange
    type(IndexMapSgnType), pointer :: idx_map

    v_exchange => get_virtual_exchange(iface_map%exchange_ids(exg_idx))

    idx_map => iface_map%exchange_maps(exg_idx)

    if (len_trim(dist_var%subcomp_name) > 0) then
      tgt_mem_path = create_mem_path(dist_var%comp_name, dist_var%subcomp_name)
    else
      tgt_mem_path = create_mem_path(dist_var%comp_name)
    end if

    src_mem_path = v_exchange%get_vrt_mem_path(dist_var%exg_var_name, '')
    call this%map_data(sol_id, &
                       dist_var%var_name, tgt_mem_path, idx_map%tgt_idx, &
                       dist_var%exg_var_name, src_mem_path, idx_map%src_idx, &
                       idx_map%sign, null(), dist_var%sync_stages)

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
                       null(), null(), stages)

  end subroutine map_data_full

  !> @brief Generic mapping between two variables in memory, using
  !< an optional sign conversion
  subroutine map_data(this, controller_id, tgt_name, tgt_path, tgt_idx, &
                      src_name, src_path, src_idx, sign_array, &
                      lookup_table, stages)
    class(MapperType) :: this
    integer(I4B) :: controller_id
    character(len=*), intent(in) :: tgt_name
    character(len=*), intent(in) :: tgt_path
    integer(I4B), dimension(:), pointer :: tgt_idx
    character(len=*), intent(in) :: src_name
    character(len=*), intent(in) :: src_path
    integer(I4B), dimension(:), pointer :: src_idx
    integer(I4B), dimension(:), pointer :: sign_array
    integer(I4B), dimension(:), pointer :: lookup_table
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
    mapped_data%lut => lookup_table
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
