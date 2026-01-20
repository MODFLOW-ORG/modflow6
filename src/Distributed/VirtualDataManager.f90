module VirtualDataManagerModule
  use KindModule, only: I4B
  use STLVecIntModule
  use VirtualDataListsModule, only: virtual_model_list, virtual_exchange_list
  use VirtualBaseModule, only: MAP_NODE_TYPE, MAP_CONN_TYPE
  use VirtualModelModule, only: VirtualModelType, get_virtual_model, &
                                get_virtual_model_from_list
  use VirtualExchangeModule, only: VirtualExchangeType, get_virtual_exchange, &
                                   get_virtual_exchange_from_list
  use VirtualSolutionModule
  use VirtualDataContainerModule
  use RouterBaseModule
  use RouterFactoryModule, only: create_router
  use ListsModule, only: basesolutionlist, baseconnectionlist
  use NumericalSolutionModule, only: NumericalSolutionType, &
                                     CastAsNumericalSolutionClass
  use NumericalModelModule, only: NumericalModelType, GetNumericalModelFromList
  use NumericalExchangeModule, only: NumericalExchangeType, &
                                     GetNumericalExchangeFromList
  use DisConnExchangeModule, only: DisConnExchangeType, &
                                   GetDisConnExchangeFromList
  use SpatialModelConnectionModule, only: SpatialModelConnectionType, &
                                          get_smc_from_list
  implicit none
  private

  type, public :: VirtualDataManagerType
    integer(I4B) :: nr_solutions
    integer(I4B), dimension(:), allocatable :: solution_ids
    type(VirtualSolutionType), dimension(:), pointer :: virtual_solutions
    class(RouterBaseType), pointer :: router
  contains
    procedure :: create => vds_create
    procedure :: init => vds_init
    procedure :: add_solution => vds_add_solution
    procedure :: activate_halo => vds_activate_halo
    procedure :: compress_halo => vds_compress_halo
    procedure :: synchronize => vds_synchronize
    procedure :: synchronize_sln => vds_synchronize_sln
    procedure :: destroy

    ! private
    procedure, private :: vds_synchronize
    procedure, private :: prepare_all
    procedure, private :: link_all
    procedure, private :: vds_synchronize_sln
    procedure, private :: prepare_sln
    procedure, private :: link_sln
    procedure, private :: count_nr_solutions
  end type VirtualDataManagerType

contains

  !> @brief Initialize the virtual data store
  subroutine vds_create(this, sim_mode)
    class(VirtualDataManagerType) :: this
    character(len=*) :: sim_mode
    ! local
    integer(I4B) :: nr_sol

    nr_sol = this%count_nr_solutions()
    allocate (this%virtual_solutions(nr_sol))
    allocate (this%solution_ids(nr_sol))

    ! we use this one as a counter:
    this%nr_solutions = 0

    ! create a router, sequential or parallel
    this%router => create_router(sim_mode, nr_sol)

  end subroutine vds_create

  !> @brief Initialize internal components
  !<
  subroutine vds_init(this)
    class(VirtualDataManagerType) :: this

    call this%router%initialize()

  end subroutine

  !> @brief Add the models and exchanges from the passed solution
  !! to the virtual data structure. This can then be used
  !< to efficiently sync only this particular solution.
  subroutine vds_add_solution(this, num_sol)
    class(VirtualDataManagerType) :: this
    class(NumericalSolutionType), pointer :: num_sol
    ! local
    integer(I4B) :: i, im, ix, ihm, ihx
    type(VirtualSolutionType), pointer :: virt_sol
    class(NumericalModelType), pointer :: num_mod
    class(DisConnExchangeType), pointer :: exg
    class(SpatialModelConnectionType), pointer :: conn
    integer(I4B) :: model_id, exg_id
    type(STLVecInt) :: model_ids, exchange_ids
    class(VirtualDataContainerType), pointer :: vdc
    logical :: found

    this%nr_solutions = this%nr_solutions + 1
    virt_sol => this%virtual_solutions(this%nr_solutions)

    call model_ids%init()
    call exchange_ids%init()

    ! build the virtual solution
    this%solution_ids(this%nr_solutions) = num_sol%id
    virt_sol%solution_id = num_sol%id
    virt_sol%numerical_solution => num_sol

    ! 1) adding all local models with a virtual model counterpart from the solution
    do im = 1, num_sol%modellist%Count()
      num_mod => GetNumericalModelFromList(num_sol%modellist, im)
      found = .false.
      do i = 1, virtual_model_list%Count()
        vdc => get_virtual_model_from_list(virtual_model_list, i)
        if (num_mod%id == vdc%id) then
          found = .true.
          exit
        end if
      end do
      if (found) then
        call model_ids%push_back(num_mod%id)
      end if
    end do

    ! 2) adding all local exchanges with a virtual exchange counterpart
    do ix = 1, num_sol%exchangelist%Count()
      exg => GetDisConnExchangeFromList(num_sol%exchangelist, ix)
      if (.not. associated(exg)) cycle ! interface model is handled separately
      found = .false.
      do i = 1, virtual_exchange_list%Count()
        vdc => get_virtual_exchange_from_list(virtual_exchange_list, i)
        if (exg%id == vdc%id) then
          found = .true.
          exit
        end if
      end do
      call exchange_ids%push_back_unique(exg%id)
    end do

    ! 3) add halo models and exchanges from interface models
    do ix = 1, num_sol%exchangelist%Count()
      conn => get_smc_from_list(num_sol%exchangelist, ix)
      if (.not. associated(conn)) cycle

      ! it's an interface model based exchanged, get
      ! halo models and halo exchanges from connection
      do ihm = 1, conn%halo_models%size
        model_id = conn%halo_models%at(ihm)
        call model_ids%push_back_unique(model_id)
      end do
      do ihx = 1, conn%halo_exchanges%size
        exg_id = conn%halo_exchanges%at(ihx)
        call exchange_ids%push_back_unique(exg_id)
      end do
    end do

    allocate (virt_sol%models(model_ids%size))
    allocate (virt_sol%exchanges(exchange_ids%size))
    allocate (virt_sol%interface_map)
    call virt_sol%interface_map%init(model_ids%size, exchange_ids%size)

    ! select virtual containers for models/exchanges
    do i = 1, model_ids%size
      vdc => get_virtual_model(model_ids%at(i))
      virt_sol%models(i)%ptr => vdc
    end do
    do i = 1, exchange_ids%size
      vdc => get_virtual_exchange(exchange_ids%at(i))
      virt_sol%exchanges(i)%ptr => vdc
    end do

    ! cleanup
    call model_ids%destroy()
    call exchange_ids%destroy()

  end subroutine vds_add_solution

  !> @brief Activates models and exchanges in the halo,
  !! i.e. the ones that have an actual chance of being used
  subroutine vds_activate_halo(this)
    class(VirtualDataManagerType) :: this
    ! local
    integer(I4B) :: im, ic, ix
    type(STLVecInt) :: halo_model_ids
    class(VirtualModelType), pointer :: vm
    class(VirtualExchangeType), pointer :: ve
    class(SpatialModelConnectionType), pointer :: conn

    call halo_model_ids%init()

    ! add halo models to list with ids (unique)
    do ic = 1, baseconnectionlist%Count()
      conn => get_smc_from_list(baseconnectionlist, ic)
      do im = 1, conn%halo_models%size
        call halo_model_ids%push_back_unique(conn%halo_models%at(im))
      end do
    end do

    ! deactivate models that are not local, and not in halo
    do im = 1, virtual_model_list%Count()
      vm => get_virtual_model_from_list(virtual_model_list, im)
      if (.not. vm%is_local) then
        if (.not. halo_model_ids%contains(vm%id)) then
          vm%is_active = .false.
        end if
      end if
    end do

    ! deactivate exchanges that are not local and outside halo
    ! (inside halo means both models are part of halo models)
    do ix = 1, virtual_exchange_list%Count()
      ve => get_virtual_exchange_from_list(virtual_exchange_list, ix)
      if (.not. ve%is_local) then
        if (.not. halo_model_ids%contains(ve%v_model1%id) .or. &
            .not. halo_model_ids%contains(ve%v_model2%id)) then
          ve%is_active = .false.
        end if
      end if
    end do

    this%router%halo_activated = .true.

    call halo_model_ids%destroy()

  end subroutine vds_activate_halo

  !> @brief Compress the halo for all solutions. This will
  !! activate the mapping tables in the virtual data items
  !< such that only relevant part of data arrays can be sync'ed
  subroutine vds_compress_halo(this)
    use ArrayHandlersModule, only: ifind
    use InputOutputModule, only: getunit
    use SimVariablesModule, only: proc_id
    use IndexMapModule
    class(VirtualDataManagerType) :: this
    ! local
    integer(I4B) :: ivm, isol, iexg, m_idx
    integer(I4B) :: outunit
    character(len=128) :: monitor_file
    type(VirtualSolutionType), pointer :: virt_sol
    class(NumericalSolutionType), pointer :: num_sol
    class(SpatialModelConnectionType), pointer :: conn
    class(VirtualDataContainerType), pointer :: vdc
    type(IndexMapType), pointer :: nmap, cmap

    ! merge the interface maps over this process
    do isol = 1, this%nr_solutions
      virt_sol => this%virtual_solutions(isol)
      num_sol => CastAsNumericalSolutionClass(virt_sol%numerical_solution)
      do iexg = 1, num_sol%exchangelist%Count()
        conn => get_smc_from_list(num_sol%exchangelist, iexg)
        if (.not. associated(conn)) cycle
        ! these are interface models, now merge their
        ! interface maps
        call virt_sol%interface_map%add(conn%interface_map)
      end do
    end do

    ! some testing
    if (.false.) then
      outunit = getunit()
      write (monitor_file, '(a,i0,a)') "iface.p", proc_id, ".log"
      open (unit=outunit, file=monitor_file)
      do isol = 1, this%nr_solutions
        write (outunit, '(a,i0,/)') "interface map for solution ", &
          this%virtual_solutions(isol)%solution_id
        virt_sol => this%virtual_solutions(isol)
        call virt_sol%interface_map%print_interface(outunit)
      end do
      close (outunit)
    end if

    ! assign reduced maps to virtual data containers
    do isol = 1, this%nr_solutions
      virt_sol => this%virtual_solutions(isol)
      do ivm = 1, size(virt_sol%models)
        vdc => virt_sol%models(ivm)%ptr
        if (.not. vdc%is_local .and. vdc%is_active) then
          m_idx = ifind(virt_sol%interface_map%model_ids, vdc%id)
          if (m_idx == -1) cycle

          nmap => virt_sol%interface_map%get_node_map(vdc%id)
          cmap => virt_sol%interface_map%get_connection_map(vdc%id)
          call vdc%set_element_map(nmap%src_idx, MAP_NODE_TYPE)
          call vdc%set_element_map(cmap%src_idx, MAP_CONN_TYPE)
        end if
      end do
    end do

  end subroutine vds_compress_halo

  !> @brief Synchronize the full virtual data store for this stage
  !<
  subroutine vds_synchronize(this, stage)
    class(VirtualDataManagerType) :: this !< this vdm
    integer(I4B) :: stage !< the stage to sync

    call this%prepare_all(stage)
    call this%link_all(stage)
    call this%router%route_all(stage)

  end subroutine vds_synchronize

  subroutine prepare_all(this, stage)
    class(VirtualDataManagerType) :: this
    integer(I4B) :: stage
    ! local
    integer(I4B) :: i
    class(VirtualDataContainerType), pointer :: vdc

    ! prepare all virtual data for this stage,
    ! cycle inactive to avoid redundant mem allocs
    do i = 1, virtual_model_list%Count()
      vdc => get_vdc_from_list(virtual_model_list, i)
      if (.not. vdc%is_active) cycle
      call vdc%prepare_stage(stage)
    end do
    do i = 1, virtual_exchange_list%Count()
      vdc => get_vdc_from_list(virtual_exchange_list, i)
      if (.not. vdc%is_active) cycle
      call vdc%prepare_stage(stage)
    end do

  end subroutine prepare_all

  subroutine link_all(this, stage)
    class(VirtualDataManagerType) :: this
    integer(I4B) :: stage
    ! local
    integer(I4B) :: i
    class(VirtualDataContainerType), pointer :: vdc

    ! link all local objects
    do i = 1, virtual_model_list%Count()
      vdc => get_vdc_from_list(virtual_model_list, i)
      call vdc%link_items(stage)
    end do
    do i = 1, virtual_exchange_list%Count()
      vdc => get_vdc_from_list(virtual_exchange_list, i)
      call vdc%link_items(stage)
    end do

  end subroutine link_all

  !> @brief Synchronize one particular solution for this stage
  !<
  subroutine vds_synchronize_sln(this, id_sln, stage)
    use ArrayHandlersModule, only: ifind
    class(VirtualDataManagerType) :: this
    integer(I4B) :: id_sln !< the id of the solution
    integer(I4B) :: stage
    ! local
    integer(I4B) :: sol_idx

    sol_idx = ifind(this%solution_ids, id_sln)
    call this%prepare_sln(this%virtual_solutions(sol_idx), stage)
    call this%link_sln(this%virtual_solutions(sol_idx), stage)
    call this%router%route_sln(this%virtual_solutions(sol_idx), stage)

  end subroutine vds_synchronize_sln

  !> @brief Force the virtual data containers (models,
  !! exchanges, etc.) to schedule their virtual data
  !< items for synchronization
  subroutine prepare_sln(this, virtual_sol, stage)
    class(VirtualDataManagerType) :: this
    type(VirtualSolutionType) :: virtual_sol
    integer(I4B) :: stage
    ! local
    integer(I4B) :: i
    class(VirtualDataContainerType), pointer :: vdc

    do i = 1, size(virtual_sol%models)
      vdc => virtual_sol%models(i)%ptr
      call vdc%prepare_stage(stage)
    end do

    do i = 1, size(virtual_sol%exchanges)
      vdc => virtual_sol%exchanges(i)%ptr
      call vdc%prepare_stage(stage)
    end do

  end subroutine prepare_sln

  !> @brief Connect virtual memory items to their
  !< sources when they are local for this stage
  subroutine link_sln(this, virtual_sol, stage)
    class(VirtualDataManagerType) :: this
    type(VirtualSolutionType) :: virtual_sol
    integer(I4B) :: stage
    ! local
    integer(I4B) :: i
    class(VirtualDataContainerType), pointer :: vdc

    do i = 1, size(virtual_sol%models)
      vdc => virtual_sol%models(i)%ptr
      call vdc%link_items(stage)
    end do

    do i = 1, size(virtual_sol%exchanges)
      vdc => virtual_sol%exchanges(i)%ptr
      call vdc%link_items(stage)
    end do

  end subroutine link_sln

  !> @brief Returns the number of Numerical Solutions
  !< in this simulation
  function count_nr_solutions(this) result(count)
    use ListsModule, only: basesolutionlist
    use NumericalSolutionModule, only: NumericalSolutionType
    class(VirtualDataManagerType) :: this
    integer(I4B) :: count
    ! local
    integer(I4B) :: isol
    class(*), pointer :: sol

    ! count nr. of numerical solutions
    count = 0
    do isol = 1, basesolutionlist%Count()
      sol => basesolutionlist%GetItem(isol)
      select type (sol)
      class is (NumericalSolutionType)
        count = count + 1
      end select
    end do

  end function count_nr_solutions

  subroutine destroy(this)
    class(VirtualDataManagerType) :: this
    ! local
    integer(I4B) :: i
    class(VirtualDataContainerType), pointer :: vdc

    do i = 1, virtual_model_list%Count()
      vdc => get_vdc_from_list(virtual_model_list, i)
      call vdc%destroy()
    end do
    call virtual_model_list%Clear(destroy=.true.)

    do i = 1, virtual_exchange_list%Count()
      vdc => get_vdc_from_list(virtual_exchange_list, i)
      call vdc%destroy()
    end do
    call virtual_exchange_list%Clear(destroy=.true.)

    do i = 1, this%nr_solutions
      deallocate (this%virtual_solutions(i)%models)
      deallocate (this%virtual_solutions(i)%exchanges)
      call this%virtual_solutions(i)%interface_map%destroy()
      deallocate (this%virtual_solutions(i)%interface_map)
    end do
    deallocate (this%virtual_solutions)

    call this%router%destroy()
    deallocate (this%router)

  end subroutine destroy

end module VirtualDataManagerModule
