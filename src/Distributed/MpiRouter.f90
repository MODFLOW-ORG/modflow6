module MpiRouterModule
  use RouterBaseModule
  use KindModule, only: I4B, LGP
  use ProfilerModule, only: g_prof
  use STLVecIntModule
  use SimVariablesModule, only: proc_id, nr_procs
  use SimStagesModule, only: STG_TO_STR, NR_SIM_STAGES
  use VirtualDataListsModule, only: virtual_model_list, &
                                    virtual_exchange_list
  use VirtualBaseModule, only: NR_VDC_ELEMENT_MAPS
  use VirtualDataContainerModule
  use VirtualExchangeModule, only: VirtualExchangeType
  use VirtualSolutionModule
  use MpiMessageBuilderModule
  use MpiMessageCacheModule
  use MpiWorldModule
  use mpi
  implicit none
  private

  public :: create_mpi_router

  type, public, extends(RouterBaseType) :: MpiRouterType
    integer(I4B), dimension(:), pointer :: model_proc_ids
    type(STLVecInt) :: senders !< the process ids to receive data from
    type(STLVecInt) :: receivers !< the process ids to send data to
    type(VdcPtrType), dimension(:), pointer :: all_models => null() !< all virtual models from the global list
    type(VdcPtrType), dimension(:), pointer :: all_exchanges => null() !< all virtual exchanges from the global list
    type(VdcPtrType), dimension(:), pointer :: rte_models => null() !< the currently active models to be routed
    type(VdcPtrType), dimension(:), pointer :: rte_exchanges => null() !< the currently active exchanges to be routed
    type(MpiMessageBuilderType) :: message_builder
    type(MpiMessageCacheType) :: msg_cache
    type(MpiWorldType), pointer :: mpi_world => null()
    integer(I4B) :: imon !< the output file unit for the mpi monitor
    logical(LGP) :: enable_monitor !< when true, log diagnostics
    integer(I4B), dimension(:, :), allocatable :: tmr_mpi_wait !< array with timer handles for MPI_Wait calls
  contains
    procedure :: initialize => mr_initialize
    procedure :: route_all => mr_route_all
    procedure :: route_sln => mr_route_sln
    procedure :: finalize => mr_finalize
    procedure :: destroy => mr_destroy
    ! private
    procedure, private :: activate
    procedure, private :: deactivate
    procedure, private :: update_senders
    procedure, private :: update_receivers
    procedure, private :: route_active
    procedure, private :: is_cached
    procedure, private :: compose_messages
    procedure, private :: load_messages
  end type MpiRouterType

contains

  !> Factory method to create MPI router
  !<
  function create_mpi_router() result(router)
    class(RouterBaseType), pointer :: router
    ! local
    class(MpiRouterType), pointer :: mpi_router

    allocate (mpi_router)
    router => mpi_router

  end function create_mpi_router

  subroutine mr_initialize(this)
    use InputOutputModule, only: getunit
    use ConstantsModule, only: LINELENGTH
    class(MpiRouterType) :: this
    ! local
    integer :: ierr
    integer(I4B) :: i
    integer(I4B) :: nr_models, nr_exchanges
    class(VirtualDataContainerType), pointer :: vdc
    character(len=LINELENGTH) :: monitor_file

    ! allocate timer handles: global + nr. solutions
    allocate (this%tmr_mpi_wait(NR_SIM_STAGES, this%nr_virt_solutions + 1))
    this%tmr_mpi_wait = -1

    ! routing over all when starting
    this%halo_activated = .false.

    ! to log or not to log
    this%enable_monitor = .false.

    ! initialize the MPI message builder and cache
    call this%message_builder%init()
    call this%msg_cache%init()

    ! get mpi world for our process
    this%mpi_world => get_mpi_world()

    ! init address list
    call this%senders%init()
    call this%receivers%init()

    ! find out where models are
    nr_models = virtual_model_list%Count()
    nr_exchanges = virtual_exchange_list%Count()
    allocate (this%model_proc_ids(nr_models))
    allocate (this%all_models(nr_models))
    allocate (this%all_exchanges(nr_exchanges))

    do i = 1, nr_models
      vdc => get_vdc_from_list(virtual_model_list, i)
      this%all_models(i)%ptr => vdc
      if (vdc%is_local) then
        this%model_proc_ids(i) = proc_id
      else
        this%model_proc_ids(i) = 0
      end if
    end do

    call MPI_Allreduce(MPI_IN_PLACE, this%model_proc_ids, nr_models, &
                       MPI_INTEGER, MPI_SUM, this%mpi_world%comm, ierr)
    call CHECK_MPI(ierr)

    ! set the process id to the models and exchanges
    do i = 1, nr_models
      vdc => get_vdc_from_list(virtual_model_list, i)
      call vdc%set_orig_rank(this%model_proc_ids(i))
    end do

    do i = 1, nr_exchanges
      vdc => get_vdc_from_list(virtual_exchange_list, i)
      this%all_exchanges(i)%ptr => vdc
      select type (vex => vdc)
      class is (VirtualExchangeType)
        call vex%set_orig_rank(vex%v_model1%orig_rank)
        if (vex%v_model1%is_local) then
          call vex%set_orig_rank(vex%v_model2%orig_rank)
        end if
      end select
    end do

    ! open log file
    if (this%enable_monitor) then
      this%imon = getunit()
      write (monitor_file, '(a,i0,a)') "mpi.p", proc_id, ".log"
      open (unit=this%imon, file=monitor_file)
      call this%message_builder%set_monitor(this%imon)

      ! write initial info
      write (this%imon, '(a,/)') ">> initialize MPI Router:"
      write (this%imon, '(2x,a,i0)') "process id: ", proc_id
      write (this%imon, '(2x,a,i0)') "nr. of processes: ", nr_procs
      write (this%imon, '(2x,a,i0)') "nr. of models: ", nr_models
      write (this%imon, '(2x,a,i0)') "nr. of exchanges: ", nr_exchanges
      write (this%imon, '(2x,2a)') "model id, processor id:"
      do i = 1, nr_models
        write (this%imon, '(4x,2i8)') i, this%model_proc_ids(i)
      end do
      write (this%imon, '(a,/)') "<< initialize done"
    end if

  end subroutine mr_initialize

  !> @brief Activate models and exchanges for routing
  !<
  subroutine activate(this, models, exchanges)
    class(MpiRouterType) :: this
    type(VdcPtrType), dimension(:), pointer :: models
    type(VdcPtrType), dimension(:), pointer :: exchanges

    this%rte_models => models
    this%rte_exchanges => exchanges
    call this%message_builder%attach_data(models, exchanges)

  end subroutine activate

  !> @brief Deactivate data after routing
  !<
  subroutine deactivate(this)
    class(MpiRouterType) :: this

    this%rte_models => null()
    this%rte_exchanges => null()
    call this%message_builder%release_data()

  end subroutine deactivate

  !> @brief This will route all remote data from the
  !! global models and exchanges over MPI, for a
  !< given stage
  subroutine mr_route_all(this, stage)
    use MemoryManagerModule, only: mem_print_detailed
    class(MpiRouterType) :: this
    integer(I4B) :: stage

    if (this%enable_monitor) then
      write (this%imon, '(/,2a)') ">> routing all: ", STG_TO_STR(stage)
    end if

    ! route all
    call this%activate(this%all_models, this%all_exchanges)
    call this%route_active(0, stage)
    call this%deactivate()

    if (this%enable_monitor) then
      write (this%imon, '(a,/)') "<< end routing all"
      !call mem_print_detailed(this%imon)
    end if

  end subroutine mr_route_all

  !> @brief This will route all remote data from models
  !! and exchanges in a particular solution over MPI,
  !< for a given stage
  subroutine mr_route_sln(this, virtual_sol, stage)
    class(MpiRouterType) :: this
    type(VirtualSolutionType) :: virtual_sol
    integer(I4B) :: stage

    if (this%enable_monitor) then
      write (this%imon, '(/,a,i0,2a)') ">> routing solution: ", &
        virtual_sol%solution_id, ", ", STG_TO_STR(stage)
    end if

    ! route for this solution
    call this%activate(virtual_sol%models, virtual_sol%exchanges)
    call this%route_active(virtual_sol%solution_id, stage)
    call this%deactivate()

    if (this%enable_monitor) then
      write (this%imon, '(a)') "<< end routing solution"
    end if

  end subroutine mr_route_sln

  !> @brief Routes the models and exchanges over MPI,
  !! either constructing the message bodies in a sequence
  !! of communication steps, or by loading from cache
  !< for efficiency.
  subroutine route_active(this, unit, stage)
    use SimModule, only: store_error
    class(MpiRouterType) :: this !< this mpi router
    integer(I4B) :: unit !< the solution id, or equal to 0 when global
    integer(I4B) :: stage !< the stage to route
    ! local
    integer(I4B) :: i
    integer(I4B) :: rnk
    integer :: ierr
    integer(kind=MPI_COUNT_KIND) :: msg_size !< need a longer int here, msg size can be > 2^31
    logical(LGP) :: from_cache
    ! mpi handles
    integer, dimension(:), allocatable :: rcv_req
    integer, dimension(:), allocatable :: snd_req
    integer, dimension(:, :), allocatable :: rcv_stat

    ! message body
    integer, dimension(:), allocatable :: body_rcv_t
    integer, dimension(:), allocatable :: body_snd_t

    ! update address list
    call this%update_senders()
    call this%update_receivers()

    ! allocate body data
    allocate (body_rcv_t(this%senders%size))
    allocate (body_snd_t(this%receivers%size))

    ! allocate handles
    allocate (rcv_req(this%senders%size))
    allocate (snd_req(this%receivers%size))
    allocate (rcv_stat(MPI_STATUS_SIZE, this%senders%size))

    ! always initialize request handles
    rcv_req = MPI_REQUEST_NULL
    snd_req = MPI_REQUEST_NULL

    if (this%enable_monitor) then
      write (this%imon, '(2x,a,*(i3))') "process ids sending data: ", &
        this%senders%get_values()
      write (this%imon, '(2x,a,*(i3))') "process ids receiving data: ", &
        this%receivers%get_values()
    end if

    ! from cache or build
    from_cache = this%is_cached(unit, stage)
    if (.not. from_cache) then
      call this%compose_messages(unit, stage, body_snd_t, body_rcv_t)
    else
      call this%load_messages(unit, stage, body_snd_t, body_rcv_t)
    end if

    if (this%enable_monitor) then
      write (this%imon, '(2x,a)') "== communicating bodies =="
    end if

    ! recv bodies
    do i = 1, this%senders%size
      rnk = this%senders%at(i)
      if (this%enable_monitor) then
        write (this%imon, '(4x,a,i0)') "receiving from process: ", rnk
      end if

      ! call extended type size function (*_x) to avoid overflow for large submodels
      call MPI_Type_size_x(body_rcv_t(i), msg_size, ierr)
      if (msg_size > 0) then
        call MPI_Irecv(MPI_BOTTOM, 1, body_rcv_t(i), rnk, stage, &
                       this%mpi_world%comm, rcv_req(i), ierr)
        call CHECK_MPI(ierr)
      end if

      if (this%enable_monitor) then
        write (this%imon, '(6x,a,i0)') "message body size: ", msg_size
      end if
    end do

    ! send bodies
    do i = 1, this%receivers%size
      rnk = this%receivers%at(i)
      if (this%enable_monitor) then
        write (this%imon, '(4x,a,i0)') "sending to process: ", rnk
      end if

      ! call extended type size function (*_x) to avoid overflow for large submodels
      call MPI_Type_size_x(body_snd_t(i), msg_size, ierr)
      if (msg_size > 0) then
        call MPI_Isend(MPI_Bottom, 1, body_snd_t(i), rnk, stage, &
                       this%mpi_world%comm, snd_req(i), ierr)
        call CHECK_MPI(ierr)
      end if

      if (this%enable_monitor) then
        write (this%imon, '(6x,a,i0)') "message body size: ", msg_size
      end if
      call flush (this%imon)
    end do

    ! wait for exchange of all messages
    call g_prof%start("MPI_WaitAll ("//trim(STG_TO_STR(stage))//")", &
                      this%tmr_mpi_wait(stage, unit + 1))
    call MPI_WaitAll(this%senders%size, rcv_req, rcv_stat, ierr)
    call g_prof%stop(this%tmr_mpi_wait(stage, unit + 1))
    call CHECK_MPI(ierr)

    deallocate (rcv_req, snd_req, rcv_stat)
    deallocate (body_rcv_t, body_snd_t)

  end subroutine route_active

  !> @brief Constructs the message bodies' MPI datatypes.
  !!
  !! Constructs the message bodies' MPI datatypes for a
  !! unit (a solution) and a given stage. This is done in
  !! a sequence of 6 steps (distributed over 3 phases):
  !!
  !! == synchronizing headers (VdcHeaderType) ==
  !!
  !!   step 1:
  !!     Receive headers from remote addresses requesting data
  !!     from virtual data containers (models, exchanges, ...) local to this process
  !!   step 2:
  !!     Send headers to remote addresses to indicate for which
  !!     virtual data containers (models, exchanges, ...) data is requested
  !!
  !! == synchronizing maps (VdcReceiverMapsType) ==
  !!
  !!   step 3:
  !!     Based on the received headers, receive element maps (which elements are
  !!     to be sent from a contiguous array) for outgoing data
  !!   step 4:
  !!     Send element maps to remote addresses to specify incoming data
  !!
  !! == construct message body data types (VirtualDataContainerType) ==
  !!
  !!   step 5:
  !!     Construct the message receive body based on the virtual
  !!     data items in the virtual data containers, and cache
  !!
  !!   step 6:
  !!     Construct the message send body, based on received header data and
  !!     and maps, from the virtual data containers, and cache
  !<
  subroutine compose_messages(this, unit, stage, body_snd_t, body_rcv_t)
    class(MpiRouterType) :: this
    integer(I4B) :: unit
    integer(I4B) :: stage
    integer, dimension(:) :: body_snd_t
    integer, dimension(:) :: body_rcv_t
    ! local
    integer(I4B) :: i, j, k
    integer(I4B) :: rnk
    integer :: ierr
    ! mpi handles
    integer, dimension(:), allocatable :: rcv_req
    integer, dimension(:), allocatable :: snd_req
    integer, dimension(:, :), allocatable :: rcv_stat
    ! message header
    integer(I4B) :: max_headers
    type(VdcHeaderType), dimension(:, :), allocatable :: headers
    integer, dimension(:), allocatable :: hdr_rcv_t
    integer, dimension(:), allocatable :: hdr_snd_t
    integer, dimension(:), allocatable :: hdr_rcv_cnt
    ! maps
    type(VdcReceiverMapsType), dimension(:, :), allocatable :: rcv_maps
    integer, dimension(:), allocatable :: map_rcv_t
    integer, dimension(:), allocatable :: map_snd_t

    ! allocate handles
    allocate (rcv_req(this%receivers%size))
    allocate (snd_req(this%senders%size))
    allocate (rcv_stat(MPI_STATUS_SIZE, this%receivers%size))

    ! init handles
    rcv_req = MPI_REQUEST_NULL
    snd_req = MPI_REQUEST_NULL

    ! allocate header data
    max_headers = size(this%rte_models) + size(this%rte_exchanges)
    allocate (hdr_rcv_t(this%receivers%size))
    allocate (hdr_snd_t(this%senders%size))
    allocate (headers(max_headers, this%receivers%size))
    allocate (hdr_rcv_cnt(this%receivers%size))

    ! allocate map data
    allocate (map_snd_t(this%senders%size))
    allocate (map_rcv_t(this%receivers%size))
    allocate (rcv_maps(max_headers, this%receivers%size)) ! for every header, we potentially need the maps

    if (this%enable_monitor) then
      write (this%imon, '(2x,a)') "== communicating headers =="
    end if

    ! first receive headers for outward data
    do i = 1, this%receivers%size
      rnk = this%receivers%at(i)
      if (this%enable_monitor) then
        write (this%imon, '(4x,a,i0)') "Ireceive header from process: ", rnk
      end if
      call this%message_builder%create_header_rcv(hdr_rcv_t(i))
      call MPI_Irecv(headers(:, i), max_headers, hdr_rcv_t(i), rnk, stage, &
                     this%mpi_world%comm, rcv_req(i), ierr)
      call CHECK_MPI(ierr)
    end do

    ! send header for incoming data
    do i = 1, this%senders%size
      rnk = this%senders%at(i)
      if (this%enable_monitor) then
        write (this%imon, '(4x,a,i0)') "send header to process: ", rnk
      end if
      call this%message_builder%create_header_snd(rnk, stage, hdr_snd_t(i))
      call MPI_Isend(MPI_BOTTOM, 1, hdr_snd_t(i), rnk, stage, &
                     this%mpi_world%comm, snd_req(i), ierr)
      call CHECK_MPI(ierr)
    end do

    ! wait for exchange of all headers
    call g_prof%start("MPI_WaitAll ("//trim(STG_TO_STR(stage))//")", &
                      this%tmr_mpi_wait(stage, unit + 1))
    call MPI_WaitAll(this%receivers%size, rcv_req, rcv_stat, ierr)
    call g_prof%stop(this%tmr_mpi_wait(stage, unit + 1))
    call CHECK_MPI(ierr)

    ! reinit handles
    rcv_req = MPI_REQUEST_NULL
    snd_req = MPI_REQUEST_NULL

    ! after WaitAll we can count incoming headers from statuses
    do i = 1, this%receivers%size
      call MPI_Get_count(rcv_stat(:, i), hdr_rcv_t(i), hdr_rcv_cnt(i), ierr)

      if (this%enable_monitor) then
        rnk = this%senders%at(i)
        write (this%imon, '(4x,a,i0)') "received headers from process: ", rnk
        write (this%imon, '(6x,a)') "expecting data for:"
        do j = 1, hdr_rcv_cnt(i)
          write (this%imon, '(6x,a,i0,a,a)') "id: ", headers(j, i)%id, &
            " type: ", trim(VDC_TYPE_TO_STR(headers(j, i)%container_type))
          write (this%imon, '(6x,a,99i6)') "map sizes: ", headers(j, i)%map_sizes
        end do
      end if
    end do

    ! clean up types
    do i = 1, this%receivers%size
      call MPI_Type_free(hdr_rcv_t(i), ierr)
    end do
    do i = 1, this%senders%size
      call MPI_Type_free(hdr_snd_t(i), ierr)
    end do

    if (this%enable_monitor) then
      write (this%imon, '(2x,a)') "== communicating maps =="
    end if

    ! allocate space for receiving maps
    do i = 1, this%receivers%size
      do j = 1, hdr_rcv_cnt(i)
        call rcv_maps(j, i)%create(headers(j, i)%map_sizes)
      end do
    end do

    ! receive maps
    do i = 1, this%receivers%size
      rnk = this%receivers%at(i)
      if (this%enable_monitor) then
        write (this%imon, '(4x,a,i0)') "Ireceive maps from process: ", rnk
      end if

      call this%message_builder%create_map_rcv(rcv_maps(:, i), hdr_rcv_cnt(i), &
                                               map_rcv_t(i))
      call MPI_Irecv(MPI_BOTTOM, 1, map_rcv_t(i), rnk, stage, &
                     this%mpi_world%comm, rcv_req(i), ierr)
      call CHECK_MPI(ierr)
    end do

    ! send maps
    do i = 1, this%senders%size
      rnk = this%senders%at(i)
      if (this%enable_monitor) then
        write (this%imon, '(4x,a,i0)') "send map to process: ", rnk
      end if

      call this%message_builder%create_map_snd(rnk, stage, map_snd_t(i))
      call MPI_Isend(MPI_BOTTOM, 1, map_snd_t(i), rnk, stage, &
                     this%mpi_world%comm, snd_req(i), ierr)
      call CHECK_MPI(ierr)
    end do

    ! wait on receiving maps
    call g_prof%start("MPI_WaitAll ("//trim(STG_TO_STR(stage))//")", &
                      this%tmr_mpi_wait(stage, unit + 1))
    call MPI_WaitAll(this%receivers%size, rcv_req, rcv_stat, ierr)
    call g_prof%stop(this%tmr_mpi_wait(stage, unit + 1))
    call CHECK_MPI(ierr)

    ! print maps
    if (this%enable_monitor) then
      do i = 1, this%receivers%size
        rnk = this%receivers%at(i)
        write (this%imon, '(4x,a,i0)') "received maps from process: ", rnk
        do j = 1, hdr_rcv_cnt(i)
          write (this%imon, '(6x,a,i0,a,a)') "id: ", headers(j, i)%id, &
            " type: ", trim(VDC_TYPE_TO_STR(headers(j, i)%container_type))
          do k = 1, NR_VDC_ELEMENT_MAPS
            write (this%imon, '(8x,i0, a,i0)') k, " nr. elements: ", &
              rcv_maps(j, i)%el_maps(k)%nr_virt_elems
            if (rcv_maps(j, i)%el_maps(k)%nr_virt_elems > 0) then
              write (this%imon, '(8x,*(i6))') &
                rcv_maps(j, i)%el_maps(k)%remote_elem_shift
            end if
          end do
        end do
      end do
    end if

    ! clean up types
    do i = 1, this%receivers%size
      call MPI_Type_free(map_rcv_t(i), ierr)
    end do
    do i = 1, this%senders%size
      call MPI_Type_free(map_snd_t(i), ierr)
    end do

    if (this%enable_monitor) then
      write (this%imon, '(2x,a)') "== composing message bodies =="
    end if

    ! construct recv bodies and cache
    do i = 1, this%senders%size
      rnk = this%senders%at(i)
      if (this%enable_monitor) then
        write (this%imon, '(4x,a,i0)') "build recv body for process: ", rnk
      end if

      call this%message_builder%create_body_rcv(rnk, stage, body_rcv_t(i))
      call this%msg_cache%put(unit, rnk, stage, MPI_BDY_RCV, body_rcv_t(i))
    end do

    ! construct send bodies and cache
    do i = 1, this%receivers%size
      rnk = this%receivers%at(i)
      if (this%enable_monitor) then
        write (this%imon, '(4x,a,i0)') "build send body for process: ", rnk
      end if

      call this%message_builder%create_body_snd( &
        rnk, stage, headers(1:hdr_rcv_cnt(i), i), &
        rcv_maps(:, i), body_snd_t(i))
      call this%msg_cache%put(unit, rnk, stage, MPI_BDY_SND, body_snd_t(i))
    end do

    ! clean up element maps
    do i = 1, this%receivers%size
      do j = 1, hdr_rcv_cnt(i)
        call rcv_maps(j, i)%destroy()
      end do
    end do

    deallocate (rcv_req, snd_req, rcv_stat)
    deallocate (hdr_rcv_t, hdr_snd_t, hdr_rcv_cnt)
    deallocate (headers)
    deallocate (map_rcv_t, map_snd_t)
    deallocate (rcv_maps)

  end subroutine compose_messages

  !> @brief Load the message body MPI datatypes from cache
  !<
  subroutine load_messages(this, unit, stage, body_snd_t, body_rcv_t)
    class(MpiRouterType) :: this
    integer(I4B) :: unit
    integer(I4B) :: stage
    integer, dimension(:), allocatable :: body_snd_t
    integer, dimension(:), allocatable :: body_rcv_t
    ! local
    integer(I4B) :: i, rnk

    if (this%enable_monitor) then
      write (this%imon, '(2x,a)') "... running from cache ..."
    end if

    do i = 1, this%receivers%size
      rnk = this%receivers%at(i)
      body_snd_t(i) = this%msg_cache%get(unit, rnk, stage, MPI_BDY_SND)
    end do
    do i = 1, this%senders%size
      rnk = this%senders%at(i)
      body_rcv_t(i) = this%msg_cache%get(unit, rnk, stage, MPI_BDY_RCV)
    end do

  end subroutine load_messages

  subroutine update_senders(this)
    class(MpiRouterType) :: this
    ! local
    integer(I4B) :: i
    class(VirtualDataContainerType), pointer :: vdc

    call this%senders%clear()

    do i = 1, size(this%rte_models)
      vdc => this%rte_models(i)%ptr
      if (.not. vdc%is_local .and. vdc%is_active) then
        call this%senders%push_back_unique(vdc%orig_rank)
      end if
    end do
    do i = 1, size(this%rte_exchanges)
      vdc => this%rte_exchanges(i)%ptr
      if (.not. vdc%is_local .and. vdc%is_active) then
        call this%senders%push_back_unique(vdc%orig_rank)
      end if
    end do

  end subroutine update_senders

  subroutine update_receivers(this)
    class(MpiRouterType) :: this
    ! local
    integer(I4B) :: i, j
    class(VirtualDataContainerType), pointer :: vdc

    call this%receivers%clear()

    if (.not. this%halo_activated) then
      ! assuming symmetry for now
      do i = 1, this%senders%size
        call this%receivers%push_back(this%senders%at(i))
      end do
    else
      ! get the receivers from the VDCs
      do i = 1, size(this%rte_models)
        vdc => this%rte_models(i)%ptr
        do j = 1, vdc%rcv_ranks%size
          call this%receivers%push_back_unique(vdc%rcv_ranks%at(j))
        end do
      end do
      do i = 1, size(this%rte_exchanges)
        vdc => this%rte_exchanges(i)%ptr
        do j = 1, vdc%rcv_ranks%size
          call this%receivers%push_back_unique(vdc%rcv_ranks%at(j))
        end do
      end do
    end if

  end subroutine update_receivers

  !> @brief Check if this stage is cached
  !<
  function is_cached(this, unit, stage) result(in_cache)
    use SimModule, only: ustop
    class(MpiRouterType) :: this
    integer(I4B) :: unit
    integer(I4B) :: stage
    logical(LGP) :: in_cache
    ! local
    integer(I4B) :: i, rnk
    integer(I4B) :: no_cache_cnt
    integer :: cached_type

    in_cache = .false.
    no_cache_cnt = 0

    do i = 1, this%receivers%size
      rnk = this%receivers%at(i)
      cached_type = this%msg_cache%get(unit, rnk, stage, MPI_BDY_SND)
      if (cached_type == NO_CACHED_VALUE) no_cache_cnt = no_cache_cnt + 1
    end do
    do i = 1, this%senders%size
      rnk = this%senders%at(i)
      cached_type = this%msg_cache%get(unit, rnk, stage, MPI_BDY_RCV)
      if (cached_type == NO_CACHED_VALUE) no_cache_cnt = no_cache_cnt + 1
    end do

    ! it should be all or nothing
    if (no_cache_cnt == this%receivers%size + this%senders%size) then
      in_cache = .false.
    else if (no_cache_cnt == 0) then
      in_cache = .true.
    else
      call ustop("Internal error: MPI message cache corrupt...")
    end if

  end function is_cached

  subroutine mr_finalize(this)
    class(MpiRouterType) :: this

    call this%msg_cache%clear()

  end subroutine mr_finalize

  subroutine mr_destroy(this)
    class(MpiRouterType) :: this

    call this%msg_cache%destroy()

    call this%senders%destroy()
    call this%receivers%destroy()

    deallocate (this%model_proc_ids)
    deallocate (this%all_models)
    deallocate (this%all_exchanges)

    deallocate (this%tmr_mpi_wait)

  end subroutine mr_destroy

end module MpiRouterModule
