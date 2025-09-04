module MpiMessageBuilderModule
  use KindModule, only: I4B, LGP
  use ConstantsModule, only: LINELENGTH
  use SimModule, only: ustop
  use MemoryTypeModule, only: MemoryType
  use STLVecIntModule
  use VirtualBaseModule
  use VirtualDataContainerModule
  use mpi
  implicit none
  private

  type, public :: VdcHeaderType
    integer(I4B) :: id
    integer(I4B) :: container_type
    integer(I4B), dimension(NR_VDC_ELEMENT_MAPS) :: map_sizes
  end type

  type, public :: VdcReceiverMapsType
    type(VdcElementMapType), dimension(NR_VDC_ELEMENT_MAPS) :: el_maps
  contains
    procedure :: create
    procedure :: destroy
  end type

  type, public :: MpiMessageBuilderType
    type(VdcPtrType), dimension(:), pointer :: vdc_models => null() !< the models to be build the message for
    type(VdcPtrType), dimension(:), pointer :: vdc_exchanges => null() !< the exchanges to be build the message for
    integer(I4B) :: imon !< the output file unit, set from outside
  contains
    procedure :: init
    procedure :: attach_data
    procedure :: release_data
    procedure :: create_header_snd
    procedure :: create_header_rcv
    procedure :: create_map_snd
    procedure :: create_map_rcv
    procedure :: create_body_rcv
    procedure :: create_body_snd
    procedure :: set_monitor
    ! private
    procedure, private :: get_vdc_from_hdr
    procedure, private :: create_vdc_snd_hdr
    procedure, private :: create_vdc_snd_map
    procedure, private :: create_vdc_snd_body
    procedure, private :: create_vdc_rcv_body
  end type

contains

  subroutine create(this, map_sizes)
    class(VdcReceiverMapsType) :: this
    integer(I4B), dimension(NR_VDC_ELEMENT_MAPS) :: map_sizes
    ! local
    integer(I4B) :: i

    do i = 1, NR_VDC_ELEMENT_MAPS
      this%el_maps(i)%nr_virt_elems = map_sizes(i)
      allocate (this%el_maps(i)%remote_elem_shift(map_sizes(i)))
    end do

  end subroutine create

  subroutine destroy(this)
    class(VdcReceiverMapsType) :: this
    ! local
    integer(I4B) :: i

    do i = 1, NR_VDC_ELEMENT_MAPS
      if (associated(this%el_maps(i)%remote_elem_shift)) then
        deallocate (this%el_maps(i)%remote_elem_shift)
      end if
    end do

  end subroutine destroy

  subroutine init(this)
    class(MpiMessageBuilderType) :: this

    this%imon = -1

  end subroutine init

  subroutine attach_data(this, vdc_models, vdc_exchanges)
    class(MpiMessageBuilderType) :: this
    type(VdcPtrType), dimension(:), pointer :: vdc_models
    type(VdcPtrType), dimension(:), pointer :: vdc_exchanges

    this%vdc_models => vdc_models
    this%vdc_exchanges => vdc_exchanges

  end subroutine attach_data

  subroutine release_data(this)
    class(MpiMessageBuilderType) :: this

    this%vdc_models => null()
    this%vdc_exchanges => null()

  end subroutine release_data

  subroutine set_monitor(this, imon)
    class(MpiMessageBuilderType) :: this
    integer(I4B) :: imon

    this%imon = imon

  end subroutine set_monitor

  !> @brief Create the header data type to send to
  !! the remote process for this particular stage.
  !! From these data, the receiver can construct the
  !< body to send back to us.
  subroutine create_header_snd(this, rank, stage, hdrs_snd_type)
    class(MpiMessageBuilderType) :: this
    integer(I4B) :: rank
    integer(I4B) :: stage
    integer, intent(out) :: hdrs_snd_type
    ! local
    integer(I4B) :: i, offset, nr_types
    class(VirtualDataContainerType), pointer :: vdc
    integer :: ierr
    type(STLVecInt) :: model_idxs, exg_idxs
    integer, dimension(:), allocatable :: blk_cnts, types
    integer(kind=MPI_ADDRESS_KIND), dimension(:), allocatable :: displs

    call model_idxs%init()
    call exg_idxs%init()

    ! determine which containers to include
    do i = 1, size(this%vdc_models)
      vdc => this%vdc_models(i)%ptr
      if (vdc%is_active .and. vdc%orig_rank == rank) then
        call model_idxs%push_back(i)
      end if
    end do
    do i = 1, size(this%vdc_exchanges)
      vdc => this%vdc_exchanges(i)%ptr
      if (vdc%is_active .and. vdc%orig_rank == rank) then
        call exg_idxs%push_back(i)
      end if
    end do

    nr_types = model_idxs%size + exg_idxs%size
    allocate (blk_cnts(nr_types))
    allocate (types(nr_types))
    allocate (displs(nr_types))

    if (this%imon > 0) then
      write (this%imon, '(6x,a,*(i3))') "create headers for models: ", &
        model_idxs%get_values()
      write (this%imon, '(6x,a,*(i3))') "create headers for exchange: ", &
        exg_idxs%get_values()
    end if

    ! loop over containers
    do i = 1, model_idxs%size
      vdc => this%vdc_models(model_idxs%at(i))%ptr
      call MPI_Get_address(vdc%id, displs(i), ierr)
      blk_cnts(i) = 1
      types(i) = this%create_vdc_snd_hdr(vdc, stage)
    end do
    offset = model_idxs%size
    do i = 1, exg_idxs%size
      vdc => this%vdc_exchanges(exg_idxs%at(i))%ptr
      call MPI_Get_address(vdc%id, displs(i + offset), ierr)
      blk_cnts(i + offset) = 1
      types(i + offset) = this%create_vdc_snd_hdr(vdc, stage)
    end do

    ! create a MPI data type for the headers to send
    call MPI_Type_create_struct(nr_types, blk_cnts, displs, types, &
                                hdrs_snd_type, ierr)
    call MPI_Type_commit(hdrs_snd_type, ierr)
    do i = 1, nr_types
      call MPI_Type_free(types(i), ierr)
    end do

    call model_idxs%destroy()
    call exg_idxs%destroy()

    deallocate (blk_cnts)
    deallocate (types)
    deallocate (displs)

  end subroutine create_header_snd

  subroutine create_header_rcv(this, hdr_rcv_type)
    class(MpiMessageBuilderType) :: this
    integer, intent(out) :: hdr_rcv_type
    ! local
    integer :: ierr

    ! this will be for one data container, the mpi recv
    ! call will accept an array of them, no need to create
    ! an overarching contiguous type...
    call MPI_Type_contiguous(NR_VDC_ELEMENT_MAPS + 2, MPI_INTEGER, &
                             hdr_rcv_type, ierr)
    call MPI_Type_commit(hdr_rcv_type, ierr)

  end subroutine create_header_rcv

  subroutine create_map_snd(this, rank, stage, map_snd_type)
    class(MpiMessageBuilderType) :: this
    integer(I4B) :: rank
    integer(I4B) :: stage
    integer, intent(out) :: map_snd_type
    ! local
    integer(I4B) :: i, offset, nr_types
    class(VirtualDataContainerType), pointer :: vdc
    integer :: ierr
    type(STLVecInt) :: model_idxs, exg_idxs
    integer, dimension(:), allocatable :: blk_cnts, types
    integer(kind=MPI_ADDRESS_KIND), dimension(:), allocatable :: displs

    call model_idxs%init()
    call exg_idxs%init()

    ! determine which containers to include,
    ! currently models + exchanges
    do i = 1, size(this%vdc_models)
      vdc => this%vdc_models(i)%ptr
      if (vdc%is_active .and. vdc%orig_rank == rank) then
        call model_idxs%push_back(i)
      end if
    end do
    do i = 1, size(this%vdc_exchanges)
      vdc => this%vdc_exchanges(i)%ptr
      if (vdc%is_active .and. vdc%orig_rank == rank) then
        call exg_idxs%push_back(i)
      end if
    end do

    nr_types = model_idxs%size + exg_idxs%size
    allocate (blk_cnts(nr_types))
    allocate (types(nr_types))
    allocate (displs(nr_types))

    if (this%imon > 0) then
      write (this%imon, '(6x,a,*(i3))') "create maps for models: ", &
        model_idxs%get_values()
      write (this%imon, '(6x,a,*(i3))') "create maps for exchange: ", &
        exg_idxs%get_values()
    end if

    ! loop over containers
    do i = 1, model_idxs%size
      vdc => this%vdc_models(model_idxs%at(i))%ptr
      call MPI_Get_address(vdc%id, displs(i), ierr)
      blk_cnts(i) = 1
      types(i) = this%create_vdc_snd_map(vdc, stage)
    end do
    offset = model_idxs%size
    do i = 1, exg_idxs%size
      vdc => this%vdc_exchanges(exg_idxs%at(i))%ptr
      call MPI_Get_address(vdc%id, displs(i + offset), ierr)
      blk_cnts(i + offset) = 1
      types(i + offset) = this%create_vdc_snd_map(vdc, stage)
    end do

    ! create a compound MPI data type for the maps
    call MPI_Type_create_struct(nr_types, blk_cnts, displs, types, &
                                map_snd_type, ierr)
    call MPI_Type_commit(map_snd_type, ierr)

    ! free the subtypes
    do i = 1, nr_types
      call MPI_Type_free(types(i), ierr)
    end do

    call model_idxs%destroy()
    call exg_idxs%destroy()

    deallocate (blk_cnts)
    deallocate (types)
    deallocate (displs)

  end subroutine create_map_snd

  subroutine create_map_rcv(this, rcv_map, nr_headers, map_rcv_type)
    class(MpiMessageBuilderType) :: this
    type(VdcReceiverMapsType), dimension(:) :: rcv_map
    integer(I4B) :: nr_headers
    integer, intent(out) :: map_rcv_type
    ! local
    integer(I4B) :: i, j, nr_elems, type_cnt
    integer :: ierr, max_nr_maps
    integer, dimension(:), allocatable :: types
    integer(kind=MPI_ADDRESS_KIND), dimension(:), allocatable :: displs
    integer, dimension(:), allocatable :: blk_cnts

    max_nr_maps = nr_headers * NR_VDC_ELEMENT_MAPS
    allocate (types(max_nr_maps))
    allocate (displs(max_nr_maps))
    allocate (blk_cnts(max_nr_maps))

    type_cnt = 0
    do i = 1, nr_headers
      do j = 1, NR_VDC_ELEMENT_MAPS
        nr_elems = rcv_map(i)%el_maps(j)%nr_virt_elems
        if (nr_elems == 0) cycle

        type_cnt = type_cnt + 1
        call MPI_Get_address(rcv_map(i)%el_maps(j)%remote_elem_shift, &
                             displs(type_cnt), ierr)
        call MPI_Type_contiguous(nr_elems, MPI_Integer, types(type_cnt), ierr)
        blk_cnts(type_cnt) = 1
      end do
    end do

    call MPI_Type_create_struct(type_cnt, blk_cnts, displs, types, &
                                map_rcv_type, ierr)
    call MPI_Type_commit(map_rcv_type, ierr)

    deallocate (types)
    deallocate (displs)
    deallocate (blk_cnts)

  end subroutine create_map_rcv

  !> @brief Create the body to receive based on the headers
  !< that have been sent
  subroutine create_body_rcv(this, rank, stage, body_rcv_type)
    class(MpiMessageBuilderType) :: this
    integer(I4B) :: rank
    integer(I4B) :: stage
    integer, intent(out) :: body_rcv_type
    ! local
    integer(I4B) :: i, nr_types, offset
    class(VirtualDataContainerType), pointer :: vdc
    type(STLVecInt) :: model_idxs, exg_idxs
    integer :: ierr
    integer, dimension(:), allocatable :: types
    integer(kind=MPI_ADDRESS_KIND), dimension(:), allocatable :: displs
    integer, dimension(:), allocatable :: blk_cnts

    call model_idxs%init()
    call exg_idxs%init()

    ! gather all containers from this rank
    do i = 1, size(this%vdc_models)
      vdc => this%vdc_models(i)%ptr
      if (vdc%is_active .and. vdc%orig_rank == rank) then
        if (this%imon > 0) then
          write (this%imon, '(6x,a,i0)') "expecting model ", vdc%id
        end if
        call model_idxs%push_back(i)
      end if
    end do
    do i = 1, size(this%vdc_exchanges)
      vdc => this%vdc_exchanges(i)%ptr
      if (vdc%is_active .and. vdc%orig_rank == rank) then
        if (this%imon > 0) then
          write (this%imon, '(6x,a,i0)') "expecting exchange ", vdc%id
        end if
        call exg_idxs%push_back(i)
      end if
    end do

    nr_types = model_idxs%size + exg_idxs%size
    allocate (types(nr_types))
    allocate (displs(nr_types))
    allocate (blk_cnts(nr_types))

    ! loop over included containers
    do i = 1, model_idxs%size
      vdc => this%vdc_models(model_idxs%at(i))%ptr
      call MPI_Get_address(vdc%id, displs(i), ierr)
      types(i) = this%create_vdc_rcv_body(vdc, rank, stage)
      blk_cnts(i) = 1
    end do
    offset = model_idxs%size
    do i = 1, exg_idxs%size
      vdc => this%vdc_exchanges(exg_idxs%at(i))%ptr
      call MPI_Get_address(vdc%id, displs(i + offset), ierr)
      blk_cnts(i + offset) = 1
      types(i + offset) = this%create_vdc_rcv_body(vdc, rank, stage)
    end do

    ! create a MPI data type for the virtual data containers to receive
    call MPI_Type_create_struct(nr_types, blk_cnts, displs, types, &
                                body_rcv_type, ierr)
    call MPI_Type_commit(body_rcv_type, ierr)
    do i = 1, nr_types
      call MPI_Type_free(types(i), ierr)
    end do

    call model_idxs%destroy()
    call exg_idxs%destroy()
    deallocate (types)
    deallocate (displs)
    deallocate (blk_cnts)

  end subroutine create_body_rcv

  !> @brief Create the body to send based on the received headers
  !<
  subroutine create_body_snd(this, rank, stage, headers, maps, body_snd_type)
    class(MpiMessageBuilderType) :: this
    integer(I4B) :: rank
    integer(I4B) :: stage
    type(VdcHeaderType), dimension(:) :: headers
    type(VdcReceiverMapsType), dimension(:) :: maps
    integer, intent(out) :: body_snd_type
    ! local
    integer(I4B) :: i, nr_headers
    class(VirtualDataContainerType), pointer :: vdc
    integer :: ierr
    integer, dimension(:), allocatable :: types
    integer(kind=MPI_ADDRESS_KIND), dimension(:), allocatable :: displs
    integer, dimension(:), allocatable :: blk_cnts

    nr_headers = size(headers)
    allocate (types(nr_headers))
    allocate (displs(nr_headers))
    allocate (blk_cnts(nr_headers))

    do i = 1, nr_headers
      vdc => this%get_vdc_from_hdr(headers(i))
      call MPI_Get_address(vdc%id, displs(i), ierr)
      types(i) = this%create_vdc_snd_body(vdc, maps(i)%el_maps, rank, stage)
      blk_cnts(i) = 1
    end do

    ! create the list of virtual data containers to receive
    call MPI_Type_create_struct(nr_headers, blk_cnts, displs, &
                                types, body_snd_type, ierr)
    call MPI_Type_commit(body_snd_type, ierr)
    do i = 1, nr_headers
      call MPI_Type_free(types(i), ierr)
    end do

    deallocate (types)
    deallocate (displs)
    deallocate (blk_cnts)

  end subroutine create_body_snd

  !> @brief Create send header for virtual data container, relative
  !< to the field ...%id
  function create_vdc_snd_hdr(this, vdc, stage) result(new_type)
    class(MpiMessageBuilderType) :: this
    class(VirtualDataContainerType) :: vdc
    integer(I4B) :: stage
    integer :: new_type ! the created MPI datatype, uncommitted
    ! local
    integer :: i, ierr
    integer, dimension(NR_VDC_ELEMENT_MAPS + 2) :: blk_cnts
    integer(kind=MPI_ADDRESS_KIND), dimension(NR_VDC_ELEMENT_MAPS + 2) :: displs
    integer, dimension(NR_VDC_ELEMENT_MAPS + 2) :: types

    call MPI_Get_address(vdc%id, displs(1), ierr)
    types(1) = MPI_INTEGER
    blk_cnts(1) = 1
    call MPI_Get_address(vdc%container_type, displs(2), ierr)
    types(2) = MPI_INTEGER
    blk_cnts(2) = 1
    do i = 1, NR_VDC_ELEMENT_MAPS
      call MPI_Get_address(vdc%element_maps(i)%nr_virt_elems, displs(i + 2), ierr)
      types(i + 2) = MPI_INTEGER
      blk_cnts(i + 2) = 1
    end do

    ! rebase to id field
    displs = displs - displs(1)
    call MPI_Type_create_struct(NR_VDC_ELEMENT_MAPS + 2, blk_cnts, &
                                displs, types, new_type, ierr)
    call MPI_Type_commit(new_type, ierr)

  end function create_vdc_snd_hdr

  !> @brief Create a MPI datatype for sending the maps
  !< with the type relative to the id field
  function create_vdc_snd_map(this, vdc, stage) result(new_type)
    class(MpiMessageBuilderType) :: this
    class(VirtualDataContainerType), pointer :: vdc
    integer(I4B) :: stage
    integer :: new_type
    ! local
    integer(I4B) :: i, type_cnt
    integer :: n_elems, ierr
    integer(kind=MPI_ADDRESS_KIND) :: offset
    integer, dimension(:), allocatable :: types
    integer(kind=MPI_ADDRESS_KIND), dimension(:), allocatable :: displs
    integer, dimension(:), allocatable :: blk_cnts

    allocate (types(NR_VDC_ELEMENT_MAPS))
    allocate (displs(NR_VDC_ELEMENT_MAPS))
    allocate (blk_cnts(NR_VDC_ELEMENT_MAPS))

    ! displ relative to id field
    call MPI_Get_address(vdc%id, offset, ierr)

    type_cnt = 0
    do i = 1, NR_VDC_ELEMENT_MAPS
      n_elems = vdc%element_maps(i)%nr_virt_elems
      if (n_elems == 0) cycle ! only non-empty maps are sent

      type_cnt = type_cnt + 1
      call MPI_Get_address(vdc%element_maps(i)%remote_elem_shift, &
                           displs(type_cnt), ierr)
      call MPI_Type_contiguous(n_elems, MPI_INTEGER, types(type_cnt), ierr)
      call MPI_Type_commit(types(type_cnt), ierr)
      blk_cnts(type_cnt) = 1
      displs(type_cnt) = displs(type_cnt) - offset
    end do

    call MPI_Type_create_struct(type_cnt, blk_cnts, displs, types, &
                                new_type, ierr)
    call MPI_Type_commit(new_type, ierr)

    do i = 1, type_cnt
      call MPI_Type_free(types(i), ierr)
    end do

    deallocate (types)
    deallocate (displs)
    deallocate (blk_cnts)

  end function create_vdc_snd_map

  function create_vdc_rcv_body(this, vdc, rank, stage) result(new_type)
    class(MpiMessageBuilderType) :: this
    class(VirtualDataContainerType), pointer :: vdc
    integer(I4B) :: rank
    integer(I4B) :: stage
    integer :: new_type
    ! local
    type(STLVecInt) :: items
    integer :: ierr
    integer(kind=MPI_ADDRESS_KIND) :: offset
    integer, dimension(:), allocatable :: types
    integer(kind=MPI_ADDRESS_KIND), dimension(:), allocatable :: displs
    integer, dimension(:), allocatable :: blk_cnts
    integer(I4B) :: i
    class(VirtualDataType), pointer :: vd

    call items%init()
    call vdc%get_recv_items(stage, rank, items)
    !if (this%imon > 0) call vdc%print_items(this%imon, items)

    allocate (types(items%size))
    allocate (displs(items%size))
    allocate (blk_cnts(items%size))

    call MPI_Get_address(vdc%id, offset, ierr)

    do i = 1, items%size
      vd => get_virtual_data_from_list(vdc%virtual_data_list, items%at(i))
      call get_mpi_datatype(this, vd, displs(i), types(i))
      blk_cnts(i) = 1
      ! rebase w.r.t. id field
      displs(i) = displs(i) - offset
    end do

    call MPI_Type_create_struct(items%size, blk_cnts, displs, &
                                types, new_type, ierr)
    call MPI_Type_commit(new_type, ierr)

    do i = 1, items%size
      vd => get_virtual_data_from_list(vdc%virtual_data_list, items%at(i))
      call free_mpi_datatype(vd, types(i))
    end do

    deallocate (types)
    deallocate (displs)
    deallocate (blk_cnts)

    call items%destroy()

  end function create_vdc_rcv_body

  function create_vdc_snd_body(this, vdc, vdc_maps, rank, stage) result(new_type)
    class(MpiMessageBuilderType) :: this
    class(VirtualDataContainerType), pointer :: vdc
    type(VdcElementMapType), dimension(:) :: vdc_maps
    integer(I4B) :: rank
    integer(I4B) :: stage
    integer :: new_type
    ! local
    type(STLVecInt) :: items
    integer :: ierr
    integer(kind=MPI_ADDRESS_KIND) :: offset
    integer, dimension(:), allocatable :: types
    integer(kind=MPI_ADDRESS_KIND), dimension(:), allocatable :: displs
    integer, dimension(:), allocatable :: blk_cnts
    integer(I4B) :: i
    class(VirtualDataType), pointer :: vd
    integer(I4B), dimension(:), pointer, contiguous :: el_map

    call items%init()
    call vdc%get_send_items(stage, rank, items)
    !if (this%imon > 0) call vdc%print_items(this%imon, items)

    allocate (types(items%size))
    allocate (displs(items%size))
    allocate (blk_cnts(items%size))

    call MPI_Get_address(vdc%id, offset, ierr)

    do i = 1, items%size
      vd => get_virtual_data_from_list(vdc%virtual_data_list, items%at(i))
      if (vd%map_type > 0) then
        el_map => vdc_maps(vd%map_type)%remote_elem_shift
      else
        el_map => null()
      end if
      call get_mpi_datatype(this, vd, displs(i), types(i), el_map)
      blk_cnts(i) = 1
      ! rebase w.r.t. id field
      displs(i) = displs(i) - offset
    end do

    call MPI_Type_create_struct(items%size, blk_cnts, displs, &
                                types, new_type, ierr)
    call MPI_Type_commit(new_type, ierr)

    do i = 1, items%size
      vd => get_virtual_data_from_list(vdc%virtual_data_list, items%at(i))
      call free_mpi_datatype(vd, types(i))
    end do

    deallocate (types)
    deallocate (displs)
    deallocate (blk_cnts)

    call items%destroy()

  end function create_vdc_snd_body

  function get_vdc_from_hdr(this, header) result(vdc)
    class(MpiMessageBuilderType) :: this
    type(VdcHeaderType) :: header
    class(VirtualDataContainerType), pointer :: vdc
    ! local
    integer(I4B) :: i

    vdc => null()
    if (header%container_type == VDC_GWFMODEL_TYPE .or. &
        header%container_type == VDC_GWTMODEL_TYPE .or. &
        header%container_type == VDC_GWEMODEL_TYPE) then
      do i = 1, size(this%vdc_models)
        vdc => this%vdc_models(i)%ptr
        if (vdc%id == header%id) return
        vdc => null()
      end do
    else if (header%container_type == VDC_GWFEXG_TYPE .or. &
             header%container_type == VDC_GWTEXG_TYPE .or. &
             header%container_type == VDC_GWEEXG_TYPE) then
      do i = 1, size(this%vdc_exchanges)
        vdc => this%vdc_exchanges(i)%ptr
        if (vdc%id == header%id) return
        vdc => null()
      end do
    end if

  end function get_vdc_from_hdr

  !> @brief Local routine to get elemental mpi data types representing
  !! the virtual data items. Types are automatically committed unless
  !< they are primitives (e.g. MPI_INTEGER)
  subroutine get_mpi_datatype(this, virtual_data, el_displ, el_type, el_map_opt)
    use SimModule, only: ustop
    class(MpiMessageBuilderType) :: this
    class(VirtualDataType), pointer :: virtual_data
    integer(kind=MPI_ADDRESS_KIND) :: el_displ
    integer :: el_type
    integer(I4B), dimension(:), pointer, contiguous, optional :: el_map_opt !< optional, and can be null
    ! local
    type(MemoryType), pointer :: mt
    integer(I4B), dimension(:), pointer, contiguous :: el_map

    el_map => null()
    if (present(el_map_opt)) el_map => el_map_opt

    if (this%imon > 0) then
      if (.not. associated(el_map)) then
        write (this%imon, '(8x,2a,i0)') virtual_data%var_name, ' all ', &
          virtual_data%virtual_mt%isize
      else
        write (this%imon, '(8x,2a,i0)') virtual_data%var_name, &
          ' with map size ', size(el_map)
      end if
    end if

    mt => virtual_data%virtual_mt

    if (associated(mt%intsclr)) then
      call get_mpitype_for_int(mt, el_displ, el_type)
    else if (associated(mt%aint1d)) then
      call get_mpitype_for_int1d(mt, el_displ, el_type, el_map)
    else if (associated(mt%dblsclr)) then
      call get_mpitype_for_dbl(mt, el_displ, el_type)
    else if (associated(mt%adbl1d)) then
      call get_mpitype_for_dbl1d(mt, el_displ, el_type, el_map)
    else if (associated(mt%adbl2d)) then
      call get_mpitype_for_dbl2d(mt, el_displ, el_type, el_map)
    else
      write (*, *) 'unsupported datatype in MPI messaging for ', &
        virtual_data%var_name, virtual_data%mem_path
      call ustop()
    end if

  end subroutine get_mpi_datatype

  !> @brief Local routine to free elemental mpi data types representing
  !! the virtual data items. This can't be done generally, because some
  !< (scalar) types are primitive and freeing them causes nasty errors...
  subroutine free_mpi_datatype(virtual_data, el_type)
    class(VirtualDataType), pointer :: virtual_data
    integer :: el_type
    ! local
    type(MemoryType), pointer :: mt
    integer :: ierr

    mt => virtual_data%virtual_mt
    if (associated(mt%intsclr)) then
      ! type is MPI_INTEGER, don't free this!
      return
    else if (associated(mt%dblsclr)) then
      ! type is MPI_DOUBLE_PRECISION, don't free this!
      return
    else if (associated(mt%logicalsclr)) then
      ! type is MPI_LOGICAL, don't free this!
      return
    else
      ! all other types are freed here
      call MPI_Type_free(el_type, ierr)
      return
    end if

  end subroutine free_mpi_datatype

  subroutine get_mpitype_for_int(mem, el_displ, el_type)
    type(MemoryType), pointer :: mem
    integer(kind=MPI_ADDRESS_KIND) :: el_displ
    integer :: el_type
    ! local
    integer :: ierr

    call MPI_Get_address(mem%intsclr, el_displ, ierr)
    el_type = MPI_INTEGER
    ! no need to commit primitive type

  end subroutine get_mpitype_for_int

  subroutine get_mpitype_for_int1d(mem, el_displ, el_type, el_map)
    type(MemoryType), pointer :: mem
    integer(kind=MPI_ADDRESS_KIND) :: el_displ
    integer :: el_type
    integer, dimension(:), pointer :: el_map
    ! local
    integer :: ierr

    ! sanity check on map
    call check_map_int1d(mem, el_map)

    call MPI_Get_address(mem%aint1d, el_displ, ierr)
    if (associated(el_map)) then
      call MPI_Type_create_indexed_block( &
        size(el_map), 1, el_map, MPI_INTEGER, el_type, ierr)
    else
      call MPI_Type_contiguous(mem%isize, MPI_INTEGER, el_type, ierr)
    end if
    call MPI_Type_commit(el_type, ierr)

  end subroutine get_mpitype_for_int1d

  subroutine get_mpitype_for_dbl(mem, el_displ, el_type)
    type(MemoryType), pointer :: mem
    integer(kind=MPI_ADDRESS_KIND) :: el_displ
    integer :: el_type
    ! local
    integer :: ierr

    call MPI_Get_address(mem%dblsclr, el_displ, ierr)
    el_type = MPI_DOUBLE_PRECISION
    ! no need to commit primitive type

  end subroutine get_mpitype_for_dbl

  subroutine get_mpitype_for_dbl1d(mem, el_displ, el_type, el_map)
    type(MemoryType), pointer :: mem
    integer(kind=MPI_ADDRESS_KIND) :: el_displ
    integer :: el_type
    integer, dimension(:), pointer :: el_map
    ! local
    integer :: ierr

    ! sanity check on map
    call check_map_dbl1d(mem, el_map)

    call MPI_Get_address(mem%adbl1d, el_displ, ierr)
    if (associated(el_map)) then
      call MPI_Type_create_indexed_block( &
        size(el_map), 1, el_map, MPI_DOUBLE_PRECISION, el_type, ierr)
    else
      call MPI_Type_contiguous(mem%isize, MPI_DOUBLE_PRECISION, el_type, ierr)
    end if
    call MPI_Type_commit(el_type, ierr)

  end subroutine get_mpitype_for_dbl1d

  subroutine get_mpitype_for_dbl2d(mem, el_displ, el_type, el_map)
    type(MemoryType), pointer :: mem
    integer(kind=MPI_ADDRESS_KIND) :: el_displ
    integer :: el_type
    integer, dimension(:), pointer :: el_map
    ! local
    integer :: ierr
    integer :: entry_type

    ! sanity check on map
    call check_map_dbl2d(mem, el_map)

    call MPI_Get_address(mem%adbl2d, el_displ, ierr)
    if (associated(el_map)) then
      call MPI_Type_contiguous( &
        size(mem%adbl2d, dim=1), MPI_DOUBLE_PRECISION, entry_type, ierr)
      call MPI_Type_create_indexed_block( &
        size(el_map), 1, el_map, entry_type, el_type, ierr)
    else
      call MPI_Type_contiguous(mem%isize, MPI_DOUBLE_PRECISION, el_type, ierr)
    end if
    call MPI_Type_commit(el_type, ierr)

  end subroutine get_mpitype_for_dbl2d

  subroutine check_map_int1d(mem, map)
    type(MemoryType), pointer :: mem
    integer, dimension(:), pointer :: map !< ZERO-based map (for creating mpi types)
    ! local
    logical(LGP) :: is_valid
    integer(I4B) :: min_idx, max_idx

    if (.not. associated(map)) return
    if (size(map) == 0) return ! nothing to check

    ! bounds check
    min_idx = minval(map) + 1
    max_idx = maxval(map) + 1
    is_valid = max_idx <= size(mem%aint1d) .and. min_idx > 0
    if (.not. is_valid) then
      write (*, '(/,4x,4a)') &
        'Error: invalid map in MPI datatype for ', &
        trim(mem%name), ' in ', trim(mem%path)
      call ustop()
    end if

  end subroutine check_map_int1d

  !> @brief Bounds check for index maps,
  !< terminates on error.
  subroutine check_map_dbl1d(mem, map)
    type(MemoryType), pointer :: mem !< memory type
    integer, dimension(:), pointer :: map !< ZERO-based map (for creating mpi types)
    ! local
    logical(LGP) :: is_valid
    integer(I4B) :: min_idx, max_idx

    if (.not. associated(map)) return
    if (size(map) == 0) return ! nothing to check

    min_idx = minval(map) + 1
    max_idx = maxval(map) + 1
    is_valid = max_idx <= size(mem%adbl1d) .and. min_idx > 0
    if (.not. is_valid) then
      write (*, '(/,4x,4a)') &
        'Error: invalid map in MPI datatype for ', &
        trim(mem%name), ' in ', trim(mem%path)
      call ustop()
    end if

  end subroutine check_map_dbl1d

  subroutine check_map_dbl2d(mem, map)
    type(MemoryType), pointer :: mem
    integer, dimension(:), pointer :: map !< ZERO-based map (for creating mpi types)
    ! local
    logical(LGP) :: is_valid
    integer(I4B) :: min_idx, max_idx

    if (.not. associated(map)) return
    if (size(map) == 0) return ! nothing to check

    min_idx = minval(map) + 1
    max_idx = maxval(map) + 1
    is_valid = max_idx <= size(mem%adbl2d, 2) .and. min_idx > 0
    if (.not. is_valid) then
      write (*, '(/,4x,4a)') &
        'Error: invalid map in MPI datatype for ', &
        trim(mem%name), ' in ', trim(mem%path)
      call ustop()
    end if

  end subroutine check_map_dbl2d

end module MpiMessageBuilderModule
