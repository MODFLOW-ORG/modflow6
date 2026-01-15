!> Refactoring issues towards parallel:
!!
!! * remove camelCase
!<
module GridConnectionModule
  use KindModule, only: I4B, DP, LGP
  use SimModule, only: ustop, store_error, store_warning
  use ConstantsModule, only: LENMEMPATH, DZERO, DPIO180, &
                             LENMODELNAME, LINELENGTH
  use CharacterStringModule
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, mem_reallocate
  use MemoryHelperModule, only: create_mem_path
  use ListModule, only: ListType, isEqualIface
  use NumericalModelModule
  use DisuModule
  use DisConnExchangeModule
  use VirtualModelModule, only: VirtualModelType, get_virtual_model, &
                                get_virtual_model_from_list
  use VirtualExchangeModule, only: VirtualExchangeType, get_virtual_exchange
  use CellWithNbrsModule
  use ConnectionsModule
  use SparseModule, only: sparsematrix
  use InterfaceMapModule
  use BaseDisModule, only: dis_transform_xy
  use CsrUtilsModule
  use STLVecIntModule
  implicit none
  private

  ! Initial nr of neighbors for sparse matrix allocation
  integer(I4B), parameter :: InitNrNeighbors = 7

  !> This class is used to construct the connections object for
  !! the interface model's spatial discretization/grid.
  !!
  !! It works as follows:
  !!
  !! 1: construct basic instance, allocate data structures
  !!    based on nr. of primary connections
  !! 2: add primary connections from exchanges
  !! 3: add secondary, tertiary, ... MODEL connections, depending
  !!    on the size of the computational stencil
  !! 4: extend the connection, creating the full data structure
  !!    and relations
  !! 5: build the connections object, from which the grid for
  !!    the interface model can be constructed
  !!
  !! A note on the different indices:
  !!
  !! We have
  !!
  !! - GLOBAL index, which technically labels the row in the solution matrix
  !! - REGIONAL index, running over all models that participate in the interface grid
  !! - LOCAL index, local to each model
  !! - INTERFACE index, numbering the cells in the interface grid
  !<
  type, public :: GridConnectionType

    character(len=LENMEMPATH) :: memoryPath
    integer(I4B) :: internalStencilDepth !< stencil size for the interior
    integer(I4B) :: exchangeStencilDepth !< stencil size at the interface
    integer(I4B) :: icondir !< icondir as defined in DIS, but reduced over all models in the interface

    class(NumericalModelType), pointer :: model => null() !< the model for which this grid connection exists
    class(DisConnExchangeType), pointer :: primaryExchange => null() !< pointer to the primary exchange for this interface

    integer(I4B), pointer :: nrOfBoundaryCells => null() !< nr of boundary cells with connection to another model
    type(CellWithNbrsType), dimension(:), pointer :: boundaryCells => null() !< cells on our side of the primary connections
    type(CellWithNbrsType), dimension(:), pointer :: connectedCells => null() !< cells on the neighbors side of the primary connection
    type(STLVecInt), pointer :: haloExchanges !< all exchanges that are potentially part of this interface

    integer(I4B), pointer :: nrOfCells => null() !< the total number of cells in the interface
    type(GlobalCellType), dimension(:), pointer :: idxToGlobal => null() !< a map from interface index to global coordinate
    integer(I4B), dimension(:), pointer, contiguous :: idxToGlobalIdx => null() !< a (flat) map from interface index to global index,
                                                                                !! stored in mem. mgr. so can be used for debugging
    type(ListType) :: regionalModels !< the models participating in the interface
    integer(I4B), dimension(:), pointer :: region_to_iface_map => null() !< (sparse) mapping from regional index to interface ixd
    integer(I4B), dimension(:), pointer :: regionalModelOffset => null() !< the new offset to compactify the range of indices
    integer(I4B), pointer :: indexCount => null() !< counts the number of cells in the interface

    type(ConnectionsType), pointer :: connections => null() !< sparse matrix with the connections
    integer(I4B), dimension(:), pointer :: connectionMask => null() !< to mask out connections from the amat coefficient calculation

    type(InterfaceMapType), pointer :: interfaceMap => null() !< defining map for the interface

  contains

    ! public
    procedure, pass(this) :: construct
    procedure, private, pass(this) :: allocateScalars
    procedure, pass(this) :: destroy
    procedure, pass(this) :: addToRegionalModels
    procedure, pass(this) :: connectPrimaryExchange
    procedure, pass(this) :: extendConnection
    procedure, pass(this) :: getDiscretization
    procedure, pass(this) :: buildInterfaceMap

    ! 'protected'
    procedure, pass(this) :: isPeriodic

    ! private routines
    procedure, private, pass(this) :: connectCell
    procedure, private, pass(this) :: buildConnections
    procedure, private, pass(this) :: addNeighbors
    procedure, private, pass(this) :: addNeighborCell
    procedure, private, pass(this) :: addRemoteNeighbors
    procedure, private, pass(this) :: get_regional_offset
    generic, private :: getInterfaceIndex => getInterfaceIndexByCell, &
      getInterfaceIndexByIndexModel
    procedure, private, pass(this) :: getInterfaceIndexByCell
    procedure, private, pass(this) :: getInterfaceIndexByIndexModel
    procedure, private, pass(this) :: registerInterfaceCells
    procedure, private, pass(this) :: addToGlobalMap
    procedure, private, pass(this) :: compressGlobalMap
    procedure, private, pass(this) :: sortInterfaceGrid
    procedure, private, pass(this) :: makePrimaryConnections
    procedure, private, pass(this) :: connectNeighborCells
    procedure, private, pass(this) :: fillConnectionDataInternal
    procedure, private, pass(this) :: fillConnectionDataFromExchanges
    procedure, private, pass(this) :: createConnectionMask
    procedure, private, pass(this) :: maskInternalConnections
    procedure, private, pass(this) :: setMaskOnConnection
  end type

contains

  !> @brief Construct the GridConnection and allocate
  !! the data structures for the primary connections
  !<
  subroutine construct(this, model, nrOfPrimaries, connectionName)
    class(GridConnectionType), intent(inout) :: this !> this instance
    class(NumericalModelType), pointer, intent(in) :: model !> the model for which the interface is constructed
    integer(I4B) :: nrOfPrimaries !> the number of primary connections between the two models
    character(len=*) :: connectionName !> the name, for memory management mostly
    ! local
    class(VirtualModelType), pointer :: v_model

    this%model => model
    this%memoryPath = create_mem_path(connectionName, 'GC')

    call this%allocateScalars()

    allocate (this%boundaryCells(nrOfPrimaries))
    allocate (this%connectedCells(nrOfPrimaries))
    allocate (this%idxToGlobal(2 * nrOfPrimaries))

    v_model => get_virtual_model(model%id)
    call this%addToRegionalModels(v_model)

    this%nrOfBoundaryCells = 0

    this%internalStencilDepth = 1
    this%exchangeStencilDepth = 1
    this%icondir = 1 ! start with assumption that vertex and cell2d data are available in all models
    this%haloExchanges => null()

  end subroutine construct

  !> @brief Make connections for the primary exchange
  !<
  subroutine connectPrimaryExchange(this, primEx)
    class(GridConnectionType) :: this !< this grid connection
    class(DisConnExchangeType), pointer :: primEx !< the primary exchange for this connection
    ! local
    integer(I4B) :: iconn

    ! store the primary exchange
    this%primaryExchange => primEx

    ! connect the cells
    do iconn = 1, primEx%nexg
      call this%connectCell(primEx%nodem1(iconn), primEx%v_model1, &
                            primEx%nodem2(iconn), primEx%v_model2)
    end do

  end subroutine connectPrimaryExchange

  !> @brief Connect neighboring cells at the interface by
  !! storing them in the boundary cell and connected cell
  !! arrays
  !<
  subroutine connectCell(this, idx1, v_model1, idx2, v_model2)
    class(GridConnectionType), intent(in) :: this !< this grid connection
    integer(I4B) :: idx1 !< local index cell 1
    class(VirtualModelType), pointer :: v_model1 !< model of cell 1
    integer(I4B) :: idx2 !< local index cell 2
    class(VirtualModelType), pointer :: v_model2 !< model of cell 2
    ! local
    type(GlobalCellType), pointer :: bnd_cell, conn_cell

    this%nrOfBoundaryCells = this%nrOfBoundaryCells + 1
    if (this%nrOfBoundaryCells > size(this%boundaryCells)) then
      write (*, *) 'Error: nr of cell connections exceeds '// &
        'capacity in grid connection, terminating...'
      call ustop()
    end if

    bnd_cell => this%boundaryCells(this%nrOfBoundaryCells)%cell
    conn_cell => this%connectedCells(this%nrOfBoundaryCells)%cell
    if (v_model1 == this%model) then
      bnd_cell%index = idx1
      bnd_cell%v_model => v_model1
      conn_cell%index = idx2
      conn_cell%v_model => v_model2
    else if (v_model2 == this%model) then
      bnd_cell%index = idx2
      bnd_cell%v_model => v_model2
      conn_cell%index = idx1
      conn_cell%v_model => v_model1
    else
      write (*, *) 'Error: unable to connect cells outside the model'
      call ustop()
    end if

  end subroutine connectCell

  !> @brief Add a model to a list of all regional models
  !<
  subroutine addToRegionalModels(this, v_model)
    class(GridConnectionType), intent(inout) :: this !< this grid connection
    class(VirtualModelType), pointer :: v_model !< the model to add to the region
    ! local
    class(*), pointer :: vm_obj

    vm_obj => v_model
    if (.not. this%regionalModels%ContainsObject(vm_obj)) then
      call this%regionalModels%Add(vm_obj)
    end if

  end subroutine addToRegionalModels

  !> @brief Extend the connection topology to deal with
  !! higher levels of connectivity (neighbors-of-neighbors, etc.)
  !!
  !! The following steps are taken:
  !! 1. Recursively add interior neighbors (own model) up to the specified depth
  !! 2. Recursively add exterior neighbors
  !! 3. Allocate a (sparse) mapping table for the region
  !! 4. Build connection object for the interface grid, and the mask
  !<
  subroutine extendConnection(this)
    class(GridConnectionType), intent(inout) :: this !< this grid connection
    ! local
    integer(I4B) :: remoteDepth, localDepth
    integer(I4B) :: icell
    integer(I4B) :: imod, regionSize, offset
    class(VirtualModelType), pointer :: v_model
    !integer(I4B), pointer :: nr_nodes

    ! we need (stencildepth-1) extra cells for the interior
    remoteDepth = this%exchangeStencilDepth
    localDepth = 2 * this%internalStencilDepth - 1
    if (localDepth < remoteDepth) then
      localDepth = remoteDepth
    end if

    ! first add the neighbors for the interior
    ! (possibly extending into other models)
    do icell = 1, this%nrOfBoundaryCells
      call this%addNeighbors(this%boundaryCells(icell), localDepth, &
                             this%connectedCells(icell)%cell, .true.)
    end do
    ! and for the exterior
    do icell = 1, this%nrOfBoundaryCells
      call this%addNeighbors(this%connectedCells(icell), remoteDepth, &
                             this%boundaryCells(icell)%cell, .false.)
    end do

    ! set up mapping for the region (models participating in interface model grid)
    allocate (this%regionalModelOffset(this%regionalModels%Count()))
    regionSize = 0
    offset = 0
    do imod = 1, this%regionalModels%Count()
      v_model => get_virtual_model_from_list(this%regionalModels, imod)
      regionSize = regionSize + v_model%dis_nodes%get()
      this%regionalModelOffset(imod) = offset
      offset = offset + v_model%dis_nodes%get()
    end do
    ! init to -1, meaning 'interface index was not assigned yet'
    allocate (this%region_to_iface_map(regionSize))
    this%region_to_iface_map = -1

    call this%buildConnections()

  end subroutine extendConnection

  !> @brief Builds a sparse matrix holding all cell connections,
  !< with new indices, and stores the mapping to the global ids
  subroutine buildConnections(this)
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    ! local
    integer(I4B) :: i, icell, iconn
    integer(I4B), dimension(:), allocatable :: nnz
    type(SparseMatrix), pointer :: sparse
    integer(I4B) :: ierror
    type(ConnectionsType), pointer :: conn
    class(VirtualModelType), pointer :: vm
    character(len=LINELENGTH) :: warnmsg

    ! Recursively generate interface cell indices, fill map to global cells,
    ! and add to region lookup table
    this%indexCount = 0
    do icell = 1, this%nrOfBoundaryCells
      call this%registerInterfaceCells(this%boundaryCells(icell))
    end do
    do icell = 1, this%nrOfBoundaryCells
      call this%registerInterfaceCells(this%connectedCells(icell))
    end do
    this%nrOfCells = this%indexCount

    ! compress lookup table
    call this%compressGlobalMap()

    ! prepare to build vertex and cell2d info, or not
    do i = 1, this%regionalModels%Count()
      vm => get_virtual_model_from_list(this%regionalModels, i)
      if (vm%dis_icondir%get() == 0) this%icondir = 0
    end do

    if (this%icondir == 0 .and. this%model%dis%icondir > 0) then
      write (warnmsg, '(3a)') 'Exchange ', trim(this%primaryExchange%name), &
        " will not use vertex and/or cell coordinate data because not all &
        &connected models have those data specified"
      call store_warning(warnmsg)
    end if

    ! sort interface indexes such that 'n > m' means 'n below m'
    call this%sortInterfaceGrid()

    ! allocate a map from interface index to global coordinates
    call mem_allocate(this%idxToGlobalIdx, this%nrOfCells, &
                      'IDXTOGLOBALIDX', this%memoryPath)

    ! create sparse data structure, to temporarily hold connections
    allocate (sparse)
    allocate (nnz(this%nrOfCells))
    nnz = InitNrNeighbors + 1
    call sparse%init(this%nrOfCells, this%nrOfCells, nnz)

    ! now (recursively) add connections to sparse, start with
    ! the primary connections (n-m from the exchange files)
    call this%makePrimaryConnections(sparse)
    ! then into own domain
    do icell = 1, this%nrOfBoundaryCells
      call this%connectNeighborCells(this%boundaryCells(icell), sparse)
    end do
    ! and same for the neighbors of connected cells
    do icell = 1, this%nrOfBoundaryCells
      call this%connectNeighborCells(this%connectedCells(icell), sparse)
    end do

    ! create connections object
    allocate (this%connections)
    conn => this%connections
    call conn%allocate_scalars(this%memoryPath)
    conn%nodes = this%nrOfCells
    conn%nja = sparse%nnz
    conn%njas = (conn%nja - conn%nodes) / 2
    call conn%allocate_arrays()
    do iconn = 1, conn%njas
      conn%anglex(iconn) = -999.
    end do

    ! fill connection from sparse
    call sparse%filliaja(conn%ia, conn%ja, ierror)
    if (ierror /= 0) then
      write (*, *) 'Error filling ia/ja in GridConnection: terminating...'
      call ustop()
    end if
    call fillisym(conn%nodes, conn%nja, conn%ia, conn%ja, conn%isym)
    call filljas(conn%nodes, conn%nja, conn%ia, conn%ja, conn%isym, conn%jas)
    call sparse%destroy()

    ! fill connection data (ihc, cl1, cl2, etc.) using data
    ! from models and exchanges
    call this%fillConnectionDataInternal()
    call this%fillConnectionDataFromExchanges()

    ! set the masks on connections
    call this%createConnectionMask()

  end subroutine buildConnections

  !< @brief Routine for finding neighbors-of-neighbors, recursively
  !<
  recursive subroutine addNeighbors(this, cellNbrs, depth, mask, interior)
    use SimModule, only: ustop
    class(GridConnectionType), intent(inout) :: this !< this grid connection
    type(CellWithNbrsType), intent(inout), target :: cellNbrs !< cell to add to
    integer(I4B), intent(inout) :: depth !< current depth (typically decreases in recursion)
    type(GlobalCellType), optional :: mask !< mask to excluded back-and-forth connection between cells
    logical(LGP) :: interior !< when true, we are adding from the exchange back into the model
    ! local
    type(GlobalCellType), pointer :: cell
    integer(I4B) :: ipos, ipos_start, ipos_end
    integer(I4B) :: nbrIdx, inbr
    integer(I4B) :: newDepth

    ! readability
    cell => cellNbrs%cell

    ! if depth == 1, then we are not adding neighbors but use
    ! the boundary and connected cell only
    if (depth < 2) then
      return
    end if
    newDepth = depth - 1

    ! find neighbors local to this cell by looping through grid connections
    ipos_start = cell%v_model%con_ia%get(cell%index) + 1
    ipos_end = cell%v_model%con_ia%get(cell%index + 1) - 1
    do ipos = ipos_start, ipos_end
      nbrIdx = cell%v_model%con_ja%get(ipos)
      call this%addNeighborCell(cellNbrs, nbrIdx, cellNbrs%cell%v_model, mask)
    end do

    ! add remote nbr using the data from the exchanges
    call this%addRemoteNeighbors(cellNbrs, mask)

    ! now find nbr-of-nbr
    do inbr = 1, cellNbrs%nrOfNbrs

      ! are we leaving the model through another exchange?
      if (interior .and. cellNbrs%cell%v_model == this%model) then
        if (.not. cellNbrs%neighbors(inbr)%cell%v_model == this%model) then
          ! decrement by 1, because the connection we are crossing is not
          ! calculated by this interface
          newDepth = newDepth - 1
        end if
      end if
      ! and add neighbors with the new depth
      call this%addNeighbors(cellNbrs%neighbors(inbr), newDepth, &
                             cellNbrs%cell, interior)
    end do

  end subroutine addNeighbors

  !> @brief Add cell neighbors across models using the stored exchange
  !! data structures
  subroutine addRemoteNeighbors(this, cellNbrs, mask)
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    type(CellWithNbrsType), intent(inout) :: cellNbrs !< cell to add to
    type(GlobalCellType), optional :: mask !< a mask to exclude back-and-forth connections
    ! local
    integer(I4B) :: ix, iexg
    class(VirtualExchangeType), pointer :: v_exchange
    class(VirtualModelType), pointer :: v_m1, v_m2

    ! loop over all exchanges
    do ix = 1, this%haloExchanges%size

      v_exchange => get_virtual_exchange(this%haloExchanges%at(ix))
      v_m1 => v_exchange%v_model1
      v_m2 => v_exchange%v_model2

      ! loop over n-m links in the exchange
      if (cellNbrs%cell%v_model == v_m1) then
        do iexg = 1, v_exchange%nexg%get()
          if (v_exchange%nodem1%get(iexg) == cellNbrs%cell%index) then
            ! we have a link, now add foreign neighbor
            call this%addNeighborCell( &
              cellNbrs, v_exchange%nodem2%get(iexg), v_m2, mask)
          end if
        end do
      end if
      ! and the reverse
      if (cellNbrs%cell%v_model == v_m2) then
        do iexg = 1, v_exchange%nexg%get()
          if (v_exchange%nodem2%get(iexg) == cellNbrs%cell%index) then
            ! we have a link, now add foreign neighbor
            call this%addNeighborCell( &
              cellNbrs, v_exchange%nodem1%get(iexg), v_m1, mask)
          end if
        end do
      end if

    end do

  end subroutine addRemoteNeighbors

  !> @brief Add neighboring cell to tree structure
  !<
  subroutine addNeighborCell(this, cellNbrs, newNbrIdx, v_nbr_model, mask)
    class(GridConnectionType), intent(in) :: this !< this grid connection instance
    type(CellWithNbrsType), intent(inout) :: cellNbrs !< the root cell which to add to
    integer(I4B), intent(in) :: newNbrIdx !< the neighboring cell's index
    class(VirtualModelType), pointer :: v_nbr_model !< the model where the new neighbor lives
    type(GlobalCellType), optional :: mask !< don't add connections to this cell (optional)

    if (present(mask)) then
      if (newNbrIdx == mask%index .and. mask%v_model == v_nbr_model) then
        return
      end if
    end if

    call cellNbrs%addNbrCell(newNbrIdx, v_nbr_model)

  end subroutine addNeighborCell

  !> @brief Recursively set interface cell indexes and
  !< add to the region-to-interface lookup table
  recursive subroutine registerInterfaceCells(this, cellWithNbrs)
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    type(CellWithNbrsType) :: cellWithNbrs !< the cell from where to start registering neighbors
    ! local
    integer(I4B) :: offset, inbr
    integer(I4B) :: regionIdx ! unique idx in the region (all connected models)
    integer(I4B) :: ifaceIdx ! unique idx in the interface grid

    offset = this%get_regional_offset(cellWithNbrs%cell%v_model)
    regionIdx = offset + cellWithNbrs%cell%index
    ifaceIdx = this%getInterfaceIndex(cellWithNbrs%cell)
    if (ifaceIdx == -1) then
      this%indexCount = this%indexCount + 1
      ifaceIdx = this%indexCount
      call this%addToGlobalMap(ifaceIdx, cellWithNbrs%cell)
      this%region_to_iface_map(regionIdx) = ifaceIdx
    end if

    ! and also for its neighbors
    do inbr = 1, cellWithNbrs%nrOfNbrs
      call this%registerInterfaceCells(cellWithNbrs%neighbors(inbr))
    end do

  end subroutine registerInterfaceCells

  !> @brief Add entry to lookup table, inflating when necessary
  !<
  subroutine addToGlobalMap(this, ifaceIdx, cell)
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    integer(I4B), intent(in) :: ifaceIdx !< unique idx in the interface grid
    type(GlobalCellType), intent(in) :: cell !< the global cell
    ! local
    integer(I4B) :: i, currentSize, newSize
    type(GlobalCellType), dimension(:), pointer :: tempMap

    ! inflate?
    currentSize = size(this%idxToGlobal)
    if (ifaceIdx > currentSize) then
      newSize = nint(1.5 * currentSize)
      allocate (tempMap(newSize))
      do i = 1, currentSize
        tempMap(i) = this%idxToGlobal(i)
      end do

      deallocate (this%idxToGlobal)
      this%idxToGlobal => tempMap
    end if

    this%idxToGlobal(ifaceIdx) = cell

  end subroutine addToGlobalMap

  !> @brief Compress lookup table to get rid of unused entries
  !<
  subroutine compressGlobalMap(this)
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    ! local
    type(GlobalCellType), dimension(:), pointer :: tempMap

    if (size(this%idxToGlobal) > this%nrOfCells) then
      allocate (tempMap(this%nrOfCells))
      tempMap(1:this%nrOfCells) = this%idxToGlobal(1:this%nrOfCells)
      deallocate (this%idxToGlobal)
      allocate (this%idxToGlobal(this%nrOfCells))
      this%idxToGlobal(1:this%nrOfCells) = tempMap(1:this%nrOfCells)
      deallocate (tempMap)
    end if

  end subroutine compressGlobalMap

  !> @brief Soft cell ids in the interface grid such that
  !< id_1 < id_2 means that cell 1 lies above cell 2
  subroutine sortInterfaceGrid(this)
    use GridSorting, only: quickSortGrid
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    ! local
    integer(I4B), dimension(:), allocatable :: newToOldIdx
    integer(I4B), dimension(:), allocatable :: oldToNewIdx
    integer(I4B) :: idxOld
    integer(I4B) :: i
    type(GlobalCellType), dimension(:), allocatable :: sortedGlobalMap
    integer(I4B), dimension(:), allocatable :: sortedRegionMap
    logical(LGP) :: z_only

    ! only use z coordinate for sorting when no cell2d info is available
    z_only = (this%icondir == 0)

    ! sort based on coordinates
    newToOldIdx = (/(i, i=1, size(this%idxToGlobal))/)
    call quickSortGrid(newToOldIdx, size(newToOldIdx), this%idxToGlobal, z_only)

    ! and invert
    allocate (oldToNewIdx(size(newToOldIdx)))
    do i = 1, size(oldToNewIdx)
      oldToNewIdx(newToOldIdx(i)) = i
    end do

    ! reorder global table
    allocate (sortedGlobalMap(size(this%idxToGlobal)))
    do i = 1, size(newToOldIdx)
      sortedGlobalMap(i) = this%idxToGlobal(newToOldIdx(i))
    end do
    do i = 1, size(newToOldIdx)
      this%idxToGlobal(i) = sortedGlobalMap(i)
    end do
    deallocate (sortedGlobalMap)

    ! reorder regional lookup table
    allocate (sortedRegionMap(size(this%region_to_iface_map)))
    do i = 1, size(sortedRegionMap)
      if (this%region_to_iface_map(i) /= -1) then
        idxOld = this%region_to_iface_map(i)
        sortedRegionMap(i) = oldToNewIdx(idxOld)
      else
        sortedRegionMap(i) = -1
      end if
    end do
    do i = 1, size(sortedRegionMap)
      this%region_to_iface_map(i) = sortedRegionMap(i)
    end do
    deallocate (sortedRegionMap)

  end subroutine sortInterfaceGrid

  !> @brief Add primary connections to the sparse data structure
  !<
  subroutine makePrimaryConnections(this, sparse)
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    type(SparseMatrix), pointer :: sparse !< the sparse data structure to hold the connections
    ! local
    integer(I4B) :: icell
    integer(I4B) :: ifaceIdx, ifaceIdxNbr

    do icell = 1, this%nrOfBoundaryCells
      ifaceIdx = this%getInterfaceIndex(this%boundaryCells(icell)%cell)
      ifaceIdxNbr = this%getInterfaceIndex(this%connectedCells(icell)%cell)

      ! add diagonals to sparse
      call sparse%addconnection(ifaceIdx, ifaceIdx, 1)
      call sparse%addconnection(ifaceIdxNbr, ifaceIdxNbr, 1)

      ! and cross terms
      call sparse%addconnection(ifaceIdx, ifaceIdxNbr, 1)
      call sparse%addconnection(ifaceIdxNbr, ifaceIdx, 1)
    end do

  end subroutine makePrimaryConnections

  !> @brief Recursively add higher order connections (from
  !! cells neighboring the primarily connected cells) to the
  !< sparse data structure
  recursive subroutine connectNeighborCells(this, cell, sparse)
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    type(CellWithNbrsType) :: cell !< the cell whose connections is to be added
    type(SparseMatrix), pointer :: sparse !< the sparse data structure to hold the connections
    ! local
    integer(I4B) :: ifaceIdx, ifaceIdxNbr ! unique idx in the interface grid
    integer(I4B) :: inbr

    ifaceIdx = this%getInterfaceIndex(cell%cell)
    do inbr = 1, cell%nrOfNbrs
      ifaceIdxNbr = this%getInterfaceIndex(cell%neighbors(inbr)%cell)

      call sparse%addconnection(ifaceIdxNbr, ifaceIdxNbr, 1)
      call sparse%addconnection(ifaceIdx, ifaceIdxNbr, 1)
      call sparse%addconnection(ifaceIdxNbr, ifaceIdx, 1)

      ! recurse
      call this%connectNeighborCells(cell%neighbors(inbr), sparse)
    end do

  end subroutine connectNeighborCells

  !> @brief Fill connection data (ihc, cl1, ...) for
  !< connections between cells within the same model.
  subroutine fillConnectionDataInternal(this)
    use ConstantsModule, only: DPI, DTWOPI
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    ! local
    type(ConnectionsType), pointer :: conn
    integer(I4B) :: n, m, ipos, isym, ipos_orig, isym_orig
    type(GlobalCellType), pointer :: ncell, mcell
    class(VirtualModelType), pointer :: v_m !< pointer to virtual model (readability)

    conn => this%connections

    do n = 1, conn%nodes
      do ipos = conn%ia(n) + 1, conn%ia(n + 1) - 1
        m = conn%ja(ipos)
        if (n > m) cycle

        isym = conn%jas(ipos)
        ncell => this%idxToGlobal(n)
        mcell => this%idxToGlobal(m)
        if (ncell%v_model == mcell%v_model) then

          ! for readability
          v_m => ncell%v_model

          ! within same model, straight copy
          ipos_orig = getCSRIndex(ncell%index, mcell%index, &
                                  v_m%con_ia%get_array(), v_m%con_ja%get_array())

          if (ipos_orig == 0) then
            ! periodic boundary conditions can add connections between cells in
            ! the same model, but they are dealt with through the exchange data
            if (this%isPeriodic(ncell%index, mcell%index)) cycle

            ! this should not be possible
            write (*, *) 'Error: cannot find cell connection in model grid'
            call ustop()
          end if

          isym_orig = v_m%con_jas%get(ipos_orig)
          conn%hwva(isym) = v_m%con_hwva%get(isym_orig)
          conn%ihc(isym) = v_m%con_ihc%get(isym_orig)
          if (ncell%index < mcell%index) then
            conn%cl1(isym) = v_m%con_cl1%get(isym_orig)
            conn%cl2(isym) = v_m%con_cl2%get(isym_orig)
            conn%anglex(isym) = v_m%con_anglex%get(isym_orig)
          else
            conn%cl1(isym) = v_m%con_cl2%get(isym_orig)
            conn%cl2(isym) = v_m%con_cl1%get(isym_orig)
            conn%anglex(isym) = mod(v_m%con_anglex%get(isym_orig) + DPI, DTWOPI)
          end if
        end if
      end do
    end do

  end subroutine fillConnectionDataInternal

  !> @brief Fill connection data (ihc, cl1, ...) for
  !< all exchanges
  subroutine fillConnectionDataFromExchanges(this)
    use ConstantsModule, only: DPIO180
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    ! local
    integer(I4B) :: inx, iexg
    integer(I4B) :: ipos, isym
    integer(I4B) :: nOffset, mOffset, nIfaceIdx, mIfaceIdx
    class(VirtualExchangeType), pointer :: v_exg
    class(VirtualModelType), pointer :: v_m1, v_m2
    type(ConnectionsType), pointer :: conn

    conn => this%connections

    do inx = 1, this%haloExchanges%size
      v_exg => get_virtual_exchange(this%haloExchanges%at(inx))

      v_m1 => v_exg%v_model1
      v_m2 => v_exg%v_model2

      if (v_exg%ianglex%get() > 0) then
        conn%ianglex = 1
      end if

      nOffset = this%get_regional_offset(v_m1)
      mOffset = this%get_regional_offset(v_m2)
      do iexg = 1, v_exg%nexg%get()
        nIfaceIdx = this%region_to_iface_map(noffset + v_exg%nodem1%get(iexg))
        mIfaceIdx = this%region_to_iface_map(moffset + v_exg%nodem2%get(iexg))
        ! not all nodes from the exchanges are part of the interface grid
        ! (think of exchanges between neighboring models, and their neighbors)
        if (nIFaceIdx == -1 .or. mIFaceIdx == -1) then
          cycle
        end if

        ipos = conn%getjaindex(nIfaceIdx, mIfaceIdx)
        ! (see prev. remark) sometimes the cells are in the interface grid,
        ! but the connection isn't. This can happen for leaf nodes of the grid.
        if (ipos == 0) then
          ! no match, safely cycle
          cycle
        end if
        isym = conn%jas(ipos)

        ! note: cl1 equals L_nm: the length from cell n to the shared
        ! face with cell m (and cl2 analogously for L_mn)
        if (nIfaceIdx < mIfaceIdx) then
          conn%cl1(isym) = v_exg%cl1%get(iexg)
          conn%cl2(isym) = v_exg%cl2%get(iexg)
          if (v_exg%ianglex%get() > 0) then
            conn%anglex(isym) = &
              v_exg%auxvar%get(v_exg%ianglex%get(), iexg) * DPIO180
          end if
        else
          conn%cl1(isym) = v_exg%cl2%get(iexg)
          conn%cl2(isym) = v_exg%cl1%get(iexg)
          if (v_exg%ianglex%get() > 0) then
            conn%anglex(isym) = mod(v_exg%auxvar%get(v_exg%ianglex%get(), iexg) &
                                    + 180.0_DP, 360.0_DP) * DPIO180
          end if
        end if
        conn%hwva(isym) = v_exg%hwva%get(iexg)
        conn%ihc(isym) = v_exg%ihc%get(iexg)

      end do
    end do

  end subroutine fillConnectionDataFromExchanges

  !> @brief Create the connection masks
  !!
  !! The level indicates the nr of connections away from
  !! the remote neighbor, the diagonal term holds the negated
  !! value of their nearest connection. We end with setting
  !< a normalized mask to the connections object: 0 or 1
  subroutine createConnectionMask(this)
    class(GridConnectionType), intent(inout) :: this !< instance of this grid connection
    ! local
    integer(I4B) :: icell, inbr, n, ipos
    integer(I4B) :: level, newMask
    type(CellWithNbrsType), pointer :: cell, nbrCell

    ! set all masks to zero to begin with
    do ipos = 1, this%connections%nja
      call this%connections%set_mask(ipos, 0)
    end do

    ! remote connections remain masked
    ! now set mask for exchange connections (level == 1)
    level = 1
    do icell = 1, this%nrOfBoundaryCells
      call this%setMaskOnConnection(this%boundaryCells(icell), &
                                    this%connectedCells(icell), level)
      ! for cross-boundary connections, we need to apply the mask to both n-m and m-n,
      ! because if the upper triangular one is disabled, its transposed (lower triangular)
      ! counter part is skipped in the NPF calculation as well.
      call this%setMaskOnConnection(this%connectedCells(icell), &
                                    this%boundaryCells(icell), level)
    end do

    ! now extend mask recursively into the internal domain (level > 1)
    do icell = 1, this%nrOfBoundaryCells
      cell => this%boundaryCells(icell)
      do inbr = 1, cell%nrOfNbrs
        nbrCell => this%boundaryCells(icell)%neighbors(inbr)
        level = 2 ! this is incremented within the recursion
        call this%maskInternalConnections(this%boundaryCells(icell), &
                                          this%boundaryCells(icell)% &
                                          neighbors(inbr), level)
      end do
    end do

    ! set normalized mask:
    ! =1 for links with connectivity <= interior stencil depth
    ! =0 otherwise
    do n = 1, this%connections%nodes
      ! set diagonals to zero
      call this%connections%set_mask(this%connections%ia(n), 0)

      do ipos = this%connections%ia(n) + 1, this%connections%ia(n + 1) - 1
        newMask = 0
        if (this%connections%mask(ipos) > 0) then
          if (this%connections%mask(ipos) < this%internalStencilDepth + 1) then
            newMask = 1
          end if
        end if
        ! set mask on off-diag
        call this%connections%set_mask(ipos, newMask)
      end do
    end do

  end subroutine createConnectionMask

  !> @brief Recursively mask connections, increasing the level as we go
  !<
  recursive subroutine maskInternalConnections(this, cell, nbrCell, level)
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    type(CellWithNbrsType), intent(inout) :: cell !< cell 1 in the connection to mask
    type(CellWithNbrsType), intent(inout) :: nbrCell !< cell 2 in the connection to mask
    integer(I4B), intent(in) :: level
    ! local
    integer(I4B) :: inbr, newLevel

    ! only set the mask for internal connections, leaving the
    ! others at 0
    if (cell%cell%v_model == this%model .and. &
        nbrCell%cell%v_model == this%model) then
      ! this will set a mask on both diagonal, and both cross terms
      call this%setMaskOnConnection(cell, nbrCell, level)
      call this%setMaskOnConnection(nbrCell, cell, level)
    end if

    ! recurse on nbrs-of-nbrs
    newLevel = level + 1
    do inbr = 1, nbrCell%nrOfNbrs
      call this%maskInternalConnections(nbrCell, &
                                        nbrCell%neighbors(inbr), &
                                        newLevel)
    end do

  end subroutine maskInternalConnections

  !> @brief Set a mask on the connection from a cell to its neighbor,
  !! (and not the transposed!) not overwriting the current level
  !< of a connection when it is smaller
  subroutine setMaskOnConnection(this, cell, nbrCell, level)
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    type(CellWithNbrsType), intent(inout) :: cell !< cell 1 in the connection
    type(CellWithNbrsType), intent(inout) :: nbrCell !< cell 2 in the connection
    integer(I4B), intent(in) :: level !< the level value to set the mask to
    ! local
    integer(I4B) :: ifaceIdx, ifaceIdxNbr
    integer(I4B) :: iposdiag, ipos
    integer(I4B) :: currentLevel

    ifaceIdx = this%getInterfaceIndex(cell%cell)
    ifaceIdxNbr = this%getInterfaceIndex(nbrCell%cell)

    ! diagonal
    iposdiag = this%connections%getjaindex(ifaceIdx, ifaceIdx)
    currentLevel = this%connections%mask(iposdiag)
    if (currentLevel == 0 .or. level < currentLevel) then
      call this%connections%set_mask(iposdiag, level)
    end if
    ! cross term
    ipos = this%connections%getjaindex(ifaceIdx, ifaceIdxNbr)
    currentLevel = this%connections%mask(ipos)
    if (currentLevel == 0 .or. level < currentLevel) then
      call this%connections%set_mask(ipos, level)
    end if

  end subroutine setMaskOnConnection

  !> @brief Get interface index from global cell
  !<
  function getInterfaceIndexByCell(this, cell) result(iface_idx)
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    type(GlobalCellType), intent(in) :: cell !< the global cell to get the interface index for
    integer(I4B) :: iface_idx !< the index in the interface model
    ! local
    integer(I4B) :: offset, region_idx

    offset = this%get_regional_offset(cell%v_model)
    region_idx = offset + cell%index
    iface_idx = this%region_to_iface_map(region_idx)

  end function getInterfaceIndexByCell

  !> @brief Get interface index from a model pointer and the local index
  !<
  function getInterfaceIndexByIndexModel(this, index, v_model) result(iface_idx)
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    integer(I4B) :: index !< the local cell index
    class(VirtualModelType), pointer :: v_model !< the cell's model
    integer(I4B) :: iface_idx !< the index in the interface model
    ! local
    integer(I4B) :: offset, region_idx

    offset = this%get_regional_offset(v_model)
    region_idx = offset + index
    iface_idx = this%region_to_iface_map(region_idx)

  end function getInterfaceIndexByIndexModel

  !> @brief Get the offset for a regional model
  !<
  function get_regional_offset(this, v_model) result(offset)
    class(GridConnectionType), intent(inout) :: this !< this grid connection instance
    class(VirtualModelType), pointer :: v_model !< the model to get the offset for
    integer(I4B) :: offset !< the index offset in the regional domain
    ! local
    integer(I4B) :: im
    class(VirtualModelType), pointer :: vm

    offset = 0
    do im = 1, this%regionalModels%Count()
      vm => get_virtual_model_from_list(this%regionalModels, im)
      if (vm == v_model) then
        offset = this%regionalModelOffset(im)
        return
      end if
    end do

  end function get_regional_offset

  !> @brief Allocate scalar data
  !<
  subroutine allocateScalars(this)
    use MemoryManagerModule, only: mem_allocate
    class(GridConnectionType) :: this !< this grid connection instance

    call mem_allocate(this%nrOfBoundaryCells, 'NRBNDCELLS', this%memoryPath)
    call mem_allocate(this%indexCount, 'IDXCOUNT', this%memoryPath)
    call mem_allocate(this%nrOfCells, 'NRCELLS', this%memoryPath)

  end subroutine allocateScalars

  !> @brief Sets the discretization (DISU) after all
  !! preprocessing by this grid connection has been done,
  !< this comes after disu_cr
  subroutine getDiscretization(this, disu)
    use ConnectionsModule
    use SparseModule, only: sparsematrix
    class(GridConnectionType) :: this !< the grid connection
    class(DisuType), pointer :: disu !< the target disu object
    ! local
    integer(I4B) :: icell, nrOfCells, idx
    class(VirtualModelType), pointer :: v_model
    real(DP) :: xglo, yglo

    ! the following is similar to dis_df
    nrOfCells = this%nrOfCells
    disu%nodes = nrOfCells
    disu%nvert = 0 ! TODO_MJR: vertex reconstruction is not supported yet
    disu%icondir = this%icondir
    disu%nodesuser = nrOfCells
    disu%nja = this%connections%nja

    call disu%allocate_arrays()

    ! these are otherwise allocated when sourced
    call disu%allocate_arrays_mem()

    ! fill data
    do icell = 1, nrOfCells
      idx = this%idxToGlobal(icell)%index
      disu%top(icell) = -huge(1.0_DP)
      disu%bot(icell) = -huge(1.0_DP)
      disu%area(icell) = -huge(1.0_DP)
    end do

    ! grid connections follow from GridConnection:
    disu%con => this%connections
    disu%njas = disu%con%njas

    if (this%icondir > 0) then
      ! copy cell x,y because available in all models
      do icell = 1, nrOfCells
        idx = this%idxToGlobal(icell)%index
        v_model => this%idxToGlobal(icell)%v_model

        ! we are merging grids with possibly (likely) different origins,
        ! transform to global coordinates:
        call dis_transform_xy(v_model%dis_xc%get(idx), &
                              v_model%dis_yc%get(idx), &
                              v_model%dis_xorigin%get(), &
                              v_model%dis_yorigin%get(), &
                              v_model%dis_angrot%get(), &
                              xglo, yglo)

        ! NB: usernodes equals internal nodes for interface
        disu%cellxy(1, icell) = xglo
        disu%xc(icell) = xglo
        disu%cellxy(2, icell) = yglo
        disu%yc(icell) = yglo
      end do
    else
      ! otherwise compress
      call mem_reallocate(disu%xc, 0, 'XC', disu%memoryPath)
      call mem_reallocate(disu%yc, 0, 'YC', disu%memoryPath)
    end if

    ! if vertices will be needed, it will look like this:
    !
    ! 1. determine total nr. of verts
    ! 2. allocate vertices list
    ! 3. create sparse
    ! 4. get vertex data per cell, add functions to base
    ! 5. add vertex (x,y) to list and connectivity to sparse
    ! 6. generate ia/ja from sparse

  end subroutine getDiscretization

  !> @brief Build interface map object for outside use
  subroutine buildInterfaceMap(this)
    use BaseModelModule, only: BaseModelType, GetBaseModelFromList
    use STLVecIntModule
    class(GridConnectionType) :: this !< this grid connection
    ! local
    integer(I4B) :: i, j, iloc, jloc
    integer(I4B) :: im, ix, mid, n
    integer(I4B) :: ipos, ipos_model
    type(STLVecInt) :: model_ids
    type(STLVecInt) :: src_idx_tmp, tgt_idx_tmp, sign_tmp
    class(VirtualExchangeType), pointer :: v_exg
    class(VirtualModelType), pointer :: vm
    class(VirtualModelType), pointer :: v_model1, v_model2
    type(InterfaceMapType), pointer :: imap
    integer(I4B), dimension(:), pointer, contiguous :: ia_ptr, ja_ptr

    allocate (this%interfaceMap)
    imap => this%interfaceMap

    ! first get the participating models
    call model_ids%init()
    do i = 1, this%nrOfCells
      call model_ids%push_back_unique(this%idxToGlobal(i)%v_model%id)
    end do

    ! initialize the map
    call imap%init(model_ids%size, this%haloExchanges%size)

    ! for each model part of this interface, ...
    do im = 1, model_ids%size
      mid = model_ids%at(im)
      imap%model_ids(im) = mid
      vm => get_virtual_model(mid)
      imap%model_names(im) = vm%name
      call src_idx_tmp%init()
      call tgt_idx_tmp%init()

      ! store the node map for this model
      do i = 1, this%nrOfCells
        if (mid == this%idxToGlobal(i)%v_model%id) then
          call src_idx_tmp%push_back(this%idxToGlobal(i)%index)
          call tgt_idx_tmp%push_back(i)
        end if
      end do

      ! and copy into interface map
      allocate (imap%node_maps(im)%src_idx(src_idx_tmp%size))
      allocate (imap%node_maps(im)%tgt_idx(tgt_idx_tmp%size))
      do i = 1, src_idx_tmp%size
        imap%node_maps(im)%src_idx(i) = src_idx_tmp%at(i)
        imap%node_maps(im)%tgt_idx(i) = tgt_idx_tmp%at(i)
      end do

      call src_idx_tmp%destroy()
      call tgt_idx_tmp%destroy()

      ! and for connections
      call src_idx_tmp%init()
      call tgt_idx_tmp%init()

      ! store the connection map for this model
      do i = 1, this%nrOfCells
        if (mid /= this%idxToGlobal(i)%v_model%id) cycle
        do ipos = this%connections%ia(i), this%connections%ia(i + 1) - 1
          j = this%connections%ja(ipos)
          if (mid /= this%idxToGlobal(j)%v_model%id) cycle

          ! i and j are now in same model (mid)
          iloc = this%idxToGlobal(i)%index
          jloc = this%idxToGlobal(j)%index
          ia_ptr => this%idxToGlobal(i)%v_model%con_ia%get_array()
          ja_ptr => this%idxToGlobal(i)%v_model%con_ja%get_array()
          ipos_model = getCSRIndex(iloc, jloc, ia_ptr, ja_ptr)
          call src_idx_tmp%push_back(ipos_model)
          call tgt_idx_tmp%push_back(ipos)
        end do
      end do

      ! copy into interface map
      allocate (imap%conn_maps(im)%src_idx(src_idx_tmp%size))
      allocate (imap%conn_maps(im)%tgt_idx(tgt_idx_tmp%size))
      do i = 1, src_idx_tmp%size
        imap%conn_maps(im)%src_idx(i) = src_idx_tmp%at(i)
        imap%conn_maps(im)%tgt_idx(i) = tgt_idx_tmp%at(i)
      end do

      call src_idx_tmp%destroy()
      call tgt_idx_tmp%destroy()

    end do

    call model_ids%destroy()

    ! for each exchange that is part of this interface
    do ix = 1, this%haloExchanges%size

      ! all exchanges in this list should have at
      ! least one relevant connection for this map
      v_exg => get_virtual_exchange(this%haloExchanges%at(ix))
      v_model1 => v_exg%v_model1
      v_model2 => v_exg%v_model2

      imap%exchange_ids(ix) = v_exg%id
      imap%exchange_names(ix) = v_exg%name

      call src_idx_tmp%init()
      call tgt_idx_tmp%init()
      call sign_tmp%init()

      do n = 1, v_exg%nexg%get()
        i = this%getInterfaceIndex(v_exg%nodem1%get(n), v_model1)
        j = this%getInterfaceIndex(v_exg%nodem2%get(n), v_model2)
        if (i == -1 .or. j == -1) cycle ! not all exchange nodes are part of the interface
        ipos = this%connections%getjaindex(i, j)
        if (ipos == 0) then
          ! this can typically happen at corner points for a larger
          ! stencil (XT3D), when both i and j are included in the
          ! interface grid as leaf nodes, but their connection is not.
          ! (c.f. 'test_gwf_ifmod_mult_exg.py')
          cycle
        end if
        call src_idx_tmp%push_back(n)
        call tgt_idx_tmp%push_back(ipos)
        call sign_tmp%push_back(1)

        ! and the reverse connection:
        call src_idx_tmp%push_back(n)
        call tgt_idx_tmp%push_back(this%connections%isym(ipos))
        call sign_tmp%push_back(-1)
      end do

      allocate (imap%exchange_maps(ix)%src_idx(src_idx_tmp%size))
      allocate (imap%exchange_maps(ix)%tgt_idx(tgt_idx_tmp%size))
      allocate (imap%exchange_maps(ix)%sign(sign_tmp%size))
      do i = 1, src_idx_tmp%size
        imap%exchange_maps(ix)%src_idx(i) = src_idx_tmp%at(i)
        imap%exchange_maps(ix)%tgt_idx(i) = tgt_idx_tmp%at(i)
        imap%exchange_maps(ix)%sign(i) = sign_tmp%at(i)
      end do

      call src_idx_tmp%destroy()
      call tgt_idx_tmp%destroy()
      call sign_tmp%destroy()

    end do

    ! set the primary exchange idx
    ! findloc cannot be used until gfortran 9...
    imap%prim_exg_idx = -1
    do i = 1, imap%nr_exchanges
      if (imap%exchange_names(i) == this%primaryExchange%name) then
        imap%prim_exg_idx = i
        exit
      end if
    end do

    ! sanity check
    ! do i = 1, interfaceMap%nr_models
    !   if (size(interfaceMap%node_map(i)%src_idx) == 0) then
    !     write(*,*) 'Error: empty node map in interface for ', &
    !                this%primaryExchange%name
    !     call ustop()
    !   end if
    !   if (size(interfaceMap%connection_map(i)%src_idx) == 0) then
    !     write(*,*) 'Error: empty connection map in interface for ', &
    !                this%primaryExchange%name
    !     call ustop()
    !   end if
    ! end do
    ! do i = 1, interfaceMap%nr_exchanges
    !   if (size(interfaceMap%exchange_map(i)%src_idx) == 0) then
    !     write(*,*) 'Error: empty exchange map in interface for ', &
    !                this%primaryExchange%name
    !     call ustop()
    !   end if
    ! end do

  end subroutine buildInterfaceMap

  !> @brief Deallocate grid connection resources
  !<
  subroutine destroy(this)
    use MemoryManagerModule, only: mem_deallocate
    class(GridConnectionType) :: this !< this grid connection instance

    call mem_deallocate(this%nrOfBoundaryCells)
    call mem_deallocate(this%indexCount)
    call mem_deallocate(this%nrOfCells)

    ! arrays
    deallocate (this%idxToGlobal)
    deallocate (this%boundaryCells)
    deallocate (this%connectedCells)

    call mem_deallocate(this%idxToGlobalIdx)

  end subroutine destroy

  !> @brief Test if the connection between nodes within
  !< the same model is periodic
  function isPeriodic(this, n, m) result(periodic)
    class(GridConnectionType), intent(in) :: this !< this grid connection instance
    integer(I4B), intent(in) :: n !< first node of the connection
    integer(I4B), intent(in) :: m !< second node of the connection
    logical :: periodic !< true when periodic
    ! local
    integer(I4B) :: icell

    periodic = .false.
    do icell = 1, this%nrOfBoundaryCells
      if (.not. this%boundaryCells(icell)%cell%v_model == &
          this%connectedCells(icell)%cell%v_model) then
        cycle
      end if

      ! one way
      if (this%boundaryCells(icell)%cell%index == n) then
        if (this%connectedCells(icell)%cell%index == m) then
          periodic = .true.
          return
        end if
      end if
      ! or the other
      if (this%boundaryCells(icell)%cell%index == m) then
        if (this%connectedCells(icell)%cell%index == n) then
          periodic = .true.
          return
        end if
      end if

    end do

  end function

end module GridConnectionModule
