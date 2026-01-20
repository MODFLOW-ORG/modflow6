module SpatialModelConnectionModule
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LINELENGTH, DZERO
  use SparseModule, only: sparsematrix
  use ConnectionsModule, only: ConnectionsType
  use CsrUtilsModule, only: getCSRIndex
  use SimModule, only: ustop, count_errors, store_error
  use NumericalModelModule, only: NumericalModelType
  use NumericalExchangeModule, only: NumericalExchangeType
  use DisConnExchangeModule, only: DisConnExchangeType
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, mem_checkin
  use MemoryHelperModule, only: create_mem_path
  use GridConnectionModule, only: GridConnectionType
  use InterfaceMapModule
  use DistVariableModule
  use VirtualDataListsModule, only: virtual_exchange_list
  use VirtualModelModule, only: VirtualModelType, get_virtual_model
  use VirtualExchangeModule, only: VirtualExchangeType, &
                                   get_virtual_exchange_from_list
  use ListModule, only: ListType
  use STLVecIntModule, only: STLVecInt
  use MatrixBaseModule
  use SparseMatrixModule

  implicit none
  private
  public :: cast_as_smc
  public :: add_smc_to_list
  public :: get_smc_from_list

  !> Class to manage spatial connection of a model to one
  !! or more models of the same type. Spatial connection here
  !! means that the model domains (spatial discretization) are
  !! adjacent and connected via DisConnExchangeType object(s).
  !! The connection itself is a Numerical Exchange as well,
  !! and part of a Numerical Solution providing the amat and rhs
  !< values for the exchange.
  type, public, extends(NumericalExchangeType) :: SpatialModelConnectionType

    class(NumericalModelType), pointer :: owner => null() !< the model whose connection this is
    class(NumericalModelType), pointer :: interface_model => null() !< the interface model
    integer(I4B), pointer :: nr_connections => null() !< total nr. of connected cells (primary)

    class(DisConnExchangeType), pointer :: prim_exchange => null() !< the exchange for which the interface model is created
    logical(LGP) :: owns_exchange !< there are two connections (in serial) for an exchange,
                                  !! one of them needs to manage/own the exchange (e.g. clean up)
    type(STLVecInt), pointer :: halo_models !< models that are potentially in the halo of this interface
    type(STLVecInt), pointer :: halo_exchanges !< exchanges that are potentially part of the halo of this interface (includes primary)
    integer(I4B), pointer :: int_stencil_depth => null() !< size of the computational stencil for the interior
                                                         !! default = 1, xt3d = 2, ...
    integer(I4B), pointer :: exg_stencil_depth => null() !< size of the computational stencil at the interface
                                                         !! default = 1, xt3d = 2, ...

    ! The following variables are equivalent to those in Numerical Solution:
    integer(I4B), pointer :: neq => null() !< nr. of equations in matrix system
    class(SparseMatrixType), pointer :: matrix => null() !< system matrix for the interface
    real(DP), dimension(:), pointer, contiguous :: rhs => null() !< rhs of interface system
    real(DP), dimension(:), pointer, contiguous :: x => null() !< dependent variable of interface system
    integer(I4B), dimension(:), pointer, contiguous :: active => null() !< cell status (c.f. ibound) of interface system

    ! these are not in the memory manager
    class(GridConnectionType), pointer :: ig_builder => null() !< facility to build the interface grid connection structure
    integer(I4B), dimension(:), pointer :: ipos_to_sln => null() !< mapping between position in the interface matrix and the solution matrix
    type(ListType) :: iface_dist_vars !< list with distributed variables for this interface
    type(InterfaceMapType), pointer :: interface_map => null() !< a map of the interface into models and exchanges

  contains

    ! public
    procedure, pass(this) :: spatialConnection_ctor
    generic :: construct => spatialConnection_ctor

    ! partly overriding NumericalExchangeType:
    procedure :: exg_df => spatialcon_df
    procedure :: exg_ar => spatialcon_ar
    procedure :: exg_ac => spatialcon_ac
    procedure :: exg_mc => spatialcon_mc
    procedure :: exg_cf => spatialcon_cf
    procedure :: exg_fc => spatialcon_fc
    procedure :: exg_da => spatialcon_da

    ! protected
    procedure, pass(this) :: spatialcon_df
    procedure, pass(this) :: spatialcon_ar
    procedure, pass(this) :: spatialcon_ac
    procedure, pass(this) :: spatialcon_cf
    procedure, pass(this) :: spatialcon_fc
    procedure, pass(this) :: spatialcon_da
    procedure, pass(this) :: spatialcon_setmodelptrs
    procedure, pass(this) :: spatialcon_connect
    procedure, pass(this) :: validateConnection
    procedure, pass(this) :: cfg_dv
    procedure, pass(this) :: createModelHalo

    ! private
    procedure, private, pass(this) :: setupGridConnection
    procedure, private, pass(this) :: getNrOfConnections
    procedure, private, pass(this) :: allocateScalars
    procedure, private, pass(this) :: allocateArrays
    procedure, private, pass(this) :: createCoefficientMatrix
    procedure, private, pass(this) :: maskOwnerConnections
    procedure, private, pass(this) :: addModelNeighbors

  end type SpatialModelConnectionType

contains ! module procedures

  !> @brief Construct the spatial connection base
  !!
  !! This constructor is typically called from a derived class.
  !<
  subroutine spatialConnection_ctor(this, model, exchange, name)
    class(SpatialModelConnectionType) :: this !< the connection
    class(NumericalModelType), intent(in), pointer :: model !< the model that owns the connection
    class(DisConnExchangeType), intent(in), pointer :: exchange !< the primary exchange from which
                                                                !! the connection is created
    character(len=*), intent(in) :: name !< the connection name (for memory management mostly)

    this%name = name
    this%memoryPath = create_mem_path(this%name)

    this%owner => model
    this%prim_exchange => exchange

    allocate (this%ig_builder)
    allocate (this%halo_models)
    allocate (this%halo_exchanges)
    allocate (this%matrix)
    call this%allocateScalars()

    this%int_stencil_depth = 1
    this%exg_stencil_depth = 1
    this%nr_connections = 0

    ! this should be set in derived ctor
    this%interface_model => null()

  end subroutine spatialConnection_ctor

  !> @brief Find all models that might participate in this interface
  !<
  subroutine createModelHalo(this)
    class(SpatialModelConnectionType) :: this !< this connection

    call this%halo_models%init()
    call this%halo_exchanges%init()

    call this%addModelNeighbors(this%owner%id, virtual_exchange_list, &
                                this%exg_stencil_depth, .true.)

  end subroutine createModelHalo

  !> @brief Add neighbors and nbrs-of-nbrs to the model tree
  !<
  recursive subroutine addModelNeighbors(this, model_id, &
                                         virtual_exchanges, &
                                         depth, is_root, mask)
    use VirtualExchangeModule, only: get_virtual_exchange
    class(SpatialModelConnectionType) :: this !< this connection
    integer(I4B) :: model_id !< the model (id) to add neighbors for
    type(ListType) :: virtual_exchanges !< list with all virtual exchanges
    integer(I4B), value :: depth !< the maximal number of exchanges between
    logical(LGP) :: is_root !< true when called for neighbor from primary exchange
    integer(I4B), optional :: mask !< don't add this one as a neighbor
    ! local
    integer(I4B) :: i, n
    class(VirtualExchangeType), pointer :: v_exg
    integer(I4B) :: neighbor_id
    integer(I4B) :: model_mask
    type(STLVecInt) :: models_at_depth !< model ids at a certain depth, to
                                       !! recurse on for nbrs-of-nbrs search

    if (.not. present(mask)) then
      model_mask = 0
    else
      model_mask = mask
    end if

    call models_at_depth%init()

    if (is_root) then
      ! first layer in the recursive search
      call models_at_depth%push_back_unique(model_id)

      ! fetch primary neighbor
      if (this%prim_exchange%v_model1%id == this%owner%id) then
        neighbor_id = this%prim_exchange%v_model2%id
      else
        neighbor_id = this%prim_exchange%v_model1%id
      end if
      ! add
      call models_at_depth%push_back_unique(neighbor_id)
      call this%halo_models%push_back_unique(neighbor_id)
      call this%halo_exchanges%push_back_unique(this%prim_exchange%id)
    else
      ! find all direct neighbors of the model and add them,
      ! avoiding duplicates
      do i = 1, virtual_exchanges%Count()
        neighbor_id = -1
        v_exg => get_virtual_exchange_from_list(virtual_exchanges, i)
        if (v_exg%v_model1%id == model_id) then
          neighbor_id = v_exg%v_model2%id
        else if (v_exg%v_model2%id == model_id) then
          neighbor_id = v_exg%v_model1%id
        end if

        ! check if there is a neighbor, and it is not masked
        ! (to prevent back-and-forth connections)
        if (neighbor_id > 0) then
          ! check if masked
          if (neighbor_id == model_mask) cycle
          call models_at_depth%push_back_unique(neighbor_id)
          call this%halo_models%push_back_unique(neighbor_id)
          call this%halo_exchanges%push_back_unique(v_exg%id)
        end if
      end do
    end if

    depth = depth - 1
    if (depth == 0) then
      ! and we're done with this branch
      call models_at_depth%destroy()
      return
    end if

    ! now recurse on the neighbors up to the specified depth
    do n = 1, models_at_depth%size
      call this%addModelNeighbors(models_at_depth%at(n), virtual_exchanges, &
                                  depth, .false., model_id)
    end do

    ! we're done with the tree
    call models_at_depth%destroy()

  end subroutine addModelNeighbors

  !> @brief Define this connection, this is where the
  !! discretization (DISU) for the interface model is
  !< created!
  subroutine spatialcon_df(this)
    class(SpatialModelConnectionType) :: this !< this connection
    ! local
    integer(I4B) :: i
    class(VirtualModelType), pointer :: v_model

    ! create the grid connection data structure
    this%nr_connections = this%getNrOfConnections()
    call this%ig_builder%construct(this%owner, &
                                   this%nr_connections, &
                                   this%name)
    this%ig_builder%internalStencilDepth = this%int_stencil_depth
    this%ig_builder%exchangeStencilDepth = this%exg_stencil_depth
    this%ig_builder%haloExchanges => this%halo_exchanges
    do i = 1, this%halo_models%size
      v_model => get_virtual_model(this%halo_models%at(i))
      call this%ig_builder%addToRegionalModels(v_model)
    end do
    call this%setupGridConnection()

    this%neq = this%ig_builder%nrOfCells
    call this%allocateArrays()

  end subroutine spatialcon_df

  !> @brief Allocate the connection,
  !<
  subroutine spatialcon_ar(this)
    class(SpatialModelConnectionType) :: this !< this connection
    ! local
    integer(I4B) :: iface_idx, glob_idx
    class(GridConnectionType), pointer :: gc

    ! fill mapping to global index (which can be
    ! done now because moffset is set in sln_df)
    gc => this%ig_builder
    do iface_idx = 1, gc%nrOfCells
      glob_idx = gc%idxToGlobal(iface_idx)%index + &
                 gc%idxToGlobal(iface_idx)%v_model%moffset%get()
      gc%idxToGlobalIdx(iface_idx) = glob_idx
    end do

  end subroutine spatialcon_ar

  !> @brief set model pointers to connection
  !<
  subroutine spatialcon_setmodelptrs(this)
    class(SpatialModelConnectionType) :: this !< this connection

    ! point x, ibound, and rhs to connection
    this%interface_model%x => this%x
    call mem_checkin(this%interface_model%x, 'X', &
                     this%interface_model%memoryPath, 'X', &
                     this%memoryPath)
    this%interface_model%rhs => this%rhs
    call mem_checkin(this%interface_model%rhs, 'RHS', &
                     this%interface_model%memoryPath, 'RHS', &
                     this%memoryPath)
    this%interface_model%ibound => this%active
    call mem_checkin(this%interface_model%ibound, 'IBOUND', &
                     this%interface_model%memoryPath, 'IBOUND', &
                     this%memoryPath)

  end subroutine spatialcon_setmodelptrs

  !> @brief map interface model connections to our sparse matrix,
  !< analogously to what happens in sln_connect.
  subroutine spatialcon_connect(this)
    class(SpatialModelConnectionType) :: this !< this connection
    ! local
    type(sparsematrix) :: sparse
    class(MatrixBaseType), pointer :: matrix_base

    call sparse%init(this%neq, this%neq, 7)
    call this%interface_model%model_ac(sparse)

    ! create amat from sparse
    call this%createCoefficientMatrix(sparse)
    call sparse%destroy()

    ! map connections
    matrix_base => this%matrix
    call this%interface_model%model_mc(matrix_base)
    call this%maskOwnerConnections()

  end subroutine spatialcon_connect

  !> @brief Mask the owner's connections
  !!
  !! Determine which connections are handled by the interface model
  !! (using the connections object in its discretization) and
  !< set their mask to zero for the owning model.
  subroutine maskOwnerConnections(this)
    use CsrUtilsModule, only: getCSRIndex
    class(SpatialModelConnectionType) :: this !< the connection
    ! local
    integer(I4B) :: ipos, n, m, nloc, mloc, csr_idx
    type(ConnectionsType), pointer :: conn

    ! set the mask on connections that are calculated by the interface model
    conn => this%interface_model%dis%con
    do n = 1, conn%nodes
      ! only for connections internal to the owning model
      if (.not. this%ig_builder%idxToGlobal(n)%v_model == this%owner) then
        cycle
      end if
      nloc = this%ig_builder%idxToGlobal(n)%index

      do ipos = conn%ia(n) + 1, conn%ia(n + 1) - 1
        m = conn%ja(ipos)
        if (.not. this%ig_builder%idxToGlobal(m)%v_model == this%owner) then
          cycle
        end if
        mloc = this%ig_builder%idxToGlobal(m)%index

        if (conn%mask(ipos) > 0) then
          ! calculated by interface model, set local model's mask to zero
          csr_idx = getCSRIndex(nloc, mloc, this%owner%ia, this%owner%ja)
          if (csr_idx == -1) then
            ! this can only happen with periodic boundary conditions,
            ! then there is no need to set the mask
            if (this%ig_builder%isPeriodic(nloc, mloc)) cycle

            write (*, *) 'Error: cannot find cell connection in global system'
            call ustop()
          end if

          if (this%owner%dis%con%mask(csr_idx) > 0) then
            call this%owner%dis%con%set_mask(csr_idx, 0)
          else
            ! edge case, this connection is already being calculated
            ! so we ignore it here. This can happen in the overlap
            ! between two different exchanges when a larger stencil
            ! (XT3D) is applied.
            call conn%set_mask(ipos, 0)
          end if
        end if
      end do
    end do

  end subroutine maskOwnerConnections

  !> @brief Add connections, handled by the interface model,
  !< to the global system's sparse
  subroutine spatialcon_ac(this, sparse)
    class(SpatialModelConnectionType) :: this !< this connection
    type(sparsematrix), intent(inout) :: sparse !< sparse matrix to store the connections
    ! local
    integer(I4B) :: n, m, ipos
    integer(I4B) :: icol_start, icol_end
    integer(I4B) :: nglo, mglo

    do n = 1, this%neq
      if (.not. this%ig_builder%idxToGlobal(n)%v_model == this%owner) then
        ! only add connections for own model to global matrix
        cycle
      end if

      nglo = this%ig_builder%idxToGlobal(n)%index + &
             this%ig_builder%idxToGlobal(n)%v_model%moffset%get()

      icol_start = this%matrix%get_first_col_pos(n)
      icol_end = this%matrix%get_last_col_pos(n)
      do ipos = icol_start, icol_end
        m = this%matrix%get_column(ipos)
        if (m == n) cycle
        mglo = this%ig_builder%idxToGlobal(m)%index + &
               this%ig_builder%idxToGlobal(m)%v_model%moffset%get()
        call sparse%addconnection(nglo, mglo, 1)
      end do

    end do

  end subroutine spatialcon_ac

  !> @brief Creates the mapping from the local system
  !< matrix to the global one
  subroutine spatialcon_mc(this, matrix_sln)
    use SimModule, only: ustop
    class(SpatialModelConnectionType) :: this !< this connection
    class(MatrixBaseType), pointer :: matrix_sln !< global matrix
    ! local
    integer(I4B) :: i, m, n, mglo, nglo, ipos, ipos_sln
    logical(LGP) :: is_owned

    allocate (this%ipos_to_sln(this%matrix%nja))
    do i = 1, this%matrix%nja
      this%ipos_to_sln(i) = -1
    end do

    do n = 1, this%neq
      is_owned = (this%ig_builder%idxToGlobal(n)%v_model == this%owner)
      if (.not. is_owned) cycle

      do ipos = this%matrix%ia(n), this%matrix%ia(n + 1) - 1
        m = this%matrix%ja(ipos)
        nglo = this%ig_builder%idxToGlobal(n)%index + &
               this%ig_builder%idxToGlobal(n)%v_model%moffset%get()
        mglo = this%ig_builder%idxToGlobal(m)%index + &
               this%ig_builder%idxToGlobal(m)%v_model%moffset%get()

        ipos_sln = matrix_sln%get_position(nglo, mglo)
        if (ipos_sln == -1) then
          ! this should not be possible
          write (*, *) 'Error: cannot find cell connection in global system'
          call ustop()
        end if
        this%ipos_to_sln(ipos) = ipos_sln

      end do
    end do

  end subroutine spatialcon_mc

  !> @brief Calculate (or adjust) matrix coefficients,
  !! in this case those which are determined or affected
  !< by the connection of a GWF model with its neighbors
  subroutine spatialcon_cf(this, kiter)
    class(SpatialModelConnectionType) :: this !< this connection
    integer(I4B), intent(in) :: kiter !< the iteration counter
    ! local
    integer(I4B) :: i

    ! reset interface system
    call this%matrix%zero_entries()
    do i = 1, this%neq
      this%rhs(i) = 0.0_DP
    end do

    ! calculate the interface model
    call this%interface_model%model_cf(kiter)

  end subroutine spatialcon_cf

  !> @brief Formulate coefficients from interface model
  !<
  subroutine spatialcon_fc(this, kiter, matrix_sln, rhs_sln, inwtflag)
    class(SpatialModelConnectionType) :: this !< this connection
    integer(I4B), intent(in) :: kiter !< the iteration counter
    class(MatrixBaseType), pointer :: matrix_sln !< the system matrix
    real(DP), dimension(:), intent(inout) :: rhs_sln !< global right-hand-side
    integer(I4B), optional, intent(in) :: inwtflag !< newton-raphson flag
    ! local
    integer(I4B) :: n, nglo
    integer(I4B) :: icol_start, icol_end, ipos
    class(MatrixBaseType), pointer :: matrix_base

    matrix_base => this%matrix
    call this%interface_model%model_fc(kiter, matrix_base, inwtflag)

    ! map back to solution matrix
    do n = 1, this%neq
      ! We only need the coefficients for our own model
      ! (i.e. rows in the matrix that belong to this%owner):
      if (.not. this%ig_builder%idxToGlobal(n)%v_model == this%owner) then
        cycle
      end if

      nglo = this%ig_builder%idxToGlobal(n)%index + &
             this%ig_builder%idxToGlobal(n)%v_model%moffset%get() - &
             matrix_sln%get_row_offset()
      rhs_sln(nglo) = rhs_sln(nglo) + this%rhs(n)

      icol_start = this%matrix%get_first_col_pos(n)
      icol_end = this%matrix%get_last_col_pos(n)
      do ipos = icol_start, icol_end
        call matrix_sln%add_value_pos(this%ipos_to_sln(ipos), &
                                      this%matrix%get_value_pos(ipos))
      end do
    end do

  end subroutine spatialcon_fc

  !> @brief Deallocation
  !<
  subroutine spatialcon_da(this)
    class(SpatialModelConnectionType) :: this !< this connection

    call mem_deallocate(this%neq)
    call mem_deallocate(this%int_stencil_depth)
    call mem_deallocate(this%exg_stencil_depth)
    call mem_deallocate(this%nr_connections)

    call mem_deallocate(this%x)
    call mem_deallocate(this%rhs)
    call mem_deallocate(this%active)

    call this%halo_models%destroy()
    call this%halo_exchanges%destroy()
    deallocate (this%halo_models)
    deallocate (this%halo_exchanges)
    call this%matrix%destroy()
    deallocate (this%matrix)

    call this%ig_builder%destroy()
    call this%iface_dist_vars%Clear(destroy=.true.)
    deallocate (this%ig_builder)
    deallocate (this%interface_map)
    deallocate (this%ipos_to_sln)

  end subroutine spatialcon_da

  !> @brief Creates the connection structure for the
  !! interface grid, starting from primary exchanges,
  !! then extending inward and outward, possibly across
  !! model boundaries.
  !<
  subroutine setupGridConnection(this)
    class(SpatialModelConnectionType) :: this !< this connection
    ! local

    ! connect cells from primary exchange
    call this%ig_builder%connectPrimaryExchange(this%prim_exchange)

    ! now scan for nbr-of-nbrs and create final data structures
    call this%ig_builder%extendConnection()

    ! construct the interface map
    call this%ig_builder%buildInterfaceMap()
    this%interface_map => this%ig_builder%interfaceMap

  end subroutine setupGridConnection

  !> @brief Allocation of scalars
  !<
  subroutine allocateScalars(this)
    use MemoryManagerModule, only: mem_allocate
    class(SpatialModelConnectionType) :: this !< this connection

    call mem_allocate(this%neq, 'NEQ', this%memoryPath)
    call mem_allocate(this%int_stencil_depth, 'INTSTDEPTH', this%memoryPath)
    call mem_allocate(this%exg_stencil_depth, 'EXGSTDEPTH', this%memoryPath)
    call mem_allocate(this%nr_connections, 'NROFCONNS', this%memoryPath)

  end subroutine allocateScalars

  !> @brief Allocation of arrays
  !<
  subroutine allocateArrays(this)
    use MemoryManagerModule, only: mem_allocate
    class(SpatialModelConnectionType) :: this !< this connection
    ! local
    integer(I4B) :: i

    call mem_allocate(this%x, this%neq, 'X', this%memoryPath)
    call mem_allocate(this%rhs, this%neq, 'RHS', this%memoryPath)
    call mem_allocate(this%active, this%neq, 'IACTIVE', this%memoryPath)

    ! c.f. NumericalSolution
    do i = 1, this%neq
      this%x(i) = DZERO
      this%active(i) = 1 ! default is active
      this%rhs(i) = DZERO
    end do

  end subroutine allocateArrays

  !> @brief Returns total nr. of primary connections
  !<
  function getNrOfConnections(this) result(nrConns)
    class(SpatialModelConnectionType) :: this !< this connection
    integer(I4B) :: nrConns
    !local

    nrConns = this%prim_exchange%nexg

  end function getNrOfConnections

  !> @brief Create connection's matrix (ia,ja,amat) from sparse
  !<
  subroutine createCoefficientMatrix(this, sparse)
    use SimModule, only: ustop
    class(SpatialModelConnectionType) :: this !< this connection
    type(sparsematrix), intent(inout) :: sparse !< the sparse matrix with the cell connections

    call sparse%sort()
    call this%matrix%init(sparse, this%memoryPath)

  end subroutine createCoefficientMatrix

  !> @brief Validate this connection
  !<
  subroutine validateConnection(this)
    class(SpatialModelConnectionType) :: this !< this connection
    ! local
    class(DisConnExchangeType), pointer :: conEx => null()
    character(len=LINELENGTH) :: errmsg

    conEx => this%prim_exchange
    if (conEx%ixt3d > 0) then
      ! if XT3D, we need these angles:
      if (conEx%v_model1%con_ianglex%get() == 0) then
        write (errmsg, '(a,a,a,a,a)') 'XT3D configured on the exchange ', &
          trim(conEx%name), ' but the discretization in model ', &
          trim(conEx%v_model1%name), ' has no ANGLDEGX specified'
        call store_error(errmsg)
      end if
      if (conEx%v_model2%con_ianglex%get() == 0) then
        write (errmsg, '(a,a,a,a,a)') 'XT3D configured on the exchange ', &
          trim(conEx%name), ' but the discretization in model ', &
          trim(conEx%v_model2%name), ' has no ANGLDEGX specified'
        call store_error(errmsg)
      end if
    end if

    if (count_errors() > 0) then
      call ustop()
    end if

  end subroutine validateConnection

  !> @brief Add a variable from the interface model to be
  !! synchronized at the configured stages by copying from
  !! the source memory in the models/exchanges that are part
  !< of this interface.
  subroutine cfg_dv(this, var_name, subcomp_name, map_type, &
                    sync_stages, exg_var_name)
    class(SpatialModelConnectionType) :: this !< this connection
    character(len=*) :: var_name !< name of variable, e.g. "K11"
    character(len=*) :: subcomp_name !< subcomponent, e.g. "NPF"
    integer(I4B) :: map_type !< type of variable map
    integer(I4B), dimension(:) :: sync_stages !< stages to sync
    character(len=*), optional :: exg_var_name !< needed for exchange variables, e.g. SIMVALS
    ! local
    type(DistVarType), pointer :: dist_var => null()
    class(*), pointer :: obj

    if (.not. present(exg_var_name)) exg_var_name = ''

    allocate (dist_var)
    dist_var%var_name = var_name
    dist_var%subcomp_name = subcomp_name
    dist_var%comp_name = this%interface_model%name
    dist_var%map_type = map_type
    dist_var%sync_stages = sync_stages
    dist_var%exg_var_name = exg_var_name

    obj => dist_var
    call this%iface_dist_vars%Add(obj)

  end subroutine cfg_dv

  !> @brief Cast to SpatialModelConnectionType
  !<
  function cast_as_smc(obj) result(res)
    implicit none
    class(*), pointer, intent(inout) :: obj !< object to be cast
    class(SpatialModelConnectionType), pointer :: res !< the instance of SpatialModelConnectionType
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    class is (SpatialModelConnectionType)
      res => obj
    end select
  end function cast_as_smc

  !> @brief Add connection to a list
  !<
  subroutine add_smc_to_list(list, conn)
    implicit none
    ! -- dummy
    type(ListType), intent(inout) :: list !< the list
    class(SpatialModelConnectionType), pointer, intent(in) :: conn !< the connection
    ! -- local
    class(*), pointer :: obj
    !
    obj => conn
    call list%Add(obj)
  end subroutine add_smc_to_list

  !> @brief Get the connection from a list
  !<
  function get_smc_from_list(list, idx) result(res)
    type(ListType), intent(inout) :: list !< the list
    integer(I4B), intent(in) :: idx !< the index of the connection
    class(SpatialModelConnectionType), pointer :: res !< the returned connection

    ! local
    class(*), pointer :: obj
    obj => list%GetItem(idx)
    res => cast_as_smc(obj)
  end function get_smc_from_list

end module SpatialModelConnectionModule
