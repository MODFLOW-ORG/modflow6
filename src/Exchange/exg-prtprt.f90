!> @brief PRT-PRT exchange module
module PrtPrtExchangeModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENPACKAGENAME, LINELENGTH, LENMEMPATH, &
                             LENMODELNAME, LENBOUNDNAME, LENAUXNAME, &
                             DZERO, DONE
  use ErrorUtilModule, only: pstop
  use ListsModule, only: basemodellist, baseexchangelist
  use SimModule, only: store_error, store_error_filename, count_errors
  use SimVariablesModule, only: errmsg
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use BaseModelModule, only: BaseModelType, GetBaseModelFromList
  use PrtModule, only: PrtModelType
  use ParticleModule, only: ParticleType, create_particle, &
                            create_particle_store
  use ParticleEventModule, only: ParticleEventType
  use MemoryHelperModule, only: create_mem_path
  use ExgPrtprtInputModule, only: ExgPrtprtParamFoundType
  use BndModule, only: BndType, GetBndFromList
  use CharacterStringModule
  use ExplicitModelModule, only: ExplicitModelType
  use MemoryManagerModule, only: mem_allocate, mem_reallocate, mem_deallocate
  use BaseDisModule, only: DisBaseType
  use DisModule, only: DisType
  use DisvModule, only: DisvType
  use GeomUtilModule, only: get_ijk, point_in_polygon
  use PrtFmiModule, only: IFLOWFACE_TOP
  use MethodModule, only: LEVEL_MODEL, LEVEL_FEATURE

  implicit none
  private

  public :: PrtPrtExchangeType
  public :: prtprt_cr, try_transfer_particle
  public :: GetPrtPrtExchangeFromList

  !> @brief Exchange that transfers particles between two PRT models.
  type, extends(BaseExchangeType) :: PrtPrtExchangeType
    character(len=LINELENGTH), pointer :: filename => null() !< name of the input file
    integer(I4B), pointer :: m1id => null() !< local index of model 1 in basemodellist
    integer(I4B), pointer :: m2id => null() !< local index of model 2 in basemodellist
    type(PrtModelType), pointer :: prtmodel1 => null() !< pointer to model 1
    type(PrtModelType), pointer :: prtmodel2 => null() !< pointer to model 2
    integer(I4B), pointer :: nexg => null() !< number of connected cell pairs
    integer(I4B), pointer :: naux => null() !< number of auxiliary variables
    character(len=LENBOUNDNAME), dimension(:), &
      pointer, contiguous :: boundname => null() !< boundnames
    character(len=LENAUXNAME), dimension(:), &
      pointer, contiguous :: auxname => null() !< vector of auxname
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: auxname_cst => null() !< copy of vector auxname that can be stored in memory manager
    real(DP), dimension(:, :), pointer, contiguous :: auxvar => null() !< array of auxiliary variable values
    integer(I4B), pointer :: iprpak => null() !< print input flag
    integer(I4B), pointer :: iprflow => null() !< print flag for cell by cell flows
    integer(I4B), pointer :: ipakcb => null() !< save flag for cell by cell flows
    integer(I4B), pointer :: inamedbound => null() !< flag to read boundnames
    integer(I4B), dimension(:), pointer, contiguous :: nodem1 => null() !< node numbers in model 1
    integer(I4B), dimension(:), pointer, contiguous :: nodem2 => null() !< node numbers in model 2
    integer(I4B), dimension(:), pointer, contiguous :: ihc => null() !< horizontal connection indicator array, size: nexg
    character(len=LENMODELNAME) :: gwfmodelname1 = '' !< name of gwf model for prt model 1
    character(len=LENMODELNAME) :: gwfmodelname2 = '' !< name of gwf model for prt model 2
    real(DP), dimension(:), pointer, contiguous :: gwfsimvals => null() !< gwf flows for exchange
    integer(I4B), pointer :: iflowface1 => null() !< index in auxvar for cell face number in model 1
    integer(I4B), pointer :: iflowface2 => null() !< index in auxvar for cell face number in model 2
    class(PrtPrtExchangeType), pointer :: self => null() !< pointer to self (for event subscriptions)
  contains
    procedure :: exg_df
    procedure :: exg_ar
    procedure :: exg_da
    procedure :: exg_ad
    procedure :: connects_model
    procedure, private :: set_model_pointers
    procedure, private :: allocate_scalars
    procedure, private :: allocate_arrays
    procedure, private :: source_options
    procedure, private :: source_dimensions
    procedure, private :: source_data
    procedure, private :: noder
    procedure, private :: get_exchange_faces
    procedure, private :: transform_particle
    procedure, private :: transfer_particle
  end type PrtPrtExchangeType

contains

  !> @brief Create a new PRT-PRT exchange and register it.
  subroutine prtprt_cr(filename, name, id, m1_id, m2_id, input_mempath)
    ! modules
    use SimVariablesModule, only: model_loc_idx
    ! dummy
    character(len=*), intent(in) :: filename !< exchange input file
    character(len=*), intent(in) :: name !< exchange name
    integer(I4B), intent(in) :: id !< exchange id
    integer(I4B), intent(in) :: m1_id !< global model 1 id
    integer(I4B), intent(in) :: m2_id !< global model 2 id
    character(len=*), intent(in) :: input_mempath !< IDM memory path for this exchange
    ! local
    class(BaseExchangeType), pointer :: baseexchange => null()
    type(PrtPrtExchangeType), pointer :: exchange => null()

    allocate (exchange)
    baseexchange => exchange
    call AddBaseExchangeToList(baseexchangelist, baseexchange)

    exchange%id = id
    exchange%name = name
    exchange%memoryPath = create_mem_path(exchange%name)
    exchange%input_mempath = input_mempath
    ! store a self-pointer for use in event subscriptions
    exchange%self => exchange

    ! allocate scalars and set defaults
    call exchange%allocate_scalars()
    exchange%filename = filename
    exchange%typename = 'PRT-PRT'

    ! NB: convert global model id to local index in basemodellist
    exchange%m1id = model_loc_idx(m1_id)
    exchange%m2id = model_loc_idx(m2_id)

    ! set model pointers
    call exchange%set_model_pointers()
  end subroutine prtprt_cr

  !> @brief Whether either of the models is connected by this exchange.
  logical function connects_model(this, model)
    class(PrtPrtExchangeType) :: this
    class(BaseModelType), pointer, intent(in) :: model

    select type (model)
    class is (PrtModelType)
      connects_model = (associated(this%prtmodel1, model) .or. &
                        associated(this%prtmodel2, model))
    end select
  end function connects_model

  !> @brief Define phase: store model pointers.
  subroutine exg_df(this)
    class(PrtPrtExchangeType) :: this
    ! local
    integer(I4B) :: iout
    class(*), pointer :: p

    iout = this%prtmodel1%iout
    call this%source_options(iout)
    call this%source_dimensions(iout)
    call this%allocate_arrays()
    call this%source_data(iout)
    !
    ! -- Validate that IFLOWFACE1 and IFLOWFACE2 are present
    if (this%iflowface1 == 0) then
      write (errmsg, '(3a)') 'PRT-PRT exchange ', trim(this%name), &
        ' requires IFLOWFACE1 as an auxiliary variable.'
      call store_error(errmsg)
    end if
    if (this%iflowface2 == 0) then
      write (errmsg, '(3a)') 'PRT-PRT exchange ', trim(this%name), &
        ' requires IFLOWFACE2 as an auxiliary variable.'
      call store_error(errmsg)
    end if
    if (count_errors() > 0) then
      call store_error_filename(this%filename)
    end if
    !
    p => this%self ! self-pointer. can't point directly to `this`
    call this%prtmodel1%handlers%subscribe(try_transfer_particle, p)
    call this%prtmodel2%handlers%subscribe(try_transfer_particle, p)
  end subroutine exg_df

  !> @brief Allocate and read: load the connected-cell list from IDM.
  subroutine exg_ar(this)
    class(PrtPrtExchangeType) :: this
  end subroutine exg_ar

  !> @brief Deallocate memory.
  subroutine exg_da(this)
    class(PrtPrtExchangeType) :: this

    call mem_deallocate(this%m1id)
    call mem_deallocate(this%m2id)
    call mem_deallocate(this%nexg)
    call mem_deallocate(this%naux)
    call mem_deallocate(this%auxname, 'AUXNAME', this%memoryPath)
    call mem_deallocate(this%auxname_cst, 'AUXNAME_CST', this%memoryPath)
    call mem_deallocate(this%iprpak)
    call mem_deallocate(this%iprflow)
    call mem_deallocate(this%ipakcb)
    call mem_deallocate(this%inamedbound)
    call mem_deallocate(this%iflowface1)
    call mem_deallocate(this%iflowface2)
    if (associated(this%nodem1)) then
      call mem_deallocate(this%nodem1, 'NODEM1', this%memoryPath)
    end if
    if (associated(this%nodem2)) then
      call mem_deallocate(this%nodem2, 'NODEM2', this%memoryPath)
    end if
    if (associated(this%ihc)) then
      call mem_deallocate(this%ihc, 'IHC', this%memoryPath)
    end if
    if (associated(this%auxvar)) then
      call mem_deallocate(this%auxvar, 'AUXVAR', this%memoryPath)
    end if
    if (associated(this%boundname)) then
      deallocate (this%boundname)
    end if
    if (associated(this%gwfsimvals)) then
      call mem_deallocate(this%gwfsimvals, 'GWFSIMVALS', this%memoryPath)
    end if
  end subroutine exg_da

  !> @brief Inject exchange boundary flows into FMI. Flows are
  !! consolidated, so that a coarse cell connected to multiple
  !! refined cells in another model will have as boundary flow
  !! the sum of flows associated with the cells connected to it.
  subroutine exg_ad(this)
    class(PrtPrtExchangeType) :: this
    integer(I4B) :: iexg, ic1, ic2, if1, if2
    real(DP) :: q

    do iexg = 1, this%nexg
      ic1 = this%nodem1(iexg)
      ic2 = this%nodem2(iexg)
      q = abs(this%gwfsimvals(iexg))
      call this%get_exchange_faces(iexg, if1, if2)
      call this%prtmodel1%fmi%add_boundary_flow(ic1, if1, -q)
      call this%prtmodel2%fmi%add_boundary_flow(ic2, if2, q)
    end do
  end subroutine

  subroutine get_exchange_faces(this, iexg, if1, if2)
    class(PrtPrtExchangeType) :: this
    integer(I4B), intent(in) :: iexg !< exchange connection index
    integer(I4B), intent(out) :: if1, if2 !< face numbers
    !
    ! -- Get face numbers from auxiliary variables
    if1 = nint(this%auxvar(this%iflowface1, iexg))
    if2 = nint(this%auxvar(this%iflowface2, iexg))
  end subroutine

  subroutine set_model_pointers(this)
    class(PrtPrtExchangeType) :: this
    class(BaseModelType), pointer :: mb => null()

    mb => GetBaseModelFromList(basemodellist, this%m1id)
    select type (mb)
    type is (PrtModelType)
      this%prtmodel1 => mb
    end select

    mb => GetBaseModelFromList(basemodellist, this%m2id)
    select type (mb)
    type is (PrtModelType)
      this%prtmodel2 => mb
    end select

    if (.not. associated(this%prtmodel1)) then
      write (errmsg, '(3a)') 'Problem with PRT-PRT exchange ', &
        trim(this%name), ': model 1 is not a PRT model.'
      call store_error(errmsg, terminate=.true.)
    end if
    if (.not. associated(this%prtmodel2)) then
      write (errmsg, '(3a)') 'Problem with PRT-PRT exchange ', &
        trim(this%name), ': model 2 is not a PRT model.'
      call store_error(errmsg, terminate=.true.)
    end if
  end subroutine set_model_pointers

  subroutine allocate_scalars(this)
    use MemoryManagerModule, only: mem_allocate
    class(PrtPrtExchangeType) :: this

    allocate (this%filename)
    this%filename = ''

    call mem_allocate(this%m1id, 'M1ID', this%memoryPath)
    call mem_allocate(this%m2id, 'M2ID', this%memoryPath)
    call mem_allocate(this%nexg, 'NEXG', this%memoryPath)
    call mem_allocate(this%naux, 'NAUX', this%memoryPath)
    call mem_allocate(this%iprpak, 'IPRPAK', this%memoryPath)
    call mem_allocate(this%iprflow, 'IPRFLOW', this%memoryPath)
    call mem_allocate(this%ipakcb, 'IPAKCB', this%memoryPath)
    call mem_allocate(this%inamedbound, 'INAMEDBOUND', this%memoryPath)
    call mem_allocate(this%auxname, LENAUXNAME, 0, &
                      'AUXNAME', this%memoryPath)
    call mem_allocate(this%auxname_cst, LENAUXNAME, 0, &
                      'AUXNAME_CST', this%memoryPath)
    call mem_allocate(this%iflowface1, 'IFLOWFACE1', this%memoryPath)
    call mem_allocate(this%iflowface2, 'IFLOWFACE2', this%memoryPath)

    this%m1id = 0
    this%m2id = 0
    this%nexg = 0
    this%naux = 0
    this%iprpak = 0
    this%iprflow = 0
    this%ipakcb = 0
    this%inamedbound = 0
    this%iflowface1 = 0
    this%iflowface2 = 0
  end subroutine allocate_scalars

  subroutine allocate_arrays(this)
    use MemoryManagerModule, only: mem_allocate
    class(PrtPrtExchangeType) :: this

    call mem_allocate(this%nodem1, this%nexg, 'NODEM1', this%memoryPath)
    call mem_allocate(this%nodem2, this%nexg, 'NODEM2', this%memoryPath)
    call mem_allocate(this%ihc, this%nexg, 'IHC', this%memoryPath)
    ! auxname array is allocated while parsing
    call mem_allocate(this%auxvar, this%naux, this%nexg, &
                      'AUXVAR', this%memoryPath)
    ! allocate boundname
    if (this%inamedbound == 1) then
      allocate (this%boundname(this%nexg))
    else
      allocate (this%boundname(1))
    end if
    this%boundname(:) = ''
  end subroutine allocate_arrays

  !> @brief Read options from IDM.
  subroutine source_options(this, iout)
    use MemoryManagerExtModule, only: mem_set_value
    use ArrayHandlersModule, only: ifind
    class(PrtPrtExchangeType) :: this
    integer(I4B), intent(in) :: iout
    ! local
    type(ExgPrtprtParamFoundType) :: found
    logical(LGP) :: found_naux
    integer(I4B) :: n, ival

    call mem_set_value(this%gwfmodelname1, 'GWFMODELNAME1', this%input_mempath, &
                       found%gwfmodelname1)
    call mem_set_value(this%gwfmodelname2, 'GWFMODELNAME2', this%input_mempath, &
                       found%gwfmodelname2)
    call mem_set_value(this%naux, 'NAUX', this%input_mempath, found_naux)
    call mem_set_value(this%ipakcb, 'IPAKCB', this%input_mempath, found%ipakcb)
    call mem_set_value(this%iprpak, 'IPRPAK', this%input_mempath, found%iprpak)
    call mem_set_value(this%iprflow, 'IPRFLOW', this%input_mempath, found%iprflow)
    call mem_set_value(this%inamedbound, 'BOUNDNAMES', this%input_mempath, &
                       found%boundnames)

    ! reallocate aux arrays if aux variables provided
    if (found_naux .and. this%naux > 0) then
      call mem_reallocate(this%auxname, LENAUXNAME, this%naux, &
                          'AUXNAME', this%memoryPath)
      call mem_reallocate(this%auxname_cst, LENAUXNAME, this%naux, &
                          'AUXNAME_CST', this%memoryPath)
      call mem_set_value(this%auxname_cst, 'AUXILIARY', this%input_mempath, &
                         found%auxiliary)
      !
      do n = 1, this%naux
        this%auxname(n) = this%auxname_cst(n)
      end do
      !
      ! -- Look for IFLOWFACE1 and IFLOWFACE2 auxiliary variables
      ival = ifind(this%auxname, 'IFLOWFACE1')
      if (ival > 0) then
        this%iflowface1 = ival
      end if
      ival = ifind(this%auxname, 'IFLOWFACE2')
      if (ival > 0) then
        this%iflowface2 = ival
      end if
    end if

    if (.not. found%gwfmodelname1) then
      write (errmsg, '(3a)') 'PRT-PRT exchange ', trim(this%name), &
        ' requires that GWFMODELNAME1 be entered in the OPTIONS block.'
      call store_error(errmsg)
    end if
    if (.not. found%gwfmodelname2) then
      write (errmsg, '(3a)') 'PRT-PRT exchange ', trim(this%name), &
        ' requires that GWFMODELNAME2 be entered in the OPTIONS block.'
      call store_error(errmsg)
    end if

    if (found%ipakcb) then
      this%ipakcb = -1
      write (iout, '(4x,a)') &
        'EXCHANGE FLOWS WILL BE SAVED TO BINARY BUDGET FILES.'
    end if

    if (found%iprpak) then
      write (iout, '(4x,a)') &
        'THE LIST OF EXCHANGES WILL BE PRINTED.'
    end if

    if (found%iprflow) then
      write (iout, '(4x,a)') &
        'EXCHANGE FLOWS WILL BE PRINTED TO LIST FILES.'
    end if

    if (found%boundnames) then
      write (iout, '(4x,a)') 'EXCHANGE BOUNDARIES HAVE NAMES IN LAST COLUMN'
    end if
  end subroutine source_options

  !> @brief Read nexg from IDM.
  subroutine source_dimensions(this, iout)
    use MemoryManagerExtModule, only: mem_set_value
    class(PrtPrtExchangeType) :: this
    integer(I4B), intent(in) :: iout
    type(ExgPrtprtParamFoundType) :: found

    call mem_set_value(this%nexg, 'NEXG', this%input_mempath, found%nexg)

    if (.not. found%nexg) then
      write (errmsg, '(3a)') 'Required NEXG not found in exchange ', &
        trim(this%name), '.'
      call store_error(errmsg, terminate=.true.)
    end if
  end subroutine source_dimensions

  !> @brief Read cellidm1/cellidm2 from IDM and map to node numbers.
  subroutine source_data(this, iout)
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    use GeomUtilModule, only: get_node
    class(PrtPrtExchangeType) :: this
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B), dimension(:, :), contiguous, pointer :: cellidm1 => null()
    integer(I4B), dimension(:, :), contiguous, pointer :: cellidm2 => null()
    integer(I4B) :: iexg, iaux
    integer(I4B), dimension(:), contiguous, pointer :: ihc
    real(DP), dimension(:, :), contiguous, pointer :: auxvar
    type(CharacterStringType), dimension(:), contiguous, pointer :: boundname

    call mem_setptr(cellidm1, 'CELLIDM1', this%input_mempath)
    call mem_setptr(cellidm2, 'CELLIDM2', this%input_mempath)
    call mem_setptr(ihc, 'IHC', this%input_mempath)
    call mem_setptr(auxvar, 'AUXVAR', this%input_mempath)
    call mem_setptr(boundname, 'BOUNDNAME', this%input_mempath)

    do iexg = 1, this%nexg
      this%nodem1(iexg) = this%noder(this%prtmodel1, cellidm1(:, iexg), iout)
      this%nodem2(iexg) = this%noder(this%prtmodel2, cellidm2(:, iexg), iout)
      this%ihc(iexg) = ihc(iexg)
      do iaux = 1, this%naux
        this%auxvar(iaux, iexg) = auxvar(iaux, iexg)
      end do
      if (this%inamedbound == 1) then
        this%boundname(iexg) = boundname(iexg)
      end if
    end do
  end subroutine source_data

  !> @brief Returns reduced node number from user specified cell id.
  function noder(this, model, cellid, iout)
    ! modules
    use GeomUtilModule, only: get_node
    ! dummy
    class(PrtPrtExchangeType) :: this !< instance of exchange object
    type(PrtModelType), pointer, intent(in) :: model
    integer(I4B), dimension(:), intent(in) :: cellid
    integer(I4B), intent(in) :: iout !< the output file unit
    integer(I4B) :: noder, node
    !
    if (model%dis%ndim == 1) then
      node = cellid(1)
    elseif (model%dis%ndim == 2) then
      node = get_node(cellid(1), 1, cellid(2), &
                      model%dis%mshape(1), 1, &
                      model%dis%mshape(2))
    else
      node = get_node(cellid(1), cellid(2), cellid(3), &
                      model%dis%mshape(1), &
                      model%dis%mshape(2), &
                      model%dis%mshape(3))
    end if
    noder = model%dis%get_nodenumber(node, 0)
  end function noder

  !> @brief Event handler for particle transfers. Checks whether the particle
  !! should be transferred in either direction through any exchange connection
  !! and, if so, conducts the transfer and return signals that it has occurred.
  function try_transfer_particle(context, particle, event) result(transferred)
    use ModelExitEventModule, only: ModelExitEventType
    use PrtPrpModule, only: PrtPrpType, ExgPrtPrpType
    ! dummy
    class(*), pointer :: context
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    logical(LGP) :: transferred
    ! local
    type(PrtPrtExchangeType), pointer :: exg
    logical(LGP) :: should
    integer(I4B) :: iexg

    transferred = .false.

    select type (context)
    type is (PrtPrtExchangeType)
      exg => context
    class default
      return
    end select

    select type (event)
    type is (ModelExitEventType)
      ! model1->model2
      should = .false.
      do iexg = 1, exg%nexg
        if ((particle%itrdomain(LEVEL_MODEL) == exg%prtmodel1%id) .and. &
            (particle%icu == exg%nodem1(iexg))) then
          should = .true.
          exit
        end if
      end do
      if (should) then
        call exg%transfer_particle(particle, iexg=iexg, from_m1=.true.)
        return
      end if

      ! model2->model1
      should = .false.
      do iexg = 1, exg%nexg
        if ((particle%itrdomain(LEVEL_MODEL) == exg%prtmodel2%id) .and. &
            (particle%icu == exg%nodem2(iexg))) then
          should = .true.
          exit
        end if
      end do
      if (should) then
        call exg%transfer_particle(particle, iexg=iexg, from_m1=.false.)
        return
      end if
    end select
  end function try_transfer_particle

  !> @brief Transfer the particle from one model to the other through
  !! the given exchange connection.
  subroutine transfer_particle(this, particle, iexg, from_m1)
    use PrtPrpModule, only: PrtPrpType, ExgPrtPrpType
    class(PrtPrtExchangeType), intent(inout) :: this
    type(ParticleType), intent(inout), pointer :: particle
    integer(I4B), intent(in) :: iexg
    logical(LGP), intent(in) :: from_m1
    ! local
    type(PrtModelType), pointer :: prt_dst
    type(ExgPrtPrpType), pointer :: exgprp
    class(BndType), pointer :: exgprp_obj
    integer(I4B) :: icu_dst

    if (from_m1) then
      prt_dst => this%prtmodel2
    else
      prt_dst => this%prtmodel1
    end if

    exgprp_obj => GetBndFromList( &
                  prt_dst%bndlist, &
                  prt_dst%bndlist%Count())

    select type (exgprp_obj)
    type is (ExgPrtPrpType)
      call this%transform_particle(particle, particle%icu, icu_dst, &
                                   nint(this%auxvar(this%iflowface1, iexg)), &
                                   from_m1=from_m1)

      particle%icu = icu_dst
      particle%advancing = .true.
      particle%itrdomain(LEVEL_MODEL) = prt_dst%id
      particle%itrdomain(LEVEL_FEATURE) = &
        prt_dst%dis%get_nodenumber(particle%icu, 1)

      exgprp => exgprp_obj
      exgprp%has_pending = .true.
      exgprp%nparticles = exgprp%nparticles + 1
      call exgprp%particles%resize(exgprp%nparticles, exgprp%memoryPath)
      call exgprp%particles%put(particle, exgprp%nparticles)
    end select
  end subroutine transfer_particle

  !> @brief Transform particle coordinates from source to destination cell,
  !! by way of intermediate transformations via cell-normalized coordinates.
  !!
  !! Transforms particle from the source model's coordinate system to coords
  !! local to the source cell, then to coords local to the destination cell,
  !! and finally to the destination model's coordinate system.
  !!
  !! If the source cell connects to multiple cells in the destination model,
  !! the destination cells' connecting faces are together assumed to coincide
  !! with the source cell's connecting face. Coordinates are then computed by
  !! reference to the extent of the cell in the 2 dimensions parallel to the
  !! connected faces, and used to determine which among potential destination
  !! cells the particle should enter.
  !<
  subroutine transform_particle(this, particle, n_src, n_dst, f_src, from_m1)
    class(PrtPrtExchangeType) :: this
    type(ParticleType), pointer :: particle
    integer(I4B), intent(in) :: n_src ! source cell
    integer(I4B), intent(out) :: n_dst ! destination cell
    integer(I4B), intent(in) :: f_src ! source face (IFLOWFACE numbering)
    logical(LGP), intent(in) :: from_m1 ! true = m1->m2, false = m2->m1
    ! local
    integer(I4B), allocatable :: iexg_candidates(:)
    integer(I4B) :: n_candidates
    integer(I4B) :: f_dst, iexg, i
    integer(I4B) :: iflowface_src, iflowface_dst
    real(DP) :: x1min, x1max, y1min, y1max, z1min, z1max
    real(DP) :: x2min, x2max, y2min, y2max, z2min, z2max
    real(DP) :: xt, yt, zt ! transformed coordinates
    real(DP) :: tx, ty, tz ! normalized parameters [0,1]
    real(DP) :: dx1, dy1, dz1, dx2, dy2, dz2
    logical :: found
    class(PrtModelType), pointer :: prt_src, prt_dst

    ! set source/destination model
    if (from_m1) then
      prt_src => this%prtmodel1
      prt_dst => this%prtmodel2
      iflowface_src = this%iflowface1
      iflowface_dst = this%iflowface2
    else
      prt_src => this%prtmodel2
      prt_dst => this%prtmodel1
      iflowface_src = this%iflowface2
      iflowface_dst = this%iflowface1
    end if

    ! find candidate connections
    call get_connections(this, n_src, f_src, iexg_candidates, &
                         n_candidates, from_m1)

    ! compute source face bounds
    call get_face_bounds(prt_src%dis, n_src, f_src, &
                         x1min, x1max, y1min, y1max, z1min, z1max)

    ! compute destination face bounds
    call get_consolidated_bounds(this, iexg_candidates, n_candidates, &
                                 x2min, x2max, y2min, y2max, z2min, z2max, &
                                 prt_dst%dis, iflowface_dst)

    ! map to cell-local coords
    dx1 = x1max - x1min
    dy1 = y1max - y1min
    dz1 = z1max - z1min
    if (dx1 > DZERO) then
      tx = (particle%x - x1min) / dx1
    else
      tx = DZERO ! constant dimension
    end if
    if (dy1 > DZERO) then
      ty = (particle%y - y1min) / dy1
    else
      ty = DZERO
    end if
    if (dz1 > DZERO) then
      tz = (particle%z - z1min) / dz1
    else
      tz = DZERO
    end if

    ! map to destination model coords
    dx2 = x2max - x2min
    dy2 = y2max - y2min
    dz2 = z2max - z2min
    if (dx2 > DZERO) then
      xt = x2min + tx * dx2
    else
      xt = x2min
    end if
    if (dy2 > DZERO) then
      yt = y2min + ty * dy2
    else
      yt = y2min
    end if
    if (dz2 > DZERO) then
      zt = z2min + tz * dz2
    else
      zt = z2min
    end if

    ! which cell is the particle in?
    found = .false.
    do i = 1, n_candidates
      iexg = iexg_candidates(i)
      if (from_m1) then
        n_dst = this%nodem2(iexg)
        f_dst = nint(this%auxvar(iflowface_dst, iexg))
      else
        n_dst = this%nodem1(iexg)
        f_dst = nint(this%auxvar(iflowface_dst, iexg))
      end if
      if (point_in_face_bounds(prt_dst%dis, n_dst, f_dst, xt, yt, zt)) then
        particle%x = xt
        particle%y = yt
        particle%z = zt
        found = .true.
        exit
      end if
    end do
    if (.not. found) then
      write (errmsg, '(a,i0,a,i0,a,3g15.6,a)') &
        'Particle transferring from cell ', n_src, ' face ', f_src, &
        ' at position (', xt, yt, zt, &
        ') not found in any destination cell'
      call pstop(1, errmsg)
    end if

    deallocate (iexg_candidates)
  end subroutine transform_particle

  !> @brief Find all exchange connections with the same source cell and face
  subroutine get_connections(this, n_src, f_src, iexg_list, n_found, from_m1)
    class(PrtPrtExchangeType) :: this
    integer(I4B), intent(in) :: n_src, f_src
    integer(I4B), allocatable, intent(out) :: iexg_list(:)
    integer(I4B), intent(out) :: n_found
    logical(LGP), intent(in) :: from_m1
    integer(I4B) :: iexg, f_test, iflowface_src
    integer(I4B), allocatable :: temp_list(:)

    ! Set up source flowface index based on direction
    if (from_m1) then
      iflowface_src = this%iflowface1
    else
      iflowface_src = this%iflowface2
    end if

    ! Allocate temporary array (max size = all connections)
    allocate (temp_list(this%nexg))
    n_found = 0

    do iexg = 1, this%nexg
      if (from_m1) then
        if (this%nodem1(iexg) == n_src) then
          f_test = nint(this%auxvar(iflowface_src, iexg))
          if (f_test == f_src) then
            n_found = n_found + 1
            temp_list(n_found) = iexg
          end if
        end if
      else
        if (this%nodem2(iexg) == n_src) then
          f_test = nint(this%auxvar(iflowface_src, iexg))
          if (f_test == f_src) then
            n_found = n_found + 1
            temp_list(n_found) = iexg
          end if
        end if
      end if
    end do

    ! Copy to right-sized output array
    allocate (iexg_list(n_found))
    iexg_list(1:n_found) = temp_list(1:n_found)
    deallocate (temp_list)
  end subroutine get_connections

  !> @brief Get consolidated bounds across multiple destination faces
  subroutine get_consolidated_bounds(this, iexg_list, n, &
                                     xmin, xmax, ymin, ymax, zmin, zmax, &
                                     dis_dst, iflowface_dst)
    class(PrtPrtExchangeType) :: this
    integer(I4B), intent(in) :: iexg_list(:), n
    real(DP), intent(out) :: xmin, xmax, ymin, ymax, zmin, zmax
    class(DisBaseType), pointer :: dis_dst
    integer(I4B), intent(in) :: iflowface_dst
    integer(I4B) :: i, iexg, n_dst, f_dst
    real(DP) :: x0, x1, y0, y1, z0, z1

    ! Initialize with extreme values
    xmin = HUGE(DONE)
    xmax = -HUGE(DONE)
    ymin = HUGE(DONE)
    ymax = -HUGE(DONE)
    zmin = HUGE(DONE)
    zmax = -HUGE(DONE)

    ! Find min/max across all destination faces
    do i = 1, n
      iexg = iexg_list(i)
      ! Destination depends on which list is being used
      ! The iexg_list comes from get_connections which already
      ! filters by source, so all entries have the same source.
      ! We need to get the destination cell/face from the opposite side.
      if (iflowface_dst == this%iflowface2) then
        n_dst = this%nodem2(iexg)
      else
        n_dst = this%nodem1(iexg)
      end if
      f_dst = nint(this%auxvar(iflowface_dst, iexg))

      call get_face_bounds(dis_dst, n_dst, f_dst, x0, x1, y0, y1, z0, z1)

      xmin = min(xmin, x0)
      xmax = max(xmax, x1)
      ymin = min(ymin, y0)
      ymax = max(ymax, y1)
      zmin = min(zmin, z0)
      zmax = max(zmax, z1)
    end do
  end subroutine get_consolidated_bounds

  !> @brief Check if point is within face bounds (polymorphic wrapper)
  function point_in_face_bounds(dis, icell, iface, x, y, z) result(inside)
    class(DisBaseType), pointer :: dis
    integer(I4B), intent(in) :: icell, iface
    real(DP), intent(in) :: x, y, z
    logical(LGP) :: inside

    select type (dis)
    type is (DisType)
      inside = point_in_face_bounds_dis(dis, icell, iface, x, y, z)
    type is (DisvType)
      inside = point_in_face_bounds_disv(dis, icell, iface, x, y, z)
    class default
      inside = .false.
    end select
  end function point_in_face_bounds

  !> @brief Check if point is within DIS face bounds
  function point_in_face_bounds_dis(dis, icell, iface, x, y, z) result(inside)
    type(DisType), pointer, intent(in) :: dis
    integer(I4B), intent(in) :: icell, iface
    real(DP), intent(in) :: x, y, z
    logical(LGP) :: inside
    real(DP) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(DP), parameter :: TOL = 1.0d-9

    call get_face_bounds_dis(dis, icell, iface, &
                             xmin, xmax, ymin, ymax, zmin, zmax)

    inside = (x >= xmin - TOL .and. x <= xmax + TOL .and. &
              y >= ymin - TOL .and. y <= ymax + TOL .and. &
              z >= zmin - TOL .and. z <= zmax + TOL)
  end function point_in_face_bounds_dis

  !> @brief Check if point is within DISV face bounds
  function point_in_face_bounds_disv(dis, icell, iface, x, y, z) result(inside)
    type(DisvType), pointer, intent(in) :: dis
    integer(I4B), intent(in) :: icell, iface
    real(DP), intent(in) :: x, y, z
    logical(LGP) :: inside
    integer(I4B) :: k, n2d, m, iv, nvert_cell, iv1, iv2
    real(DP) :: ztop, zbot, x1, y1, x2, y2
    real(DP), allocatable :: poly(:, :)
    real(DP), parameter :: TOL = 1.0d-9

    n2d = mod(icell - 1, dis%ncpl) + 1
    k = (icell - 1) / dis%ncpl + 1

    if (k == 1) then
      ztop = dis%top1d(n2d)
    else
      ztop = dis%bot2d(n2d, k - 1)
    end if
    zbot = dis%bot2d(n2d, k)

    inside = .false.

    if (iface == -1 .or. iface == -2) then
      ! Top or bottom face - use 2D polygon containment
      if ((iface == -1 .and. abs(z - ztop) < TOL) .or. &
          (iface == -2 .and. abs(z - zbot) < TOL)) then
        ! Build polygon from vertices
        nvert_cell = dis%iavert(n2d + 1) - dis%iavert(n2d)
        allocate (poly(2, nvert_cell))
        do m = 1, nvert_cell
          iv = dis%javert(dis%iavert(n2d) + m - 1)
          poly(1, m) = dis%vertices(1, iv)
          poly(2, m) = dis%vertices(2, iv)
        end do
        inside = point_in_polygon(x, y, poly)
        deallocate (poly)
      end if

    else
      ! Lateral face - check z bounds and point on edge
      if (z < zbot - TOL .or. z > ztop + TOL) return

      ! Check if (x,y) is on the edge
      nvert_cell = dis%iavert(n2d + 1) - dis%iavert(n2d)
      iv1 = dis%javert(dis%iavert(n2d) + iface - 1)
      iv2 = dis%javert(dis%iavert(n2d) + mod(iface, nvert_cell))

      x1 = dis%vertices(1, iv1)
      y1 = dis%vertices(2, iv1)
      x2 = dis%vertices(1, iv2)
      y2 = dis%vertices(2, iv2)

      inside = point_on_line_segment(x, y, x1, y1, x2, y2, TOL)
    end if
  end function point_in_face_bounds_disv

  !> @brief Get face bounds (polymorphic wrapper)
  subroutine get_face_bounds(dis, icell, iface, xmin, xmax, ymin, ymax, &
                             zmin, zmax)
    class(DisBaseType), pointer :: dis
    integer(I4B), intent(in) :: icell, iface
    real(DP), intent(out) :: xmin, xmax, ymin, ymax, zmin, zmax

    select type (dis)
    type is (DisType)
      call get_face_bounds_dis(dis, icell, iface, xmin, xmax, &
                               ymin, ymax, zmin, zmax)
    type is (DisvType)
      call get_face_bounds_disv(dis, icell, iface, xmin, xmax, &
                                ymin, ymax, zmin, zmax)
    class default
      call store_error('Unsupported discretization type for PRT-PRT exchange')
      call store_error_filename('')
    end select
  end subroutine get_face_bounds

  !> @brief Get face bounds for DISV (vertex) grid
  subroutine get_face_bounds_disv(dis, icell, iface, xmin, xmax, ymin, ymax, &
                                  zmin, zmax)
    type(DisvType), pointer, intent(in) :: dis
    integer(I4B), intent(in) :: icell, iface
    real(DP), intent(out) :: xmin, xmax, ymin, ymax, zmin, zmax
    integer(I4B) :: k, n2d, m, iv, nvert_cell, iv1, iv2
    real(DP) :: ztop, zbot

    ! Get layer and 2D cell index
    n2d = mod(icell - 1, dis%ncpl) + 1
    k = (icell - 1) / dis%ncpl + 1

    ! Z bounds
    if (k == 1) then
      ztop = dis%top1d(n2d)
    else
      ztop = dis%bot2d(n2d, k - 1)
    end if
    zbot = dis%bot2d(n2d, k)

    if (iface == -1 .or. iface == -2) then
      ! Top or bottom face - horizontal polygon
      ! Get bounding box of polygon vertices
      xmin = HUGE(DONE)
      xmax = -HUGE(DONE)
      ymin = HUGE(DONE)
      ymax = -HUGE(DONE)

      do m = dis%iavert(n2d), dis%iavert(n2d + 1) - 1
        iv = dis%javert(m)
        xmin = min(xmin, dis%vertices(1, iv))
        xmax = max(xmax, dis%vertices(1, iv))
        ymin = min(ymin, dis%vertices(2, iv))
        ymax = max(ymax, dis%vertices(2, iv))
      end do

      if (iface == -1) then
        zmin = ztop
        zmax = ztop
      else
        zmin = zbot
        zmax = zbot
      end if

    else
      ! Lateral face - vertical rectangle along an edge
      nvert_cell = dis%iavert(n2d + 1) - dis%iavert(n2d)
      iv1 = dis%javert(dis%iavert(n2d) + iface - 1)
      iv2 = dis%javert(dis%iavert(n2d) + mod(iface, nvert_cell))

      xmin = min(dis%vertices(1, iv1), dis%vertices(1, iv2))
      xmax = max(dis%vertices(1, iv1), dis%vertices(1, iv2))
      ymin = min(dis%vertices(2, iv1), dis%vertices(2, iv2))
      ymax = max(dis%vertices(2, iv1), dis%vertices(2, iv2))
      zmin = zbot
      zmax = ztop
    end if
  end subroutine get_face_bounds_disv

  !> @brief Get face bounds for DIS (structured) grid
  subroutine get_face_bounds_dis(dis, icell, iface, xmin, xmax, ymin, ymax, &
                                 zmin, zmax)
    type(DisType), pointer, intent(in) :: dis
    integer(I4B), intent(in) :: icell, iface
    real(DP), intent(out) :: xmin, xmax, ymin, ymax, zmin, zmax
    integer(I4B) :: k, i, j
    real(DP) :: ztop, zbot

    call get_ijk(icell, dis%nrow, dis%ncol, dis%nlay, i, j, k)

    ! X bounds
    xmin = dis%xorigin
    if (j > 1) xmin = xmin + sum(dis%delr(1:j - 1))
    xmax = xmin + dis%delr(j)

    ! Y bounds (rows numbered north to south)
    ymax = dis%yorigin + sum(dis%delc)
    if (i > 1) ymax = ymax - sum(dis%delc(1:i - 1))
    ymin = ymax - dis%delc(i)

    ! Z bounds
    if (k == 1) then
      ztop = dis%top2d(j, i)
    else
      ztop = dis%bot3d(j, i, k - 1)
    end if
    zbot = dis%bot3d(j, i, k)
    zmin = zbot
    zmax = ztop

    ! Constrain to face (set min=max for constant dimension)
    select case (iface)
    case (1) ! West face
      xmax = xmin
    case (3) ! East face
      xmin = xmax
    case (4) ! South face
      ymax = ymin
    case (2) ! North face
      ymin = ymax
    case (-2) ! Bottom face
      zmax = zbot
    case (-1) ! Top face
      zmin = ztop
    end select

  end subroutine get_face_bounds_dis

  !> @brief Check if a point lies on a line segment
  function point_on_line_segment(x, y, x1, y1, x2, y2, tol) result(on_seg)
    real(DP), intent(in) :: x, y, x1, y1, x2, y2, tol
    logical(LGP) :: on_seg
    real(DP) :: dx, dy, dx_seg, dy_seg, cross, dot, len_sq

    dx = x - x1
    dy = y - y1
    dx_seg = x2 - x1
    dy_seg = y2 - y1

    ! Cross product (zero if collinear)
    cross = abs(dx * dy_seg - dy * dx_seg)

    ! Dot product and length check (within segment)
    dot = dx * dx_seg + dy * dy_seg
    len_sq = dx_seg**2 + dy_seg**2

    on_seg = .false.
    if (cross < tol .and. dot >= -tol .and. dot <= len_sq + tol) then
      on_seg = .true.
    end if
  end function point_on_line_segment

  !> @brief Get PRT-PRT exchange from list
  function GetPrtPrtExchangeFromList(list, idx) result(res)
    use ListModule, only: ListType
    implicit none
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: idx
    class(PrtPrtExchangeType), pointer :: res
    class(*), pointer :: obj

    obj => list%GetItem(idx)
    res => CastAsPrtPrtExchange(obj)
  end function GetPrtPrtExchangeFromList

  !> @brief Cast object as PRT-PRT exchange
  function CastAsPrtPrtExchange(obj) result(res)
    implicit none
    class(*), pointer, intent(inout) :: obj
    class(PrtPrtExchangeType), pointer :: res

    res => null()
    if (.not. associated(obj)) return

    select type (obj)
    class is (PrtPrtExchangeType)
      res => obj
    end select
  end function CastAsPrtPrtExchange

end module PrtPrtExchangeModule
