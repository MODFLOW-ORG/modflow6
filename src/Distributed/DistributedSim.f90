module DistributedSimModule
  use KindModule, only: I4B, LGP
  use SimVariablesModule, only: idm_context, simulation_mode, nr_procs, proc_id, &
                                errmsg, warnmsg
  use ConstantsModule, only: LENMEMPATH, LENMODELNAME, LINELENGTH, LENPACKAGETYPE
  use ArrayHandlersModule, only: ifind
  use CharacterStringModule, only: CharacterStringType
  use SimModule, only: store_error, store_error_filename, count_errors, &
                       store_warning
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, mem_setptr, &
                                 mem_print_detailed, get_isize
  use MemoryManagerExtModule, only: mem_set_value
  use MemoryHelperModule, only: create_mem_path

  implicit none
  private

  public :: DistributedSimType
  public :: get_dsim

  type :: DistributedSimType
    character(len=LENMEMPATH) :: memory_path
    integer(I4B), pointer :: nr_models !< the total (global) number of models, equals the length of the model block in mfsim.nam
    integer(I4B), dimension(:), pointer, contiguous :: load_mask => null() !< mask for loading models and exchanges, 1 when active on this processor, else 0
    integer(I4B), dimension(:), pointer, contiguous :: model_ranks => null() !< load balance: model rank (0,...,nr_procs-1) per global model id (array index)
    logical(LGP), pointer :: print_ptable !< when true, the partition table is printed to file
  contains
    procedure :: create
    procedure :: get_load_mask
    procedure :: get_load_balance
    procedure :: destroy
    ! private
    procedure, private :: create_load_mask
    procedure, private :: set_load_balance_from_input
    procedure, private :: set_load_balance_default
    procedure, private :: validate_load_balance
    procedure, private :: print_load_balance
  end type

  ! singleton, private member
  class(DistributedSimType), private, pointer :: dist_sim => null()

contains

  !> @brief Get pointer to the distributed simulation object
  !<
  function get_dsim() result(ds)
    class(DistributedSimType), pointer :: ds

    if (.not. associated(dist_sim)) then
      allocate (dist_sim)
      call dist_sim%create()
    end if
    ds => dist_sim

  end function get_dsim

  !> Create the distributed simulation object from the simulation input ctx
  !<
  subroutine create(this)
    class(DistributedSimType) :: this
    ! local
    character(len=LENMEMPATH) :: input_mempath
    integer(I4B), pointer :: nmod

    this%memory_path = create_mem_path(component='SIM')

    input_mempath = create_mem_path(component='SIM', context=idm_context)
    call mem_setptr(nmod, 'NUMMODELS', input_mempath)

    call mem_allocate(this%nr_models, 'NUMMODELS', this%memory_path)
    this%nr_models = nmod

    call mem_allocate(this%print_ptable, 'PRINT_PTABLE', this%memory_path)
    this%print_ptable = .false.

  end subroutine create

  !> @brief Return pointer to the load mask for models
  !!
  !! Get a load mask to determine which models
  !! should be loaded by idm on this process. This is in
  !! sync with models create. The mask array is allocated
  !! with its size equal to the global number of models.
  !! It is returned as (1, 1, 0, 0, ... 0) with each entry
  !! being a load mask for the model at the corresponding
  !< location in the 'MNAME' array of the IDM.
  function get_load_mask(this) result(load_mask)
    class(DistributedSimType) :: this
    integer(I4B), dimension(:), pointer :: load_mask

    if (.not. associated(this%load_mask)) then
      call this%create_load_mask()
    end if
    load_mask => this%load_mask

  end function get_load_mask

  !> @brief Create a load mask for IDM from the load balance array
  !<
  subroutine create_load_mask(this)
    class(DistributedSimType) :: this
    ! local
    integer(I4B), dimension(:), pointer :: model_ranks => null() !< the load balance
    integer(I4B) :: m_id !< model id

    call mem_allocate(this%load_mask, this%nr_models, 'LOADMASK', &
                      this%memory_path)
    this%load_mask = 0

    ! get load balance (probably the first call, so creates it)
    model_ranks => this%get_load_balance()

    ! set mask from balance
    do m_id = 1, this%nr_models
      if (model_ranks(m_id) == proc_id) then
        this%load_mask(m_id) = 1
      else
        this%load_mask(m_id) = 0
      end if
    end do

  end subroutine create_load_mask

  !> @brief Get the model load balance for the simulation
  !<
  function get_load_balance(this) result(mranks)
    use SimVariablesModule, only: iout
    use UtlHpcInputModule, only: UtlHpcParamFoundType
    class(DistributedSimType) :: this !< this distributed sim instance
    integer(I4B), dimension(:), pointer :: mranks !< the load balance: array of ranks per model id
    ! local
    integer(I4B) :: isize
    logical(LGP) :: hpc6_present, partitions_present
    character(len=LENMEMPATH) :: simnam_mempath, hpc_mempath
    type(UtlHpcParamFoundType) :: found

    ! if load balance available, return here:
    if (associated(this%model_ranks)) then
      mranks => this%model_ranks
      return
    end if

    call mem_allocate(this%model_ranks, this%nr_models, 'MODELRANKS', &
                      this%memory_path)

    ! check for optional HPC file
    simnam_mempath = create_mem_path('SIM', 'NAM', idm_context)
    call get_isize('HPC6_FILENAME', simnam_mempath, isize)
    hpc6_present = isize > 0

    ! handle serial case
    if (simulation_mode == 'SEQUENTIAL') then
      if (hpc6_present) then
        write (warnmsg, *) "Ignoring PARTITIONS block in HPC file when "// &
          "running a serial process"
        call store_warning(warnmsg)
      end if

      ! single process, everything on cpu 0:
      this%model_ranks = 0
      mranks => this%model_ranks
      return
    end if

    ! continue for PARALLEL mode only:
    write (iout, '(/1x,a)') 'PROCESSING HPC DATA'

    hpc_mempath = create_mem_path('UTL', 'HPC', idm_context)
    ! source optional print input flag
    call mem_set_value(this%print_ptable, 'PRINT_TABLE', hpc_mempath, &
                       found%print_table)
    ! check if optional partition block exists
    call get_isize('MNAME', hpc_mempath, isize)
    partitions_present = isize > 0

    ! fill model ranks (i.e. the load balance)
    if (partitions_present) then
      ! set balance from HPC file
      call this%set_load_balance_from_input()
      call this%validate_load_balance()
      write (iout, '(1x,a)') 'Read partition data from HPC file'
    else
      ! no HPC file present, set balance with default algorithm
      call this%set_load_balance_default()
      write (iout, '(1x,a)') 'Generate default partition data'
    end if

    mranks => this%model_ranks

    ! print to listing file
    if (this%print_ptable) then
      call this%print_load_balance()
    end if

    write (iout, '(1x,a)') 'END OF HPC DATA'

  end function get_load_balance

  !> @brief Load load balance from the input configuration
  !<
  subroutine set_load_balance_from_input(this)
    class(DistributedSimType) :: this !< this distributed sim instance
    ! local
    character(len=LENMEMPATH) :: simnam_mempath, hpc_mempath
    character(len=LENMODELNAME) :: model_name
    type(CharacterStringType), dimension(:), contiguous, pointer :: mnames !< model names (all) from the simulation nam file
    type(CharacterStringType), dimension(:), contiguous, pointer :: mnames_hpc !< model names in the hpc file
    integer(I4B), dimension(:), contiguous, pointer :: mranks_hpc !< rank numbers in the hpc file
    integer(I4B) :: i, model_idx
    integer(I4B) :: target_rank
    integer(I4B), dimension(:), allocatable :: rank_used
    type(CharacterStringType), dimension(:), contiguous, pointer :: hpc_names !< helper array to get the hpc filename
    character(len=LINELENGTH) :: hpc_filename !< the HPC option file

    ! set to uninitialized
    this%model_ranks = -1

    ! from IDM
    simnam_mempath = create_mem_path('SIM', 'NAM', idm_context)
    hpc_mempath = create_mem_path('UTL', 'HPC', idm_context)
    call mem_setptr(mnames, 'MNAME', simnam_mempath)
    call mem_setptr(mnames_hpc, 'MNAME', hpc_mempath)
    call mem_setptr(mranks_hpc, 'MRANK', hpc_mempath)
    call mem_setptr(hpc_names, 'HPC6_FILENAME', simnam_mempath)

    ! FILEIN options give an array, so take the first:
    hpc_filename = hpc_names(1)

    ! check: valid model names
    do i = 1, size(mnames_hpc)
      if (ifind(mnames, mnames_hpc(i)) == -1) then
        model_name = mnames_hpc(i)
        write (errmsg, *) "HPC input error: undefined model name (", &
          trim(model_name), ")"
        call store_error(errmsg)
      end if
    end do
    ! check: valid ranks
    do i = 1, size(mranks_hpc)
      target_rank = mranks_hpc(i)
      if (target_rank < 0 .or. target_rank > nr_procs - 1) then
        model_name = mnames_hpc(i)
        write (errmsg, '(a,i0,2a)') "HPC input error: invalid target rank (", &
          target_rank, ") for model ", trim(model_name)
        call store_error(errmsg)
      end if
    end do
    if (count_errors() > 0) then
      call store_error_filename(hpc_filename)
    end if

    ! construct rank array
    do i = 1, size(mnames_hpc)
      model_idx = ifind(mnames, mnames_hpc(i))
      this%model_ranks(model_idx) = mranks_hpc(i)
    end do

    ! check: all models acquired rank
    do i = 1, size(this%model_ranks)
      if (this%model_ranks(i) == -1) then
        model_name = mnames(i)
        write (errmsg, '(2a)') "HPC input error: no target rank for model ", &
          trim(model_name)
        call store_error(errmsg)
      end if
    end do
    if (count_errors() > 0) then
      call store_error_filename(hpc_filename)
    end if

    ! check: no idle ranks
    allocate (rank_used(nr_procs))
    rank_used = 0
    do i = 1, size(this%model_ranks)
      if (this%model_ranks(i) >= 0 .and. this%model_ranks(i) < nr_procs) then
        rank_used(this%model_ranks(i) + 1) = 1
      end if
    end do
    do i = 1, size(rank_used)
      if (rank_used(i) == 0) then
        write (errmsg, '(a,i0,a)') "HPC input error: rank ", i - 1, &
          " has no models assigned"
        call store_error(errmsg)
      end if
    end do
    deallocate (rank_used)
    if (count_errors() > 0) then
      call store_error_filename(hpc_filename)
    end if

  end subroutine set_load_balance_from_input

  !> @brief Distribute the models over the available
  !! processes in a parallel run. Expects an array sized
  !< to the number of models in the global simulation
  subroutine set_load_balance_default(this)
    class(DistributedSimType) :: this !< this distributed sim. instance
    ! local
    integer(I4B) :: im, imm, ie, ip, cnt
    integer(I4B) :: nr_models, nr_gwf_models
    integer(I4B) :: nr_exchanges
    integer(I4B) :: min_per_proc, nr_left
    integer(I4B) :: rank
    integer(I4B), dimension(:), allocatable :: nr_models_proc
    character(len=LENPACKAGETYPE) :: model_type_str
    character(len=LENMEMPATH) :: input_mempath
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mtypes !< model types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mnames !< model names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: etypes !< exg types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: emnames_a !< model a names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: emnames_b !< model b names

    this%model_ranks = 0

    ! load IDM data
    input_mempath = create_mem_path('SIM', 'NAM', idm_context)
    call mem_setptr(mtypes, 'MTYPE', input_mempath)
    call mem_setptr(mnames, 'MNAME', input_mempath)
    call mem_setptr(etypes, 'EXGTYPE', input_mempath)
    call mem_setptr(emnames_a, 'EXGMNAMEA', input_mempath)
    call mem_setptr(emnames_b, 'EXGMNAMEB', input_mempath)

    ! count flow models
    nr_models = size(mnames)
    nr_gwf_models = 0
    do im = 1, nr_models
      if (mtypes(im) == 'GWF6') then
        nr_gwf_models = nr_gwf_models + 1
      end if

      if (mtypes(im) == 'GWF6' .or. &
          mtypes(im) == 'GWT6' .or. &
          mtypes(im) == 'GWE6') then
        cycle
      end if

      model_type_str = mtypes(im)
      write (errmsg, *) 'Model type ', model_type_str, &
        ' not supported in parallel mode.'
      call store_error(errmsg, terminate=.true.)
    end do

    ! calculate nr of flow models for each rank
    allocate (nr_models_proc(nr_procs))
    min_per_proc = nr_gwf_models / nr_procs
    nr_left = nr_gwf_models - nr_procs * min_per_proc
    cnt = 1
    do ip = 1, nr_procs
      rank = ip - 1
      nr_models_proc(ip) = min_per_proc
      if (rank < nr_left) then
        nr_models_proc(ip) = nr_models_proc(ip) + 1
      end if
    end do

    ! assign ranks for flow models
    rank = 0
    do im = 1, nr_models
      if (mtypes(im) == 'GWF6') then
        if (nr_models_proc(rank + 1) == 0) then
          rank = rank + 1
        end if
        this%model_ranks(im) = rank
        nr_models_proc(rank + 1) = nr_models_proc(rank + 1) - 1
      end if
    end do

    ! match other models to flow
    nr_exchanges = size(etypes)
    do im = 1, nr_models
      if (mtypes(im) == 'GWT6') then

        ! find match
        do ie = 1, nr_exchanges
          if (etypes(ie) == 'GWF6-GWT6' .and. mnames(im) == emnames_b(ie)) then
            rank = 0
            do imm = 1, nr_models
              if (mnames(imm) == emnames_a(ie)) then
                rank = this%model_ranks(imm)
                exit
              end if
            end do
            this%model_ranks(im) = rank
            exit
          end if
        end do

      else if (mtypes(im) == 'GWE6') then
        do ie = 1, nr_exchanges
          if (etypes(ie) == 'GWF6-GWE6' .and. mnames(im) == emnames_b(ie)) then
            rank = 0
            do imm = 1, nr_models
              if (mnames(imm) == emnames_a(ie)) then
                rank = this%model_ranks(imm)
                exit
              end if
            end do
            this%model_ranks(im) = rank
            exit
          end if
        end do

      else
        cycle ! e.g., for a flow model
      end if
    end do

    ! cleanup
    deallocate (nr_models_proc)

  end subroutine set_load_balance_default

  !> @brief Check validity of load balance configuration
  !<
  subroutine validate_load_balance(this)
    class(DistributedSimType) :: this
    ! local
    character(len=LENMEMPATH) :: input_mempath
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mtypes !< model types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mnames !< model names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: etypes !< exg types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: emnames_a !< model a names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: emnames_b !< model b names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: hpc_names !< helper array to get hpc filename
    integer(I4B) :: ie
    integer(I4B) :: idx_a, idx_b
    integer(I4B) :: rank_a, rank_b
    integer(I4B) :: nr_exchanges
    character(len=LINELENGTH) :: hpc_filename !< the HPC option file
    character(len=LENMODELNAME) :: name_a, name_b
    character(len=LINELENGTH) :: exg_type

    ! load IDM data
    input_mempath = create_mem_path('SIM', 'NAM', idm_context)
    call mem_setptr(mtypes, 'MTYPE', input_mempath)
    call mem_setptr(mnames, 'MNAME', input_mempath)
    call mem_setptr(etypes, 'EXGTYPE', input_mempath)
    call mem_setptr(emnames_a, 'EXGMNAMEA', input_mempath)
    call mem_setptr(emnames_b, 'EXGMNAMEB', input_mempath)
    call mem_setptr(hpc_names, 'HPC6_FILENAME', input_mempath)

    ! FILEIN options give an array, so take the first:
    hpc_filename = hpc_names(1)

    nr_exchanges = size(etypes)

    ! loop over exchanges
    do ie = 1, nr_exchanges
      if (etypes(ie) == 'GWF6-GWT6' .or. etypes(ie) == 'GWF6-GWE6') then
        idx_a = ifind(mnames, emnames_a(ie))
        idx_b = ifind(mnames, emnames_b(ie))
        rank_a = this%model_ranks(idx_a)
        rank_b = this%model_ranks(idx_b)
        if (rank_a /= rank_b) then
          name_a = emnames_a(ie)
          name_b = emnames_b(ie)
          exg_type = etypes(ie)
          write (errmsg, '(7a)') "HPC input error: models ", &
            trim(name_a), " and ", trim(name_b), " with a ", &
            trim(exg_type), " coupling have to be assigned to the same rank"
          call store_error(errmsg)
        end if
      end if
    end do

    if (count_errors() > 0) then
      call store_error_filename(hpc_filename)
    end if

  end subroutine validate_load_balance

  !> @brief Print the load balance table to the listing file
  !<
  subroutine print_load_balance(this)
    use TableModule, only: TableType, table_cr
    use ConstantsModule, only: TABLEFT, TABCENTER
    use SimVariablesModule, only: iout, proc_id
    class(DistributedSimType) :: this
    ! local
    type(TableType), pointer :: inputtab => null()
    character(len=LINELENGTH) :: tag, term
    character(len=LENMEMPATH) :: input_mempath
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mtypes !< model types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mnames !< model names
    integer(I4B) :: im, nr_models

    input_mempath = create_mem_path('SIM', 'NAM', idm_context)

    call mem_setptr(mtypes, 'MTYPE', input_mempath)
    call mem_setptr(mnames, 'MNAME', input_mempath)

    ! setup table
    nr_models = size(mnames)
    call table_cr(inputtab, 'HPC', 'HPC PARTITION DATA')
    call inputtab%table_df(nr_models, 5, iout)

    ! add columns
    tag = 'ID'
    call inputtab%initialize_column(tag, 8, alignment=TABLEFT)
    tag = 'NAME'
    call inputtab%initialize_column(tag, LENMODELNAME + 4, alignment=TABLEFT)
    tag = 'TYPE'
    call inputtab%initialize_column(tag, 8, alignment=TABLEFT)
    tag = 'RANK'
    call inputtab%initialize_column(tag, 8, alignment=TABLEFT)
    tag = 'LOCAL'
    call inputtab%initialize_column(tag, 8, alignment=TABLEFT)

    do im = 1, nr_models
      call inputtab%add_term(im)
      term = mnames(im)
      call inputtab%add_term(term)
      term = mtypes(im)
      call inputtab%add_term(term)
      call inputtab%add_term(this%model_ranks(im))
      term = ''
      if (this%model_ranks(im) == proc_id) term = 'X'
      call inputtab%add_term(term)
    end do

    ! deallocate
    call inputtab%table_da()
    deallocate (inputtab)

  end subroutine print_load_balance

  !> @brief clean up
  !<
  subroutine destroy(this)
    class(DistributedSimType) :: this

    if (associated(this%load_mask)) then
      call mem_deallocate(this%load_mask)
      call mem_deallocate(this%model_ranks)
    end if

    call mem_deallocate(this%nr_models)
    call mem_deallocate(this%print_ptable)

    ! delete singleton instance
    if (associated(dist_sim)) deallocate (dist_sim)

  end subroutine destroy

end module DistributedSimModule
