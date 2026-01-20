module VirtualBaseModule
  use KindModule, only: I4B, DP, LGP
  use ListModule
  use ConstantsModule, only: LENVARNAME, LENMEMPATH, LENCOMPONENTNAME
  use MemoryTypeModule, only: MemoryType
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, &
                                 get_from_memorystore
  implicit none
  private

  public :: get_virtual_data_from_list

  integer(I4B), public, parameter :: MAP_ALL_TYPE = 0
  integer(I4B), public, parameter :: MAP_NODE_TYPE = 1
  integer(I4B), public, parameter :: MAP_CONN_TYPE = 2
  integer(I4B), public, parameter :: NR_VDC_ELEMENT_MAPS = 2

  !> This is a generic data structure to virtualize pieces
  !! of memory in 2 distinct ways:
  !!
  !!  1) Virtualize remote memory
  !!  This concerns memory residing on another process.
  !!  Typically, these pieces are subsets of certain model
  !!  and exchange data and lookup tables are kept with the
  !!  data to manage their mapping. The stage(s) at which
  !!  to synchronize the virtual memory is stored as well.
  !!
  !!  2) Virtualize local memory
  !!  In this case no virtual memory item is created, no
  !!  lookup tables and synchronization are necessary.
  !!  The virtual memory item will be pointed to the
  !!  original memory location at the requested
  !!  synchronization stage.
  !<
  type, abstract, public :: VirtualDataType
    logical(LGP) :: is_remote = .false. !< is remote memory, when true (default is false)
    character(len=LENVARNAME) :: var_name !< variable name
    character(len=LENCOMPONENTNAME) :: subcmp_name !< subcomponent name, e.g. package name
    character(len=LENMEMPATH) :: mem_path !< memory path
    integer(I4B), dimension(:), allocatable :: sync_stages !< stage(s) at which to synchronize
    integer(I4B) :: map_type !< the type of map
    logical(LGP) :: is_reduced !< when true, the discontinuous remote data is compressed
                               !! into contiguous virtual memory
    integer(I4B), dimension(:), &
      pointer, contiguous :: remote_elem_shift => null() !< contiguous list with 0-based remote indexes
                                !! (this is important for creating mpi data types)
    integer(I4B), dimension(:), &
      pointer, contiguous :: remote_to_virtual => null() !< sparse list which maps remote index to virtual
    type(MemoryType), pointer :: virtual_mt => null()
  contains
    procedure(vm_allocate_if), deferred :: vm_allocate
    procedure(vm_deallocate_if), deferred :: vm_deallocate
    procedure :: base => vm_to_base
    procedure :: check_stage => vm_check_stage
    procedure :: link => vm_link
    procedure :: get_element_map
  end type

  type, public, extends(VirtualDataType) :: VirtualIntType
    integer(I4B), private, pointer :: intsclr => null()
  contains
    procedure :: vm_allocate => vm_allocate_int
    procedure :: vm_deallocate => vm_deallocate_int
    procedure :: get => get_int
  end type

  type, public, extends(VirtualDataType) :: VirtualInt1dType
    integer(I4B), dimension(:), pointer, contiguous :: int1d => null()
  contains
    procedure :: vm_allocate => vm_allocate_int1d
    procedure :: vm_deallocate => vm_deallocate_int1d
    procedure :: get => get_int1d
    procedure :: get_array => get_array_int1d
  end type

  type, public, extends(VirtualDataType) :: VirtualDblType
    real(DP), private, pointer :: dblsclr => null()
  contains
    procedure :: vm_allocate => vm_allocate_dbl
    procedure :: vm_deallocate => vm_deallocate_dbl
    procedure :: get => get_dbl
  end type

  type, public, extends(VirtualDataType) :: VirtualDbl1dType
    real(DP), dimension(:), pointer, contiguous :: dbl1d => null()
  contains
    procedure :: vm_allocate => vm_allocate_dbl1d
    procedure :: vm_deallocate => vm_deallocate_dbl1d
    procedure :: get => get_dbl1d
    procedure :: get_array => get_array_dbl1d
  end type

  type, public, extends(VirtualDataType) :: VirtualDbl2dType
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d => null()
  contains
    procedure :: vm_allocate => vm_allocate_dbl2D
    procedure :: vm_deallocate => vm_deallocate_dbl2D
    procedure :: get => get_dbl2d
    procedure :: get_array => get_array_dbl2d
  end type

  ! etc...
  abstract interface
    subroutine vm_allocate_if(this, var_name, mem_path, shape)
      import VirtualDataType, I4B
      class(VirtualDataType) :: this
      character(len=*) :: var_name
      character(len=*) :: mem_path
      integer(I4B), dimension(:) :: shape
    end subroutine vm_allocate_if
    subroutine vm_deallocate_if(this)
      import VirtualDataType
      class(VirtualDataType) :: this
    end subroutine vm_deallocate_if
  end interface

contains

  function vm_to_base(this) result(base_ptr)
    class(VirtualDataType), target :: this
    class(VirtualDataType), pointer :: base_ptr

    base_ptr => this

  end function vm_to_base

  !> @brief Check if this data item requires syncing
  !< for this particular stage
  function vm_check_stage(this, stage) result(has_stage)
    use ArrayHandlersModule, only: ifind
    class(VirtualDataType), target :: this
    integer(I4B) :: stage, stg_idx
    logical(LGP) :: has_stage

    has_stage = .false.
    if (allocated(this%sync_stages)) then
      stg_idx = ifind(this%sync_stages, stage)
      has_stage = (stg_idx > 0)
    end if

  end function vm_check_stage

  subroutine vm_link(this)
    class(VirtualDataType), target :: this
    ! local
    logical(LGP) :: found

    call get_from_memorystore(this%var_name, this%mem_path, &
                              this%virtual_mt, found)

  end subroutine vm_link

  !> @brief Return array with offsets for elements
  !< mapped in this virtual data item
  function get_element_map(this) result(el_map)
    class(VirtualDataType), target :: this
    integer(I4B), dimension(:), pointer, contiguous :: el_map

    el_map => null()
    if (this%map_type > 0) then
      el_map => this%remote_elem_shift
    end if

  end function get_element_map

  subroutine vm_allocate_int(this, var_name, mem_path, shape)
    class(VirtualIntType) :: this
    character(len=*) :: var_name
    character(len=*) :: mem_path
    integer(I4B), dimension(:) :: shape

    call mem_allocate(this%intsclr, var_name, mem_path)

  end subroutine vm_allocate_int

  subroutine vm_deallocate_int(this)
    class(VirtualIntType) :: this

    if (this%is_remote) call mem_deallocate(this%intsclr)

  end subroutine vm_deallocate_int

  subroutine vm_allocate_int1d(this, var_name, mem_path, shape)
    class(VirtualInt1dType) :: this
    character(len=*) :: var_name
    character(len=*) :: mem_path
    integer(I4B), dimension(:) :: shape

    call mem_allocate(this%int1d, shape(1), var_name, mem_path)

  end subroutine vm_allocate_int1d

  subroutine vm_deallocate_int1d(this)
    class(VirtualInt1dType) :: this

    if (this%is_remote) call mem_deallocate(this%int1d)

  end subroutine vm_deallocate_int1d

  subroutine vm_allocate_dbl(this, var_name, mem_path, shape)
    class(VirtualDblType) :: this
    character(len=*) :: var_name
    character(len=*) :: mem_path
    integer(I4B), dimension(:) :: shape

    call mem_allocate(this%dblsclr, var_name, mem_path)

  end subroutine vm_allocate_dbl

  subroutine vm_deallocate_dbl(this)
    class(VirtualDblType) :: this

    if (this%is_remote) call mem_deallocate(this%dblsclr)

  end subroutine vm_deallocate_dbl

  subroutine vm_allocate_dbl1d(this, var_name, mem_path, shape)
    class(VirtualDbl1dType) :: this
    character(len=*) :: var_name
    character(len=*) :: mem_path
    integer(I4B), dimension(:) :: shape

    call mem_allocate(this%dbl1d, shape(1), var_name, mem_path)

  end subroutine vm_allocate_dbl1d

  subroutine vm_deallocate_dbl1d(this)
    class(VirtualDbl1dType) :: this

    if (this%is_remote) call mem_deallocate(this%dbl1d)

  end subroutine vm_deallocate_dbl1d

  subroutine vm_allocate_dbl2d(this, var_name, mem_path, shape)
    class(VirtualDbl2dType) :: this
    character(len=*) :: var_name
    character(len=*) :: mem_path
    integer(I4B), dimension(:) :: shape

    call mem_allocate(this%dbl2d, shape(1), shape(2), var_name, mem_path)

  end subroutine vm_allocate_dbl2d

  subroutine vm_deallocate_dbl2d(this)
    class(VirtualDbl2dType) :: this

    if (this%is_remote) call mem_deallocate(this%dbl2d)

  end subroutine vm_deallocate_dbl2d

  function get_int(this) result(val)
    class(VirtualIntType) :: this
    integer(I4B) :: val

    val = this%virtual_mt%intsclr

  end function get_int

  function get_int1d(this, i_rmt) result(val)
    class(VirtualInt1dType) :: this
    integer(I4B) :: i_rmt
    integer(I4B) :: val
    ! local
    integer(I4B) :: i_vrt

    if (this%is_reduced) then
      i_vrt = this%remote_to_virtual(i_rmt)
    else
      i_vrt = i_rmt
    end if
    val = this%virtual_mt%aint1d(i_vrt)

  end function get_int1d

  function get_array_int1d(this) result(array)
    class(VirtualInt1dType) :: this
    integer(I4B), dimension(:), pointer, contiguous :: array

    array => this%virtual_mt%aint1d

  end function get_array_int1d

  function get_dbl(this) result(val)
    class(VirtualDblType) :: this
    real(DP) :: val

    val = this%virtual_mt%dblsclr

  end function get_dbl

  function get_dbl1d(this, i_rmt) result(val)
    class(VirtualDbl1dType) :: this
    integer(I4B) :: i_rmt
    real(DP) :: val
    ! local
    integer(I4B) :: i_vrt

    if (this%is_reduced) then
      i_vrt = this%remote_to_virtual(i_rmt)
    else
      i_vrt = i_rmt
    end if
    val = this%virtual_mt%adbl1d(i_vrt)

  end function get_dbl1d

  function get_array_dbl1d(this) result(array)
    class(VirtualDbl1dType) :: this
    real(DP), dimension(:), pointer, contiguous :: array

    array => this%virtual_mt%adbl1d

  end function get_array_dbl1d

  function get_dbl2d(this, j_cmp, i_rmt) result(val)
    class(VirtualDbl2dType) :: this
    integer(I4B) :: j_cmp
    integer(I4B) :: i_rmt
    real(DP) :: val
    ! local
    integer(I4B) :: i_vrt

    if (this%is_reduced) then
      i_vrt = this%remote_to_virtual(i_rmt)
    else
      i_vrt = i_rmt
    end if
    val = this%virtual_mt%adbl2d(j_cmp, i_vrt)

  end function get_dbl2d

  function get_array_dbl2d(this) result(array)
    class(VirtualDbl2dType) :: this
    real(DP), dimension(:, :), pointer, contiguous :: array

    array => this%virtual_mt%adbl2d

  end function get_array_dbl2d

  function get_virtual_data_from_list(list, idx) result(vd)
    type(ListType) :: list
    integer(I4B) :: idx
    class(VirtualDataType), pointer :: vd
    ! local
    class(*), pointer :: obj_ptr

    vd => null()
    obj_ptr => list%GetItem(idx)
    select type (obj_ptr)
    class is (VirtualDataType)
      vd => obj_ptr
    end select

  end function get_virtual_data_from_list

end module VirtualBaseModule
