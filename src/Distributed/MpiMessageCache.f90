module MpiMessageCacheModule
  use KindModule, only: I4B
  use SimStagesModule, only: NR_SIM_STAGES
  use ListModule
  use STLVecIntModule
  use MpiUnitCacheModule
  implicit none
  private

  ! the message types for caching during a simulation stage:
  integer(I4B), public, parameter :: MPI_BDY_RCV = 1 !< receiving data (body) from ranks
  integer(I4B), public, parameter :: MPI_BDY_SND = 2 !< sending data (body) to ranks
  integer(I4B), public, parameter :: NR_MSG_TYPES = 2 !< the total number of message types to be cached

  ! expose this from the unit cache module
  public :: NO_CACHED_VALUE

  !> @brief Facility to cache the constructed MPI datatypes.
  !! This will avoid having to construct them over and over
  !! again for the communication inside the timestep loop.
  !! This class deals with separate caches for different
  !! units (solutions or global) and for different types of
  !< messages within the communication stage.
  type, public :: MpiMessageCacheType
    type(STLVecInt) :: cached_ids !< a vector with ids for the cached units (solution ids)
    type(ListType) :: unit_caches !< a list with caches per unit
  contains
    procedure :: init => mmc_init
    procedure :: get => mmc_get
    procedure :: put => mmc_put
    procedure :: clear => mmc_clear
    procedure :: destroy => mmc_destroy
  end type MpiMessageCacheType

contains

  !< @brief Initialize the MPI type cache system.
  !<
  subroutine mmc_init(this)
    class(MpiMessageCacheType) :: this !< the message cache

    call this%cached_ids%init()

  end subroutine mmc_init

  !< @brief Get the cached mpi datatype for the given
  !! unit, rank, stage, and message element. Returns
  !< NO_CACHED_VALUE when not in cache.
  function mmc_get(this, unit, rank, stage, msg_id) result(mpi_type)
    class(MpiMessageCacheType) :: this !< the message cache
    integer(I4B) :: unit !< the unit (solution or global)
    integer(I4B) :: rank !< the rank of the MPI process to communicate with
    integer(I4B) :: stage !< the simulation stage at which the message is sent
    integer(I4B) :: msg_id !< the message type as an integer between 1 and NR_MSG_TYPES (see above for predefined values)
    integer :: mpi_type !< the resulting mpi datatype
    ! local
    integer(I4B) :: unit_idx
    class(*), pointer :: obj_ptr

    mpi_type = NO_CACHED_VALUE

    unit_idx = this%cached_ids%get_index(unit)
    if (unit_idx == -1) return ! not cached

    obj_ptr => this%unit_caches%GetItem(unit_idx)
    select type (obj_ptr)
    class is (MpiUnitCacheType)
      mpi_type = obj_ptr%get_cached(rank, stage, msg_id)
    end select

  end function mmc_get

  !> @brief Put the mpi datatype for this particular unit,
  !! rank, and stage in cache. The datatype should be
  !< committed to the type database externally.
  subroutine mmc_put(this, unit, rank, stage, msg_id, mpi_type)
    class(MpiMessageCacheType) :: this !< the message cache
    integer(I4B) :: unit !< the unit (solution or global)
    integer(I4B) :: rank !< the rank of the MPI process to communicate with
    integer(I4B) :: stage !< the simulation stage at which the message is sent
    integer(I4B) :: msg_id !< the message type as an integer between 1 and NR_MSG_TYPES (see above for predefined values)
    integer :: mpi_type !< the mpi datatype to cache
    ! local
    integer(I4B) :: unit_idx
    type(MpiUnitCacheType), pointer :: new_cache
    class(*), pointer :: obj_ptr

    unit_idx = this%cached_ids%get_index(unit)
    if (unit_idx == -1) then
      ! add to vector with cached unit ids
      call this%cached_ids%push_back(unit)
      ! create and add unit cache
      allocate (new_cache)
      call new_cache%init(NR_SIM_STAGES, NR_MSG_TYPES)
      obj_ptr => new_cache
      call this%unit_caches%Add(obj_ptr)
      unit_idx = this%cached_ids%size
    end if

    ! get the cache for this unit
    obj_ptr => this%unit_caches%GetItem(unit_idx)
    select type (obj_ptr)
    class is (MpiUnitCacheType)
      call obj_ptr%cache(rank, stage, msg_id, mpi_type)
    end select

  end subroutine mmc_put

  !< @brief Clear the MPI type cache system
  !<
  subroutine mmc_clear(this)
    class(MpiMessageCacheType) :: this !< the message cache
    ! local
    integer(I4B) :: i
    class(*), pointer :: obj_ptr

    ! clear caches
    do i = 1, this%cached_ids%size
      obj_ptr => this%unit_caches%GetItem(i)
      select type (obj_ptr)
      class is (MpiUnitCacheType)
        call obj_ptr%clear()
      end select
    end do

  end subroutine mmc_clear

  !< @brief Destroy the MPI type cache system.
  !<
  subroutine mmc_destroy(this)
    class(MpiMessageCacheType) :: this !< the message cache
    ! local
    integer(I4B) :: i
    class(*), pointer :: obj_ptr

    ! clear caches
    do i = 1, this%cached_ids%size
      obj_ptr => this%unit_caches%GetItem(i)
      select type (obj_ptr)
      class is (MpiUnitCacheType)
        call obj_ptr%destroy()
      end select
    end do
    call this%unit_caches%Clear(destroy=.true.)

    call this%cached_ids%destroy()

  end subroutine mmc_destroy

end module
