module CircularBufferModule

  use KindModule, only: DP, I4B
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use SimModule, only: store_error

  implicit none
  private
  public :: CircularBufferType

  type :: CircularBufferType
    private
    integer(I4B) :: n_capacity !< Maximum number of elements in the buffer
    integer(I4B) :: n_size = 0 !< Current number of elements in the buffer
    integer(I4B) :: head_index = 0 !< Slot index of the back (newest) element; 0 when empty
    real(DP), pointer, contiguous :: data(:, :) !< The buffer array

  contains
    procedure :: push_back
    procedure :: pop_back
    generic :: front => front_impl, front_at_impl
    generic :: back => back_impl, back_at_impl
    procedure :: size => get_size
    procedure :: capacity => get_capacity
    procedure :: is_empty
    procedure :: is_full
    final :: destructor

    procedure, private :: front_impl
    procedure, private :: front_at_impl
    procedure, private :: back_impl
    procedure, private :: back_at_impl
  end type CircularBufferType

  interface CircularBufferType
    module procedure constructor
  end interface CircularBufferType

contains

  function constructor(capacity, neq, name, mem_path) result(buffer)
    type(CircularBufferType) :: buffer
    ! -- dummy
    integer(I4B), intent(in) :: capacity !< maximum number of elements
    integer(I4B), intent(in) :: neq !< number of equations (size of each element)
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored

    call mem_allocate(buffer%data, neq, capacity, name, mem_path)
    buffer%n_capacity = capacity

  end function constructor

  subroutine destructor(this)
    ! -- dummy
    type(CircularBufferType), intent(inout) :: this

    call mem_deallocate(this%data)
  end subroutine destructor

  !> @brief Insert an element at the back (newest end) of the buffer.
  !!
  !! If the buffer is full the front (oldest) element is overwritten.
  !<
  subroutine push_back(this, element)
    ! -- dummy
    class(CircularBufferType), target :: this
    real(DP), dimension(:), intent(in) :: element
    ! -- local
    integer(I4B) :: insert_index

    if (this%n_size < this%n_capacity) then
      this%n_size = this%n_size + 1
    end if

    insert_index = mod(this%head_index, this%n_capacity) + 1
    this%data(:, insert_index) = element
    this%head_index = insert_index

  end subroutine push_back

  !> @brief Remove the element at the back (newest end) of the buffer.
  !!
  !! Error stops if the buffer is empty.
  !<
  subroutine pop_back(this)
    ! -- dummy
    class(CircularBufferType), target :: this

    if (this%n_size == 0) then
      call store_error("Cannot pop_back from an empty CircularBufferType")
      return
    end if

    this%n_size = this%n_size - 1
    if (this%n_size > 0) then
      this%head_index = mod(this%head_index - 2 + this%n_capacity, &
                            this%n_capacity) + 1
    else
      this%head_index = 0
    end if

  end subroutine pop_back

  !> @brief Return a pointer to the front (oldest) element.
  !!
  !! Error stops if the buffer is empty.
  !<
  function front_impl(this) result(element)
    ! -- dummy
    class(CircularBufferType), target :: this
    ! -- return
    real(DP), pointer, dimension(:) :: element

    element => this%front_at_impl(1)

  end function front_impl

  !> @brief Return a pointer to element at 1-based index n counted from the front.
  !!
  !! n=1 is the oldest (front) element; n=size() is the newest (back).
  !! Error stops if n is out of bounds.
  !<
  function front_at_impl(this, n) result(element)
    ! -- dummy
    class(CircularBufferType), target :: this
    integer(I4B), intent(in) :: n !< 1-based index from front (1=oldest, size()=newest)
    ! -- return
    real(DP), pointer, dimension(:) :: element
    ! -- local
    integer(I4B) :: index

    if (n < 1 .or. n > this%n_size) then
      call store_error("Index out of bounds in CircularBufferType%front")
      element => null()
      return
    end if

    index = mod(this%head_index - this%n_size + n - 1 + this%n_capacity, &
                this%n_capacity) + 1
    element => this%data(:, index)

  end function front_at_impl

  !> @brief Return a pointer to the back (newest) element.
  !!
  !! Error stops if the buffer is empty.
  !<
  function back_impl(this) result(element)
    ! -- dummy
    class(CircularBufferType), target :: this
    ! -- return
    real(DP), pointer, dimension(:) :: element

    element => this%back_at_impl(1)

  end function back_impl

  !> @brief Return a pointer to element at 1-based index n counted from the back.
  !!
  !! n=1 is the newest (back) element; n=size() is the oldest (front).
  !! Error stops if n is out of bounds.
  !<
  function back_at_impl(this, n) result(element)
    ! -- dummy
    class(CircularBufferType), target :: this
    integer(I4B), intent(in) :: n !< 1-based index from back (1=newest, size()=oldest)
    ! -- return
    real(DP), pointer, dimension(:) :: element
    ! -- local
    integer(I4B) :: index

    if (n < 1 .or. n > this%n_size) then
      call store_error("Index out of bounds in CircularBufferType%back")
      element => null()
      return
    end if

    index = mod(this%head_index - n + this%n_capacity, this%n_capacity) + 1
    element => this%data(:, index)

  end function back_at_impl

  !> @brief Return the number of elements currently stored in the buffer.
  !<
  function get_size(this) result(n)
    ! -- dummy
    class(CircularBufferType), target :: this
    ! -- return
    integer(I4B) :: n

    n = this%n_size

  end function get_size

  !> @brief Return the maximum number of elements the buffer can hold.
  !<
  function get_capacity(this) result(n)
    ! -- dummy
    class(CircularBufferType), target :: this
    ! -- return
    integer(I4B) :: n

    n = this%n_capacity

  end function get_capacity

  !> @brief Return .true. if the buffer contains no elements.
  !<
  function is_empty(this) result(empty)
    ! -- dummy
    class(CircularBufferType), target :: this
    ! -- return
    logical :: empty

    empty = (this%n_size == 0)

  end function is_empty

  !> @brief Return .true. if the number of elements equals the capacity.
  !<
  function is_full(this) result(full)
    ! -- dummy
    class(CircularBufferType), target :: this
    ! -- return
    logical :: full

    full = (this%n_size == this%n_capacity)

  end function is_full

end module CircularBufferModule
