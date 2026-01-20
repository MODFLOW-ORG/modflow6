module STLVecIntModule
  use KindModule, only: I4B, LGP
  use SimModule, only: ustop
  use ArrayHandlersModule, only: ExpandArray
  implicit none
  private
  public :: STLVecInt

  integer(I4B), parameter :: defaultInitialCapacity = 4

  ! This is a dynamic vector type for integers
  type :: STLVecInt
    integer(I4B), private, allocatable :: values(:) !< the internal array for storage
    integer(I4B) :: size !< the number of elements
    integer(I4B), private :: capacity !< the reserved storage
  contains
    procedure, pass(this) :: init !< allocate memory, init size and capacity
    procedure, pass(this) :: push_back !< adds an element at the end of the vector
    procedure, pass(this) :: push_back_unique !< adds an element at the end of the vector, if not present yet
    procedure, pass(this) :: pop !< removes the last element
    procedure, pass(this) :: add_array !< adds elements of array at the end of the vector
    procedure, pass(this) :: add_array_unique !< adds elements of array at the end of the vector, if not present yet
    procedure, pass(this) :: at !< random access, unsafe, no bounds checking
    procedure, pass(this) :: at_safe !< random access with bounds checking
    procedure, pass(this) :: set !< set value at index, no bounds checking
    procedure, pass(this) :: clear !< empties the vector, leaves memory unchanged
    procedure, pass(this) :: shrink_to_fit !< reduces the allocated memory to fit the actual vector size
    procedure, pass(this) :: destroy !< deletes the memory
    procedure, pass(this) :: contains !< true when element already present
    procedure, pass(this) :: get_index !< return index of first occurrence of value in array, -1 when not present
    procedure, pass(this) :: get_values !< returns a copy of the values
    ! private
    procedure, private, pass(this) :: expand
  end type STLVecInt

contains ! module routines

  subroutine init(this, capacity)
    class(STLVecInt), intent(inout) :: this
    integer(I4B), intent(in), optional :: capacity ! the initial capacity, when given

    if (present(capacity)) then
      this%capacity = capacity
    else
      this%capacity = defaultInitialCapacity
    end if

    allocate (this%values(this%capacity))
    this%size = 0

  end subroutine init

  subroutine push_back(this, newValue)
    class(STLVecInt), intent(inout) :: this
    integer(I4B) :: newValue
    ! check capacity
    if (this%size + 1 > this%capacity) then
      call this%expand()
    end if

    this%size = this%size + 1
    this%values(this%size) = newValue

  end subroutine push_back

  subroutine push_back_unique(this, newValue)
    class(STLVecInt), intent(inout) :: this
    integer(I4B) :: newValue

    if (.not. this%contains(newValue)) then
      call this%push_back(newValue)
    end if

  end subroutine push_back_unique

  subroutine pop(this)
    class(STLVecInt), intent(inout) :: this

    if (this%size > 0) then
      this%size = this%size - 1
    else
      write (*, *) 'STLVecInt exception: cannot pop from an empty array'
      call ustop()
    end if

  end subroutine

  subroutine add_array(this, array)
    class(STLVecInt), intent(inout) :: this
    integer(I4B), dimension(:), pointer :: array
    ! local
    integer(I4B) :: i

    do i = 1, size(array)
      call this%push_back(array(i))
    end do

  end subroutine add_array

  subroutine add_array_unique(this, array)
    class(STLVecInt), intent(inout) :: this
    integer(I4B), dimension(:), pointer :: array
    ! local
    integer(I4B) :: i

    do i = 1, size(array)
      if (.not. this%contains(array(i))) then
        call this%push_back(array(i))
      end if
    end do

  end subroutine add_array_unique

  function at(this, idx) result(value)
    class(STLVecInt), intent(in) :: this
    integer(I4B), intent(in) :: idx
    integer(I4B) :: value

    value = this%values(idx)

  end function at

  function at_safe(this, idx) result(value)
    class(STLVecInt), intent(inout) :: this
    integer(I4B), intent(in) :: idx
    integer(I4B) :: value

    if (idx > this%size) then
      write (*, *) 'STLVecInt exception: access out of bounds, index ', idx, &
        ' exceeds actual size (', this%size, ')'
      call ustop()
    end if
    value = this%at(idx)

  end function at_safe

  subroutine set(this, idx, value)
    class(STLVecInt), intent(inout) :: this
    integer(I4B), intent(in) :: idx
    integer(I4B) :: value

    this%values(idx) = value

  end subroutine set

  subroutine clear(this)
    class(STLVecInt), intent(inout) :: this

    ! really, this is all there is to it...
    this%size = 0

  end subroutine clear

  subroutine shrink_to_fit(this)
    class(STLVecInt), intent(inout) :: this
    ! local
    integer(I4B), allocatable :: tempValues(:)
    integer(I4B) :: i, newSize

    if (this%size == this%capacity) then
      return
    end if

    ! store temp
    newSize = this%size
    allocate (tempValues(newSize))
    do i = 1, newSize
      tempValues(i) = this%values(i)
    end do

    ! reinit
    call this%destroy()
    call this%init(newSize)

    ! copy back
    do i = 1, newSize
      call this%push_back(tempValues(i))
    end do

  end subroutine shrink_to_fit

  subroutine destroy(this)
    class(STLVecInt), intent(inout) :: this

    if (allocated(this%values)) then
      deallocate (this%values)
      this%size = 0
      this%capacity = 0
    else
      write (*, *) 'STLVecInt exception: cannot delete an unallocated array'
      call ustop()
    end if

  end subroutine destroy

  ! expand the array with the given strategy, at
  ! least by 1
  subroutine expand(this)
    class(STLVecInt), intent(inout) :: this
    integer(I4B) :: increment

    ! expansion strategy
    increment = this%capacity / 2 + 1
    call ExpandArray(this%values, increment)
    this%capacity = this%capacity + increment

  end subroutine expand

  ! check if the element is already present
  function contains(this, val) result(res)
    class(STLVecInt), intent(inout) :: this
    integer(I4B) :: val
    logical(LGP) :: res
    ! local
    integer(I4B) :: i

    res = .false.
    do i = 1, this%size
      if (this%at(i) == val) then
        res = .true.
        return
      end if
    end do

  end function contains

  !> @brief Return index of first occurrence,
  !< returns -1 when not present
  function get_index(this, val) result(idx)
    class(STLVecInt), intent(inout) :: this
    integer(I4B) :: val
    integer(I4B) :: idx
    ! local
    integer(I4B) :: i

    idx = -1
    do i = 1, this%size
      if (this%at(i) == val) then
        idx = i
        return
      end if
    end do

  end function get_index

  function get_values(this) result(values)
    class(STLVecInt), intent(in) :: this
    integer(I4B), dimension(:), allocatable :: values

    values = this%values(1:this%size)

  end function get_values

end module STLVecIntModule
