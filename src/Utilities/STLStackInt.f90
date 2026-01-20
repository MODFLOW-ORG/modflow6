module STLStackIntModule
  use KindModule, only: I4B, LGP
  use STLVecIntModule
  use SimModule, only: ustop
  implicit none
  private
  public :: STLStackInt

  !> @brief A derived type representing a stack of integers.
  !!
  !! This type provides a stack data structure specifically for integers.
  !! It includes methods for typical stack operations such as push, pop,
  !< and checking if the stack is empty.
  type :: STLStackInt
    type(STLVecInt), private :: stack !< the internal stack
  contains
    procedure, pass(this) :: init !< allocate memory, init size and capacity
    procedure, pass(this) :: destroy !< deletes the memory
    procedure, pass(this) :: push !< adds an element at the end of the vector
    procedure, pass(this) :: pop !< removes the last element
    procedure, pass(this) :: top !< returns the last element (without removing it)
    procedure, pass(this) :: size !< returns the size of the stack
  end type STLStackInt

contains ! module routines

  subroutine init(this, capacity)
    class(STLStackInt), intent(inout) :: this
    integer(I4B), intent(in), optional :: capacity ! the initial capacity, when given

    ! init the vector
    if (present(capacity)) then
      call this%stack%init(capacity)
    else
      call this%stack%init()
    end if

  end subroutine init

  subroutine push(this, newValue)
    class(STLStackInt), intent(inout) :: this
    integer(I4B) :: newValue

    call this%stack%push_back(newValue)

  end subroutine push

  subroutine pop(this)
    class(STLStackInt), intent(inout) :: this

    if (this%stack%size == 0) then
      write (*, *) 'STLStackInt exception: cannot pop an empty stack'
      call ustop()
    end if
    this%stack%size = this%stack%size - 1

  end subroutine pop

  function top(this) result(top_value)
    class(STLStackInt), intent(in) :: this
    integer(I4B) :: top_value

    if (this%stack%size == 0) then
      write (*, *) 'STLStackInt exception: cannot get top of an empty stack'
      call ustop()
    end if
    top_value = this%stack%at(this%stack%size)

  end function top

  function size(this) result(size_value)
    class(STLStackInt), intent(in) :: this
    integer(I4B) :: size_value

    size_value = this%stack%size

  end function size

  subroutine destroy(this)
    class(STLStackInt), intent(inout) :: this

    call this%stack%destroy()

  end subroutine destroy

end module STLStackIntModule
