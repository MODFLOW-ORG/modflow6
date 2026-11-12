module VectorBaseModule
  use KindModule, only: I4B, DP, LGP
  implicit none
  private

  type, public, abstract :: VectorBaseType
    logical(LGP) :: is_mem_managed
  contains
    procedure(create_mm_if), deferred :: create_mm
    procedure(destroy_if), deferred :: destroy
    procedure(get_array_if), deferred :: get_array
    procedure(get_ownership_range_if), deferred :: get_ownership_range
    procedure(get_size_if), deferred :: get_size
    procedure(get_value_local_if), deferred :: get_value_local
    procedure(zero_entries_if), deferred :: zero_entries
    procedure(set_value_local_if), deferred :: set_value_local
    procedure(axpy_if), deferred :: axpy
    procedure(norm2_if), deferred :: norm2
    procedure(print_if), deferred :: print
  end type VectorBaseType

  abstract interface
    subroutine create_mm_if(this, n, name, mem_path)
      import VectorBaseType, I4B
      class(VectorBaseType) :: this !< this vector
      integer(I4B) :: n !< the nr. of elements in the (local) vector
      character(len=*) :: name !< the variable name (for access through memory manager)
      character(len=*) :: mem_path !< memory path for storing the underlying memory items
    end subroutine create_mm_if
    subroutine destroy_if(this)
      import VectorBaseType
      class(VectorBaseType) :: this !< this vector
    end subroutine destroy_if
    function get_array_if(this) result(array)
      import VectorBaseType, DP
      class(VectorBaseType) :: this !< this vector
      real(DP), dimension(:), pointer, contiguous :: array
    end function get_array_if
    subroutine get_ownership_range_if(this, start, end)
      import VectorBaseType, I4B
      class(VectorBaseType) :: this !< this vector
      integer(I4B) :: start, end !< local range in global numbering
    end subroutine
    !> @return the global vector size
    function get_size_if(this) result(size)
      import VectorBaseType, I4B
      class(VectorBaseType) :: this !< this vector
      integer(I4B) :: size
    end function
    !> @return the value at the index
    function get_value_local_if(this, idx) result(val)
      import VectorBaseType, I4B, DP
      class(VectorBaseType) :: this !< this vector
      integer(I4B) :: idx !< index in local numbering
      real(DP) :: val
    end function
    subroutine zero_entries_if(this)
      import VectorBaseType
      class(VectorBaseType) :: this
    end subroutine
    subroutine set_value_local_if(this, idx, val)
      import VectorBaseType, I4B, DP
      class(VectorBaseType) :: this !< this vector
      integer(I4B) :: idx !< index in local numbering
      real(DP) :: val !< the value to set
    end subroutine
    subroutine axpy_if(this, alpha, vec_x)
      import VectorBaseType, DP
      class(VectorBaseType) :: this !< this vector y => alpha*x + y
      real(DP) :: alpha !< the factor
      class(VectorBaseType), pointer :: vec_x !< the vector to add
    end subroutine
    !> @return the calculated 2-norm
    function norm2_if(this) result(n2)
      import VectorBaseType, DP
      class(VectorBaseType) :: this !< this vector
      real(DP) :: n2
    end function
    subroutine print_if(this)
      import VectorBaseType
      class(VectorBaseType) :: this
    end subroutine print_if
  end interface

end module VectorBaseModule
