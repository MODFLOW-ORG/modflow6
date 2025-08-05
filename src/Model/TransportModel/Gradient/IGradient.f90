module IGradient
  use KindModule, only: DP, I4B

  implicit none
  private

  public :: IGradientType

  !> @brief Abstract interface for cell-based gradient computation.
  !!
  !! This module defines the abstract type `IGradientType`, which provides a deferred
  !! interface for computing the gradient of a scalar field at a given cell.
  !! Any concrete gradient implementation must extend this type and implement the `get` method.
  !<
  type, abstract :: IGradientType
  contains
    procedure(get_if), deferred :: get
  end type IGradientType

  abstract interface

    function get_if(this, n, c) result(grad_c)
      ! -- import
      import IGradientType
      import DP, I4B
      ! -- dummy
      class(IGradientType), target :: this
      integer(I4B), intent(in) :: n
      real(DP), dimension(:), intent(in) :: c
      !-- return
      real(DP), dimension(3) :: grad_c
    end function

  end interface

contains

end module IGradient
