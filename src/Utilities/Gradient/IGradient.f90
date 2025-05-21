module IGradient
  use KindModule, only: DP, I4B

  implicit none
  private

  public :: IGradientType

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

    ! subroutine compute_cell_gradient_if(this)
    !   import IGradientType
    !   class(IGradientType) :: this
    ! end subroutine

  end interface

contains

end module IGradient
