module CachedGradientModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE

  Use IGradient
  use BaseDisModule, only: DisBaseType

  implicit none
  private

  public :: CachedGradientType

  !> @brief Decorator that adds caching to any gradient computation implementation
  !!
  !! This class wraps any IGradientType implementation and provides caching functionality
  !! using the Decorator pattern. When set_field is called, it pre-computes gradients
  !! for all nodes and stores them in memory for fast O(1) retrieval. This trades memory
  !! for speed when gradients are accessed multiple times for the same scalar field.
  !!
  !! The class takes ownership of the wrapped gradient object via move semantics and
  !! provides the same interface as any other gradient implementation. This allows it
  !! to be used transparently in place of the original gradient object.
  !!
  !! Usage pattern:
  !! 1. Create a base gradient implementation (e.g., LeastSquaresGradientType)
  !! 2. Wrap it with CachedGradientType using the constructor
  !! 3. Call set_field() once to pre-compute all gradients
  !! 4. Call get() multiple times for fast cached lookups
  !!
  !! @note The wrapped gradient object is moved (not copied) during construction
  !!       for efficient memory management.
  !<
  type, extends(IGradientType) :: CachedGradientType
    private
    class(DisBaseType), pointer :: dis
    class(IGradientType), allocatable :: gradient
    real(DP), dimension(:, :), allocatable :: cached_gradients ! gradients at nodes
  contains
    procedure :: get
    procedure :: set_field
    final :: destructor
  end type CachedGradientType

  interface CachedGradientType
    module procedure Constructor
  end interface CachedGradientType

contains

  function constructor(gradient, dis) result(cached_gradient)
    ! --dummy
    class(IGradientType), allocatable, intent(inout) :: gradient
    class(DisBaseType), pointer, intent(in) :: dis
    !-- return
    type(CachedGradientType) :: cached_gradient

    cached_gradient%dis => dis

    call move_alloc(gradient, cached_gradient%gradient) ! Take ownership
    allocate (cached_gradient%cached_gradients(dis%nodes, 3))

  end function constructor

  subroutine destructor(this)
    ! -- dummy
    type(CachedGradientType), intent(inout) :: this

    deallocate (this%cached_gradients)
  end subroutine destructor

  function get(this, n) result(grad_c)
    ! -- dummy
    class(CachedGradientType), target :: this
    integer(I4B), intent(in) :: n
    !-- return
    real(DP), dimension(3) :: grad_c

    grad_c = this%cached_gradients(n, :)
  end function get

  subroutine set_field(this, phi)
    ! -- dummy
    class(CachedGradientType), target :: this
    real(DP), dimension(:), pointer, intent(in) :: phi
    ! -- local
    integer(I4B) :: n

    call this%gradient%set_field(phi)
    do n = 1, this%dis%nodes
      this%cached_gradients(n, :) = this%gradient%get(n)
    end do

  end subroutine set_field

end module CachedGradientModule
