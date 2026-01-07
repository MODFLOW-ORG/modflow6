module LocalCellExtremaModule
  use KindModule, only: DP, I4B
  use BaseDisModule, only: DisBaseType

  implicit none
  private

  public :: LocalCellExtremaType

  !> @brief Computes and caches local extrema for each cell and its connected neighbors
  !!
  !! This class computes the minimum and maximum values within the local
  !! stencil of each cell (the cell itself plus all directly connected neighboring cells).
  !! The extrema are computed once when the scalar field is set and then cached for
  !! fast retrieval. This is particularly useful for TVD limiters and slope
  !! limiting algorithms that need to enforce monotonicity constraints.
  !!
  !! The local extrema computation follows the connectivity pattern defined by the
  !! discretization object, examining all cells that share a face with the target cell.
  !! This creates a computational stencil that includes the central cell and its
  !! immediate neighbors.
  !!
  !! @note The get_min() and get_max() methods return pointers for zero-copy performance
  !!       in performance-critical loops. Callers should treat returned values as read-only.
  !<
  type :: LocalCellExtremaType
    private
    class(DisBaseType), pointer :: dis
    real(DP), dimension(:), allocatable :: min
    real(DP), dimension(:), allocatable :: max
  contains
    procedure :: set_field
    procedure :: get_min
    procedure :: get_max
    final :: destructor

    procedure, private :: compute_local_extrema
    procedure, private :: find_local_extrema
  end type LocalCellExtremaType

  interface LocalCellExtremaType
    module procedure constructor
  end interface LocalCellExtremaType

contains
  function constructor(dis) result(local_extrema)
    ! -- return
    type(LocalCellExtremaType) :: local_extrema
    ! -- dummy
    class(DisBaseType), pointer, intent(in) :: dis

    local_extrema%dis => dis

    allocate (local_extrema%min(dis%nodes))
    allocate (local_extrema%max(dis%nodes))
  end function constructor

  subroutine destructor(this)
    ! -- dummy
    type(LocalCellExtremaType), intent(inout) :: this

    deallocate (this%min)
    deallocate (this%max)
  end subroutine destructor

  subroutine set_field(this, phi)
    ! -- dummy
    class(LocalCellExtremaType), target :: this
    real(DP), intent(in), dimension(:), pointer :: phi

    call this%compute_local_extrema(phi)
  end subroutine set_field

  function get_min(this, n) result(min_val)
    ! -- dummy
    class(LocalCellExtremaType), target :: this
    integer(I4B), intent(in) :: n ! Node index
    !-- return
    real(DP), pointer :: min_val

    min_val => this%min(n)
  end function get_min

  function get_max(this, n) result(max_val)
    ! -- dummy
    class(LocalCellExtremaType), target :: this
    integer(I4B), intent(in) :: n ! Node index
    !-- return
    real(DP), pointer :: max_val

    max_val => this%max(n)
  end function get_max

  subroutine compute_local_extrema(this, phi)
    ! -- dummy
    class(LocalCellExtremaType), target :: this
    real(DP), intent(in), dimension(:) :: phi
    ! -- local
    integer(I4B) :: n
    real(DP) :: min_phi, max_phi

    do n = 1, this%dis%nodes
      call this%find_local_extrema(n, phi, min_phi, max_phi)
      this%min(n) = min_phi
      this%max(n) = max_phi
    end do
  end subroutine compute_local_extrema

  subroutine find_local_extrema(this, n, phi, min_phi, max_phi)
    ! -- dummy
    class(LocalCellExtremaType), target :: this
    integer(I4B), intent(in) :: n
    real(DP), intent(in), dimension(:) :: phi
    real(DP), intent(out) :: min_phi, max_phi
    ! -- local
    integer(I4B) :: ipos, m

    min_phi = phi(n)
    max_phi = phi(n)
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)
      min_phi = min(min_phi, phi(m))
      max_phi = max(max_phi, phi(m))
    end do
  end subroutine

end module LocalCellExtremaModule
