module LeastSquaresGradientModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE

  Use IGradient
  use BaseDisModule, only: DisBaseType
  use PseudoInverseModule, only: pinv
  use DisUtilsModule, only: number_connected_faces, node_distance

  implicit none
  private

  public :: LeastSquaresGradientType

  type Array2D
    real(DP), dimension(:, :), allocatable :: data
  end type Array2D

  !> @brief Weighted least-squares gradient method for structured and unstructured grids.
  !!
  !! This class implements a least-squares gradient reconstruction for use on both structured and unstructured grids.
  !! For each cell, it precomputes and caches a gradient reconstruction matrix using the Moore-Penrose pseudoinverse,
  !! based on the geometry and connectivity of the mesh. The operator is created once during initialization
  !! and can then be efficiently applied to any scalar field to compute the gradient in each cell.
  !! The gradient can then be computed by multiplying the reconstruction matrix with the difference vector.
  !! ∇ɸ = R * ∑(ɸ_i - ɸ_up), where i are the neighboring cells.
  !!
  !! - The gradient operator is constructed using normalized direction vectors between cell centers,
  !!   scaled by the inverse of the distance.
  !! - The least-squares approach ensures robust gradients even for irregular or rank-deficient stencils.
  !! - The operator is cached for each cell, so gradient computation is efficient for repeated queries.
  !! - The class provides a `get` method to compute the gradient for any cell and scalar field.
  !!
  !! @note Boundary cells are not handled in a special manner. This may impact the quality of the gradient
  !!       near boundaries, especially if a cell does not have enough neighbors (fewer than three in 3D).
  !<
  type, extends(IGradientType) :: LeastSquaresGradientType
    class(DisBaseType), pointer :: dis
    type(Array2D), allocatable, dimension(:) :: R ! Gradient reconstruction matrix
  contains
    procedure :: get

    procedure, private :: compute_cell_gradient
    procedure, private :: create_gradient_reconstruction_matrix
  end type LeastSquaresGradientType

  interface LeastSquaresGradientType
    module procedure Constructor
  end interface LeastSquaresGradientType

contains
  function constructor(dis) Result(gradient)
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    !-- return
    type(LeastSquaresGradientType) :: gradient
    ! -- local
    integer(I4B) :: n, nodes

    gradient%dis => dis
    nodes = dis%nodes

    ! -- Compute the gradient rec
    nodes = dis%nodes
    allocate (gradient%R(dis%nodes))
    do n = 1, nodes
      gradient%R(n)%data = gradient%create_gradient_reconstruction_matrix(n)
    end do
  end function constructor

  function create_gradient_reconstruction_matrix(this, n) result(R)
    ! -- dummy
    class(LeastSquaresGradientType) :: this
    integer(I4B), intent(in) :: n ! Cell index for which to create the operator
    real(DP), dimension(:, :), allocatable :: R ! The resulting gradient reconstruction matrix (3 x number_connections)
    ! -- local
    integer(I4B) :: number_connections ! Number of connected neighboring cells
    integer(I4B) :: ipos, local_pos, m ! Loop indices and neighbor cell index
    real(DP) :: length ! Distance between cell centers
    real(DP), dimension(3) :: dnm ! Vector from cell n to neighbor m
    real(DP), dimension(:, :), allocatable :: d ! Matrix of normalized direction vectors (number_connections x 3)
    real(DP), dimension(:, :), allocatable :: grad_scale ! Diagonal scaling matrix (number_connections x number_connections),
    ! where each diagonal entry is the inverse of the distance between

    number_connections = number_connected_faces(this%dis, n)

    allocate (d(number_connections, 3))
    allocate (R(3, number_connections))
    allocate (grad_scale(number_connections, number_connections))

    grad_scale = 0
    d = 0

    ! Assemble the distance matrix
    ! Handle the internal connections
    local_pos = 1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)

      dnm = node_distance(this%dis, n, m)
      length = norm2(dnm)

      d(local_pos, :) = dnm / length
      grad_scale(local_pos, local_pos) = 1.0_dp / length

      local_pos = local_pos + 1
    end do

    ! Compute the gradient reconstructions matrix
    R = matmul(pinv(d), grad_scale)

  end function create_gradient_reconstruction_matrix

  function get(this, n, c) result(grad_c)
    ! -- dummy
    class(LeastSquaresGradientType), target :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: c
    !-- return
    real(DP), dimension(3) :: grad_c

    grad_c = this%compute_cell_gradient(n, c)
  end function get

  function compute_cell_gradient(this, n, phi_new) result(grad_c)
    ! -- return
    real(DP), dimension(3) :: grad_c
    ! -- dummy
    class(LeastSquaresGradientType), target :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: phi_new
    ! -- local
    real(DP), dimension(:, :), pointer :: R
    integer(I4B) :: ipos, local_pos
    integer(I4B) :: number_connections

    integer(I4B) :: m
    real(DP), dimension(:), allocatable :: dc

    ! Assemble the concentration difference vector
    number_connections = number_connected_faces(this%dis, n)
    allocate (dc(number_connections))
    local_pos = 1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)
      dc(local_pos) = phi_new(m) - phi_new(n)
      local_pos = local_pos + 1
    end do

    ! Compute the cells gradient
    R => this%R(n)%data
    grad_c = matmul(R, dc)

  end function compute_cell_gradient

end module LeastSquaresGradientModule
