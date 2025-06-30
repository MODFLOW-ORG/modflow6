module LeastSquaredGradientModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE

  Use IGradient
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use PseudoInverseModule, only: pinv
  use DisUtilsModule, only: number_connected_faces, node_distance

  implicit none
  private

  public :: LeastSquaredGradientType

  type Array2D
    real(DP), dimension(:, :), allocatable :: data
  end type Array2D

  !> @brief Weighted least-squares gradient method for structured and unstructured grids.
  !!
  !! This class implements a least-squares gradient reconstruction for use on both structured and unstructured grids.
  !! For each cell, it precomputes and caches a gradient operator using the Moore-Penrose pseudoinverse,
  !! based on the geometry and connectivity of the mesh. The operator is created once during initialization
  !! and can then be efficiently applied to any scalar field to compute the gradient in each cell.
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
  type, extends(IGradientType) :: LeastSquaredGradientType
    class(DisBaseType), pointer :: dis
    type(TspFmiType), pointer :: fmi
    type(Array2D), allocatable, dimension(:) :: grad_op
  contains
    procedure :: get

    procedure, private :: compute_cell_gradient
    procedure, private :: create_grad_operator
  end type LeastSquaredGradientType

  interface LeastSquaredGradientType
    module procedure Constructor
  end interface LeastSquaredGradientType

contains
  function constructor(dis, fmi) Result(gradient)
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    type(TspFmiType), pointer, intent(in) :: fmi
    !-- return
    type(LeastSquaredGradientType) :: gradient
    ! -- local
    integer(I4B) :: n, nodes

    gradient%dis => dis
    gradient%fmi => fmi

    nodes = dis%nodes

    ! -- Compute the gradient operator
    nodes = dis%nodes
    allocate (gradient%grad_op(dis%nodes))
    do n = 1, nodes
      gradient%grad_op(n)%data = gradient%create_grad_operator(n)
    end do
  end function constructor

  function create_grad_operator(this, n) result(grad_op)
    ! -- dummy
    class(LeastSquaredGradientType) :: this
    integer(I4B), intent(in) :: n ! Cell index for which to create the operator
    real(DP), dimension(:, :), allocatable :: grad_op ! The resulting gradient operator (3 x number_connections)
    ! -- local
    integer(I4B) :: number_connections ! Number of connected neighboring cells
    integer(I4B) :: ipos, local_pos, m ! Loop indices and neighbor cell index
    real(DP) :: length ! Distance between cell centers
    real(DP), dimension(3) :: dnm ! Vector from cell n to neighbor m
    real(DP), dimension(:, :), allocatable :: d ! Matrix of normalized direction vectors (number_connections x 3)
    real(DP), dimension(:, :), allocatable :: d_trans ! Transpose of d (3 x number_connections)
    real(DP), dimension(:, :), allocatable :: grad_scale ! Diagonal scaling matrix (number_connections x number_connections),
    ! where each diagonal entry is the inverse of the distance between
    ! the cell center and its neighbor.
    real(DP), dimension(3, 3) :: g
    real(DP), dimension(3, 3) :: g_inv

    number_connections = number_connected_faces(this%dis, n)

    allocate (d(number_connections, 3))
    allocate (d_trans(3, number_connections))
    allocate (grad_op(3, number_connections))
    allocate (grad_scale(number_connections, number_connections))

    grad_scale = 0
    d = 0

    ! Assemble the distance matrix
    ! Handle the internal connections
    local_pos = 1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)

      dnm = node_distance(this%dis, this%fmi, n, m)
      length = norm2(dnm)

      d(local_pos, :) = dnm / length
      grad_scale(local_pos, local_pos) = 1.0_dp / length

      local_pos = local_pos + 1
    end do

    d_trans = transpose(d)

    ! Compute the G and inverse G matrices
    g = matmul(d_trans, d)
    g_inv = pinv(g)

    ! Compute the gradient operator
    grad_op = matmul(matmul(g_inv, d_trans), grad_scale)

  end function create_grad_operator

  function get(this, n, c) result(grad_c)
    ! -- dummy
    class(LeastSquaredGradientType), target :: this
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
    class(LeastSquaredGradientType), target :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: phi_new
    ! -- local
    real(DP), dimension(:, :), pointer :: grad_op
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
    grad_op => this%grad_op(n)%data
    grad_c = matmul(grad_op, dc)

  end function compute_cell_gradient

end module LeastSquaredGradientModule
