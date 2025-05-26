module LeastSquaredGradientBoundaryModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE

  Use IGradient
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use BoundaryFacesModule, only: BoundaryFacesType
  use SVDModule, only: pinv
  use DisInfoModule, only: number_connected_faces, number_faces

  implicit none
  private

  public :: LeastSquaredGradientBoundaryType

  type Array2D
    real(DP), dimension(:, :), allocatable :: data
  end type Array2D

  type, extends(IGradientType) :: LeastSquaredGradientBoundaryType
    private
    class(DisBaseType), pointer :: dis
    type(TspFmiType), pointer :: fmi
    type(Array2D), allocatable, dimension(:) :: grad_op
    type(BoundaryFacesType), allocatable :: boundary_faces
  contains
    procedure :: get

    procedure, private :: compute_cell_gradient
    procedure, private :: create_grad_operator
    procedure, private :: node_distance
  end type LeastSquaredGradientBoundaryType

  interface LeastSquaredGradientBoundaryType
    module procedure Constructor
  end interface LeastSquaredGradientBoundaryType

contains
  function constructor(dis, fmi) Result(gradient)
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    type(TspFmiType), pointer, intent(in) :: fmi
    !-- return
    type(LeastSquaredGradientBoundaryType) :: gradient
    ! -- local
    integer(I4B) :: n, nodes

    gradient%dis => dis
    gradient%fmi => fmi

    nodes = dis%nodes

    ! -- Create boundary Cells
    gradient%boundary_faces = BoundaryFacesType(dis)

    ! -- Compute the gradient operator
    nodes = dis%nodes
    allocate (gradient%grad_op(dis%nodes))
    do n = 1, nodes
      gradient%grad_op(n)%data = gradient%create_grad_operator(n)
    end do
  end function constructor

  function create_grad_operator(this, n) result(grad_op)
    ! -- dummy
    class(LeastSquaredGradientBoundaryType) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:, :), allocatable :: grad_op
    ! -- local
    integer(I4B) :: number_connections
    integer(I4B) :: ipos, local_pos, m
    real(DP) :: length
    real(DP), dimension(3) :: dnm
    real(DP), dimension(:, :), allocatable :: d
    real(DP), dimension(:, :), allocatable :: d_trans
    real(DP), dimension(:, :), allocatable :: grad_scale
    real(DP), dimension(3, 3) :: g
    real(DP), dimension(3, 3) :: g_inv
    integer(I4B) :: number_boundaries, number_sides

    number_connections = number_connected_faces(this%dis, n)
    number_boundaries = this%boundary_faces%ia(n + 1) - this%boundary_faces%ia(n)
    number_sides = number_connections + number_boundaries

    allocate (d(number_sides, 3))
    allocate (d_trans(3, number_sides))
    allocate (grad_op(3, number_sides))
    allocate (grad_scale(number_sides, number_connections))

    grad_scale = 0
    d = 0

    ! Assemble the distance matrix
    ! Handle the internal connections
    local_pos = 1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)

      dnm = this%node_distance(n, m)
      length = norm2(dnm)

      d(local_pos, :) = dnm / length
      grad_scale(local_pos, local_pos) = 1.0_dp / length

      local_pos = local_pos + 1
    end do

    ! Handle the boundary cells
    do ipos = this%boundary_faces%ia(n), this%boundary_faces%ia(n + 1) - 1
      d(local_pos, :) = this%boundary_faces%get_normal(ipos)
      local_pos = local_pos + 1
    end do

    d_trans = transpose(d)

    ! Compute the G and inverse G matrices
    g = matmul(d_trans, d)
    g_inv = pinv(g)

    ! Compute the gradient operator
    grad_op = matmul(matmul(g_inv, d_trans), grad_scale)

  end function create_grad_operator

  function node_distance(this, n, m) result(d)
    ! -- return
    real(DP), dimension(3) :: d
    ! -- dummy
    class(LeastSquaredGradientBoundaryType) :: this
    integer(I4B), intent(in) :: n, m
    ! -- local
    real(DP) :: x_dir, y_dir, z_dir, length
    real(DP) :: satn, satm
    integer(I4B) :: ipos, isympos, ihc
    real(DP), dimension(3) :: xn, xm

    isympos = -1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      if (this%dis%con%ja(ipos) == m) then
        isympos = this%dis%con%jas(ipos)
        exit
      end if
    end do

    if (isympos == -1) then
      ! -- if the connection is not found, then return the distance between the two nodes
      xn(1) = this%dis%xc(n)
      xn(2) = this%dis%yc(n)
      xn(3) = (this%dis%top(n) + this%dis%bot(n)) / 2.0_dp

      xm(1) = this%dis%xc(m)
      xm(2) = this%dis%yc(m)
      xm(3) = (this%dis%top(m) + this%dis%bot(m)) / 2.0_dp

      d = xm - xn
      return
    end if

    ihc = this%dis%con%ihc(isympos)
    if (associated(this%fmi%gwfsat)) then
      satn = this%fmi%gwfsat(n)
      satm = this%fmi%gwfsat(m)
    else
      satn = DONE
      satm = DONE
    end if

    call this%dis%connection_vector(n, m, .true., satn, satm, ihc, x_dir, &
                                    y_dir, z_dir, length)
    d(1) = x_dir * length
    d(2) = y_dir * length
    d(3) = z_dir * length

  end function node_distance

  function get(this, n, c) result(grad_c)
    ! -- dummy
    class(LeastSquaredGradientBoundaryType), target :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: c
    !-- return
    real(DP), dimension(3) :: grad_c

    grad_c = this%compute_cell_gradient(n, c)
  end function get

  function compute_cell_gradient(this, n, cnew) result(grad_c)
    ! -- return
    real(DP), dimension(3) :: grad_c
    ! -- dummy
    class(LeastSquaredGradientBoundaryType), target :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: cnew
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
      dc(local_pos) = cnew(m) - cnew(n)
      local_pos = local_pos + 1
    end do

    ! Compute the cells gradient
    grad_op => this%grad_op(n)%data
    grad_c = matmul(grad_op, dc)

  end function compute_cell_gradient

end module LeastSquaredGradientBoundaryModule
