module GreenGaussGradientModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE

  Use IGradient
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType

  use DisInfoModule, only: number_connected_faces
  use BoundaryFacesModule, only: BoundaryFacesType

  implicit none
  private

  public :: GreenGaussGradientType

  type Array2D
    real(DP), dimension(:, :), allocatable :: data
  end type Array2D

  type, extends(IGradientType) :: GreenGaussGradientType
    private
    class(DisBaseType), pointer :: dis
    type(TspFmiType), pointer :: fmi
    type(Array2D), allocatable, dimension(:) :: grad_op
    type(BoundaryFacesType), allocatable :: boundary_faces
  contains
    procedure :: get

    procedure, private :: create_grad_operator
    procedure, private :: compute_cell_gradient
  end type GreenGaussGradientType

  interface GreenGaussGradientType
    module procedure Constructor
  end interface GreenGaussGradientType

contains

  function constructor(dis, fmi) Result(gradient)
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    type(TspFmiType), pointer, intent(in) :: fmi
    !-- return
    type(GreenGaussGradientType) :: gradient
    ! -- local
    integer(I4B) :: n, nodes

    gradient%dis => dis
    gradient%fmi => fmi

    nodes = dis%nodes
    nodes = dis%nodes
    allocate (gradient%grad_op(dis%nodes))

    ! -- Create boundary Cells
    gradient%boundary_faces = BoundaryFacesType(dis)

    ! -- Compute the gradient operator
    do n = 1, nodes
      gradient%grad_op(n)%data = gradient%create_grad_operator(n)
    end do

  end function constructor

  function create_grad_operator(this, n) result(grad_op)
    ! -- dummy
    class(GreenGaussGradientType) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:, :), allocatable :: grad_op
    ! -- local
    integer(I4B) :: isympos, m, ihc, ipos, local_pos
    real(DP) :: cl1, cl2, height, area, volume
    real(DP), dimension(3) :: normal
    real(DP), dimension(:, :), allocatable :: normal_matrix
    real(DP), dimension(:, :), allocatable :: geometric_coefficients
    integer(I4B) :: number_connections, number_boundaries, number_sides

    number_connections = number_connected_faces(this%dis, n)
    number_boundaries = this%boundary_faces%ia(n + 1) - this%boundary_faces%ia(n)
    number_sides = number_connections + number_boundaries

    allocate (normal_matrix(3, number_sides))
    allocate (geometric_coefficients(number_sides, number_sides + 1))
    allocate (grad_op(3, number_sides + 1))

    geometric_coefficients = 0.0_dp

    height = this%dis%top(n) - this%dis%bot(n)
    volume = this%dis%get_cell_volume(n, this%dis%top(n))

    ! Handle the internal connections
    local_pos = 1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)
      isympos = this%dis%con%jas(ipos)
      ihc = this%dis%con%ihc(isympos)

      call this%dis%connection_normal( &
        n, m, ihc, normal(1), normal(2), normal(3), ipos)
      normal_matrix(:, local_pos) = normal

      area = this%dis%con%hwva(isympos)
      if (ihc > 0) then
        area = area * height
      end if

      cl1 = this%dis%con%cl1(isympos)
      cl2 = this%dis%con%cl2(isympos)
      geometric_coefficients(local_pos, 1) = cl2 / (cl1 + cl2) * area
      geometric_coefficients(local_pos, local_pos + 1) = cl1 / (cl1 + cl2) * area

      local_pos = local_pos + 1
    end do

    ! Handle the boundary cells
    do ipos = this%boundary_faces%ia(n), this%boundary_faces%ia(n + 1) - 1
      normal_matrix(:, local_pos) = this%boundary_faces%get_normal(ipos)

      area = this%boundary_faces%get_area(ipos)
      cl1 = this%boundary_faces%get_cl1(ipos)
      cl2 = cl1

      geometric_coefficients(local_pos, 1) = cl2 / (cl1 + cl2) * area
      geometric_coefficients(local_pos, local_pos + 1) = cl1 / (cl1 + cl2) * area

      local_pos = local_pos + 1
    end do

    grad_op = matmul(normal_matrix, geometric_coefficients) / volume

  end function create_grad_operator

  function get(this, n, c) result(grad_c)
    ! -- dummy
    class(GreenGaussGradientType), target :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: c
    !-- return
    real(DP), dimension(3) :: grad_c

    grad_c = this%compute_cell_gradient(n, c)
  end function get

  function compute_cell_gradient(this, n, cnew) result(grad_c)
    ! -- dummy
    class(GreenGaussGradientType), target :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: cnew
    !-- return
    real(DP), dimension(3) :: grad_c
    ! -- local
    real(DP), dimension(:, :), pointer :: grad_op
    real(DP), dimension(:), allocatable :: c
    integer(I4B) :: m, local_pos, ipos
    integer(I4B) :: number_connections, number_boundaries, number_sides

    number_connections = number_connected_faces(this%dis, n)
    number_boundaries = this%boundary_faces%ia(n + 1) - this%boundary_faces%ia(n)
    number_sides = number_connections + number_boundaries
    allocate (c(number_sides + 1))

    c(1) = cnew(n)
    local_pos = 2
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)
      c(local_pos) = cnew(m)
      local_pos = local_pos + 1
    end do

    do ipos = this%boundary_faces%ia(n), this%boundary_faces%ia(n + 1) - 1
      c(local_pos) = cnew(n)
      local_pos = local_pos + 1
    end do

    grad_op => this%grad_op(n)%data
    grad_c = matmul(grad_op, c)

  end function compute_cell_gradient

end module GreenGaussGradientModule
