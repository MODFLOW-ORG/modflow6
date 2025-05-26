module ExtendedLeastSquaredGradientModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE

  Use IGradient
  use BaseDisModule, only: DisBaseType
  use MathUtilModule, only: eye
  use SVDModule, only: pinv
  use TspFmiModule, only: TspFmiType

  implicit none
  private

  public :: ExtendedLeastSquaredGradientType

  type Array1D
    integer(I4B), dimension(:), allocatable :: data
  end type Array1D

  type Array2D
    real(DP), dimension(:, :), allocatable :: data
  end type Array2D

  type, extends(IGradientType) :: ExtendedLeastSquaredGradientType
    private
    class(DisBaseType), pointer :: dis
    type(TspFmiType), pointer :: fmi
    type(Array2D), allocatable, dimension(:) :: grad_op
    type(Array1D), allocatable, dimension(:) :: connected_cells
  contains
    procedure :: get

    procedure, private :: compute_cell_gradient
    procedure, private :: find_connected_cells
    procedure, private :: create_grad_operator
    procedure, private :: node_distance
  end type ExtendedLeastSquaredGradientType

  interface ExtendedLeastSquaredGradientType
    module procedure Constructor
  end interface ExtendedLeastSquaredGradientType

contains
  function constructor(dis, fmi) Result(gradient)
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    type(TspFmiType), pointer, intent(in) :: fmi
    !-- return
    type(ExtendedLeastSquaredGradientType) :: gradient
    ! -- local
    integer(I4B) :: n, nodes

    gradient%dis => dis
    gradient%fmi => fmi

    nodes = dis%nodes

    allocate (gradient%connected_cells(dis%nodes))
    do n = 1, nodes
      gradient%connected_cells(n)%data = gradient%find_connected_cells(n)
    end do

    allocate (gradient%grad_op(dis%nodes))
    do n = 1, nodes
      gradient%grad_op(n)%data = gradient%create_grad_operator(n)
    end do

  end function constructor

  function find_connected_cells(this, n) result(res)
    use PtrHashTableModule, only: PtrHashTableType
    use IteratorModule, only: IteratorType
    ! -- dummy
    class(ExtendedLeastSquaredGradientType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), allocatable :: res(:)
    ! -- local
    integer(I4B) :: ipos, ipos2, m, mm, local_pos
    character(1) :: c
    class(*), pointer :: m_ptr, mm_ptr
    integer(I4B), pointer :: i_ptr
    type(PtrHashTableType) :: connected_cells
    class(IteratorType), allocatable :: itr

    ! Find 1st and 2nd degree connected cells (Direct enighbours and neighbours of neighbours)
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)

      m_ptr => this%dis%con%ja(ipos)
      c = char(this%dis%con%ja(ipos))
      if (.not. connected_cells%contains(c)) then
        call connected_cells%add(c, m_ptr)
      end if

      do ipos2 = this%dis%con%ia(m) + 1, this%dis%con%ia(m + 1) - 1
        mm = this%dis%con%ja(ipos2)

        mm_ptr => this%dis%con%ja(ipos2)
        c = char(this%dis%con%ja(ipos2))
        if (.not. connected_cells%contains(c) .and. mm /= n) then
          call connected_cells%add(c, mm_ptr)
        end if
      end do
    end do

    ! Store the results in an array
    allocate (res(connected_cells%count()))
    local_pos = 1
    allocate (itr, source=connected_cells%iterator())
    do while (itr%has_next())
      call itr%next()

      select type (val => itr%value())
      type is (integer(I4B))
        i_ptr => val
      end select
      res(local_pos) = i_ptr
      local_pos = local_pos + 1
    end do

  end function find_connected_cells

  function create_grad_operator(this, n) result(grad_op)
    use MathUtilModule, only: eye
    ! -- dummy
    class(ExtendedLeastSquaredGradientType) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:, :), allocatable :: grad_op
    ! -- local
    integer(I4B) :: ipos, m
    real(DP) :: length
    real(DP), dimension(3) :: dnm
    real(DP), dimension(:, :), allocatable :: d
    real(DP), dimension(:, :), allocatable :: d_trans
    real(DP), dimension(:, :), allocatable :: grad_scale
    real(DP), dimension(:, :), allocatable :: W
    real(DP), dimension(3, 3) :: g
    real(DP), dimension(3, 3) :: g_inv
    integer(I4B), allocatable :: connected_cells(:)
    integer(I4B) :: num_connected_cells

    connected_cells = this%connected_cells(n)%data
    num_connected_cells = size(connected_cells)

    allocate (d(num_connected_cells, 3))
    allocate (d_trans(3, num_connected_cells))
    allocate (grad_op(3, num_connected_cells))
    allocate (grad_scale(num_connected_cells, num_connected_cells))
    allocate (W(num_connected_cells, num_connected_cells))

    grad_scale = 0
    d = 0
    W = eye(num_connected_cells)

    ! Assemble the distance matrix
    do ipos = 1, num_connected_cells
      m = connected_cells(ipos)

      dnm = this%node_distance(n, m)
      length = norm2(dnm)

      d(ipos, :) = dnm / length
      grad_scale(ipos, ipos) = 1.0_dp / length
    end do

    d_trans = transpose(d)

    ! Compute the G and inverse G matrices
    g = matmul(d_trans, matmul(W, d))
    g_inv = pinv(g)

    ! Compute the gradient operator
    grad_op = matmul(matmul(matmul(g_inv, d_trans), W), grad_scale)

  end function create_grad_operator

  function node_distance(this, n, m) result(d)
    ! -- return
    real(DP), dimension(3) :: d
    ! -- dummy
    class(ExtendedLeastSquaredGradientType) :: this
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
    class(ExtendedLeastSquaredGradientType), target :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: c
    !-- return
    real(DP), dimension(3) :: grad_c

    grad_c = this%compute_cell_gradient(n, c)
  end function get

  function compute_cell_gradient(this, n, cnew) result(grad_c)
    ! -- dummy
    class(ExtendedLeastSquaredGradientType), target :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: cnew
    !-- return
    real(DP), dimension(3) :: grad_c
    ! -- local
    real(DP), dimension(:, :), pointer :: grad_op
    integer(I4B) :: ipos

    integer(I4B) :: m
    real(DP), dimension(:), allocatable :: dc

    integer(I4B), allocatable :: connected_cells(:)
    integer(I4B) :: num_connected_cells

    connected_cells = this%connected_cells(n)%data
    num_connected_cells = size(connected_cells)

    ! Assemble the concentration difference vector
    allocate (dc(num_connected_cells))
    do ipos = 1, num_connected_cells
      m = connected_cells(ipos)
      dc(ipos) = cnew(m) - cnew(n)
    end do

    ! Compute the cells gradient
    grad_op => this%grad_op(n)%data
    grad_c = matmul(grad_op, dc)

  end function compute_cell_gradient

end module
