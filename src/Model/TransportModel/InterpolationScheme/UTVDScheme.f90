module UTVDSchemeModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, DSAME, DHALF
  use InterpolationSchemeInterfaceModule, only: InterpolationSchemeInterface, &
                                                CoefficientsType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use IGradient, only: IGradientType

  use DisUtilsModule, only: node_distance
  use LocalCellExtremaModule, only: LocalCellExtremaType

  implicit none
  private

  public :: UTVDSchemeType

  !> @brief Total Variation Diminishing (TVD) interpolation scheme.
  !!
  !! This class implements a high-resolution, TVD interpolation scheme for use in transport modeling.
  !! It extends a generic interpolation scheme interface and supports multiple TVD limiters (van Leer,
  !! Koren, Superbee, van Albada, etc.) for controlling numerical diffusion and oscillations.
  !! The default limiter is van Leer, but others can be selected by changing the `limiter_id` member.
  !!
  !! The scheme uses a combination of low-order upwind and high-order limited terms to compute
  !! face concentrations. The high-order term is constructed using a gradient-based virtual node
  !! value, following the approach described by Darwish et al. An additional TVD clamp is applied
  !! to the virtual node value to enforce TVD compliance, especially on grids where the original
  !! method may not guarantee monotonicity.
  !!
  !! - Supports both structured and unstructured grids via polymorphic discretization and gradient objects.
  !! - The limiter can be selected via the `limiter_id` member (default is van Leer).
  !! - The method `find_local_extrema` finds the minimum and maximum values among the current cell and its neighbors,
  !!   which is used to enforce the TVD condition.
  !! - The `compute` method calculates the face coefficients for the transport equation.
  !<
  type, extends(InterpolationSchemeInterface) :: UTVDSchemeType
    private
    class(DisBaseType), pointer :: dis
    type(TspFmiType), pointer :: fmi
    class(IGradientType), pointer :: gradient

    real(DP), dimension(:), pointer :: phi
    type(LocalCellExtremaType), allocatable :: min_max_phi ! local minimum values at nodes
    integer(I4B) :: limiter_id = 2 ! default to van Leer limiter
    real(DP), dimension(:, :), allocatable :: cached_node_distance ! distance vectors
  contains
    procedure :: compute
    procedure :: set_field
    final :: destructor

    procedure, private :: limiter
    procedure, private :: compute_node_distance
  end type UTVDSchemeType

  interface UTVDSchemeType
    module procedure constructor
  end interface UTVDSchemeType

contains

  function constructor(dis, fmi, gradient) &
    result(interpolation_scheme)
    ! -- return
    type(UTVDSchemeType) :: interpolation_scheme
    ! -- dummy
    class(DisBaseType), pointer, intent(in) :: dis
    type(TspFmiType), pointer, intent(in) :: fmi
    class(IGradientType), allocatable, intent(in), target :: gradient

    interpolation_scheme%dis => dis
    interpolation_scheme%fmi => fmi

    interpolation_scheme%gradient => gradient
    interpolation_scheme%min_max_phi = LocalCellExtremaType(dis)

    allocate (interpolation_scheme%cached_node_distance(dis%njas, 3))
    call compute_node_distance(interpolation_scheme)

  end function constructor

  subroutine destructor(this)
    ! -- dummy
    type(UTVDSchemeType), intent(inout) :: this

    deallocate (this%cached_node_distance)
  end subroutine destructor

  !> @brief Set the scalar field for which interpolation will be computed
  !!
  !! This method establishes a pointer to the scalar field data and updates
  !! any dependent cached data (gradients and local extrema) to ensure
  !! subsequent interpolation computations use the current field values.
  !!
  !<
  subroutine set_field(this, phi)
    ! -- dummy
    class(UTVDSchemeType), target :: this
    real(DP), intent(in), dimension(:), pointer :: phi

    this%phi => phi
    call this%gradient%set_field(phi)
    call this%min_max_phi%set_field(phi)
  end subroutine set_field

  function compute(this, n, m, iposnm) result(phi_face)
    !-- return
    type(CoefficientsType), target :: phi_face ! Output: coefficients for the face between cells n and m
    ! -- dummy
    class(UTVDSchemeType), target :: this
    integer(I4B), intent(in) :: n ! Index of the first cell
    integer(I4B), intent(in) :: m ! Index of the second cell
    integer(I4B), intent(in) :: iposnm ! Position in the connectivity array for the n-m connection
    ! -- local
    integer(I4B) :: iup, idn, isympos ! Indices for upwind, downwind, and symmetric position in connectivity
    real(DP) :: qnm ! Flow rate across the face between n and m
    real(DP), pointer :: coef_up, coef_dn ! Pointers to upwind and downwind coefficients in phi_face
    real(DP), dimension(3) :: grad_c ! gradient at upwind cell
    real(DP), dimension(3) :: dnm ! vector from upwind to downwind cell
    real(DP) :: smooth ! Smoothness indicator for limiter i.e. ratio of gradients
    real(DP) :: alimiter ! Value of the TVD limiter
    real(DP) :: cl1, cl2 ! Connection lengths from upwind and downwind cells to the face
    real(DP) :: relative_distance ! Relative distance factor for high-order term
    real(DP) :: c_virtual ! Virtual node concentration (Darwish method)
    real(DP), pointer :: min_phi, max_phi ! Local minimum and maximum among cell and neighbors

    isympos = this%dis%con%jas(iposnm)
    qnm = this%fmi%gwfflowja(iposnm)
    !
    ! -- Find upstream node
    if (qnm > DZERO) then
      ! -- positive flow into n means m is upstream
      iup = m
      idn = n

      cl1 = this%dis%con%cl2(isympos)
      cl2 = this%dis%con%cl1(isympos)

      coef_up => phi_face%c_m
      coef_dn => phi_face%c_n
    else
      iup = n
      idn = m

      cl1 = this%dis%con%cl1(isympos)
      cl2 = this%dis%con%cl2(isympos)

      coef_up => phi_face%c_n
      coef_dn => phi_face%c_m
    end if
    !
    ! Determine direction of distance vector from upwind to downwind cell
    ! The cached_node_distance always stores vector from lower-numbered node to higher-numbered node.
    ! Since we need dnm to point from upwind (iup) to downwind (idn), we must adjust the sign:
    ! - If iup > idn: the cached vector points from idn to iup, so we negate it to get iup to idn
    ! - If iup < idn: the cached vector already points from iup to idn, so use it as-is
    if (iup > idn) then
      dnm = -this%cached_node_distance(isympos, :)
    else
      dnm = this%cached_node_distance(isympos, :)
    end if
    !
    ! -- Add low order terms
    coef_up = DONE
    !
    ! -- Add high order terms
    !
    ! -- Return if straddled cells have same value
    if (abs(this%phi(idn) - this%phi(iup)) < DSAME) return
    !
    ! -- Compute cell concentration gradient
    grad_c = this%gradient%get(iup)
    !
    ! Darwish's method to compute virtual node concentration
    c_virtual = this%phi(idn) - 2.0_dp * (dot_product(grad_c, dnm))
    !
    ! Enforce local TVD condition.
    ! This is done by limiting the virtual concentration to the range of
    ! the max and min concentration of the neighbouring cells.
    min_phi => this%min_max_phi%get_min(iup)
    max_phi => this%min_max_phi%get_max(iup)

    if (c_virtual > max_phi) then
      c_virtual = max_phi
    end if

    if (c_virtual < max(min_phi, DZERO)) then
      c_virtual = max(min_phi, DZERO)
    end if
    !
    ! -- Compute smoothness factor
    smooth = (this%phi(iup) - c_virtual) / (this%phi(idn) - this%phi(iup))
    !
    ! -- Compute limiter
    alimiter = this%limiter(smooth)

    ! High order term is:
    relative_distance = cl1 / (cl1 + cl2)
    phi_face%rhs = -relative_distance * alimiter * (this%phi(idn) - this%phi(iup))

    ! Alternative way of writing the high order term by adding it to the
    ! coefficients matrix. The equation to be added is:
    ! high_order = cl1 / (cl1 + cl2) * alimiter * qnm * (phi(idn) - phi(iup))
    ! This is split into two parts:
    ! coef_up = coef_up - relative_distance * alimiter
    ! coef_dn = coef_dn + relative_distance * alimiter

  end function compute

  function limiter(this, r) result(theta)
    ! -- return
    real(DP) :: theta ! limited slope
    ! -- dummy
    class(UTVDSchemeType) :: this
    real(DP) :: r ! ratio of successive gradients

    select case (this%limiter_id)
    case (2) ! van Leer
      theta = max(0.0_dp, min((r + dabs(r)) / (1.0_dp + dabs(r)), 2.0_dp))
    case (3) ! Koren
      theta = max(0.0_dp, min(2.0_dp * r, &
                              1.0_dp / 3.0_dp + 2.0_dp / 3.0_dp * r, 2.0_dp))
    case (4) ! Superbee
      theta = max(0.0_dp, min(2.0_dp * r, 1.0_dp), min(r, 2.0_dp))
    case (5) ! van Albada
      theta = max(0.0_dp, (r * r + r) / (r * r + 1.0_dp))
    case (6) ! Koren modified
      theta = max(0.0_dp, min(4.0_dp * r * r + r, &
                              1.0_dp / 3.0_dp + 2.0_dp / 3.0_dp * r, 2.0_dp))
    case default
      theta = DZERO
    end select
  end function

  subroutine compute_node_distance(this)
    ! -- dummy
    class(UTVDSchemeType), target :: this
    ! -- local
    integer(I4B) :: n, m, ipos, isympos

    this%cached_node_distance = 0.0_dp
    do n = 1, this%dis%nodes
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        if (m <= n) cycle

        isympos = this%dis%con%jas(ipos)
        this%cached_node_distance(isympos, :) = node_distance(this%dis, n, m)
      end do
    end do

  end subroutine compute_node_distance

end module UTVDSchemeModule
