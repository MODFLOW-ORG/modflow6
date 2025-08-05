module UTVDSchemeModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, DSAME, DHALF
  use InterpolationSchemeInterfaceModule, only: InterpolationSchemeInterface, &
                                                CoefficientsType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use IGradient, only: IGradientType
  use DisUtilsModule, only: node_distance

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
    integer(I4B) :: limiter_id = 2 ! default to van Leer limiter
  contains
    procedure :: compute

    procedure, private :: find_local_extrema
    procedure, private :: limiter
  end type UTVDSchemeType

  interface UTVDSchemeType
    module procedure constructor
  end interface UTVDSchemeType

contains
  function constructor(dis, fmi, gradient) &
    result(interpolation_scheme)
    ! -- return
    type(UTVDSchemeType) :: interpolation_scheme
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    type(TspFmiType), pointer, intent(in) :: fmi
    class(IGradientType), allocatable, target, intent(in) :: gradient

    interpolation_scheme%dis => dis
    interpolation_scheme%fmi => fmi
    interpolation_scheme%gradient => gradient

  end function constructor

  function compute(this, n, m, iposnm, phi) result(phi_face)
    !-- return
    type(CoefficientsType), target :: phi_face ! Output: coefficients for the face between cells n and m
    ! -- dummy
    class(UTVDSchemeType), target :: this
    integer(I4B), intent(in) :: n ! Index of the first cell
    integer(I4B), intent(in) :: m ! Index of the second cell
    integer(I4B), intent(in) :: iposnm ! Position in the connectivity array for the n-m connection
    real(DP), intent(in), dimension(:) :: phi ! Array of scalar values (e.g., concentrations) at all cells
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
    real(DP) :: min_phi, max_phi ! Local minimum and maximum among cell and neighbors

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
    ! -- Add low order terms
    coef_up = DONE
    !
    ! -- Add high order terms
    !
    ! -- Return if straddled cells have same value
    if (abs(phi(idn) - phi(iup)) < DSAME) return
    !
    ! -- Compute cell concentration gradient
    grad_c = this%gradient%get(iup, phi)
    !
    ! Darwish's method to compute virtual node concentration
    dnm = node_distance(this%dis, iup, idn)
    c_virtual = phi(idn) - 2.0_dp * (dot_product(grad_c, dnm))
    !
    ! Enforce local TVD condition.
    ! This is done by limiting the virtual concentration to the range of
    ! the max and min concentration of the neighbouring cells.
    call find_local_extrema(this, iup, phi, min_phi, max_phi)

    if (c_virtual > max_phi) then
      c_virtual = max_phi
    end if

    if (c_virtual < max(min_phi, DZERO)) then
      c_virtual = max(min_phi, DZERO)
    end if
    !
    ! -- Compute smoothness factor
    smooth = (phi(iup) - c_virtual) / (phi(idn) - phi(iup))
    !
    ! -- Compute limiter
    alimiter = this%limiter(smooth)

    ! High order term is:
    relative_distance = cl1 / (cl1 + cl2)
    phi_face%rhs = -relative_distance * alimiter * (phi(idn) - phi(iup))

    ! Alternative way of writing the high order term by adding it to the
    ! coefficients matrix. The equation to be added is:
    ! high_order = cl1 / (cl1 + cl2) * alimiter * qnm * (phi(idn) - phi(iup))
    ! This is split into two parts:
    ! coef_up = coef_up - relative_distance * alimiter
    ! coef_dn = coef_dn + relative_distance * alimiter

  end function compute

  subroutine find_local_extrema(this, n, phi, min_phi, max_phi)
    ! -- dummy
    class(UTVDSchemeType) :: this
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

end module UTVDSchemeModule
