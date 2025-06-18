module TVDSchemeModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, DSAME, DHALF
  use IInterpolationSchemeModule, only: IInterpolationSchemeType, CoefficientsType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use IGradient, only: IGradientType
  use DisInfoModule, only: node_distance
  use TspSourceInfoProviderModule, only: TspSourceInfoProviderType

  implicit none
  private

  public :: TVDSchemeType

  type, extends(IInterpolationSchemeType) :: TVDSchemeType
    private
    class(DisBaseType), pointer :: dis
    type(TspFmiType), pointer :: fmi
    class(IGradientType), pointer :: gradient
    class(TspSourceInfoProviderType), pointer :: source_info_provider => null()
    integer(I4B) :: limiter_id = 2 ! default to van Leer limiter
  contains
    procedure :: compute

    procedure, private :: limiter
  end type TVDSchemeType

  interface TVDSchemeType
    module procedure Constructor
  end interface TVDSchemeType

contains
  function Constructor(dis, fmi, gradient, source_info_provider) &
    result(interpolation_scheme)
    ! -- return
    type(TVDSchemeType) :: interpolation_scheme
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    type(TspFmiType), pointer, intent(in) :: fmi
    class(IGradientType), allocatable, target, intent(in) :: gradient
    class(TspSourceInfoProviderType), pointer, intent(in) :: source_info_provider

    interpolation_scheme%dis => dis
    interpolation_scheme%fmi => fmi
    interpolation_scheme%gradient => gradient
    interpolation_scheme%source_info_provider => source_info_provider

  end function Constructor

  function compute(this, n, m, iposnm, phi) result(phi_face)
    !-- return
    type(CoefficientsType), target :: phi_face
    ! -- dummy
    class(TVDSchemeType), target :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: iposnm
    real(DP), intent(in), dimension(:) :: phi
    ! -- local
    integer(I4B) :: iup, idn, isympos
    real(DP) :: qnm
    real(DP), pointer :: coef_up, coef_dn
    real(DP), dimension(3) :: grad_c, dnm
    real(DP) :: smooth, alimiter
    real(DP) :: cl1, cl2, relative_distance
    real(DP) :: c_virtual, c_source
    real(DP) :: min_phi, max_phi
    integer(I4B) :: ipos, mm

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
    dnm = node_distance(this%dis, this%fmi, iup, idn)
    c_virtual = phi(idn) - 2.0_dp * (dot_product(grad_c, dnm))
    !
    ! If there is a source term, we may need to limit the virtual concentration
    c_source = DZERO
    if (associated(this%source_info_provider)) then
      if (this%source_info_provider%has_source(iup)) then
        c_source = this%source_info_provider%get_source(iup)
      end if
    end if

    ! Enforce local TVD condition.
    ! This is done by limiting the virtual concentration to the range of
    ! the max and min concentration of the neighbouring cells.
    min_phi = phi(iup)
    max_phi = phi(iup)
    do ipos = this%dis%con%ia(iup) + 1, this%dis%con%ia(iup + 1) - 1
      mm = this%dis%con%ja(ipos)
      min_phi = min(min_phi, phi(mm))
      max_phi = max(max_phi, phi(mm))
    end do
    !
    ! Apply the local TVD condition
    if (c_virtual > max(max_phi, c_source)) then
      c_virtual = max(max_phi, c_source)
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
    ! cl1 / (cl1 + cl2) * alimiter * qnm * (phi(idn) - phi(iup))
    ! This is split into two parts:
    relative_distance = cl1 / (cl1 + cl2)
    coef_up = coef_up - relative_distance * alimiter
    coef_dn = coef_dn + relative_distance * alimiter

    ! Alternative way of writing the high order term by adding it to the rhs:
    ! phi_face%rhs = -relative_distance * alimiter * (phi(idn) - phi(iup))

  end function compute

  function limiter(this, r) result(theta)
    ! -- return
    real(DP) :: theta ! limited slope
    ! -- dummy
    class(TVDSchemeType) :: this
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
    CASE DEFAULT
      theta = DZERO
    end select

  end function

end module TVDSchemeModule
