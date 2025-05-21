module TVDSchemeModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, DSAME, DHALF
  use IInterpolationSchemeModule, only: IInterpolationSchemeType, CoefficientsType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use IGradient, only: IGradientType

  implicit none
  private

  public :: TVDSchemeType

  type, extends(IInterpolationSchemeType) :: TVDSchemeType
    private
    class(DisBaseType), pointer :: dis
    type(TspFmiType), pointer :: fmi
    class(IGradientType), pointer :: gradient
    integer(I4B) :: limiter_id = 2 ! default to van Leer limiter
  contains
    procedure :: compute

    procedure, private :: node_distance
    procedure, private :: limiter
  end type TVDSchemeType

  interface TVDSchemeType
    module procedure Constructor
  end interface TVDSchemeType

contains
  function Constructor(dis, fmi, gradient) result(interpolation_scheme)
    ! -- return
    type(TVDSchemeType) :: interpolation_scheme
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    type(TspFmiType), pointer, intent(in) :: fmi
    class(IGradientType), allocatable, target, intent(in) :: gradient

    interpolation_scheme%dis => dis
    interpolation_scheme%fmi => fmi
    interpolation_scheme%gradient => gradient

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
    integer(I4B) :: iup, idn, isympos, ihc
    real(DP) :: qnm
    real(DP), pointer :: coef_up, coef_dn
    real(DP), dimension(3) :: grad_c, dnm, normal
    real(DP) :: smooth, alimiter
    real(DP) :: cl1, cl2, c_virtual, relative_distance

    isympos = this%dis%con%jas(iposnm)
    qnm = this%fmi%gwfflowja(iposnm)

    ihc = this%dis%con%ihc(isympos)
    call this%dis%connection_normal( &
      n, m, ihc, normal(1), normal(2), normal(3), iposnm)

    !
    ! -- Find upstream node
    if (qnm > DZERO) then
      ! -- positive flow into n means m is upstream
      iup = m
      idn = n

      cl1 = this%dis%con%cl2(isympos)
      cl2 = this%dis%con%cl1(isympos)

      normal = -normal

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
    ! -- Compute smoothness factor
    dnm = this%node_distance(iup, idn)
    ! dnm = normal * (cl1 + cl2)
    smooth = 2.0_dp * (dot_product(grad_c, dnm)) / &
             (phi(idn) - phi(iup)) - 1.0_dp
    !
    ! -- Correct smoothness factor to prevent negative concentration
    ! c_virtual = phi(iup) - smooth * (phi(idn) - phi(iup))
    ! if (c_virtual < DZERO) then
    !   smooth = phi(iup) / (phi(idn) - phi(iup))
    ! end if
    !
    ! -- Compute limiter
    alimiter = this%limiter(smooth)

    ! High order term is:
    ! DHALF * alimiter * qnm * (phi(idn) - phi(iup))
    ! This is split into two parts:
    ! relative_distance = cl1 / (cl1 + cl2)
    relative_distance = DHALF
    ! coef_up = coef_up - relative_distance * alimiter
    ! coef_dn = coef_dn + relative_distance * alimiter

    phi_face%rhs = -relative_distance * alimiter * (phi(idn) - phi(iup))

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

  function node_distance(this, n, m) result(d)
    ! -- return
    real(DP), dimension(3) :: d
    ! -- dummy
    class(TVDSchemeType) :: this
    integer(I4B), intent(in) :: n, m
    ! -- local
    real(DP) :: x_dir, y_dir, z_dir, length
    real(DP) :: satn, satm
    integer(I4B) :: ipos, isympos, ihc

    isympos = -1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      if (this%dis%con%ja(ipos) == m) then
        isympos = this%dis%con%jas(ipos)
        exit
      end if
    end do

    if (isympos == -1) then
      ! handle exception
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

end module TVDSchemeModule
