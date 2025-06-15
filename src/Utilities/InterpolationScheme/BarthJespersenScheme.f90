module BarthJespersenSchemeModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE, DSAME, DPI, DEM6
  use IInterpolationSchemeModule, only: IInterpolationSchemeType, CoefficientsType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use IGradient, only: IGradientType

  implicit none
  private

  public :: BarthJespersenSchemeType

  type, extends(IInterpolationSchemeType) :: BarthJespersenSchemeType
    private
    class(DisBaseType), pointer :: dis
    type(TspFmiType), pointer :: fmi
    class(IGradientType), pointer :: gradient
    real(DP), dimension(:, :), allocatable :: face_centers_con
  contains
    procedure :: compute

    procedure, private :: get_cell_position
    procedure, private :: face_center
    procedure, private :: face_centers
    procedure, private :: face_limited
    procedure, private :: compute_gradient_limiter
    procedure, private :: create_face_centers
    procedure, private :: find_connecting_face_center

  end type BarthJespersenSchemeType

  interface BarthJespersenSchemeType
    module procedure Constructor
  end interface BarthJespersenSchemeType

contains

  function Constructor(dis, fmi, gradient) result(interpolation_scheme)
    ! -- return
    type(BarthJespersenSchemeType) :: interpolation_scheme
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    type(TspFmiType), pointer, intent(in) :: fmi
    class(IGradientType), allocatable, target, intent(in) :: gradient

    interpolation_scheme%dis => dis
    interpolation_scheme%fmi => fmi
    interpolation_scheme%gradient => gradient

    call interpolation_scheme%create_face_centers()

  end function Constructor

  function compute(this, n, m, iposnm, phi) result(phi_face)
    !-- return
    type(CoefficientsType), target :: phi_face
    ! -- dummy
    class(BarthJespersenSchemeType), target :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: iposnm
    real(DP), intent(in), dimension(:) :: phi
    ! -- local
    integer(I4B) :: iup, idn, isympos
    real(DP) :: qnm
    real(DP), pointer :: coef_up, coef_dn
    real(DP), dimension(3) :: grad_c, xup, xf
    real(DP) :: grad_lim_cup

    isympos = this%dis%con%jas(iposnm)
    qnm = this%fmi%gwfflowja(iposnm)

    !
    ! -- Find upstream node
    if (qnm > DZERO) then
      ! -- positive flow into n means m is upstream
      iup = m
      idn = n

      coef_up => phi_face%c_m
      coef_dn => phi_face%c_n
    else
      iup = n
      idn = m

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

    xup = this%get_cell_position(iup)
    xf = this%face_center(iup, idn)

    grad_c = this%gradient%get(iup, phi)
    grad_lim_cup = this%compute_gradient_limiter(iup, grad_c, phi)

    phi_face%rhs = -grad_lim_cup * dot_product(grad_c, xf - xup)

  end function compute

  function get_cell_position(this, n) result(xn)
    ! -- return
    real(DP), dimension(3) :: xn
    ! -- dummy
    class(BarthJespersenSchemeType) :: this
    integer(I4B), intent(in) :: n
    ! -- local

    xn(1) = this%dis%xc(n)
    xn(2) = this%dis%yc(n)
    xn(3) = (this%dis%top(n) + this%dis%bot(n)) / 2.0_dp

  end function get_cell_position

  function compute_gradient_limiter(this, n, grad_c, cnew) result(psi)
    ! -- return
    real(DP) :: psi
    ! -- dummy
    class(BarthJespersenSchemeType) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: grad_c
    real(DP), intent(in), dimension(:) :: cnew
    ! -- local
    integer(I4B) :: ipos, m
    real(DP) :: min_c, max_c
    real(DP) :: del_plus_max, del_plus_min

    min_c = cnew(n)
    max_c = cnew(n)
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)
      min_c = min(min_c, cnew(m))
      max_c = max(max_c, cnew(m))
    end do

    del_plus_max = max_c - cnew(n)
    del_plus_min = min_c - cnew(n)

    psi = this%face_limited(n, del_plus_max, del_plus_min, grad_c)

  end function compute_gradient_limiter

  function face_limited(this, n, del_plus_max, del_plus_min, grad_c) result(psi)

    ! -- return
    real(DP) :: psi
    ! -- dummy
    class(BarthJespersenSchemeType) :: this
    integer(I4B), intent(in) :: n
    real(DP), intent(in) :: del_plus_max
    real(DP), intent(in) :: del_plus_min
    real(DP), dimension(3), intent(in) :: grad_c
    ! -- local
    integer(I4B) :: ipos, m
    real(DP), dimension(3) :: xn, xf, dnf
    real(DP) :: del_min, psi_m, height
    real(DP) :: K, dx, eps2, del_plus, top

    xn = this%get_cell_position(n)
    psi = 1.0_dp

    K = 0.01_dp

    top = this%dis%top(n)
    height = this%dis%top(n) - this%dis%bot(n)
    dx = 2.0_dp * sqrt(this%dis%get_cell_volume(n, top) / (DPI * height))
    eps2 = (K * dx)**3_dp

    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)

      xf = this%face_center(n, m)
      dnf = xf - xn

      del_min = dot_product(grad_c, dnf)

      if (del_min > 0.0_dp) then
        del_plus = del_plus_max
      elseif (del_min < 0.0_dp) then
        del_plus = del_plus_min
      else
        cycle
      end if

      psi_m = &
        1.0_dp / del_min * &
        ((del_plus**2_dp + eps2) * del_min + 2_dp * del_min**2_dp * del_plus) / &
        (del_plus**2_dp + 2_dp * del_min**2_dp + del_plus * del_min + eps2)

      psi = min(psi, psi_m)
    end do

  end function face_limited

  function face_center(this, n, m) result(xf)
    ! -- return
    real(DP), dimension(3) :: xf
    ! -- dummy
    class(BarthJespersenSchemeType) :: this
    integer(I4B), intent(in) :: n, m
    ! -- local
    integer(I4B) :: ipos, isympos

    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      if (m == this%dis%con%ja(ipos)) exit
    end do

    isympos = this%dis%con%jas(ipos)

    xf = this%face_centers_con(isympos, :)

  end function face_center

  subroutine create_face_centers(this)
    ! -- dummy
    class(BarthJespersenSchemeType), target :: this
    ! -- local
    integer(I4B) :: n, ipos, m, isympos, ihc
    real(DP) :: cl1
    real(DP), dimension(3) :: xn, normal

    allocate (this%face_centers_con(this%dis%con%njas, 3))

    do n = 1, this%dis%nodes
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        isympos = this%dis%con%jas(ipos)
        ihc = this%dis%con%ihc(isympos)

        if (ihc == 0) then
          xn = this%get_cell_position(n)
          cl1 = this%dis%con%cl1(isympos)
          call this%dis%connection_normal( &
            n, m, ihc, normal(1), normal(2), normal(3), ipos)
          this%face_centers_con(isympos, :) = xn + normal * cl1
        else
          this%face_centers_con(isympos, :) = &
            this%find_connecting_face_center(n, m)
        end if

      end do
    end do

  end subroutine create_face_centers

  function find_connecting_face_center(this, n, m) result(xf)
    ! -- return
    real(DP), dimension(3) :: xf
    ! -- dummy
    class(BarthJespersenSchemeType), target :: this
    integer(I4B), intent(in) :: n, m
    ! -- local
    real(DP), dimension(:, :), allocatable :: polyverts_n, polyverts_m
    integer(I4B) :: ipos, ipos2

    polyverts_n = this%face_centers(n)
    polyverts_m = this%face_centers(m)

    do ipos = 1, size(polyverts_n, 1)
      do ipos2 = 1, size(polyverts_m, 1)
        if (is_same(polyverts_n(ipos, :), polyverts_m(ipos2, :))) then
          xf = polyverts_n(ipos, :)
          return
        end if
      end do
    end do
  end function

  function face_centers(this, n) result(xf)
    ! -- return
    real(DP), dimension(:, :), allocatable :: xf
    ! -- dummy
    class(BarthJespersenSchemeType) :: this
    integer(I4B), intent(in) :: n
    ! -- local
    integer(I4B) :: iedge
    real(DP), dimension(3) :: xn, v1, v2
    real(DP), dimension(:, :), allocatable :: polyverts
    integer(I4B) :: number_of_edges

    call this%dis%get_polyverts(n, polyverts, .true.)
    number_of_edges = size(polyverts, 2) - 1

    allocate (xf(number_of_edges, 3))

    xn = this%get_cell_position(n)
    do iedge = 1, number_of_edges
      v1(1:2) = polyverts(:, iedge)
      v2(1:2) = polyverts(:, iedge + 1)
      V1(3) = xn(3)
      V2(3) = xn(3)

      xf(iedge, :) = (v1 + v2) / 2.0_dp
    end do

  end function face_centers

  function is_same(A, B) result(res)
    ! -- return
    logical :: res
    ! -- dummy
    real(DP), dimension(3), intent(in) :: A, B

    res = (abs(A(1) - B(1)) < DEM6) .and. &
          (abs(A(2) - B(2)) < DEM6)

  end function is_same

end module BarthJespersenSchemeModule
