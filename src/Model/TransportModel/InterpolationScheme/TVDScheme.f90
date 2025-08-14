module TVDSchemeModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, DPREC, DHALF, DTWO
  use InterpolationSchemeInterfaceModule, only: InterpolationSchemeInterface, &
                                                CoefficientsType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType

  implicit none
  private

  public :: TVDSchemeType

  type, extends(InterpolationSchemeInterface) :: TVDSchemeType
    private
    class(DisBaseType), pointer :: dis
    type(TspFmiType), pointer :: fmi
    real(DP), dimension(:), pointer :: phi
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
  contains
    procedure :: compute
    procedure :: set_field
  end type TVDSchemeType

  interface TVDSchemeType
    module procedure constructor
  end interface TVDSchemeType

contains
  function constructor(dis, fmi, ibound) result(interpolation_scheme)
    ! -- return
    type(TVDSchemeType) :: interpolation_scheme
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    type(TspFmiType), pointer, intent(in) :: fmi
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: ibound

    interpolation_scheme%dis => dis
    interpolation_scheme%fmi => fmi
    interpolation_scheme%ibound => ibound

  end function constructor

  !> @brief Set the scalar field for which interpolation will be computed
  !!
  !! This method establishes a pointer to the scalar field data for
  !! subsequent TVD interpolation computations.
  !<
  subroutine set_field(this, phi)
    ! -- dummy
    class(TVDSchemeType), target :: this
    real(DP), intent(in), dimension(:), pointer :: phi

    this%phi => phi
  end subroutine set_field

  function compute(this, n, m, iposnm) result(phi_face)
    !-- return
    type(CoefficientsType), target :: phi_face
    ! -- dummy
    class(TVDSchemeType), target :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: iposnm
    ! -- local
    integer(I4B) :: ipos, isympos, iup, idn, i2up, j
    real(DP) :: qnm, qmax, qupj, elupdn, elup2up
    real(DP) :: smooth, cdiff, alimiter
    real(DP), pointer :: coef_up, coef_dn
    !
    ! -- Find upstream node
    isympos = this%dis%con%jas(iposnm)
    qnm = this%fmi%gwfflowja(iposnm)
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
    elupdn = this%dis%con%cl1(isympos) + this%dis%con%cl2(isympos)
    !
    ! -- Add low order terms
    coef_up = DONE
    !
    ! -- Add high order terms
    !
    ! -- Find second node upstream to iup
    i2up = 0
    qmax = DZERO
    do ipos = this%dis%con%ia(iup) + 1, this%dis%con%ia(iup + 1) - 1
      j = this%dis%con%ja(ipos)
      if (this%ibound(j) == 0) cycle
      qupj = this%fmi%gwfflowja(ipos)
      isympos = this%dis%con%jas(ipos)
      if (qupj > qmax) then
        qmax = qupj
        i2up = j
        elup2up = this%dis%con%cl1(isympos) + this%dis%con%cl2(isympos)
      end if
    end do
    !
    ! -- Calculate flux limiting term
    if (i2up > 0) then
      smooth = DZERO
      cdiff = ABS(this%phi(idn) - this%phi(iup))
      if (cdiff > DPREC) then
        smooth = (this%phi(iup) - this%phi(i2up)) / elup2up * &
                 elupdn / (this%phi(idn) - this%phi(iup))
      end if
      if (smooth > DZERO) then
        alimiter = DTWO * smooth / (DONE + smooth)
        phi_face%rhs = -DHALF * alimiter * (this%phi(idn) - this%phi(iup))
      end if
    end if

  end function compute

end module TVDSchemeModule
