module UzrStorageModule
  use KindModule, only: I4B, LGP, DP
  use ConstantsModule, only: DONE, DTWO, DHALF, DZERO, DPREC, LENVARNAME
  use MatrixBaseModule, only: MatrixBaseType
  use BaseDisModule, only: DisBaseType
  use GwfStoModule, only: GwfStoType
  use GwfStoExtModule, only: GwfStoExtType
  use UzrSoilModelModule, only: SoilModelType
  implicit none
  private

  public :: UzrStorageType

  type, extends(GwfStoExtType) :: UzrStorageType
    class(SoilModelType), pointer :: soil_model => null()

    integer(I4B), pointer :: storage_scheme => null() !< points to UZR data
    integer(I4B), pointer, dimension(:), contiguous :: uzr_iunsat => null() !< points to UZR data
    class(DisBaseType), pointer :: gwf_dis => null() !< points to UZR data
    type(GwfStoType), pointer :: gwf_sto => null() !< points to UZR data
  contains
    procedure :: initialize
    procedure :: is_active => uft_is_active
    procedure :: fc => uft_fc
    procedure :: fn => uft_fn
    procedure :: cq => uft_cq
    procedure :: destroy
    ! private
    procedure, private :: calculate_coeffs
    procedure, private :: calculate_coeffs_nwt
  end type UzrStorageType

contains

  subroutine initialize(this, iunsat, scheme, soil_model, dis, sto)
    class(UzrStorageType), intent(inout) :: this
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: iunsat
    integer(I4B), pointer, intent(in) :: scheme
    class(SoilModelType), pointer, intent(in) :: soil_model
    class(DisBaseType), pointer, intent(in) :: dis
    type(GwfStoType), pointer, intent(in) :: sto

    this%uzr_iunsat => iunsat

    this%storage_scheme => scheme
    this%soil_model => soil_model

    this%gwf_dis => dis
    this%gwf_sto => sto

  end subroutine initialize

  function uft_is_active(this, n) result(is_active)
    class(UzrStorageType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    logical(LGP) :: is_active

    is_active = .false.
    if (this%uzr_iunsat(n) == 1) then
      is_active = .true.
    end if

  end function uft_is_active

  subroutine uft_fc(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
    use GwfStorageUtilsModule, only: SsCapacity
    class(UzrStorageType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: h_old
    real(DP), dimension(:), intent(in) :: h_new
    ! local
    real(DP), dimension(4) :: coeffs !< the coefficients to add to the linear system
    integer(I4B) :: idiag !< the position of the diagonal element

    if (this%gwf_sto%inewton == 0) then
      call this%calculate_coeffs(n, h_old, h_new, coeffs)
    else
      call this%calculate_coeffs_nwt(n, h_old, h_new, coeffs)
    end if

    idiag = this%gwf_dis%con%ia(n)
    call matrix_sln%add_value_pos(idxglo(idiag), coeffs(1))
    rhs(n) = rhs(n) + coeffs(2)
    call matrix_sln%add_value_pos(idxglo(idiag), coeffs(3))
    rhs(n) = rhs(n) + coeffs(4)

  end subroutine uft_fc

  !> @brief Calculate the matrix and RHS coefficients for formulation and cq
  !<
  subroutine calculate_coeffs(this, n, h_old, h_new, coeffs)
    use GwfStorageUtilsModule, only: SsCapacity
    class(UzrStorageType), intent(inout) :: this !< this instance
    integer(I4B), intent(in) :: n !< the node number
    real(DP), dimension(:), intent(in) :: h_old !< the old head
    real(DP), dimension(:), intent(in) :: h_new !< the new head
    real(DP), dimension(4), intent(inout) :: coeffs !< the coefficients: aterm_1, rhs_1, aterm_2, rhs_2
    ! local
    real(DP) :: sc1 !< specific storage capacity
    real(DP) :: top !< the top of the cell
    real(DP) :: bot !< the bottom of the cell
    real(DP) :: area !< area of the cell
    real(DP) :: thk !< cell thickness
    real(DP) :: z !< the nodal elevation for n
    real(DP) :: s1 !< the current saturation
    real(DP) :: s0 !< the saturation at the previous time level
    real(DP) :: Cm !< the current moisture capacity
    real(DP) :: dsdh_lim !< slope
    real(DP) :: h1 !< the current hydraulic head
    real(DP) :: h0 !< the hydraulic head at the previous time level
    real(DP) :: psi1 !< the current saturation
    real(DP) :: psi0 !< the saturation at the previous time level
    real(DP) :: phi !< the cell porosity

    coeffs(:) = DZERO

    top = this%gwf_dis%top(n)
    bot = this%gwf_dis%bot(n)
    thk = top - bot
    area = this%gwf_dis%area(n)

    h1 = h_new(n)
    h0 = h_old(n)
    z = DHALF * (top + bot)
    psi1 = h1 - z
    psi0 = h0 - z

    phi = this%soil_model%porosity(n)
    s1 = this%soil_model%saturation(psi1, n)
    s0 = this%soil_model%saturation(psi0, n)
    Cm = this%soil_model%capacity(psi1, n)
    dsdh_lim = Cm / phi

    sc1 = SsCapacity(this%gwf_sto%istor_coef, top, bot, area, this%gwf_sto%ss(n))

    ! specific storage
    call get_specific_storage_terms(s1, h0, sc1, coeffs(1), coeffs(2))

    select case (this%storage_scheme)
    case (0, 1)
      ! chord slope
      call get_unsat_storage_terms_CS(s1, s0, h1, h0, dsdh_lim, &
                                      phi, z, area, thk, &
                                      coeffs(3), coeffs(4))
    case (2)
      ! modified picard
      call get_unsat_storage_terms_MP(s1, s0, h1, h0, Cm, &
                                      phi, z, area, thk, &
                                      coeffs(3), coeffs(4))
    end select

  end subroutine calculate_coeffs

  !> @brief Calculate the matrix and RHS coefficients for newton formulation
  !<
  subroutine calculate_coeffs_nwt(this, n, h_old, h_new, coeffs)
    use TdisModule, only: delt
    use GwfStorageUtilsModule, only: SsCapacity
    class(UzrStorageType), intent(inout) :: this !< this instance
    integer(I4B), intent(in) :: n !< the node number
    real(DP), dimension(:), intent(in) :: h_old !< the old head
    real(DP), dimension(:), intent(in) :: h_new !< the new head
    real(DP), dimension(4), intent(inout) :: coeffs !< the coefficients: aterm_1, rhs_1, aterm_2, rhs_2
    ! local
    real(DP) :: sc1 !< specific storage capacity
    real(DP) :: top !< the top of the cell
    real(DP) :: bot !< the bottom of the cell
    real(DP) :: area !< area of the cell
    real(DP) :: thk !< cell thickness
    real(DP) :: z !< the nodal elevation for n
    real(DP) :: s1 !< the current saturation
    real(DP) :: s0 !< the saturation at the previous time level
    real(DP) :: Cm !< the current moisture capacity
    real(DP) :: dsdh_lim !< slope
    real(DP) :: h1 !< the current head
    real(DP) :: h0 !< the head at the previous time level
    real(DP) :: psi1 !< the current saturation
    real(DP) :: psi0 !< the saturation at the previous time level
    real(DP) :: phi !< the cell porosity
    real(DP) :: Qss, Qus, dQss_dh, dQus_dh !< the storage rate and derivative

    coeffs(:) = DZERO

    top = this%gwf_dis%top(n)
    bot = this%gwf_dis%bot(n)
    thk = top - bot
    area = this%gwf_dis%area(n)

    h1 = h_new(n)
    h0 = h_old(n)
    z = DHALF * (top + bot)
    psi1 = h1 - z
    psi0 = h0 - z

    phi = this%soil_model%porosity(n)
    s1 = this%soil_model%saturation(psi1, n)
    s0 = this%soil_model%saturation(psi0, n)
    Cm = this%soil_model%capacity(psi1, n)
    dsdh_lim = Cm / phi

    sc1 = SsCapacity(this%gwf_sto%istor_coef, top, bot, area, this%gwf_sto%ss(n))

    Qss = -sc1 * s1 * (h1 - h0) / delt
    Qus = -phi * area * thk * (s1 - s0) / delt

    dQss_dh = -sc1 * s1 / delt - (sc1 * (h1 - h0) / delt) * dsdh_lim
    dQus_dh = -phi * area * thk * dsdh_lim / delt

    coeffs(1) = dQss_dh
    coeffs(2) = -Qss + dQss_dh * h1
    coeffs(3) = dQus_dh
    coeffs(4) = -Qus + dQus_dh * h1

  end subroutine calculate_coeffs_nwt

  subroutine uft_fn(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
    use GwfStorageUtilsModule, only: SsCapacity
    class(UzrStorageType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: h_old
    real(DP), dimension(:), intent(in) :: h_new

    ! TODO_UZR

  end subroutine uft_fn

  subroutine uft_cq(this, n, flowja, h_new, h_old)
    class(UzrStorageType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(inout) :: flowja
    real(DP), dimension(:), intent(in) :: h_new
    real(DP), dimension(:), intent(in) :: h_old
    ! local
    real(DP), dimension(4) :: coeffs !< the linear system coefficients to calculate the flux from
    integer(I4B) :: idiag !< the position of the diagonal element
    real(DP) :: flow_ss !< the specific storage rate for node n
    real(DP) :: flow_sy !< the (unsaturated) specific yield rate for node n

    if (this%gwf_sto%inewton == 0) then
      call this%calculate_coeffs(n, h_old, h_new, coeffs)
    else
      call this%calculate_coeffs_nwt(n, h_old, h_new, coeffs)
    end if

    ! q_n = A_nn * h_n - rhs_n
    flow_ss = coeffs(1) * h_new(n) - coeffs(2)
    flow_sy = coeffs(3) * h_new(n) - coeffs(4)

    ! update flowja
    idiag = this%gwf_dis%con%ia(n)
    flowja(idiag) = flowja(idiag) + flow_ss + flow_sy

    ! store rate
    this%gwf_sto%strgss(n) = flow_ss
    this%gwf_sto%strgsy(n) = flow_sy

  end subroutine uft_cq

  subroutine destroy(this)
    class(UzrStorageType) :: this
  end subroutine destroy

  subroutine get_specific_storage_terms(s_new, h_old, sc1, aterm, rhsterm)
    use TdisModule, only: delt
    real(DP), intent(in) :: s_new
    real(DP), intent(in) :: h_old
    real(DP), intent(in) :: sc1
    real(DP), intent(inout) :: aterm
    real(DP), intent(inout) :: rhsterm

    ! specific storage
    aterm = -sc1 * s_new / delt
    rhsterm = aterm * h_old

  end subroutine get_specific_storage_terms

  subroutine get_unsat_storage_terms_CS(s_new, s_old, h_new, h_old, dsdh_lim, &
                                        phi, z, area, thk, aterm, rhsterm)
    use TdisModule, only: delt
    real(DP), intent(in) :: s_new
    real(DP), intent(in) :: s_old
    real(DP), intent(in) :: h_new
    real(DP), intent(in) :: h_old
    real(DP), intent(in) :: dsdh_lim
    real(DP), intent(in) :: phi
    real(DP), intent(in) :: z
    real(DP), intent(in) :: area
    real(DP), intent(in) :: thk
    real(DP), intent(inout) :: aterm
    real(DP), intent(inout) :: rhsterm
    ! local
    real(DP) :: slope

    ! chord sloping
    if (abs(h_new - h_old) > DPREC) then
      slope = (s_new - s_old) / (h_new - h_old)
    else
      slope = dsdh_lim
    end if

    aterm = -phi * area * thk * slope / delt
    rhsterm = aterm * h_old

  end subroutine get_unsat_storage_terms_CS

  subroutine get_unsat_storage_terms_MP(s_new, s_old, h_new, h_old, Cm, phi, &
                                        z_n, area, thk, aterm, rhsterm)
    use TdisModule, only: delt
    real(DP), intent(in) :: s_new
    real(DP), intent(in) :: s_old
    real(DP), intent(in) :: h_new
    real(DP), intent(in) :: h_old
    real(DP), intent(in) :: Cm
    real(DP), intent(in) :: phi
    real(DP), intent(in) :: z_n
    real(DP), intent(in) :: area
    real(DP), intent(in) :: thk
    real(DP), intent(inout) :: aterm
    real(DP), intent(inout) :: rhsterm
    ! local
    real(DP) :: rhs_ds, rhs_dh, a_dh

    ! part 1: change in s
    rhs_ds = phi * area * thk * (s_new - s_old) / delt

    ! part 2: change in h
    a_dh = -area * thk * Cm / delt
    rhs_dh = a_dh * h_new

    aterm = a_dh
    rhsterm = rhs_ds + rhs_dh

  end subroutine get_unsat_storage_terms_MP

end module UzrStorageModule
