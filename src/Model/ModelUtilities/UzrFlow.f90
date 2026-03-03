module UzrFlowModule
  use KindModule, only: I4B, LGP, DP
  use ConstantsModule, only: DONE, DTWO, DHALF, DZERO, LENVARNAME
  use MatrixBaseModule, only: MatrixBaseType
  use BaseDisModule, only: DisBaseType
  use GwfNpfModule, only: GwfNpfType
  use GwfNpfExtModule, only: GwfNpfExtType
  use UzrSoilModelModule, only: SoilModelType
  implicit none
  private

  public :: UzrFlowType
  public :: kr_averaging_name

  character(len=LENVARNAME), parameter, dimension(3) :: kr_averaging_name = &
      &[character(len=LENVARNAME) :: 'GEOMETRIC', 'ARITHMETIC', 'UPSTREAM']

  enum, bind(C)
    enumerator :: KR_GEOMETRIC = 1 !< Geometric mean of relative permeability
    enumerator :: KR_ARITHMETIC = 2 !< Arithmetic mean of relative permeability
    enumerator :: KR_UPSTREAM = 3 !< Upstream relative permeability
  end enum

  type, extends(GwfNpfExtType) :: UzrFlowType
    integer(I4B), pointer, dimension(:), contiguous :: iunsat => null() !< see UZR
    integer(I4B), pointer :: kr_averaging => null() !< see UZR
    real(DP), pointer, dimension(:), contiguous :: krel => null() !< pointer to NPF k_r
    class(SoilModelType), pointer :: soil_model => null() !< soil model used to get relative permeability
    class(DisBaseType), pointer :: gwf_dis => null()
    type(GwfNpfType), pointer :: gwf_npf => null()
  contains
    procedure :: initialize
    procedure :: is_active => uft_is_active
    procedure :: cf => uft_cf
    procedure :: fc => uft_fc
    procedure :: fn => uft_fn
    procedure :: cq => uft_cq
    procedure :: destroy
    ! private
    procedure, private :: calculate_coeffs
    procedure, private :: calculate_coeffs_nwt
  end type UzrFlowType

contains

  subroutine initialize(this, iunsat, kr_avg, soil_model, dis, npf)
    class(UzrFlowType), intent(inout) :: this
    integer(I4B), pointer, dimension(:), contiguous, intent(in) :: iunsat
    integer(I4B), pointer :: kr_avg
    class(SoilModelType), pointer :: soil_model
    class(DisBaseType), pointer :: dis
    type(GwfNpfType), pointer, intent(in) :: npf

    this%iunsat => iunsat
    this%kr_averaging => kr_avg
    this%soil_model => soil_model

    this%gwf_dis => dis
    this%gwf_npf => npf

    this%krel => this%gwf_npf%krel

  end subroutine initialize

  function uft_is_active(this, n, m) result(is_active)
    class(UzrFlowType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    logical(LGP) :: is_active

    is_active = .false.
    if (this%iunsat(n) == 1 .or. this%iunsat(m) == 1) then
      is_active = .true.
    end if

  end function uft_is_active

  subroutine uft_cf(this, kiter, n)
    class(UzrFlowType), intent(inout) :: this
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: n
    ! local
    real(DP) :: z_n, psi

    ! calculate k_r for node n
    z_n = DHALF * (this%gwf_dis%bot(n) + this%gwf_dis%top(n))
    psi = this%gwf_npf%hnew(n) - z_n
    this%krel(n) = this%soil_model%krelative(psi, n)

  end subroutine uft_cf

  subroutine uft_fc(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)
    class(UzrFlowType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: ipos
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: hnew
    ! local
    integer(I4B) :: idiag !< diagonal position
    integer(I4B) :: isym !< transposed connection
    integer(I4B) :: isymcon !< position of reverse connection m-n
    real(DP), dimension(3) :: coeffs !< the linear system coefficients for the flow

    if (this%gwf_npf%inewton == 0) then

      ! calculate system coefficients
      call this%calculate_coeffs(n, m, ipos, hnew, coeffs)

      ! Fill row n
      idiag = this%gwf_dis%con%ia(n)
      call matrix_sln%add_value_pos(idxglo(idiag), coeffs(1))
      call matrix_sln%add_value_pos(idxglo(ipos), coeffs(2))

      ! Fill row m
      isymcon = this%gwf_dis%con%isym(ipos)
      idiag = this%gwf_dis%con%ia(m)
      call matrix_sln%add_value_pos(idxglo(idiag), coeffs(1))
      call matrix_sln%add_value_pos(idxglo(isymcon), coeffs(2))
    else

      ! calculate system coefficients for newton formulation for n
      call this%calculate_coeffs_nwt(n, m, ipos, hnew, coeffs)

      ! Fill row n and add to RHS
      idiag = this%gwf_dis%con%ia(n)
      call matrix_sln%add_value_pos(idxglo(idiag), coeffs(1))
      call matrix_sln%add_value_pos(idxglo(ipos), coeffs(2))
      rhs(n) = rhs(n) + coeffs(3)

      ! calculate system coefficients for newton formulation for transposed: m
      isym = this%gwf_dis%con%isym(ipos)
      call this%calculate_coeffs_nwt(m, n, isym, hnew, coeffs)

      ! Fill row m and add to RHS
      idiag = this%gwf_dis%con%ia(m)
      call matrix_sln%add_value_pos(idxglo(idiag), coeffs(1))
      call matrix_sln%add_value_pos(idxglo(isym), coeffs(2))
      rhs(m) = rhs(m) + coeffs(3)

    end if

  end subroutine uft_fc

  subroutine calculate_coeffs(this, n, m, ipos, hnew, coeffs)
    class(UzrFlowType), intent(inout) :: this !, this instance
    integer(I4B), intent(in) :: n !< node n
    integer(I4B), intent(in) :: m !< node m
    integer(I4B), intent(in) :: ipos !< index in ja array for connection n-m
    real(DP), dimension(:), intent(in) :: hnew !< the new head
    real(DP), dimension(3), intent(inout) :: coeffs !< the coefficients for the linear system: A_nn, A_nm
    ! local
    real(DP) :: sat_cond !< conductance at full saturation
    real(DP) :: cond !< conductance at current saturation
    real(DP) :: kr_avg !< weighted rel. permeability between nodes

    coeffs(:) = DZERO

    sat_cond = this%gwf_npf%condsat(this%gwf_dis%con%jas(ipos))

    ! averaging of k_r
    kr_avg = kr_averaging(this%krel(n), this%krel(m), hnew(n), hnew(m), &
                          this%kr_averaging)

    ! calculate unsaturated conductance
    cond = kr_avg * sat_cond

    ! coefficients for row n
    coeffs(1) = -cond ! nn
    coeffs(2) = cond ! nm

  end subroutine calculate_coeffs

  subroutine calculate_coeffs_nwt(this, n, m, ipos, hnew, coeffs)
    class(UzrFlowType), intent(inout) :: this !, this instance
    integer(I4B), intent(in) :: n !< node n
    integer(I4B), intent(in) :: m !< node m
    integer(I4B), intent(in) :: ipos !< index in ja array for connection n-m
    real(DP), dimension(:), intent(in) :: hnew !< the new head
    real(DP), dimension(3), intent(inout) :: coeffs !< the coefficients for the linear system: A_nn, A_nm, rhs_n
    ! local
    integer(I4B) :: iup !< the upstream node
    real(DP) :: sat_cond !< conductance at full saturation
    real(DP) :: z_up !< the nodal elevation for n
    real(DP) :: psi !< the pressure head
    real(DP) :: kr_up !< upstream weighted rel. permeability between nodes
    real(DP) :: shift !< the shift for the numerical derivative
    real(DP) :: dkrdh_n, dkrdh_m !< numerical derivative of upstream weighted kr
    real(DP) :: Q_nm, dQdh_n, dQdh_m !< newton rates and derivatives

    coeffs(:) = DZERO

    sat_cond = this%gwf_npf%condsat(this%gwf_dis%con%jas(ipos))

    ! determine upstream
    if (hnew(n) > hnew(m)) then
      iup = n
    else
      iup = m
    end if

    ! upstream k_relative
    kr_up = this%krel(iup)

    ! dkr/dh
    z_up = DHALF * (this%gwf_dis%bot(iup) + this%gwf_dis%top(iup))
    psi = hnew(iup) - z_up
    shift = 1.0e-6
    if (iup == m) then
      dkrdh_n = DZERO
      dkrdh_m = (this%soil_model%krelative(psi + shift, iup) - &
                 this%soil_model%krelative(psi - shift, iup)) / (2.0 * shift)
    else
      dkrdh_n = (this%soil_model%krelative(psi + shift, iup) - &
                 this%soil_model%krelative(psi - shift, iup)) / (2.0 * shift)
      dkrdh_m = DZERO
    end if

    Q_nm = kr_up * sat_cond * (hnew(m) - hnew(n))
    dQdh_n = -kr_up * sat_cond
    dQdh_m = kr_up * sat_cond
    if (n == iup) then
      dQdh_n = dQdh_n + dkrdh_n * sat_cond * (hnew(m) - hnew(n))
    else
      dQdh_m = dQdh_m + dkrdh_m * sat_cond * (hnew(m) - hnew(n))
    end if

    coeffs(1) = dQdh_n ! diagonal
    coeffs(2) = dQdh_m ! off-diagonal
    coeffs(3) = -Q_nm + dQdh_n * hnew(n) + dQdh_m * hnew(m) ! RHS

  end subroutine calculate_coeffs_nwt

  subroutine uft_fn(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)
    class(UzrFlowType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: ipos
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: hnew

    ! TODO_UZR

  end subroutine uft_fn

  subroutine uft_cq(this, n, m, ipos, flowja, h_new)
    class(UzrFlowType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: ipos
    real(DP), dimension(:), intent(inout) :: flowja
    real(DP), dimension(:), intent(in) :: h_new
    ! local
    real(DP), dimension(3) :: coeffs !< the linear system coefficients: A_nn, A_nm, rhs_n
    real(DP) :: flow_nm !< the flow rate into node n from m

    if (this%gwf_npf%inewton == 0) then
      call this%calculate_coeffs(n, m, ipos, h_new, coeffs)
    else
      call this%calculate_coeffs_nwt(n, m, ipos, h_new, coeffs)
    end if

    flow_nm = coeffs(2) * h_new(m) + coeffs(1) * h_new(n) - coeffs(3)
    flowja(ipos) = flow_nm
    flowja(this%gwf_dis%con%isym(ipos)) = -flow_nm

  end subroutine uft_cq

  subroutine destroy(this)
    class(UzrFlowType) :: this

    this%gwf_dis => null()
    this%gwf_npf => null()
    this%iunsat => null()

  end subroutine destroy

  pure function kr_averaging(kr_n, kr_m, h_n, h_m, iavg) result(kr_avg)
    real(DP), intent(in) :: kr_n !< kr for node n
    real(DP), intent(in) :: kr_m !< kr for node m
    real(DP), intent(in) :: h_n !< h for node n
    real(DP), intent(in) :: h_m !< h for node m
    integer(I4B), intent(in) :: iavg !< averaging method
    real(DP) :: kr_avg !< averaged kr

    select case (iavg)
    case (KR_GEOMETRIC)
      kr_avg = sqrt(kr_n * kr_m)
    case (KR_ARITHMETIC)
      kr_avg = DHALF * (kr_n + kr_m)
    case (KR_UPSTREAM)
      if (h_n > h_m) then
        kr_avg = kr_n
      else
        kr_avg = kr_m
      end if
    case default
      kr_avg = sqrt(kr_n * kr_m)
    end select

  end function kr_averaging

end module UzrFlowModule
