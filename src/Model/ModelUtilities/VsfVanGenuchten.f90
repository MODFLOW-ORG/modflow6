module VsfVanGenuchtenModule
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE, DTWO, DHALF
  use SimModule, only: store_error, count_errors, ustop
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use MemoryManagerExtModule, only: mem_set_value
  use VsfSoilModelModule, only: SoilModelType, DKR_SHIFT
  use GwfVsfInputModule, only: GwfVsfParamFoundType
  implicit none
  private

  public :: VanGenuchtenModelType

  ! TODO_VSF: document these parameters based on the tech ref
  type, extends(SoilModelType) :: VanGenuchtenModelType
    real(DP), dimension(:), pointer, contiguous :: alpha => null()
    real(DP), dimension(:), pointer, contiguous :: n => null()
    real(DP), dimension(:), pointer, contiguous :: hae => null() !< air-entry pressure head (<= 0), 0: none
    real(DP), pointer :: delta => null() !< width of the kr smoothing interval below saturation, 0: none
    real(DP), dimension(:), pointer, contiguous :: sc => null() !< unscaled effective saturation at the air-entry head
    real(DP), dimension(:), pointer, contiguous :: krscale => null() !< unscaled kr at the air-entry head
    real(DP), dimension(:), pointer, contiguous :: krd => null() !< unscaled kr at -delta
    real(DP), dimension(:), pointer, contiguous :: dkrd => null() !< unscaled dkr/dpsi at -delta
  contains
    procedure :: load => load_vangenuchten
    procedure :: saturation => saturation_vangenuchten
    procedure :: capacity => capacity_vangenuchten
    procedure :: krelative => krelative_vangenuchten
    procedure :: kr_dkrdh => kr_dkrdh_vangenuchten
    procedure :: sat_capacity => sat_capacity_vangenuchten
    procedure :: destroy => destroy_vangenuchten
    procedure, private :: kr_raw
    procedure, private :: kr_dkr_raw
    procedure, private :: se_dse_raw
  end type VanGenuchtenModelType

contains

  !> @brief Load the Van Genuchten model from input data
  !<
  subroutine load_vangenuchten(this, porosity, sat_res)
    class(VanGenuchtenModelType), intent(inout) :: this !< this instance
    real(DP), dimension(:), pointer, contiguous :: porosity !< porosity
    real(DP), dimension(:), pointer, contiguous :: sat_res !< residual saturation
    ! local
    type(GwfVsfParamFoundType) :: found
    integer(I4B) :: i, nae, nsm
    real(DP) :: se, dse, kr, dkr, f, m

    this%porosity => porosity
    this%sat_res => sat_res

    call mem_allocate(this%alpha, this%nodes, 'ALPHAVGN', this%mem_path)
    call mem_allocate(this%n, this%nodes, 'NVGN', this%mem_path)
    call mem_allocate(this%hae, this%nodes, 'HAEVGN', this%mem_path)
    call mem_allocate(this%delta, 'KRSMOOTHING', this%mem_path)
    this%delta = DZERO

    ! load alpha and n
    call mem_set_value(this%alpha, 'ALPHAVGN', this%input_mem_path, &
                       this%map, found%alphavgn)
    call mem_set_value(this%n, 'NVGN', this%input_mem_path, &
                       this%map, found%nvgn)
    if (.not. found%alphavgn) then
      call store_error('Van Genuchten soil data (ALPHAVGN) &
        &not specified in input', .false.)
    end if
    if (.not. found%nvgn) then
      call store_error('Van Genuchten soil data (NVGN) &
        &not specified in input', .false.)
    end if

    ! optional air-entry head and kr smoothing
    do i = 1, this%nodes
      this%hae(i) = DZERO
    end do
    call mem_set_value(this%hae, 'HAEVGN', this%input_mem_path, &
                       this%map, found%haevgn)
    call mem_set_value(this%delta, 'KR_SMOOTHING', this%input_mem_path, &
                       found%kr_smoothing)
    if (any(this%hae > DZERO)) then
      call store_error('Van Genuchten air-entry head (HAEVGN) must not be &
        &positive', .false.)
    end if
    if (this%delta < DZERO) then
      call store_error('KR_SMOOTHING must not be negative', .false.)
    end if
    if (count_errors() > 0) call ustop()

    ! precompute the scaling at the air-entry head and the values at the
    ! start of the smoothing interval
    nae = 0
    if (any(this%hae < DZERO)) nae = this%nodes
    nsm = 0
    if (this%delta > DZERO) nsm = this%nodes
    call mem_allocate(this%sc, nae, 'SCVGN', this%mem_path)
    call mem_allocate(this%krscale, nae, 'KRSCALEVGN', this%mem_path)
    call mem_allocate(this%krd, nsm, 'KRDVGN', this%mem_path)
    call mem_allocate(this%dkrd, nsm, 'DKRDVGN', this%mem_path)
    do i = 1, nae
      if (this%hae(i) < DZERO) then
        call this%se_dse_raw(this%hae(i), i, se, dse)
        m = DONE - DONE / this%n(i)
        f = (DONE - se**(DONE / m))**m
        this%sc(i) = se
        this%krscale(i) = sqrt(se) * (DONE - f)**2
      else
        this%sc(i) = DONE
        this%krscale(i) = DONE
      end if
    end do
    do i = 1, nsm
      call this%kr_dkr_raw(-this%delta, i, kr, dkr)
      this%krd(i) = kr
      this%dkrd(i) = dkr
    end do

  end subroutine load_vangenuchten

  function saturation_vangenuchten(this, psi, i) result(s)
    class(VanGenuchtenModelType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: s !< the saturation
    ! local
    real(DP) :: cm

    call this%sat_capacity(psi, i, s, cm)

  end function saturation_vangenuchten

  function capacity_vangenuchten(this, psi, i) result(Cm)
    class(VanGenuchtenModelType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: Cm !< the moisture capacity
    ! local
    real(DP) :: s

    call this%sat_capacity(psi, i, s, Cm)

  end function capacity_vangenuchten

  function krelative_vangenuchten(this, psi, i) result(kr)
    class(VanGenuchtenModelType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: kr !< the relative permeability
    ! local
    real(DP) :: dkr

    call this%kr_dkrdh(psi, i, kr, dkr)

  end function krelative_vangenuchten

  !> @brief Mualem-Van Genuchten relative permeability, without
  !< regularisation
  function kr_raw(this, psi, i) result(kr)
    class(VanGenuchtenModelType), intent(inout) :: this
    real(DP), intent(in) :: psi !< pressure head
    integer(I4B), intent(in) :: i !< node number
    real(DP) :: kr
    ! local
    real(DP) :: nv, m, lx, l1x, se, f

    if (psi >= DZERO) then
      kr = DONE
      return
    end if
    nv = this%n(i)
    m = DONE - DONE / nv
    lx = nv * log(-this%alpha(i) * psi)
    l1x = log(DONE + exp(lx))
    se = exp(-m * l1x)
    f = exp(m * (lx - l1x))
    kr = sqrt(se) * (DONE - f)**2

  end function kr_raw

  !> @brief Mualem-Van Genuchten relative permeability and its derivative,
  !! without regularisation
  !!
  !! With x = (alpha |psi|)^n: se = (1 + x)^-m and
  !! (1 - se^(1/m))^m = (x / (1 + x))^m, so kr and dkr/dpsi need only
  !! three exponentials and two logarithms. Close to saturation, where
  !! dkr/dpsi becomes singular for n < 2, a central difference is used.
  !<
  subroutine kr_dkr_raw(this, psi, i, kr, dkr)
    class(VanGenuchtenModelType), intent(inout) :: this
    real(DP), intent(in) :: psi !< pressure head
    integer(I4B), intent(in) :: i !< node number
    real(DP), intent(out) :: kr !< relative permeability
    real(DP), intent(out) :: dkr !< d(kr)/d(psi)
    ! local
    real(DP) :: nv, m, a, lx, x, l1x, se, sqse, f, dse

    if (psi > -DKR_SHIFT) then
      kr = this%kr_raw(psi, i)
      dkr = (this%kr_raw(psi + DKR_SHIFT, i) - &
             this%kr_raw(psi - DKR_SHIFT, i)) / (DTWO * DKR_SHIFT)
      return
    end if

    nv = this%n(i)
    m = DONE - DONE / nv
    a = -this%alpha(i) * psi
    lx = nv * log(a)
    x = exp(lx)
    l1x = log(DONE + x)
    se = exp(-m * l1x)
    f = exp(m * (lx - l1x))
    sqse = sqrt(se)
    kr = sqse * (DONE - f)**2
    ! d(se)/d(psi) = m n alpha (alpha |psi|)^(n-1) (1 + x)^-(m+1)
    dse = m * nv * this%alpha(i) * (x / a) * se / (DONE + x)
    dkr = dse * (DHALF * (DONE - f)**2 + DTWO * (DONE - f) * f / x) / sqse

  end subroutine kr_dkr_raw

  !> @brief Van Genuchten effective saturation and its derivative, without
  !< regularisation
  subroutine se_dse_raw(this, psi, i, se, dse)
    class(VanGenuchtenModelType), intent(inout) :: this
    real(DP), intent(in) :: psi !< pressure head
    integer(I4B), intent(in) :: i !< node number
    real(DP), intent(out) :: se !< effective saturation
    real(DP), intent(out) :: dse !< d(se)/d(psi)
    ! local
    real(DP) :: nv, m, a, x

    if (psi >= DZERO) then
      se = DONE
      dse = DZERO
      return
    end if
    nv = this%n(i)
    m = DONE - DONE / nv
    a = -this%alpha(i) * psi
    x = exp(nv * log(a))
    se = exp(-m * log(DONE + x))
    dse = m * nv * this%alpha(i) * (x / a) * se / (DONE + x)

  end subroutine se_dse_raw

  !> @brief Fused relative permeability and its derivative
  !!
  !! With an air-entry head hae < 0 (Vogel et al., 2001) the soil is
  !! saturated for psi >= hae and kr is the Mualem integral rescaled to 1
  !! at hae. Otherwise, with a smoothing width delta > 0, kr on (-delta, 0)
  !! is a cubic matching kr and its slope at -delta and reaching 1 with
  !! zero slope at saturation.
  !<
  subroutine kr_dkrdh_vangenuchten(this, psi, i, kr, dkr)
    class(VanGenuchtenModelType), intent(inout) :: this
    real(DP), intent(in) :: psi !< pressure head
    integer(I4B), intent(in) :: i !< node number
    real(DP), intent(out) :: kr !< relative permeability
    real(DP), intent(out) :: dkr !< d(kr)/d(psi)
    ! local
    real(DP) :: t, t2, t3, m0

    if (this%hae(i) < DZERO) then
      if (psi >= this%hae(i)) then
        kr = DONE
        dkr = DZERO
      else
        call this%kr_dkr_raw(psi, i, kr, dkr)
        kr = kr / this%krscale(i)
        dkr = dkr / this%krscale(i)
      end if
    else if (this%delta > DZERO .and. psi > -this%delta) then
      if (psi >= DZERO) then
        kr = DONE
        dkr = DZERO
      else
        ! cubic Hermite on t in [0, 1): kr(-delta), slope, then 1 and 0
        t = (psi + this%delta) / this%delta
        t2 = t * t
        t3 = t2 * t
        m0 = this%dkrd(i) * this%delta
        kr = (2.0_DP * t3 - 3.0_DP * t2 + DONE) * this%krd(i) + &
             (t3 - DTWO * t2 + t) * m0 + (-2.0_DP * t3 + 3.0_DP * t2)
        dkr = ((6.0_DP * t2 - 6.0_DP * t) * this%krd(i) + &
               (3.0_DP * t2 - 4.0_DP * t + DONE) * m0 + &
               (-6.0_DP * t2 + 6.0_DP * t)) / this%delta
      end if
    else
      call this%kr_dkr_raw(psi, i, kr, dkr)
    end if

  end subroutine kr_dkrdh_vangenuchten

  !> @brief Fused saturation and moisture capacity
  !<
  subroutine sat_capacity_vangenuchten(this, psi, i, s, cm)
    class(VanGenuchtenModelType), intent(inout) :: this
    real(DP), intent(in) :: psi !< pressure head
    integer(I4B), intent(in) :: i !< node number
    real(DP), intent(out) :: s !< saturation
    real(DP), intent(out) :: cm !< moisture capacity
    ! local
    real(DP) :: se, dse

    if (this%hae(i) < DZERO) then
      if (psi >= this%hae(i)) then
        se = DONE
        dse = DZERO
      else
        call this%se_dse_raw(psi, i, se, dse)
        se = se / this%sc(i)
        dse = dse / this%sc(i)
      end if
    else
      call this%se_dse_raw(psi, i, se, dse)
    end if
    s = (DONE - this%sat_res(i)) * se + this%sat_res(i)
    cm = this%porosity(i) * (DONE - this%sat_res(i)) * dse

  end subroutine sat_capacity_vangenuchten

  !> @brief clean up
  !<
  subroutine destroy_vangenuchten(this)
    class(VanGenuchtenModelType), intent(inout) :: this !< this instance

    this%porosity => null()
    this%sat_res => null()

    if (associated(this%alpha)) then
      call mem_deallocate(this%alpha)
      call mem_deallocate(this%n)
      call mem_deallocate(this%hae)
      call mem_deallocate(this%delta)
      call mem_deallocate(this%sc)
      call mem_deallocate(this%krscale)
      call mem_deallocate(this%krd)
      call mem_deallocate(this%dkrd)
    end if

    call this%destroy_base()

  end subroutine destroy_vangenuchten

end module VsfVanGenuchtenModule
