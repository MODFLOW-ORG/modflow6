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
  contains
    procedure :: load => load_vangenuchten
    procedure :: saturation => saturation_vangenuchten
    procedure :: capacity => capacity_vangenuchten
    procedure :: krelative => krelative_vangenuchten
    procedure :: kr_dkrdh => kr_dkrdh_vangenuchten
    procedure :: sat_capacity => sat_capacity_vangenuchten
    procedure :: destroy => destroy_vangenuchten
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

    this%porosity => porosity
    this%sat_res => sat_res

    call mem_allocate(this%alpha, this%nodes, 'ALPHAVGN', this%mem_path)
    call mem_allocate(this%n, this%nodes, 'NVGN', this%mem_path)

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
    if (count_errors() > 0) call ustop()

  end subroutine load_vangenuchten

  function saturation_vangenuchten(this, psi, i) result(s)
    class(VanGenuchtenModelType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: s !< the saturation
    ! local
    real(DP) :: s_eff
    real(DP) :: m

    m = DONE - DONE / this%n(i)

    if (psi > DZERO) then
      s = DONE
    else
      s_eff = DONE / ((DONE + ((-this%alpha(i) * psi)**this%n(i)))**m)
      s = (DONE - this%sat_res(i)) * s_eff + this%sat_res(i)
    end if

  end function saturation_vangenuchten

  function capacity_vangenuchten(this, psi, i) result(Cm)
    class(VanGenuchtenModelType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: Cm !< the saturation
    ! local
    real(DP) :: num, denom
    real(DP) :: m

    m = DONE - DONE / this%n(i)

    if (psi > DZERO) then
      Cm = DZERO
    else
      num = m * this%n(i) * this%alpha(i) * &
            ((-this%alpha(i) * psi)**(this%n(i) - DONE))
      denom = (DONE + (-this%alpha(i) * psi)**this%n(i))**(m + DONE)
      Cm = this%porosity(i) * (DONE - this%sat_res(i)) * num / denom
    end if

  end function capacity_vangenuchten

  function krelative_vangenuchten(this, psi, i) result(kr)
    class(VanGenuchtenModelType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: kr !< the saturation
    ! local
    real(DP) :: se !< effective saturation
    real(DP) :: term
    real(DP) :: m

    m = DONE - DONE / this%n(i)

    if (psi > DZERO) then
      kr = DONE
    else
      se = (this%saturation(psi, i) - this%sat_res(i)) / (DONE - this%sat_res(i))
      term = DONE - se**(DONE / m)
      kr = sqrt(se) * ((DONE - term**m)**DTWO)
      ! TODO_VSF: formalize this:
      if (.false.) then
        kr = kr_smoothing(kr, se, 0.999_DP, 0.0001_DP)
      end if
    end if

  end function krelative_vangenuchten

  !> @brief Fused relative permeability and its analytical derivative
  !!
  !! With x = (alpha |psi|)^n: se = (1 + x)^-m and
  !! (1 - se^(1/m))^m = (x / (1 + x))^m, so kr and dkr/dpsi need only
  !! three exponentials and two logarithms. Close to saturation, where
  !! dkr/dpsi becomes singular for n < 2, the central difference of the
  !! default implementation is kept.
  !<
  subroutine kr_dkrdh_vangenuchten(this, psi, i, kr, dkr)
    class(VanGenuchtenModelType), intent(inout) :: this
    real(DP), intent(in) :: psi !< pressure head
    integer(I4B), intent(in) :: i !< node number
    real(DP), intent(out) :: kr !< relative permeability
    real(DP), intent(out) :: dkr !< d(kr)/d(psi)
    ! local
    real(DP) :: nv, m, a, lx, x, l1x, se, sqse, f, dse

    if (psi > -DKR_SHIFT) then
      kr = this%krelative(psi, i)
      dkr = (this%krelative(psi + DKR_SHIFT, i) - &
             this%krelative(psi - DKR_SHIFT, i)) / (DTWO * DKR_SHIFT)
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
    real(DP) :: nv, m, a, x, se

    if (psi >= DZERO) then
      s = DONE
      cm = DZERO
      return
    end if

    nv = this%n(i)
    m = DONE - DONE / nv
    a = -this%alpha(i) * psi
    x = exp(nv * log(a))
    se = exp(-m * log(DONE + x))
    s = (DONE - this%sat_res(i)) * se + this%sat_res(i)
    cm = this%porosity(i) * (DONE - this%sat_res(i)) * &
         m * nv * this%alpha(i) * (x / a) * se / (DONE + x)

  end subroutine sat_capacity_vangenuchten

  function kr_smoothing(kr, s, s0, omega) result(kr_smooth)
    real(DP) :: kr
    real(DP) :: s
    real(DP) :: s0
    real(DP) :: omega
    real(DP) :: kr_smooth
    ! local
    real(DP) :: t !< transition function

    t = DHALF + DHALF * tanh((s - s0) / omega)
    kr_smooth = (DONE - t) * kr + t * DONE

  end function kr_smoothing

  !> @brief clean up
  !<
  subroutine destroy_vangenuchten(this)
    class(VanGenuchtenModelType), intent(inout) :: this !< this instance

    this%porosity => null()
    this%sat_res => null()

    if (associated(this%alpha)) then
      call mem_deallocate(this%alpha)
      call mem_deallocate(this%n)
    end if

    call this%destroy_base()

  end subroutine destroy_vangenuchten

end module VsfVanGenuchtenModule
