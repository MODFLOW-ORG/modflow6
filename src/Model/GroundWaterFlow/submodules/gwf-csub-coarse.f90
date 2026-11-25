!> @brief CSUB package coarse-grained storage
!<
submodule(GwfCsubModule) GwfCsubModuleCoarse
contains

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_cg_fc(this, node, tled, area, hcell, hcellold, hcof, &
                               rhs)
    ! -- dummy variables
    class(GwfCsubType) :: this
    integer(I4B), intent(in) :: node !< cell node number
    real(DP), intent(in) :: tled !< recripicol of the time step length
    real(DP), intent(in) :: area !< horizontal cell area
    real(DP), intent(in) :: hcell !< current head
    real(DP), intent(in) :: hcellold !< previous head
    real(DP), intent(inout) :: hcof !< coarse-grained A matrix entry
    real(DP), intent(inout) :: rhs !< coarse-grained right-hand side entry
    ! -- local variables
    real(DP) :: top
    real(DP) :: bot
    real(DP) :: tthk
    real(DP) :: snold
    real(DP) :: snnew
    real(DP) :: hbar
    real(DP) :: sske
    real(DP) :: rho1
    !
    ! -- initialize variables
    rhs = DZERO
    hcof = DZERO
    !
    ! -- aquifer elevations and thickness
    top = this%dis%top(node)
    bot = this%dis%bot(node)
    tthk = this%cg_thickini(node)
    !
    ! -- calculate hcof and rhs terms if coarse-grained materials present
    if (tthk > DZERO) then
      !
      ! -- calculate aquifer saturation
      call this%csub_calc_sat(node, hcell, hcellold, snnew, snold)
      !
      ! -- calculate corrected head (hbar)
      hbar = sQuadratic0sp(hcell, bot, this%satomega)
      !
      ! -- storage coefficients
      call this%csub_cg_calc_sske(node, sske, hcell)
      rho1 = sske * area * tthk * tled
      !
      ! -- update sk and ske
      this%cg_ske(node) = sske * tthk * snold
      this%cg_sk(node) = sske * tthk * snnew
      !
      ! -- calculate hcof and rhs term
      hcof = -rho1 * snnew
      rhs = rho1 * snold * this%cg_es0(node) - &
            rho1 * snnew * (this%cg_gs(node) + bot)
      !
      ! -- calculate and apply the flow correction term
      rhs = rhs - rho1 * snnew * (hcell - hbar)
    end if
  end subroutine csub_cg_fc

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_cg_fn(this, node, tled, area, hcell, hcof, rhs)
    ! -- dummy variables
    class(GwfCsubType) :: this
    integer(I4B), intent(in) :: node !< node number
    real(DP), intent(in) :: tled !< reciprocal of the time step length
    real(DP), intent(in) :: area !< horizontal cell area
    real(DP), intent(in) :: hcell !< current head in cell
    real(DP), intent(inout) :: hcof !< coarse-grained A matrix entry
    real(DP), intent(inout) :: rhs !< coarse-grained right-hand side entry
    ! -- local variables
    real(DP) :: top
    real(DP) :: bot
    real(DP) :: tthk
    real(DP) :: snnew
    real(DP) :: snold
    real(DP) :: satderv
    real(DP) :: hbar
    real(DP) :: hbarderv
    real(DP) :: sske
    real(DP) :: rho1
    !
    ! -- initialize variables
    rhs = DZERO
    hcof = DZERO
    !
    ! -- aquifer elevations and thickness
    top = this%dis%top(node)
    bot = this%dis%bot(node)
    tthk = this%cg_thickini(node)
    !
    ! -- calculate newton terms if coarse-grained materials present
    if (tthk > DZERO) then
      !
      ! -- calculate aquifer saturation - only need snnew
      call this%csub_calc_sat(node, hcell, top, snnew, snold)
      !
      ! -- calculate saturation derivative
      satderv = this%csub_calc_sat_derivative(node, hcell)
      !
      ! -- calculate corrected head (hbar)
      hbar = sQuadratic0sp(hcell, bot, this%satomega)
      !
      ! -- calculate the derivative of the hbar functions
      hbarderv = sQuadratic0spDerivative(hcell, bot, this%satomega)
      !
      ! -- storage coefficients
      call this%csub_cg_calc_sske(node, sske, hcell)
      rho1 = sske * area * tthk * tled
      !
      ! -- calculate hcof term
      hcof = rho1 * snnew * (DONE - hbarderv) + &
             rho1 * (this%cg_gs(node) - hbar + bot) * satderv
      !
      ! -- Add additional term if using lagged effective stress
      if (this%ieslag /= 0) then
        hcof = hcof - rho1 * this%cg_es0(node) * satderv
      end if
      !
      ! -- calculate rhs term
      rhs = hcof * hcell
    end if
  end subroutine csub_cg_fn

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_cg_calc_sske(this, n, sske, hcell)
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    integer(I4B), intent(in) :: n !< cell node number
    real(DP), intent(inout) :: sske !< coarse grained Sske
    real(DP), intent(in) :: hcell !< current head in cell
    ! -- local variables
    real(DP) :: top
    real(DP) :: bot
    real(DP) :: hbar
    real(DP) :: znode
    real(DP) :: es
    real(DP) :: es0
    real(DP) :: theta
    real(DP) :: f
    real(DP) :: f0
    !
    ! -- initialize variables
    sske = DZERO
    !
    ! -- calculate factor for the head-based case
    if (this%lhead_based .EQV. .TRUE.) then
      f = DONE
      f0 = DONE
      !
      ! -- calculate factor for the effective stress case
    else
      top = this%dis%top(n)
      bot = this%dis%bot(n)
      !
      ! -- calculate corrected head (hbar)
      hbar = sQuadratic0sp(hcell, bot, this%satomega)
      !
      ! -- calculate znode
      znode = this%csub_calc_znode(top, bot, hbar)
      !
      ! -- calculate effective stress and theta
      es = this%cg_es(n)
      es0 = this%cg_es0(n)
      theta = this%cg_thetaini(n)
      !
      ! -- calculate the compression index factors for the delay
      !    node relative to the center of the cell based on the
      !    current and previous head
      call this%csub_calc_sfacts(n, bot, znode, theta, es, es0, &
                                 this%cg_gs(n), f)
    end if
    sske = f * this%cg_ske_cr(n)
  end subroutine csub_cg_calc_sske

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_cg_calc_comp(this, node, hcell, hcellold, comp)
    ! -- dummy variables
    class(GwfCsubType) :: this
    integer(I4B), intent(in) :: node !< cell node number
    real(DP), intent(in) :: hcell !< current head in cell
    real(DP), intent(in) :: hcellold !< previous head in cell
    real(DP), intent(inout) :: comp !< coarse-grained compaction
    ! -- local variables
    real(DP) :: area
    real(DP) :: tled
    real(DP) :: hcof
    real(DP) :: rhs
    !
    ! -- initialize variables
    area = DONE
    tled = DONE
    !
    ! -- calculate terms
    call this%csub_cg_fc(node, tled, area, hcell, hcellold, hcof, rhs)
    !
    ! - calculate compaction
    comp = hcof * hcell - rhs
  end subroutine csub_cg_calc_comp

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_cg_update(this, node)
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    integer(I4B), intent(in) :: node !< cell node number
    ! -- local variables
    character(len=20) :: cellid
    real(DP) :: comp
    real(DP) :: thick
    real(DP) :: theta
    !
    ! -- update thickness and theta
    comp = this%cg_tcomp(node) + this%cg_comp(node)
    call this%dis%noder_to_string(node, cellid)
    if (ABS(comp) > DZERO) then
      thick = this%cg_thickini(node)
      theta = this%cg_thetaini(node)
      call compact_layer(comp, thick, theta)
      if (thick <= DZERO) then
        write (errmsg, '(a,1x,a,1x,a,g0,a)') &
          'Adjusted thickness for cell', trim(adjustl(cellid)), &
          'is less than or equal to 0 (', thick, ').'
        call store_error(errmsg)
      end if
      if (theta <= DZERO) then
        write (errmsg, '(a,1x,a,1x,a,g0,a)') &
          'Adjusted theta for cell', trim(adjustl(cellid)), &
          'is less than or equal to 0 (', theta, ').'
        call store_error(errmsg)
      end if
      this%cg_thick(node) = thick
      this%cg_theta(node) = theta
    end if
  end subroutine csub_cg_update

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_cg_wcomp_fc(this, node, tled, area, hcell, hcellold, &
                                     hcof, rhs)
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    integer(I4B), intent(in) :: node !< cell node number
    real(DP), intent(in) :: tled !< reciprocal of the time step length
    real(DP), intent(in) :: area !< horizontal cell area
    real(DP), intent(in) :: hcell !< current head in cell
    real(DP), intent(in) :: hcellold !< previous head in cell
    real(DP), intent(inout) :: hcof !< coarse-grained A matrix entry
    real(DP), intent(inout) :: rhs !< coarse-grained right-hand side entry
    ! -- local variables
    real(DP) :: top
    real(DP) :: bot
    real(DP) :: tthk
    real(DP) :: tthk0
    real(DP) :: snold
    real(DP) :: snnew
    real(DP) :: wc
    real(DP) :: wc0
    !
    ! -- initialize variables
    rhs = DZERO
    hcof = DZERO
    !
    ! -- aquifer elevations and thickness
    top = this%dis%top(node)
    bot = this%dis%bot(node)
    tthk = this%cg_thick(node)
    tthk0 = this%cg_thick0(node)
    !
    ! -- aquifer saturation
    call this%csub_calc_sat(node, hcell, hcellold, snnew, snold)
    !
    ! -- storage coefficients
    wc0 = this%brg * area * tthk0 * this%cg_theta0(node) * tled
    wc = this%brg * area * tthk * this%cg_theta(node) * tled
    !
    ! -- calculate hcof term
    hcof = -wc * snnew
    !
    ! -- calculate rhs term
    rhs = -wc0 * snold * hcellold
  end subroutine csub_cg_wcomp_fc

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_cg_wcomp_fn(this, node, tled, area, hcell, hcellold, &
                                     hcof, rhs)
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    integer(I4B), intent(in) :: node !< cell node number
    real(DP), intent(in) :: tled !< reciprocal of the time step length
    real(DP), intent(in) :: area !< horizontal cell area
    real(DP), intent(in) :: hcell !< current head in cell
    real(DP), intent(in) :: hcellold !< previous head in cell
    real(DP), intent(inout) :: hcof !< coarse-grained A matrix entry
    real(DP), intent(inout) :: rhs !< coarse-grained right-hand side entry
    ! -- local variables
    real(DP) :: top
    real(DP) :: bot
    real(DP) :: tthk
    real(DP) :: tthk0
    real(DP) :: satderv
    real(DP) :: f
    real(DP) :: wc
    real(DP) :: wc0
    !
    ! -- initialize variables
    rhs = DZERO
    hcof = DZERO
    !
    ! -- aquifer elevations and thickness
    top = this%dis%top(node)
    bot = this%dis%bot(node)
    tthk = this%cg_thick(node)
    !
    ! -- calculate saturation derivative
    satderv = this%csub_calc_sat_derivative(node, hcell)
    !
    ! -- calculate water compressibility factor
    f = this%brg * area * tled
    !
    ! -- water compressibility coefficient
    wc = f * tthk * this%cg_theta(node)
    !
    ! -- calculate hcof term
    hcof = -wc * hcell * satderv
    !
    ! -- Add additional term if using lagged effective stress
    if (this%ieslag /= 0) then
      tthk0 = this%cg_thick0(node)
      wc0 = f * tthk0 * this%cg_theta0(node)
      hcof = hcof + wc * hcellold * satderv
    end if
    !
    ! -- calculate rhs term
    rhs = hcof * hcell
  end subroutine csub_cg_wcomp_fn

end submodule GwfCsubModuleCoarse
