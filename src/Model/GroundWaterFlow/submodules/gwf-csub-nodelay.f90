!> @brief CSUB package no-delay interbed storage
!<
submodule(GwfCsubModule) GwfCsubModuleNodelay
contains

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_nodelay_update(this, i)
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    integer(I4B), intent(in) :: i
    ! -- local variables
    real(DP) :: comp
    real(DP) :: thick
    real(DP) :: theta
    !
    ! -- update thickness and theta
    comp = this%tcomp(i) + this%comp(i)
    if (ABS(comp) > DZERO) then
      thick = this%thickini(i)
      theta = this%thetaini(i)
      call compact_layer(comp, thick, theta)
      if (thick <= DZERO) then
        write (errmsg, '(a,1x,i0,1x,a,g0,a)') &
          'Adjusted thickness for no-delay interbed', i, &
          'is less than or equal to 0 (', thick, ').'
        call store_error(errmsg)
      end if
      if (theta <= DZERO) then
        write (errmsg, '(a,1x,i0,1x,a,g0,a)') &
          'Adjusted theta for no-delay interbed', i, &
          'is less than or equal to 0 (', theta, ').'
        call store_error(errmsg)
      end if
      this%thick(i) = thick
      this%theta(i) = theta
    end if
  end subroutine csub_nodelay_update

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_nodelay_fc(this, ib, hcell, hcellold, rho1, rho2, &
                                    rhs, argtled)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy variables
    class(GwfCsubType) :: this
    integer(I4B), intent(in) :: ib !< interbed number
    real(DP), intent(in) :: hcell !< current head in the cell
    real(DP), intent(in) :: hcellold !< previous head in the cell
    real(DP), intent(inout) :: rho1 !< current storage coefficient value using Sske
    real(DP), intent(inout) :: rho2 !< current storage coefficient value based on Ssk
    real(DP), intent(inout) :: rhs !< no-delay interbed contribution to the right-hand side
    real(DP), intent(in), optional :: argtled !< optional reciprocal of the time step length
    ! -- local variables
    integer(I4B) :: node
    real(DP) :: tled
    real(DP) :: top
    real(DP) :: bot
    real(DP) :: thick
    real(DP) :: hbar
    real(DP) :: znode
    real(DP) :: snold
    real(DP) :: snnew
    real(DP) :: sto_fac
    real(DP) :: sto_fac0
    real(DP) :: area
    real(DP) :: theta
    real(DP) :: es
    real(DP) :: es0
    real(DP) :: f
    real(DP) :: f0
    real(DP) :: rcorr
    !
    ! -- process optional variables
    if (present(argtled)) then
      tled = argtled
    else
      tled = DONE / delt
    end if
    node = this%nodelist(ib)
    area = this%dis%get_area(node)
    bot = this%dis%bot(node)
    top = this%dis%top(node)
    thick = this%thickini(ib)
    !
    ! -- calculate corrected head (hbar)
    hbar = sQuadratic0sp(hcell, bot, this%satomega)
    !
    ! -- set iconvert
    this%iconvert(ib) = 0
    !
    ! -- aquifer saturation
    call this%csub_calc_sat(node, hcell, hcellold, snnew, snold)
    if (this%lhead_based .EQV. .TRUE.) then
      f = DONE
      f0 = DONE
    else
      znode = this%csub_calc_znode(top, bot, hbar)
      es = this%cg_es(node)
      es0 = this%cg_es0(node)
      theta = this%thetaini(ib)
      !
      ! -- calculate the compression index factors for the delay
      !    node relative to the center of the cell based on the
      !    current and previous head
      call this%csub_calc_sfacts(node, bot, znode, theta, es, es0, &
                                 this%cg_gs(node), f)
    end if
    sto_fac = tled * snnew * thick * f
    sto_fac0 = tled * snold * thick * f
    !
    ! -- calculate rho1 and rho2
    rho1 = this%rci(ib) * sto_fac0
    rho2 = this%rci(ib) * sto_fac
    if (this%cg_es(node) > this%pcs(ib)) then
      this%iconvert(ib) = 1
      rho2 = this%ci(ib) * sto_fac
    end if
    !
    ! -- calculate correction term
    rcorr = rho2 * (hcell - hbar)
    !
    ! -- fill right-hand side
    if (this%ielastic(ib) /= 0) then
      rhs = rho1 * this%cg_es0(node) - &
            rho2 * (this%cg_gs(node) + bot) - &
            rcorr
    else
      rhs = -rho2 * (this%cg_gs(node) + bot) + &
            (this%pcs(ib) * (rho2 - rho1)) + &
            (rho1 * this%cg_es0(node)) - &
            rcorr
    end if
    !
    ! -- save ske and sk
    this%ske(ib) = rho1
    this%sk(ib) = rho2
  end subroutine csub_nodelay_fc

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_nodelay_calc_comp(this, ib, hcell, hcellold, comp, &
                                           rho1, rho2)
    ! -- dummy variables
    class(GwfCsubType) :: this
    integer(I4B), intent(in) :: ib !< interbed number
    real(DP), intent(in) :: hcell !< current head for the cell
    real(DP), intent(in) :: hcellold !< previous head for the cell
    real(DP), intent(inout) :: comp !< no-delay interbed compaction
    real(DP), intent(inout) :: rho1 !< current storage coefficient based on Sske
    real(DP), intent(inout) :: rho2 !< current storage coefficient based on Ssk
    ! -- local variables
    integer(I4B) :: node
    real(DP) :: es
    real(DP) :: es0
    real(DP) :: pcs
    real(DP) :: tled
    real(DP) :: rhs
    !
    ! -- initialize variables
    node = this%nodelist(ib)
    tled = DONE
    es = this%cg_es(node)
    es0 = this%cg_es0(node)
    pcs = this%pcs(ib)
    !
    ! -- calculate no-delay interbed rho1 and rho2
    call this%csub_nodelay_fc(ib, hcell, hcellold, rho1, rho2, rhs, argtled=tled)
    !
    ! -- calculate no-delay interbed compaction
    if (this%ielastic(ib) /= 0) then
      comp = rho2 * es - rho1 * es0
    else
      comp = -pcs * (rho2 - rho1) - (rho1 * es0) + (rho2 * es)
    end if
  end subroutine csub_nodelay_calc_comp

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_interbed_fc(this, ib, node, area, hcell, hcellold, &
                                     hcof, rhs)
    ! -- dummy variables
    class(GwfCsubType) :: this
    integer(I4B), intent(in) :: ib !< interbed number
    integer(I4B), intent(in) :: node !< cell node number
    real(DP), intent(in) :: area !< horizontal cell area
    real(DP), intent(in) :: hcell !< current head in cell
    real(DP), intent(in) :: hcellold !< previous head in cell
    real(DP), intent(inout) :: hcof !< interbed A matrix entry
    real(DP), intent(inout) :: rhs !< interbed right-hand side
    ! -- local variables
    real(DP) :: snnew
    real(DP) :: snold
    real(DP) :: comp
    real(DP) :: compi
    real(DP) :: compe
    real(DP) :: rho1
    real(DP) :: rho2
    real(DP) :: f
    !
    ! -- initialize variables
    rhs = DZERO
    hcof = DZERO
    comp = DZERO
    compi = DZERO
    compe = DZERO
    !
    ! -- skip inactive and constant head cells
    if (this%ibound(node) > 0) then
      if (this%idelay(ib) == 0) then
        !
        ! -- update material properties
        if (this%iupdatematprop /= 0) then
          if (this%ieslag == 0) then
            !
            ! -- calculate compaction
            call this%csub_nodelay_calc_comp(ib, hcell, hcellold, comp, &
                                             rho1, rho2)
            this%comp(ib) = comp
            !
            ! -- update thickness and void ratio
            call this%csub_nodelay_update(ib)
          end if
        end if
        !
        ! -- calculate no-delay interbed rho1 and rho2
        call this%csub_nodelay_fc(ib, hcell, hcellold, rho1, hcof, rhs)
        f = area
      else
        !
        ! -- calculate cell saturation
        call this%csub_calc_sat(node, hcell, hcellold, snnew, snold)
        !
        ! -- update material properties
        if (this%iupdatematprop /= 0) then
          if (this%ieslag == 0) then
            !
            ! -- calculate compaction
            call this%csub_delay_calc_comp(ib, hcell, hcellold, &
                                           comp, compi, compe)
            this%comp(ib) = comp
            !
            ! -- update thickness and void ratio
            call this%csub_delay_update(ib)
          end if
        end if
        !
        ! -- calculate delay interbed hcof and rhs
        call this%csub_delay_sln(ib, hcell)
        call this%csub_delay_fc(ib, hcof, rhs)
        f = area * this%rnb(ib)
      end if
      rhs = rhs * f
      hcof = -hcof * f
    end if
  end subroutine csub_interbed_fc

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_interbed_fn(this, ib, node, hcell, hcellold, hcof, rhs)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy variables
    class(GwfCsubType) :: this
    integer(I4B), intent(in) :: ib !< interbed number
    integer(I4B), intent(in) :: node !< cell node number
    real(DP), intent(in) :: hcell !< current head in a cell
    real(DP), intent(in) :: hcellold !< previous head in a cell
    real(DP), intent(inout) :: hcof !< interbed A matrix entry
    real(DP), intent(inout) :: rhs !< interbed right-hand side entry
    ! -- local variables
    integer(I4B) :: idelay
    real(DP) :: hcofn
    real(DP) :: rhsn
    real(DP) :: top
    real(DP) :: bot
    real(DP) :: tled
    real(DP) :: tthk
    real(DP) :: snnew
    real(DP) :: snold
    real(DP) :: f
    real(DP) :: satderv
    real(DP) :: hbar
    real(DP) :: hbarderv
    real(DP) :: rho1
    real(DP) :: rho2
    !
    ! -- initialize variables
    rhs = DZERO
    hcof = DZERO
    rhsn = DZERO
    hcofn = DZERO
    satderv = DZERO
    idelay = this%idelay(ib)
    top = this%dis%top(node)
    bot = this%dis%bot(node)
    !
    ! -- skip inactive and constant head cells
    if (this%ibound(node) > 0) then
      tled = DONE / delt
      tthk = this%thickini(ib)
      !
      ! -- calculate cell saturation
      call this%csub_calc_sat(node, hcell, hcellold, snnew, snold)
      !
      ! -- no-delay interbeds
      if (idelay == 0) then
        !
        ! -- initialize factor
        f = DONE
        !
        ! -- calculate the saturation derivative
        satderv = this%csub_calc_sat_derivative(node, hcell)
        !
        ! -- calculate corrected head (hbar)
        hbar = sQuadratic0sp(hcell, bot, this%satomega)
        !
        ! -- calculate the derivative of the hbar functions
        hbarderv = sQuadratic0spDerivative(hcell, bot, this%satomega)
        !
        ! -- calculate storage coefficient
        call this%csub_nodelay_fc(ib, hcell, hcellold, rho1, rho2, rhsn)
        !
        ! -- calculate hcofn term
        hcofn = rho2 * (DONE - hbarderv) * snnew + &
                rho2 * (this%cg_gs(node) - hbar + bot) * satderv
        if (this%ielastic(ib) == 0) then
          hcofn = hcofn - rho2 * this%pcs(ib) * satderv
        end if
        !
        ! -- Add additional term if using lagged effective stress
        if (this%ieslag /= 0) then
          if (this%ielastic(ib) /= 0) then
            hcofn = hcofn - rho1 * this%cg_es0(node) * satderv
          else
            hcofn = hcofn - rho1 * (this%pcs(ib) - this%cg_es0(node)) * satderv
          end if
        end if
      end if
    end if
  end subroutine csub_interbed_fn

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_nodelay_wcomp_fc(this, ib, node, tled, area, hcell, &
                                          hcellold, hcof, rhs)
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    integer(I4B), intent(in) :: ib !< interbed number
    integer(I4B), intent(in) :: node !< cell node number
    real(DP), intent(in) :: tled !< reciprocal of time step length
    real(DP), intent(in) :: area !< horizontal cell area
    real(DP), intent(in) :: hcell !< current head in cell
    real(DP), intent(in) :: hcellold !< previous head in cell
    real(DP), intent(inout) :: hcof !< no-delay A matrix entry
    real(DP), intent(inout) :: rhs !< no-delay right-hand side entry
    ! -- local variables
    real(DP) :: top
    real(DP) :: bot
    real(DP) :: snold
    real(DP) :: snnew
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
    !
    ! -- calculate cell saturation
    call this%csub_calc_sat(node, hcell, hcellold, snnew, snold)
    !
    !
    f = this%brg * area * tled
    wc0 = f * this%theta0(ib) * this%thick0(ib)
    wc = f * this%theta(ib) * this%thick(ib)
    hcof = -wc * snnew
    rhs = -wc0 * snold * hcellold
  end subroutine csub_nodelay_wcomp_fc

  ! -- documented on the GwfCsubType type-bound procedure in gwf-csub.f90
  module subroutine csub_nodelay_wcomp_fn(this, ib, node, tled, area, hcell, &
                                          hcellold, hcof, rhs)
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    integer(I4B), intent(in) :: ib !< interbed number
    integer(I4B), intent(in) :: node !< cell node number
    real(DP), intent(in) :: tled !< reciprocal of time step length
    real(DP), intent(in) :: area !< horizontal cell area
    real(DP), intent(in) :: hcell !< current head in cell
    real(DP), intent(in) :: hcellold !< previous head in cell
    real(DP), intent(inout) :: hcof !< no-delay A matrix entry
    real(DP), intent(inout) :: rhs !< no-delay right-hand side entry
    ! -- local variables
    real(DP) :: top
    real(DP) :: bot
    real(DP) :: f
    real(DP) :: wc
    real(DP) :: wc0
    real(DP) :: satderv
    !
    ! -- initialize variables
    rhs = DZERO
    hcof = DZERO
    !
    ! -- aquifer elevations and thickness
    top = this%dis%top(node)
    bot = this%dis%bot(node)
    !
    !
    f = this%brg * area * tled
    !
    ! -- calculate saturation derivative
    satderv = this%csub_calc_sat_derivative(node, hcell)
    !
    ! -- calculate the current water compressibility factor
    wc = f * this%theta(ib) * this%thick(ib)
    !
    ! -- calculate derivative term
    hcof = -wc * hcell * satderv
    !
    ! -- Add additional term if using lagged effective stress
    if (this%ieslag /= 0) then
      wc0 = f * this%theta0(ib) * this%thick0(ib)
      hcof = hcof + wc0 * hcellold * satderv
    end if
    !
    ! -- set rhs
    rhs = hcof * hcell
  end subroutine csub_nodelay_wcomp_fn

end submodule GwfCsubModuleNodelay
