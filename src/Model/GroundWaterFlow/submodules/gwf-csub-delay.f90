!> @brief CSUB package delay interbeds
!<
submodule(GwfCsubModule) GwfCsubModuleDelay
  use LinearAlgebraUtilsModule, only: solve_tridiagonal

contains

  !> @brief Check delay interbed head
  !!
  !! Method to determine if the delay interbed head in any delay cell
  !! in a non-convertible gwf cell is less than the top of each delay
  !! interbed cell.
  !!
  !<
  module procedure csub_delay_head_check
  ! -- local variables
  integer(I4B) :: iviolate
  integer(I4B) :: idelay
  integer(I4B) :: node
  integer(I4B) :: n
  real(DP) :: z
  real(DP) :: h
  real(DP) :: dzhalf
  real(DP) :: ztop
  !
  ! -- initialize variables
  iviolate = 0
  idelay = this%idelay(ib)
  node = this%nodelist(ib)
  !
  ! -- evaluate every delay cell
  idelaycells: do n = 1, this%ndelaycells
    z = this%dbz(n, idelay)
    h = this%dbh(n, idelay)
    dzhalf = DHALF * this%dbdzini(1, idelay)
    !
    ! -- non-convertible cell
    if (this%stoiconv(node) == 0) then
      ztop = z + dzhalf
      if (h < ztop) then
        this%idb_nconv_count(1) = this%idb_nconv_count(1) + 1
        iviolate = 1
      end if
    end if
    !
    ! -- terminate the loop
    if (iviolate > 0) then
      exit idelaycells
    end if
  end do idelaycells
  end procedure csub_delay_head_check

  !> @brief Solve delay interbed continuity equation
  !!
  !! Method to calculate solve the delay interbed continuity equation for a
  !! delay interbed. The method encapsulates the non-linear loop and calls the
  !! linear solution.
  !!
  !<
  module procedure csub_delay_sln
                                                 !! if the maximum head change variable
                                                 !! in a delay bed should be updated
  ! -- local variables
  logical(LGP) :: lupdate
  integer(I4B) :: n
  integer(I4B) :: icnvg
  integer(I4B) :: iter
  integer(I4B) :: idelay
  real(DP) :: dh
  real(DP) :: dhmax
  real(DP) :: dhmax0
  real(DP), parameter :: dclose = DHUNDRED * DPREC
  !
  ! -- initialize variables
  if (present(update)) then
    lupdate = update
  else
    lupdate = .true.
  end if
  !
  ! -- calculate geostatic and effective stress for each delay bed cell
  call this%csub_delay_calc_stress(ib, hcell)
  !
  ! -- terminate if the aquifer head is below the top of delay interbeds
  if (count_errors() > 0) then
    call store_error_filename(this%input_fname)
  end if
  !
  ! -- solve for delay bed heads
  if (this%thickini(ib) > DZERO) then
    icnvg = 0
    iter = 0
    idelay = this%idelay(ib)
    do
      iter = iter + 1
      !
      ! -- assemble coefficients
      call this%csub_delay_assemble(ib, hcell)
      !
      ! -- solve for head change in delay interbed cells
      call solve_tridiagonal(this%ndelaycells, &
                             this%dbal, this%dbad, this%dbau, &
                             this%dbrhs, this%dbdh, this%dbaw)
      !
      ! -- calculate maximum head change and update delay bed heads
      dhmax = DZERO
      do n = 1, this%ndelaycells
        dh = this%dbdh(n) - this%dbh(n, idelay)
        if (abs(dh) > abs(dhmax)) then
          dhmax = dh
          if (lupdate) then
            this%dbdhmax(idelay) = dhmax
          end if
        end if
        ! -- update delay bed heads
        this%dbh(n, idelay) = this%dbdh(n)
      end do
      !
      ! -- update delay bed stresses
      call this%csub_delay_calc_stress(ib, hcell)
      !
      ! -- check delay bed convergence
      if (abs(dhmax) < dclose) then
        icnvg = 1
      else if (iter /= 1) then
        if (abs(dhmax) - abs(dhmax0) < DPREC) then
          icnvg = 1
        end if
      end if
      if (icnvg == 1) then
        exit
      end if
      dhmax0 = dhmax
    end do
  end if
  end procedure csub_delay_sln

  !> @brief Calculate delay interbed znode and z relative to interbed center
  !!
  !! Method to calculate the initial center of each delay interbed cell,
  !! assuming the delay bed head is equal to the top of the delay interbed.
  !! The method also calculates the distance of the center of each delay
  !! bed cell from the center of the delay interbed (z_offset) that is used
  !! to calculate average skeletal specific storage values for a delay interbed
  !! centered on the center of the saturated thickness for a cell.
  !!
  !<
  module procedure csub_delay_init_zcell
  ! -- local variables
  integer(I4B) :: n
  integer(I4B) :: node
  integer(I4B) :: idelay
  real(DP) :: bot
  real(DP) :: top
  real(DP) :: hbar
  real(DP) :: znode
  real(DP) :: dzz
  real(DP) :: z
  real(DP) :: zr
  real(DP) :: b
  real(DP) :: dz
  !
  ! -- initialize variables
  idelay = this%idelay(ib)
  node = this%nodelist(ib)
  b = this%thickini(ib)
  bot = this%dis%bot(node)
  top = bot + b
  hbar = top
  !
  ! -- calculate znode based on assumption that the delay bed bottom
  !    is equal to the cell bottom
  znode = this%csub_calc_znode(top, bot, hbar)
  dz = DHALF * this%dbdzini(1, idelay)
  dzz = DHALF * b
  z = znode + dzz
  zr = dzz
  !
  ! -- calculate z and z relative to znode for each delay
  !    interbed node
  do n = 1, this%ndelaycells
    ! z of node relative to bottom of cell
    z = z - dz
    this%dbz(n, idelay) = z
    z = z - dz
    ! z relative to znode
    zr = zr - dz
    if (ABS(zr) < dz) then
      zr = DZERO
    end if
    this%dbrelz(n, idelay) = zr
    zr = zr - dz
  end do
  end procedure csub_delay_init_zcell

  !> @brief Calculate delay interbed stress values
  !!
  !! Method to calculate the geostatic and effective stress in delay interbed
  !! cells using the passed the current head value in a cell.
  !!
  !<
  module procedure csub_delay_calc_stress
  ! -- local variables
  integer(I4B) :: n
  integer(I4B) :: idelay
  integer(I4B) :: node
  real(DP) :: sigma
  real(DP) :: topaq
  real(DP) :: botaq
  real(DP) :: dzhalf
  real(DP) :: sadd
  real(DP) :: sgm
  real(DP) :: sgs
  real(DP) :: h
  real(DP) :: hbar
  real(DP) :: z
  real(DP) :: top
  real(DP) :: bot
  real(DP) :: phead
  !
  ! -- initialize variables
  idelay = this%idelay(ib)
  node = this%nodelist(ib)
  sigma = this%cg_gs(node)
  topaq = this%dis%top(node)
  botaq = this%dis%bot(node)
  dzhalf = DHALF * this%dbdzini(1, idelay)
  top = this%dbz(1, idelay) + dzhalf
  !
  ! -- calculate corrected head (hbar)
  hbar = sQuadratic0sp(hcell, botaq, this%satomega)
  !
  ! -- calculate the geostatic load in the cell at the top of the interbed.
  sgm = this%sgm(node)
  sgs = this%sgs(node)
  if (hcell < top) then
    sadd = ((top - hbar) * sgm) + ((hbar - botaq) * sgs)
  else
    sadd = (top - botaq) * sgs
  end if
  sigma = sigma - sadd
  !
  ! -- calculate geostatic and effective stress for each interbed node.
  do n = 1, this%ndelaycells
    h = this%dbh(n, idelay)
    !
    ! -- geostatic calculated at the bottom of the delay cell
    z = this%dbz(n, idelay)
    top = z + dzhalf
    bot = z - dzhalf
    !
    ! -- calculate corrected head (hbar)
    hbar = sQuadratic0sp(h, bot, this%satomega)
    !
    ! -- geostatic stress calculation
    if (h < top) then
      sadd = ((top - hbar) * sgm) + ((hbar - bot) * sgs)
    else
      sadd = (top - bot) * sgs
    end if
    sigma = sigma + sadd
    phead = hbar - bot
    this%dbgeo(n, idelay) = sigma
    this%dbes(n, idelay) = sigma - phead
  end do
  end procedure csub_delay_calc_stress

  !> @brief Calculate delay interbed cell storage coefficients
  !<
  module procedure csub_delay_calc_ssksske
  ! -- local variables
  integer(I4B) :: idelay
  integer(I4B) :: ielastic
  integer(I4B) :: node
  real(DP) :: topcell
  real(DP) :: botcell
  real(DP) :: hbarcell
  real(DP) :: zcell
  real(DP) :: zcenter
  real(DP) :: dzhalf
  real(DP) :: top
  real(DP) :: bot
  real(DP) :: h
  real(DP) :: hbar
  real(DP) :: znode
  real(DP) :: zbot
  real(DP) :: es
  real(DP) :: es0
  real(DP) :: theta
  real(DP) :: f
  real(DP) :: f0
  real(DP) :: pcs
  real(DP) :: estop
  real(DP) :: w
  real(DP) :: dwde
  !
  ! -- initialize variables
  sske = DZERO
  ssk = DZERO
  idelay = this%idelay(ib)
  ielastic = this%ielastic(ib)
  !
  ! -- calculate factor for the head-based case
  if (this%lhead_based .EQV. .TRUE.) then
    f = DONE
    f0 = f
    !
    ! -- calculate factor for the effective stress case
  else
    node = this%nodelist(ib)
    theta = this%dbthetaini(n, idelay)
    !
    ! -- set top and bottom of layer
    topcell = this%dis%top(node)
    botcell = this%dis%bot(node)
    !
    ! -- calculate corrected head for the cell (hbarcell)
    hbarcell = sQuadratic0sp(hcell, botcell, this%satomega)
    !
    ! -- set location of delay node relative to the center
    !    of the cell based on current head
    zcell = this%csub_calc_znode(topcell, botcell, hbarcell)
    !
    ! -- set variables for delay interbed zcell calculations
    zcenter = zcell + this%dbrelz(n, idelay)
    dzhalf = DHALF * this%dbdzini(1, idelay)
    top = zcenter + dzhalf
    bot = zcenter - dzhalf
    h = this%dbh(n, idelay)
    !
    ! -- calculate corrected head for the delay interbed cell (hbar)
    hbar = sQuadratic0sp(h, bot, this%satomega)
    !
    ! -- calculate the center of the saturated portion of the
    !    delay interbed cell
    znode = this%csub_calc_znode(top, bot, hbar)
    !
    ! -- set reference point for bottom of delay interbed cell that is used to
    !    scale the effective stress at the bottom of the delay interbed cell
    zbot = this%dbz(n, idelay) - dzhalf
    !
    ! -- set the effective stress
    es = this%dbes(n, idelay)
    es0 = this%dbes0(n, idelay)
    !
    ! -- calculate the compression index factors for the delay
    !    node relative to the center of the cell based on the
    !    current and previous head
    call this%csub_calc_sfacts(node, zbot, znode, theta, es, es0, &
                               this%dbgeo(n, idelay), f)
  end if
  this%idbconvert(n, idelay) = 0
  sske = f * this%rci(ib)
  ssk = f * this%rci(ib)
  if (present(dsskde)) dsskde = DZERO
  ! -- wfac is the inelastic fraction of the storage change used to split the
  !    reported elastic/inelastic budget; 0 while elastic, 1 once inelastic,
  !    and equal to the smoothing weight w across the transition window
  if (present(wfac)) wfac = DZERO
  if (ielastic == 0) then
    es = this%dbes(n, idelay)
    pcs = this%dbpcs(n, idelay)
    ! -- require pcs > DZERO so the smoothing window (pcsomega * pcs) is
    !    positive, as sQuadraticSaturation and its derivative need
    if (this%pcsomega > DZERO .and. pcs > DZERO) then
      ! -- blend elastic (rci) -> inelastic (ci) skeletal storage over a
      !    window of pcsomega * pcs above pcs; w runs 0 (elastic) to 1 (inelastic)
      estop = pcs + this%pcsomega * pcs
      w = sQuadraticSaturation(estop, pcs, es)
      ssk = f * (this%rci(ib) + w * (this%ci(ib) - this%rci(ib)))
      if (w > DHALF) this%idbconvert(n, idelay) = 1
      if (present(wfac)) wfac = w
      if (present(dsskde)) then
        dwde = sQuadraticSaturationDerivative(estop, pcs, es)
        dsskde = f * (this%ci(ib) - this%rci(ib)) * dwde
      end if
    else
      ! -- original hard elastic<->inelastic switch
      if (es > pcs) then
        this%idbconvert(n, idelay) = 1
        ssk = f * this%ci(ib)
        if (present(wfac)) wfac = DONE
      end if
    end if
  end if
  end procedure csub_delay_calc_ssksske

  !> @brief Assemble delay interbed coefficients
  !!
  !! Method to assemble matrix and right-hand side coefficients for a delay
  !! interbed. The method calls the appropriate standard or Newton-Raphson
  !! assembly routines and fills all of the entries for a delay interbed.
  !!
  !<
  module procedure csub_delay_assemble
  ! -- local variables
  integer(I4B) :: n
  real(DP) :: aii
  real(DP) :: au
  real(DP) :: al
  real(DP) :: r
  !
  ! -- calculate matrix terms for each delay bed cell
  do n = 1, this%ndelaycells
    !
    ! -- assemble terms
    if (this%inewton == 0) then
      call this%csub_delay_assemble_fc(ib, n, hcell, aii, au, al, r)
    else
      call this%csub_delay_assemble_fn(ib, n, hcell, aii, au, al, r)
    end if
    !
    ! -- add terms
    this%dbal(n) = al
    this%dbau(n) = au
    this%dbad(n) = aii
    this%dbrhs(n) = r
  end do
  end procedure csub_delay_assemble

  !> @brief Assemble delay interbed standard formulation coefficients
  !!
  !! Method to assemble standard formulation matrix and right-hand side
  !! coefficients for a delay interbed.
  !!
  !<
  module procedure csub_delay_assemble_fc
  ! -- modules
  use TdisModule, only: delt
  ! -- local variables
  integer(I4B) :: node
  integer(I4B) :: idelay
  integer(I4B) :: ielastic
  real(DP) :: dzini
  real(DP) :: dzhalf
  real(DP) :: c
  real(DP) :: c2
  real(DP) :: c3
  real(DP) :: tled
  real(DP) :: wcf
  real(DP) :: smult
  real(DP) :: sske
  real(DP) :: ssk
  real(DP) :: z
  real(DP) :: ztop
  real(DP) :: zbot
  real(DP) :: dz
  real(DP) :: dz0
  real(DP) :: theta
  real(DP) :: theta0
  real(DP) :: dsn
  real(DP) :: dsn0
  real(DP) :: gs
  real(DP) :: es0
  real(DP) :: pcs
  real(DP) :: wc
  real(DP) :: wc0
  real(DP) :: h
  real(DP) :: h0
  real(DP) :: hbar
  !
  ! -- initialize accumulators
  aii = DZERO
  au = DZERO
  al = DZERO
  r = DZERO
  !
  ! -- initialize local variables
  idelay = this%idelay(ib)
  ielastic = this%ielastic(ib)
  node = this%nodelist(ib)
  dzini = this%dbdzini(1, idelay)
  dzhalf = DHALF * dzini
  tled = DONE / delt
  c = this%kv(ib) / dzini
  c2 = DTWO * c
  c3 = DTHREE * c
  !
  ! -- add qdb terms
  aii = aii - c2
  !
  ! -- top or bottom cell
  if (n == 1 .or. n == this%ndelaycells) then
    aii = aii - c
    r = r - c2 * hcell
  end if
  !
  ! -- lower qdb term
  if (n > 1) then
    al = c
  end if
  !
  ! -- upper qdb term
  if (n < this%ndelaycells) then
    au = c
  end if
  !
  ! -- current and previous delay cell states
  z = this%dbz(n, idelay)
  ztop = z + dzhalf
  zbot = z - dzhalf
  h = this%dbh(n, idelay)
  h0 = this%dbh0(n, idelay)
  dz = this%dbdz(n, idelay)
  dz0 = this%dbdz0(n, idelay)
  theta = this%dbtheta(n, idelay)
  theta0 = this%dbtheta0(n, idelay)
  !
  ! -- calculate corrected head (hbar)
  hbar = sQuadratic0sp(h, zbot, this%satomega)
  !
  ! -- calculate saturation
  call this%csub_delay_calc_sat(node, idelay, n, h, h0, dsn, dsn0)
  !
  ! -- calculate ssk and sske
  call this%csub_delay_calc_ssksske(ib, n, hcell, ssk, sske)
  !
  ! -- calculate and add storage terms
  smult = dzini * tled
  gs = this%dbgeo(n, idelay)
  es0 = this%dbes0(n, idelay)
  pcs = this%dbpcs(n, idelay)
  aii = aii - smult * dsn * ssk
  if (ielastic /= 0) then
    r = r - smult * &
        (dsn * ssk * (gs + zbot) - dsn0 * sske * es0)
  else
    r = r - smult * &
        (dsn * ssk * (gs + zbot - pcs) + dsn0 * sske * (pcs - es0))
  end if
  !
  ! -- add storage correction term
  r = r + smult * dsn * ssk * (h - hbar)
  !
  ! -- add water compressibility terms
  wcf = this%brg * tled
  wc = dz * wcf * theta
  wc0 = dz0 * wcf * theta0
  aii = aii - dsn * wc
  r = r - dsn0 * wc0 * h0
  end procedure csub_delay_assemble_fc

  !> @brief Assemble delay interbed Newton-Raphson formulation coefficients
  !!
  !! Method to assemble Newton-Raphson formulation matrix and right-hand side
  !! coefficients for a delay interbed.
  !!
  !<
  module procedure csub_delay_assemble_fn
  ! -- modules
  use TdisModule, only: delt
  ! -- local variables
  integer(I4B) :: node
  integer(I4B) :: idelay
  integer(I4B) :: ielastic
  real(DP) :: dzini
  real(DP) :: dzhalf
  real(DP) :: c
  real(DP) :: c2
  real(DP) :: c3
  real(DP) :: tled
  real(DP) :: wcf
  real(DP) :: smult
  real(DP) :: sske
  real(DP) :: ssk
  real(DP) :: z
  real(DP) :: ztop
  real(DP) :: zbot
  real(DP) :: dz
  real(DP) :: dz0
  real(DP) :: theta
  real(DP) :: theta0
  real(DP) :: dsn
  real(DP) :: dsn0
  real(DP) :: dsnderv
  real(DP) :: wc
  real(DP) :: wc0
  real(DP) :: h
  real(DP) :: h0
  real(DP) :: hbar
  real(DP) :: hbarderv
  real(DP) :: gs
  real(DP) :: es0
  real(DP) :: pcs
  real(DP) :: qsto
  real(DP) :: stoderv
  real(DP) :: dsskde
  real(DP) :: qwc
  real(DP) :: wcderv
  !
  ! -- initialize accumulators
  aii = DZERO
  au = DZERO
  al = DZERO
  r = DZERO
  !
  ! -- initialize local variables
  idelay = this%idelay(ib)
  ielastic = this%ielastic(ib)
  node = this%nodelist(ib)
  dzini = this%dbdzini(1, idelay)
  dzhalf = DHALF * dzini
  tled = DONE / delt
  c = this%kv(ib) / dzini
  c2 = DTWO * c
  c3 = DTHREE * c
  !
  ! -- add qdb terms
  aii = aii - c2
  !
  ! -- top or bottom cell
  if (n == 1 .or. n == this%ndelaycells) then
    aii = aii - c
    r = r - c2 * hcell
  end if
  !
  ! -- lower qdb term
  if (n > 1) then
    al = c
  end if
  !
  ! -- upper qdb term
  if (n < this%ndelaycells) then
    au = c
  end if
  !
  ! -- current and previous delay cell states
  z = this%dbz(n, idelay)
  ztop = z + dzhalf
  zbot = z - dzhalf
  h = this%dbh(n, idelay)
  h0 = this%dbh0(n, idelay)
  dz = this%dbdz(n, idelay)
  dz0 = this%dbdz0(n, idelay)
  theta = this%dbtheta(n, idelay)
  theta0 = this%dbtheta0(n, idelay)
  !
  ! -- calculate corrected head (hbar)
  hbar = sQuadratic0sp(h, zbot, this%satomega)
  !
  ! -- calculate the derivative of the hbar functions
  hbarderv = sQuadratic0spDerivative(h, zbot, this%satomega)
  !
  ! -- calculate saturation
  call this%csub_delay_calc_sat(node, idelay, n, h, h0, dsn, dsn0)
  !
  ! -- calculate the derivative of the saturation
  dsnderv = this%csub_delay_calc_sat_derivative(node, idelay, n, hcell)
  !
  ! -- calculate ssk, sske, and the smoothing derivative dsskde
  call this%csub_delay_calc_ssksske(ib, n, hcell, ssk, sske, dsskde)
  !
  ! -- calculate storage terms
  smult = dzini * tled
  gs = this%dbgeo(n, idelay)
  es0 = this%dbes0(n, idelay)
  pcs = this%dbpcs(n, idelay)
  if (ielastic /= 0) then
    qsto = smult * (dsn * ssk * (gs - hbar + zbot) - dsn0 * sske * es0)
    stoderv = -smult * dsn * ssk * hbarderv + &
              smult * ssk * (gs - hbar + zbot) * dsnderv
  else
    qsto = smult * (dsn * ssk * (gs - hbar + zbot - pcs) + &
                    dsn0 * sske * (pcs - es0))
    stoderv = -smult * dsn * ssk * hbarderv + &
              smult * ssk * (gs - hbar + zbot - pcs) * dsnderv
    ! -- derivative of ssk through the smoothed preconsolidation switch.
    !    es = gs - hbar + zbot  =>  d(es)/d(h) = -hbarderv
    stoderv = stoderv - &
              smult * dsn * dsskde * hbarderv * (gs - hbar + zbot - pcs)
  end if
  !
  ! -- Add additional term if using lagged effective stress
  if (this%ieslag /= 0) then
    if (ielastic /= 0) then
      stoderv = stoderv - smult * sske * es0 * dsnderv
    else
      stoderv = stoderv + smult * sske * (pcs - es0) * dsnderv
    end if
  end if
  !
  ! -- add newton-raphson storage terms
  aii = aii + stoderv
  r = r - qsto + stoderv * h
  !
  ! -- add water compressibility terms
  wcf = this%brg * tled
  wc = dz * wcf * theta
  wc0 = dz0 * wcf * theta0
  qwc = dsn0 * wc0 * h0 - dsn * wc * h
  wcderv = -dsn * wc - wc * h * dsnderv
  !
  ! -- Add additional term if using lagged effective stress
  if (this%ieslag /= 0) then
    wcderv = wcderv + wc0 * h0 * dsnderv
  end if
  !
  ! -- add newton-raphson water compressibility terms
  aii = aii + wcderv
  r = r - qwc + wcderv * h
  end procedure csub_delay_assemble_fn

  !> @brief Calculate delay interbed saturation
  !!
  !! Method to calculate the saturation in a delay interbed cell.
  !<
  module procedure csub_delay_calc_sat
  ! -- local variables
  real(DP) :: dzhalf
  real(DP) :: top
  real(DP) :: bot
  !
  ! -- calculate delay interbed cell saturation
  if (this%stoiconv(node) /= 0) then
    dzhalf = DHALF * this%dbdzini(n, idelay)
    top = this%dbz(n, idelay) + dzhalf
    bot = this%dbz(n, idelay) - dzhalf
    snnew = sQuadraticSaturation(top, bot, hcell, this%satomega)
    snold = sQuadraticSaturation(top, bot, hcellold, this%satomega)
  else
    snnew = DONE
    snold = DONE
  end if
  if (this%ieslag /= 0) then
    snold = snnew
  end if
  end procedure csub_delay_calc_sat

  !> @brief Calculate the delay interbed cell saturation derivative
  !!
  !! Function to calculate the derivative of the saturation with
  !! respect to the current head in delay interbed cell n.
  !!
  !! @return      satderv              derivative of saturation
  !<
  module procedure csub_delay_calc_sat_derivative
  ! -- local variables
  real(DP) :: dzhalf
  real(DP) :: top
  real(DP) :: bot

  if (this%stoiconv(node) /= 0) then
    dzhalf = DHALF * this%dbdzini(n, idelay)
    top = this%dbz(n, idelay) + dzhalf
    bot = this%dbz(n, idelay) - dzhalf
    satderv = sQuadraticSaturationDerivative(top, bot, hcell, this%satomega)
  else
    satderv = DZERO
  end if
  end procedure csub_delay_calc_sat_derivative

  !> @brief Calculate delay interbed storage change
  !!
  !! Method to calculate the storage change in a delay interbed.
  !<
  module procedure csub_delay_calc_dstor
  ! -- local variables
  integer(I4B) :: idelay
  integer(I4B) :: ielastic
  integer(I4B) :: node
  integer(I4B) :: n
  real(DP) :: sske
  real(DP) :: ssk
  real(DP) :: wfac
  real(DP) :: fmult
  real(DP) :: v1
  real(DP) :: v2
  real(DP) :: ske
  real(DP) :: sk
  real(DP) :: z
  real(DP) :: zbot
  real(DP) :: h
  real(DP) :: h0
  real(DP) :: dsn
  real(DP) :: dsn0
  real(DP) :: hbar
  real(DP) :: dzhalf
  !
  ! -- initialize variables
  idelay = this%idelay(ib)
  ielastic = this%ielastic(ib)
  node = this%nodelist(ib)
  stoe = DZERO
  stoi = DZERO
  ske = DZERO
  sk = DZERO
  !
  !
  if (this%thickini(ib) > DZERO) then
    fmult = this%dbdzini(1, idelay)
    dzhalf = DHALF * this%dbdzini(1, idelay)
    do n = 1, this%ndelaycells
      call this%csub_delay_calc_ssksske(ib, n, hcell, ssk, sske, wfac=wfac)
      z = this%dbz(n, idelay)
      zbot = z - dzhalf
      h = this%dbh(n, idelay)
      h0 = this%dbh0(n, idelay)
      call this%csub_delay_calc_sat(node, idelay, n, h, h0, dsn, dsn0)
      hbar = sQuadratic0sp(h, zbot, this%satomega)
      if (ielastic /= 0) then
        v1 = dsn * ssk * (this%dbgeo(n, idelay) - hbar + zbot) - &
             dsn0 * sske * this%dbes0(n, idelay)
        v2 = DZERO
      else
        v1 = dsn * ssk * (this%dbgeo(n, idelay) - hbar + zbot - &
                          this%dbpcs(n, idelay))
        v2 = dsn0 * sske * (this%dbpcs(n, idelay) - this%dbes0(n, idelay))
      end if
      !
      ! -- split the storage change into inelastic and elastic components
      !    weighted by wfac so the reported budget blends across the smoothed
      !    transition; wfac is 0/1 for the hard switch, reproducing the split
      stoi = stoi + wfac * v1 * fmult
      stoe = stoe + ((DONE - wfac) * v1 + v2) * fmult
      !
      ! calculate inelastic and elastic storativity
      ske = ske + sske * fmult
      sk = sk + ssk * fmult
    end do
  end if
  !
  ! -- save ske and sk
  this%ske(ib) = ske
  this%sk(ib) = sk
  end procedure csub_delay_calc_dstor

  !> @brief Calculate delay interbed water compressibility
  !!
  !! Method to calculate the change in water compressibility in a delay interbed.
  !<
  module procedure csub_delay_calc_wcomp
  ! -- modules
  use TdisModule, only: delt
  ! -- local variables
  integer(I4B) :: idelay
  integer(I4B) :: node
  integer(I4B) :: n
  real(DP) :: tled
  real(DP) :: h
  real(DP) :: h0
  real(DP) :: dz
  real(DP) :: dz0
  real(DP) :: dsn
  real(DP) :: dsn0
  real(DP) :: wc
  real(DP) :: wc0
  real(DP) :: v
  !
  ! -- initialize variables
  dwc = DZERO
  !
  !
  if (this%thickini(ib) > DZERO) then
    idelay = this%idelay(ib)
    node = this%nodelist(ib)
    tled = DONE / delt
    do n = 1, this%ndelaycells
      h = this%dbh(n, idelay)
      h0 = this%dbh0(n, idelay)
      dz = this%dbdz(n, idelay)
      dz0 = this%dbdz0(n, idelay)
      call this%csub_delay_calc_sat(node, idelay, n, h, h0, dsn, dsn0)
      wc = dz * this%brg * this%dbtheta(n, idelay)
      wc0 = dz0 * this%brg * this%dbtheta0(n, idelay)
      v = dsn0 * wc0 * h0 - dsn * wc * h
      dwc = dwc + v * tled
    end do
  end if
  end procedure csub_delay_calc_wcomp

  !> @brief Calculate delay interbed compaction
  !!
  !! Method to calculate the compaction in a delay interbed.
  !<
  module procedure csub_delay_calc_comp
  ! -- local variables
  integer(I4B) :: idelay
  integer(I4B) :: ielastic
  integer(I4B) :: node
  integer(I4B) :: n
  real(DP) :: snnew
  real(DP) :: snold
  real(DP) :: sske
  real(DP) :: ssk
  real(DP) :: wfac
  real(DP) :: fmult
  real(DP) :: h
  real(DP) :: h0
  real(DP) :: dsn
  real(DP) :: dsn0
  real(DP) :: v
  real(DP) :: v1
  real(DP) :: v2
  !
  ! -- initialize variables
  idelay = this%idelay(ib)
  ielastic = this%ielastic(ib)
  node = this%nodelist(ib)
  comp = DZERO
  compi = DZERO
  compe = DZERO
  !
  ! -- calculate cell saturation
  call this%csub_calc_sat(node, hcell, hcellold, snnew, snold)
  !
  ! -- calculate compaction
  if (this%thickini(ib) > DZERO) then
    fmult = this%dbdzini(1, idelay)
    do n = 1, this%ndelaycells
      h = this%dbh(n, idelay)
      h0 = this%dbh0(n, idelay)
      call this%csub_delay_calc_sat(node, idelay, n, h, h0, dsn, dsn0)
      call this%csub_delay_calc_ssksske(ib, n, hcell, ssk, sske, wfac=wfac)
      if (ielastic /= 0) then
        v1 = dsn * ssk * this%dbes(n, idelay) - sske * this%dbes0(n, idelay)
        v2 = DZERO
      else
        v1 = dsn * ssk * (this%dbes(n, idelay) - this%dbpcs(n, idelay))
        v2 = dsn0 * sske * (this%dbpcs(n, idelay) - this%dbes0(n, idelay))
      end if
      v = (v1 + v2) * fmult
      comp = comp + v
      !
      ! -- save compaction data
      this%dbcomp(n, idelay) = v * snnew
      !
      ! -- split compaction into inelastic and elastic components weighted by
      !    wfac so the reported budget blends across the smoothed transition;
      !    wfac is 0/1 for the hard switch, reproducing the original split
      compi = compi + wfac * v1 * fmult
      compe = compe + ((DONE - wfac) * v1 + v2) * fmult
    end do
  end if
  !
  ! -- fill compaction
  comp = comp * this%rnb(ib)
  compi = compi * this%rnb(ib)
  compe = compe * this%rnb(ib)
  end procedure csub_delay_calc_comp

  !> @brief Update delay interbed material properties
  !!
  !! Method to update the thickness and porosity of each delay interbed cell.
  !!
  !<
  module procedure csub_delay_update
  ! -- local variables
  integer(I4B) :: idelay
  integer(I4B) :: n
  real(DP) :: comp
  real(DP) :: thick
  real(DP) :: theta
  real(DP) :: tthick
  real(DP) :: wtheta
  !
  ! -- initialize variables
  idelay = this%idelay(ib)
  comp = DZERO
  tthick = DZERO
  wtheta = DZERO
  !
  !
  do n = 1, this%ndelaycells
    !
    ! -- initialize compaction for delay cell
    comp = this%dbtcomp(n, idelay) + this%dbcomp(n, idelay)
    !
    ! -- scale compaction by rnb to get the compaction for
    !    the interbed system (as opposed to the full system)
    comp = comp / this%rnb(ib)
    !
    ! -- update thickness and theta
    if (ABS(comp) > DZERO) then
      thick = this%dbdzini(n, idelay)
      theta = this%dbthetaini(n, idelay)
      call this%csub_adj_matprop(comp, thick, theta)
      if (thick <= DZERO) then
        write (errmsg, '(2(a,i0),a,g0,a)') &
          'Adjusted thickness for delay interbed (', ib, &
          ') cell (', n, ') is less than or equal to 0 (', thick, ').'
        call store_error(errmsg)
      end if
      if (theta <= DZERO) then
        write (errmsg, '(2(a,i0),a,g0,a)') &
          'Adjusted theta for delay interbed (', ib, &
          ') cell (', n, 'is less than or equal to 0 (', theta, ').'
        call store_error(errmsg)
      end if
      this%dbdz(n, idelay) = thick
      this%dbtheta(n, idelay) = theta
      tthick = tthick + thick
      wtheta = wtheta + thick * theta
    else
      thick = this%dbdz(n, idelay)
      theta = this%dbtheta(n, idelay)
      tthick = tthick + thick
      wtheta = wtheta + thick * theta
    end if
  end do
  !
  ! -- calculate thickness weighted theta and save thickness and weighted
  !    theta values for delay interbed
  if (tthick > DZERO) then
    wtheta = wtheta / tthick
  else
    tthick = DZERO
    wtheta = DZERO
  end if
  this%thick(ib) = tthick
  this%theta(ib) = wtheta
  end procedure csub_delay_update

  !> @brief Calculate delay interbed contribution to the cell
  !!
  !! Method to calculate the coefficients to calculate the delay interbed
  !! contribution to a cell. The product of hcof* h - rhs equals the
  !! delay contribution to the cell
  !<
  module procedure csub_delay_fc
  ! -- local variables
  integer(I4B) :: idelay
  real(DP) :: c1
  real(DP) :: c2
  !
  ! -- initialize variables
  idelay = this%idelay(ib)
  hcof = DZERO
  rhs = DZERO
  if (this%thickini(ib) > DZERO) then
    ! -- calculate terms for gwf matrix
    c1 = DTWO * this%kv(ib) / this%dbdzini(1, idelay)
    rhs = -c1 * this%dbh(1, idelay)
    c2 = DTWO * &
         this%kv(ib) / this%dbdzini(this%ndelaycells, idelay)
    rhs = rhs - c2 * this%dbh(this%ndelaycells, idelay)
    hcof = c1 + c2
  end if
  end procedure csub_delay_fc

  !> @brief Calculate the flow from delay interbed top or bottom
  !!
  !! Function to calculate the flow from across the top or bottom of
  !! a delay interbed.
  !!
  !! @return  q  flow across the top or bottom of a delay interbed
  !<
  module procedure csub_calc_delay_flow
  ! -- local variables
  integer(I4B) :: idelay
  real(DP) :: c
  !
  ! -- calculate flow between delay interbed and GWF
  idelay = this%idelay(ib)
  c = DTWO * this%kv(ib) / this%dbdzini(n, idelay)
  q = c * (hcell - this%dbh(n, idelay))
  end procedure csub_calc_delay_flow

end submodule GwfCsubModuleDelay
