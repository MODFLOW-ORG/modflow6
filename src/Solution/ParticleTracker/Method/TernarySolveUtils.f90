module TernarySolveUtils

  use KindModule, only: I4B, DP, LGP
  use GeomUtilModule, only: skew
  use ConstantsModule, only: DPREC, DZERO, DONE
  use MathUtilModule, only: f1d, zero_ch, zero_br
  use ErrorUtilModule, only: pstop

  private
  public :: traverse_triangle
  public :: canonical
  public :: get_w
  public :: solve_coefs
  public :: step_analytical
  public :: find_exit_bary
  public :: get_t_alpt
  public :: get_bet_outflow_bary
  public :: get_bet_soln_limits
  public :: soln_brent
  public :: soln_chand
  public :: alpfun

  ! global data
  real(DP) ca1, ca2, ca3, cb1, cb2 !< Analytical solution coefficients
  real(DP) waa, wab, wba, wbb !< Elements of the "velocity matrix," W
  real(DP) :: cv0(2), cv1(2), cv2(2) !< "Canonical" velocity components at corners of triangular subcell
  integer(I4B) icase !< Case index for analytical solution

contains

  !> @brief Traverse triangular cell
  subroutine traverse_triangle(isolv, tol, texit, &
                               alpexit, betexit, &
                               itrifaceenter, itrifaceexit, &
                               alp1, bet1, alp2, bet2, alpi, beti)
    ! dummy
    integer(I4B), intent(in) :: isolv !< solution method
    real(DP), intent(in) :: tol !< solution tolerance
    real(DP), intent(out) :: texit !< time particle exits the cell
    real(DP) :: alpexit
    real(DP) :: betexit !< alpha and beta coefficients
    integer(I4B) :: itrifaceenter
    integer(I4B) :: itrifaceexit !< entry and exit faces
    real(DP) :: alp1
    real(DP) :: bet1
    real(DP) :: alp2
    real(DP) :: bet2
    real(DP) :: alpi
    real(DP) :: beti !< alpha and beta coefficients
    ! local
    real(DP) :: texit0
    real(DP) :: alpexit0
    real(DP) :: betexit0
    real(DP) :: texit1
    real(DP) :: alpexit1
    real(DP) :: betexit1
    real(DP) :: texit2
    real(DP) :: alpexit2
    real(DP) :: betexit2

    ! Compute elements of matrix W
    call get_w(alp1, bet1, alp2, bet2, waa, wab, wba, wbb)

    ! Determine alpha and beta analytically as functions of time
    call solve_coefs(alpi, beti)

    ! Compute exit time (travel time to exit) and exit location
    call find_exit_bary(isolv, 0, itrifaceenter, &
                        alpi, beti, tol, &
                        texit0, alpexit0, betexit0)
    call find_exit_bary(isolv, 1, itrifaceenter, &
                        alpi, beti, tol, &
                        texit1, alpexit1, betexit1)
    call find_exit_bary(isolv, 2, itrifaceenter, &
                        alpi, beti, tol, &
                        texit2, alpexit2, betexit2)
    texit = min(texit0, texit1, texit2)

    ! Note that while the numbering of triangle faces is generally zero-based
    ! (0, 1, 2), itrifaceexit, which gets passed out, is one-based (1, 2, 3).
    if (texit .eq. texit0) then
      alpexit = alpexit0
      betexit = betexit0
      itrifaceexit = 1
    else if (texit .eq. texit1) then
      alpexit = alpexit1
      betexit = betexit1
      itrifaceexit = 2
    else if (texit .eq. texit2) then
      alpexit = alpexit2
      betexit = betexit2
      itrifaceexit = 3
    end if
    if (texit .eq. huge(1d0)) itrifaceexit = 0

  end subroutine

  !> @brief Set coordinates to "canonical" configuration
  subroutine canonical(x0, y0, x1, y1, x2, y2, &
                       v0x, v0y, v1x, v1y, v2x, v2y, &
                       xi, yi, &
                       rxx, rxy, ryx, ryy, &
                       sxx, sxy, syy, &
                       alp0, bet0, alp1, bet1, alp2, bet2, alpi, beti)
    ! dummy
    real(DP) :: x0
    real(DP) :: y0
    real(DP) :: x1
    real(DP) :: y1
    real(DP) :: x2
    real(DP) :: y2
    real(DP) :: v0x
    real(DP) :: v0y
    real(DP) :: v1x
    real(DP) :: v1y
    real(DP) :: v2x
    real(DP) :: v2y
    real(DP) :: xi
    real(DP) :: yi
    real(DP) :: rxx
    real(DP) :: rxy
    real(DP) :: ryx
    real(DP) :: ryy !< rotation matrix
    real(DP), intent(inout) :: sxx, sxy, syy !< skew matrix entries (top left, top right, bottom right)
    real(DP) :: alp0
    real(DP) :: bet0
    real(DP) :: alp1
    real(DP) :: bet1
    real(DP) :: alp2
    real(DP) :: bet2
    real(DP) :: alpi
    real(DP) :: beti !< alpha and beta coefficients
    ! local
    real(DP) :: baselen
    real(DP) :: oobaselen
    real(DP) :: sinomega
    real(DP) :: cosomega
    real(DP) :: x1diff
    real(DP) :: y1diff
    real(DP) :: x2diff
    real(DP) :: y2diff
    real(DP) :: xidiff
    real(DP) :: yidiff
    real(DP) :: rot(2, 2), res(2)

    ! Translate and rotate coordinates to "canonical" configuration
    x1diff = x1 - x0
    y1diff = y1 - y0
    x2diff = x2 - x0
    y2diff = y2 - y0
    baselen = dsqrt(x1diff * x1diff + y1diff * y1diff)
    oobaselen = DONE / baselen
    cosomega = x1diff * oobaselen
    sinomega = y1diff * oobaselen
    rxx = cosomega
    rxy = sinomega
    ryx = -sinomega
    ryy = cosomega
    alp0 = DZERO
    bet0 = DZERO
    alp1 = baselen
    bet1 = DZERO

    rot = reshape((/rxx, ryx, rxy, ryy/), shape(rot))
    res = matmul(rot, (/x2diff, y2diff/))
    alp2 = res(1)
    bet2 = res(2)

    cv0 = matmul(rot, (/v0x, v0y/))
    cv1 = matmul(rot, (/v1x, v1y/))
    cv2 = matmul(rot, (/v2x, v2y/))

    xidiff = xi - x0
    yidiff = yi - y0
    res = matmul(rot, (/xidiff, yidiff/))
    alpi = res(1)
    beti = res(2)

    sxx = DONE / alp1
    syy = DONE / bet2
    sxy = -alp2 * sxx * syy
    alp1 = DONE
    alp2 = DZERO
    bet2 = DONE
    cv0 = skew(cv0, (/sxx, sxy, syy/))
    cv1 = skew(cv1, (/sxx, sxy, syy/))
    cv2 = skew(cv2, (/sxx, sxy, syy/))
    res = (/alpi, beti/)
    res = skew(res, (/sxx, sxy, syy/))
    alpi = res(1)
    beti = res(2)

  end subroutine

  !> @brief Compute elements of W matrix
  subroutine get_w( &
    alp1, bet1, alp2, bet2, &
    waa, wab, wba, wbb)
    ! dummy
    real(DP) :: alp1
    real(DP) :: bet1
    real(DP) :: alp2
    real(DP) :: bet2 !< triangle face points
    real(DP) :: waa
    real(DP) :: wab
    real(DP) :: wba
    real(DP) :: wbb !< w matrix
    ! local
    real(DP) :: v1alpdiff
    real(DP) :: v2alpdiff
    real(DP) :: v2betdiff

    ! Note: wab is the "alpha,beta" entry in matrix W
    !    and the alpha component of the w^(beta) vector
    v1alpdiff = cv1(1) - cv0(1)
    v2alpdiff = cv2(1) - cv0(1)
    v2betdiff = cv2(2) - cv0(2)
    waa = v1alpdiff
    wab = v2alpdiff
    wba = DZERO
    wbb = v2betdiff

  end subroutine

  !> @brief Compute analytical solution coefficients depending on case
  subroutine solve_coefs(alpi, beti)
    ! dummy
    real(DP) :: alpi
    real(DP) :: beti
    ! local
    real(DP) :: zerotol
    real(DP) :: wratv
    real(DP) :: acoef
    real(DP) :: bcoef
    real(DP) :: afact
    real(DP) :: bfact
    real(DP) :: vfact
    real(DP) :: oowaa

    zerotol = 1d-10 ! todo AMP: consider tolerance
    if (dabs(wbb) .gt. zerotol) then
      wratv = (wab / wbb) * cv0(2)
      acoef = cv0(1) - wratv
      bcoef = wratv + wab * beti
      afact = acoef / waa
      vfact = cv0(2) / wbb
      ! Coefs for beta do not depend on whether waa = 0 or not
      cb1 = -vfact ! const term in beta
      cb2 = vfact + beti ! coef for e(wbb*t) term in beta
      ! Coefs for alpha
      if (dabs(waa) .gt. zerotol) then
        ! Case waa <> 0, wbb <> 0
        if (dabs(wbb - waa) .gt. zerotol) then
          ! Subcase wbb <> waa
          bfact = bcoef / (wbb - waa)
          ca1 = -afact ! const term in alpha
          ca2 = alpi + afact - bfact ! coef for exp(waa*t) term in alpha
          ca3 = bfact ! coef for exp(wbb*t) term in alpha
          icase = 1
        else
          ! Subcase wbb = waa
          ca1 = -afact ! const term in alpha
          ca2 = alpi + afact ! coef for exp(waa*t) term in alpha
          ca3 = bcoef ! coef for t*exp(waa*t) term in alpha
          icase = -1
        end if
      else
        ! Case waa = 0, wbb <> 0
        bfact = bcoef / wbb
        ca1 = alpi - bfact ! const term in alpha
        ca2 = acoef ! coef for t term in alpha
        ca3 = bfact ! coef for exp(wbb*t) term in alpha
        icase = 2
      end if
    else
      ! Coefs for beta do not depend on whether waa = 0 or not
      cb1 = beti ! const term in beta
      cb2 = cv0(2) ! coef for t term in beta
      if (dabs(waa) .gt. zerotol) then
        ! Case waa <> 0, wbb = 0
        oowaa = 1d0 / waa
        vfact = (wab * oowaa) * cv0(2)
        ca1 = -oowaa * (cv0(1) + wab * beti + vfact) ! const term in alpha
        ca2 = -vfact ! coef for t term in alpha
        ca3 = alpi - ca1 ! coef for exp(waa*t) term in alpha
        icase = 3
      else
        ! Case waa = 0, wbb = 0
        ca1 = alpi ! const term in alpha
        ca2 = cv0(1) + wab * beti ! coef for t term in alpha
        ca3 = 5d-1 * wab * cv0(2) ! coef for t^2 term in alpha
        icase = 4
      end if
    end if

  end subroutine

  !> @brief Step (evaluate) analytically depending on case
  subroutine step_analytical(t, alp, bet)
    ! dummy
    real(DP), intent(in) :: t
    real(DP) :: alp
    real(DP) :: bet

    if (icase .eq. 1) then
      alp = ca1 + ca2 * dexp(waa * t) + ca3 * dexp(wbb * t)
      bet = cb1 + cb2 * dexp(wbb * t)
    else if (icase .eq. -1) then
      alp = ca1 + (ca2 + ca3 * t) * dexp(waa * t)
      bet = cb1 + cb2 * dexp(wbb * t)
    else if (icase .eq. 2) then
      alp = ca1 + ca2 * t + ca3 * dexp(wbb * t)
      bet = cb1 + cb2 * dexp(wbb * t)
    else if (icase .eq. 3) then
      alp = ca1 + ca2 * t + ca3 * dexp(waa * t)
      bet = cb1 + cb2 * t
    else if (icase .eq. 4) then
      alp = ca1 + (ca2 + ca3 * t) * t
      bet = cb1 + cb2 * t
    end if

  end subroutine

  !> @brief Find the exit time and location in barycentric coordinates.
  subroutine find_exit_bary(isolv, itriface, itrifaceenter, &
                            alpi, beti, tol, &
                            texit, alpexit, betexit)
    ! dummy
    integer(I4B) :: isolv
    integer(I4B) :: itriface
    integer(I4B) :: itrifaceenter
    real(DP) :: alpi
    real(DP) :: beti
    real(DP) :: tol
    real(DP) :: texit
    real(DP) :: alpexit
    real(DP) :: betexit
    ! local
    real(DP) :: alplo
    real(DP) :: alphi
    real(DP) :: alpt
    real(DP) :: alplim
    real(DP) :: fax
    real(DP) :: fbx
    real(DP) :: t
    real(DP) :: tlo
    real(DP) :: thi
    real(DP) :: v0alpstar
    real(DP) :: valpi
    real(DP) :: v1n
    real(DP) :: v2n
    real(DP) :: vbeti
    real(DP) :: betlo
    real(DP) :: bethi
    real(DP) :: betsollo
    real(DP) :: betsolhi
    real(DP) :: betoutlo
    real(DP) :: betouthi
    integer(I4B) :: ibettrend
    real(DP) :: zerotol

    zerotol = 1d-10 ! todo AMP: consider tolerance

    ! Use iterative scheme or numerical integration indicated by isolv.
    if (itriface .eq. 0) then
      ! Checking for exit on canonical face 0 (beta = 0)
      if (itrifaceenter .eq. 0) then
        ! Entrance face, so no exit. (Normal velocity is uniform along face 0,
        ! so it cannot be both an entrance and an exit.)
        texit = huge(1d0)
      else
        ! Not the entrance face, so check for outflow
        if (cv0(2) .ge. DZERO) then ! check beta velocity along beta=0 face
          ! Inflow or no flow, so no exit
          texit = huge(DONE)
        else
          ! Outflow, so check beta-velocity at the initial location,
          ! recognizing that it will never change sign along the
          ! trajectory (and will not be blocked from zero by an asymptote)
          vbeti = cv0(2) + wbb * beti
          if (vbeti .ge. DZERO) then
            ! Can't exit along beta = 0
            texit = huge(DONE)
          else
            ! get alpt and check it
            call get_t_alpt(DZERO, t, alpt)
            if ((alpt .ge. DZERO) .and. (alpt .le. DONE)) then
              ! alpt within the edge, so exit found
              texit = t
              alpexit = alpt
              betexit = DZERO
            else
              ! alpt not within the edge, so not an exit
              texit = huge(DONE)
            end if
          end if
        end if
      end if
      ! End canonical face 0 (beta = 0)
    else
      ! Checking for exit on canonical face 1 (gamma = 0.) or 2 (alpha = 0.)
      if (itriface .eq. 1) then
        ! Normal velocities (gamma components) at ends of canonical face 1
        v1n = -cv1(1) - cv1(2)
        v2n = -cv2(1) - cv2(2)
      else
        ! Normal velocities (alpha components) at ends of canonical face 2
        v1n = cv0(1)
        v2n = cv2(1)
      end if
      if ((v1n .ge. DZERO) .and. (v2n .ge. DZERO)) then
        ! No outflow at vn1 and vn2 corners; no outflow interval, so no exit.
        texit = huge(DONE)
      else
        ! Find outflow interval
        call get_bet_outflow_bary(v1n, v2n, betoutlo, betouthi)
        ! Find trend of and limits on beta from beta{t} solution
        call get_bet_soln_limits(beti, betsollo, betsolhi, ibettrend)
        ! Look for exit
        if (ibettrend .eq. 0) then
          ! Beta is constant, so check if it's within the outflow interval;
          ! if not, no exit; if so, solve for t and alpha
          if ((beti .gt. betouthi) .or. (beti .lt. betoutlo)) then
            texit = huge(1d0)
          else
            ! Check alpha-velocity at the initial location,
            ! recognizing that it will never change sign along the
            ! trajectory (and will not be blocked from zero by an asymptote)
            ! in this special case
            v0alpstar = cv0(1) + wab * beti
            valpi = v0alpstar + waa * alpi
            if ((itriface .eq. 1) .and. (valpi .le. DZERO)) then
              ! Can't exit along gamma = 0.
              texit = huge(DONE)
            else if ((itriface .eq. 2) .and. (valpi .ge. DZERO)) then
              ! Can't exit along alpha = 0.
              texit = huge(DONE)
            else
              ! get exit
              if (itriface .eq. 1) then
                alpexit = DONE - beti
              else
                alpexit = DZERO
              end if
              betexit = beti
              if (dabs(waa) > zerotol) then
                alplim = -v0alpstar / waa
                texit = dlog(alpexit - alplim / (alpi - alplim)) / waa
              else
                texit = (alpexit - alpi) / v0alpstar
              end if
            end if
          end if
          ! End constant-beta case
        else
          ! Beta varies along trajectory; combine outflow and soln limits on beta
          bethi = min(betouthi, betsolhi)
          betlo = max(betoutlo, betsollo)
          if (betlo .gt. bethi) then
            ! If bounds on bet leave no feasible interval, no exit
            texit = huge(DONE)
          else
            ! Check sign of function value at beta bounds
            call get_t_alpt(bethi, thi, alphi)
            call get_t_alpt(betlo, tlo, alplo)
            if (itriface .eq. 1) then
              fax = DONE - betlo - alplo
              fbx = DONE - bethi - alphi
            else
              fax = alplo
              fbx = alphi
            end if
            if (fax * fbx .gt. DZERO) then
              ! Root not bracketed; no exit
              texit = huge(DONE)
            else
              if (isolv .eq. 1) then
                ! Use Brent's method with initial bounds on beta of betlo and bethi,
                ! assuming they bound the root
                call soln_brent(itriface, betlo, bethi, tol, texit, &
                                alpexit, betexit)
              else if (isolv .eq. 2) then
                ! Use Chandrupatla's method with initial bounds on beta of betlo and bethi,
                ! assuming they bound the root
                call soln_chand(itriface, betlo, bethi, tol, texit, &
                                alpexit, betexit)
              else
                call pstop(1, "Invalid isolv, expected one of: 1, 2")
              end if
            end if
          end if
          ! End variable-beta case
        end if
      end if
      ! End canonical face 1 (gamma = 0.) or 2 (alpha = 0.)
    end if

    if (texit .ne. huge(DONE) .and. texit .lt. DZERO) then
      call pstop(1, "texit is negative (unexpected)") ! shouldn't get here
    end if

  end subroutine

  !> @brief Brent's method applied to canonical face 1 (gamma = 0)
  function fbary1(bet) result(fb)
    ! dummy
    real(DP), intent(in) :: bet
    real(DP) :: fb
    ! local
    real(DP) :: t
    real(DP) :: alpt

    ! Evaluate gamma{t{beta}} = 1. - alpha{t{beta}} - beta
    call get_t_alpt(bet, t, alpt)
    fb = DONE - alpt - bet
  end function

  !> @brief Brent's method applied to canonical face 2 (alpha = 0)
  function fbary2(bet) result(fb)
    ! dummy
    real(DP), intent(in) :: bet
    real(DP) :: fb
    ! local
    real(DP) :: t
    real(DP) :: alpt

    ! Evaluate alpha{t{beta}}
    call get_t_alpt(bet, t, alpt)
    fb = alpt
  end function

  !> @brief Given beta evaluate t and alpha depending on case
  subroutine get_t_alpt(bet, t, alp)
    ! dummy
    real(DP), intent(in) :: bet
    real(DP) :: t
    real(DP) :: alp
    ! local
    real(DP) :: term
    real(DP) :: zerotol
    real(DP) :: waat
    real(DP) :: coef
    real(DP) :: waatlim
    real(DP) :: ewaatlim

    ! assumes cb2<>0, wbb<>0
    zerotol = 1d-10 ! todo AMP: consider tolerance
    term = (bet - cb1) / cb2
    waatlim = 50.0_DP
    ewaatlim = 5d+21 ! approx. e^waatlim

    if (icase .eq. 1) then
      term = max(term, zerotol)
      t = dlog(term) / wbb
      waat = waa * t
      if (waat > waatlim) then
        alp = sign(ewaatlim, ca2)
      else
        alp = ca1 + ca2 * dexp(waat) + ca3 * term
      end if
    else if (icase .eq. -1) then
      term = max(term, zerotol)
      t = dlog(term) / wbb
      waat = waa * t
      coef = ca2 + ca3 * t
      if (waat > waatlim) then
        alp = sign(ewaatlim, coef)
      else
        alp = ca1 + coef * dexp(waat)
      end if
    else if (icase .eq. 2) then
      term = max(term, zerotol)
      t = dlog(term) / wbb
      alp = ca1 + ca2 * t + ca3 * term
    else if (icase .eq. 3) then
      t = term
      waat = waa * t
      if (waat > waatlim) then
        alp = sign(ewaatlim, ca3)
      else
        alp = ca1 + ca2 * t + ca3 * dexp(waat)
      end if
    else if (icase .eq. 4) then
      t = term
      alp = ca1 + (ca2 + ca3 * t) * t
    end if

  end subroutine

  !> @brief Find outflow interval
  subroutine get_bet_outflow_bary(vn1, vn2, betoutlo, betouthi)
    ! dummy
    real(DP) :: vn1
    real(DP) :: vn2
    real(DP) :: betoutlo
    real(DP) :: betouthi
    ! local
    real(DP) :: vndiff

    vndiff = vn2 - vn1
    if (vn1 .lt. DZERO) then
      ! Outflow at vn1 corner
      betoutlo = DZERO
      if (vn2 .le. DZERO) then
        ! Outflow along entire edge (except possibly no-flow right at vn2 corner)
        betouthi = DONE
      else
        ! Outflow along part of edge
        betouthi = -vn1 / vndiff
      end if
    else
      ! Outflow at vn2 corner
      betouthi = DONE
      if (vn1 .le. DZERO) then
        ! Outflow along entire edge (except possibly no-flow right at vn1 corner)
        betoutlo = DZERO
      else
        ! Outflow along part of edge
        betoutlo = -vn1 / vndiff
      end if
    end if

  end subroutine

  !> @brief Find trend of and limits on beta from beta{t} solution
  subroutine get_bet_soln_limits(beti, betsollo, betsolhi, ibettrend)
    ! dummy
    real(DP), intent(in) :: beti
    real(DP) :: betsollo
    real(DP) :: betsolhi
    integer(I4B), intent(inout) :: ibettrend
    ! local
    real(DP) :: betlim

    if (icase > 2) then
      if (cv0(2) .gt. DZERO) then
        betsolhi = huge(DONE)
        betsollo = beti
        ibettrend = 1
      else if (cv0(2) .lt. DZERO) then
        betsolhi = beti
        betsollo = -huge(DONE)
        ibettrend = -1
      else
        betsolhi = beti
        betsollo = beti
        ibettrend = 0
      end if
    else if (wbb .gt. DZERO) then
      betlim = -cv0(2) / wbb
      if (beti .gt. betlim) then
        betsolhi = huge(DONE)
        betsollo = beti
        ibettrend = 1
      else if (beti .lt. betlim) then
        betsolhi = beti
        betsollo = -huge(DONE)
        ibettrend = -1
      else
        betsolhi = beti
        betsollo = beti
        ibettrend = 0
      end if
    else if (wbb .lt. DZERO) then
      betlim = -cv0(2) / wbb
      if (beti .gt. betlim) then
        betsolhi = beti
        betsollo = betlim
        ibettrend = -1
      else if (beti .lt. betlim) then
        betsolhi = betlim
        betsollo = beti
        ibettrend = 1
      else
        betsolhi = beti
        betsollo = beti
        ibettrend = 0
      end if
    end if

  end subroutine

  !> @brief Use Brent's method with initial bounds on beta of betlo and bethi
  subroutine soln_brent(itriface, betlo, bethi, tol, &
                        texit, alpexit, betexit)
    ! dummy
    integer(I4B), intent(in) :: itriface
    real(DP) :: betlo
    real(DP) :: bethi
    real(DP), intent(in) :: tol
    real(DP) :: texit
    real(DP) :: alpexit
    real(DP) :: betexit
    ! local
    real(DP) :: blo
    real(DP) :: bhi
    procedure(f1d), pointer :: f

    ! assuming betlo and bethi bracket the root
    blo = betlo
    bhi = bethi
    if (itriface .eq. 1) then
      f => fbary1
      betexit = zero_br(blo, bhi, f, tol)
    else
      f => fbary2
      betexit = zero_br(blo, bhi, f, tol)
    end if
    call get_t_alpt(betexit, texit, alpexit)

  end subroutine

  !> @brief Use Chandrupatla's method with initial bounds on beta of betlo and bethi
  subroutine soln_chand(itriface, betlo, bethi, tol, &
                        texit, alpexit, betexit)
    ! dummy
    integer(I4B), intent(in) :: itriface
    real(DP) :: betlo
    real(DP) :: bethi
    real(DP), intent(in) :: tol
    real(DP) :: texit
    real(DP) :: alpexit
    real(DP) :: betexit
    ! local
    real(DP) :: blo
    real(DP) :: bhi
    procedure(f1d), pointer :: f

    ! note: assuming betlo and bethi bracket the root
    blo = betlo
    bhi = bethi
    if (itriface .eq. 1) then
      f => fbary1
      betexit = zero_ch(blo, bhi, f, tol)
    else
      f => fbary2
      betexit = zero_ch(blo, bhi, f, tol)
    end if
    call get_t_alpt(betexit, texit, alpexit)

  end subroutine

end module TernarySolveUtils
