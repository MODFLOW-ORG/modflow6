module TestSmoothing
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DHALF, DONE, DTWO
  use testdrive, only: check, error_type, new_unittest, unittest_type
  use SmoothingModule, only: sSCurve, sCubicLinear, sRamp, sCubic, sLinear, &
                             sQuadratic, sChSmooth, sLinearSaturation, &
                             sCubicSaturation, sQuadraticSaturation, &
                             svanGenuchtenSaturation, &
                             sQuadraticSaturationDerivative, &
                             sQSaturation, sQSaturationDerivative, &
                             sSlope, sSlopeDerivative, sQuadratic0sp, &
                             sQuadratic0spDerivative, sQuadraticSlope, &
                             sQuadraticSlopeDerivative
  implicit none
  private
  public :: collect_smoothing

  real(DP), parameter :: tol = 1.0e-9_DP
  real(DP), parameter :: hfd = 1.0e-6_DP

contains

  subroutine collect_smoothing(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("sRamp_endpoints", test_sramp_endpoints), &
                new_unittest("sRamp_exact_outside_transition", &
                             test_sramp_exact_outside), &
                new_unittest("sRamp_c1_continuity", test_sramp_continuity), &
                new_unittest("sRamp_derivative", test_sramp_derivative), &
                new_unittest("sRamp_monotonic", test_sramp_monotonic), &
                new_unittest("sSCurve", test_sscurve), &
                new_unittest("sCubicLinear", test_scubiclinear), &
                new_unittest("sCubic", test_scubic), &
                new_unittest("sLinear", test_slinear), &
                new_unittest("sQuadratic", test_squadratic), &
                new_unittest("sChSmooth", test_schsmooth), &
                new_unittest("sLinearSaturation", test_slinearsat), &
                new_unittest("sCubicSaturation", test_scubicsat), &
                new_unittest("sQuadraticSaturation", test_squadsat), &
                new_unittest("sQuadraticSaturationDerivative", &
                             test_squadsat_deriv), &
                new_unittest("svanGenuchtenSaturation", test_svg_sat), &
                new_unittest("sQSaturation", test_sqsat), &
                new_unittest("sQSaturationDerivative", test_sqsat_deriv), &
                new_unittest("sQuadratic0sp", test_squad0sp), &
                new_unittest("sQuadratic0spDerivative", test_squad0sp_deriv), &
                new_unittest("sQuadraticSlope", test_squadslope), &
                new_unittest("sQuadraticSlopeDerivative", &
                             test_squadslope_deriv), &
                new_unittest("sSlope", test_sslope), &
                new_unittest("sSlopeDerivative", test_sslope_deriv) &
                ]
  end subroutine collect_smoothing

  !> @brief relative + absolute closeness test for finite-difference checks
  logical function fd_close(d, fd)
    real(DP), intent(in) :: d
    real(DP), intent(in) :: fd
    fd_close = abs(d - fd) <= 1.0e-6_DP + 1.0e-4_DP * abs(d)
  end function fd_close

  ! -------------------------------------------------------------------------
  ! sRamp: smoothed ramp, value ~ max(x, 0)
  ! -------------------------------------------------------------------------

  subroutine test_sramp_endpoints(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: w, y, dydx
    w = 0.5_DP
    call sRamp(-1.0_DP, w, dydx, y)
    call check(error, y == DZERO .and. dydx == DZERO, "x<0")
    if (allocated(error)) return
    call sRamp(DZERO, w, dydx, y)
    call check(error, y == DZERO .and. dydx == DZERO, "x=0")
    if (allocated(error)) return
    call sRamp(2.0_DP, w, dydx, y)
    call check(error, abs(y - 2.0_DP) < tol .and. abs(dydx - DONE) < tol, "x>w")
    if (allocated(error)) return
    call sRamp(w, w, dydx, y)
    call check(error, abs(y - w) < tol .and. abs(dydx - DONE) < tol, "x=w")
  end subroutine test_sramp_endpoints

  subroutine test_sramp_exact_outside(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: w, x, y, dydx
    integer(I4B) :: i
    w = 0.25_DP
    do i = -5, 5
      x = real(i, DP)
      if (abs(x) < w) cycle
      call sRamp(x, w, dydx, y)
      call check(error, abs(y - max(x, DZERO)) < tol, "max(x,0)")
      if (allocated(error)) return
    end do
  end subroutine test_sramp_exact_outside

  subroutine test_sramp_continuity(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: w, ylo, yhi, dlo, dhi
    w = 0.5_DP
    call sRamp(-hfd, w, dlo, ylo)
    call sRamp(hfd, w, dhi, yhi)
    call check(error, abs(yhi - ylo) < 1.0e-5_DP, "value at 0")
    if (allocated(error)) return
    call check(error, abs(dhi - dlo) < 1.0e-5_DP, "slope at 0")
    if (allocated(error)) return
    call sRamp(w - hfd, w, dlo, ylo)
    call sRamp(w + hfd, w, dhi, yhi)
    call check(error, abs(yhi - ylo) < 1.0e-5_DP, "value at w")
    if (allocated(error)) return
    call check(error, abs(dhi - dlo) < 1.0e-5_DP, "slope at w")
  end subroutine test_sramp_continuity

  subroutine test_sramp_derivative(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: w, x, y, dydx, yp, ym, dummy, fd
    integer(I4B) :: i
    w = 1.0_DP
    do i = 1, 9
      x = 0.1_DP * real(i, DP)
      call sRamp(x, w, dydx, y)
      call sRamp(x + hfd, w, dummy, yp)
      call sRamp(x - hfd, w, dummy, ym)
      fd = (yp - ym) / (DTWO * hfd)
      call check(error, fd_close(dydx, fd), "dydx vs fd")
      if (allocated(error)) return
    end do
  end subroutine test_sramp_derivative

  subroutine test_sramp_monotonic(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: w, x, y, dydx, yprev
    integer(I4B) :: i
    w = 0.5_DP
    yprev = -DONE
    do i = -20, 20
      x = 0.1_DP * real(i, DP)
      call sRamp(x, w, dydx, y)
      call check(error, y >= yprev - tol, "non-decreasing")
      if (allocated(error)) return
      call check(error, dydx >= -tol, "slope >= 0")
      if (allocated(error)) return
      yprev = y
    end do
  end subroutine test_sramp_monotonic

  ! -------------------------------------------------------------------------
  ! 0 -> 1 smoothing subroutines, signature (x, range, dydx, y)
  ! -------------------------------------------------------------------------

  subroutine test_sscurve(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: w, x, y, dydx, yp, ym, dd
    integer(I4B) :: i
    w = 1.0_DP
    ! endpoints
    call sSCurve(-1.0_DP, w, dydx, y)
    call check(error, y == DZERO .and. dydx == DZERO, "x<0")
    if (allocated(error)) return
    call sSCurve(2.0_DP, w, dydx, y)
    call check(error, abs(y - DONE) < tol .and. dydx == DZERO, "x>range")
    if (allocated(error)) return
    ! bounds and derivative
    do i = 1, 9
      x = 0.1_DP * real(i, DP)
      call sSCurve(x, w, dydx, y)
      call check(error, y >= -tol .and. y <= DONE + tol, "in [0,1]")
      if (allocated(error)) return
      call sSCurve(x + hfd, w, dd, yp)
      call sSCurve(x - hfd, w, dd, ym)
      call check(error, fd_close(dydx, (yp - ym) / (DTWO * hfd)), "dydx")
      if (allocated(error)) return
    end do
  end subroutine test_sscurve

  subroutine test_scubiclinear(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: w, x, y, dydx, yp, ym, dd
    integer(I4B) :: i
    w = 1.0_DP
    call sCubicLinear(-1.0_DP, w, dydx, y)
    call check(error, y == DZERO .and. dydx == DZERO, "x<0")
    if (allocated(error)) return
    call sCubicLinear(2.0_DP, w, dydx, y)
    call check(error, abs(y - DONE) < tol .and. dydx == DZERO, "x>range")
    if (allocated(error)) return
    do i = 1, 9
      x = 0.1_DP * real(i, DP)
      call sCubicLinear(x, w, dydx, y)
      call check(error, y >= -tol .and. y <= DONE + tol, "in [0,1]")
      if (allocated(error)) return
      call sCubicLinear(x + hfd, w, dd, yp)
      call sCubicLinear(x - hfd, w, dd, ym)
      call check(error, fd_close(dydx, (yp - ym) / (DTWO * hfd)), "dydx")
      if (allocated(error)) return
    end do
  end subroutine test_scubiclinear

  subroutine test_scubic(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: x, w, y, dydx, xp, xm, rr, yp, ym, dd
    integer(I4B) :: i
    ! x >= range -> 1, slope 0
    x = 2.0_DP
    w = 1.0_DP
    call sCubic(x, w, dydx, y)
    call check(error, abs(y - DONE) < tol .and. dydx == DZERO, "x>range")
    if (allocated(error)) return
    ! interior: in [0,1] and slope matches finite difference
    do i = 2, 8
      x = 0.1_DP * real(i, DP)
      w = 1.0_DP
      call sCubic(x, w, dydx, y)
      call check(error, y >= -tol .and. y <= DONE + tol, "in [0,1]")
      if (allocated(error)) return
      xp = 0.1_DP * real(i, DP) + hfd
      rr = 1.0_DP
      call sCubic(xp, rr, dd, yp)
      xm = 0.1_DP * real(i, DP) - hfd
      rr = 1.0_DP
      call sCubic(xm, rr, dd, ym)
      call check(error, fd_close(dydx, (yp - ym) / (DTWO * hfd)), "dydx")
      if (allocated(error)) return
    end do
  end subroutine test_scubic

  subroutine test_slinear(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: x, w, y, dydx, xp, xm, rr, yp, ym, dd
    integer(I4B) :: i
    x = 2.0_DP
    w = 1.0_DP
    call sLinear(x, w, dydx, y)
    call check(error, abs(y - DONE) < tol .and. dydx == DZERO, "x>range")
    if (allocated(error)) return
    do i = 2, 8
      x = 0.1_DP * real(i, DP)
      w = 1.0_DP
      call sLinear(x, w, dydx, y)
      call check(error, y >= -tol .and. y <= DONE + tol, "in [0,1]")
      if (allocated(error)) return
      xp = 0.1_DP * real(i, DP) + hfd
      rr = 1.0_DP
      call sLinear(xp, rr, dd, yp)
      xm = 0.1_DP * real(i, DP) - hfd
      rr = 1.0_DP
      call sLinear(xm, rr, dd, ym)
      call check(error, fd_close(dydx, (yp - ym) / (DTWO * hfd)), "dydx")
      if (allocated(error)) return
    end do
  end subroutine test_slinear

  subroutine test_squadratic(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: x, w, y, dydx, xp, xm, rr, yp, ym, dd
    integer(I4B) :: i
    x = 2.0_DP
    w = 1.0_DP
    call sQuadratic(x, w, dydx, y)
    call check(error, abs(y - DONE) < tol .and. dydx == DZERO, "x>range")
    if (allocated(error)) return
    do i = 2, 8
      x = 0.1_DP * real(i, DP)
      w = 1.0_DP
      call sQuadratic(x, w, dydx, y)
      call check(error, y >= -tol .and. y <= DONE + tol, "in [0,1]")
      if (allocated(error)) return
      xp = 0.1_DP * real(i, DP) + hfd
      rr = 1.0_DP
      call sQuadratic(xp, rr, dd, yp)
      xm = 0.1_DP * real(i, DP) - hfd
      rr = 1.0_DP
      call sQuadratic(xm, rr, dd, ym)
      call check(error, fd_close(dydx, (yp - ym) / (DTWO * hfd)), "dydx")
      if (allocated(error)) return
    end do
  end subroutine test_squadratic

  ! -------------------------------------------------------------------------
  ! sChSmooth: smooths channel variables over a tiny range near zero depth
  ! -------------------------------------------------------------------------

  subroutine test_schsmooth(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: d, smooth, dwdh, sp, sm, dd, h, fd
    integer(I4B) :: i
    h = 1.0e-9_DP
    ! d <= 0 -> 0, d >= 1e-5 -> 1
    call sChSmooth(-1.0_DP, smooth, dwdh)
    call check(error, smooth == DZERO .and. dwdh == DZERO, "d<0")
    if (allocated(error)) return
    call sChSmooth(1.0_DP, smooth, dwdh)
    call check(error, abs(smooth - DONE) < tol .and. dwdh == DZERO, "d>range")
    if (allocated(error)) return
    ! in the smoothing band: bounded and slope matches finite difference
    do i = 1, 9
      d = 1.0e-6_DP * real(i, DP)
      call sChSmooth(d, smooth, dwdh)
      call check(error, smooth >= -tol .and. smooth <= DONE + tol, "in [0,1]")
      if (allocated(error)) return
      call sChSmooth(d + h, sp, dd)
      call sChSmooth(d - h, sm, dd)
      fd = (sp - sm) / (DTWO * h)
      call check(error, fd_close(dwdh, fd), "dwdh vs fd")
      if (allocated(error)) return
    end do
  end subroutine test_schsmooth

  ! -------------------------------------------------------------------------
  ! saturation functions, signature (top, bot, x[, ...]) -> [0, 1]
  ! -------------------------------------------------------------------------

  subroutine test_slinearsat(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: top, bot, x, y
    integer(I4B) :: i
    top = 10.0_DP
    bot = 2.0_DP
    call check(error, sLinearSaturation(top, bot, 1.0_DP) == DZERO, "below bot")
    if (allocated(error)) return
    call check(error, abs(sLinearSaturation(top, bot, 11.0_DP) - DONE) < tol, &
               "above top")
    if (allocated(error)) return
    ! midpoint is 0.5, and values stay in [0, 1]
    call check(error, abs(sLinearSaturation(top, bot, 6.0_DP) - DHALF) < tol, &
               "midpoint")
    if (allocated(error)) return
    do i = 2, 10
      x = real(i, DP)
      y = sLinearSaturation(top, bot, x)
      call check(error, y >= -tol .and. y <= DONE + tol, "in [0,1]")
      if (allocated(error)) return
    end do
  end subroutine test_slinearsat

  subroutine test_scubicsat(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: top, bot, x, y
    integer(I4B) :: i
    top = 10.0_DP
    bot = 0.0_DP
    call check(error, sCubicSaturation(top, bot, -1.0_DP) == DZERO, "below bot")
    if (allocated(error)) return
    call check(error, abs(sCubicSaturation(top, bot, 11.0_DP) - DONE) < tol, &
               "above top")
    if (allocated(error)) return
    do i = 0, 10
      x = real(i, DP)
      y = sCubicSaturation(top, bot, x)
      call check(error, y >= -tol .and. y <= DONE + tol, "in [0,1]")
      if (allocated(error)) return
    end do
  end subroutine test_scubicsat

  subroutine test_squadsat(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: top, bot, x, y
    integer(I4B) :: i
    top = 10.0_DP
    bot = 0.0_DP
    call check(error, abs(sQuadraticSaturation(top, bot, -1.0_DP)) < tol, &
               "below bot")
    if (allocated(error)) return
    call check(error, abs(sQuadraticSaturation(top, bot, 11.0_DP) - DONE) < tol, &
               "above top")
    if (allocated(error)) return
    do i = 0, 10
      x = real(i, DP)
      y = sQuadraticSaturation(top, bot, x)
      call check(error, y >= -tol .and. y <= DONE + tol, "in [0,1]")
      if (allocated(error)) return
    end do
  end subroutine test_squadsat

  subroutine test_squadsat_deriv(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: top, bot, eps, x, d, fd, yp, ym
    integer(I4B) :: i
    real(DP), dimension(3) :: xs
    top = 10.0_DP
    bot = 0.0_DP
    eps = 0.1_DP
    xs = [0.5_DP, 5.0_DP, 9.5_DP]
    do i = 1, size(xs)
      x = xs(i)
      d = sQuadraticSaturationDerivative(top, bot, x, eps)
      yp = sQuadraticSaturation(top, bot, x + hfd, eps)
      ym = sQuadraticSaturation(top, bot, x - hfd, eps)
      fd = (yp - ym) / (DTWO * hfd)
      call check(error, fd_close(d, fd), "deriv vs fd")
      if (allocated(error)) return
    end do
  end subroutine test_squadsat_deriv

  subroutine test_svg_sat(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: top, bot, alpha, beta, sr, x, y
    integer(I4B) :: i
    top = 1.0_DP
    bot = 0.0_DP
    alpha = 1.0_DP
    beta = 2.0_DP
    sr = 0.1_DP
    ! at zero/negative capillary pressure (x >= b/2) the function returns 0
    call check(error, svanGenuchtenSaturation(top, bot, 0.6_DP, alpha, beta, &
                                              sr) == DZERO, "pc<=0")
    if (allocated(error)) return
    ! elsewhere it stays within [0, 1]
    do i = 0, 4
      x = 0.1_DP * real(i, DP)
      y = svanGenuchtenSaturation(top, bot, x, alpha, beta, sr)
      call check(error, y >= -tol .and. y <= DONE + tol, "in [0,1]")
      if (allocated(error)) return
    end do
  end subroutine test_svg_sat

  ! -------------------------------------------------------------------------
  ! sQSaturation pair, signature (top, bot, x[, c1, c2])
  ! -------------------------------------------------------------------------

  subroutine test_sqsat(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: top, bot, x, y
    integer(I4B) :: i
    top = 10.0_DP
    bot = 0.0_DP
    call check(error, sQSaturation(top, bot, -1.0_DP) == DZERO, "below bot")
    if (allocated(error)) return
    call check(error, abs(sQSaturation(top, bot, 11.0_DP) - DONE) < tol, &
               "above top")
    if (allocated(error)) return
    do i = 0, 10
      x = real(i, DP)
      y = sQSaturation(top, bot, x)
      call check(error, y >= -tol .and. y <= DONE + tol, "in [0,1]")
      if (allocated(error)) return
    end do
  end subroutine test_sqsat

  subroutine test_sqsat_deriv(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: top, bot, x, d, fd, yp, ym
    integer(I4B) :: i
    top = 10.0_DP
    bot = 0.0_DP
    do i = 1, 9
      x = real(i, DP)
      d = sQSaturationDerivative(top, bot, x)
      yp = sQSaturation(top, bot, x + hfd)
      ym = sQSaturation(top, bot, x - hfd)
      fd = (yp - ym) / (DTWO * hfd)
      call check(error, fd_close(d, fd), "deriv vs fd")
      if (allocated(error)) return
    end do
  end subroutine test_sqsat_deriv

  ! -------------------------------------------------------------------------
  ! sQuadratic0sp pair, signature (x, xi[, tomega]): smooths max(xi, x)
  ! -------------------------------------------------------------------------

  subroutine test_squad0sp(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: xi, omega, x
    integer(I4B) :: i
    xi = 1.0_DP
    omega = 0.5_DP
    ! well left of xi -> xi; well right of xi -> x
    call check(error, abs(sQuadratic0sp(0.0_DP, xi, omega) - xi) < tol, "left")
    if (allocated(error)) return
    call check(error, abs(sQuadratic0sp(2.0_DP, xi, omega) - 2.0_DP) < tol, &
               "right")
    if (allocated(error)) return
    ! always >= xi (smoothed maximum)
    do i = 0, 20
      x = 0.1_DP * real(i, DP)
      call check(error, sQuadratic0sp(x, xi, omega) >= xi - tol, ">= xi")
      if (allocated(error)) return
    end do
  end subroutine test_squad0sp

  subroutine test_squad0sp_deriv(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: xi, omega, x, d, fd, yp, ym
    integer(I4B) :: i
    xi = 1.0_DP
    omega = 0.5_DP
    do i = -4, 4
      x = xi + 0.05_DP * real(i, DP)
      d = sQuadratic0spDerivative(x, xi, omega)
      yp = sQuadratic0sp(x + hfd, xi, omega)
      ym = sQuadratic0sp(x - hfd, xi, omega)
      fd = (yp - ym) / (DTWO * hfd)
      call check(error, fd_close(d, fd), "deriv vs fd")
      if (allocated(error)) return
    end do
  end subroutine test_squad0sp_deriv

  ! -------------------------------------------------------------------------
  ! sQuadraticSlope pair: smooths a change in slope from sm to sp at xi
  ! -------------------------------------------------------------------------

  subroutine test_squadslope(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: xi, yi, sm, sp, omega, eps, ylo, yhi
    xi = 2.0_DP
    yi = 1.0_DP
    sm = 0.5_DP
    sp = -1.5_DP
    omega = 0.5_DP
    eps = DHALF * omega
    ! left/right limbs are the lines yi + slope*dx
    call check(error, abs(sQuadraticSlope(0.0_DP, xi, yi, sm, sp, omega) - &
                          (yi + sm * (0.0_DP - xi))) < tol, "left limb")
    if (allocated(error)) return
    call check(error, abs(sQuadraticSlope(4.0_DP, xi, yi, sm, sp, omega) - &
                          (yi + sp * (4.0_DP - xi))) < tol, "right limb")
    if (allocated(error)) return
    ! value is continuous where the quadratic meets each limb (xi +/- eps)
    ylo = sQuadraticSlope(xi - eps - hfd, xi, yi, sm, sp, omega)
    yhi = sQuadraticSlope(xi - eps + hfd, xi, yi, sm, sp, omega)
    call check(error, abs(yhi - ylo) < 1.0e-5_DP, "continuous at xi-eps")
    if (allocated(error)) return
    ylo = sQuadraticSlope(xi + eps - hfd, xi, yi, sm, sp, omega)
    yhi = sQuadraticSlope(xi + eps + hfd, xi, yi, sm, sp, omega)
    call check(error, abs(yhi - ylo) < 1.0e-5_DP, "continuous at xi+eps")
  end subroutine test_squadslope

  subroutine test_squadslope_deriv(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: xi, yi, sm, sp, omega, x, d, fd, yp, ym
    integer(I4B) :: i
    xi = 2.0_DP
    yi = 1.0_DP
    sm = 0.5_DP
    sp = -1.5_DP
    omega = 0.5_DP
    do i = -4, 4
      x = xi + 0.05_DP * real(i, DP)
      d = sQuadraticSlopeDerivative(x, xi, sm, sp, omega)
      yp = sQuadraticSlope(x + hfd, xi, yi, sm, sp, omega)
      ym = sQuadraticSlope(x - hfd, xi, yi, sm, sp, omega)
      fd = (yp - ym) / (DTWO * hfd)
      call check(error, fd_close(d, fd), "deriv vs fd")
      if (allocated(error)) return
    end do
  end subroutine test_squadslope_deriv

  ! -------------------------------------------------------------------------
  ! sSlope pair: smooths a change in slope from sm to sp at xi
  ! -------------------------------------------------------------------------

  subroutine test_sslope(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: xi, yi, sm, sp, taa, x, d, fd, yp, ym
    integer(I4B) :: i
    xi = 2.0_DP
    yi = 1.0_DP
    sm = 0.5_DP
    sp = -1.5_DP
    ! limbs (default tiny smoothing) are the lines yi + slope*dx
    call check(error, abs(sSlope(0.0_DP, xi, yi, sm, sp) - &
                          (yi + sm * (0.0_DP - xi))) < 1.0e-6_DP, "left limb")
    if (allocated(error)) return
    call check(error, abs(sSlope(4.0_DP, xi, yi, sm, sp) - &
                          (yi + sp * (4.0_DP - xi))) < 1.0e-6_DP, "right limb")
    if (allocated(error)) return
    call check(error, abs(sSlope(xi, xi, yi, sm, sp) - yi) < 1.0e-6_DP, "at xi")
    if (allocated(error)) return
    ! sSlopeDerivative is the derivative of sSlope across the smoothing window
    taa = 0.2_DP
    do i = -5, 5
      x = xi + 0.1_DP * real(i, DP)
      d = sSlopeDerivative(x, xi, sm, sp, taa)
      yp = sSlope(x + hfd, xi, yi, sm, sp, taa)
      ym = sSlope(x - hfd, xi, yi, sm, sp, taa)
      fd = (yp - ym) / (DTWO * hfd)
      call check(error, fd_close(d, fd), "deriv vs fd")
      if (allocated(error)) return
    end do
  end subroutine test_sslope

  subroutine test_sslope_deriv(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: xi, sm, sp, d
    xi = 0.0_DP
    sm = 2.0_DP
    sp = -1.0_DP
    ! far left -> sm
    d = sSlopeDerivative(-1.0_DP, xi, sm, sp)
    call check(error, abs(d - sm) < 1.0e-6_DP, "left limit sm")
    if (allocated(error)) return
    ! far right -> sp
    d = sSlopeDerivative(1.0_DP, xi, sm, sp)
    call check(error, abs(d - sp) < 1.0e-6_DP, "right limit sp")
    if (allocated(error)) return
    ! at xi -> average of sm and sp
    d = sSlopeDerivative(xi, xi, sm, sp)
    call check(error, abs(d - DHALF * (sm + sp)) < tol, "midpoint avg")
  end subroutine test_sslope_deriv

end module TestSmoothing
