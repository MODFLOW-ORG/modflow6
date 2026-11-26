module TestCompactionUtil
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: DZERO, DONE
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       to_string, unittest_type
  use MathUtilModule, only: is_close
  use CompactionUtilModule, only: porosity_to_void_ratio, &
                                  void_ratio_to_porosity, compact_layer, &
                                  adjust_eff_stress, storage_factor, &
                                  storage_newton_term, wcomp_newton_term
  use SmoothingModule, only: sQuadraticSaturation, &
                             sQuadraticSaturationDerivative, sQuadratic0sp, &
                             sQuadratic0spDerivative
  implicit none
  private
  public :: collect_compactionutil

  ! -- cell and smoothing used by the Newton-Raphson tests
  real(DP), parameter :: TOP = 10.0_DP
  real(DP), parameter :: BOT = DZERO
  real(DP), parameter :: OMEGA = 0.1_DP
  real(DP), parameter :: DH = 1.0e-4_DP
  ! -- heads that are dry, in the smoothed intervals, partially saturated,
  !    and above the top
  real(DP), dimension(6), parameter :: HEADS = &
                                       [-0.5_DP, 0.02_DP, 0.5_DP, 5.0_DP, &
                                        9.5_DP, 12.0_DP]

contains

  subroutine collect_compactionutil(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("void_ratio_porosity", test_void_ratio_porosity), &
                new_unittest("compact_layer_solids", test_compact_layer_solids), &
                new_unittest("compact_layer_zero_thickness", &
                             test_compact_layer_zero_thickness), &
                new_unittest("adjust_eff_stress", test_adjust_eff_stress), &
                new_unittest("storage_factor", test_storage_factor), &
                new_unittest("storage_newton_term", test_storage_newton_term), &
                new_unittest("wcomp_newton_term", test_wcomp_newton_term) &
                ]
  end subroutine collect_compactionutil

  !> @brief Check known values and the round trip
  !<
  subroutine test_void_ratio_porosity(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(4), parameter :: thetas = &
                                         [0.05_DP, 0.25_DP, 0.4_DP, 0.6_DP]
    real(DP) :: theta
    integer(I4B) :: i

    call check(error, is_close(porosity_to_void_ratio(0.25_DP), &
                               DONE / 3.0_DP))
    if (allocated(error)) return
    call check(error, is_close(void_ratio_to_porosity(DONE), 0.5_DP))
    if (allocated(error)) return

    do i = 1, size(thetas)
      theta = void_ratio_to_porosity(porosity_to_void_ratio(thetas(i)))
      call check(error, is_close(theta, thetas(i), rtol=1.0e-14_DP), &
                 "theta = "//to_string(theta))
      if (allocated(error)) return
    end do
  end subroutine test_void_ratio_porosity

  !> @brief Solid thickness is unchanged by compaction and expansion
  !<
  subroutine test_compact_layer_solids(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(3), parameter :: comps = [0.5_DP, -0.2_DP, DZERO]
    real(DP) :: thick
    real(DP) :: theta
    real(DP) :: solids0
    real(DP) :: solids1
    integer(I4B) :: i

    do i = 1, size(comps)
      thick = 10.0_DP
      theta = 0.4_DP
      solids0 = thick / (DONE + porosity_to_void_ratio(theta))
      call compact_layer(comps(i), thick, theta)
      solids1 = thick / (DONE + porosity_to_void_ratio(theta))
      call check(error, is_close(thick, 10.0_DP - comps(i)), &
                 "thick = "//to_string(thick))
      if (allocated(error)) return
      call check(error, is_close(solids1, solids0, rtol=1.0e-12_DP), &
                 "solids "//to_string(solids0)//" -> "//to_string(solids1))
      if (allocated(error)) return
    end do

    ! compaction lowers the porosity
    thick = 10.0_DP
    theta = 0.4_DP
    call compact_layer(0.5_DP, thick, theta)
    call check(error, theta < 0.4_DP, "theta = "//to_string(theta))
  end subroutine test_compact_layer_solids

  !> @brief A layer with no thickness keeps its porosity
  !<
  subroutine test_compact_layer_zero_thickness(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: thick
    real(DP) :: theta

    thick = DZERO
    theta = 0.3_DP
    call compact_layer(0.1_DP, thick, theta)
    call check(error, is_close(theta, 0.3_DP, rtol=1.0e-14_DP), &
               "theta = "//to_string(theta))
    if (allocated(error)) return
    call check(error, is_close(thick, -0.1_DP), "thick = "//to_string(thick))
  end subroutine test_compact_layer_zero_thickness

  !> @brief Effective stress changes linearly with elevation
  !<
  subroutine test_adjust_eff_stress(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: es

    es = adjust_eff_stress(100.0_DP, 2.0_DP, 5.0_DP, 1.8_DP)
    call check(error, is_close(es, 97.6_DP), "es = "//to_string(es))
    if (allocated(error)) return

    es = adjust_eff_stress(100.0_DP, 2.0_DP, 2.0_DP, 1.8_DP)
    call check(error, is_close(es, 100.0_DP), "es = "//to_string(es))
  end subroutine test_adjust_eff_stress

  !> @brief Check a known value and zero effective stress
  !<
  subroutine test_storage_factor(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: fact

    fact = storage_factor(50.0_DP, 0.2_DP)
    call check(error, is_close(fact, 0.016_DP), "fact = "//to_string(fact))
    if (allocated(error)) return

    fact = storage_factor(DZERO, 0.2_DP)
    call check(error, fact == DZERO, "fact = "//to_string(fact))
  end subroutine test_storage_factor

  !> @brief Newton storage term matches a finite-difference derivative
  !!
  !! The storage flow is formed as hcof h - rhs from the Picard terms
  !! (csub_cg_fc and csub_nodelay_fc) for elastic and inelastic storage,
  !! with and without the lagged saturation.
  !<
  subroutine test_storage_newton_term(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), parameter :: rho = 2.0_DP
    real(DP), parameter :: rho0 = 1.5_DP
    real(DP), parameter :: gs = 30.0_DP
    real(DP), parameter :: es0 = 22.0_DP
    real(DP), parameter :: h0 = 6.0_DP
    real(DP), dimension(2), parameter :: pcss = [DZERO, 25.0_DP]
    real(DP) :: h
    real(DP) :: sn
    real(DP) :: es
    real(DP) :: dqdh
    real(DP) :: hcof
    real(DP) :: expected
    logical(LGP) :: lagged
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: k

    do k = 1, 2
      lagged = (k == 2)
      do j = 1, size(pcss)
        do i = 1, size(HEADS)
          h = HEADS(i)
          sn = sQuadraticSaturation(TOP, BOT, h, OMEGA)
          es = gs + BOT - sQuadratic0sp(h, BOT, OMEGA)
          dqdh = (storage_flow(h + DH) - storage_flow(h - DH)) / (2.0_DP * DH)
          expected = dqdh + rho * sn
          hcof = storage_newton_term(rho, rho0, sn, &
                                     sQuadraticSaturationDerivative(TOP, BOT, &
                                                                    h, OMEGA), &
                                     sQuadratic0spDerivative(h, BOT, OMEGA), &
                                     es, es0, pcss(j), lagged)
          call check(error, is_close(hcof, expected, rtol=1.0e-6_DP, &
                                     atol=1.0e-8_DP), &
                     "h = "//to_string(h)//" pcs = "//to_string(pcss(j))// &
                     " hcof = "//to_string(hcof)//" expected = "// &
                     to_string(expected))
          if (allocated(error)) return
        end do
      end do
    end do

  contains

    !> @brief Storage flow, hcof h - rhs, from the Picard terms
    !<
    function storage_flow(hh) result(q)
      real(DP), intent(in) :: hh
      real(DP) :: q
      real(DP) :: snh
      real(DP) :: sn0
      real(DP) :: hbar
      real(DP) :: rhs

      snh = sQuadraticSaturation(TOP, BOT, hh, OMEGA)
      if (lagged) then
        sn0 = snh
      else
        sn0 = sQuadraticSaturation(TOP, BOT, h0, OMEGA)
      end if
      hbar = sQuadratic0sp(hh, BOT, OMEGA)
      rhs = -rho * snh * (gs + BOT) + pcss(j) * (rho * snh - rho0 * sn0) + &
            rho0 * sn0 * es0 - rho * snh * (hh - hbar)
      q = -rho * snh * hh - rhs
    end function storage_flow

  end subroutine test_storage_newton_term

  !> @brief Newton water compressibility term matches a finite difference
  !<
  subroutine test_wcomp_newton_term(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), parameter :: wc = 0.3_DP
    real(DP), parameter :: wc0 = 0.25_DP
    real(DP), parameter :: h0 = 6.0_DP
    real(DP) :: h
    real(DP) :: sn
    real(DP) :: dqdh
    real(DP) :: hcof
    real(DP) :: expected
    logical(LGP) :: lagged
    integer(I4B) :: i
    integer(I4B) :: k

    do k = 1, 2
      lagged = (k == 2)
      do i = 1, size(HEADS)
        h = HEADS(i)
        sn = sQuadraticSaturation(TOP, BOT, h, OMEGA)
        dqdh = (wcomp_flow(h + DH) - wcomp_flow(h - DH)) / (2.0_DP * DH)
        expected = dqdh + wc * sn
        hcof = wcomp_newton_term(wc, wc0, h, h0, &
                                 sQuadraticSaturationDerivative(TOP, BOT, h, &
                                                                OMEGA), &
                                 lagged)
        call check(error, is_close(hcof, expected, rtol=1.0e-6_DP, &
                                   atol=1.0e-8_DP), &
                   "h = "//to_string(h)//" hcof = "//to_string(hcof)// &
                   " expected = "//to_string(expected))
        if (allocated(error)) return
      end do
    end do

  contains

    !> @brief Water compressibility flow, hcof h - rhs, from the Picard terms
    !<
    function wcomp_flow(hh) result(q)
      real(DP), intent(in) :: hh
      real(DP) :: q
      real(DP) :: snh
      real(DP) :: sn0

      snh = sQuadraticSaturation(TOP, BOT, hh, OMEGA)
      if (lagged) then
        sn0 = snh
      else
        sn0 = sQuadraticSaturation(TOP, BOT, h0, OMEGA)
      end if
      q = -wc * snh * hh + wc0 * sn0 * h0
    end function wcomp_flow

  end subroutine test_wcomp_newton_term

end module TestCompactionUtil
