module TestUltimateScheme
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DHALF, DONE
  use testdrive, only: check, error_type, new_unittest, unittest_type
  use UltimateSchemeModule, only: universal_limiter, multidim_bounds
  implicit none
  private
  public :: collect_ultimatescheme

  real(DP), parameter :: tol = 1.0e-12_DP
  ! sample values for the parameter sweeps
  real(DP), dimension(5), parameter :: phis = &
                                       [0.0_DP, 0.2_DP, 0.5_DP, 0.8_DP, 1.0_DP]
  real(DP), dimension(4), parameter :: courants = &
                                       [0.1_DP, 0.25_DP, 0.5_DP, 1.0_DP]

contains

  subroutine collect_ultimatescheme(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("universal_limiter_courant_one", &
                             test_universal_courant_one), &
                new_unittest("universal_limiter_extremum", &
                             test_universal_extremum), &
                new_unittest("universal_limiter_bounds", &
                             test_universal_bounds), &
                new_unittest("universal_limiter_mirror", &
                             test_universal_mirror), &
                new_unittest("multidim_bounds_contain_upwind", &
                             test_multidim_contain_upwind), &
                new_unittest("multidim_bounds_match_1d", &
                             test_multidim_match_1d), &
                new_unittest("multidim_bounds_no_outflow", &
                             test_multidim_no_outflow), &
                new_unittest("multidim_bounds_keep_cell_bounded", &
                             test_multidim_keep_bounded) &
                ]
  end subroutine collect_ultimatescheme

  !> @brief At a Courant number of one a monotone profile gives the upwind value
  subroutine test_universal_courant_one(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: i
    real(DP) :: phi

    do i = 1, size(phis)
      ! decreasing and increasing profiles
      phi = universal_limiter(phis(i), 1.0_DP, 0.6_DP, 0.2_DP, DONE)
      call check(error, abs(phi - 0.6_DP) < tol, "decreasing profile")
      if (allocated(error)) return
      phi = universal_limiter(phis(i), 0.2_DP, 0.6_DP, 1.0_DP, DONE)
      call check(error, abs(phi - 0.6_DP) < tol, "increasing profile")
      if (allocated(error)) return
    end do
  end subroutine test_universal_courant_one

  !> @brief A local extremum upstream reverts to the upwind value
  subroutine test_universal_extremum(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: i, j
    real(DP) :: phi

    do i = 1, size(phis)
      do j = 1, size(courants)
        phi = universal_limiter(phis(i), 0.2_DP, 0.8_DP, 0.3_DP, courants(j))
        call check(error, abs(phi - 0.8_DP) < tol, "maximum")
        if (allocated(error)) return
        phi = universal_limiter(phis(i), 0.7_DP, 0.1_DP, 0.4_DP, courants(j))
        call check(error, abs(phi - 0.1_DP) < tol, "minimum")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_universal_extremum

  !> @brief A decreasing profile clamps the face value to its monotone range
  subroutine test_universal_bounds(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: phi

    ! with a Courant number of one half the lower bound is
    ! max(0.2, 1 + (0.6 - 1) / 0.5) = 0.2 and the upper bound is 0.6
    phi = universal_limiter(0.0_DP, 1.0_DP, 0.6_DP, 0.2_DP, DHALF)
    call check(error, abs(phi - 0.2_DP) < tol, "below the range")
    if (allocated(error)) return
    phi = universal_limiter(0.9_DP, 1.0_DP, 0.6_DP, 0.2_DP, DHALF)
    call check(error, abs(phi - 0.6_DP) < tol, "above the range")
    if (allocated(error)) return
    phi = universal_limiter(0.4_DP, 1.0_DP, 0.6_DP, 0.2_DP, DHALF)
    call check(error, abs(phi - 0.4_DP) < tol, "within the range")
    if (allocated(error)) return
    ! a steeper upstream gradient raises the lower bound to
    ! 1 + (0.8 - 1) / 0.25 = 0.2 against a downstream value of 0
    phi = universal_limiter(0.0_DP, 1.0_DP, 0.8_DP, 0.0_DP, 0.25_DP)
    call check(error, abs(phi - 0.2_DP) < tol, "courant bound")
  end subroutine test_universal_bounds

  !> @brief Reflecting every value about one half reflects the result
  subroutine test_universal_mirror(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: i, j, k, l, ic
    real(DP) :: phi, phimirror

    do i = 1, size(phis)
      do j = 1, size(phis)
        do k = 1, size(phis)
          do l = 1, size(phis)
            do ic = 1, size(courants)
              phi = universal_limiter(phis(i), phis(j), phis(k), phis(l), &
                                      courants(ic))
              phimirror = universal_limiter(DONE - phis(i), DONE - phis(j), &
                                            DONE - phis(k), DONE - phis(l), &
                                            courants(ic))
              call check(error, abs(phi + phimirror - DONE) < tol, "mirror")
              if (allocated(error)) return
            end do
          end do
        end do
      end do
    end do
  end subroutine test_universal_mirror

  !> @brief The bounds always admit the upwind value
  subroutine test_multidim_contain_upwind(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: i, j, k, ic, jc
    real(DP) :: phimin, phimax, lower, upper

    do i = 1, size(phis)
      do j = 1, size(phis)
        do k = 1, size(phis)
          phimin = min(phis(i), phis(j))
          phimax = max(phis(i), phis(j))
          do ic = 1, size(courants)
            do jc = 1, size(courants)
              ! unequal inflow and outflow as well as equal
              call multidim_bounds(phis(k), min(phimin, phis(k)), &
                                   max(phimax, phis(k)), phimin, phimax, &
                                   courants(ic), courants(jc), lower, upper)
              call check(error, lower <= phis(k) + tol .and. &
                         upper >= phis(k) - tol, "upwind value excluded")
              if (allocated(error)) return
            end do
          end do
        end do
      end do
    end do
  end subroutine test_multidim_contain_upwind

  !> @brief With one inflow and one outflow face the bounds match the 1D limiter
  subroutine test_multidim_match_1d(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: ic
    real(DP) :: c, lower, upper, phi

    ! a decreasing profile 1, 0.8, 0 steep enough that the Courant number sets
    ! the lower bound: the upstream range is [0.8, 1] and the downstream range
    ! [0, 0.8]
    do ic = 1, size(courants)
      c = courants(ic)
      call multidim_bounds(0.8_DP, 0.8_DP, 1.0_DP, 0.0_DP, 0.8_DP, c, c, &
                           lower, upper)
      phi = universal_limiter(DZERO, 1.0_DP, 0.8_DP, 0.0_DP, c)
      call check(error, abs(lower - phi) < tol, "lower bound")
      if (allocated(error)) return
      phi = universal_limiter(DONE, 1.0_DP, 0.8_DP, 0.0_DP, c)
      call check(error, abs(upper - phi) < tol, "upper bound")
      if (allocated(error)) return
    end do
  end subroutine test_multidim_match_1d

  !> @brief Without outflow only the downstream range applies
  subroutine test_multidim_no_outflow(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: lower, upper

    call multidim_bounds(0.5_DP, 0.2_DP, 0.9_DP, 0.3_DP, 0.7_DP, DHALF, DZERO, &
                         lower, upper)
    call check(error, abs(lower - 0.3_DP) < tol, "lower bound")
    if (allocated(error)) return
    call check(error, abs(upper - 0.7_DP) < tol, "upper bound")
  end subroutine test_multidim_no_outflow

  !> @brief Faces at the bounds keep the explicit update within the cell range
  !!
  !! For equal inflow and outflow, the update is
  !! phi + c * (mean inflow face value) - c * (mean outflow face value); the
  !! worst cases pair inflow at one end of the range with outflow at the
  !! opposite bound.
  !<
  subroutine test_multidim_keep_bounded(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: i, j, k, ic
    real(DP) :: phi, phimin, phimax, c, lower, upper, phinew

    do i = 1, size(phis)
      do j = 1, size(phis)
        do k = 1, size(phis)
          phi = phis(k)
          phimin = min(phis(i), phis(j), phi)
          phimax = max(phis(i), phis(j), phi)
          do ic = 1, size(courants)
            c = courants(ic)
            ! a downstream range wide enough not to bind
            call multidim_bounds(phi, phimin, phimax, phimin - DONE, &
                                 phimax + DONE, c, c, lower, upper)
            phinew = phi + c * phimin - c * upper
            call check(error, phinew >= phimin - tol, "undershoot")
            if (allocated(error)) return
            phinew = phi + c * phimax - c * lower
            call check(error, phinew <= phimax + tol, "overshoot")
            if (allocated(error)) return
          end do
        end do
      end do
    end do
  end subroutine test_multidim_keep_bounded

end module TestUltimateScheme
