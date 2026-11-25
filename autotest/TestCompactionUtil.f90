module TestCompactionUtil
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       to_string, unittest_type
  use MathUtilModule, only: is_close
  use CompactionUtilModule, only: porosity_to_void_ratio, &
                                  void_ratio_to_porosity, compact_layer, &
                                  adjust_eff_stress, storage_factor
  implicit none
  private
  public :: collect_compactionutil

contains

  subroutine collect_compactionutil(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("void_ratio_porosity", test_void_ratio_porosity), &
                new_unittest("compact_layer_solids", test_compact_layer_solids), &
                new_unittest("compact_layer_zero_thickness", &
                             test_compact_layer_zero_thickness), &
                new_unittest("adjust_eff_stress", test_adjust_eff_stress), &
                new_unittest("storage_factor", test_storage_factor) &
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

end module TestCompactionUtil
