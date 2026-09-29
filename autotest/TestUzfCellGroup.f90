module TestUzfCellGroup
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DHALF, DONE, DTWO
  use MemoryHelperModule, only: create_mem_path
  use UzfCellGroupModule, only: UzfCellGroupType
  use testdrive, only: error_type, unittest_type, new_unittest, check

  implicit none
  private
  public :: collect_uzfcellgroup

  integer(I4B), parameter :: NCELLS = 3
  integer(I4B), parameter :: NWAV = 20
  integer(I4B), parameter :: ICELL = 2

  ! -- a three wave train for ICELL, deepest wave first
  real(DP), parameter :: THETA_RES = 0.05_DP
  real(DP), parameter :: THETA_SAT = 0.35_DP
  real(DP), dimension(3), parameter :: DEPTH = [10.0_DP, 6.0_DP, 2.0_DP]
  real(DP), dimension(3), parameter :: THETA = [0.20_DP, 0.15_DP, 0.10_DP]

contains

  subroutine collect_uzfcellgroup(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("store_load_roundtrip", test_store_load_roundtrip), &
                new_unittest("store_load_leaves_other_cells", &
                             test_store_load_other_cells), &
                new_unittest("shift_waves_up", test_shift_waves_up), &
                new_unittest("shift_waves_down", test_shift_waves_down), &
                new_unittest("unsat_stor", test_unsat_stor), &
                new_unittest("unsat_stor_partial", test_unsat_stor_partial), &
                new_unittest("get_wcnew", test_get_wcnew), &
                new_unittest("water_content_at_depth", test_wc_at_depth), &
                new_unittest("caph", test_caph) &
                ]
  end subroutine collect_uzfcellgroup

  !> @brief Build a cell group with a known three wave train in ICELL
  !<
  subroutine setup(uzf, tag)
    type(UzfCellGroupType), intent(inout) :: uzf
    character(len=*), intent(in) :: tag
    integer(I4B) :: j
    integer(I4B) :: icell

    call uzf%init(NCELLS, NWAV, create_mem_path('TESTUZF', tag))
    do icell = 1, NCELLS
      uzf%theta_res(icell) = THETA_RES
      uzf%theta_sat(icell) = THETA_SAT
      uzf%bc_eps(icell) = 4.0_DP
      uzf%vks(icell) = 1.0_DP
      uzf%air_entry(icell) = 0.5_DP
      uzf%celtop(icell) = 100.0_DP
      uzf%celbot(icell) = 80.0_DP
      uzf%water_table(icell) = 90.0_DP
      uzf%nwaves(icell) = 3
      do j = 1, 3
        uzf%wave_depth(j, icell) = DEPTH(j)
        uzf%wave_theta(j, icell) = THETA(j)
        uzf%wave_flux(j, icell) = DZERO
        uzf%wave_speed(j, icell) = DZERO
      end do
    end do
  end subroutine setup

  !> @brief Storing then restoring a wave train reproduces it exactly
  !<
  subroutine test_store_load_roundtrip(error)
    type(error_type), allocatable, intent(out) :: error
    type(UzfCellGroupType) :: uzf
    integer(I4B) :: j

    call setup(uzf, 'ROUNDTRIP')
    call uzf%store_waves(ICELL, uzf%wavsav)

    ! -- overwrite the train, as a trial solution would
    uzf%nwaves(ICELL) = 1
    do j = 1, 3
      uzf%wave_depth(j, ICELL) = -DONE
      uzf%wave_theta(j, ICELL) = -DONE
      uzf%wave_flux(j, ICELL) = -DONE
      uzf%wave_speed(j, ICELL) = -DONE
    end do

    call uzf%load_waves(ICELL, uzf%wavsav)

    call check(error, uzf%nwaves(ICELL), 3)
    if (allocated(error)) return
    do j = 1, 3
      call check(error, uzf%wave_depth(j, ICELL), DEPTH(j))
      if (allocated(error)) return
      call check(error, uzf%wave_theta(j, ICELL), THETA(j))
      if (allocated(error)) return
    end do
  end subroutine test_store_load_roundtrip

  !> @brief Storing and restoring one cell does not touch its neighbors
  !<
  subroutine test_store_load_other_cells(error)
    type(error_type), allocatable, intent(out) :: error
    type(UzfCellGroupType) :: uzf
    integer(I4B) :: j

    call setup(uzf, 'OTHERCELLS')
    uzf%wave_theta(1, 1) = 0.25_DP
    uzf%wave_theta(1, 3) = 0.27_DP

    call uzf%store_waves(ICELL, uzf%wavsav)
    do j = 1, 3
      uzf%wave_theta(j, ICELL) = DZERO
    end do
    call uzf%load_waves(ICELL, uzf%wavsav)

    call check(error, uzf%wave_theta(1, 1), 0.25_DP)
    if (allocated(error)) return
    call check(error, uzf%wave_theta(1, 3), 0.27_DP)
  end subroutine test_store_load_other_cells

  !> @brief Shifting by -1 opens position 1 for a new wave at the surface
  !<
  subroutine test_shift_waves_up(error)
    type(error_type), allocatable, intent(out) :: error
    type(UzfCellGroupType) :: uzf

    call setup(uzf, 'SHIFTUP')
    ! -- move waves 1..3 up to 2..4, counting down so the source stays ahead
    call uzf%shift_waves(ICELL, -1, 4, 2, -1)

    call check(error, uzf%wave_depth(2, ICELL), DEPTH(1))
    if (allocated(error)) return
    call check(error, uzf%wave_depth(3, ICELL), DEPTH(2))
    if (allocated(error)) return
    call check(error, uzf%wave_depth(4, ICELL), DEPTH(3))
    if (allocated(error)) return
    call check(error, uzf%wave_theta(2, ICELL), THETA(1))
  end subroutine test_shift_waves_up

  !> @brief Shifting by +1 drops the deepest wave off the bottom
  !<
  subroutine test_shift_waves_down(error)
    type(error_type), allocatable, intent(out) :: error
    type(UzfCellGroupType) :: uzf

    call setup(uzf, 'SHIFTDOWN')
    ! -- move waves 2..3 down to 1..2, counting up so the source stays ahead
    call uzf%shift_waves(ICELL, 1, 1, 2, 1)

    call check(error, uzf%wave_depth(1, ICELL), DEPTH(2))
    if (allocated(error)) return
    call check(error, uzf%wave_depth(2, ICELL), DEPTH(3))
    if (allocated(error)) return
    call check(error, uzf%wave_theta(1, ICELL), THETA(2))
  end subroutine test_shift_waves_down

  !> @brief Mobile water over the full train is the sum of the wave segments
  !<
  subroutine test_unsat_stor(error)
    type(error_type), allocatable, intent(out) :: error
    type(UzfCellGroupType) :: uzf
    real(DP) :: d
    real(DP) :: expected

    call setup(uzf, 'UNSATSTOR')
    ! -- surface to 2 at THETA(3), 2 to 6 at THETA(2), 6 to 10 at THETA(1)
    expected = (THETA(3) - THETA_RES) * (DEPTH(3) - DZERO) + &
               (THETA(2) - THETA_RES) * (DEPTH(2) - DEPTH(3)) + &
               (THETA(1) - THETA_RES) * (DEPTH(1) - DEPTH(2))
    d = DEPTH(1)

    call check(error, uzf%unsat_stor(ICELL, d), expected, thr=1.0e-12_DP)
  end subroutine test_unsat_stor

  !> @brief Mobile water above a depth inside the train stops at that depth
  !<
  subroutine test_unsat_stor_partial(error)
    type(error_type), allocatable, intent(out) :: error
    type(UzfCellGroupType) :: uzf
    real(DP) :: d
    real(DP) :: expected

    call setup(uzf, 'UNSATPART')
    ! -- depth 4 lies in the segment carrying THETA(2)
    expected = (THETA(3) - THETA_RES) * (DEPTH(3) - DZERO) + &
               (THETA(2) - THETA_RES) * (4.0_DP - DEPTH(3))
    d = 4.0_DP

    call check(error, uzf%unsat_stor(ICELL, d), expected, thr=1.0e-12_DP)
  end subroutine test_unsat_stor_partial

  !> @brief Cell water content is the depth-averaged wave content
  !<
  subroutine test_get_wcnew(error)
    type(error_type), allocatable, intent(out) :: error
    type(UzfCellGroupType) :: uzf
    real(DP) :: d
    real(DP) :: thk
    real(DP) :: expected

    call setup(uzf, 'GETWCNEW')
    thk = uzf%celtop(ICELL) - uzf%water_table(ICELL)
    d = thk
    expected = uzf%unsat_stor(ICELL, d) / thk + THETA_RES

    call check(error, uzf%get_wcnew(ICELL), expected, thr=1.0e-12_DP)
  end subroutine test_get_wcnew

  !> @brief Water content at a depth is the content of the wave spanning it
  !<
  subroutine test_wc_at_depth(error)
    type(error_type), allocatable, intent(out) :: error
    type(UzfCellGroupType) :: uzf

    call setup(uzf, 'WCATDEPTH')
    ! -- depth 4 lies in the segment carrying THETA(2)
    call check(error, uzf%get_water_content_at_depth(ICELL, 4.0_DP), THETA(2), &
               thr=1.0e-9_DP)
    if (allocated(error)) return
    ! -- below the water table the cell is saturated
    call check(error, uzf%get_water_content_at_depth(ICELL, 15.0_DP), THETA_SAT)
  end subroutine test_wc_at_depth

  !> @brief Capillary head follows the Brooks-Corey retention curve
  !<
  subroutine test_caph(error)
    type(error_type), allocatable, intent(out) :: error
    type(UzfCellGroupType) :: uzf
    real(DP) :: star
    real(DP) :: lambda
    real(DP) :: expected
    real(DP) :: tho

    call setup(uzf, 'CAPH')
    tho = 0.20_DP
    star = (tho - THETA_RES) / (THETA_SAT - THETA_RES)
    lambda = DTWO / (uzf%bc_eps(ICELL) - 3.0_DP)
    expected = uzf%air_entry(ICELL) * star**(-DONE / lambda)

    call check(error, uzf%caph(ICELL, tho), expected, thr=1.0e-12_DP)
    if (allocated(error)) return
    ! -- at saturation the head is the air entry potential
    call check(error, uzf%caph(ICELL, THETA_SAT), uzf%air_entry(ICELL), &
               thr=1.0e-12_DP)
    if (allocated(error)) return
    ! -- above saturation there is no capillary head
    call check(error, uzf%caph(ICELL, THETA_SAT + DONE), DZERO)
  end subroutine test_caph

end module TestUzfCellGroup
