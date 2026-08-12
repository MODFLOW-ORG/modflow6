module TestStrandedMass

  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE, DHALF, DTWO
  use MathUtilModule, only: is_close
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       to_string, unittest_type
  use StrandedMassModule, only: strand_rate_sorbed, retained_volume, &
                                return_fraction, decay_amount

  implicit none
  private
  public :: collect_strandedmass

contains

  subroutine collect_strandedmass(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("retained_no_drainage", &
                             test_retained_no_drainage), &
                new_unittest("retained_volume", test_retained_volume), &
                new_unittest("strand_rate_linear", test_strand_rate_linear), &
                new_unittest("return_fraction_saturated", &
                             test_return_fraction_saturated), &
                new_unittest("return_fraction_dry", test_return_fraction_dry), &
                new_unittest("return_fraction_full", &
                             test_return_fraction_full), &
                new_unittest("return_empty_reservoir", &
                             test_return_empty_reservoir), &
                new_unittest("decay_amount", test_decay_amount), &
                new_unittest("strand_return_inverse", &
                             test_strand_return_inverse) &
                ]
  end subroutine collect_strandedmass

  !> @brief No drainage retains nothing
  !<
  subroutine test_retained_no_drainage(error)
    type(error_type), allocatable, intent(out) :: error

    call check(error, retained_volume(DZERO, 100.0_DP, 0.3_DP, DZERO) == DZERO)
    if (allocated(error)) return
    call check(error, strand_rate_sorbed(DZERO, 100.0_DP, 0.9_DP, 1600.0_DP, &
                                         0.2_DP, 1.0_DP) == DZERO)
    if (allocated(error)) return
  end subroutine test_retained_no_drainage

  !> @brief Water retained against drainage is what the flow model does not release
  !<
  subroutine test_retained_volume(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: ds, vcell, thetam, sy, released, thickness

    ds = 0.25_DP
    thickness = 40.0_DP
    vcell = 100.0_DP * thickness
    thetam = 0.3_DP

    ! the flow model releases the drainable part of the pore space, so the
    ! water that stays behind is the rest of it
    sy = 0.15_DP
    released = sy * ds * thickness * (vcell / thickness)
    call check(error, is_close(retained_volume(ds, vcell, thetam, released), &
                               (thetam - sy) * ds * vcell))
    if (allocated(error)) return

    ! a specific yield equal to the porosity retains nothing
    sy = thetam
    released = sy * ds * vcell
    call check(error, retained_volume(ds, vcell, thetam, released) == DZERO)
    if (allocated(error)) return

    ! a specific yield greater than the porosity cannot retain a negative volume
    released = 2.0_DP * thetam * ds * vcell
    call check(error, retained_volume(ds, vcell, thetam, released) == DZERO)
    if (allocated(error)) return
  end subroutine test_retained_volume

  !> @brief Sorbed mass stranded by drainage, for a linear isotherm
  !<
  subroutine test_strand_rate_linear(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: ds, vcell, conc, volfracm, rhob, kd, delt
    real(DP) :: rate, expected

    ds = 0.25_DP
    vcell = 100.0_DP
    conc = 3.0_DP
    volfracm = 0.9_DP
    rhob = 1600.0_DP
    kd = 1.0e-4_DP
    delt = 2.0_DP

    ! for a linear isotherm the sorbed concentration is kd * C
    rate = strand_rate_sorbed(ds, vcell, volfracm, rhob, kd * conc, delt)
    expected = ds * vcell * volfracm * rhob * kd * conc / delt
    call check(error, is_close(rate, expected))
    if (allocated(error)) return
  end subroutine test_strand_rate_linear

  !> @brief A cell that has stranded nothing returns nothing
  !<
  subroutine test_return_fraction_saturated(error)
    type(error_type), allocatable, intent(out) :: error

    call check(error, return_fraction(DZERO, DZERO) == DZERO)
    if (allocated(error)) return
    call check(error, return_fraction(0.1_DP, DZERO) == DZERO)
    if (allocated(error)) return
  end subroutine test_return_fraction_saturated

  !> @brief Rewetting part of what drained returns that share
  !<
  subroutine test_return_fraction_dry(error)
    type(error_type), allocatable, intent(out) :: error

    call check(error, is_close(return_fraction(0.2_DP, 0.5_DP), 0.4_DP))
    if (allocated(error)) return
  end subroutine test_return_fraction_dry

  !> @brief Rewetting all of what drained returns the whole reservoir
  !<
  subroutine test_return_fraction_full(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: held

    held = 0.7_DP
    call check(error, return_fraction(held, held) == DONE)
    if (allocated(error)) return

    ! rewetting more than was drained cannot return more than everything
    call check(error, return_fraction(DONE, held) == DONE)
    if (allocated(error)) return
  end subroutine test_return_fraction_full

  !> @brief An empty reservoir returns zero, never a negative mass
  !<
  subroutine test_return_empty_reservoir(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: mass, ret

    mass = DZERO
    ret = mass * return_fraction(0.2_DP, 0.5_DP)
    call check(error, ret == DZERO)
    if (allocated(error)) return

    ! a falling water table asks for no return at all
    call check(error, return_fraction(-0.2_DP, 0.5_DP) == DZERO)
    if (allocated(error)) return
  end subroutine test_return_empty_reservoir

  !> @brief First-order decay over one step
  !<
  subroutine test_decay_amount(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: mass, lambda, delt

    mass = 10.0_DP
    lambda = 0.1_DP
    delt = 2.0_DP
    call check(error, is_close(decay_amount(mass, lambda, delt), &
                               mass * (DONE - exp(-lambda * delt))))
    if (allocated(error)) return

    ! decay off, or an empty reservoir, loses nothing
    call check(error, decay_amount(mass, DZERO, delt) == DZERO)
    if (allocated(error)) return
    call check(error, decay_amount(DZERO, lambda, delt) == DZERO)
    if (allocated(error)) return

    ! a negative rate is production, so the reservoir gains mass
    call check(error, decay_amount(mass, -lambda, delt) < DZERO)
    if (allocated(error)) return

    ! an empty reservoir cannot produce mass
    call check(error, decay_amount(DZERO, -lambda, delt) == DZERO)
    if (allocated(error)) return
  end subroutine test_decay_amount

  !> @brief Stranding and returning are exact inverses over a full cycle
  !!
  !! This is the algebraic core of conservation. The reservoir is taken to be
  !! spread uniformly through the unsaturated part of the cell, so a cell that
  !! drains from saturated and rewets to saturated recovers all of it. A partial
  !! cycle recovers only the share of the unsaturated part that rewets, because
  !! mass left in the part that stays dry stays there.
  !<
  subroutine test_strand_return_inverse(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: ds, vcell, conc, volfracm, rhob, kd, delt, thetam, released
    real(DP) :: stranded, returned

    vcell = 250.0_DP
    conc = 7.5_DP
    volfracm = 0.85_DP
    rhob = 1750.0_DP
    kd = 5.0e-5_DP
    delt = 3.0_DP
    thetam = 0.3_DP
    ds = 0.4_DP
    released = 0.15_DP * ds * vcell

    ! drain, holding the solute of the retained water and all of the sorbate
    stranded = retained_volume(ds, vcell, thetam, released) * conc + &
               strand_rate_sorbed(ds, vcell, volfracm, rhob, kd * conc, delt) * &
               delt

    ! rewetting all of what drained returns the whole reservoir
    returned = stranded * return_fraction(ds, ds)
    call check(error, is_close(returned, stranded))
    if (allocated(error)) return

    ! a partial rewetting returns the corresponding share of what drained
    returned = stranded * return_fraction(0.3_DP, ds)
    call check(error, is_close(returned, stranded * 0.3_DP / ds))
    if (allocated(error)) return
  end subroutine test_strand_return_inverse

end module TestStrandedMass
