module TestStrandedMass

  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE, DHALF, DTWO, DNODATA
  use MathUtilModule, only: is_close
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       to_string, unittest_type
  use StrandedMassModule, only: StrandedMassType, strand_rate_sorbed, &
                                retained_volume, return_fraction, decay_amount, &
                                flush_fraction, deposit_fraction

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
                new_unittest("return_fraction_small", &
                             test_return_fraction_small), &
                new_unittest("return_empty_reservoir", &
                             test_return_empty_reservoir), &
                new_unittest("decay_amount", test_decay_amount), &
                new_unittest("strand_return_inverse", &
                             test_strand_return_inverse), &
                new_unittest("advance_restore", test_advance_restore), &
                new_unittest("flush_fraction", test_flush_fraction) &
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

  !> @brief Mass held from a very small drainage returns in full
  !!
  !! Rewetting by as much as the cell drained returns all of it, however
  !! small the drained fraction, so no mass is left that cannot return.
  !<
  subroutine test_return_fraction_small(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: held

    held = 1.0e-9_DP
    call check(error, return_fraction(held, held) == DONE)
    if (allocated(error)) return
    call check(error, return_fraction(2.0_DP * held, held) == DONE)
    if (allocated(error)) return
    call check(error, is_close(return_fraction(0.5_DP * held, held), 0.5_DP))
    if (allocated(error)) return
  end subroutine test_return_fraction_small

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

  !> @brief A repeated time step starts from the reservoirs the last one ended with
  !<
  subroutine test_advance_restore(error)
    type(error_type), allocatable, intent(out) :: error
    type(StrandedMassType) :: strand

    call strand%init('TESTSTRAND', 'MST-STRND', 2)
    !
    ! -- the end of a time step; the second cell has no new saturation
    strand%stranded_aqueous(1) = 1.0_DP
    strand%stranded_sorbed(1) = 2.0_DP
    strand%held_aq(1) = 0.3_DP
    strand%held_srb(1) = 0.4_DP
    strand%stranded_et(1) = 3.0_DP
    strand%held_et(1) = 0.5_DP
    strand%sat_new(1) = 0.6_DP
    strand%sat_old(2) = 0.8_DP
    call strand%advance()
    call check(error, strand%stranded_aqueous0(1) == 1.0_DP .and. &
               strand%stranded_sorbed0(1) == 2.0_DP .and. &
               strand%held_aq0(1) == 0.3_DP .and. &
               strand%held_srb0(1) == 0.4_DP, &
               'advance did not store the reservoirs')
    if (allocated(error)) return
    call check(error, strand%sat_old(1) == 0.6_DP .and. &
               strand%sat_old(2) == 0.8_DP .and. &
               strand%sat_new(2) == DNODATA, &
               'advance did not keep the saturation of each cell')
    if (allocated(error)) return
    !
    ! -- a time step that fails changes the reservoirs and is repeated
    strand%stranded_aqueous(1) = 5.0_DP
    strand%stranded_sorbed(1) = 6.0_DP
    strand%held_aq(1) = 0.7_DP
    strand%held_srb(1) = 0.8_DP
    strand%stranded_et(1) = 7.0_DP
    strand%held_et(1) = 0.9_DP
    call strand%restore()
    call check(error, strand%stranded_aqueous(1) == 1.0_DP .and. &
               strand%stranded_sorbed(1) == 2.0_DP .and. &
               strand%held_aq(1) == 0.3_DP .and. &
               strand%held_srb(1) == 0.4_DP .and. &
               strand%stranded_et(1) == 3.0_DP .and. &
               strand%held_et(1) == 0.5_DP, &
               'restore did not return the reservoirs to the start of the step')
    if (allocated(error)) return
    call check(error, strand%sat_old(1) == 0.6_DP, &
               'restore changed the saturation of the previous time step')
    call strand%da()
  end subroutine test_advance_restore

  !> @brief Recharge flushes a deposit as a well-mixed volume
  !<
  subroutine test_flush_fraction(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: x, f, g

    ! -- without recharge nothing is flushed and all of a deposit stays
    call check(error, flush_fraction(DZERO, 10.0_DP, 150.0_DP) == DZERO)
    if (allocated(error)) return
    call check(error, deposit_fraction(DZERO, 10.0_DP, 150.0_DP) == DONE)
    if (allocated(error)) return
    !
    ! -- a cell with no drained pore volume keeps nothing
    call check(error, flush_fraction(5.0_DP, 10.0_DP, DZERO) == DONE)
    if (allocated(error)) return
    call check(error, deposit_fraction(5.0_DP, 10.0_DP, DZERO) == DZERO)
    if (allocated(error)) return
    !
    ! -- at steady state the deposit retained over a time step, g x, equals
    !    the share of the deposit flushed, f, for any length of time step
    x = 5.0_DP * 10.0_DP / 150.0_DP
    f = flush_fraction(5.0_DP, 10.0_DP, 150.0_DP)
    g = deposit_fraction(5.0_DP, 10.0_DP, 150.0_DP)
    call check(error, is_close(f, DONE - exp(-x)), &
               'flush fraction '//to_string(f))
    if (allocated(error)) return
    call check(error, is_close(g * x, f), &
               'deposit fraction '//to_string(g))
    if (allocated(error)) return
    !
    ! -- a very small flush keeps nearly all of the deposit, without the
    !    loss of precision of the full expression
    g = deposit_fraction(1.0e-12_DP, 1.0_DP, 1.0_DP)
    call check(error, is_close(g, DONE - 0.5e-12_DP), &
               'small deposit fraction '//to_string(g))
  end subroutine test_flush_fraction

end module TestStrandedMass
