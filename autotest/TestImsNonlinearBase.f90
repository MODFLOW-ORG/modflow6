module TestImsNonlinearBase
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE, DEM20
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       to_string, unittest_type
  use ImsNonlinearBaseModule, only: ims_nl_underrelax, ims_nl_maxval, &
                                    ims_nl_calcdx, ims_nl_dxmax, &
                                    ims_nl_backtrack_flag, &
                                    ims_nl_apply_backtrack, &
                                    ims_nl_has_converged, &
                                    ims_nl_nur_has_converged, &
                                    ims_nl_residual, ims_nl_l2norm
  use VectorBaseModule, only: VectorBaseType
  use SparseMatrixModule, only: SparseMatrixType
  use SparseModule, only: sparsematrix
  implicit none
  private
  public :: collect_imsnonlinearbase

  real(DP), parameter :: tol = 1.0e-9_DP

contains

  subroutine collect_imsnonlinearbase(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("underrelax_simple", test_underrelax_simple), &
                new_unittest("underrelax_simple_inactive", &
                             test_underrelax_simple_inactive), &
                new_unittest("underrelax_cooley_kiter1", &
                             test_underrelax_cooley_kiter1), &
                new_unittest("underrelax_cooley_relax", &
                             test_underrelax_cooley_relax), &
                new_unittest("underrelax_dbd_kiter1", &
                             test_underrelax_dbd_kiter1), &
                new_unittest("underrelax_dbd_flipflop", &
                             test_underrelax_dbd_flipflop), &
                new_unittest("maxval_normalized", test_maxval_normalized), &
                new_unittest("maxval_zero_denom", test_maxval_zero_denom), &
                new_unittest("calcdx", test_calcdx), &
                new_unittest("dxmax", test_dxmax), &
                new_unittest("dxmax_signed", test_dxmax_signed), &
                new_unittest("backtrack_flag_on", test_backtrack_flag_on), &
                new_unittest("backtrack_flag_off", test_backtrack_flag_off), &
                new_unittest("apply_backtrack", test_apply_backtrack), &
                new_unittest("has_converged", test_has_converged), &
                new_unittest("nur_has_converged", test_nur_has_converged), &
                new_unittest("residual", test_residual), &
                new_unittest("l2norm", test_l2norm) &
                ]
  end subroutine collect_imsnonlinearbase

  !> @brief Build a 3x3 diagonal system A = diag(2,3,5), x = b = [1,1,1].
  !! Helper (not a registered test).
  !<
  subroutine build_diag_system(matrix, vec_x, vec_rhs, mem_path)
    type(SparseMatrixType), intent(inout) :: matrix
    class(VectorBaseType), pointer :: vec_x
    class(VectorBaseType), pointer :: vec_rhs
    character(len=*), intent(in) :: mem_path
    type(sparsematrix) :: sparse
    integer(I4B) :: i

    call sparse%init(3, 3, 1)
    do i = 1, 3
      call sparse%addconnection(i, i, 1)
    end do
    call matrix%init(sparse, mem_path)
    call sparse%destroy()

    call matrix%set_diag_value(1, 2.0_DP)
    call matrix%set_diag_value(2, 3.0_DP)
    call matrix%set_diag_value(3, 5.0_DP)

    vec_x => matrix%create_vec(3)
    vec_rhs => matrix%create_vec(3)
    do i = 1, 3
      call vec_x%set_value_local(i, 1.0_DP)
      call vec_rhs%set_value_local(i, 1.0_DP)
    end do
  end subroutine build_diag_system

  !> @brief Release a system built by build_diag_system. Helper.
  !<
  subroutine free_diag_system(matrix, vec_x, vec_rhs)
    type(SparseMatrixType), intent(inout) :: matrix
    class(VectorBaseType), pointer :: vec_x
    class(VectorBaseType), pointer :: vec_rhs

    call vec_x%destroy()
    deallocate (vec_x)
    call vec_rhs%destroy()
    deallocate (vec_rhs)
    call matrix%destroy()
  end subroutine free_diag_system

  !> @brief SIMPLE dampening (nonmeth == 1): x = xtemp + gamma * (x - xtemp)
  !<
  subroutine test_underrelax_simple(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: neq = 2
    integer(I4B) :: active(neq)
    real(DP) :: x(neq), xtemp(neq), dxold(neq), wsave(neq), hchold(neq)
    real(DP) :: deold(neq)
    real(DP) :: bigch_store, bigchold, relaxold
    real(DP) :: gamma, theta, akappa, amomentum

    gamma = 0.5_DP; theta = 0.7_DP; akappa = 0.1_DP; amomentum = 0.0_DP
    bigch_store = DZERO; bigchold = DZERO; relaxold = DZERO
    active = [1, 1]
    xtemp = [1.0_DP, 2.0_DP]
    x = [2.0_DP, 4.0_DP]
    dxold = DZERO; wsave = DZERO; hchold = DZERO; deold = DZERO

    call ims_nl_underrelax(1, 1, DONE, neq, active, x, xtemp, &
                           gamma, theta, akappa, amomentum, &
                           bigch_store, bigchold, relaxold, &
                           dxold, wsave, hchold, deold)

    ! node 1: delx = 1.0, x = 1.0 + 0.5*1.0 = 1.5
    call check(error, abs(x(1) - 1.5_DP) < tol, "x(1) expected 1.5")
    if (allocated(error)) return
    ! node 2: delx = 2.0, x = 2.0 + 0.5*2.0 = 3.0
    call check(error, abs(x(2) - 3.0_DP) < tol, "x(2) expected 3.0")
    if (allocated(error)) return
    call check(error, abs(dxold(1) - 1.0_DP) < tol, "dxold(1) expected 1.0")
    if (allocated(error)) return
    call check(error, abs(dxold(2) - 2.0_DP) < tol, "dxold(2) expected 2.0")
  end subroutine test_underrelax_simple

  !> @brief SIMPLE dampening leaves inactive nodes untouched
  !<
  subroutine test_underrelax_simple_inactive(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: neq = 2
    integer(I4B) :: active(neq)
    real(DP) :: x(neq), xtemp(neq), dxold(neq), wsave(neq), hchold(neq)
    real(DP) :: deold(neq)
    real(DP) :: bigch_store, bigchold, relaxold
    real(DP) :: gamma, theta, akappa, amomentum

    gamma = 0.5_DP; theta = 0.7_DP; akappa = 0.1_DP; amomentum = 0.0_DP
    bigch_store = DZERO; bigchold = DZERO; relaxold = DZERO
    active = [1, 0]
    xtemp = [1.0_DP, 2.0_DP]
    x = [2.0_DP, 6.0_DP]
    dxold = [-999.0_DP, -999.0_DP]
    wsave = DZERO; hchold = DZERO; deold = DZERO

    call ims_nl_underrelax(1, 1, DONE, neq, active, x, xtemp, &
                           gamma, theta, akappa, amomentum, &
                           bigch_store, bigchold, relaxold, &
                           dxold, wsave, hchold, deold)

    ! active node updated
    call check(error, abs(x(1) - 1.5_DP) < tol, "x(1) expected 1.5")
    if (allocated(error)) return
    ! inactive node unchanged (x and dxold both left alone)
    call check(error, abs(x(2) - 6.0_DP) < tol, "x(2) inactive unchanged")
    if (allocated(error)) return
    call check(error, abs(dxold(2) - (-999.0_DP)) < tol, &
               "dxold(2) inactive unchanged")
  end subroutine test_underrelax_simple_inactive

  !> @brief COOLEY (nonmeth == 2) first iteration: relax = 1, no x update,
  !! state seeded from bigch
  !<
  subroutine test_underrelax_cooley_kiter1(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: neq = 1
    integer(I4B) :: active(neq)
    real(DP) :: x(neq), xtemp(neq), dxold(neq), wsave(neq), hchold(neq)
    real(DP) :: deold(neq)
    real(DP) :: bigch_store, bigchold, relaxold
    real(DP) :: gamma, theta, akappa, amomentum

    gamma = 0.5_DP; theta = 0.7_DP; akappa = 0.1_DP; amomentum = 0.0_DP
    bigch_store = DZERO; bigchold = DZERO; relaxold = DZERO
    active = [1]
    xtemp = [0.0_DP]
    x = [10.0_DP]
    dxold = DZERO; wsave = DZERO; hchold = DZERO; deold = DZERO

    call ims_nl_underrelax(2, 1, 5.0_DP, neq, active, x, xtemp, &
                           gamma, theta, akappa, amomentum, &
                           bigch_store, bigchold, relaxold, &
                           dxold, wsave, hchold, deold)

    ! kiter==1: relax = 1.0 so x is not updated
    call check(error, abs(x(1) - 10.0_DP) < tol, "x(1) unchanged (relax=1)")
    if (allocated(error)) return
    call check(error, abs(relaxold - 1.0_DP) < tol, "relaxold expected 1.0")
    if (allocated(error)) return
    call check(error, abs(bigch_store - 5.0_DP) < tol, "bigch_store expected 5.0")
    if (allocated(error)) return
    ! bigchold = (1-gamma)*bigch + gamma*bigch = bigch = 5.0
    call check(error, abs(bigchold - 5.0_DP) < tol, "bigchold expected 5.0")
  end subroutine test_underrelax_cooley_kiter1

  !> @brief COOLEY subsequent iteration computes relaxation and updates x.
  !! State: bigchold=5, relaxold=1; bigch=-3, gamma=0.5.
  !! es = -3/5 = -0.6 (>= -1) => relax = (3-0.6)/(3+0.6) = 2/3
  !<
  subroutine test_underrelax_cooley_relax(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: neq = 1
    integer(I4B) :: active(neq)
    real(DP) :: x(neq), xtemp(neq), dxold(neq), wsave(neq), hchold(neq)
    real(DP) :: deold(neq)
    real(DP) :: bigch_store, bigchold, relaxold
    real(DP) :: gamma, theta, akappa, amomentum
    real(DP) :: relax_expect

    gamma = 0.5_DP; theta = 0.7_DP; akappa = 0.1_DP; amomentum = 0.0_DP
    bigch_store = DZERO; bigchold = 5.0_DP; relaxold = 1.0_DP
    active = [1]
    xtemp = [0.0_DP]
    x = [9.0_DP]
    dxold = DZERO; wsave = DZERO; hchold = DZERO; deold = DZERO

    call ims_nl_underrelax(2, 2, -3.0_DP, neq, active, x, xtemp, &
                           gamma, theta, akappa, amomentum, &
                           bigch_store, bigchold, relaxold, &
                           dxold, wsave, hchold, deold)

    relax_expect = 2.0_DP / 3.0_DP
    call check(error, abs(relaxold - relax_expect) < tol, "relaxold = 2/3")
    if (allocated(error)) return
    ! x = xtemp + relax*delx = 0 + (2/3)*9 = 6.0
    call check(error, abs(x(1) - 6.0_DP) < tol, "x(1) expected 6.0")
    if (allocated(error)) return
    call check(error, abs(dxold(1) - 9.0_DP) < tol, "dxold(1) expected 9.0")
    if (allocated(error)) return
    ! bigchold = (1-0.5)*(-3) + 0.5*5 = -1.5 + 2.5 = 1.0
    call check(error, abs(bigchold - 1.0_DP) < tol, "bigchold expected 1.0")
  end subroutine test_underrelax_cooley_relax

  !> @brief DBD (nonmeth == 3) first iteration seeds wsave/hchold/deold.
  !! akappa=0.1 => ww=min(1, 1+0.1)=1; x unchanged; state recorded.
  !<
  subroutine test_underrelax_dbd_kiter1(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: neq = 1
    integer(I4B) :: active(neq)
    real(DP) :: x(neq), xtemp(neq), dxold(neq), wsave(neq), hchold(neq)
    real(DP) :: deold(neq)
    real(DP) :: bigch_store, bigchold, relaxold
    real(DP) :: gamma, theta, akappa, amomentum

    gamma = 0.5_DP; theta = 0.7_DP; akappa = 0.1_DP; amomentum = 0.0_DP
    bigch_store = DZERO; bigchold = DZERO; relaxold = DZERO
    active = [1]
    xtemp = [1.0_DP]
    x = [3.0_DP]
    dxold = DZERO; wsave = DZERO; hchold = DZERO; deold = DZERO

    call ims_nl_underrelax(3, 1, DONE, neq, active, x, xtemp, &
                           gamma, theta, akappa, amomentum, &
                           bigch_store, bigchold, relaxold, &
                           dxold, wsave, hchold, deold)

    ! delx = 2.0; ww capped at 1.0; x = xtemp + 2.0*1.0 = 3.0 (unchanged)
    call check(error, abs(x(1) - 3.0_DP) < tol, "x(1) expected 3.0")
    if (allocated(error)) return
    call check(error, abs(wsave(1) - 1.0_DP) < tol, "wsave(1) expected 1.0")
    if (allocated(error)) return
    ! kiter==1: hchold = delx = 2.0
    call check(error, abs(hchold(1) - 2.0_DP) < tol, "hchold(1) expected 2.0")
    if (allocated(error)) return
    call check(error, abs(deold(1) - 2.0_DP) < tol, "deold(1) expected 2.0")
    if (allocated(error)) return
    call check(error, abs(dxold(1) - 2.0_DP) < tol, "dxold(1) expected 2.0")
  end subroutine test_underrelax_dbd_kiter1

  !> @brief DBD flip-flop: opposite-sign change decreases the relaxation
  !! factor via theta. State from a prior iteration: wsave=1, hchold=2, deold=2.
  !<
  subroutine test_underrelax_dbd_flipflop(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: neq = 1
    integer(I4B) :: active(neq)
    real(DP) :: x(neq), xtemp(neq), dxold(neq), wsave(neq), hchold(neq)
    real(DP) :: deold(neq)
    real(DP) :: bigch_store, bigchold, relaxold
    real(DP) :: gamma, theta, akappa, amomentum

    gamma = 0.5_DP; theta = 0.7_DP; akappa = 0.1_DP; amomentum = 0.0_DP
    bigch_store = DZERO; bigchold = DZERO; relaxold = DZERO
    active = [1]
    xtemp = [3.0_DP]
    x = [1.0_DP]
    dxold = DZERO
    wsave = [1.0_DP]
    hchold = [2.0_DP]
    deold = [2.0_DP]

    call ims_nl_underrelax(3, 2, -2.0_DP, neq, active, x, xtemp, &
                           gamma, theta, akappa, amomentum, &
                           bigch_store, bigchold, relaxold, &
                           dxold, wsave, hchold, deold)

    ! delx = -2.0; deold*delx = -4 < 0 (flip-flop) => ww = theta*1 = 0.7
    call check(error, abs(wsave(1) - 0.7_DP) < tol, "wsave(1) expected 0.7")
    if (allocated(error)) return
    ! kiter>1: hchold = (1-0.5)*(-2) + 0.5*2 = 0.0
    call check(error, abs(hchold(1) - 0.0_DP) < tol, "hchold(1) expected 0.0")
    if (allocated(error)) return
    ! amom = 0 (kiter<=4); delx = -2.0*0.7 + 0 = -1.4; x = 3.0 - 1.4 = 1.6
    call check(error, abs(x(1) - 1.6_DP) < tol, "x(1) expected 1.6")
    if (allocated(error)) return
    call check(error, abs(deold(1) - (-2.0_DP)) < tol, "deold(1) expected -2.0")
  end subroutine test_underrelax_dbd_flipflop

  !> @brief maxval uses a normalized comparison to pick the largest magnitude
  !<
  subroutine test_maxval_normalized(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: vmax

    ! vmax starts at 1.0; -3.0 has larger magnitude => vmax = -3.0;
    ! 2.0 has smaller magnitude than 3.0 => vmax unchanged (signed result kept)
    vmax = DZERO
    call ims_nl_maxval(3, [1.0_DP, -3.0_DP, 2.0_DP], vmax)
    call check(error, abs(vmax - (-3.0_DP)) < tol, "vmax expected -3.0")
  end subroutine test_maxval_normalized

  !> @brief maxval handles a zero running value without dividing by zero
  !<
  subroutine test_maxval_zero_denom(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: vmax

    ! first value 0 => denom falls back to DPREC, so any nonzero later value wins
    vmax = DZERO
    call ims_nl_maxval(2, [0.0_DP, 0.5_DP], vmax)
    call check(error, abs(vmax - 0.5_DP) < tol, "vmax expected 0.5")
  end subroutine test_maxval_zero_denom

  !> @brief calcdx returns x - xtemp for active cells, zero for inactive
  !<
  subroutine test_calcdx(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: neq = 3
    real(DP) :: dx(neq)

    dx = -999.0_DP
    call ims_nl_calcdx(neq, [1, 0, 1], [5.0_DP, 5.0_DP, 5.0_DP], &
                       [1.0_DP, 1.0_DP, 1.0_DP], dx)
    call check(error, abs(dx(1) - 4.0_DP) < tol, "dx(1) expected 4.0")
    if (allocated(error)) return
    call check(error, abs(dx(2) - 0.0_DP) < tol, "dx(2) inactive => 0.0")
    if (allocated(error)) return
    call check(error, abs(dx(3) - 4.0_DP) < tol, "dx(3) expected 4.0")
  end subroutine test_calcdx

  !> @brief dxmax returns the largest-magnitude change and its location,
  !! excluding inactive cells
  !<
  subroutine test_dxmax(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: neq = 3
    real(DP) :: dxmax
    integer(I4B) :: loc

    dxmax = DZERO
    loc = 0
    ! node 3 has the biggest change but is inactive => excluded
    call ims_nl_dxmax(neq, [1, 1, 0], [2.0_DP, 5.0_DP, 100.0_DP], &
                      [1.0_DP, 1.0_DP, 1.0_DP], dxmax, loc)
    call check(error, abs(dxmax - 4.0_DP) < tol, "dxmax expected 4.0")
    if (allocated(error)) return
    call check(error, loc == 2, "loc expected 2")
  end subroutine test_dxmax

  !> @brief dxmax preserves the sign of the maximum change
  !<
  subroutine test_dxmax_signed(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: neq = 2
    real(DP) :: dxmax
    integer(I4B) :: loc

    dxmax = DZERO
    loc = 0
    call ims_nl_dxmax(neq, [1, 1], [1.0_DP, 1.0_DP], &
                      [5.0_DP, 1.0_DP], dxmax, loc)
    call check(error, abs(dxmax - (-4.0_DP)) < tol, "dxmax expected -4.0")
    if (allocated(error)) return
    call check(error, loc == 1, "loc expected 1")
  end subroutine test_dxmax_signed

  !> @brief backtracking flag set when reduced max change exceeds closure;
  !! inactive cells excluded from the max
  !<
  subroutine test_backtrack_flag_on(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: bt_flag

    ! active max change = 9 (node 3 inactive, ignored); 0.5*9 = 4.5 >= 1.0
    bt_flag = ims_nl_backtrack_flag(3, [1, 1, 0], &
                                    [2.0_DP, 10.0_DP, 1000.0_DP], &
                                    [1.0_DP, 1.0_DP, 1.0_DP], 0.5_DP, 1.0_DP)
    call check(error, bt_flag == 1, "bt_flag expected 1")
  end subroutine test_backtrack_flag_on

  !> @brief backtracking flag stays off when reduced max change is below closure
  !<
  subroutine test_backtrack_flag_off(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: bt_flag

    ! 0.05*9 = 0.45 < 1.0
    bt_flag = ims_nl_backtrack_flag(3, [1, 1, 0], &
                                    [2.0_DP, 10.0_DP, 1000.0_DP], &
                                    [1.0_DP, 1.0_DP, 1.0_DP], 0.05_DP, 1.0_DP)
    call check(error, bt_flag == 0, "bt_flag expected 0")
  end subroutine test_backtrack_flag_off

  !> @brief apply_backtrack scales the change by breduc for active cells only
  !<
  subroutine test_apply_backtrack(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: neq = 2
    real(DP) :: x(neq)

    x = [10.0_DP, 10.0_DP]
    call ims_nl_apply_backtrack(neq, [1, 0], x, [0.0_DP, 0.0_DP], 0.25_DP)
    ! active: x = 0 + 0.25*10 = 2.5
    call check(error, abs(x(1) - 2.5_DP) < tol, "x(1) expected 2.5")
    if (allocated(error)) return
    ! inactive: unchanged
    call check(error, abs(x(2) - 10.0_DP) < tol, "x(2) inactive unchanged")
  end subroutine test_apply_backtrack

  !> @brief convergence when |max_dvc| <= dvclose (sign-independent)
  !<
  subroutine test_has_converged(error)
    type(error_type), allocatable, intent(out) :: error

    call check(error, ims_nl_has_converged(0.0005_DP, 0.001_DP), &
               "0.0005 within 0.001 => converged")
    if (allocated(error)) return
    call check(error, ims_nl_has_converged(-0.0005_DP, 0.001_DP), &
               "-0.0005 within 0.001 => converged")
    if (allocated(error)) return
    call check(error, ims_nl_has_converged(0.001_DP, 0.001_DP), &
               "boundary equal => converged")
    if (allocated(error)) return
    call check(error,.not. ims_nl_has_converged(0.01_DP, 0.001_DP), &
               "0.01 exceeds 0.001 => not converged")
  end subroutine test_has_converged

  !> @brief NUR convergence requires BOTH changes within closure
  !<
  subroutine test_nur_has_converged(error)
    type(error_type), allocatable, intent(out) :: error

    call check(error, &
               ims_nl_nur_has_converged(0.0005_DP, 0.0005_DP, 0.001_DP), &
               "both within => converged")
    if (allocated(error)) return
    call check(error, &
               .not. ims_nl_nur_has_converged(0.0005_DP, 0.01_DP, 0.001_DP), &
               "hncg exceeds => not converged")
    if (allocated(error)) return
    call check(error, &
               .not. ims_nl_nur_has_converged(0.01_DP, 0.0005_DP, 0.001_DP), &
               "dxold_max exceeds => not converged")
  end subroutine test_nur_has_converged

  !> @brief residual r = A*x - b, zeroed for inactive cells.
  !! A = diag(2,3,5), x = b = [1,1,1] => r = [1,2,4]
  !<
  subroutine test_residual(error)
    type(error_type), allocatable, intent(out) :: error
    type(SparseMatrixType) :: matrix
    class(VectorBaseType), pointer :: vec_x, vec_rhs, vec_resid
    real(DP) :: r_all(3), r_inact(3)
    integer(I4B) :: i

    call build_diag_system(matrix, vec_x, vec_rhs, 'TEST/RESID')
    vec_resid => matrix%create_vec(3)

    ! all active: r = diag(2,3,5)*[1,1,1] - [1,1,1] = [1,2,4]
    call ims_nl_residual(matrix, vec_x, vec_rhs, 3, [1, 1, 1], vec_resid)
    do i = 1, 3
      r_all(i) = vec_resid%get_value_local(i)
    end do

    ! middle cell inactive: r(2) forced to 0, others unchanged
    call ims_nl_residual(matrix, vec_x, vec_rhs, 3, [1, 0, 1], vec_resid)
    do i = 1, 3
      r_inact(i) = vec_resid%get_value_local(i)
    end do

    ! free resources before asserting so a failure cannot leak
    call vec_resid%destroy()
    deallocate (vec_resid)
    call free_diag_system(matrix, vec_x, vec_rhs)

    call check(error, abs(r_all(1) - 1.0_DP) < tol, "r(1) expected 1")
    if (allocated(error)) return
    call check(error, abs(r_all(2) - 2.0_DP) < tol, "r(2) expected 2")
    if (allocated(error)) return
    call check(error, abs(r_all(3) - 4.0_DP) < tol, "r(3) expected 4")
    if (allocated(error)) return
    call check(error, abs(r_inact(2) - 0.0_DP) < tol, "r(2) inactive => 0")
    if (allocated(error)) return
    call check(error, abs(r_inact(1) - 1.0_DP) < tol, "r(1) unchanged")
    if (allocated(error)) return
    call check(error, abs(r_inact(3) - 4.0_DP) < tol, "r(3) unchanged")
  end subroutine test_residual

  !> @brief l2norm = || A*x - b ||_2. For r = [1,2,4], ||r|| = sqrt(21)
  !<
  subroutine test_l2norm(error)
    type(error_type), allocatable, intent(out) :: error
    type(SparseMatrixType) :: matrix
    class(VectorBaseType), pointer :: vec_x, vec_rhs
    real(DP) :: l2

    call build_diag_system(matrix, vec_x, vec_rhs, 'TEST/L2NORM')
    l2 = -1.0_DP
    call ims_nl_l2norm(matrix, vec_x, vec_rhs, 3, [1, 1, 1], l2)
    call free_diag_system(matrix, vec_x, vec_rhs)

    call check(error, abs(l2 - sqrt(21.0_DP)) < tol, "l2norm expected sqrt(21)")
  end subroutine test_l2norm

end module TestImsNonlinearBase
