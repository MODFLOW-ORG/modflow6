!> @brief Unit tests for IMS linear solver preconditioners
!!
!! Tests ILU0, MILU0, ILUT, MILUT, and AMG preconditioners on a small
!! 1D Poisson tridiagonal system. For the ILU family the factorization
!! is exact on a tridiagonal matrix (zero fill), so the solve residual
!! should be at machine-precision level. The AMG tests exercise the
!! V-cycle with the ILU0 and ILU0_ALL smoothers, and verify
!! convergence in Richardson iterations.
!<
module TestImsLinearMethods
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use IMSLinearBaseModule, only: ims_base_pccrs, ims_base_pcilu0, &
                                 ims_base_ilu0a, ims_calc_pcdims
  use ImsLinearAmgModule, only: ImsAmgDataType
  use ImsLinearSettingsModule, only: SMOOTHER_ILU0, SMOOTHER_ILU0_ALL

  implicit none
  private
  public :: collect_imslinearmethods

  integer(I4B), parameter :: N_TEST = 20 !< unknowns in the test system

  !> Explicit interfaces for the SPARSKIT2 routines (no module wrapper)
  interface
    subroutine ilut(n, a, ja, ia, lfil, droptol, alu, jlu, ju, iwk, &
                    w, jw, ierr, relax, izero, delta)
      integer(kind=4) :: n, lfil, iwk, ierr, izero
      real(kind=8) :: droptol, relax, delta
      real(kind=8) :: a(*), alu(*), w(*)
      integer(kind=4) :: ja(*), ia(*), jlu(*), ju(*), jw(*)
    end subroutine ilut
    subroutine lusol(n, y, x, alu, jlu, ju)
      integer(kind=4) :: n
      real(kind=8) :: y(*), x(*), alu(*)
      integer(kind=4) :: jlu(*), ju(*)
    end subroutine lusol
  end interface

contains

  subroutine collect_imslinearmethods(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("ilu0_residual", test_ilu0_residual), &
                new_unittest("milu0_residual", test_milu0_residual), &
                new_unittest("ilut_residual", test_ilut_residual), &
                new_unittest("milut_residual", test_milut_residual), &
                new_unittest("amg_fixed_omega", test_amg_fixed_omega), &
                new_unittest("amg_adaptive_omega", test_amg_adaptive_omega), &
                new_unittest("amg_hierarchy_levels", &
                             test_amg_hierarchy_levels), &
                new_unittest("amg_nsmooth_1", test_amg_nsmooth_1), &
                new_unittest("amg_nonsymmetric_soc", &
                             test_amg_nonsymmetric_soc), &
                new_unittest("amg_sthresh_convergence", &
                             test_amg_sthresh_convergence), &
                new_unittest("amg_sthresh_hierarchy", &
                             test_amg_sthresh_hierarchy), &
                new_unittest("amg_ilu0all_converges", &
                             test_amg_ilu0all_converges) &
                ]
  end subroutine collect_imslinearmethods

  ! ---------------------------------------------------------------------------
  ! Helpers
  ! ---------------------------------------------------------------------------

  !> @brief Build a 1D heterogeneous diffusion tridiagonal CSR matrix of size n
  !!
  !! Builds the matrix for a 1D diffusion problem, -d/dx(K du/dx) = f, with
  !! the value fixed to zero at both ends. The conductance K at each cell
  !! interface alternates between a high value (10) and a low value (0.1),
  !! so every pair of neighbouring interfaces differs by 100:1 -- a hard
  !! case for the solver.
  !!
  !! @param[in]  n       number of unknowns
  !! @param[out] ia      CSR row pointers (n+1)
  !! @param[out] ja      CSR column indices (3n-2)
  !! @param[out] a       CSR values (3n-2)
  !! @param[out] b       right-hand side vector b = A * x_exact
  !! @param[out] x_exact manufactured exact solution x(i) = sin(i*pi/(n+1))
  !<
  subroutine build_poisson1d(n, ia, ja, a, b, x_exact)
    integer(I4B), intent(in) :: n
    integer(I4B), intent(out) :: ia(n + 1)
    integer(I4B), allocatable, intent(out) :: ja(:)
    real(DP), allocatable, intent(out) :: a(:)
    real(DP), intent(out) :: b(n), x_exact(n)
    integer(I4B) :: nja, i, pos
    real(DP) :: pi
    real(DP) :: kf(n - 1) !< conductance K at each interface, i = 1 .. n-1

    ! alternate between high (10) and low (0.1) so every large-K interface
    ! is bordered by a small-K interface and vice versa
    do i = 1, n - 1
      if (mod(i, 2) == 1) then
        kf(i) = 10.0_DP
      else
        kf(i) = 0.1_DP
      end if
    end do

    nja = 3 * n - 2
    allocate (ja(nja), a(nja))

    ia(1) = 1
    ia(2) = 3
    do i = 3, n
      ia(i) = ia(i - 1) + 3
    end do
    ia(n + 1) = nja + 1

    pos = 1
    ! row 1: left boundary - diagonal first (MODFLOW 6 CSR: diagonal at ia(n))
    ja(pos) = 1; a(pos) = kf(1); pos = pos + 1
    ja(pos) = 2; a(pos) = -kf(1); pos = pos + 1
    ! interior rows: diagonal first, then off-diagonals
    do i = 2, n - 1
      ja(pos) = i; a(pos) = kf(i - 1) + kf(i); pos = pos + 1
      ja(pos) = i - 1; a(pos) = -kf(i - 1); pos = pos + 1
      ja(pos) = i + 1; a(pos) = -kf(i); pos = pos + 1
    end do
    ! row n: right boundary - diagonal first
    ja(pos) = n; a(pos) = kf(n - 1); pos = pos + 1
    ja(pos) = n - 1; a(pos) = -kf(n - 1)

    pi = acos(-1.0_DP)
    do i = 1, n
      x_exact(i) = sin(pi * real(i, DP) / real(n + 1, DP))
    end do

    ! b = A * x_exact (manufactured RHS)
    b(1) = kf(1) * x_exact(1) - kf(1) * x_exact(2)
    do i = 2, n - 1
      b(i) = -kf(i - 1) * x_exact(i - 1) + &
             (kf(i - 1) + kf(i)) * x_exact(i) - &
             kf(i) * x_exact(i + 1)
    end do
    b(n) = -kf(n - 1) * x_exact(n - 1) + kf(n - 1) * x_exact(n)
  end subroutine build_poisson1d

  !> @brief Compute the relative residual ||b - A*x|| / ||b||
  !!
  !! @param[in] n   number of unknowns
  !! @param[in] nja number of nonzeros
  !! @param[in] ia  CSR row pointers
  !! @param[in] ja  CSR column indices
  !! @param[in] a   CSR values
  !! @param[in] b   right-hand side
  !! @param[in] x   approximate solution
  !<
  real(DP) function rel_residual(n, nja, ia, ja, a, b, x) result(rr)
    integer(I4B), intent(in) :: n, nja
    integer(I4B), intent(in) :: ia(n + 1), ja(nja)
    real(DP), intent(in) :: a(nja), b(n), x(n)
    real(DP) :: ax(n)
    integer(I4B) :: i, j
    real(DP) :: b_sq, r_sq

    do i = 1, n
      ax(i) = DZERO
      do j = ia(i), ia(i + 1) - 1
        ax(i) = ax(i) + a(j) * x(ja(j))
      end do
    end do
    b_sq = DZERO
    r_sq = DZERO
    do i = 1, n
      b_sq = b_sq + b(i) * b(i)
      r_sq = r_sq + (b(i) - ax(i))**2
    end do
    if (b_sq > DZERO) then
      rr = sqrt(r_sq / b_sq)
    else
      rr = sqrt(r_sq)
    end if
  end function rel_residual

  !> @brief Build a 1D upwind convection-diffusion CSR matrix of size n
  !!
  !! Builds the matrix for a 1D problem with both diffusion (D=1) and
  !! advection (v=5), with the value fixed to zero at both ends. Advection
  !! uses upwind differencing, which makes the matrix strongly
  !! nonsymmetric: a_{i,i-1} = -(D+v) = -6 but a_{i-1,i} = -D = -1. The
  !! ratio v/D = 5 is the cell Peclet number (advection vs diffusion).
  !<
  subroutine build_convdiff1d(n, ia, ja, a, b, x_exact)
    integer(I4B), intent(in) :: n
    integer(I4B), intent(out) :: ia(n + 1)
    integer(I4B), allocatable, intent(out) :: ja(:)
    real(DP), allocatable, intent(out) :: a(:)
    real(DP), intent(out) :: b(n), x_exact(n)
    real(DP), parameter :: D_COEF = 1.0_DP, V_ADV = 5.0_DP
    integer(I4B) :: nja, i, pos
    real(DP) :: pi
    real(DP) :: a_diag, a_lower, a_upper

    a_diag = 2.0_DP * D_COEF + V_ADV
    a_lower = -(D_COEF + V_ADV)
    a_upper = -D_COEF

    nja = 3 * n - 2
    allocate (ja(nja), a(nja))

    ia(1) = 1
    ia(2) = 3
    do i = 3, n
      ia(i) = ia(i - 1) + 3
    end do
    ia(n + 1) = nja + 1

    pi = acos(-1.0_DP)
    do i = 1, n
      x_exact(i) = sin(pi * real(i, DP) / real(n + 1, DP))
    end do

    pos = 1
    ! row 1: diagonal first, then upper off-diagonal (no lower neighbor)
    ja(pos) = 1; a(pos) = D_COEF + V_ADV; pos = pos + 1
    ja(pos) = 2; a(pos) = a_upper; pos = pos + 1
    ! interior rows: diagonal first, then lower, then upper off-diagonal
    do i = 2, n - 1
      ja(pos) = i; a(pos) = a_diag; pos = pos + 1
      ja(pos) = i - 1; a(pos) = a_lower; pos = pos + 1
      ja(pos) = i + 1; a(pos) = a_upper; pos = pos + 1
    end do
    ! row n: diagonal first, then lower off-diagonal (no upper neighbor)
    ja(pos) = n; a(pos) = a_diag; pos = pos + 1
    ja(pos) = n - 1; a(pos) = a_lower

    ! b = A * x_exact (manufactured RHS)
    b(1) = (D_COEF + V_ADV) * x_exact(1) + a_upper * x_exact(2)
    do i = 2, n - 1
      b(i) = a_lower * x_exact(i - 1) + a_diag * x_exact(i) &
             + a_upper * x_exact(i + 1)
    end do
    b(n) = a_lower * x_exact(n - 1) + a_diag * x_exact(n)
  end subroutine build_convdiff1d

  ! ---------------------------------------------------------------------------
  ! Tests
  ! ---------------------------------------------------------------------------

  !> @brief ILU0 factorization and solve of a tridiagonal system
  !!
  !! ILU0 adds no new nonzeros (zero fill), so on a tridiagonal matrix it
  !! is an exact factorization and solving A*z = b is accurate to near
  !! machine precision.
  !<
  subroutine test_ilu0_residual(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: ia(N_TEST + 1), iapc(N_TEST + 1)
    integer(I4B), allocatable :: ja(:), japc(:)
    real(DP), allocatable :: a(:), apc(:)
    real(DP) :: b(N_TEST), x_exact(N_TEST), z(N_TEST), w(N_TEST)
    integer(I4B) :: iw(N_TEST), ipcflag, nja, niapc, njapc, njlu, njw, nwlu
    integer(I4B) :: i

    call build_poisson1d(N_TEST, ia, ja, a, b, x_exact)
    nja = size(ja)

    call ims_calc_pcdims(N_TEST, nja, ia, 0, 1, niapc, njapc, njlu, njw, nwlu)
    allocate (japc(njapc), apc(njapc))
    do i = 1, njapc
      japc(i) = 0
      apc(i) = DZERO
    end do

    call ims_base_pccrs(N_TEST, nja, ia, ja, iapc, japc)
    ipcflag = 0
    call ims_base_pcilu0(nja, N_TEST, a, ia, ja, apc, iapc, japc, &
                         iw, w, 0.0_DP, ipcflag, 0.0_DP)
    call check(error, ipcflag == 0, &
               "ILU0 factorization should succeed without pivot fix")
    if (allocated(error)) return

    call ims_base_ilu0a(nja, N_TEST, apc, iapc, japc, b, z)

    call check(error, rel_residual(N_TEST, nja, ia, ja, a, b, z) < 1.0e-12_DP, &
               "ILU0 solve of tridiagonal should be "// &
               "exact to near machine precision")
  end subroutine test_ilu0_residual

  !> @brief MILU0 factorization and solve of a tridiagonal system
  !!
  !! For a tridiagonal matrix there are no dropped off-diagonal terms, so
  !! the MILU0 relaxation has no effect and the result is identical to ILU0.
  !<
  subroutine test_milu0_residual(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: ia(N_TEST + 1), iapc(N_TEST + 1)
    integer(I4B), allocatable :: ja(:), japc(:)
    real(DP), allocatable :: a(:), apc(:)
    real(DP) :: b(N_TEST), x_exact(N_TEST), z(N_TEST), w(N_TEST)
    integer(I4B) :: iw(N_TEST), ipcflag, nja, niapc, njapc, njlu, njw, nwlu
    integer(I4B) :: i

    call build_poisson1d(N_TEST, ia, ja, a, b, x_exact)
    nja = size(ja)

    call ims_calc_pcdims(N_TEST, nja, ia, 0, 2, niapc, njapc, njlu, njw, nwlu)
    allocate (japc(njapc), apc(njapc))
    do i = 1, njapc
      japc(i) = 0
      apc(i) = DZERO
    end do

    call ims_base_pccrs(N_TEST, nja, ia, ja, iapc, japc)
    ipcflag = 0
    call ims_base_pcilu0(nja, N_TEST, a, ia, ja, apc, iapc, japc, &
                         iw, w, 1.0_DP, ipcflag, 0.0_DP)
    call check(error, ipcflag == 0, &
               "MILU0 factorization should succeed without pivot fix")
    if (allocated(error)) return

    call ims_base_ilu0a(nja, N_TEST, apc, iapc, japc, b, z)

    call check(error, rel_residual(N_TEST, nja, ia, ja, a, b, z) < 1.0e-12_DP, &
               "MILU0 solve of tridiagonal (no dropped terms) should be exact")
  end subroutine test_milu0_residual

  !> @brief ILUT factorization and solve of a tridiagonal system
  !!
  !! ILUT keeps only the largest `lfil` entries per row. On a tridiagonal
  !! matrix each row has just one upper off-diagonal, and MODFLOW's ILUT
  !! keeps lfil-1 of these, so lfil=2 is needed to keep that one entry and
  !! factor the matrix exactly.
  !<
  subroutine test_ilut_residual(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: ia(N_TEST + 1)
    integer(I4B), allocatable :: ja(:), jlu(:), ju(:), jw(:)
    real(DP), allocatable :: a(:), alu(:), w_lu(:)
    real(DP) :: b(N_TEST), x_exact(N_TEST), z(N_TEST)
    integer(I4B) :: nja, niapc, njapc, njlu, njw, nwlu
    integer(I4B) :: ierr, izero, i

    call build_poisson1d(N_TEST, ia, ja, a, b, x_exact)
    nja = size(ja)

    call ims_calc_pcdims(N_TEST, nja, ia, 2, 3, niapc, njapc, njlu, njw, nwlu)
    allocate (alu(njapc), jlu(njlu), ju(N_TEST), jw(njw), w_lu(nwlu))
    do i = 1, njapc
      alu(i) = DZERO
    end do
    do i = 1, njlu
      jlu(i) = 0
    end do
    do i = 1, N_TEST
      ju(i) = 0
    end do
    do i = 1, njw
      jw(i) = 0
    end do
    do i = 1, nwlu
      w_lu(i) = DZERO
    end do

    ierr = 0; izero = 0
    call ilut(N_TEST, a, ja, ia, 2, 0.0_DP, alu, jlu, ju, njapc, w_lu, jw, &
              ierr, 0.0_DP, izero, 0.0_DP)
    call check(error, ierr == 0, "ILUT factorization should succeed")
    if (allocated(error)) return

    call lusol(N_TEST, b, z, alu, jlu, ju)

    call check(error, rel_residual(N_TEST, nja, ia, ja, a, b, z) < 1.0e-12_DP, &
               "ILUT(lfil=2) solve of tridiagonal should be exact")
  end subroutine test_ilut_residual

  !> @brief MILUT factorization and solve of a tridiagonal system
  !!
  !! MILUT is ILUT with relaxation (relax=1.0) that redistributes dropped
  !! terms onto the diagonal. With lfil=2 and droptol=0 on a tridiagonal
  !! there are no dropped terms and the factorization is exact.
  !<
  subroutine test_milut_residual(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: ia(N_TEST + 1)
    integer(I4B), allocatable :: ja(:), jlu(:), ju(:), jw(:)
    real(DP), allocatable :: a(:), alu(:), w_lu(:)
    real(DP) :: b(N_TEST), x_exact(N_TEST), z(N_TEST)
    integer(I4B) :: nja, niapc, njapc, njlu, njw, nwlu
    integer(I4B) :: ierr, izero, i

    call build_poisson1d(N_TEST, ia, ja, a, b, x_exact)
    nja = size(ja)

    call ims_calc_pcdims(N_TEST, nja, ia, 2, 4, niapc, njapc, njlu, njw, nwlu)
    allocate (alu(njapc), jlu(njlu), ju(N_TEST), jw(njw), w_lu(nwlu))
    do i = 1, njapc
      alu(i) = DZERO
    end do
    do i = 1, njlu
      jlu(i) = 0
    end do
    do i = 1, N_TEST
      ju(i) = 0
    end do
    do i = 1, njw
      jw(i) = 0
    end do
    do i = 1, nwlu
      w_lu(i) = DZERO
    end do

    ierr = 0; izero = 0
    call ilut(N_TEST, a, ja, ia, 2, 0.0_DP, alu, jlu, ju, njapc, w_lu, jw, &
              ierr, 1.0_DP, izero, 0.0_DP)
    call check(error, ierr == 0, "MILUT factorization should succeed")
    if (allocated(error)) return

    call lusol(N_TEST, b, z, alu, jlu, ju)

    call check(error, rel_residual(N_TEST, nja, ia, ja, a, b, z) < 1.0e-12_DP, &
               "MILUT(lfil=2) solve of tridiagonal "// &
               "(no dropped terms) should be exact")
  end subroutine test_milut_residual

  !> @brief Run N_ITER AMG Richardson iterations and return relative residual
  !!
  !! Defaults to a V-cycle with the ILU0 smoother. An optional argument lets
  !! the smoother be changed for the ILU0_ALL test.
  !<
  real(DP) function amg_vcycle_residual(ia, ja, a, omega, &
                                        nsmooth, n_iter, sthresh, &
                                        smoother_type_opt) result(rr)
    integer(I4B), intent(in) :: ia(N_TEST + 1)
    integer(I4B), intent(in) :: ja(:)
    real(DP), intent(in) :: a(:)
    real(DP), intent(in) :: omega !< damping factor; 0.0 = adaptive optimal
    integer(I4B), intent(in) :: nsmooth !< smoother sweeps per level
    integer(I4B), intent(in) :: n_iter
    real(DP), intent(in), optional :: sthresh !< aggregation strength threshold
    integer(I4B), intent(in), optional :: smoother_type_opt !< SMOOTHER_ILU0 etc.
    real(DP) :: b(N_TEST), x_exact(N_TEST)
    real(DP) :: x(N_TEST), r(N_TEST), z(N_TEST)
    integer(I4B) :: nja, k, i, j
    real(DP) :: ax_i, pi, sthresh_
    integer(I4B) :: smoother_type_
    type(ImsAmgDataType) :: amg

    nja = size(ja)
    pi = acos(-1.0_DP)
    do i = 1, N_TEST
      x_exact(i) = sin(pi * real(i, DP) / real(N_TEST + 1, DP))
    end do
    ! compute b = A * x_exact from the passed matrix (works for any K field)
    do i = 1, N_TEST
      b(i) = DZERO
      do j = ia(i), ia(i + 1) - 1
        b(i) = b(i) + a(j) * x_exact(ja(j))
      end do
    end do

    sthresh_ = 0.0_DP
    if (present(sthresh)) sthresh_ = sthresh
    smoother_type_ = SMOOTHER_ILU0
    if (present(smoother_type_opt)) smoother_type_ = smoother_type_opt
    call amg%setup(N_TEST, nja, ia, ja, a, max_levels=5, &
                   nsmooth=nsmooth, omega=omega, sthresh=sthresh_, &
                   smoother_type=smoother_type_)

    do i = 1, N_TEST
      x(i) = DZERO
    end do
    do k = 1, n_iter
      do i = 1, N_TEST
        ax_i = DZERO
        do j = ia(i), ia(i + 1) - 1
          ax_i = ax_i + a(j) * x(ja(j))
        end do
        r(i) = b(i) - ax_i
      end do
      call amg%apply(N_TEST, r, z)
      do i = 1, N_TEST
        x(i) = x(i) + z(i)
      end do
    end do

    call amg%destroy()
    rr = rel_residual(N_TEST, nja, ia, ja, a, b, x)
  end function amg_vcycle_residual

  !> @brief AMG with fixed undamped Gauss-Seidel (omega=1) converges
  !!
  !! Uses RELAXATION_FACTOR = 1.0 (undamped GS) and verifies that 15
  !! V-cycle Richardson iterations reduce the relative residual below 1e-6.
  !<
  subroutine test_amg_fixed_omega(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: ia(N_TEST + 1)
    integer(I4B), allocatable :: ja(:)
    real(DP), allocatable :: a(:)
    real(DP) :: b(N_TEST), x_exact(N_TEST)
    integer(I4B) :: nja

    call build_poisson1d(N_TEST, ia, ja, a, b, x_exact)
    nja = size(ja)

    call check(error, &
               amg_vcycle_residual(ia, ja, a, omega=1.0_DP, &
                                   nsmooth=4, n_iter=15) < 1.0e-6_DP, &
               "AMG fixed omega=1: relative residual "// &
               "should fall below 1e-6 in 15 cycles")
  end subroutine test_amg_fixed_omega

  !> @brief AMG with per-sweep optimal adaptive omega (omega=0) converges
  !!
  !! Uses RELAXATION_FACTOR = 0.0 which activates per-application optimal
  !! damping: alpha = (r,z)/(z,Az) capped at 1.0.  Verifies convergence to
  !! below 1e-6 in at most 15 V-cycles on the 1D Poisson system.
  !<
  subroutine test_amg_adaptive_omega(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: ia(N_TEST + 1)
    integer(I4B), allocatable :: ja(:)
    real(DP), allocatable :: a(:)
    real(DP) :: b(N_TEST), x_exact(N_TEST)
    integer(I4B) :: nja

    call build_poisson1d(N_TEST, ia, ja, a, b, x_exact)
    nja = size(ja)

    call check(error, &
               amg_vcycle_residual(ia, ja, a, omega=0.0_DP, &
                                   nsmooth=4, n_iter=15) < 1.0e-6_DP, &
               "AMG adaptive omega=0: relative residual "// &
               "should fall below 1e-6 in 15 cycles")
  end subroutine test_amg_adaptive_omega

  !> @brief AMG hierarchy has at least 2 levels for the N_TEST system
  !!
  !! For N_TEST=20 the solver pairs up the 20 unknowns into ~10 coarse
  !! nodes, giving a two-level hierarchy. Verifies nlevels >= 2 and that
  !! the coarsest level is strictly smaller than the fine level.
  !<
  subroutine test_amg_hierarchy_levels(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: ia(N_TEST + 1)
    integer(I4B), allocatable :: ja(:)
    real(DP), allocatable :: a(:)
    real(DP) :: b(N_TEST), x_exact(N_TEST)
    integer(I4B) :: nja
    type(ImsAmgDataType) :: amg

    call build_poisson1d(N_TEST, ia, ja, a, b, x_exact)
    nja = size(ja)

    call amg%setup(N_TEST, nja, ia, ja, a, max_levels=5, nsmooth=2, &
                   omega=1.0_DP, sthresh=0.0_DP, &
                   smoother_type=SMOOTHER_ILU0)

    call check(error, amg%nlevels >= 2, &
               "AMG on N_TEST=20 should build at least 2 hierarchy levels")
    if (allocated(error)) then
      call amg%destroy()
      return
    end if

    call check(error, amg%levels(amg%nlevels - 1)%neq < N_TEST, &
               "Coarsest level should have fewer unknowns than fine level")

    call amg%destroy()
  end subroutine test_amg_hierarchy_levels

  !> @brief AMG with NUMBER_SMOOTHING_ITERATIONS = 1 converges
  !!
  !! nsmooth=1 is the minimum allowed value. Convergence is slower than
  !! the default (nsmooth=2), so 30 V-cycles are allowed. Verifies the
  !! nsmooth parameter is wired through and the solver still converges.
  !<
  subroutine test_amg_nsmooth_1(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: ia(N_TEST + 1)
    integer(I4B), allocatable :: ja(:)
    real(DP), allocatable :: a(:)
    real(DP) :: b(N_TEST), x_exact(N_TEST)
    integer(I4B) :: nja

    call build_poisson1d(N_TEST, ia, ja, a, b, x_exact)
    nja = size(ja)

    call check(error, &
               amg_vcycle_residual(ia, ja, a, omega=0.0_DP, &
                                   nsmooth=1, n_iter=30) < 1.0e-4_DP, &
               "AMG nsmooth=1: should converge below "// &
               "1e-4 in 30 cycles")
  end subroutine test_amg_nsmooth_1

  !> @brief AMG converges on a nonsymmetric convection-diffusion matrix
  !!
  !! The matrix comes from a 1D convection-diffusion problem (D=1, v=5,
  !! Peclet number 5) with upwind differencing, so it is strongly
  !! nonsymmetric: a_{i,i-1} = -6 but a_{i-1,i} = -1. Aggregation uses a
  !! symmetric strength-of-connection measure, (|a_{ij}|+|a_{ji}|)/2 = 3.5,
  !! so aggregates stay balanced and do not pick up a direction from the
  !! advection term. Verifies convergence below 1e-4 in 30 adaptive
  !! V-cycles.
  !<
  subroutine test_amg_nonsymmetric_soc(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: ia(N_TEST + 1)
    integer(I4B), allocatable :: ja(:)
    real(DP), allocatable :: a(:)
    real(DP) :: b(N_TEST), x_exact(N_TEST)
    integer(I4B) :: nja

    call build_convdiff1d(N_TEST, ia, ja, a, b, x_exact)
    nja = size(ja)

    call check(error, &
               amg_vcycle_residual(ia, ja, a, omega=0.0_DP, &
                                   nsmooth=4, n_iter=30) < 1.0e-4_DP, &
               "AMG nonsymmetric (Pe=5): should converge "// &
               "below 1e-4 in 30 cycles with symmetric SOC")
  end subroutine test_amg_nonsymmetric_soc

  !> @brief AMG with a positive aggregation strength threshold still converges
  !!
  !! The 1D heterogeneous Poisson matrix has 100:1 conductance contrast between
  !! adjacent interfaces (alternating K=10 and K=0.1).  With sthresh=0.25, only
  !! connections whose combined strength is at least 25 percent of the
  !! strongest connection at that node are eligible for aggregation.  This
  !! filters the weak K=0.1 connections, producing layer-aligned aggregates.
  !! Verifies that convergence still occurs within 30 V-cycles.
  !<
  subroutine test_amg_sthresh_convergence(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: ia(N_TEST + 1)
    integer(I4B), allocatable :: ja(:)
    real(DP), allocatable :: a(:)
    real(DP) :: b(N_TEST), x_exact(N_TEST)
    integer(I4B) :: nja

    call build_poisson1d(N_TEST, ia, ja, a, b, x_exact)
    nja = size(ja)

    call check(error, &
               amg_vcycle_residual(ia, ja, a, omega=0.0_DP, &
                                   nsmooth=4, n_iter=30, &
                                   sthresh=0.25_DP) < 1.0e-4_DP, &
               "AMG sthresh=0.25 on heterogeneous K: should converge "// &
               "below 1e-4 in 30 cycles")
  end subroutine test_amg_sthresh_convergence

  !> @brief AMG with a positive strength threshold still builds a coarse
  !! hierarchy
  !!
  !! With sthresh=0.25, the weakest connections in the heterogeneous Poisson
  !! matrix are excluded from aggregation decisions.  The hierarchy must still
  !! coarsen to at least two levels and the coarsest level must be strictly
  !! smaller than the fine level, confirming that the threshold does not
  !! prevent coarsening entirely.
  !<
  subroutine test_amg_sthresh_hierarchy(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: ia(N_TEST + 1)
    integer(I4B), allocatable :: ja(:)
    real(DP), allocatable :: a(:)
    real(DP) :: b(N_TEST), x_exact(N_TEST)
    integer(I4B) :: nja
    type(ImsAmgDataType) :: amg

    call build_poisson1d(N_TEST, ia, ja, a, b, x_exact)
    nja = size(ja)

    call amg%setup(N_TEST, nja, ia, ja, a, max_levels=5, nsmooth=2, &
                   omega=1.0_DP, sthresh=0.25_DP, &
                   smoother_type=SMOOTHER_ILU0)

    call check(error, amg%nlevels >= 2, &
               "AMG sthresh=0.25: should still build "// &
               "at least 2 hierarchy levels")
    if (allocated(error)) then
      call amg%destroy()
      return
    end if

    call check(error, amg%levels(amg%nlevels - 1)%neq < N_TEST, &
               "AMG sthresh=0.25: coarsest level should have fewer unknowns "// &
               "than fine level")

    call amg%destroy()
  end subroutine test_amg_sthresh_hierarchy

  !> @brief AMG ILU0_ALL smoother (SMOOTHER_TYPE ILU0_ALL) converges
  !!
  !! Activates SMOOTHER_ILU0_ALL, which extends ILU(0) correction sweeps
  !! to every level of the AMG hierarchy rather than just the finest level.
  !! Verifies convergence below 1e-6 in 15 Richardson iterations on the
  !! 1D heterogeneous Poisson system.
  !<
  subroutine test_amg_ilu0all_converges(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B) :: ia(N_TEST + 1)
    integer(I4B), allocatable :: ja(:)
    real(DP), allocatable :: a(:)
    real(DP) :: b(N_TEST), x_exact(N_TEST)
    integer(I4B) :: nja

    call build_poisson1d(N_TEST, ia, ja, a, b, x_exact)
    nja = size(ja)

    call check(error, &
               amg_vcycle_residual(ia, ja, a, omega=0.0_DP, &
                                   nsmooth=2, n_iter=15, &
                                   smoother_type_opt=SMOOTHER_ILU0_ALL) &
               < 1.0e-6_DP, &
               "AMG ILU0_ALL smoother: relative residual "// &
               "should fall below 1e-6 in 15 cycles")
  end subroutine test_amg_ilu0all_converges

end module TestImsLinearMethods
