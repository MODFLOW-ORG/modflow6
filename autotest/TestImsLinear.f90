module TestImsLinear
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       to_string, unittest_type
  use ConvergenceSummaryModule, only: ConvergenceSummaryType
  use ImsLinearSettingsModule, only: IPC_ILU0, IPC_MILU0
  use IMSLinearBaseModule, only: ims_base_cg, ims_base_bcgs, &
                                 ims_base_scale, ims_base_cnvg_scale, &
                                 ims_base_pccrs, ims_base_pcilu0, &
                                 ims_base_residual, ims_base_pcu, &
                                 ims_base_ilu0a
  implicit none
  private
  public :: collect_imslinear

  ! 4 by 4 five-point grid
  integer(I4B), parameter :: NC = 4
  integer(I4B), parameter :: NEQ = NC * NC
  integer(I4B), parameter :: NJA = NEQ + 4 * NC * (NC - 1)
  real(DP), parameter :: RTOL = 1.0d-8

contains

  subroutine collect_imslinear(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("cnvg_scale_reorder", test_cnvg_scale_reorder), &
                new_unittest("cg_diagonal", test_cg_diagonal), &
                new_unittest("cg_l2norm", test_cg_l2norm), &
                new_unittest("bcgs_diagonal", test_bcgs_diagonal), &
                new_unittest("bcgs_l2norm", test_bcgs_l2norm), &
                new_unittest("pcilu0_bad_pivot", test_pcilu0_bad_pivot), &
                new_unittest("pcu_bad_pivot", test_pcu_bad_pivot) &
                ]
  end subroutine collect_imslinear

  !> @brief Build a grid with conductances from 1e3 to 1e5, diagonal first
  !<
  subroutine setup_grid(ia, ja, a, b)
    integer(I4B), dimension(NEQ + 1), intent(out) :: ia
    integer(I4B), dimension(NJA), intent(out) :: ja
    real(DP), dimension(NJA), intent(out) :: a
    real(DP), dimension(NEQ), intent(out) :: b
    integer(I4B) :: i, j, n, m, k, idiag
    real(DP) :: c

    k = 0
    do n = 1, NEQ
      i = (n - 1) / NC + 1
      j = n - (i - 1) * NC
      k = k + 1
      ia(n) = k
      idiag = k
      ja(k) = n
      a(k) = DONE
      do m = 1, NEQ
        if (m == n) cycle
        if (.not. ((abs(m - n) == NC) .or. &
                   (abs(m - n) == 1 .and. (m - 1) / NC == (n - 1) / NC))) &
          cycle
        c = 10.0d0**(3.0d0 + 2.0d0 * modulo(3 * min(m, n) + max(m, n), 7) &
                     / 6.0d0)
        k = k + 1
        ja(k) = m
        a(k) = -c
        a(idiag) = a(idiag) + c
      end do
      b(n) = 1.0d3 * real(n, DP)
    end do
    ia(NEQ + 1) = k + 1
  end subroutine setup_grid

  !> @brief Return the largest absolute value of a vector
  !<
  function absmax(v) result(vmax)
    real(DP), dimension(:), intent(in) :: v
    real(DP) :: vmax
    integer(I4B) :: n

    vmax = DZERO
    do n = 1, size(v)
      vmax = max(vmax, abs(v(n)))
    end do
  end function absmax

  !> @brief Run one inner iteration and compare the reported head change
  !! and residual with model unit values from the unscaled system
  !<
  subroutine check_one_iteration(error, ilinmeth, iscl)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), intent(in) :: ilinmeth
    integer(I4B), intent(in) :: iscl
    integer(I4B), dimension(NEQ + 1) :: ia, iapc
    integer(I4B), dimension(NJA) :: ja, japc
    integer(I4B), dimension(NEQ) :: iw
    integer(I4B), dimension(1) :: jlu, lorder
    integer(I4B), dimension(2) :: convmodstart
    integer(I4B), dimension(:), pointer :: bounds
    real(DP), dimension(NJA) :: a, a0, apc
    real(DP), dimension(NEQ) :: b, b0, x, d, p, q, z, t, v, dhat, phat, qhat
    real(DP), dimension(NEQ) :: w, dscale, dscale2, rscale, hscale
    real(DP), dimension(NEQ) :: x0, xm, dx, r
    character(len=31), dimension(1) :: caccel
    type(ConvergenceSummaryType), pointer :: summary
    integer(I4B) :: n, icnvg, innerit, ipcflag
    real(DP) :: l2norm0, dvexp, rexp

    call setup_grid(ia, ja, a0, b0)
    do n = 1, NJA
      a(n) = a0(n)
    end do
    do n = 1, NEQ
      b(n) = b0(n)
      x(n) = DZERO
    end do

    ! scale, then set the factors and the preconditioner
    call ims_base_scale(0, iscl, NEQ, NJA, ia, ja, a, x, b, dscale, dscale2)
    lorder(1) = 0
    call ims_base_cnvg_scale(NEQ, 0, dscale, dscale2, lorder, rscale, hscale)
    call ims_base_pccrs(NEQ, NJA, ia, ja, iapc, japc)
    ipcflag = 0
    call ims_base_pcilu0(NJA, NEQ, a, ia, ja, apc, iapc, japc, iw, w, &
                         DZERO, ipcflag, DZERO)
    call ims_base_residual(NEQ, NJA, x, b, d, a, ia, ja)
    l2norm0 = DONE
    do n = 1, NEQ
      x0(n) = x(n) * dscale2(n)
    end do

    allocate (summary, bounds(2))
    bounds(1) = 1
    bounds(2) = NEQ
    ! a memory path for each case, since each init allocates new entries
    call summary%init(1, bounds, &
                      'IMSTEST'//to_string(ilinmeth)//to_string(iscl))
    convmodstart(1) = 1
    convmodstart(2) = NEQ + 1
    icnvg = 0
    jlu(1) = 0
    if (ilinmeth == 1) then
      call ims_base_cg(icnvg, 1, innerit, NEQ, NJA, NEQ, NJA, IPC_ILU0, 0, 0, &
                       iscl, rscale, hscale, 1.0d-30, 1.0d-30, l2norm0, &
                       DONE, ia, ja, a, iapc, japc, apc, x, b, d, p, q, z, &
                       1, iw, jlu, 1, 1, convmodstart, caccel, summary)
    else
      call ims_base_bcgs(icnvg, 1, innerit, NEQ, NJA, NEQ, NJA, IPC_ILU0, 0, &
                         0, iscl, rscale, hscale, 1.0d-30, 1.0d-30, l2norm0, &
                         DONE, ia, ja, a, iapc, japc, apc, x, b, d, p, q, &
                         t, v, dhat, phat, qhat, 1, iw, jlu, 1, 1, &
                         convmodstart, caccel, summary)
    end if

    ! head change and residual of the unscaled system
    do n = 1, NEQ
      xm(n) = x(n) * dscale2(n)
      dx(n) = xm(n) - x0(n)
    end do
    call ims_base_residual(NEQ, NJA, xm, b0, r, a0, ia, ja)
    dvexp = absmax(dx)
    rexp = absmax(r)

    call check(error, abs(abs(summary%dvmax(1)) - dvexp) <= RTOL * dvexp, &
               'head change '//to_string(abs(summary%dvmax(1)))// &
               ' expected '//to_string(dvexp))
    if (allocated(error)) return
    call check(error, abs(abs(summary%rmax(1)) - rexp) <= RTOL * rexp, &
               'residual '//to_string(abs(summary%rmax(1)))// &
               ' expected '//to_string(rexp))

    call summary%destroy()
    deallocate (summary, bounds)
  end subroutine check_one_iteration

  !> @brief Factors follow the reordering of the system
  !<
  subroutine test_cnvg_scale_reorder(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: N4 = 4
    integer(I4B), dimension(N4) :: lorder, lorder0
    real(DP), dimension(N4) :: dscale, dscale2, rscale, hscale, node
    integer(I4B) :: n, k

    lorder0 = [3, 1, 4, 2]
    do n = 1, N4
      lorder(n) = lorder0(n)
      dscale(n) = real(n, DP)
      dscale2(n) = 1.0d1 * real(n, DP)
      node(n) = real(n, DP)
    end do
    call ims_base_cnvg_scale(N4, 1, dscale, dscale2, lorder, rscale, hscale)

    ! the node at each position of the reordered system
    call dvperm(N4, node, lorder)
    do n = 1, N4
      k = nint(node(n))
      call check(error, rscale(n) == DONE / dscale(k), &
                 'rscale at '//to_string(n))
      if (allocated(error)) return
      call check(error, hscale(n) == dscale2(k), 'hscale at '//to_string(n))
      if (allocated(error)) return
      call check(error, lorder(n) == lorder0(n), 'lorder at '//to_string(n))
      if (allocated(error)) return
    end do
  end subroutine test_cnvg_scale_reorder

  subroutine test_cg_diagonal(error)
    type(error_type), allocatable, intent(out) :: error
    call check_one_iteration(error, 1, 1)
  end subroutine test_cg_diagonal

  subroutine test_cg_l2norm(error)
    type(error_type), allocatable, intent(out) :: error
    call check_one_iteration(error, 1, 2)
  end subroutine test_cg_l2norm

  subroutine test_bcgs_diagonal(error)
    type(error_type), allocatable, intent(out) :: error
    call check_one_iteration(error, 2, 1)
  end subroutine test_bcgs_diagonal

  subroutine test_bcgs_l2norm(error)
    type(error_type), allocatable, intent(out) :: error
    call check_one_iteration(error, 2, 2)
  end subroutine test_bcgs_l2norm

  !> @brief Build a 3 by 3 matrix whose MILU0 pivot changes sign
  !<
  subroutine setup_bad_pivot(ia, ja, a, iapc, japc)
    integer(I4B), dimension(4), intent(out) :: ia, iapc
    integer(I4B), dimension(7), intent(out) :: ja, japc
    real(DP), dimension(7), intent(out) :: a
    real(DP), parameter :: c = 0.95d0

    ia = [1, 4, 6, 8]
    ja = [1, 2, 3, 2, 1, 3, 1]
    a = [DONE, c, c, DONE, c, DONE, c]
    call ims_base_pccrs(3, 7, ia, ja, iapc, japc)
  end subroutine setup_bad_pivot

  !> @brief A bad pivot is reported to the caller
  !<
  subroutine test_pcilu0_bad_pivot(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), dimension(4) :: ia, iapc
    integer(I4B), dimension(7) :: ja, japc
    integer(I4B), dimension(3) :: iw
    real(DP), dimension(7) :: a, apc
    real(DP), dimension(3) :: w
    integer(I4B) :: ipcflag

    call setup_bad_pivot(ia, ja, a, iapc, japc)
    ipcflag = 0
    call ims_base_pcilu0(7, 3, a, ia, ja, apc, iapc, japc, iw, w, &
                         0.97d0, ipcflag, DZERO)
    call check(error, ipcflag == 1, 'ipcflag '//to_string(ipcflag))
  end subroutine test_pcilu0_bad_pivot

  !> @brief The retries end with a complete, finite factorization
  !<
  subroutine test_pcu_bad_pivot(error)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), dimension(4) :: ia, iapc
    integer(I4B), dimension(7) :: ja, japc
    integer(I4B), dimension(3) :: iw
    integer(I4B), dimension(1) :: jlu, jw
    real(DP), dimension(7) :: a, apc
    real(DP), dimension(3) :: w, r, d
    real(DP), dimension(1) :: wlu
    real(DP), parameter :: SENTINEL = -9.0d30
    integer(I4B) :: iout, n

    call setup_bad_pivot(ia, ja, a, iapc, japc)
    do n = 1, 7
      apc(n) = SENTINEL
    end do
    open (newunit=iout, status='scratch')
    call ims_base_pcu(iout, 7, 3, 3, 7, IPC_MILU0, 0.97d0, a, ia, ja, &
                      apc, iapc, japc, iw, w, 0, DZERO, 1, 1, 1, jlu, jw, wlu)
    close (iout)
    do n = 1, 7
      call check(error, apc(n) /= SENTINEL, &
                 'apc not set at '//to_string(n))
      if (allocated(error)) return
    end do

    do n = 1, 3
      r(n) = DONE
    end do
    call ims_base_ilu0a(7, 3, apc, iapc, japc, r, d)
    do n = 1, 3
      call check(error, ieee_is_finite(d(n)), 'd not finite at '//to_string(n))
      if (allocated(error)) return
    end do
  end subroutine test_pcu_bad_pivot

end module TestImsLinear
