module TestImsLinear
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       to_string, unittest_type
  use ConvergenceSummaryModule, only: ConvergenceSummaryType
  use ImsLinearSettingsModule, only: IPC_ILU0
  use IMSLinearBaseModule, only: ims_base_cg, ims_base_bcgs, &
                                 ims_base_scale, ims_base_cnvg_scale, &
                                 ims_base_pccrs, ims_base_pcilu0, &
                                 ims_base_residual, ims_base_calc_order
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
                new_unittest("cg_diagonal_rcm", test_cg_diagonal_rcm), &
                new_unittest("bcgs_l2norm_rcm", test_bcgs_l2norm_rcm) &
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
  subroutine check_one_iteration(error, ilinmeth, iscl, iord)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), intent(in) :: ilinmeth
    integer(I4B), intent(in) :: iscl
    integer(I4B), intent(in) :: iord
    integer(I4B), dimension(NEQ + 1) :: ia, iao, iapc
    integer(I4B), dimension(NJA) :: ja, jao, japc
    integer(I4B), dimension(NEQ) :: iw, lorder, iorder, idum
    integer(I4B), dimension(1) :: jlu
    integer(I4B), dimension(2) :: convmodstart
    integer(I4B), dimension(:), pointer :: bounds
    real(DP), dimension(NJA) :: a, a0, ao, apc
    real(DP), dimension(NEQ) :: b, b0, x, d, p, q, z, t, v, dhat, phat, qhat
    real(DP), dimension(NEQ) :: w, dscale, dscale2, rscale, hscale
    real(DP), dimension(NEQ) :: x0, xm, dx, r
    character(len=31), dimension(1) :: caccel
    type(ConvergenceSummaryType), pointer :: summary
    integer(I4B) :: n, icnvg, innerit, ipcflag, nmoved
    real(DP) :: l2norm0, dvexp, rexp

    call setup_grid(ia, ja, a0, b0)
    do n = 1, NJA
      a(n) = a0(n)
    end do
    do n = 1, NEQ
      b(n) = b0(n)
      x(n) = DZERO
    end do

    ! scale and reorder as imslinear_ap does, then set the factors and the
    ! preconditioner
    call ims_base_scale(0, iscl, NEQ, NJA, ia, ja, a, x, b, dscale, dscale2)
    do n = 1, NEQ
      x0(n) = x(n) * dscale2(n)
    end do
    if (iord /= 0) then
      call ims_base_calc_order(iord, NEQ, NJA, ia, ja, lorder, iorder)
      nmoved = 0
      do n = 1, NEQ
        if (lorder(n) /= n) nmoved = nmoved + 1
      end do
      call check(error, nmoved > 0, 'reordering is the identity')
      if (allocated(error)) return
      call dperm(NEQ, a, ja, ia, ao, jao, iao, lorder, idum, 1)
      call dvperm(NEQ, x, lorder)
      call dvperm(NEQ, b, lorder)
    else
      do n = 1, NEQ + 1
        iao(n) = ia(n)
      end do
      do n = 1, NJA
        jao(n) = ja(n)
        ao(n) = a(n)
      end do
    end if
    call ims_base_cnvg_scale(NEQ, iord, dscale, dscale2, lorder, rscale, &
                             hscale)
    call ims_base_pccrs(NEQ, NJA, iao, jao, iapc, japc)
    ipcflag = 0
    call ims_base_pcilu0(NJA, NEQ, ao, iao, jao, apc, iapc, japc, iw, w, &
                         DZERO, ipcflag, DZERO)
    call ims_base_residual(NEQ, NJA, x, b, d, ao, iao, jao)
    l2norm0 = DONE

    allocate (summary, bounds(2))
    bounds(1) = 1
    bounds(2) = NEQ
    ! a memory path for each case, since each init allocates new entries
    call summary%init(1, bounds, 'IMSTEST'//to_string(ilinmeth)// &
                      to_string(iscl)//to_string(iord))
    convmodstart(1) = 1
    convmodstart(2) = NEQ + 1
    icnvg = 0
    jlu(1) = 0
    if (ilinmeth == 1) then
      call ims_base_cg(icnvg, 1, innerit, NEQ, NJA, NEQ, NJA, IPC_ILU0, 0, 0, &
                       iscl, rscale, hscale, 1.0d-30, 1.0d-30, l2norm0, &
                       DONE, iao, jao, ao, iapc, japc, apc, x, b, d, p, q, z, &
                       1, iw, jlu, 1, 1, convmodstart, caccel, summary)
    else
      call ims_base_bcgs(icnvg, 1, innerit, NEQ, NJA, NEQ, NJA, IPC_ILU0, 0, &
                         0, iscl, rscale, hscale, 1.0d-30, 1.0d-30, l2norm0, &
                         DONE, iao, jao, ao, iapc, japc, apc, x, b, d, p, q, &
                         t, v, dhat, phat, qhat, 1, iw, jlu, 1, 1, &
                         convmodstart, caccel, summary)
    end if

    ! head change and residual of the unscaled system in model order
    if (iord /= 0) call dvperm(NEQ, x, iorder)
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
    if (.not. allocated(error)) &
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
    call check_one_iteration(error, 1, 1, 0)
  end subroutine test_cg_diagonal

  subroutine test_cg_l2norm(error)
    type(error_type), allocatable, intent(out) :: error
    call check_one_iteration(error, 1, 2, 0)
  end subroutine test_cg_l2norm

  subroutine test_bcgs_diagonal(error)
    type(error_type), allocatable, intent(out) :: error
    call check_one_iteration(error, 2, 1, 0)
  end subroutine test_bcgs_diagonal

  subroutine test_bcgs_l2norm(error)
    type(error_type), allocatable, intent(out) :: error
    call check_one_iteration(error, 2, 2, 0)
  end subroutine test_bcgs_l2norm

  subroutine test_cg_diagonal_rcm(error)
    type(error_type), allocatable, intent(out) :: error
    call check_one_iteration(error, 1, 1, 1)
  end subroutine test_cg_diagonal_rcm

  subroutine test_bcgs_l2norm_rcm(error)
    type(error_type), allocatable, intent(out) :: error
    call check_one_iteration(error, 2, 2, 1)
  end subroutine test_bcgs_l2norm_rcm

end module TestImsLinear
