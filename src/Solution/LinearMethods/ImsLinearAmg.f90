!> @brief Algebraic Multi-Grid (AMG) preconditioner for the IMS linear solver
!!
!! This module implements a simple aggregation-based AMG preconditioner
!! used as an alternative to the ILU-based preconditioners in the IMS solver.
!<
module ImsLinearAmgModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DHALF, DONE, DHUNDRED, DEP20, LINELENGTH
  use SimModule, only: store_error
  use IMSLinearMisc, only: ims_base_pccrs, ims_base_pcilu0, ims_base_ilu0a
  use ImsLinearSettingsModule, only: SMOOTHER_ILU0, SMOOTHER_ILU0_ALL
  use ProfilerModule, only: g_prof, LEN_SECTION_TITLE

  implicit none
  private

  public :: ImsAmgDataType, ims_amg_setup, ims_amg_apply, &
            ims_amg_da, ims_amg_update, ims_amg_summary

  !> maximum number of smoother sweeps used for the coarsest-level solve
  integer(I4B), parameter :: COARSE_SOLVE_ITERS = 50

  !> per-level decay applied to the strength threshold
  !!
  !! Galerkin coarsening densifies the operator, so a fixed threshold measured
  !! against the largest entry in a row rejects a growing share of connections
  !! and stalls the coarsening. Relaxing the threshold as levels coarsen keeps
  !! the aggregation productive.
  real(DP), parameter :: STHRESH_DECAY = DHALF

  !> @brief Storage for one level in the AMG hierarchy
  !<
  type :: AmgLevelType
    integer(I4B) :: neq = 0 !< number of equations at this level
    integer(I4B) :: nja = 0 !< number of non-zeros at this level
    integer(I4B) :: nc = 0 !< number of coarse-grid aggregates
    integer(I4B), allocatable :: ia(:) !< CSR row pointers
    integer(I4B), allocatable :: ja(:) !< CSR column pointers
    integer(I4B), allocatable :: agg(:) !< aggregate assignment (fine to coarse)
    real(DP), allocatable :: a(:) !< coefficient values
    real(DP), allocatable :: diag(:) !< inverse of diagonal entries
    real(DP), allocatable :: r(:) !< right-hand side work vector
    real(DP), allocatable :: e(:) !< error/solution work vector
    ! ILU(0) smoother data (finest level only)
    integer(I4B), allocatable :: iapc(:) !< ILU0 CSR row pointers (size neq+1)
    integer(I4B), allocatable :: japc(:) !< ILU0 CSR column pointers (size nja)
    real(DP), allocatable :: apc(:) !< ILU0 factorization values (size nja)
    real(DP), allocatable :: work(:) !< smoother work vector (size neq)
    ! aggregation diagnostics (set by build_aggregation)
    real(DP) :: sthresh = DZERO !< strength threshold applied at this level
    integer(I4B) :: nsingle = 0 !< aggregates holding one cell
    integer(I4B) :: ncand = 0 !< unaggregated neighbours considered
    integer(I4B) :: nweak = 0 !< candidates rejected by the strength threshold
    ! per-level profiler handle/title (used only under PROFILE_OPTION DETAIL)
    integer(I4B) :: itmr_smooth = -1 !< "AMG smooth L<n>" timer handle
    character(len=LEN_SECTION_TITLE) :: tmr_title = "" !< section title for this level
  end type AmgLevelType

  !> @brief AMG preconditioner data type
  !<
  type :: ImsAmgDataType
    type(AmgLevelType), allocatable :: levels(:) !< AMG level hierarchy
    integer(I4B) :: nlevels = 0 !< number of active levels
    integer(I4B) :: nsmooth = 2 !< number of smoother iterations
    real(DP) :: omega = DONE !< relaxation factor (<=0 means adaptive)
    integer(I4B) :: smoother_type = SMOOTHER_ILU0_ALL !< smoother type enum
    integer(I4B) :: overflow_count = 0 !< outer iterations with overflow-scale residuals
    logical :: overflow_this_outer = .false. !< overflow seen in current outer iteration
    ! profiler section handles (-1 until the first g_prof%start creates them)
    integer(I4B) :: itmr_setup = -1 !< "AMG setup" timer handle
    integer(I4B) :: itmr_apply = -1 !< "AMG apply" timer handle
    integer(I4B) :: itmr_cvals = -1 !< "AMG coarse values" timer handle
    integer(I4B) :: itmr_ilu = -1 !< "AMG ILU0 factor" timer handle
  contains
    procedure :: setup => amg_setup
    procedure :: apply => amg_apply
    procedure :: destroy => amg_destroy
  end type ImsAmgDataType

contains

  ! ---------------------------------------------------------------------------
  ! Plain subroutine interface (ILU-style, used by ImsLinear and ImsLinearBase)
  ! ---------------------------------------------------------------------------

  !> @brief Set up the AMG preconditioner hierarchy (plain subroutine interface)
  !!
  !! Builds the AMG hierarchy from the fine-level matrix using greedy pairwise
  !! aggregation. The actual number of levels may be fewer than max_levels
  !! if the coarse grid becomes too small.
  !<
  subroutine ims_amg_setup(amg, neq, nja, ia, ja, a, max_levels, nsmooth, &
                           omega, sthresh, smoother_type)
    type(ImsAmgDataType), intent(inout) :: amg !< AMG preconditioner
    integer(I4B), intent(in) :: neq !< number of equations (fine level)
    integer(I4B), intent(in) :: nja !< number of non-zeros (fine level)
    integer(I4B), dimension(neq + 1), intent(in) :: ia !< CSR row pointers
    integer(I4B), dimension(nja), intent(in) :: ja !< CSR column pointers
    real(DP), dimension(nja), intent(in) :: a !< coefficient matrix
    integer(I4B), intent(in) :: max_levels !< maximum number of levels
    integer(I4B), intent(in) :: nsmooth !< number of smoother iterations
    real(DP), intent(in) :: omega !< relaxation factor
    real(DP), intent(in) :: sthresh !< strength-of-connection threshold
    integer(I4B), intent(in) :: smoother_type !< SMOOTHER_ILU0 or SMOOTHER_ILU0_ALL
    ! -- local
    integer(I4B) :: nlevs
    integer(I4B) :: l
    integer(I4B) :: nc
    integer(I4B) :: n

    ! If hierarchy already built with matching sparsity, only update values.
    if (amg%nlevels > 0 .and. allocated(amg%levels)) then
      if (amg%levels(0)%neq == neq .and. amg%levels(0)%nja == nja) then
        call ims_amg_update(amg, nja, a)
        amg%nsmooth = nsmooth
        amg%omega = omega
        amg%smoother_type = smoother_type
        return
      end if
    end if

    call ims_amg_da(amg)

    nlevs = max(1, max_levels)
    allocate (amg%levels(0:nlevs - 1))
    amg%nlevels = 1

    amg%levels(0)%neq = neq
    amg%levels(0)%nja = nja
    allocate (amg%levels(0)%ia(neq + 1))
    allocate (amg%levels(0)%ja(nja))
    allocate (amg%levels(0)%a(nja))
    allocate (amg%levels(0)%agg(neq))
    do n = 1, neq + 1
      amg%levels(0)%ia(n) = ia(n)
    end do
    do n = 1, nja
      amg%levels(0)%ja(n) = ja(n)
      amg%levels(0)%a(n) = a(n)
    end do
    do n = 1, neq
      amg%levels(0)%agg(n) = 0
    end do
    amg%levels(0)%nc = 0

    do l = 0, nlevs - 2
      if (amg%levels(l)%neq <= 2) exit
      call build_aggregation(amg%levels(l), sthresh * STHRESH_DECAY**l)
      nc = amg%levels(l)%nc
      if (nc <= 1) exit
      call build_coarse_matrix(amg%levels(l), amg%levels(l + 1))
      allocate (amg%levels(l + 1)%agg(nc))
      do n = 1, nc
        amg%levels(l + 1)%agg(n) = 0
      end do
      amg%levels(l + 1)%nc = 0
      amg%nlevels = l + 2
    end do

    do l = 0, amg%nlevels - 1
      call alloc_work(amg%levels(l))
      write (amg%levels(l)%tmr_title, '(a,i0)') "AMG smooth L", l
    end do

    amg%nsmooth = nsmooth
    amg%omega = omega
    amg%smoother_type = smoother_type

    if (smoother_type == SMOOTHER_ILU0) then
      call build_ilu0(amg%levels(0), .false.)
    else
      do l = 0, amg%nlevels - 1
        call build_ilu0(amg%levels(l), l > 0)
      end do
    end if

  end subroutine ims_amg_setup

  !> @brief Write the AMG hierarchy to the listing file
  !!
  !! Reports per-level size and coarsening statistics along with the grid and
  !! operator complexities, which together indicate whether the aggregation is
  !! coarsening effectively and what one cycle costs relative to one fine-level
  !! matrix-vector product.
  !<
  subroutine ims_amg_summary(amg, iout, sthresh)
    ! -- modules
    use ConstantsModule, only: TABCENTER, TABRIGHT
    use TableModule, only: TableType, table_cr
    ! -- dummy
    type(ImsAmgDataType), intent(in) :: amg !< AMG preconditioner
    integer(I4B), intent(in) :: iout !< listing file unit
    real(DP), intent(in) :: sthresh !< strength-of-connection threshold
    ! -- local
    type(TableType), pointer :: tbl => null()
    character(len=LINELENGTH) :: tag
    integer(I4B) :: l
    integer(I4B) :: neq_tot, nja_tot
    real(DP) :: ratio, fsingle, fweak
    ! -- formats
02010 format(/, ' GRID COMPLEXITY (sum rows / fine rows)         =', F10.3, /, &
            ' OPERATOR COMPLEXITY (sum nonzeros / fine nnz)  =', F10.3, /, &
            ' STRENGTH THRESHOLD (finest level)              =', E15.5, //)

    if (iout <= 0) return
    if (amg%nlevels <= 0 .or. .not. allocated(amg%levels)) return

    tag = 'AMG HIERARCHY'
    call table_cr(tbl, 'AMG', tag)
    call tbl%table_df(amg%nlevels, 8, iout)
    tag = 'LEVEL'
    call tbl%initialize_column(tag, 7, alignment=TABCENTER)
    tag = 'ROWS'
    call tbl%initialize_column(tag, 12, alignment=TABRIGHT)
    tag = 'NONZEROS'
    call tbl%initialize_column(tag, 14, alignment=TABRIGHT)
    tag = 'NONZEROS PER ROW'
    call tbl%initialize_column(tag, 12, alignment=TABRIGHT)
    tag = 'STRENGTH THRESHOLD'
    call tbl%initialize_column(tag, 12, alignment=TABRIGHT)
    tag = 'COARSENING RATIO'
    call tbl%initialize_column(tag, 12, alignment=TABRIGHT)
    tag = 'SINGLETON AGGREGATES, IN PERCENT'
    call tbl%initialize_column(tag, 12, alignment=TABRIGHT)
    tag = 'WEAK CONNECTIONS, IN PERCENT'
    call tbl%initialize_column(tag, 12, alignment=TABRIGHT)

    neq_tot = 0
    nja_tot = 0
    do l = 0, amg%nlevels - 1
      neq_tot = neq_tot + amg%levels(l)%neq
      nja_tot = nja_tot + amg%levels(l)%nja
      ! coarsening ratio and diagnostics are undefined on the coarsest level,
      ! which is solved rather than aggregated
      if (l < amg%nlevels - 1 .and. amg%levels(l)%nc > 0) then
        ratio = real(amg%levels(l)%neq, DP) / real(amg%levels(l)%nc, DP)
        fsingle = DHUNDRED * real(amg%levels(l)%nsingle, DP) / &
                  real(amg%levels(l)%nc, DP)
      else
        ratio = DZERO
        fsingle = DZERO
      end if
      if (amg%levels(l)%ncand > 0) then
        fweak = DHUNDRED * real(amg%levels(l)%nweak, DP) / &
                real(amg%levels(l)%ncand, DP)
      else
        fweak = DZERO
      end if
      call tbl%add_term(l)
      call tbl%add_term(amg%levels(l)%neq)
      call tbl%add_term(amg%levels(l)%nja)
      call tbl%add_term(real(amg%levels(l)%nja, DP) / &
                        real(max(1, amg%levels(l)%neq), DP))
      call tbl%add_term(amg%levels(l)%sthresh)
      call tbl%add_term(ratio)
      call tbl%add_term(fsingle)
      call tbl%add_term(fweak)
    end do
    call tbl%table_da()
    deallocate (tbl)
    nullify (tbl)

    write (iout, 2010) &
      real(neq_tot, DP) / real(max(1, amg%levels(0)%neq), DP), &
      real(nja_tot, DP) / real(max(1, amg%levels(0)%nja), DP), &
      sthresh

  end subroutine ims_amg_summary

  !> @brief Apply the AMG preconditioner: z = M^{-1} r
  !! (plain subroutine interface)
  !<
  subroutine ims_amg_apply(amg, neq, r_in, z_out)
    type(ImsAmgDataType), intent(inout) :: amg !< AMG preconditioner
    integer(I4B), intent(in) :: neq !< number of equations
    real(DP), dimension(neq), intent(in) :: r_in !< input residual vector
    real(DP), dimension(neq), intent(inout) :: z_out !< output precond. vector
    ! -- local
    integer(I4B) :: n, k
    real(DP) :: rnorm
    real(DP) :: znorm
    real(DP) :: rz_dot, zaz_dot
    real(DP) :: alpha
    real(DP) :: az_val

    if (amg%nlevels == 0) then
      do n = 1, neq
        z_out(n) = r_in(n)
      end do
      return
    end if

    ! Guard: an overflow-scale residual would corrupt the ILU factorization,
    ! so return a zero update. amg_note_overflow counts these per outer
    ! iteration and stops the run if they persist.
    rnorm = DZERO
    do n = 1, neq
      if (abs(r_in(n)) > rnorm) rnorm = abs(r_in(n))
    end do
    if (rnorm > DEP20) then
      call amg_note_overflow(amg, 'input residual norm', rnorm)
      do n = 1, neq
        z_out(n) = DZERO
      end do
      return
    end if

    do n = 1, neq
      amg%levels(0)%r(n) = r_in(n)
      amg%levels(0)%e(n) = DZERO
    end do

    call amg_vcycle(amg%levels, amg%nlevels, amg%nsmooth)

    znorm = DZERO
    do n = 1, neq
      if (abs(amg%levels(0)%e(n)) > znorm) znorm = abs(amg%levels(0)%e(n))
    end do
    if (znorm > DEP20) then
      call amg_note_overflow(amg, 'output correction norm', znorm)
      do n = 1, neq
        z_out(n) = DZERO
      end do
      return
    end if

    if (amg%omega <= DZERO) then
      ! -- Choose alpha to make the error as small as possible:
      ! -- alpha = (r, z) / (z, Az). Only correct when A is symmetric
      ! -- positive definite (SPD), as it is for the CG solver.
      rz_dot = DZERO
      zaz_dot = DZERO
      do n = 1, neq
        rz_dot = rz_dot + r_in(n) * amg%levels(0)%e(n)
      end do
      do n = 1, neq
        az_val = DZERO
        do k = amg%levels(0)%ia(n), amg%levels(0)%ia(n + 1) - 1
          az_val = az_val + amg%levels(0)%a(k) * &
                   amg%levels(0)%e(amg%levels(0)%ja(k))
        end do
        zaz_dot = zaz_dot + amg%levels(0)%e(n) * az_val
      end do
      if (zaz_dot > DZERO) then
        alpha = rz_dot / zaz_dot
        if (alpha > DONE) alpha = DONE
      else
        alpha = DONE
      end if
      do n = 1, neq
        z_out(n) = alpha * amg%levels(0)%e(n)
      end do
    else
      do n = 1, neq
        z_out(n) = amg%omega * amg%levels(0)%e(n)
      end do
    end if

  end subroutine ims_amg_apply

  !> @brief Deallocate the AMG preconditioner hierarchy
  !! (plain subroutine interface)
  !<
  subroutine ims_amg_da(amg)
    type(ImsAmgDataType), intent(inout) :: amg !< AMG preconditioner
    ! -- local
    integer(I4B) :: l

    if (allocated(amg%levels)) then
      do l = lbound(amg%levels, 1), ubound(amg%levels, 1)
        if (allocated(amg%levels(l)%ia)) deallocate (amg%levels(l)%ia)
        if (allocated(amg%levels(l)%ja)) deallocate (amg%levels(l)%ja)
        if (allocated(amg%levels(l)%a)) deallocate (amg%levels(l)%a)
        if (allocated(amg%levels(l)%agg)) deallocate (amg%levels(l)%agg)
        if (allocated(amg%levels(l)%diag)) deallocate (amg%levels(l)%diag)
        if (allocated(amg%levels(l)%r)) deallocate (amg%levels(l)%r)
        if (allocated(amg%levels(l)%e)) deallocate (amg%levels(l)%e)
        if (allocated(amg%levels(l)%iapc)) deallocate (amg%levels(l)%iapc)
        if (allocated(amg%levels(l)%japc)) deallocate (amg%levels(l)%japc)
        if (allocated(amg%levels(l)%apc)) deallocate (amg%levels(l)%apc)
        if (allocated(amg%levels(l)%work)) deallocate (amg%levels(l)%work)
      end do
      deallocate (amg%levels)
    end if
    amg%nlevels = 0

  end subroutine ims_amg_da

  ! ---------------------------------------------------------------------------
  ! Type-bound wrappers (delegate to the plain subroutines above)
  ! ---------------------------------------------------------------------------

  subroutine amg_setup(this, neq, nja, ia, ja, a, max_levels, nsmooth, &
                       omega, sthresh, smoother_type)
    class(ImsAmgDataType), intent(inout) :: this
    integer(I4B), intent(in) :: neq
    integer(I4B), intent(in) :: nja
    integer(I4B), dimension(neq + 1), intent(in) :: ia
    integer(I4B), dimension(nja), intent(in) :: ja
    real(DP), dimension(nja), intent(in) :: a
    integer(I4B), intent(in) :: max_levels
    integer(I4B), intent(in) :: nsmooth
    real(DP), intent(in) :: omega
    real(DP), intent(in) :: sthresh
    integer(I4B), intent(in) :: smoother_type
    call ims_amg_setup(this, neq, nja, ia, ja, a, max_levels, nsmooth, &
                       omega, sthresh, smoother_type)
  end subroutine amg_setup

  subroutine amg_apply(this, neq, r_in, z_out)
    class(ImsAmgDataType), intent(inout) :: this
    integer(I4B), intent(in) :: neq
    real(DP), dimension(neq), intent(in) :: r_in
    real(DP), dimension(neq), intent(inout) :: z_out
    call ims_amg_apply(this, neq, r_in, z_out)
  end subroutine amg_apply

  subroutine amg_destroy(this)
    class(ImsAmgDataType), intent(inout) :: this
    call ims_amg_da(this)
  end subroutine amg_destroy

  ! ---------------------------------------------------------------------------
  ! Private helper subroutines
  ! ---------------------------------------------------------------------------

  !> @brief Record an overflow-scale norm and stop the run if they persist
  !!
  !! Called when a residual or correction norm exceeds the overflow threshold.
  !! Counts at most one overflow per outer iteration and, once OVERFLOW_MAX of
  !! them have occurred, stops the run with guidance on likely causes.
  !<
  subroutine amg_note_overflow(amg, label, norm)
    type(ImsAmgDataType), intent(inout) :: amg !< AMG preconditioner
    character(len=*), intent(in) :: label !< name of the norm (for the message)
    real(DP), intent(in) :: norm !< the overflow-scale norm value
    ! -- local
    integer(I4B), parameter :: OVERFLOW_MAX = 5
    character(len=LINELENGTH) :: errmsg

    ! only count the first overflow seen in a given outer iteration
    if (amg%overflow_this_outer) return
    amg%overflow_this_outer = .true.
    amg%overflow_count = amg%overflow_count + 1
    if (amg%overflow_count < OVERFLOW_MAX) return

    write (errmsg, '(a,es12.4,a)') &
      'AMG preconditioner: '//trim(label)//' = ', &
      norm, ' exceeds the overflow threshold (1.0e+20).'
    call store_error( &
      trim(errmsg)// &
      ' The Newton iterate has diverged beyond recoverable'// &
      ' range. Possible causes and remedies:'// &
      ' (1) NO_PTC in the IMS OPTIONS block -- PTC is'// &
      ' required for steady-state stress periods in'// &
      ' strongly heterogeneous models; remove NO_PTC or'// &
      ' provide a near-converged initial condition.'// &
      ' (2) icelltype=1 on topmost-active (water table)'// &
      ' cells -- set icelltype=0 to prevent conductance'// &
      ' oscillations between outer iterations.'// &
      ' (3) Low K at topmost-active cells -- consider'// &
      ' setting minimum K values to prevent near-zero'// &
      ' conductance values.', terminate=.TRUE.)

  end subroutine amg_note_overflow

  !> @brief Update AMG matrix values without rebuilding the hierarchy structure
  !!
  !! Called when the sparsity pattern is unchanged (same neq/nja) but matrix
  !! values have changed (Newton update). Copies new fine-level values and
  !! recomputes coarse-level values and diagonal inverses at every level.
  !<
  subroutine ims_amg_update(amg, nja, a)
    type(ImsAmgDataType), intent(inout) :: amg !< AMG preconditioner
    integer(I4B), intent(in) :: nja !< number of fine-level non-zeros
    real(DP), dimension(nja), intent(in) :: a !< fine-level coefficient matrix
    ! -- local
    integer(I4B) :: l, k

    amg%overflow_this_outer = .false.
    do k = 1, nja
      amg%levels(0)%a(k) = a(k)
    end do

    call g_prof%start("AMG coarse values", amg%itmr_cvals)
    do l = 0, amg%nlevels - 2
      call build_coarse_values(amg%levels(l), amg%levels(l + 1))
    end do

    do l = 0, amg%nlevels - 1
      call update_diag(amg%levels(l))
    end do
    call g_prof%stop(amg%itmr_cvals)

    call g_prof%start("AMG ILU0 factor", amg%itmr_ilu)
    if (amg%smoother_type == SMOOTHER_ILU0) then
      call build_ilu0(amg%levels(0), .false.)
    else
      do l = 0, amg%nlevels - 1
        call build_ilu0(amg%levels(l), l > 0)
      end do
    end if
    call g_prof%stop(amg%itmr_ilu)

  end subroutine ims_amg_update

  !> @brief Build the inverse aggregation map (coarse aggregate -> fine cells)
  !!
  !! For each coarse aggregate ic, inv_ja(inv_ia(ic) .. inv_ia(ic + 1) - 1)
  !! lists the fine cells assigned to it. This is a counting sort of lev_f%agg,
  !! shared by build_coarse_matrix and build_coarse_values.
  !<
  subroutine build_inverse_aggregation(lev_f, inv_ia, inv_ja)
    type(AmgLevelType), intent(in) :: lev_f !< fine level (agg, neq, nc)
    integer(I4B), allocatable, intent(out) :: inv_ia(:) !< aggregate row pointers
    integer(I4B), allocatable, intent(out) :: inv_ja(:) !< fine cells, by aggregate
    ! -- local
    integer(I4B) :: i, ic, nc
    integer(I4B), allocatable :: cnt(:)

    nc = lev_f%nc
    allocate (inv_ia(nc + 1), inv_ja(lev_f%neq), cnt(nc))

    ! count fine cells per aggregate, then turn counts into CSR row pointers
    do i = 1, nc + 1
      inv_ia(i) = 0
    end do
    do i = 1, lev_f%neq
      inv_ia(lev_f%agg(i) + 1) = inv_ia(lev_f%agg(i) + 1) + 1
    end do
    inv_ia(1) = 1
    do ic = 1, nc
      inv_ia(ic + 1) = inv_ia(ic + 1) + inv_ia(ic)
    end do

    ! scatter each fine cell into its aggregate's slot
    do i = 1, nc
      cnt(i) = 0
    end do
    do i = 1, lev_f%neq
      ic = lev_f%agg(i)
      inv_ja(inv_ia(ic) + cnt(ic)) = i
      cnt(ic) = cnt(ic) + 1
    end do
    deallocate (cnt)

  end subroutine build_inverse_aggregation

  !> @brief Recompute coarse-level matrix values using the existing sparsity
  !!
  !! Performs only pass 2 of build_coarse_matrix: accumulates fine-level
  !! values into the already-allocated coarse CSR arrays lev_c%a.
  !<
  subroutine build_coarse_values(lev_f, lev_c)
    type(AmgLevelType), intent(in) :: lev_f !< fine level (ia, ja, a, agg)
    type(AmgLevelType), intent(inout) :: lev_c !< coarse level (ia, ja pre-built;
                                                 !! a updated)
    ! -- local
    integer(I4B) :: i, k, fi, ic, jc, nc, ipos, row_len
    integer(I4B), allocatable :: inv_ia(:), inv_ja(:)
    integer(I4B), allocatable :: mark(:), cols(:)
    real(DP), allocatable :: acols(:)

    nc = lev_f%nc
    call build_inverse_aggregation(lev_f, inv_ia, inv_ja)

    allocate (mark(nc), cols(nc), acols(nc))
    do i = 1, nc
      mark(i) = 0
      acols(i) = DZERO
    end do

    do ic = 1, nc
      row_len = 0
      do fi = inv_ia(ic), inv_ia(ic + 1) - 1
        i = inv_ja(fi)
        do k = lev_f%ia(i), lev_f%ia(i + 1) - 1
          jc = lev_f%agg(lev_f%ja(k))
          if (mark(jc) == 0) then
            row_len = row_len + 1
            mark(jc) = row_len
            cols(row_len) = jc
          end if
          acols(jc) = acols(jc) + lev_f%a(k)
        end do
      end do
      ipos = lev_c%ia(ic)
      do k = 1, row_len
        jc = cols(k)
        lev_c%a(ipos + k - 1) = acols(jc)
        mark(jc) = 0
        acols(jc) = DZERO
      end do
    end do

    deallocate (inv_ia, inv_ja, mark, cols, acols)

  end subroutine build_coarse_values

  !> @brief Recompute inverse diagonal for one AMG level
  !<
  subroutine update_diag(lev)
    type(AmgLevelType), intent(inout) :: lev !< AMG level
    ! -- local
    integer(I4B) :: i, k
    real(DP) :: dval, offsum

    do i = 1, lev%neq
      dval = DZERO
      offsum = DZERO
      do k = lev%ia(i), lev%ia(i + 1) - 1
        if (lev%ja(k) == i) then
          dval = lev%a(k)
        else
          offsum = offsum + abs(lev%a(k))
        end if
      end do
      ! -- l1 diagonal: a row that is not diagonally dominant is given the
      !    absolute row sum instead, which keeps the smoother from amplifying
      !    the error. Dominant rows, which is the usual case, are unchanged.
      if (abs(dval) < offsum) then
        if (dval < DZERO) then
          dval = -offsum
        else
          dval = offsum
        end if
      end if
      if (abs(dval) > DZERO) then
        lev%diag(i) = DONE / dval
      else
        lev%diag(i) = DONE
      end if
    end do

  end subroutine update_diag

  !> @brief Build greedy pairwise aggregation for one AMG level
  !<
  subroutine build_aggregation(lev, sthresh)
    type(AmgLevelType), intent(inout) :: lev !< AMG level
    real(DP), intent(in) :: sthresh !< strength-of-connection threshold
    ! -- local
    integer(I4B) :: i, j, k, kji
    integer(I4B) :: nc
    integer(I4B) :: best_j
    integer(I4B) :: i0, i1, j0, j1
    real(DP) :: maxoff_i
    real(DP) :: str, str_ji
    real(DP) :: best_str

    do i = 1, lev%neq
      lev%agg(i) = 0
    end do
    nc = 0
    lev%sthresh = sthresh
    lev%nsingle = 0
    lev%ncand = 0
    lev%nweak = 0

    do i = 1, lev%neq
      if (lev%agg(i) /= 0) cycle
      maxoff_i = DZERO
      i0 = lev%ia(i)
      i1 = lev%ia(i + 1) - 1
      do k = i0, i1
        j = lev%ja(k)
        if (j /= i) then
          ! symmetric strength: average of |a_ij| and |a_ji|
          str_ji = DZERO
          j0 = lev%ia(j)
          j1 = lev%ia(j + 1) - 1
          do kji = j0, j1
            if (lev%ja(kji) == i) then
              str_ji = abs(lev%a(kji))
              exit
            end if
          end do
          str = DHALF * (abs(lev%a(k)) + str_ji)
          if (str > maxoff_i) maxoff_i = str
        end if
      end do
      best_j = 0
      best_str = DZERO
      do k = i0, i1
        j = lev%ja(k)
        if (j == i) cycle
        if (lev%agg(j) /= 0) cycle
        str_ji = DZERO
        j0 = lev%ia(j)
        j1 = lev%ia(j + 1) - 1
        do kji = j0, j1
          if (lev%ja(kji) == i) then
            str_ji = abs(lev%a(kji))
            exit
          end if
        end do
        str = DHALF * (abs(lev%a(k)) + str_ji)
        lev%ncand = lev%ncand + 1
        if (sthresh > DZERO .and. str < sthresh * maxoff_i) then
          lev%nweak = lev%nweak + 1
          cycle
        end if
        if (str > best_str) then
          best_j = j
          best_str = str
        end if
      end do
      nc = nc + 1
      lev%agg(i) = nc
      if (best_j > 0) then
        lev%agg(best_j) = nc
      else
        lev%nsingle = lev%nsingle + 1
      end if
    end do

    lev%nc = nc

  end subroutine build_aggregation

  !> @brief Build the coarse-grid matrix using the aggregation from lev_f
  !!
  !! Builds the coarse matrix with a sparse "mark array" instead of a dense
  !! nc x nc scratch matrix. This cuts memory from O(nc^2) to O(nc + neq_f)
  !! and run time from O(nc^2) to O(nja_f). For a 12,000-cell model at the
  !! first coarsening level (nc ~ 6,000) it avoids a 288 MB dense allocation.
  !! It also stays practical on unstructured DISV and DISU grids, where nc^2
  !! can be far larger than the actual number of coarse connections.
  !<
  subroutine build_coarse_matrix(lev_f, lev_c)
    type(AmgLevelType), intent(inout) :: lev_f !< fine level
    type(AmgLevelType), intent(inout) :: lev_c !< coarse level (output)
    ! -- local
    integer(I4B) :: i, k, fi, ic, jc, nc, nja_c, ipos, row_len
    integer(I4B), allocatable :: inv_ia(:), inv_ja(:)
    integer(I4B), allocatable :: mark(:), cols(:)
    real(DP), allocatable :: acols(:)

    nc = lev_f%nc
    lev_c%neq = nc
    call build_inverse_aggregation(lev_f, inv_ia, inv_ja)

    ! Sparse accumulation using mark(jc) to detect new coarse connections.
    ! acols(jc) accumulates the coarse entry A_c(ic,jc) while processing row ic.
    allocate (mark(nc), cols(nc), acols(nc))
    do i = 1, nc
      mark(i) = 0
      acols(i) = DZERO
    end do

    ! Pass 1: determine coarse sparsity pattern (ia) and count nja_c.
    allocate (lev_c%ia(nc + 1))
    lev_c%ia(1) = 1
    nja_c = 0
    do ic = 1, nc
      row_len = 0
      do fi = inv_ia(ic), inv_ia(ic + 1) - 1
        i = inv_ja(fi)
        do k = lev_f%ia(i), lev_f%ia(i + 1) - 1
          jc = lev_f%agg(lev_f%ja(k))
          if (mark(jc) == 0) then
            row_len = row_len + 1
            mark(jc) = row_len
            cols(row_len) = jc
          end if
        end do
      end do
      nja_c = nja_c + row_len
      lev_c%ia(ic + 1) = nja_c + 1
      do k = 1, row_len
        mark(cols(k)) = 0
      end do
    end do

    allocate (lev_c%ja(nja_c), lev_c%a(nja_c))
    lev_c%nja = nja_c

    ! Pass 2: accumulate values and write to CSR.
    do ic = 1, nc
      row_len = 0
      do fi = inv_ia(ic), inv_ia(ic + 1) - 1
        i = inv_ja(fi)
        do k = lev_f%ia(i), lev_f%ia(i + 1) - 1
          jc = lev_f%agg(lev_f%ja(k))
          if (mark(jc) == 0) then
            row_len = row_len + 1
            mark(jc) = row_len
            cols(row_len) = jc
          end if
          acols(jc) = acols(jc) + lev_f%a(k)
        end do
      end do
      ipos = lev_c%ia(ic)
      do k = 1, row_len
        jc = cols(k)
        lev_c%ja(ipos + k - 1) = jc
        lev_c%a(ipos + k - 1) = acols(jc)
        mark(jc) = 0
        acols(jc) = DZERO
      end do
    end do

    deallocate (inv_ia, inv_ja, mark, cols, acols)

  end subroutine build_coarse_matrix

  !> @brief Allocate work vectors and compute inverse diagonal for one level
  !<
  subroutine alloc_work(lev)
    type(AmgLevelType), intent(inout) :: lev !< AMG level
    ! -- local
    integer(I4B) :: i

    allocate (lev%r(lev%neq))
    allocate (lev%e(lev%neq))
    do i = 1, lev%neq
      lev%r(i) = DZERO
      lev%e(i) = DZERO
    end do

    if (.not. allocated(lev%diag)) then
      allocate (lev%diag(lev%neq))
    end if
    call update_diag(lev)

  end subroutine alloc_work

  !> @brief Forward Gauss-Seidel smoother sweep
  !<
  subroutine gs_forward(lev, niter)
    type(AmgLevelType), intent(inout) :: lev !< AMG level
    integer(I4B), intent(in) :: niter !< number of sweeps
    ! -- local
    integer(I4B) :: iter, i, k
    real(DP) :: delta

    do iter = 1, niter
      do i = 1, lev%neq
        delta = lev%r(i)
        do k = lev%ia(i), lev%ia(i + 1) - 1
          delta = delta - lev%a(k) * lev%e(lev%ja(k))
        end do
        lev%e(i) = lev%e(i) + lev%diag(i) * delta
      end do
    end do

  end subroutine gs_forward

  !> @brief Backward Gauss-Seidel smoother sweep
  !<
  subroutine gs_backward(lev, niter)
    type(AmgLevelType), intent(inout) :: lev !< AMG level
    integer(I4B), intent(in) :: niter !< number of sweeps
    ! -- local
    integer(I4B) :: iter, i, k
    real(DP) :: delta

    do iter = 1, niter
      do i = lev%neq, 1, -1
        delta = lev%r(i)
        do k = lev%ia(i), lev%ia(i + 1) - 1
          delta = delta - lev%a(k) * lev%e(lev%ja(k))
        end do
        lev%e(i) = lev%e(i) + lev%diag(i) * delta
      end do
    end do

  end subroutine gs_backward

  !> @brief Downward (pre-smooth) sweep: ILU(0) if built, else forward GS
  !<
  subroutine smooth_down(lev, niter)
    type(AmgLevelType), intent(inout) :: lev !< AMG level
    integer(I4B), intent(in) :: niter !< number of sweeps
    if (allocated(lev%apc)) then
      call ilu0_smooth(lev, niter)
    else
      call gs_forward(lev, niter)
    end if
  end subroutine smooth_down

  !> @brief Upward (post-smooth) sweep: ILU(0) if built, else backward GS
  !<
  subroutine smooth_up(lev, niter)
    type(AmgLevelType), intent(inout) :: lev !< AMG level
    integer(I4B), intent(in) :: niter !< number of sweeps
    if (allocated(lev%apc)) then
      call ilu0_smooth(lev, niter)
    else
      call gs_backward(lev, niter)
    end if
  end subroutine smooth_up

  !> @brief Restrict the fine-level residual to the coarse level
  !!
  !! Forms the fine residual r_f - A_f e_f, sums it into the coarse right-hand
  !! side by aggregate, and zeros the coarse error ready for the next solve.
  !<
  subroutine restrict_residual(lev_f, lev_c)
    type(AmgLevelType), intent(in) :: lev_f !< fine level
    type(AmgLevelType), intent(inout) :: lev_c !< coarse level
    ! -- local
    integer(I4B) :: i, k, c
    real(DP) :: res_i

    do i = 1, lev_c%neq
      lev_c%r(i) = DZERO
    end do

    do i = 1, lev_f%neq
      res_i = lev_f%r(i)
      do k = lev_f%ia(i), lev_f%ia(i + 1) - 1
        res_i = res_i - lev_f%a(k) * lev_f%e(lev_f%ja(k))
      end do
      c = lev_f%agg(i)
      lev_c%r(c) = lev_c%r(c) + res_i
    end do

    do i = 1, lev_c%neq
      lev_c%e(i) = DZERO
    end do
  end subroutine restrict_residual

  !> @brief Add the coarse-level correction back onto the fine-level error
  !<
  subroutine prolong_add(lev_f, lev_c)
    type(AmgLevelType), intent(inout) :: lev_f !< fine level
    type(AmgLevelType), intent(in) :: lev_c !< coarse level
    ! -- local
    integer(I4B) :: i

    do i = 1, lev_f%neq
      lev_f%e(i) = lev_f%e(i) + lev_c%e(lev_f%agg(i))
    end do
  end subroutine prolong_add

  !> @brief Perform one AMG V-cycle
  !!
  !! Non-recursive implementation. Descends from the fine level to the
  !! coarsest level smoothing and restricting, then ascends with corrections.
  !<
  subroutine amg_vcycle(levels, nlevels, nsmooth)
    type(AmgLevelType), intent(inout) :: levels(0:) !< AMG level array
    integer(I4B), intent(in) :: nlevels !< number of active levels
    integer(I4B), intent(in) :: nsmooth !< smoother iterations
    ! -- local
    integer(I4B) :: l

    ! Descend: pre-smooth then restrict the residual to the next coarser level
    do l = 0, nlevels - 2
      call g_prof%start(levels(l)%tmr_title, levels(l)%itmr_smooth)
      call smooth_down(levels(l), nsmooth)
      call g_prof%stop(levels(l)%itmr_smooth)
      call restrict_residual(levels(l), levels(l + 1))
    end do

    ! Approximate solve on the coarsest level
    call g_prof%start(levels(nlevels - 1)%tmr_title, &
                      levels(nlevels - 1)%itmr_smooth)
    call coarse_solve(levels(nlevels - 1))
    call g_prof%stop(levels(nlevels - 1)%itmr_smooth)

    ! Ascend: prolong the coarse correction then post-smooth
    do l = nlevels - 2, 0, -1
      call prolong_add(levels(l), levels(l + 1))
      call g_prof%start(levels(l)%tmr_title, levels(l)%itmr_smooth)
      call smooth_up(levels(l), nsmooth)
      call g_prof%stop(levels(l)%itmr_smooth)
    end do

  end subroutine amg_vcycle

  !> @brief Euclidean norm of the residual on one AMG level
  !<
  function level_resnorm(lev) result(rnorm)
    type(AmgLevelType), intent(in) :: lev !< AMG level
    real(DP) :: rnorm !< residual norm
    ! -- local
    integer(I4B) :: i, k
    real(DP) :: res

    rnorm = DZERO
    do i = 1, lev%neq
      res = lev%r(i)
      do k = lev%ia(i), lev%ia(i + 1) - 1
        res = res - lev%a(k) * lev%e(lev%ja(k))
      end do
      rnorm = rnorm + res * res
    end do
    rnorm = sqrt(rnorm)

  end function level_resnorm

  !> @brief Solve the coarsest level with the smoother
  !!
  !! The correction is discarded if the residual grows, so the sweep count is
  !! unchanged on a healthy problem. Gauss-Seidel and ILU(0)
  !! only converge for a diagonally dominant or definite matrix, and a coarse
  !! matrix formed from an anisotropic Newton matrix need not be either, so
  !! sweeping a fixed number of times can amplify the error instead of
  !! reducing it.
  !<
  subroutine coarse_solve(lev)
    type(AmgLevelType), intent(inout) :: lev !< coarsest AMG level
    ! -- local
    integer(I4B) :: it, i
    real(DP) :: rnorm, rnorm0

    rnorm0 = level_resnorm(lev)
    if (rnorm0 <= DZERO) return

    do it = 1, COARSE_SOLVE_ITERS
      call smooth_down(lev, 1)
      rnorm = level_resnorm(lev)
      ! -- the comparison is negated so a NaN residual also exits
      if (.not. (rnorm <= rnorm0)) then
        do i = 1, lev%neq
          lev%e(i) = DZERO
        end do
        exit
      end if
    end do

  end subroutine coarse_solve

  !> @brief Build or update the ILU(0) factorization for the finest AMG level
  !<
  subroutine build_ilu0(lev, use_l1)
    type(AmgLevelType), intent(inout) :: lev !< AMG level
    logical, intent(in) :: use_l1 !< factor an l1-modified matrix
    ! -- local
    integer(I4B), allocatable :: iw(:)
    real(DP), allocatable :: w(:)
    integer(I4B), allocatable :: dpos(:)
    real(DP), allocatable :: dsave(:)
    integer(I4B) :: ipcflag
    integer(I4B) :: i, k, n
    real(DP) :: amax

    ! Guard: skip factorization if matrix contains overflow-scale values;
    ! the V-cycle output check in ims_amg_apply will catch the bad result.
    amax = DZERO
    do n = 1, lev%nja
      if (abs(lev%a(n)) > amax) amax = abs(lev%a(n))
    end do
    if (amax > DEP20) return

    if (.not. allocated(lev%iapc)) then
      allocate (lev%iapc(lev%neq + 1))
      allocate (lev%japc(lev%nja))
      allocate (lev%apc(lev%nja))
      allocate (lev%work(lev%neq))
      do n = 1, lev%neq + 1
        lev%iapc(n) = 0
      end do
      do n = 1, lev%nja
        lev%japc(n) = 0
        lev%apc(n) = DZERO
      end do
      do n = 1, lev%neq
        lev%work(n) = DZERO
      end do
      call ims_base_pccrs(lev%neq, lev%nja, lev%ia, lev%ja, lev%iapc, lev%japc)
    end if

    allocate (iw(lev%neq), w(lev%neq))
    ipcflag = 0
    ! -- factor an l1-modified matrix on the coarser levels, where Galerkin
    !    coarsening of a Newton matrix can leave rows that are not diagonally
    !    dominant and whose factorization amplifies the error. lev%a is
    !    restored afterwards because the residual calculations need it.
    if (use_l1) then
      allocate (dpos(lev%neq), dsave(lev%neq))
      do i = 1, lev%neq
        dpos(i) = 0
        do k = lev%ia(i), lev%ia(i + 1) - 1
          if (lev%ja(k) == i) then
            dpos(i) = k
            dsave(i) = lev%a(k)
            lev%a(k) = DONE / lev%diag(i)
            exit
          end if
        end do
      end do
    end if

    call ims_base_pcilu0(lev%nja, lev%neq, lev%a, lev%ia, lev%ja, &
                         lev%apc, lev%iapc, lev%japc, iw, w, &
                         DZERO, ipcflag, DZERO)
    deallocate (iw, w)

    if (use_l1) then
      do i = 1, lev%neq
        if (dpos(i) > 0) lev%a(dpos(i)) = dsave(i)
      end do
      deallocate (dpos, dsave)
    end if

  end subroutine build_ilu0

  !> @brief Apply ILU(0) correction sweeps as the AMG smoother
  !! at the finest level
  !!
  !! Each sweep computes the residual, applies the ILU(0) solve, and adds the
  !! result to the current error estimate.
  !<
  subroutine ilu0_smooth(lev, niter)
    type(AmgLevelType), intent(inout) :: lev !< finest AMG level
    integer(I4B), intent(in) :: niter !< number of correction sweeps
    ! -- local
    integer(I4B) :: iter, i, k

    do iter = 1, niter
      do i = 1, lev%neq
        lev%work(i) = lev%r(i)
        do k = lev%ia(i), lev%ia(i + 1) - 1
          lev%work(i) = lev%work(i) - lev%a(k) * lev%e(lev%ja(k))
        end do
      end do
      call ims_base_ilu0a(lev%nja, lev%neq, lev%apc, lev%iapc, lev%japc, &
                          lev%work, lev%work)
      do i = 1, lev%neq
        lev%e(i) = lev%e(i) + lev%work(i)
      end do
    end do

  end subroutine ilu0_smooth

end module ImsLinearAmgModule
