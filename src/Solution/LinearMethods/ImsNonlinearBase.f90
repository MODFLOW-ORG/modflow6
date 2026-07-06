
!> @brief This module contains the IMS nonlinear (outer) accelerator subroutines
!!
!! This module contains the stateless IMS nonlinear (outer-iteration)
!! subroutines used by a MODFLOW 6 solution. The routines are free
!! subroutines that take explicit arguments (mirroring IMSLinearBaseModule)
!! so that they can be reused outside of NumericalSolutionType and unit
!! tested directly. All computations are local; in parallel the caller is
!! responsible for supplying already-reduced quantities (e.g. bigch).
!<
module ImsNonlinearBaseModule
  ! -- modules
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DHALF, DONE, DTHREE, DEM20, DPREC
  use MatrixBaseModule, only: MatrixBaseType
  use VectorBaseModule, only: VectorBaseType

  implicit none
  private

  public :: ims_nl_underrelax
  public :: ims_nl_maxval
  public :: ims_nl_calcdx
  public :: ims_nl_dxmax
  public :: ims_nl_backtrack_flag
  public :: ims_nl_apply_backtrack
  public :: ims_nl_has_converged
  public :: ims_nl_nur_has_converged
  public :: ims_nl_residual
  public :: ims_nl_l2norm

contains

  !> @ brief Under-relaxation
  !!
  !!  Under relax using the simple, cooley, or delta-bar-delta methods.
  !!  This is the local (stateless) body of NumericalSolutionType%sln_underrelax;
  !!  the persistent under-relaxation state is passed in/out explicitly.
  !<
  subroutine ims_nl_underrelax(nonmeth, kiter, bigch, neq, active, x, xtemp, &
                               gamma, theta, akappa, amomentum, &
                               bigch_store, bigchold, relaxold, &
                               dxold, wsave, hchold, deold)
    ! -- dummy variables
    integer(I4B), intent(in) :: nonmeth !< under-relaxation method used
    integer(I4B), intent(in) :: kiter !< Picard iteration number
    real(DP), intent(in) :: bigch !< maximum dependent-variable change
    integer(I4B), intent(in) :: neq !< number of equations
    integer(I4B), dimension(neq), intent(in) :: active !< active cell flag (1)
    real(DP), dimension(neq), intent(inout) :: x !< current dependent-variable
    real(DP), dimension(neq), intent(in) :: xtemp !< previous dependent-variable
    real(DP), intent(in) :: gamma !< under-relaxation gamma
    real(DP), intent(in) :: theta !< under-relaxation theta
    real(DP), intent(in) :: akappa !< under-relaxation kappa
    real(DP), intent(in) :: amomentum !< under-relaxation momentum term
    real(DP), intent(inout) :: bigch_store !< stored maximum dependent-variable change
    real(DP), intent(inout) :: bigchold !< cooley under-relaxation weight
    real(DP), intent(inout) :: relaxold !< previous relaxation factor
    real(DP), dimension(neq), intent(inout) :: dxold !< previous dependent-variable change
    real(DP), dimension(neq), intent(inout) :: wsave !< DBD sign-change factor
    real(DP), dimension(neq), intent(inout) :: hchold !< DBD weighted dependent-variable change
    real(DP), dimension(neq), intent(inout) :: deold !< DBD dependent-variable change variable
    ! -- local variables
    integer(I4B) :: n
    real(DP) :: ww
    real(DP) :: delx
    real(DP) :: relax
    real(DP) :: es
    real(DP) :: aes
    real(DP) :: amom
    !
    ! -- option for using simple dampening (as done by MODFLOW-2005 PCG)
    if (nonmeth == 1) then
      do n = 1, neq
        !
        ! -- skip inactive nodes
        if (active(n) < 1) cycle
        !
        ! -- compute step-size (delta x)
        delx = x(n) - xtemp(n)
        dxold(n) = delx

        ! -- dampen dependent variable solution
        x(n) = xtemp(n) + gamma * delx
      end do
      !
      ! -- option for using cooley underrelaxation
    else if (nonmeth == 2) then
      !
      ! -- set bigch
      bigch_store = bigch
      !
      ! -- initialize values for first iteration
      if (kiter == 1) then
        relax = DONE
        relaxold = DONE
        bigchold = bigch
      else
        !
        ! -- compute relaxation factor
        es = bigch_store / (bigchold * relaxold)
        aes = abs(es)
        if (es < -DONE) then
          relax = dhalf / aes
        else
          relax = (DTHREE + es) / (DTHREE + aes)
        end if
      end if
      relaxold = relax
      !
      ! -- modify cooley to use weighted average of past changes
      bigchold = (DONE - gamma) * bigch_store + gamma * &
                 bigchold
      !
      ! -- compute new dependent variable after under-relaxation
      if (relax < DONE) then
        do n = 1, neq
          !
          ! -- skip inactive nodes
          if (active(n) < 1) cycle
          !
          ! -- update dependent variable
          delx = x(n) - xtemp(n)
          dxold(n) = delx
          x(n) = xtemp(n) + relax * delx
        end do
      end if
      !
      ! -- option for using delta-bar-delta scheme to under-relax for all equations
    else if (nonmeth == 3) then
      do n = 1, neq
        !
        ! -- skip inactive nodes
        if (active(n) < 1) cycle
        !
        ! -- compute step-size (delta x) and initialize d-b-d parameters
        delx = x(n) - xtemp(n)
        !
        ! -- initialize values for first iteration
        if (kiter == 1) then
          wsave(n) = DONE
          hchold(n) = DEM20
          deold(n) = DZERO
        end if
        !
        ! -- compute new relaxation term as per delta-bar-delta
        ww = wsave(n)
        !
        ! for flip-flop condition, decrease factor
        if (deold(n) * delx < DZERO) then
          ww = theta * wsave(n)
          ! -- when change is of same sign, increase factor
        else
          ww = wsave(n) + akappa
        end if
        if (ww > DONE) ww = DONE
        wsave(n) = ww
        !
        ! -- compute weighted average of past changes in hchold
        if (kiter == 1) then
          hchold(n) = delx
        else
          hchold(n) = (DONE - gamma) * delx + &
                      gamma * hchold(n)
        end if
        !
        ! -- store slope (change) term for next iteration
        deold(n) = delx
        dxold(n) = delx
        !
        ! -- compute accepted step-size and new dependent variable
        amom = DZERO
        if (kiter > 4) amom = amomentum
        delx = delx * ww + amom * hchold(n)
        x(n) = xtemp(n) + delx
      end do
      !
    end if
  end subroutine ims_nl_underrelax

  !> @ brief Determine maximum value
  !!
  !!  Determine the maximum value in a vector using a normalized comparison.
  !!  This is the local body of NumericalSolutionType%sln_maxval.
  !<
  subroutine ims_nl_maxval(nsize, v, vmax)
    ! -- dummy variables
    integer(I4B), intent(in) :: nsize !< length of vector
    real(DP), dimension(nsize), intent(in) :: v !< input vector
    real(DP), intent(inout) :: vmax !< maximum value
    ! -- local variables
    integer(I4B) :: n
    real(DP) :: d
    real(DP) :: denom
    real(DP) :: dnorm
    !
    ! -- determine maximum value
    vmax = v(1)
    do n = 2, nsize
      d = v(n)
      denom = abs(vmax)
      if (denom == DZERO) then
        denom = DPREC
      end if
      !
      ! -- calculate normalized value
      dnorm = abs(d) / denom
      if (dnorm > DONE) then
        vmax = d
      end if
    end do
  end subroutine ims_nl_maxval

  !> @ brief Calculate dependent-variable change
  !!
  !!  Calculate the dependent-variable change for every cell. This is the
  !!  local body of NumericalSolutionType%sln_calcdx.
  !<
  subroutine ims_nl_calcdx(neq, active, x, xtemp, dx)
    ! -- dummy variables
    integer(I4B), intent(in) :: neq !< number of equations
    integer(I4B), dimension(neq), intent(in) :: active !< active cell flag (1)
    real(DP), dimension(neq), intent(in) :: x !< current dependent-variable
    real(DP), dimension(neq), intent(in) :: xtemp !< previous dependent-variable
    real(DP), dimension(neq), intent(inout) :: dx !< dependent-variable change
    ! -- local
    integer(I4B) :: n
    !
    ! -- calculate dependent-variable change
    do n = 1, neq
      ! -- skip inactive nodes
      if (active(n) < 1) then
        dx(n) = DZERO
      else
        dx(n) = x(n) - xtemp(n)
      end if
    end do
  end subroutine ims_nl_calcdx

  !> @ brief Determine maximum dependent-variable change
  !!
  !!  Determine the maximum dependent-variable change and its location.
  !!  This is the local body of NumericalSolutionType%sln_get_dxmax.
  !<
  subroutine ims_nl_dxmax(neq, active, x, xtemp, dxmax, loc)
    ! -- dummy variables
    integer(I4B), intent(in) :: neq !< number of equations
    integer(I4B), dimension(neq), intent(in) :: active !< active cell flag (1)
    real(DP), dimension(neq), intent(in) :: x !< current dependent-variable
    real(DP), dimension(neq), intent(in) :: xtemp !< previous dependent-variable
    real(DP), intent(inout) :: dxmax !< maximum dependent-variable change
    integer(I4B), intent(inout) :: loc !< location of the maximum change
    ! -- local variables
    integer(I4B) :: nb
    real(DP) :: bigch
    real(DP) :: abigch
    integer(I4B) :: n
    real(DP) :: hdif
    real(DP) :: ahdif
    !
    ! -- determine the maximum change
    nb = 0
    bigch = DZERO
    abigch = DZERO
    do n = 1, neq
      if (active(n) < 1) cycle
      hdif = x(n) - xtemp(n)
      ahdif = abs(hdif)
      if (ahdif > abigch) then
        bigch = hdif
        abigch = ahdif
        nb = n
      end if
    end do
    !
    !-----store maximum change value and location
    dxmax = bigch
    loc = nb
  end subroutine ims_nl_dxmax

  !> @ brief Determine the backtracking flag
  !!
  !!  Return 1 when the maximum dependent-variable change scaled by the
  !!  reduction factor still exceeds the closure criterion. This is the
  !!  local body of NumericalSolutionType%get_backtracking_flag.
  !<
  function ims_nl_backtrack_flag(neq, active, x, xtemp, breduc, dvclose) &
    result(bt_flag)
    ! -- dummy variables
    integer(I4B), intent(in) :: neq !< number of equations
    integer(I4B), dimension(neq), intent(in) :: active !< active cell flag (1)
    real(DP), dimension(neq), intent(in) :: x !< current dependent-variable
    real(DP), dimension(neq), intent(in) :: xtemp !< previous dependent-variable
    real(DP), intent(in) :: breduc !< backtracking reduction factor
    real(DP), intent(in) :: dvclose !< dependent-variable closure criterion
    ! -- return
    integer(I4B) :: bt_flag !< (1) backtracking performed (0) not performed
    ! -- local
    integer(I4B) :: n
    real(DP) :: dx
    real(DP) :: dx_abs
    real(DP) :: dx_abs_max

    ! default is off
    bt_flag = 0

    ! find max. change
    dx_abs_max = 0.0
    do n = 1, neq
      if (active(n) < 1) cycle
      dx = x(n) - xtemp(n)
      dx_abs = abs(dx)
      if (dx_abs > dx_abs_max) dx_abs_max = dx_abs
    end do

    ! if backtracking, set flag
    if (breduc * dx_abs_max >= dvclose) then
      bt_flag = 1
    end if

  end function ims_nl_backtrack_flag

  !> @ brief Update x with backtracking
  !!
  !!  Scale the dependent-variable change by the reduction factor. This is
  !!  the local body of NumericalSolutionType%apply_backtracking.
  !<
  subroutine ims_nl_apply_backtrack(neq, active, x, xtemp, breduc)
    ! -- dummy variables
    integer(I4B), intent(in) :: neq !< number of equations
    integer(I4B), dimension(neq), intent(in) :: active !< active cell flag (1)
    real(DP), dimension(neq), intent(inout) :: x !< current dependent-variable
    real(DP), dimension(neq), intent(in) :: xtemp !< previous dependent-variable
    real(DP), intent(in) :: breduc !< backtracking reduction factor
    ! -- local
    integer(I4B) :: n
    real(DP) :: delx

    do n = 1, neq
      if (active(n) < 1) cycle
      delx = breduc * (x(n) - xtemp(n))
      x(n) = xtemp(n) + delx
    end do

  end subroutine ims_nl_apply_backtrack

  !> @ brief Convergence check
  !!
  !!  True when the maximum dependent-variable change is within the closure
  !!  criterion. This is the local body of
  !!  NumericalSolutionType%sln_has_converged.
  !<
  function ims_nl_has_converged(max_dvc, dvclose) result(has_converged)
    ! -- dummy variables
    real(DP), intent(in) :: max_dvc !< the maximum dependent variable change
    real(DP), intent(in) :: dvclose !< dependent-variable closure criterion
    ! -- return
    logical(LGP) :: has_converged !< True, when converged

    has_converged = .false.
    if (abs(max_dvc) <= dvclose) then
      has_converged = .true.
    end if

  end function ims_nl_has_converged

  !> @ brief Convergence check after Newton under-relaxation
  !!
  !!  True when both the unrelaxed maximum change and the largest change at
  !!  the end of the Picard iteration are within the closure criterion. This
  !!  is the local body of NumericalSolutionType%sln_nur_has_converged.
  !<
  function ims_nl_nur_has_converged(dxold_max, hncg, dvclose) &
    result(has_converged)
    ! -- dummy variables
    real(DP), intent(in) :: dxold_max !< maximum change for unrelaxed cells
    real(DP), intent(in) :: hncg !< largest change at end of Picard iteration
    real(DP), intent(in) :: dvclose !< dependent-variable closure criterion
    ! -- return
    logical(LGP) :: has_converged !< True, when converged

    has_converged = .false.
    if (abs(dxold_max) <= dvclose .and. &
        abs(hncg) <= dvclose) then
      has_converged = .true.
    end if

  end function ims_nl_nur_has_converged

  !> @ brief Calculate the residual vector r = A*x - b
  !!
  !!  Compute the residual and zero it for inactive cells. This is the local
  !!  body of NumericalSolutionType%sln_calc_residual. The caller supplies the
  !!  system matrix, the x and rhs vectors, and a residual vector to fill.
  !<
  subroutine ims_nl_residual(matrix, vec_x, vec_rhs, neq, active, vec_resid)
    ! -- dummy variables
    class(MatrixBaseType) :: matrix !< the system matrix A
    class(VectorBaseType), pointer :: vec_x !< the dependent-variable vector x
    class(VectorBaseType), pointer :: vec_rhs !< the right-hand side vector b
    integer(I4B), intent(in) :: neq !< number of equations
    integer(I4B), dimension(neq), intent(in) :: active !< active cell flag (1)
    class(VectorBaseType), pointer :: vec_resid !< the residual vector to fill
    ! -- local
    integer(I4B) :: n

    call matrix%multiply(vec_x, vec_resid) ! r = A*x

    call vec_resid%axpy(-1.0_DP, vec_rhs) ! r = r - b

    do n = 1, neq
      if (active(n) < 1) then
        call vec_resid%set_value_local(n, 0.0_DP) ! r_i = 0 if inactive
      end if
    end do

  end subroutine ims_nl_residual

  !> @ brief Calculate the L-2 norm of the residual for all active cells
  !!
  !!  L2norm = || A*x - b ||_2 with inactive cells zeroed. This is the local
  !!  body of NumericalSolutionType%sln_l2norm; a temporary residual vector is
  !!  created from the matrix and released here.
  !<
  subroutine ims_nl_l2norm(matrix, vec_x, vec_rhs, neq, active, l2norm)
    ! -- dummy variables
    class(MatrixBaseType) :: matrix !< the system matrix A
    class(VectorBaseType), pointer :: vec_x !< the dependent-variable vector x
    class(VectorBaseType), pointer :: vec_rhs !< the right-hand side vector b
    integer(I4B), intent(in) :: neq !< number of equations
    integer(I4B), dimension(neq), intent(in) :: active !< active cell flag (1)
    real(DP), intent(inout) :: l2norm !< calculated L-2 norm
    ! -- local
    class(VectorBaseType), pointer :: vec_resid

    ! calc. residual vector
    vec_resid => matrix%create_vec(neq)
    call ims_nl_residual(matrix, vec_x, vec_rhs, neq, active, vec_resid)

    ! 2-norm
    l2norm = vec_resid%norm2()

    ! clean up temp. vector
    call vec_resid%destroy()
    deallocate (vec_resid)

  end subroutine ims_nl_l2norm

end module ImsNonlinearBaseModule
