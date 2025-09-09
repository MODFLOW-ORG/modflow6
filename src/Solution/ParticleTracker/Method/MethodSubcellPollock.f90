module MethodSubcellPollockModule
  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use MethodModule, only: LEVEL_SUBFEATURE
  use MethodSubcellModule, only: MethodSubcellType
  use SubcellRectModule, only: SubcellRectType, create_subcell_rect
  use ParticleModule, only: ParticleType
  use PrtFmiModule, only: PrtFmiType
  use BaseDisModule, only: DisBaseType
  use CellModule, only: CellType
  use ConstantsModule, only: DZERO, DONE
  use DomainModule, only: DomainType
  use SubcellModule, only: SubcellType
  use ListModule, only: ListType
  use ExitSolutionModule, only: ExitSolutionType, &
                                LinearExitSolutionType, &
                                OK_EXIT, OK_EXIT_CONSTANT, &
                                NO_EXIT_STATIONARY, NO_EXIT_NO_OUTFLOW
  implicit none
  private
  public :: MethodSubcellPollockType
  public :: create_method_subcell_pollock
  public :: calculate_dt

  !> @brief Rectangular subcell tracking method
  type, extends(MethodSubcellType) :: MethodSubcellPollockType
    private
    real(DP), allocatable, public :: qextl1(:), qextl2(:), qintl(:) !< external and internal subcell flows
    type(LinearExitSolutionType), public :: exit_solutions(3) !< candidate exit solutions
  contains
    procedure, public :: find_exits
    procedure, public :: pick_exit
    procedure, public :: apply => apply_msp
    procedure, public :: deallocate
    procedure, private :: track_subcell
  end type MethodSubcellPollockType

contains

  !> @brief Create a new Pollock's subcell method
  subroutine create_method_subcell_pollock(method)
    ! dummy
    type(MethodSubcellPollockType), pointer :: method
    ! local
    type(SubcellRectType), pointer :: subcell

    allocate (method)
    call create_subcell_rect(subcell)
    method%subcell => subcell
    method%name => method%subcell%type
    method%delegates = .false.
  end subroutine create_method_subcell_pollock

  !> @brief Deallocate the Pollock's subcell method
  subroutine deallocate (this)
    class(MethodSubcellPollockType), intent(inout) :: this
    deallocate (this%name)
  end subroutine deallocate

  !> @brief Apply Pollock's method to a rectangular subcell
  subroutine apply_msp(this, particle, tmax)
    ! dummy
    class(MethodSubcellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local
    real(DP) :: x_origin
    real(DP) :: y_origin
    real(DP) :: z_origin
    real(DP) :: sinrot
    real(DP) :: cosrot

    select type (subcell => this%subcell)
    type is (SubcellRectType)
      ! Transform particle position into local subcell coordinates,
      ! track particle across subcell, convert back to model coords
      ! (sinrot and cosrot should be 0 and 1, respectively, i.e. no
      ! rotation, also no z translation; only x and y translations)
      x_origin = subcell%xOrigin
      y_origin = subcell%yOrigin
      z_origin = subcell%zOrigin
      sinrot = subcell%sinrot
      cosrot = subcell%cosrot
      call particle%transform(x_origin, y_origin)
      call this%track_subcell(subcell, particle, tmax)
      call particle%transform(x_origin, y_origin, invert=.true.)
    end select
  end subroutine apply_msp

  !> @brief Track a particle across a rectangular subcell using Pollock's method
  !!
  !! This subroutine consists partly of code written by
  !! David W. Pollock of the USGS for MODPATH 7. PRT's
  !! authors take responsibility for its application in
  !! this context and for any modifications or errors.
  !<
  subroutine track_subcell(this, subcell, particle, tmax)
    use TdisModule, only: endofsimulation
    use ParticleModule, only: ACTIVE, TERM_NO_EXITS_SUB, TERM_TIMEOUT
    use ParticleEventModule, only: TIMESTEP, FEATEXIT
    ! dummy
    class(MethodSubcellPollockType), intent(inout) :: this
    class(SubcellRectType), intent(in) :: subcell
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local
    real(DP) :: dt, dtexit, texit
    real(DP) :: t, x, y, z
    real(DP) :: t0, x0, y0, z0
    integer(I4B) :: i, exit_face, exit_soln
    type(LinearExitSolutionType) :: exit_x, exit_y, exit_z

    t0 = particle%ttrack
    x0 = particle%x / subcell%dx
    y0 = particle%y / subcell%dy
    z0 = particle%z / subcell%dz

    ! Find exit solution in each direction
    call this%find_exits(particle, subcell)

    exit_x = this%exit_solutions(1)
    exit_y = this%exit_solutions(2)
    exit_z = this%exit_solutions(3)

    ! Subcell has no exit face, terminate the particle
    ! TODO: consider ramifications
    if (all([this%exit_solutions%status] == NO_EXIT_NO_OUTFLOW)) then
      call this%terminate(particle, status=TERM_NO_EXITS_SUB)
      return
    end if

    ! If particle stationary and extended tracking is on, terminate here if it's
    ! the last timestep. TODO: temporary solution, consider where to catch this?
    ! Should we really have to special case this here? We do diverge from MP7 in
    ! guaranteeing that every particle terminates at the end of the simulation..
    ! ideally that would be handled at a higher scope but with extended tracking
    ! tmax is not the end of the simulation, it's just a wildly high upper bound.
    if (all([this%exit_solutions%status] == NO_EXIT_STATIONARY) .and. &
        particle%extend .and. endofsimulation) then
      call this%terminate(particle, status=TERM_TIMEOUT)
      return
    end if

    ! Pick exit solution, face, travel time, and time
    exit_soln = this%pick_exit(particle)
    if (exit_soln == 0) then
      exit_face = 0
      dtexit = 1.0d+30
    else
      exit_face = this%exit_solutions(exit_soln)%iboundary
      dtexit = this%exit_solutions(exit_soln)%dt
    end if
    texit = particle%ttrack + dtexit

    ! Select user tracking times to solve. If this is the first time step
    ! of the simulation, include all times before it begins; if it is the
    ! last time step include all times after it ends only if the 'extend'
    ! option is on, otherwise times in this period and time step only.
    call this%tracktimes%advance()
    if (this%tracktimes%any()) then
      do i = this%tracktimes%selection(1), this%tracktimes%selection(2)
        t = this%tracktimes%times(i)
        if (t < particle%ttrack) cycle
        if (t >= texit .or. t >= tmax) exit
        dt = t - t0
        x = new_x(exit_x%v, exit_x%dvdx, subcell%vx1, subcell%vx2, &
                  dt, x0, subcell%dx, exit_x%status == 1)
        y = new_x(exit_y%v, exit_y%dvdx, subcell%vy1, subcell%vy2, &
                  dt, y0, subcell%dy, exit_y%status == 1)
        z = new_x(exit_z%v, exit_z%dvdx, subcell%vz1, subcell%vz2, &
                  dt, z0, subcell%dz, exit_z%status == 1)
        particle%x = x * subcell%dx
        particle%y = y * subcell%dy
        particle%z = z * subcell%dz
        particle%ttrack = t
        particle%istatus = ACTIVE
        call this%usertime(particle)
      end do
    end if

    if (texit .gt. tmax) then
      ! The computed exit time is greater than the maximum time, so set
      ! final time for particle trajectory equal to maximum time and
      ! calculate particle location at that final time.
      t = tmax
      dt = t - t0
      x = new_x(exit_x%v, exit_x%dvdx, subcell%vx1, subcell%vx2, &
                dt, x0, subcell%dx, exit_x%status == 1)
      y = new_x(exit_y%v, exit_y%dvdx, subcell%vy1, subcell%vy2, &
                dt, y0, subcell%dy, exit_y%status == 1)
      z = new_x(exit_z%v, exit_z%dvdx, subcell%vz1, subcell%vz2, &
                dt, z0, subcell%dz, exit_z%status == 1)
      exit_face = 0
      particle%istatus = ACTIVE
      particle%advancing = .false.

      ! Set final particle location in local (unscaled) subcell coordinates,
      ! final time for particle trajectory
      particle%x = x * subcell%dx
      particle%y = y * subcell%dy
      particle%z = z * subcell%dz
      particle%ttrack = t
      particle%iboundary(LEVEL_SUBFEATURE) = exit_face

      ! Save particle track record
      call this%timestep(particle)
    else
      ! The computed exit time is less than or equal to the maximum time,
      ! so set final time for particle trajectory equal to exit time and
      ! calculate exit location.
      t = texit
      dt = dtexit
      if ((exit_face .eq. 1) .or. (exit_face .eq. 2)) then
        x = DZERO
        y = new_x(exit_y%v, exit_y%dvdx, subcell%vy1, subcell%vy2, &
                  dt, y0, subcell%dy, exit_y%status == 1)
        z = new_x(exit_z%v, exit_z%dvdx, subcell%vz1, subcell%vz2, &
                  dt, z0, subcell%dz, exit_z%status == 1)
        if (exit_face .eq. 2) x = DONE
      else if ((exit_face .eq. 3) .or. (exit_face .eq. 4)) then
        x = new_x(exit_x%v, exit_x%dvdx, subcell%vx1, subcell%vx2, dt, &
                  x0, subcell%dx, exit_x%status == 1)
        y = DZERO
        z = new_x(exit_z%v, exit_z%dvdx, subcell%vz1, subcell%vz2, dt, &
                  z0, subcell%dz, exit_z%status == 1)
        if (exit_face .eq. 4) y = DONE
      else if ((exit_face .eq. 5) .or. (exit_face .eq. 6)) then
        x = new_x(exit_x%v, exit_x%dvdx, subcell%vx1, subcell%vx2, &
                  dt, x0, subcell%dx, exit_x%status == 1)
        y = new_x(exit_y%v, exit_y%dvdx, subcell%vy1, subcell%vy2, &
                  dt, y0, subcell%dy, exit_y%status == 1)
        z = DZERO
        if (exit_face .eq. 6) z = DONE
      else
        print *, "programmer error, invalid exit face", exit_face
        call pstop(1)
      end if

      ! Set final particle location in local (unscaled) subcell coordinates,
      ! final time for particle trajectory, and exit face
      particle%x = x * subcell%dx
      particle%y = y * subcell%dy
      particle%z = z * subcell%dz
      particle%ttrack = t
      particle%iboundary(LEVEL_SUBFEATURE) = exit_face

      ! Save particle track record
      call this%subcellexit(particle)
    end if

  end subroutine track_subcell

  !> @brief Pick the exit solution with the shortest travel time
  function pick_exit(this, particle) result(exit_soln)
    class(MethodSubcellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer(I4B) :: exit_soln
    ! local
    real(DP) :: dtmin

    exit_soln = 0
    dtmin = 1.0d+30

    if (this%exit_solutions(1)%status < 2) then
      exit_soln = 1 ! x
      dtmin = this%exit_solutions(1)%dt
    end if
    if (this%exit_solutions(2)%status < 2 .and. &
        this%exit_solutions(2)%dt < dtmin) then
      exit_soln = 2 ! y
      dtmin = this%exit_solutions(2)%dt
    end if
    if (this%exit_solutions(3)%status < 2 .and. &
        this%exit_solutions(3)%dt < dtmin) then
      exit_soln = 3 ! z
    end if

  end function pick_exit

  !> @brief Compute candidate exit solutions
  subroutine find_exits(this, particle, domain)
    class(MethodSubcellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(DomainType), intent(in) :: domain
    ! local
    real(DP) :: x0, y0, z0

    select type (domain)
    type is (SubcellRectType)
      ! Initial particle location in scaled subcell coordinates
      x0 = particle%x / domain%dx
      y0 = particle%y / domain%dy
      z0 = particle%z / domain%dz

      ! Calculate exit solutions for each coordinate direction
      this%exit_solutions = [ &
                            find_exit(domain%vx1, domain%vx2, domain%dx, x0), &
                            find_exit(domain%vy1, domain%vy2, domain%dy, y0), &
                            find_exit(domain%vz1, domain%vz2, domain%dz, z0) &
                            ]

      ! Set exit faces
      if (this%exit_solutions(1)%v < DZERO) then
        this%exit_solutions(1)%iboundary = 1
      else if (this%exit_solutions(1)%v > DZERO) then
        this%exit_solutions(1)%iboundary = 2
      end if
      if (this%exit_solutions(2)%v < DZERO) then
        this%exit_solutions(2)%iboundary = 3
      else if (this%exit_solutions(2)%v > DZERO) then
        this%exit_solutions(2)%iboundary = 4
      end if
      if (this%exit_solutions(3)%v < DZERO) then
        this%exit_solutions(3)%iboundary = 5
      else if (this%exit_solutions(3)%v > DZERO) then
        this%exit_solutions(3)%iboundary = 6
      end if
    end select
  end subroutine find_exits

  !> @brief Find an exit solution for one dimension
  function find_exit(v1, v2, dx, xL) result(solution)
    ! dummy
    real(DP), intent(in) :: v1
    real(DP), intent(in) :: v2
    real(DP), intent(in) :: dx
    real(DP), intent(in) :: xL
    type(LinearExitSolutionType) :: solution

    solution = LinearExitSolutionType()
    solution%status = calculate_dt(v1, v2, dx, xL, &
                                   solution%v, solution%dvdx, solution%dt)
  end function find_exit

  !> @brief Calculate particle travel time to exit and exit status.
  !!
  !! This subroutine consists partly of code written by and/or adapted from
  !! David W. Pollock of the USGS for MODPATH 7. The authors of the present
  !! code are responsible for its appropriate application in this context
  !! and for any modifications or errors.
  !<
  function calculate_dt(v1, v2, dx, xL, v, dvdx, dt) result(status)
    ! dummy
    real(DP) :: v1
    real(DP) :: v2
    real(DP) :: dx
    real(DP) :: xL
    real(DP) :: v
    real(DP) :: dvdx
    real(DP) :: dt
    ! result
    integer(I4B) :: status
    ! local
    real(DP) :: v2a
    real(DP) :: v1a
    real(DP) :: dv
    real(DP) :: dva
    real(DP) :: vv
    real(DP) :: vvv
    real(DP) :: zro
    real(DP) :: zrom
    real(DP) :: x
    real(DP) :: tol
    real(DP) :: vr1
    real(DP) :: vr2
    real(DP) :: vr
    real(DP) :: v1v2
    logical(LGP) :: noOutflow

    ! Initialize variables.
    status = -1
    dt = 1.0d+20
    v2a = v2
    if (v2a .lt. DZERO) v2a = -v2a
    v1a = v1
    if (v1a .lt. DZERO) v1a = -v1a
    dv = v2 - v1
    dva = dv
    if (dva .lt. DZERO) dva = -dva

    ! Check for a uniform zero velocity in this direction.
    ! If so, set status = 2 and return (dt = 1.0d+20).
    tol = 1.0d-15
    if ((v2a .lt. tol) .and. (v1a .lt. tol)) then
      v = DZERO
      dvdx = DZERO
      status = NO_EXIT_STATIONARY
      return
    end if

    ! Check for uniform non-zero velocity in this direction.
    ! If so, set compute dt using the constant velocity,
    ! set status = 1 and return.
    vv = v1a
    if (v2a .gt. vv) vv = v2a
    vvv = dva / vv
    if (vvv .lt. 1.0d-4) then
      zro = tol
      zrom = -zro
      v = v1
      x = xL * dx
      if (v1 .gt. zro) dt = (dx - x) / v1
      if (v1 .lt. zrom) dt = -x / v1
      dvdx = DZERO
      status = OK_EXIT_CONSTANT
      return
    end if

    ! Velocity has a linear variation.
    ! Compute velocity corresponding to particle position.
    dvdx = dv / dx
    v = (DONE - xL) * v1 + xL * v2

    ! If flow is into the cell from both sides there is no outflow.
    ! In that case, set status = 3 and return.
    noOutflow = .true.
    if (v1 .lt. DZERO) noOutflow = .false.
    if (v2 .gt. DZERO) noOutflow = .false.
    if (noOutflow) then
      status = NO_EXIT_NO_OUTFLOW
      return
    end if

    ! If there is a divide in the cell for this flow direction, check to
    ! see if the particle is located exactly on the divide. If it is, move
    ! it very slightly to get it off the divide. This avoids possible
    ! numerical problems related to stagnation points.
    if ((v1 .le. DZERO) .and. (v2 .ge. DZERO)) then
      if (abs(v) .le. DZERO) then
        v = 1.0d-20
        if (v2 .le. DZERO) v = -v
      end if
    end if

    ! If there is a flow divide, this check finds out what side of the
    ! divide the particle is on and sets the value of vr appropriately
    ! to reflect that location.
    vr1 = v1 / v
    vr2 = v2 / v
    vr = vr1
    if (vr .le. DZERO) then
      vr = vr2
    end if

    ! If the product v1*v2 > 0, the velocity is in the same direction
    ! throughout the cell (i.e. no flow divide). If so, set the value
    ! of vr to reflect the appropriate direction.
    v1v2 = v1 * v2
    if (v1v2 .gt. DZERO) then
      if (v .gt. DZERO) vr = vr2
      if (v .lt. DZERO) vr = vr1
    end if

    ! Check if vr is (very close to) zero.
    ! If so, set status = 2 and return (dt = 1.0d+20).
    if (dabs(vr) .lt. 1.0d-10) then
      v = DZERO
      dvdx = DZERO
      status = NO_EXIT_STATIONARY
      return
    end if

    ! Compute travel time to exit face. Return with status = 0.
    dt = log(vr) / dvdx
    status = OK_EXIT

  end function calculate_dt

  !> @brief Update a cell-local coordinate based on a time increment.
  !!
  !! This subroutine consists partly or entirely of code written by
  !! David W. Pollock of the USGS for MODPATH 7. The authors of the present
  !! code are responsible for its appropriate application in this context
  !! and for any modifications or errors.
  !<
  pure function new_x(v, dvdx, v1, v2, dt, x, dx, velocity_profile) result(newx)
    ! dummy
    real(DP), intent(in) :: v
    real(DP), intent(in) :: dvdx
    real(DP), intent(in) :: v1
    real(DP), intent(in) :: v2
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: x
    real(DP), intent(in) :: dx
    logical(LGP), intent(in), optional :: velocity_profile
    ! result
    real(DP) :: newx
    logical(LGP) :: lprofile

    ! process optional arguments
    if (present(velocity_profile)) then
      lprofile = velocity_profile
    else
      lprofile = .false.
    end if

    ! recompute coordinate
    newx = x
    if (lprofile) then
      newx = newx + (v1 * dt / dx)
    else if (v .ne. DZERO) then
      newx = newx + (v * (exp(dvdx * dt) - DONE) / dvdx / dx)
    end if

    ! clamp to [0, 1]
    if (newx .lt. DZERO) newx = DZERO
    if (newx .gt. DONE) newx = DONE

  end function new_x

end module MethodSubcellPollockModule
