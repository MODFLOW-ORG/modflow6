module MethodSubcellPollockModule
  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use MethodModule, only: MethodType
  use SubcellRectModule, only: SubcellRectType, create_subcell_rect
  use ParticleModule, only: ParticleType
  use PrtFmiModule, only: PrtFmiType
  use BaseDisModule, only: DisBaseType
  use CellModule, only: CellType
  use ConstantsModule, only: DZERO, DONE
  implicit none
  private
  public :: MethodSubcellPollockType
  public :: create_method_subcell_pollock
  public :: calculate_dt

  !> @brief Rectangular subcell tracking method
  type, extends(MethodType) :: MethodSubcellPollockType
    private
    real(DP), allocatable, public :: qextl1(:), qextl2(:), qintl(:) !< external and internal subcell flows
  contains
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
    real(DP) :: xOrigin
    real(DP) :: yOrigin
    real(DP) :: zOrigin
    real(DP) :: sinrot
    real(DP) :: cosrot

    select type (subcell => this%subcell)
    type is (SubcellRectType)
      ! Transform particle position into local subcell coordinates,
      !    track particle across subcell, convert back to model coords
      !    (sinrot and cosrot should be 0 and 1, respectively, i.e. no
      !    rotation, also no z translation; only x and y translations)
      xOrigin = subcell%xOrigin
      yOrigin = subcell%yOrigin
      zOrigin = subcell%zOrigin
      sinrot = subcell%sinrot
      cosrot = subcell%cosrot
      call particle%transform(xOrigin, yOrigin)
      call this%track_subcell(subcell, particle, tmax)
      call particle%transform(xOrigin, yOrigin, invert=.true.)
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
    use ParticleModule, only: ACTIVE, TERM_NO_EXITS_SUB, LVL_SUBCELL
    ! dummy
    class(MethodSubcellPollockType), intent(inout) :: this
    class(SubcellRectType), intent(in) :: subcell
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local
    real(DP) :: vx
    real(DP) :: dvxdx
    real(DP) :: vy
    real(DP) :: dvydy
    real(DP) :: vz
    real(DP) :: dvzdz
    real(DP) :: dtexitx
    real(DP) :: dtexity
    real(DP) :: dtexitz
    real(DP) :: dtexit
    real(DP) :: texit
    real(DP) :: dt
    real(DP) :: t
    real(DP) :: t0
    real(DP) :: x
    real(DP) :: y
    real(DP) :: z
    integer(I4B) :: statusVX
    integer(I4B) :: statusVY
    integer(I4B) :: statusVZ
    integer(I4B) :: i
    real(DP) :: initialX
    real(DP) :: initialY
    real(DP) :: initialZ
    integer(I4B) :: exitFace
    integer(I4B) :: reason

    reason = -1

    ! Initial particle location in scaled subcell coordinates
    initialX = particle%x / subcell%dx
    initialY = particle%y / subcell%dy
    initialZ = particle%z / subcell%dz

    ! Compute time of travel to each possible exit face
    statusVX = calculate_dt(subcell%vx1, subcell%vx2, subcell%dx, &
                            initialX, vx, dvxdx, dtexitx)
    statusVY = calculate_dt(subcell%vy1, subcell%vy2, subcell%dy, &
                            initialY, vy, dvydy, dtexity)
    statusVZ = calculate_dt(subcell%vz1, subcell%vz2, subcell%dz, &
                            initialZ, vz, dvzdz, dtexitz)

    !  Subcell has no exit face, terminate the particle
    !  todo: after initial release, consider ramifications
    if ((statusVX .eq. 3) .and. (statusVY .eq. 3) .and. (statusVZ .eq. 3)) then
      particle%istatus = TERM_NO_EXITS_SUB
      particle%advancing = .false.
      call this%save(particle, reason=3)
      return
    end if

    ! Determine (earliest) exit face and corresponding travel time to exit
    exitFace = 0
    dtexit = 1.0d+30
    if ((statusVX .lt. 2) .or. (statusVY .lt. 2) .or. (statusVZ .lt. 2)) then
      ! Consider x-oriented faces
      dtexit = dtexitx
      if (vx .lt. DZERO) then
        exitFace = 1
      else if (vx .gt. 0) then
        exitFace = 2
      end if
      ! Consider y-oriented faces
      if (dtexity .lt. dtexit) then
        dtexit = dtexity
        if (vy .lt. DZERO) then
          exitFace = 3
        else if (vy .gt. DZERO) then
          exitFace = 4
        end if
      end if
      ! Consider z-oriented faces
      if (dtexitz .lt. dtexit) then
        dtexit = dtexitz
        if (vz .lt. DZERO) then
          exitFace = 5
        else if (vz .gt. DZERO) then
          exitFace = 6
        end if
      end if
    else
    end if

    ! Compute exit time
    texit = particle%ttrack + dtexit
    t0 = particle%ttrack

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
        x = new_x(vx, dvxdx, subcell%vx1, subcell%vx2, &
                  dt, initialX, subcell%dx, statusVX == 1)
        y = new_x(vy, dvydy, subcell%vy1, subcell%vy2, &
                  dt, initialY, subcell%dy, statusVY == 1)
        z = new_x(vz, dvzdz, subcell%vz1, subcell%vz2, &
                  dt, initialZ, subcell%dz, statusVZ == 1)
        particle%x = x * subcell%dx
        particle%y = y * subcell%dy
        particle%z = z * subcell%dz
        particle%ttrack = t
        particle%istatus = ACTIVE
        call this%save(particle, reason=5)
      end do
    end if

    if (texit .gt. tmax) then
      ! The computed exit time is greater than the maximum time, so set
      ! final time for particle trajectory equal to maximum time and
      ! calculate particle location at that final time.
      t = tmax
      dt = t - t0
      x = new_x(vx, dvxdx, subcell%vx1, subcell%vx2, &
                dt, initialX, subcell%dx, statusVX == 1)
      y = new_x(vy, dvydy, subcell%vy1, subcell%vy2, &
                dt, initialY, subcell%dy, statusVY == 1)
      z = new_x(vz, dvzdz, subcell%vz1, subcell%vz2, &
                dt, initialZ, subcell%dz, statusVZ == 1)
      exitFace = 0
      particle%istatus = ACTIVE
      particle%advancing = .false.
      reason = 2 ! timestep end
    else
      ! The computed exit time is less than or equal to the maximum time,
      ! so set final time for particle trajectory equal to exit time and
      ! calculate exit location.
      t = texit
      dt = dtexit
      if ((exitFace .eq. 1) .or. (exitFace .eq. 2)) then
        x = DZERO
        y = new_x(vy, dvydy, subcell%vy1, subcell%vy2, &
                  dt, initialY, subcell%dy, statusVY == 1)
        z = new_x(vz, dvzdz, subcell%vz1, subcell%vz2, &
                  dt, initialZ, subcell%dz, statusVZ == 1)
        if (exitFace .eq. 2) x = DONE
      else if ((exitFace .eq. 3) .or. (exitFace .eq. 4)) then
        x = new_x(vx, dvxdx, subcell%vx1, subcell%vx2, dt, &
                  initialX, subcell%dx, statusVX == 1)
        y = DZERO
        z = new_x(vz, dvzdz, subcell%vz1, subcell%vz2, dt, &
                  initialZ, subcell%dz, statusVZ == 1)
        if (exitFace .eq. 4) y = DONE
      else if ((exitFace .eq. 5) .or. (exitFace .eq. 6)) then
        x = new_x(vx, dvxdx, subcell%vx1, subcell%vx2, &
                  dt, initialX, subcell%dx, statusVX == 1)
        y = new_x(vy, dvydy, subcell%vy1, subcell%vy2, &
                  dt, initialY, subcell%dy, statusVY == 1)
        z = DZERO
        if (exitFace .eq. 6) z = DONE
      else
        print *, "programmer error, invalid exit face", exitFace
        call pstop(1)
      end if
      reason = 1 ! cell transition
    end if

    ! Set final particle location in local (unscaled) subcell coordinates,
    ! final time for particle trajectory, and exit face
    particle%x = x * subcell%dx
    particle%y = y * subcell%dy
    particle%z = z * subcell%dz
    particle%ttrack = t
    particle%iboundary(LVL_SUBCELL) = exitFace

    ! Save particle track record
    if (reason >= 0) call this%save(particle, reason=reason)

  end subroutine track_subcell

  !> @brief Calculate particle travel time to exit and exit status.
  !!
  !! This subroutine consists partly or entirely of code written by
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
      status = 2
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
      status = 1
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
      status = 3
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
      status = 2
      return
    end if

    ! Compute travel time to exit face. Return with status = 0.
    dt = log(vr) / dvdx
    status = 0

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
