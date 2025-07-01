module MethodCellModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DONE, DZERO
  use MethodModule, only: MethodType
  use ParticleModule, only: ParticleType
  use CellDefnModule, only: CellDefnType
  implicit none

  private
  public :: MethodCellType

  type, abstract, extends(MethodType) :: MethodCellType
  contains
    procedure, public :: check
  end type MethodCellType

contains

  !> @brief Check reporting/terminating conditions before tracking
  !! the particle across the cell.
  !!
  !! Check a number of conditions determining whether to continue
  !! tracking the particle or terminate it, as well as whether to
  !! record any output data as per selected reporting conditions.
  !<
  subroutine check(this, particle, cell_defn, tmax)
    ! modules
    use TdisModule, only: endofsimulation, totimc, totim
    use ParticleModule, only: TERM_WEAKSINK, TERM_NO_EXITS, &
                              TERM_STOPZONE, TERM_INACTIVE
    use ParticleEventModule, only: CELLEXIT, TERMINATE, &
                                   TIMESTEP, WEAKSINK, USERTIME
    ! dummy
    class(MethodCellType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(CellDefnType), pointer, intent(inout) :: cell_defn
    real(DP), intent(in) :: tmax
    ! local
    logical(LGP) :: dry_cell, dry_particle, no_exit_face, stop_zone, weak_sink
    integer(I4B) :: i
    real(DP) :: t, ttrackmax

    dry_cell = this%fmi%ibdgwfsat0(cell_defn%icell) == 0
    dry_particle = particle%z > cell_defn%top
    no_exit_face = cell_defn%inoexitface > 0
    stop_zone = cell_defn%izone > 0 .and. particle%istopzone == cell_defn%izone
    weak_sink = cell_defn%iweaksink > 0

    particle%izone = cell_defn%izone
    if (stop_zone) then
      call this%events%terminate(particle, status=TERM_STOPZONE)
      return
    end if

    if (no_exit_face .and. .not. dry_cell) then
      call this%events%terminate(particle, status=TERM_NO_EXITS)
      return
    end if

    if (weak_sink) then
      if (particle%istopweaksink > 0) then
        call this%events%terminate(particle, status=TERM_WEAKSINK)
        return
      else
        call this%events%weaksink(particle)
      end if
    end if

    if (dry_cell) then
      if (particle%idrymeth == 0) then
        ! drop to cell bottom. handled by pass
        ! to bottom method, nothing to do here
        no_exit_face = .false.
      else if (particle%idrymeth == 1) then
        ! stop
        call this%events%terminate(particle, status=TERM_INACTIVE)
        return
      else if (particle%idrymeth == 2) then
        ! stay
        particle%advancing = .false.
        no_exit_face = .false.

        ! we might report tracking times
        ! out of order here, but we want
        ! the particle termination event
        ! (if this is the last time step)
        ! to have the maximum tracking t,
        ! so we need to keep tabs on it.
        ttrackmax = totim

        ! update tracking time to time
        ! step end time and save record
        particle%ttrack = totim
        call this%events%timestep(particle)

        ! record user tracking times
        call this%tracktimes%advance()
        if (this%tracktimes%any()) then
          do i = this%tracktimes%selection(1), this%tracktimes%selection(2)
            t = this%tracktimes%times(i)
            if (t < totimc) cycle
            if (t >= tmax) exit
            particle%ttrack = t
            call this%events%usertime(particle)
            if (t > ttrackmax) ttrackmax = t
          end do
        end if

        ! terminate if last period/step
        if (endofsimulation) then
          particle%ttrack = ttrackmax
          call this%events%terminate(particle, status=TERM_NO_EXITS)
          return
        end if
      end if
    else if (dry_particle .and. this%name /= "passtobottom") then
      if (particle%idrymeth == 0) then
        ! drop to water table
        particle%z = cell_defn%top
        call this%events%cellexit(particle)
      else if (particle%idrymeth == 1) then
        ! stop
        call this%events%terminate(particle, status=TERM_INACTIVE)
        return
      else if (particle%idrymeth == 2) then
        ! stay
        particle%advancing = .false.
        no_exit_face = .false.

        ! we might report tracking times
        ! out of order here, but we want
        ! the particle termination event
        ! (if this is the last time step)
        ! to have the maximum tracking t,
        ! so we need to keep tabs on it.
        ttrackmax = totim

        ! update tracking time to time
        ! step end time and save record
        particle%ttrack = totim
        call this%events%timestep(particle)

        ! record user tracking times
        call this%tracktimes%advance()
        if (this%tracktimes%any()) then
          do i = this%tracktimes%selection(1), this%tracktimes%selection(2)
            t = this%tracktimes%times(i)
            if (t < totimc) cycle
            if (t >= tmax) exit
            particle%ttrack = t
            call this%events%usertime(particle)
            if (t > ttrackmax) ttrackmax = t
          end do
        end if
      end if
    end if

    if (no_exit_face) then
      particle%advancing = .false.
      particle%istatus = TERM_NO_EXITS
      call this%events%terminate(particle)
      return
    end if

  end subroutine check

end module MethodCellModule
