module MethodCellModule

  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use ConstantsModule, only: DONE, DZERO
  use MethodModule, only: MethodType, LEVEL_FEATURE
  use ParticleModule, only: ParticleType, ACTIVE, TERM_NO_EXITS, TERM_BOUNDARY
  use ParticleEventModule, only: ParticleEventType
  use CellExitEventModule, only: CellExitEventType
  use CellDefnModule, only: CellDefnType
  use IteratorModule, only: IteratorType
  implicit none

  private
  public :: MethodCellType

  type, abstract, extends(MethodType) :: MethodCellType
  contains
    procedure, public :: assess
    procedure, public :: cellexit
    procedure, public :: forms_cycle
    procedure, public :: store_event
    procedure, public :: get_level
    procedure, public :: try_pass
  end type MethodCellType

contains
  !> @brief Try passing the particle to the next subdomain.
  subroutine try_pass(this, particle, nextlevel, advancing)
    class(MethodCellType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer(I4B) :: nextlevel, ic, iface
    logical(LGP) :: advancing, on_face, on_top_face, partly_sat

    if (.not. particle%advancing) then
      advancing = .false.
      particle%iboundary = 0
      return
    end if

    call this%pass(particle)

    ic = particle%itrdomain(LEVEL_FEATURE)
    iface = particle%iboundary(LEVEL_FEATURE) - 1 ! cell is closed
    on_face = iface >= 0
    on_top_face = this%fmi%max_faces == iface
    partly_sat = this%fmi%gwfsat(this%cell%defn%icell) < DONE

    ! if at top and the cell is partially saturated,
    ! we're at the water table, not the cell top, so
    ! no exit event. if cell top is a boundary face,
    ! terminate. if the dry tracking method is stop,
    ! terminate. otherwise, leave the particle where
    ! it is- keep tracking it on the next time step.
    if (on_top_face .and. partly_sat) then
      advancing = .false.
      particle%advancing = .false.
      if (this%fmi%is_net_out_boundary_face(ic, iface) .or. &
          particle%idrymeth == 1) & ! dry_tracking_method stop
        call this%terminate(particle, status=TERM_BOUNDARY)
      return
    end if

    if (on_face) then
      advancing = .false.
      call this%cellexit(particle)
    end if
  end subroutine try_pass

  !> @brief Check reporting/terminating conditions before tracking
  !! the particle across the cell.
  !!
  !! Check a number of conditions determining whether to continue
  !! tracking the particle or terminate it, as well as whether to
  !! record any output data as per selected reporting conditions.
  !<
  subroutine assess(this, particle, cell_defn, tmax)
    ! modules
    use TdisModule, only: endofsimulation, totimc, totim
    use ParticleModule, only: TERM_WEAKSINK, TERM_NO_EXITS, &
                              TERM_STOPZONE, TERM_INACTIVE
    use ParticleEventModule, only: TERMINATE, TIMESTEP, WEAKSINK, USERTIME
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
      call this%terminate(particle, status=TERM_STOPZONE)
      return
    end if

    if (no_exit_face .and. .not. dry_cell) then
      call this%terminate(particle, status=TERM_NO_EXITS)
      return
    end if

    if (weak_sink) then
      if (particle%istopweaksink > 0) then
        call this%terminate(particle, status=TERM_WEAKSINK)
        return
      else
        call this%weaksink(particle)
      end if
    end if

    if (dry_cell) then
      if (particle%idrymeth == 0) then
        ! drop to cell bottom. handled by pass
        ! to bottom method, nothing to do here
        no_exit_face = .false.
      else if (particle%idrymeth == 1) then
        ! stop
        call this%terminate(particle, status=TERM_INACTIVE)
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
        call this%timestep(particle)

        ! record user tracking times
        call this%tracktimes%advance()
        if (this%tracktimes%any()) then
          do i = this%tracktimes%selection(1), this%tracktimes%selection(2)
            t = this%tracktimes%times(i)
            if (t < totimc) cycle
            if (t >= tmax) exit
            particle%ttrack = t
            call this%usertime(particle)
            if (t > ttrackmax) ttrackmax = t
          end do
        end if

        ! terminate if last period/step
        if (endofsimulation) then
          particle%ttrack = ttrackmax
          call this%terminate(particle, status=TERM_NO_EXITS)
          return
        end if
      end if
    else if (dry_particle .and. this%name /= "passtobottom") then
      if (particle%idrymeth == 0) then
        ! drop to water table
        particle%z = cell_defn%top
      else if (particle%idrymeth == 1) then
        ! stop
        call this%terminate(particle, status=TERM_INACTIVE)
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
        call this%timestep(particle)

        ! record user tracking times
        call this%tracktimes%advance()
        if (this%tracktimes%any()) then
          do i = this%tracktimes%selection(1), this%tracktimes%selection(2)
            t = this%tracktimes%times(i)
            if (t < totimc) cycle
            if (t >= tmax) exit
            particle%ttrack = t
            call this%usertime(particle)
            if (t > ttrackmax) ttrackmax = t
          end do
        end if
      end if
    end if

    if (no_exit_face) then
      particle%advancing = .false.
      particle%istatus = TERM_NO_EXITS
      call this%terminate(particle)
      return
    end if

  end subroutine assess

  !> @brief Particle exits a cell.
  subroutine cellexit(this, particle)
    class(MethodCellType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer :: event
    ! local
    integer(I4B) :: i, nhist
    class(*), pointer :: prev

    allocate (CellExitEventType :: event)
    select type (event)
    type is (CellExitEventType)
      event%exit_face = particle%iboundary(LEVEL_FEATURE)
    end select
    call this%events%dispatch(particle, event)
    if (particle%icycwin == 0) then
      deallocate (event)
      return
    end if
    if (this%forms_cycle(particle, event)) then
      ! print event history
      print *, "Cyclic pathline detected"
      nhist = particle%history%Count()
      do i = 1, nhist
        prev => particle%history%GetItem(i)
        select type (prev)
        class is (ParticleEventType)
          print *, "Back ", nhist - i + 1, ": ", prev%get_text()
        end select
      end do
      print *, "Current :", event%get_text()
      call pstop(1, 'Cyclic pathline detected, aborting')
    else
      call this%store_event(particle, event)
    end if
  end subroutine cellexit

  !> @brief Check if the event forms a cycle in the particle path.
  function forms_cycle(this, particle, event) result(found_cycle)
    ! dummy
    class(MethodCellType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    ! local
    class(IteratorType), allocatable :: itr
    logical(LGP) :: found_cycle

    found_cycle = .false.
    select type (event)
    type is (CellExitEventType)
      itr = particle%history%Iterator()
      do while (itr%has_next())
        call itr%next()
        select type (prev => itr%value())
        class is (CellExitEventType)
          if (event%icu == prev%icu .and. &
              event%ilay == prev%ilay .and. &
              event%izone == prev%izone .and. &
              event%exit_face == prev%exit_face .and. &
              event%exit_face /= 0) then
            found_cycle = .true.
            exit
          end if
        end select
      end do
    end select
  end function forms_cycle

  !> @brief Save the event in the particle's history.
  !! Acts like a queue, the oldest event is removed
  !! when the event count exceeds the maximum size.
  subroutine store_event(this, particle, event)
    ! dummy
    class(MethodCellType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    ! local
    class(*), pointer :: p

    select type (event)
    type is (CellExitEventType)
      p => event
      call particle%history%Add(p)
      if (particle%history%Count() > particle%icycwin) &
        call particle%history%RemoveNode(1, .true.)
    end select
  end subroutine store_event

  !> @brief Get the cell method's level.
  function get_level(this) result(level)
    class(MethodCellType), intent(in) :: this
    integer(I4B) :: level
    level = LEVEL_FEATURE
  end function get_level

end module MethodCellModule
