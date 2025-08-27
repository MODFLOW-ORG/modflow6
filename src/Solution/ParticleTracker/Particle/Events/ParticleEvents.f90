module ParticleEventsModule
  use KindModule, only: DP, I4B, LGP
  use ListModule, only: ListType
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType
  implicit none

  private

  type, public, abstract :: ParticleEventConsumerType
  contains
    procedure(handle_event), deferred :: handle_event
  end type ParticleEventConsumerType

  type, public :: ParticleEventDispatcherType
    type(ListType) :: consumers
  contains
    procedure, public :: subscribe
    procedure, public :: dispatch
    procedure :: destroy
  end type ParticleEventDispatcherType

  abstract interface
    subroutine handle_event(this, particle, event)
      import ParticleEventConsumerType, ParticleType, ParticleEventType
      class(ParticleEventConsumerType), intent(inout) :: this
      type(ParticleType), pointer, intent(in) :: particle
      class(ParticleEventType), pointer, intent(in) :: event
    end subroutine handle_event
  end interface

contains
  !> @brief Subscribe a consumer to the dispatcher.
  subroutine subscribe(this, consumer)
    class(ParticleEventDispatcherType), intent(inout) :: this
    class(ParticleEventConsumerType), target, intent(inout) :: consumer
    class(*), pointer :: p
    p => consumer
    call this%consumers%Add(p)
  end subroutine subscribe

  !> @brief Dispatch an event.
  subroutine dispatch(this, particle, event)
    use TdisModule, only: kper, kstp, totimc
    ! dummy
    class(ParticleEventDispatcherType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(inout) :: event
    ! local
    integer(I4B) :: i, per, stp
    real(DP) :: x, y, z
    class(*), pointer :: p

    ! If tracking time falls exactly on a boundary between time steps,
    ! report the previous time step for this datum. This is to follow
    ! MP7's behavior, and because the particle will have been tracked
    ! up to this instant under the previous time step's conditions, so
    ! the time step we're about to start shouldn't get "credit" for it.
    per = kper
    stp = kstp
    if (particle%ttrack == totimc .and. (per > 1 .or. stp > 1)) then
      if (stp > 1) then
        stp = stp - 1
      else if (per > 1) then
        per = per - 1
        stp = 1
      end if
    end if

    ! Convert to model coordinates if we need to
    x = particle%x
    y = particle%y
    z = particle%z
    call particle%get_model_coords(x, y, z)

    event%kper = per
    event%kstp = stp
    event%imdl = particle%imdl
    event%iprp = particle%iprp
    event%irpt = particle%irpt
    event%ilay = particle%ilay
    event%icu = particle%icu
    event%izone = particle%izone
    event%trelease = particle%trelease
    event%ttrack = particle%ttrack
    event%x = x
    event%y = y
    event%z = z
    event%istatus = particle%istatus

    do i = 1, this%consumers%Count()
      p => this%consumers%GetItem(i)
      select type (consumer => p)
      class is (ParticleEventConsumerType)
        call consumer%handle_event(particle, event)
      end select
    end do
  end subroutine dispatch

  !> @brief Destroy the dispatcher.
  subroutine destroy(this)
    class(ParticleEventDispatcherType), intent(inout) :: this
    call this%consumers%Clear()
  end subroutine destroy

end module ParticleEventsModule
