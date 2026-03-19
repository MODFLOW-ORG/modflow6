module ParticleEventsModule
  use KindModule, only: DP, I4B, LGP
  use ListModule, only: ListType
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType
  implicit none

  private
  public :: handle_event

  !> @brief Subscription to particle events: a procedure to
  !! handle the event with an unlimited pointer for storing
  !! arbitrary context the handling procedure may reference.
  type, public :: ParticleEventSubscription
    procedure(handle_event), pointer, nopass :: handler
    class(*), pointer :: context
  end type

  !> @brief Dispatcher for particle events. Consumers subscribe
  !! handlers to the dispatcher. Events may be dispatched, with
  !! the first handler to claim the event stopping propagation,
  !! or broadcast, with all handlers receiving the event.
  type, public :: ParticleEventDispatcherType
    type(ListType) :: subscriptions
  contains
    procedure, public :: subscribe
    procedure, public :: dispatch
    procedure, public :: broadcast
    procedure :: prep_event
    procedure :: destroy
  end type ParticleEventDispatcherType

  abstract interface
    !> @brief Event handler interface. Handlers may signal
    !! to the dispatching caller whether they have handled
    !! the event, but the signal is ignored for broadcasts.
    logical function handle_event(context, particle, event)
      import :: ParticleType, ParticleEventType
      class(*), pointer :: context
      type(ParticleType), pointer, intent(inout) :: particle
      class(ParticleEventType), pointer, intent(in) :: event
    end function
  end interface

contains
  !> @brief Add a subscription to the dispatcher.
  subroutine subscribe(this, handler, context)
    class(ParticleEventDispatcherType), intent(inout) :: this
    procedure(handle_event) :: handler
    class(*), pointer :: context
    ! local
    type(ParticleEventSubscription), pointer :: subscription
    class(*), pointer :: p

    allocate (subscription)
    subscription%handler => handler
    subscription%context => context
    p => subscription
    call this%subscriptions%Add(p)
  end subroutine subscribe

  !> @brief Prepare an event for dispatching, loading it with
  !! the current state of the particle. For internal use only.
  subroutine prep_event(this, particle, event)
    use TdisModule, only: kper, kstp
    ! dummy
    class(ParticleEventDispatcherType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    ! local
    real(DP) :: x, y, z

    x = particle%x
    y = particle%y
    z = particle%z
    ! switch to model coordinates if needed
    call particle%get_model_coords(x, y, z)

    event%kper = kper
    event%kstp = kstp
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
  end subroutine prep_event

  !> @brief Dispatch an event for claiming. The first
  !! subscriber to claim the event stops propagation.
  subroutine dispatch(this, particle, event)
    class(ParticleEventDispatcherType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    ! local
    logical(LGP) :: handled
    integer(I4B) :: i
    class(*), pointer :: p

    call this%prep_event(particle, event)

    do i = 1, this%subscriptions%Count()
      p => this%subscriptions%GetItem(i)
      select type (subscription => p)
      type is (ParticleEventSubscription)
        handled = subscription%handler( &
                  subscription%context, &
                  particle, &
                  event)
        if (handled) exit
      end select
    end do
  end subroutine dispatch

  !> @brief Broadcast an event to all subscribers.
  !! All get the event regardless of return value.
  subroutine broadcast(this, particle, event)
    class(ParticleEventDispatcherType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    ! local
    logical(LGP) :: handled
    integer(I4B) :: i
    class(*), pointer :: p

    call this%prep_event(particle, event)

    do i = 1, this%subscriptions%Count()
      p => this%subscriptions%GetItem(i)
      select type (subscription => p)
      type is (ParticleEventSubscription)
        handled = subscription%handler(subscription%context, particle, event)
      end select
    end do
  end subroutine broadcast

  !> @brief Destroy the dispatcher.
  subroutine destroy(this)
    class(ParticleEventDispatcherType), intent(inout) :: this
    call this%subscriptions%Clear(destroy=.true.)
  end subroutine destroy

end module ParticleEventsModule
