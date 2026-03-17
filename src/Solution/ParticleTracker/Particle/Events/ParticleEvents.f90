module ParticleEventsModule
  use KindModule, only: DP, I4B, LGP
  use ListModule, only: ListType
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType
  implicit none

  private
  public :: event_callback

  ! type, public, abstract :: ParticleEventConsumerType
  ! contains
  !   procedure(handle_event), deferred :: handle_event
  ! end type ParticleEventConsumerType

  type, public :: ParticleEventSubscription
    procedure(event_callback), pointer, nopass :: callback
    class(*), pointer :: context
  end type

  type, public :: ParticleEventDispatcherType
    type(ListType) :: subscriptions
  contains
    procedure, public :: subscribe
    procedure, public :: dispatch
    procedure :: destroy
  end type ParticleEventDispatcherType

  abstract interface
    ! subroutine handle_event(this, particle, event)
    !   import ParticleEventConsumerType, ParticleType, ParticleEventType
    !   class(ParticleEventConsumerType), intent(inout) :: this
    !   type(ParticleType), pointer, intent(in) :: particle
    !   class(ParticleEventType), pointer, intent(in) :: event
    ! end subroutine handle_event

    subroutine event_callback(context, particle, event)
      import :: ParticleType, ParticleEventType
      class(*), pointer :: context
      type(ParticleType), pointer, intent(inout) :: particle
      class(ParticleEventType), pointer, intent(in) :: event
    end subroutine
  end interface

contains
  !> @brief Subscribe a consumer to the dispatcher.
  subroutine subscribe(this, callback, context)
    ! dummy
    class(ParticleEventDispatcherType), intent(inout) :: this
    procedure(event_callback) :: callback
    class(*), pointer :: context
    ! local
    type(ParticleEventSubscription), pointer :: subscription
    class(*), pointer :: p

    allocate (subscription)
    subscription%callback => callback
    subscription%context => context
    p => subscription
    call this%subscriptions%Add(p)
  end subroutine subscribe

  !> @brief Dispatch an event.
  subroutine dispatch(this, particle, event)
    use TdisModule, only: kper, kstp
    ! dummy
    class(ParticleEventDispatcherType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    ! local
    integer(I4B) :: i
    real(DP) :: x, y, z
    class(*), pointer :: p

    ! Convert to model coordinates if we need to
    x = particle%x
    y = particle%y
    z = particle%z
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

    do i = 1, this%subscriptions%Count()
      p => this%subscriptions%GetItem(i)
      select type (subscription => p)
      type is (ParticleEventSubscription)
        call subscription%callback(subscription%context, particle, event)
      end select
    end do
  end subroutine dispatch

  !> @brief Destroy the dispatcher.
  subroutine destroy(this)
    class(ParticleEventDispatcherType), intent(inout) :: this
    call this%subscriptions%Clear(destroy=.true.)
  end subroutine destroy

end module ParticleEventsModule
