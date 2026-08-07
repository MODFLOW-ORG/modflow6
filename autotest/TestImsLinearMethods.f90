module TestImsLinearMethods
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: DZERO
  use testdrive, only: check, error_type, new_unittest, unittest_type
  use ImsLinearSettingsModule, only: resolve_ipc, ImsLinearSettingsType, &
                                     IPC_ILU0, IPC_MILU0, IPC_ILUT, IPC_MILUT
  use ImsLinearPeriodModule, only: ImsLinearPeriodType
  implicit none
  private
  public :: collect_imslinearmethods

contains

  subroutine collect_imslinearmethods(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("resolve_ipc_ilu_family", test_resolve_ipc), &
                new_unittest("period_utility_inactive", test_period_inactive) &
                ]
  end subroutine collect_imslinearmethods

  !> @brief resolve_ipc maps the ILU controls to the expected preconditioner
  !<
  subroutine test_resolve_ipc(error)
    type(error_type), allocatable, intent(out) :: error
    !
    ! -- level == 0 -> ILU0, or MILU0 when relax > 0
    call check(error, resolve_ipc(0, DZERO) == IPC_ILU0)
    if (allocated(error)) return
    call check(error, resolve_ipc(0, 0.97_DP) == IPC_MILU0)
    if (allocated(error)) return
    !
    ! -- level > 0 -> ILUT, or MILUT when relax > 0
    call check(error, resolve_ipc(5, DZERO) == IPC_ILUT)
    if (allocated(error)) return
    call check(error, resolve_ipc(5, 0.97_DP) == IPC_MILUT)
  end subroutine test_resolve_ipc

  !> @brief A disabled period-settings utility is inactive and changes nothing
  !<
  subroutine test_period_inactive(error)
    type(error_type), allocatable, intent(out) :: error
    type(ImsLinearPeriodType) :: period
    type(ImsLinearSettingsType) :: settings
    logical(LGP) :: changed
    !
    ! -- inunit == 0 disables the utility; read_period must be a no-op
    call period%init(0, 0)
    call check(error,.not. period%active)
    if (allocated(error)) return
    call period%read_period(settings, 1, .true., .true., changed)
    call check(error,.not. changed)
    if (allocated(error)) return
    call period%da()
  end subroutine test_period_inactive

end module TestImsLinearMethods
