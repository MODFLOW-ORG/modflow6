!> @brief Tests of the XMI, run on the simulation in the working directory
!!
!! The suite steps the two-model simulation in .mf6minsim through the XMI
!! and is skipped when the working directory has no mfsim.nam. A simulation
!! cannot be initialized twice in one process, so the whole lifecycle is one
!! test.
!<
module TestXmi
  use testdrive, only: error_type, unittest_type, new_unittest, check, &
                       skip_test
  use iso_c_binding, only: c_int, c_double, c_char, c_loc
  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENMEMADDRESS
  use mf6bmiError, only: BMI_SUCCESS
  use mf6bmi, only: bmi_initialize, bmi_finalize, get_current_time, &
                    get_time_step, get_value_at_indices
  use mf6xmi, only: xmi_get_version, get_var_address, &
                    xmi_get_subcomponent_count, xmi_prepare_time_step, &
                    xmi_prepare_solve, xmi_solve, xmi_finalize_solve, &
                    xmi_finalize_time_step
  use TestBmi, only: c_string, f_string
  implicit none
  private
  public :: collect_xmi

contains

  subroutine collect_xmi(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [new_unittest("lifecycle", test_lifecycle)]
  end subroutine collect_xmi

  !> @brief Step the simulation through the XMI and check the result
  subroutine test_lifecycle(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: max_iter = 100
    ! -- BMI_LENVERSION characters
    character(kind=c_char, len=1) :: c_version(256)
    character(kind=c_char, len=1) :: c_address(LENMEMADDRESS + 1)
    integer(c_int), target :: inds(1)
    integer(c_int) :: count, converged
    real(c_double), target :: head(1)
    real(c_double) :: dt, time
    integer(I4B) :: i, kiter
    logical :: exists

    inquire (file='mfsim.nam', exist=exists)
    if (.not. exists) then
      call skip_test(error, "no simulation in the working directory")
      return
    end if

    call check(error, bmi_initialize() == BMI_SUCCESS, "initialize")
    if (allocated(error)) return

    run: block
      call check(error, xmi_get_version(c_version) == BMI_SUCCESS, &
                 "get_version")
      if (allocated(error)) exit run
      call check(error, len(f_string(c_version)) > 0, "empty version")
      if (allocated(error)) exit run

      call check(error, get_var_address(c_string('LEFTMODEL'), &
                                        c_string('NPF'), c_string('K11'), &
                                        c_address) == BMI_SUCCESS, &
                 "get_var_address")
      if (allocated(error)) exit run
      call check(error, f_string(c_address) == 'LEFTMODEL/NPF/K11', &
                 "address '"//f_string(c_address)//"'")
      if (allocated(error)) exit run

      ! -- both models are in one solution
      call check(error, xmi_get_subcomponent_count(count) == BMI_SUCCESS, &
                 "get_subcomponent_count")
      if (allocated(error)) exit run
      call check(error, count == 1, "subcomponent count")
      if (allocated(error)) exit run

      ! -- one time step, solved outer iteration by outer iteration
      call check(error, get_time_step(dt) == BMI_SUCCESS, "get_time_step")
      if (allocated(error)) exit run
      call check(error, xmi_prepare_time_step(dt) == BMI_SUCCESS, &
                 "prepare_time_step")
      if (allocated(error)) exit run
      call check(error, xmi_prepare_solve(1) == BMI_SUCCESS, "prepare_solve")
      if (allocated(error)) exit run
      converged = 0
      do kiter = 1, max_iter
        call check(error, xmi_solve(1, converged) == BMI_SUCCESS, "solve")
        if (allocated(error)) exit run
        if (converged == 1) exit
      end do
      call check(error, converged == 1, "the solution did not converge")
      if (allocated(error)) exit run
      call check(error, xmi_finalize_solve(1) == BMI_SUCCESS, "finalize_solve")
      if (allocated(error)) exit run
      call check(error, xmi_finalize_time_step() == BMI_SUCCESS, &
                 "finalize_time_step")
      if (allocated(error)) exit run
      call check(error, get_current_time(time) == BMI_SUCCESS, &
                 "get_current_time")
      if (allocated(error)) exit run
      call check(error, time == 1.0_DP, "current time after the time step")
      if (allocated(error)) exit run

      ! -- the constant heads hold the two ends of the section
      inds(1) = 0
      call check(error, get_value_at_indices(c_string('LEFTMODEL/X'), &
                                             c_loc(head), inds, 1) &
                 == BMI_SUCCESS, "get_value_at_indices of the left head")
      if (allocated(error)) exit run
      call check(error, head(1) == 1.0_DP, "head at the left end")
      if (allocated(error)) exit run
      inds(1) = 4
      call check(error, get_value_at_indices(c_string('RIGHTMODEL/X'), &
                                             c_loc(head), inds, 1) &
                 == BMI_SUCCESS, "get_value_at_indices of the right head")
      if (allocated(error)) exit run
      call check(error, head(1) == 10.0_DP, "head at the right end")
      if (allocated(error)) exit run
    end block run

    if (allocated(error)) then
      i = bmi_finalize()
    else
      call check(error, bmi_finalize() == BMI_SUCCESS, "finalize")
    end if

  end subroutine test_lifecycle

end module TestXmi
