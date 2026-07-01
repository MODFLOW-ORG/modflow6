module mf6api
  use Mf6CoreModule, only: Mf6PrepareRetryLoop, &
                           Mf6StartRetry, &
                           Mf6FinishRetry
  use mf6bmiError
  use iso_c_binding, only: c_bool, c_int
  use KindModule, only: LGP
  implicit none

contains

  !> @brief Prepare the retry loop
  !<
  function api_prepare_retryloop() result(bmi_status) &
    bind(C, name="prepare_retryloop")
    !DIR$ ATTRIBUTES DLLEXPORT :: api_prepare_retryloop
    ! -- dummy variables
    integer(kind=c_int) :: bmi_status !< BMI status code

    call Mf6PrepareRetryLoop()
    bmi_status = BMI_SUCCESS

  end function api_prepare_retryloop

  !> @brief Call this before the solve
  !<
  function api_start_retry() result(bmi_status) &
    bind(C, name="start_retry")
    !DIR$ ATTRIBUTES DLLEXPORT :: api_start_retry
    ! -- dummy variables
    integer(kind=c_int) :: bmi_status !< BMI status code

    call Mf6StartRetry()
    bmi_status = BMI_SUCCESS

  end function api_start_retry

  !> @brief Call this after the solve, when true the retry
  !< loop should be exited and the timestep is finished
  function api_finish_retry(finish_retry) result(bmi_status) &
    bind(C, name="finish_retry")
    !DIR$ ATTRIBUTES DLLEXPORT :: api_finish_retry
    ! -- dummy variables
    integer(kind=c_int) :: bmi_status !< BMI status code
    logical(kind=c_bool), intent(out) :: finish_retry

    finish_retry = Mf6FinishRetry()

    bmi_status = BMI_SUCCESS

  end function api_finish_retry

end module mf6api
