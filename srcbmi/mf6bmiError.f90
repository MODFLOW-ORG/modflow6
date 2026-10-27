!> @brief Detailed error information for the BMI
!!
!! This module contains error codes and detailed error
!! messages (as format strings) for the BMI/XMI.
!<
module mf6bmiError
  use iso_c_binding, only: c_char, c_int, C_NULL_CHAR
  use KindModule, only: I4B
  use ConstantsModule, only: MAXCHARLEN
  use SimVariablesModule, only: istdout

  integer, parameter :: BMI_FAILURE = 1 !< BMI status code for failure (taken from bmi.f90, CSDMS)
  integer, parameter :: BMI_SUCCESS = 0 !< BMI status code for success (taken from bmi.f90, CSDMS)

  integer(I4B), parameter :: LENERRMESSAGE = 1024 !< max length for the error message
  integer(c_int), bind(C, name="BMI_LENERRMESSAGE") :: BMI_LENERRMESSAGE = &
                                                       LENERRMESSAGE + 1 !< max. length for the (exported) C-style error message
  !DIR$ ATTRIBUTES DLLEXPORT :: BMI_LENERRMESSAGE

  character(len=LENERRMESSAGE) :: bmi_last_error = 'No BMI error reported' !< module variable containing the last error as a Fortran string

  character(len=*), parameter :: fmt_general_err = & !< General bmi error, args: context/detail
                                 "('BMI Error, ', a)"
  character(len=*), parameter :: fmt_unknown_var = & !< Variable unknown, args: variable name, memory path
                                 "('BMI Error, unknown variable: ', a, ' at ', a)"
  character(len=*), parameter :: fmt_invalid_var = & !< Invalid variable address, args: variable address
                                 "('BMI Error, invalid address string: ', a)"
  character(len=*), parameter :: fmt_unsupported_rank = & !< Unsupported rank, args: variable name
                                 "('BMI Error, unsupported rank for variable: &
                                 &', a)"
  character(len=*), parameter :: fmt_unsupported_type = & !< Unsupported type, args: variable name
                                 "('BMI Error, unsupported type for variable: &
                                 &', a)"
  character(len=*), parameter :: fmt_invalid_mem_access = & !< Invalid memory access, args: variable name
                                 "('Fatal BMI Error, invalid access of memory &
                                 &for variable: ', a)"
  character(len=*), parameter :: fmt_fail_cvg_sol = & !< Solution failed to converge, args: detail
                                 "('BMI Error, Numerical Solution ', i3, &
                                 &' failed to converge')"

contains

  !> @brief Sets the last BMI error message and copies
    !! it to an exported C-string
  !<
  subroutine report_bmi_error(err_msg)
    ! -- dummy variables
    character(len=*), intent(in) :: err_msg !< the error message
    bmi_last_error = err_msg
    write (istdout, *) trim(err_msg)
  end subroutine report_bmi_error

  !> @brief Get the last error in the BMI as a character array
    !! with size BMI_LENERRMESSAGE
  !<
  function get_last_bmi_error(c_error) result(bmi_status) &
    bind(C, name="get_last_bmi_error")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_last_bmi_error
    ! -- dummy variables
    character(kind=c_char, len=1), intent(out) :: c_error(BMI_LENERRMESSAGE) !< C style char array with error
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    integer(I4B) :: i, length

    length = len(trim(bmi_last_error))
    do i = 1, length
      c_error(i) = bmi_last_error(i:i)
    end do
    c_error(length + 1) = C_NULL_CHAR

    bmi_status = BMI_SUCCESS
  end function get_last_bmi_error

end module mf6bmiError
