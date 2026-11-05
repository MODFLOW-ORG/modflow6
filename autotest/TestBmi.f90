!> @brief Tests of the BMI, run on the simulation in the working directory
!!
!! The suite runs the two-model simulation in .mf6minsim and is skipped when
!! the working directory has no mfsim.nam. A simulation cannot be initialized
!! twice in one process, so the whole lifecycle is one test.
!<
module TestBmi
  use testdrive, only: error_type, unittest_type, new_unittest, check, &
                       skip_test
  use iso_c_binding, only: c_double, c_char, c_ptr, c_loc
  use KindModule, only: DP, I4B
  use mf6bmiError, only: BMI_SUCCESS
  use mf6bmiUtil, only: string_to_char_array
  use mf6bmi, only: bmi_initialize, bmi_finalize, get_value
  use MemoryManagerModule, only: mem_setptr
  use MemoryAliasModule, only: mem_alias_deprecated
  use SimModule, only: count_warnings
  implicit none
  private
  public :: collect_bmi, c_string

contains

  subroutine collect_bmi(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [new_unittest("lifecycle", test_lifecycle)]
  end subroutine collect_bmi

  !> @brief Convert a Fortran string to a nul-terminated C string
  function c_string(str) result(c_str)
    character(len=*), intent(in) :: str
    character(kind=c_char, len=1) :: c_str(len(str) + 1)
    c_str = string_to_char_array(str, len(str))
  end function c_string

  !> @brief Run the simulation through the BMI
  subroutine test_lifecycle(error)
    type(error_type), allocatable, intent(out) :: error
    real(c_double), target :: x(5)
    real(c_double), target :: x_old(5)
    real(DP), dimension(:), pointer, contiguous :: xptr
    type(c_ptr) :: ptr
    integer(I4B) :: i
    integer(I4B) :: nwarn
    logical :: exists

    inquire (file='mfsim.nam', exist=exists)
    if (.not. exists) then
      call skip_test(error, "no simulation in the working directory")
      return
    end if

    call check(error, bmi_initialize() == BMI_SUCCESS, "initialize")
    if (allocated(error)) return

    run: block
      ! -- a deprecated alias resolves and is reported once, the current
      !    name is not reported
      call mem_setptr(xptr, 'X', 'LEFTMODEL')
      call mem_alias_deprecated(xptr, 'X_OLD', 'LEFTMODEL', 'X', 'LEFTMODEL', &
                                '6.9.0')
      nwarn = count_warnings()
      ptr = c_loc(x)
      call check(error, get_value(c_string('LEFTMODEL/X'), ptr) &
                 == BMI_SUCCESS, "get_value of X")
      if (allocated(error)) exit run
      call check(error, count_warnings() == nwarn, "X reported as deprecated")
      if (allocated(error)) exit run
      ptr = c_loc(x_old)
      do i = 1, 2
        call check(error, get_value(c_string('LEFTMODEL/X_OLD'), ptr) &
                   == BMI_SUCCESS, "get_value of X_OLD")
        if (allocated(error)) exit run
      end do
      call check(error, count_warnings() == nwarn + 1, &
                 "X_OLD not reported once as deprecated")
      if (allocated(error)) exit run
      do i = 1, 5
        call check(error, x_old(i) == x(i), "X_OLD differs from X")
        if (allocated(error)) exit run
      end do
    end block run

    if (allocated(error)) then
      i = bmi_finalize()
    else
      call check(error, bmi_finalize() == BMI_SUCCESS, "finalize")
    end if

  end subroutine test_lifecycle

end module TestBmi
