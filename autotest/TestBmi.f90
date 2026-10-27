!> @brief Tests of the BMI, run on the simulation in the working directory
!!
!! The suite runs the two-model simulation in .mf6minsim (one stress period
!! of 1 day with one time step) and is skipped when the working directory has
!! no mfsim.nam. A simulation cannot be initialized twice in one process, so
!! the whole lifecycle is one test.
!<
module TestBmi
  use testdrive, only: error_type, unittest_type, new_unittest, check, &
                       skip_test
  use iso_c_binding, only: c_int, c_double, c_char, c_ptr, c_loc
  use KindModule, only: DP, I4B
  use mf6bmiError, only: BMI_SUCCESS, BMI_FAILURE
  use mf6bmiUtil, only: LENUNITS, LENVARLOCATION, string_to_char_array, &
                        char_array_to_string, strlen
  use mf6bmi, only: bmi_initialize, bmi_update_until, bmi_finalize, &
                    get_current_time, get_end_time, get_time_units, &
                    get_var_units, get_var_location, get_value, &
                    get_value_at_indices, set_value_at_indices
  use MemoryManagerModule, only: mem_setptr
  use MemoryAliasModule, only: mem_alias_deprecated
  use SimModule, only: count_warnings
  implicit none
  private
  public :: collect_bmi, c_string, f_string

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

  !> @brief Convert a nul-terminated C string to a Fortran string
  function f_string(c_str) result(str)
    character(kind=c_char, len=1), intent(in) :: c_str(:)
    character(len=:), allocatable :: str
    str = char_array_to_string(c_str, strlen(c_str, size(c_str)))
  end function f_string

  !> @brief Run the simulation through the BMI and check each new function
  subroutine test_lifecycle(error)
    type(error_type), allocatable, intent(out) :: error
    character(kind=c_char, len=1) :: c_units(LENUNITS + 1)
    character(kind=c_char, len=1) :: c_location(LENVARLOCATION + 1)
    integer(c_int), target :: inds(3), int_vals(1)
    real(c_double), target :: x(5), k(5), spdis(15), vals(3), new_k(2)
    real(c_double), target :: x_old(5)
    real(DP), dimension(:), pointer, contiguous :: xptr
    real(c_double), dimension(5), parameter :: &
      k_expected = [1.0_DP, 2.0_DP, 1.0_DP, 3.0_DP, 1.0_DP]
    real(c_double) :: time
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
      ! -- TIME_UNITS days in the TDIS Package
      call check(error, get_time_units(c_units) == BMI_SUCCESS, &
                 "get_time_units")
      if (allocated(error)) exit run
      call check(error, f_string(c_units) == 'd', &
                 "time units '"//f_string(c_units)//"', expected 'd'")
      if (allocated(error)) exit run

      call check(error, get_var_units(c_string('LEFTMODEL/X'), c_units) &
                 == BMI_SUCCESS, "get_var_units")
      if (allocated(error)) exit run
      call check(error, f_string(c_units) == 'unknown', "variable units")
      if (allocated(error)) exit run
      call check(error, get_var_location(c_string('LEFTMODEL/X'), c_location) &
                 == BMI_SUCCESS, "get_var_location")
      if (allocated(error)) exit run
      call check(error, f_string(c_location) == 'node', "variable location")
      if (allocated(error)) exit run
      call check(error, get_var_units(c_string('LEFTMODEL/NOPE'), c_units) &
                 == BMI_FAILURE, "get_var_units of an unknown variable")
      if (allocated(error)) exit run

      ! -- the simulation ends at 1 day, so 2 days cannot be reached
      call check(error, bmi_update_until(2.0_DP) == BMI_FAILURE, &
                 "update_until after the end time")
      if (allocated(error)) exit run
      call check(error, get_current_time(time) == BMI_SUCCESS, &
                 "get_current_time")
      if (allocated(error)) exit run
      call check(error, time == 0.0_DP, "time advanced by a failed update_until")
      if (allocated(error)) exit run
      call check(error, bmi_update_until(1.0_DP) == BMI_SUCCESS, &
                 "update_until the end time")
      if (allocated(error)) exit run
      call check(error, get_current_time(time) == BMI_SUCCESS, &
                 "get_current_time")
      if (allocated(error)) exit run
      call check(error, time == 1.0_DP, "current time after update_until")
      if (allocated(error)) exit run
      call check(error, bmi_update_until(0.5_DP) == BMI_FAILURE, &
                 "update_until before the current time")
      if (allocated(error)) exit run

      ! -- a rank-1 double precision variable
      ptr = c_loc(x)
      call check(error, get_value(c_string('LEFTMODEL/X'), ptr) &
                 == BMI_SUCCESS, "get_value of X")
      if (allocated(error)) exit run
      inds = [4, 0, 2]
      call check(error, get_value_at_indices(c_string('LEFTMODEL/X'), &
                                             c_loc(vals), inds, 3) &
                 == BMI_SUCCESS, "get_value_at_indices of X")
      if (allocated(error)) exit run
      do i = 1, 3
        call check(error, vals(i) == x(inds(i) + 1), &
                   "get_value_at_indices of X differs from get_value")
        if (allocated(error)) exit run
      end do

      ! -- a rank-2 variable, Fortran shape (3, 5), indexed in C order
      ptr = c_loc(spdis)
      call check(error, get_value(c_string('LEFTMODEL/NPF/SPDIS'), ptr) &
                 == BMI_SUCCESS, "get_value of SPDIS")
      if (allocated(error)) exit run
      inds = [14, 7, 0]
      call check(error, get_value_at_indices(c_string('LEFTMODEL/NPF/SPDIS'), &
                                             c_loc(vals), inds, 3) &
                 == BMI_SUCCESS, "get_value_at_indices of SPDIS")
      if (allocated(error)) exit run
      do i = 1, 3
        call check(error, vals(i) == spdis(inds(i) + 1), &
                   "get_value_at_indices of SPDIS differs from get_value")
        if (allocated(error)) exit run
      end do

      ! -- an integer variable: the model shape is (1, 1, 5)
      inds(1) = 2
      call check(error, get_value_at_indices(c_string('LEFTMODEL/DIS/MSHAPE'), &
                                             c_loc(int_vals), inds, 1) &
                 == BMI_SUCCESS, "get_value_at_indices of MSHAPE")
      if (allocated(error)) exit run
      call check(error, int_vals(1) == 5, "MSHAPE")
      if (allocated(error)) exit run

      ! -- an index past the end of the variable
      inds(1) = 5
      call check(error, get_value_at_indices(c_string('LEFTMODEL/X'), &
                                             c_loc(vals), inds, 1) &
                 == BMI_FAILURE, "get_value_at_indices out of range")
      if (allocated(error)) exit run

      ! -- set two of the five horizontal hydraulic conductivities
      inds(1:2) = [1, 3]
      new_k = [2.0_DP, 3.0_DP]
      call check(error, set_value_at_indices(c_string('LEFTMODEL/NPF/K11'), &
                                             inds, 2, c_loc(new_k)) &
                 == BMI_SUCCESS, "set_value_at_indices of K11")
      if (allocated(error)) exit run
      ptr = c_loc(k)
      call check(error, get_value(c_string('LEFTMODEL/NPF/K11'), ptr) &
                 == BMI_SUCCESS, "get_value of K11")
      if (allocated(error)) exit run
      do i = 1, 5
        call check(error, k(i) == k_expected(i), &
                   "K11 after set_value_at_indices")
        if (allocated(error)) exit run
      end do

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
