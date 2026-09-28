module TestBudget
  use, intrinsic :: ieee_arithmetic, only: ieee_next_after
  use KindModule, only: I4B, DP
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       unittest_type
  use BudgetModule, only: value_to_string
  implicit none
  private
  public :: collect_budget

  real(DP), parameter :: big = 9.99999d11
  real(DP), parameter :: small = 0.1d0

contains

  subroutine collect_budget(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("value_to_string_fixed_point", &
                             test_value_to_string_fixed_point), &
                new_unittest("value_to_string_two_digit_exponent", &
                             test_value_to_string_two_digit_exponent), &
                new_unittest("value_to_string_three_digit_exponent", &
                             test_value_to_string_three_digit_exponent) &
                ]
  end subroutine collect_budget

  subroutine check_writes_e(error, val)
    type(error_type), allocatable, intent(out) :: error
    real(DP), intent(in) :: val
    character(len=17) :: string
    real(DP) :: readval
    integer(I4B) :: ios

    call value_to_string(val, string, big, small)
    if (index(string, 'E') == 0) then
      call test_failed(error, "missing E character: '"//string//"'")
      return
    end if
    read (string, *, iostat=ios) readval
    call check(error, ios == 0, "could not read back: '"//string//"'")
    if (allocated(error)) return
    call check(error, abs(readval - val) <= 1.d-4 * abs(val), &
               "value mismatch: '"//string//"'")
  end subroutine check_writes_e

  subroutine test_value_to_string_fixed_point(error)
    type(error_type), allocatable, intent(out) :: error
    character(len=17) :: string

    call value_to_string(0.d0, string, big, small)
    call check(error, string == '           0.0000')
    if (allocated(error)) return
    call value_to_string(-123.45678d0, string, big, small)
    call check(error, string == '        -123.4568')
  end subroutine test_value_to_string_fixed_point

  subroutine test_value_to_string_two_digit_exponent(error)
    type(error_type), allocatable, intent(out) :: error
    character(len=17) :: string

    call value_to_string(1.5d-20, string, big, small)
    call check(error, string == '       1.5000E-20', "'"//string//"'")
    if (allocated(error)) return
    call value_to_string(-1.5d20, string, big, small)
    call check(error, string == '      -1.5000E+20', "'"//string//"'")
    if (allocated(error)) return
    call value_to_string(1.d-99, string, big, small)
    call check(error, string == '       1.0000E-99', "'"//string//"'")
    if (allocated(error)) return
    call value_to_string(9.99994d99, string, big, small)
    call check(error, string == '       9.9999E+99', "'"//string//"'")
  end subroutine test_value_to_string_two_digit_exponent

  subroutine test_value_to_string_three_digit_exponent(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: vals(11)
    integer(I4B) :: i

    vals = [8.5159d-100, 5.d-101, 1.d-100, &
            ieee_next_after(1.d-99, 0.d0), &
            9.99995d99, ieee_next_after(9.99995d99, 0.d0), &
            ieee_next_after(9.99995d99, 1.d100), 9.99996d99, &
            1.d100, 1.d-300, 1.d300]
    do i = 1, size(vals)
      call check_writes_e(error, vals(i))
      if (allocated(error)) return
      call check_writes_e(error, -vals(i))
      if (allocated(error)) return
    end do
  end subroutine test_value_to_string_three_digit_exponent

end module TestBudget
