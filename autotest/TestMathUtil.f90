module TestMathUtil
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DNODATA, DZERO, DONE
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       to_string, unittest_type
  use MathUtilModule, only: f1d, is_close, mod_offset, &
                            zero_ch, zero_br, &
                            get_perturbation, &
                            arange, linspace
  implicit none
  private
  public :: collect_mathutil

contains

  subroutine collect_mathutil(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("is_close_symmetric", test_is_close_symmetric), &
                new_unittest("is_close_symmetric_near_0", &
                             test_is_close_symmetric_near_0), &
                new_unittest("mod_offset", &
                             test_mod_offset), &
                new_unittest("zero_ch", &
                             test_zero_ch), &
                new_unittest("zero_br", &
                             test_zero_br), &
                new_unittest("get_perturbation", &
                             test_get_perturbation), &
                new_unittest("arange", test_arange), &
                new_unittest("linspace", test_linspace) &
                ]
  end subroutine collect_mathutil

  subroutine test_mod_offset(error)
    type(error_type), allocatable, intent(out) :: error

    ! with no offset specified, should behave just like mod
    call check(error, mod_offset(2, 2) == 0)
    call check(error, mod_offset(2, 3) == 2)
    call check(error, mod_offset(2.0_DP, 2.0_DP) == 0.0_DP)
    call check(error, mod_offset(2.0_DP, 3.0_DP) == 2.0_DP)

    ! with offset d specified, if the result x = a mod n falls
    ! between 0 and n - 1, the new result x = a mod_d n falls
    ! between d and d + n - 1.
    call check(error, mod_offset(2, 3, -2) == -1)
    call check(error, mod_offset(2, 3, -1) == -1)
    call check(error, mod_offset(2, 3, 0) == 2)
    call check(error, mod_offset(2, 3, 1) == 2)
    call check(error, mod_offset(2, 3, 2) == 2)
    call check(error, mod_offset(2, 3, 3) == 5)
    call check(error, mod_offset(2, 3, 4) == 5)
    call check(error, mod_offset(2.0_DP, 3.0_DP, -1.0_DP) == -1.0_DP)
    call check(error, mod_offset(2.0_DP, 3.0_DP, 2.0_DP) == 2.0_DP)
    call check(error, mod_offset(2.0_DP, 3.0_DP, 3.0_DP) == 5.0_DP)
  end subroutine test_mod_offset

  subroutine test_is_close_symmetric(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: a, b, rtol

    ! exact match
    a = 1.0_DP
    b = 1.0_DP
    call check(error, is_close(a, b), &
               "exp eq: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", eps=default")
    if (allocated(error)) return

    ! mismatch with default rtol
    b = 1.0001_DP
    call check(error, (.not. (is_close(a, b))), &
               "exp ne: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", eps=default")
    if (allocated(error)) return

    ! inexact match with large rtol
    rtol = 1d-2
    call check(error, is_close(a, b, rtol=rtol), &
               "exp eq: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", rtol="//to_string(rtol))
    if (allocated(error)) return

    ! mismatch when we reduce rtol
    rtol = 0.5d-5
    call check(error, (.not. is_close(a, b, rtol=rtol)), &
               "exp ne: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", rtol="//to_string(rtol))
    if (allocated(error)) return

    ! +/-0
    call check(error, is_close(0.0_DP, -0.0_DP), &
               "exp ne: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", eps=default")

    ! DNODATA
    call check(error, (.not. is_close(0.0_DP, DNODATA)), &
               "exp ne: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", eps=default")
    call check(error, is_close(DNODATA, DNODATA), &
               "exp ne: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", eps=default")
    call check(error, (.not. is_close(DNODATA, DNODATA / 10)), &
               "exp ne: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", eps=default")
    call check(error, (.not. is_close(DNODATA, DNODATA * 10)), &
               "exp ne: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", eps=default")

  end subroutine test_is_close_symmetric

  subroutine test_is_close_symmetric_near_0(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: a, b, rtol, atol

    a = 0.0_DP
    b = 0.0_DP
    call check(error, is_close(a, b), &
               "exp eq: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", rtol=default")
    if (allocated(error)) return

    a = DZERO
    b = DZERO
    call check(error, is_close(a, b), &
               "exp eq: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", rtol=default")
    if (allocated(error)) return

    b = 1d-4
    call check(error, (.not. is_close(a, b)), &
               "exp eq: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", rtol=default")
    if (allocated(error)) return

    rtol = 0.999_DP
    call check(error, &
               ! expect failure, see above
               (.not. is_close(a, b, rtol=rtol)), &
               "exp eq: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", rtol="//to_string(rtol))
    if (allocated(error)) return

    ! absolute comparison is appropriate when a and/or b are near or equal to 0
    b = 1d-4
    atol = 1d-3
    call check(error, is_close(a, b, atol=atol), &
               "exp eq: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", atol="//to_string(atol))
    if (allocated(error)) return

    ! make sure the absolute tolerance is applied
    b = 1d-4
    atol = 1d-5
    call check(error, (.not. is_close(a, b, atol=atol)), &
               "exp eq: a="//to_string(a)// &
               ", b="//to_string(b)// &
               ", atol="//to_string(atol))
    if (allocated(error)) return

  end subroutine test_is_close_symmetric_near_0

  pure function sine(bet) result(s)
    real(DP), intent(in) :: bet
    real(DP) :: s
    s = sin(bet)
  end function sine

  subroutine test_zero_ch(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), parameter :: pi = 4 * atan(1.0_DP)
    real(DP) :: z
    procedure(f1d), pointer :: f

    f => sine

    z = zero_ch(-1.0_DP, 1.0_DP, f, 0.001_DP)
    call check(error, is_close(z, 0.0_DP, atol=1d-6), &
               'expected 0, got: '//to_string(z))
    if (allocated(error)) return

    z = zero_ch(-4.0_DP, -1.0_DP, f, 0.001_DP)
    call check(error, is_close(z, -pi, atol=1d-6), &
               'expected -pi, got: '//to_string(z))
    if (allocated(error)) return

    z = zero_ch(1.0_DP, 4.0_DP, f, 0.001_DP)
    call check(error, is_close(z, pi, atol=1d-6), &
               'expected pi, got: '//to_string(z))
    if (allocated(error)) return
  end subroutine test_zero_ch

  subroutine test_zero_br(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), parameter :: pi = 4 * atan(1.0_DP)
    real(DP) :: z
    procedure(f1d), pointer :: f

    f => sine

    z = zero_br(-1.0_DP, 1.0_DP, f, 0.001_DP)
    call check(error, is_close(z, 0.0_DP, atol=1d-6), &
               'expected 0, got: '//to_string(z))
    if (allocated(error)) return

    z = zero_br(-4.0_DP, -1.0_DP, f, 0.001_DP)
    call check(error, is_close(z, -pi, atol=1d-6), &
               'expected -pi, got: '//to_string(z))
    if (allocated(error)) return

    z = zero_br(1.0_DP, 4.0_DP, f, 0.001_DP)
    call check(error, is_close(z, pi, atol=1d-6), &
               'expected pi, got: '//to_string(z))
    if (allocated(error)) return
  end subroutine test_zero_br

  subroutine test_get_perturbation(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP) :: x, eps
    real(DP) :: v1, v2

    ! test perturbation calculation for pos and negative x
    v1 = get_perturbation(1.D0)
    v2 = -get_perturbation(-1.D0)
    call check(error, &
               is_close(v1, v2, atol=1d-12), &
               'expected '//to_string(v1)//' got: '//to_string(v2))
    if (allocated(error)) return

    ! test derivative calculation for sin(x) where x=1
    x = 1.d0
    eps = get_perturbation(x)
    v1 = (sin(x + eps) - sin(x)) / eps
    v2 = cos(x)
    call check(error, &
               is_close(v1, v2, atol=1d-5), &
               'expected '//to_string(v1)//' got: '//to_string(v2))
    if (allocated(error)) return

    ! test derivative calculation for sin(x) where x=0
    x = 0.d0
    eps = get_perturbation(x)
    v1 = (sin(x + eps) - sin(x)) / eps
    v2 = cos(x)
    call check(error, &
               is_close(v1, v2, atol=1d-5), &
               'expected '//to_string(v1)//' got: '//to_string(v2))
    if (allocated(error)) return

    ! test derivative calculation for x ** 2
    x = 1.d6
    eps = get_perturbation(x)
    v1 = ((x + eps)**2 - x**2) / eps
    v2 = 2 * x
    call check(error, &
               is_close(v1, v2, atol=1d-1), &
               'expected '//to_string(v1)//' got: '//to_string(v2))
    if (allocated(error)) return

  end subroutine test_get_perturbation

  subroutine test_arange(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), allocatable :: a(:)
    integer(I4B) :: i

    a = arange(DZERO, 10.0_DP, DONE)
    call check(error, size(a) == 10, "wrong size: "//to_string(size(a)))
    if (allocated(error)) return
    do i = 1, 10
      call check(error, is_close(a(i), real(i - 1, DP)))
      if (allocated(error)) return
    end do

    deallocate (a)

    a = arange(DZERO, DONE, 0.1_DP)
    call check(error, size(a) == 10, "wrong size: "//to_string(size(a)))
    if (allocated(error)) return
    do i = 1, 10
      call check(error, is_close(a(i), real(i - 1, DP) / 10.0_DP))
      if (allocated(error)) return
    end do

  end subroutine test_arange

  subroutine test_linspace(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), allocatable :: a(:)
    integer(I4B) :: i

    a = linspace(DONE, 10.0_DP, 10)
    call check(error, size(a) == 10)
    if (allocated(error)) return
    do i = 1, 10
      call check(error, is_close(a(i), real(i, DP)))
      if (allocated(error)) return
    end do

  end subroutine test_linspace

end module TestMathUtil
