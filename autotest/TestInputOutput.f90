module TestInputOutput
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use KindModule, only: I4B, DP
  use ConstantsModule, only: LINELENGTH
  use InputOutputModule, only: urword
  implicit none
  private
  public :: collect_inputoutput

contains

  subroutine collect_inputoutput(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("urword", test_urword) &
                ]
  end subroutine collect_inputoutput

  subroutine test_urword(error)
    type(error_type), allocatable, intent(out) :: error
    character(len=LINELENGTH) :: line
    integer(I4B) :: icol
    integer(I4B) :: istart
    integer(I4B) :: istop
    integer(I4B) :: ncode
    integer(I4B) :: n
    integer(I4B) :: iout
    integer(I4B) :: in
    real(DP) :: r

    ! parse integer
    line = "77"
    icol = 1
    ncode = 2
    call URWORD(line, icol, istart, istop, ncode, n, r, iout, in)
    call check(error, n == 77)
    if (allocated(error)) return

    ! parse float
    line = "1.0"
    icol = 1
    ncode = 3
    call URWORD(line, icol, istart, istop, ncode, n, r, iout, in)
    call check(error, r == 1.D0)
    if (allocated(error)) return

    ! parse string
    line = "mymodel"
    icol = 1
    ncode = 0
    call URWORD(line, icol, istart, istop, ncode, n, r, iout, in)
    call check(error, line(istart:istop) == "mymodel")
    if (allocated(error)) return

    ! parse string and convert to caps
    line = "mymodel"
    icol = 1
    ncode = 1
    call URWORD(line, icol, istart, istop, ncode, n, r, iout, in)
    call check(error, line(istart:istop) == "MYMODEL")
    if (allocated(error)) return

    ! parse float with more than 30 characters
    line = "1000000000000000000000000000000000000000"
    icol = 1
    ncode = 3
    call URWORD(line, icol, istart, istop, ncode, n, r, iout, in)
    call check(error, r == 1D39)
    if (allocated(error)) return

    ! parse empty line into a zero float value
    line = "       "
    icol = 1
    ncode = 3
    r = -999.D0
    call URWORD(line, icol, istart, istop, ncode, n, r, iout, in)
    call check(error, r == 0.D0)
    if (allocated(error)) return

    ! parse empty line into a zero integer value
    line = "       "
    icol = 1
    ncode = 2
    n = -999
    call URWORD(line, icol, istart, istop, ncode, n, r, iout, in)
    call check(error, n == 0)
    if (allocated(error)) return

  end subroutine test_urword

end module TestInputOutput
