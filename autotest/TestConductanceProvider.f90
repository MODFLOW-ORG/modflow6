module TestConductanceProvider
  use KindModule, only: I4B, DP
  use testdrive, only: check, error_type, new_unittest, unittest_type
  use ConstantsModule, only: DZERO, DONE, DTWO, DHALF
  use MathUtilModule, only: is_close
  use HGeoUtilModule, only: hyeff
  use ConductanceProviderModule, only: ConductanceProviderType
  implicit none
  private
  public :: collect_conductanceprovider

  !> Minimal provider used to exercise the abstract contract. eff_hy encodes
  !> its inputs into the result so the test can confirm they pass through.
  type, extends(ConductanceProviderType) :: MockProvider
  contains
    procedure :: eff_hy => mock_eff_hy
  end type MockProvider

contains

  subroutine collect_conductanceprovider(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("provider_dispatch", test_provider_dispatch), &
                new_unittest("hyeff_axis_aligned", test_hyeff_axis_aligned), &
                new_unittest("hyeff_diagonal", test_hyeff_diagonal) &
                ]
  end subroutine collect_conductanceprovider

  function mock_eff_hy(this, n, ihc, vg) result(hy)
    class(MockProvider), intent(in) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: ihc
    real(DP), dimension(3), intent(in) :: vg
    real(DP) :: hy
    hy = real(n, DP) + real(ihc, DP) + vg(1)
  end function mock_eff_hy

  !> The abstract type is callable through a base-class pointer and forwards
  !> n, ihc, and vg to the concrete implementation.
  subroutine test_provider_dispatch(error)
    type(error_type), allocatable, intent(out) :: error
    type(MockProvider), target :: mock
    class(ConductanceProviderType), pointer :: provider
    real(DP) :: hy

    provider => mock
    hy = provider%eff_hy(7, 1, [DHALF, DZERO, DZERO])
    call check(error, is_close(hy, 7.0_DP + 1.0_DP + 0.5_DP))
  end subroutine test_provider_dispatch

  !> Along a principal axis (no rotation) the effective K is that axis' K.
  subroutine test_hyeff_axis_aligned(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), parameter :: k11 = 10.0_DP, k22 = 2.0_DP, k33 = 5.0_DP
    real(DP) :: hy

    hy = hyeff(k11, k22, k33, DZERO, DZERO, DZERO, DONE, DZERO, DZERO, 1)
    call check(error, is_close(hy, k11))
    if (allocated(error)) return
    hy = hyeff(k11, k22, k33, DZERO, DZERO, DZERO, DZERO, DONE, DZERO, 1)
    call check(error, is_close(hy, k22))
    if (allocated(error)) return
    hy = hyeff(k11, k22, k33, DZERO, DZERO, DZERO, DZERO, DZERO, DONE, 1)
    call check(error, is_close(hy, k33))
  end subroutine test_hyeff_axis_aligned

  !> A 45-degree in-plane direction gives the arithmetic mean of k11 and k22.
  subroutine test_hyeff_diagonal(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), parameter :: k11 = 10.0_DP, k22 = 2.0_DP, k33 = 5.0_DP
    real(DP) :: c, hy

    c = sqrt(DHALF)
    hy = hyeff(k11, k22, k33, DZERO, DZERO, DZERO, c, c, DZERO, 1)
    call check(error, is_close(hy, DHALF * (k11 + k22)))
  end subroutine test_hyeff_diagonal

end module TestConductanceProvider
