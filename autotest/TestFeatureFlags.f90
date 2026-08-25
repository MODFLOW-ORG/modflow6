module TestFeatureFlags
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use FeatureFlagsModule, only: is_release_mode
  use VersionModule, only: IDEVELOPMODE

  implicit none
  private
  public :: collect_feature_flags

contains

  subroutine collect_feature_flags(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("is_release_mode", test_is_release_mode) &
                ]
  end subroutine collect_feature_flags

  subroutine test_is_release_mode(error)
    type(error_type), allocatable, intent(out) :: error
    call check(error, is_release_mode() .eqv. (IDEVELOPMODE == 0))
  end subroutine test_is_release_mode

end module TestFeatureFlags
