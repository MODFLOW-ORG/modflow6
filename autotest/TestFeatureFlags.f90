module TestFeatureFlags
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use FeatureFlagsModule, only: developmode
  use ConstantsModule, only: LINELENGTH
  use VersionModule, only: IDEVELOPMODE

  implicit none
  private
  public :: collect_feature_flags

contains

  subroutine collect_feature_flags(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                ! expect failure if in release mode, otherwise pass
                new_unittest("developmode", test_developmode, &
                             should_fail=(IDEVELOPMODE == 0)) &
                ]
  end subroutine collect_feature_flags

  subroutine test_developmode(error)
    type(error_type), allocatable, intent(out) :: error
    character(len=LINELENGTH) :: errmsg
    call developmode(errmsg)
  end subroutine test_developmode

end module TestFeatureFlags
