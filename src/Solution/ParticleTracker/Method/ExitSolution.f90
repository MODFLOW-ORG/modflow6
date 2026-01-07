module ExitSolutionModule

  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: DZERO

  implicit none

  !> @brief Exit status codes
  enum, bind(C)
    enumerator :: OK_EXIT = 0 !< exit found using velocity interpolation
    ! below only used for linear solution
    enumerator :: OK_EXIT_CONSTANT = 1 !< exit found, constant velocity
    enumerator :: NO_EXIT_STATIONARY = 2 !< no exit, zero velocity
    enumerator :: NO_EXIT_NO_OUTFLOW = 3 !< no exit, no outflow
  end enum

  !> @brief Base type for exit solutions
  type :: ExitSolutionType
    integer(I4B) :: status = -1 !< domain exit status code
    integer(I4B) :: iboundary = 0 !< boundary number
    real(DP) :: dt = 1.0d+20 !< time to exit
  end type ExitSolutionType

  !> @brief Linear velocity interpolation exit solution
  type, extends(ExitSolutionType) :: LinearExitSolutionType
    real(DP) :: v = DZERO !< particle velocity
    real(DP) :: dvdx = DZERO !< velocity gradient
  end type LinearExitSolutionType

end module ExitSolutionModule
