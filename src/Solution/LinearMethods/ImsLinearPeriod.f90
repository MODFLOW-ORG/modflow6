!> @brief Stress-period-varying IMS linear settings
!!
!! Optional utility, enabled from the IMS OPTIONS block, that overrides selected
!! linear settings on a per-stress-period basis using keyword PERIOD blocks with
!! sticky-until-changed semantics. The implementation is in the ims-linear-period
!! submodule.
!<
module ImsLinearPeriodModule
  use KindModule, only: DP, I4B, LGP
  use BlockParserModule, only: BlockParserType
  use ImsLinearSettingsModule, only: ImsLinearSettingsType
  implicit none
  private
  public :: ImsLinearPeriodType

  !> @brief Reader/state for period-varying IMS linear settings
  !<
  type :: ImsLinearPeriodType
    logical(LGP) :: active = .false. !< period data supplied and in use
    integer(I4B) :: inunit = 0 !< period data file unit (0 = disabled)
    integer(I4B) :: iout = 0 !< listing file unit
    integer(I4B) :: ionper = 0 !< stress period of the next PERIOD block to read
    integer(I4B) :: lastonper = 0 !< previous ionper (increasing-period check)
    type(BlockParserType) :: parser !< parser bound to the period data file
  contains
    procedure :: init => imslinperiod_init
    procedure :: read_period => imslinperiod_read_period
    procedure :: da => imslinperiod_da
  end type ImsLinearPeriodType

  interface
    !> @brief Initialize the period-settings utility from an open file unit
    !<
    module subroutine imslinperiod_init(this, inunit, iout)
      class(ImsLinearPeriodType), intent(inout) :: this !< period-settings utility
      integer(I4B), intent(in) :: inunit !< period file unit (0 = disabled)
      integer(I4B), intent(in) :: iout !< listing file unit
    end subroutine imslinperiod_init

    !> @brief Apply this stress period's linear-setting overrides
    !!
    !! Read the PERIOD block for stress period kper, if present, and mutate the
    !! running settings in place (sticky-until-changed). allow_tol/allow_precond
    !! gate which categories are honored; a disabled keyword is an error.
    !<
    module subroutine imslinperiod_read_period(this, settings, kper, &
                                               allow_tol, allow_precond, changed)
      class(ImsLinearPeriodType), intent(inout) :: this !< period-settings utility
      type(ImsLinearSettingsType), intent(inout) :: settings !< running linear settings
      integer(I4B), intent(in) :: kper !< current stress period
      logical(LGP), intent(in) :: allow_tol !< inner tolerances/maximum may vary
      logical(LGP), intent(in) :: allow_precond !< preconditioner settings may vary
      logical(LGP), intent(out) :: changed !< true if any setting was changed
    end subroutine imslinperiod_read_period

    !> @brief Deallocate the period-settings utility
    !<
    module subroutine imslinperiod_da(this)
      class(ImsLinearPeriodType), intent(inout) :: this !< period-settings utility
    end subroutine imslinperiod_da
  end interface

end module ImsLinearPeriodModule
