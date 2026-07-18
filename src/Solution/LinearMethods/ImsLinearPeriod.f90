!> @brief Stress-period-varying IMS linear settings
!!
!! Optional utility, enabled from the IMS OPTIONS block, that overrides selected
!! linear settings (inner closure tolerances, inner maximum, and the
!! preconditioner) on a per-stress-period basis. The period file uses keyword
!! PERIOD blocks drawn from the LINEAR-block vocabulary with sticky-until-changed
!! semantics: a setting persists to later periods until a subsequent PERIOD block
!! overrides it, and periods without a block inherit the running settings.
!!
!! The implementation lives in the ims-linear-period submodule so the core IMS
!! source files stay small.
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
    !! Reads the PERIOD block for stress period kper (if present) and mutates the
    !! running settings in place (sticky-until-changed). Sets changed to .true.
    !! when any setting was modified so the caller can reconfigure the solver.
    !! The allow_tol/allow_precond flags gate which setting categories the active
    !! solver mode honors; a keyword whose category is disabled is rejected with
    !! an error rather than silently ignored.
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
