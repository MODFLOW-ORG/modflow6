module GwfNpfOptionsModule
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DONE

  implicit none
  private

  !> Data structure and helper methods for passing NPF options
  !! into npf_df, as an alternative to reading those from file
  !<
  type, public :: GwfNpfOptionsType
    integer(I4B) :: ithickstrt !< same as npf variable
    integer(I4B) :: ihighcellsat !< same as npf variable
    integer(I4B) :: iperched !< same as npf variable
    integer(I4B) :: ivarcv !< same as npf variable
    integer(I4B) :: idewatcv !< same as npf variable
    integer(I4B) :: irewet !< same as npf variable
    real(DP) :: wetfct !< same as npf variable
    integer(I4B) :: iwetit !< same as npf variable
    integer(I4B) :: ihdwet !< same as npf variable
  contains
    procedure, pass(this) :: construct
    procedure, pass(this) :: destroy
  end type GwfNpfOptionsType

contains

  !> @brief construct the input options, set variables to their defaults
  !<
  subroutine construct(this)
    class(GwfNpfOptionsType), intent(inout) :: this !< the NPF options, as in the input OPTIONS block

    this%ithickstrt = 0
    this%ihighcellsat = 0
    this%iperched = 0
    this%ivarcv = 0
    this%idewatcv = 0
    this%irewet = 0
    this%wetfct = DONE
    this%iwetit = 1
    this%ihdwet = 0

  end subroutine construct

  !> @brief cleans up
  !<
  subroutine destroy(this)
    class(GwfNpfOptionsType), intent(inout) :: this !< the NPF options, as in the input OPTIONS block

    ! nothing to be done here for now...

  end subroutine destroy

end module GwfNpfOptionsModule
