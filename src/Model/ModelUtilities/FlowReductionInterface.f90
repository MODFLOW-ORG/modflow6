module FlowReductionInterfaceModule
  use KindModule, only: DP, I4B

  implicit none
  private
  public :: FlowReductionType

  !> @brief Interface for supplying a flow-reduction interval.
  !!
  !! A stress package that reduces a withdrawal, such as the Well Package
  !! with AUTO_FLOW_REDUCE, tapers the rate over an interval of cell
  !! elevations. By default that interval spans the cell, so the rate falls
  !! away as the cell empties.
  !!
  !! Where only part of a cell holds the fluid being withdrawn, the taper
  !! should be referenced to that part instead. A package that knows how the
  !! cell is divided implements this interface and narrows the interval; the
  !! stress package reads it without knowing which package provided it, or
  !! whether one did at all.
  !<
  type, abstract :: FlowReductionType
  contains
    procedure(interval_if), deferred :: interval !< Narrow the interval for node n
  end type FlowReductionType

  abstract interface
    !> @brief Narrow the flow-reduction interval for a node
    !!
    !! The bounds arrive holding the interval the caller would otherwise use,
    !! normally the cell bottom and top. An implementation raises bt, lowers
    !! tp, or both, to the part of the cell its fluid occupies. Leaving them
    !! unchanged means the whole cell is available.
    !<
    subroutine interval_if(this, n, bt, tp)
      ! -- import
      import :: FlowReductionType, DP, I4B
      ! -- dummy
      class(FlowReductionType), intent(in) :: this !< flow reduction provider
      integer(I4B), intent(in) :: n !< reduced node number
      real(DP), intent(inout) :: bt !< bottom of the interval
      real(DP), intent(inout) :: tp !< top of the interval
    end subroutine interval_if
  end interface

end module FlowReductionInterfaceModule
