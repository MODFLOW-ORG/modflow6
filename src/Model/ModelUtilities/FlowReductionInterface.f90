module FlowReductionInterfaceModule
  use KindModule, only: DP, I4B

  implicit none
  private
  public :: FlowReductionType

  !> @brief Interface for supplying a flow-reduction factor.
  !!
  !! A stress package that reduces a withdrawal, such as the Well Package
  !! with AUTO_FLOW_REDUCE, multiplies the requested rate by a factor that
  !! falls from one to zero as the fluid runs out. On its own the package
  !! judges that from the cell: the factor tapers over a length measured up
  !! from the cell bottom, and only a convertible cell can run dry.
  !!
  !! A flow formulation may know better how much of the fluid being
  !! withdrawn the cell holds: a sharp interface knows the elevation range
  !! the fluid occupies, and an unsaturated-flow formulation knows its
  !! saturation. Such a formulation implements this interface and returns the
  !! factor itself, with its derivative, so the stress package applies it
  !! without knowing which formulation provided it, or whether one did at all.
  !<
  type, abstract :: FlowReductionType
  contains
    procedure(reduce_if), deferred :: reduce !< Factor for node n
  end type FlowReductionType

  abstract interface
    !> @brief Flow-reduction factor for a node, and its derivative
    !!
    !! The taper length is what the stress package was asked to taper over,
    !! resolved to a length: the per-boundary auxiliary value or the package
    !! value, a fraction of the cell thickness already converted. An
    !! implementation decides how much of its fluid the cell holds, tapers
    !! the factor over that length from where the fluid runs out, and may
    !! ignore the length if its fluid never runs out at a definite elevation.
    !! The factor is one where nothing limits the withdrawal and zero where
    !! there is nothing to withdraw; the derivative is with respect to the
    !! head of the cell in this model.
    !<
    subroutine reduce_if(this, n, h, taper, f, dfdh)
      ! -- import
      import :: FlowReductionType, DP, I4B
      ! -- dummy
      class(FlowReductionType), intent(in) :: this !< flow reduction provider
      integer(I4B), intent(in) :: n !< reduced node number
      real(DP), intent(in) :: h !< head in the cell
      real(DP), intent(in) :: taper !< length to taper the factor over
      real(DP), intent(out) :: f !< factor to apply to the requested rate
      real(DP), intent(out) :: dfdh !< derivative of the factor with respect to head
    end subroutine reduce_if
  end interface

end module FlowReductionInterfaceModule
