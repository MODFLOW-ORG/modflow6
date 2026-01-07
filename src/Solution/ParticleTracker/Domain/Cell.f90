module CellModule

  use CellDefnModule, only: CellDefnType
  use DomainModule, only: DomainType
  use KindModule, only: I4B
  implicit none
  private
  public :: CellType

  !> @brief Base type for grid cells of a concrete type. Contains
  !! a cell-definition which is information shared by cell types.
  type, abstract, extends(DomainType) :: CellType
    type(CellDefnType), pointer :: defn => null() ! cell defn
  end type CellType

end module CellModule
