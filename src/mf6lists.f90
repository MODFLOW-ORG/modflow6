module ListsModule
  ! -- Contains lists of base types BaseModelType,
  !    BaseSolutionType, SolutionGroupType, and
  !    BaseExchangeType for use by any MF6 module.

  use KindModule, only: DP, I4B
  use ListModule, only: ListType

  implicit none
  private
  public :: basemodellist, basesolutionlist, solutiongrouplist, &
            baseexchangelist, baseconnectionlist
  public :: lists_da

  ! -- list of all models in simulation
  type(ListType) :: basemodellist

  ! -- list of all solutions in simulation
  type(ListType) :: basesolutionlist

  ! -- list of all solutiongroups in simulation
  type(ListType) :: solutiongrouplist

  ! -- list of all exchanges in simulation
  type(ListType) :: baseexchangelist

  ! -- list of all connections in simulation
  type(ListType) :: baseconnectionlist

contains

  !< @brief deallocate
  subroutine lists_da()
    !
    call basemodellist%Clear()
    call basesolutionlist%Clear()
    call solutiongrouplist%Clear()
    call baseexchangelist%Clear()
    call baseconnectionlist%Clear()
  end subroutine lists_da

end module ListsModule
