module NumericalExchangeModule

  use KindModule, only: DP, I4B
  use BaseModelModule, only: BaseModelType
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use NumericalModelModule, only: NumericalModelType
  use ListModule, only: ListType
  use MatrixBaseModule

  implicit none

  private
  public :: NumericalExchangeType, &
            AddNumericalExchangeToList, GetNumericalExchangeFromList

  type, extends(BaseExchangeType) :: NumericalExchangeType
  contains

    procedure :: exg_df
    procedure :: exg_ac
    procedure :: exg_mc
    procedure :: exg_ar
    !procedure :: exg_rp (not needed yet; base exg_rp does nothing)
    procedure :: exg_cf
    procedure :: exg_fc
    procedure :: exg_cc
    procedure :: exg_cq
    procedure :: exg_bd
    procedure :: exg_ot
    procedure :: exg_da
    procedure :: get_iasym
  end type NumericalExchangeType

contains

  !> @brief Define the exchange
  subroutine exg_df(this)
    class(NumericalExchangeType) :: this
  end subroutine exg_df

  !> @brief If an implicit exchange then add connections to sparse
  subroutine exg_ac(this, sparse)
    use SparseModule, only: sparsematrix
    class(NumericalExchangeType) :: this
    type(sparsematrix), intent(inout) :: sparse
  end subroutine exg_ac

  !> @brief Map the connections in the global matrix
  subroutine exg_mc(this, matrix_sln)
    class(NumericalExchangeType) :: this
    class(MatrixBaseType), pointer :: matrix_sln
  end subroutine exg_mc

  !> @brief Allocate and read
  subroutine exg_ar(this)
    class(NumericalExchangeType) :: this
  end subroutine exg_ar

  !> @brief Calculate conductance, and for explicit exchanges, set the
  !! conductance in the boundary package
  subroutine exg_cf(this, kiter)
    class(NumericalExchangeType) :: this
    integer(I4B), intent(in) :: kiter
  end subroutine exg_cf

  !> @brief Fill the matrix
  subroutine exg_fc(this, kiter, matrix_sln, rhs_sln, inwtflag)
    class(NumericalExchangeType) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs_sln
    integer(I4B), optional, intent(in) :: inwtflag
  end subroutine exg_fc

  !> @brief Additional convergence check
  subroutine exg_cc(this, icnvg)
    class(NumericalExchangeType) :: this
    integer(I4B), intent(inout) :: icnvg
  end subroutine exg_cc

  !> @brief Calculate flow
  subroutine exg_cq(this, icnvg, isuppress_output, isolnid)
    class(NumericalExchangeType) :: this
    integer(I4B), intent(inout) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    integer(I4B), intent(in) :: isolnid
  end subroutine exg_cq

  !> @brief Exchange budget
  subroutine exg_bd(this, icnvg, isuppress_output, isolnid)
    class(NumericalExchangeType) :: this
    integer(I4B), intent(inout) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    integer(I4B), intent(in) :: isolnid
  end subroutine exg_bd

  !> @brief Output
  subroutine exg_ot(this)
    class(NumericalExchangeType) :: this
  end subroutine exg_ot

  !> @brief Deallocate memory
  subroutine exg_da(this)
    class(NumericalExchangeType) :: this
  end subroutine exg_da

  function get_iasym(this) result(iasym)
    class(NumericalExchangeType) :: this
    integer(I4B) :: iasym
    iasym = 0
  end function get_iasym

  function CastAsNumericalExchangeClass(obj) result(res)
    implicit none
    class(*), pointer, intent(inout) :: obj
    class(NumericalExchangeType), pointer :: res

    res => null()
    if (.not. associated(obj)) return

    select type (obj)
    class is (NumericalExchangeType)
      res => obj
    end select
  end function CastAsNumericalExchangeClass

  !> @brief Add numerical exchange to a list
  subroutine AddNumericalExchangeToList(list, exchange)
    implicit none
    type(ListType), intent(inout) :: list
    class(NumericalExchangeType), pointer, intent(in) :: exchange
    class(*), pointer :: obj

    obj => exchange
    call list%Add(obj)
  end subroutine AddNumericalExchangeToList

  !> @brief Retrieve a specific numerical exchange from a list
  function GetNumericalExchangeFromList(list, idx) result(res)
    implicit none
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: idx
    class(NumericalExchangeType), pointer :: res
    class(*), pointer :: obj

    obj => list%GetItem(idx)
    res => CastAsNumericalExchangeClass(obj)
  end function GetNumericalExchangeFromList

end module NumericalExchangeModule
