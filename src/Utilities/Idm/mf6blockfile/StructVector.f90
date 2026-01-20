!> @brief This module contains the StructVectorModule
!!
!! This module contains a generic type for storing
!! different types of vectors.
!!
!<
module StructVectorModule

  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: DNODATA, LENMEMPATH, LENVARNAME, LINELENGTH, &
                             LENTIMESERIESNAME
  use ListModule, only: ListType
  use InputDefinitionModule, only: InputParamDefinitionType
  use CharacterStringModule, only: CharacterStringType
  use STLVecIntModule, only: STLVecInt

  implicit none
  private
  public :: StructVectorType, TSStringLocType

  !> @brief derived type which describes time series string field
  !<
  type :: TSStringLocType
    integer(I4B) :: structarray_col !< global SA column index
    integer(I4B) :: col !< SV column (1 if 1d array)
    integer(I4B) :: row !< SV row
    character(LINELENGTH) :: token !< TS string token
  contains
  end type TSStringLocType

  !> @brief derived type for generic vector
  !!
  !! This derived type is used in the StructArrayType to
  !! store any type of vector.
  !!
  !<
  type StructVectorType
    type(InputParamDefinitionType), pointer :: idt !< input definition
    ! SA vector attributes
    integer(I4B) :: memtype = 0 !< SA memtype
    integer(I4B) :: icol = 0 !< SA column
    integer(I4B) :: size = 0 !< size of array
    ! Data pointers
    integer(I4B), dimension(:), pointer, contiguous :: int1d => null()
    integer(I4B), dimension(:, :), pointer, contiguous :: int2d => null()
    real(DP), dimension(:), pointer, contiguous :: dbl1d => null()
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      charstr1d => null()
    type(STLVecInt), pointer :: intvector => null()
    ! Shape data pointers
    integer(I4B), pointer :: intshape => null()
    integer(I4B), dimension(:), pointer, contiguous :: intvector_shape => null()
    ! TimeSeries strings
    type(ListType) :: ts_strlocs
  contains
    procedure :: clear => sv_clear
    procedure :: read_token => sv_read_token
    procedure :: add_ts_strloc => sv_add_ts_strloc
    procedure :: get_ts_strloc => sv_get_ts_strloc
  end type StructVectorType

contains

  function sv_read_token(this, token, structarray_col, col, row) result(val)
    class(StructVectorType) :: this
    character(len=*), intent(in) :: token
    integer(I4B), intent(in) :: structarray_col
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    real(DP) :: val
    integer(I4B) :: istat
    real(DP) :: r
    ! initialize
    val = DNODATA
    read (token, *, iostat=istat) r
    if (istat == 0) then
      val = r
    else
      call this%add_ts_strloc(token, structarray_col, col, row)
    end if
  end function sv_read_token

  subroutine sv_add_ts_strloc(this, token, structarray_col, col, row)
    class(StructVectorType) :: this
    character(len=*), intent(in) :: token
    integer(I4B), intent(in) :: structarray_col
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    class(TSStringLocType), pointer :: str_field
    class(*), pointer :: obj
    allocate (str_field)
    str_field%structarray_col = structarray_col
    str_field%col = col
    str_field%row = row
    str_field%token = token
    obj => str_field
    call this%ts_strlocs%Add(obj)
  end subroutine sv_add_ts_strloc

  function sv_get_ts_strloc(this, idx) result(res)
    class(StructVectorType) :: this
    integer(I4B), intent(in) :: idx !< package number
    class(TSStringLocType), pointer :: res
    class(*), pointer :: obj
    ! initialize res
    res => null()
    ! get the package from the list
    obj => this%ts_strlocs%GetItem(idx)
    if (associated(obj)) then
      select type (obj)
      class is (TSStringLocType)
        res => obj
      end select
    end if
  end function sv_get_ts_strloc

  !> @brief
  !<
  subroutine sv_clear(this)
    class(StructVectorType) :: this
    class(TSStringLocType), pointer :: ts_strloc
    integer(I4B) :: n
    do n = 1, this%ts_strlocs%Count()
      ts_strloc => this%get_ts_strloc(n)
      deallocate (ts_strloc)
      nullify (ts_strloc)
    end do
    call this%ts_strlocs%Clear()
  end subroutine sv_clear

end module StructVectorModule
