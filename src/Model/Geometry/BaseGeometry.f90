module BaseGeometryModule

  use KindModule, only: DP, I4B

  implicit none
  private
  public BaseGeometryType

  integer(I4B), parameter :: GEONAMELEN = 20

  type :: BaseGeometryType
    character(len=20) :: geo_type = 'UNDEFINED'
    integer(I4B) :: id = 0
    character(len=GEONAMELEN) :: name = ''

  contains

    procedure :: area_sat
    procedure :: perimeter_sat
    procedure :: area_wet
    procedure :: perimeter_wet
    procedure :: set_attribute
    procedure :: print_attributes
  end type BaseGeometryType

contains

  function area_sat(this)
    ! -- return
    real(DP) :: area_sat
    ! -- dummy
    class(BaseGeometryType) :: this
    !
    area_sat = 0.d0
  end function area_sat

  function perimeter_sat(this)
    ! -- return
    real(DP) :: perimeter_sat
    ! -- dummy
    class(BaseGeometryType) :: this
    !
    perimeter_sat = 0.d0
  end function perimeter_sat

  function area_wet(this, depth)
    ! -- return
    real(DP) :: area_wet
    ! -- dummy
    class(BaseGeometryType) :: this
    real(DP), intent(in) :: depth
    !
    area_wet = 0.d0
  end function area_wet

  function perimeter_wet(this, depth)
    ! -- return
    real(DP) :: perimeter_wet
    ! -- dummy
    class(BaseGeometryType) :: this
    real(DP), intent(in) :: depth
    !
    perimeter_wet = 0.d0
  end function perimeter_wet

  subroutine set_attribute(this, line)
    ! -- dummy
    class(BaseGeometryType) :: this
    character(len=*), intent(inout) :: line
  end subroutine set_attribute

  !> @brief Print the attributes for this object
  !<
  subroutine print_attributes(this, iout)
    ! -- dummy
    class(BaseGeometryType) :: this
    ! -- local
    integer(I4B), intent(in) :: iout
    ! -- formats
    character(len=*), parameter :: fmtid = "(4x,a,i0)"
    character(len=*), parameter :: fmtnm = "(4x,a,a)"
    !
    write (iout, fmtid) 'ID = ', this%id
    write (iout, fmtnm) 'NAME = ', trim(adjustl(this%name))
    write (iout, fmtnm) 'GEOMETRY TYPE = ', trim(adjustl(this%geo_type))
  end subroutine print_attributes

end module BaseGeometryModule
