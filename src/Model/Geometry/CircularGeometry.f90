module CircularGeometryModule
  use KindModule, only: DP, I4B
  use BaseGeometryModule, only: BaseGeometryType
  use ConstantsModule, only: DZERO

  implicit none

  private
  public :: CircularGeometryType

  type, extends(BaseGeometryType) :: CircularGeometryType
    real(DP) :: radius = DZERO

  contains

    procedure :: area_sat
    procedure :: perimeter_sat
    procedure :: area_wet
    procedure :: perimeter_wet
    procedure :: set_attribute
    procedure :: print_attributes
  end type CircularGeometryType

contains

  !> @brief Return area as if geometry is fully saturated
  !<
  function area_sat(this)
    ! -- modules
    use ConstantsModule, only: DTWO, DPI
    ! -- return
    real(DP) :: area_sat
    ! -- dummy
    class(CircularGeometryType) :: this
    !
    ! -- Calculate area
    area_sat = DPI * this%radius**DTWO
  end function area_sat

  !> @brief Return perimeter as if geometry is fully saturated
  !<
  function perimeter_sat(this)
    ! -- modules
    use ConstantsModule, only: DTWO, DPI
    ! -- return
    real(DP) :: perimeter_sat
    ! -- dummy
    class(CircularGeometryType) :: this
    !
    ! -- Calculate area
    perimeter_sat = DTWO * DPI * this%radius
  end function perimeter_sat

  !> @brief Return wetted area
  !<
  function area_wet(this, depth)
    ! -- modules
    use ConstantsModule, only: DTWO, DPI, DZERO
    ! -- return
    real(DP) :: area_wet
    ! -- dummy
    class(CircularGeometryType) :: this
    real(DP), intent(in) :: depth
    !
    ! -- Calculate area
    if (depth <= DZERO) then
      area_wet = DZERO
    elseif (depth <= this%radius) then
      area_wet = this%radius * this%radius * &
                 acos((this%radius - depth) / this%radius) - &
                 (this%radius - depth) * &
                 sqrt(this%radius * this%radius - (this%radius - depth)**DTWO)
    elseif (depth <= DTWO * this%radius) then
      area_wet = this%radius * this%radius * &
                 (DPI - acos((depth - this%radius) / this%radius)) - &
                 (this%radius - depth) * &
                 sqrt(this%radius * this%radius - (this%radius - depth)**DTWO)
    else
      area_wet = DPI * this%radius * this%radius
    end if
  end function area_wet

  !> @brief Return wetted perimeter
  !<
  function perimeter_wet(this, depth)
    ! -- modules
    use ConstantsModule, only: DTWO, DPI
    ! -- return
    real(DP) :: perimeter_wet
    ! -- dummy
    class(CircularGeometryType) :: this
    real(DP), intent(in) :: depth
    !
    ! -- Calculate area
    if (depth <= DZERO) then
      perimeter_wet = DZERO
    elseif (depth <= this%radius) then
      perimeter_wet = DTWO * this%radius * acos((this%radius - depth) / &
                                                this%radius)
    elseif (depth <= DTWO * this%radius) then
      perimeter_wet = DTWO * this%radius * (DPI - acos((depth - this%radius) / &
                                                       this%radius))
    else
      perimeter_wet = DTWO * DPI * this%radius
    end if
  end function perimeter_wet

  !> @brief Set a parameter for this circular object
  !<
  subroutine set_attribute(this, line)
    ! -- module
    use InputOutputModule, only: urword
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error, count_errors
    ! -- dummy
    class(CircularGeometryType) :: this
    character(len=LINELENGTH) :: errmsg
    character(len=*), intent(inout) :: line
    ! -- local
    integer(I4B) :: lloc, istart, istop, ival
    real(DP) :: rval
    !
    ! -- should change this and set id if uninitialized or store it
    lloc = 1
    call urword(line, lloc, istart, istop, 2, ival, rval, 0, 0)
    this%id = ival

    ! -- Parse the attribute
    call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
    select case (line(istart:istop))
    case ('NAME')
      call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)
      this%name = line(istart:istop)
    case ('RADIUS')
      call urword(line, lloc, istart, istop, 3, ival, rval, 0, 0)
      this%radius = rval
    case default
      write (errmsg, '(a,a)') &
        'Unknown circular geometry attribute: ', line(istart:istop)
      call store_error(errmsg, terminate=.TRUE.)
    end select
  end subroutine set_attribute

  !> @brief Print the attributes for this object
  !<
  subroutine print_attributes(this, iout)
    ! -- dummy
    class(CircularGeometryType) :: this
    ! -- local
    integer(I4B), intent(in) :: iout
    ! -- formats
    character(len=*), parameter :: fmtnm = "(4x,a,a)"
    character(len=*), parameter :: fmttd = "(4x,a,1(1PG15.6))"
    !
    ! -- call parent to print parent attributes
    call this%BaseGeometryType%print_attributes(iout)
    !
    ! -- Print specifics of this geometry type
    write (iout, fmttd) 'RADIUS = ', this%radius
    write (iout, fmttd) 'SATURATED AREA = ', this%area_sat()
    write (iout, fmttd) 'SATURATED WETTED PERIMETER = ', this%perimeter_sat()
  end subroutine print_attributes

end module CircularGeometryModule
