module UzrDevSoilModule
  use UzrSoilModelModule, only: SoilModelType
  use KindModule, only: I4B, DP
  use ConstantsModule, only: LINELENGTH, LENVARNAME, DONE, DZERO, DPREC
  use SimModule, only: store_error
  implicit none
  private

  public :: UzrDevSoilType

  type, extends(SoilModelType) :: UzrDevSoilType
    character(len=LINELENGTH) :: dev_soil_param
  contains
    procedure :: load => load_custom
    procedure :: saturation => saturation_custom
    procedure :: capacity => capacity_custom
    procedure :: krelative => krelative_custom
    procedure :: destroy => destroy_custom
  end type UzrDevSoilType

contains

  subroutine load_custom(this, porosity, sat_res)
    use MemoryManagerExtModule, only: mem_set_value
    use GwfUzrInputModule, only: GwfUzrParamFoundType
    class(UzrDevSoilType), intent(inout) :: this !< this instance
    real(DP), dimension(:), pointer, contiguous :: porosity !< porosity
    real(DP), dimension(:), pointer, contiguous :: sat_res !< residual saturation
    ! local
    type(GwfUzrParamFoundType) :: found

    this%porosity => porosity
    this%sat_res => sat_res

    this%dev_soil_param = ""
    call mem_set_value(this%dev_soil_param, 'DEV_SOIL_PARAM', &
                       this%input_mem_path, found%dev_soil_param)
    if (.not. found%dev_soil_param) then
      call store_error('Soil parameter set not specified in input', .true.)
    end if

  end subroutine load_custom

  function saturation_custom(this, psi, i) result(s)
    class(UzrDevSoilType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: s !< the saturation
    ! local
    real(DP) :: s_eff

    if (psi > -DPREC) then
      s = DONE
      return
    end if

    select case (this%dev_soil_param)
    case ("HVK_YOLO_1977")
      s_eff = 739.0 / (739.0 + (log(-psi))**4.0)
      s = (DONE - this%sat_res(i)) * s_eff + this%sat_res(i)
    case ("HVK_SAND_1977")
      s_eff = 1.611e+06 / (1.611e+06 + (-psi)**3.96)
      s = (DONE - this%sat_res(i)) * s_eff + this%sat_res(i)
    case default
      call store_error('Invalid soil model: '//trim(this%dev_soil_param), .true.)
    end select

  end function saturation_custom

  function capacity_custom(this, psi, i) result(Cm)
    class(UzrDevSoilType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: Cm !< the saturation
    ! local
    real(DP) :: num, denom

    if (psi > -DPREC) then
      Cm = DZERO
      return
    end if

    select case (this%dev_soil_param)
    case ("HVK_YOLO_1977")
      num = (739 * 4) * (this%porosity(i) - this%porosity(i) * &
                         this%sat_res(i)) * (log(-psi))**3.0
      denom = psi * (739 + (log(-psi))**4)**2
      Cm = num / denom
    case ("HVK_SAND_1977")
      num = (3.96 * 1.611e+06) * &
            (this%porosity(i) - this%porosity(i) * this%sat_res(i)) * (-psi)**2.96
      denom = psi * (1.611e+06 + (-psi)**3.96)**2
      Cm = num / denom
    case default
      call store_error('Invalid soil model: '//trim(this%dev_soil_param), .true.)
    end select

  end function capacity_custom

  function krelative_custom(this, psi, i) result(kr)
    class(UzrDevSoilType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: kr !< the saturation

    if (psi > DZERO) then
      kr = DONE
      return
    end if

    select case (this%dev_soil_param)
    case ("HVK_YOLO_1977")
      kr = 124.6 / (124.6 + (-psi)**1.77)
    case ("HVK_SAND_1977")
      kr = 1.175e+06 / (1.175e+06 + (-psi)**4.74)
    case default
      call store_error('Invalid soil model: '//trim(this%dev_soil_param), .true.)
    end select

  end function krelative_custom

!> @brief clean up
!<
  subroutine destroy_custom(this)
    class(UzrDevSoilType), intent(inout) :: this !< this instance

    this%porosity => null()
    this%sat_res => null()

    call this%destroy_base()

  end subroutine destroy_custom

end module UzrDevSoilModule
