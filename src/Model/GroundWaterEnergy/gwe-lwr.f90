!> @brief This module contains methods for calculating the longwave radiation heat flux
!!
!! This module contains the methods used to calculate the longwave radiation
!! heat flux for surface-water boundaries, like streams and lakes.  In its
!! current form, this class acts like a package to a package, similar to the
!! TVK package (or utiliity) that can be invoked from the NPF package.  Once
!! this package is completed in its prototyped form, it will likely be moved
!! around.
!<

module LongwaveModule

  use ConstantsModule, only: LINELENGTH, LENMEMPATH, DZERO, LENVARNAME, &
                             DSTEFANBOLTZMANN, DCTOK, DFOUR
  use KindModule, only: I4B, DP
  use MemoryManagerModule, only: mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use TimeSeriesLinkModule, only: TimeSeriesLinkType
  use TimeSeriesManagerModule, only: TimeSeriesManagerType, tsmanager_cr
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use PbstBaseModule !, only: PbstBaseType, pbstbase_da

  implicit none

  public :: LwrType
  public :: lwr_cr

  character(len=16) :: text = '          LWR'

  type, extends(PbstBaseType) :: LwrType

    real(DP), pointer :: lwrefl => null() !< longwave reflectance of water surface
    real(DP), pointer :: emissw => null() !< emissivity of water
    real(DP), pointer :: emissr => null() !< emissivity of riparian canopy
    real(DP), dimension(:), pointer, contiguous :: shd => null() !< shade fraction
    real(DP), dimension(:), pointer, contiguous :: atmc => null() !< atmospheric composition adjustment
    real(DP), dimension(:), pointer, contiguous :: tatm => null() !< temperature of the atmosphere
    real(DP), dimension(:), pointer, contiguous :: rh => null() !< relative humidity
    real(DP), dimension(:), pointer, contiguous :: ea => null() !< ambient vapor pressure of the atmosphere
    real(DP), dimension(:), pointer, contiguous :: es => null() !< saturation vapor pressure at air temperature
    real(DP), dimension(:), pointer, contiguous :: ew => null() !< saturation vapor pressure at water temperature

  contains

    procedure :: da => lwr_da
    procedure :: pbst_ar => lwr_ar_set_pointers
    procedure, public :: lwr_cq
    procedure, private :: calc_lwr

  end type LwrType

contains

  !> @brief Create a new LwrType object
  !!
  !! Create a new longwave radiation flux (LwrType) object. Initially for use with
  !! the SFE package.
  !<
  subroutine lwr_cr(this, name_model, inunit, iout, ncv)
    ! -- dummy
    type(LwrType), pointer, intent(out) :: this
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    integer(I4B), target, intent(in) :: ncv
    !
    allocate (this)
    call this%init(name_model, 'LWR', 'LWR', inunit, iout, ncv)
    this%text = text
  end subroutine lwr_cr

  !> @brief Announce package and set pointers to variables
  !!
  !! Announce package version and set array and variable pointers from the ABC
  !! package for access by LWR.
  !<
  subroutine lwr_ar_set_pointers(this)
    ! -- dummy
    class(LwrType) :: this
    ! -- local
    character(len=LENMEMPATH) :: abcMemoryPath
    !
    ! -- print a message noting that the LWR utility is active
    write (this%iout, '(a)') &
      'LWR -- LONGWAVE RADIATION WILL BE INCLUDED IN THE ATMOSPHERIC '// &
      'BOUNDARY CONDITIONS FOR THE STREAMFLOW ENERGY TRANSPORT PACKAGE'
    !
    ! -- set pointers to variables hosted in the ABC package
    abcMemoryPath = create_mem_path(this%name_model, 'ABC')
    call mem_setptr(this%lwrefl, 'LWREFL', abcMemoryPath)
    call mem_setptr(this%shd, 'SHD', abcMemoryPath)
    call mem_setptr(this%atmc, 'ATMC', abcMemoryPath)
    call mem_setptr(this%tatm, 'TATM', abcMemoryPath)
    call mem_setptr(this%emissw, 'EMISSW', abcMemoryPath)
    call mem_setptr(this%emissr, 'EMISSR', abcMemoryPath)
    call mem_setptr(this%rh, 'RH', abcMemoryPath)
    call mem_setptr(this%rh, 'EA', abcMemoryPath)
    call mem_setptr(this%rh, 'ES', abcMemoryPath)
    call mem_setptr(this%rh, 'EW', abcMemoryPath)
    !
    ! -- set the riparian canopy emissivity constant
    this%emissr = 0.97_DP
    !
    ! -- create time series manager
    call tsmanager_cr(this%tsmanager, this%iout, &
                      removeTsLinksOnCompletion=.true., &
                      extendTsToEndOfSimulation=.true.)
  end subroutine lwr_ar_set_pointers

  !> @brief Calculate Longwave Radiation Heat Flux
  !!
  !! Calculate and return the longwave radiation heat flux for one reach
  !<
  subroutine lwr_cq(this, ifno, tstrm, lwrflx)
    ! -- dummy
    class(LwrType), intent(inout) :: this
    integer(I4B), intent(in) :: ifno !< stream reach integer id
    real(DP), intent(in) :: tstrm !< temperature of the stream reach
    real(DP), intent(inout) :: lwrflx !< calculated longwave radiation heat flux amount
    ! -- local
    real(DP) :: sat_vap_ta
    real(DP) :: amb_vap_atm
    real(DP) :: emissa
    real(DP) :: emisss
    real(DP) :: lwratm
    real(DP) :: lwrstrm
    !
    ! -- intermediate calculations
    !
    ! -- atmospheric emissivity (A.14)
    emissa = this%epsa(this%ea(ifno), this%tatm(ifno) + DCTOK, this%atmc(ifno))
    !
    ! -- shade-altered above-channel emissivity [Eq. 3, Fogg et al. (2023)]
    emisss = this%epss(this%shd(ifno), emissa, this%emissr)
    !
    ! -- long wave radiation transmitted from the atmosphere to the water surface (A.12)
    lwratm = this%calc_lwr(emisss, this%tatm(ifno))
    !
    ! -- long wave radiation transmitted from water surface to the atmosphere (A.13)
    lwrstrm = this%calc_lwr(-this%emissw, tstrm)
    !
    ! -- longwave radiation heat flux
    lwrflx = lwratm * (1 - this%lwrefl) + lwrstrm
  end subroutine lwr_cq

  !> @brief Deallocate package memory
  !!
  !! Deallocate TVK package scalars and arrays.
  !<
  subroutine lwr_da(this)
    ! -- dummy
    class(LwrType) :: this
    !
    ! -- deallocate time series
    nullify (this%shd)
    nullify (this%lwrefl)
    nullify (this%tatm)
    nullify (this%emissw)
    nullify (this%emissr)
    nullify (this%atmc)
    nullify (this%rh)
    nullify (this%ea)
    nullify (this%es)
    nullify (this%ew)
    !
    ! -- deallocate parent
    call pbstbase_da(this)
  end subroutine lwr_da

  !> @brief Calculate longwave radiation
  !!
  !! Function for calculating incoming or outgoing longwave radiation
  !<
  function calc_lwr(this, eps, temp) result(lwr)
    ! -- dummy
    class(LwrType) :: this
    real(DP) :: eps !< epsilon, representing either emissivity of the atmosphere or the shade-weighted emissivity of the atm
    real(DP) :: temp !< temperature, representing either the temperature of the stream or the atmosphere
    ! -- return
    real(DP) :: lwr !< longwave radiation
    !
    ! -- generalized equation for calculating longwave radiation
    lwr = eps * DSTEFANBOLTZMANN * temp**DFOUR
  end function calc_lwr

end module LongwaveModule
