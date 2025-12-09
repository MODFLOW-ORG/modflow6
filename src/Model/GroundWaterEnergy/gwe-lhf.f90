!> @brief This module contains methods for calculating latent heat flux
!!
!! This module contains the methods used to calculate the latent heat flux
!! for surface-water boundaries, like streams and lakes.  In its current form,
!! this class acts like a package to a package, similar to the TVK package that
!! can be invoked from the NPF package.  Once this package is completed in its
!! prototyped form, it will likely be moved around.
!<

module LatHeatModule
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, DZERO, LENVARNAME, DCTOK
  use KindModule, only: I4B, DP
  use MemoryManagerModule, only: mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use TimeSeriesLinkModule, only: TimeSeriesLinkType
  use TimeSeriesManagerModule, only: TimeSeriesManagerType, tsmanager_cr
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use PbstBaseModule, only: PbstBaseType, pbstbase_da

  implicit none

  public :: LhfType
  public :: lhf_cr

  character(len=16) :: text = '          LHF'

  type, extends(PbstBaseType) :: LhfType

    real(DP), pointer :: wfslope => null() !< wind function slope
    real(DP), pointer :: wfint => null() !< wind function intercept
    real(DP), dimension(:), pointer, contiguous :: wspd => null() !< wind speed
    real(DP), dimension(:), pointer, contiguous :: tatm => null() !< temperature of the atmosphere
    real(DP), dimension(:), pointer, contiguous :: ea => null() !< saturation vapor pressure at air temperature
    real(DP), dimension(:), pointer, contiguous :: ew => null() !< saturation vapor pressure at water temperature

  contains

    procedure :: da => lhf_da
    procedure :: pbst_ar => lhf_ar_set_pointers
    procedure, public :: lhf_cq

  end type LhfType

contains

  !> @brief Create a new LhfType object
  !!
  !! Create a new latent heat flux (LhfType) object. Initially for use with
  !! the SFE-ABC package.
  !<
  subroutine lhf_cr(this, name_model, inunit, iout, ncv)
    ! -- dummy
    type(LhfType), pointer, intent(out) :: this
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    integer(I4B), target, intent(in) :: ncv
    !
    allocate (this)
    call this%init(name_model, 'LHF', 'LHF', inunit, iout, ncv)
    this%text = text
    !
  end subroutine lhf_cr

  !> @brief Announce package and set pointers to variables
  !!
  !! Announce package version and set array and variable pointers from the ABC
  !! package for access by LHF.
  !<
  subroutine lhf_ar_set_pointers(this)
    ! -- dummy
    class(LhfType) :: this
    ! -- local
    character(len=LENMEMPATH) :: abcMemoryPath
    !
    ! -- print a message noting that the SHF utility is active
    write (this%iout, '(a)') &
      'LHF -- LATENT HEAT WILL BE INCLUDED IN THE ATMOSPHERIC BOUNDARY '// &
      'CONDITIONS FOR THE STREAMFLOW ENERGY TRANSPORT PACKAGE'
    !
    ! -- set pointers to variables hosted in the ABC package
    abcMemoryPath = create_mem_path(this%name_model, 'ABC')
    call mem_setptr(this%wfslope, 'WFSLOPE', abcMemoryPath)
    call mem_setptr(this%wfint, 'WFINT', abcMemoryPath)
    call mem_setptr(this%wspd, 'WSPD', abcMemoryPath)
    call mem_setptr(this%tatm, 'TATM', abcMemoryPath)
    call mem_setptr(this%ea, 'EA', abcMemoryPath)
    call mem_setptr(this%ew, 'EW', abcMemoryPath)
    !
    ! -- create time series manager
    call tsmanager_cr(this%tsmanager, this%iout, &
                      removeTsLinksOnCompletion=.true., &
                      extendTsToEndOfSimulation=.true.)
  end subroutine lhf_ar_set_pointers

  !> @brief Calculate Latent Heat Flux
  !!
  !! Calculate and return the latent heat flux for one reach
  !<
  subroutine lhf_cq(this, ifno, tstrm, rhow, lhflx)
    ! -- dummy
    class(LhfType), intent(inout) :: this
    integer(I4B), intent(in) :: ifno !< stream reach integer id
    real(DP), intent(in) :: tstrm !< temperature of the stream reach
    real(DP), intent(in) :: rhow !< density of water
    real(DP), intent(inout) :: lhflx !< calculated latent heat flux amount
    ! -- local
    real(DP) :: l !< latent heat vaporization
    real(DP) :: sat_vap_tw
    real(DP) :: sat_vap_ta
    real(DP) :: amb_vap_atm
    real(DP) :: evap
    !
    ! -- calculate latent heat of vaporization (water temperature dependent) Eq. A.16
    l = 2499.64_DP - (2.51_DP * tstrm) ! tstrm must be in degrees C for now
    !
    ! -- mass-transfer method for calculating evap rate (A.17)
    evap = (this%wfint + this%wfslope * this%wspd(ifno)) * &
           (this%ew(ifno) - this%ea(ifno))
    !
    ! -- calculate latent heat flux (A.15)
    lhflx = evap * l * rhow
  end subroutine lhf_cq

  !> @brief Deallocate package memory
  !!
  !! Deallocate TVK package scalars and arrays.
  !<
  subroutine lhf_da(this)
    ! -- dummy
    class(LhfType) :: this
    !
    ! -- nullify pointers to other package variables
    nullify (this%wfint)
    nullify (this%wfslope)
    nullify (this%wspd)
    nullify (this%tatm)
    nullify (this%ea)
    nullify (this%ew)
    !
    ! -- deallocate parent
    call pbstbase_da(this)
  end subroutine lhf_da

end module LatHeatModule
