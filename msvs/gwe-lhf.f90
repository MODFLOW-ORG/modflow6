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
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use PbstBaseModule, only: PbstBaseType, pbstbase_da

  implicit none
  
  private

  public :: LhfType
  public :: lhf_cr

  character(len=16) :: text = '          LHF'

  type, extends(PbstBaseType) :: LhfType
    real(DP), pointer :: wfslope => null() !< ABC wind function slope
    real(DP), pointer :: wfint => null() !< ABC wind function intercept 
    real(DP), dimension(:), pointer, contiguous :: wspd => null() !< ABC wind speed
    real(DP), dimension(:), pointer, contiguous :: tatm => null() !< ABC temperature of the atmosphere
    real(DP), dimension(:), pointer, contiguous :: rh => null() !< ABC relative humidity

  contains

    procedure :: da => lhf_da
    procedure :: ar_set_pointers => lhf_ar_set_pointers
    procedure :: read_option => lhf_read_option
    procedure :: get_pointer_to_value => lhf_get_pointer_to_value
    !procedure :: pbst_options => lhf_options
    !procedure :: subpck_set_stressperiod => lhf_set_stressperiod
    !procedure :: pbst_allocate_arrays => lhf_allocate_arrays
    !procedure, private :: lhf_allocate_scalars
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
    ! -- formats
    character(len=*), parameter :: fmtlhf = &
      "(1x,/1x,'LHF -- LATENT HEAT FLUX PACKAGE, VERSION 1, 08/19/2025', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! -- Print a message identifying the LHF package
    write (this%iout, fmtlhf) this%inunit
    !
    ! -- Set pointers to other package variables
    ! -- ABC
    abcMemoryPath = create_mem_path(this%name_model, 'ABC')
    call mem_setptr(this%wfslope, 'WFSLOPE', abcMemoryPath)
    call mem_setptr(this%wfint, 'WFINT', abcMemoryPath)
    call mem_setptr(this%wspd, 'WSPD', abcMemoryPath)
    call mem_setptr(this%tatm, 'TATM', abcMemoryPath)
    call mem_setptr(this%rh, 'RH', abcMemoryPath)
   
  end subroutine lhf_ar_set_pointers
  
  !> @brief Get an array value pointer given a variable name and node index
  !!
  !! Return a pointer to the given node's value in the appropriate ABC array
  !! based on the given variable name string.
  !<
  function lhf_get_pointer_to_value(this, n, varName) result(bndElem)
    ! -- dummy
    class(LhfType) :: this
    integer(I4B), intent(in) :: n
    character(len=*), intent(in) :: varName
    ! -- return
    real(DP), pointer :: bndElem
    !
    select case (varName)
    case ('TATM')
      bndElem => this%tatm(n)
    case ('WSPD')
      bndElem => this%wspd(n)
    case ('RH')
      bndElem => this%rh(n)  
    case default
      bndElem => null()
    end select
  end function lhf_get_pointer_to_value


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
    real(DP) :: latent_heat_vap
    real(DP) :: sat_vap_tw
    real(DP) :: sat_vap_ta
    real(DP) :: amb_vap_atm
    real(DP) :: evap_rate
    !
    ! -- calculate latent heat flux using mass transfer equation
    !! EQ BASED ON TEMPERATURES IN CELSIUS
    ! -- temperature dependent latent heat of vaporization:
    latent_heat_vap = (2499.64 - (2.51*tstrm))*1000 ! tstrm in C
    
    sat_vap_tw = 6.1275*exp(17.2693882*((tstrm)/(tstrm + DCTOK - 35.86)))
    sat_vap_ta = 6.1275*exp(17.2693882*((this%tatm(ifno))/(this%tatm(ifno) + DCTOK - 35.86)))
    
    amb_vap_atm = (this%rh(ifno)/100)*sat_vap_ta
    
    evap_rate = (this%wfint + this%wfslope*this%wspd(ifno))*(sat_vap_tw - amb_vap_atm)
    
    lhflx = evap_rate * latent_heat_vap * rhow
  end subroutine lhf_cq

  !> @brief Deallocate package memory
  !!
  !! Deallocate TVK package scalars and arrays.
  !<
  subroutine lhf_da(this)
    ! -- dummy
    class(LhfType) :: this
    !  
    ! -- Nullify pointers to other package variables
    nullify (this%wfint)
    nullify (this%wfslope)
    nullify (this%wspd)
    nullify (this%tatm)
    nullify (this%rh)
    !
    ! -- Deallocate parent
    call pbstbase_da(this)
  end subroutine lhf_da

  !> @brief Read a LHF-specific option from the OPTIONS block
  !!
  !! Process a single LHF-specific option. Used when reading the OPTIONS block
  !! of the LHF package input file.
  !<
  function lhf_read_option(this, keyword) result(success)
    ! -- dummy
    class(LhfType) :: this
    character(len=*), intent(in) :: keyword
    ! -- return
    logical :: success
    !
    ! -- There are no LHF-specific options, so just return false
    success = .false.
  end function lhf_read_option

end module LatHeatModule