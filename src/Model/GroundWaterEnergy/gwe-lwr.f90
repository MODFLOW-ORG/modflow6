!> @brief This module contains methods for calculating the longwave radiation heat flux
!!
!! This module contains the methods used to calculate the longwave radiation heat flux
!! for surface-water boundaries, like streams and lakes.  In its current form,
!! this class acts like a package to a package, similar to the TVK package that
!! can be invoked from the NPF package.  Once this package is completed in its
!! prototyped form, it will likely be moved around.
!<

module LongwaveModule
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, DZERO, LENVARNAME, DSTEFANBOLTZMANN, DCTOK
  use KindModule, only: I4B, DP
  use MemoryManagerModule, only: mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use TimeSeriesLinkModule, only: TimeSeriesLinkType
  use TimeSeriesManagerModule, only: TimeSeriesManagerType, tsmanager_cr
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use PbstBaseModule, only: PbstBaseType, pbstbase_da

  implicit none

  private

  public :: LwrType
  public :: lwr_cr

  character(len=16) :: text = '          LWR'

  type, extends(PbstBaseType) :: LwrType

    real(DP), pointer :: lwrefl => null() !< ABC longwave reflectance of water surface 
    real(DP), pointer :: emissw => null() !< ABC emissivity of water
    real(DP), pointer :: emissr => null() !< ABC emissivity of riparian canopy
    real(DP), dimension(:), pointer, contiguous :: shd => null() !< ABC shade fraction
    real(DP), dimension(:), pointer, contiguous :: atmc => null() !< ABC atmospheric composition adjustment
    real(DP), dimension(:), pointer, contiguous :: tatm => null() !< ABC temperature of the atmosphere
    real(DP), dimension(:), pointer, contiguous :: rh => null() !< ABC relative humidity

  contains

    procedure :: da => lwr_da
    procedure :: pbst_ar => lwr_ar_set_pointers
    procedure :: read_option => lwr_read_option
    !procedure :: get_pointer_to_value => lwr_get_pointer_to_value
    procedure, public :: lwr_cq

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
    !
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
    ! -- formats
    !character(len=*), parameter :: fmtlwr = &
    !
    ! -- Print a message noting that the LWR utility is active
    write (this%iout, '(a)') &
      'LWR -- LONGWAVE RADIATION WILL BE INCLUDED IN THE ATMOSPHERIC '// &
      'BOUNDARY CONDITIONS FOR THE STREAMFLOW ENERGY TRANSPORT PACKAGE'
    !
    ! -- Set pointers to other package variables
    ! -- ABC
    abcMemoryPath = create_mem_path(this%name_model, 'ABC')
    call mem_setptr(this%lwrefl, 'LWREFL', abcMemoryPath)
    call mem_setptr(this%shd, 'SHD', abcMemoryPath)
    call mem_setptr(this%atmc, 'ATMC', abcMemoryPath)
    call mem_setptr(this%tatm, 'TATM', abcMemoryPath)
    call mem_setptr(this%emissw, 'EMISSW', abcMemoryPath)
    call mem_setptr(this%emissr, 'EMISSR', abcMemoryPath)
    call mem_setptr(this%rh, 'RH', abcMemoryPath)
    write(*,*) "Here where set ptr happens"
    !
    ! -- Set the riparian canopy emissivity constant
    this%emissr = 0.97_DP
    !
    ! -- create time series manager
    call tsmanager_cr(this%tsmanager, this%iout, &
                      removeTsLinksOnCompletion=.true., &
                      extendTsToEndOfSimulation=.true.)
  end subroutine lwr_ar_set_pointers
  
  !> @brief Get an array value pointer given a variable name and node index
  !!
  !! Return a pointer to the given node's value in the appropriate ABC array
  !! based on the given variable name string.
  !<
  !function lwr_get_pointer_to_value(this, n, varName) result(bndElem)
  !  ! -- dummy
  !  class(LwrType) :: this
  !  integer(I4B), intent(in) :: n
  !  character(len=*), intent(in) :: varName
  !  ! -- return
  !  real(DP), pointer :: bndElem
  !  
  !  select case(varName)
  !  case ('TATM')
  !    bndElem => this%tatm(n)
  !  case ('ATMC')
  !    bndElem => this%atmc(n)
  !  case ('RH')
  !    bndElem => this%rh(n)  
  !  case default
  !    bndElem => null()
  !  end select
  !end function
  
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
    write(*,*) "value of tatm is ",this%tatm(ifno)
    sat_vap_ta = 6.1275_DP * exp(17.2693882_DP * (this%tatm(ifno) / (this%tatm(ifno) + DCTOK - 35.86_DP)))
    amb_vap_atm = (this%rh(ifno) / 100.0_DP) * sat_vap_ta
    emissa = 1.24_DP * this%atmc(ifno) * ((amb_vap_atm/this%tatm(ifno))**(1/7))
    !
    emisss = (1 - this%shd(ifno)) * emissa + this%shd(ifno) * this%emissr
    !
    lwratm = emisss * DSTEFANBOLTZMANN * this%tatm(ifno)**4
    lwrstrm = -this%emissw * DSTEFANBOLTZMANN * tstrm**4
    !
    ! -- calculate longwave radiation heat flux
    lwrflx = lwratm * (1 - this%lwrefl) + lwrstrm
  end subroutine lwr_cq

  !> @brief Deallocate package memory
  !!
  !! Deallocate TVK package scalars and arrays.
  !<
  subroutine lwr_da(this)
    ! -- modules
    !use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(LwrType) :: this
    !
    ! -- Deallocate time series
    nullify (this%shd)
    nullify (this%lwrefl)
    nullify (this%tatm)
    nullify (this%emissw)
    nullify (this%emissr)
    nullify (this%atmc)
    nullify (this%rh)
    !
    ! -- Deallocate parent
    call pbstbase_da(this)
  end subroutine lwr_da

  !> @brief Read a LWR-specific option from the OPTIONS block
  !!
  !! Process a single LWR-specific option. Used when reading the OPTIONS block
  !! of the LWR package input file.
  !<
  function lwr_read_option(this, keyword) result(success)
    ! -- dummy
    class(LwrType) :: this
    character(len=*), intent(in) :: keyword
    ! -- return
    logical :: success
    !
    ! -- There are no LWR-specific options, so just return false
    success = .false.
  end function lwr_read_option

end module LongwaveModule
