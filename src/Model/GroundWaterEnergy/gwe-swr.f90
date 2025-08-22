!> @brief This module contains methods for calculating the shortwave radiation heat flux
!!
!! This module contains the methods used to calculate the shortwave radiation heat flux
!! for surface-water boundaries, like streams and lakes.  In its current form,
!! this class acts like a package to a package, similar to the TVK package that
!! can be invoked from the NPF package.  Once this package is completed in its
!! prototyped form, it will likely be moved around.
!<

module ShortwaveModule
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, DZERO, LENVARNAME
  use KindModule, only: I4B, DP
  use MemoryManagerModule, only: mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use PbstBaseModule, only: PbstBaseType, pbstbase_da

  implicit none

  private

  public :: SwrType
  public :: swr_cr

  character(len=16) :: text = '          SWR'

  type, extends(PbstBaseType) :: SwrType

    real(DP), dimension(:), pointer, contiguous :: solr => null() !< solar radiation
    real(DP), dimension(:), pointer, contiguous :: shd => null() !< shade fraction
    real(DP), dimension(:), pointer, contiguous :: swrefl => null() !< shortwave reflectance of water surface

  contains

    procedure :: da => swr_da
    procedure :: read_option => swr_read_option
    procedure :: ar_set_pointers => swr_ar_set_pointers
    procedure :: get_pointer_to_value => swr_get_pointer_to_value
    !procedure :: pbst_options => swr_options
    !procedure :: subpck_set_stressperiod => swr_set_stressperiod
    !procedure :: pbst_allocate_arrays => swr_allocate_arrays
    !procedure, private :: swr_allocate_scalars
    procedure, public :: swr_cq

  end type SwrType

contains

  !> @brief Create a new SwrType object
  !!
  !! Create a new shortwave radiation flux (SwrType) object. Initially for use with
  !! the SFE package.
  !<
  subroutine swr_cr(this, name_model, inunit, iout, ncv)
    ! -- dummy
    type(SwrType), pointer, intent(out) :: this
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    integer(I4B), target, intent(in) :: ncv
    !
    allocate (this)
    call this%init(name_model, 'SWR', 'SWR', inunit, iout, ncv)
    this%text = text
    !
    ! -- allocate scalars
    !call swr%swr_allocate_scalars()
  end subroutine swr_cr

   !> @brief Announce package and set pointers to variables
  !!
  !! Announce package version and set array and variable pointers from the ABC
  !! package for access by SWR.
  !<
  subroutine swr_ar_set_pointers(this)
    ! -- dummy
    class(SwrType) :: this
    ! -- local
    character(len=LENMEMPATH) :: abcMemoryPath
    ! -- formats
    character(len=*), parameter :: fmtswr = &
      "(1x,/1x,'SWR -- SHORTWAVE RADIATION PACKAGE, VERSION 1, 05/01/2025', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! -- Print a message identifying the SWR package
    write (this%iout, fmtswr) this%inunit
    !
    ! -- Set pointers to other package variables
    ! -- ABC
    abcMemoryPath = create_mem_path(this%name_model, 'ABC')
    call mem_setptr(this%solr, 'SOLR', abcMemoryPath)
    call mem_setptr(this%shd, 'SHD', abcMemoryPath)
    call mem_setptr(this%swrefl, 'SWREFL', abcMemoryPath)
   
  end subroutine swr_ar_set_pointers
  
  !!> @brief Get an array value pointer given a variable name and node index
  !!!
  !!! Return a pointer to the given node's value in the appropriate ABC array
  !!! based on the given variable name string.
  !!<
  function swr_get_pointer_to_value(this, n, varName) result(bndElem)
      ! -- dummy
      class(SwrType) :: this
      integer(I4B), intent(in) :: n
      character(len=*), intent(in) :: varName
      ! -- return
      real(DP), pointer :: bndElem
      
      select case(varName)
          case default
          bndElem => null()
      end select
    end function
  
  !function swr_get_pointer_to_value(this, n, varName) result(bndElem)
  !  ! -- dummy
  !  class(SwrType) :: this
  !  integer(I4B), intent(in) :: n
  !  character(len=*), intent(in) :: varName
  !  ! -- return
  !  real(DP), pointer :: bndElem
  !  !
  !  select case (varName)
  !  case ('TATM')
  !    bndElem => this%tatm(n)
  !  case ('WSPD')
  !    bndElem => this%wspd(n)
  !  case default
  !    bndElem => null()
  !  end select
  !end function shf_get_pointer_to_value

  !> @brief Allocate scalars specific to the streamflow energy transport (SFE)
  !! package.
  !! NO SCALARS IN THIS YET
  !<
  !subroutine swr_allocate_scalars(this)
  ! -- modules
  !use MemoryManagerModule, only: mem_allocate
  ! -- dummy
  !class(SwrType) :: this
  !
  ! -- allocate
  !call mem_allocate(this%solr, 'SOLR', this%memoryPath)
  !call mem_allocate(this%shd, 'SHD', this%memoryPath)
  !call mem_allocate(this%swrefl, 'SWREFL', this%memoryPath)
  !
  ! -- initialize to default values
  !this%solr = 600 ! W/m3
  !this%shd = 0.25 ! dimensionless fraction
  !this%swrefl = 0.03 ! dimensionless fraction
  !end subroutine swr_allocate_scalars

  !> @brief Allocate arrays specific to the sensible heat flux (SWR) package
  !<
  !subroutine swr_allocate_arrays(this)
  !  ! -- modules
  !  use MemoryManagerModule, only: mem_allocate
  !  ! -- dummy
  !  class(SwrType), intent(inout) :: this
  !  ! -- local
  !  integer(I4B) :: n
  !  !
  !  ! -- time series
  !  call mem_allocate(this%solr, this%ncv, 'SOLR', this%memoryPath)
  !  call mem_allocate(this%shd, this%ncv, 'SHD', this%memoryPath)
  !  call mem_allocate(this%swrefl, this%ncv, 'SWREFL', this%memoryPath)
  !  !
  !  ! -- initialize
  !  do n = 1, this%ncv
  !    this%solr(n) = DZERO
  !    this%shd(n) = DZERO
  !    this%swrefl(n) = DZERO
  !  end do
  !end subroutine

  !> @brief Set options specific to the SwrType
  !!
  ! NOT USED RIGHT NOW BECAUSE NO CONSTANTS
  !! This routine overrides PbstBaseType%bnd_options
  !<
  !subroutine swr_options(this, option, found)
  ! -- dummy
  !  class(SwrType), intent(inout) :: this
  !  character(len=*), intent(inout) :: option
  !  logical, intent(inout) :: found
  !
  !  found = .true.
  !  select case (option)
  !  case ('DENSITY_AIR')
  !    this%rhoa = this%parser%GetDouble()
  !    if (this%rhoa <= 0.0) then
  !      write (errmsg, '(a)') 'Specified value for the density of &
  !        &the atmosphere must be greater than 0.0.'
  !      call store_error(errmsg)
  !      call this%parser%StoreErrorUnit()
  !    else
  !      write (this%iout, '(4x,a,1pg15.6)') &
  !        "The density of the atmosphere has been set to: ", this%rhoa
  !    end if
  !  case ('HEAT_CAPACITY_AIR')
  !    this%cpa = this%parser%GetDouble()
  !    if (this%cpa <= 0.0) then
  !      write (errmsg, '(a)') 'Specified value for the heat capacity of &
  !        &the atmosphere must be greater than 0.0.'
  !      call store_error(errmsg)
  !      call this%parser%StoreErrorUnit()
  !    else
  !      write (this%iout, '(4x,a,1pg15.6)') &
  !        "The heat capacity of the atmosphere has been set to: ", this%cpa
  !    end if
  !  case ('DRAG_COEFFICIENT')
  !    this%cd = this%parser%GetDouble()
  !    if (this%cd <= 0.0) then
  !      write (errmsg, '(a)') 'Specified value for the drag coefficient &
  !        &must be greater than 0.0.'
  !      call store_error(errmsg)
  !      call this%parser%StoreErrorUnit()
  !    else
  !      write (this%iout, '(4x,a,1pg15.6)') &
  !        "The heat capacity of the atmosphere has been set to: ", this%cpa
  !    end if
  !  case default
  !    write (errmsg, '(a,a)') 'Unknown SHF option: ', trim(option)
  !    call store_error(errmsg)
  !    call this%parser%StoreErrorUnit()
  !  end select
  ! end subroutine swr_options

  !> @brief Calculate Shortwave Radiation Heat Flux
  !!
  !! Calculate and return the shortwave radiation heat flux for one reach
  !<
  subroutine swr_cq(this, ifno, swrflx)
    ! -- dummy
    class(SwrType), intent(inout) :: this
    integer(I4B), intent(in) :: ifno !< stream reach integer id
    !real(DP), intent(in) :: tstrm !< temperature of the stream reach
    real(DP), intent(inout) :: swrflx !< calculated shortwave radiation heat flux amount
    ! -- local
    !real(DP) :: swr_const
    !
    ! -- calculate shortwave radiation heat flux (version: user input of sol rad data)
    swrflx = (1 - this%shd(ifno)) * (1 - this%swrefl(ifno)) * this%solr(ifno)
  end subroutine swr_cq

  !> @brief Deallocate package memory
  !!
  !! Deallocate TVK package scalars and arrays.
  !<
  subroutine swr_da(this)
    ! -- modules
    !use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(SwrType) :: this
    !
    ! -- Nullify pointers to other package variables
    !call mem_deallocate(this%rhoa)
    !call mem_deallocate(this%cpa)
    !call mem_deallocate(this%cd)
    !
    ! -- Deallocate time series
    nullify (this%shd)
    nullify (this%swrefl)
    nullify (this%solr)
    !
    ! -- Deallocate parent
    call pbstbase_da(this)
  end subroutine swr_da

  !> @brief Read a SWR-specific option from the OPTIONS block
  !!
  !! Process a single SWR-specific option. Used when reading the OPTIONS block
  !! of the SWR package input file.
  !<
  function swr_read_option(this, keyword) result(success)
    ! -- dummy
    class(SwrType) :: this
    character(len=*), intent(in) :: keyword
    ! -- return
    logical :: success
    !
    ! -- There are no SWR-specific options, so just return false
    success = .false.
  end function swr_read_option

  !> @brief Set the stress period attributes based on the keyword
  !<
!  subroutine swr_set_stressperiod(this, itemno, keyword, found)
!    ! -- module
!    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
!    ! -- dummy
!    class(SwrType), intent(inout) :: this
!    integer(I4B), intent(in) :: itemno
!    character(len=*), intent(in) :: keyword
!    logical, intent(inout) :: found
!    ! -- local
!    character(len=LINELENGTH) :: text
!    integer(I4B) :: ierr
!    integer(I4B) :: jj
!    real(DP), pointer :: bndElem => null()
!    !
!    ! <shd> SHADE
!    ! <swrefl> REFLECTANCE OF SHORTWAVE RADIATION OFF WATER SURFACE
!    ! <solr> SOLAR RADIATION
!    !
!    found = .true.
!    select case (keyword)
!    case ('SHD')
!      ierr = this%pbst_check_valid(itemno)
!      if (ierr /= 0) then
!        goto 999
!      end if
!      call this%parser%GetString(text)
!      jj = 1
!      bndElem => this%shd(itemno)
!      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
!                                         this%packName, 'BND', this%tsManager, &
!                                         this%iprpak, 'SHD')
!    case ('SWREFL')
!      ierr = this%pbst_check_valid(itemno)
!      if (ierr /= 0) then
!        goto 999
!      end if
!      call this%parser%GetString(text)
!      jj = 1
!      bndElem => this%swrefl(itemno)
!      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
!                                         this%packName, 'BND', this%tsManager, &
!                                         this%iprpak, 'SWREFL')
!    case ('SOLR')
!      ierr = this%pbst_check_valid(itemno)
!      if (ierr /= 0) then
!        goto 999
!      end if
!      call this%parser%GetString(text)
!      jj = 1
!      bndElem => this%solr(itemno)
!      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
!                                         this%packName, 'BND', this%tsManager, &
!                                         this%iprpak, 'SOLR')
!    case default
!      !
!      ! -- Keyword not recognized so return to caller with found = .false.
!      found = .false.
!    end select
!    !
!999 continue
!  end subroutine swr_set_stressperiod

end module ShortwaveModule
