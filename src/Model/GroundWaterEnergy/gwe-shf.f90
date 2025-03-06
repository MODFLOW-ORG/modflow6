!> @brief This module contains methods for calculating sensible heat flux
!!
!! This module contains the methods used to calculate the sensible heat flux
!! for surface-water boundaries, like streams and lakes.  In its current form,
!! this class acts like a package to a package, similar to the TVK package that
!! can be invoked from the NPF package.  Once this package is completed in its
!! prototyped form, it will likely be moved around.
!<

module SensHeatModule
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, DZERO, LENVARNAME
  use KindModule, only: I4B, DP
  use MemoryManagerModule, only: mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use PbstBaseModule, only: PbstBaseType, pbstbase_da

  implicit none
  
  private

  public :: ShfType
  public :: shf_cr

  character(len=16) :: text = '          SHF'

  type, extends(PbstBaseType) :: ShfType
    real(DP), pointer :: rhoa => null() !< ABC density of air
    real(DP), pointer :: cpa => null() !< ABC heat capacity of air
    real(DP), pointer :: cd => null() !< ABC drag coefficient
    real(DP), dimension(:), pointer, contiguous :: wspd => null() !< ABC wind speed
    real(DP), dimension(:), pointer, contiguous :: tatm => null() !< ABC temperature of the atmosphere

  contains

    procedure :: da => shf_da
    procedure :: ar_set_pointers => shf_ar_set_pointers
    procedure :: read_option => shf_read_option
    procedure :: get_pointer_to_value => shf_get_pointer_to_value
    !procedure :: pbst_options => shf_options
    !procedure :: subpck_set_stressperiod => shf_set_stressperiod
    !procedure :: pbst_allocate_arrays => shf_allocate_arrays
    !procedure, private :: shf_allocate_scalars
    procedure, public :: shf_cq

  end type ShfType

contains

  !> @brief Create a new ShfType object
  !!
  !! Create a new sensible heat flux (ShfType) object. Initially for use with
  !! the SFE-ABC package.
  !<
  subroutine shf_cr(this, name_model, inunit, iout, ncv)
    ! -- dummy
    type(ShfType), pointer, intent(out) :: this
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    integer(I4B), target, intent(in) :: ncv
    !
    allocate (this)
    call this%init(name_model, 'SHF', 'SHF', inunit, iout, ncv)
    this%text = text
    !
    ! -- allocate scalars
    !call shf%shf_allocate_scalars()
  end subroutine shf_cr

  !> @brief Announce package and set pointers to variables
  !!
  !! Announce package version and set array and variable pointers from the ABC
  !! package for access by SHF.
  !<
  subroutine shf_ar_set_pointers(this)
    ! -- dummy
    class(ShfType) :: this
    ! -- local
    character(len=LENMEMPATH) :: abcMemoryPath
    ! -- formats
    character(len=*), parameter :: fmtshf = &
      "(1x,/1x,'SHF -- SENSIBLE HEAT FLUX PACKAGE, VERSION 1, 03/12/2025', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! -- Print a message identifying the SHF package
    write (this%iout, fmtshf) this%inunit
    !
    ! -- Set pointers to other package variables
    ! -- ABC
    abcMemoryPath = create_mem_path(this%name_model, 'ABC')
    call mem_setptr(this%rhoa, 'RHOA', abcMemoryPath)
    call mem_setptr(this%cpa, 'CPA', abcMemoryPath)
    call mem_setptr(this%cd, 'CD', abcMemoryPath)
    call mem_setptr(this%wspd, 'WSPD', abcMemoryPath)
    call mem_setptr(this%tatm, 'TATM', abcMemoryPath)
   
  end subroutine shf_ar_set_pointers
  
  !> @brief Get an array value pointer given a variable name and node index
  !!
  !! Return a pointer to the given node's value in the appropriate ABC array
  !! based on the given variable name string.
  !<
  function shf_get_pointer_to_value(this, n, varName) result(bndElem)
    ! -- dummy
    class(ShfType) :: this
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
    case default
      bndElem => null()
    end select
  end function shf_get_pointer_to_value

  
  !> @brief Allocate scalars specific to the streamflow energy transport (SFE)
  !! package.
  !<
  !subroutine shf_allocate_scalars(this)
  !  ! -- modules
  !  use MemoryManagerModule, only: mem_allocate
  !  ! -- dummy
  !  class(ShfType) :: this
  !  !
  !  ! -- allocate
  !  call mem_allocate(this%rhoa, 'RHOA', this%memoryPath)
  !  call mem_allocate(this%cpa, 'CPA', this%memoryPath)
  !  call mem_allocate(this%cd, 'CD', this%memoryPath)
  !  !
  !  ! -- initialize to default values
  !  this%rhoa = 1.225 ! kg/m3
  !  this%cpa = 717.0 ! J/kg/C
  !  this%cd = 0.002 ! unitless
  !end subroutine shf_allocate_scalars

  !> @brief Allocate arrays specific to the sensible heat flux (SHF) package
  !<
  !subroutine shf_allocate_arrays(this)
  !  ! -- modules
  !  use MemoryManagerModule, only: mem_allocate
  !  ! -- dummy
  !  class(ShfType), intent(inout) :: this
  !  ! -- local
  !  integer(I4B) :: n
  !  !
  !  ! -- time series
  !  call mem_allocate(this%wspd, this%ncv, 'WSPD', this%memoryPath)
  !  call mem_allocate(this%tatm, this%ncv, 'TATM', this%memoryPath)
  !  !
  !  ! -- initialize
  !  do n = 1, this%ncv
  !    this%wspd(n) = DZERO
  !    this%tatm(n) = DZERO
  !  end do
  !end subroutine

  !> @brief Set options specific to the ShfType
  !!
  !! This routine overrides PbstBaseType%bnd_options
  !<
<<<<<<< HEAD
  !subroutine shf_options(this, option, found)
  !  ! -- dummy
  !  class(ShfType), intent(inout) :: this
  !  character(len=*), intent(inout) :: option
  !  logical, intent(inout) :: found
  !  !
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
  !end subroutine shf_options

  !> @brief Calculate Sensible Heat Flux
  !!
  !! Calculate and return the sensible heat flux for one reach
  !<
  subroutine shf_cq(this, ifno, tstrm, shflx)
    ! -- dummy
    class(ShfType), intent(inout) :: this
    integer(I4B), intent(in) :: ifno !< stream reach integer id
    real(DP), intent(in) :: tstrm !< temperature of the stream reach
    real(DP), intent(inout) :: shflx !< calculated sensible heat flux amount
    ! -- local
    real(DP) :: shf_const
    !
    ! -- calculate sensible heat flux using HGS equation
    shf_const = this%cd * this%cpa * this%rhoa
    shflx = shf_const * this%wspd(ifno) * (this%tatm(ifno) - tstrm)
  end subroutine shf_cq

=======
  subroutine shf_options(this, option, found)
    ! -- dummy
    class(ShfType), intent(inout) :: this
    character(len=*), intent(inout) :: option
    logical, intent(inout) :: found
    !
    found = .true.
    select case (option)
    case ('DENSITY_AIR')
      this%rhoa = this%parser%GetDouble()
      if (this%rhoa <= 0.0) then
        write (errmsg, '(a)') 'Specified value for the density of &
          &the atmosphere must be greater than 0.0.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      else
        write (this%iout, '(4x,a,1pg15.6)') &
          "The density of the atmosphere has been set to: ", this%rhoa
      end if
    case ('HEAT_CAPACITY_AIR')
      this%cpa = this%parser%GetDouble()
      if (this%cpa <= 0.0) then
        write (errmsg, '(a)') 'Specified value for the heat capacity of &
          &the atmosphere must be greater than 0.0.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      else
        write (this%iout, '(4x,a,1pg15.6)') &
          "The heat capacity of the atmosphere has been set to: ", this%cpa
      end if
    case ('DRAG_COEFFICIENT')
      this%cd = this%parser%GetDouble()
      if (this%cd <= 0.0) then
        write (errmsg, '(a)') 'Specified value for the drag coefficient &
          &must be greater than 0.0.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      else
        write (this%iout, '(4x,a,1pg15.6)') &
          "The heat capacity of the atmosphere has been set to: ", this%cpa
      end if
    case default
      write (errmsg, '(a,a)') 'Unknown SHF option: ', trim(option)
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end select
  end subroutine shf_options
  
>>>>>>> 4baf8469ab (read shf options block; will need to circle back and test TS functionality later)
  !> @brief Deallocate package memory
  !!
  !! Deallocate TVK package scalars and arrays.
  !<
  subroutine shf_da(this)
  !  ! -- modules
  !  use MemoryManagerModule, only: mem_deallocate
  !  ! -- dummy
    class(ShfType) :: this
  !  !
  !  ! -- Nullify pointers to other package variables
  !  call mem_deallocate(this%rhoa)
  !  call mem_deallocate(this%cpa)
  !  call mem_deallocate(this%cd)
  !  !
  !  ! -- Deallocate time series
  !  call mem_deallocate(this%wspd)
  !  call mem_deallocate(this%tatm)
    ! -- Nullify pointers to other package variables
    nullify (this%rhoa)
    nullify (this%cpa)
    nullify (this%cd)
    nullify (this%wspd)
    nullify (this%tatm)
  !  !
  !  ! -- Deallocate parent
    call pbstbase_da(this)
  end subroutine shf_da

  !> @brief Read a SHF-specific option from the OPTIONS block
  !!
  !! Process a single SHF-specific option. Used when reading the OPTIONS block
  !! of the SHF package input file.
  !<
  function shf_read_option(this, keyword) result(success)
    ! -- dummy
    class(ShfType) :: this
    character(len=*), intent(in) :: keyword
    ! -- return
    logical :: success
    !
    ! -- There are no SHF-specific options, so just return false
    success = .false.
  end function shf_read_option

  !> @brief Set the stress period attributes based on the keyword
  !<
!  subroutine shf_set_stressperiod(this, itemno, keyword, found)
!    ! -- module
!    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
!    ! -- dummy
!    class(ShfType), intent(inout) :: this
!    integer(I4B), intent(in) :: itemno
!    character(len=*), intent(in) :: keyword
!    logical, intent(inout) :: found
!    ! -- local
!    character(len=LINELENGTH) :: text
!    integer(I4B) :: ierr
!    integer(I4B) :: jj
!    real(DP), pointer :: bndElem => null()
!    !
!    ! <wspd> WIND SPEED
!    ! <tatm> TEMPERATURE OF THE ATMOSPHERE
!    !
!    found = .true.
!    select case (keyword)
!    case ('WSPD')
!      ierr = this%pbst_check_valid(itemno)
!      if (ierr /= 0) then
!        goto 999
!      end if
!      call this%parser%GetString(text)
!      jj = 1
!      bndElem => this%wspd(itemno)
!      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
!                                         this%packName, 'BND', this%tsManager, &
!                                         this%iprpak, 'WSPD')
!    case ('TATM')
!      ierr = this%pbst_check_valid(itemno)
!      if (ierr /= 0) then
!        goto 999
!      end if
!      call this%parser%GetString(text)
!      jj = 1
!      bndElem => this%tatm(itemno)
!      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
!                                         this%packName, 'BND', this%tsManager, &
!                                         this%iprpak, 'TATM')
!    case default
!      !
!      ! -- Keyword not recognized so return to caller with found = .false.
!      found = .false.
!    end select
!    !
!999 continue
!  end subroutine shf_set_stressperiod

end module SensHeatModule
