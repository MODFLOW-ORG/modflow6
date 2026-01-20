!> @brief This module contains the time-varying storage package methods
!!
!! This module contains the methods used to allow storage parameters in the
!! STO package (specific storage and specific yield) to be varied throughout
!! a simulation.
!!
!<
module TvsModule
  use BaseDisModule, only: DisBaseType
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, DZERO
  use KindModule, only: I4B, DP
  use MemoryManagerModule, only: mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use TvBaseModule, only: TvBaseType, tvbase_da

  implicit none

  private

  public :: TvsType
  public :: tvs_cr

  type, extends(TvBaseType) :: TvsType
    integer(I4B), pointer :: integratechanges => null() !< STO flag indicating if mid-simulation ss and sy changes should be integrated via an additional matrix formulation term
    integer(I4B), pointer :: iusesy => null() !< STO flag set if any cell is convertible (0, 1)
    real(DP), dimension(:), pointer, contiguous :: ss => null() !< STO specific storage or storage coefficient
    real(DP), dimension(:), pointer, contiguous :: sy => null() !< STO specific yield

  contains

    procedure :: da => tvs_da
    procedure :: ar_set_pointers => tvs_ar_set_pointers
    procedure :: read_option => tvs_read_option
    procedure :: get_pointer_to_value => tvs_get_pointer_to_value
    procedure :: set_changed_at => tvs_set_changed_at
    procedure :: reset_change_flags => tvs_reset_change_flags
    procedure :: validate_change => tvs_validate_change
  end type TvsType

contains

  !> @brief Create a new TvsType object
  !!
  !! Create a new time-varying storage (TVS) object.
  !<
  subroutine tvs_cr(tvs, name_model, inunit, iout)
    ! -- dummy
    type(TvsType), pointer, intent(out) :: tvs
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    !
    allocate (tvs)
    call tvs%init(name_model, 'TVS', 'TVS', inunit, iout)
  end subroutine tvs_cr

  !> @brief Announce package and set pointers to variables
  !!
  !! Announce package version, set array and variable pointers from the STO
  !! package for access by TVS, and enable storage change integration.
  !<
  subroutine tvs_ar_set_pointers(this)
    ! -- dummy
    class(TvsType) :: this
    ! -- local
    character(len=LENMEMPATH) :: stoMemoryPath
    ! -- formats
    character(len=*), parameter :: fmttvs = &
      "(1x,/1x,'TVS -- TIME-VARYING S PACKAGE, VERSION 1, 08/18/2021', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! -- Print a message identifying the TVS package
    write (this%iout, fmttvs) this%inunit
    !
    ! -- Set pointers to other package variables
    ! -- STO
    stoMemoryPath = create_mem_path(this%name_model, 'STO')
    call mem_setptr(this%integratechanges, 'INTEGRATECHANGES', stoMemoryPath)
    call mem_setptr(this%iusesy, 'IUSESY', stoMemoryPath)
    call mem_setptr(this%ss, 'SS', stoMemoryPath)
    call mem_setptr(this%sy, 'SY', stoMemoryPath)
    !
    ! -- Instruct STO to integrate storage changes, since TVS is active
    this%integratechanges = 1
  end subroutine tvs_ar_set_pointers

  !> @brief Read a TVS-specific option from the OPTIONS block
  !!
  !! Process a single TVS-specific option. Used when reading the OPTIONS block
  !! of the TVS package input file.
  !<
  function tvs_read_option(this, keyword) result(success)
    ! -- dummy
    class(TvsType) :: this
    character(len=*), intent(in) :: keyword
    ! -- return
    logical :: success
    ! -- formats
    character(len=*), parameter :: fmtdsci = &
      "(4X, 'DISABLE_STORAGE_CHANGE_INTEGRATION OPTION:', /, 1X, &
      &'Storage derivative terms will not be added to STO matrix formulation')"
    !
    select case (keyword)
    case ('DISABLE_STORAGE_CHANGE_INTEGRATION')
      success = .true.
      this%integratechanges = 0
      write (this%iout, fmtdsci)
    case default
      success = .false.
    end select
  end function tvs_read_option

  !> @brief Get an array value pointer given a variable name and node index
  !!
  !! Return a pointer to the given node's value in the appropriate STO array
  !! based on the given variable name string.
  !<
  function tvs_get_pointer_to_value(this, n, varName) result(bndElem)
    ! -- dummy
    class(TvsType) :: this
    integer(I4B), intent(in) :: n
    character(len=*), intent(in) :: varName
    ! -- return
    real(DP), pointer :: bndElem
    !
    select case (varName)
    case ('SS')
      bndElem => this%ss(n)
    case ('SY')
      bndElem => this%sy(n)
    case default
      bndElem => null()
    end select
  end function tvs_get_pointer_to_value

  !> @brief Mark property changes as having occurred at (kper, kstp)
  !!
  !! Deferred procedure implementation called by the TvBaseType code when a
  !! property value change occurs at (kper, kstp).
  !<
  subroutine tvs_set_changed_at(this, kper, kstp)
    ! -- dummy
    class(TvsType) :: this
    integer(I4B), intent(in) :: kper
    integer(I4B), intent(in) :: kstp
    !
    ! -- No need to record TVS/STO changes, as no other packages cache
    ! -- Ss or Sy values
  end subroutine tvs_set_changed_at

  !> @brief Clear all per-node change flags
  !!
  !! Deferred procedure implementation called by the TvBaseType code when a
  !! new time step commences, indicating that any previously set per-node
  !! property value change flags should be reset.
  !<
  subroutine tvs_reset_change_flags(this)
    ! -- dummy
    class(TvsType) :: this
    !
    ! -- No need to record TVS/STO changes, as no other packages cache
    ! -- Ss or Sy values
  end subroutine tvs_reset_change_flags

  !> @brief Check that a given property value is valid
  !!
  !! Deferred procedure implementation called by the TvBaseType code after a
  !! property value change occurs. Check if the property value of the given
  !! variable at the given node is invalid, and log an error if so.
  !<
  subroutine tvs_validate_change(this, n, varName)
    ! -- dummy
    class(TvsType) :: this
    integer(I4B), intent(in) :: n
    character(len=*), intent(in) :: varName
    ! -- local
    character(len=LINELENGTH) :: cellstr
    ! -- formats
    character(len=*), parameter :: fmtserr = &
      "(1x, a, ' changed storage property ', a, ' is < 0 for cell ', a,' ', &
      &1pg15.6)"
    character(len=*), parameter :: fmtsyerr = &
      "(1x, a, ' cannot change ', a ,' for cell ', a, ' because SY is unused &
      &in this model (all ICONVERT flags are 0).')"
    !
    ! -- Check the changed value is ok and convert to storage capacity
    if (varName == 'SS') then
      if (this%ss(n) < DZERO) then
        call this%dis%noder_to_string(n, cellstr)
        write (errmsg, fmtserr) trim(adjustl(this%packName)), 'SS', &
          trim(cellstr), this%ss(n)
        call store_error(errmsg)
      end if
    elseif (varName == 'SY') then
      if (this%iusesy /= 1) then
        call this%dis%noder_to_string(n, cellstr)
        write (errmsg, fmtsyerr) trim(adjustl(this%packName)), 'SY', &
          trim(cellstr)
        call store_error(errmsg)
      elseif (this%sy(n) < DZERO) then
        call this%dis%noder_to_string(n, cellstr)
        write (errmsg, fmtserr) trim(adjustl(this%packName)), 'SY', &
          trim(cellstr), this%sy(n)
        call store_error(errmsg)
      end if
    end if
  end subroutine tvs_validate_change

  !> @brief Deallocate package memory
  !!
  !! Deallocate TVS package scalars and arrays.
  !<
  subroutine tvs_da(this)
    ! -- dummy
    class(TvsType) :: this
    !
    ! -- Nullify pointers to other package variables
    nullify (this%integratechanges)
    nullify (this%iusesy)
    nullify (this%ss)
    nullify (this%sy)
    !
    ! -- Deallocate parent
    call tvbase_da(this)
  end subroutine tvs_da

end module TvsModule
