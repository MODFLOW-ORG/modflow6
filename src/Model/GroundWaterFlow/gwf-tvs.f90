!> @brief This module contains the time-varying storage package methods
!!
!! This module contains the methods used to allow storage parameters in the
!! STO package (specific storage and specific yield) to be varied throughout
!! a simulation.
!!
!<
module TvsModule
  use BaseDisModule, only: DisBaseType
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, DZERO, DNODATA
  use KindModule, only: I4B, DP
  use MemoryManagerModule, only: mem_setptr, get_isize
  use MemoryHelperModule, only: create_mem_path
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use TdisModule, only: kper
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
    real(DP), dimension(:), pointer, contiguous :: ss_src => null() !< input SS values
    real(DP), dimension(:), pointer, contiguous :: sy_src => null() !< input SY values

  contains

    procedure :: da => tvs_da
    procedure :: ar_set_pointers => tvs_ar_set_pointers
    procedure :: source_package_options => tvs_source_package_options
    procedure :: apply_row_changes => tvs_apply_row_changes
    procedure :: set_changed_at => tvs_set_changed_at
    procedure :: reset_change_flags => tvs_reset_change_flags
    procedure :: validate_change => tvs_validate_change
  end type TvsType

contains

  !> @brief Create a new TvsType object
  !!
  !! Create a new time-varying storage (TVS) object.
  !<
  subroutine tvs_cr(tvs, name_model, mempath, inunit, iout)
    ! -- dummy
    type(TvsType), pointer, intent(out) :: tvs
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    !
    allocate (tvs)
    call tvs%init(name_model, 'TVS', 'TVS', mempath, inunit, iout)
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
      &' INPUT READ FROM MEMPATH ', A, //)"
    !
    write (this%iout, fmttvs) this%input_mempath
    !
    stoMemoryPath = create_mem_path(this%name_model, 'STO')
    call mem_setptr(this%integratechanges, 'INTEGRATECHANGES', stoMemoryPath)
    call mem_setptr(this%iusesy, 'IUSESY', stoMemoryPath)
    call mem_setptr(this%ss, 'SS', stoMemoryPath)
    call mem_setptr(this%sy, 'SY', stoMemoryPath)
    !
    ! -- Instruct STO to integrate storage changes, since TVS is active
    this%integratechanges = 1
    !
    ! -- set input mempath pointers
    call mem_setptr(this%ss_src, 'SS', this%input_mempath)
    call mem_setptr(this%sy_src, 'SY', this%input_mempath)
  end subroutine tvs_ar_set_pointers

  !> @brief Source TVS-specific options from the input memory path.
  !<
  subroutine tvs_source_package_options(this)
    ! -- dummy
    class(TvsType) :: this
    ! -- locals
    integer(I4B) :: isize
    ! -- formats
    character(len=*), parameter :: fmtdsci = &
      "(4X, 'DISABLE_STORAGE_CHANGE_INTEGRATION OPTION:', /, 6X, &
      &'Storage derivative terms will not be added to STO matrix formulation')"
    !
    ! -- DISABLE_STORAGE_CHANGE_INTEGRATION is a keyword; check via get_isize
    call get_isize('DISABLE_SC_INT', this%input_mempath, isize)
    if (isize > 0) then
      this%integratechanges = 0
      write (this%iout, fmtdsci)
    end if
  end subroutine tvs_source_package_options

  !> @brief Apply input SS/SY column changes for period-data row n to node.
  !<
  subroutine tvs_apply_row_changes(this, n, node)
    ! -- dummy
    class(TvsType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: node
    ! -- local
    character(len=LINELENGTH) :: cellstr
    ! -- formats
    character(len=*), parameter :: fmtvalchg = &
      "(a, ' package: Setting ', a, ' value for cell ', a, ' at start of &
      &stress period ', i0, ' = ', g12.5)"
    !
    if (this%ss_src(n) /= DNODATA) then
      this%ss(node) = this%ss_src(n)
      call this%validate_change(node, 'SS')
      if (this%iprpak /= 0) then
        call this%dis%noder_to_string(node, cellstr)
        write (this%iout, fmtvalchg) &
          trim(adjustl(this%packName)), 'SS', trim(cellstr), kper, this%ss(node)
      end if
    end if
    !
    if (this%sy_src(n) /= DNODATA) then
      this%sy(node) = this%sy_src(n)
      call this%validate_change(node, 'SY')
      if (this%iprpak /= 0) then
        call this%dis%noder_to_string(node, cellstr)
        write (this%iout, fmtvalchg) &
          trim(adjustl(this%packName)), 'SY', trim(cellstr), kper, this%sy(node)
      end if
    end if
  end subroutine tvs_apply_row_changes

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
    nullify (this%integratechanges)
    nullify (this%iusesy)
    nullify (this%ss)
    nullify (this%sy)
    nullify (this%ss_src)
    nullify (this%sy_src)
    call tvbase_da(this)
  end subroutine tvs_da

end module TvsModule
