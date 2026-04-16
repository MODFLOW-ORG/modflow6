!> @brief This module contains time-varying conductivity package methods
!!
!! This module contains the methods used to allow hydraulic conductivity
!! parameters in the NPF package (K11, K22, K33) to be varied throughout
!! a simulation.
!!
!<
module TvkModule
  use BaseDisModule, only: DisBaseType
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, DZERO, DNODATA
  use KindModule, only: I4B, DP
  use MemoryManagerModule, only: mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use TdisModule, only: kper
  use TvBaseModule, only: TvBaseType, tvbase_da

  implicit none

  private

  public :: TvkType
  public :: tvk_cr

  type, extends(TvBaseType) :: TvkType
    integer(I4B), pointer :: ik22overk => null() !< NPF flag that k22 is specified as anisotropy ratio
    integer(I4B), pointer :: ik33overk => null() !< NPF flag that k33 is specified as anisotropy ratio
    real(DP), dimension(:), pointer, contiguous :: k11 => null() !< NPF hydraulic conductivity; if anisotropic, then this is Kx prior to rotation
    real(DP), dimension(:), pointer, contiguous :: k22 => null() !< NPF hydraulic conductivity; if specified then this is Ky prior to rotation
    real(DP), dimension(:), pointer, contiguous :: k33 => null() !< NPF hydraulic conductivity; if specified then this is Kz prior to rotation
    integer(I4B), pointer :: kchangeper => null() !< NPF last stress period in which any node K (or K22, or K33) values were changed (0 if unchanged from start of simulation)
    integer(I4B), pointer :: kchangestp => null() !< NPF last time step in which any node K (or K22, or K33) values were changed (0 if unchanged from start of simulation)
    integer(I4B), dimension(:), pointer, contiguous :: nodekchange => null() !< NPF grid array of flags indicating for each node whether its K (or K22, or K33) value changed (1) at (kchangeper, kchangestp) or not (0)
    real(DP), dimension(:), pointer, contiguous :: k11_src => null() !< input K values
    real(DP), dimension(:), pointer, contiguous :: k22_src => null() !< input K22 values
    real(DP), dimension(:), pointer, contiguous :: k33_src => null() !< input K33 values

  contains

    procedure :: da => tvk_da
    procedure :: ar_set_pointers => tvk_ar_set_pointers
    procedure :: apply_row_changes => tvk_apply_row_changes
    procedure :: set_changed_at => tvk_set_changed_at
    procedure :: reset_change_flags => tvk_reset_change_flags
    procedure :: validate_change => tvk_validate_change
  end type TvkType

contains

  !> @brief Create a new TvkType object
  !!
  !! Create a new time-varying conductivity (TvkType) object.
  !<
  subroutine tvk_cr(tvk, name_model, mempath, inunit, iout)
    ! -- dummy
    type(TvkType), pointer, intent(out) :: tvk
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    !
    allocate (tvk)
    call tvk%init(name_model, 'TVK', 'TVK', mempath, inunit, iout)
  end subroutine tvk_cr

  !> @brief Announce package and set pointers to variables
  !!
  !! Announce package version and set array and variable pointers from the NPF
  !! package for access by TVK.
  !<
  subroutine tvk_ar_set_pointers(this)
    ! -- dummy
    class(TvkType) :: this
    ! -- local
    character(len=LENMEMPATH) :: npfMemoryPath
    ! -- formats
    character(len=*), parameter :: fmttvk = &
      "(1x,/1x,'TVK -- TIME-VARYING K PACKAGE, VERSION 1, 08/18/2021', &
      &' INPUT READ FROM MEMPATH ', A, //)"
    !
    write (this%iout, fmttvk) this%input_mempath
    !
    npfMemoryPath = create_mem_path(this%name_model, 'NPF')
    call mem_setptr(this%ik22overk, 'IK22OVERK', npfMemoryPath)
    call mem_setptr(this%ik33overk, 'IK33OVERK', npfMemoryPath)
    call mem_setptr(this%k11, 'K11', npfMemoryPath)
    call mem_setptr(this%k22, 'K22', npfMemoryPath)
    call mem_setptr(this%k33, 'K33', npfMemoryPath)
    call mem_setptr(this%kchangeper, 'KCHANGEPER', npfMemoryPath)
    call mem_setptr(this%kchangestp, 'KCHANGESTP', npfMemoryPath)
    call mem_setptr(this%nodekchange, 'NODEKCHANGE', npfMemoryPath)
    !
    ! -- set input mempath pointers
    call mem_setptr(this%k11_src, 'K', this%input_mempath)
    call mem_setptr(this%k22_src, 'K22', this%input_mempath)
    call mem_setptr(this%k33_src, 'K33', this%input_mempath)
  end subroutine tvk_ar_set_pointers

  !> @brief Apply input K/K22/K33 column changes for period-data row n to node.
  !<
  subroutine tvk_apply_row_changes(this, n, node)
    ! -- dummy
    class(TvkType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: node
    ! -- local
    character(len=LINELENGTH) :: cellstr
    ! -- formats
    character(len=*), parameter :: fmtvalchg = &
      "(a, ' package: Setting ', a, ' value for cell ', a, ' at start of &
      &stress period ', i0, ' = ', g12.5)"
    !
    ! -- K is processed before K22/K33 so that validate_change can use
    ! -- the already-updated k11 value when ik22overk/ik33overk are set.
    if (this%k11_src(n) /= DNODATA) then
      this%k11(node) = this%k11_src(n)
      call this%validate_change(node, 'K')
      if (this%iprpak /= 0) then
        call this%dis%noder_to_string(node, cellstr)
        write (this%iout, fmtvalchg) &
          trim(adjustl(this%packName)), 'K', trim(cellstr), kper, this%k11(node)
      end if
    end if
    !
    if (this%k22_src(n) /= DNODATA) then
      this%k22(node) = this%k22_src(n)
      call this%validate_change(node, 'K22')
      if (this%iprpak /= 0) then
        call this%dis%noder_to_string(node, cellstr)
        write (this%iout, fmtvalchg) &
          trim(adjustl(this%packName)), 'K22', trim(cellstr), kper, &
          this%k22(node)
      end if
    end if
    !
    if (this%k33_src(n) /= DNODATA) then
      this%k33(node) = this%k33_src(n)
      call this%validate_change(node, 'K33')
      if (this%iprpak /= 0) then
        call this%dis%noder_to_string(node, cellstr)
        write (this%iout, fmtvalchg) &
          trim(adjustl(this%packName)), 'K33', trim(cellstr), kper, &
          this%k33(node)
      end if
    end if
  end subroutine tvk_apply_row_changes

  !> @brief Mark property changes as having occurred at (kper, kstp)
  !!
  !! Deferred procedure implementation called by the TvBaseType code when a
  !! property value change occurs at (kper, kstp).
  !<
  subroutine tvk_set_changed_at(this, kper, kstp)
    ! -- dummy
    class(TvkType) :: this
    integer(I4B), intent(in) :: kper
    integer(I4B), intent(in) :: kstp
    !
    this%kchangeper = kper
    this%kchangestp = kstp
  end subroutine tvk_set_changed_at

  !> @brief Clear all per-node change flags
  !!
  !! Deferred procedure implementation called by the TvBaseType code when a
  !! new time step commences, indicating that any previously set per-node
  !! property value change flags should be reset.
  !<
  subroutine tvk_reset_change_flags(this)
    ! -- dummy variables
    class(TvkType) :: this
    ! -- local variables
    integer(I4B) :: i
    !
    do i = 1, this%dis%nodes
      this%nodekchange(i) = 0
    end do
  end subroutine tvk_reset_change_flags

  !> @brief Check that a given property value is valid
  !!
  !! Deferred procedure implementation called by the TvBaseType code after a
  !! property value change occurs. Check if the property value of the given
  !! variable at the given node is invalid, and log an error if so. Update
  !! K22 and K33 values appropriately when specified as anisotropy.
  !<
  subroutine tvk_validate_change(this, n, varName)
    ! -- dummy
    class(TvkType) :: this
    integer(I4B), intent(in) :: n
    character(len=*), intent(in) :: varName
    ! -- local
    character(len=LINELENGTH) :: cellstr
    ! -- formats
    character(len=*), parameter :: fmtkerr = &
      "(1x, a, ' changed hydraulic property ', a, ' is <= 0 for cell ', a, &
      &' ', 1pg15.6)"
    !
    ! -- Mark the node as being changed this time step
    this%nodekchange(n) = 1
    !
    ! -- Check the changed value is ok
    if (varName == 'K') then
      if (this%k11(n) <= DZERO) then
        call this%dis%noder_to_string(n, cellstr)
        write (errmsg, fmtkerr) &
          trim(adjustl(this%packName)), 'K', trim(cellstr), this%k11(n)
        call store_error(errmsg)
      end if
    elseif (varName == 'K22') then
      if (this%ik22overk == 1) then
        this%k22(n) = this%k22(n) * this%k11(n)
      end if
      if (this%k22(n) <= DZERO) then
        call this%dis%noder_to_string(n, cellstr)
        write (errmsg, fmtkerr) &
          trim(adjustl(this%packName)), 'K22', trim(cellstr), this%k22(n)
        call store_error(errmsg)
      end if
    elseif (varName == 'K33') then
      if (this%ik33overk == 1) then
        this%k33(n) = this%k33(n) * this%k11(n)
      end if
      if (this%k33(n) <= DZERO) then
        call this%dis%noder_to_string(n, cellstr)
        write (errmsg, fmtkerr) &
          trim(adjustl(this%packName)), 'K33', trim(cellstr), this%k33(n)
        call store_error(errmsg)
      end if
    end if
  end subroutine tvk_validate_change

  !> @brief Deallocate package memory
  !!
  !! Deallocate TVK package scalars and arrays.
  !<
  subroutine tvk_da(this)
    ! -- dummy
    class(TvkType) :: this
    !
    nullify (this%ik22overk)
    nullify (this%ik33overk)
    nullify (this%k11)
    nullify (this%k22)
    nullify (this%k33)
    nullify (this%kchangeper)
    nullify (this%kchangestp)
    nullify (this%nodekchange)
    nullify (this%k11_src)
    nullify (this%k22_src)
    nullify (this%k33_src)
    call tvbase_da(this)
  end subroutine tvk_da

end module TvkModule
