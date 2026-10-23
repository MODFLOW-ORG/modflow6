!> @brief Previous variable names checked in as deprecated aliases
!!
!! A package that renames or moves a variable checks in the previous name as
!! an alias of the variable, so a program that addresses it through the API
!! still resolves, and records the replacement for the API to report.
!<
module MemoryAliasModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENMEMADDRESS
  use MemoryHelperModule, only: create_mem_address
  use MemoryManagerModule, only: mem_allocate, mem_checkin, mem_reassignptr
  use HashTableModule, only: HashTableType, hash_table_cr, hash_table_da
  use CharacterStringModule, only: CharacterStringType
  use ArrayHandlersModule, only: ExpandArray

  implicit none
  private
  public :: mem_alias_deprecated
  public :: get_deprecated_alias
  public :: mem_alias_da

  integer(I4B), parameter :: LENVERSION = 16 !< maximum length of a version string

  interface mem_alias_deprecated
    module procedure &
      alias_int, &
      alias_logical1d, &
      alias_int1d, &
      alias_int2d, &
      alias_dbl1d, &
      alias_dbl2d, &
      alias_charstr1d
  end interface mem_alias_deprecated

  type(HashTableType), pointer :: alias_index => null() !< alias address to record number
  integer(I4B) :: nalias = 0 !< number of recorded aliases
  character(len=LENMEMADDRESS), dimension(:), allocatable :: replacements !< address of the replacement variable
  character(len=LENVERSION), dimension(:), allocatable :: versions !< version the alias was deprecated in

contains

  !> @brief Check in a deprecated name for an integer scalar
  !<
  subroutine alias_int(sclr, name, mem_path, name2, mem_path2, version)
    integer(I4B), pointer, intent(in) :: sclr !< the replacement variable
    character(len=*), intent(in) :: name !< deprecated variable name
    character(len=*), intent(in) :: mem_path !< path of the deprecated name
    character(len=*), intent(in) :: name2 !< replacement variable name
    character(len=*), intent(in) :: mem_path2 !< path of the replacement variable
    character(len=*), intent(in) :: version !< version the name was deprecated in
    ! -- local
    integer(I4B), pointer :: alias => null()
    !
    if (.not. associated(sclr)) return
    !
    ! -- a scalar has no checkin, so point a new entry at the replacement
    call mem_allocate(alias, name, mem_path)
    call mem_reassignptr(alias, name, mem_path, name2, mem_path2)
    call record_alias(name, mem_path, name2, mem_path2, version)
  end subroutine alias_int

  !> @brief Check in a deprecated name for a 1d logical array
  !<
  subroutine alias_logical1d(alog, name, mem_path, name2, mem_path2, version)
    logical(LGP), dimension(:), pointer, contiguous, intent(in) :: alog !< the replacement variable
    character(len=*), intent(in) :: name !< deprecated variable name
    character(len=*), intent(in) :: mem_path !< path of the deprecated name
    character(len=*), intent(in) :: name2 !< replacement variable name
    character(len=*), intent(in) :: mem_path2 !< path of the replacement variable
    character(len=*), intent(in) :: version !< version the name was deprecated in
    !
    if (.not. associated(alog)) return
    call mem_checkin(alog, name, mem_path, name2, mem_path2)
    call record_alias(name, mem_path, name2, mem_path2, version)
  end subroutine alias_logical1d

  !> @brief Check in a deprecated name for a 1d integer array
  !<
  subroutine alias_int1d(aint, name, mem_path, name2, mem_path2, version)
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: aint !< the replacement variable
    character(len=*), intent(in) :: name !< deprecated variable name
    character(len=*), intent(in) :: mem_path !< path of the deprecated name
    character(len=*), intent(in) :: name2 !< replacement variable name
    character(len=*), intent(in) :: mem_path2 !< path of the replacement variable
    character(len=*), intent(in) :: version !< version the name was deprecated in
    !
    if (.not. associated(aint)) return
    call mem_checkin(aint, name, mem_path, name2, mem_path2)
    call record_alias(name, mem_path, name2, mem_path2, version)
  end subroutine alias_int1d

  !> @brief Check in a deprecated name for a 2d integer array
  !<
  subroutine alias_int2d(aint, name, mem_path, name2, mem_path2, version)
    integer(I4B), dimension(:, :), pointer, contiguous, intent(inout) :: aint !< the replacement variable
    character(len=*), intent(in) :: name !< deprecated variable name
    character(len=*), intent(in) :: mem_path !< path of the deprecated name
    character(len=*), intent(in) :: name2 !< replacement variable name
    character(len=*), intent(in) :: mem_path2 !< path of the replacement variable
    character(len=*), intent(in) :: version !< version the name was deprecated in
    !
    if (.not. associated(aint)) return
    call mem_checkin(aint, name, mem_path, name2, mem_path2)
    call record_alias(name, mem_path, name2, mem_path2, version)
  end subroutine alias_int2d

  !> @brief Check in a deprecated name for a 1d double precision array
  !<
  subroutine alias_dbl1d(adbl, name, mem_path, name2, mem_path2, version)
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: adbl !< the replacement variable
    character(len=*), intent(in) :: name !< deprecated variable name
    character(len=*), intent(in) :: mem_path !< path of the deprecated name
    character(len=*), intent(in) :: name2 !< replacement variable name
    character(len=*), intent(in) :: mem_path2 !< path of the replacement variable
    character(len=*), intent(in) :: version !< version the name was deprecated in
    !
    if (.not. associated(adbl)) return
    call mem_checkin(adbl, name, mem_path, name2, mem_path2)
    call record_alias(name, mem_path, name2, mem_path2, version)
  end subroutine alias_dbl1d

  !> @brief Check in a deprecated name for a 2d double precision array
  !<
  subroutine alias_dbl2d(adbl, name, mem_path, name2, mem_path2, version)
    real(DP), dimension(:, :), pointer, contiguous, intent(inout) :: adbl !< the replacement variable
    character(len=*), intent(in) :: name !< deprecated variable name
    character(len=*), intent(in) :: mem_path !< path of the deprecated name
    character(len=*), intent(in) :: name2 !< replacement variable name
    character(len=*), intent(in) :: mem_path2 !< path of the replacement variable
    character(len=*), intent(in) :: version !< version the name was deprecated in
    !
    if (.not. associated(adbl)) return
    call mem_checkin(adbl, name, mem_path, name2, mem_path2)
    call record_alias(name, mem_path, name2, mem_path2, version)
  end subroutine alias_dbl2d

  !> @brief Check in a deprecated name for a 1d character string array
  !<
  subroutine alias_charstr1d(acharstr1d, ilen, name, mem_path, name2, &
                             mem_path2, version)
    type(CharacterStringType), dimension(:), &
      pointer, contiguous, intent(inout) :: acharstr1d !< the replacement variable
    integer(I4B), intent(in) :: ilen !< string length
    character(len=*), intent(in) :: name !< deprecated variable name
    character(len=*), intent(in) :: mem_path !< path of the deprecated name
    character(len=*), intent(in) :: name2 !< replacement variable name
    character(len=*), intent(in) :: mem_path2 !< path of the replacement variable
    character(len=*), intent(in) :: version !< version the name was deprecated in
    !
    if (.not. associated(acharstr1d)) return
    call mem_checkin(acharstr1d, ilen, name, mem_path, name2, mem_path2)
    call record_alias(name, mem_path, name2, mem_path2, version)
  end subroutine alias_charstr1d

  !> @brief Record the replacement and version of a deprecated name
  !<
  subroutine record_alias(name, mem_path, name2, mem_path2, version)
    character(len=*), intent(in) :: name !< deprecated variable name
    character(len=*), intent(in) :: mem_path !< path of the deprecated name
    character(len=*), intent(in) :: name2 !< replacement variable name
    character(len=*), intent(in) :: mem_path2 !< path of the replacement variable
    character(len=*), intent(in) :: version !< version the name was deprecated in
    !
    if (.not. associated(alias_index)) call hash_table_cr(alias_index)
    nalias = nalias + 1
    call ExpandArray(replacements)
    call ExpandArray(versions)
    replacements(nalias) = create_mem_address(mem_path2, name2)
    versions(nalias) = version
    call alias_index%add(create_mem_address(mem_path, name), nalias)
  end subroutine record_alias

  !> @brief Get the replacement and version of a deprecated name
  !<
  subroutine get_deprecated_alias(name, mem_path, found, replacement, version)
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path of the variable
    logical(LGP), intent(out) :: found !< true when the name is deprecated
    character(len=LENMEMADDRESS), intent(out) :: replacement !< address of the replacement variable
    character(len=*), intent(out) :: version !< version the name was deprecated in
    ! -- local
    integer(I4B) :: i
    !
    found = .false.
    replacement = ''
    version = ''
    if (.not. associated(alias_index)) return
    i = alias_index%get(create_mem_address(mem_path, name))
    if (i == 0) return
    found = .true.
    replacement = replacements(i)
    version = versions(i)
  end subroutine get_deprecated_alias

  !> @brief Clear the recorded aliases
  !!
  !! The alias entries themselves are deallocated with the memory manager.
  !<
  subroutine mem_alias_da()
    if (associated(alias_index)) call hash_table_da(alias_index)
    if (allocated(replacements)) deallocate (replacements)
    if (allocated(versions)) deallocate (versions)
    nalias = 0
  end subroutine mem_alias_da

end module MemoryAliasModule
