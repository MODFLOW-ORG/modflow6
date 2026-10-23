module TestMemoryAlias
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LENMEMADDRESS
  use MemoryManagerModule, only: mem_allocate, mem_setptr, &
                                 get_from_memorystore
  use MemoryTypeModule, only: MemoryType
  use CharacterStringModule, only: CharacterStringType
  use MemoryAliasModule, only: mem_alias_deprecated, get_deprecated_alias, &
                               mem_alias_da
  use testdrive, only: error_type, unittest_type, new_unittest, check

  implicit none
  private
  public :: collect_memoryalias

  character(len=*), parameter :: MEMPATH = 'TESTALIAS'
  character(len=*), parameter :: VERSION = '6.9.0'

contains

  subroutine collect_memoryalias(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("alias_arrays", test_alias_arrays), &
                new_unittest("alias_scalar", test_alias_scalar), &
                new_unittest("alias_unassociated", test_alias_unassociated), &
                new_unittest("alias_other_path", test_alias_other_path), &
                new_unittest("alias_da", test_alias_da) &
                ]
  end subroutine collect_memoryalias

  !> @brief Check that a name is a deprecated alias of a replacement
  !<
  subroutine check_recorded(error, name, mem_path, replacement)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: mem_path
    character(len=*), intent(in) :: replacement
    logical(LGP) :: found
    character(len=LENMEMADDRESS) :: address
    character(len=16) :: ver

    call get_deprecated_alias(name, mem_path, found, address, ver)
    call check(error, found, name//" not recorded as deprecated")
    if (allocated(error)) return
    call check(error, address == replacement, &
               name//" replacement '"//trim(address)//"'")
    if (allocated(error)) return
    call check(error, ver == VERSION, name//" version '"//trim(ver)//"'")
  end subroutine check_recorded

  !> @brief Alias each array rank and type and share its memory
  !<
  subroutine test_alias_arrays(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), dimension(:), pointer, contiguous :: ia1, ib1
    integer(I4B), dimension(:, :), pointer, contiguous :: ia2, ib2
    real(DP), dimension(:), pointer, contiguous :: da1, db1
    real(DP), dimension(:, :), pointer, contiguous :: da2, db2
    logical(LGP), dimension(:), pointer, contiguous :: la1, lb1
    type(CharacterStringType), dimension(:), pointer, contiguous :: ca1, cb1
    integer(I4B) :: n

    call mem_allocate(la1, 3, 'NEWLOG1D', MEMPATH)
    call mem_allocate(ia1, 3, 'NEWINT1D', MEMPATH)
    call mem_allocate(ia2, 2, 3, 'NEWINT2D', MEMPATH)
    call mem_allocate(da1, 3, 'NEWDBL1D', MEMPATH)
    call mem_allocate(da2, 2, 3, 'NEWDBL2D', MEMPATH)
    call mem_allocate(ca1, 8, 3, 'NEWSTR1D', MEMPATH)
    do n = 1, 3
      la1(n) = .true.
      ia1(n) = n
      da1(n) = real(n, DP)
    end do

    call mem_alias_deprecated(la1, 'OLDLOG1D', MEMPATH, 'NEWLOG1D', MEMPATH, &
                              VERSION)
    call mem_alias_deprecated(ia1, 'OLDINT1D', MEMPATH, 'NEWINT1D', MEMPATH, &
                              VERSION)
    call mem_alias_deprecated(ia2, 'OLDINT2D', MEMPATH, 'NEWINT2D', MEMPATH, &
                              VERSION)
    call mem_alias_deprecated(da1, 'OLDDBL1D', MEMPATH, 'NEWDBL1D', MEMPATH, &
                              VERSION)
    call mem_alias_deprecated(da2, 'OLDDBL2D', MEMPATH, 'NEWDBL2D', MEMPATH, &
                              VERSION)
    call mem_alias_deprecated(ca1, 8, 'OLDSTR1D', MEMPATH, 'NEWSTR1D', &
                              MEMPATH, VERSION)

    ! -- the previous names resolve to the same memory
    call mem_setptr(lb1, 'OLDLOG1D', MEMPATH)
    call mem_setptr(ib1, 'OLDINT1D', MEMPATH)
    call mem_setptr(ib2, 'OLDINT2D', MEMPATH)
    call mem_setptr(db1, 'OLDDBL1D', MEMPATH)
    call mem_setptr(db2, 'OLDDBL2D', MEMPATH)
    call mem_setptr(cb1, 'OLDSTR1D', MEMPATH)
    call check(error, associated(lb1, la1), "OLDLOG1D shares memory")
    if (allocated(error)) return
    call check(error, associated(ib1, ia1), "OLDINT1D shares memory")
    if (allocated(error)) return
    call check(error, associated(ib2, ia2), "OLDINT2D shares memory")
    if (allocated(error)) return
    call check(error, associated(db1, da1), "OLDDBL1D shares memory")
    if (allocated(error)) return
    call check(error, associated(db2, da2), "OLDDBL2D shares memory")
    if (allocated(error)) return
    call check(error, associated(cb1, ca1), "OLDSTR1D shares memory")
    if (allocated(error)) return

    ! -- each previous name is recorded with its replacement
    call check_recorded(error, 'OLDINT1D', MEMPATH, MEMPATH//'/NEWINT1D')
    if (allocated(error)) return
    call check_recorded(error, 'OLDDBL2D', MEMPATH, MEMPATH//'/NEWDBL2D')
    if (allocated(error)) return
    call check_recorded(error, 'OLDLOG1D', MEMPATH, MEMPATH//'/NEWLOG1D')
    if (allocated(error)) return
    call check_recorded(error, 'OLDSTR1D', MEMPATH, MEMPATH//'/NEWSTR1D')
    if (allocated(error)) return

    ! -- the replacement itself is not deprecated
    call check_not_recorded(error, 'NEWINT1D', MEMPATH)
  end subroutine test_alias_arrays

  !> @brief Alias an integer scalar and share its memory
  !<
  subroutine test_alias_scalar(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), pointer :: inew, iold

    call mem_allocate(inew, 'NEWSCALAR', MEMPATH)
    inew = 7
    call mem_alias_deprecated(inew, 'OLDSCALAR', MEMPATH, 'NEWSCALAR', &
                              MEMPATH, VERSION)

    call mem_setptr(iold, 'OLDSCALAR', MEMPATH)
    call check(error, associated(iold, inew), "OLDSCALAR shares memory")
    if (allocated(error)) return
    iold = 11
    call check(error, inew == 11, "write through OLDSCALAR not seen")
    if (allocated(error)) return
    call check_recorded(error, 'OLDSCALAR', MEMPATH, MEMPATH//'/NEWSCALAR')
  end subroutine test_alias_scalar

  !> @brief Skip an alias when the replacement does not exist
  !<
  subroutine test_alias_unassociated(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(:), pointer, contiguous :: absent => null()
    type(MemoryType), pointer :: mt
    logical(LGP) :: found

    call mem_alias_deprecated(absent, 'OLDABSENT', MEMPATH, 'NEWABSENT', &
                              MEMPATH, VERSION)
    call get_from_memorystore('OLDABSENT', MEMPATH, mt, found, check=.false.)
    call check(error,.not. found, "OLDABSENT checked in")
    if (allocated(error)) return
    call check_not_recorded(error, 'OLDABSENT', MEMPATH)
  end subroutine test_alias_unassociated

  !> @brief Alias a variable stored at a different memory path
  !<
  subroutine test_alias_other_path(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(:), pointer, contiguous :: anew, aold
    character(len=*), parameter :: input_path = '__INPUT__/TESTALIAS'

    call mem_allocate(anew, 2, 'MOVED', input_path)
    call mem_alias_deprecated(anew, 'MOVED_OLD', MEMPATH, 'MOVED', &
                              input_path, VERSION)
    call mem_setptr(aold, 'MOVED_OLD', MEMPATH)
    call check(error, associated(aold, anew), "MOVED_OLD shares memory")
    if (allocated(error)) return
    call check_recorded(error, 'MOVED_OLD', MEMPATH, input_path//'/MOVED')
  end subroutine test_alias_other_path

  !> @brief Clear the recorded aliases
  !<
  subroutine test_alias_da(error)
    type(error_type), allocatable, intent(out) :: error
    real(DP), dimension(:), pointer, contiguous :: anew

    call mem_allocate(anew, 2, 'CLEARED', MEMPATH)
    call mem_alias_deprecated(anew, 'CLEARED_OLD', MEMPATH, 'CLEARED', &
                              MEMPATH, VERSION)
    call check_recorded(error, 'CLEARED_OLD', MEMPATH, MEMPATH//'/CLEARED')
    if (allocated(error)) return
    call mem_alias_da()
    call check_not_recorded(error, 'CLEARED_OLD', MEMPATH)
  end subroutine test_alias_da

  !> @brief Check that a name is not recorded as deprecated
  !<
  subroutine check_not_recorded(error, name, mem_path)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: mem_path
    logical(LGP) :: found
    character(len=LENMEMADDRESS) :: address
    character(len=16) :: ver

    call get_deprecated_alias(name, mem_path, found, address, ver)
    call check(error,.not. found, name//" recorded as deprecated")
  end subroutine check_not_recorded

end module TestMemoryAlias
