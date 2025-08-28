module MemoryTypeModule

  use KindModule, only: DP, LGP, I4B
  use ConstantsModule, only: LENMEMPATH, LENMEMADDRESS, LENTIMESERIESNAME, &
                             LENVARNAME, MAXMEMRANK, LENMEMTYPE, &
                             TABSTRING, TABINTEGER, &
                             TABCENTER, TABLEFT, TABRIGHT
  use CharacterStringModule, only: CharacterStringType
  use TableModule, only: TableType
  use MemoryHelperModule, only: create_mem_address

  implicit none
  private
  public :: MemoryType

  type MemoryType
    character(len=LENVARNAME) :: name !< name of the array
    character(len=LENVARNAME) :: mastername = 'none' !< name of the master array
    character(len=LENMEMPATH) :: path !< path to memory object
    character(len=LENMEMPATH) :: masterPath = 'none' !< path to master memory object
    character(len=LENMEMTYPE) :: memtype !< type (LOGICAL or INTEGER or DOUBLE)
    integer(I4B) :: id !< id, not used
    integer(I4B) :: nrealloc = 0 !< number of times reallocated
    integer(I4B) :: isize = -1 !< size of the array, equal to the number of array elements; 1 for scalars
    integer(I4B) :: element_size = 0 !< byte size of an element; string length
    integer(I4B) :: set_handler_idx = 0 !< index of side effect handler for external access
    logical(LGP) :: master = .true. !< master copy, others point to this one
    character(len=:), pointer :: strsclr => null() !< pointer to the character string
    logical(LGP), pointer :: logicalsclr => null() !< pointer to the logical
    integer(I4B), pointer :: intsclr => null() !< pointer to the integer
    real(DP), pointer :: dblsclr => null() !< pointer to the double
    ! The 1d character string array is handled differently than the other arrays due to a bug in gfortran 11.3 and 12.1.
    ! Due to this bug the length of the string is not stored in the array descriptor. With a segmentation fault as a result
    ! on deallocation.
    class(*), dimension(:), pointer, contiguous :: astr1d => null() !< pointer to the 1d character string array
    logical(LGP), dimension(:), pointer, contiguous :: alogical1d => null() !< pointer to 1d logical array
    integer(I4B), dimension(:), pointer, contiguous :: aint1d => null() !< pointer to 1d integer array
    integer(I4B), dimension(:, :), pointer, contiguous :: aint2d => null() !< pointer to 2d integer array
    integer(I4B), dimension(:, :, :), pointer, contiguous :: aint3d => null() !< pointer to 3d integer array
    real(DP), dimension(:), pointer, contiguous :: adbl1d => null() !< pointer to 1d double array
    real(DP), dimension(:, :), pointer, contiguous :: adbl2d => null() !< pointer to 2d double array
    real(DP), dimension(:, :, :), pointer, contiguous :: adbl3d => null() !< pointer to 3d double array
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      acharstr1d => null() !< pointer to the 1d character string array
  contains
    procedure :: table_entry
    procedure :: mt_associated
    procedure :: mt_deallocate
  end type

contains

  subroutine table_entry(this, memtab)
    ! -- dummy
    class(MemoryType) :: this
    type(TableType), intent(inout) :: memtab
    ! -- local
    character(len=16) :: cmem
    character(len=LENMEMADDRESS) :: cptr
    integer(I4B) :: ipos
    ! -- formats
    !
    ! -- determine memory type
    ipos = index(this%memtype, ' (')
    if (ipos < 1) then
      ipos = 16
    else
      ipos = min(16, ipos - 1)
    end if
    cmem = this%memtype(1:ipos)
    !
    ! -- Set pointer string
    cptr = '--'
    if (.not. this%master) then
      cptr = create_mem_address(this%masterPath, this%mastername)
    end if
    !
    ! -- write data to the table
    call memtab%add_term(this%path)
    call memtab%add_term(this%name)
    call memtab%add_term(cmem)
    call memtab%add_term(this%isize)
    call memtab%add_term(cptr)
  end subroutine table_entry

  function mt_associated(this) result(al)
    class(MemoryType) :: this
    logical :: al
    al = .false.
    if (associated(this%strsclr)) al = .true.
    if (associated(this%logicalsclr)) al = .true.
    if (associated(this%intsclr)) al = .true.
    if (associated(this%dblsclr)) al = .true.
    if (associated(this%alogical1d)) al = .true.
    if (associated(this%astr1d)) al = .true.
    if (associated(this%aint1d)) al = .true.
    if (associated(this%aint2d)) al = .true.
    if (associated(this%aint3d)) al = .true.
    if (associated(this%adbl1d)) al = .true.
    if (associated(this%adbl2d)) al = .true.
    if (associated(this%adbl3d)) al = .true.
    if (associated(this%acharstr1d)) al = .true.
  end function mt_associated

  subroutine mt_deallocate(this)
    use iso_c_binding, only: c_loc, c_ptr, c_null_ptr, c_f_pointer
    class(MemoryType) :: this
    integer(I4B) :: n
    type(c_ptr) :: cptr

    character(len=1), dimension(:), pointer :: astr1d

    if (associated(this%strsclr)) then
      if (this%master) deallocate (this%strsclr)
      nullify (this%strsclr)
    end if

    if (associated(this%logicalsclr)) then
      if (this%master) deallocate (this%logicalsclr)
      nullify (this%logicalsclr)
    end if

    if (associated(this%intsclr)) then
      if (this%master) deallocate (this%intsclr)
      nullify (this%intsclr)
    end if

    if (associated(this%dblsclr)) then
      if (this%master) deallocate (this%dblsclr)
      nullify (this%dblsclr)
    end if

    ! Handle the dealloction of the 1d character string array differently due to a bug in gfortran 11.3, 12.1, 13.1 and 13.2.
    ! Due to a bug in the gfortran compiler we can't use a deferred length character variable
    ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=106317
    !
    ! We use a c_ptr to cast the pointer to a string array with a length of 1. The actual length of the array is
    ! computed by the actual length of the string multiplied by the array size.
    ! So we go from the actual character(len=element_size), dimension(isize) to a character(len=1), dimension(isize*element_size).

    if (associated(this%astr1d)) then
      select type (item => this%astr1d)
      type is (character(*))
        cptr = c_loc(item)
      class default
        cptr = c_null_ptr
      end select

      call c_f_pointer(cptr, astr1d, [this%isize * this%element_size])

#if __GFORTRAN__  && ((__GNUC__ < 13) || (__GNUC__ == 13 && __GNUC_MINOR__ < 3))
      if (this%master) deallocate (astr1d)
#else
      if (this%master) deallocate (this%astr1d)
#endif

      nullify (this%astr1d)
    end if

    if (associated(this%alogical1d)) then
      if (this%master) deallocate (this%alogical1d)
      nullify (this%alogical1d)
    end if

    if (associated(this%aint1d)) then
      if (this%master) deallocate (this%aint1d)
      nullify (this%aint1d)
    end if

    if (associated(this%aint2d)) then
      if (this%master) deallocate (this%aint2d)
      nullify (this%aint2d)
    end if

    if (associated(this%aint3d)) then
      if (this%master) deallocate (this%aint3d)
      nullify (this%aint3d)
    end if

    if (associated(this%adbl1d)) then
      if (this%master) deallocate (this%adbl1d)
      nullify (this%adbl1d)
    end if

    if (associated(this%adbl2d)) then
      if (this%master) deallocate (this%adbl2d)
      nullify (this%adbl2d)
    end if

    if (associated(this%adbl3d)) then
      if (this%master) deallocate (this%adbl3d)
      nullify (this%adbl3d)
    end if

    if (associated(this%acharstr1d)) then
      if (this%master) then
        do n = 1, size(this%acharstr1d)
          call this%acharstr1d(n)%destroy()
        end do
        deallocate (this%acharstr1d)
      end if
      nullify (this%acharstr1d)
    end if
  end subroutine mt_deallocate

end module MemoryTypeModule
