module MemoryManagerModule

  use KindModule, only: DP, LGP, I4B, I8B
  use ConstantsModule, only: DZERO, DONE, &
                             DEM3, DEM6, DEM9, DEP3, DEP6, DEP9, &
                             LENMEMPATH, LENMEMSEPARATOR, LENVARNAME, &
                             LENMEMADDRESS, LENCOMPONENTNAME, &
                             LENMEMTYPE, LINELENGTH, &
                             TABSTRING, TABUCSTRING, &
                             TABINTEGER, TABREAL, TABCENTER, TABLEFT, &
                             TABRIGHT
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, count_errors
  use MemoryTypeModule, only: MemoryType
  use MemoryStoreModule, only: MemoryStoreType
  use MemoryContainerIteratorModule, only: MemoryContainerIteratorType
  use MemoryHelperModule, only: mem_check_length, split_mem_path, &
                                strip_context_mem_path, get_mem_path_context
  use TableModule, only: TableType, table_cr
  use CharacterStringModule, only: CharacterStringType

  implicit none
  private
  public :: mem_allocate
  public :: mem_checkin
  public :: mem_reallocate
  public :: mem_setptr
  public :: mem_copyptr
  public :: mem_reassignptr
  public :: mem_deallocate
  public :: mem_write_usage
  public :: mem_da
  public :: mem_set_print_option
  public :: get_from_memorystore

  public :: get_mem_type
  public :: get_mem_rank
  public :: get_mem_elem_size
  public :: get_mem_shape
  public :: get_isize
  public :: copy_dbl1d

  public :: memorystore
  public :: mem_print_detailed

  type(MemoryStoreType) :: memorystore
  type(TableType), pointer :: memtab => null()
  integer(I8B) :: nvalues_alogical = 0
  integer(I8B) :: nvalues_astr = 0
  integer(I8B) :: nvalues_aint = 0
  integer(I8B) :: nvalues_adbl = 0
  integer(I4B) :: iprmem = 0

  interface mem_allocate
    module procedure &
      allocate_logical, &
      allocate_logical1d, &
      allocate_str, &
      allocate_str1d, &
      allocate_int, &
      allocate_int1d, &
      allocate_int2d, &
      allocate_int3d, &
      allocate_dbl, &
      allocate_dbl1d, &
      allocate_dbl2d, &
      allocate_dbl3d, &
      allocate_charstr1d
  end interface mem_allocate

  interface mem_checkin
    module procedure &
      checkin_logical1d, &
      checkin_int1d, &
      checkin_int2d, &
      checkin_dbl1d, &
      checkin_dbl2d, &
      checkin_charstr1d
  end interface mem_checkin

  interface mem_reallocate
    module procedure &
      reallocate_logical1d, &
      reallocate_int1d, &
      reallocate_int2d, &
      reallocate_dbl1d, &
      reallocate_dbl2d, &
      reallocate_str1d, &
      reallocate_charstr1d
  end interface mem_reallocate

  interface mem_setptr
    module procedure &
      setptr_logical, &
      setptr_logical1d, &
      setptr_int, &
      setptr_int1d, &
      setptr_int2d, &
      setptr_int3d, &
      setptr_dbl, &
      setptr_dbl1d, &
      setptr_dbl2d, &
      setptr_dbl3d, &
      setptr_str, &
      setptr_str1d, &
      setptr_charstr1d
  end interface mem_setptr

  interface mem_copyptr
    module procedure &
      copyptr_logical1d, &
      copyptr_int1d, &
      copyptr_int2d, &
      copyptr_dbl1d, &
      copyptr_dbl2d
  end interface mem_copyptr

  interface mem_reassignptr
    module procedure &
      reassignptr_int, &
      reassignptr_logical1d, &
      reassignptr_int1d, &
      reassignptr_int2d, &
      reassignptr_dbl1d, &
      reassignptr_dbl2d
  end interface mem_reassignptr

  interface mem_deallocate
    module procedure &
      deallocate_logical, &
      deallocate_logical1d, &
      deallocate_str, &
      deallocate_str1d, &
      deallocate_charstr1d, &
      deallocate_int, &
      deallocate_int1d, &
      deallocate_int2d, &
      deallocate_int3d, &
      deallocate_dbl, &
      deallocate_dbl1d, &
      deallocate_dbl2d, &
      deallocate_dbl3d
  end interface mem_deallocate

contains

  !> @ brief Get the variable memory type
  !!
  !! Returns any of 'LOGICAL', 'INTEGER', 'DOUBLE', 'STRING'.
  !! returns 'UNKNOWN' when the variable is not found.
  !<
  subroutine get_mem_type(name, mem_path, var_type)
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where the variable is stored
    character(len=LENMEMTYPE), intent(out) :: var_type !< memory type
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    mt => null()
    var_type = 'UNKNOWN'
    call get_from_memorystore(name, mem_path, mt, found)
    if (found) then
      var_type = mt%memtype
    end if
  end subroutine get_mem_type

  !> @ brief Get the variable rank
  !!
  !! Returns rank = -1 when not found.
  !<
  subroutine get_mem_rank(name, mem_path, rank)
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< mem_path
    integer(I4B), intent(out) :: rank !< rank
    ! -- local
    type(MemoryType), pointer :: mt => null()
    logical(LGP) :: found
    ! -- code
    !
    ! -- initialize rank to a value to communicate failure
    rank = -1
    !
    ! -- get the entry from the memory manager
    call get_from_memorystore(name, mem_path, mt, found)
    !
    ! -- set rank
    if (found) then
      if (associated(mt%logicalsclr)) rank = 0
      if (associated(mt%intsclr)) rank = 0
      if (associated(mt%dblsclr)) rank = 0
      if (associated(mt%alogical1d)) rank = 1
      if (associated(mt%aint1d)) rank = 1
      if (associated(mt%aint2d)) rank = 2
      if (associated(mt%aint3d)) rank = 3
      if (associated(mt%adbl1d)) rank = 1
      if (associated(mt%adbl2d)) rank = 2
      if (associated(mt%adbl3d)) rank = 3
      if (associated(mt%strsclr)) rank = 0
      if (associated(mt%astr1d)) rank = 1
      if (associated(mt%acharstr1d)) rank = 1
    end if
  end subroutine get_mem_rank

  !> @ brief Get the memory size of a single element of the stored variable
  !!
  !! Memory size in bytes, returns size = -1 when not found. This is
  !< also string length.
  subroutine get_mem_elem_size(name, mem_path, size)
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where the variable is stored
    integer(I4B), intent(out) :: size !< size of the variable in bytes
    ! -- local
    type(MemoryType), pointer :: mt => null()
    logical(LGP) :: found
    ! -- code
    !
    ! -- initialize size to a value to communicate failure
    size = -1
    !
    ! -- get the entry from the memory manager
    call get_from_memorystore(name, mem_path, mt, found)
    !
    ! -- set memory size
    if (found) then
      size = mt%element_size
    end if
  end subroutine get_mem_elem_size

  !> @ brief Get the variable memory shape
  !!
  !! Returns an integer array with the shape (Fortran ordering),
  !! and set shape(1) = -1 when not found.
  !<
  subroutine get_mem_shape(name, mem_path, mem_shape)
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where the variable is stored
    integer(I4B), dimension(:), intent(out) :: mem_shape !< shape of the variable
    ! -- local
    type(MemoryType), pointer :: mt => null()
    logical(LGP) :: found
    ! -- code
    !
    ! -- get the entry from the memory manager
    call get_from_memorystore(name, mem_path, mt, found)
    !
    ! -- set shape
    if (found) then
      if (associated(mt%logicalsclr)) mem_shape = shape(mt%logicalsclr)
      if (associated(mt%intsclr)) mem_shape = shape(mt%logicalsclr)
      if (associated(mt%dblsclr)) mem_shape = shape(mt%dblsclr)
      if (associated(mt%alogical1d)) mem_shape = shape(mt%alogical1d)
      if (associated(mt%aint1d)) mem_shape = shape(mt%aint1d)
      if (associated(mt%aint2d)) mem_shape = shape(mt%aint2d)
      if (associated(mt%aint3d)) mem_shape = shape(mt%aint3d)
      if (associated(mt%adbl1d)) mem_shape = shape(mt%adbl1d)
      if (associated(mt%adbl2d)) mem_shape = shape(mt%adbl2d)
      if (associated(mt%adbl3d)) mem_shape = shape(mt%adbl3d)
      if (associated(mt%strsclr)) mem_shape = shape(mt%strsclr)
      if (associated(mt%astr1d)) mem_shape = shape(mt%astr1d)
      if (associated(mt%acharstr1d)) mem_shape = shape(mt%acharstr1d)
      ! -- to communicate failure
    else
      mem_shape(1) = -1
    end if
  end subroutine get_mem_shape

  !> @ brief Get the number of elements for this variable
  !!
  !! Returns with isize = -1 when not found.
  !! Return 1 for scalars.
  !<
  subroutine get_isize(name, mem_path, isize)
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where the variable is stored
    integer(I4B), intent(out) :: isize !< number of elements (flattened)
    ! -- local
    type(MemoryType), pointer :: mt => null()
    logical(LGP) :: found
    logical(LGP) :: terminate
    ! -- code
    !
    ! -- initialize isize to a value to communicate failure
    isize = -1
    !
    ! -- don't exit program if variable not found
    terminate = .false.
    !
    ! -- get the entry from the memory manager
    call get_from_memorystore(name, mem_path, mt, found, terminate)
    !
    ! -- set isize
    if (found) then
      isize = mt%isize
    end if
  end subroutine get_isize

  !> @ brief Get a memory type entry from the memory list
  !!
  !! Default value for @par check is .true. which means that this
  !! routine will kill the program when the memory entry cannot be found.
  !<
  subroutine get_from_memorystore(name, mem_path, mt, found, check)
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where the variable is stored
    type(MemoryType), pointer, intent(inout) :: mt !< memory type entry
    logical(LGP), intent(out) :: found !< set to .true. when found
    logical(LGP), intent(in), optional :: check !< to suppress aborting the program when not found,
                                                !! set check = .false.
    ! -- local
    logical(LGP) check_opt
    ! -- code
    mt => memorystore%get(name, mem_path)
    found = associated(mt)

    check_opt = .true.
    if (present(check)) then
      check_opt = check
    end if
    if (check_opt) then
      if (.not. found) then
        errmsg = "Programming error in memory manager. Variable '"// &
                 trim(name)//"' in '"//trim(mem_path)//"' cannot be "// &
                 "assigned because it does not exist in memory manager."
        call store_error(errmsg, terminate=.TRUE.)
      end if
    end if
  end subroutine get_from_memorystore

  !> @brief Issue allocation error message and stop program execution
  !<
  subroutine allocate_error(varname, mem_path, istat, isize)
    character(len=*), intent(in) :: varname !< variable name
    character(len=*), intent(in) :: mem_path !< path where the variable is stored
    integer(I4B), intent(in) :: istat !< status code
    integer(I4B), intent(in) :: isize !< size of allocation
    ! -- local
    character(len=20) :: csize
    character(len=20) :: cstat
    ! -- code
    !
    ! -- initialize character variables
    write (csize, '(i0)') isize
    write (cstat, '(i0)') istat
    !
    ! -- create error message
    errmsg = "Error trying to allocate memory. Path '"//trim(mem_path)// &
             "' variable name '"//trim(varname)//"' size '"//trim(csize)// &
             "'. Error message is '"//trim(adjustl(errmsg))// &
             "'. Status code is "//trim(cstat)//'.'
    !
    ! -- store error and stop program execution
    call store_error(errmsg, terminate=.TRUE.)
  end subroutine allocate_error

  !> @brief Allocate a logical scalar
  !<
  subroutine allocate_logical(sclr, name, mem_path)
    logical(LGP), pointer, intent(inout) :: sclr !< variable for allocation
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where the variable is stored
    ! -- local
    integer(I4B) :: istat
    type(MemoryType), pointer :: mt
    ! -- code
    !
    ! -- check variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- allocate the logical scalar
    allocate (sclr, stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, 1)
    end if
    !
    ! -- update counter
    nvalues_alogical = nvalues_alogical + 1
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%logicalsclr => sclr
    mt%element_size = LGP
    mt%isize = 1
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a)") 'LOGICAL'
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine allocate_logical

  !> @brief Allocate a 1-dimensional logical array
  !<
  subroutine allocate_logical1d(alog, nrow, name, mem_path)
    logical(LGP), dimension(:), pointer, contiguous, intent(inout) :: alog !< variable for allocation
    integer(I4B), intent(in) :: nrow !< number of rows
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! --local
    type(MemoryType), pointer :: mt
    integer(I4B) :: istat
    integer(I4B) :: isize
    ! -- code
    !
    ! -- check variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- set isize
    isize = nrow
    !
    ! -- allocate logical array
    allocate (alog(nrow), stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    !
    ! -- update counter
    nvalues_alogical = nvalues_alogical + isize
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%alogical1d => alog
    mt%element_size = LGP
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' (',i0,')')") 'LOGICAL', isize
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine allocate_logical1d

  !> @brief Allocate a character string
  !<
  subroutine allocate_str(sclr, ilen, name, mem_path)
    integer(I4B), intent(in) :: ilen !< string length
    character(len=ilen), pointer, intent(inout) :: sclr !< variable for allocation
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where the variable is stored
    ! -- local
    integer(I4B) :: istat
    type(MemoryType), pointer :: mt
    ! -- format
    ! -- code
    !
    ! -- make sure ilen is greater than 0
    if (ilen < 1) then
      errmsg = 'Programming error in allocate_str. ILEN must be greater than 0.'
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- check variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- allocate string
    allocate (character(len=ilen) :: sclr, stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, 1)
    end if
    !
    ! -- set sclr to a empty string
    sclr = ' '
    !
    ! -- update counter
    nvalues_astr = nvalues_astr + ilen
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%strsclr => sclr
    mt%element_size = ilen
    mt%isize = 1
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' LEN=',i0)") 'STRING', ilen
    !
    ! -- add defined length string to the memory manager list
    call memorystore%add(mt)
  end subroutine allocate_str

  !> @brief Allocate a 1-dimensional defined length string array
  !<
  subroutine allocate_str1d(astr1d, ilen, nrow, name, mem_path)
    integer(I4B), intent(in) :: ilen !< string length
    character(len=ilen), dimension(:), &
      pointer, contiguous, intent(inout) :: astr1d !< variable for allocation
    integer(I4B), intent(in) :: nrow !< number of strings in array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where the variable is stored
    ! -- local variables
    type(MemoryType), pointer :: mt
    character(len=ilen) :: string
    integer(I4B) :: n
    integer(I4B) :: istat
    integer(I4B) :: isize
    ! -- code
    !
    ! -- initialize string
    string = ''
    !
    ! -- make sure ilen is greater than 0
    if (ilen < 1) then
      errmsg = 'Programming error in allocate_str1d. '// &
               'ILEN must be greater than 0.'
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- check variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- calculate isize
    isize = nrow
    !
    ! -- allocate defined length string array
    allocate (character(len=ilen) :: astr1d(nrow), stat=istat, errmsg=errmsg)
    !
    ! -- check for error condition
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    !
    ! -- fill deferred length string with empty string
    do n = 1, nrow
      astr1d(n) = string
    end do
    !
    ! -- update counter
    nvalues_astr = nvalues_astr + isize
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%astr1d => astr1d
    mt%element_size = ilen
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' LEN=',i0,' (',i0,')')") 'STRING', ilen, nrow
    !
    ! -- add deferred length character array to the memory manager list
    call memorystore%add(mt)
  end subroutine allocate_str1d

  !> @brief Allocate a 1-dimensional array of deferred-length CharacterStringType
  !<
  subroutine allocate_charstr1d(acharstr1d, ilen, nrow, name, mem_path)
    type(CharacterStringType), dimension(:), &
      pointer, contiguous, intent(inout) :: acharstr1d !< variable for allocation
    integer(I4B), intent(in) :: ilen !< string length
    integer(I4B), intent(in) :: nrow !< number of strings in array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where the variable is stored
    ! -- local variables
    character(len=ilen) :: string
    type(MemoryType), pointer :: mt
    integer(I4B) :: n
    integer(I4B) :: istat
    integer(I4B) :: isize
    ! -- code
    !
    ! -- initialize string
    string = ''
    !
    ! -- check variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- calculate isize
    isize = nrow
    !
    ! -- allocate deferred length string array
    allocate (acharstr1d(nrow), stat=istat, errmsg=errmsg)
    !
    ! -- check for error condition
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    !
    ! -- fill deferred length string with empty string
    do n = 1, nrow
      acharstr1d(n) = string
    end do
    !
    ! -- update counter
    nvalues_astr = nvalues_astr + isize
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%acharstr1d => acharstr1d
    mt%element_size = ilen
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' LEN=',i0,' (',i0,')')") 'STRING', ilen, nrow
    !
    ! -- add deferred length character array to the memory manager list
    call memorystore%add(mt)
  end subroutine allocate_charstr1d

  !> @brief Allocate a integer scalar
  !<
  subroutine allocate_int(sclr, name, mem_path)
    integer(I4B), pointer, intent(inout) :: sclr !< variable for allocation
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where the variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    integer(I4B) :: istat
    ! -- code
    !
    ! -- check variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- allocate integer scalar
    allocate (sclr, stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, 1)
    end if
    !
    ! -- update counter
    nvalues_aint = nvalues_aint + 1
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%intsclr => sclr
    mt%element_size = I4B
    mt%isize = 1
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a)") 'INTEGER'
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine allocate_int

  !> @brief Allocate a 1-dimensional integer array
  !<
  subroutine allocate_int1d(aint, nrow, name, mem_path)
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: aint !< variable for allocation
    integer(I4B), intent(in) :: nrow !< integer array number of rows
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! --local
    type(MemoryType), pointer :: mt
    integer(I4B) :: istat
    integer(I4B) :: isize
    ! -- code
    !
    ! -- check variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- set isize
    isize = nrow
    !
    ! -- allocate integer array
    allocate (aint(nrow), stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    !
    ! -- update counter
    nvalues_aint = nvalues_aint + isize
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%aint1d => aint
    mt%element_size = I4B
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' (',i0,')')") 'INTEGER', isize
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine allocate_int1d

  !> @brief Allocate a 2-dimensional integer array
  !<
  subroutine allocate_int2d(aint, ncol, nrow, name, mem_path)
    integer(I4B), dimension(:, :), pointer, contiguous, intent(inout) :: aint !< variable for allocation
    integer(I4B), intent(in) :: ncol !< number of columns
    integer(I4B), intent(in) :: nrow !< number of rows
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    integer(I4B) :: istat
    integer(I4B) :: isize
    ! -- code
    !
    ! -- check the variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- set isize
    isize = ncol * nrow
    !
    ! -- allocate the integer array
    allocate (aint(ncol, nrow), stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    !
    ! -- update the counter
    nvalues_aint = nvalues_aint + isize
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%aint2d => aint
    mt%element_size = I4B
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' (',i0,',',i0,')')") 'INTEGER', ncol, nrow
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine allocate_int2d

  !> @brief Allocate a 3-dimensional integer array
  !<
  subroutine allocate_int3d(aint, ncol, nrow, nlay, name, mem_path)
    integer(I4B), dimension(:, :, :), pointer, contiguous, intent(inout) :: aint !< variable for allocation
    integer(I4B), intent(in) :: ncol !< number of columns
    integer(I4B), intent(in) :: nrow !< number of rows
    integer(I4B), intent(in) :: nlay !< number of layers
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    integer(I4B) :: istat
    integer(I4B) :: isize
    ! -- code
    !
    ! -- check variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- set isize
    isize = ncol * nrow * nlay
    !
    ! -- allocate integer array
    allocate (aint(ncol, nrow, nlay), stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    !
    ! -- update counter
    nvalues_aint = nvalues_aint + isize
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%aint3d => aint
    mt%element_size = I4B
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' (',i0,',',i0,',',i0,')')") 'INTEGER', ncol, &
      nrow, nlay
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine allocate_int3d

  !> @brief Allocate a real scalar
  !<
  subroutine allocate_dbl(sclr, name, mem_path)
    real(DP), pointer, intent(inout) :: sclr !< variable for allocation
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    integer(I4B) :: istat
    ! -- code
    !
    ! -- check variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- allocate real scalar
    allocate (sclr, stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, 1)
    end if
    !
    ! -- update counter
    nvalues_aint = nvalues_aint + 1
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%dblsclr => sclr
    mt%element_size = DP
    mt%isize = 1
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a)") 'DOUBLE'
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine allocate_dbl

  !> @brief Allocate a 1-dimensional real array
  !<
  subroutine allocate_dbl1d(adbl, nrow, name, mem_path)
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: adbl !< variable for allocation
    integer(I4B), intent(in) :: nrow !< number of rows
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    integer(I4B) :: istat
    integer(I4B) :: isize
    ! -- code
    !
    ! -- check the variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- set isize
    isize = nrow
    !
    ! -- allocate the real array
    allocate (adbl(nrow), stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    !
    ! -- update counter
    nvalues_adbl = nvalues_adbl + isize
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%adbl1d => adbl
    mt%element_size = DP
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' (',i0,')')") 'DOUBLE', isize
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine allocate_dbl1d

  !> @brief Allocate a 2-dimensional real array
  !<
  subroutine allocate_dbl2d(adbl, ncol, nrow, name, mem_path)
    real(DP), dimension(:, :), pointer, contiguous, intent(inout) :: adbl !< variable for allocation
    integer(I4B), intent(in) :: ncol !< number of columns
    integer(I4B), intent(in) :: nrow !< number of rows
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    integer(I4B) :: istat
    integer(I4B) :: isize
    ! -- code
    !
    ! -- check the variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- set isize
    isize = ncol * nrow
    !
    ! -- allocate the real array
    allocate (adbl(ncol, nrow), stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    !
    ! -- update counter
    nvalues_adbl = nvalues_adbl + isize
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%adbl2d => adbl
    mt%element_size = DP
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' (',i0,',',i0,')')") 'DOUBLE', ncol, nrow
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine allocate_dbl2d

  !> @brief Allocate a 3-dimensional real array
  !<
  subroutine allocate_dbl3d(adbl, ncol, nrow, nlay, name, mem_path)
    real(DP), dimension(:, :, :), pointer, contiguous, intent(inout) :: adbl !< variable for allocation
    integer(I4B), intent(in) :: ncol !< number of columns
    integer(I4B), intent(in) :: nrow !< number of rows
    integer(I4B), intent(in) :: nlay !< number of layers
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    integer(I4B) :: istat
    integer(I4B) :: isize
    ! -- code
    !
    ! -- check the variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- set isize
    isize = ncol * nrow * nlay
    !
    ! -- allocate the real array
    allocate (adbl(ncol, nrow, nlay), stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    !
    ! -- update the counter
    nvalues_adbl = nvalues_adbl + isize
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%adbl3d => adbl
    mt%element_size = DP
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' (',i0,',',i0,',',i0,')')") 'DOUBLE', ncol, &
      nrow, nlay
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine allocate_dbl3d

  !> @brief Check in an existing 1d logical array with a new address (name + path)
  !<
  subroutine checkin_logical1d(alog, name, mem_path, name2, mem_path2)
    logical(LGP), dimension(:), pointer, contiguous, intent(in) :: alog !< the existing array
    character(len=*), intent(in) :: name !< new variable name
    character(len=*), intent(in) :: mem_path !< new path where variable is stored
    character(len=*), intent(in) :: name2 !< existing variable name
    character(len=*), intent(in) :: mem_path2 !< existing path where variable is stored
    ! --local
    type(MemoryType), pointer :: mt
    integer(I4B) :: isize
    ! -- code
    !
    ! -- check variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- set isize
    isize = size(alog)
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%alogical1d => alog
    mt%element_size = LGP
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' (',i0,')')") 'LOGICAL', isize
    !
    ! -- set master information
    mt%master = .false.
    mt%mastername = name2
    mt%masterPath = mem_path2
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine checkin_logical1d

  !> @brief Check in an existing 1d integer array with a new address (name + path)
  !<
  subroutine checkin_int1d(aint, name, mem_path, name2, mem_path2)
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: aint !< the existing array
    character(len=*), intent(in) :: name !< new variable name
    character(len=*), intent(in) :: mem_path !< new path where variable is stored
    character(len=*), intent(in) :: name2 !< existing variable name
    character(len=*), intent(in) :: mem_path2 !< existing path where variable is stored
    ! --local
    type(MemoryType), pointer :: mt
    integer(I4B) :: isize
    ! -- code
    !
    ! -- check variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- set isize
    isize = size(aint)
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%aint1d => aint
    mt%element_size = I4B
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' (',i0,')')") 'INTEGER', isize
    !
    ! -- set master information
    mt%master = .false.
    mt%mastername = name2
    mt%masterPath = mem_path2
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine checkin_int1d

  !> @brief Check in an existing 2d integer array with a new address (name + path)
  !<
  subroutine checkin_int2d(aint2d, name, mem_path, name2, mem_path2)
    integer(I4B), dimension(:, :), pointer, contiguous, intent(inout) :: aint2d !< the existing 2d array
    character(len=*), intent(in) :: name !< new variable name
    character(len=*), intent(in) :: mem_path !< new path where variable is stored
    character(len=*), intent(in) :: name2 !< existing variable name
    character(len=*), intent(in) :: mem_path2 !< existing path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    integer(I4B) :: ncol, nrow, isize
    ! -- code
    !
    ! -- check the variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- set isize
    ncol = size(aint2d, dim=1)
    nrow = size(aint2d, dim=2)
    isize = ncol * nrow
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%aint2d => aint2d
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' (',i0,',',i0,')')") 'INTEGER', ncol, nrow
    !
    ! -- set master information
    mt%master = .false.
    mt%mastername = name2
    mt%masterPath = mem_path2
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine checkin_int2d

  !> @brief Check in an existing 1d double precision array with a new address (name + path)
  !<
  subroutine checkin_dbl1d(adbl, name, mem_path, name2, mem_path2)
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: adbl !< the existing array
    character(len=*), intent(in) :: name !< new variable name
    character(len=*), intent(in) :: mem_path !< new path where variable is stored
    character(len=*), intent(in) :: name2 !< existing variable name
    character(len=*), intent(in) :: mem_path2 !< existing path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    integer(I4B) :: isize
    ! -- code
    !
    ! -- check the variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- set isize
    isize = size(adbl)
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%adbl1d => adbl
    mt%element_size = DP
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' (',i0,')')") 'DOUBLE', isize
    !
    ! -- set master information
    mt%master = .false.
    mt%mastername = name2
    mt%masterPath = mem_path2
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine checkin_dbl1d

  !> @brief Check in an existing 2d double precision array with a new address (name + path)
  !<
  subroutine checkin_dbl2d(adbl2d, name, mem_path, name2, mem_path2)
    real(DP), dimension(:, :), pointer, contiguous, intent(inout) :: adbl2d !< the existing 2d array
    character(len=*), intent(in) :: name !< new variable name
    character(len=*), intent(in) :: mem_path !< new path where variable is stored
    character(len=*), intent(in) :: name2 !< existing variable name
    character(len=*), intent(in) :: mem_path2 !< existing path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    integer(I4B) :: ncol, nrow, isize
    ! -- code
    !
    ! -- check the variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- set isize
    ncol = size(adbl2d, dim=1)
    nrow = size(adbl2d, dim=2)
    isize = ncol * nrow
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%adbl2d => adbl2d
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' (',i0,',',i0,')')") 'DOUBLE', ncol, nrow
    !
    ! -- set master information
    mt%master = .false.
    mt%mastername = name2
    mt%masterPath = mem_path2
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine checkin_dbl2d

  !> @brief Check in an existing 1d CharacterStringType array with a new address (name + path)
  !<
  subroutine checkin_charstr1d(acharstr1d, ilen, name, mem_path, name2, mem_path2)
    type(CharacterStringType), dimension(:), &
      pointer, contiguous, intent(inout) :: acharstr1d !< the existing array
    integer(I4B), intent(in) :: ilen
    character(len=*), intent(in) :: name !< new variable name
    character(len=*), intent(in) :: mem_path !< new path where variable is stored
    character(len=*), intent(in) :: name2 !< existing variable name
    character(len=*), intent(in) :: mem_path2 !< existing path where variable is stored
    ! --local
    type(MemoryType), pointer :: mt
    integer(I4B) :: isize
    ! -- code
    !
    ! -- check variable name length
    call mem_check_length(name, LENVARNAME, "variable")
    !
    ! -- set isize
    isize = size(acharstr1d)
    !
    ! -- allocate memory type
    allocate (mt)
    !
    ! -- set memory type
    mt%acharstr1d => acharstr1d
    mt%element_size = ilen
    mt%isize = isize
    mt%name = name
    mt%path = mem_path
    write (mt%memtype, "(a,' LEN=',i0,' (',i0,')')") 'STRING', ilen, isize
    !
    ! -- set master information
    mt%master = .false.
    mt%mastername = name2
    mt%masterPath = mem_path2
    !
    ! -- add memory type to the memory list
    call memorystore%add(mt)
  end subroutine checkin_charstr1d

  !> @brief Reallocate a 1-dimensional defined length string array
  !<
  subroutine reallocate_str1d(astr, ilen, nrow, name, mem_path)
    integer(I4B), intent(in) :: ilen !< string length
    integer(I4B), intent(in) :: nrow !< number of rows
    character(len=ilen), dimension(:), pointer, contiguous, intent(inout) :: astr !< the reallocated string array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    character(len=ilen), dimension(:), allocatable :: astrtemp
    integer(I4B) :: istat
    integer(I4B) :: isize
    integer(I4B) :: isize_old
    integer(I4B) :: nrow_old
    integer(I4B) :: n
    !
    ! -- Find and assign mt
    call get_from_memorystore(name, mem_path, mt, found)
    !
    ! -- reallocate astr1d
    if (found) then
      isize_old = mt%isize
      if (isize_old > 0) then
        nrow_old = size(astr)
      else
        nrow_old = 0
      end if
      !
      ! -- calculate isize
      isize = nrow
      !
      ! -- allocate astrtemp
      allocate (astrtemp(nrow), stat=istat, errmsg=errmsg)
      if (istat /= 0) then
        call allocate_error(name, mem_path, istat, isize)
      end if
      !
      ! -- copy existing values
      do n = 1, nrow_old
        astrtemp(n) = astr(n)
      end do
      !
      ! -- fill new values with missing values
      do n = nrow_old + 1, nrow
        astrtemp(n) = ''
      end do
      !
      ! -- deallocate mt pointer, repoint, recalculate isize
      deallocate (astr)
      !
      ! -- allocate astr1d
      allocate (astr(nrow), stat=istat, errmsg=errmsg)
      if (istat /= 0) then
        call allocate_error(name, mem_path, istat, isize)
      end if
      !
      ! -- fill the reallocate character array
      do n = 1, nrow
        astr(n) = astrtemp(n)
      end do
      !
      ! -- deallocate temporary storage
      deallocate (astrtemp)
      !
      ! -- reset memory manager values
      mt%astr1d => astr
      mt%element_size = ilen
      mt%isize = isize
      mt%nrealloc = mt%nrealloc + 1
      mt%master = .true.
      nvalues_astr = nvalues_astr + isize - isize_old
      write (mt%memtype, "(a,' LEN=',i0,' (',i0,')')") 'STRING', ilen, nrow
    else
      errmsg = "Programming error, variable '"//trim(name)//"' from '"// &
               trim(mem_path)//"' is not defined in the memory manager. Use "// &
               "mem_allocate instead."
      call store_error(errmsg, terminate=.TRUE.)
    end if
  end subroutine reallocate_str1d

  !> @brief Reallocate a 1-dimensional deferred length string array
  !<
  subroutine reallocate_charstr1d(acharstr1d, ilen, nrow, name, mem_path)
    type(CharacterStringType), dimension(:), pointer, contiguous, &
      intent(inout) :: acharstr1d !< the reallocated charstring array
    integer(I4B), intent(in) :: ilen !< string length
    integer(I4B), intent(in) :: nrow !< number of rows
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    type(CharacterStringType), dimension(:), allocatable :: astrtemp
    character(len=ilen) :: string
    integer(I4B) :: istat
    integer(I4B) :: isize
    integer(I4B) :: isize_old
    integer(I4B) :: nrow_old
    integer(I4B) :: n
    !
    ! -- Initialize string
    string = ''
    !
    ! -- Find and assign mt
    call get_from_memorystore(name, mem_path, mt, found)
    !
    ! -- reallocate astr1d
    if (found) then
      isize_old = mt%isize
      if (isize_old > 0) then
        nrow_old = size(acharstr1d)
      else
        nrow_old = 0
      end if
      !
      ! -- calculate isize
      isize = nrow
      !
      ! -- allocate astrtemp
      allocate (astrtemp(nrow), stat=istat, errmsg=errmsg)
      if (istat /= 0) then
        call allocate_error(name, mem_path, istat, isize)
      end if
      !
      ! -- copy existing values
      do n = 1, nrow_old
        astrtemp(n) = acharstr1d(n)
        call acharstr1d(n)%destroy()
      end do
      !
      ! -- fill new values with missing values
      do n = nrow_old + 1, nrow
        astrtemp(n) = string
      end do
      !
      ! -- deallocate mt pointer, repoint, recalculate isize
      deallocate (acharstr1d)
      !
      ! -- allocate astr1d
      allocate (acharstr1d(nrow), stat=istat, errmsg=errmsg)
      if (istat /= 0) then
        call allocate_error(name, mem_path, istat, isize)
      end if
      !
      ! -- fill the reallocated character array
      do n = 1, nrow
        acharstr1d(n) = astrtemp(n)
        call astrtemp(n)%destroy()
      end do
      !
      ! -- deallocate temporary storage
      deallocate (astrtemp)
      !
      ! -- reset memory manager values
      mt%acharstr1d => acharstr1d
      mt%element_size = ilen
      mt%isize = isize
      mt%nrealloc = mt%nrealloc + 1
      mt%master = .true.
      nvalues_astr = nvalues_astr + isize - isize_old
      write (mt%memtype, "(a,' LEN=',i0,' (',i0,')')") 'STRING', ilen, nrow
    else
      errmsg = "Programming error, variable '"//trim(name)//"' from '"// &
               trim(mem_path)//"' is not defined in the memory manager. Use "// &
               "mem_allocate instead."
      call store_error(errmsg, terminate=.TRUE.)
    end if
  end subroutine reallocate_charstr1d

  !> @brief Reallocate a 1-dimensional logical array
  !<
  subroutine reallocate_logical1d(alog, nrow, name, mem_path)
    logical(LGP), dimension(:), pointer, contiguous, intent(inout) :: alog !< the reallocated logical array
    integer(I4B), intent(in) :: nrow !< number of rows
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    integer(I4B) :: istat
    integer(I4B) :: isize
    integer(I4B) :: i
    integer(I4B) :: isizeold
    integer(I4B) :: ifill
    ! -- code
    !
    ! -- Find and assign mt
    call get_from_memorystore(name, mem_path, mt, found)
    !
    ! -- Allocate aint and then refill
    isize = nrow
    isizeold = size(mt%alogical1d)
    ifill = min(isizeold, isize)
    allocate (alog(nrow), stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    do i = 1, ifill
      alog(i) = mt%alogical1d(i)
    end do
    !
    ! -- deallocate mt pointer, repoint, recalculate isize
    deallocate (mt%alogical1d)
    mt%alogical1d => alog
    mt%element_size = LGP
    mt%isize = isize
    mt%nrealloc = mt%nrealloc + 1
    mt%master = .true.
    nvalues_alogical = nvalues_alogical + isize - isizeold
  end subroutine reallocate_logical1d

  !> @brief Reallocate a 1-dimensional integer array
  !<
  subroutine reallocate_int1d(aint, nrow, name, mem_path)
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: aint !< the reallocated integer array
    integer(I4B), intent(in) :: nrow !< number of rows
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    integer(I4B) :: istat
    integer(I4B) :: isize
    integer(I4B) :: i
    integer(I4B) :: isizeold
    integer(I4B) :: ifill
    ! -- code
    !
    ! -- Find and assign mt
    call get_from_memorystore(name, mem_path, mt, found)
    !
    ! -- Allocate aint and then refill
    isize = nrow
    isizeold = size(mt%aint1d)
    ifill = min(isizeold, isize)
    allocate (aint(nrow), stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    do i = 1, ifill
      aint(i) = mt%aint1d(i)
    end do
    !
    ! -- deallocate mt pointer, repoint, recalculate isize
    deallocate (mt%aint1d)
    mt%aint1d => aint
    mt%element_size = I4B
    mt%isize = isize
    mt%nrealloc = mt%nrealloc + 1
    mt%master = .true.
    nvalues_aint = nvalues_aint + isize - isizeold
  end subroutine reallocate_int1d

  !> @brief Reallocate a 2-dimensional integer array
  !<
  subroutine reallocate_int2d(aint, ncol, nrow, name, mem_path)
    integer(I4B), dimension(:, :), pointer, contiguous, intent(inout) :: aint !< the reallocated 2d integer array
    integer(I4B), intent(in) :: ncol !< number of columns
    integer(I4B), intent(in) :: nrow !< number of rows
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    integer(I4B) :: istat
    integer(I4B), dimension(2) :: ishape
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: isize
    integer(I4B) :: isizeold
    ! -- code
    !
    ! -- Find and assign mt
    call get_from_memorystore(name, mem_path, mt, found)
    !
    ! -- Allocate aint and then refill
    ishape = shape(mt%aint2d)
    isize = nrow * ncol
    isizeold = ishape(1) * ishape(2)
    allocate (aint(ncol, nrow), stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    do i = 1, ishape(2)
      do j = 1, ishape(1)
        aint(j, i) = mt%aint2d(j, i)
      end do
    end do
    !
    ! -- deallocate mt pointer, repoint, recalculate isize
    deallocate (mt%aint2d)
    mt%aint2d => aint
    mt%element_size = I4B
    mt%isize = isize
    mt%nrealloc = mt%nrealloc + 1
    mt%master = .true.
    nvalues_aint = nvalues_aint + isize - isizeold
    write (mt%memtype, "(a,' (',i0,',',i0,')')") 'INTEGER', ncol, nrow
  end subroutine reallocate_int2d

  !> @brief Reallocate a 1-dimensional real array
  !<
  subroutine reallocate_dbl1d(adbl, nrow, name, mem_path)
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: adbl !< the reallocated 1d real array
    integer(I4B), intent(in) :: nrow !< number of rows
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    integer(I4B) :: istat
    integer(I4B) :: isize
    integer(I4B) :: i
    integer(I4B) :: isizeold
    integer(I4B) :: ifill
    logical(LGP) :: found
    ! -- code
    !
    ! -- Find and assign mt
    call get_from_memorystore(name, mem_path, mt, found)
    !
    ! -- Allocate adbl and then refill
    isize = nrow
    isizeold = size(mt%adbl1d)
    ifill = min(isizeold, isize)
    allocate (adbl(nrow), stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    do i = 1, ifill
      adbl(i) = mt%adbl1d(i)
    end do
    !
    ! -- deallocate mt pointer, repoint, recalculate isize
    deallocate (mt%adbl1d)
    mt%adbl1d => adbl
    mt%element_size = DP
    mt%isize = isize
    mt%nrealloc = mt%nrealloc + 1
    mt%master = .true.
    nvalues_adbl = nvalues_adbl + isize - isizeold
    write (mt%memtype, "(a,' (',i0,')')") 'DOUBLE', isize
  end subroutine reallocate_dbl1d

  !> @brief Reallocate a 2-dimensional real array
  !<
  subroutine reallocate_dbl2d(adbl, ncol, nrow, name, mem_path)
    real(DP), dimension(:, :), pointer, contiguous, intent(inout) :: adbl !< the reallocated 2d real array
    integer(I4B), intent(in) :: ncol !< number of columns
    integer(I4B), intent(in) :: nrow !< number of rows
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    integer(I4B) :: istat
    integer(I4B), dimension(2) :: ishape
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: isize
    integer(I4B) :: isizeold
    ! -- code
    !
    ! -- Find and assign mt
    call get_from_memorystore(name, mem_path, mt, found)
    !
    ! -- Allocate adbl and then refill
    ishape = shape(mt%adbl2d)
    isize = nrow * ncol
    isizeold = ishape(1) * ishape(2)
    allocate (adbl(ncol, nrow), stat=istat, errmsg=errmsg)
    if (istat /= 0) then
      call allocate_error(name, mem_path, istat, isize)
    end if
    do i = 1, ishape(2)
      do j = 1, ishape(1)
        adbl(j, i) = mt%adbl2d(j, i)
      end do
    end do
    !
    ! -- deallocate mt pointer, repoint, recalculate isize
    deallocate (mt%adbl2d)
    mt%adbl2d => adbl
    mt%element_size = DP
    mt%isize = isize
    mt%nrealloc = mt%nrealloc + 1
    mt%master = .true.
    nvalues_adbl = nvalues_adbl + isize - isizeold
    write (mt%memtype, "(a,' (',i0,',',i0,')')") 'DOUBLE', ncol, nrow
  end subroutine reallocate_dbl2d

  !> @brief Set pointer to a logical scalar
  !<
  subroutine setptr_logical(sclr, name, mem_path)
    logical(LGP), pointer, intent(inout) :: sclr !< pointer to logical scalar
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    sclr => mt%logicalsclr
  end subroutine setptr_logical

  !> @brief Set pointer to integer scalar
  !<
  subroutine setptr_int(sclr, name, mem_path)
    integer(I4B), pointer, intent(inout) :: sclr !< pointer to integer scalar
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    sclr => mt%intsclr
  end subroutine setptr_int

  !> @brief Set pointer to 1d logical array
  !<
  subroutine setptr_logical1d(alog, name, mem_path)
    logical(LGP), dimension(:), pointer, contiguous, intent(inout) :: alog !< pointer to 1d logical array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    alog => mt%alogical1d
  end subroutine setptr_logical1d

  !> @brief Set pointer to 1d integer array
  !<
  subroutine setptr_int1d(aint, name, mem_path)
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: aint !< pointer to 1d integer array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    aint => mt%aint1d
  end subroutine setptr_int1d

  !> @brief Set pointer to 2d integer array
  !<
  subroutine setptr_int2d(aint, name, mem_path)
    integer(I4B), dimension(:, :), pointer, contiguous, intent(inout) :: aint !< pointer to 2d integer array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    aint => mt%aint2d
  end subroutine setptr_int2d

  !> @brief Set pointer to 3d integer array
  !<
  subroutine setptr_int3d(aint, name, mem_path)
    integer(I4B), dimension(:, :, :), pointer, contiguous, intent(inout) :: aint !< pointer to 3d integer array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    aint => mt%aint3d
  end subroutine setptr_int3d

  !> @brief Set pointer to a real scalar
  !<
  subroutine setptr_dbl(sclr, name, mem_path)
    real(DP), pointer, intent(inout) :: sclr !< pointer to a real scalar
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    sclr => mt%dblsclr
  end subroutine setptr_dbl

  !> @brief Set pointer to a 1d real array
  !<
  subroutine setptr_dbl1d(adbl, name, mem_path)
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: adbl !< pointer to 1d real array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    adbl => mt%adbl1d
  end subroutine setptr_dbl1d

  !> @brief Set pointer to a 2d real array
  !<
  subroutine setptr_dbl2d(adbl, name, mem_path)
    real(DP), dimension(:, :), pointer, contiguous, intent(inout) :: adbl !< pointer to 2d real array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    adbl => mt%adbl2d
  end subroutine setptr_dbl2d

  !> @brief Set pointer to a 3d real array
  !<
  subroutine setptr_dbl3d(adbl, name, mem_path)
    real(DP), dimension(:, :, :), pointer, contiguous, intent(inout) :: adbl !< pointer to 3d real array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    adbl => mt%adbl3d
  end subroutine setptr_dbl3d

  !> @brief Set pointer to a string (scalar)
  !<
  subroutine setptr_str(asrt, name, mem_path)
    character(len=:), pointer :: asrt !< pointer to the character string
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    asrt => mt%strsclr
  end subroutine setptr_str

  !> @brief Set pointer to a fixed-length string array
  !<
  subroutine setptr_str1d(astr1d, name, mem_path)
    character(len=:), dimension(:), &
      pointer, contiguous, intent(inout) :: astr1d !< pointer to the string array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    select type (item => mt%astr1d)
    type is (character(*))
      astr1d => item
    class default
      astr1d => null()
    end select
  end subroutine setptr_str1d

  !> @brief Set pointer to an array of CharacterStringType
  !<
  subroutine setptr_charstr1d(acharstr1d, name, mem_path)
    type(CharacterStringType), dimension(:), pointer, contiguous, &
      intent(inout) :: acharstr1d !< the reallocated charstring array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    acharstr1d => mt%acharstr1d
  end subroutine setptr_charstr1d

  !> @brief Make a copy of a 1-dimensional logical array
  !<
  subroutine copyptr_logical1d(alog, name, mem_path, mem_path_copy)
    logical(LGP), dimension(:), pointer, contiguous, intent(inout) :: alog !< returned copy of 1d logical array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    character(len=*), intent(in), optional :: mem_path_copy !< optional path where the copy will be stored,
                                                            !! if passed then the copy is added to the
                                                            !! memory manager
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    integer(I4B) :: n
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    alog => null()
    ! -- check the copy into the memory manager
    if (present(mem_path_copy)) then
      call allocate_logical1d(alog, size(mt%alogical1d), mt%name, mem_path_copy)
      ! -- create a local copy
    else
      allocate (alog(size(mt%alogical1d)))
    end if
    do n = 1, size(mt%alogical1d)
      alog(n) = mt%alogical1d(n)
    end do
  end subroutine copyptr_logical1d

  !> @brief Make a copy of a 1-dimensional integer array
  !<
  subroutine copyptr_int1d(aint, name, mem_path, mem_path_copy)
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: aint !< returned copy of 1d integer array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    character(len=*), intent(in), optional :: mem_path_copy !< optional path where the copy will be stored,
                                                            !! if passed then the copy is added to the
                                                            !! memory manager
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    integer(I4B) :: n
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    aint => null()
    ! -- check the copy into the memory manager
    if (present(mem_path_copy)) then
      call allocate_int1d(aint, size(mt%aint1d), mt%name, mem_path_copy)
      ! -- create a local copy
    else
      allocate (aint(size(mt%aint1d)))
    end if
    do n = 1, size(mt%aint1d)
      aint(n) = mt%aint1d(n)
    end do
  end subroutine copyptr_int1d

  !> @brief Make a copy of a 2-dimensional integer array
  !<
  subroutine copyptr_int2d(aint, name, mem_path, mem_path_copy)
    integer(I4B), dimension(:, :), pointer, contiguous, intent(inout) :: aint !< returned copy of 2d integer array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    character(len=*), intent(in), optional :: mem_path_copy !< optional path where the copy will be stored,
                                                            !! if passed then the copy is added to the
                                                            !! memory manager
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: ncol
    integer(I4B) :: nrow
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    aint => null()
    ncol = size(mt%aint2d, dim=1)
    nrow = size(mt%aint2d, dim=2)
    ! -- check the copy into the memory manager
    if (present(mem_path_copy)) then
      call allocate_int2d(aint, ncol, nrow, mt%name, mem_path_copy)
      ! -- create a local copy
    else
      allocate (aint(ncol, nrow))
    end if
    do i = 1, nrow
      do j = 1, ncol
        aint(j, i) = mt%aint2d(j, i)
      end do
    end do
  end subroutine copyptr_int2d

  !> @brief Make a copy of a 1-dimensional real array
  !<
  subroutine copyptr_dbl1d(adbl, name, mem_path, mem_path_copy)
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: adbl !< returned copy of 1d real array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    character(len=*), intent(in), optional :: mem_path_copy !< optional path where the copy will be stored,
                                                            !! if passed then the copy is added to the
                                                            !! memory manager
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    integer(I4B) :: n
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    adbl => null()
    ! -- check the copy into the memory manager
    if (present(mem_path_copy)) then
      call allocate_dbl1d(adbl, size(mt%adbl1d), mt%name, mem_path_copy)
      ! -- create a local copy
    else
      allocate (adbl(size(mt%adbl1d)))
    end if
    do n = 1, size(mt%adbl1d)
      adbl(n) = mt%adbl1d(n)
    end do
  end subroutine copyptr_dbl1d

  !> @brief Make a copy of a 2-dimensional real array
  !<
  subroutine copyptr_dbl2d(adbl, name, mem_path, mem_path_copy)
    real(DP), dimension(:, :), pointer, contiguous, intent(inout) :: adbl !< returned copy of 2d real array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    character(len=*), intent(in), optional :: mem_path_copy !< optional path where the copy will be stored,
                                                            !! if passed then the copy is added to the
                                                            !! memory manager
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: ncol
    integer(I4B) :: nrow
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    adbl => null()
    ncol = size(mt%adbl2d, dim=1)
    nrow = size(mt%adbl2d, dim=2)
    ! -- check the copy into the memory manager
    if (present(mem_path_copy)) then
      call allocate_dbl2d(adbl, ncol, nrow, mt%name, mem_path_copy)
      ! -- create a local copy
    else
      allocate (adbl(ncol, nrow))
    end if
    do i = 1, nrow
      do j = 1, ncol
        adbl(j, i) = mt%adbl2d(j, i)
      end do
    end do
  end subroutine copyptr_dbl2d

  !> @brief Copy values from a 1-dimensional real array in the memory
  !< manager to a passed 1-dimensional real array
  subroutine copy_dbl1d(adbl, name, mem_path)
    real(DP), dimension(:), intent(inout) :: adbl !< target array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    logical(LGP) :: found
    integer(I4B) :: n
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    do n = 1, size(mt%adbl1d)
      adbl(n) = mt%adbl1d(n)
    end do
  end subroutine copy_dbl1d

  !> @brief Set the pointer for an integer scalar to
  !< a target array already stored in the memory manager
  subroutine reassignptr_int(sclr, name, mem_path, name_target, mem_path_target)
    integer(I4B), pointer, intent(inout) :: sclr !< pointer to integer scalar
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    character(len=*), intent(in) :: name_target !< name of target variable
    character(len=*), intent(in) :: mem_path_target !< path where target variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    type(MemoryType), pointer :: mt2
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    call get_from_memorystore(name_target, mem_path_target, mt2, found)
    if (associated(sclr)) then
      nvalues_aint = nvalues_aint - 1
      deallocate (sclr)
    end if
    sclr => mt2%intsclr
    mt%intsclr => sclr
    mt%element_size = I4B
    mt%isize = 1
    write (mt%memtype, "(a,' (',i0,')')") 'INTEGER', mt%isize
    !
    ! -- set master information
    mt%master = .false.
    mt%mastername = name_target
    mt%masterPath = mem_path_target
  end subroutine reassignptr_int

  !> @brief Set the pointer for a 1-dimensional logical array to
  !< a target array already stored in the memory manager
  subroutine reassignptr_logical1d(alog, name, mem_path, name_target, &
                                   mem_path_target)
    logical(LGP), dimension(:), pointer, contiguous, intent(inout) :: alog !< array pointer
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    character(len=*), intent(in) :: name_target !< name of target variable
    character(len=*), intent(in) :: mem_path_target !< path where target variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    type(MemoryType), pointer :: mt2
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    call get_from_memorystore(name_target, mem_path_target, mt2, found)
    if (size(alog) > 0) then
      nvalues_alogical = nvalues_alogical - size(alog)
      deallocate (alog)
    end if
    alog => mt2%alogical1d
    mt%alogical1d => alog
    mt%element_size = LGP
    mt%isize = size(alog)
    write (mt%memtype, "(a,' (',i0,')')") 'LOGICAL', mt%isize
    !
    ! -- set master information
    mt%master = .false.
    mt%mastername = name_target
    mt%masterPath = mem_path_target
  end subroutine reassignptr_logical1d

  !> @brief Set the pointer for a 1-dimensional integer array to
  !< a target array already stored in the memory manager
  subroutine reassignptr_int1d(aint, name, mem_path, name_target, mem_path_target)
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: aint !< pointer to 1d integer array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    character(len=*), intent(in) :: name_target !< name of target variable
    character(len=*), intent(in) :: mem_path_target !< path where target variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    type(MemoryType), pointer :: mt2
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    call get_from_memorystore(name_target, mem_path_target, mt2, found)
    if (size(aint) > 0) then
      nvalues_aint = nvalues_aint - size(aint)
      deallocate (aint)
    end if
    aint => mt2%aint1d
    mt%aint1d => aint
    mt%element_size = I4B
    mt%isize = size(aint)
    write (mt%memtype, "(a,' (',i0,')')") 'INTEGER', mt%isize
    !
    ! -- set master information
    mt%master = .false.
    mt%mastername = name_target
    mt%masterPath = mem_path_target
  end subroutine reassignptr_int1d

  !> @brief Set the pointer for a 2-dimensional integer array to
  !< a target array already stored in the memory manager
  subroutine reassignptr_int2d(aint, name, mem_path, name_target, mem_path_target)
    integer(I4B), dimension(:, :), pointer, contiguous, intent(inout) :: aint !< pointer to 2d integer array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    character(len=*), intent(in) :: name_target !< name of target variable
    character(len=*), intent(in) :: mem_path_target !< path where target variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    type(MemoryType), pointer :: mt2
    logical(LGP) :: found
    integer(I4B) :: ncol
    integer(I4B) :: nrow
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    call get_from_memorystore(name_target, mem_path_target, mt2, found)
    if (size(aint) > 0) then
      nvalues_aint = nvalues_aint - size(aint)
      deallocate (aint)
    end if
    aint => mt2%aint2d
    mt%aint2d => aint
    mt%element_size = I4B
    mt%isize = size(aint)
    ncol = size(aint, dim=1)
    nrow = size(aint, dim=2)
    write (mt%memtype, "(a,' (',i0,',',i0,')')") 'INTEGER', ncol, nrow
    !
    ! -- set master information
    mt%master = .false.
    mt%mastername = name_target
    mt%masterPath = mem_path_target
  end subroutine reassignptr_int2d

  !> @brief Set the pointer for a 1-dimensional real array to
  !< a target array already stored in the memory manager
  subroutine reassignptr_dbl1d(adbl, name, mem_path, name_target, mem_path_target)
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: adbl !< pointer to 1d real array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    character(len=*), intent(in) :: name_target !< name of target variable
    character(len=*), intent(in) :: mem_path_target !< path where target variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    type(MemoryType), pointer :: mt2
    logical(LGP) :: found
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    call get_from_memorystore(name_target, mem_path_target, mt2, found)
    if (size(adbl) > 0) then
      nvalues_adbl = nvalues_adbl - size(adbl)
      deallocate (adbl)
    end if
    adbl => mt2%adbl1d
    mt%adbl1d => adbl
    mt%element_size = DP
    mt%isize = size(adbl)
    write (mt%memtype, "(a,' (',i0,')')") 'DOUBLE', mt%isize
    !
    ! -- set master information
    mt%master = .false.
    mt%mastername = name_target
    mt%masterPath = mem_path_target
  end subroutine reassignptr_dbl1d

  !> @brief Set the pointer for a 2-dimensional real array to
  !< a target array already stored in the memory manager
  subroutine reassignptr_dbl2d(adbl, name, mem_path, name_target, mem_path_target)
    real(DP), dimension(:, :), pointer, contiguous, intent(inout) :: adbl !< pointer to 2d real array
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    character(len=*), intent(in) :: name_target !< name of target variable
    character(len=*), intent(in) :: mem_path_target !< path where target variable is stored
    ! -- local
    type(MemoryType), pointer :: mt
    type(MemoryType), pointer :: mt2
    logical(LGP) :: found
    integer(I4B) :: ncol
    integer(I4b) :: nrow
    ! -- code
    call get_from_memorystore(name, mem_path, mt, found)
    call get_from_memorystore(name_target, mem_path_target, mt2, found)
    if (size(adbl) > 0) then
      nvalues_adbl = nvalues_adbl - size(adbl)
      deallocate (adbl)
    end if
    adbl => mt2%adbl2d
    mt%adbl2d => adbl
    mt%element_size = DP
    mt%isize = size(adbl)
    ncol = size(adbl, dim=1)
    nrow = size(adbl, dim=2)
    write (mt%memtype, "(a,' (',i0,',',i0,')')") 'DOUBLE', ncol, nrow
    !
    ! -- set master information
    mt%master = .false.
    mt%mastername = name_target
    mt%masterPath = mem_path_target
  end subroutine reassignptr_dbl2d

  !> @brief DEPRECATED. The memory manager will handle the deallocation of the pointer.
  !<
  subroutine deallocate_str(sclr, name, mem_path)
    character(len=*), pointer, intent(inout) :: sclr !< pointer to string
    character(len=*), intent(in), optional :: name !< variable name
    character(len=*), intent(in), optional :: mem_path !< path where variable is stored
    ! -- code
    return
  end subroutine deallocate_str

  !> @brief Deallocate an array of defined-length character strings
  !!
  !<
  subroutine deallocate_str1d(astr1d, name, mem_path)
    character(len=*), dimension(:), pointer, contiguous, intent(inout) :: astr1d !< array of strings
    character(len=*), optional, intent(in) :: name !< variable name
    character(len=*), optional, intent(in) :: mem_path !< path where variable is stored
    ! -- code
    return

  end subroutine deallocate_str1d

  !> @brief DEPRECATED. The memory manager will handle the deallocation of the pointer.
  !!
  !<
  subroutine deallocate_charstr1d(astr1d, name, mem_path)
    type(CharacterStringType), dimension(:), pointer, contiguous, &
      intent(inout) :: astr1d !< array of strings
    character(len=*), optional, intent(in) :: name !< variable name
    character(len=*), optional, intent(in) :: mem_path !< path where variable is stored
    ! -- code
    return
  end subroutine deallocate_charstr1d

  !> @brief DEPRECATED. The memory manager will handle the deallocation of the pointer.
  !<
  subroutine deallocate_logical(sclr)
    logical(LGP), pointer, intent(inout) :: sclr !< logical scalar to deallocate
    ! -- code
    return
  end subroutine deallocate_logical

  !> @brief DEPRECATED. The memory manager will handle the deallocation of the pointer.
  !<
  subroutine deallocate_int(sclr)
    integer(I4B), pointer, intent(inout) :: sclr !< integer variable to deallocate
    ! -- code
    return
  end subroutine deallocate_int

  !> @brief DEPRECATED. The memory manager will handle the deallocation of the pointer.
  !<
  subroutine deallocate_dbl(sclr)
    real(DP), pointer, intent(inout) :: sclr !< real variable to deallocate
    ! -- code
    return
  end subroutine deallocate_dbl

  !> @brief DEPRECATED. The memory manager will handle the deallocation of the pointer.
  !<
  subroutine deallocate_logical1d(alog, name, mem_path)
    logical(LGP), dimension(:), pointer, contiguous, intent(inout) :: alog !< 1d logical array to deallocate
    character(len=*), optional :: name !< variable name
    character(len=*), optional :: mem_path !< path where variable is stored
    ! -- code
    return
  end subroutine deallocate_logical1d

  !> @brief DEPRECATED. The memory manager will handle the deallocation of the pointer.
  !<
  subroutine deallocate_int1d(aint, name, mem_path)
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: aint !< 1d integer array to deallocate
    character(len=*), optional :: name !< variable name
    character(len=*), optional :: mem_path !< path where variable is stored
    ! -- code
    return
  end subroutine deallocate_int1d

  !> @brief DEPRECATED. The memory manager will handle the deallocation of the pointer.
  !<
  subroutine deallocate_int2d(aint, name, mem_path)
    integer(I4B), dimension(:, :), pointer, contiguous, intent(inout) :: aint !< 2d integer array to deallocate
    character(len=*), optional :: name !< variable name
    character(len=*), optional :: mem_path !< path where variable is stored
    ! -- code
    return
  end subroutine deallocate_int2d

  !> @brief DEPRECATED. The memory manager will handle the deallocation of the pointer.
  !<
  subroutine deallocate_int3d(aint, name, mem_path)
    integer(I4B), dimension(:, :, :), pointer, contiguous, intent(inout) :: aint !< 3d integer array to deallocate
    character(len=*), optional :: name !< variable name
    character(len=*), optional :: mem_path !< path where variable is stored
    ! -- code
    return
  end subroutine deallocate_int3d

  !> @brief DEPRECATED. The memory manager will handle the deallocation of the pointer.
  !<
  subroutine deallocate_dbl1d(adbl, name, mem_path)
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: adbl !< 1d real array to deallocate
    character(len=*), optional :: name !< variable name
    character(len=*), optional :: mem_path !< path where variable is stored
    ! -- code
    return
  end subroutine deallocate_dbl1d

  !> @brief DEPRECATED. The memory manager will handle the deallocation of the pointer.
  !<
  subroutine deallocate_dbl2d(adbl, name, mem_path)
    real(DP), dimension(:, :), pointer, contiguous, intent(inout) :: adbl !< 2d real array to deallocate
    character(len=*), optional :: name !< variable name
    character(len=*), optional :: mem_path !< path where variable is stored
    ! -- code
    return
  end subroutine deallocate_dbl2d

  !> @brief DEPRECATED. The memory manager will handle the deallocation of the pointer.
  !<
  subroutine deallocate_dbl3d(adbl, name, mem_path)
    real(DP), dimension(:, :, :), pointer, contiguous, intent(inout) :: adbl !< 3d real array to deallocate
    character(len=*), optional :: name !< variable name
    character(len=*), optional :: mem_path !< path where variable is stored
    ! -- code
    return
  end subroutine deallocate_dbl3d

  !> @brief Set the memory print option
  !<
  subroutine mem_set_print_option(iout, keyword, error_msg)
    integer(I4B), intent(in) :: iout !< unit number for mfsim.lst
    character(len=*), intent(in) :: keyword !< memory print option
    character(len=*), intent(inout) :: error_msg !< returned error message if keyword is not valid option
    ! -- local
    ! -- format
    ! -- code
    select case (keyword)
    case ('NONE')
      iprmem = 0
      write (iout, '(4x, a)') &
        'LIMITED MEMORY INFORMATION WILL BE WRITTEN.'
    case ('SUMMARY')
      iprmem = 1
      write (iout, '(4x, a)') &
        'A SUMMARY OF SIMULATION MEMORY INFORMATION WILL BE WRITTEN.'
    case ('ALL')
      iprmem = 2
      write (iout, '(4x, a)') &
        'ALL SIMULATION MEMORY INFORMATION WILL BE WRITTEN.'
    case default
      error_msg = "Unknown memory print option '"//trim(keyword)//"."
    end select
  end subroutine mem_set_print_option

  !> @brief Create a table if memory_print_option is 'SUMMARY'
  !<
  subroutine mem_summary_table(iout, nrows, cunits)
    integer(I4B), intent(in) :: iout !< unit number for mfsim.lst
    integer(I4B), intent(in) :: nrows !< number of table rows
    character(len=*), intent(in) :: cunits !< memory units (bytes, kilobytes, megabytes, or gigabytes)
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    integer(I4B) :: nterms
    ! -- formats
    ! -- code
    nterms = 6
    !
    ! -- set up table title
    title = 'SUMMARY INFORMATION ON VARIABLES STORED IN THE MEMORY MANAGER, '// &
            'IN '//trim(cunits)
    !
    ! -- set up stage tableobj
    call table_cr(memtab, 'MEM SUM', title)
    call memtab%table_df(nrows, nterms, iout)
    !
    ! -- data type
    text = 'COMPONENT'
    call memtab%initialize_column(text, 20, alignment=TABLEFT)
    !
    ! -- memory allocated for characters
    text = 'CHARACTER'
    call memtab%initialize_column(text, 15, alignment=TABCENTER)
    !
    ! -- memory allocated for logical
    text = 'LOGICAL'
    call memtab%initialize_column(text, 15, alignment=TABCENTER)
    !
    ! -- memory allocated for integers
    text = 'INTEGER'
    call memtab%initialize_column(text, 15, alignment=TABCENTER)
    !
    ! -- memory allocated for reals
    text = 'REAL'
    call memtab%initialize_column(text, 15, alignment=TABCENTER)
    !
    ! -- total memory allocated
    text = 'TOTAL'
    call memtab%initialize_column(text, 15, alignment=TABCENTER)
  end subroutine mem_summary_table

  !> @brief Create a table if memory_print_option is 'ALL'
  !<
  subroutine mem_detailed_table(iout, nrows)
    integer(I4B), intent(in) :: iout !< unit number for mfsim.lst
    integer(I4B), intent(in) :: nrows !< number of table rows
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    integer(I4B) :: nterms
    ! -- formats
    ! -- code
    nterms = 5
    !
    ! -- set up table title
    title = 'DETAILED INFORMATION ON VARIABLES STORED IN THE MEMORY MANAGER'
    !
    ! -- set up stage tableobj
    call table_cr(memtab, 'MEM DET', title)
    call memtab%table_df(nrows, nterms, iout)
    !
    ! -- origin
    text = 'ORIGIN'
    call memtab%initialize_column(text, LENMEMPATH, alignment=TABLEFT)
    !
    ! -- variable
    text = 'VARIABLE NAME'
    call memtab%initialize_column(text, LENVARNAME, alignment=TABLEFT)
    !
    ! -- data type
    text = 'DATA TYPE'
    call memtab%initialize_column(text, 16, alignment=TABLEFT)
    !
    ! -- size
    text = 'NUMBER OF ITEMS'
    call memtab%initialize_column(text, 20, alignment=TABRIGHT)
    !
    ! -- is it a pointer
    text = 'ASSOCIATED VARIABLE'
    call memtab%initialize_column(text, LENMEMADDRESS, alignment=TABLEFT)
  end subroutine mem_detailed_table

  !> @brief Write a row for the memory_print_option 'SUMMARY' table
  !<
  subroutine mem_summary_line(component, rchars, rlog, rint, rreal, bytes)
    character(len=*), intent(in) :: component !< character defining the program component (e.g. solution)
    real(DP), intent(in) :: rchars !< allocated size of characters (in common units)
    real(DP), intent(in) :: rlog !< allocated size of logical (in common units)
    real(DP), intent(in) :: rint !< allocated size of integer variables (in common units)
    real(DP), intent(in) :: rreal !< allocated size of real variables (in common units)
    real(DP), intent(in) :: bytes !< total allocated memory in memory manager (in common units)
    ! -- formats
    ! -- code
    !
    ! -- write line
    call memtab%add_term(component)
    call memtab%add_term(rchars)
    call memtab%add_term(rlog)
    call memtab%add_term(rint)
    call memtab%add_term(rreal)
    call memtab%add_term(bytes)
  end subroutine mem_summary_line

  !> @brief Determine appropriate memory unit and conversion factor
  !<
  subroutine mem_units(bytes, fact, cunits)
    ! -- dummy
    real(DP), intent(in) :: bytes !< total nr. of bytes
    real(DP), intent(inout) :: fact !< conversion factor
    character(len=*), intent(inout) :: cunits !< string with memory unit
    ! -- local
    ! -- formats
    ! -- code
    !
    ! -- initialize factor and unit string
    cunits = 'UNKNOWN'
    fact = DONE
    !
    ! -- set factor
    if (bytes < DEP3) then
      fact = DONE
      cunits = 'BYTES'
    else if (bytes < DEP6) then
      fact = DEM3
      cunits = 'KILOBYTES'
    else if (bytes < DEP9) then
      fact = DEM6
      cunits = 'MEGABYTES'
    else
      fact = DEM9
      cunits = 'GIGABYTES'
    end if
  end subroutine mem_units

  !> @brief Create and fill a table with the total allocated memory
  !< in the memory manager
  subroutine mem_summary_total(iout, bytes)
    integer(I4B), intent(in) :: iout !< unit number for mfsim.lst
    real(DP), intent(in) :: bytes !< total number of bytes allocated in the memory manager
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    character(LEN=10) :: cunits
    integer(I4B) :: nterms
    integer(I4B) :: nrows
    real(DP) :: fact
    real(DP) :: smb
    ! -- formats
    ! -- code
    !
    ! -- calculate factor and memory units
    call mem_units(bytes, fact, cunits)
    !
    ! -- set table terms
    nterms = 2
    nrows = 6
    !
    ! -- set up table title
    title = 'MEMORY MANAGER TOTAL STORAGE BY DATA TYPE, IN '//trim(cunits)
    !
    ! -- set up stage tableobj
    call table_cr(memtab, 'MEM TOT', title)
    call memtab%table_df(nrows, nterms, iout)
    !
    ! -- data type
    text = 'DATA TYPE'
    call memtab%initialize_column(text, 15, alignment=TABLEFT)
    !
    ! -- number of values
    text = 'ALLOCATED MEMORY'
    call memtab%initialize_column(text, 15, alignment=TABCENTER)
    !
    ! -- write data
    !
    ! -- characters
    smb = real(nvalues_astr, DP) * fact
    call memtab%add_term('Character')
    call memtab%add_term(smb)
    !
    ! -- logicals
    smb = real(nvalues_alogical * LGP, DP) * fact
    call memtab%add_term('Logical')
    call memtab%add_term(smb)
    !
    ! -- integers
    smb = real(nvalues_aint * I4B, DP) * fact
    call memtab%add_term('Integer')
    call memtab%add_term(smb)
    !
    ! -- reals
    smb = real(nvalues_adbl * DP, DP) * fact
    call memtab%add_term('Real')
    call memtab%add_term(smb)
    !
    ! -- total memory usage
    call memtab%print_separator()
    smb = bytes * fact
    call memtab%add_term('Total')
    call memtab%add_term(smb)
    !
    ! -- Virtual memory
    smb = calc_virtual_mem() * fact
    call memtab%add_term('Virtual')
    call memtab%add_term(smb)
    !
    ! -- deallocate table
    call mem_cleanup_table()
  end subroutine mem_summary_total

  !> @brief Generic function to clean a memory manager table
  !<
  subroutine mem_cleanup_table()
    ! -- local
    ! -- formats
    ! -- code
    call memtab%table_da()
    deallocate (memtab)
    nullify (memtab)
  end subroutine mem_cleanup_table

  !> @brief Write memory manager memory usage based on the
  !! user-specified memory_print_option
  !!
  !! The total memory usage by data types (int, real, etc.)
  !! is written for every simulation.
  !<
  subroutine mem_write_usage(iout)
    integer(I4B), intent(in) :: iout !< unit number for mfsim.lst
    ! -- local
    class(MemoryType), pointer :: mt
    character(len=LENMEMADDRESS), allocatable, dimension(:) :: cunique
    ! character(len=LENMEMPATH) :: mem_path
    character(len=LENMEMPATH) :: context
    character(len=LENCOMPONENTNAME) :: component
    character(len=LENCOMPONENTNAME) :: subcomponent
    character(len=LENMEMADDRESS) :: context_component
    character(LEN=10) :: cunits
    type(MemoryContainerIteratorType), allocatable :: itr
    integer(I4B) :: icomp
    integer(I4B) :: ilen
    integer(I8B) :: nchars
    integer(I8B) :: nlog
    integer(I8B) :: nint
    integer(I8B) :: nreal
    real(DP) :: simbytes
    real(DP) :: fact
    real(DP) :: rchars
    real(DP) :: rlog
    real(DP) :: rint
    real(DP) :: rreal
    real(DP) :: bytes
    ! -- formats
    ! -- code
    !
    ! -- Calculate simulation memory allocation
    simbytes = (nvalues_astr + &
                nvalues_alogical * LGP + &
                nvalues_aint * I4B + &
                nvalues_adbl * DP)
    simbytes = real(simbytes, DP)
    !
    ! -- calculate factor and memory units
    call mem_units(simbytes, fact, cunits)
    !
    ! -- Write summary table for simulation components
    if (iprmem == 1) then
      !
      ! -- Find unique names of simulation components
      call mem_unique_origins(cunique)
      call mem_summary_table(iout, size(cunique), cunits)
      do icomp = 1, size(cunique)
        nchars = 0
        nlog = 0
        nint = 0
        nreal = 0
        bytes = DZERO
        ilen = len_trim(cunique(icomp))
        itr = memorystore%iterator()
        do while (itr%has_next())
          call itr%next()
          mt => itr%value()
          call split_mem_path(mt%path, component, subcomponent)
          context = get_mem_path_context(mt%path)
          context_component = trim(context)//component
          if (cunique(icomp) /= context_component(1:ilen)) cycle
          if (.not. mt%master) cycle
          if (mt%memtype(1:6) == 'STRING') then
            nchars = nchars + mt%isize * mt%element_size
          else if (mt%memtype(1:7) == 'LOGICAL') then
            nlog = nlog + mt%isize
          else if (mt%memtype(1:7) == 'INTEGER') then
            nint = nint + mt%isize
          else if (mt%memtype(1:6) == 'DOUBLE') then
            nreal = nreal + mt%isize
          end if
        end do
        !
        ! -- calculate size of each data type in bytes
        rchars = real(nchars, DP) * fact
        rlog = real(nlog * LGP, DP) * fact
        rint = real(nint * I4B, DP) * fact
        rreal = real(nreal * DP, DP) * fact
        !
        ! -- calculate total storage in bytes
        bytes = rchars + rlog + rint + rreal
        !
        ! -- write data
        call mem_summary_line(cunique(icomp), rchars, rlog, rint, rreal, bytes)
      end do
      call mem_cleanup_table()
    end if
    !
    ! -- Write table with all variables for iprmem == 2
    if (iprmem == 2) then
      call mem_print_detailed(iout)
    end if
    !
    ! -- Write total memory allocation
    call mem_summary_total(iout, simbytes)
  end subroutine mem_write_usage

  subroutine mem_print_detailed(iout)
    integer(I4B) :: iout
    ! local
    class(MemoryType), pointer :: mt
    type(MemoryContainerIteratorType), allocatable :: itr

    call mem_detailed_table(iout, memorystore%count())
    itr = memorystore%iterator()
    do while (itr%has_next())
      call itr%next()
      mt => itr%value()
      call mt%table_entry(memtab)
    end do
    call mem_cleanup_table()

  end subroutine mem_print_detailed

  !> @brief Sum up virtual memory, i.e. memory
  !< that is owned by other processes
  function calc_virtual_mem() result(vmem_size)
    real(DP) :: vmem_size
    ! local
    type(MemoryContainerIteratorType), allocatable :: itr
    type(MemoryType), pointer :: mt

    vmem_size = DZERO
    itr = memorystore%iterator()
    do while (itr%has_next())
      call itr%next()
      mt => itr%value()
      if (index(mt%path, "__P") == 1) then
        vmem_size = mt%element_size * mt%isize + vmem_size
      end if
    end do

  end function calc_virtual_mem

  !> @brief Deallocate memory in the memory manager
  !<
  subroutine mem_da()
    ! -- modules
    use VersionModule, only: IDEVELOPMODE
    ! -- local
    class(MemoryType), pointer :: mt
    type(MemoryContainerIteratorType), allocatable :: itr
    ! -- code
    itr = memorystore%iterator()
    do while (itr%has_next())
      call itr%next()
      mt => itr%value()
      call mt%mt_deallocate()
      if (IDEVELOPMODE == 1) call mem_da_check(mt)
      deallocate (mt)
    end do

    call memorystore%clear()
    if (count_errors() > 0) then
      call store_error('Could not clear memory list.', terminate=.TRUE.)
    end if
  end subroutine mem_da

  subroutine mem_da_check(mt)
    ! -- modules
    use InputOutputModule, only: UPCASE
    ! -- dummy
    class(MemoryType), pointer :: mt
    ! -- local
    character(len=LINELENGTH) :: error_msg
    character(len=LENVARNAME) :: ucname
    !
    ! -- check if memory has been deallocated
    if (mt%mt_associated() .and. mt%element_size == -1) then
      error_msg = trim(adjustl(mt%path))//' '// &
                  trim(adjustl(mt%name))//' has invalid element size'
      call store_error(trim(error_msg))
    end if
    !
    ! -- check if memory has been deallocated
    if (mt%mt_associated() .and. mt%isize > 0) then
      error_msg = trim(adjustl(mt%path))//' '// &
                  trim(adjustl(mt%name))//' not deallocated'
      call store_error(trim(error_msg))
    end if
    !
    ! -- check case of varname
    ucname = mt%name
    call UPCASE(ucname)
    if (mt%name /= ucname) then
      error_msg = trim(adjustl(mt%path))//' '// &
                  trim(adjustl(mt%name))//' not upper case'
      call store_error(trim(error_msg))
    end if
  end subroutine mem_da_check

  !> @brief Create a array with unique first components from all memory paths.
  !! Only the first component of the memory path is evaluated.
  !<
  subroutine mem_unique_origins(cunique)
    ! -- modules
    use ArrayHandlersModule, only: ExpandArray, ifind
    ! -- dummy
    character(len=LENMEMADDRESS), allocatable, dimension(:), intent(inout) :: &
      cunique !< array with unique first components
    ! -- local
    class(MemoryType), pointer :: mt
    character(len=LENMEMPATH) :: context
    character(len=LENCOMPONENTNAME) :: component
    character(len=LENCOMPONENTNAME) :: subcomponent
    character(len=LENMEMADDRESS) :: context_component
    type(MemoryContainerIteratorType), allocatable :: itr
    integer(I4B) :: ipa
    ! -- code
    !
    ! -- initialize cunique
    allocate (cunique(0))
    !
    ! -- find unique origins
    itr = memorystore%iterator()
    do while (itr%has_next())
      call itr%next()
      mt => itr%value()
      call split_mem_path(mt%path, component, subcomponent)
      context = get_mem_path_context(mt%path)
      context_component = trim(context)//component
      ipa = ifind(cunique, context_component)
      if (ipa < 1) then
        call ExpandArray(cunique, 1)
        cunique(size(cunique)) = context_component
      end if
    end do
  end subroutine mem_unique_origins

end module MemoryManagerModule
