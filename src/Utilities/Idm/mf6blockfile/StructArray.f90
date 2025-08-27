!> @brief This module contains the StructArrayModule
!!
!! This module contains the routines for reading a
!! structured list, which consists of a separate vector
!! for each column in the list.
!!
!<
module StructArrayModule

  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: DZERO, IZERO, DNODATA, &
                             LINELENGTH, LENMEMPATH, LENVARNAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use StructVectorModule, only: StructVectorType
  use InputDefinitionModule, only: InputParamDefinitionType
  use MemoryManagerModule, only: mem_allocate, mem_reallocate, mem_setptr
  use CharacterStringModule, only: CharacterStringType
  use STLVecIntModule, only: STLVecInt
  use IdmLoggerModule, only: idm_log_var
  use BlockParserModule, only: BlockParserType
  use ModflowInputModule, only: ModflowInputType

  implicit none
  private
  public :: StructArrayType
  public :: constructStructArray, destructStructArray

  !> @brief type for structured array
  !!
  !! This type is used to read and store a list
  !! that consists of multiple one-dimensional
  !! vectors.
  !!
  !<
  type StructArrayType
    integer(I4B) :: ncol
    integer(I4B) :: nrow
    integer(I4B) :: blocknum
    logical(LGP) :: deferred_shape = .false.
    integer(I4B) :: deferred_size_init = 5
    character(len=LENMEMPATH) :: mempath
    character(len=LENMEMPATH) :: component_mempath
    type(StructVectorType), dimension(:), allocatable :: struct_vectors
    integer(I4B), dimension(:), allocatable :: startidx
    integer(I4B), dimension(:), allocatable :: numcols
    type(ModflowInputType) :: mf6_input
  contains
    procedure :: mem_create_vector
    procedure :: count
    procedure :: get
    procedure :: allocate_int_type
    procedure :: allocate_dbl_type
    procedure :: allocate_charstr_type
    procedure :: allocate_int1d_type
    procedure :: allocate_dbl1d_type
    procedure :: read_param
    procedure :: read_from_parser
    procedure :: read_from_binary
    procedure :: memload_vectors
    procedure :: load_deferred_vector
    procedure :: log_structarray_vars
    procedure :: check_reallocate

  end type StructArrayType

contains

  !> @brief constructor for a struct_array
  !<
  function constructStructArray(mf6_input, ncol, nrow, blocknum, mempath, &
                                component_mempath) result(struct_array)
    type(ModflowInputType), intent(in) :: mf6_input
    integer(I4B), intent(in) :: ncol !< number of columns in the StructArrayType
    integer(I4B), intent(in) :: nrow !< number of rows in the StructArrayType
    integer(I4B), intent(in) :: blocknum !< valid block number or 0
    character(len=*), intent(in) :: mempath !< memory path for storing the vector
    character(len=*), intent(in) :: component_mempath
    type(StructArrayType), pointer :: struct_array !< new StructArrayType

    ! allocate StructArrayType
    allocate (struct_array)

    ! set description of input
    struct_array%mf6_input = mf6_input

    ! set number of arrays
    struct_array%ncol = ncol

    ! set rows if known or set deferred
    struct_array%nrow = nrow
    if (struct_array%nrow == -1) then
      struct_array%nrow = 0
      struct_array%deferred_shape = .true.
    end if

    ! set blocknum
    if (blocknum > 0) then
      struct_array%blocknum = blocknum
    else
      struct_array%blocknum = 0
    end if

    ! set mempath
    struct_array%mempath = mempath
    struct_array%component_mempath = component_mempath

    ! allocate StructVectorType objects
    allocate (struct_array%struct_vectors(ncol))
    allocate (struct_array%startidx(ncol))
    allocate (struct_array%numcols(ncol))
  end function constructStructArray

  !> @brief destructor for a struct_array
  !<
  subroutine destructStructArray(struct_array)
    type(StructArrayType), pointer, intent(inout) :: struct_array !< StructArrayType to destroy
    deallocate (struct_array%struct_vectors)
    deallocate (struct_array%startidx)
    deallocate (struct_array%numcols)
    deallocate (struct_array)
    nullify (struct_array)
  end subroutine destructStructArray

  !> @brief create new vector in StructArrayType
  !<
  subroutine mem_create_vector(this, icol, idt)
    class(StructArrayType) :: this !< StructArrayType
    integer(I4B), intent(in) :: icol !< column to create
    type(InputParamDefinitionType), pointer :: idt
    type(StructVectorType) :: sv
    integer(I4B) :: numcol

    ! initialize
    numcol = 1
    sv%idt => idt
    sv%icol = icol

    ! set size
    if (this%deferred_shape) then
      sv%size = this%deferred_size_init
    else
      sv%size = this%nrow
    end if

    ! allocate array memory for StructVectorType
    select case (idt%datatype)
    case ('INTEGER')
      call this%allocate_int_type(sv)
    case ('DOUBLE')
      call this%allocate_dbl_type(sv)
    case ('STRING', 'KEYWORD')
      call this%allocate_charstr_type(sv)
    case ('INTEGER1D')
      call this%allocate_int1d_type(sv)
      if (sv%memtype == 5) then
        numcol = sv%intshape
      end if
    case ('DOUBLE1D')
      call this%allocate_dbl1d_type(sv)
      numcol = sv%intshape
    case default
      errmsg = 'IDM unimplemented. StructArray::mem_create_vector &
               &type='//trim(idt%datatype)
      call store_error(errmsg, .true.)
    end select

    ! set the object in the Struct Array
    this%struct_vectors(icol) = sv
    this%numcols(icol) = numcol
    if (icol == 1) then
      this%startidx(icol) = 1
    else
      this%startidx(icol) = this%startidx(icol - 1) + this%numcols(icol - 1)
    end if
  end subroutine mem_create_vector

  function count(this)
    class(StructArrayType) :: this !< StructArrayType
    integer(I4B) :: count
    count = size(this%struct_vectors)
  end function count

  subroutine set_pointer(sv, sv_target)
    type(StructVectorType), pointer :: sv
    type(StructVectorType), target :: sv_target
    sv => sv_target
  end subroutine set_pointer

  function get(this, idx) result(sv)
    class(StructArrayType) :: this !< StructArrayType
    integer(I4B), intent(in) :: idx
    type(StructVectorType), pointer :: sv
    call set_pointer(sv, this%struct_vectors(idx))
  end function get

  !> @brief allocate integer input type
  !<
  subroutine allocate_int_type(this, sv)
    class(StructArrayType) :: this !< StructArrayType
    type(StructVectorType), intent(inout) :: sv
    integer(I4B), dimension(:), pointer, contiguous :: int1d
    integer(I4B) :: j, nrow

    if (this%deferred_shape) then
      ! shape not known, allocate locally
      nrow = this%deferred_size_init
      allocate (int1d(this%deferred_size_init))
    else
      ! shape known, allocate in managed memory
      nrow = this%nrow
      call mem_allocate(int1d, this%nrow, sv%idt%mf6varname, this%mempath)
    end if

    ! initialize vector values
    do j = 1, nrow
      int1d(j) = IZERO
    end do

    sv%memtype = 1
    sv%int1d => int1d
  end subroutine allocate_int_type

  !> @brief allocate double input type
  !<
  subroutine allocate_dbl_type(this, sv)
    class(StructArrayType) :: this !< StructArrayType
    type(StructVectorType), intent(inout) :: sv
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    integer(I4B) :: j, nrow

    if (this%deferred_shape) then
      ! shape not known, allocate locally
      nrow = this%deferred_size_init
      allocate (dbl1d(this%deferred_size_init))
    else
      ! shape known, allocate in managed memory
      nrow = this%nrow
      call mem_allocate(dbl1d, this%nrow, sv%idt%mf6varname, this%mempath)
    end if

    ! initialize
    do j = 1, nrow
      dbl1d(j) = DZERO
    end do

    sv%memtype = 2
    sv%dbl1d => dbl1d
  end subroutine allocate_dbl_type

  !> @brief allocate charstr input type
  !<
  subroutine allocate_charstr_type(this, sv)
    class(StructArrayType) :: this !< StructArrayType
    type(StructVectorType), intent(inout) :: sv
    type(CharacterStringType), dimension(:), pointer, contiguous :: charstr1d
    integer(I4B) :: j

    if (this%deferred_shape) then
      allocate (charstr1d(this%deferred_size_init))
    else
      call mem_allocate(charstr1d, LINELENGTH, this%nrow, &
                        sv%idt%mf6varname, this%mempath)
    end if

    do j = 1, this%nrow
      charstr1d(j) = ''
    end do

    sv%memtype = 3
    sv%charstr1d => charstr1d
  end subroutine allocate_charstr_type

  !> @brief allocate int1d input type
  !<
  subroutine allocate_int1d_type(this, sv)
    use ConstantsModule, only: LENMODELNAME
    use MemoryHelperModule, only: create_mem_path
    use SimVariablesModule, only: idm_context
    class(StructArrayType) :: this !< StructArrayType
    type(StructVectorType), intent(inout) :: sv
    integer(I4B), dimension(:, :), pointer, contiguous :: int2d
    type(STLVecInt), pointer :: intvector
    integer(I4B), pointer :: ncelldim, exgid
    character(len=LENMEMPATH) :: input_mempath
    character(len=LENMODELNAME) :: mname
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: charstr1d
    integer(I4B) :: nrow, n, m

    if (sv%idt%shape == 'NCELLDIM') then
      ! if EXCHANGE set to NCELLDIM of appropriate model
      if (this%mf6_input%component_type == 'EXG') then
        ! set pointer to EXGID
        call mem_setptr(exgid, 'EXGID', this%mf6_input%mempath)
        ! set pointer to appropriate exchange model array
        input_mempath = create_mem_path('SIM', 'NAM', idm_context)
        if (sv%idt%tagname == 'CELLIDM1') then
          call mem_setptr(charstr1d, 'EXGMNAMEA', input_mempath)
        else if (sv%idt%tagname == 'CELLIDM2') then
          call mem_setptr(charstr1d, 'EXGMNAMEB', input_mempath)
        end if

        ! set the model name
        mname = charstr1d(exgid)

        ! set ncelldim pointer
        input_mempath = create_mem_path(component=mname, context=idm_context)
        call mem_setptr(ncelldim, sv%idt%shape, input_mempath)
      else
        call mem_setptr(ncelldim, sv%idt%shape, this%component_mempath)
      end if

      if (this%deferred_shape) then
        ! shape not known, allocate locally
        nrow = this%deferred_size_init
        allocate (int2d(ncelldim, this%deferred_size_init))
      else
        ! shape known, allocate in managed memory
        nrow = this%nrow
        call mem_allocate(int2d, ncelldim, this%nrow, &
                          sv%idt%mf6varname, this%mempath)
      end if

      ! initialize
      do m = 1, nrow
        do n = 1, ncelldim
          int2d(n, m) = IZERO
        end do
      end do

      sv%memtype = 5
      sv%int2d => int2d
      sv%intshape => ncelldim
    else
      ! allocate intvector object
      allocate (intvector)
      ! initialize STLVecInt
      call intvector%init()
      sv%memtype = 4
      sv%intvector => intvector
      sv%size = -1
      ! set pointer to dynamic shape
      call mem_setptr(sv%intvector_shape, sv%idt%shape, this%mempath)
    end if
  end subroutine allocate_int1d_type

  !> @brief allocate dbl1d input type
  !<
  subroutine allocate_dbl1d_type(this, sv)
    use MemoryManagerModule, only: get_isize
    class(StructArrayType) :: this !< StructArrayType
    type(StructVectorType), intent(inout) :: sv
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B), pointer :: naux, nseg, nseg_1
    integer(I4B) :: nseg1_isize, n, m

    if (sv%idt%shape == 'NAUX') then
      call mem_setptr(naux, sv%idt%shape, this%mempath)
      call mem_allocate(dbl2d, naux, this%nrow, sv%idt%mf6varname, this%mempath)

      ! initialize
      do m = 1, this%nrow
        do n = 1, naux
          dbl2d(n, m) = DZERO
        end do
      end do

      sv%memtype = 6
      sv%dbl2d => dbl2d
      sv%intshape => naux
    else if (sv%idt%shape == 'NSEG-1') then
      call mem_setptr(nseg, 'NSEG', this%mempath)
      call get_isize('NSEG_1', this%mempath, nseg1_isize)

      if (nseg1_isize < 0) then
        call mem_allocate(nseg_1, 'NSEG_1', this%mempath)
        nseg_1 = nseg - 1
      else
        call mem_setptr(nseg_1, 'NSEG_1', this%mempath)
      end if

      ! allocate
      call mem_allocate(dbl2d, nseg_1, this%nrow, sv%idt%mf6varname, this%mempath)

      ! initialize
      do m = 1, this%nrow
        do n = 1, nseg_1
          dbl2d(n, m) = DZERO
        end do
      end do

      sv%memtype = 6
      sv%dbl2d => dbl2d
      sv%intshape => nseg_1
    else
      errmsg = 'IDM unimplemented. StructArray::allocate_dbl1d_type &
               & unsupported shape "'//trim(sv%idt%shape)//'".'
      call store_error(errmsg, terminate=.TRUE.)
    end if
  end subroutine allocate_dbl1d_type

  subroutine load_deferred_vector(this, icol)
    use MemoryManagerModule, only: get_isize
    class(StructArrayType) :: this !< StructArrayType
    integer(I4B), intent(in) :: icol
    integer(I4B) :: i, j, isize
    integer(I4B), dimension(:), pointer, contiguous :: p_int1d
    integer(I4B), dimension(:, :), pointer, contiguous :: p_int2d
    real(DP), dimension(:), pointer, contiguous :: p_dbl1d
    type(CharacterStringType), dimension(:), pointer, contiguous :: p_charstr1d
    character(len=LENVARNAME) :: varname
    logical(LGP) :: overwrite

    overwrite = .true.
    if (this%struct_vectors(icol)%idt%blockname == 'SOLUTIONGROUP') &
      overwrite = .false.

    ! set varname
    varname = this%struct_vectors(icol)%idt%mf6varname
    ! check if already mem managed variable
    call get_isize(varname, this%mempath, isize)

    ! allocate and load based on memtype
    select case (this%struct_vectors(icol)%memtype)
    case (1) ! memtype integer
      if (isize > -1) then
        ! variable exists, reallocate and append
        call mem_setptr(p_int1d, varname, this%mempath)

        if (overwrite) then
          ! overwrite existing array
          if (this%nrow > isize) then
            ! reallocate
            call mem_reallocate(p_int1d, this%nrow, varname, this%mempath)
          end if

          ! write new data
          do i = 1, this%nrow
            p_int1d(i) = this%struct_vectors(icol)%int1d(i)
          end do

          if (isize > this%nrow) then
            ! initialize excess space
            do i = this%nrow + 1, isize
              p_int1d(i) = IZERO
            end do
          end if
        else
          ! reallocate to new size
          call mem_reallocate(p_int1d, this%nrow + isize, varname, this%mempath)

          ! write new data after existing
          do i = 1, this%nrow
            p_int1d(isize + i) = this%struct_vectors(icol)%int1d(i)
          end do
        end if
      else
        ! allocate memory manager vector
        call mem_allocate(p_int1d, this%nrow, varname, this%mempath)

        ! load local vector to managed memory
        do i = 1, this%nrow
          p_int1d(i) = this%struct_vectors(icol)%int1d(i)
        end do
      end if

      ! deallocate local memory
      deallocate (this%struct_vectors(icol)%int1d)

      ! update structvector
      this%struct_vectors(icol)%int1d => p_int1d
      this%struct_vectors(icol)%size = this%nrow
    case (2) ! memtype real
      if (isize > -1) then
        call mem_setptr(p_dbl1d, varname, this%mempath)

        if (overwrite) then
          if (this%nrow > isize) then
            call mem_reallocate(p_dbl1d, this%nrow, varname, this%mempath)
          end if

          do i = 1, this%nrow
            p_dbl1d(i) = this%struct_vectors(icol)%dbl1d(i)
          end do

          if (isize > this%nrow) then
            do i = this%nrow + 1, isize
              p_dbl1d(i) = DZERO
            end do
          end if
        else
          call mem_reallocate(p_dbl1d, this%nrow + isize, varname, &
                              this%mempath)
          do i = 1, this%nrow
            p_dbl1d(isize + i) = this%struct_vectors(icol)%dbl1d(i)
          end do
        end if
      else
        call mem_allocate(p_dbl1d, this%nrow, varname, this%mempath)

        do i = 1, this%nrow
          p_dbl1d(i) = this%struct_vectors(icol)%dbl1d(i)
        end do
      end if

      deallocate (this%struct_vectors(icol)%dbl1d)

      this%struct_vectors(icol)%dbl1d => p_dbl1d
      this%struct_vectors(icol)%size = this%nrow
      !
    case (3) ! memtype charstring
      if (isize > -1) then
        call mem_setptr(p_charstr1d, varname, this%mempath)

        if (overwrite) then
          if (this%nrow > isize) then
            call mem_reallocate(p_charstr1d, LINELENGTH, this%nrow, varname, &
                                this%mempath)
          end if

          do i = 1, this%nrow
            p_charstr1d(i) = this%struct_vectors(icol)%charstr1d(i)
          end do

          if (isize > this%nrow) then
            do i = this%nrow + 1, isize
              p_charstr1d(i) = ''
            end do
          end if
        else
          call mem_reallocate(p_charstr1d, LINELENGTH, this%nrow + isize, &
                              varname, this%mempath)
          do i = 1, this%nrow
            p_charstr1d(isize + i) = this%struct_vectors(icol)%charstr1d(i)
          end do
        end if
      else
        call mem_allocate(p_charstr1d, LINELENGTH, this%nrow, varname, &
                          this%mempath)
        do i = 1, this%nrow
          p_charstr1d(i) = this%struct_vectors(icol)%charstr1d(i)
          call this%struct_vectors(icol)%charstr1d(i)%destroy()
        end do
      end if

      deallocate (this%struct_vectors(icol)%charstr1d)

      this%struct_vectors(icol)%charstr1d => p_charstr1d
      this%struct_vectors(icol)%size = this%nrow
    case (4) ! memtype intvector
      errmsg = 'StructArray::load_deferred_vector &
               &intvector reallocate unimplemented.'
      call store_error(errmsg, terminate=.TRUE.)
    case (5)
      if (isize > -1) then
        errmsg = 'StructArray::load_deferred_vector &
                 &int2d reallocate unimplemented.'
        call store_error(errmsg, terminate=.TRUE.)
      else
        call mem_allocate(p_int2d, this%struct_vectors(icol)%intshape, &
                          this%nrow, varname, this%mempath)
        do i = 1, this%nrow
          do j = 1, this%struct_vectors(icol)%intshape
            p_int2d(j, i) = this%struct_vectors(icol)%int2d(j, i)
          end do
        end do
      end if

      deallocate (this%struct_vectors(icol)%int2d)

      this%struct_vectors(icol)%int2d => p_int2d
      this%struct_vectors(icol)%size = this%nrow
    case (6)
      errmsg = 'StructArray::load_deferred_vector &
               &dbl2d reallocate unimplemented.'
      call store_error(errmsg, terminate=.TRUE.)
    case default
    end select
  end subroutine load_deferred_vector

  !> @brief load deferred vectors into managed memory
  !<
  subroutine memload_vectors(this)
    class(StructArrayType) :: this !< StructArrayType
    integer(I4B) :: icol, j
    integer(I4B), dimension(:), pointer, contiguous :: p_intvector
    character(len=LENVARNAME) :: varname

    do icol = 1, this%ncol
      ! set varname
      varname = this%struct_vectors(icol)%idt%mf6varname

      if (this%struct_vectors(icol)%memtype == 4) then
        ! intvectors always need to be loaded
        ! size intvector to number of values read
        call this%struct_vectors(icol)%intvector%shrink_to_fit()

        ! allocate memory manager vector
        call mem_allocate(p_intvector, &
                          this%struct_vectors(icol)%intvector%size, &
                          varname, this%mempath)

        ! load local vector to managed memory
        do j = 1, this%struct_vectors(icol)%intvector%size
          p_intvector(j) = this%struct_vectors(icol)%intvector%at(j)
        end do

        ! cleanup local memory
        call this%struct_vectors(icol)%intvector%destroy()
        deallocate (this%struct_vectors(icol)%intvector)
        nullify (this%struct_vectors(icol)%intvector_shape)
      else if (this%deferred_shape) then
        ! load as shape wasn't known
        call this%load_deferred_vector(icol)
      end if
    end do
  end subroutine memload_vectors

  !> @brief log information about the StructArrayType
  !<
  subroutine log_structarray_vars(this, iout)
    class(StructArrayType) :: this !< StructArrayType
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B) :: j
    integer(I4B), dimension(:), pointer, contiguous :: int1d

    ! idm variable logging
    do j = 1, this%ncol
      ! log based on memtype
      select case (this%struct_vectors(j)%memtype)
      case (1) ! memtype integer
        call idm_log_var(this%struct_vectors(j)%int1d, &
                         this%struct_vectors(j)%idt%tagname, &
                         this%mempath, iout)
      case (2) ! memtype real
        if (this%struct_vectors(j)%ts_strlocs%count() > 0) then
          call idm_log_var(this%struct_vectors(j)%idt%tagname, &
                           this%mempath, iout, .false.)
        else
          call idm_log_var(this%struct_vectors(j)%dbl1d, &
                           this%struct_vectors(j)%idt%tagname, &
                           this%mempath, iout)
        end if
      case (4) ! memtype intvector
        call mem_setptr(int1d, this%struct_vectors(j)%idt%mf6varname, &
                        this%mempath)
        call idm_log_var(int1d, this%struct_vectors(j)%idt%tagname, &
                         this%mempath, iout)
      case (5) ! memtype int2d
        call idm_log_var(this%struct_vectors(j)%int2d, &
                         this%struct_vectors(j)%idt%tagname, &
                         this%mempath, iout)
      case (6) ! memtype dbl2d
        if (this%struct_vectors(j)%ts_strlocs%count() > 0) then
          call idm_log_var(this%struct_vectors(j)%idt%tagname, &
                           this%mempath, iout, .false.)
        else
          call idm_log_var(this%struct_vectors(j)%dbl2d, &
                           this%struct_vectors(j)%idt%tagname, &
                           this%mempath, iout)
        end if
      end select
    end do
  end subroutine log_structarray_vars

  !> @brief reallocate local memory for deferred vectors if necessary
  !<
  subroutine check_reallocate(this)
    class(StructArrayType) :: this !< StructArrayType
    integer(I4B) :: i, j, k, newsize
    integer(I4B), dimension(:), pointer, contiguous :: p_int1d
    integer(I4B), dimension(:, :), pointer, contiguous :: p_int2d
    real(DP), dimension(:), pointer, contiguous :: p_dbl1d
    type(CharacterStringType), dimension(:), pointer, contiguous :: p_charstr1d
    integer(I4B) :: reallocate_mult

    ! set growth rate
    reallocate_mult = 2

    do j = 1, this%ncol
      ! reallocate based on memtype
      select case (this%struct_vectors(j)%memtype)
      case (1) ! memtype integer
        ! check if more space needed
        if (this%nrow > this%struct_vectors(j)%size) then
          ! calculate new size
          newsize = this%struct_vectors(j)%size * reallocate_mult
          ! allocate new vector
          allocate (p_int1d(newsize))

          ! copy from old to new
          do i = 1, this%struct_vectors(j)%size
            p_int1d(i) = this%struct_vectors(j)%int1d(i)
          end do

          ! deallocate old vector
          deallocate (this%struct_vectors(j)%int1d)

          ! update struct array object
          this%struct_vectors(j)%int1d => p_int1d
          this%struct_vectors(j)%size = newsize
        end if
      case (2) ! memtype real
        if (this%nrow > this%struct_vectors(j)%size) then
          newsize = this%struct_vectors(j)%size * reallocate_mult
          allocate (p_dbl1d(newsize))

          do i = 1, this%struct_vectors(j)%size
            p_dbl1d(i) = this%struct_vectors(j)%dbl1d(i)
          end do

          deallocate (this%struct_vectors(j)%dbl1d)

          this%struct_vectors(j)%dbl1d => p_dbl1d
          this%struct_vectors(j)%size = newsize
        end if
        !
      case (3) ! memtype charstring
        if (this%nrow > this%struct_vectors(j)%size) then
          newsize = this%struct_vectors(j)%size * reallocate_mult
          allocate (p_charstr1d(newsize))

          do i = 1, this%struct_vectors(j)%size
            p_charstr1d(i) = this%struct_vectors(j)%charstr1d(i)
            call this%struct_vectors(j)%charstr1d(i)%destroy()
          end do

          deallocate (this%struct_vectors(j)%charstr1d)

          this%struct_vectors(j)%charstr1d => p_charstr1d
          this%struct_vectors(j)%size = newsize
        end if
      case (5)
        if (this%nrow > this%struct_vectors(j)%size) then
          newsize = this%struct_vectors(j)%size * reallocate_mult
          allocate (p_int2d(this%struct_vectors(j)%intshape, newsize))

          do i = 1, this%struct_vectors(j)%size
            do k = 1, this%struct_vectors(j)%intshape
              p_int2d(k, i) = this%struct_vectors(j)%int2d(k, i)
            end do
          end do

          deallocate (this%struct_vectors(j)%int2d)

          this%struct_vectors(j)%int2d => p_int2d
          this%struct_vectors(j)%size = newsize
        end if
        ! TODO: case (6)
      case default
        errmsg = 'IDM unimplemented. StructArray::check_reallocate &
                 &unsupported memtype.'
        call store_error(errmsg, terminate=.TRUE.)
      end select
    end do
  end subroutine check_reallocate

  subroutine read_param(this, parser, sv_col, irow, timeseries, iout, auxcol)
    class(StructArrayType) :: this !< StructArrayType
    type(BlockParserType), intent(inout) :: parser !< block parser to read from
    integer(I4B), intent(in) :: sv_col
    integer(I4B), intent(in) :: irow
    logical(LGP), intent(in) :: timeseries
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B), optional, intent(in) :: auxcol
    integer(I4B) :: n, intval, numval, icol
    character(len=LINELENGTH) :: str
    character(len=:), allocatable :: line
    logical(LGP) :: preserve_case

    select case (this%struct_vectors(sv_col)%memtype)
    case (1) ! memtype integer
      ! if reloadable block and first col, store blocknum
      if (sv_col == 1 .and. this%blocknum > 0) then
        ! store blocknum
        this%struct_vectors(sv_col)%int1d(irow) = this%blocknum
      else
        ! read and store int
        this%struct_vectors(sv_col)%int1d(irow) = parser%GetInteger()
      end if
    case (2) ! memtype real
      if (this%struct_vectors(sv_col)%idt%timeseries .and. timeseries) then
        call parser%GetString(str)
        if (present(auxcol)) then
          icol = auxcol
        else
          icol = 1
        end if
        this%struct_vectors(sv_col)%dbl1d(irow) = &
          this%struct_vectors(sv_col)%read_token(str, this%startidx(sv_col), &
                                                 icol, irow)
      else
        this%struct_vectors(sv_col)%dbl1d(irow) = parser%GetDouble()
      end if
    case (3) ! memtype charstring
      if (this%struct_vectors(sv_col)%idt%shape /= '') then
        ! if last column with any shape, store rest of line
        if (sv_col == this%ncol) then
          call parser%GetRemainingLine(line)
          this%struct_vectors(sv_col)%charstr1d(irow) = line
          deallocate (line)
        end if
      else
        ! read string token
        preserve_case = (.not. this%struct_vectors(sv_col)%idt%preserve_case)
        call parser%GetString(str, preserve_case)
        this%struct_vectors(sv_col)%charstr1d(irow) = str
      end if
    case (4) ! memtype intvector
      ! get shape for this row
      numval = this%struct_vectors(sv_col)%intvector_shape(irow)
      ! read and store row values
      do n = 1, numval
        intval = parser%GetInteger()
        call this%struct_vectors(sv_col)%intvector%push_back(intval)
      end do
    case (5) ! memtype int2d
      ! read and store row values
      do n = 1, this%struct_vectors(sv_col)%intshape
        this%struct_vectors(sv_col)%int2d(n, irow) = parser%GetInteger()
      end do
    case (6) ! memtype dbl2d
      ! read and store row values
      do n = 1, this%struct_vectors(sv_col)%intshape
        if (this%struct_vectors(sv_col)%idt%timeseries .and. timeseries) then
          call parser%GetString(str)
          icol = this%startidx(sv_col) + n - 1
          this%struct_vectors(sv_col)%dbl2d(n, irow) = &
            this%struct_vectors(sv_col)%read_token(str, icol, n, irow)
        else
          this%struct_vectors(sv_col)%dbl2d(n, irow) = parser%GetDouble()
        end if
      end do
    end select
  end subroutine read_param

  !> @brief read from the block parser to fill the StructArrayType
  !<
  function read_from_parser(this, parser, timeseries, iout) result(irow)
    class(StructArrayType) :: this !< StructArrayType
    type(BlockParserType) :: parser !< block parser to read from
    logical(LGP), intent(in) :: timeseries
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B) :: irow, j
    logical(LGP) :: endOfBlock

    ! initialize index irow
    irow = 0

    ! read entire block
    do
      ! read next line
      call parser%GetNextLine(endOfBlock)
      if (endOfBlock) then
        ! no more lines
        exit
      else if (this%deferred_shape) then
        ! shape unknown, track lines read
        this%nrow = this%nrow + 1
        ! check and update memory allocation
        call this%check_reallocate()
      end if
      ! update irow index
      irow = irow + 1
      ! handle line reads by column memtype
      do j = 1, this%ncol
        call this%read_param(parser, j, irow, timeseries, iout)
      end do
    end do
    ! if deferred shape vectors were read, load to input path
    call this%memload_vectors()
    ! log loaded variables
    if (iout > 0) then
      call this%log_structarray_vars(iout)
    end if
  end function read_from_parser

  !> @brief read from binary input to fill the StructArrayType
  !<
  function read_from_binary(this, inunit, iout) result(irow)
    class(StructArrayType) :: this !< StructArrayType
    integer(I4B), intent(in) :: inunit !< unit number for binary input
    integer(I4B), intent(in) :: iout !< unit number for output
    integer(I4B) :: irow, ierr
    integer(I4B) :: j, k
    integer(I4B) :: intval, numval
    character(len=LINELENGTH) :: fname
    character(len=*), parameter :: fmtlsterronly = &
      "('Error reading LIST from file: ',&
      &1x,a,1x,' on UNIT: ',I0)"

    ! set error and exit if deferred shape
    if (this%deferred_shape) then
      errmsg = 'IDM unimplemented. StructArray::read_from_binary deferred shape &
               &not supported for binary inputs.'
      call store_error(errmsg, terminate=.TRUE.)
    end if
    ! initialize
    irow = 0
    ierr = 0
    readloop: do
      ! update irow index
      irow = irow + 1
      ! handle line reads by column memtype
      do j = 1, this%ncol
        select case (this%struct_vectors(j)%memtype)
        case (1) ! memtype integer
          read (inunit, iostat=ierr) this%struct_vectors(j)%int1d(irow)
        case (2) ! memtype real
          read (inunit, iostat=ierr) this%struct_vectors(j)%dbl1d(irow)
        case (3) ! memtype charstring
          errmsg = 'List style binary inputs not supported &
                   &for text columns, tag='// &
                   trim(this%struct_vectors(j)%idt%tagname)//'.'
          call store_error(errmsg, terminate=.TRUE.)
        case (4) ! memtype intvector
          ! get shape for this row
          numval = this%struct_vectors(j)%intvector_shape(irow)
          ! read and store row values
          do k = 1, numval
            if (ierr == 0) then
              read (inunit, iostat=ierr) intval
              call this%struct_vectors(j)%intvector%push_back(intval)
            end if
          end do
        case (5) ! memtype int2d
          ! read and store row values
          do k = 1, this%struct_vectors(j)%intshape
            if (ierr == 0) then
              read (inunit, iostat=ierr) this%struct_vectors(j)%int2d(k, irow)
            end if
          end do
        case (6) ! memtype dbl2d
          do k = 1, this%struct_vectors(j)%intshape
            if (ierr == 0) then
              read (inunit, iostat=ierr) this%struct_vectors(j)%dbl2d(k, irow)
            end if
          end do
        end select

        ! handle error cases
        select case (ierr)
        case (0)
          ! no error
        case (:-1)
          ! End of block was encountered
          irow = irow - 1
          exit readloop
        case (1:)
          ! Error
          inquire (unit=inunit, name=fname)
          write (errmsg, fmtlsterronly) trim(adjustl(fname)), inunit
          call store_error(errmsg, terminate=.TRUE.)
        case default
        end select
      end do
      if (irow == this%nrow) exit readloop
    end do readloop

    ! Stop if errors were detected
    !if (count_errors() > 0) then
    !  call store_error_unit(inunit)
    !end if

    ! if deferred shape vectors were read, load to input path
    call this%memload_vectors()

    ! log loaded variables
    if (iout > 0) then
      call this%log_structarray_vars(iout)
    end if
  end function read_from_binary

end module StructArrayModule
