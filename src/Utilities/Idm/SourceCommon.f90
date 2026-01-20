!> @brief This module contains the SourceCommonModule
!!
!! This module contains source independent input
!! processing helper routines.
!!
!<
module SourceCommonModule

  use KindModule, only: DP, I4B, LGP
  use SimVariablesModule, only: errmsg
  use ConstantsModule, only: LINELENGTH, LENPACKAGETYPE, LENPACKAGENAME, &
                             LENCOMPONENTNAME
  use SimModule, only: store_error, store_error_filename

  implicit none
  private
  public :: package_source_type
  public :: idm_component_type, idm_subcomponent_type, idm_subcomponent_name
  public :: set_model_shape
  public :: get_shape_from_string
  public :: get_layered_shape
  public :: file_ext
  public :: ifind_charstr
  public :: filein_fname
  public :: inlen_check

contains

  !> @brief source identifier from model namfile FNAME array
  !!
  !! Return the source type for a package listed in the
  !! model nam file packages block FNAME field.
  !!
  !<
  function package_source_type(sourcename) result(sourcetype)
    use InputOutputModule, only: upcase
    character(len=*), intent(in) :: sourcename
    character(len=LENPACKAGENAME) :: sourcetype
    character(len=LENPACKAGENAME) :: ext
    ext = file_ext(sourcename)
    select case (ext)
    case default
      sourcetype = 'MF6FILE'
    end select
  end function package_source_type

  !> @brief component from package or model type
  !!
  !! Return the component type typically derived from package file type,
  !! i.e. return GWF when input is GWF6. This function checks the
  !! resultant component type and throws a terminating error if not
  !! supported by IDM in some capacity.
  !!
  !<
  function idm_component_type(component) result(component_type)
    use IdmDfnSelectorModule, only: idm_component
    character(len=*), intent(in) :: component
    character(len=LENCOMPONENTNAME) :: component_type
    integer(I4B) :: i, ilen, idx

    ! initialize
    component_type = ''
    idx = 0

    ilen = len_trim(component)
    do i = 1, ilen
      if (component(i:i) == '6' .or. component(i:i) == '-') then
      else
        idx = idx + 1
        component_type(idx:idx) = component(i:i)
      end if
    end do

    if (.not. idm_component(component_type)) then
      write (errmsg, '(a)') &
        'IDP input error, unrecognized component: "'//trim(component)//'"'
      call store_error(errmsg, .true.)
    end if
  end function idm_component_type

  !> @brief component from package or model type
  !!
  !! Return the subcomponent type typically derived from package file type,
  !! i.e. return CHD when input is CHD6. Note this function is called on
  !! file types that are both idm integrated and not and should not set
  !! an error based on this difference.
  !!
  !<
  function idm_subcomponent_type(component, subcomponent) &
    result(subcomponent_type)
    character(len=*), intent(in) :: component !< component, e.g. GWF6
    character(len=*), intent(in) :: subcomponent !< subcomponent, e.g. CHD6
    character(len=LENCOMPONENTNAME) :: subcomponent_type
    character(len=LENCOMPONENTNAME) :: component_type
    integer(I4B) :: i, ilen, idx

    ! initialize
    subcomponent_type = ''
    idx = 0

    ! verify component
    component_type = idm_component_type(component)

    ilen = len_trim(subcomponent)
    do i = 1, ilen
      if (subcomponent(i:i) == '6' .or. subcomponent(i:i) == '-') then
      else
        idx = idx + 1
        subcomponent_type(idx:idx) = subcomponent(i:i)
      end if
    end do
  end function idm_subcomponent_type

  !> @brief model package subcomponent name
  !!
  !! Return the IDM component name, which is the package type for
  !! base packages and the package name for multi package (i.e.
  !! stress) types.
  !!
  !<
  function idm_subcomponent_name(component_type, subcomponent_type, sc_name) &
    result(subcomponent_name)
    use IdmDfnSelectorModule, only: idm_multi_package
    character(len=*), intent(in) :: component_type
    character(len=*), intent(in) :: subcomponent_type
    character(len=*), intent(in) :: sc_name
    character(len=LENPACKAGENAME) :: subcomponent_name
    subcomponent_name = ''
    if (idm_multi_package(component_type, subcomponent_type)) then
      subcomponent_name = sc_name
    else
      subcomponent_name = subcomponent_type
    end if
  end function idm_subcomponent_name

  !> @brief input file extension
  !!
  !! Return a file extension, or an empty string if
  !! not identified.
  !!
  !<
  function file_ext(filename) result(ext)
    use IdmDfnSelectorModule, only: idm_multi_package
    character(len=*), intent(in) :: filename
    character(len=LENPACKAGETYPE) :: ext
    integer(I4B) :: idx
    ! initialize
    ext = ''
    idx = 0
    ! identify '.' character position from back of string
    idx = index(filename, '.', back=.true.)
    if (idx > 0) then
      ext = filename(idx + 1:len_trim(filename))
    end if
  end function file_ext

  subroutine get_shape_from_string(shape_string, array_shape, memoryPath)
    use InputOutputModule, only: parseline
    use MemoryManagerModule, only: mem_setptr
    character(len=*), intent(in) :: shape_string
    integer(I4B), dimension(:), allocatable, intent(inout) :: array_shape
    character(len=*), intent(in) :: memoryPath !< memorypath to put loaded information
    integer(I4B) :: ndim
    integer(I4B) :: i
    integer(I4B), pointer :: int_ptr
    character(len=16), dimension(:), allocatable :: array_shape_string
    character(len=:), allocatable :: shape_string_copy

    ! parse the string into multiple words
    shape_string_copy = trim(shape_string)//' '
    call ParseLine(shape_string_copy, ndim, array_shape_string)
    allocate (array_shape(ndim))

    ! find shape in memory manager and put into array_shape
    do i = 1, ndim
      call mem_setptr(int_ptr, array_shape_string(i), memoryPath)
      array_shape(i) = int_ptr
    end do
  end subroutine get_shape_from_string

  subroutine get_layered_shape(mshape, nlay, layer_shape)
    integer(I4B), dimension(:), intent(in) :: mshape
    integer(I4B), intent(out) :: nlay
    integer(I4B), dimension(:), allocatable, intent(out) :: layer_shape
    integer(I4B) :: ndim

    ndim = size(mshape)
    nlay = 0

    if (ndim == 1) then ! disu
      nlay = 1
      allocate (layer_shape(1))
      layer_shape(1) = mshape(1)
    else if (ndim == 2) then ! disv
      nlay = mshape(1)
      allocate (layer_shape(1))
      layer_shape(1) = mshape(2)
    else if (ndim == 3) then ! disu
      nlay = mshape(1)
      allocate (layer_shape(2))
      layer_shape(1) = mshape(3) ! ncol
      layer_shape(2) = mshape(2) ! nrow
    end if
  end subroutine get_layered_shape

  !> @brief routine for setting the model shape
  !!
  !! The model shape must be set in the memory manager because
  !! individual packages need to know the shape of the arrays
  !! to read.
  !!
  !<
  subroutine set_model_shape(ftype, fname, model_mempath, dis_mempath, &
                             model_shape)
    use ConstantsModule, only: DISUNDEF, DIS, DISV, DISU, DIS2D, DISV1D
    use MemoryManagerModule, only: mem_allocate, mem_setptr, get_isize
    character(len=*), intent(in) :: ftype
    character(len=*), intent(in) :: fname
    character(len=*), intent(in) :: model_mempath
    character(len=*), intent(in) :: dis_mempath
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: model_shape
    integer(I4B), pointer :: ndim1
    integer(I4B), pointer :: ndim2
    integer(I4B), pointer :: ndim3
    integer(I4B), pointer :: ncelldim
    integer(I4B), pointer :: distype
    integer(I4B) :: dim1_size, dim2_size, dim3_size, dis_type

    ! initialize dis_type
    dis_type = DISUNDEF

    ! allocate and set model shape in model input context
    select case (ftype)
    case ('DIS6')
      ! set dis_type
      dis_type = DIS
      call get_isize('NLAY', dis_mempath, dim1_size)
      call get_isize('NROW', dis_mempath, dim2_size)
      call get_isize('NCOL', dis_mempath, dim3_size)

      if (dim1_size <= 0) then
        write (errmsg, '(a)') &
          'Required input dimension "NLAY" not found.'
        call store_error(errmsg)
      end if

      if (dim2_size <= 0) then
        write (errmsg, '(a)') &
          'Required input dimension "NROW" not found.'
        call store_error(errmsg)
      end if

      if (dim3_size <= 0) then
        write (errmsg, '(a)') &
          'Required input dimension "NCOL" not found.'
        call store_error(errmsg)
      end if

      if (dim1_size >= 1 .and. dim2_size >= 1 .and. dim3_size >= 1) then
        call mem_allocate(model_shape, 3, 'MODEL_SHAPE', model_mempath)
        call mem_setptr(ndim1, 'NLAY', dis_mempath)
        call mem_setptr(ndim2, 'NROW', dis_mempath)
        call mem_setptr(ndim3, 'NCOL', dis_mempath)
        model_shape = [ndim1, ndim2, ndim3]
      else
        call store_error_filename(fname)
      end if
    case ('DIS2D6')
      ! set dis_type
      dis_type = DIS2D
      call get_isize('NROW', dis_mempath, dim1_size)
      call get_isize('NCOL', dis_mempath, dim2_size)

      if (dim1_size <= 0) then
        write (errmsg, '(a)') &
          'Required input dimension "NROW" not found.'
        call store_error(errmsg)
      end if

      if (dim2_size <= 0) then
        write (errmsg, '(a)') &
          'Required input dimension "NCOL" not found.'
        call store_error(errmsg)
      end if

      if (dim1_size >= 1 .and. dim2_size >= 1) then
        call mem_allocate(model_shape, 2, 'MODEL_SHAPE', model_mempath)
        call mem_setptr(ndim1, 'NROW', dis_mempath)
        call mem_setptr(ndim2, 'NCOL', dis_mempath)
        model_shape = [ndim1, ndim2]
      else
        call store_error_filename(fname)
      end if
    case ('DISV6')
      ! set dis_type
      dis_type = DISV
      call get_isize('NLAY', dis_mempath, dim1_size)
      call get_isize('NCPL', dis_mempath, dim2_size)

      if (dim1_size <= 0) then
        write (errmsg, '(a)') &
          'Required input dimension "NLAY" not found.'
        call store_error(errmsg)
      end if

      if (dim2_size <= 0) then
        write (errmsg, '(a)') &
          'Required input dimension "NCPL" not found.'
        call store_error(errmsg)
      end if

      if (dim1_size >= 1 .and. dim2_size >= 1) then
        call mem_allocate(model_shape, 2, 'MODEL_SHAPE', model_mempath)
        call mem_setptr(ndim1, 'NLAY', dis_mempath)
        call mem_setptr(ndim2, 'NCPL', dis_mempath)
        model_shape = [ndim1, ndim2]
      else
        call store_error_filename(fname)
      end if
    case ('DISV2D6')
      call get_isize('NODES', dis_mempath, dim1_size)

      if (dim1_size <= 0) then
        write (errmsg, '(a)') &
          'Required input dimension "NODES" not found.'
        call store_error(errmsg)
      end if

      if (dim1_size >= 1) then
        call mem_allocate(model_shape, 1, 'MODEL_SHAPE', model_mempath)
        call mem_setptr(ndim1, 'NODES', dis_mempath)
        model_shape = [ndim1]
      else
        call store_error_filename(fname)
      end if
    case ('DISU6', 'DISV1D6')
      ! set dis_type
      if (ftype == 'DISU6') then
        dis_type = DISU
      else if (ftype == 'DISV1D6') then
        dis_type = DISV1D
      end if

      call get_isize('NODES', dis_mempath, dim1_size)

      if (dim1_size <= 0) then
        write (errmsg, '(a)') &
          'Required input dimension "NODES" not found.'
        call store_error(errmsg)
        call store_error_filename(fname)
      end if

      call mem_allocate(model_shape, 1, 'MODEL_SHAPE', model_mempath)
      call mem_setptr(ndim1, 'NODES', dis_mempath)
      model_shape = [ndim1]
    case default
      errmsg = 'Unknown discretization type.  IDM cannot set shape for "' &
               //trim(ftype)//"'"
      call store_error(errmsg)
      call store_error_filename(fname)
    end select

    ! allocate and set ncelldim in model input context
    call mem_allocate(ncelldim, 'NCELLDIM', model_mempath)
    ncelldim = size(model_shape)

    ! allocate and set distype in model input context
    call mem_allocate(distype, 'DISENUM', model_mempath)
    distype = dis_type
  end subroutine set_model_shape

  function ifind_charstr(array, str)
    use CharacterStringModule, only: CharacterStringType
    implicit none
    integer(I4B) :: ifind_charstr
    type(CharacterStringType), dimension(:), intent(in) :: array
    character(len=*) :: str
    character(len=LINELENGTH) :: compare_str
    integer(I4B) :: i
    ifind_charstr = -1
    findloop: do i = 1, size(array)
      compare_str = array(i)
      if (compare_str == str) then
        ifind_charstr = i
        exit findloop
      end if
    end do findloop
  end function ifind_charstr

  !> @brief enforce and set a single input filename provided via FILEIN keyword
  !!
  !! Set a FILEIN filename provided via an OPTIONS block.
  !! Only use this function if a maximum of one FILEIN file name
  !! string is expected.
  !!
  !! Return true if single FILEIN file name found and set, return
  !! false if FILEIN tag not found.
  !!
  !<
  function filein_fname(filename, tagname, input_mempath, input_fname) &
    result(found)
    use MemoryManagerModule, only: mem_setptr, get_isize
    use CharacterStringModule, only: CharacterStringType
    character(len=*), intent(inout) :: filename
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: input_mempath
    character(len=*), intent(in) :: input_fname
    logical(LGP) :: found
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: fnames
    integer(I4B) :: isize

    ! initialize
    found = .false.
    filename = ''

    call get_isize(tagname, input_mempath, isize)

    if (isize > 0) then
      if (isize /= 1) then
        errmsg = 'Multiple FILEIN keywords detected for tag "'//trim(tagname)// &
                 '" in OPTIONS block. Only one entry allowed.'
        call store_error(errmsg)
        call store_error_filename(input_fname)
      end if

      call mem_setptr(fnames, tagname, input_mempath)
      filename = fnames(1)
      found = .true.
    end if
  end function filein_fname

  !> @brief store an error for input exceeding internal name length
  !<
  subroutine inlen_check(input_name, mf6_name, maxlen, name_type)
    use CharacterStringModule, only: CharacterStringType
    type(CharacterStringType), intent(in) :: input_name
    character(len=*), intent(inout) :: mf6_name
    integer(I4B), intent(in) :: maxlen
    character(len=*), intent(in) :: name_type
    character(len=LINELENGTH) :: input_str
    integer(I4B) :: ilen

    ! initialize
    mf6_name = ''
    input_str = input_name
    ilen = len_trim(input_str)
    if (ilen > maxlen) then
      write (errmsg, '(a,i0,a)') &
        'Input name "'//trim(input_str)//'" exceeds maximum allowed length (', &
        maxlen, ') for '//trim(name_type)//'.'
      call store_error(errmsg)
    end if

    ! set truncated name
    mf6_name = trim(input_str)
  end subroutine inlen_check

end module SourceCommonModule
