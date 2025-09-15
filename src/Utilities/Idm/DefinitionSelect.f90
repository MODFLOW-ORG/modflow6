!> @brief This module contains the DefinitionSelectModule
!!
!! This module contains the routines for getting parameter
!! definitions, aggregate definitions, and block definitions
!! for the different package types.
!!
!<
module DefinitionSelectModule

  use KindModule, only: I4B
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType

  implicit none
  private
  public :: get_param_definition_type
  public :: get_aggregate_definition_type
  public :: split_record_dfn_tag1, split_record_dfn_tag2
  public :: idt_parse_rectype
  public :: idt_datatype
  public :: idt_default

contains

  !> @brief allocate and set RECARRAY, KEYSTRING or RECORD param list
  !<
  subroutine idt_parse_rectype(idt, cols, ncol)
    use ConstantsModule, only: LINELENGTH
    use InputOutputModule, only: parseline
    type(InputParamDefinitionType), pointer, intent(in) :: idt
    character(len=LINELENGTH), dimension(:), allocatable, &
      intent(inout) :: cols
    integer(I4B), intent(inout) :: ncol
    character(len=:), allocatable :: parse_str
    character(len=LINELENGTH), dimension(:), allocatable :: param_cols
    integer(I4B) :: param_ncol, n

    ! initialize
    if (allocated(cols)) deallocate (cols)
    ncol = 0

    ! split definition
    parse_str = trim(idt%datatype)//' '
    call parseline(parse_str, param_ncol, param_cols)

    if (param_ncol > 1) then
      if (param_cols(1) == 'RECARRAY' .or. &
          param_cols(1) == 'KEYSTRING' .or. &
          param_cols(1) == 'RECORD') then
        ! exclude 1st column
        allocate (cols(param_ncol - 1))
        do n = 2, param_ncol
          cols(n - 1) = param_cols(n)
        end do
        ! set ncol
        ncol = param_ncol - 1
      end if
    end if

    ! cleanup
    if (allocated(param_cols)) deallocate (param_cols)
    if (allocated(parse_str)) deallocate (parse_str)
  end subroutine idt_parse_rectype

  !> @brief return input definition type datatype
  !<
  function idt_datatype(idt) result(datatype)
    use ConstantsModule, only: LINELENGTH
    type(InputParamDefinitionType), pointer, intent(in) :: idt
    character(len=LINELENGTH) :: datatype
    if (idt%datatype(1:9) == 'KEYSTRING') then
      datatype = 'KEYSTRING'
    else if (idt%datatype(1:8) == 'RECARRAY') then
      datatype = 'RECARRAY'
    else if (idt%datatype(1:6) == 'RECORD') then
      datatype = 'RECORD'
    else
      datatype = idt%datatype
    end if
  end function idt_datatype

  !> @brief Return parameter definition
  !<
  function get_param_definition_type(input_definition_types, &
                                     component_type, subcomponent_type, &
                                     blockname, tagname, filename) &
    result(idt)
    type(InputParamDefinitionType), dimension(:), intent(in), target :: &
      input_definition_types
    character(len=*), intent(in) :: component_type !< component type, such as GWF or GWT
    character(len=*), intent(in) :: subcomponent_type !< subcomponent type, such as DIS or NPF
    character(len=*), intent(in) :: blockname !< name of the block
    character(len=*), intent(in) :: tagname !< name of the input tag
    character(len=*), intent(in) :: filename !< input filename
    type(InputParamDefinitionType), pointer :: idt !< corresponding InputParameterDefinitionType for this tag
    type(InputParamDefinitionType), pointer :: tmp_ptr
    integer(I4B) :: i

    nullify (idt)
    do i = 1, size(input_definition_types)
      tmp_ptr => input_definition_types(i)
      if (tmp_ptr%component_type == component_type .and. &
          tmp_ptr%subcomponent_type == subcomponent_type .and. &
          tmp_ptr%blockname == blockname .and. &
          tmp_ptr%tagname == tagname) then
        idt => input_definition_types(i)
        exit
      end if
    end do

    if (.not. associated(idt)) then
      write (errmsg, '(a,a,a,a,a)') &
        'Input file tag not found: "', trim(tagname), &
        '" in block "', trim(blockname), &
        '".'
      call store_error(errmsg)
      call store_error_filename(filename)
    end if
  end function get_param_definition_type

  !> @brief Return aggregate definition
  !<
  function get_aggregate_definition_type(input_definition_types, component_type, &
                                         subcomponent_type, blockname) result(idt)
    type(InputParamDefinitionType), dimension(:), intent(in), target :: &
      input_definition_types
    character(len=*), intent(in) :: component_type !< component type, such as GWF or GWT
    character(len=*), intent(in) :: subcomponent_type !< subcomponent type, such as DIS or NPF
    character(len=*), intent(in) :: blockname !< name of the block
    type(InputParamDefinitionType), pointer :: idt !< corresponding InputParameterDefinitionType for this block
    type(InputParamDefinitionType), pointer :: tmp_ptr
    integer(I4B) :: i

    nullify (idt)
    do i = 1, size(input_definition_types)
      tmp_ptr => input_definition_types(i)
      if (tmp_ptr%component_type == component_type .and. &
          tmp_ptr%subcomponent_type == subcomponent_type .and. &
          tmp_ptr%blockname == blockname) then
        idt => input_definition_types(i)
        exit
      end if
    end do

    if (.not. associated(idt)) then
      write (errmsg, '(a,a,a,a,a,a,a)') &
        'Idm aggregate definition not found: ', trim(blockname), &
        '. Component="', trim(component_type), &
        '", subcomponent="', trim(subcomponent_type), '".'
      call store_error(errmsg, .true.)
    end if
  end function get_aggregate_definition_type

  !> @brief Return aggregate definition
  !!
  !! Split a component RECORD datatype definition whose second element matches
  !! tagname into an array of character tokens
  !<
  subroutine split_record_dfn_tag1(input_definition_types, component_type, &
                                   subcomponent_type, tagname, nwords, words)
    use InputOutputModule, only: parseline
    type(InputParamDefinitionType), dimension(:), intent(in), target :: &
      input_definition_types
    character(len=*), intent(in) :: component_type !< component type, such as GWF or GWT
    character(len=*), intent(in) :: subcomponent_type !< subcomponent type, such as DIS or NPF
    character(len=*), intent(in) :: tagname !< name of the input tag
    integer(I4B), intent(inout) :: nwords
    character(len=40), dimension(:), allocatable, intent(inout) :: words
    type(InputParamDefinitionType), pointer :: tmp_ptr
    integer(I4B) :: i
    character(len=:), allocatable :: parse_str

    ! initialize to deallocated
    if (allocated(words)) deallocate (words)

    ! return all tokens of multi-record type that matches the first
    ! tag following the expected first token "RECORD"
    do i = 1, size(input_definition_types)

      ! initialize
      nwords = 0

      ! set ptr to current definition
      tmp_ptr => input_definition_types(i)

      ! match for definition to split
      if (tmp_ptr%component_type == component_type .and. &
          tmp_ptr%subcomponent_type == subcomponent_type .and. &
          idt_datatype(tmp_ptr) == 'RECORD') then

        ! set split string
        parse_str = trim(input_definition_types(i)%datatype)//' '

        ! split
        call parseline(parse_str, nwords, words)

        ! check for match and manage memory
        if (nwords >= 2) then
          if (words(1) == 'RECORD' .and. words(2) == tagname) then
            exit
          end if
        end if

        ! deallocate
        if (allocated(parse_str)) deallocate (parse_str)
        if (allocated(words)) deallocate (words)
      end if
    end do
  end subroutine split_record_dfn_tag1

  !> @brief Return aggregate definition
  !!
  !! Split a component RECORD datatype definition whose second and third
  !! elements match tagnames into an array of character tokens
  !<
  subroutine split_record_dfn_tag2(input_definition_types, component_type, &
                                   subcomponent_type, tagname, tag2, nwords, &
                                   words)
    use InputOutputModule, only: parseline
    type(InputParamDefinitionType), dimension(:), intent(in), target :: &
      input_definition_types
    character(len=*), intent(in) :: component_type !< component type, such as GWF or GWT
    character(len=*), intent(in) :: subcomponent_type !< subcomponent type, such as DIS or NPF
    character(len=*), intent(in) :: tagname !< name of the input tag
    character(len=*), intent(in) :: tag2
    integer(I4B), intent(inout) :: nwords
    character(len=40), dimension(:), allocatable, intent(inout) :: words
    type(InputParamDefinitionType), pointer :: tmp_ptr
    integer(I4B) :: i
    character(len=:), allocatable :: parse_str

    ! initialize to deallocated
    if (allocated(words)) deallocate (words)

    ! return all tokens of multi-record type that matches the first
    ! tag following the expected first token "RECORD"
    do i = 1, size(input_definition_types)

      ! initialize
      nwords = 0

      ! set ptr to current definition
      tmp_ptr => input_definition_types(i)

      ! match for definition to split
      if (tmp_ptr%component_type == component_type .and. &
          tmp_ptr%subcomponent_type == subcomponent_type .and. &
          idt_datatype(tmp_ptr) == 'RECORD') then

        ! set split string
        parse_str = trim(input_definition_types(i)%datatype)//' '

        ! split
        call parseline(parse_str, nwords, words)

        ! check for match and manage memory
        if (nwords >= 2) then
          if (words(1) == 'RECORD' .and. &
              words(2) == tagname .and. &
              words(3) == tag2) then
            exit
          end if
        end if

        ! deallocate
        if (allocated(parse_str)) deallocate (parse_str)
        if (allocated(words)) deallocate (words)
      end if
    end do
  end subroutine split_record_dfn_tag2

  !> @brief return allocated input definition type
  !<
  function idt_default(component_type, subcomponent_type, blockname, tagname, &
                       mf6varname, datatype) result(idt)
    ! -- modules
    ! -- dummy
    character(len=*), intent(in) :: component_type
    character(len=*), intent(in) :: subcomponent_type
    character(len=*), intent(in) :: blockname
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: mf6varname
    character(len=*), intent(in) :: datatype
    ! -- result
    type(InputParamDefinitionType), pointer :: idt

    allocate (idt)

    idt%component_type = trim(component_type)
    idt%subcomponent_type = trim(subcomponent_type)
    idt%blockname = trim(blockname)
    idt%tagname = trim(tagname)
    idt%mf6varname = trim(mf6varname)
    idt%datatype = trim(datatype)
    idt%shape = ''
    idt%required = .true.
    idt%prerelease = .false.
    idt%in_record = .false.
    idt%preserve_case = .false.
    idt%layered = .false.
    idt%timeseries = .false.
  end function idt_default

end module DefinitionSelectModule
