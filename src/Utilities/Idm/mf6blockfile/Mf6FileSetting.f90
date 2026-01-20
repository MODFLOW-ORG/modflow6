!> @brief This module contains the Mf6FileSettingLoadModule
!!
!! This module contains the routines for reading a dfn setting
!! with multiple tokens into a single CharacterStringType array
!! using the StructArrayType.
!!
!! Preceding param columns can be configured per package (e.g. OC)
!!
!<
module Mf6FileSettingLoadModule

  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LINELENGTH
  use InputDefinitionModule, only: InputParamDefinitionType
  use ModflowInputModule, only: ModflowInputType
  use AsciiInputLoadTypeModule, only: AsciiDynamicPkgLoadBaseType
  use StructArrayModule, only: StructArrayType, constructStructArray, &
                               destructStructArray
  use BlockParserModule, only: BlockParserType
  use LoadContextModule, only: LoadContextType

  implicit none
  private
  public :: SettingLoadType

  !> @brief Pointer type for read state variable
  !<
  type IdtPtrType
    type(InputParamDefinitionType), pointer :: idt
  end type IdtPtrType

  !> @brief Setting package loader
  !!
  !<
  type, extends(AsciiDynamicPkgLoadBaseType) :: SettingLoadType
    type(StructArrayType), pointer :: structarray => null() !< struct array list based load type
    type(IdtPtrType), dimension(:), allocatable :: idts !< idts for struct array input cols
    type(LoadContextType) :: ctx !< input load context
  contains
    procedure :: ainit => settingload_init
    procedure :: df
    procedure :: rp
    procedure :: reset
    procedure :: set_leading_idts
    procedure :: set_idts
    procedure :: set_oc_idts
    procedure :: destroy
  end type SettingLoadType

contains

  subroutine settingload_init(this, mf6_input, component_name, &
                              component_input_name, input_name, &
                              iperblock, parser, iout)
    use DefinitionSelectModule, only: idt_default
    use LoadMf6FileModule, only: LoadMf6FileType
    class(SettingLoadType), intent(inout) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), intent(in) :: component_name
    character(len=*), intent(in) :: component_input_name
    character(len=*), intent(in) :: input_name
    integer(I4B), intent(in) :: iperblock
    type(BlockParserType), pointer, intent(inout) :: parser
    integer(I4B), intent(in) :: iout
    type(LoadMf6FileType) :: loader
    integer(I4B) :: icol, numset

    ! init loader
    call this%DynamicPkgLoadType%init(mf6_input, component_name, &
                                      component_input_name, input_name, &
                                      iperblock, iout)
    ! initialize static loader
    call loader%load(parser, mf6_input, this%nc_vars, this%input_name, iout)

    ! create and allocate load context
    call this%ctx%init(mf6_input)
    call this%ctx%allocate_arrays()

    ! allocate idt arrays and set idts for leading cols
    numset = this%set_leading_idts()
    icol = numset + 1

    ! set setting idt
    this%idts(icol)%idt => &
      idt_default(mf6_input%component_type, mf6_input%subcomponent_type, &
                  'PERIOD', 'SETTING', 'SETTING', 'STRING')
    ! force all setting tokens to be read
    this%idts(icol)%idt%shape = 'LINELENGTH'
  end subroutine settingload_init

  subroutine df(this)
    class(SettingLoadType), intent(inout) :: this
  end subroutine df

  subroutine rp(this, parser)
    class(SettingLoadType), intent(inout) :: this
    type(BlockParserType), pointer, intent(inout) :: parser

    ! recreate structarray for period load
    call this%reset()

    ! read from ascii
    this%ctx%nbound = &
      this%structarray%read_from_parser(parser, .false., this%iout)
  end subroutine rp

  subroutine reset(this)
    class(SettingLoadType), intent(inout) :: this
    integer(I4B) :: icol

    if (associated(this%structarray)) then
      ! destroy the structured array reader
      call destructStructArray(this%structarray)
    end if

    ! construct and set up the struct array object
    this%structarray => constructStructArray(this%mf6_input, size(this%idts), &
                                             -1, 0, this%mf6_input%mempath, &
                                             this%mf6_input%component_mempath)
    ! set up struct array
    do icol = 1, size(this%idts)
      ! allocate variable in memory manager
      call this%structarray%mem_create_vector(icol, this%idts(icol)%idt)
    end do
  end subroutine reset

  function set_leading_idts(this) result(nset)
    class(SettingLoadType), intent(inout) :: this
    integer(I4B) :: nset

    ! This loader is intended to be generic for LIST based dynamic input that
    ! includes any number of columns preceding a final keystring type
    ! setting column. This routine handles exception cases. OC, for example,
    ! which could adhere to this format however changes to the dfn to make it
    ! so would introduce breaking changes to existing FloPy3 user scripts.

    ! set leading idts
    select case (this%mf6_input%subcomponent_type)
    case ('OC')
      nset = this%set_oc_idts()
    case default
      nset = this%set_idts()
    end select
  end function set_leading_idts

  function set_idts(this) result(nset)
    use SimModule, only: store_error, store_error_filename
    use InputDefinitionModule, only: InputParamDefinitionType
    use DefinitionSelectModule, only: idt_parse_rectype, &
                                      get_aggregate_definition_type
    class(SettingLoadType), intent(inout) :: this
    integer(I4B) :: nset
    type(InputParamDefinitionType), pointer :: aidt, idt
    character(len=LINELENGTH), dimension(:), allocatable :: cols
    integer(I4B) :: nparam, iparam, ip, ilen

    nset = 0

    aidt => &
      get_aggregate_definition_type(this%mf6_input%aggregate_dfns, &
                                    this%mf6_input%component_type, &
                                    this%mf6_input%subcomponent_type, &
                                    'PERIOD')

    call idt_parse_rectype(aidt, cols, nparam)

    ! allocate idts
    ilen = len_trim(cols(nparam))
    if (cols(nparam) (ilen - 6:ilen) == 'SETTING') then
      allocate (this%idts(nparam))
      nset = nparam - 1
    else
      call store_error('Internal IDM error: trailing setting param not found')
      call store_error_filename(this%input_name)
    end if

    ! set leading idts
    do iparam = 1, nset
      do ip = 1, size(this%mf6_input%param_dfns)
        idt => this%mf6_input%param_dfns(ip)
        if (idt%tagname == cols(iparam)) then
          this%idts(iparam)%idt => idt
        end if
      end do
    end do

    if (allocated(cols)) deallocate (cols)
    return
  end function set_idts

  function set_oc_idts(this) result(nset)
    use InputDefinitionModule, only: InputParamDefinitionType
    use DefinitionSelectModule, only: idt_default
    class(SettingLoadType), intent(inout) :: this
    integer(I4B) :: nset
    type(InputParamDefinitionType), pointer :: idt, idt_ocaction
    integer(I4B) :: nparam, ip

    ! generalize first kw field to a string
    idt_ocaction => &
      idt_default(this%mf6_input%component_type, &
                  this%mf6_input%subcomponent_type, &
                  'PERIOD', 'OCACTION', 'OCACTION', 'STRING')

    ! set 2 leading params
    nset = 2

    ! allocate for 3 including ocsetting
    nparam = nset + 1
    allocate (this%idts(nparam))

    ! first generalized idt
    this%idts(1)%idt => idt_ocaction

    ! set rtype idt
    do ip = 1, size(this%mf6_input%param_dfns)
      idt => this%mf6_input%param_dfns(ip)
      if (idt%tagname == 'RTYPE') then
        this%idts(2)%idt => idt
      end if
    end do
  end function set_oc_idts

  subroutine destroy(this)
    use MemoryManagerModule, only: mem_deallocate
    class(SettingLoadType), intent(inout) :: this
    integer(I4B) :: icol

    if (associated(this%structarray)) then
      ! destroy the structured array reader
      call destructStructArray(this%structarray)
    end if

    do icol = 1, size(this%idts)
      !deallocate (this%idts(icol)%idt)
      nullify (this%idts(icol)%idt)
    end do

    if (allocated(this%idts)) deallocate (this%idts)

    call this%DynamicPkgLoadType%destroy()
  end subroutine destroy

end module Mf6FileSettingLoadModule
