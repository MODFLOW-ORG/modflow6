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
    procedure :: ainit => oc_init
    procedure :: df
    procedure :: rp
    procedure :: reset
    procedure :: destroy
  end type SettingLoadType

contains

  subroutine oc_init(this, mf6_input, component_name, component_input_name, &
                     input_name, iperblock, parser, iout)
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
    integer(I4B) :: nparam, icol

    ! init loader
    call this%DynamicPkgLoadType%init(mf6_input, component_name, &
                                      component_input_name, input_name, &
                                      iperblock, iout)
    ! initialize static loader
    call loader%load(parser, mf6_input, this%nc_vars, this%input_name, iout)

    ! create and allocate load context
    call this%ctx%init(mf6_input)
    call this%ctx%allocate_arrays()

    ! set nparam
    nparam = 1
    select case (mf6_input%subcomponent_type)
    case ('OC')
      nparam = nparam + 2
    case default
    end select

    ! allocate idts
    allocate (this%idts(nparam))

    ! set leading idts
    select case (mf6_input%subcomponent_type)
    case ('OC')
      this%idts(1)%idt => &
        idt_default(mf6_input%component_type, mf6_input%subcomponent_type, &
                    'PERIOD', 'OCACTION', 'OCACTION', 'STRING')
      this%idts(2)%idt => &
        idt_default(mf6_input%component_type, mf6_input%subcomponent_type, &
                    'PERIOD', 'RTYPE', 'RTYPE', 'STRING')
      icol = 3
    case default
      icol = 1
    end select

    ! set setting idt
    this%idts(icol)%idt => &
      idt_default(mf6_input%component_type, mf6_input%subcomponent_type, &
                  'PERIOD', 'SETTING', 'SETTING', 'STRING')
    ! force all setting tokens to be read
    this%idts(icol)%idt%shape = 'LINELENGTH'
  end subroutine oc_init

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

  subroutine destroy(this)
    use MemoryManagerModule, only: mem_deallocate
    class(SettingLoadType), intent(inout) :: this
    integer(I4B) :: icol

    if (associated(this%structarray)) then
      ! destroy the structured array reader
      call destructStructArray(this%structarray)
    end if

    do icol = 1, size(this%idts)
      ! allocate variable in memory manager
      deallocate (this%idts(icol)%idt)
      nullify (this%idts(icol)%idt)
    end do

    call this%DynamicPkgLoadType%destroy()
  end subroutine destroy

end module Mf6FileSettingLoadModule
