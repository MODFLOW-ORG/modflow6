!> @brief This module contains the Mf6FileStoInputModule
!!
!! This module contains the routines for reading STO period block input
!!
!<
module Mf6FileStoInputModule

  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LINELENGTH
  use InputDefinitionModule, only: InputParamDefinitionType
  use MemoryManagerModule, only: mem_setptr, mem_allocate
  use ModflowInputModule, only: ModflowInputType
  use AsciiInputLoadTypeModule, only: AsciiDynamicPkgLoadBaseType

  implicit none
  private
  public :: StoInputType

  !> @brief STO package loader
  !!
  !<
  type, extends(AsciiDynamicPkgLoadBaseType) :: StoInputType
    character(len=LINELENGTH), pointer :: storage => null()
  contains
    procedure :: ainit => sto_init
    procedure :: rp => sto_rp
    procedure :: destroy => sto_destroy
  end type StoInputType

contains

  subroutine sto_init(this, mf6_input, component_name, component_input_name, &
                      input_name, iperblock, parser, iout)
    use MemoryManagerExtModule, only: mem_set_value
    use BlockParserModule, only: BlockParserType
    use LoadMf6FileModule, only: LoadMf6FileType
    class(StoInputType), intent(inout) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), intent(in) :: component_name
    character(len=*), intent(in) :: component_input_name
    character(len=*), intent(in) :: input_name
    integer(I4B), intent(in) :: iperblock
    type(BlockParserType), pointer, intent(inout) :: parser
    integer(I4B), intent(in) :: iout
    type(LoadMf6FileType) :: loader

    ! init loader
    call this%DynamicPkgLoadType%init(mf6_input, component_name, &
                                      component_input_name, input_name, &
                                      iperblock, iout)
    ! initialize static loader
    call loader%load(parser, mf6_input, this%nc_vars, this%input_name, iout)

    ! allocate storage string
    call mem_allocate(this%storage, LINELENGTH, 'STORAGE', this%mf6_input%mempath)

    ! initialize storage to TRANSIENT (model iss=0)
    this%storage = 'TRANSIENT'
  end subroutine sto_init

  subroutine sto_rp(this, parser)
    use BlockParserModule, only: BlockParserType
    use DefinitionSelectModule, only: get_param_definition_type
    use IdmLoggerModule, only: idm_log_header, idm_log_close, idm_log_var
    class(StoInputType), intent(inout) :: this
    type(BlockParserType), pointer, intent(inout) :: parser
    character(len=LINELENGTH) :: tagname
    type(InputParamDefinitionType), pointer :: idt
    logical(LGP) :: endOfBlock

    ! read next line
    call parser%GetNextLine(endOfBlock)

    ! return if no input
    if (endOfBlock) return

    ! read the tag
    call parser%GetStringCaps(tagname)

    ! verify tag is supported
    idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                     this%mf6_input%component_type, &
                                     this%mf6_input%subcomponent_type, &
                                     'PERIOD', tagname, this%input_name)
    ! set storage
    this%storage = idt%tagname

    ! only one input line is expected, terminate block
    call parser%terminateblock()

    ! log lst file header
    call idm_log_header(this%mf6_input%component_name, &
                        this%mf6_input%subcomponent_name, this%iout)

    call idm_log_var(this%storage, tagname, this%mf6_input%mempath, this%iout)

    ! close logging statement
    call idm_log_close(this%mf6_input%component_name, &
                       this%mf6_input%subcomponent_name, this%iout)
  end subroutine sto_rp

  subroutine sto_destroy(this)
    use MemoryManagerModule, only: mem_deallocate
    class(StoInputType), intent(inout) :: this
    call mem_deallocate(this%storage, 'STORAGE', this%mf6_input%mempath)
    call this%DynamicPkgLoadType%destroy()
  end subroutine sto_destroy

end module Mf6FileStoInputModule
