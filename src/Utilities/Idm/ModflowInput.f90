!> @brief This module contains the ModflowInputModule
!!
!! This module contains a helper object and function
!! for accessing the ModflowInput, which is a
!! description of the structure of a modflow input
!! file.
!!
!<
module ModflowInputModule

  use KindModule, only: I4B, LGP
  use ConstantsModule, only: LENMEMPATH, LENCOMPONENTNAME, &
                             LENPACKAGETYPE
  use MemoryHelperModule, only: create_mem_path
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  use IdmDfnSelectorModule, only: block_definitions, &
                                  aggregate_definitions, &
                                  param_definitions
  use SimVariablesModule, only: idm_context

  implicit none
  private
  public :: ModflowInputType, getModflowInput

  !> @brief derived type for storing input definition for a file
  !!
  !! This derived type contains the information needed to read
  !! a specific modflow input file, including block definitions,
  !! aggregate definitions (structarrays), and individual
  !! parameter definitions.
  !!
  !<
  type ModflowInputType
    character(len=LENCOMPONENTNAME) :: pkgtype
    character(len=LENCOMPONENTNAME) :: component_type
    character(len=LENCOMPONENTNAME) :: subcomponent_type
    character(len=LENCOMPONENTNAME) :: component_name
    character(len=LENCOMPONENTNAME) :: subcomponent_name
    character(len=LENMEMPATH) :: mempath
    character(len=LENMEMPATH) :: component_mempath
    character(len=LENCOMPONENTNAME) :: load_scope
    type(InputBlockDefinitionType), dimension(:), pointer :: block_dfns
    type(InputParamDefinitionType), dimension(:), pointer :: aggregate_dfns
    type(InputParamDefinitionType), dimension(:), pointer :: param_dfns
  end type ModflowInputType

contains

  !> @brief function to return ModflowInputType
  !<
  function getModflowInput(pkgtype, component_type, subcomponent_type, &
                           component_name, subcomponent_name, &
                           load_scope, filename) &
    result(mf6_input)
    character(len=*), intent(in) :: pkgtype !< package type to load, such as DIS6, DISV6, NPF6
    character(len=*), intent(in) :: component_type !< component type, such as GWF or GWT
    character(len=*), intent(in) :: subcomponent_type !< subcomponent type, such as DIS or NPF
    character(len=*), intent(in) :: component_name !< component name, such as MYGWFMODEL
    character(len=*), intent(in) :: subcomponent_name !< subcomponent name, such as MYWELLPACKAGE
    character(len=*), intent(in) :: load_scope !< e.g. SIM, MODEL, EXCHANGE
    character(len=*), optional, intent(in) :: filename !< optional name of package input file
    type(ModflowInputType) :: mf6_input
    character(len=LENPACKAGETYPE) :: dfn_subcomponent_type

    ! set subcomponent type
    if (present(filename)) then
      dfn_subcomponent_type = update_sc_type(pkgtype, filename, component_type, &
                                             subcomponent_type)
    else
      dfn_subcomponent_type = trim(subcomponent_type)
    end if

    ! set input attributes
    mf6_input%pkgtype = trim(pkgtype)
    mf6_input%component_type = trim(component_type)
    mf6_input%subcomponent_type = trim(dfn_subcomponent_type)
    mf6_input%component_name = trim(component_name)
    mf6_input%subcomponent_name = trim(subcomponent_name)
    mf6_input%load_scope = trim(load_scope)

    ! set mempaths
    mf6_input%mempath = create_mem_path(component_name, subcomponent_name, &
                                        idm_context)
    mf6_input%component_mempath = create_mem_path(component=component_name, &
                                                  context=idm_context)

    ! set input definitions
    mf6_input%block_dfns => block_definitions(mf6_input%component_type, &
                                              mf6_input%subcomponent_type)
    mf6_input%aggregate_dfns => aggregate_definitions(mf6_input%component_type, &
                                                      mf6_input%subcomponent_type)
    mf6_input%param_dfns => param_definitions(mf6_input%component_type, &
                                              mf6_input%subcomponent_type)
  end function getModflowInput

  function update_sc_type(filetype, filename, component_type, subcomponent_type) &
    result(sc_type)
    character(len=*), intent(in) :: component_type
    character(len=*), intent(in) :: subcomponent_type
    character(len=*), intent(in) :: filetype
    character(len=*), intent(in) :: filename
    character(len=LENPACKAGETYPE) :: sc_type
    sc_type = subcomponent_type
    select case (subcomponent_type)
    case ('CHD', 'DRN', 'EVT', 'GHB', 'RCH', 'RIV', 'SCP', 'WEL')
      sc_type = readarray(filetype, filename, component_type, &
                          subcomponent_type)
    case default
    end select
  end function update_sc_type

  function readarray(filetype, filename, component_type, subcomponent_type) &
    result(sc_type)
    use ConstantsModule, only: LINELENGTH
    use InputOutputModule, only: openfile, getunit
    use BlockParserModule, only: BlockParserType
    character(len=*), intent(in) :: component_type
    character(len=*), intent(in) :: subcomponent_type
    character(len=*), intent(in) :: filetype
    character(len=*), intent(in) :: filename
    character(len=LENPACKAGETYPE) :: sc_type
    type(BlockParserType) :: parser
    integer(I4B) :: ierr, inunit
    logical(LGP) :: isfound
    logical(LGP) :: endOfBlock
    character(len=LINELENGTH) :: keyword

    sc_type = subcomponent_type
    inunit = getunit()
    call openfile(inunit, 0, trim(adjustl(filename)), filetype, &
                  'FORMATTED', 'SEQUENTIAL', 'OLD')
    call parser%Initialize(inunit, 0)

    ! get options block
    call parser%GetBlock('OPTIONS', isfound, ierr, &
                         supportOpenClose=.true., blockRequired=.false.)

    ! parse options block if detected
    if (isfound) then
      do
        call parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call parser%GetStringCaps(keyword)
        select case (keyword)
        case ('READASARRAYS')
          write (sc_type, '(a)') trim(subcomponent_type)//'A'
        case ('READARRAYGRID')
          write (sc_type, '(a)') trim(subcomponent_type)//'G'
        case default
        end select
      end do
    end if

    call parser%clear()
  end function readarray

end module ModflowInputModule
