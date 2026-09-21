! ** Do Not Modify! MODFLOW 6 system generated file. **
module IdmSimDfnSelectorModule

  use ConstantsModule, only: LENVARNAME
  use SimModule, only: store_error
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  use SimNamInputModule
  use SimTdisInputModule

  implicit none
  private
  public :: sim_param_definitions
  public :: sim_aggregate_definitions
  public :: sim_block_definitions
  public :: sim_idm_multi_package
  public :: sim_idm_is_advanced
  public :: sim_idm_subpackages
  public :: sim_idm_integrated

contains

  subroutine set_param_pointer(input_dfn, input_dfn_target)
    type(InputParamDefinitionType), dimension(:), pointer :: input_dfn
    type(InputParamDefinitionType), dimension(:), target :: input_dfn_target
    input_dfn => input_dfn_target
  end subroutine set_param_pointer

  subroutine set_block_pointer(input_dfn, input_dfn_target)
    type(InputBlockDefinitionType), dimension(:), pointer :: input_dfn
    type(InputBlockDefinitionType), dimension(:), target :: input_dfn_target
    input_dfn => input_dfn_target
  end subroutine set_block_pointer

  subroutine set_subpkg_pointer(subpkg_list, subpkg_list_target)
    character(len=16), dimension(:), pointer :: subpkg_list
    character(len=16), dimension(:), target :: subpkg_list_target
    subpkg_list => subpkg_list_target
  end subroutine set_subpkg_pointer

  function sim_param_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, sim_nam_param_definitions)
    case ('TDIS')
      call set_param_pointer(input_definition, sim_tdis_param_definitions)
    case default
    end select
    return
  end function sim_param_definitions

  function sim_aggregate_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, sim_nam_aggregate_definitions)
    case ('TDIS')
      call set_param_pointer(input_definition, sim_tdis_aggregate_definitions)
    case default
    end select
    return
  end function sim_aggregate_definitions

  function sim_block_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputBlockDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_block_pointer(input_definition, sim_nam_block_definitions)
    case ('TDIS')
      call set_block_pointer(input_definition, sim_tdis_block_definitions)
    case default
    end select
    return
  end function sim_block_definitions

  function sim_idm_multi_package(subcomponent) result(multi_package)
    character(len=*), intent(in) :: subcomponent
    logical :: multi_package
    select case (subcomponent)
    case ('NAM')
      multi_package = sim_nam_multi_package
    case ('TDIS')
      multi_package = sim_tdis_multi_package
    case default
      call store_error('Idm selector subcomponent not found; '//&
                       &'component="SIM"'//&
                       &', subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function sim_idm_multi_package

  function sim_idm_is_advanced(subcomponent) result(is_advanced)
    character(len=*), intent(in) :: subcomponent
    logical :: is_advanced
    select case (subcomponent)
    case ('NAM')
      is_advanced = sim_nam_is_advanced
    case ('TDIS')
      is_advanced = sim_tdis_is_advanced
    case default
      call store_error('Idm selector subcomponent not found; '//&
                       &'component="SIM"'//&
                       &', subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function sim_idm_is_advanced

  function sim_idm_subpackages(subcomponent) result(subpackages)
    character(len=*), intent(in) :: subcomponent
    character(len=16), dimension(:), pointer :: subpackages
    select case (subcomponent)
    case ('NAM')
      call set_subpkg_pointer(subpackages, sim_nam_subpackages)
    case ('TDIS')
      call set_subpkg_pointer(subpackages, sim_tdis_subpackages)
    case default
    end select
    return
  end function sim_idm_subpackages

  function sim_idm_integrated(subcomponent) result(integrated)
    character(len=*), intent(in) :: subcomponent
    logical :: integrated
    integrated = .false.
    select case (subcomponent)
    case ('NAM')
      integrated = .true.
    case ('TDIS')
      integrated = .true.
    case default
    end select
    return
  end function sim_idm_integrated

end module IdmSimDfnSelectorModule
