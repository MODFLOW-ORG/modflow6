! ** Do Not Modify! MODFLOW 6 system generated file. **
module IdmUtlDfnSelectorModule

  use ConstantsModule, only: LENVARNAME
  use SimModule, only: store_error
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  use UtlHpcInputModule
  use UtlNcfInputModule
  use UtlObsInputModule

  implicit none
  private
  public :: utl_param_definitions
  public :: utl_aggregate_definitions
  public :: utl_block_definitions
  public :: utl_idm_multi_package
  public :: utl_idm_subpackages
  public :: utl_idm_integrated

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

  function utl_param_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('HPC')
      call set_param_pointer(input_definition, utl_hpc_param_definitions)
    case ('NCF')
      call set_param_pointer(input_definition, utl_ncf_param_definitions)
    case ('OBS')
      call set_param_pointer(input_definition, utl_obs_param_definitions)
    case default
    end select
    return
  end function utl_param_definitions

  function utl_aggregate_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('HPC')
      call set_param_pointer(input_definition, utl_hpc_aggregate_definitions)
    case ('NCF')
      call set_param_pointer(input_definition, utl_ncf_aggregate_definitions)
    case ('OBS')
      call set_param_pointer(input_definition, utl_obs_aggregate_definitions)
    case default
    end select
    return
  end function utl_aggregate_definitions

  function utl_block_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputBlockDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('HPC')
      call set_block_pointer(input_definition, utl_hpc_block_definitions)
    case ('NCF')
      call set_block_pointer(input_definition, utl_ncf_block_definitions)
    case ('OBS')
      call set_block_pointer(input_definition, utl_obs_block_definitions)
    case default
    end select
    return
  end function utl_block_definitions

  function utl_idm_multi_package(subcomponent) result(multi_package)
    character(len=*), intent(in) :: subcomponent
    logical :: multi_package
    select case (subcomponent)
    case ('HPC')
      multi_package = utl_hpc_multi_package
    case ('NCF')
      multi_package = utl_ncf_multi_package
    case ('OBS')
      multi_package = utl_obs_multi_package
    case default
      call store_error('Idm selector subcomponent not found; '//&
                       &'component="UTL"'//&
                       &', subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function utl_idm_multi_package

  function utl_idm_subpackages(subcomponent) result(subpackages)
    character(len=*), intent(in) :: subcomponent
    character(len=16), dimension(:), pointer :: subpackages
    select case (subcomponent)
    case ('HPC')
      call set_subpkg_pointer(subpackages, utl_hpc_subpackages)
    case ('NCF')
      call set_subpkg_pointer(subpackages, utl_ncf_subpackages)
    case ('OBS')
      call set_subpkg_pointer(subpackages, utl_obs_subpackages)
    case default
    end select
    return
  end function utl_idm_subpackages

  function utl_idm_integrated(subcomponent) result(integrated)
    character(len=*), intent(in) :: subcomponent
    logical :: integrated
    integrated = .false.
    select case (subcomponent)
    case ('HPC')
      integrated = .true.
    case ('NCF')
      integrated = .true.
    case ('OBS')
      integrated = .true.
    case default
    end select
    return
  end function utl_idm_integrated

end module IdmUtlDfnSelectorModule
