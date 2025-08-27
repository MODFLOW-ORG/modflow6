! ** Do Not Modify! MODFLOW 6 system generated file. **
module IdmGweDfnSelectorModule

  use ConstantsModule, only: LENVARNAME
  use SimModule, only: store_error
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  use GweNamInputModule
  use GweAdvInputModule
  use GweEslInputModule
  use GweIcInputModule
  use GweCtpInputModule
  use GweCndInputModule
  use GweDisvInputModule
  use GweDisuInputModule
  use GweDisInputModule
  use GweEstInputModule
  use GweFmiInputModule
  use GweSsmInputModule

  implicit none
  private
  public :: gwe_param_definitions
  public :: gwe_aggregate_definitions
  public :: gwe_block_definitions
  public :: gwe_idm_multi_package
  public :: gwe_idm_subpackages
  public :: gwe_idm_integrated

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

  function gwe_param_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, gwe_nam_param_definitions)
    case ('ADV')
      call set_param_pointer(input_definition, gwe_adv_param_definitions)
    case ('ESL')
      call set_param_pointer(input_definition, gwe_esl_param_definitions)
    case ('IC')
      call set_param_pointer(input_definition, gwe_ic_param_definitions)
    case ('CTP')
      call set_param_pointer(input_definition, gwe_ctp_param_definitions)
    case ('CND')
      call set_param_pointer(input_definition, gwe_cnd_param_definitions)
    case ('DISV')
      call set_param_pointer(input_definition, gwe_disv_param_definitions)
    case ('DISU')
      call set_param_pointer(input_definition, gwe_disu_param_definitions)
    case ('DIS')
      call set_param_pointer(input_definition, gwe_dis_param_definitions)
    case ('EST')
      call set_param_pointer(input_definition, gwe_est_param_definitions)
    case ('FMI')
      call set_param_pointer(input_definition, gwe_fmi_param_definitions)
    case ('SSM')
      call set_param_pointer(input_definition, gwe_ssm_param_definitions)
    case default
    end select
    return
  end function gwe_param_definitions

  function gwe_aggregate_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, gwe_nam_aggregate_definitions)
    case ('ADV')
      call set_param_pointer(input_definition, gwe_adv_aggregate_definitions)
    case ('ESL')
      call set_param_pointer(input_definition, gwe_esl_aggregate_definitions)
    case ('IC')
      call set_param_pointer(input_definition, gwe_ic_aggregate_definitions)
    case ('CTP')
      call set_param_pointer(input_definition, gwe_ctp_aggregate_definitions)
    case ('CND')
      call set_param_pointer(input_definition, gwe_cnd_aggregate_definitions)
    case ('DISV')
      call set_param_pointer(input_definition, gwe_disv_aggregate_definitions)
    case ('DISU')
      call set_param_pointer(input_definition, gwe_disu_aggregate_definitions)
    case ('DIS')
      call set_param_pointer(input_definition, gwe_dis_aggregate_definitions)
    case ('EST')
      call set_param_pointer(input_definition, gwe_est_aggregate_definitions)
    case ('FMI')
      call set_param_pointer(input_definition, gwe_fmi_aggregate_definitions)
    case ('SSM')
      call set_param_pointer(input_definition, gwe_ssm_aggregate_definitions)
    case default
    end select
    return
  end function gwe_aggregate_definitions

  function gwe_block_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputBlockDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_block_pointer(input_definition, gwe_nam_block_definitions)
    case ('ADV')
      call set_block_pointer(input_definition, gwe_adv_block_definitions)
    case ('ESL')
      call set_block_pointer(input_definition, gwe_esl_block_definitions)
    case ('IC')
      call set_block_pointer(input_definition, gwe_ic_block_definitions)
    case ('CTP')
      call set_block_pointer(input_definition, gwe_ctp_block_definitions)
    case ('CND')
      call set_block_pointer(input_definition, gwe_cnd_block_definitions)
    case ('DISV')
      call set_block_pointer(input_definition, gwe_disv_block_definitions)
    case ('DISU')
      call set_block_pointer(input_definition, gwe_disu_block_definitions)
    case ('DIS')
      call set_block_pointer(input_definition, gwe_dis_block_definitions)
    case ('EST')
      call set_block_pointer(input_definition, gwe_est_block_definitions)
    case ('FMI')
      call set_block_pointer(input_definition, gwe_fmi_block_definitions)
    case ('SSM')
      call set_block_pointer(input_definition, gwe_ssm_block_definitions)
    case default
    end select
    return
  end function gwe_block_definitions

  function gwe_idm_multi_package(subcomponent) result(multi_package)
    character(len=*), intent(in) :: subcomponent
    logical :: multi_package
    select case (subcomponent)
    case ('NAM')
      multi_package = gwe_nam_multi_package
    case ('ADV')
      multi_package = gwe_adv_multi_package
    case ('ESL')
      multi_package = gwe_esl_multi_package
    case ('IC')
      multi_package = gwe_ic_multi_package
    case ('CTP')
      multi_package = gwe_ctp_multi_package
    case ('CND')
      multi_package = gwe_cnd_multi_package
    case ('DISV')
      multi_package = gwe_disv_multi_package
    case ('DISU')
      multi_package = gwe_disu_multi_package
    case ('DIS')
      multi_package = gwe_dis_multi_package
    case ('EST')
      multi_package = gwe_est_multi_package
    case ('FMI')
      multi_package = gwe_fmi_multi_package
    case ('SSM')
      multi_package = gwe_ssm_multi_package
    case default
      call store_error('Idm selector subcomponent not found; '//&
                       &'component="GWE"'//&
                       &', subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function gwe_idm_multi_package

  function gwe_idm_subpackages(subcomponent) result(subpackages)
    character(len=*), intent(in) :: subcomponent
    character(len=16), dimension(:), pointer :: subpackages
    select case (subcomponent)
    case ('NAM')
      call set_subpkg_pointer(subpackages, gwe_nam_subpackages)
    case ('ADV')
      call set_subpkg_pointer(subpackages, gwe_adv_subpackages)
    case ('ESL')
      call set_subpkg_pointer(subpackages, gwe_esl_subpackages)
    case ('IC')
      call set_subpkg_pointer(subpackages, gwe_ic_subpackages)
    case ('CTP')
      call set_subpkg_pointer(subpackages, gwe_ctp_subpackages)
    case ('CND')
      call set_subpkg_pointer(subpackages, gwe_cnd_subpackages)
    case ('DISV')
      call set_subpkg_pointer(subpackages, gwe_disv_subpackages)
    case ('DISU')
      call set_subpkg_pointer(subpackages, gwe_disu_subpackages)
    case ('DIS')
      call set_subpkg_pointer(subpackages, gwe_dis_subpackages)
    case ('EST')
      call set_subpkg_pointer(subpackages, gwe_est_subpackages)
    case ('FMI')
      call set_subpkg_pointer(subpackages, gwe_fmi_subpackages)
    case ('SSM')
      call set_subpkg_pointer(subpackages, gwe_ssm_subpackages)
    case default
    end select
    return
  end function gwe_idm_subpackages

  function gwe_idm_integrated(subcomponent) result(integrated)
    character(len=*), intent(in) :: subcomponent
    logical :: integrated
    integrated = .false.
    select case (subcomponent)
    case ('NAM')
      integrated = .true.
    case ('ADV')
      integrated = .true.
    case ('ESL')
      integrated = .true.
    case ('IC')
      integrated = .true.
    case ('CTP')
      integrated = .true.
    case ('CND')
      integrated = .true.
    case ('DISV')
      integrated = .true.
    case ('DISU')
      integrated = .true.
    case ('DIS')
      integrated = .true.
    case ('EST')
      integrated = .true.
    case ('FMI')
      integrated = .true.
    case ('SSM')
      integrated = .true.
    case default
    end select
    return
  end function gwe_idm_integrated

end module IdmGweDfnSelectorModule
