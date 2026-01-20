! ** Do Not Modify! MODFLOW 6 system generated file. **
module IdmDfnSelectorModule

  use ConstantsModule, only: LENVARNAME
  use SimModule, only: store_error
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  use IdmSimDfnSelectorModule
  use IdmGwfDfnSelectorModule
  use IdmGwtDfnSelectorModule
  use IdmGweDfnSelectorModule
  use IdmSwfDfnSelectorModule
  use IdmChfDfnSelectorModule
  use IdmOlfDfnSelectorModule
  use IdmPrtDfnSelectorModule
  use IdmExgDfnSelectorModule
  use IdmUtlDfnSelectorModule

  implicit none
  private
  public :: param_definitions
  public :: aggregate_definitions
  public :: block_definitions
  public :: idm_multi_package
  public :: idm_subpackages
  public :: idm_integrated
  public :: idm_component

contains

  function param_definitions(component, subcomponent) result(input_definition)
    character(len=*), intent(in) :: component
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (component)
    case ('SIM')
      input_definition => sim_param_definitions(subcomponent)
    case ('GWF')
      input_definition => gwf_param_definitions(subcomponent)
    case ('GWT')
      input_definition => gwt_param_definitions(subcomponent)
    case ('GWE')
      input_definition => gwe_param_definitions(subcomponent)
    case ('SWF')
      input_definition => swf_param_definitions(subcomponent)
    case ('CHF')
      input_definition => chf_param_definitions(subcomponent)
    case ('OLF')
      input_definition => olf_param_definitions(subcomponent)
    case ('PRT')
      input_definition => prt_param_definitions(subcomponent)
    case ('EXG')
      input_definition => exg_param_definitions(subcomponent)
    case ('UTL')
      input_definition => utl_param_definitions(subcomponent)
    case default
    end select
    return
  end function param_definitions

  function aggregate_definitions(component, subcomponent) result(input_definition)
    character(len=*), intent(in) :: component
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (component)
    case ('SIM')
      input_definition => sim_aggregate_definitions(subcomponent)
    case ('GWF')
      input_definition => gwf_aggregate_definitions(subcomponent)
    case ('GWT')
      input_definition => gwt_aggregate_definitions(subcomponent)
    case ('GWE')
      input_definition => gwe_aggregate_definitions(subcomponent)
    case ('SWF')
      input_definition => swf_aggregate_definitions(subcomponent)
    case ('CHF')
      input_definition => chf_aggregate_definitions(subcomponent)
    case ('OLF')
      input_definition => olf_aggregate_definitions(subcomponent)
    case ('PRT')
      input_definition => prt_aggregate_definitions(subcomponent)
    case ('EXG')
      input_definition => exg_aggregate_definitions(subcomponent)
    case ('UTL')
      input_definition => utl_aggregate_definitions(subcomponent)
    case default
    end select
    return
  end function aggregate_definitions

  function block_definitions(component, subcomponent) result(input_definition)
    character(len=*), intent(in) :: component
    character(len=*), intent(in) :: subcomponent
    type(InputBlockDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (component)
    case ('SIM')
      input_definition => sim_block_definitions(subcomponent)
    case ('GWF')
      input_definition => gwf_block_definitions(subcomponent)
    case ('GWT')
      input_definition => gwt_block_definitions(subcomponent)
    case ('GWE')
      input_definition => gwe_block_definitions(subcomponent)
    case ('SWF')
      input_definition => swf_block_definitions(subcomponent)
    case ('CHF')
      input_definition => chf_block_definitions(subcomponent)
    case ('OLF')
      input_definition => olf_block_definitions(subcomponent)
    case ('PRT')
      input_definition => prt_block_definitions(subcomponent)
    case ('EXG')
      input_definition => exg_block_definitions(subcomponent)
    case ('UTL')
      input_definition => utl_block_definitions(subcomponent)
    case default
    end select
    return
  end function block_definitions

  function idm_multi_package(component, subcomponent) result(multi_package)
    character(len=*), intent(in) :: component
    character(len=*), intent(in) :: subcomponent
    logical :: multi_package
    select case (component)
    case ('SIM')
      multi_package = sim_idm_multi_package(subcomponent)
    case ('GWF')
      multi_package = gwf_idm_multi_package(subcomponent)
    case ('GWT')
      multi_package = gwt_idm_multi_package(subcomponent)
    case ('GWE')
      multi_package = gwe_idm_multi_package(subcomponent)
    case ('SWF')
      multi_package = swf_idm_multi_package(subcomponent)
    case ('CHF')
      multi_package = chf_idm_multi_package(subcomponent)
    case ('OLF')
      multi_package = olf_idm_multi_package(subcomponent)
    case ('PRT')
      multi_package = prt_idm_multi_package(subcomponent)
    case ('EXG')
      multi_package = exg_idm_multi_package(subcomponent)
    case ('UTL')
      multi_package = utl_idm_multi_package(subcomponent)
    case default
      call store_error('Idm selector component not found; '//&
                       &'component="'//trim(component)//&
                       &'", subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function idm_multi_package

  function idm_subpackages(component, subcomponent) result(subpackages)
    character(len=*), intent(in) :: component
    character(len=*), intent(in) :: subcomponent
    character(len=16), dimension(:), pointer :: subpackages
    select case (component)
    case ('SIM')
      subpackages => sim_idm_subpackages(subcomponent)
    case ('GWF')
      subpackages => gwf_idm_subpackages(subcomponent)
    case ('GWT')
      subpackages => gwt_idm_subpackages(subcomponent)
    case ('GWE')
      subpackages => gwe_idm_subpackages(subcomponent)
    case ('SWF')
      subpackages => swf_idm_subpackages(subcomponent)
    case ('CHF')
      subpackages => chf_idm_subpackages(subcomponent)
    case ('OLF')
      subpackages => olf_idm_subpackages(subcomponent)
    case ('PRT')
      subpackages => prt_idm_subpackages(subcomponent)
    case ('EXG')
      subpackages => exg_idm_subpackages(subcomponent)
    case ('UTL')
      subpackages => utl_idm_subpackages(subcomponent)
    case default
      call store_error('Idm selector component not found; '//&
                       &'component="'//trim(component)//&
                       &'", subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function idm_subpackages

  function idm_integrated(component, subcomponent) result(integrated)
    character(len=*), intent(in) :: component
    character(len=*), intent(in) :: subcomponent
    logical :: integrated
    integrated = .false.
    select case (component)
    case ('SIM')
      integrated = sim_idm_integrated(subcomponent)
    case ('GWF')
      integrated = gwf_idm_integrated(subcomponent)
    case ('GWT')
      integrated = gwt_idm_integrated(subcomponent)
    case ('GWE')
      integrated = gwe_idm_integrated(subcomponent)
    case ('SWF')
      integrated = swf_idm_integrated(subcomponent)
    case ('CHF')
      integrated = chf_idm_integrated(subcomponent)
    case ('OLF')
      integrated = olf_idm_integrated(subcomponent)
    case ('PRT')
      integrated = prt_idm_integrated(subcomponent)
    case ('EXG')
      integrated = exg_idm_integrated(subcomponent)
    case ('UTL')
      integrated = utl_idm_integrated(subcomponent)
    case default
    end select
    return
  end function idm_integrated

  function idm_component(component) result(integrated)
    character(len=*), intent(in) :: component
    logical :: integrated
    integrated = .false.
    select case (component)
    case ('SIM')
      integrated = .true.
    case ('GWF')
      integrated = .true.
    case ('GWT')
      integrated = .true.
    case ('GWE')
      integrated = .true.
    case ('SWF')
      integrated = .true.
    case ('CHF')
      integrated = .true.
    case ('OLF')
      integrated = .true.
    case ('PRT')
      integrated = .true.
    case ('EXG')
      integrated = .true.
    case ('UTL')
      integrated = .true.
    case default
    end select
    return
  end function idm_component

end module IdmDfnSelectorModule
