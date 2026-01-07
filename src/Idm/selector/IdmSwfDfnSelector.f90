! ** Do Not Modify! MODFLOW 6 system generated file. **
module IdmSwfDfnSelectorModule

  use ConstantsModule, only: LENVARNAME
  use SimModule, only: store_error
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  use SwfNamInputModule
  use SwfDisv1DInputModule
  use SwfDis2DInputModule
  use SwfDisv2DInputModule
  use SwfCxsInputModule
  use SwfDfwInputModule
  use SwfIcInputModule
  use SwfOcInputModule
  use SwfCdbInputModule
  use SwfChdInputModule
  use SwfFlwInputModule
  use SwfPcpInputModule
  use SwfEvpInputModule
  use SwfStoInputModule
  use SwfZdgInputModule

  implicit none
  private
  public :: swf_param_definitions
  public :: swf_aggregate_definitions
  public :: swf_block_definitions
  public :: swf_idm_multi_package
  public :: swf_idm_subpackages
  public :: swf_idm_integrated

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

  function swf_param_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, swf_nam_param_definitions)
    case ('DISV1D')
      call set_param_pointer(input_definition, swf_disv1d_param_definitions)
    case ('DIS2D')
      call set_param_pointer(input_definition, swf_dis2d_param_definitions)
    case ('DISV2D')
      call set_param_pointer(input_definition, swf_disv2d_param_definitions)
    case ('CXS')
      call set_param_pointer(input_definition, swf_cxs_param_definitions)
    case ('DFW')
      call set_param_pointer(input_definition, swf_dfw_param_definitions)
    case ('IC')
      call set_param_pointer(input_definition, swf_ic_param_definitions)
    case ('OC')
      call set_param_pointer(input_definition, swf_oc_param_definitions)
    case ('CDB')
      call set_param_pointer(input_definition, swf_cdb_param_definitions)
    case ('CHD')
      call set_param_pointer(input_definition, swf_chd_param_definitions)
    case ('FLW')
      call set_param_pointer(input_definition, swf_flw_param_definitions)
    case ('PCP')
      call set_param_pointer(input_definition, swf_pcp_param_definitions)
    case ('EVP')
      call set_param_pointer(input_definition, swf_evp_param_definitions)
    case ('STO')
      call set_param_pointer(input_definition, swf_sto_param_definitions)
    case ('ZDG')
      call set_param_pointer(input_definition, swf_zdg_param_definitions)
    case default
    end select
    return
  end function swf_param_definitions

  function swf_aggregate_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, swf_nam_aggregate_definitions)
    case ('DISV1D')
      call set_param_pointer(input_definition, swf_disv1d_aggregate_definitions)
    case ('DIS2D')
      call set_param_pointer(input_definition, swf_dis2d_aggregate_definitions)
    case ('DISV2D')
      call set_param_pointer(input_definition, swf_disv2d_aggregate_definitions)
    case ('CXS')
      call set_param_pointer(input_definition, swf_cxs_aggregate_definitions)
    case ('DFW')
      call set_param_pointer(input_definition, swf_dfw_aggregate_definitions)
    case ('IC')
      call set_param_pointer(input_definition, swf_ic_aggregate_definitions)
    case ('OC')
      call set_param_pointer(input_definition, swf_oc_aggregate_definitions)
    case ('CDB')
      call set_param_pointer(input_definition, swf_cdb_aggregate_definitions)
    case ('CHD')
      call set_param_pointer(input_definition, swf_chd_aggregate_definitions)
    case ('FLW')
      call set_param_pointer(input_definition, swf_flw_aggregate_definitions)
    case ('PCP')
      call set_param_pointer(input_definition, swf_pcp_aggregate_definitions)
    case ('EVP')
      call set_param_pointer(input_definition, swf_evp_aggregate_definitions)
    case ('STO')
      call set_param_pointer(input_definition, swf_sto_aggregate_definitions)
    case ('ZDG')
      call set_param_pointer(input_definition, swf_zdg_aggregate_definitions)
    case default
    end select
    return
  end function swf_aggregate_definitions

  function swf_block_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputBlockDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_block_pointer(input_definition, swf_nam_block_definitions)
    case ('DISV1D')
      call set_block_pointer(input_definition, swf_disv1d_block_definitions)
    case ('DIS2D')
      call set_block_pointer(input_definition, swf_dis2d_block_definitions)
    case ('DISV2D')
      call set_block_pointer(input_definition, swf_disv2d_block_definitions)
    case ('CXS')
      call set_block_pointer(input_definition, swf_cxs_block_definitions)
    case ('DFW')
      call set_block_pointer(input_definition, swf_dfw_block_definitions)
    case ('IC')
      call set_block_pointer(input_definition, swf_ic_block_definitions)
    case ('OC')
      call set_block_pointer(input_definition, swf_oc_block_definitions)
    case ('CDB')
      call set_block_pointer(input_definition, swf_cdb_block_definitions)
    case ('CHD')
      call set_block_pointer(input_definition, swf_chd_block_definitions)
    case ('FLW')
      call set_block_pointer(input_definition, swf_flw_block_definitions)
    case ('PCP')
      call set_block_pointer(input_definition, swf_pcp_block_definitions)
    case ('EVP')
      call set_block_pointer(input_definition, swf_evp_block_definitions)
    case ('STO')
      call set_block_pointer(input_definition, swf_sto_block_definitions)
    case ('ZDG')
      call set_block_pointer(input_definition, swf_zdg_block_definitions)
    case default
    end select
    return
  end function swf_block_definitions

  function swf_idm_multi_package(subcomponent) result(multi_package)
    character(len=*), intent(in) :: subcomponent
    logical :: multi_package
    select case (subcomponent)
    case ('NAM')
      multi_package = swf_nam_multi_package
    case ('DISV1D')
      multi_package = swf_disv1d_multi_package
    case ('DIS2D')
      multi_package = swf_dis2d_multi_package
    case ('DISV2D')
      multi_package = swf_disv2d_multi_package
    case ('CXS')
      multi_package = swf_cxs_multi_package
    case ('DFW')
      multi_package = swf_dfw_multi_package
    case ('IC')
      multi_package = swf_ic_multi_package
    case ('OC')
      multi_package = swf_oc_multi_package
    case ('CDB')
      multi_package = swf_cdb_multi_package
    case ('CHD')
      multi_package = swf_chd_multi_package
    case ('FLW')
      multi_package = swf_flw_multi_package
    case ('PCP')
      multi_package = swf_pcp_multi_package
    case ('EVP')
      multi_package = swf_evp_multi_package
    case ('STO')
      multi_package = swf_sto_multi_package
    case ('ZDG')
      multi_package = swf_zdg_multi_package
    case default
      call store_error('Idm selector subcomponent not found; '//&
                       &'component="SWF"'//&
                       &', subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function swf_idm_multi_package

  function swf_idm_subpackages(subcomponent) result(subpackages)
    character(len=*), intent(in) :: subcomponent
    character(len=16), dimension(:), pointer :: subpackages
    select case (subcomponent)
    case ('NAM')
      call set_subpkg_pointer(subpackages, swf_nam_subpackages)
    case ('DISV1D')
      call set_subpkg_pointer(subpackages, swf_disv1d_subpackages)
    case ('DIS2D')
      call set_subpkg_pointer(subpackages, swf_dis2d_subpackages)
    case ('DISV2D')
      call set_subpkg_pointer(subpackages, swf_disv2d_subpackages)
    case ('CXS')
      call set_subpkg_pointer(subpackages, swf_cxs_subpackages)
    case ('DFW')
      call set_subpkg_pointer(subpackages, swf_dfw_subpackages)
    case ('IC')
      call set_subpkg_pointer(subpackages, swf_ic_subpackages)
    case ('OC')
      call set_subpkg_pointer(subpackages, swf_oc_subpackages)
    case ('CDB')
      call set_subpkg_pointer(subpackages, swf_cdb_subpackages)
    case ('CHD')
      call set_subpkg_pointer(subpackages, swf_chd_subpackages)
    case ('FLW')
      call set_subpkg_pointer(subpackages, swf_flw_subpackages)
    case ('PCP')
      call set_subpkg_pointer(subpackages, swf_pcp_subpackages)
    case ('EVP')
      call set_subpkg_pointer(subpackages, swf_evp_subpackages)
    case ('STO')
      call set_subpkg_pointer(subpackages, swf_sto_subpackages)
    case ('ZDG')
      call set_subpkg_pointer(subpackages, swf_zdg_subpackages)
    case default
    end select
    return
  end function swf_idm_subpackages

  function swf_idm_integrated(subcomponent) result(integrated)
    character(len=*), intent(in) :: subcomponent
    logical :: integrated
    integrated = .false.
    select case (subcomponent)
    case ('NAM')
      integrated = .true.
    case ('DISV1D')
      integrated = .true.
    case ('DIS2D')
      integrated = .true.
    case ('DISV2D')
      integrated = .true.
    case ('CXS')
      integrated = .true.
    case ('DFW')
      integrated = .true.
    case ('IC')
      integrated = .true.
    case ('OC')
      integrated = .true.
    case ('CDB')
      integrated = .true.
    case ('CHD')
      integrated = .true.
    case ('FLW')
      integrated = .true.
    case ('PCP')
      integrated = .true.
    case ('EVP')
      integrated = .true.
    case ('STO')
      integrated = .true.
    case ('ZDG')
      integrated = .true.
    case default
    end select
    return
  end function swf_idm_integrated

end module IdmSwfDfnSelectorModule
