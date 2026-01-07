! ** Do Not Modify! MODFLOW 6 system generated file. **
module IdmChfDfnSelectorModule

  use ConstantsModule, only: LENVARNAME
  use SimModule, only: store_error
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  use ChfNamInputModule
  use ChfDisv1DInputModule
  use ChfCxsInputModule
  use ChfDfwInputModule
  use ChfIcInputModule
  use ChfOcInputModule
  use ChfCdbInputModule
  use ChfChdInputModule
  use ChfFlwInputModule
  use ChfPcpInputModule
  use ChfEvpInputModule
  use ChfStoInputModule
  use ChfZdgInputModule

  implicit none
  private
  public :: chf_param_definitions
  public :: chf_aggregate_definitions
  public :: chf_block_definitions
  public :: chf_idm_multi_package
  public :: chf_idm_subpackages
  public :: chf_idm_integrated

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

  function chf_param_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, chf_nam_param_definitions)
    case ('DISV1D')
      call set_param_pointer(input_definition, chf_disv1d_param_definitions)
    case ('CXS')
      call set_param_pointer(input_definition, chf_cxs_param_definitions)
    case ('DFW')
      call set_param_pointer(input_definition, chf_dfw_param_definitions)
    case ('IC')
      call set_param_pointer(input_definition, chf_ic_param_definitions)
    case ('OC')
      call set_param_pointer(input_definition, chf_oc_param_definitions)
    case ('CDB')
      call set_param_pointer(input_definition, chf_cdb_param_definitions)
    case ('CHD')
      call set_param_pointer(input_definition, chf_chd_param_definitions)
    case ('FLW')
      call set_param_pointer(input_definition, chf_flw_param_definitions)
    case ('PCP')
      call set_param_pointer(input_definition, chf_pcp_param_definitions)
    case ('EVP')
      call set_param_pointer(input_definition, chf_evp_param_definitions)
    case ('STO')
      call set_param_pointer(input_definition, chf_sto_param_definitions)
    case ('ZDG')
      call set_param_pointer(input_definition, chf_zdg_param_definitions)
    case default
    end select
    return
  end function chf_param_definitions

  function chf_aggregate_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, chf_nam_aggregate_definitions)
    case ('DISV1D')
      call set_param_pointer(input_definition, chf_disv1d_aggregate_definitions)
    case ('CXS')
      call set_param_pointer(input_definition, chf_cxs_aggregate_definitions)
    case ('DFW')
      call set_param_pointer(input_definition, chf_dfw_aggregate_definitions)
    case ('IC')
      call set_param_pointer(input_definition, chf_ic_aggregate_definitions)
    case ('OC')
      call set_param_pointer(input_definition, chf_oc_aggregate_definitions)
    case ('CDB')
      call set_param_pointer(input_definition, chf_cdb_aggregate_definitions)
    case ('CHD')
      call set_param_pointer(input_definition, chf_chd_aggregate_definitions)
    case ('FLW')
      call set_param_pointer(input_definition, chf_flw_aggregate_definitions)
    case ('PCP')
      call set_param_pointer(input_definition, chf_pcp_aggregate_definitions)
    case ('EVP')
      call set_param_pointer(input_definition, chf_evp_aggregate_definitions)
    case ('STO')
      call set_param_pointer(input_definition, chf_sto_aggregate_definitions)
    case ('ZDG')
      call set_param_pointer(input_definition, chf_zdg_aggregate_definitions)
    case default
    end select
    return
  end function chf_aggregate_definitions

  function chf_block_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputBlockDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_block_pointer(input_definition, chf_nam_block_definitions)
    case ('DISV1D')
      call set_block_pointer(input_definition, chf_disv1d_block_definitions)
    case ('CXS')
      call set_block_pointer(input_definition, chf_cxs_block_definitions)
    case ('DFW')
      call set_block_pointer(input_definition, chf_dfw_block_definitions)
    case ('IC')
      call set_block_pointer(input_definition, chf_ic_block_definitions)
    case ('OC')
      call set_block_pointer(input_definition, chf_oc_block_definitions)
    case ('CDB')
      call set_block_pointer(input_definition, chf_cdb_block_definitions)
    case ('CHD')
      call set_block_pointer(input_definition, chf_chd_block_definitions)
    case ('FLW')
      call set_block_pointer(input_definition, chf_flw_block_definitions)
    case ('PCP')
      call set_block_pointer(input_definition, chf_pcp_block_definitions)
    case ('EVP')
      call set_block_pointer(input_definition, chf_evp_block_definitions)
    case ('STO')
      call set_block_pointer(input_definition, chf_sto_block_definitions)
    case ('ZDG')
      call set_block_pointer(input_definition, chf_zdg_block_definitions)
    case default
    end select
    return
  end function chf_block_definitions

  function chf_idm_multi_package(subcomponent) result(multi_package)
    character(len=*), intent(in) :: subcomponent
    logical :: multi_package
    select case (subcomponent)
    case ('NAM')
      multi_package = chf_nam_multi_package
    case ('DISV1D')
      multi_package = chf_disv1d_multi_package
    case ('CXS')
      multi_package = chf_cxs_multi_package
    case ('DFW')
      multi_package = chf_dfw_multi_package
    case ('IC')
      multi_package = chf_ic_multi_package
    case ('OC')
      multi_package = chf_oc_multi_package
    case ('CDB')
      multi_package = chf_cdb_multi_package
    case ('CHD')
      multi_package = chf_chd_multi_package
    case ('FLW')
      multi_package = chf_flw_multi_package
    case ('PCP')
      multi_package = chf_pcp_multi_package
    case ('EVP')
      multi_package = chf_evp_multi_package
    case ('STO')
      multi_package = chf_sto_multi_package
    case ('ZDG')
      multi_package = chf_zdg_multi_package
    case default
      call store_error('Idm selector subcomponent not found; '//&
                       &'component="CHF"'//&
                       &', subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function chf_idm_multi_package

  function chf_idm_subpackages(subcomponent) result(subpackages)
    character(len=*), intent(in) :: subcomponent
    character(len=16), dimension(:), pointer :: subpackages
    select case (subcomponent)
    case ('NAM')
      call set_subpkg_pointer(subpackages, chf_nam_subpackages)
    case ('DISV1D')
      call set_subpkg_pointer(subpackages, chf_disv1d_subpackages)
    case ('CXS')
      call set_subpkg_pointer(subpackages, chf_cxs_subpackages)
    case ('DFW')
      call set_subpkg_pointer(subpackages, chf_dfw_subpackages)
    case ('IC')
      call set_subpkg_pointer(subpackages, chf_ic_subpackages)
    case ('OC')
      call set_subpkg_pointer(subpackages, chf_oc_subpackages)
    case ('CDB')
      call set_subpkg_pointer(subpackages, chf_cdb_subpackages)
    case ('CHD')
      call set_subpkg_pointer(subpackages, chf_chd_subpackages)
    case ('FLW')
      call set_subpkg_pointer(subpackages, chf_flw_subpackages)
    case ('PCP')
      call set_subpkg_pointer(subpackages, chf_pcp_subpackages)
    case ('EVP')
      call set_subpkg_pointer(subpackages, chf_evp_subpackages)
    case ('STO')
      call set_subpkg_pointer(subpackages, chf_sto_subpackages)
    case ('ZDG')
      call set_subpkg_pointer(subpackages, chf_zdg_subpackages)
    case default
    end select
    return
  end function chf_idm_subpackages

  function chf_idm_integrated(subcomponent) result(integrated)
    character(len=*), intent(in) :: subcomponent
    logical :: integrated
    integrated = .false.
    select case (subcomponent)
    case ('NAM')
      integrated = .true.
    case ('DISV1D')
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
  end function chf_idm_integrated

end module IdmChfDfnSelectorModule
