! ** Do Not Modify! MODFLOW 6 system generated file. **
module IdmOlfDfnSelectorModule

  use ConstantsModule, only: LENVARNAME
  use SimModule, only: store_error
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  use OlfNamInputModule
  use OlfDis2DInputModule
  use OlfDisv2DInputModule
  use OlfDfwInputModule
  use OlfIcInputModule
  use OlfOcInputModule
  use OlfCdbInputModule
  use OlfChdInputModule
  use OlfFlwInputModule
  use OlfPcpInputModule
  use OlfEvpInputModule
  use OlfStoInputModule
  use OlfZdgInputModule

  implicit none
  private
  public :: olf_param_definitions
  public :: olf_aggregate_definitions
  public :: olf_block_definitions
  public :: olf_idm_multi_package
  public :: olf_idm_subpackages
  public :: olf_idm_integrated

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

  function olf_param_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, olf_nam_param_definitions)
    case ('DIS2D')
      call set_param_pointer(input_definition, olf_dis2d_param_definitions)
    case ('DISV2D')
      call set_param_pointer(input_definition, olf_disv2d_param_definitions)
    case ('DFW')
      call set_param_pointer(input_definition, olf_dfw_param_definitions)
    case ('IC')
      call set_param_pointer(input_definition, olf_ic_param_definitions)
    case ('OC')
      call set_param_pointer(input_definition, olf_oc_param_definitions)
    case ('CDB')
      call set_param_pointer(input_definition, olf_cdb_param_definitions)
    case ('CHD')
      call set_param_pointer(input_definition, olf_chd_param_definitions)
    case ('FLW')
      call set_param_pointer(input_definition, olf_flw_param_definitions)
    case ('PCP')
      call set_param_pointer(input_definition, olf_pcp_param_definitions)
    case ('EVP')
      call set_param_pointer(input_definition, olf_evp_param_definitions)
    case ('STO')
      call set_param_pointer(input_definition, olf_sto_param_definitions)
    case ('ZDG')
      call set_param_pointer(input_definition, olf_zdg_param_definitions)
    case default
    end select
    return
  end function olf_param_definitions

  function olf_aggregate_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, olf_nam_aggregate_definitions)
    case ('DIS2D')
      call set_param_pointer(input_definition, olf_dis2d_aggregate_definitions)
    case ('DISV2D')
      call set_param_pointer(input_definition, olf_disv2d_aggregate_definitions)
    case ('DFW')
      call set_param_pointer(input_definition, olf_dfw_aggregate_definitions)
    case ('IC')
      call set_param_pointer(input_definition, olf_ic_aggregate_definitions)
    case ('OC')
      call set_param_pointer(input_definition, olf_oc_aggregate_definitions)
    case ('CDB')
      call set_param_pointer(input_definition, olf_cdb_aggregate_definitions)
    case ('CHD')
      call set_param_pointer(input_definition, olf_chd_aggregate_definitions)
    case ('FLW')
      call set_param_pointer(input_definition, olf_flw_aggregate_definitions)
    case ('PCP')
      call set_param_pointer(input_definition, olf_pcp_aggregate_definitions)
    case ('EVP')
      call set_param_pointer(input_definition, olf_evp_aggregate_definitions)
    case ('STO')
      call set_param_pointer(input_definition, olf_sto_aggregate_definitions)
    case ('ZDG')
      call set_param_pointer(input_definition, olf_zdg_aggregate_definitions)
    case default
    end select
    return
  end function olf_aggregate_definitions

  function olf_block_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputBlockDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_block_pointer(input_definition, olf_nam_block_definitions)
    case ('DIS2D')
      call set_block_pointer(input_definition, olf_dis2d_block_definitions)
    case ('DISV2D')
      call set_block_pointer(input_definition, olf_disv2d_block_definitions)
    case ('DFW')
      call set_block_pointer(input_definition, olf_dfw_block_definitions)
    case ('IC')
      call set_block_pointer(input_definition, olf_ic_block_definitions)
    case ('OC')
      call set_block_pointer(input_definition, olf_oc_block_definitions)
    case ('CDB')
      call set_block_pointer(input_definition, olf_cdb_block_definitions)
    case ('CHD')
      call set_block_pointer(input_definition, olf_chd_block_definitions)
    case ('FLW')
      call set_block_pointer(input_definition, olf_flw_block_definitions)
    case ('PCP')
      call set_block_pointer(input_definition, olf_pcp_block_definitions)
    case ('EVP')
      call set_block_pointer(input_definition, olf_evp_block_definitions)
    case ('STO')
      call set_block_pointer(input_definition, olf_sto_block_definitions)
    case ('ZDG')
      call set_block_pointer(input_definition, olf_zdg_block_definitions)
    case default
    end select
    return
  end function olf_block_definitions

  function olf_idm_multi_package(subcomponent) result(multi_package)
    character(len=*), intent(in) :: subcomponent
    logical :: multi_package
    select case (subcomponent)
    case ('NAM')
      multi_package = olf_nam_multi_package
    case ('DIS2D')
      multi_package = olf_dis2d_multi_package
    case ('DISV2D')
      multi_package = olf_disv2d_multi_package
    case ('DFW')
      multi_package = olf_dfw_multi_package
    case ('IC')
      multi_package = olf_ic_multi_package
    case ('OC')
      multi_package = olf_oc_multi_package
    case ('CDB')
      multi_package = olf_cdb_multi_package
    case ('CHD')
      multi_package = olf_chd_multi_package
    case ('FLW')
      multi_package = olf_flw_multi_package
    case ('PCP')
      multi_package = olf_pcp_multi_package
    case ('EVP')
      multi_package = olf_evp_multi_package
    case ('STO')
      multi_package = olf_sto_multi_package
    case ('ZDG')
      multi_package = olf_zdg_multi_package
    case default
      call store_error('Idm selector subcomponent not found; '//&
                       &'component="OLF"'//&
                       &', subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function olf_idm_multi_package

  function olf_idm_subpackages(subcomponent) result(subpackages)
    character(len=*), intent(in) :: subcomponent
    character(len=16), dimension(:), pointer :: subpackages
    select case (subcomponent)
    case ('NAM')
      call set_subpkg_pointer(subpackages, olf_nam_subpackages)
    case ('DIS2D')
      call set_subpkg_pointer(subpackages, olf_dis2d_subpackages)
    case ('DISV2D')
      call set_subpkg_pointer(subpackages, olf_disv2d_subpackages)
    case ('DFW')
      call set_subpkg_pointer(subpackages, olf_dfw_subpackages)
    case ('IC')
      call set_subpkg_pointer(subpackages, olf_ic_subpackages)
    case ('OC')
      call set_subpkg_pointer(subpackages, olf_oc_subpackages)
    case ('CDB')
      call set_subpkg_pointer(subpackages, olf_cdb_subpackages)
    case ('CHD')
      call set_subpkg_pointer(subpackages, olf_chd_subpackages)
    case ('FLW')
      call set_subpkg_pointer(subpackages, olf_flw_subpackages)
    case ('PCP')
      call set_subpkg_pointer(subpackages, olf_pcp_subpackages)
    case ('EVP')
      call set_subpkg_pointer(subpackages, olf_evp_subpackages)
    case ('STO')
      call set_subpkg_pointer(subpackages, olf_sto_subpackages)
    case ('ZDG')
      call set_subpkg_pointer(subpackages, olf_zdg_subpackages)
    case default
    end select
    return
  end function olf_idm_subpackages

  function olf_idm_integrated(subcomponent) result(integrated)
    character(len=*), intent(in) :: subcomponent
    logical :: integrated
    integrated = .false.
    select case (subcomponent)
    case ('NAM')
      integrated = .true.
    case ('DIS2D')
      integrated = .true.
    case ('DISV2D')
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
  end function olf_idm_integrated

end module IdmOlfDfnSelectorModule
