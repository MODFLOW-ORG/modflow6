! ** Do Not Modify! MODFLOW 6 system generated file. **
module IdmGwfDfnSelectorModule

  use ConstantsModule, only: LENVARNAME
  use SimModule, only: store_error
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  use GwfNamInputModule
  use GwfChdInputModule
  use GwfCsubInputModule
  use GwfDisInputModule
  use GwfDisuInputModule
  use GwfDisvInputModule
  use GwfDrnInputModule
  use GwfDrngInputModule
  use GwfEvtInputModule
  use GwfEvtaInputModule
  use GwfGhbInputModule
  use GwfGhbgInputModule
  use GwfIcInputModule
  use GwfNpfInputModule
  use GwfRchInputModule
  use GwfRchaInputModule
  use GwfRivInputModule
  use GwfStoInputModule
  use GwfWelInputModule

  implicit none
  private
  public :: gwf_param_definitions
  public :: gwf_aggregate_definitions
  public :: gwf_block_definitions
  public :: gwf_idm_multi_package
  public :: gwf_idm_subpackages
  public :: gwf_idm_integrated

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

  function gwf_param_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, gwf_nam_param_definitions)
    case ('CHD')
      call set_param_pointer(input_definition, gwf_chd_param_definitions)
    case ('CSUB')
      call set_param_pointer(input_definition, gwf_csub_param_definitions)
    case ('DIS')
      call set_param_pointer(input_definition, gwf_dis_param_definitions)
    case ('DISU')
      call set_param_pointer(input_definition, gwf_disu_param_definitions)
    case ('DISV')
      call set_param_pointer(input_definition, gwf_disv_param_definitions)
    case ('DRN')
      call set_param_pointer(input_definition, gwf_drn_param_definitions)
    case ('DRNG')
      call set_param_pointer(input_definition, gwf_drng_param_definitions)
    case ('EVT')
      call set_param_pointer(input_definition, gwf_evt_param_definitions)
    case ('EVTA')
      call set_param_pointer(input_definition, gwf_evta_param_definitions)
    case ('GHB')
      call set_param_pointer(input_definition, gwf_ghb_param_definitions)
    case ('GHBG')
      call set_param_pointer(input_definition, gwf_ghbg_param_definitions)
    case ('IC')
      call set_param_pointer(input_definition, gwf_ic_param_definitions)
    case ('NPF')
      call set_param_pointer(input_definition, gwf_npf_param_definitions)
    case ('RCH')
      call set_param_pointer(input_definition, gwf_rch_param_definitions)
    case ('RCHA')
      call set_param_pointer(input_definition, gwf_rcha_param_definitions)
    case ('RIV')
      call set_param_pointer(input_definition, gwf_riv_param_definitions)
    case ('STO')
      call set_param_pointer(input_definition, gwf_sto_param_definitions)
    case ('WEL')
      call set_param_pointer(input_definition, gwf_wel_param_definitions)
    case default
    end select
    return
  end function gwf_param_definitions

  function gwf_aggregate_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, gwf_nam_aggregate_definitions)
    case ('CHD')
      call set_param_pointer(input_definition, gwf_chd_aggregate_definitions)
    case ('CSUB')
      call set_param_pointer(input_definition, gwf_csub_aggregate_definitions)
    case ('DIS')
      call set_param_pointer(input_definition, gwf_dis_aggregate_definitions)
    case ('DISU')
      call set_param_pointer(input_definition, gwf_disu_aggregate_definitions)
    case ('DISV')
      call set_param_pointer(input_definition, gwf_disv_aggregate_definitions)
    case ('DRN')
      call set_param_pointer(input_definition, gwf_drn_aggregate_definitions)
    case ('DRNG')
      call set_param_pointer(input_definition, gwf_drng_aggregate_definitions)
    case ('EVT')
      call set_param_pointer(input_definition, gwf_evt_aggregate_definitions)
    case ('EVTA')
      call set_param_pointer(input_definition, gwf_evta_aggregate_definitions)
    case ('GHB')
      call set_param_pointer(input_definition, gwf_ghb_aggregate_definitions)
    case ('GHBG')
      call set_param_pointer(input_definition, gwf_ghbg_aggregate_definitions)
    case ('IC')
      call set_param_pointer(input_definition, gwf_ic_aggregate_definitions)
    case ('NPF')
      call set_param_pointer(input_definition, gwf_npf_aggregate_definitions)
    case ('RCH')
      call set_param_pointer(input_definition, gwf_rch_aggregate_definitions)
    case ('RCHA')
      call set_param_pointer(input_definition, gwf_rcha_aggregate_definitions)
    case ('RIV')
      call set_param_pointer(input_definition, gwf_riv_aggregate_definitions)
    case ('STO')
      call set_param_pointer(input_definition, gwf_sto_aggregate_definitions)
    case ('WEL')
      call set_param_pointer(input_definition, gwf_wel_aggregate_definitions)
    case default
    end select
    return
  end function gwf_aggregate_definitions

  function gwf_block_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputBlockDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_block_pointer(input_definition, gwf_nam_block_definitions)
    case ('CHD')
      call set_block_pointer(input_definition, gwf_chd_block_definitions)
    case ('CSUB')
      call set_block_pointer(input_definition, gwf_csub_block_definitions)
    case ('DIS')
      call set_block_pointer(input_definition, gwf_dis_block_definitions)
    case ('DISU')
      call set_block_pointer(input_definition, gwf_disu_block_definitions)
    case ('DISV')
      call set_block_pointer(input_definition, gwf_disv_block_definitions)
    case ('DRN')
      call set_block_pointer(input_definition, gwf_drn_block_definitions)
    case ('DRNG')
      call set_block_pointer(input_definition, gwf_drng_block_definitions)
    case ('EVT')
      call set_block_pointer(input_definition, gwf_evt_block_definitions)
    case ('EVTA')
      call set_block_pointer(input_definition, gwf_evta_block_definitions)
    case ('GHB')
      call set_block_pointer(input_definition, gwf_ghb_block_definitions)
    case ('GHBG')
      call set_block_pointer(input_definition, gwf_ghbg_block_definitions)
    case ('IC')
      call set_block_pointer(input_definition, gwf_ic_block_definitions)
    case ('NPF')
      call set_block_pointer(input_definition, gwf_npf_block_definitions)
    case ('RCH')
      call set_block_pointer(input_definition, gwf_rch_block_definitions)
    case ('RCHA')
      call set_block_pointer(input_definition, gwf_rcha_block_definitions)
    case ('RIV')
      call set_block_pointer(input_definition, gwf_riv_block_definitions)
    case ('STO')
      call set_block_pointer(input_definition, gwf_sto_block_definitions)
    case ('WEL')
      call set_block_pointer(input_definition, gwf_wel_block_definitions)
    case default
    end select
    return
  end function gwf_block_definitions

  function gwf_idm_multi_package(subcomponent) result(multi_package)
    character(len=*), intent(in) :: subcomponent
    logical :: multi_package
    select case (subcomponent)
    case ('NAM')
      multi_package = gwf_nam_multi_package
    case ('CHD')
      multi_package = gwf_chd_multi_package
    case ('CSUB')
      multi_package = gwf_csub_multi_package
    case ('DIS')
      multi_package = gwf_dis_multi_package
    case ('DISU')
      multi_package = gwf_disu_multi_package
    case ('DISV')
      multi_package = gwf_disv_multi_package
    case ('DRN')
      multi_package = gwf_drn_multi_package
    case ('DRNG')
      multi_package = gwf_drng_multi_package
    case ('EVT')
      multi_package = gwf_evt_multi_package
    case ('EVTA')
      multi_package = gwf_evta_multi_package
    case ('GHB')
      multi_package = gwf_ghb_multi_package
    case ('GHBG')
      multi_package = gwf_ghbg_multi_package
    case ('IC')
      multi_package = gwf_ic_multi_package
    case ('NPF')
      multi_package = gwf_npf_multi_package
    case ('RCH')
      multi_package = gwf_rch_multi_package
    case ('RCHA')
      multi_package = gwf_rcha_multi_package
    case ('RIV')
      multi_package = gwf_riv_multi_package
    case ('STO')
      multi_package = gwf_sto_multi_package
    case ('WEL')
      multi_package = gwf_wel_multi_package
    case default
      call store_error('Idm selector subcomponent not found; '//&
                       &'component="GWF"'//&
                       &', subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function gwf_idm_multi_package

  function gwf_idm_subpackages(subcomponent) result(subpackages)
    character(len=*), intent(in) :: subcomponent
    character(len=16), dimension(:), pointer :: subpackages
    select case (subcomponent)
    case ('NAM')
      call set_subpkg_pointer(subpackages, gwf_nam_subpackages)
    case ('CHD')
      call set_subpkg_pointer(subpackages, gwf_chd_subpackages)
    case ('CSUB')
      call set_subpkg_pointer(subpackages, gwf_csub_subpackages)
    case ('DIS')
      call set_subpkg_pointer(subpackages, gwf_dis_subpackages)
    case ('DISU')
      call set_subpkg_pointer(subpackages, gwf_disu_subpackages)
    case ('DISV')
      call set_subpkg_pointer(subpackages, gwf_disv_subpackages)
    case ('DRN')
      call set_subpkg_pointer(subpackages, gwf_drn_subpackages)
    case ('DRNG')
      call set_subpkg_pointer(subpackages, gwf_drng_subpackages)
    case ('EVT')
      call set_subpkg_pointer(subpackages, gwf_evt_subpackages)
    case ('EVTA')
      call set_subpkg_pointer(subpackages, gwf_evta_subpackages)
    case ('GHB')
      call set_subpkg_pointer(subpackages, gwf_ghb_subpackages)
    case ('GHBG')
      call set_subpkg_pointer(subpackages, gwf_ghbg_subpackages)
    case ('IC')
      call set_subpkg_pointer(subpackages, gwf_ic_subpackages)
    case ('NPF')
      call set_subpkg_pointer(subpackages, gwf_npf_subpackages)
    case ('RCH')
      call set_subpkg_pointer(subpackages, gwf_rch_subpackages)
    case ('RCHA')
      call set_subpkg_pointer(subpackages, gwf_rcha_subpackages)
    case ('RIV')
      call set_subpkg_pointer(subpackages, gwf_riv_subpackages)
    case ('STO')
      call set_subpkg_pointer(subpackages, gwf_sto_subpackages)
    case ('WEL')
      call set_subpkg_pointer(subpackages, gwf_wel_subpackages)
    case default
    end select
    return
  end function gwf_idm_subpackages

  function gwf_idm_integrated(subcomponent) result(integrated)
    character(len=*), intent(in) :: subcomponent
    logical :: integrated
    integrated = .false.
    select case (subcomponent)
    case ('NAM')
      integrated = .true.
    case ('CHD')
      integrated = .true.
    case ('CSUB')
      integrated = .true.
    case ('DIS')
      integrated = .true.
    case ('DISU')
      integrated = .true.
    case ('DISV')
      integrated = .true.
    case ('DRN')
      integrated = .true.
    case ('DRNG')
      integrated = .true.
    case ('EVT')
      integrated = .true.
    case ('EVTA')
      integrated = .true.
    case ('GHB')
      integrated = .true.
    case ('GHBG')
      integrated = .true.
    case ('IC')
      integrated = .true.
    case ('NPF')
      integrated = .true.
    case ('RCH')
      integrated = .true.
    case ('RCHA')
      integrated = .true.
    case ('RIV')
      integrated = .true.
    case ('STO')
      integrated = .true.
    case ('WEL')
      integrated = .true.
    case default
    end select
    return
  end function gwf_idm_integrated

end module IdmGwfDfnSelectorModule
