! ** Do Not Modify! MODFLOW 6 system generated file. **
module IdmGwtDfnSelectorModule

  use ConstantsModule, only: LENVARNAME
  use SimModule, only: store_error
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  use GwtNamInputModule
  use GwtDisInputModule
  use GwtDisuInputModule
  use GwtDisvInputModule
  use GwtDspInputModule
  use GwtCncInputModule
  use GwtIcInputModule
  use GwtIstInputModule
  use GwtMstInputModule

  implicit none
  private
  public :: gwt_param_definitions
  public :: gwt_aggregate_definitions
  public :: gwt_block_definitions
  public :: gwt_idm_multi_package
  public :: gwt_idm_subpackages
  public :: gwt_idm_integrated

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

  function gwt_param_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, gwt_nam_param_definitions)
    case ('DIS')
      call set_param_pointer(input_definition, gwt_dis_param_definitions)
    case ('DISU')
      call set_param_pointer(input_definition, gwt_disu_param_definitions)
    case ('DISV')
      call set_param_pointer(input_definition, gwt_disv_param_definitions)
    case ('DSP')
      call set_param_pointer(input_definition, gwt_dsp_param_definitions)
    case ('CNC')
      call set_param_pointer(input_definition, gwt_cnc_param_definitions)
    case ('IC')
      call set_param_pointer(input_definition, gwt_ic_param_definitions)
    case ('IST')
      call set_param_pointer(input_definition, gwt_ist_param_definitions)
    case ('MST')
      call set_param_pointer(input_definition, gwt_mst_param_definitions)
    case default
    end select
    return
  end function gwt_param_definitions

  function gwt_aggregate_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_param_pointer(input_definition, gwt_nam_aggregate_definitions)
    case ('DIS')
      call set_param_pointer(input_definition, gwt_dis_aggregate_definitions)
    case ('DISU')
      call set_param_pointer(input_definition, gwt_disu_aggregate_definitions)
    case ('DISV')
      call set_param_pointer(input_definition, gwt_disv_aggregate_definitions)
    case ('DSP')
      call set_param_pointer(input_definition, gwt_dsp_aggregate_definitions)
    case ('CNC')
      call set_param_pointer(input_definition, gwt_cnc_aggregate_definitions)
    case ('IC')
      call set_param_pointer(input_definition, gwt_ic_aggregate_definitions)
    case ('IST')
      call set_param_pointer(input_definition, gwt_ist_aggregate_definitions)
    case ('MST')
      call set_param_pointer(input_definition, gwt_mst_aggregate_definitions)
    case default
    end select
    return
  end function gwt_aggregate_definitions

  function gwt_block_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputBlockDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('NAM')
      call set_block_pointer(input_definition, gwt_nam_block_definitions)
    case ('DIS')
      call set_block_pointer(input_definition, gwt_dis_block_definitions)
    case ('DISU')
      call set_block_pointer(input_definition, gwt_disu_block_definitions)
    case ('DISV')
      call set_block_pointer(input_definition, gwt_disv_block_definitions)
    case ('DSP')
      call set_block_pointer(input_definition, gwt_dsp_block_definitions)
    case ('CNC')
      call set_block_pointer(input_definition, gwt_cnc_block_definitions)
    case ('IC')
      call set_block_pointer(input_definition, gwt_ic_block_definitions)
    case ('IST')
      call set_block_pointer(input_definition, gwt_ist_block_definitions)
    case ('MST')
      call set_block_pointer(input_definition, gwt_mst_block_definitions)
    case default
    end select
    return
  end function gwt_block_definitions

  function gwt_idm_multi_package(subcomponent) result(multi_package)
    character(len=*), intent(in) :: subcomponent
    logical :: multi_package
    select case (subcomponent)
    case ('NAM')
      multi_package = gwt_nam_multi_package
    case ('DIS')
      multi_package = gwt_dis_multi_package
    case ('DISU')
      multi_package = gwt_disu_multi_package
    case ('DISV')
      multi_package = gwt_disv_multi_package
    case ('DSP')
      multi_package = gwt_dsp_multi_package
    case ('CNC')
      multi_package = gwt_cnc_multi_package
    case ('IC')
      multi_package = gwt_ic_multi_package
    case ('IST')
      multi_package = gwt_ist_multi_package
    case ('MST')
      multi_package = gwt_mst_multi_package
    case default
      call store_error('Idm selector subcomponent not found; '//&
                       &'component="GWT"'//&
                       &', subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function gwt_idm_multi_package

  function gwt_idm_subpackages(subcomponent) result(subpackages)
    character(len=*), intent(in) :: subcomponent
    character(len=16), dimension(:), pointer :: subpackages
    select case (subcomponent)
    case ('NAM')
      call set_subpkg_pointer(subpackages, gwt_nam_subpackages)
    case ('DIS')
      call set_subpkg_pointer(subpackages, gwt_dis_subpackages)
    case ('DISU')
      call set_subpkg_pointer(subpackages, gwt_disu_subpackages)
    case ('DISV')
      call set_subpkg_pointer(subpackages, gwt_disv_subpackages)
    case ('DSP')
      call set_subpkg_pointer(subpackages, gwt_dsp_subpackages)
    case ('CNC')
      call set_subpkg_pointer(subpackages, gwt_cnc_subpackages)
    case ('IC')
      call set_subpkg_pointer(subpackages, gwt_ic_subpackages)
    case ('IST')
      call set_subpkg_pointer(subpackages, gwt_ist_subpackages)
    case ('MST')
      call set_subpkg_pointer(subpackages, gwt_mst_subpackages)
    case default
    end select
    return
  end function gwt_idm_subpackages

  function gwt_idm_integrated(subcomponent) result(integrated)
    character(len=*), intent(in) :: subcomponent
    logical :: integrated
    integrated = .false.
    select case (subcomponent)
    case ('NAM')
      integrated = .true.
    case ('DIS')
      integrated = .true.
    case ('DISU')
      integrated = .true.
    case ('DISV')
      integrated = .true.
    case ('DSP')
      integrated = .true.
    case ('CNC')
      integrated = .true.
    case ('IC')
      integrated = .true.
    case ('IST')
      integrated = .true.
    case ('MST')
      integrated = .true.
    case default
    end select
    return
  end function gwt_idm_integrated

end module IdmGwtDfnSelectorModule
