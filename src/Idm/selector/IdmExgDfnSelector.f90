! ** Do Not Modify! MODFLOW 6 system generated file. **
module IdmExgDfnSelectorModule

  use ConstantsModule, only: LENVARNAME
  use SimModule, only: store_error
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  use ExgChfgwfInputModule
  use ExgGwfgwfInputModule
  use ExgGwfgwtInputModule
  use ExgGwtgwtInputModule
  use ExgGwfgweInputModule
  use ExgGwegweInputModule
  use ExgGwfprtInputModule
  use ExgOlfgwfInputModule

  implicit none
  private
  public :: exg_param_definitions
  public :: exg_aggregate_definitions
  public :: exg_block_definitions
  public :: exg_idm_multi_package
  public :: exg_idm_subpackages
  public :: exg_idm_integrated

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

  function exg_param_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('CHFGWF')
      call set_param_pointer(input_definition, exg_chfgwf_param_definitions)
    case ('GWFGWF')
      call set_param_pointer(input_definition, exg_gwfgwf_param_definitions)
    case ('GWFGWT')
      call set_param_pointer(input_definition, exg_gwfgwt_param_definitions)
    case ('GWTGWT')
      call set_param_pointer(input_definition, exg_gwtgwt_param_definitions)
    case ('GWFGWE')
      call set_param_pointer(input_definition, exg_gwfgwe_param_definitions)
    case ('GWEGWE')
      call set_param_pointer(input_definition, exg_gwegwe_param_definitions)
    case ('GWFPRT')
      call set_param_pointer(input_definition, exg_gwfprt_param_definitions)
    case ('OLFGWF')
      call set_param_pointer(input_definition, exg_olfgwf_param_definitions)
    case default
    end select
    return
  end function exg_param_definitions

  function exg_aggregate_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputParamDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('CHFGWF')
      call set_param_pointer(input_definition, exg_chfgwf_aggregate_definitions)
    case ('GWFGWF')
      call set_param_pointer(input_definition, exg_gwfgwf_aggregate_definitions)
    case ('GWFGWT')
      call set_param_pointer(input_definition, exg_gwfgwt_aggregate_definitions)
    case ('GWTGWT')
      call set_param_pointer(input_definition, exg_gwtgwt_aggregate_definitions)
    case ('GWFGWE')
      call set_param_pointer(input_definition, exg_gwfgwe_aggregate_definitions)
    case ('GWEGWE')
      call set_param_pointer(input_definition, exg_gwegwe_aggregate_definitions)
    case ('GWFPRT')
      call set_param_pointer(input_definition, exg_gwfprt_aggregate_definitions)
    case ('OLFGWF')
      call set_param_pointer(input_definition, exg_olfgwf_aggregate_definitions)
    case default
    end select
    return
  end function exg_aggregate_definitions

  function exg_block_definitions(subcomponent) result(input_definition)
    character(len=*), intent(in) :: subcomponent
    type(InputBlockDefinitionType), dimension(:), pointer :: input_definition
    nullify (input_definition)
    select case (subcomponent)
    case ('CHFGWF')
      call set_block_pointer(input_definition, exg_chfgwf_block_definitions)
    case ('GWFGWF')
      call set_block_pointer(input_definition, exg_gwfgwf_block_definitions)
    case ('GWFGWT')
      call set_block_pointer(input_definition, exg_gwfgwt_block_definitions)
    case ('GWTGWT')
      call set_block_pointer(input_definition, exg_gwtgwt_block_definitions)
    case ('GWFGWE')
      call set_block_pointer(input_definition, exg_gwfgwe_block_definitions)
    case ('GWEGWE')
      call set_block_pointer(input_definition, exg_gwegwe_block_definitions)
    case ('GWFPRT')
      call set_block_pointer(input_definition, exg_gwfprt_block_definitions)
    case ('OLFGWF')
      call set_block_pointer(input_definition, exg_olfgwf_block_definitions)
    case default
    end select
    return
  end function exg_block_definitions

  function exg_idm_multi_package(subcomponent) result(multi_package)
    character(len=*), intent(in) :: subcomponent
    logical :: multi_package
    select case (subcomponent)
    case ('CHFGWF')
      multi_package = exg_chfgwf_multi_package
    case ('GWFGWF')
      multi_package = exg_gwfgwf_multi_package
    case ('GWFGWT')
      multi_package = exg_gwfgwt_multi_package
    case ('GWTGWT')
      multi_package = exg_gwtgwt_multi_package
    case ('GWFGWE')
      multi_package = exg_gwfgwe_multi_package
    case ('GWEGWE')
      multi_package = exg_gwegwe_multi_package
    case ('GWFPRT')
      multi_package = exg_gwfprt_multi_package
    case ('OLFGWF')
      multi_package = exg_olfgwf_multi_package
    case default
      call store_error('Idm selector subcomponent not found; '//&
                       &'component="EXG"'//&
                       &', subcomponent="'//trim(subcomponent)//'".', .true.)
    end select
    return
  end function exg_idm_multi_package

  function exg_idm_subpackages(subcomponent) result(subpackages)
    character(len=*), intent(in) :: subcomponent
    character(len=16), dimension(:), pointer :: subpackages
    select case (subcomponent)
    case ('CHFGWF')
      call set_subpkg_pointer(subpackages, exg_chfgwf_subpackages)
    case ('GWFGWF')
      call set_subpkg_pointer(subpackages, exg_gwfgwf_subpackages)
    case ('GWFGWT')
      call set_subpkg_pointer(subpackages, exg_gwfgwt_subpackages)
    case ('GWTGWT')
      call set_subpkg_pointer(subpackages, exg_gwtgwt_subpackages)
    case ('GWFGWE')
      call set_subpkg_pointer(subpackages, exg_gwfgwe_subpackages)
    case ('GWEGWE')
      call set_subpkg_pointer(subpackages, exg_gwegwe_subpackages)
    case ('GWFPRT')
      call set_subpkg_pointer(subpackages, exg_gwfprt_subpackages)
    case ('OLFGWF')
      call set_subpkg_pointer(subpackages, exg_olfgwf_subpackages)
    case default
    end select
    return
  end function exg_idm_subpackages

  function exg_idm_integrated(subcomponent) result(integrated)
    character(len=*), intent(in) :: subcomponent
    logical :: integrated
    integrated = .false.
    select case (subcomponent)
    case ('CHFGWF')
      integrated = .true.
    case ('GWFGWF')
      integrated = .true.
    case ('GWFGWT')
      integrated = .true.
    case ('GWTGWT')
      integrated = .true.
    case ('GWFGWE')
      integrated = .true.
    case ('GWEGWE')
      integrated = .true.
    case ('GWFPRT')
      integrated = .true.
    case ('OLFGWF')
      integrated = .true.
    case default
    end select
    return
  end function exg_idm_integrated

end module IdmExgDfnSelectorModule
