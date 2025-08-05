!> @brief This module contains the IdmLoadModule
!!
!! This module contains routines for managing static
!! and dynamic input loading for supported sources.
!!
!<
module IdmLoadModule

  use KindModule, only: DP, I4B, LGP
  use SimVariablesModule, only: errmsg
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, LENMODELNAME, &
                             LENEXCHANGENAME, LENCOMPONENTNAME
  use SimModule, only: store_error, count_errors, store_error_filename
  use ListModule, only: ListType
  use InputLoadTypeModule, only: StaticPkgLoadBaseType, &
                                 DynamicPkgLoadBaseType, &
                                 ModelDynamicPkgsType, &
                                 model_inputs
  use InputDefinitionModule, only: InputParamDefinitionType
  use ModflowInputModule, only: ModflowInputType, getModflowInput

  implicit none
  private
  public :: simnam_load
  public :: simtdis_load
  public :: load_models
  public :: load_exchanges
  public :: idm_df
  public :: idm_rp
  public :: idm_ad
  public :: idm_da

contains

  !> @brief advance package dynamic data for period steps
  !<
  subroutine idm_df()
    use InputLoadTypeModule, only: GetDynamicModelFromList
    class(ModelDynamicPkgsType), pointer :: model_dynamic_input
    integer(I4B) :: n
    do n = 1, model_inputs%Count()
      model_dynamic_input => GetDynamicModelFromList(model_inputs, n)
      call model_dynamic_input%df()
    end do
  end subroutine idm_df

  !> @brief load package dynamic data for period
  !<
  subroutine idm_rp()
    use InputLoadTypeModule, only: GetDynamicModelFromList
    class(ModelDynamicPkgsType), pointer :: model_dynamic_input
    integer(I4B) :: n
    do n = 1, model_inputs%Count()
      model_dynamic_input => GetDynamicModelFromList(model_inputs, n)
      call model_dynamic_input%rp()
    end do
  end subroutine idm_rp

  !> @brief advance package dynamic data for period steps
  !<
  subroutine idm_ad()
    use InputLoadTypeModule, only: GetDynamicModelFromList
    class(ModelDynamicPkgsType), pointer :: model_dynamic_input
    integer(I4B) :: n
    do n = 1, model_inputs%Count()
      model_dynamic_input => GetDynamicModelFromList(model_inputs, n)
      call model_dynamic_input%ad()
    end do
  end subroutine idm_ad

  !> @brief idm deallocate routine
  !<
  subroutine idm_da(iout)
    use SimVariablesModule, only: idm_context
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path, split_mem_path
    use MemoryManagerExtModule, only: memorystore_remove
    use CharacterStringModule, only: CharacterStringType
    integer(I4B), intent(in) :: iout
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mempaths
    character(len=LENCOMPONENTNAME) :: exg_comp, exg_subcomp
    character(len=LENMEMPATH) :: input_mempath, mempath
    integer(I4B) :: n

    ! deallocate dynamic loaders
    call dynamic_da(iout)

    ! deallocate EXG mempaths
    input_mempath = create_mem_path('SIM', 'NAM', idm_context)
    call mem_setptr(mempaths, 'EXGMEMPATHS', input_mempath)
    do n = 1, size(mempaths)
      mempath = mempaths(n)
      if (mempath /= '') then
        call split_mem_path(mempath, exg_comp, exg_subcomp)
        call memorystore_remove(exg_comp, exg_subcomp, idm_context)
      end if
    end do

    ! deallocate input context SIM paths
    call memorystore_remove('UTL', 'HPC', idm_context)
    call memorystore_remove('SIM', 'TDIS', idm_context)
    call memorystore_remove('SIM', 'NAM', idm_context)
    call memorystore_remove(component='SIM', context=idm_context)
  end subroutine idm_da

  !> @brief load an integrated model package from supported source
  !<
  recursive subroutine input_load(component_type, subcomponent_type, modelname, &
                                  pkgname, pkgtype, filename, modelfname, &
                                  nc_vars, iout)
    use ModelPackageInputsModule, only: ModelPackageInputsType
    use NCFileVarsModule, only: NCFileVarsType
    use SourceLoadModule, only: create_input_loader
    character(len=*), intent(in) :: component_type
    character(len=*), intent(in) :: subcomponent_type
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: pkgtype
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: modelfname
    type(NCFileVarsType), pointer, intent(in) :: nc_vars
    integer(I4B), intent(in) :: iout
    class(StaticPkgLoadBaseType), pointer :: static_loader
    class(DynamicPkgLoadBaseType), pointer :: dynamic_loader
    class(ModelDynamicPkgsType), pointer :: dynamic_model
    integer(I4B) :: n

    ! create model package loader
    static_loader => &
      create_input_loader(component_type, subcomponent_type, modelname, pkgname, &
                          pkgtype, 'MODEL', filename, modelfname, nc_vars)

    ! load static input and set dynamic loader
    dynamic_loader => static_loader%load(iout)

    ! set pointer to model dynamic packages list
    dynamic_model => &
      dynamic_models(static_loader%mf6_input%component_type, modelname, &
                     static_loader%component_input_name, nc_vars%nc_fname, &
                     nc_vars%ncid, iout)

    if (associated(dynamic_loader)) then
      ! add dynamic pkg loader to list
      call dynamic_model%add(dynamic_loader)
    end if

    ! create subpackage list
    call static_loader%create_subpkg_list()

    ! load idm integrated subpackages
    do n = 1, static_loader%subpkg_list%pnum
      ! load subpackage
      call input_load(static_loader%subpkg_list%component_types(n), &
                      static_loader%subpkg_list%subcomponent_types(n), &
                      static_loader%mf6_input%component_name, &
                      static_loader%subpkg_list%subcomponent_types(n), &
                      static_loader%subpkg_list%pkgtypes(n), &
                      static_loader%subpkg_list%filenames(n), &
                      modelfname, nc_vars, iout)
    end do

    ! cleanup
    call static_loader%destroy()
    deallocate (static_loader)
  end subroutine input_load

  !> @brief load integrated model package files
  !<
  subroutine load_model_pkgs(model_pkg_inputs, iout)
    use ModelPackageInputsModule, only: ModelPackageInputsType
    use NCFileVarsModule, only: NCFileVarsType
    use SourceLoadModule, only: open_source_file, netcdf_context
    use IdmDfnSelectorModule, only: idm_integrated
    type(ModelPackageInputsType), intent(inout) :: model_pkg_inputs
    integer(i4B), intent(in) :: iout
    type(NCFileVarsType), pointer :: nc_vars
    integer(I4B) :: itype, ipkg

    nc_vars => netcdf_context(model_pkg_inputs%modeltype, &
                              model_pkg_inputs%component_type, &
                              model_pkg_inputs%modelname, &
                              model_pkg_inputs%modelfname, iout)
    ! load package instances by type
    do itype = 1, size(model_pkg_inputs%pkglist)
      ! load package instances
      do ipkg = 1, model_pkg_inputs%pkglist(itype)%pnum
        if (idm_integrated(model_pkg_inputs%component_type, &
                           model_pkg_inputs%pkglist(itype)%subcomponent_type)) &
          then
          ! only load if model pkg can read from input context
          call input_load(model_pkg_inputs%component_type, &
                          model_pkg_inputs%pkglist(itype)%subcomponent_type, &
                          model_pkg_inputs%modelname, &
                          model_pkg_inputs%pkglist(itype)%pkgnames(ipkg), &
                          model_pkg_inputs%pkglist(itype)%pkgtype, &
                          model_pkg_inputs%pkglist(itype)%filenames(ipkg), &
                          model_pkg_inputs%modelfname, nc_vars, iout)
        else
          ! open input file for package parser
          model_pkg_inputs%pkglist(itype)%inunits(ipkg) = &
            open_source_file(model_pkg_inputs%pkglist(itype)%pkgtype, &
                             model_pkg_inputs%pkglist(itype)%filenames(ipkg), &
                             model_pkg_inputs%modelfname, iout)
        end if
      end do
    end do

    ! cleanup
    call nc_vars%destroy()
    deallocate (nc_vars)
    nullify (nc_vars)
  end subroutine load_model_pkgs

  !> @brief load model namfiles and model package files
  !<
  subroutine load_models(iout)
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: mem_setptr
    use CharacterStringModule, only: CharacterStringType
    use DistributedSimModule, only: DistributedSimType, get_dsim
    use SimVariablesModule, only: idm_context, simfile
    use ModelPackageInputsModule, only: ModelPackageInputsType
    use SourceCommonModule, only: idm_component_type, inlen_check
    use SourceLoadModule, only: load_modelnam
    integer(I4B), intent(in) :: iout
    type(DistributedSimType), pointer :: ds
    integer(I4B), dimension(:), pointer :: model_loadmask
    character(len=LENMEMPATH) :: input_mempath
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mtypes !< model types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mfnames !< model file names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mnames !< model names
    character(len=LINELENGTH) :: mtype, mfname
    character(len=LENMODELNAME) :: mname
    type(ModelPackageInputsType), allocatable :: model_pkg_inputs
    integer(I4B) :: n

    ! get model mask
    ds => get_dsim()
    model_loadmask => ds%get_load_mask()

    ! set input memory path
    input_mempath = create_mem_path('SIM', 'NAM', idm_context)

    ! set pointers to input context model attribute arrays
    call mem_setptr(mtypes, 'MTYPE', input_mempath)
    call mem_setptr(mfnames, 'MFNAME', input_mempath)
    call mem_setptr(mnames, 'MNAME', input_mempath)

    do n = 1, size(mtypes)
      ! attributes for this model
      mtype = mtypes(n)
      mfname = mfnames(n)
      call inlen_check(mnames(n), mname, LENMODELNAME, 'MODELNAME')

      ! terminate if errors were detected
      if (count_errors() > 0) then
        call store_error_filename(simfile)
      end if

      ! load specified model inputs
      if (model_loadmask(n) > 0) then
        ! load model nam file
        call load_modelnam(mtype, mfname, mname, iout)
        ! create description of model packages
        allocate (model_pkg_inputs)
        call model_pkg_inputs%init(mtype, mfname, mname, iout)
        ! load packages
        call load_model_pkgs(model_pkg_inputs, iout)
        ! publish pkg info to input context
        call model_pkg_inputs%memload()
        ! cleanup
        call model_pkg_inputs%destroy()
        deallocate (model_pkg_inputs)
      end if
    end do
  end subroutine load_models

  !> @brief load exchange files
  !<
  subroutine load_exchanges(iout)
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: mem_setptr, mem_allocate, &
                                   mem_deallocate, get_isize
    use CharacterStringModule, only: CharacterStringType
    use SimVariablesModule, only: idm_context, simfile
    use DistributedSimModule, only: DistributedSimType, get_dsim
    use SourceCommonModule, only: idm_subcomponent_type, ifind_charstr, &
                                  inlen_check
    use SourceLoadModule, only: create_input_loader, remote_model_ndim
    integer(I4B), intent(in) :: iout
    type(DistributedSimType), pointer :: ds
    integer(I4B), dimension(:), pointer :: model_loadmask
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: etypes !< exg types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: efiles !< exg file names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: emnames_a !< model a names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: emnames_b !< model b names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: emempaths !< exg mempaths
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mtypes !< model types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mfnames !< model file names
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mnames !< model names
    character(len=LENMEMPATH) :: input_mempath, mempath
    integer(I4B), pointer :: exgid, ncelldim
    character(len=LINELENGTH) :: exgtype, efname, mfname
    character(len=LENMODELNAME) :: mname1, mname2, mname
    character(len=LENCOMPONENTNAME) :: sc_type, sc_name, mtype
    class(StaticPkgLoadBaseType), pointer :: static_loader
    class(DynamicPkgLoadBaseType), pointer :: dynamic_loader
    integer(I4B) :: n, m1_idx, m2_idx, irem, isize

    ! get model mask
    ds => get_dsim()
    model_loadmask => ds%get_load_mask()

    ! set input memory path
    input_mempath = create_mem_path('SIM', 'NAM', idm_context)

    ! set pointers to input context exg and model attribute arrays
    call mem_setptr(etypes, 'EXGTYPE', input_mempath)
    call mem_setptr(efiles, 'EXGFILE', input_mempath)
    call mem_setptr(emnames_a, 'EXGMNAMEA', input_mempath)
    call mem_setptr(emnames_b, 'EXGMNAMEB', input_mempath)
    call mem_setptr(mtypes, 'MTYPE', input_mempath)
    call mem_setptr(mfnames, 'MFNAME', input_mempath)
    call mem_setptr(mnames, 'MNAME', input_mempath)

    ! allocate mempaths array for exchanges
    call mem_allocate(emempaths, LENMEMPATH, size(etypes), 'EXGMEMPATHS', &
                      input_mempath)

    ! load exchanges for local models
    do n = 1, size(etypes)
      ! attributes for this exchange
      exgtype = etypes(n)
      efname = efiles(n)
      call inlen_check(emnames_a(n), mname1, LENMODELNAME, 'MODELNAME')
      call inlen_check(emnames_b(n), mname2, LENMODELNAME, 'MODELNAME')

      ! initialize mempath as no path
      emempaths(n) = ''
      irem = 0

      ! set indexes for exchange model names
      m1_idx = ifind_charstr(mnames, mname1)
      m2_idx = ifind_charstr(mnames, mname2)

      if (m1_idx <= 0 .or. m2_idx <= 0) then
        errmsg = 'Exchange has invalid (unrecognized) model name(s):'
        if (m1_idx <= 0) errmsg = trim(errmsg)//' '//trim(mname1)
        if (m2_idx <= 0) errmsg = trim(errmsg)//' '//trim(mname2)
        call store_error(errmsg)
      end if

      ! terminate if errors were detected
      if (count_errors() > 0) then
        call store_error_filename(simfile)
      end if

      ! load the exchange input if either model local
      if (model_loadmask(m1_idx) > 0 .or. model_loadmask(m2_idx) > 0) then
        ! set index if either model is remote
        if (model_loadmask(m1_idx) == 0) then
          irem = m1_idx
        else if (model_loadmask(m2_idx) == 0) then
          irem = m2_idx
        end if

        ! allocate and set remote model NCELLDIM
        if (irem > 0) then
          mtype = mtypes(irem)
          mfname = mfnames(irem)
          mname = mnames(irem)
          mempath = create_mem_path(component=mname, context=idm_context)
          call get_isize('NCELLDIM', mempath, isize)
          if (isize < 0) then
            call mem_allocate(ncelldim, 'NCELLDIM', mempath)
            ncelldim = remote_model_ndim(mtype, mfname)
          else
            call mem_setptr(ncelldim, 'NCELLDIM', mempath)
          end if
        else
          nullify (ncelldim)
        end if

        ! set subcomponent strings
        sc_type = trim(idm_subcomponent_type('EXG', exgtype))
        write (sc_name, '(a,i0)') trim(sc_type)//'_', n

        ! create and set exchange mempath
        mempath = create_mem_path('EXG', sc_name, idm_context)
        emempaths(n) = mempath

        ! allocate and set exgid
        call mem_allocate(exgid, 'EXGID', mempath)
        exgid = n

        ! create exchange loader
        static_loader => create_input_loader('EXG', sc_type, 'EXG', sc_name, &
                                             exgtype, 'SIM', efname, simfile)
        ! load static input
        dynamic_loader => static_loader%load(iout)

        if (associated(dynamic_loader)) then
          errmsg = 'IDM unimplemented. Dynamic Exchanges not supported.'
          call store_error(errmsg)
          call store_error_filename(efname)
        else
          call static_loader%destroy()
          deallocate (static_loader)
        end if
      end if
    end do

    ! clean up temporary NCELLDIM for remote models
    do n = 1, size(mnames)
      if (model_loadmask(n) == 0) then
        mname = mnames(n)
        mempath = create_mem_path(component=mname, context=idm_context)
        call get_isize('NCELLDIM', mempath, isize)
        if (isize > 0) then
          call mem_setptr(ncelldim, 'NCELLDIM', mempath)
          call mem_deallocate(ncelldim)
        end if
      end if
    end do
  end subroutine load_exchanges

  !> @brief MODFLOW 6 mfsim.nam input load routine
  !<
  subroutine simnam_load(paramlog)
    use SourceLoadModule, only: load_simnam
    integer(I4B), intent(inout) :: paramlog
    ! load sim nam file
    call load_simnam()
    ! allocate any unallocated simnam params
    call simnam_allocate()
    ! read and set input parameter logging keyword
    paramlog = input_param_log()
    ! memload summary info
    call simnam_load_dim()
  end subroutine simnam_load

  !> @brief MODFLOW 6 tdis input load routine
  !<
  subroutine simtdis_load()
    use SourceLoadModule, only: load_simtdis
    ! load sim tdis file
    call load_simtdis()
  end subroutine simtdis_load

  !> @brief retrieve list of model dynamic loaders
  !<
  function dynamic_models(modeltype, modelname, modelfname, nc_fname, &
                          ncid, iout) result(model_dynamic_input)
    use InputLoadTypeModule, only: AddDynamicModelToList, GetDynamicModelFromList
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: modelfname
    character(len=*), intent(in) :: nc_fname
    integer(I4B), intent(in) :: ncid
    integer(I4B), intent(in) :: iout
    class(ModelDynamicPkgsType), pointer :: model_dynamic_input
    class(ModelDynamicPkgsType), pointer :: temp
    integer(I4B) :: id

    ! initialize
    nullify (model_dynamic_input)

    ! assign model loader object if found
    do id = 1, model_inputs%Count()
      temp => GetDynamicModelFromList(model_inputs, id)
      if (temp%modelname == modelname) then
        model_dynamic_input => temp
        exit
      end if
    end do

    ! create if not found
    if (.not. associated(model_dynamic_input)) then
      allocate (model_dynamic_input)
      call model_dynamic_input%init(modeltype, modelname, modelfname, &
                                    nc_fname, ncid, iout)
      call AddDynamicModelToList(model_inputs, model_dynamic_input)
    end if
  end function dynamic_models

  !> @brief deallocate all model dynamic loader collections
  !<
  subroutine dynamic_da(iout)
    use InputLoadTypeModule, only: GetDynamicModelFromList
    use SourceLoadModule, only: nc_close
    integer(I4B), intent(in) :: iout
    class(ModelDynamicPkgsType), pointer :: model_dynamic_input
    integer(I4B) :: n
    do n = 1, model_inputs%Count()
      model_dynamic_input => GetDynamicModelFromList(model_inputs, n)
      call nc_close(model_dynamic_input%ncid, model_dynamic_input%nc_fname)
      call model_dynamic_input%destroy()
      deallocate (model_dynamic_input)
      nullify (model_dynamic_input)
    end do
    call model_inputs%Clear()
  end subroutine dynamic_da

  !> @brief return sim input context PRINT_INPUT value
  !<
  function input_param_log() result(paramlog)
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: mem_setptr
    use SimVariablesModule, only: idm_context
    character(len=LENMEMPATH) :: simnam_mempath
    integer(I4B) :: paramlog
    integer(I4B), pointer :: p
    ! read and set input value of PRINT_INPUT
    simnam_mempath = create_mem_path('SIM', 'NAM', idm_context)
    call mem_setptr(p, 'PRINT_INPUT', simnam_mempath)
    paramlog = p
  end function input_param_log

  !> @brief load simulation summary info to input context
  !<
  subroutine simnam_load_dim()
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    use SimVariablesModule, only: idm_context
    use CharacterStringModule, only: CharacterStringType
    character(len=LENMEMPATH) :: sim_mempath, simnam_mempath
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mtypes !< model types
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: etypes !< model types
    integer(I4B), pointer :: nummodels
    integer(I4B), pointer :: numexchanges

    ! initialize
    nullify (nummodels)
    nullify (numexchanges)

    ! set memory paths
    sim_mempath = create_mem_path(component='SIM', context=idm_context)
    simnam_mempath = create_mem_path('SIM', 'NAM', idm_context)

    ! set pointers to loaded simnam arrays
    call mem_setptr(mtypes, 'MTYPE', simnam_mempath)
    call mem_setptr(etypes, 'EXGTYPE', simnam_mempath)

    ! allocate variables
    call mem_allocate(nummodels, 'NUMMODELS', sim_mempath)
    call mem_allocate(numexchanges, 'NUMEXCHANGES', sim_mempath)

    ! set values
    nummodels = size(mtypes)
    numexchanges = size(etypes)
  end subroutine simnam_load_dim

  !> @brief set sim nam input context default integer value
  !<
  subroutine allocate_simnam_int(input_mempath, idt)
    use MemoryManagerModule, only: mem_allocate
    use SimVariablesModule, only: isimcontinue, isimcheck, simfile
    character(len=LENMEMPATH), intent(in) :: input_mempath
    type(InputParamDefinitionType), pointer, intent(in) :: idt
    integer(I4B), pointer :: intvar

    ! allocate and set default
    call mem_allocate(intvar, idt%mf6varname, input_mempath)

    select case (idt%mf6varname)
    case ('CONTINUE')
      intvar = isimcontinue
    case ('NOCHECK')
      intvar = isimcheck
    case ('MAXERRORS')
      intvar = 1000 !< MessageType max_message
    case ('MXITER')
      intvar = 1
    case ('PRINT_INPUT')
      intvar = 0
    case default
      write (errmsg, '(a,a)') &
        'Idm SIMNAM Load default value setting '&
        &'is unhandled for this variable: ', &
        trim(idt%mf6varname)
      call store_error(errmsg)
      call store_error_filename(simfile)
    end select
  end subroutine allocate_simnam_int

  !> @brief MODFLOW 6 mfsim.nam parameter allocate and set
  !<
  subroutine allocate_simnam_param(input_mempath, idt)
    use SimVariablesModule, only: simfile
    use MemoryManagerModule, only: mem_allocate
    use CharacterStringModule, only: CharacterStringType
    use DefinitionSelectModule, only: idt_datatype
    character(len=LENMEMPATH), intent(in) :: input_mempath
    type(InputParamDefinitionType), pointer, intent(in) :: idt
    character(len=LINELENGTH), pointer :: cstr
    type(CharacterStringType), dimension(:), &
      pointer, contiguous :: acharstr1d

    select case (idt_datatype(idt))
    case ('KEYWORD', 'INTEGER')
      if (idt%in_record) then
        ! no-op
      else
        ! allocate and set default
        call allocate_simnam_int(input_mempath, idt)
      end if
    case ('STRING')
      ! did this param originate from sim namfile RECARRAY type
      if (idt%in_record) then
        ! allocate 0 size CharacterStringType array
        call mem_allocate(acharstr1d, LINELENGTH, 0, idt%mf6varname, &
                          input_mempath)
      else
        ! allocate empty string
        call mem_allocate(cstr, LINELENGTH, idt%mf6varname, input_mempath)
        cstr = ''
      end if
    case ('RECORD')
      ! no-op
    case default
      write (errmsg, '(a,a)') &
        'IdmLoad allocate simnam param unhandled datatype: ', &
        trim(idt%datatype)
      call store_error(errmsg)
      call store_error_filename(simfile)
    end select
  end subroutine allocate_simnam_param

  !> @brief MODFLOW 6 mfsim.nam input context parameter allocation
  !<
  subroutine simnam_allocate()
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: get_isize, mem_allocate
    use SimVariablesModule, only: idm_context
    use CharacterStringModule, only: CharacterStringType
    character(len=LENMEMPATH) :: input_mempath
    type(ModflowInputType) :: mf6_input
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: iparam, isize

    ! set memory path
    input_mempath = create_mem_path('SIM', 'NAM', idm_context)
    ! create description of input
    mf6_input = getModflowInput('NAM6', 'SIM', 'NAM', 'SIM', 'NAM', 'ROOT')

    ! allocate sim namfile parameters if not in input context
    do iparam = 1, size(mf6_input%param_dfns)
      ! assign param definition pointer
      idt => mf6_input%param_dfns(iparam)
      ! check if variable is already allocated
      call get_isize(idt%mf6varname, input_mempath, isize)
      if (isize < 0) then
        ! allocate and set parameter
        call allocate_simnam_param(input_mempath, idt)
      end if
    end do
  end subroutine simnam_allocate

end module IdmLoadModule
