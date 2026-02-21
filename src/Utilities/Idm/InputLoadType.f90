!> @brief This module contains the InputLoadTypeModule
!!
!! This module defines types that support generic IDM
!! static and dynamic input loading.
!!
!<
module InputLoadTypeModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENCOMPONENTNAME, LENMODELNAME, &
                             LENMEMPATH, LENVARNAME, LENFTYPE
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use ModflowInputModule, only: ModflowInputType
  use ListModule, only: ListType
  use InputDefinitionModule, only: InputParamDefinitionType
  use NCFileVarsModule, only: NCPackageVarsType

  implicit none
  private
  public :: StaticPkgLoadBaseType
  public :: DynamicPkgLoadBaseType
  public :: ModelDynamicPkgsType
  public :: AddDynamicModelToList, GetDynamicModelFromList
  public :: StaticPkgLoadType, DynamicPkgLoadType
  public :: SubPackageListType
  public :: model_inputs

  !> @brief type representing package subpackage list
  type :: SubPackageListType
    character(len=LENCOMPONENTNAME), dimension(:), allocatable :: pkgtypes
    character(len=LENCOMPONENTNAME), dimension(:), allocatable :: component_types
    character(len=LENCOMPONENTNAME), dimension(:), &
      allocatable :: subcomponent_types
    character(len=LENCOMPONENTNAME), dimension(:), &
      allocatable :: subcomponent_names
    character(len=LINELENGTH), dimension(:), allocatable :: filenames
    character(len=LENCOMPONENTNAME) :: component_type
    character(len=LENCOMPONENTNAME) :: component_name
    integer(I4B) :: pnum
  contains
    procedure :: create => subpkg_create
    procedure :: add => subpkg_add
    procedure :: set_names => subpkg_names
    procedure :: destroy => subpkg_destroy
  end type SubPackageListType

  !> @brief Static loader type
  !!
  !! This type is a base concrete type for a static input loader
  !!
  !<
  type StaticPkgLoadType
    type(ModflowInputType) :: mf6_input !< description of modflow6 input
    type(NCPackageVarsType), pointer :: nc_vars => null()
    character(len=LENCOMPONENTNAME) :: component_name !< name of component
    character(len=LINELENGTH) :: component_input_name !< component input name, e.g. model name file
    character(len=LINELENGTH) :: input_name !< input name, e.g. package *.chd file
    integer(I4B) :: iperblock !< index of period block on block definition list
    type(SubPackageListType) :: subpkg_list !< list of input subpackages
  contains
    procedure :: init => static_init
    procedure :: create_subpkg_list
    procedure :: destroy => static_destroy
  end type StaticPkgLoadType

  !> @brief Base abstract type for static input loader
  !!
  !! IDM sources should extend and implement this type
  !!
  !<
  type, abstract, extends(StaticPkgLoadType) :: StaticPkgLoadBaseType
  contains
    procedure(load_if), deferred :: load
  end type StaticPkgLoadBaseType

  !> @brief Dynamic loader type
  !!
  !! This type is a base concrete type for a dynamic (period) input loader
  !!
  !<
  type :: DynamicPkgLoadType
    type(ModflowInputType) :: mf6_input !< description of modflow6 input
    type(NCPackageVarsType), pointer :: nc_vars => null()
    character(len=LENCOMPONENTNAME) :: component_name !< name of component
    character(len=LINELENGTH) :: component_input_name !< component input name, e.g. model name file
    character(len=LINELENGTH) :: input_name !< input name, e.g. package *.chd file
    character(len=LINELENGTH), dimension(:), allocatable :: param_names !< dynamic param tagnames
    logical(LGP) :: readasarrays !< readasarrays style input package
    logical(LGP) :: readarraygrid !< readarraygrid style input package
    logical(LGP) :: has_setting !< period block contains setting keystring param
    integer(I4B) :: iperblock !< index of period block on block definition list
    integer(I4B) :: iout !< inunit number for logging
    integer(I4B) :: nparam !< number of in scope params
  contains
    procedure :: init => dynamic_init
    procedure :: df => dynamic_df
    procedure :: ad => dynamic_ad
    procedure :: destroy => dynamic_destroy
  end type DynamicPkgLoadType

  !> @brief Base abstract type for dynamic input loader
  !!
  !! IDM sources should extend and implement this type
  !!
  !<
  type, abstract, extends(DynamicPkgLoadType) :: DynamicPkgLoadBaseType
  contains
    procedure(period_load_if), deferred :: rp
  end type DynamicPkgLoadBaseType

  !> @brief load interfaces for source static and dynamic types
  !<
  abstract interface
    function load_if(this, iout) result(dynamic_loader)
      import StaticPkgLoadBaseType, DynamicPkgLoadBaseType, I4B
      class(StaticPkgLoadBaseType), intent(inout) :: this
      integer(I4B), intent(in) :: iout
      class(DynamicPkgLoadBaseType), pointer :: dynamic_loader
    end function load_if
    subroutine period_load_if(this)
      import DynamicPkgLoadBaseType, I4B
      class(DynamicPkgLoadBaseType), intent(inout) :: this
    end subroutine
  end interface

  !> @brief type for storing a dynamic package load list
  !!
  !! This type is used to store a list of package
  !! dynamic load types for a model
  !!
  !<
  type :: ModelDynamicPkgsType
    character(len=LENCOMPONENTNAME) :: modeltype !< type of model
    character(len=LENMODELNAME) :: modelname !< name of model
    character(len=LINELENGTH) :: modelfname !< name of model input file
    type(ListType) :: pkglist !< model package list
    character(len=LINELENGTH) :: nc_fname !< name of model netcdf input
    integer(I4B) :: ncid !< netcdf file handle
    integer(I4B) :: iout
  contains
    procedure :: init => dynamicpkgs_init
    procedure :: add => dynamicpkgs_add
    procedure :: get => dynamicpkgs_get
    procedure :: rp => dynamicpkgs_rp
    procedure :: df => dynamicpkgs_df
    procedure :: ad => dynamicpkgs_ad
    procedure :: size => dynamicpkgs_size
    procedure :: destroy => dynamicpkgs_destroy
  end type ModelDynamicPkgsType

  type(ListType) :: model_inputs

contains

  !> @brief create a new package type
  !<
  subroutine subpkg_create(this, component_type, component_name)
    class(SubPackageListType) :: this
    character(len=*), intent(in) :: component_type
    character(len=*), intent(in) :: component_name

    ! initialize
    this%pnum = 0
    this%component_type = component_type
    this%component_name = component_name

    ! allocate arrays
    allocate (this%pkgtypes(0))
    allocate (this%component_types(0))
    allocate (this%subcomponent_types(0))
    allocate (this%subcomponent_names(0))
    allocate (this%filenames(0))
  end subroutine subpkg_create

  !> @brief create a new package type
  !<
  subroutine subpkg_add(this, pkgtype, component_type, subcomponent_type, &
                        filename)
    use ArrayHandlersModule, only: expandarray
    class(SubPackageListType) :: this
    character(len=*), intent(in) :: pkgtype
    character(len=*), intent(in) :: component_type
    character(len=*), intent(in) :: subcomponent_type
    character(len=*), intent(in) :: filename

    ! reallocate
    call expandarray(this%pkgtypes)
    call expandarray(this%component_types)
    call expandarray(this%subcomponent_types)
    call expandarray(this%subcomponent_names)
    call expandarray(this%filenames)

    ! add new package instance
    this%pnum = this%pnum + 1
    this%pkgtypes(this%pnum) = pkgtype
    this%component_types(this%pnum) = component_type
    this%subcomponent_types(this%pnum) = subcomponent_type
    this%subcomponent_names(this%pnum) = ''
    this%filenames(this%pnum) = filename
  end subroutine subpkg_add

  !> @brief set subpackage mempaths
  !<
  subroutine subpkg_names(this, sc_type, sc_name, sc_mempath, modelfname)
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    use ArrayHandlersModule, only: expandarray
    use SimVariablesModule, only: idm_context
    use CharacterStringModule, only: CharacterStringType
    use ModelPackageInputModule, only: multi_package_type
    use IdmDfnSelectorModule, only: idm_multi_package, idm_integrated
    use SourceCommonModule, only: idm_utl_type, idm_component_type
    class(SubPackageListType) :: this
    character(len=*), intent(in) :: sc_type !< parent subcomponent type
    character(len=*), intent(in) :: sc_name !< parent subcomponent name
    character(len=*), intent(in) :: sc_mempath !< parent mempath
    character(len=*), intent(in) :: modelfname !< model name file
    character(len=LINELENGTH), dimension(:), allocatable :: ptypes
    integer(I4B), dimension(:), allocatable :: nptypes
    type(CharacterStringType), dimension(:), &
      pointer, contiguous :: mempaths
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pnames, ftypes
    character(len=LINELENGTH), pointer :: input_fname
    character(len=LENVARNAME) :: mempath_tag
    character(len=LENVARNAME) :: pname_prefix, pname, ptype, sctype
    character(len=LENMEMPATH) :: mempath, model_mempath
    integer(I4B) :: n, m, idx
    logical(LGP) :: multi

    ! no subpackages to load
    if (size(this%pkgtypes) == 0) return

    ! assume utility types do not support subpackages
    if (idm_utl_type(this%component_type, sc_type)) return

    ! initialize
    pname_prefix = ''
    allocate (ptypes(0))
    allocate (nptypes(0))

    ! reset updated readasarrays, readarraygrid types
    select case (sc_type)
    case ('EVTA', 'RCHA', &
          'RIVG', 'CHDG', &
          'WELG', 'DRNG', &
          'GHBG')
      sctype = sc_type(1:3)
    case default
      sctype = sc_type
    end select

    ! set package type
    ptype = trim(sctype)//'6'

    ! determine if multi instance package
    if (idm_integrated(this%component_type, sctype)) then
      multi = idm_multi_package(this%component_type, sctype)
    else
      multi = multi_package_type(this%component_type, sctype, ptype)
    end if

    ! set package name prefix based on parent package
    if (this%component_type == 'EXG') then
      ! no-op
    else
      if (multi) then
        idx = 0
        model_mempath = create_mem_path(this%component_name, 'NAM', idm_context)
        call mem_setptr(pnames, 'PNAME', model_mempath)
        call mem_setptr(ftypes, 'FTYPE', model_mempath)

        ! identify instance number of the package type to distinguish in name
        do n = 1, size(pnames)
          if (ftypes(n) == ptype) then
            idx = idx + 1
            ! multi packages must match on type and name
            pname = pnames(n)
            if (pname == '') then
              write (pname, '(a,i0)') trim(sctype)//'-', idx
            end if
            if (pname == sc_name) then
              write (pname_prefix, '(a,i0,a)') trim(sctype), idx, '-'
              exit
            end if
          end if
        end do

        if (pname_prefix == '') then
          errmsg = &
            'Internal IDM error: subpackage load cannot identify &
            &package "'//trim(sc_name)//'" in model name file &
            &packages block.'
          call store_error(errmsg)
          call store_error_filename(modelfname)
        end if

      else
        write (pname_prefix, '(2a)') trim(sctype), '-'
      end if
    end if

    ! count number of each package type
    do n = 1, size(this%pkgtypes)
      idx = 0
      do m = 1, size(ptypes)
        if (ptypes(m) == this%pkgtypes(n)) then
          nptypes(m) = nptypes(m) + 1
          idx = m
        end if
      end do
      if (idx == 0) then
        call expandarray(ptypes)
        call expandarray(nptypes)
        ptypes(size(ptypes)) = this%pkgtypes(n)
        nptypes(size(ptypes)) = 1
      end if
    end do

    ! set subpackage names and mempaths
    do n = 1, size(ptypes)
      idx = 0
      mempath_tag = trim(ptypes(n))//'_MEMPATH'
      call mem_allocate(mempaths, LENMEMPATH, nptypes(n), &
                        mempath_tag, sc_mempath)
      do m = 1, size(this%pkgtypes)
        if (this%pkgtypes(m) == ptypes(n)) then
          idx = idx + 1
          write (this%subcomponent_names(m), '(a,i0)') &
            trim(pname_prefix)//trim(this%subcomponent_types(m)), idx
          mempath = create_mem_path(this%component_name, &
                                    this%subcomponent_names(m), &
                                    idm_context)
          mempaths(idx) = mempath
          call mem_allocate(input_fname, LINELENGTH, 'INPUT_FNAME', mempath)
          input_fname = trim(this%filenames(m))
        end if
      end do
    end do

    ! cleanup
    deallocate (ptypes)
    deallocate (nptypes)
  end subroutine subpkg_names

  !> @brief create a new package type
  !<
  subroutine subpkg_destroy(this)
    class(SubPackageListType) :: this
    ! allocate arrays
    deallocate (this%pkgtypes)
    deallocate (this%component_types)
    deallocate (this%subcomponent_types)
    deallocate (this%subcomponent_names)
    deallocate (this%filenames)
  end subroutine subpkg_destroy

  !> @brief initialize static package loader
  !!
  !<
  subroutine static_init(this, mf6_input, component_name, component_input_name, &
                         input_name)
    class(StaticPkgLoadType), intent(inout) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), intent(in) :: component_name
    character(len=*), intent(in) :: component_input_name
    character(len=*), intent(in) :: input_name
    integer(I4B) :: iblock

    this%mf6_input = mf6_input
    this%component_name = component_name
    this%component_input_name = component_input_name
    this%input_name = input_name
    this%iperblock = 0

    ! create subpackage list
    call this%subpkg_list%create(this%mf6_input%component_type, &
                                 this%mf6_input%component_name)

    ! identify period block definition
    do iblock = 1, size(mf6_input%block_dfns)
      if (mf6_input%block_dfns(iblock)%blockname == 'PERIOD') then
        this%iperblock = iblock
        exit
      end if
    end do
  end subroutine static_init

  !> @brief create the subpackage list
  !!
  !<
  subroutine create_subpkg_list(this)
    use IdmDfnSelectorModule, only: idm_subpackages, idm_integrated, &
                                    idm_multi_package
    use SourceCommonModule, only: filein_fname
    use MemoryManagerModule, only: mem_setptr, get_isize
    use ArrayHandlersModule, only: expandarray
    use CharacterStringModule, only: CharacterStringType
    class(StaticPkgLoadType), intent(inout) :: this
    character(len=16), dimension(:), pointer :: subpkgs
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: fnames
    character(len=LINELENGTH) :: tag, fname, pkgtype
    character(len=LENFTYPE) :: c_type, sc_type
    character(len=16) :: subpkg
    integer(I4B) :: idx, n, m, isize

    ! set pointer to package (idm integrated) subpackage list
    subpkgs => idm_subpackages(this%mf6_input%component_type, &
                               this%mf6_input%subcomponent_type)

    do n = 1, size(subpkgs)
      ! check for input matching this supported subpackage
      subpkg = subpkgs(n)
      idx = index(subpkg, '-')

      if (idx > 0) then
        ! split string in component/subcomponent types
        c_type = subpkg(1:idx - 1)
        sc_type = subpkg(idx + 1:len_trim(subpkg))

        if (idm_integrated(c_type, sc_type)) then
          ! construct FILEIN filename tag
          pkgtype = trim(sc_type)//'6'
          tag = trim(pkgtype)//'_FILENAME'
          call get_isize(tag, this%mf6_input%mempath, isize)
          if (isize > 0) then

            call mem_setptr(fnames, tag, this%mf6_input%mempath)
            do m = 1, size(fnames)
              fname = fnames(m)
              call this%subpkg_list%add(pkgtype, c_type, sc_type, fname)
            end do
          end if
        else
          errmsg = 'Identified subpackage is not IDM integrated. Remove dfn &
                   &subpackage tagline for package "'//trim(subpkg)//'".'
          call store_error(errmsg)
          call store_error_filename(this%input_name)
        end if
      end if
    end do

    call this%subpkg_list%set_names(this%mf6_input%subcomponent_type, &
                                    this%mf6_input%subcomponent_name, &
                                    this%mf6_input%mempath, &
                                    this%component_input_name)
  end subroutine create_subpkg_list

  subroutine static_destroy(this)
    class(StaticPkgLoadType), intent(inout) :: this
    call this%subpkg_list%destroy()
    if (associated(this%nc_vars)) then
      call this%nc_vars%destroy()
      deallocate (this%nc_vars)
      nullify (this%nc_vars)
    end if
  end subroutine static_destroy

  !> @brief initialize dynamic package loader
  !!
  !! Any managed memory pointed to from model/package context
  !! must be allocated when dynamic loader is initialized.
  !!
  !<
  subroutine dynamic_init(this, mf6_input, component_name, component_input_name, &
                          input_name, iperblock, iout)
    use SimVariablesModule, only: errmsg
    use InputDefinitionModule, only: InputParamDefinitionType
    use DefinitionSelectModule, only: idt_datatype
    class(DynamicPkgLoadType), intent(inout) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), intent(in) :: component_name
    character(len=*), intent(in) :: component_input_name
    character(len=*), intent(in) :: input_name
    integer(I4B), intent(in) :: iperblock
    integer(I4B), intent(in) :: iout
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: iparam, ilen

    this%mf6_input = mf6_input
    this%component_name = component_name
    this%component_input_name = component_input_name
    this%input_name = input_name
    this%readasarrays = .false.
    this%readarraygrid = .false.
    this%has_setting = .false.
    this%iperblock = iperblock
    this%nparam = 0
    this%iout = iout
    nullify (idt)

    ! throw error and exit if not found
    if (this%iperblock == 0) then
      write (errmsg, '(a,a)') &
        'Programming error. (IDM) PERIOD block not found in '&
        &'dynamic package input block dfns: ', &
        trim(mf6_input%subcomponent_name)
      call store_error(errmsg)
      call store_error_filename(this%input_name)
    end if

    ! set readasarrays and readarraygrid
    if (mf6_input%block_dfns(iperblock)%aggregate) then
      ! no-op, list based input
    else
      do iparam = 1, size(mf6_input%param_dfns)
        idt => mf6_input%param_dfns(iparam)
        if (idt%blockname == 'OPTIONS') then
          select case (idt%tagname)
          case ('READASARRAYS')
            this%readasarrays = .true.
          case ('READARRAYGRID')
            this%readarraygrid = .true.
          case default
            ! no-op
          end select
        end if
      end do
    end if

    ! determine if has setting type
    do iparam = 1, size(mf6_input%param_dfns)
      idt => mf6_input%param_dfns(iparam)
      if (idt%blockname == 'PERIOD') then
        if (idt_datatype(idt) == 'KEYSTRING') then
          ilen = len_trim(idt%tagname)
          if (idt%tagname(ilen - 6:ilen) == 'SETTING') then
            this%has_setting = .true.
          end if
        end if
      end if
    end do
  end subroutine dynamic_init

  !> @brief dynamic package loader define
  !!
  !<
  subroutine dynamic_df(this)
    class(DynamicPkgLoadType), intent(inout) :: this
    ! override in derived type
  end subroutine dynamic_df

  !> @brief dynamic package loader advance
  !!
  !<
  subroutine dynamic_ad(this)
    class(DynamicPkgLoadType), intent(inout) :: this
    ! override in derived type
  end subroutine dynamic_ad

  !> @brief dynamic package loader destroy
  !!
  !<
  subroutine dynamic_destroy(this)
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    class(DynamicPkgLoadType), intent(inout) :: this

    ! clean up netcdf variables structure
    if (associated(this%nc_vars)) then
      call this%nc_vars%destroy()
      deallocate (this%nc_vars)
      nullify (this%nc_vars)
    end if

    ! deallocate package static and dynamic input context
    call memorystore_remove(this%mf6_input%component_name, &
                            this%mf6_input%subcomponent_name, &
                            idm_context)
  end subroutine dynamic_destroy

  !> @brief model dynamic packages init
  !!
  !<
  subroutine dynamicpkgs_init(this, modeltype, modelname, modelfname, nc_fname, &
                              ncid, iout)
    class(ModelDynamicPkgsType), intent(inout) :: this
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: modelfname
    character(len=*), intent(in) :: nc_fname
    integer(I4B), intent(in) :: ncid
    integer(I4B), intent(in) :: iout
    this%modeltype = modeltype
    this%modelname = modelname
    this%modelfname = modelfname
    this%nc_fname = nc_fname
    this%ncid = ncid
    this%iout = iout
  end subroutine dynamicpkgs_init

  !> @brief add package to model dynamic packages list
  !!
  !<
  subroutine dynamicpkgs_add(this, dynamic_pkg)
    class(ModelDynamicPkgsType), intent(inout) :: this
    class(DynamicPkgLoadBaseType), pointer, intent(inout) :: dynamic_pkg
    class(*), pointer :: obj
    obj => dynamic_pkg
    call this%pkglist%add(obj)
  end subroutine dynamicpkgs_add

  !> @brief retrieve package from model dynamic packages list
  !!
  !<
  function dynamicpkgs_get(this, idx) result(res)
    class(ModelDynamicPkgsType), intent(inout) :: this
    integer(I4B), intent(in) :: idx
    class(DynamicPkgLoadBaseType), pointer :: res
    class(*), pointer :: obj
    nullify (res)
    obj => this%pkglist%GetItem(idx)
    if (associated(obj)) then
      select type (obj)
      class is (DynamicPkgLoadBaseType)
        res => obj
      end select
    end if
  end function dynamicpkgs_get

  !> @brief read and prepare model dynamic packages
  !!
  !<
  subroutine dynamicpkgs_rp(this)
    use IdmLoggerModule, only: idm_log_period_header, idm_log_period_close
    class(ModelDynamicPkgsType), intent(inout) :: this
    class(DynamicPkgLoadBaseType), pointer :: dynamic_pkg
    integer(I4B) :: n
    call idm_log_period_header(this%modelname, this%iout)
    do n = 1, this%pkglist%Count()
      dynamic_pkg => this%get(n)
      call dynamic_pkg%rp()
    end do
    call idm_log_period_close(this%iout)
  end subroutine dynamicpkgs_rp

  !> @brief define model dynamic packages
  !!
  !<
  subroutine dynamicpkgs_df(this)
    class(ModelDynamicPkgsType), intent(inout) :: this
    class(DynamicPkgLoadBaseType), pointer :: dynamic_pkg
    integer(I4B) :: n
    do n = 1, this%pkglist%Count()
      dynamic_pkg => this%get(n)
      call dynamic_pkg%df()
    end do
  end subroutine dynamicpkgs_df

  !> @brief advance model dynamic packages
  !!
  !<
  subroutine dynamicpkgs_ad(this)
    class(ModelDynamicPkgsType), intent(inout) :: this
    class(DynamicPkgLoadBaseType), pointer :: dynamic_pkg
    integer(I4B) :: n
    do n = 1, this%pkglist%Count()
      dynamic_pkg => this%get(n)
      call dynamic_pkg%ad()
    end do
  end subroutine dynamicpkgs_ad

  !> @brief get size of model dynamic packages list
  !!
  !<
  function dynamicpkgs_size(this) result(size)
    class(ModelDynamicPkgsType), intent(inout) :: this
    integer(I4B) :: size
    size = this%pkglist%Count()
  end function dynamicpkgs_size

  !> @brief destroy model dynamic packages object
  !!
  !<
  subroutine dynamicpkgs_destroy(this)
    class(ModelDynamicPkgsType), intent(inout) :: this
    class(DynamicPkgLoadBaseType), pointer :: dynamic_pkg
    integer(I4B) :: n
    ! destroy dynamic loaders
    do n = 1, this%pkglist%Count()
      dynamic_pkg => this%get(n)
      call dynamic_pkg%destroy()
      deallocate (dynamic_pkg)
      nullify (dynamic_pkg)
    end do
    call this%pkglist%Clear()
  end subroutine dynamicpkgs_destroy

  !> @brief add model dynamic packages object to list
  !!
  !<
  subroutine AddDynamicModelToList(list, model_dynamic)
    type(ListType), intent(inout) :: list !< package list
    class(ModelDynamicPkgsType), pointer, intent(inout) :: model_dynamic
    class(*), pointer :: obj
    obj => model_dynamic
    call list%Add(obj)
  end subroutine AddDynamicModelToList

  !> @brief get model dynamic packages object from list
  !!
  !<
  function GetDynamicModelFromList(list, idx) result(res)
    type(ListType), intent(inout) :: list !< spd list
    integer(I4B), intent(in) :: idx !< package number
    class(ModelDynamicPkgsType), pointer :: res
    class(*), pointer :: obj
    ! initialize res
    nullify (res)
    ! get the object from the list
    obj => list%GetItem(idx)
    if (associated(obj)) then
      select type (obj)
      class is (ModelDynamicPkgsType)
        res => obj
      end select
    end if
  end function GetDynamicModelFromList

end module InputLoadTypeModule
