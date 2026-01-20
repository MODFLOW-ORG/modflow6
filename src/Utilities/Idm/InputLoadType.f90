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
  public :: model_inputs

  !> @brief type representing package subpackage list
  type :: SubPackageListType
    character(len=LENCOMPONENTNAME), dimension(:), allocatable :: pkgtypes
    character(len=LENCOMPONENTNAME), dimension(:), allocatable :: component_types
    character(len=LENCOMPONENTNAME), dimension(:), &
      allocatable :: subcomponent_types
    character(len=LINELENGTH), dimension(:), allocatable :: filenames
    character(len=LENMEMPATH) :: mempath
    character(len=LENCOMPONENTNAME) :: component_name
    integer(I4B) :: pnum
  contains
    procedure :: create => subpkg_create
    procedure :: add => subpkg_add
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
  subroutine subpkg_create(this, mempath, component_name)
    class(SubPackageListType) :: this
    character(len=*), intent(in) :: mempath
    character(len=*), intent(in) :: component_name

    ! initialize
    this%pnum = 0
    this%mempath = mempath
    this%component_name = component_name

    ! allocate arrays
    allocate (this%pkgtypes(0))
    allocate (this%component_types(0))
    allocate (this%subcomponent_types(0))
    allocate (this%filenames(0))
  end subroutine subpkg_create

  !> @brief create a new package type
  !<
  subroutine subpkg_add(this, pkgtype, component_type, subcomponent_type, &
                        tagname, filename)
    use ArrayHandlersModule, only: expandarray
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: mem_allocate
    use SimVariablesModule, only: idm_context
    class(SubPackageListType) :: this
    character(len=*), intent(in) :: pkgtype
    character(len=*), intent(in) :: component_type
    character(len=*), intent(in) :: subcomponent_type
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: filename
    character(len=LENVARNAME) :: mempath_tag
    character(len=LENMEMPATH), pointer :: subpkg_mempath
    character(len=LINELENGTH), pointer :: input_fname
    integer(I4B) :: idx, trimlen

    ! reallocate
    call expandarray(this%pkgtypes)
    call expandarray(this%component_types)
    call expandarray(this%subcomponent_types)
    call expandarray(this%filenames)

    ! add new package instance
    this%pnum = this%pnum + 1
    this%pkgtypes(this%pnum) = pkgtype
    this%component_types(this%pnum) = component_type
    this%subcomponent_types(this%pnum) = subcomponent_type
    this%filenames(this%pnum) = filename

    ! initialize mempath tag
    mempath_tag = tagname
    trimlen = len_trim(tagname)
    idx = 0

    ! create mempath tagname
    idx = index(tagname, '_')
    if (idx > 0) then
      if (tagname(idx + 1:trimlen) == 'FILENAME') then
        write (mempath_tag, '(a)') tagname(1:idx)//'MEMPATH'
      end if
    end if

    ! allocate mempath variable for subpackage
    call mem_allocate(subpkg_mempath, LENMEMPATH, mempath_tag, &
                      this%mempath)

    ! create and set the mempath
    subpkg_mempath = &
      create_mem_path(this%component_name, &
                      subcomponent_type, idm_context)

    ! allocate and initialize filename for subpackage
    call mem_allocate(input_fname, LINELENGTH, 'INPUT_FNAME', subpkg_mempath)
    input_fname = filename
  end subroutine subpkg_add

  !> @brief create a new package type
  !<
  subroutine subpkg_destroy(this)
    class(SubPackageListType) :: this
    ! allocate arrays
    deallocate (this%pkgtypes)
    deallocate (this%component_types)
    deallocate (this%subcomponent_types)
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
    call this%subpkg_list%create(this%mf6_input%mempath, &
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
    class(StaticPkgLoadType), intent(inout) :: this
    character(len=16), dimension(:), pointer :: subpkgs
    character(len=LINELENGTH) :: tag, fname, pkgtype
    character(len=LENFTYPE) :: c_type, sc_type
    character(len=16) :: subpkg
    integer(I4B) :: idx, n

    ! set pointer to package (idm integrated) subpackage list
    subpkgs => idm_subpackages(this%mf6_input%component_type, &
                               this%mf6_input%subcomponent_type)

    ! check if tag matches subpackage
    do n = 1, size(subpkgs)
      subpkg = subpkgs(n)
      idx = index(subpkg, '-')
      ! split sp string into component/subcomponent
      if (idx > 0) then
        ! split string in component/subcomponent types
        c_type = subpkg(1:idx - 1)
        sc_type = subpkg(idx + 1:len_trim(subpkg))
        if (idm_integrated(c_type, sc_type)) then
          ! set pkgtype and input filename tag
          pkgtype = trim(sc_type)//'6'
          tag = trim(pkgtype)//'_FILENAME'
          ! support single instance of each subpackage
          if (idm_multi_package(c_type, sc_type)) then
            errmsg = 'Multi-instance subpackages not supported. Remove dfn &
                     &subpackage tagline for package "'//trim(subpkg)//'".'
            call store_error(errmsg)
            call store_error_filename(this%input_name)
          else
            if (filein_fname(fname, tag, this%mf6_input%mempath, &
                             this%input_name)) then
              call this%subpkg_list%add(pkgtype, c_type, sc_type, &
                                        trim(tag), trim(fname))
            end if
          end if
        else
          errmsg = 'Identified subpackage is not IDM integrated. Remove dfn &
                   &subpackage tagline for package "'//trim(subpkg)//'".'
          call store_error(errmsg)
          call store_error_filename(this%input_name)
        end if
      end if
    end do
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
      ! no-op
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
