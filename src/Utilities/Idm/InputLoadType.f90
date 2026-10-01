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
    logical(LGP) :: has_keystring !< period block uses keystring-based dispatch
    integer(I4B) :: iperblock !< index of period block on block definition list
    integer(I4B) :: iout !< inunit number for logging
    integer(I4B) :: nparam !< number of in scope params
  contains
    procedure :: init => dynamic_init
    procedure :: df => dynamic_df
    procedure :: ts_advance => dynamic_ad
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
    procedure :: ts_advance => dynamicpkgs_ad
    procedure :: size => dynamicpkgs_size
    procedure :: destroy => dynamicpkgs_destroy
  end type ModelDynamicPkgsType

  type(ListType) :: model_inputs

contains

  !> @brief initialize a SubPackageListType object
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

  !> @brief append one subpackage file instance to the list
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

  !> @brief Assign subpackage names and mempaths for IDM-integrated subpackages.
  !<
  subroutine subpkg_names(this, parent_sctype, parent_scname, &
                          parent_mempath, modelfname)
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: mem_allocate
    use SimVariablesModule, only: idm_context
    use CharacterStringModule, only: CharacterStringType
    use SourceCommonModule, only: idm_utl_type
    class(SubPackageListType) :: this
    character(len=*), intent(in) :: parent_sctype
    character(len=*), intent(in) :: parent_scname
    character(len=*), intent(in) :: parent_mempath
    character(len=*), intent(in) :: modelfname
    character(len=LINELENGTH), dimension(:), allocatable :: subptypes
    integer(I4B), dimension(:), allocatable :: nsubptypes
    type(CharacterStringType), dimension(:), pointer, contiguous :: mempaths
    character(len=LINELENGTH), pointer :: input_fname
    character(len=LENVARNAME) :: mempath_key
    character(len=LENVARNAME) :: subpkg_prefix
    character(len=LENMEMPATH) :: mempath
    integer(I4B) :: subpkg_inst, n, m

    ! nothing to do if no subpackages were added
    if (size(this%pkgtypes) == 0) return

    ! UTL packages do not themselves have subpackages
    if (idm_utl_type(this%component_type, parent_sctype)) return

    ! build subpkg_prefix from the parent package identity
    subpkg_prefix = build_subpkg_prefix(this%component_type, &
                                        this%component_name, parent_sctype, &
                                        parent_scname, modelfname)

    ! deduplicate pkgtypes into unique list with per-type counts
    call deduplicate_pkgtypes(this%pkgtypes, subptypes, nsubptypes)

    ! allocate mempath arrays for each subpackage type, create and
    ! store the memory paths for package side access.
    do n = 1, size(subptypes)
      subpkg_inst = 0
      mempath_key = trim(subptypes(n))//'_MEMPATH'
      call mem_allocate(mempaths, LENMEMPATH, nsubptypes(n), &
                        mempath_key, parent_mempath)
      do m = 1, size(this%pkgtypes)
        if (this%pkgtypes(m) == subptypes(n)) then
          subpkg_inst = subpkg_inst + 1
          ! set the subpackage name
          write (this%subcomponent_names(m), '(a,i0)') &
            trim(subpkg_prefix)//trim(this%subcomponent_types(m)), subpkg_inst
          ! create and set mempath
          mempath = create_mem_path(this%component_name, &
                                    this%subcomponent_names(m), &
                                    idm_context)
          mempaths(subpkg_inst) = mempath
          ! create and set INPUT_FNAME string in each new memory path.
          call mem_allocate(input_fname, LINELENGTH, 'INPUT_FNAME', mempath)
          input_fname = trim(this%filenames(m))
        end if
      end do
    end do

    deallocate (subptypes)
    deallocate (nsubptypes)
  end subroutine subpkg_names

  !> @brief Build the subpackage name prefix for the given parent package.
  !!
  !! For single-instance parents (e.g. NPF), prefix is '<TYPE>-'.
  !! For multi-instance parents (e.g. WEL), prefix includes the instance
  !! number: '<TYPE><N>-' (e.g. 'WEL1-').  EXG packages return ''.
  !<
  function build_subpkg_prefix(component_type, component_name, &
                               parent_sctype, parent_scname, &
                               modelfname) result(subpkg_prefix)
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: mem_setptr
    use SimVariablesModule, only: idm_context
    use CharacterStringModule, only: CharacterStringType
    use ModelPackageInputModule, only: multi_package_type
    use IdmDfnSelectorModule, only: idm_multi_package, idm_integrated
    use SourceCommonModule, only: idm_pkg_instance_name
    character(len=*), intent(in) :: component_type
    character(len=*), intent(in) :: component_name
    character(len=*), intent(in) :: parent_sctype
    character(len=*), intent(in) :: parent_scname
    character(len=*), intent(in) :: modelfname
    character(len=LENVARNAME) :: subpkg_prefix
    type(CharacterStringType), dimension(:), contiguous, pointer :: pnames, ftypes
    character(len=LENVARNAME) :: parent_type, parent_ftype, parent_name
    character(len=LENMEMPATH) :: model_mempath
    integer(I4B) :: parent_inst, n
    logical(LGP) :: multi

    subpkg_prefix = ''

    ! EXG (exchange) packages have no model NAM file and don't need a prefix
    if (component_type == 'EXG') return

    ! resolve definition names to the namefile packages block type name
    select case (parent_sctype)
    case ('EVTA', 'RCHA', 'SPCA', 'RIVG', 'CHDG', 'WELG', 'DRNG', 'GHBG')
      parent_type = parent_sctype(1:3)
    case default
      parent_type = parent_sctype
    end select

    ! build the filetype string used to match FTYPE in the NAM packages block
    parent_ftype = trim(parent_type)//'6'

    ! determine if multi-package type
    if (idm_integrated(component_type, parent_type)) then
      multi = idm_multi_package(component_type, parent_type)
    else
      multi = multi_package_type(component_type, parent_type, parent_ftype)
    end if

    if (multi) then
      ! identify instance number of this package type in the namefile packages
      ! block and use to set subpackage prefix
      model_mempath = create_mem_path(component_name, 'NAM', idm_context)
      call mem_setptr(pnames, 'PNAME', model_mempath)
      call mem_setptr(ftypes, 'FTYPE', model_mempath)

      parent_inst = 0
      do n = 1, size(pnames)
        if (ftypes(n) == parent_ftype) then
          parent_inst = parent_inst + 1
          parent_name = pnames(n)
          if (parent_name == '') &
            parent_name = idm_pkg_instance_name(parent_type, parent_inst)
          if (parent_name == parent_scname) then
            write (subpkg_prefix, '(a,i0,a)') trim(parent_type), parent_inst, '-'
            exit
          end if
        end if
      end do

      if (subpkg_prefix == '') then
        errmsg = &
          'Internal IDM error: subpackage load cannot identify &
          &package "'//trim(parent_scname)//'" in model name file &
          &packages block.'
        call store_error(errmsg)
        call store_error_filename(modelfname)
      end if
    else
      ! single-instance parent: prefix is '<TYPE>-', e.g. 'NPF-'
      write (subpkg_prefix, '(2a)') trim(parent_type), '-'
    end if
  end function build_subpkg_prefix

  !> @brief Deduplicate pkgtypes into unique entries with counts (run-length encoding).
  !!
  !! INVARIANT: pkgtypes entries for the same type must be contiguous.
  !<
  subroutine deduplicate_pkgtypes(pkgtypes, subptypes, nsubptypes)
    use ArrayHandlersModule, only: expandarray
    character(len=LENCOMPONENTNAME), intent(in) :: pkgtypes(:)
    character(len=LINELENGTH), allocatable, intent(out) :: subptypes(:)
    integer(I4B), allocatable, intent(out) :: nsubptypes(:)
    character(len=LENCOMPONENTNAME) :: prev
    integer(I4B) :: n, ntype

    allocate (subptypes(0))
    allocate (nsubptypes(0))
    prev = ''
    ntype = 0
    do n = 1, size(pkgtypes)
      if (pkgtypes(n) /= prev) then
        ntype = ntype + 1
        prev = pkgtypes(n)
        call expandarray(subptypes)
        call expandarray(nsubptypes)
        subptypes(ntype) = prev
        nsubptypes(ntype) = 1
      else
        nsubptypes(ntype) = nsubptypes(ntype) + 1
      end if
    end do
  end subroutine deduplicate_pkgtypes

  !> @brief destroy a SubPackageListType object
  !<
  subroutine subpkg_destroy(this)
    class(SubPackageListType) :: this
    ! deallocate arrays
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
    use IdmDfnSelectorModule, only: idm_subpackages, idm_integrated
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

    ! check each subpackage type this package supports
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
            ! add all input files of this type to subpackage type list
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

    ! create subpackage names and use to store mempaths in memory manager
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
    use LoadContextModule, only: is_keystring_period
    class(DynamicPkgLoadType), intent(inout) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), intent(in) :: component_name
    character(len=*), intent(in) :: component_input_name
    character(len=*), intent(in) :: input_name
    integer(I4B), intent(in) :: iperblock
    integer(I4B), intent(in) :: iout
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: iparam

    this%mf6_input = mf6_input
    this%component_name = component_name
    this%component_input_name = component_input_name
    this%input_name = input_name
    this%readasarrays = .false.
    this%readarraygrid = .false.
    this%has_keystring = .false.
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

    ! detect keystring packages
    if (mf6_input%block_dfns(iperblock)%aggregate) then
      if (is_keystring_period(mf6_input)) this%has_keystring = .true.
    end if
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
      call dynamic_pkg%ts_advance()
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
