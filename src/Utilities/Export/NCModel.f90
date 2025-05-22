!> @brief This module contains the NCModelExportModule
!!
!! This module defines a model export and base type for
!! supported netcdf files and is not dependent on
!! netcdf libraries.
!!
!<
module NCModelExportModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENCOMPONENTNAME, LENMODELNAME, &
                             LENMEMPATH, LENBIGLINE, LENVARNAME, MVALIDATE, &
                             DIS, DISU, DISV
  use SimVariablesModule, only: isim_mode, idm_context, errmsg
  use SimModule, only: store_error, store_error_filename
  use InputLoadTypeModule, only: ModelDynamicPkgsType
  use ModflowInputModule, only: ModflowInputType
  use BoundInputContextModule, only: ReadStateVarType, rsv_name
  use ListModule, only: ListType

  implicit none
  private
  public :: NCBaseModelExportType, NCModelExportType
  public :: NCExportAnnotation
  public :: ExportPackageType
  public :: NETCDF_UNDEF, NETCDF_STRUCTURED, NETCDF_MESH2D
  public :: export_longname

  !> @brief netcdf export types enumerator
  !<
  ENUM, BIND(C)
    ENUMERATOR :: NETCDF_UNDEF = 0 !< undefined netcdf export type
    ENUMERATOR :: NETCDF_STRUCTURED = 1 !< netcdf structrured export
    ENUMERATOR :: NETCDF_MESH2D = 2 !< netcdf ugrid layered mesh export
  END ENUM

  type :: ExportPackageType
    type(ModflowInputType) :: mf6_input !< description of modflow6 input
    character(len=LINELENGTH), dimension(:), allocatable :: param_names !< dynamic param tagnames
    type(ReadStateVarType), dimension(:), allocatable :: param_reads !< param read states
    integer(I4B), dimension(:), pointer, contiguous :: mshape => null() !< model shape
    integer(I4B), pointer :: iper !< most recent package rp load
    integer(I4B) :: iper_export !< most recent period of netcdf package export
    integer(I4B) :: nparam !< number of in scope params
  contains
    procedure :: init => epkg_init
    procedure :: destroy => epkg_destroy
  end type ExportPackageType

  !> @brief netcdf export attribute annotations
  !<
  type :: NCExportAnnotation
    character(len=LINELENGTH) :: title !< file scoped title attribute
    character(len=LINELENGTH) :: model !< file scoped model attribute
    character(len=LINELENGTH) :: grid !< grid type
    character(len=LINELENGTH) :: history !< file scoped history attribute
    character(len=LINELENGTH) :: source !< file scoped source attribute
    character(len=LINELENGTH) :: conventions !< file scoped conventions attribute
    character(len=LINELENGTH) :: stdname !< dependent variable standard name
    character(len=LINELENGTH) :: longname !< dependent variable long name
  contains
    procedure :: set
  end type NCExportAnnotation

  !> @brief base class for an export model
  !<
  type :: NCModelExportType
    type(ListType) :: pkglist
    character(len=LENMODELNAME) :: modelname !< name of model
    character(len=LENCOMPONENTNAME) :: modeltype !< type of model
    character(len=LINELENGTH) :: modelfname !< name of model input file
    character(len=LINELENGTH) :: nc_fname !< name of netcdf export file
    character(len=LINELENGTH) :: gridmap_name !< name of grid mapping variable
    character(len=LINELENGTH) :: mesh_name = 'mesh' !< name of mesh container variable
    character(len=LENMEMPATH) :: dis_mempath !< discretization input mempath
    character(len=LENMEMPATH) :: ncf_mempath !< netcdf utility package input mempath
    character(len=LENBIGLINE) :: wkt !< wkt user string
    character(len=LINELENGTH) :: datetime !< export file creation time
    character(len=LINELENGTH) :: xname !< dependent variable name
    character(len=LINELENGTH) :: lenunits !< unidata udunits length units
    type(NCExportAnnotation) :: annotation !< export file annotation
    real(DP), dimension(:), pointer, contiguous :: x !< dependent variable pointer
    integer(I4B) :: disenum !< type of discretization
    integer(I4B) :: ncid !< netcdf file descriptor
    integer(I4B) :: stepcnt !< simulation step count
    integer(I4B) :: totnstp !< simulation total number of steps
    integer(I4B), pointer :: deflate !< variable deflate level
    integer(I4B), pointer :: shuffle !< variable shuffle filter
    integer(I4B), pointer :: input_attr !< assign variable input attr
    integer(I4B), pointer :: chunk_time !< chunking parameter for time dimension
    integer(I4B) :: iout !< lst file descriptor
    logical(LGP) :: chunking_active !< have chunking parameters been provided
  contains
    procedure :: init => export_init
    procedure :: get => export_get
    procedure :: input_attribute
    procedure :: destroy => export_destroy
  end type NCModelExportType

  !> @brief abstract type for model netcdf export type
  !<
  type, abstract, extends(NCModelExportType) :: NCBaseModelExportType
  contains
    procedure :: export_input
    procedure(model_define), deferred :: df
    procedure(model_step), deferred :: step
    procedure(package_export), deferred :: package_step
    procedure(package_export_ilayer), deferred :: package_step_ilayer
  end type NCBaseModelExportType

  !> @brief abstract interfaces for model netcdf export type
  !<
  abstract interface
    subroutine model_define(this)
      import NCBaseModelExportType
      class(NCBaseModelExportType), intent(inout) :: this
    end subroutine
    subroutine model_step(this)
      import NCBaseModelExportType
      class(NCBaseModelExportType), intent(inout) :: this
    end subroutine
    subroutine package_export(this, export_pkg)
      import NCBaseModelExportType, ExportPackageType
      class(NCBaseModelExportType), intent(inout) :: this
      class(ExportPackageType), pointer, intent(in) :: export_pkg
    end subroutine
    subroutine package_export_ilayer(this, export_pkg, ilayer_varname, &
                                     ilayer)
      import NCBaseModelExportType, ExportPackageType, I4B
      class(NCBaseModelExportType), intent(inout) :: this
      class(ExportPackageType), pointer, intent(in) :: export_pkg
      character(len=*), intent(in) :: ilayer_varname
      integer(I4B), intent(in) :: ilayer
    end subroutine
  end interface

contains

  !> @brief initialize dynamic package export object
  !<
  subroutine epkg_init(this, mf6_input, mshape, param_names, &
                       nparam)
    use SimVariablesModule, only: idm_context
    use MemoryManagerModule, only: mem_setptr
    use MemoryManagerExtModule, only: mem_set_value
    use MemoryHelperModule, only: create_mem_path
    class(ExportPackageType), intent(inout) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: mshape !< model shape
    character(len=LINELENGTH), dimension(:), allocatable, &
      intent(in) :: param_names
    integer(I4B), intent(in) :: nparam
    integer(I4B) :: n
    character(len=LENVARNAME) :: rs_varname
    character(len=LENMEMPATH) :: input_mempath
    integer(I4B), pointer :: rsvar

    this%mf6_input = mf6_input
    this%mshape => mshape
    this%nparam = nparam
    this%iper_export = 0

    input_mempath = create_mem_path(component=mf6_input%component_name, &
                                    subcomponent=mf6_input%subcomponent_name, &
                                    context=idm_context)

    ! allocate param arrays
    allocate (this%param_names(nparam))
    allocate (this%param_reads(nparam))

    ! set param arrays
    do n = 1, nparam
      this%param_names(n) = param_names(n)
      rs_varname = rsv_name(param_names(n))
      call mem_setptr(rsvar, rs_varname, mf6_input%mempath)
      this%param_reads(n)%invar => rsvar
    end do

    ! set pointer to loaded input period
    call mem_setptr(this%iper, 'IPER', mf6_input%mempath)
  end subroutine epkg_init

  !> @brief destroy dynamic package export object
  !<
  subroutine epkg_destroy(this)
    use InputDefinitionModule, only: InputParamDefinitionType
    class(ExportPackageType), intent(inout) :: this
    if (allocated(this%param_names)) deallocate (this%param_names)
  end subroutine epkg_destroy

  !> @brief set netcdf file scoped attributes
  !<
  subroutine set(this, modelname, modeltype, modelfname, nctype)
    use VersionModule, only: VERSION
    class(NCExportAnnotation), intent(inout) :: this
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelfname
    integer(I4B), intent(in) :: nctype
    character(len=LINELENGTH) :: fullname
    integer :: values(8)

    this%title = ''
    this%model = ''
    this%grid = ''
    this%history = ''
    this%source = ''
    this%conventions = ''
    this%stdname = ''
    this%longname = ''

    ! set file conventions
    this%conventions = 'CF-1.11'
    if (nctype == NETCDF_MESH2D) this%conventions = &
      trim(this%conventions)//' UGRID-1.0'

    ! set model specific attributes
    select case (modeltype)
    case ('GWF')
      fullname = 'Groundwater Flow'
      this%title = trim(modelname)//' hydraulic head'
      this%longname = 'head'
    case ('GWT')
      fullname = 'Groundwater Transport'
      this%title = trim(modelname)//' concentration'
      this%longname = 'concentration'
    case ('GWE')
      fullname = 'Groundwater Energy'
      this%title = trim(modelname)//' temperature'
      this%longname = 'temperature'
    case default
      errmsg = trim(modeltype)//' models not supported for NetCDF export.'
      call store_error(errmsg)
      call store_error_filename(modelfname)
    end select

    if (isim_mode == MVALIDATE) then
      this%title = trim(this%title)//' array input'
    end if

    ! set export type
    if (nctype == NETCDF_MESH2D) then
      this%grid = 'LAYERED MESH'
    else if (nctype == NETCDF_STRUCTURED) then
      this%grid = 'STRUCTURED'
    end if

    ! model description string
    this%model = trim(modelname)//': MODFLOW 6 '//trim(fullname)// &
                 ' ('//trim(modeltype)//') model'

    ! modflow6 version string
    this%source = 'MODFLOW 6 '//trim(adjustl(VERSION))

    ! create timestamp
    call date_and_time(values=values)
    write (this%history, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
      'first created ', values(1), '/', values(2), '/', values(3), ' ', &
      values(5), ':', values(6), ':', values(7), '.', values(8)
  end subroutine set

  !> @brief initialization of model netcdf export
  !<
  subroutine export_init(this, modelname, modeltype, modelfname, nc_fname, &
                         disenum, nctype, iout)
    use TdisModule, only: datetime0, nstp
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerExtModule, only: mem_set_value
    use InputOutputModule, only: lowcase
    use UtlNcfInputModule, only: UtlNcfParamFoundType
    class(NCModelExportType), intent(inout) :: this
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelfname
    character(len=*), intent(in) :: nc_fname
    integer(I4B), intent(in) :: disenum
    integer(I4B), intent(in) :: nctype
    integer(I4B), intent(in) :: iout
    character(len=LENMEMPATH) :: model_mempath
    type(UtlNcfParamFoundType) :: ncf_found
    logical(LGP) :: found_mempath

    ! allocate
    allocate (this%deflate)
    allocate (this%shuffle)
    allocate (this%input_attr)
    allocate (this%chunk_time)

    ! initialize
    this%modelname = modelname
    this%modeltype = modeltype
    this%modelfname = modelfname
    this%nc_fname = nc_fname
    this%gridmap_name = ''
    this%ncf_mempath = ''
    this%wkt = ''
    this%datetime = ''
    this%xname = ''
    this%lenunits = ''
    this%disenum = disenum
    this%ncid = 0
    this%stepcnt = 0
    this%totnstp = 0
    this%deflate = -1
    this%shuffle = 0
    this%input_attr = 1
    this%chunk_time = -1
    this%iout = iout
    this%chunking_active = .false.

    ! set file scoped attributes
    call this%annotation%set(modelname, modeltype, modelfname, nctype)

    ! set dependent variable basename
    select case (modeltype)
    case ('GWF')
      this%xname = 'head'
    case ('GWT')
      this%xname = 'concentration'
    case ('GWE')
      this%xname = 'temperature'
    case default
      errmsg = trim(modeltype)//' models not supported for NetCDF export.'
      call store_error(errmsg)
      call store_error_filename(modelfname)
    end select

    ! set discretization input mempath
    if (disenum == DIS) then
      this%dis_mempath = create_mem_path(modelname, 'DIS', idm_context)
    else if (disenum == DISU) then
      this%dis_mempath = create_mem_path(modelname, 'DISU', idm_context)
    else if (disenum == DISV) then
      this%dis_mempath = create_mem_path(modelname, 'DISV', idm_context)
    end if

    ! set dependent variable pointer
    model_mempath = create_mem_path(component=modelname)
    call mem_setptr(this%x, 'X', model_mempath)

    ! set ncf_mempath if provided
    call mem_set_value(this%ncf_mempath, 'NCF6_MEMPATH', this%dis_mempath, &
                       found_mempath)

    if (found_mempath) then
      call mem_set_value(this%wkt, 'WKT', this%ncf_mempath, &
                         ncf_found%wkt)
      call mem_set_value(this%deflate, 'DEFLATE', this%ncf_mempath, &
                         ncf_found%deflate)
      call mem_set_value(this%shuffle, 'SHUFFLE', this%ncf_mempath, &
                         ncf_found%shuffle)
      call mem_set_value(this%input_attr, 'ATTR_OFF', this%ncf_mempath, &
                         ncf_found%attr_off)
      call mem_set_value(this%chunk_time, 'CHUNK_TIME', this%ncf_mempath, &
                         ncf_found%chunk_time)
    end if

    if (ncf_found%wkt) then
      this%gridmap_name = 'projection'
    end if

    ! ATTR_OFF turns off modflow 6 input attributes
    if (ncf_found%attr_off) then
      this%input_attr = 0
    end if

    ! set datetime string
    if (datetime0 /= '') then
      this%datetime = 'days since '//trim(datetime0)
    else
      ! January 1, 1970 at 00:00:00 UTC
      this%datetime = 'days since 1970-01-01T00:00:00'
    end if

    ! set total nstp
    this%totnstp = sum(nstp)
  end subroutine export_init

  !> @brief retrieve dynamic export object from package list
  !<
  function export_get(this, idx) result(res)
    use ListModule, only: ListType
    class(NCModelExportType), intent(inout) :: this
    integer(I4B), intent(in) :: idx
    class(ExportPackageType), pointer :: res
    class(*), pointer :: obj
    nullify (res)
    obj => this%pkglist%GetItem(idx)
    if (associated(obj)) then
      select type (obj)
      class is (ExportPackageType)
        res => obj
      end select
    end if
  end function export_get

  !> @brief build modflow6_input attribute string
  !<
  function input_attribute(this, pkgname, idt) result(attr)
    use InputOutputModule, only: lowcase
    use MemoryHelperModule, only: memPathSeparator
    use InputDefinitionModule, only: InputParamDefinitionType
    class(NCModelExportType), intent(inout) :: this
    character(len=*), intent(in) :: pkgname
    type(InputParamDefinitionType), pointer, intent(in) :: idt
    character(len=LINELENGTH) :: attr
    attr = ''
    if (this%input_attr > 0) then
      attr = trim(this%modelname)//memPathSeparator//trim(pkgname)// &
             memPathSeparator//trim(idt%mf6varname)
    end if
  end function input_attribute

  !> @brief build netcdf variable longname
  !<
  function export_longname(longname, pkgname, tagname, layer, iper) result(lname)
    use InputOutputModule, only: lowcase
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    integer(I4B), intent(in) :: layer
    integer(I4B), optional, intent(in) :: iper
    character(len=LINELENGTH) :: lname
    character(len=LINELENGTH) :: pname, vname
    pname = pkgname
    vname = tagname
    call lowcase(pname)
    call lowcase(vname)
    if (longname == '') then
      lname = trim(pname)//' '//trim(vname)
    else
      lname = longname
    end if
    if (layer > 0) then
      write (lname, '(a,i0)') trim(lname)//' layer=', layer
    end if
    if (present(iper)) then
      if (iper > 0) then
        write (lname, '(a,i0)') trim(lname)//' period=', iper
      end if
    end if
  end function export_longname

  !> @brief netcdf dynamic package period export
  !<
  subroutine export_input(this)
    use TdisModule, only: kper
    use ArrayHandlersModule, only: ifind
    class(NCBaseModelExportType), intent(inout) :: this
    integer(I4B) :: idx, ilayer
    class(ExportPackageType), pointer :: export_pkg
    character(len=LENVARNAME) :: ilayer_varname

    do idx = 1, this%pkglist%Count()
      export_pkg => this%get(idx)
      ! last loaded data is not current period
      if (export_pkg%iper /= kper) cycle
      ! period input already exported
      if (export_pkg%iper_export >= export_pkg%iper) cycle
      ! set exported iper
      export_pkg%iper_export = export_pkg%iper

      ! initialize ilayer
      ilayer = 0

      ! set expected ilayer index variable name
      ilayer_varname = 'I'//trim(export_pkg%mf6_input%subcomponent_type(1:3))

      ! is ilayer variable in param name list
      ilayer = ifind(export_pkg%param_names, ilayer_varname)

      ! layer index variable is required to be first defined in period block
      if (ilayer == 1) then
        call this%package_step_ilayer(export_pkg, ilayer_varname, ilayer)
      else
        call this%package_step(export_pkg)
      end if
    end do
  end subroutine export_input

  !> @brief destroy model netcdf export object
  !<
  subroutine export_destroy(this)
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    class(NCModelExportType), intent(inout) :: this
    ! override in derived class
    deallocate (this%deflate)
    deallocate (this%shuffle)
    deallocate (this%input_attr)
    deallocate (this%chunk_time)
    ! Deallocate idm memory
    if (this%ncf_mempath /= '') then
      call memorystore_remove(this%modelname, 'NCF', idm_context)
    end if
  end subroutine export_destroy

end module NCModelExportModule
