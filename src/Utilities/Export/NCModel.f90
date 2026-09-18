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
                             DIS, DISU, DISV, DPIO180, DZERO
  use SimVariablesModule, only: isim_mode, idm_context, errmsg
  use SimModule, only: store_error, store_error_filename
  use InputLoadTypeModule, only: ModelDynamicPkgsType
  use ModflowInputModule, only: ModflowInputType
  use LoadContextModule, only: ReadStateVarType, rsv_name
  use ListModule, only: ListType

  implicit none
  private
  public :: NCBaseModelExportType, NCModelExportType
  public :: NCExportAnnotation
  public :: ExportPackageType
  public :: NETCDF_UNDEF, NETCDF_STRUCTURED, NETCDF_MESH2D
  public :: export_longname, export_varname
  public :: wkt_to_cf_gridmapping, wrap_rotated_crs
  public :: wkt_transverse_mercator_params
  public :: wkt_lambert_conformal_conic_params
  public :: wkt_albers_conical_equal_area_params

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
    integer(I4B), dimension(:, :), allocatable :: varids_param
    integer(I4B), dimension(:, :), allocatable :: varids_aux
    integer(I4B), dimension(:), pointer, contiguous :: mshape => null() !< model shape
    integer(I4B), pointer :: iper !< most recent package rp load
    integer(I4B) :: eper !< most recent period of netcdf package export
    integer(I4B) :: nparam !< number of in scope params
    integer(I4B) :: naux !< number of auxiliary variables
  contains
    procedure :: init => epkg_init
    procedure :: destroy => epkg_destroy
  end type ExportPackageType

  !> @brief netcdf export attribute annotations
  !<
  type :: NCExportAnnotation
    character(len=LINELENGTH) :: title !< file scoped title attribute
    character(len=LINELENGTH) :: model !< file scoped model attribute
    character(len=LINELENGTH) :: mesh !< mesh type
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
    character(len=LENBIGLINE) :: wkt !< WKT1 (OGC 01-009) user string
    character(len=LENBIGLINE) :: crs_wkt !< WKT2 (ISO 19162:2019) user string; falls back to wkt if empty
    character(len=LINELENGTH) :: datetime !< export file creation time
    character(len=LINELENGTH) :: xname !< dependent variable name
    character(len=LINELENGTH) :: lenunits !< unidata udunits length units
    type(NCExportAnnotation) :: annotation !< export file annotation
    real(DP), dimension(:), pointer, contiguous :: x !< dependent variable pointer
    integer(I4B) :: disenum !< type of discretization
    integer(I4B) :: ncid !< netcdf file descriptor
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
    procedure :: istp
    procedure :: destroy => export_destroy
  end type NCModelExportType

  !> @brief abstract type for model netcdf export type
  !<
  type, abstract, extends(NCModelExportType) :: NCBaseModelExportType
  contains
    procedure :: export_input
    procedure(model_define), deferred :: df
    procedure(package_export), deferred :: export_df
    procedure(model_step), deferred :: step
    procedure(package_export), deferred :: package_step
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
  subroutine epkg_init(this, mf6_input, mshape, naux, param_names, &
                       nparam)
    use SimVariablesModule, only: idm_context
    use MemoryManagerModule, only: mem_setptr
    use MemoryManagerExtModule, only: mem_set_value
    use MemoryHelperModule, only: create_mem_path
    class(ExportPackageType), intent(inout) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: mshape !< model shape
    integer(I4B), intent(in) :: naux
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
    this%naux = naux
    this%eper = 0

    input_mempath = create_mem_path(component=mf6_input%component_name, &
                                    subcomponent=mf6_input%subcomponent_name, &
                                    context=idm_context)

    ! allocate param arrays
    allocate (this%param_names(nparam))
    allocate (this%param_reads(nparam))
    allocate (this%varids_param(nparam, mshape(1)))
    allocate (this%varids_aux(naux, mshape(1)))

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
  subroutine set(this, modelname, modeltype, modelfname, nctype, disenum)
    use VersionModule, only: FULLVERSION
    use InputOutputModule, only: lowcase
    class(NCExportAnnotation), intent(inout) :: this
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelfname
    integer(I4B), intent(in) :: nctype
    integer(I4B), intent(in) :: disenum
    integer :: values(8)

    this%title = ''
    this%model = ''
    this%mesh = ''
    this%grid = ''
    this%history = ''
    this%source = ''
    this%conventions = ''
    this%stdname = ''
    this%longname = ''

    ! set file conventions
    this%conventions = 'CF-1.13'
    if (nctype == NETCDF_MESH2D) this%conventions = &
      trim(this%conventions)//' UGRID-1.0'

    ! set model specific attributes
    select case (modeltype)
    case ('GWF')
      this%title = trim(modelname)//' hydraulic head'
      this%longname = 'head'
    case ('GWT')
      this%title = trim(modelname)//' concentration'
      this%longname = 'concentration'
    case ('GWE')
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

    ! set mesh type
    if (nctype == NETCDF_MESH2D) then
      this%mesh = 'layered'
    end if

    ! set grid type
    if (disenum == DIS) then
      this%grid = 'structured'
    else if (disenum == DISV) then
      this%grid = 'vertex'
    end if

    ! model description string (lowercase, no version suffix, matching CF convention)
    this%model = trim(modeltype)//': '//trim(modelname)
    call lowcase(this%model)

    ! modflow6 version string
    this%source = 'MODFLOW 6 '//trim(adjustl(FULLVERSION))

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
    use TdisModule, only: datetime0, nper, nstp
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerExtModule, only: mem_set_value
    use InputOutputModule, only: lowcase
    use SourceCommonModule, only: filein_fname
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
    this%crs_wkt = ''
    this%datetime = ''
    this%xname = ''
    this%lenunits = ''
    this%disenum = disenum
    this%ncid = 0
    this%totnstp = 0
    this%deflate = -1
    this%shuffle = 0
    this%input_attr = 1
    this%chunk_time = -1
    this%iout = iout
    this%chunking_active = .false.

    ! set file scoped attributes
    call this%annotation%set(modelname, modeltype, modelfname, nctype, disenum)

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
    if (filein_fname(this%ncf_mempath, 'NCF6_MEMPATH', this%dis_mempath, &
                     modelfname)) then
      call mem_set_value(this%wkt, 'WKT', this%ncf_mempath, &
                         ncf_found%wkt)
      call mem_set_value(this%crs_wkt, 'CRS_WKT', this%ncf_mempath, &
                         ncf_found%crs_wkt)
      call mem_set_value(this%deflate, 'DEFLATE', this%ncf_mempath, &
                         ncf_found%deflate)
      call mem_set_value(this%shuffle, 'SHUFFLE', this%ncf_mempath, &
                         ncf_found%shuffle)
      call mem_set_value(this%input_attr, 'ATTR_OFF', this%ncf_mempath, &
                         ncf_found%attr_off)
      call mem_set_value(this%chunk_time, 'CHUNK_TIME', this%ncf_mempath, &
                         ncf_found%chunk_time)
    end if

    if (ncf_found%wkt .or. ncf_found%crs_wkt) then
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
    if (isim_mode == MVALIDATE) then
      this%totnstp = nper
    else
      this%totnstp = sum(nstp)
    end if
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

  !> @brief build modflow_input attribute string
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
             memPathSeparator//trim(idt%tagname)
      call lowcase(attr)
    end if
  end function input_attribute

  !> @brief step index for timeseries data
  !<
  function istp(this)
    use TdisModule, only: kstp, kper, nstp
    class(NCModelExportType), intent(inout) :: this
    integer(I4B) :: n, istp
    istp = kstp
    if (kper > 1) then
      do n = 1, kper - 1
        istp = istp + nstp(n)
      end do
    end if
  end function istp

  !> @brief build netcdf variable name
  !<
  function export_varname(pkgname, tagname, mempath, layer, iaux) &
    result(varname)
    use MemoryManagerModule, only: mem_setptr
    use CharacterStringModule, only: CharacterStringType
    use InputOutputModule, only: lowcase
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: mempath
    integer(I4B), optional, intent(in) :: layer
    integer(I4B), optional, intent(in) :: iaux
    character(len=LINELENGTH) :: varname
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: auxnames
    character(len=LINELENGTH) :: pname, vname
    vname = tagname
    pname = pkgname

    if (present(iaux)) then
      if (iaux > 0) then
        if (tagname == 'AUX') then
          ! reset vname to auxiliary variable name
          call mem_setptr(auxnames, 'AUXILIARY', mempath)
          vname = auxnames(iaux)
        end if
      end if
    end if

    call lowcase(vname)
    call lowcase(pname)
    varname = trim(pname)//'_'//trim(vname)

    if (present(layer)) then
      if (layer > 0) then
        !write (varname, '(a,i0)') trim(varname)//'_L', layer
        write (varname, '(a,i0)') trim(varname)//'_l', layer
      end if
    end if
  end function export_varname

  !> @brief build netcdf variable longname
  !<
  function export_longname(longname, pkgname, tagname, mempath, layer, iaux, &
                           component_type, subcomponent_type) result(lname)
    use MemoryManagerModule, only: mem_setptr
    use CharacterStringModule, only: CharacterStringType
    use InputOutputModule, only: lowcase
    use IdmDfnSelectorModule, only: idm_multi_package
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: mempath
    integer(I4B), optional, intent(in) :: layer
    integer(I4B), optional, intent(in) :: iaux
    character(len=*), optional, intent(in) :: component_type
    character(len=*), optional, intent(in) :: subcomponent_type
    character(len=LINELENGTH) :: lname
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: auxnames
    character(len=LINELENGTH) :: pname, vname, auxname
    pname = pkgname
    vname = tagname
    call lowcase(pname)
    call lowcase(vname)
    if (longname == '') then
      lname = trim(pname)//' '//trim(vname)
    else
      lname = longname
      if (present(component_type) .and. present(subcomponent_type)) then
        if (idm_multi_package(component_type, subcomponent_type)) then
          lname = trim(pname)//' '//trim(lname)
        end if
      end if
    end if

    if (present(iaux)) then
      if (iaux > 0) then
        if (tagname == 'AUX') then
          ! reset vname to auxiliary variable name
          call mem_setptr(auxnames, 'AUXILIARY', mempath)
          auxname = auxnames(iaux)
          call lowcase(auxname)
          lname = trim(lname)//' '//trim(auxname)
        end if
      end if
    end if

    if (present(layer)) then
      if (layer > 0) then
        write (lname, '(a,i0)') trim(lname)//' layer ', layer
      end if
    end if
  end function export_longname

  !> @brief netcdf dynamic package period export
  !<
  subroutine export_input(this)
    use TdisModule, only: kper
    class(NCBaseModelExportType), intent(inout) :: this
    integer(I4B) :: idx
    class(ExportPackageType), pointer :: export_pkg
    do idx = 1, this%pkglist%Count()
      export_pkg => this%get(idx)
      ! period input already exported
      if (export_pkg%eper >= kper) cycle
      ! update export package
      call this%package_step(export_pkg)
      ! update exported iper
      export_pkg%eper = kper
    end do
  end subroutine export_input

  !> @brief map a WKT string to a CF-1.13 grid_mapping_name
  !!
  !! Accepts either WKT1 (OGC 01-009, PROJECTION["name"]) or WKT2
  !! (ISO 19162:2019, METHOD["name"]) format.  WKT1 is tried first;
  !! if PROJECTION[ is absent the WKT2 METHOD[ path is attempted.
  !! Returns '' for geographic CRS or unrecognised projection names.
  !! Only common groundwater projections included.  Matching is
  !! case-insensitive.
  !<
  function wkt_to_cf_gridmapping(wkt) result(gmname)
    use InputOutputModule, only: upcase
    character(len=*), intent(in) :: wkt
    character(len=LINELENGTH) :: gmname
    character(len=LINELENGTH) :: proj_name
    character(len=LENBIGLINE) :: wkt_upper
    integer :: istart, iend

    gmname = ''
    proj_name = ''

    ! uppercase copy for case-insensitive keyword search
    wkt_upper = wkt
    call upcase(wkt_upper)

    ! --- WKT1 path: PROJECTION["name"] ---
    istart = index(wkt_upper, 'PROJECTION[')
    if (istart /= 0) then
      istart = istart + len('PROJECTION[')
      do while (istart <= len(wkt))
        if (wkt(istart:istart) == '"') exit
        istart = istart + 1
      end do
      if (istart > len(wkt)) return
      istart = istart + 1
      iend = index(wkt(istart:), '"')
      if (iend == 0) return
      iend = istart + iend - 2
      if (iend < istart) return
      proj_name = wkt(istart:iend)
      call upcase(proj_name)
      select case (trim(proj_name))
      case ('TRANSVERSE_MERCATOR')
        gmname = 'transverse_mercator'
      case ('LAMBERT_CONFORMAL_CONIC_2SP', 'LAMBERT_CONFORMAL_CONIC_1SP')
        gmname = 'lambert_conformal_conic'
      case ('ALBERS_CONIC_EQUAL_AREA')
        gmname = 'albers_conical_equal_area'
      case ('MERCATOR_1SP', 'MERCATOR_2SP')
        gmname = 'mercator'
      case ('POLAR_STEREOGRAPHIC')
        gmname = 'polar_stereographic'
      end select
      return
    end if

    ! --- WKT2 path: METHOD["name"] ---
    istart = index(wkt_upper, 'METHOD[')
    if (istart == 0) return
    istart = istart + len('METHOD[')
    do while (istart <= len(wkt))
      if (wkt(istart:istart) == '"') exit
      istart = istart + 1
    end do
    if (istart > len(wkt)) return
    istart = istart + 1
    iend = index(wkt(istart:), '"')
    if (iend == 0) return
    iend = istart + iend - 2
    if (iend < istart) return
    proj_name = wkt(istart:iend)
    call upcase(proj_name)
    select case (trim(proj_name))
    case ('TRANSVERSE MERCATOR')
      gmname = 'transverse_mercator'
    case ('LAMBERT CONIC CONFORMAL (2SP)', 'LAMBERT CONIC CONFORMAL (1SP)')
      gmname = 'lambert_conformal_conic'
    case ('ALBERS EQUAL AREA')
      gmname = 'albers_conical_equal_area'
    case ('MERCATOR (VARIANT A)', 'MERCATOR (VARIANT B)')
      gmname = 'mercator'
    case ('POLAR STEREOGRAPHIC (VARIANT A)', 'POLAR STEREOGRAPHIC (VARIANT B)')
      gmname = 'polar_stereographic'
    end select
  end function wkt_to_cf_gridmapping

  !> @brief find the first number following a quoted parameter name
  !!
  !! Matches e.g. PARAMETER["central_meridian",-117] (WKT1) or
  !! PARAMETER["Longitude of natural origin",-117,...] (WKT2).
  !! Best-effort only -- MF6 has no CRS parsing library.  found=.false.
  !! if name is not present or no number follows it.
  !<
  function wkt_find_param(wkt, name, found) result(value)
    use InputOutputModule, only: upcase
    character(len=*), intent(in) :: wkt
    character(len=*), intent(in) :: name
    logical(LGP), intent(out) :: found
    real(DP) :: value
    character(len=LENBIGLINE) :: wkt_upper
    character(len=LINELENGTH) :: target_upper
    integer :: ipos, istart, iend, ios, wkt_len

    found = .false.
    value = DZERO
    wkt_len = len_trim(wkt)
    target_upper = '"'//trim(name)//'"'
    wkt_upper = wkt
    call upcase(wkt_upper)
    call upcase(target_upper)
    ipos = index(wkt_upper(1:wkt_len), trim(target_upper))
    if (ipos == 0) return

    istart = ipos + len_trim(target_upper)
    do while (istart <= wkt_len)
      if (wkt(istart:istart) == ',') exit
      istart = istart + 1
    end do
    istart = istart + 1
    iend = istart
    do while (iend <= wkt_len)
      if (wkt(iend:iend) == ',' .or. wkt(iend:iend) == ']') exit
      iend = iend + 1
    end do
    if (iend <= istart) return
    read (wkt(istart:iend - 1), *, iostat=ios) value
    if (ios == 0) found = .true.
  end function wkt_find_param

  !> @brief find semi-major axis and inverse flattening from a WKT string
  !!
  !! Matches SPHEROID["name",A,B,...] (WKT1) or ELLIPSOID["name",A,B,...]
  !! (WKT2).  Best-effort only.  found=.false. if neither keyword is
  !! present or the two numbers cannot be read.
  !<
  subroutine wkt_find_ellipsoid(wkt, semi_major_axis, inverse_flattening, found)
    use InputOutputModule, only: upcase
    character(len=*), intent(in) :: wkt
    real(DP), intent(out) :: semi_major_axis
    real(DP), intent(out) :: inverse_flattening
    logical(LGP), intent(out) :: found
    character(len=LENBIGLINE) :: wkt_upper
    integer :: ipos, i, iend, ios, wkt_len, nquote

    found = .false.
    semi_major_axis = DZERO
    inverse_flattening = DZERO
    wkt_len = len_trim(wkt)
    wkt_upper = wkt
    call upcase(wkt_upper)

    ipos = index(wkt_upper(1:wkt_len), 'SPHEROID[')
    if (ipos == 0) ipos = index(wkt_upper(1:wkt_len), 'ELLIPSOID[')
    if (ipos == 0) return

    ! skip past the quoted ellipsoid name (its closing quote)
    i = ipos
    nquote = 0
    do while (i <= wkt_len)
      if (wkt(i:i) == '"') nquote = nquote + 1
      if (nquote == 2) exit
      i = i + 1
    end do
    if (nquote < 2) return
    i = i + 1
    if (wkt(i:i) == ',') i = i + 1

    iend = i
    do while (iend <= wkt_len)
      if (wkt(iend:iend) == ',') exit
      iend = iend + 1
    end do
    if (iend <= i) return
    read (wkt(i:iend - 1), *, iostat=ios) semi_major_axis
    if (ios /= 0) return

    i = iend + 1
    iend = i
    do while (iend <= wkt_len)
      if (wkt(iend:iend) == ',' .or. wkt(iend:iend) == ']') exit
      iend = iend + 1
    end do
    if (iend <= i) return
    read (wkt(i:iend - 1), *, iostat=ios) inverse_flattening
    if (ios == 0) found = .true.
  end subroutine wkt_find_ellipsoid

  !> @brief extract CF transverse_mercator grid_mapping parameters
  !!
  !! Best-effort extraction from a WKT1 or WKT2 string -- MF6 has no
  !! CRS parsing library.  Tries WKT1 parameter names first, falling
  !! back to WKT2 names.  found=.false. if any parameter cannot be
  !! located.
  !<
  subroutine wkt_transverse_mercator_params(wkt, longitude_of_central_meridian, &
                                            latitude_of_projection_origin, &
                                            scale_factor_at_central_meridian, &
                                            false_easting, false_northing, &
                                            semi_major_axis, inverse_flattening, &
                                            found)
    character(len=*), intent(in) :: wkt
    real(DP), intent(out) :: longitude_of_central_meridian
    real(DP), intent(out) :: latitude_of_projection_origin
    real(DP), intent(out) :: scale_factor_at_central_meridian
    real(DP), intent(out) :: false_easting
    real(DP), intent(out) :: false_northing
    real(DP), intent(out) :: semi_major_axis
    real(DP), intent(out) :: inverse_flattening
    logical(LGP), intent(out) :: found
    logical(LGP) :: f1, f2, f3, f4, f5, f6

    longitude_of_central_meridian = wkt_find_param(wkt, 'central_meridian', f1)
    if (.not. f1) longitude_of_central_meridian = &
      wkt_find_param(wkt, 'Longitude of natural origin', f1)

    latitude_of_projection_origin = wkt_find_param(wkt, 'latitude_of_origin', f2)
    if (.not. f2) latitude_of_projection_origin = &
      wkt_find_param(wkt, 'Latitude of natural origin', f2)

    scale_factor_at_central_meridian = wkt_find_param(wkt, 'scale_factor', f3)
    if (.not. f3) scale_factor_at_central_meridian = &
      wkt_find_param(wkt, 'Scale factor at natural origin', f3)

    false_easting = wkt_find_param(wkt, 'false_easting', f4)
    if (.not. f4) false_easting = wkt_find_param(wkt, 'False easting', f4)

    false_northing = wkt_find_param(wkt, 'false_northing', f5)
    if (.not. f5) false_northing = wkt_find_param(wkt, 'False northing', f5)

    call wkt_find_ellipsoid(wkt, semi_major_axis, inverse_flattening, f6)

    found = f1 .and. f2 .and. f3 .and. f4 .and. f5 .and. f6
  end subroutine wkt_transverse_mercator_params

  !> @brief extract CF lambert_conformal_conic grid_mapping parameters
  !!
  !! Best-effort extraction from a WKT1 or WKT2 string -- MF6 has no
  !! CRS parsing library.  2SP (two standard parallels) only: CF's
  !! lambert_conformal_conic has no scale_factor attribute, so a 1SP
  !! (scale-factor-based) WKT cannot be exactly represented without
  !! trigonometric conversion MF6 has no library for.  is_2sp=.false.
  !! if a second standard parallel is not present (i.e. likely 1SP);
  !! found=.false. if is_2sp but any other parameter cannot be located.
  !<
  subroutine wkt_lambert_conformal_conic_params(wkt, standard_parallel_1, &
                                                standard_parallel_2, &
                                                longitude_of_central_meridian, &
                                                latitude_of_projection_origin, &
                                                false_easting, false_northing, &
                                                semi_major_axis, &
                                                inverse_flattening, is_2sp, found)
    character(len=*), intent(in) :: wkt
    real(DP), intent(out) :: standard_parallel_1
    real(DP), intent(out) :: standard_parallel_2
    real(DP), intent(out) :: longitude_of_central_meridian
    real(DP), intent(out) :: latitude_of_projection_origin
    real(DP), intent(out) :: false_easting
    real(DP), intent(out) :: false_northing
    real(DP), intent(out) :: semi_major_axis
    real(DP), intent(out) :: inverse_flattening
    logical(LGP), intent(out) :: is_2sp
    logical(LGP), intent(out) :: found
    logical(LGP) :: f1, f2, f3, f4, f5, f6, f7

    standard_parallel_2 = wkt_find_param(wkt, 'standard_parallel_2', f2)
    if (.not. f2) standard_parallel_2 = &
      wkt_find_param(wkt, 'Latitude of 2nd standard parallel', f2)
    is_2sp = f2
    if (.not. is_2sp) then
      found = .false.
      return
    end if

    standard_parallel_1 = wkt_find_param(wkt, 'standard_parallel_1', f1)
    if (.not. f1) standard_parallel_1 = &
      wkt_find_param(wkt, 'Latitude of 1st standard parallel', f1)

    longitude_of_central_meridian = wkt_find_param(wkt, 'central_meridian', f3)
    if (.not. f3) longitude_of_central_meridian = &
      wkt_find_param(wkt, 'Longitude of false origin', f3)

    latitude_of_projection_origin = wkt_find_param(wkt, 'latitude_of_origin', f4)
    if (.not. f4) latitude_of_projection_origin = &
      wkt_find_param(wkt, 'Latitude of false origin', f4)

    false_easting = wkt_find_param(wkt, 'false_easting', f5)
    if (.not. f5) false_easting = &
      wkt_find_param(wkt, 'Easting at false origin', f5)

    false_northing = wkt_find_param(wkt, 'false_northing', f6)
    if (.not. f6) false_northing = &
      wkt_find_param(wkt, 'Northing at false origin', f6)

    call wkt_find_ellipsoid(wkt, semi_major_axis, inverse_flattening, f7)

    found = f1 .and. f3 .and. f4 .and. f5 .and. f6 .and. f7
  end subroutine wkt_lambert_conformal_conic_params

  !> @brief extract CF albers_conical_equal_area grid_mapping parameters
  !!
  !! Best-effort extraction from a WKT1 or WKT2 string -- MF6 has no
  !! CRS parsing library.  found=.false. if any parameter cannot be
  !! located.
  !<
  subroutine wkt_albers_conical_equal_area_params(wkt, standard_parallel_1, &
                                                  standard_parallel_2, &
                                                  longitude_of_central_meridian, &
                                                  latitude_of_projection_origin, &
                                                  false_easting, &
                                                  false_northing, &
                                                  semi_major_axis, &
                                                  inverse_flattening, found)
    character(len=*), intent(in) :: wkt
    real(DP), intent(out) :: standard_parallel_1
    real(DP), intent(out) :: standard_parallel_2
    real(DP), intent(out) :: longitude_of_central_meridian
    real(DP), intent(out) :: latitude_of_projection_origin
    real(DP), intent(out) :: false_easting
    real(DP), intent(out) :: false_northing
    real(DP), intent(out) :: semi_major_axis
    real(DP), intent(out) :: inverse_flattening
    logical(LGP), intent(out) :: found
    logical(LGP) :: f1, f2, f3, f4, f5, f6, f7

    standard_parallel_1 = wkt_find_param(wkt, 'standard_parallel_1', f1)
    if (.not. f1) standard_parallel_1 = &
      wkt_find_param(wkt, 'Latitude of 1st standard parallel', f1)

    standard_parallel_2 = wkt_find_param(wkt, 'standard_parallel_2', f2)
    if (.not. f2) standard_parallel_2 = &
      wkt_find_param(wkt, 'Latitude of 2nd standard parallel', f2)

    longitude_of_central_meridian = wkt_find_param(wkt, 'longitude_of_center', f3)
    if (.not. f3) longitude_of_central_meridian = &
      wkt_find_param(wkt, 'Longitude of false origin', f3)

    latitude_of_projection_origin = wkt_find_param(wkt, 'latitude_of_center', f4)
    if (.not. f4) latitude_of_projection_origin = &
      wkt_find_param(wkt, 'Latitude of false origin', f4)

    false_easting = wkt_find_param(wkt, 'false_easting', f5)
    if (.not. f5) false_easting = &
      wkt_find_param(wkt, 'Easting at false origin', f5)

    false_northing = wkt_find_param(wkt, 'false_northing', f6)
    if (.not. f6) false_northing = &
      wkt_find_param(wkt, 'Northing at false origin', f6)

    call wkt_find_ellipsoid(wkt, semi_major_axis, inverse_flattening, f7)

    found = f1 .and. f2 .and. f3 .and. f4 .and. f5 .and. f6 .and. f7
  end subroutine wkt_albers_conical_equal_area_params

  !> @brief wrap a WKT2 PROJCRS in a derived CRS encoding grid rotation
  !!
  !! Builds a DerivedProjectedCRS wrapping the given WKT2 PROJCRS in an
  !! EPSG:9624 (affine parametric transformation) deriving conversion
  !! parameterized from xorigin/yorigin/angrot (degrees). Per ISO 19111,
  !! a deriving conversion is directed base->derived, so the parameters
  !! encode the world->local (inverse) rotation; a CRS-aware consumer
  !! applies the inverse to resolve true position from local coordinates.
  !! Returns '' if wkt is not recognized as a WKT2 PROJCRS -- a
  !! best-effort structural check only, as MF6 has no CRS parsing
  !! library and does not validate the WKT beyond this.
  !<
  function wrap_rotated_crs(wkt, xorigin, yorigin, angrot) result(derived_wkt)
    use InputOutputModule, only: upcase
    character(len=*), intent(in) :: wkt
    real(DP), intent(in) :: xorigin, yorigin, angrot
    character(len=LENBIGLINE) :: derived_wkt
    character(len=LENBIGLINE) :: wkt_trim, wkt_upper, base_body, conversion
    character(len=30) :: a0s, a1s, a2s, b0s, b1s, b2s
    integer :: nstart, i, depth, closures, iend, wkt_len
    real(DP) :: ang, a0, a1, a2, b0, b1, b2

    derived_wkt = ''

    ! best-effort check that wkt is a WKT2 PROJCRS
    wkt_trim = trim(adjustl(wkt))
    wkt_len = len_trim(wkt_trim)
    wkt_upper = wkt_trim
    call upcase(wkt_upper)
    if (wkt_len < 9 .or. wkt_upper(1:8) /= 'PROJCRS[') return

    ! locate the end of PROJCRS's 2nd top-level child (the mandatory
    ! BASEGEOGCRS/BASEGEODCRS and CONVERSION elements, in that order),
    ! by tracking bracket depth from just inside the opening '['
    nstart = 9
    depth = 0
    closures = 0
    iend = 0
    do i = nstart, wkt_len
      if (wkt_trim(i:i) == '[') then
        depth = depth + 1
      else if (wkt_trim(i:i) == ']') then
        depth = depth - 1
        if (depth == 0) then
          closures = closures + 1
          if (closures == 2) then
            iend = i
            exit
          end if
        end if
      end if
    end do
    if (iend == 0) return
    base_body = wkt_trim(nstart:iend)

    ! EPSG:9624 affine parameters, base(world) -> derived(local grid)
    ang = angrot * DPIO180
    a0 = -(xorigin * cos(ang) + yorigin * sin(ang))
    a1 = cos(ang)
    a2 = sin(ang)
    b0 = xorigin * sin(ang) - yorigin * cos(ang)
    b1 = -sin(ang)
    b2 = cos(ang)
    write (a0s, '(es22.15)') a0
    write (a1s, '(es22.15)') a1
    write (a2s, '(es22.15)') a2
    write (b0s, '(es22.15)') b0
    write (b1s, '(es22.15)') b1
    write (b2s, '(es22.15)') b2

    conversion = 'DERIVINGCONVERSION["MODFLOW 6 grid rotation",'// &
                 'METHOD["Affine parametric transformation",'// &
                 'ID["EPSG",9624]],'// &
                 'PARAMETER["A0",'//trim(adjustl(a0s))// &
                 ',LENGTHUNIT["metre",1]],'// &
                 'PARAMETER["A1",'//trim(adjustl(a1s))// &
                 ',SCALEUNIT["unity",1]],'// &
                 'PARAMETER["A2",'//trim(adjustl(a2s))// &
                 ',SCALEUNIT["unity",1]],'// &
                 'PARAMETER["B0",'//trim(adjustl(b0s))// &
                 ',LENGTHUNIT["metre",1]],'// &
                 'PARAMETER["B1",'//trim(adjustl(b1s))// &
                 ',SCALEUNIT["unity",1]],'// &
                 'PARAMETER["B2",'//trim(adjustl(b2s))// &
                 ',SCALEUNIT["unity",1]]]'

    derived_wkt = 'DERIVEDPROJCRS["MODFLOW 6 rotated grid CRS",'// &
                  'BASEPROJCRS['//trim(base_body)//'],'// &
                  trim(conversion)//','// &
                  'CS[Cartesian,2],'// &
                  'AXIS["easting (X)",east,ORDER[1],LENGTHUNIT["metre",1]],'// &
                  'AXIS["northing (Y)",north,ORDER[2],LENGTHUNIT["metre",1]]]'
  end function wrap_rotated_crs

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
