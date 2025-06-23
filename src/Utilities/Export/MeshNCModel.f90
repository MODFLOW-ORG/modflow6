!> @brief This module contains the MeshModelModule
!!
!! This module defines a base class for UGRID based
!! (mesh) model netcdf exports. It is dependent on
!! external netcdf libraries.
!<
module MeshModelModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENCOMPONENTNAME, LENMEMPATH, &
                             DNODATA, DHNOFLO
  use SimVariablesModule, only: errmsg, warnmsg
  use SimModule, only: store_error, store_warning, store_error_filename
  use MemoryManagerModule, only: mem_setptr
  use InputDefinitionModule, only: InputParamDefinitionType
  use CharacterStringModule, only: CharacterStringType
  use NCModelExportModule, only: export_longname, export_varname, &
                                 NCBaseModelExportType
  use NetCDFCommonModule, only: nf_verify
  use netcdf

  implicit none
  private
  public :: MeshNCDimIdType, MeshNCVarIdType
  public :: Mesh2dModelType
  public :: ncvar_chunk
  public :: ncvar_deflate
  public :: ncvar_gridmap
  public :: ncvar_mf6attr

  !> @brief type for storing model export dimension ids
  !<
  type :: MeshNCDimIdType
    integer(I4B) :: nmesh_node !< number of nodes in mesh
    integer(I4B) :: nmesh_face !< number of faces in mesh
    integer(I4B) :: max_nmesh_face_nodes !< max number of nodes in a single face
    integer(I4B) :: nlay !< number of layers
    integer(I4B) :: time !< number of steps
  contains
  end type MeshNCDimIdType

  !> @brief type for storing model export variable ids
  !<
  type :: MeshNCVarIdType
    integer(I4B) :: mesh !< mesh container variable
    integer(I4B) :: mesh_node_x !< mesh nodes x array
    integer(I4B) :: mesh_node_y !< mesh nodes y array
    integer(I4B) :: mesh_face_x !< mesh faces x location array
    integer(I4B) :: mesh_face_y !< mesh faces y location array
    integer(I4B) :: mesh_face_xbnds !< mesh faces 2D x bounds array
    integer(I4B) :: mesh_face_ybnds !< mesh faces 2D y bounds array
    integer(I4B) :: mesh_face_nodes !< mesh faces 2D nodes array
    integer(I4B) :: time !< time coordinate variable
    integer(I4B), dimension(:), allocatable :: export !< in scope layer export
    integer(I4B), dimension(:), allocatable :: dependent !< layered dependent variables array
  contains
  end type MeshNCVarIdType

  !> @brief base ugrid netcdf export type
  !<
  type, abstract, extends(NCBaseModelExportType) :: MeshModelType
    type(MeshNCDimIdType) :: dim_ids !< dimension ids
    type(MeshNCVarIdType) :: var_ids !< variable ids
    integer(I4B) :: nlay !< number of layers
    integer(I4B), pointer :: chunk_face !< chunking parameter for face dimension
  contains
    procedure :: mesh_init
    procedure :: mesh_destroy
    procedure :: df_export
    procedure :: export_df
    procedure :: create_timeseries
    procedure :: add_global_att
    procedure(nc_array_export_if), deferred :: export_input_array
    procedure :: export_input_arrays
    procedure :: add_pkg_data
    procedure :: define_dependent
    procedure :: define_gridmap
  end type MeshModelType

  !> @brief abstract interfaces for derived ugrid netcd export types
  !<
  abstract interface
    subroutine nc_array_export_if(this, pkgtype, pkgname, mempath, idt)
      import MeshModelType, InputParamDefinitionType
      class(MeshModelType), intent(inout) :: this
      character(len=*), intent(in) :: pkgtype
      character(len=*), intent(in) :: pkgname
      character(len=*), intent(in) :: mempath
      type(InputParamDefinitionType), pointer, intent(in) :: idt
    end subroutine
  end interface

  type, abstract, extends(MeshModelType) :: Mesh2dModelType
  contains
    procedure :: create_mesh
  end type Mesh2dModelType

contains

  !> @brief initialize
  !<
  subroutine mesh_init(this, modelname, modeltype, modelfname, nc_fname, &
                       disenum, nctype, lenuni, iout)
    use MemoryManagerExtModule, only: mem_set_value
    class(MeshModelType), intent(inout) :: this
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelfname
    character(len=*), intent(in) :: nc_fname
    integer(I4B), intent(in) :: disenum
    integer(I4B), intent(in) :: nctype
    integer(I4B), intent(in) :: lenuni
    integer(I4B), intent(in) :: iout
    logical(LGP) :: found

    ! initialize base class
    call this%NCModelExportType%init(modelname, modeltype, modelfname, nc_fname, &
                                     disenum, nctype, iout)

    ! allocate and initialize
    allocate (this%chunk_face)
    this%chunk_face = -1

    ! update values from input context
    if (this%ncf_mempath /= '') then
      call mem_set_value(this%chunk_face, 'CHUNK_FACE', this%ncf_mempath, found)
    end if

    if (this%chunk_time > 0 .and. this%chunk_face > 0) then
      this%chunking_active = .true.
    else if (this%chunk_time > 0 .or. this%chunk_face > 0) then
      this%chunk_face = -1
      this%chunk_time = -1
      write (warnmsg, '(a)') 'Ignoring user provided NetCDF chunking parameter. &
        &Define chunk_time and chunk_face input parameters to see an effect in &
        &file "'//trim(nc_fname)//'".'
      call store_warning(warnmsg)
    end if

    if (lenuni == 1) then
      this%lenunits = 'ft'
    else
      this%lenunits = 'm'
    end if

    ! create the netcdf file
    call nf_verify(nf90_create(this%nc_fname, &
                               IOR(NF90_CLOBBER, NF90_NETCDF4), this%ncid), &
                   this%nc_fname)
  end subroutine mesh_init

  !> @brief initialize
  !<
  subroutine mesh_destroy(this)
    use MemoryManagerExtModule, only: mem_set_value
    class(MeshModelType), intent(inout) :: this
    call nf_verify(nf90_close(this%ncid), this%nc_fname)
    deallocate (this%chunk_face)
    nullify (this%chunk_face)
  end subroutine mesh_destroy

  !> @brief define timeseries input variables
  !<
  subroutine df_export(this)
    use NCModelExportModule, only: ExportPackageType
    class(MeshModelType), intent(inout) :: this
    class(ExportPackageType), pointer :: export_pkg
    integer(I4B) :: idx
    do idx = 1, this%pkglist%Count()
      export_pkg => this%get(idx)
      call this%export_df(export_pkg)
    end do
  end subroutine df_export

  !> @brief define export package
  !<
  subroutine export_df(this, export_pkg)
    use NCModelExportModule, only: ExportPackageType
    use DefinitionSelectModule, only: get_param_definition_type
    class(MeshModelType), intent(inout) :: this
    class(ExportPackageType), pointer, intent(in) :: export_pkg
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: iparam, iaux, layer

    ! export defined period input
    do iparam = 1, export_pkg%nparam
      ! initialize
      iaux = 0
      layer = 0
      ! set input definition
      idt => &
        get_param_definition_type(export_pkg%mf6_input%param_dfns, &
                                  export_pkg%mf6_input%component_type, &
                                  export_pkg%mf6_input%subcomponent_type, &
                                  'PERIOD', export_pkg%param_names(iparam), '')

      select case (idt%shape)
      case ('NCPL')
        call this%create_timeseries(idt, iparam, iaux, layer, export_pkg)
      case ('NODES')
        do layer = 1, this%nlay
          call this%create_timeseries(idt, iparam, iaux, layer, export_pkg)
        end do
      case ('NAUX NCPL')
        do iaux = 1, export_pkg%naux
          call this%create_timeseries(idt, iparam, iaux, layer, export_pkg)
        end do
      case ('NAUX NODES')
        do iaux = 1, export_pkg%naux
          do layer = 1, this%nlay
            call this%create_timeseries(idt, iparam, iaux, layer, export_pkg)
          end do
        end do
      case default
      end select
    end do
  end subroutine export_df

  !> @brief create timeseries export variable
  !<
  subroutine create_timeseries(this, idt, iparam, iaux, layer, export_pkg)
    use ConstantsModule, only: DNODATA
    use NCModelExportModule, only: ExportPackageType
    class(MeshModelType), intent(inout) :: this
    type(InputParamDefinitionType), pointer, intent(in) :: idt
    integer(I4B), intent(in) :: iparam
    integer(I4B), intent(in) :: iaux
    integer(I4B), intent(in) :: layer
    class(ExportPackageType), pointer, intent(in) :: export_pkg
    character(len=LINELENGTH) :: varname, longname, nc_tag
    integer(I4B) :: varid

    ! set variable input tag
    nc_tag = this%input_attribute(export_pkg%mf6_input%subcomponent_name, &
                                  idt)

    ! set names
    varname = export_varname(export_pkg%mf6_input%subcomponent_name, &
                             idt%tagname, export_pkg%mf6_input%mempath, &
                             layer=layer, iaux=iaux)
    longname = export_longname(idt%longname, &
                               export_pkg%mf6_input%subcomponent_name, &
                               idt%tagname, export_pkg%mf6_input%mempath, &
                               layer=layer, iaux=iaux)

    ! create the netcdf dependent layer variable
    select case (idt%datatype)
    case ('DOUBLE1D', 'DOUBLE2D')
      call nf_verify(nf90_def_var(this%ncid, varname, NF90_DOUBLE, &
                                  (/this%dim_ids%nmesh_face, &
                                    this%dim_ids%time/), &
                                  varid), &
                     this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, varid, &
                                  '_FillValue', (/DNODATA/)), &
                     this%nc_fname)
    case ('INTEGER1D')
      call nf_verify(nf90_def_var(this%ncid, varname, NF90_INT, &
                                  (/this%dim_ids%nmesh_face, &
                                    this%dim_ids%time/), &
                                  varid), &
                     this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, varid, &
                                  '_FillValue', (/NF90_FILL_INT/)), &
                     this%nc_fname)
    end select

    ! apply chunking parameters
    if (this%chunking_active) then
      call nf_verify(nf90_def_var_chunking(this%ncid, &
                                           varid, &
                                           NF90_CHUNKED, &
                                           (/this%chunk_face, &
                                             this%chunk_time/)), &
                     this%nc_fname)
    end if

    ! deflate and shuffle
    call ncvar_deflate(this%ncid, varid, this%deflate, &
                       this%shuffle, this%nc_fname)

    ! assign variable attributes
    call nf_verify(nf90_put_att(this%ncid, varid, &
                                'units', this%lenunits), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, varid, &
                                'long_name', longname), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, varid, &
                                'mesh', this%mesh_name), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, varid, &
                                'location', 'face'), this%nc_fname)

    ! add grid mapping and mf6 attr
    call ncvar_gridmap(this%ncid, varid, &
                       this%gridmap_name, this%nc_fname)
    call ncvar_mf6attr(this%ncid, varid, layer, iaux, nc_tag, this%nc_fname)

    ! store variable id
    if (idt%tagname == 'AUX') then
      if (layer > 0) then
        export_pkg%varids_aux(iaux, layer) = varid
      else
        export_pkg%varids_aux(iaux, 1) = varid
      end if
    else
      if (layer > 0) then
        export_pkg%varids_param(iparam, layer) = varid
      else
        export_pkg%varids_param(iparam, 1) = varid
      end if
    end if
  end subroutine create_timeseries

  !> @brief create file (group) attributes
  !<
  subroutine add_global_att(this)
    class(MeshModelType), intent(inout) :: this
    ! file scoped title
    call nf_verify(nf90_put_att(this%ncid, NF90_GLOBAL, 'title', &
                                this%annotation%title), this%nc_fname)
    ! source (MODFLOW 6)
    call nf_verify(nf90_put_att(this%ncid, NF90_GLOBAL, 'source', &
                                this%annotation%source), this%nc_fname)
    ! grid type (MODFLOW 6)
    call nf_verify(nf90_put_att(this%ncid, NF90_GLOBAL, 'modflow_grid', &
                                this%annotation%grid), this%nc_fname)
    ! mesh type (MODFLOW 6)
    if (this%annotation%mesh /= '') then
      call nf_verify(nf90_put_att(this%ncid, NF90_GLOBAL, 'mesh', &
                                  this%annotation%mesh), this%nc_fname)

    end if
    ! MODFLOW 6 model type
    call nf_verify(nf90_put_att(this%ncid, NF90_GLOBAL, 'modflow_model', &
                                this%annotation%model), this%nc_fname)
    ! generation datetime
    call nf_verify(nf90_put_att(this%ncid, NF90_GLOBAL, 'history', &
                                this%annotation%history), this%nc_fname)
    ! supported conventions
    call nf_verify(nf90_put_att(this%ncid, NF90_GLOBAL, 'Conventions', &
                                this%annotation%conventions), &
                   this%nc_fname)
  end subroutine add_global_att

  !> @brief write package gridded input data
  !<
  subroutine export_input_arrays(this, pkgtype, pkgname, mempath, param_dfns)
    use MemoryManagerModule, only: get_isize
    class(MeshModelType), intent(inout) :: this
    character(len=*), intent(in) :: pkgtype
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: mempath
    type(InputParamDefinitionType), dimension(:), pointer, &
      intent(in) :: param_dfns
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: iparam, isize
    ! export griddata block parameters
    do iparam = 1, size(param_dfns)
      ! assign param definition pointer
      idt => param_dfns(iparam)
      ! for now only griddata is exported
      if (idt%blockname == 'GRIDDATA') then
        ! veriy variable is allocated
        call get_isize(idt%mf6varname, mempath, isize)
        if (isize > 0) then
          call this%export_input_array(pkgtype, pkgname, mempath, idt)
        end if
      end if
    end do
  end subroutine export_input_arrays

  !> @brief determine packages to write gridded input
  !<
  subroutine add_pkg_data(this)
    use SimVariablesModule, only: idm_context
    use MemoryManagerExtModule, only: mem_set_value
    use MemoryHelperModule, only: create_mem_path
    use IdmDfnSelectorModule, only: param_definitions
    use SourceCommonModule, only: idm_subcomponent_type
    use IdmDfnSelectorModule, only: idm_multi_package
    class(MeshModelType), intent(inout) :: this
    character(LENCOMPONENTNAME) :: ptype, pname, pkgtype
    character(len=LENMEMPATH) :: input_mempath
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgtypes => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgnames => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mempaths => null()
    type(InputParamDefinitionType), dimension(:), pointer :: param_dfns
    character(len=LENMEMPATH) :: mempath
    integer(I4B) :: n
    integer(I4B), pointer :: export_arrays
    logical(LGP) :: found

    input_mempath = create_mem_path(component=this%modelname, context=idm_context)

    ! set pointers to model path package info
    call mem_setptr(pkgtypes, 'PKGTYPES', input_mempath)
    call mem_setptr(pkgnames, 'PKGNAMES', input_mempath)
    call mem_setptr(mempaths, 'MEMPATHS', input_mempath)

    allocate (export_arrays)

    do n = 1, size(mempaths)
      ! initialize export_arrays
      export_arrays = 0

      ! set package attributes
      mempath = mempaths(n)
      pname = pkgnames(n)
      ptype = pkgtypes(n)

      ! export input arrays
      if (mempath /= '') then
        ! update export
        call mem_set_value(export_arrays, 'EXPORT_NC', mempath, found)
        if (export_arrays > 0) then
          pkgtype = idm_subcomponent_type(this%modeltype, ptype)
          param_dfns => param_definitions(this%modeltype, pkgtype)
          call this%export_input_arrays(ptype, pname, mempath, param_dfns)
        end if
      end if
    end do

    ! cleanup
    deallocate (export_arrays)
  end subroutine add_pkg_data

  !> @brief create the model layer dependent variables
  !<
  subroutine define_dependent(this)
    class(MeshModelType), intent(inout) :: this
    character(len=LINELENGTH) :: varname, longname
    integer(I4B) :: k

    ! create a dependent variable for each layer
    do k = 1, this%nlay
      ! initialize names
      varname = ''
      longname = ''

      ! set layer variable and longnames
      write (varname, '(a,i0)') trim(this%xname)//'_l', k
      write (longname, '(a,i0,a)') trim(this%annotation%longname)// &
        ' (layer ', k, ')'

      ! create the netcdf dependent layer variable
      call nf_verify(nf90_def_var(this%ncid, varname, NF90_DOUBLE, &
                                  (/this%dim_ids%nmesh_face, &
                                    this%dim_ids%time/), &
                                  this%var_ids%dependent(k)), &
                     this%nc_fname)

      ! apply chunking parameters
      if (this%chunking_active) then
        call nf_verify(nf90_def_var_chunking(this%ncid, &
                                             this%var_ids%dependent(k), &
                                             NF90_CHUNKED, &
                                             (/this%chunk_face, &
                                               this%chunk_time/)), &
                       this%nc_fname)
      end if

      ! deflate and shuffle
      call ncvar_deflate(this%ncid, this%var_ids%dependent(k), this%deflate, &
                         this%shuffle, this%nc_fname)

      ! assign variable attributes
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%dependent(k), &
                                  'units', this%lenunits), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%dependent(k), &
                                  'standard_name', this%annotation%stdname), &
                     this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%dependent(k), &
                                  'long_name', longname), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%dependent(k), &
                                  '_FillValue', (/DHNOFLO/)), &
                     this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%dependent(k), &
                                  'mesh', this%mesh_name), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%dependent(k), &
                                  'location', 'face'), this%nc_fname)

      ! add grid mapping
      call ncvar_gridmap(this%ncid, this%var_ids%dependent(k), &
                         this%gridmap_name, this%nc_fname)
    end do
  end subroutine define_dependent

  !> @brief create the file grid mapping container variable
  !<
  subroutine define_gridmap(this)
    class(MeshModelType), intent(inout) :: this
    integer(I4B) :: var_id

    ! was projection info provided
    if (this%wkt /= '') then
      ! create projection variable
      call nf_verify(nf90_redef(this%ncid), this%nc_fname)
      call nf_verify(nf90_def_var(this%ncid, this%gridmap_name, NF90_INT, &
                                  var_id), this%nc_fname)
      ! cf-conventions prefers 'crs_wkt'
      !call nf_verify(nf90_put_att(this%ncid, var_id, 'crs_wkt', this%wkt), &
      !               this%nc_fname)
      ! QGIS recognizes 'wkt'
      call nf_verify(nf90_put_att(this%ncid, var_id, 'wkt', this%wkt), &
                     this%nc_fname)
      call nf_verify(nf90_enddef(this%ncid), this%nc_fname)
      call nf_verify(nf90_put_var(this%ncid, var_id, 1), &
                     this%nc_fname)
    end if
  end subroutine define_gridmap

  !> @brief create the file mesh container variable
  !<
  subroutine create_mesh(this)
    class(Mesh2dModelType), intent(inout) :: this

    ! create mesh container variable
    call nf_verify(nf90_def_var(this%ncid, this%mesh_name, NF90_INT, &
                                this%var_ids%mesh), this%nc_fname)

    ! assign container variable attributes
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh, 'cf_role', &
                                'mesh_topology'), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh, 'long_name', &
                                '2D mesh topology'), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh, &
                                'topology_dimension', 2), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh, 'face_dimension', &
                                'nmesh_face'), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh, &
                                'node_coordinates', 'mesh_node_x mesh_node_y'), &
                   this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh, &
                                'face_coordinates', 'mesh_face_x mesh_face_y'), &
                   this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh, &
                                'face_node_connectivity', 'mesh_face_nodes'), &
                   this%nc_fname)

    ! create mesh x node (mesh vertex) variable
    call nf_verify(nf90_def_var(this%ncid, 'mesh_node_x', NF90_DOUBLE, &
                                (/this%dim_ids%nmesh_node/), &
                                this%var_ids%mesh_node_x), this%nc_fname)

    ! assign mesh x node variable attributes
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_node_x, &
                                'units', this%lenunits), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_node_x, &
                                'standard_name', 'projection_x_coordinate'), &
                   this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_node_x, &
                                'long_name', 'Easting'), this%nc_fname)

    if (this%wkt /= '') then
      ! associate with projection
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_node_x, &
                                  'grid_mapping', this%gridmap_name), &
                     this%nc_fname)
    end if

    ! create mesh y node (mesh vertex) variable
    call nf_verify(nf90_def_var(this%ncid, 'mesh_node_y', NF90_DOUBLE, &
                                (/this%dim_ids%nmesh_node/), &
                                this%var_ids%mesh_node_y), this%nc_fname)

    ! assign mesh y variable attributes
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_node_y, &
                                'units', this%lenunits), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_node_y, &
                                'standard_name', 'projection_y_coordinate'), &
                   this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_node_y, &
                                'long_name', 'Northing'), this%nc_fname)

    if (this%wkt /= '') then
      ! associate with projection
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_node_y, &
                                  'grid_mapping', this%gridmap_name), &
                     this%nc_fname)
    end if

    ! create mesh x face (cell vertex) variable
    call nf_verify(nf90_def_var(this%ncid, 'mesh_face_x', NF90_DOUBLE, &
                                (/this%dim_ids%nmesh_face/), &
                                this%var_ids%mesh_face_x), this%nc_fname)

    ! assign mesh x face variable attributes
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_x, &
                                'units', this%lenunits), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_x, &
                                'standard_name', 'projection_x_coordinate'), &
                   this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_x, &
                                'long_name', 'Easting'), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_x, 'bounds', &
                                'mesh_face_xbnds'), this%nc_fname)
    if (this%wkt /= '') then
      ! associate with projection
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_x, &
                                  'grid_mapping', this%gridmap_name), &
                     this%nc_fname)
    end if

    ! create mesh x cell bounds variable
    call nf_verify(nf90_def_var(this%ncid, 'mesh_face_xbnds', NF90_DOUBLE, &
                                (/this%dim_ids%max_nmesh_face_nodes, &
                                  this%dim_ids%nmesh_face/), &
                                this%var_ids%mesh_face_xbnds), &
                   this%nc_fname)

    ! create mesh y face (cell vertex) variable
    call nf_verify(nf90_def_var(this%ncid, 'mesh_face_y', NF90_DOUBLE, &
                                (/this%dim_ids%nmesh_face/), &
                                this%var_ids%mesh_face_y), this%nc_fname)

    ! assign mesh y face variable attributes
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_y, &
                                'units', this%lenunits), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_y, &
                                'standard_name', 'projection_y_coordinate'), &
                   this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_y, &
                                'long_name', 'Northing'), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_y, 'bounds', &
                                'mesh_face_ybnds'), this%nc_fname)

    if (this%wkt /= '') then
      ! associate with projection
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_y, &
                                  'grid_mapping', this%gridmap_name), &
                     this%nc_fname)
    end if

    ! create mesh y cell bounds variable
    call nf_verify(nf90_def_var(this%ncid, 'mesh_face_ybnds', NF90_DOUBLE, &
                                (/this%dim_ids%max_nmesh_face_nodes, &
                                  this%dim_ids%nmesh_face/), &
                                this%var_ids%mesh_face_ybnds), &
                   this%nc_fname)

    ! create mesh face nodes variable
    call nf_verify(nf90_def_var(this%ncid, 'mesh_face_nodes', NF90_INT, &
                                (/this%dim_ids%max_nmesh_face_nodes, &
                                  this%dim_ids%nmesh_face/), &
                                this%var_ids%mesh_face_nodes), &
                   this%nc_fname)

    ! assign variable attributes
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_nodes, &
                                'cf_role', 'face_node_connectivity'), &
                   this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_nodes, &
                                'long_name', &
                                'Vertices bounding cell (counterclockwise)'), &
                   this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_nodes, &
                                '_FillValue', (/NF90_FILL_INT/)), &
                   this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%mesh_face_nodes, &
                                'start_index', 1), this%nc_fname)
  end subroutine create_mesh

  !> @brief define variable chunking
  !<
  subroutine ncvar_chunk(ncid, varid, chunk_face, nc_fname)
    integer(I4B), intent(in) :: ncid
    integer(I4B), intent(in) :: varid
    integer(I4B), intent(in) :: chunk_face
    character(len=*), intent(in) :: nc_fname
    if (chunk_face > 0) then
      call nf_verify(nf90_def_var_chunking(ncid, varid, NF90_CHUNKED, &
                                           (/chunk_face/)), nc_fname)
    end if
  end subroutine ncvar_chunk

  !> @brief define variable compression
  !<
  subroutine ncvar_deflate(ncid, varid, deflate, shuffle, nc_fname)
    integer(I4B), intent(in) :: ncid
    integer(I4B), intent(in) :: varid
    integer(I4B), intent(in) :: deflate
    integer(I4B), intent(in) :: shuffle
    character(len=*), intent(in) :: nc_fname
    if (deflate >= 0) then
      call nf_verify(nf90_def_var_deflate(ncid, varid, shuffle=shuffle, &
                                          deflate=1, deflate_level=deflate), &
                     nc_fname)
    end if
  end subroutine ncvar_deflate

  !> @brief put variable gridmap attributes
  !<
  subroutine ncvar_gridmap(ncid, varid, gridmap_name, nc_fname)
    integer(I4B), intent(in) :: ncid
    integer(I4B), intent(in) :: varid
    character(len=*), intent(in) :: gridmap_name
    character(len=*), intent(in) :: nc_fname
    if (gridmap_name /= '') then
      call nf_verify(nf90_put_att(ncid, varid, 'coordinates', &
                                  'mesh_face_x mesh_face_y'), nc_fname)
      call nf_verify(nf90_put_att(ncid, varid, 'grid_mapping', &
                                  gridmap_name), nc_fname)
    end if
  end subroutine ncvar_gridmap

  !> @brief put variable internal attributes
  !<
  subroutine ncvar_mf6attr(ncid, varid, layer, iaux, nc_tag, nc_fname)
    integer(I4B), intent(in) :: ncid
    integer(I4B), intent(in) :: varid
    integer(I4B), intent(in) :: layer
    integer(I4B), intent(in) :: iaux
    character(len=*), intent(in) :: nc_tag
    character(len=*), intent(in) :: nc_fname
    if (nc_tag /= '') then
      call nf_verify(nf90_put_att(ncid, varid, 'modflow_input', &
                                  nc_tag), nc_fname)
      if (layer > 0) then
        call nf_verify(nf90_put_att(ncid, varid, 'layer', &
                                    layer), nc_fname)
      end if
      if (iaux > 0) then
        call nf_verify(nf90_put_att(ncid, varid, 'modflow_iaux', &
                                    iaux), nc_fname)
      end if
    end if
  end subroutine ncvar_mf6attr

end module MeshModelModule
