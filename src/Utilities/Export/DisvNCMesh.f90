!> @brief This module contains the MeshDisvModelModule
!!
!! This module defines UGRID layered mesh compliant netcdf
!! export type for DISV models. It is dependent on netcdf
!! libraries.
!!
!<
module MeshDisvModelModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENBIGLINE, LENCOMPONENTNAME, &
                             LENMEMPATH
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use MemoryManagerModule, only: mem_setptr
  use InputDefinitionModule, only: InputParamDefinitionType
  use CharacterStringModule, only: CharacterStringType
  use MeshModelModule, only: Mesh2dModelType, MeshNCDimIdType, MeshNCVarIdType, &
                             ncvar_chunk, ncvar_deflate, ncvar_gridmap, &
                             ncvar_mf6attr, export_varname
  use NCModelExportModule, only: export_longname
  use DisvModule, only: DisvType
  use NetCDFCommonModule, only: nf_verify
  use netcdf

  implicit none
  private
  public :: Mesh2dDisvExportType

  ! UGRID layered mesh DISV
  type, extends(Mesh2dModelType) :: Mesh2dDisvExportType
    type(DisvType), pointer :: disv => null() !< pointer to model disv package
  contains
    procedure :: init => disv_export_init
    procedure :: destroy => disv_export_destroy
    procedure :: df
    procedure :: step
    procedure :: export_input_array
    procedure :: package_step_ilayer
    procedure :: package_step
    procedure :: export_layer_2d
    procedure :: define_dim
    procedure :: add_mesh_data
  end type Mesh2dDisvExportType

contains

  !> @brief netcdf export disv init
  !<
  subroutine disv_export_init(this, modelname, modeltype, modelfname, nc_fname, &
                              disenum, nctype, iout)
    use ArrayHandlersModule, only: expandarray
    class(Mesh2dDisvExportType), intent(inout) :: this
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelfname
    character(len=*), intent(in) :: nc_fname
    integer(I4B), intent(in) :: disenum
    integer(I4B), intent(in) :: nctype
    integer(I4B), intent(in) :: iout

    ! set nlay
    this%nlay = this%disv%nlay

    ! allocate var_id arrays
    allocate (this%var_ids%dependent(this%nlay))

    ! initialize base class
    call this%mesh_init(modelname, modeltype, modelfname, nc_fname, disenum, &
                        nctype, this%disv%lenuni, iout)
  end subroutine disv_export_init

  !> @brief netcdf export disv destroy
  !<
  subroutine disv_export_destroy(this)
    class(Mesh2dDisvExportType), intent(inout) :: this
    deallocate (this%var_ids%dependent)
    ! destroy base class
    call this%mesh_destroy()
    call this%NCModelExportType%destroy()
  end subroutine disv_export_destroy

  !> @brief netcdf export define
  !<
  subroutine df(this)
    use ConstantsModule, only: MVALIDATE
    use SimVariablesModule, only: isim_mode
    class(Mesh2dDisvExportType), intent(inout) :: this
    ! put root group file scope attributes
    call this%add_global_att()
    ! define root group dimensions and coordinate variables
    call this%define_dim()
    ! define mesh variables
    call this%create_mesh()
    if (isim_mode /= MVALIDATE) then
      ! define the dependent variable
      call this%define_dependent()
    end if
    ! exit define mode
    call nf_verify(nf90_enddef(this%ncid), this%nc_fname)
    ! create mesh
    call this%add_mesh_data()
    ! define and set package input griddata
    call this%add_pkg_data()
    ! define and set gridmap variable
    call this%define_gridmap()
    ! synchronize file
    call nf_verify(nf90_sync(this%ncid), this%nc_fname)
  end subroutine df

  !> @brief netcdf export step
  !<
  subroutine step(this)
    use ConstantsModule, only: DHNOFLO
    use TdisModule, only: totim
    class(Mesh2dDisvExportType), intent(inout) :: this
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    integer(I4B) :: n, k, nvals
    integer(I4B), dimension(2) :: dis_shape
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d

    ! initialize
    nullify (dbl1d)
    nullify (dbl2d)

    ! increment step
    this%stepcnt = this%stepcnt + 1

    dis_shape(1) = this%disv%ncpl
    dis_shape(2) = this%disv%nlay

    nvals = product(dis_shape)

    ! add data to dependent variable
    if (size(this%disv%nodeuser) < &
        size(this%disv%nodereduced)) then
      ! allocate nodereduced size 1d array
      allocate (dbl1d(size(this%disv%nodereduced)))

      ! initialize DHNOFLO for non-active cells
      dbl1d = DHNOFLO

      ! update active cells
      do n = 1, size(this%disv%nodereduced)
        if (this%disv%nodereduced(n) > 0) then
          dbl1d(n) = this%x(this%disv%nodereduced(n))
        end if
      end do

      dbl2d(1:dis_shape(1), 1:dis_shape(2)) => dbl1d(1:nvals)
    else
      dbl2d(1:dis_shape(1), 1:dis_shape(2)) => this%x(1:nvals)
    end if

    do k = 1, this%disv%nlay
      ! extend array with step data
      call nf_verify(nf90_put_var(this%ncid, &
                                  this%var_ids%dependent(k), dbl2d(:, k), &
                                  start=(/1, this%stepcnt/), &
                                  count=(/this%disv%ncpl, 1/)), &
                     this%nc_fname)
    end do

    ! write to time coordinate variable
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%time, &
                                totim, start=(/this%stepcnt/)), &
                   this%nc_fname)

    ! update file
    call nf_verify(nf90_sync(this%ncid), this%nc_fname)

    ! cleanup
    if (associated(dbl1d)) deallocate (dbl1d)
    nullify (dbl1d)
    nullify (dbl2d)
  end subroutine step

  !> @brief netcdf export package dynamic input with ilayer index variable
  !<
  subroutine package_step_ilayer(this, export_pkg, ilayer_varname, ilayer)
    use ConstantsModule, only: DNODATA, DZERO
    use TdisModule, only: kper
    use DefinitionSelectModule, only: get_param_definition_type
    use NCModelExportModule, only: ExportPackageType
    class(Mesh2dDisvExportType), intent(inout) :: this
    class(ExportPackageType), pointer, intent(in) :: export_pkg
    character(len=*), intent(in) :: ilayer_varname
    integer(I4B), intent(in) :: ilayer
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B), dimension(:), pointer, contiguous :: int1d
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B), dimension(:), pointer, contiguous :: ialayer
    real(DP), dimension(:), contiguous, pointer :: dbl1d_ptr
    character(len=LINELENGTH) :: nc_varname, input_attr
    integer(I4B) :: n, iparam, nvals
    logical(LGP) :: ilayer_read

    ! initialize
    nullify (ialayer)
    ilayer_read = .false.

    ! set pointer to ilayer variable
    call mem_setptr(ialayer, export_pkg%param_names(ilayer), &
                    export_pkg%mf6_input%mempath)

    ! check if layer index variable was read
    if (export_pkg%param_reads(ilayer)%invar == 1) then
      ilayer_read = .true.
    end if

    ! export defined period input
    do iparam = 1, export_pkg%nparam
      ! check if variable was read this period
      if (export_pkg%param_reads(iparam)%invar < 1) cycle

      ! set input definition
      idt => &
        get_param_definition_type(export_pkg%mf6_input%param_dfns, &
                                  export_pkg%mf6_input%component_type, &
                                  export_pkg%mf6_input%subcomponent_type, &
                                  'PERIOD', export_pkg%param_names(iparam), '')

      ! set variable name and input string
      nc_varname = trim(export_pkg%mf6_input%subcomponent_name)//'_'// &
                   trim(idt%mf6varname)
      input_attr = this%input_attribute(export_pkg%mf6_input%subcomponent_name, &
                                        idt)

      ! export arrays
      select case (idt%datatype)
      case ('INTEGER1D')
        call mem_setptr(int1d, idt%mf6varname, export_pkg%mf6_input%mempath)
        call nc_export_int1d(this%ncid, this%dim_ids, this%var_ids, this%disv, &
                             int1d, nc_varname, &
                             export_pkg%mf6_input%subcomponent_name, &
                             idt%tagname, this%gridmap_name, idt%shape, &
                             idt%longname, input_attr, this%deflate, &
                             this%shuffle, this%chunk_face, kper, this%nc_fname)
      case ('DOUBLE1D')
        call mem_setptr(dbl1d, idt%mf6varname, export_pkg%mf6_input%mempath)
        call this%export_layer_2d(export_pkg, idt, ilayer_read, ialayer, &
                                  dbl1d, nc_varname, input_attr)
      case ('DOUBLE2D')
        call mem_setptr(dbl2d, idt%mf6varname, export_pkg%mf6_input%mempath)
        nvals = this%disv%ncpl
        do n = 1, size(dbl2d, dim=1) !naux
          dbl1d_ptr(1:nvals) => dbl2d(n, :)
          if (all(dbl1d_ptr == DZERO)) then
          else
            call this%export_layer_2d(export_pkg, idt, ilayer_read, ialayer, &
                                      dbl1d_ptr, nc_varname, input_attr, n)
          end if
        end do
      case default
        errmsg = 'EXPORT ilayer unsupported datatype='//trim(idt%datatype)
        call store_error(errmsg, .true.)
      end select
    end do

    ! synchronize file
    call nf_verify(nf90_sync(this%ncid), this%nc_fname)
  end subroutine package_step_ilayer

  !> @brief netcdf export package dynamic input
  !<
  subroutine package_step(this, export_pkg)
    use NCModelExportModule, only: ExportPackageType
    class(Mesh2dDisvExportType), intent(inout) :: this
    class(ExportPackageType), pointer, intent(in) :: export_pkg
    errmsg = 'NetCDF period export not supported for model='// &
             trim(this%modelname)//', package='// &
             trim(export_pkg%mf6_input%subcomponent_name)
    call store_error(errmsg, .true.)
    ! synchronize file
    call nf_verify(nf90_sync(this%ncid), this%nc_fname)
  end subroutine package_step

  !> @brief export layer variable as full grid
  !<
  subroutine export_layer_2d(this, export_pkg, idt, ilayer_read, ialayer, &
                             dbl1d, nc_varname, input_attr, iaux)
    use ConstantsModule, only: DNODATA, DZERO
    use NCModelExportModule, only: ExportPackageType
    class(Mesh2dDisvExportType), intent(inout) :: this
    class(ExportPackageType), pointer, intent(in) :: export_pkg
    type(InputParamDefinitionType), pointer, intent(in) :: idt
    logical(LGP), intent(in) :: ilayer_read
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: ialayer
    real(DP), dimension(:), pointer, contiguous, intent(in) :: dbl1d
    character(len=*), intent(in) :: nc_varname
    character(len=*), intent(in) :: input_attr
    integer(I4B), optional, intent(in) :: iaux
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B) :: n, j, k, idxaux

    ! initialize
    idxaux = 0
    if (present(iaux)) then
      idxaux = iaux
    end if

    allocate (dbl2d(export_pkg%mshape(2), export_pkg%mshape(1)))

    if (ilayer_read) then
      do k = 1, size(dbl2d, dim=2)
        n = 0
        do j = 1, size(dbl2d, dim=1)
          n = n + 1
          if (ialayer(n) == k) then
            dbl2d(j, k) = dbl1d(n)
          else
            dbl2d(j, k) = DNODATA
          end if
        end do
      end do
    else
      dbl2d = DNODATA
      dbl2d(:, 1) = dbl1d(:)
    end if

    call nc_export_dbl2d(this%ncid, this%dim_ids, this%var_ids, this%disv, &
                         dbl2d, nc_varname, &
                         export_pkg%mf6_input%subcomponent_name, idt%tagname, &
                         this%gridmap_name, idt%shape, idt%longname, input_attr, &
                         this%deflate, this%shuffle, this%chunk_face, &
                         export_pkg%iper, idxaux, this%nc_fname)

    deallocate (dbl2d)
  end subroutine export_layer_2d

  !> @brief netcdf export an input array
  !<
  subroutine export_input_array(this, pkgtype, pkgname, mempath, idt)
    class(Mesh2dDisvExportType), intent(inout) :: this
    character(len=*), intent(in) :: pkgtype
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: mempath
    type(InputParamDefinitionType), pointer, intent(in) :: idt
    integer(I4B), dimension(:), pointer, contiguous :: int1d
    integer(I4B), dimension(:, :), pointer, contiguous :: int2d
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    character(len=LINELENGTH) :: nc_varname, input_attr
    integer(I4B) :: iper, iaux

    iper = 0
    iaux = 0

    ! set package base name
    nc_varname = trim(pkgname)//'_'//trim(idt%mf6varname)
    ! put input attributes
    input_attr = this%input_attribute(pkgname, idt)

    select case (idt%datatype)
    case ('INTEGER1D')
      call mem_setptr(int1d, idt%mf6varname, mempath)
      call nc_export_int1d(this%ncid, this%dim_ids, this%var_ids, this%disv, &
                           int1d, nc_varname, pkgname, idt%tagname, &
                           this%gridmap_name, idt%shape, idt%longname, &
                           input_attr, this%deflate, this%shuffle, &
                           this%chunk_face, iper, this%nc_fname)
    case ('INTEGER2D')
      call mem_setptr(int2d, idt%mf6varname, mempath)
      call nc_export_int2d(this%ncid, this%dim_ids, this%var_ids, this%disv, &
                           int2d, nc_varname, pkgname, idt%tagname, &
                           this%gridmap_name, idt%shape, idt%longname, &
                           input_attr, this%deflate, this%shuffle, &
                           this%chunk_face, this%nc_fname)
    case ('DOUBLE1D')
      call mem_setptr(dbl1d, idt%mf6varname, mempath)
      call nc_export_dbl1d(this%ncid, this%dim_ids, this%var_ids, this%disv, &
                           dbl1d, nc_varname, pkgname, idt%tagname, &
                           this%gridmap_name, idt%shape, idt%longname, &
                           input_attr, this%deflate, this%shuffle, &
                           this%chunk_face, this%nc_fname)
    case ('DOUBLE2D')
      call mem_setptr(dbl2d, idt%mf6varname, mempath)
      call nc_export_dbl2d(this%ncid, this%dim_ids, this%var_ids, this%disv, &
                           dbl2d, nc_varname, pkgname, idt%tagname, &
                           this%gridmap_name, idt%shape, idt%longname, &
                           input_attr, this%deflate, this%shuffle, &
                           this%chunk_face, iper, iaux, this%nc_fname)
    case default
      ! no-op, no other datatypes exported
    end select
  end subroutine export_input_array

  !> @brief netcdf export define dimensions
  !<
  subroutine define_dim(this)
    use ConstantsModule, only: MVALIDATE
    use SimVariablesModule, only: isim_mode
    class(Mesh2dDisvExportType), intent(inout) :: this
    integer(I4B), dimension(:), contiguous, pointer :: ncvert

    ! set pointers to input context
    call mem_setptr(ncvert, 'NCVERT', this%dis_mempath)

    ! time
    if (isim_mode /= MVALIDATE) then
      call nf_verify(nf90_def_dim(this%ncid, 'time', this%totnstp, &
                                  this%dim_ids%time), this%nc_fname)
      call nf_verify(nf90_def_var(this%ncid, 'time', NF90_DOUBLE, &
                                  this%dim_ids%time, this%var_ids%time), &
                     this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%time, 'calendar', &
                                  'standard'), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%time, 'units', &
                                  this%datetime), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%time, 'axis', 'T'), &
                     this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%time, 'standard_name', &
                                  'time'), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%time, 'long_name', &
                                  'time'), this%nc_fname)
    end if

    ! mesh
    call nf_verify(nf90_def_dim(this%ncid, 'nmesh_node', this%disv%nvert, &
                                this%dim_ids%nmesh_node), this%nc_fname)
    call nf_verify(nf90_def_dim(this%ncid, 'nmesh_face', this%disv%ncpl, &
                                this%dim_ids%nmesh_face), this%nc_fname)
    call nf_verify(nf90_def_dim(this%ncid, 'max_nmesh_face_nodes', &
                                maxval(ncvert), &
                                this%dim_ids%max_nmesh_face_nodes), &
                   this%nc_fname)

    ! ncpl, nlay
    call nf_verify(nf90_def_dim(this%ncid, 'nlay', this%disv%nlay, &
                                this%dim_ids%nlay), this%nc_fname)
  end subroutine define_dim

  !> @brief netcdf export add mesh information
  !<
  subroutine add_mesh_data(this)
    use BaseDisModule, only: dis_transform_xy
    class(Mesh2dDisvExportType), intent(inout) :: this
    integer(I4B), dimension(:), contiguous, pointer :: icell2d => null()
    integer(I4B), dimension(:), contiguous, pointer :: ncvert => null()
    integer(I4B), dimension(:), contiguous, pointer :: icvert => null()
    real(DP), dimension(:), contiguous, pointer :: cell_x => null()
    real(DP), dimension(:), contiguous, pointer :: cell_y => null()
    real(DP), dimension(:), contiguous, pointer :: vert_x => null()
    real(DP), dimension(:), contiguous, pointer :: vert_y => null()
    real(DP), dimension(:), contiguous, pointer :: cell_xt => null()
    real(DP), dimension(:), contiguous, pointer :: cell_yt => null()
    real(DP), dimension(:), contiguous, pointer :: vert_xt => null()
    real(DP), dimension(:), contiguous, pointer :: vert_yt => null()
    real(DP) :: x_transform, y_transform
    integer(I4B) :: n, m, idx, cnt, iv, maxvert
    integer(I4B), dimension(:), allocatable :: verts
    real(DP), dimension(:), allocatable :: bnds
    integer(I4B) :: istop

    ! set pointers to input context
    call mem_setptr(icell2d, 'ICELL2D', this%dis_mempath)
    call mem_setptr(ncvert, 'NCVERT', this%dis_mempath)
    call mem_setptr(icvert, 'ICVERT', this%dis_mempath)
    call mem_setptr(cell_x, 'XC', this%dis_mempath)
    call mem_setptr(cell_y, 'YC', this%dis_mempath)
    call mem_setptr(vert_x, 'XV', this%dis_mempath)
    call mem_setptr(vert_y, 'YV', this%dis_mempath)

    ! allocate x, y transform arrays
    allocate (cell_xt(size(cell_x)))
    allocate (cell_yt(size(cell_y)))
    allocate (vert_xt(size(vert_x)))
    allocate (vert_yt(size(vert_y)))

    ! set mesh container variable value to 1
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh, 1), &
                   this%nc_fname)

    ! transform vert x and y
    do n = 1, size(vert_x)
      call dis_transform_xy(vert_x(n), vert_y(n), &
                            this%disv%xorigin, &
                            this%disv%yorigin, &
                            this%disv%angrot, &
                            x_transform, y_transform)
      vert_xt(n) = x_transform
      vert_yt(n) = y_transform
    end do

    ! transform cell x and y
    do n = 1, size(cell_x)
      call dis_transform_xy(cell_x(n), cell_y(n), &
                            this%disv%xorigin, &
                            this%disv%yorigin, &
                            this%disv%angrot, &
                            x_transform, y_transform)
      cell_xt(n) = x_transform
      cell_yt(n) = y_transform
    end do

    ! write node_x and node_y arrays to netcdf file
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_node_x, &
                                vert_xt), this%nc_fname)
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_node_y, &
                                vert_yt), this%nc_fname)

    ! write face_x and face_y arrays to netcdf file
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_face_x, &
                                cell_xt), this%nc_fname)
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_face_y, &
                                cell_yt), this%nc_fname)

    ! initialize max vertices required to define cell
    maxvert = maxval(ncvert)

    ! allocate temporary arrays
    allocate (verts(maxvert))
    allocate (bnds(maxvert))

    ! set face nodes array
    cnt = 0
    do n = 1, size(ncvert)
      verts = NF90_FILL_INT
      idx = cnt + ncvert(n)
      iv = 0
      istop = cnt + 1
      do m = idx, istop, -1
        cnt = cnt + 1
        iv = iv + 1
        verts(iv) = icvert(m)
      end do

      ! write face nodes array to netcdf file
      call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_face_nodes, &
                                  verts, start=(/1, n/), &
                                  count=(/maxvert, 1/)), &
                     this%nc_fname)

      ! set face y bounds array
      bnds = NF90_FILL_DOUBLE
      do m = 1, size(bnds)
        if (verts(m) /= NF90_FILL_INT) then
          bnds(m) = vert_yt(verts(m))
        end if
        ! write face y bounds array to netcdf file
        call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_face_ybnds, &
                                    bnds, start=(/1, n/), &
                                    count=(/maxvert, 1/)), &
                       this%nc_fname)
      end do

      ! set face x bounds array
      bnds = NF90_FILL_DOUBLE
      do m = 1, size(bnds)
        if (verts(m) /= NF90_FILL_INT) then
          bnds(m) = vert_xt(verts(m))
        end if
        ! write face x bounds array to netcdf file
        call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_face_xbnds, &
                                    bnds, start=(/1, n/), &
                                    count=(/maxvert, 1/)), &
                       this%nc_fname)
      end do
    end do

    ! cleanup
    deallocate (bnds)
    deallocate (verts)
    deallocate (cell_xt)
    deallocate (cell_yt)
    deallocate (vert_xt)
    deallocate (vert_yt)
  end subroutine add_mesh_data

  !> @brief netcdf export 1D integer array
  !<
  subroutine nc_export_int1d(ncid, dim_ids, var_ids, dis, p_mem, nc_varname, &
                             pkgname, tagname, gridmap_name, shapestr, longname, &
                             nc_tag, deflate, shuffle, chunk_face, iper, nc_fname)
    integer(I4B), intent(in) :: ncid
    type(MeshNCDimIdType), intent(inout) :: dim_ids
    type(MeshNCVarIdType), intent(inout) :: var_ids
    type(DisvType), pointer, intent(in) :: dis
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: p_mem
    character(len=*), intent(in) :: nc_varname
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: gridmap_name
    character(len=*), intent(in) :: shapestr
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: nc_tag
    integer(I4B), intent(in) :: deflate
    integer(I4B), intent(in) :: shuffle
    integer(I4B), intent(in) :: chunk_face
    integer(I4B), intent(in) :: iper
    character(len=*), intent(in) :: nc_fname
    integer(I4B), dimension(2) :: dis_shape
    integer(I4B), dimension(:, :), pointer, contiguous :: int2d
    integer(I4B) :: axis_sz, nvals, k
    integer(I4B), dimension(:), allocatable :: var_id
    character(len=LINELENGTH) :: longname_l, varname_l

    if (shapestr == 'NCPL') then
      ! set names
      varname_l = export_varname(nc_varname)
      longname_l = export_longname(longname, pkgname, tagname, layer=0, iper=iper)

      allocate (var_id(1))
      axis_sz = dim_ids%nmesh_face

      ! reenter define mode and create variable
      call nf_verify(nf90_redef(ncid), nc_fname)
      call nf_verify(nf90_def_var(ncid, varname_l, NF90_INT, &
                                  (/axis_sz/), var_id(1)), &
                     nc_fname)

      ! apply chunking parameters
      call ncvar_chunk(ncid, var_id(1), chunk_face, nc_fname)
      ! deflate and shuffle
      call ncvar_deflate(ncid, var_id(1), deflate, shuffle, nc_fname)

      ! put attr
      call nf_verify(nf90_put_att(ncid, var_id(1), '_FillValue', &
                                  (/NF90_FILL_INT/)), nc_fname)
      call nf_verify(nf90_put_att(ncid, var_id(1), 'long_name', &
                                  longname_l), nc_fname)

      ! add grid mapping and mf6 attr
      call ncvar_gridmap(ncid, var_id(1), gridmap_name, nc_fname)
      call ncvar_mf6attr(ncid, var_id(1), 0, iper, 0, nc_tag, nc_fname)

      ! exit define mode and write data
      call nf_verify(nf90_enddef(ncid), nc_fname)
      call nf_verify(nf90_put_var(ncid, var_id(1), p_mem), &
                     nc_fname)

    else
      allocate (var_id(dis%nlay))

      ! reenter define mode and create variable
      call nf_verify(nf90_redef(ncid), nc_fname)
      do k = 1, dis%nlay
        ! set names
        varname_l = export_varname(nc_varname, layer=k, iper=iper)
        longname_l = export_longname(longname, pkgname, tagname, layer=k, &
                                     iper=iper)

        call nf_verify(nf90_def_var(ncid, varname_l, NF90_INT, &
                                    (/dim_ids%nmesh_face/), var_id(k)), &
                       nc_fname)

        ! apply chunking parameters
        call ncvar_chunk(ncid, var_id(k), chunk_face, nc_fname)
        ! defalte and shuffle
        call ncvar_deflate(ncid, var_id(k), deflate, shuffle, nc_fname)

        ! put attr
        call nf_verify(nf90_put_att(ncid, var_id(k), '_FillValue', &
                                    (/NF90_FILL_INT/)), nc_fname)
        call nf_verify(nf90_put_att(ncid, var_id(k), 'long_name', &
                                    longname_l), nc_fname)

        ! add grid mapping and mf6 attr
        call ncvar_gridmap(ncid, var_id(k), gridmap_name, nc_fname)
        call ncvar_mf6attr(ncid, var_id(k), k, iper, 0, nc_tag, nc_fname)
      end do

      ! reshape input
      dis_shape(1) = dis%ncpl
      dis_shape(2) = dis%nlay
      nvals = product(dis_shape)
      int2d(1:dis_shape(1), 1:dis_shape(2)) => p_mem(1:nvals)

      ! exit define mode and write data
      call nf_verify(nf90_enddef(ncid), nc_fname)
      do k = 1, dis%nlay
        call nf_verify(nf90_put_var(ncid, var_id(k), int2d(:, k)), nc_fname)
      end do

      ! cleanup
      deallocate (var_id)
    end if
  end subroutine nc_export_int1d

  !> @brief netcdf export 2D integer array
  !<
  subroutine nc_export_int2d(ncid, dim_ids, var_ids, disv, p_mem, nc_varname, &
                             pkgname, tagname, gridmap_name, shapestr, longname, &
                             nc_tag, deflate, shuffle, chunk_face, nc_fname)
    integer(I4B), intent(in) :: ncid
    type(MeshNCDimIdType), intent(inout) :: dim_ids
    type(MeshNCVarIdType), intent(inout) :: var_ids
    type(DisvType), pointer, intent(in) :: disv
    integer(I4B), dimension(:, :), pointer, contiguous, intent(in) :: p_mem
    character(len=*), intent(in) :: nc_varname
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: gridmap_name
    character(len=*), intent(in) :: shapestr
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: nc_tag
    integer(I4B), intent(in) :: deflate
    integer(I4B), intent(in) :: shuffle
    integer(I4B), intent(in) :: chunk_face
    character(len=*), intent(in) :: nc_fname
    integer(I4B), dimension(:), allocatable :: var_id
    character(len=LINELENGTH) :: longname_l, varname_l
    integer(I4B) :: k

    allocate (var_id(disv%nlay))

    ! reenter define mode and create variable
    call nf_verify(nf90_redef(ncid), nc_fname)
    do k = 1, disv%nlay
      ! set names
      varname_l = export_varname(nc_varname, layer=k)
      longname_l = export_longname(longname, pkgname, tagname, k)

      call nf_verify(nf90_def_var(ncid, varname_l, NF90_INT, &
                                  (/dim_ids%nmesh_face/), var_id(k)), &
                     nc_fname)

      ! apply chunking parameters
      call ncvar_chunk(ncid, var_id(k), chunk_face, nc_fname)
      ! deflate and shuffle
      call ncvar_deflate(ncid, var_id(k), deflate, shuffle, nc_fname)

      ! put attr
      call nf_verify(nf90_put_att(ncid, var_id(k), '_FillValue', &
                                  (/NF90_FILL_INT/)), nc_fname)
      call nf_verify(nf90_put_att(ncid, var_id(k), 'long_name', &
                                  longname_l), nc_fname)

      ! add grid mapping and mf6 attr
      call ncvar_gridmap(ncid, var_id(k), gridmap_name, nc_fname)
      call ncvar_mf6attr(ncid, var_id(k), k, 0, 0, nc_tag, nc_fname)
    end do

    ! exit define mode and write data
    call nf_verify(nf90_enddef(ncid), nc_fname)
    do k = 1, disv%nlay
      call nf_verify(nf90_put_var(ncid, var_id(k), p_mem(:, k)), nc_fname)
    end do

    deallocate (var_id)
  end subroutine nc_export_int2d

  !> @brief netcdf export 1D double array
  !<
  subroutine nc_export_dbl1d(ncid, dim_ids, var_ids, dis, p_mem, nc_varname, &
                             pkgname, tagname, gridmap_name, shapestr, longname, &
                             nc_tag, deflate, shuffle, chunk_face, nc_fname)
    integer(I4B), intent(in) :: ncid
    type(MeshNCDimIdType), intent(inout) :: dim_ids
    type(MeshNCVarIdType), intent(inout) :: var_ids
    type(DisvType), pointer, intent(in) :: dis
    real(DP), dimension(:), pointer, contiguous, intent(in) :: p_mem
    character(len=*), intent(in) :: nc_varname
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: gridmap_name
    character(len=*), intent(in) :: shapestr
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: nc_tag
    integer(I4B), intent(in) :: deflate
    integer(I4B), intent(in) :: shuffle
    integer(I4B), intent(in) :: chunk_face
    character(len=*), intent(in) :: nc_fname
    integer(I4B), dimension(2) :: dis_shape
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    integer(I4B) :: axis_sz, nvals, k
    integer(I4B), dimension(:), allocatable :: var_id
    character(len=LINELENGTH) :: longname_l, varname_l

    if (shapestr == 'NCPL') then
      ! set names
      varname_l = export_varname(nc_varname)
      longname_l = export_longname(longname, pkgname, tagname, 0)

      allocate (var_id(1))
      axis_sz = dim_ids%nmesh_face

      ! reenter define mode and create variable
      call nf_verify(nf90_redef(ncid), nc_fname)
      call nf_verify(nf90_def_var(ncid, varname_l, NF90_DOUBLE, &
                                  (/axis_sz/), var_id(1)), &
                     nc_fname)

      ! apply chunking parameters
      call ncvar_chunk(ncid, var_id(1), chunk_face, nc_fname)
      ! deflate and shuffle
      call ncvar_deflate(ncid, var_id(1), deflate, shuffle, nc_fname)

      ! put attr
      call nf_verify(nf90_put_att(ncid, var_id(1), '_FillValue', &
                                  (/NF90_FILL_DOUBLE/)), nc_fname)
      call nf_verify(nf90_put_att(ncid, var_id(1), 'long_name', &
                                  longname_l), nc_fname)

      ! add grid mapping and mf6 attr
      call ncvar_gridmap(ncid, var_id(1), gridmap_name, nc_fname)
      call ncvar_mf6attr(ncid, var_id(1), 0, 0, 0, nc_tag, nc_fname)

      ! exit define mode and write data
      call nf_verify(nf90_enddef(ncid), nc_fname)
      call nf_verify(nf90_put_var(ncid, var_id(1), p_mem), &
                     nc_fname)

    else
      allocate (var_id(dis%nlay))

      ! reenter define mode and create variable
      call nf_verify(nf90_redef(ncid), nc_fname)
      do k = 1, dis%nlay
        ! set names
        varname_l = export_varname(nc_varname, layer=k)
        longname_l = export_longname(longname, pkgname, tagname, k)

        call nf_verify(nf90_def_var(ncid, varname_l, NF90_DOUBLE, &
                                    (/dim_ids%nmesh_face/), var_id(k)), &
                       nc_fname)

        ! apply chunking parameters
        call ncvar_chunk(ncid, var_id(k), chunk_face, nc_fname)
        ! deflate and shuffle
        call ncvar_deflate(ncid, var_id(k), deflate, shuffle, nc_fname)

        ! put attr
        call nf_verify(nf90_put_att(ncid, var_id(k), '_FillValue', &
                                    (/NF90_FILL_DOUBLE/)), nc_fname)
        call nf_verify(nf90_put_att(ncid, var_id(k), 'long_name', &
                                    longname_l), nc_fname)

        ! add grid mapping and mf6 attr
        call ncvar_gridmap(ncid, var_id(k), gridmap_name, nc_fname)
        call ncvar_mf6attr(ncid, var_id(k), k, 0, 0, nc_tag, nc_fname)
      end do

      ! reshape input
      dis_shape(1) = dis%ncpl
      dis_shape(2) = dis%nlay
      nvals = product(dis_shape)
      dbl2d(1:dis_shape(1), 1:dis_shape(2)) => p_mem(1:nvals)

      ! exit define mode and write data
      call nf_verify(nf90_enddef(ncid), nc_fname)
      do k = 1, dis%nlay
        call nf_verify(nf90_put_var(ncid, var_id(k), dbl2d(:, k)), nc_fname)
      end do

      ! cleanup
      deallocate (var_id)
    end if
  end subroutine nc_export_dbl1d

  !> @brief netcdf export 2D double array
  !<
  subroutine nc_export_dbl2d(ncid, dim_ids, var_ids, disv, p_mem, nc_varname, &
                             pkgname, tagname, gridmap_name, shapestr, longname, &
                             nc_tag, deflate, shuffle, chunk_face, iper, iaux, &
                             nc_fname)
    use ConstantsModule, only: DNODATA
    integer(I4B), intent(in) :: ncid
    type(MeshNCDimIdType), intent(inout) :: dim_ids
    type(MeshNCVarIdType), intent(inout) :: var_ids
    type(DisvType), pointer, intent(in) :: disv
    real(DP), dimension(:, :), pointer, contiguous, intent(in) :: p_mem
    character(len=*), intent(in) :: nc_varname
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: gridmap_name
    character(len=*), intent(in) :: shapestr
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: nc_tag
    integer(I4B), intent(in) :: deflate
    integer(I4B), intent(in) :: shuffle
    integer(I4B), intent(in) :: chunk_face
    integer(I4B), intent(in) :: iper
    integer(I4B), intent(in) :: iaux
    character(len=*), intent(in) :: nc_fname
    integer(I4B), dimension(:), allocatable :: var_id
    character(len=LINELENGTH) :: longname_l, varname_l
    integer(I4B) :: k
    real(DP) :: fill_value

    if (iper > 0) then
      fill_value = DNODATA
    else
      fill_value = NF90_FILL_DOUBLE
    end if

    allocate (var_id(disv%nlay))

    ! reenter define mode and create variable
    call nf_verify(nf90_redef(ncid), nc_fname)
    do k = 1, disv%nlay
      ! set names
      varname_l = export_varname(nc_varname, layer=k, iper=iper, iaux=iaux)
      longname_l = export_longname(longname, pkgname, tagname, layer=k, iper=iper)

      call nf_verify(nf90_def_var(ncid, varname_l, NF90_DOUBLE, &
                                  (/dim_ids%nmesh_face/), var_id(k)), &
                     nc_fname)

      ! apply chunking parameters
      call ncvar_chunk(ncid, var_id(k), chunk_face, nc_fname)
      ! deflate and shuffle
      call ncvar_deflate(ncid, var_id(k), deflate, shuffle, nc_fname)

      ! put attr
      call nf_verify(nf90_put_att(ncid, var_id(k), '_FillValue', &
                                  (/fill_value/)), nc_fname)
      call nf_verify(nf90_put_att(ncid, var_id(k), 'long_name', &
                                  longname_l), nc_fname)

      ! add grid mapping and mf6 attr
      call ncvar_gridmap(ncid, var_id(k), gridmap_name, nc_fname)
      call ncvar_mf6attr(ncid, var_id(k), k, iper, iaux, nc_tag, nc_fname)
    end do

    ! exit define mode and write data
    call nf_verify(nf90_enddef(ncid), nc_fname)
    do k = 1, disv%nlay
      call nf_verify(nf90_put_var(ncid, var_id(k), p_mem(:, k)), nc_fname)
    end do

    deallocate (var_id)
  end subroutine nc_export_dbl2d

end module MeshDisvModelModule
