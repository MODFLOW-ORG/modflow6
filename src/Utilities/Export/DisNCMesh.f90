!> @brief This module contains the MeshDisModelModule
!!
!! This module defines UGRID layered mesh compliant netcdf
!! export type for DIS models. It is dependent on netcdf
!! libraries.
!!
!<
module MeshDisModelModule

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
  use DisModule, only: DisType
  use NetCDFCommonModule, only: nf_verify
  use netcdf

  implicit none
  private
  public :: Mesh2dDisExportType

  ! UGRID layered mesh (ULM) DIS
  type, extends(Mesh2dModelType) :: Mesh2dDisExportType
    type(DisType), pointer :: dis => null() !< pointer to model dis package
    integer(I4B) :: x_dim !< ncol dimension id
    integer(I4B) :: y_dim !< nrow dimension id
  contains
    procedure :: init => dis_export_init
    procedure :: destroy => dis_export_destroy
    procedure :: df
    procedure :: step
    procedure :: export_input_array
    procedure :: package_step_ilayer
    procedure :: package_step
    procedure :: export_layer_3d
    procedure :: define_dim
    procedure :: add_mesh_data
  end type Mesh2dDisExportType

contains

  !> @brief netcdf export dis init
  !<
  subroutine dis_export_init(this, modelname, modeltype, modelfname, nc_fname, &
                             disenum, nctype, iout)
    use ArrayHandlersModule, only: expandarray
    class(Mesh2dDisExportType), intent(inout) :: this
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelfname
    character(len=*), intent(in) :: nc_fname
    integer(I4B), intent(in) :: disenum
    integer(I4B), intent(in) :: nctype
    integer(I4B), intent(in) :: iout

    ! set nlay
    this%nlay = this%dis%nlay

    ! allocate var_id arrays
    allocate (this%var_ids%dependent(this%nlay))

    ! initialize base class
    call this%mesh_init(modelname, modeltype, modelfname, nc_fname, disenum, &
                        nctype, this%dis%lenuni, iout)
  end subroutine dis_export_init

  !> @brief netcdf export dis destroy
  !<
  subroutine dis_export_destroy(this)
    class(Mesh2dDisExportType), intent(inout) :: this
    deallocate (this%var_ids%dependent)
    ! destroy base class
    call this%mesh_destroy()
    call this%NCModelExportType%destroy()
  end subroutine dis_export_destroy

  !> @brief netcdf export define
  !<
  subroutine df(this)
    use ConstantsModule, only: MVALIDATE
    use SimVariablesModule, only: isim_mode
    class(Mesh2dDisExportType), intent(inout) :: this
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
    class(Mesh2dDisExportType), intent(inout) :: this
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    integer(I4B) :: n, k, nvals
    integer(I4B), dimension(2) :: dis_shape
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d

    ! initialize
    nullify (dbl1d)
    nullify (dbl2d)

    ! increment step
    this%stepcnt = this%stepcnt + 1

    dis_shape(1) = this%dis%ncol * this%dis%nrow
    dis_shape(2) = this%dis%nlay

    nvals = product(dis_shape)

    ! add data to dependent variable
    if (size(this%dis%nodeuser) < &
        size(this%dis%nodereduced)) then
      ! allocate nodereduced size 1d array
      allocate (dbl1d(size(this%dis%nodereduced)))

      ! initialize DHNOFLO for non-active cells
      dbl1d = DHNOFLO

      ! update active cells
      do n = 1, size(this%dis%nodereduced)
        if (this%dis%nodereduced(n) > 0) then
          dbl1d(n) = this%x(this%dis%nodereduced(n))
        end if
      end do

      dbl2d(1:dis_shape(1), 1:dis_shape(2)) => dbl1d(1:nvals)
    else
      dbl2d(1:dis_shape(1), 1:dis_shape(2)) => this%x(1:nvals)
    end if

    do k = 1, this%dis%nlay
      ! extend array with step data
      call nf_verify(nf90_put_var(this%ncid, &
                                  this%var_ids%dependent(k), dbl2d(:, k), &
                                  start=(/1, this%stepcnt/), &
                                  count=(/(this%dis%ncol * this%dis%nrow), 1/)), &
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
    class(Mesh2dDisExportType), intent(inout) :: this
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
        call nc_export_int1d(this%ncid, this%dim_ids, this%x_dim, this%y_dim, &
                             this%var_ids, this%dis, int1d, nc_varname, &
                             export_pkg%mf6_input%subcomponent_name, &
                             idt%tagname, this%gridmap_name, idt%shape, &
                             idt%longname, input_attr, this%deflate, &
                             this%shuffle, this%chunk_face, kper, this%nc_fname)
      case ('DOUBLE1D')
        call mem_setptr(dbl1d, idt%mf6varname, export_pkg%mf6_input%mempath)
        call this%export_layer_3d(export_pkg, idt, ilayer_read, ialayer, &
                                  dbl1d, nc_varname, input_attr)
      case ('DOUBLE2D')
        call mem_setptr(dbl2d, idt%mf6varname, export_pkg%mf6_input%mempath)
        nvals = this%dis%ncol * this%dis%nrow

        do n = 1, size(dbl2d, dim=1) !naux
          dbl1d_ptr(1:nvals) => dbl2d(n, :)
          if (all(dbl1d_ptr == DZERO)) then
          else
            call this%export_layer_3d(export_pkg, idt, ilayer_read, ialayer, &
                                      dbl1d_ptr, nc_varname, input_attr, n)
          end if
        end do
      case default
        errmsg = 'EXPORT ilayaer unsupported datatype='//trim(idt%datatype)
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
    class(Mesh2dDisExportType), intent(inout) :: this
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
  subroutine export_layer_3d(this, export_pkg, idt, ilayer_read, ialayer, &
                             dbl1d, nc_varname, input_attr, iaux)
    use ConstantsModule, only: DNODATA, DZERO
    use NCModelExportModule, only: ExportPackageType
    class(Mesh2dDisExportType), intent(inout) :: this
    class(ExportPackageType), pointer, intent(in) :: export_pkg
    type(InputParamDefinitionType), pointer, intent(in) :: idt
    logical(LGP), intent(in) :: ilayer_read
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: ialayer
    real(DP), dimension(:), pointer, contiguous, intent(in) :: dbl1d
    character(len=*), intent(in) :: nc_varname
    character(len=*), intent(in) :: input_attr
    integer(I4B), optional, intent(in) :: iaux
    real(DP), dimension(:, :, :), pointer, contiguous :: dbl3d
    integer(I4B) :: n, i, j, k, nvals, idxaux
    real(DP), dimension(:, :), contiguous, pointer :: dbl2d_ptr

    ! initialize
    idxaux = 0
    if (present(iaux)) then
      idxaux = iaux
    end if

    allocate (dbl3d(export_pkg%mshape(3), export_pkg%mshape(2), &
                    export_pkg%mshape(1)))

    if (ilayer_read) then
      do k = 1, size(dbl3d, dim=3)
        n = 0
        do i = 1, size(dbl3d, dim=2)
          do j = 1, size(dbl3d, dim=1)
            n = n + 1
            if (ialayer(n) == k) then
              dbl3d(j, i, k) = dbl1d(n)
            else
              dbl3d(j, i, k) = DNODATA
            end if
          end do
        end do
      end do
    else
      dbl3d = DNODATA
      nvals = export_pkg%mshape(3) * export_pkg%mshape(2)
      dbl2d_ptr(1:export_pkg%mshape(3), 1:export_pkg%mshape(2)) => dbl1d(1:nvals)
      dbl3d(:, :, 1) = dbl2d_ptr(:, :)
    end if

    call nc_export_dbl3d(this%ncid, this%dim_ids, this%var_ids, this%dis, dbl3d, &
                         nc_varname, export_pkg%mf6_input%subcomponent_name, &
                         idt%tagname, this%gridmap_name, idt%shape, &
                         idt%longname, input_attr, this%deflate, this%shuffle, &
                         this%chunk_face, export_pkg%iper, idxaux, this%nc_fname)

    deallocate (dbl3d)
  end subroutine export_layer_3d

  !> @brief netcdf export an input array
  !<
  subroutine export_input_array(this, pkgtype, pkgname, mempath, idt)
    class(Mesh2dDisExportType), intent(inout) :: this
    character(len=*), intent(in) :: pkgtype
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: mempath
    type(InputParamDefinitionType), pointer, intent(in) :: idt
    integer(I4B), dimension(:), pointer, contiguous :: int1d
    integer(I4B), dimension(:, :), pointer, contiguous :: int2d
    integer(I4B), dimension(:, :, :), pointer, contiguous :: int3d
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    real(DP), dimension(:, :, :), pointer, contiguous :: dbl3d
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
      call nc_export_int1d(this%ncid, this%dim_ids, this%x_dim, this%y_dim, &
                           this%var_ids, this%dis, int1d, nc_varname, pkgname, &
                           idt%tagname, this%gridmap_name, idt%shape, &
                           idt%longname, input_attr, this%deflate, this%shuffle, &
                           this%chunk_face, iper, this%nc_fname)
    case ('INTEGER2D')
      call mem_setptr(int2d, idt%mf6varname, mempath)
      call nc_export_int2d(this%ncid, this%dim_ids, this%var_ids, this%dis, &
                           int2d, nc_varname, pkgname, idt%tagname, &
                           this%gridmap_name, idt%shape, idt%longname, &
                           input_attr, this%deflate, this%shuffle, &
                           this%chunk_face, this%nc_fname)
    case ('INTEGER3D')
      call mem_setptr(int3d, idt%mf6varname, mempath)
      call nc_export_int3d(this%ncid, this%dim_ids, this%var_ids, this%dis, &
                           int3d, nc_varname, pkgname, idt%tagname, &
                           this%gridmap_name, idt%shape, idt%longname, &
                           input_attr, this%deflate, this%shuffle, &
                           this%chunk_face, this%nc_fname)
    case ('DOUBLE1D')
      call mem_setptr(dbl1d, idt%mf6varname, mempath)
      call nc_export_dbl1d(this%ncid, this%dim_ids, this%x_dim, this%y_dim, &
                           this%var_ids, this%dis, dbl1d, nc_varname, pkgname, &
                           idt%tagname, this%gridmap_name, idt%shape, &
                           idt%longname, input_attr, this%deflate, this%shuffle, &
                           this%chunk_face, this%nc_fname)
    case ('DOUBLE2D')
      call mem_setptr(dbl2d, idt%mf6varname, mempath)
      call nc_export_dbl2d(this%ncid, this%dim_ids, this%var_ids, this%dis, &
                           dbl2d, nc_varname, pkgname, idt%tagname, &
                           this%gridmap_name, idt%shape, idt%longname, &
                           input_attr, this%deflate, this%shuffle, &
                           this%chunk_face, this%nc_fname)
    case ('DOUBLE3D')
      call mem_setptr(dbl3d, idt%mf6varname, mempath)
      call nc_export_dbl3d(this%ncid, this%dim_ids, this%var_ids, this%dis, &
                           dbl3d, nc_varname, pkgname, idt%tagname, &
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
    class(Mesh2dDisExportType), intent(inout) :: this

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
    call nf_verify(nf90_def_dim(this%ncid, 'nmesh_node', &
                                ((this%dis%ncol + 1) * (this%dis%nrow + 1)), &
                                this%dim_ids%nmesh_node), this%nc_fname)
    call nf_verify(nf90_def_dim(this%ncid, 'nmesh_face', &
                                (this%dis%ncol * this%dis%nrow), &
                                this%dim_ids%nmesh_face), this%nc_fname)
    call nf_verify(nf90_def_dim(this%ncid, 'max_nmesh_face_nodes', 4, &
                                this%dim_ids%max_nmesh_face_nodes), &
                   this%nc_fname)

    ! x, y, nlay
    call nf_verify(nf90_def_dim(this%ncid, 'nlay', this%dis%nlay, &
                                this%dim_ids%nlay), this%nc_fname)
    call nf_verify(nf90_def_dim(this%ncid, 'x', this%dis%ncol, &
                                this%x_dim), this%nc_fname)
    call nf_verify(nf90_def_dim(this%ncid, 'y', this%dis%nrow, &
                                this%y_dim), this%nc_fname)
  end subroutine define_dim

  !> @brief netcdf export add mesh information
  !<
  subroutine add_mesh_data(this)
    use BaseDisModule, only: dis_transform_xy
    class(Mesh2dDisExportType), intent(inout) :: this
    integer(I4B) :: cnt, maxvert, m
    integer(I4B), dimension(:), allocatable :: verts
    real(DP), dimension(:), allocatable :: bnds
    integer(I4B) :: i, j
    real(DP) :: x, y, x_transform, y_transform
    real(DP), dimension(:), allocatable :: node_x, node_y
    real(DP), dimension(:), allocatable :: cell_x, cell_y

    ! initialize max vertices required to define cell
    maxvert = 4

    ! set mesh container variable value to 1
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh, 1), &
                   this%nc_fname)

    ! allocate temporary arrays
    allocate (verts(maxvert))
    allocate (bnds(maxvert))
    allocate (node_x(((this%dis%ncol + 1) * (this%dis%nrow + 1))))
    allocate (node_y(((this%dis%ncol + 1) * (this%dis%nrow + 1))))
    allocate (cell_x((this%dis%ncol * this%dis%nrow)))
    allocate (cell_y((this%dis%ncol * this%dis%nrow)))

    ! set node_x and node_y arrays
    cnt = 0
    node_x = NF90_FILL_DOUBLE
    node_y = NF90_FILL_DOUBLE
    y = sum(this%dis%delc)
    do j = this%dis%nrow, 0, -1
      x = 0
      do i = this%dis%ncol, 0, -1
        cnt = cnt + 1
        call dis_transform_xy(x, y, &
                              this%dis%xorigin, &
                              this%dis%yorigin, &
                              this%dis%angrot, &
                              x_transform, y_transform)
        node_x(cnt) = x_transform
        node_y(cnt) = y_transform
        if (i > 0) x = x + this%dis%delr(i)
      end do
      if (j > 0) y = y - this%dis%delc(j)
    end do

    ! write node_x and node_y arrays to netcdf file
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_node_x, node_x), &
                   this%nc_fname)
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_node_y, node_y), &
                   this%nc_fname)

    ! set cell_x and cell_y arrays
    cnt = 1
    cell_x = NF90_FILL_DOUBLE
    cell_y = NF90_FILL_DOUBLE
    do j = 1, this%dis%nrow
      y = this%dis%celly(j)
      do i = 1, this%dis%ncol
        x = this%dis%cellx(i)
        call dis_transform_xy(x, y, &
                              this%dis%xorigin, &
                              this%dis%yorigin, &
                              this%dis%angrot, &
                              x_transform, y_transform)
        cell_x(cnt) = x_transform
        cell_y(cnt) = y_transform
        cnt = cnt + 1
      end do
    end do

    ! write face_x and face_y arrays to netcdf file
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_face_x, cell_x), &
                   this%nc_fname)
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_face_y, cell_y), &
                   this%nc_fname)

    ! set face nodes array
    cnt = 0
    do i = 1, this%dis%nrow
      do j = 1, this%dis%ncol
        cnt = cnt + 1
        verts = NF90_FILL_INT
        verts(1) = cnt + this%dis%ncol + i
        verts(2) = cnt + this%dis%ncol + i + 1
        if (i > 1) then
          verts(3) = cnt + i
          verts(4) = cnt + i - 1
        else
          verts(3) = cnt + 1
          verts(4) = cnt
        end if

        ! write face nodes array to netcdf file
        call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_face_nodes, &
                                    verts, start=(/1, cnt/), &
                                    count=(/maxvert, 1/)), &
                       this%nc_fname)

        ! set face y bounds array
        bnds = NF90_FILL_DOUBLE
        do m = 1, size(bnds)
          if (verts(m) /= NF90_FILL_INT) then
            bnds(m) = node_y(verts(m))
          end if
          ! write face y bounds array to netcdf file
          call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_face_ybnds, &
                                      bnds, start=(/1, cnt/), &
                                      count=(/maxvert, 1/)), &
                         this%nc_fname)
        end do

        ! set face x bounds array
        bnds = NF90_FILL_DOUBLE
        do m = 1, size(bnds)
          if (verts(m) /= NF90_FILL_INT) then
            bnds(m) = node_x(verts(m))
          end if
          ! write face x bounds array to netcdf file
          call nf_verify(nf90_put_var(this%ncid, this%var_ids%mesh_face_xbnds, &
                                      bnds, start=(/1, cnt/), &
                                      count=(/maxvert, 1/)), &
                         this%nc_fname)
        end do
      end do
    end do

    ! cleanup
    deallocate (bnds)
    deallocate (verts)
    deallocate (node_x)
    deallocate (node_y)
    deallocate (cell_x)
    deallocate (cell_y)
  end subroutine add_mesh_data

  !> @brief netcdf export 1D integer
  !<
  subroutine nc_export_int1d(ncid, dim_ids, x_dim, y_dim, var_ids, dis, p_mem, &
                             nc_varname, pkgname, tagname, gridmap_name, &
                             shapestr, longname, nc_tag, deflate, shuffle, &
                             chunk_face, iper, nc_fname)
    integer(I4B), intent(in) :: ncid
    type(MeshNCDimIdType), intent(inout) :: dim_ids
    integer(I4B), intent(in) :: x_dim
    integer(I4B), intent(in) :: y_dim
    type(MeshNCVarIdType), intent(inout) :: var_ids
    type(DisType), pointer, intent(in) :: dis
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
    integer(I4B), dimension(3) :: dis_shape
    integer(I4B), dimension(1) :: layer_shape
    integer(I4B), dimension(:, :, :), pointer, contiguous :: int3d
    integer(I4B), dimension(:), pointer, contiguous :: int1d
    integer(I4B) :: axis_dim, nvals, k
    integer(I4B), dimension(:), allocatable :: var_id
    character(len=LINELENGTH) :: longname_l, varname_l

    if (shapestr == 'NROW' .or. &
        shapestr == 'NCOL' .or. &
        shapestr == 'NCPL') then

      select case (shapestr)
      case ('NROW')
        axis_dim = y_dim
      case ('NCOL')
        axis_dim = x_dim
      case ('NCPL')
        axis_dim = dim_ids%nmesh_face
      end select

      ! set names
      varname_l = export_varname(nc_varname, layer=0, iper=iper)
      longname_l = export_longname(longname, pkgname, tagname, layer=0, iper=iper)

      allocate (var_id(1))

      ! reenter define mode and create variable
      call nf_verify(nf90_redef(ncid), nc_fname)
      call nf_verify(nf90_def_var(ncid, varname_l, NF90_INT, &
                                  (/axis_dim/), var_id(1)), &
                     nc_fname)

      ! NROW/NCOL shapes use default chunking
      call ncvar_deflate(ncid, var_id(1), deflate, shuffle, nc_fname)

      ! put attr
      call nf_verify(nf90_put_att(ncid, var_id(1), '_FillValue', &
                                  (/NF90_FILL_INT/)), nc_fname)
      call nf_verify(nf90_put_att(ncid, var_id(1), 'long_name', &
                                  longname_l), nc_fname)

      ! add mf6 attr
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
      dis_shape(1) = dis%ncol
      dis_shape(2) = dis%nrow
      dis_shape(3) = dis%nlay
      nvals = product(dis_shape)
      int3d(1:dis_shape(1), 1:dis_shape(2), 1:dis_shape(3)) => p_mem(1:nvals)

      ! exit define mode and write data
      call nf_verify(nf90_enddef(ncid), nc_fname)
      layer_shape(1) = dis%nrow * dis%ncol
      do k = 1, dis%nlay
        int1d(1:layer_shape(1)) => int3d(:, :, k)
        call nf_verify(nf90_put_var(ncid, var_id(k), int1d), nc_fname)
      end do

      ! cleanup
      deallocate (var_id)
    end if
  end subroutine nc_export_int1d

  !> @brief netcdf export 2D integer
  !<
  subroutine nc_export_int2d(ncid, dim_ids, var_ids, dis, p_mem, nc_varname, &
                             pkgname, tagname, gridmap_name, shapestr, longname, &
                             nc_tag, deflate, shuffle, chunk_face, nc_fname)
    integer(I4B), intent(in) :: ncid
    type(MeshNCDimIdType), intent(inout) :: dim_ids
    type(MeshNCVarIdType), intent(inout) :: var_ids
    type(DisType), pointer, intent(in) :: dis
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
    integer(I4B) :: var_id
    integer(I4B), dimension(:), pointer, contiguous :: int1d
    integer(I4B), dimension(1) :: layer_shape
    character(len=LINELENGTH) :: longname_l, varname_l

    ! set names
    varname_l = export_varname(nc_varname)
    longname_l = export_longname(longname, pkgname, tagname, 0)

    ! reenter define mode and create variable
    call nf_verify(nf90_redef(ncid), nc_fname)
    call nf_verify(nf90_def_var(ncid, varname_l, NF90_INT, &
                                (/dim_ids%nmesh_face/), var_id), &
                   nc_fname)

    ! apply chunking parameters
    call ncvar_chunk(ncid, var_id, chunk_face, nc_fname)
    ! deflate and shuffle
    call ncvar_deflate(ncid, var_id, deflate, shuffle, nc_fname)

    ! put attr
    call nf_verify(nf90_put_att(ncid, var_id, '_FillValue', &
                                (/NF90_FILL_INT/)), nc_fname)
    call nf_verify(nf90_put_att(ncid, var_id, 'long_name', &
                                longname_l), nc_fname)

    ! add grid mapping and mf6 attr
    call ncvar_gridmap(ncid, var_id, gridmap_name, nc_fname)
    call ncvar_mf6attr(ncid, var_id, 0, 0, 0, nc_tag, nc_fname)

    ! exit define mode and write data
    call nf_verify(nf90_enddef(ncid), nc_fname)
    layer_shape(1) = dis%nrow * dis%ncol
    int1d(1:layer_shape(1)) => p_mem
    call nf_verify(nf90_put_var(ncid, var_id, int1d), nc_fname)
  end subroutine nc_export_int2d

  !> @brief netcdf export 3D integer
  !<
  subroutine nc_export_int3d(ncid, dim_ids, var_ids, dis, p_mem, nc_varname, &
                             pkgname, tagname, gridmap_name, shapestr, longname, &
                             nc_tag, deflate, shuffle, chunk_face, nc_fname)
    integer(I4B), intent(in) :: ncid
    type(MeshNCDimIdType), intent(inout) :: dim_ids
    type(MeshNCVarIdType), intent(inout) :: var_ids
    type(DisType), pointer, intent(in) :: dis
    integer(I4B), dimension(:, :, :), pointer, contiguous, intent(in) :: p_mem
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
    integer(I4B), dimension(:), pointer, contiguous :: int1d
    character(len=LINELENGTH) :: longname_l, varname_l
    integer(I4B), dimension(1) :: layer_shape
    integer(I4B) :: k

    allocate (var_id(dis%nlay))

    ! reenter define mode and create variable
    call nf_verify(nf90_redef(ncid), nc_fname)
    do k = 1, dis%nlay
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
    layer_shape(1) = dis%nrow * dis%ncol
    do k = 1, dis%nlay
      int1d(1:layer_shape(1)) => p_mem(:, :, k)
      call nf_verify(nf90_put_var(ncid, var_id(k), int1d), nc_fname)
    end do

    ! cleanup
    deallocate (var_id)
  end subroutine nc_export_int3d

  !> @brief netcdf export 1D double
  !<
  subroutine nc_export_dbl1d(ncid, dim_ids, x_dim, y_dim, var_ids, dis, p_mem, &
                             nc_varname, pkgname, tagname, gridmap_name, &
                             shapestr, longname, nc_tag, deflate, shuffle, &
                             chunk_face, nc_fname)
    integer(I4B), intent(in) :: ncid
    type(MeshNCDimIdType), intent(inout) :: dim_ids
    integer(I4B), intent(in) :: x_dim
    integer(I4B), intent(in) :: y_dim
    type(MeshNCVarIdType), intent(inout) :: var_ids
    type(DisType), pointer, intent(in) :: dis
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
    integer(I4B), dimension(3) :: dis_shape
    integer(I4B), dimension(1) :: layer_shape
    real(DP), dimension(:, :, :), pointer, contiguous :: dbl3d
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    integer(I4B) :: axis_dim, nvals, k
    integer(I4B), dimension(:), allocatable :: var_id
    character(len=LINELENGTH) :: longname_l, varname_l

    if (shapestr == 'NROW' .or. &
        shapestr == 'NCOL') then ! .or. &
      !shapestr == 'NCPL') then

      select case (shapestr)
      case ('NROW')
        axis_dim = y_dim
      case ('NCOL')
        axis_dim = x_dim
        !case ('NCPL')
        !  axis_dim = dim_ids%nmesh_face
      end select

      ! set names
      varname_l = export_varname(nc_varname)
      longname_l = export_longname(longname, pkgname, tagname, 0)

      allocate (var_id(1))

      ! reenter define mode and create variable
      call nf_verify(nf90_redef(ncid), nc_fname)
      call nf_verify(nf90_def_var(ncid, varname_l, NF90_DOUBLE, &
                                  (/axis_dim/), var_id(1)), &
                     nc_fname)

      ! NROW/NCOL shapes use default chunking
      call ncvar_deflate(ncid, var_id(1), deflate, shuffle, nc_fname)

      ! put attr
      call nf_verify(nf90_put_att(ncid, var_id(1), '_FillValue', &
                                  (/NF90_FILL_DOUBLE/)), nc_fname)
      call nf_verify(nf90_put_att(ncid, var_id(1), 'long_name', &
                                  longname_l), nc_fname)

      ! add mf6 attr
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
        ! defalte and shuffle
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
      dis_shape(1) = dis%ncol
      dis_shape(2) = dis%nrow
      dis_shape(3) = dis%nlay
      nvals = product(dis_shape)
      dbl3d(1:dis_shape(1), 1:dis_shape(2), 1:dis_shape(3)) => p_mem(1:nvals)

      ! exit define mode and write data
      call nf_verify(nf90_enddef(ncid), nc_fname)
      layer_shape(1) = dis%nrow * dis%ncol
      do k = 1, dis%nlay
        dbl1d(1:layer_shape(1)) => dbl3d(:, :, k)
        call nf_verify(nf90_put_var(ncid, var_id(k), dbl1d), nc_fname)
      end do

      ! cleanup
      deallocate (var_id)
    end if
  end subroutine nc_export_dbl1d

  !> @brief netcdf export 2D double
  !<
  subroutine nc_export_dbl2d(ncid, dim_ids, var_ids, dis, p_mem, nc_varname, &
                             pkgname, tagname, gridmap_name, shapestr, longname, &
                             nc_tag, deflate, shuffle, chunk_face, nc_fname)
    integer(I4B), intent(in) :: ncid
    type(MeshNCDimIdType), intent(inout) :: dim_ids
    type(MeshNCVarIdType), intent(inout) :: var_ids
    type(DisType), pointer, intent(in) :: dis
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
    character(len=*), intent(in) :: nc_fname
    integer(I4B) :: var_id
    character(len=LINELENGTH) :: longname_l, varname_l
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    integer(I4B), dimension(1) :: layer_shape

    ! set names
    varname_l = export_varname(nc_varname)
    longname_l = export_longname(longname, pkgname, tagname, 0)

    ! reenter define mode and create variable
    call nf_verify(nf90_redef(ncid), nc_fname)
    call nf_verify(nf90_def_var(ncid, varname_l, NF90_DOUBLE, &
                                (/dim_ids%nmesh_face/), var_id), &
                   nc_fname)

    ! apply chunking parameters
    call ncvar_chunk(ncid, var_id, chunk_face, nc_fname)
    ! deflate and shuffle
    call ncvar_deflate(ncid, var_id, deflate, shuffle, nc_fname)

    ! put attr
    call nf_verify(nf90_put_att(ncid, var_id, '_FillValue', &
                                (/NF90_FILL_DOUBLE/)), nc_fname)
    call nf_verify(nf90_put_att(ncid, var_id, 'long_name', &
                                longname_l), nc_fname)

    ! add grid mapping and mf6 attr
    call ncvar_gridmap(ncid, var_id, gridmap_name, nc_fname)
    call ncvar_mf6attr(ncid, var_id, 0, 0, 0, nc_tag, nc_fname)

    ! exit define mode and write data
    call nf_verify(nf90_enddef(ncid), nc_fname)
    layer_shape(1) = dis%nrow * dis%ncol
    dbl1d(1:layer_shape(1)) => p_mem
    call nf_verify(nf90_put_var(ncid, var_id, dbl1d), nc_fname)
  end subroutine nc_export_dbl2d

  !> @brief netcdf export 3D double
  !<
  subroutine nc_export_dbl3d(ncid, dim_ids, var_ids, dis, p_mem, nc_varname, &
                             pkgname, tagname, gridmap_name, shapestr, longname, &
                             nc_tag, deflate, shuffle, chunk_face, iper, iaux, &
                             nc_fname)
    use ConstantsModule, only: DNODATA
    integer(I4B), intent(in) :: ncid
    type(MeshNCDimIdType), intent(inout) :: dim_ids
    type(MeshNCVarIdType), intent(inout) :: var_ids
    type(DisType), pointer, intent(in) :: dis
    real(DP), dimension(:, :, :), pointer, contiguous, intent(in) :: p_mem
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
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    character(len=LINELENGTH) :: longname_l, varname_l
    integer(I4B), dimension(1) :: layer_shape
    integer(I4B) :: k
    real(DP) :: fill_value

    if (iper > 0) then
      fill_value = DNODATA
    else
      fill_value = NF90_FILL_DOUBLE
    end if

    allocate (var_id(dis%nlay))

    ! reenter define mode and create variable
    call nf_verify(nf90_redef(ncid), nc_fname)
    do k = 1, dis%nlay
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
      !end if
    end do

    ! exit define mode and write data
    call nf_verify(nf90_enddef(ncid), nc_fname)
    layer_shape(1) = dis%nrow * dis%ncol
    do k = 1, dis%nlay
      dbl1d(1:layer_shape(1)) => p_mem(:, :, k)
      call nf_verify(nf90_put_var(ncid, var_id(k), dbl1d), nc_fname)
    end do

    ! cleanup
    deallocate (var_id)
  end subroutine nc_export_dbl3d

end module MeshDisModelModule
