!> @brief This module contains the DisNCStructuredModule
!!
!! This module defines a STRUCTURED (non-ugrid) netcdf
!! export type for DIS models. It is dependent on netcdf
!! libraries.
!!
!<
module DisNCStructuredModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENBIGLINE, LENCOMPONENTNAME, &
                             LENMEMPATH, LENVARNAME, DNODATA, DZERO
  use SimVariablesModule, only: errmsg, warnmsg
  use SimModule, only: store_error, store_warning, store_error_filename
  use MemoryManagerModule, only: mem_setptr
  use InputDefinitionModule, only: InputParamDefinitionType
  use CharacterStringModule, only: CharacterStringType
  use NCModelExportModule, only: NCBaseModelExportType, export_longname
  use DisModule, only: DisType
  use NetCDFCommonModule, only: nf_verify
  use netcdf

  implicit none
  private
  public :: DisNCStructuredType

  type :: StructuredNCDimIdType
    integer(I4B) :: x !< number of columns
    integer(I4B) :: y !< number of rows
    integer(I4B) :: z !< number of layers
    integer(I4B) :: ncpl !< number of cells in layer
    integer(I4B) :: time !< number of steps
    integer(I4B) :: bnd !< number in boundary
  contains
  end type StructuredNCDimIdType

  type :: StructuredNCVarIdType
    integer(I4B) :: x !< x coordinate variable
    integer(I4B) :: y !< y coordinate variable
    integer(I4B) :: z !< z coordinate variable
    integer(I4B) :: time !< time coordinate variable
    integer(I4B) :: dependent !< dependent variable
    integer(I4B) :: x_bnds !< x boundaries 2D array
    integer(I4B) :: y_bnds !< y boundaries 2D array
    integer(I4B) :: z_bnds !< z boundaries 2D array
    integer(I4B) :: latitude !< latitude 2D array
    integer(I4B) :: longitude !< longitude 2D array
  contains
  end type StructuredNCVarIdType

  type, extends(NCBaseModelExportType) :: DisNCStructuredType
    type(StructuredNCDimIdType) :: dim_ids !< structured dimension ids type
    type(StructuredNCVarIdType) :: var_ids !< structured variable ids type
    type(DisType), pointer :: dis => null() !< pointer to model dis package
    integer(I4B) :: nlay !< number of layers
    real(DP), dimension(:), pointer, contiguous :: latitude => null() !< lat input array pointer
    real(DP), dimension(:), pointer, contiguous :: longitude => null() !< lon input array pointer
    integer(I4B), pointer :: chunk_z !< chunking parameter for z dimension
    integer(I4B), pointer :: chunk_y !< chunking parameter for y dimension
    integer(I4B), pointer :: chunk_x !< chunking parameter for x dimension
    integer(I4B), dimension(:), allocatable :: layers !< layers array
    logical(LGP) :: latlon !< are lat and lon arrays to be written to netcdf file
  contains
    procedure :: init => dis_export_init
    procedure :: destroy => dis_export_destroy
    procedure :: df
    procedure :: step
    procedure :: export_input_array
    procedure :: export_input_arrays
    procedure :: package_step_ilayer
    procedure :: package_step
    procedure :: export_layer_3d
    procedure :: add_pkg_data
    procedure :: add_global_att
    procedure :: define_dim
    procedure :: define_dependent
    procedure :: define_gridmap
    procedure :: define_geocoords
    procedure :: add_proj_data
    procedure :: add_grid_data
  end type DisNCStructuredType

  interface nc_export_array
    module procedure nc_export_int1d, nc_export_int2d, &
      nc_export_int3d, nc_export_dbl1d, &
      nc_export_dbl2d, nc_export_dbl3d
  end interface nc_export_array

contains

  !> @brief netcdf export dis init
  !<
  subroutine dis_export_init(this, modelname, modeltype, modelfname, nc_fname, &
                             disenum, nctype, iout)
    use MemoryManagerModule, only: get_isize
    use MemoryManagerExtModule, only: mem_set_value
    class(DisNCStructuredType), intent(inout) :: this
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelfname
    character(len=*), intent(in) :: nc_fname
    integer(I4B), intent(in) :: disenum
    integer(I4B), intent(in) :: nctype
    integer(I4B), intent(in) :: iout
    integer(I4B) :: k, latsz, lonsz
    logical(LGP) :: found

    ! set nlay
    this%nlay = this%dis%nlay

    ! allocate
    allocate (this%chunk_z)
    allocate (this%chunk_y)
    allocate (this%chunk_x)
    allocate (this%layers(this%nlay))

    ! initialize
    this%chunk_z = -1
    this%chunk_y = -1
    this%chunk_x = -1
    do k = 1, this%nlay
      this%layers(k) = k
    end do

    this%latlon = .false.

    ! initialize base class
    call this%NCModelExportType%init(modelname, modeltype, modelfname, nc_fname, &
                                     disenum, nctype, iout)

    ! update values from input context
    if (this%ncf_mempath /= '') then
      call mem_set_value(this%chunk_z, 'CHUNK_Z', this%ncf_mempath, found)
      call mem_set_value(this%chunk_y, 'CHUNK_Y', this%ncf_mempath, found)
      call mem_set_value(this%chunk_x, 'CHUNK_X', this%ncf_mempath, found)

      if (this%chunk_time > 0 .and. this%chunk_z > 0 .and. &
          this%chunk_y > 0 .and. this%chunk_x > 0) then
        this%chunking_active = .true.
      else if (this%chunk_time > 0 .or. this%chunk_z > 0 .or. &
               this%chunk_y > 0 .or. this%chunk_x > 0) then
        this%chunk_time = -1
        this%chunk_z = -1
        this%chunk_y = -1
        this%chunk_x = -1
        write (warnmsg, '(a)') 'Ignoring user provided NetCDF chunking &
          &parameters. Define chunk_time, chunk_x, chunk_y and chunk_z input &
          &parameters to see an effect in file "'//trim(nc_fname)//'".'
        call store_warning(warnmsg)
      end if

      call get_isize('LATITUDE', this%ncf_mempath, latsz)
      call get_isize('LONGITUDE', this%ncf_mempath, lonsz)

      if (latsz > 0 .and. lonsz > 0) then
        this%latlon = .true.
        if (this%wkt /= '') then
          write (warnmsg, '(a)') 'Ignoring user provided NetCDF wkt parameter &
            &as longitude and latitude arrays have been provided. &
            &Applies to file "'//trim(nc_fname)//'".'
          call store_warning(warnmsg)
          this%wkt = ''
          this%gridmap_name = ''
        end if
        call mem_setptr(this%latitude, 'LATITUDE', this%ncf_mempath)
        call mem_setptr(this%longitude, 'LONGITUDE', this%ncf_mempath)
      end if

      if (this%wkt /= '') then
        if (this%dis%angrot /= DZERO) then
          write (warnmsg, '(a)') 'WKT parameter set with structured rotated &
            &grid. Projected coordinates will have grid local values. &
            &Applies to file "'//trim(nc_fname)//'".'
          call store_warning(warnmsg)
        end if
      end if
    end if

    if (this%dis%lenuni == 1) then
      this%lenunits = 'ft'
    else
      this%lenunits = 'm'
    end if

    ! create the netcdf file
    call nf_verify(nf90_create(this%nc_fname, &
                               IOR(NF90_CLOBBER, NF90_NETCDF4), this%ncid), &
                   this%nc_fname)
  end subroutine dis_export_init

  !> @brief netcdf export dis destroy
  !<
  subroutine dis_export_destroy(this)
    class(DisNCStructuredType), intent(inout) :: this
    call nf_verify(nf90_close(this%ncid), this%nc_fname)
    deallocate (this%chunk_z)
    deallocate (this%chunk_y)
    deallocate (this%chunk_x)
    deallocate (this%layers)
    nullify (this%chunk_z)
    nullify (this%chunk_y)
    nullify (this%chunk_x)
    ! destroy base class
    call this%NCModelExportType%destroy()
  end subroutine dis_export_destroy

  !> @brief netcdf export define
  !<
  subroutine df(this)
    use ConstantsModule, only: MVALIDATE
    use SimVariablesModule, only: isim_mode
    class(DisNCStructuredType), intent(inout) :: this
    ! put root group file scope attributes
    call this%add_global_att()
    ! define root group dimensions and coordinate variables
    call this%define_dim()
    ! define grid projection variables
    call this%define_geocoords()
    if (isim_mode /= MVALIDATE) then
      ! define the dependent variable
      call this%define_dependent()
    end if
    ! exit define mode
    call nf_verify(nf90_enddef(this%ncid), this%nc_fname)
    ! add data locations
    call this%add_grid_data()
    ! add projection data
    call this%add_proj_data()
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
    class(DisNCStructuredType), intent(inout) :: this
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    integer(I4B) :: n

    this%stepcnt = this%stepcnt + 1

    if (size(this%dis%nodeuser) < &
        size(this%dis%nodereduced)) then
      allocate (dbl1d(size(this%dis%nodereduced)))
      dbl1d = DHNOFLO
      do n = 1, size(this%dis%nodereduced)
        if (this%dis%nodereduced(n) > 0) then
          dbl1d(n) = this%x(this%dis%nodereduced(n))
        end if
      end do
      ! write step data to dependent variable
      call nf_verify(nf90_put_var(this%ncid, &
                                  this%var_ids%dependent, dbl1d, &
                                  start=(/1, 1, 1, this%stepcnt/), &
                                  count=(/this%dis%ncol, &
                                          this%dis%nrow, &
                                          this%dis%nlay, 1/)), &
                     this%nc_fname)
      deallocate (dbl1d)
    else
      ! write step data to dependent variable
      call nf_verify(nf90_put_var(this%ncid, &
                                  this%var_ids%dependent, this%x, &
                                  start=(/1, 1, 1, this%stepcnt/), &
                                  count=(/this%dis%ncol, &
                                          this%dis%nrow, &
                                          this%dis%nlay, 1/)), &
                     this%nc_fname)
    end if

    ! write to time coordinate variable
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%time, &
                                totim, start=(/this%stepcnt/)), &
                   this%nc_fname)

    ! synchronize file
    call nf_verify(nf90_sync(this%ncid), this%nc_fname)
  end subroutine step

  !> @brief netcdf export an input array
  !<
  subroutine export_input_array(this, pkgtype, pkgname, mempath, idt)
    class(DisNCStructuredType), intent(inout) :: this
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

    ! initialize
    iper = 0
    iaux = 0

    ! set variable name and input attribute string
    nc_varname = export_varname(pkgname, idt)
    input_attr = this%input_attribute(pkgname, idt)

    select case (idt%datatype)
    case ('INTEGER1D')
      call mem_setptr(int1d, idt%mf6varname, mempath)
      call nc_export_array(this%ncid, this%dim_ids, this%var_ids, this%dis, &
                           int1d, nc_varname, pkgname, idt%tagname, &
                           idt%shape, idt%longname, input_attr, &
                           this%gridmap_name, this%latlon, this%deflate, &
                           this%shuffle, this%chunk_z, this%chunk_y, &
                           this%chunk_x, iper, this%nc_fname)
    case ('INTEGER2D')
      call mem_setptr(int2d, idt%mf6varname, mempath)
      call nc_export_array(this%ncid, this%dim_ids, this%var_ids, this%dis, &
                           int2d, nc_varname, pkgname, idt%tagname, &
                           idt%shape, idt%longname, input_attr, &
                           this%gridmap_name, this%latlon, this%deflate, &
                           this%shuffle, this%chunk_z, this%chunk_y, &
                           this%chunk_x, this%nc_fname)
    case ('INTEGER3D')
      call mem_setptr(int3d, idt%mf6varname, mempath)
      call nc_export_array(this%ncid, this%dim_ids, this%var_ids, this%dis, &
                           int3d, nc_varname, pkgname, idt%tagname, &
                           idt%shape, idt%longname, input_attr, &
                           this%gridmap_name, this%latlon, this%deflate, &
                           this%shuffle, this%chunk_z, this%chunk_y, &
                           this%chunk_x, this%nc_fname)
    case ('DOUBLE1D')
      call mem_setptr(dbl1d, idt%mf6varname, mempath)
      call nc_export_array(this%ncid, this%dim_ids, this%var_ids, this%dis, &
                           dbl1d, nc_varname, pkgname, idt%tagname, &
                           idt%shape, idt%longname, input_attr, &
                           this%gridmap_name, this%latlon, this%deflate, &
                           this%shuffle, this%chunk_z, this%chunk_y, &
                           this%chunk_x, iper, this%nc_fname)
    case ('DOUBLE2D')
      call mem_setptr(dbl2d, idt%mf6varname, mempath)
      call nc_export_array(this%ncid, this%dim_ids, this%var_ids, this%dis, &
                           dbl2d, nc_varname, pkgname, idt%tagname, &
                           idt%shape, idt%longname, input_attr, &
                           this%gridmap_name, this%latlon, this%deflate, &
                           this%shuffle, this%chunk_z, this%chunk_y, &
                           this%chunk_x, this%nc_fname)
    case ('DOUBLE3D')
      call mem_setptr(dbl3d, idt%mf6varname, mempath)
      call nc_export_array(this%ncid, this%dim_ids, this%var_ids, this%dis, &
                           dbl3d, nc_varname, pkgname, idt%tagname, &
                           idt%shape, idt%longname, input_attr, &
                           this%gridmap_name, this%latlon, this%deflate, &
                           this%shuffle, this%chunk_z, this%chunk_y, &
                           this%chunk_x, iper, iaux, this%nc_fname)
    case default
      ! no-op, no other datatypes exported
    end select
  end subroutine export_input_array

  !> @brief write package gridded input data
  !<
  subroutine export_input_arrays(this, pkgtype, pkgname, mempath, param_dfns)
    use MemoryManagerModule, only: get_isize
    class(DisNCStructuredType), intent(inout) :: this
    character(len=*), intent(in) :: pkgtype
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: mempath
    type(InputParamDefinitionType), dimension(:), pointer, &
      intent(in) :: param_dfns
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: iparam, isize
    do iparam = 1, size(param_dfns)
      ! assign param definition pointer
      idt => param_dfns(iparam)
      ! for now only griddata is exported
      if (idt%blockname == 'GRIDDATA') then
        ! check if variable is already allocated
        call get_isize(idt%mf6varname, mempath, isize)
        if (isize > 0) then
          call this%export_input_array(pkgtype, pkgname, mempath, idt)
        end if
      end if
    end do
  end subroutine export_input_arrays

  !> @brief netcdf export package dynamic input with ilayer index variable
  !<
  subroutine package_step_ilayer(this, export_pkg, ilayer_varname, ilayer)
    use TdisModule, only: kper
    use NCModelExportModule, only: ExportPackageType
    use DefinitionSelectModule, only: get_param_definition_type
    class(DisNCStructuredType), intent(inout) :: this
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
                                  'PERIOD', export_pkg%param_names(iparam), &
                                  this%nc_fname)
      ! set variable name and input attrs
      nc_varname = export_varname(export_pkg%mf6_input%subcomponent_name, idt, &
                                  iper=kper)
      input_attr = this%input_attribute(export_pkg%mf6_input%subcomponent_name, &
                                        idt)
      ! export arrays
      select case (idt%datatype)
      case ('INTEGER1D')
        call mem_setptr(int1d, idt%mf6varname, export_pkg%mf6_input%mempath)
        call nc_export_array(this%ncid, this%dim_ids, this%var_ids, this%dis, &
                             int1d, nc_varname, &
                             export_pkg%mf6_input%subcomponent_name, &
                             idt%tagname, idt%shape, idt%longname, input_attr, &
                             this%gridmap_name, this%latlon, this%deflate, &
                             this%shuffle, this%chunk_z, this%chunk_y, &
                             this%chunk_x, export_pkg%iper, this%nc_fname)
      case ('DOUBLE1D')
        call mem_setptr(dbl1d, idt%mf6varname, export_pkg%mf6_input%mempath)
        call this%export_layer_3d(export_pkg, idt, ilayer_read, ialayer, &
                                  dbl1d, nc_varname, input_attr)
      case ('DOUBLE2D')
        call mem_setptr(dbl2d, idt%mf6varname, export_pkg%mf6_input%mempath)
        nvals = this%dis%ncol * this%dis%nrow
        do n = 1, size(dbl2d, dim=1) ! naux
          dbl1d_ptr(1:nvals) => dbl2d(n, :)
          if (all(dbl1d_ptr == DZERO)) then
          else
            call this%export_layer_3d(export_pkg, idt, ilayer_read, ialayer, &
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
    use TdisModule, only: kper
    use NCModelExportModule, only: ExportPackageType
    use DefinitionSelectModule, only: get_param_definition_type
    class(DisNCStructuredType), intent(inout) :: this
    class(ExportPackageType), pointer, intent(in) :: export_pkg
    integer(I4B), dimension(:), pointer, contiguous :: int1d
    real(DP), dimension(:), pointer, contiguous :: dbl1d
    type(InputParamDefinitionType), pointer :: idt
    character(len=LINELENGTH) :: nc_varname, input_attr
    integer(I4B) :: iparam

    do iparam = 1, export_pkg%nparam
      ! set input definition
      idt => get_param_definition_type(export_pkg%mf6_input%param_dfns, &
                                       export_pkg%mf6_input%component_type, &
                                       export_pkg%mf6_input%subcomponent_type, &
                                       'PERIOD', export_pkg%param_names(iparam), &
                                       this%nc_fname)

      ! set variable name and input attribute string
      nc_varname = export_varname(export_pkg%mf6_input%subcomponent_name, idt, &
                                  iper=kper)
      input_attr = this%input_attribute(export_pkg%mf6_input%subcomponent_name, &
                                        idt)

      ! export arrays
      select case (idt%datatype)
      case ('INTEGER1D')
        call mem_setptr(int1d, export_pkg%param_names(iparam), &
                        export_pkg%mf6_input%mempath)
        call nc_export_array(this%ncid, this%dim_ids, this%var_ids, this%dis, &
                             int1d, nc_varname, &
                             export_pkg%mf6_input%subcomponent_name, &
                             idt%tagname, idt%shape, idt%longname, input_attr, &
                             this%gridmap_name, this%latlon, this%deflate, &
                             this%shuffle, this%chunk_z, this%chunk_y, &
                             this%chunk_x, kper, this%nc_fname)
      case ('DOUBLE1D')
        call mem_setptr(dbl1d, export_pkg%param_names(iparam), &
                        export_pkg%mf6_input%mempath)
        call nc_export_array(this%ncid, this%dim_ids, this%var_ids, this%dis, &
                             dbl1d, nc_varname, &
                             export_pkg%mf6_input%subcomponent_name, &
                             idt%tagname, idt%shape, idt%longname, input_attr, &
                             this%gridmap_name, this%latlon, this%deflate, &
                             this%shuffle, this%chunk_z, this%chunk_y, &
                             this%chunk_x, kper, this%nc_fname)
      case default
        errmsg = 'EXPORT unsupported datatype='//trim(idt%datatype)
        call store_error(errmsg, .true.)
      end select
    end do

    ! synchronize file
    call nf_verify(nf90_sync(this%ncid), this%nc_fname)
  end subroutine package_step

  !> @brief export layer variable as full grid
  !<
  subroutine export_layer_3d(this, export_pkg, idt, ilayer_read, ialayer, &
                             dbl1d, nc_varname, input_attr, iaux)
    use TdisModule, only: kper
    use NCModelExportModule, only: ExportPackageType
    class(DisNCStructuredType), intent(inout) :: this
    class(ExportPackageType), pointer, intent(in) :: export_pkg
    type(InputParamDefinitionType), pointer, intent(in) :: idt
    logical(LGP), intent(in) :: ilayer_read
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: ialayer
    real(DP), dimension(:), pointer, contiguous, intent(in) :: dbl1d
    character(len=*), intent(inout) :: nc_varname
    character(len=*), intent(in) :: input_attr
    integer(I4B), optional, intent(in) :: iaux
    real(DP), dimension(:, :, :), pointer, contiguous :: dbl3d
    integer(I4B) :: n, i, j, k, nvals, idxaux
    real(DP), dimension(:, :), contiguous, pointer :: dbl2d_ptr

    ! initialize
    idxaux = 0
    if (present(iaux)) then
      nc_varname = export_varname(export_pkg%mf6_input%subcomponent_name, &
                                  idt, iper=kper, iaux=iaux)
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

    call nc_export_array(this%ncid, this%dim_ids, this%var_ids, this%dis, dbl3d, &
                         nc_varname, export_pkg%mf6_input%subcomponent_name, &
                         idt%tagname, idt%shape, idt%longname, input_attr, &
                         this%gridmap_name, this%latlon, this%deflate, &
                         this%shuffle, this%chunk_z, this%chunk_y, this%chunk_x, &
                         export_pkg%iper, idxaux, this%nc_fname)

    deallocate (dbl3d)
  end subroutine export_layer_3d

  !> @brief determine packages to write gridded input
  !<
  subroutine add_pkg_data(this)
    use SimVariablesModule, only: idm_context
    use MemoryManagerExtModule, only: mem_set_value
    use MemoryHelperModule, only: create_mem_path
    use IdmDfnSelectorModule, only: param_definitions
    use SourceCommonModule, only: idm_subcomponent_type
    class(DisNCStructuredType), intent(inout) :: this
    character(LENCOMPONENTNAME) :: ptype, pname, pkgtype
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgtypes => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgnames => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mempaths => null()
    type(InputParamDefinitionType), dimension(:), pointer :: param_dfns
    character(len=LENMEMPATH) :: input_mempath, mempath
    integer(I4B) :: n
    integer(I4B), pointer :: export_arrays
    logical(LGP) :: found

    input_mempath = create_mem_path(component=this%modelname, context=idm_context)

    ! set pointers to model path package info
    call mem_setptr(pkgtypes, 'PKGTYPES', input_mempath)
    call mem_setptr(pkgnames, 'PKGNAMES', input_mempath)
    call mem_setptr(mempaths, 'MEMPATHS', input_mempath)

    do n = 1, size(mempaths)
      ! allocate export_arrays
      allocate (export_arrays)
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

      ! cleanup
      deallocate (export_arrays)
    end do
  end subroutine add_pkg_data

  !> @brief create file (group) attributes
  !<
  subroutine add_global_att(this)
    class(DisNCStructuredType), intent(inout) :: this
    ! file scoped title
    call nf_verify(nf90_put_att(this%ncid, NF90_GLOBAL, 'title', &
                                this%annotation%title), this%nc_fname)
    ! source (MODFLOW 6)
    call nf_verify(nf90_put_att(this%ncid, NF90_GLOBAL, 'source', &
                                this%annotation%source), this%nc_fname)
    ! export type (MODFLOW 6)
    call nf_verify(nf90_put_att(this%ncid, NF90_GLOBAL, 'modflow_grid', &
                                this%annotation%grid), this%nc_fname)
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

  !> @brief netcdf export define dimensions
  !<
  subroutine define_dim(this)
    use ConstantsModule, only: MVALIDATE
    use SimVariablesModule, only: isim_mode
    class(DisNCStructuredType), intent(inout) :: this

    ! bound dim
    call nf_verify(nf90_def_dim(this%ncid, 'bnd', 2, this%dim_ids%bnd), &
                   this%nc_fname)

    ! Time
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
      !call nf_verify(nf90_put_att(ncid, var_ids%time, 'bounds', 'time_bnds'), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%time, 'standard_name', &
                                  'time'), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%time, 'long_name', &
                                  'time'), this%nc_fname)
    end if

    ! Z dimension
    call nf_verify(nf90_def_dim(this%ncid, 'z', this%dis%nlay, this%dim_ids%z), &
                   this%nc_fname)
    call nf_verify(nf90_def_var(this%ncid, 'z', NF90_DOUBLE, this%dim_ids%z, &
                                this%var_ids%z), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%z, 'units', 'layer'), &
                   this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%z, 'long_name', &
                                'layer number'), this%nc_fname)
    !call nf_verify(nf90_put_att(this%ncid, this%var_ids%z, 'bounds', 'z_bnds'), &
    !               this%nc_fname)
    !call nf_verify(nf90_def_var(this%ncid, 'z_bnds', NF90_DOUBLE, &
    !                            (/this%dim_ids%bnd, this%dim_ids%z/), &
    !                            this%var_ids%z_bnds), this%nc_fname)
    !call nf_verify(nf90_put_var(this%ncid, this%var_ids%z_bnds, &
    !                            this%elev_bnds), this%nc_fname)

    ! Y dimension
    call nf_verify(nf90_def_dim(this%ncid, 'y', this%dis%nrow, this%dim_ids%y), &
                   this%nc_fname)
    call nf_verify(nf90_def_var(this%ncid, 'y', NF90_DOUBLE, this%dim_ids%y, &
                                this%var_ids%y), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%y, 'units', &
                                this%lenunits), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%y, 'axis', 'Y'), &
                   this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%y, 'standard_name', &
                                'projection_y_coordinate'), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%y, 'long_name', &
                                'Northing'), this%nc_fname)
    if (this%wkt /= '') then
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%y, 'grid_mapping', &
                                  this%gridmap_name), this%nc_fname)
    end if
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%y, 'bounds', 'y_bnds'), &
                   this%nc_fname)
    call nf_verify(nf90_def_var(this%ncid, 'y_bnds', NF90_DOUBLE, &
                                (/this%dim_ids%bnd, this%dim_ids%y/), &
                                this%var_ids%y_bnds), this%nc_fname)

    ! X dimension
    call nf_verify(nf90_def_dim(this%ncid, 'x', this%dis%ncol, this%dim_ids%x), &
                   this%nc_fname)
    call nf_verify(nf90_def_var(this%ncid, 'x', NF90_DOUBLE, this%dim_ids%x, &
                                this%var_ids%x), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%x, 'units', &
                                this%lenunits), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%x, 'axis', 'X'), &
                   this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%x, 'standard_name', &
                                'projection_x_coordinate'), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%x, 'long_name', &
                                'Easting'), this%nc_fname)
    if (this%wkt /= '') then
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%x, 'grid_mapping', &
                                  this%gridmap_name), this%nc_fname)
    end if
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%x, 'bounds', 'x_bnds'), &
                   this%nc_fname)
    call nf_verify(nf90_def_var(this%ncid, 'x_bnds', NF90_DOUBLE, &
                                (/this%dim_ids%bnd, this%dim_ids%x/), &
                                this%var_ids%x_bnds), this%nc_fname)

    ! NCPL dimension
    call nf_verify(nf90_def_dim(this%ncid, 'ncpl', &
                                this%dis%ncol * this%dis%nrow, &
                                this%dim_ids%ncpl), this%nc_fname)
  end subroutine define_dim

  !> @brief create the model layer dependent variables
  !<
  subroutine define_dependent(this)
    use ConstantsModule, only: DHNOFLO
    class(DisNCStructuredType), intent(inout) :: this

    call nf_verify(nf90_def_var(this%ncid, this%xname, NF90_DOUBLE, &
                                (/this%dim_ids%x, this%dim_ids%y, &
                                  this%dim_ids%z, this%dim_ids%time/), &
                                this%var_ids%dependent), &
                   this%nc_fname)

    ! apply chunking parameters
    if (this%chunking_active) then
      call nf_verify(nf90_def_var_chunking(this%ncid, &
                                           this%var_ids%dependent, &
                                           NF90_CHUNKED, &
                                           (/this%chunk_x, this%chunk_y, &
                                             this%chunk_z, this%chunk_time/)), &
                     this%nc_fname)
    end if

    ! deflate and shuffle
    call ncvar_deflate(this%ncid, this%var_ids%dependent, this%deflate, &
                       this%shuffle, this%nc_fname)

    ! put attr
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%dependent, &
                                'units', this%lenunits), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%dependent, &
                                'standard_name', this%annotation%stdname), &
                   this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%dependent, 'long_name', &
                                this%annotation%longname), this%nc_fname)
    call nf_verify(nf90_put_att(this%ncid, this%var_ids%dependent, '_FillValue', &
                                (/DHNOFLO/)), this%nc_fname)

    ! add grid mapping
    call ncvar_gridmap(this%ncid, this%var_ids%dependent, this%gridmap_name, &
                       this%latlon, this%nc_fname)
  end subroutine define_dependent

  !> @brief create the file grid mapping container variable
  !<
  subroutine define_gridmap(this)
    class(DisNCStructuredType), intent(inout) :: this
    integer(I4B) :: var_id
    if (this%wkt /= '') then
      call nf_verify(nf90_redef(this%ncid), this%nc_fname)
      call nf_verify(nf90_def_var(this%ncid, this%gridmap_name, NF90_INT, &
                                  var_id), this%nc_fname)
      ! TODO: consider variants epsg_code, spatial_ref, esri_pe_string, wkt, etc
      call nf_verify(nf90_put_att(this%ncid, var_id, 'crs_wkt', this%wkt), &
                     this%nc_fname)
      call nf_verify(nf90_enddef(this%ncid), this%nc_fname)
      call nf_verify(nf90_put_var(this%ncid, var_id, 1), &
                     this%nc_fname)
    end if
  end subroutine define_gridmap

  !> @brief define grid projection variables
  !<
  subroutine define_geocoords(this)
    class(DisNCStructuredType), intent(inout) :: this
    if (this%latlon) then
      ! lat
      call nf_verify(nf90_def_var(this%ncid, 'lat', NF90_DOUBLE, &
                                  (/this%dim_ids%x, this%dim_ids%y/), &
                                  this%var_ids%latitude), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%latitude, &
                                  'units', 'degrees_north'), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%latitude, &
                                  'standard_name', 'latitude'), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%latitude, &
                                  'long_name', 'latitude'), this%nc_fname)

      ! lon
      call nf_verify(nf90_def_var(this%ncid, 'lon', NF90_DOUBLE, &
                                  (/this%dim_ids%x, this%dim_ids%y/), &
                                  this%var_ids%longitude), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%longitude, &
                                  'units', 'degrees_east'), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%longitude, &
                                  'standard_name', 'longitude'), this%nc_fname)
      call nf_verify(nf90_put_att(this%ncid, this%var_ids%longitude, &
                                  'long_name', 'longitude'), this%nc_fname)
    end if
  end subroutine define_geocoords

  !> @brief add grid projection data
  !<
  subroutine add_proj_data(this)
    class(DisNCStructuredType), intent(inout) :: this
    if (this%latlon) then
      ! lat
      call nf_verify(nf90_put_var(this%ncid, this%var_ids%latitude, &
                                  this%latitude, start=(/1, 1/), &
                                  count=(/this%dis%ncol, this%dis%nrow/)), &
                     this%nc_fname)

      ! lon
      call nf_verify(nf90_put_var(this%ncid, this%var_ids%longitude, &
                                  this%longitude, start=(/1, 1/), &
                                  count=(/this%dis%ncol, this%dis%nrow/)), &
                     this%nc_fname)
    end if
  end subroutine add_proj_data

  !> @brief add grid coordinates
  !<
  subroutine add_grid_data(this)
    class(DisNCStructuredType), intent(inout) :: this
    integer(I4B) :: ibnd, n !, k, i, j
    real(DP), dimension(:, :), pointer, contiguous :: dbl2d
    real(DP), dimension(:), allocatable :: x, y
    real(DP) :: xoff, yoff

    if (this%dis%angrot /= DZERO) then
      xoff = DZERO
      yoff = DZERO
    else
      xoff = this%dis%xorigin
      yoff = this%dis%yorigin
    end if

    allocate (x(size(this%dis%cellx)))
    allocate (y(size(this%dis%celly)))

    do n = 1, size(this%dis%cellx)
      x(n) = this%dis%cellx(n) + xoff
    end do

    do n = 1, size(this%dis%celly)
      y(n) = this%dis%celly(n) + yoff
    end do

    call nf_verify(nf90_put_var(this%ncid, this%var_ids%x, x), &
                   this%nc_fname)
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%y, y), &
                   this%nc_fname)
    ! TODO see cf-conventions 4.3.3. Parametric Vertical Coordinate
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%z, this%layers), &
                   this%nc_fname)

    deallocate (x)
    deallocate (y)

    ! bounds x
    allocate (dbl2d(2, size(this%dis%cellx)))
    ibnd = 1
    do n = 1, size(this%dis%cellx)
      if (ibnd == 1) then
        dbl2d(1, ibnd) = xoff
        dbl2d(2, ibnd) = xoff + this%dis%delr(ibnd)
      else
        dbl2d(1, ibnd) = dbl2d(1, ibnd - 1) + this%dis%delr(ibnd)
        dbl2d(2, ibnd) = dbl2d(2, ibnd - 1) + this%dis%delr(ibnd)
      end if
      ibnd = ibnd + 1
    end do
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%x_bnds, dbl2d), &
                   this%nc_fname)
    deallocate (dbl2d)

    ! bounds y
    allocate (dbl2d(2, size(this%dis%celly)))
    ibnd = 1
    do n = size(this%dis%celly), 1, -1
      if (ibnd == 1) then
        dbl2d(1, ibnd) = yoff + sum(this%dis%delc) - this%dis%delc(n)
        dbl2d(2, ibnd) = yoff + sum(this%dis%delc)
      else
        dbl2d(1, ibnd) = dbl2d(1, ibnd - 1) - this%dis%delc(n)
        dbl2d(2, ibnd) = dbl2d(2, ibnd - 1) - this%dis%delc(n)
      end if
      ibnd = ibnd + 1
    end do
    call nf_verify(nf90_put_var(this%ncid, this%var_ids%y_bnds, dbl2d), &
                   this%nc_fname)
    deallocate (dbl2d)
  end subroutine add_grid_data

  !> @brief define 2d variable chunking
  !<
  subroutine ncvar_chunk2d(ncid, varid, chunk_x, chunk_y, nc_fname)
    integer(I4B), intent(in) :: ncid
    integer(I4B), intent(in) :: varid
    integer(I4B), intent(in) :: chunk_x
    integer(I4B), intent(in) :: chunk_y
    character(len=*), intent(in) :: nc_fname
    if (chunk_y > 0 .and. chunk_x > 0) then
      call nf_verify(nf90_def_var_chunking(ncid, varid, NF90_CHUNKED, &
                                           (/chunk_x, chunk_y/)), nc_fname)
    end if
  end subroutine ncvar_chunk2d

  !> @brief define 3d variable chunking
  !<
  subroutine ncvar_chunk3d(ncid, varid, chunk_x, chunk_y, chunk_z, nc_fname)
    integer(I4B), intent(in) :: ncid
    integer(I4B), intent(in) :: varid
    integer(I4B), intent(in) :: chunk_x
    integer(I4B), intent(in) :: chunk_y
    integer(I4B), intent(in) :: chunk_z
    character(len=*), intent(in) :: nc_fname
    if (chunk_z > 0 .and. chunk_y > 0 .and. chunk_x > 0) then
      call nf_verify(nf90_def_var_chunking(ncid, varid, NF90_CHUNKED, &
                                           (/chunk_x, chunk_y, chunk_z/)), &
                     nc_fname)
    end if
  end subroutine ncvar_chunk3d

  !> @brief define variable compression
  !<
  subroutine ncvar_deflate(ncid, varid, deflate, shuffle, nc_fname)
    integer(I4B), intent(in) :: ncid
    integer(I4B), intent(in) :: varid
    integer(I4B), intent(in) :: deflate
    integer(I4B), intent(in) :: shuffle
    character(len=*), intent(in) :: nc_fname
    ! deflate and shuffle
    if (deflate >= 0) then
      call nf_verify(nf90_def_var_deflate(ncid, varid, shuffle=shuffle, &
                                          deflate=1, deflate_level=deflate), &
                     nc_fname)
    end if
  end subroutine ncvar_deflate

  !> @brief put variable gridmap attributes
  !<
  subroutine ncvar_gridmap(ncid, varid, gridmap_name, latlon, nc_fname)
    integer(I4B), intent(in) :: ncid
    integer(I4B), intent(in) :: varid
    character(len=*), intent(in) :: gridmap_name
    logical(LGP), intent(in) :: latlon
    character(len=*), intent(in) :: nc_fname
    if (gridmap_name /= '') then
      call nf_verify(nf90_put_att(ncid, varid, 'coordinates', 'x y'), &
                     nc_fname)
      call nf_verify(nf90_put_att(ncid, varid, 'grid_mapping', gridmap_name), &
                     nc_fname)
    else if (latlon) then
      call nf_verify(nf90_put_att(ncid, varid, 'coordinates', 'lon lat'), &
                     nc_fname)
    end if
  end subroutine ncvar_gridmap

  !> @brief put variable internal modflow6 attributes
  !<
  subroutine ncvar_mf6attr(ncid, varid, iper, iaux, nc_tag, nc_fname)
    integer(I4B), intent(in) :: ncid
    integer(I4B), intent(in) :: varid
    integer(I4B), intent(in) :: iper
    integer(I4B), intent(in) :: iaux
    character(len=*), intent(in) :: nc_tag
    character(len=*), intent(in) :: nc_fname
    if (nc_tag /= '') then
      call nf_verify(nf90_put_att(ncid, varid, 'modflow6_input', &
                                  nc_tag), nc_fname)
      if (iper > 0) then
        call nf_verify(nf90_put_att(ncid, varid, 'modflow6_iper', &
                                    iper), nc_fname)
      end if
      if (iaux > 0) then
        call nf_verify(nf90_put_att(ncid, varid, 'modflow6_iaux', &
                                    iaux), nc_fname)
      end if
    end if
  end subroutine ncvar_mf6attr

  !> @brief netcdf export 1D integer
  !<
  subroutine nc_export_int1d(ncid, dim_ids, var_ids, dis, p_mem, nc_varname, &
                             pkgname, tagname, shapestr, longname, nc_tag, &
                             gridmap_name, latlon, deflate, shuffle, chunk_z, &
                             chunk_y, chunk_x, iper, nc_fname)
    integer(I4B), intent(in) :: ncid
    type(StructuredNCDimIdType), intent(inout) :: dim_ids
    type(StructuredNCVarIdType), intent(inout) :: var_ids
    type(DisType), pointer, intent(in) :: dis
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: p_mem
    character(len=*), intent(in) :: nc_varname
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: shapestr
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: nc_tag
    character(len=*), intent(in) :: gridmap_name
    logical(LGP), intent(in) :: latlon
    integer(I4B), intent(in) :: deflate
    integer(I4B), intent(in) :: shuffle
    integer(I4B), intent(in) :: chunk_z
    integer(I4B), intent(in) :: chunk_y
    integer(I4B), intent(in) :: chunk_x
    integer(I4B), intent(in) :: iper
    character(len=*), intent(in) :: nc_fname
    integer(I4B) :: var_id, axis_sz
    character(len=LINELENGTH) :: longname_l

    if (shapestr == 'NROW' .or. &
        shapestr == 'NCOL' .or. &
        shapestr == 'NCPL') then

      select case (shapestr)
      case ('NROW')
        axis_sz = dim_ids%y
      case ('NCOL')
        axis_sz = dim_ids%x
      case ('NCPL')
        axis_sz = dim_ids%ncpl
      end select

      longname_l = export_longname(longname, pkgname, tagname, layer=0, iper=iper)

      ! reenter define mode and create variable
      call nf_verify(nf90_redef(ncid), nc_fname)
      call nf_verify(nf90_def_var(ncid, nc_varname, NF90_INT, &
                                  (/axis_sz/), var_id), &
                     nc_fname)

      ! NROW/NCOL shapes use default chunking
      call ncvar_deflate(ncid, var_id, deflate, shuffle, nc_fname)

      ! put attr
      call nf_verify(nf90_put_att(ncid, var_id, '_FillValue', &
                                  (/NF90_FILL_INT/)), nc_fname)
      call nf_verify(nf90_put_att(ncid, var_id, 'long_name', &
                                  longname_l), nc_fname)

      ! add mf6 attr
      call ncvar_mf6attr(ncid, var_id, iper, 0, nc_tag, nc_fname)

      ! exit define mode and write data
      call nf_verify(nf90_enddef(ncid), nc_fname)
      call nf_verify(nf90_put_var(ncid, var_id, p_mem), &
                     nc_fname)

    else
      ! reenter define mode and create variable
      call nf_verify(nf90_redef(ncid), nc_fname)
      call nf_verify(nf90_def_var(ncid, nc_varname, NF90_INT, &
                                  (/dim_ids%x, dim_ids%y, dim_ids%z/), var_id), &
                     nc_fname)

      ! apply chunking parameters
      call ncvar_chunk3d(ncid, var_id, chunk_x, chunk_y, chunk_z, nc_fname)
      ! deflate and shuffle
      call ncvar_deflate(ncid, var_id, deflate, shuffle, nc_fname)

      ! put attr
      call nf_verify(nf90_put_att(ncid, var_id, '_FillValue', &
                                  (/NF90_FILL_INT/)), nc_fname)
      call nf_verify(nf90_put_att(ncid, var_id, 'long_name', &
                                  longname), nc_fname)

      ! add grid mapping and mf6 attr
      call ncvar_gridmap(ncid, var_id, gridmap_name, latlon, nc_fname)
      call ncvar_mf6attr(ncid, var_id, 0, 0, nc_tag, nc_fname)

      ! exit define mode and write data
      call nf_verify(nf90_enddef(ncid), nc_fname)
      call nf_verify(nf90_put_var(ncid, var_id, p_mem, start=(/1, 1, 1/), &
                                  count=(/dis%ncol, dis%nrow, dis%nlay/)), &
                     nc_fname)
    end if
  end subroutine nc_export_int1d

  !> @brief netcdf export 2D integer
  !<
  subroutine nc_export_int2d(ncid, dim_ids, var_ids, dis, p_mem, nc_varname, &
                             pkgname, tagname, shapestr, longname, nc_tag, &
                             gridmap_name, latlon, deflate, shuffle, chunk_z, &
                             chunk_y, chunk_x, nc_fname)
    integer(I4B), intent(in) :: ncid
    type(StructuredNCDimIdType), intent(inout) :: dim_ids
    type(StructuredNCVarIdType), intent(inout) :: var_ids
    type(DisType), pointer, intent(in) :: dis
    integer(I4B), dimension(:, :), pointer, contiguous, intent(in) :: p_mem
    character(len=*), intent(in) :: nc_varname
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: shapestr
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: nc_tag
    character(len=*), intent(in) :: gridmap_name
    logical(LGP), intent(in) :: latlon
    integer(I4B), intent(in) :: deflate
    integer(I4B), intent(in) :: shuffle
    integer(I4B), intent(in) :: chunk_z
    integer(I4B), intent(in) :: chunk_y
    integer(I4B), intent(in) :: chunk_x
    character(len=*), intent(in) :: nc_fname
    integer(I4B) :: var_id

    ! reenter define mode and create variable
    call nf_verify(nf90_redef(ncid), nc_fname)
    call nf_verify(nf90_def_var(ncid, nc_varname, NF90_INT, &
                                (/dim_ids%x, dim_ids%y/), var_id), &
                   nc_fname)

    ! apply chunking parameters
    call ncvar_chunk2d(ncid, var_id, chunk_x, chunk_y, nc_fname)
    ! deflate and shuffle
    call ncvar_deflate(ncid, var_id, deflate, shuffle, nc_fname)

    ! put attr
    call nf_verify(nf90_put_att(ncid, var_id, '_FillValue', &
                                (/NF90_FILL_INT/)), nc_fname)
    call nf_verify(nf90_put_att(ncid, var_id, 'long_name', &
                                longname), nc_fname)

    ! add grid mapping and mf6 attr
    call ncvar_gridmap(ncid, var_id, gridmap_name, latlon, nc_fname)
    call ncvar_mf6attr(ncid, var_id, 0, 0, nc_tag, nc_fname)

    ! exit define mode and write data
    call nf_verify(nf90_enddef(ncid), nc_fname)
    call nf_verify(nf90_put_var(ncid, var_id, p_mem, start=(/1, 1/), &
                                count=(/dis%ncol, dis%nrow/)), &
                   nc_fname)
  end subroutine nc_export_int2d

  !> @brief netcdf export 3D integer
  !<
  subroutine nc_export_int3d(ncid, dim_ids, var_ids, dis, p_mem, nc_varname, &
                             pkgname, tagname, shapestr, longname, nc_tag, &
                             gridmap_name, latlon, deflate, shuffle, chunk_z, &
                             chunk_y, chunk_x, nc_fname)
    integer(I4B), intent(in) :: ncid
    type(StructuredNCDimIdType), intent(inout) :: dim_ids
    type(StructuredNCVarIdType), intent(inout) :: var_ids
    type(DisType), pointer, intent(in) :: dis
    integer(I4B), dimension(:, :, :), pointer, contiguous, intent(in) :: p_mem
    character(len=*), intent(in) :: nc_varname
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: shapestr
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: nc_tag
    character(len=*), intent(in) :: gridmap_name
    logical(LGP), intent(in) :: latlon
    integer(I4B), intent(in) :: deflate
    integer(I4B), intent(in) :: shuffle
    integer(I4B), intent(in) :: chunk_z
    integer(I4B), intent(in) :: chunk_y
    integer(I4B), intent(in) :: chunk_x
    character(len=*), intent(in) :: nc_fname
    integer(I4B) :: var_id

    ! reenter define mode and create variable
    call nf_verify(nf90_redef(ncid), nc_fname)
    call nf_verify(nf90_def_var(ncid, nc_varname, NF90_INT, &
                                (/dim_ids%x, dim_ids%y, dim_ids%z/), var_id), &
                   nc_fname)

    ! apply chunking parameters
    call ncvar_chunk3d(ncid, var_id, chunk_x, chunk_y, chunk_z, nc_fname)
    ! deflate and shuffle
    call ncvar_deflate(ncid, var_id, deflate, shuffle, nc_fname)

    ! put attr
    call nf_verify(nf90_put_att(ncid, var_id, '_FillValue', &
                                (/NF90_FILL_INT/)), nc_fname)
    call nf_verify(nf90_put_att(ncid, var_id, 'long_name', &
                                longname), nc_fname)

    ! add grid mapping and mf6 attr
    call ncvar_gridmap(ncid, var_id, gridmap_name, latlon, nc_fname)
    call ncvar_mf6attr(ncid, var_id, 0, 0, nc_tag, nc_fname)

    ! exit define mode and write data
    call nf_verify(nf90_enddef(ncid), nc_fname)
    call nf_verify(nf90_put_var(ncid, var_id, p_mem, start=(/1, 1, 1/), &
                                count=(/dis%ncol, dis%nrow, dis%nlay/)), &
                   nc_fname)
  end subroutine nc_export_int3d

  !> @brief netcdf export 1D double
  !<
  subroutine nc_export_dbl1d(ncid, dim_ids, var_ids, dis, p_mem, nc_varname, &
                             pkgname, tagname, shapestr, longname, nc_tag, &
                             gridmap_name, latlon, deflate, shuffle, chunk_z, &
                             chunk_y, chunk_x, iper, nc_fname)
    use ConstantsModule, only: DNODATA
    integer(I4B), intent(in) :: ncid
    type(StructuredNCDimIdType), intent(inout) :: dim_ids
    type(StructuredNCVarIdType), intent(inout) :: var_ids
    type(DisType), pointer, intent(in) :: dis
    real(DP), dimension(:), pointer, contiguous, intent(in) :: p_mem
    character(len=*), intent(in) :: nc_varname
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: shapestr
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: nc_tag
    character(len=*), intent(in) :: gridmap_name
    logical(LGP), intent(in) :: latlon
    integer(I4B), intent(in) :: deflate
    integer(I4B), intent(in) :: shuffle
    integer(I4B), intent(in) :: chunk_z
    integer(I4B), intent(in) :: chunk_y
    integer(I4B), intent(in) :: chunk_x
    integer(I4B), intent(in) :: iper
    character(len=*), intent(in) :: nc_fname
    integer(I4B) :: var_id, axis_sz
    real(DP) :: fill_value
    character(len=LINELENGTH) :: longname_l

    if (shapestr == 'NROW' .or. &
        shapestr == 'NCOL' .or. &
        shapestr == 'NCPL') then

      select case (shapestr)
      case ('NROW')
        axis_sz = dim_ids%y
      case ('NCOL')
        axis_sz = dim_ids%x
      case ('NCPL')
        axis_sz = dim_ids%ncpl
      end select

      ! reenter define mode and create variable
      call nf_verify(nf90_redef(ncid), nc_fname)
      call nf_verify(nf90_def_var(ncid, nc_varname, NF90_DOUBLE, &
                                  (/axis_sz/), var_id), &
                     nc_fname)

      ! NROW/NCOL shapes use default chunking
      call ncvar_deflate(ncid, var_id, deflate, shuffle, nc_fname)

      ! put attr
      call nf_verify(nf90_put_att(ncid, var_id, '_FillValue', &
                                  (/NF90_FILL_DOUBLE/)), nc_fname)
      call nf_verify(nf90_put_att(ncid, var_id, 'long_name', &
                                  longname), nc_fname)

      ! add mf6 attr
      call ncvar_mf6attr(ncid, var_id, 0, 0, nc_tag, nc_fname)

      ! exit define mode and write data
      call nf_verify(nf90_enddef(ncid), nc_fname)
      call nf_verify(nf90_put_var(ncid, var_id, p_mem), &
                     nc_fname)

    else
      if (iper > 0) then
        fill_value = DNODATA
      else
        fill_value = NF90_FILL_DOUBLE
      end if

      longname_l = export_longname(longname, pkgname, tagname, layer=0, iper=iper)

      ! reenter define mode and create variable
      call nf_verify(nf90_redef(ncid), nc_fname)
      call nf_verify(nf90_def_var(ncid, nc_varname, NF90_DOUBLE, &
                                  (/dim_ids%x, dim_ids%y, dim_ids%z/), var_id), &
                     nc_fname)

      ! apply chunking parameters
      call ncvar_chunk3d(ncid, var_id, chunk_x, chunk_y, chunk_z, nc_fname)
      ! deflate and shuffle
      call ncvar_deflate(ncid, var_id, deflate, shuffle, nc_fname)

      ! put attr
      call nf_verify(nf90_put_att(ncid, var_id, '_FillValue', &
                                  (/fill_value/)), nc_fname)
      call nf_verify(nf90_put_att(ncid, var_id, 'long_name', &
                                  longname_l), nc_fname)

      ! add grid mapping and mf6 attr
      call ncvar_gridmap(ncid, var_id, gridmap_name, latlon, nc_fname)
      call ncvar_mf6attr(ncid, var_id, iper, 0, nc_tag, nc_fname)

      ! exit define mode and write data
      call nf_verify(nf90_enddef(ncid), nc_fname)
      call nf_verify(nf90_put_var(ncid, var_id, p_mem, start=(/1, 1, 1/), &
                                  count=(/dis%ncol, dis%nrow, dis%nlay/)), &
                     nc_fname)
    end if
  end subroutine nc_export_dbl1d

  !> @brief netcdf export 2D double
  !<
  subroutine nc_export_dbl2d(ncid, dim_ids, var_ids, dis, p_mem, nc_varname, &
                             pkgname, tagname, shapestr, longname, nc_tag, &
                             gridmap_name, latlon, deflate, shuffle, chunk_z, &
                             chunk_y, chunk_x, nc_fname)
    integer(I4B), intent(in) :: ncid
    type(StructuredNCDimIdType), intent(inout) :: dim_ids
    type(StructuredNCVarIdType), intent(inout) :: var_ids
    type(DisType), pointer, intent(in) :: dis
    real(DP), dimension(:, :), pointer, contiguous, intent(in) :: p_mem
    character(len=*), intent(in) :: nc_varname
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: shapestr
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: nc_tag
    character(len=*), intent(in) :: gridmap_name
    logical(LGP), intent(in) :: latlon
    integer(I4B), intent(in) :: deflate
    integer(I4B), intent(in) :: shuffle
    integer(I4B), intent(in) :: chunk_z
    integer(I4B), intent(in) :: chunk_y
    integer(I4B), intent(in) :: chunk_x
    character(len=*), intent(in) :: nc_fname
    integer(I4B) :: var_id

    ! reenter define mode and create variable
    call nf_verify(nf90_redef(ncid), nc_fname)
    call nf_verify(nf90_def_var(ncid, nc_varname, NF90_DOUBLE, &
                                (/dim_ids%x, dim_ids%y/), var_id), &
                   nc_fname)

    ! apply chunking parameters
    call ncvar_chunk2d(ncid, var_id, chunk_x, chunk_y, nc_fname)
    ! deflate and shuffle
    call ncvar_deflate(ncid, var_id, deflate, shuffle, nc_fname)

    ! put attr
    call nf_verify(nf90_put_att(ncid, var_id, '_FillValue', &
                                (/NF90_FILL_DOUBLE/)), nc_fname)
    call nf_verify(nf90_put_att(ncid, var_id, 'long_name', &
                                longname), nc_fname)

    ! add grid mapping and mf6 attr
    call ncvar_gridmap(ncid, var_id, gridmap_name, latlon, nc_fname)
    call ncvar_mf6attr(ncid, var_id, 0, 0, nc_tag, nc_fname)

    ! exit define mode and write data
    call nf_verify(nf90_enddef(ncid), nc_fname)
    call nf_verify(nf90_put_var(ncid, var_id, p_mem, start=(/1, 1/), &
                                count=(/dis%ncol, dis%nrow/)), &
                   nc_fname)
  end subroutine nc_export_dbl2d

  !> @brief netcdf export 3D double
  !<
  subroutine nc_export_dbl3d(ncid, dim_ids, var_ids, dis, p_mem, nc_varname, &
                             pkgname, tagname, shapestr, longname, nc_tag, &
                             gridmap_name, latlon, deflate, shuffle, chunk_z, &
                             chunk_y, chunk_x, iper, iaux, nc_fname)
    use ConstantsModule, only: DNODATA
    integer(I4B), intent(in) :: ncid
    type(StructuredNCDimIdType), intent(inout) :: dim_ids
    type(StructuredNCVarIdType), intent(inout) :: var_ids
    type(DisType), pointer, intent(in) :: dis
    real(DP), dimension(:, :, :), pointer, contiguous, intent(in) :: p_mem
    character(len=*), intent(in) :: nc_varname
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: shapestr
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: nc_tag
    character(len=*), intent(in) :: gridmap_name
    logical(LGP), intent(in) :: latlon
    integer(I4B), intent(in) :: deflate
    integer(I4B), intent(in) :: shuffle
    integer(I4B), intent(in) :: chunk_z
    integer(I4B), intent(in) :: chunk_y
    integer(I4B), intent(in) :: chunk_x
    integer(I4B), intent(in) :: iper
    integer(I4B), intent(in) :: iaux
    character(len=*), intent(in) :: nc_fname
    integer(I4B) :: var_id
    real(DP) :: fill_value
    character(len=LINELENGTH) :: longname_l

    if (iper > 0) then
      fill_value = DNODATA
    else
      fill_value = NF90_FILL_DOUBLE
    end if

    longname_l = export_longname(longname, pkgname, tagname, layer=0, iper=iper)

    ! reenter define mode and create variable
    call nf_verify(nf90_redef(ncid), nc_fname)
    call nf_verify(nf90_def_var(ncid, nc_varname, NF90_DOUBLE, &
                                (/dim_ids%x, dim_ids%y, dim_ids%z/), var_id), &
                   nc_fname)

    ! apply chunking parameters
    call ncvar_chunk3d(ncid, var_id, chunk_x, chunk_y, chunk_z, nc_fname)
    ! deflate and shuffle
    call ncvar_deflate(ncid, var_id, deflate, shuffle, nc_fname)

    ! put attr
    call nf_verify(nf90_put_att(ncid, var_id, '_FillValue', &
                                (/fill_value/)), nc_fname)
    call nf_verify(nf90_put_att(ncid, var_id, 'long_name', &
                                longname_l), nc_fname)

    ! add grid mapping and mf6 attr
    call ncvar_gridmap(ncid, var_id, gridmap_name, latlon, nc_fname)
    call ncvar_mf6attr(ncid, var_id, iper, iaux, nc_tag, nc_fname)

    ! exit define mode and write data
    call nf_verify(nf90_enddef(ncid), nc_fname)
    call nf_verify(nf90_put_var(ncid, var_id, p_mem, start=(/1, 1, 1/), &
                                count=(/dis%ncol, dis%nrow, dis%nlay/)), &
                   nc_fname)
  end subroutine nc_export_dbl3d

  !> @brief build netcdf variable name
  !<
  function export_varname(pkgname, idt, iper, iaux) result(varname)
    use InputOutputModule, only: lowcase
    character(len=*), intent(in) :: pkgname
    type(InputParamDefinitionType), pointer, intent(in) :: idt
    integer(I4B), optional, intent(in) :: iper
    integer(I4B), optional, intent(in) :: iaux
    character(len=LINELENGTH) :: varname
    character(len=LINELENGTH) :: pname, vname
    pname = pkgname
    vname = idt%mf6varname
    call lowcase(pname)
    call lowcase(vname)
    if (present(iper)) then
      if (present(iaux)) then
        write (varname, '(a,i0,a,i0)') trim(pname)//'_'//trim(vname)// &
          '_p', iper, 'a', iaux
      else
        write (varname, '(a,i0)') trim(pname)//'_'//trim(vname)//'_p', iper
      end if
    else
      varname = trim(pname)//'_'//trim(vname)
    end if
  end function export_varname

end module DisNCStructuredModule
