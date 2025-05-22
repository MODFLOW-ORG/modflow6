!> @brief This module contains the NCContextBuildModule
!!
!! Read NetCDF input file and add modflow6 input variables
!1 information to internal NCFileVarsType structure.
!!
!<
module NCContextBuildModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENCOMPONENTNAME
  use SimModule, only: store_error, store_error_filename
  use SimVariablesModule, only: errmsg
  use NCFileVarsModule, only: NCFileVarsType
  use NetCDFCommonModule, only: nf_verify, NETCDF_ATTR_STRLEN
  use netcdf

  implicit none
  private
  public :: open_ncfile
  public :: create_netcdf_context

contains

  !> @brief open netcdf file
  !<
  function open_ncfile(nc_fname, iout) result(ncid)
    use MemoryManagerExtModule, only: mem_set_value
    use NetCDFCommonModule, only: nc_fopen
    character(len=*) :: nc_fname
    integer(I4B) :: iout
    integer(I4B) :: ncid
    logical(LGP) :: exists

    ! initialize
    ncid = 0

    ! check if NETCDF file exists
    inquire (file=nc_fname, exist=exists)
    if (.not. exists) then
      write (errmsg, '(a,a,a)') 'Specified NetCDF input file does &
        &not exist [file=', trim(nc_fname), '].'
      call store_error(errmsg, .true.)
    end if

    ! open
    ncid = nc_fopen(nc_fname, iout)
  end function open_ncfile

  !> @brief add a package input variable to nc_vars structure
  !<
  subroutine add_package_var(modeltype, modelname, nc_vars, input_name, varid, &
                             iout)
    use InputOutputModule, only: lowcase, upcase
    use MemoryHelperModule, only: split_mem_address, split_mem_path
    use SourceCommonModule, only: idm_subcomponent_type
    use SourceCommonModule, only: idm_subcomponent_name
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelname
    type(NCFileVarsType), intent(inout) :: nc_vars
    character(len=*), intent(in) :: input_name
    integer(I4B), intent(in) :: varid
    integer(I4B), intent(in) :: iout
    character(len=NETCDF_ATTR_STRLEN) :: input_str
    character(len=LENCOMPONENTNAME) :: c_name, sc_name
    character(len=LINELENGTH) :: mempath, varname
    integer(I4B) :: layer, period, iaux, mf6_layer, mf6_period, mf6_iaux
    logical(LGP) :: success

    ! initialize
    layer = -1
    period = -1
    iaux = -1
    varname = ''
    c_name = ''
    sc_name = ''

    ! process mf6_input attribute
    if (nf90_get_att(nc_vars%ncid, varid, 'modflow6_input', &
                     input_str) == NF90_NOERR) then
      ! mf6_input should provide a memory address
      call split_mem_address(input_str, mempath, varname, success)

      if (success) then
        ! split the mempath
        call split_mem_path(mempath, c_name, sc_name)
        ! set read tokens to upper case
        call upcase(varname)
        call upcase(c_name)
        call upcase(sc_name)
        ! check for optional layer attribute
        if (nf90_get_att(nc_vars%ncid, varid, &
                         'modflow6_layer', mf6_layer) == NF90_NOERR) then
          layer = mf6_layer
        end if

        ! check for optional period attribute
        if (nf90_get_att(nc_vars%ncid, varid, &
                         'modflow6_iper', mf6_period) == NF90_NOERR) then
          period = mf6_period
        end if

        ! check for optional period attribute
        if (nf90_get_att(nc_vars%ncid, varid, &
                         'modflow6_iaux', mf6_iaux) == NF90_NOERR) then
          iaux = mf6_iaux
        end if

        ! add the variable to netcdf description
        call nc_vars%add(sc_name, varname, layer, period, iaux, varid)
      else
        errmsg = 'NetCDF variable invalid modflow6_input attribute: "'// &
                 trim(input_str)//'".'
        call store_error(errmsg)
        call store_error_filename(nc_vars%nc_fname)
      end if
    end if
  end subroutine add_package_var

  !> @brief verify global attribute modflow_grid is present and return value
  !<
  function verify_global_attr(modeltype, modelname, input_name, nc_fname, ncid) &
    result(grid)
    use InputOutputModule, only: lowcase, upcase
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: input_name
    character(len=*), intent(in) :: nc_fname
    integer(I4B), intent(in) :: ncid
    character(len=NETCDF_ATTR_STRLEN) :: grid

    ! initialize grid
    grid = ''

    ! verify expected mf6_modeltype file attribute
    if (nf90_get_att(ncid, NF90_GLOBAL, "modflow_grid", &
                     grid) == NF90_NOERR) then
      ! set grid to upper case
      call upcase(grid)
    else
      errmsg = 'NetCDF input file global attribute "modflow_grid" not found.'
      call store_error(errmsg)
      call store_error_filename(nc_fname)
    end if
  end function verify_global_attr

  !> @brief create internal description of modflow6 input variables in netcdf file
  !<
  subroutine create_netcdf_context(modeltype, modelname, input_name, &
                                   nc_vars, nc_fname, ncid, iout)
    use InputOutputModule, only: lowcase, upcase
    character(len=*), intent(in) :: modeltype
    character(len=*), intent(in) :: modelname
    character(len=*), intent(in) :: input_name
    type(NCFileVarsType), intent(inout) :: nc_vars
    character(len=*), intent(in) :: nc_fname
    integer(I4B), intent(in) :: ncid
    integer(I4B), intent(in) :: iout
    integer(I4B) :: ndim, nvar, nattr, unlimDimID
    integer(I4B), dimension(:), allocatable :: varids
    character(len=LINELENGTH) :: grid
    integer(I4B) :: iparam

    ! check global attributes
    grid = verify_global_attr(modeltype, modelname, input_name, nc_fname, ncid)

    ! initialize netcdf input structure
    call nc_vars%init(modelname, nc_fname, ncid, grid)

    ! inquire for root dataset info
    call nf_verify(nf90_inquire(ncid, ndim, nvar, nattr, unlimdimid), &
                   nc_vars%nc_fname)

    ! allocate and set varids
    allocate (varids(nvar))
    call nf_verify(nf90_inq_varids(ncid, nvar, varids), nc_vars%nc_fname)
    do iparam = 1, nvar
      ! validate and add netcdf file input variable
      call add_package_var(modeltype, modelname, nc_vars, input_name, &
                           varids(iparam), iout)
    end do

    ! cleanup
    deallocate (varids)
  end subroutine create_netcdf_context

end module NCContextBuildModule
