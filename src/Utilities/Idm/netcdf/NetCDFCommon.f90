!> @brief This module contains the NetCDFCommonModule
!!
!! Common NetCDF interfaces and constants
!!
!<
module NetCDFCommonModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use netcdf

  implicit none
  private
  public :: nc_fopen, nc_fclose
  public :: NETCDF_MAX_DIM
  public :: NETCDF_ATTR_STRLEN
  public :: nf_verify

  integer(I4B), parameter :: NETCDF_MAX_DIM = 6
  integer(I4B), parameter :: NETCDF_ATTR_STRLEN = 80

contains

  !> @brief Open netcdf file
  !<
  function nc_fopen(nc_fname, iout) result(ncid)
    character(len=*), intent(in) :: nc_fname
    integer(I4B), intent(in) :: iout
    integer(I4B) :: ncid
    ! initialize
    ncid = -1
    ! open netcdf file
    call nf_verify(nf90_open(nc_fname, NF90_NOWRITE, ncid), nc_fname)
  end function nc_fopen

  !> @brief Close netcdf file
  !<
  subroutine nc_fclose(ncid, nc_fname)
    integer(I4B), intent(in) :: ncid
    character(len=*), intent(in) :: nc_fname
    ! close netcdf file
    call nf_verify(nf90_close(ncid), nc_fname)
  end subroutine nc_fclose

  !> @brief error check a netcdf-fortran interface call
  !<
  subroutine nf_verify(res, nc_fname)
    integer(I4B), intent(in) :: res
    character(len=*), intent(in) :: nc_fname
    character(len=LINELENGTH) :: errstr

    ! strings are set for a subset of errors
    ! but the exit status will always be reported
    if (res /= NF90_NOERR) then
      !
      select case (res)
      case (-33) ! (NC_EBADID)
        errstr = 'Not a netcdf id'
      case (NF90_EINVAL) ! (-36)
        errstr = 'Invalid Argument'
      case (NF90_EPERM) ! (-37)
        errstr = 'Write to read only'
      case (-38) ! (NC_ENOTINDEFINE)
        errstr = 'Operation not allowed in data mode'
      case (-39) ! (NC_EINDEFINE)
        errstr = 'Operation not allowed in define mode'
      case (NF90_EINVALCOORDS) ! (-40)
        errstr = 'Index exceeds dimension bound'
      case (NF90_ENAMEINUSE) ! (-42)
        errstr = 'String match to name in use'
      case (NF90_ENOTATT) ! (-43)
        errstr = 'Attribute not found'
      case (-45) ! (NC_EBADTYPE)
        errstr = 'Not a netcdf data type'
      case (NF90_EBADDIM) ! (-46)
        errstr = 'Invalid dimension id or name'
      case (NF90_ENOTVAR) ! (-49)
        errstr = 'Variable not found'
      case (NF90_ENOTNC) ! (-51)
        errstr = 'Not a netcdf file'
      case (NF90_ECHAR) ! (-56)
        errstr = 'Attempt to convert between text & numbers'
      case (NF90_EEDGE) ! (-57)
        errstr = 'Edge+start exceeds dimension bound'
      case (NF90_ESTRIDE) ! (-58)
        errstr = 'Illegal stride'
      case (NF90_EBADNAME) ! (-59)
        errstr = 'Attribute or variable name contains illegal characters'
      case (-127) ! (NC_EBADCHUNK)
        errstr = 'Bad chunksize.'
      case default
        errstr = ''
      end select

      if (errstr /= '') then
        write (errmsg, '(a,a,a,i0,a)') 'NetCDF library error [error="', &
          trim(errstr), '", exit code=', res, '].'
      else
        write (errmsg, '(a,i0,a)') 'NetCDF library error [exit code=', &
          res, '].'
      end if

      call store_error(errmsg)
      call store_error_filename(nc_fname)
    end if
  end subroutine nf_verify

end module NetCDFCommonModule
