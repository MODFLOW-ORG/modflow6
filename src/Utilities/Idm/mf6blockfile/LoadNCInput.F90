!> @brief This module contains the LoadNCInputModule
!!
!! Access interfaces for netcdf reads optioned from a
!! traditional modflow 6 input. These interfaces throw
!! an error when netcdf libraries aren't compiled in.
!<
module LoadNCInputModule

  use KindModule, only: DP, I4B, LGP
  use ModflowInputModule, only: ModflowInputType
  use InputDefinitionModule, only: InputParamDefinitionType
  use NCFileVarsModule, only: NCPackageVarsType
#if defined(__WITH_NETCDF__)
  use NCArrayReaderModule, only: netcdf_array_load
#endif

  implicit none
  private
  public :: netcdf_read_array

  interface netcdf_read_array
    module procedure nc_read_int1d, nc_read_int2d, &
      nc_read_int3d, nc_read_dbl1d, nc_read_dbl2d, &
      nc_read_dbl3d
  end interface netcdf_read_array

contains

#if !defined(__WITH_NETCDF__)
  !> @brief Set an error and exit if NetCDF libraries aren't available
  !<
  subroutine error_and_exit(tagname, input_fname)
    use SimVariablesModule, only: errmsg
    use SimModule, only: store_error, store_error_filename
    character(len=*), intent(in) :: tagname
    character(len=*), intent(in) :: input_fname
    write (errmsg, '(a)') &
      'Cannot load NetCDF array for variable "'//trim(tagname)// &
      '". NetCDF libraries are not available. Ensure Extended MODFLOW 6 in use.'
    call store_error(errmsg)
    call store_error_filename(input_fname)
  end subroutine error_and_exit
#endif

  !> @brief Read a NetCDF integer 1d array
  !<
  subroutine nc_read_int1d(int1d, mshape, idt, mf6_input, nc_vars, input_fname, &
                           iout, kper)
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: int1d
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    type(ModflowInputType), intent(in) :: mf6_input
    type(NCPackageVarsType), pointer, intent(in) :: nc_vars
    character(len=*), intent(in) :: input_fname
    integer(I4B), intent(in) :: iout
    integer(I4B), optional, intent(in) :: kper
#if defined(__WITH_NETCDF__)
    call netcdf_array_load(int1d, mshape, idt, mf6_input, nc_vars, input_fname, &
                           iout, kper)
#else
    call error_and_exit(idt%tagname, input_fname)
#endif
  end subroutine nc_read_int1d

  !> @brief Read a NetCDF integer 2d array
  !<
  subroutine nc_read_int2d(int2d, mshape, idt, mf6_input, nc_vars, input_fname, &
                           iout)
    integer(I4B), dimension(:, :), pointer, contiguous, intent(in) :: int2d
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    type(ModflowInputType), intent(in) :: mf6_input
    type(NCPackageVarsType), pointer, intent(in) :: nc_vars
    character(len=*), intent(in) :: input_fname
    integer(I4B), intent(in) :: iout
#if defined(__WITH_NETCDF__)
    call netcdf_array_load(int2d, mshape, idt, mf6_input, nc_vars, input_fname, &
                           iout)
#else
    call error_and_exit(idt%tagname, input_fname)
#endif
  end subroutine nc_read_int2d

  !> @brief Read a NetCDF integer 3d array
  !<
  subroutine nc_read_int3d(int3d, mshape, idt, mf6_input, nc_vars, input_fname, &
                           iout)
    integer(I4B), dimension(:, :, :), pointer, contiguous, intent(in) :: int3d
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    type(ModflowInputType), intent(in) :: mf6_input
    type(NCPackageVarsType), pointer, intent(in) :: nc_vars
    character(len=*), intent(in) :: input_fname
    integer(I4B), intent(in) :: iout
#if defined(__WITH_NETCDF__)
    call netcdf_array_load(int3d, mshape, idt, mf6_input, nc_vars, input_fname, &
                           iout)
#else
    call error_and_exit(idt%tagname, input_fname)
#endif
  end subroutine nc_read_int3d

  !> @brief Read a NetCDF double 1d array
  !<
  subroutine nc_read_dbl1d(dbl1d, mshape, idt, mf6_input, nc_vars, input_fname, &
                           iout, kper, iaux)
    real(DP), dimension(:), pointer, contiguous, intent(in) :: dbl1d
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    type(ModflowInputType), intent(in) :: mf6_input
    type(NCPackageVarsType), pointer, intent(in) :: nc_vars
    character(len=*), intent(in) :: input_fname
    integer(I4B), intent(in) :: iout
    integer(I4B), optional, intent(in) :: kper
    integer(I4B), optional, intent(in) :: iaux
#if defined(__WITH_NETCDF__)
    call netcdf_array_load(dbl1d, mshape, idt, mf6_input, nc_vars, input_fname, &
                           iout, kper=kper, iaux=iaux)
#else
    call error_and_exit(idt%tagname, input_fname)
#endif
  end subroutine nc_read_dbl1d

  !> @brief Read a NetCDF double 2d array
  !<
  subroutine nc_read_dbl2d(dbl2d, mshape, idt, mf6_input, nc_vars, input_fname, &
                           iout)
    real(DP), dimension(:, :), pointer, contiguous, intent(in) :: dbl2d
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    type(ModflowInputType), intent(in) :: mf6_input
    type(NCPackageVarsType), pointer, intent(in) :: nc_vars
    character(len=*), intent(in) :: input_fname
    integer(I4B), intent(in) :: iout
#if defined(__WITH_NETCDF__)
    call netcdf_array_load(dbl2d, mshape, idt, mf6_input, nc_vars, input_fname, &
                           iout)
#else
    call error_and_exit(idt%tagname, input_fname)
#endif
  end subroutine nc_read_dbl2d

  !> @brief Read a NetCDF double 3d array
  !<
  subroutine nc_read_dbl3d(dbl3d, mshape, idt, mf6_input, nc_vars, input_fname, &
                           iout)
    real(DP), dimension(:, :, :), pointer, contiguous, intent(in) :: dbl3d
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: mshape !< model shape
    type(InputParamDefinitionType), intent(in) :: idt !< input data type object describing this record
    type(ModflowInputType), intent(in) :: mf6_input
    type(NCPackageVarsType), pointer, intent(in) :: nc_vars
    character(len=*), intent(in) :: input_fname
    integer(I4B), intent(in) :: iout
#if defined(__WITH_NETCDF__)
    call netcdf_array_load(dbl3d, mshape, idt, mf6_input, nc_vars, input_fname, &
                           iout)
#else
    call error_and_exit(idt%tagname, input_fname)
#endif
  end subroutine nc_read_dbl3d

end module LoadNCInputModule
