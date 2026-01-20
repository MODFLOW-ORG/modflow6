!> @brief This module contains the MODFLOW 6 BMI
!!
!! This BMI interface matches the CSDMS standard, with a few modifications:
!!
!! - This interface will build into a shared library that can be called from other
!!   executables and scripts, not necessarily written in Fortran. Therefore we have
!!   omitted the type-boundness of the routines, since they cannot have the
!!   bind(C,"...") attribute.
!! - MODFLOW has internal data arrays with rank > 1 that we would like to expose. An
!!   example would be access to data in the BOUND array of GWF boundary packages (BndType).
!!   The get_value_ptr calls below support this, returning a C-style pointer to the arrays
!!   and methods have been added to query the variable's rank and shape.
!!
!! Note on style: BMI apparently uses underscores, we use underscores in some
!! places but camelcase in other. Since this is a dedicated BMI interface module,
!! we'll use underscores here as well.
!<
module mf6bmi
  use mf6bmiUtil
  use mf6bmiError
  use Mf6CoreModule
  use TdisModule, only: kper, kstp
  use iso_c_binding, only: c_int, c_char, c_double, C_NULL_CHAR, c_loc, c_ptr, &
                           c_f_pointer
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENMEMPATH, LENVARNAME
  use CharacterStringModule
  use MemoryManagerModule, only: mem_setptr, get_mem_elem_size, get_isize, &
                                 get_mem_rank, get_mem_shape, get_mem_type, &
                                 memorystore
  use MemoryContainerIteratorModule, only: MemoryContainerIteratorType
  use MemoryTypeModule, only: MemoryType
  use MemoryHelperModule, only: create_mem_address
  use SimVariablesModule, only: simstdout, istdout
  use InputOutputModule, only: getunit
  implicit none

  integer(c_int), bind(C, name="ISTDOUTTOFILE") :: istdout_to_file = 1 !< output control: =0 to screen, >0 to file
  !DIR$ ATTRIBUTES DLLEXPORT :: istdout_to_file

contains

  function bmi_get_component_name(name) result(bmi_status) &
    bind(C, name="get_component_name")
    !DIR$ ATTRIBUTES DLLEXPORT :: bmi_get_component_name
    ! -- dummy variables
    character(kind=c_char), intent(inout) :: name(BMI_LENCOMPONENTNAME)
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables

    name = string_to_char_array('MODFLOW 6', 9)
    bmi_status = BMI_SUCCESS

  end function bmi_get_component_name

  !> @brief Initialize the computational core
  !!
  !! It is required to have the MODFLOW 6 configuration file 'mfsim.nam'
  !! available in the current working directory when calling this routine.
  !!
  !! NOTE: initialization should always be matched with a call to
  !! bmi_finalize(). However, we currently do not support the reinitialization
  !! of a model in the same memory space... You would have to create a new
  !! process for that.
  !<
  function bmi_initialize() result(bmi_status) bind(C, name="initialize")
    !DIR$ ATTRIBUTES DLLEXPORT :: bmi_initialize
    ! -- dummy variables
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables

    if (istdout_to_file > 0) then
      ! -- open stdout file mfsim.stdout
      istdout = getunit()
      !
      ! -- set STDOUT to a physical file unit
      open (unit=istdout, file=simstdout)
    end if
    !
    ! -- initialize MODFLOW 6
    call Mf6Initialize()

    bmi_status = BMI_SUCCESS

  end function bmi_initialize

  !> @brief Perform a computational time step
  !!
  !! It will prepare the timestep, call the calculation routine
  !! on all the solution groups in the simulation, and finalize
  !! the timestep by printing out diagnostics and writing output.
  !! It can be called in succession to perform multiple steps up
  !! to the simulation's end time is reached.
  !<
  function bmi_update() result(bmi_status) bind(C, name="update")
    !DIR$ ATTRIBUTES DLLEXPORT :: bmi_update
    ! -- dummy variables
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    logical :: hasConverged

    hasConverged = Mf6Update()

    bmi_status = BMI_SUCCESS
  end function bmi_update

  !> @brief Clean up the initialized simulation
  !!
  !! Performs teardown tasks for the initialized simulation, this
  !! call should match the call to bmi_initialize()
  !<
  function bmi_finalize() result(bmi_status) bind(C, name="finalize")
    !DIR$ ATTRIBUTES DLLEXPORT :: bmi_finalize
    ! -- modules
    use SimVariablesModule, only: iforcestop
    ! -- dummy variables
    integer(kind=c_int) :: bmi_status !< BMI status code

    ! we don't want a full stop() here, this disables it:
    iforcestop = 0
    call Mf6Finalize()

    bmi_status = BMI_SUCCESS

  end function bmi_finalize

  !> @brief Get the start time of the simulation
  !!
  !! As MODFLOW currently does not have internal time, this will be
  !! returning 0.0 for now. New version...
  !<
  function get_start_time(start_time) result(bmi_status) &
    bind(C, name="get_start_time")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_start_time
    ! -- dummy variables
    real(kind=c_double), intent(out) :: start_time !< start time
    integer(kind=c_int) :: bmi_status !< BMI status code

    start_time = 0.0_DP
    bmi_status = BMI_SUCCESS

  end function get_start_time

  !> @brief Get the end time of the simulation
  !!
  !! As MODFLOW does currently does not have internal time, this will be
  !! equal to the total runtime.
  !<
  function get_end_time(end_time) result(bmi_status) bind(C, name="get_end_time")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_end_time
    ! -- modules
    use TdisModule, only: totalsimtime
    ! -- dummy variables
    real(kind=c_double), intent(out) :: end_time !< end time
    integer(kind=c_int) :: bmi_status !< BMI status code

    end_time = totalsimtime
    bmi_status = BMI_SUCCESS

  end function get_end_time

  !> @brief Get the current time of the simulation
  !!
  !! As MODFLOW currently does not have internal time, this will be
  !! equal to the time passed w.r.t. the start time of the simulation.
  !<
  function get_current_time(current_time) result(bmi_status) &
    bind(C, name="get_current_time")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_current_time
    ! -- modules
    use TdisModule, only: totim
    ! -- dummy variables
    real(kind=c_double), intent(out) :: current_time !< current time
    integer(kind=c_int) :: bmi_status !< BMI status code

    current_time = totim
    bmi_status = BMI_SUCCESS

  end function get_current_time

  !> @brief Get the time step for the simulation
  !!
  !! Note that the returned value may vary between and within stress periods,
  !! depending on your time discretization settings in the TDIS package.
  !<
  function get_time_step(time_step) result(bmi_status) &
    bind(C, name="get_time_step")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_time_step
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy variables
    real(kind=c_double), intent(out) :: time_step !< current time step
    integer(kind=c_int) :: bmi_status !< BMI status code

    time_step = delt
    bmi_status = BMI_SUCCESS

  end function get_time_step

  !> @brief Get the number of input variables in the simulation
  !!
  !! This concerns all variables stored in the memory manager
  !<
  function get_input_item_count(count) result(bmi_status) &
    bind(C, name="get_input_item_count")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_input_item_count
    ! -- dummy variables
    integer(kind=c_int), intent(out) :: count !< the number of input variables
    integer(kind=c_int) :: bmi_status !< BMI status code

    count = memorystore%count()

    bmi_status = BMI_SUCCESS

  end function get_input_item_count

  !> @brief Get the number of output variables in the simulation
  !!
  !! This concerns all variables stored in the memory manager
  !<
  function get_output_item_count(count) result(bmi_status) &
    bind(C, name="get_output_item_count")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_output_item_count
    ! -- dummy variables
    integer(kind=c_int), intent(out) :: count !< the number of output variables
    integer(kind=c_int) :: bmi_status !< BMI status code

    count = memorystore%count()

    bmi_status = BMI_SUCCESS

  end function get_output_item_count

  !> @brief Returns all input variables in the simulation
  !!
  !! This functions returns the full address for all variables in the
  !! memory manager
  !!
  !! The array @p c_names should be pre-allocated of proper size:
  !!
  !! size = BMI_LENVARADDRESS *  get_input_item_count()
  !!
  !! The strings will be written contiguously with stride equal to
  !! BMI_LENVARADDRESS and nul-terminated where the trailing spaces start:
  !!
  !! c_names = 'variable_address_1\x00 ... variable_address_2\x00 ... ' etc.
  !<
  function get_input_var_names(c_names) result(bmi_status) &
    bind(C, name="get_input_var_names")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_input_var_names
    ! -- dummy variables
    character(kind=c_char, len=1), intent(inout) :: c_names(*) !< array with memory paths for input variables
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    integer(I4B) :: start, i
    type(MemoryContainerIteratorType), allocatable :: itr
    type(MemoryType), pointer :: mt => null()
    character(len=LENMEMADDRESS) :: var_address

    start = 1
    itr = memorystore%iterator()
    do while (itr%has_next())
      call itr%next()
      mt => itr%value()
      var_address = create_mem_address(mt%path, mt%name)
      do i = 1, len(trim(var_address))
        c_names(start + i - 1) = var_address(i:i)
      end do
      c_names(start + i) = c_null_char
      start = start + BMI_LENVARADDRESS
    end do

    bmi_status = BMI_SUCCESS

  end function get_input_var_names

  !> @brief Returns all output variables in the simulation
  !!
  !! This function works analogously to get_input_var_names(),
  !! and currently returns the same set of memory variables,
  !! which is all of them!
  !<
  function get_output_var_names(c_names) result(bmi_status) &
    bind(C, name="get_output_var_names")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_output_var_names
    ! -- dummy variables
    character(kind=c_char, len=1), intent(inout) :: c_names(*) !< array with memory paths for output variables
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    integer(I4B) :: start, i
    type(MemoryContainerIteratorType), allocatable :: itr
    type(MemoryType), pointer :: mt => null()
    character(len=LENMEMADDRESS) :: var_address

    start = 1
    itr = memorystore%iterator()
    do while (itr%has_next())
      call itr%next()
      mt => itr%value()
      var_address = create_mem_address(mt%path, mt%name)
      do i = 1, len(trim(var_address))
        c_names(start + i - 1) = var_address(i:i)
      end do
      c_names(start + i) = c_null_char
      start = start + BMI_LENVARADDRESS
    end do

    bmi_status = BMI_SUCCESS

  end function get_output_var_names

  !> @brief Get the size (in bytes) of a single element of a variable
  !<
  function get_var_itemsize(c_var_address, var_size) result(bmi_status) &
    bind(C, name="get_var_itemsize")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_var_itemsize
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    integer(kind=c_int), intent(out) :: var_size !< size of the element in bytes
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    call get_mem_elem_size(var_name, mem_path, var_size)
    if (var_size == -1) bmi_status = BMI_FAILURE

  end function get_var_itemsize

  !> @brief Get size of the variable, in bytes
  !<
  function get_var_nbytes(c_var_address, var_nbytes) result(bmi_status) &
    bind(C, name="get_var_nbytes")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_var_nbytes
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    integer(kind=c_int), intent(out) :: var_nbytes !< size in bytes
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    integer(I4B) :: var_size, isize
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    call get_mem_elem_size(var_name, mem_path, var_size)
    if (var_size == -1) bmi_status = BMI_FAILURE
    call get_isize(var_name, mem_path, isize)
    if (isize == -1) bmi_status = BMI_FAILURE

    var_nbytes = var_size * isize

  end function get_var_nbytes

  !> @brief Copy the double precision values of a variable into the array
  !!
  !! The copied variable us located at @p c_var_address. The caller should
  !! provide @p c_arr_ptr pointing to an array of the proper shape (the
  !! BMI function get_var_shape() can be used to create it). Multi-dimensional
  !! arrays are supported.
  !<
  function get_value_double(c_var_address, c_arr_ptr) result(bmi_status) &
    bind(C, name="get_value_double")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_value_double
    ! -- modules
    use MemorySetHandlerModule, only: on_memory_set
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    type(c_ptr), intent(in) :: c_arr_ptr !< pointer to the double precision array
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid
    integer(I4B) :: rank
    real(DP), pointer :: src_ptr, tgt_ptr
    real(DP), dimension(:), pointer, contiguous :: src1D_ptr, tgt1D_ptr
    real(DP), dimension(:, :), pointer, contiguous :: src2D_ptr, tgt2D_ptr
    real(DP), dimension(:, :, :), pointer, contiguous :: src3D_ptr, tgt3D_ptr
    integer(I4B) :: i, j, k

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    ! convert pointer and copy data from memory manager into
    ! the passed array, using loops to avoid stack overflow
    rank = -1
    call get_mem_rank(var_name, mem_path, rank)

    if (rank == 0) then
      call mem_setptr(src_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, tgt_ptr)
      tgt_ptr = src_ptr
    else if (rank == 1) then
      call mem_setptr(src1D_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, tgt1D_ptr, shape(src1D_ptr))
      do i = 1, size(tgt1D_ptr)
        tgt1D_ptr(i) = src1D_ptr(i)
      end do
    else if (rank == 2) then
      call mem_setptr(src2D_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, tgt2D_ptr, shape(src2D_ptr))
      do j = 1, size(tgt2D_ptr, 2)
        do i = 1, size(tgt2D_ptr, 1)
          tgt2D_ptr(i, j) = src2D_ptr(i, j)
        end do
      end do
    else if (rank == 3) then
      call mem_setptr(src3D_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, tgt3D_ptr, shape(src3D_ptr))
      do k = 1, size(tgt3D_ptr, 3)
        do j = 1, size(tgt3D_ptr, 2)
          do i = 1, size(tgt3D_ptr, 1)
            tgt3D_ptr(i, j, k) = src3D_ptr(i, j, k)
          end do
        end do
      end do
    else
      write (bmi_last_error, fmt_unsupported_rank) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

  end function get_value_double

  !> @brief Copy the integer values of a variable into the array
  !!
  !! The copied variable is located at @p c_var_address. The caller should
  !! provide @p c_arr_ptr pointing to an array of the proper shape (the
  !! BMI function get_var_shape() can be used to create it). Multi-dimensional
  !! arrays are supported.
  !<
  function get_value_int(c_var_address, c_arr_ptr) result(bmi_status) &
    bind(C, name="get_value_int")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_value_int
    ! -- modules
    use MemorySetHandlerModule, only: on_memory_set
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    type(c_ptr), intent(in) :: c_arr_ptr !< pointer to the integer array
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid
    integer(I4B) :: rank
    integer(I4B), pointer :: src_ptr, tgt_ptr
    integer(I4B), dimension(:), pointer, contiguous :: src1D_ptr, tgt1D_ptr
    integer(I4B), dimension(:, :), pointer, contiguous :: src2D_ptr, tgt2D_ptr
    integer(I4B), dimension(:, :, :), pointer, contiguous :: src3D_ptr, tgt3D_ptr
    integer(I4B) :: i, j, k

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    ! convert pointer and copy data from memory manager into
    ! the passed array, using loops to avoid stack overflow
    rank = -1
    call get_mem_rank(var_name, mem_path, rank)

    if (rank == 0) then
      call mem_setptr(src_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, tgt_ptr)
      tgt_ptr = src_ptr
    else if (rank == 1) then
      call mem_setptr(src1D_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, tgt1D_ptr, shape(src1D_ptr))
      do i = 1, size(tgt1D_ptr)
        tgt1D_ptr(i) = src1D_ptr(i)
      end do
    else if (rank == 2) then
      call mem_setptr(src2D_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, tgt2D_ptr, shape(src2D_ptr))
      do j = 1, size(tgt2D_ptr, 2)
        do i = 1, size(tgt2D_ptr, 1)
          tgt2D_ptr(i, j) = src2D_ptr(i, j)
        end do
      end do
    else if (rank == 3) then
      call mem_setptr(src3D_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, tgt3D_ptr, shape(src3D_ptr))
      do k = 1, size(tgt3D_ptr, 3)
        do j = 1, size(tgt3D_ptr, 2)
          do i = 1, size(tgt3D_ptr, 1)
            tgt3D_ptr(i, j, k) = src3D_ptr(i, j, k)
          end do
        end do
      end do
    else
      write (bmi_last_error, fmt_unsupported_rank) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

  end function get_value_int

  !> @brief Copy the logical scalar value into the array
  !!
  !! The copied variable us located at @p c_var_address. The caller should
  !! provide @p c_arr_ptr pointing to a scalar array with rank=0.
  !<
  function get_value_bool(c_var_address, c_arr_ptr) result(bmi_status) &
    bind(C, name="get_value_bool")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_value_bool
    ! -- modules
    use MemorySetHandlerModule, only: on_memory_set
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    type(c_ptr), intent(in) :: c_arr_ptr !< pointer to the logical array
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid
    integer(I4B) :: rank
    logical(LGP), pointer :: src_ptr, tgt_ptr

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    rank = -1
    call get_mem_rank(var_name, mem_path, rank)

    ! convert pointer
    if (rank == 0) then
      call mem_setptr(src_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, tgt_ptr)
      tgt_ptr = src_ptr
    else
      write (bmi_last_error, fmt_unsupported_rank) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

  end function get_value_bool

  !> @brief Copy the string(s) of a variable into the array
  !!
  !! The copied variable is located at @p c_var_address. The caller should
  !! provide @p c_arr_ptr pointing to an array of the proper shape (the
  !! BMI function get_var_shape() can be used to create it). For strings
  !! currently scalars and 1d arrays (of CharacterStringType) are
  !< supported
  function get_value_string(c_var_address, c_arr_ptr) result(bmi_status) &
    bind(C, name="get_value_string")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_value_string
    ! -- modules
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    type(c_ptr), intent(in) :: c_arr_ptr !< pointer to the string array
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid
    integer(I4B) :: rank
    character(len=:), pointer :: srcstr
    character(kind=c_char), pointer :: tgtstr(:)
    type(CharacterStringType), dimension(:), pointer, contiguous :: srccharstr1d
    character(kind=c_char), pointer :: tgtstr1d(:, :)
    character(:), allocatable :: tempstr
    integer(I4B) :: i, ilen, isize

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    ! single string, or array of strings (CharacterStringType)
    rank = -1
    call get_mem_rank(var_name, mem_path, rank)

    if (rank == 0) then
      ! a string scalar
      call mem_setptr(srcstr, var_name, mem_path)
      call get_mem_elem_size(var_name, mem_path, ilen)
      call c_f_pointer(c_arr_ptr, tgtstr, shape=[ilen + 1])

      tgtstr(1:len(srcstr) + 1) = string_to_char_array(srcstr, len(srcstr))

    else if (rank == 1) then
      ! an array of strings
      call mem_setptr(srccharstr1d, var_name, mem_path)
      if (.not. associated(srccharstr1d)) then
        write (bmi_last_error, fmt_general_err) 'string type not supported in API'
        call report_bmi_error(bmi_last_error)
        bmi_status = BMI_FAILURE
        return
      end if

      ! create fortran pointer to C data array
      call get_isize(var_name, mem_path, isize)
      call get_mem_elem_size(var_name, mem_path, ilen)
      call c_f_pointer(c_arr_ptr, tgtstr1d, shape=[ilen + 1, isize])

      ! allocate work array to handle CharacterStringType,
      ! and copy the strings
      allocate (character(ilen) :: tempstr)
      do i = 1, isize
        tempstr = srccharstr1d(i)
        tgtstr1d(1:ilen + 1, i) = string_to_char_array(tempstr, ilen)
      end do
      deallocate (tempstr)
    else
      write (bmi_last_error, fmt_unsupported_rank) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

  end function get_value_string

  !> @brief Copy the value of a variable into the array
  !!
  !! The copied variable is located at @p c_var_address. The caller should
  !! provide @p c_arr_ptr pointing to an array of the proper shape (the
  !! BMI function get_var_shape() can be used to create it). Multi-dimensional
  !! arrays are supported.
  !<
  function get_value(c_var_address, c_arr_ptr) result(bmi_status) &
    bind(C, name="get_value")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_value
    ! -- modules
    use ConstantsModule, only: LENMEMTYPE
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    type(c_ptr), intent(inout) :: c_arr_ptr !< pointer to the array
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENMEMTYPE) :: mem_type
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    call get_mem_type(var_name, mem_path, mem_type)

    if (index(mem_type, "DOUBLE") /= 0) then
      bmi_status = get_value_double(c_var_address, c_arr_ptr)
    else if (index(mem_type, "INTEGER") /= 0) then
      bmi_status = get_value_int(c_var_address, c_arr_ptr)
    else if (index(mem_type, "LOGICAL") /= 0) then
      bmi_status = get_value_bool(c_var_address, c_arr_ptr)
    else if (index(mem_type, "STRING") /= 0) then
      bmi_status = get_value_string(c_var_address, c_arr_ptr)
    else
      write (bmi_last_error, fmt_unsupported_type) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

  end function get_value

  !> @brief Get a pointer to an array
  !!
  !! The array is located at @p c_var_address. There is no copying of data involved.
  !! Multi-dimensional arrays are supported and the get_var_rank() function
  !! can be used to get the variable's dimensionality, and get_var_shape() for
  !! its shape.
  !<
  function get_value_ptr(c_var_address, c_arr_ptr) result(bmi_status) &
    bind(C, name="get_value_ptr")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_value_ptr
    ! -- modules
    use ConstantsModule, only: LENMEMTYPE
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    type(c_ptr), intent(inout) :: c_arr_ptr !< pointer to the array
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENMEMTYPE) :: mem_type
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    call get_mem_type(var_name, mem_path, mem_type)

    if (index(mem_type, "DOUBLE") /= 0) then
      bmi_status = get_value_ptr_double(c_var_address, c_arr_ptr)
    else if (index(mem_type, "INTEGER") /= 0) then
      bmi_status = get_value_ptr_int(c_var_address, c_arr_ptr)
    else if (index(mem_type, "LOGICAL") /= 0) then
      bmi_status = get_value_ptr_bool(c_var_address, c_arr_ptr)
    else
      write (bmi_last_error, fmt_unsupported_type) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

  end function get_value_ptr

  !> @brief Get a pointer to the array of double precision numbers
  !!
  !! The array is located at @p c_var_address. There is no copying of data involved.
  !! Multi-dimensional arrays are supported and the get_var_rank() function
  !! can be used to get the variable's dimensionality, and get_var_shape() for
  !! its shape.
  !<
  function get_value_ptr_double(c_var_address, c_arr_ptr) result(bmi_status) &
    bind(C, name="get_value_ptr_double")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_value_ptr_double
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    type(c_ptr), intent(inout) :: c_arr_ptr !< pointer to the array
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid
    real(DP), pointer :: scalar_ptr
    real(DP), dimension(:), pointer, contiguous :: array_ptr
    real(DP), dimension(:, :), pointer, contiguous :: array2D_ptr
    real(DP), dimension(:, :, :), pointer, contiguous :: array3D_ptr
    integer(I4B) :: rank

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    rank = -1
    call get_mem_rank(var_name, mem_path, rank)
    if (rank == 0) then
      call mem_setptr(scalar_ptr, var_name, mem_path)
      c_arr_ptr = c_loc(scalar_ptr)
    else if (rank == 1) then
      call mem_setptr(array_ptr, var_name, mem_path)
      c_arr_ptr = c_loc(array_ptr)
    else if (rank == 2) then
      call mem_setptr(array2D_ptr, var_name, mem_path)
      c_arr_ptr = c_loc(array2D_ptr)
    else if (rank == 3) then
      call mem_setptr(array3D_ptr, var_name, mem_path)
      c_arr_ptr = c_loc(array3D_ptr)
    else
      write (bmi_last_error, fmt_unsupported_rank) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

  end function get_value_ptr_double

  !> @brief Get a pointer to the array of integer numbers
  !!
  !! The array is located at @p c_var_address. There is no copying of data involved.
  !! Multi-dimensional arrays are supported and the get_var_rank() function
  !! can be used to get the variable's dimensionality.
  !<
  function get_value_ptr_int(c_var_address, c_arr_ptr) result(bmi_status) &
    bind(C, name="get_value_ptr_int")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_value_ptr_int
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    type(c_ptr), intent(inout) :: c_arr_ptr !< pointer to the array
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid
    integer(I4B) :: rank
    integer(I4B), pointer :: scalar_ptr
    integer(I4B), dimension(:), pointer, contiguous :: array_ptr
    integer(I4B), dimension(:, :), pointer, contiguous :: array2D_ptr
    integer(I4B), dimension(:, :, :), pointer, contiguous :: array3D_ptr

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    rank = -1
    call get_mem_rank(var_name, mem_path, rank)

    if (rank == 0) then
      call mem_setptr(scalar_ptr, var_name, mem_path)
      c_arr_ptr = c_loc(scalar_ptr)
    else if (rank == 1) then
      call mem_setptr(array_ptr, var_name, mem_path)
      c_arr_ptr = c_loc(array_ptr)
    else if (rank == 2) then
      call mem_setptr(array2D_ptr, var_name, mem_path)
      c_arr_ptr = c_loc(array2D_ptr)
    else if (rank == 3) then
      call mem_setptr(array3D_ptr, var_name, mem_path)
      c_arr_ptr = c_loc(array3D_ptr)
    else
      write (bmi_last_error, fmt_unsupported_rank) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

  end function get_value_ptr_int

  !> @brief Get a pointer to the logical scalar value
  !!
  !! Only scalar values (with rank=0) are supported.
  !<
  function get_value_ptr_bool(c_var_address, c_arr_ptr) result(bmi_status) &
    bind(C, name="get_value_ptr_bool")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_value_ptr_bool
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    type(c_ptr), intent(inout) :: c_arr_ptr !< pointer to the array
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid
    logical(LGP), pointer :: scalar_ptr
    integer(I4B) :: rank

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    rank = -1
    call get_mem_rank(var_name, mem_path, rank)
    if (rank == 0) then
      call mem_setptr(scalar_ptr, var_name, mem_path)
      c_arr_ptr = c_loc(scalar_ptr)
    else
      write (bmi_last_error, fmt_unsupported_rank) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

  end function get_value_ptr_bool

  !> @brief Set new values for a given variable
  !!
  !! The array pointed to by @p c_arr_ptr can have rank equal to 0, 1, or 2
  !! and should have a C-style layout, which is particularly important for
  !! rank > 1.
  !<
  function set_value(c_var_address, c_arr_ptr) result(bmi_status) &
    bind(C, name="set_value")
    !DIR$ ATTRIBUTES DLLEXPORT :: set_value
    ! -- modules
    use ConstantsModule, only: LENMEMTYPE
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    type(c_ptr), intent(inout) :: c_arr_ptr !< pointer to the array
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENMEMTYPE) :: mem_type
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    call get_mem_type(var_name, mem_path, mem_type)

    if (index(mem_type, "DOUBLE") /= 0) then
      bmi_status = set_value_double(c_var_address, c_arr_ptr)
    else if (index(mem_type, "INTEGER") /= 0) then
      bmi_status = set_value_int(c_var_address, c_arr_ptr)
    else if (index(mem_type, "LOGICAL") /= 0) then
      bmi_status = set_value_bool(c_var_address, c_arr_ptr)
    else
      write (bmi_last_error, fmt_unsupported_type) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

  end function set_value

  !> @brief Set new values for a variable of type double
  !!
  !! The array pointed to by @p c_arr_ptr can have rank equal to 0, 1, or 2
  !! and should have a C-style layout, which is particularly important for
  !! rank > 1.
  !<
  function set_value_double(c_var_address, c_arr_ptr) result(bmi_status) &
    bind(C, name="set_value_double")
    !DIR$ ATTRIBUTES DLLEXPORT :: set_value_double
    ! -- modules
    use MemorySetHandlerModule, only: on_memory_set
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    type(c_ptr), intent(in) :: c_arr_ptr !< pointer to the double precision array
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid
    integer(I4B) :: rank
    real(DP), pointer :: src_ptr, tgt_ptr
    real(DP), dimension(:), pointer, contiguous :: src1D_ptr, tgt1D_ptr
    real(DP), dimension(:, :), pointer, contiguous :: src2D_ptr, tgt2D_ptr
    integer(I4B) :: i, j
    integer(I4B) :: status

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    ! convert pointer and copy, using loops to avoid stack overflow
    rank = -1
    call get_mem_rank(var_name, mem_path, rank)

    if (rank == 0) then
      call mem_setptr(tgt_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, src_ptr)
      tgt_ptr = src_ptr
    else if (rank == 1) then
      call mem_setptr(tgt1D_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, src1D_ptr, shape(tgt1D_ptr))
      do i = 1, size(tgt1D_ptr)
        tgt1D_ptr(i) = src1D_ptr(i)
      end do
    else if (rank == 2) then
      call mem_setptr(tgt2D_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, src2D_ptr, shape(tgt2D_ptr))
      do j = 1, size(tgt2D_ptr, 2)
        do i = 1, size(tgt2D_ptr, 1)
          tgt2D_ptr(i, j) = src2D_ptr(i, j)
        end do
      end do
    else
      write (bmi_last_error, fmt_unsupported_rank) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

    ! trigger event:
    call on_memory_set(var_name, mem_path, status)
    if (status /= 0) then
      ! something went terribly wrong here, aborting
      write (bmi_last_error, fmt_invalid_mem_access) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

  end function set_value_double

  !> @brief Set new values for a variable of type integer
  !!
  !! The array pointed to by @p c_arr_ptr can have rank equal to 0, 1, or 2.
  !<
  function set_value_int(c_var_address, c_arr_ptr) result(bmi_status) &
    bind(C, name="set_value_int")
    !DIR$ ATTRIBUTES DLLEXPORT :: set_value_int
    ! -- modules
    use MemorySetHandlerModule, only: on_memory_set
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    type(c_ptr), intent(in) :: c_arr_ptr !< pointer to the integer array
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid
    integer(I4B) :: rank
    integer(I4B), pointer :: src_ptr, tgt_ptr
    integer(I4B), dimension(:), pointer, contiguous :: src1D_ptr, tgt1D_ptr
    integer(I4B), dimension(:, :), pointer, contiguous :: src2D_ptr, tgt2D_ptr
    integer(I4B) :: i, j
    integer(I4B) :: status

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    ! convert pointer and copy, using loops to avoid stack overflow
    rank = -1
    call get_mem_rank(var_name, mem_path, rank)

    if (rank == 0) then
      call mem_setptr(tgt_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, src_ptr)
      tgt_ptr = src_ptr
    else if (rank == 1) then
      call mem_setptr(tgt1D_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, src1D_ptr, shape(tgt1D_ptr))
      do i = 1, size(tgt1D_ptr)
        tgt1D_ptr(i) = src1D_ptr(i)
      end do
    else if (rank == 2) then
      call mem_setptr(tgt2D_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, src2D_ptr, shape(tgt2D_ptr))
      do j = 1, size(tgt2D_ptr, 2)
        do i = 1, size(tgt2D_ptr, 1)
          tgt2D_ptr(i, j) = src2D_ptr(i, j)
        end do
      end do
    else
      write (bmi_last_error, fmt_unsupported_rank) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

    ! trigger event:
    call on_memory_set(var_name, mem_path, status)
    if (status /= 0) then
      ! something went terribly wrong here, aborting
      write (bmi_last_error, fmt_invalid_mem_access) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

  end function set_value_int

  !> @brief Set new value for a logical scalar variable
  !!
  !! The array pointed to by @p c_arr_ptr must have a rank equal to 0.
  !<
  function set_value_bool(c_var_address, c_arr_ptr) result(bmi_status) &
    bind(C, name="set_value_bool")
    !DIR$ ATTRIBUTES DLLEXPORT :: set_value_bool
    ! -- modules
    use MemorySetHandlerModule, only: on_memory_set
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    type(c_ptr), intent(in) :: c_arr_ptr !< pointer to the logical array
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid
    integer(I4B) :: rank
    logical(LGP), pointer :: src_ptr, tgt_ptr
    integer(I4B) :: status

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    rank = -1
    call get_mem_rank(var_name, mem_path, rank)

    ! convert pointer
    if (rank == 0) then
      call mem_setptr(tgt_ptr, var_name, mem_path)
      call c_f_pointer(c_arr_ptr, src_ptr)
      tgt_ptr = src_ptr
    else
      write (bmi_last_error, fmt_unsupported_rank) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

    ! trigger event:
    call on_memory_set(var_name, mem_path, status)
    if (status /= 0) then
      ! something went terribly wrong here, aborting
      write (bmi_last_error, fmt_invalid_mem_access) trim(var_name)
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
      return
    end if

  end function set_value_bool

  !> @brief Get the variable type as a string
  !!
  !! The type returned is that of a single element.
  !! When the variable cannot be found, the string
  !< 'UNKNOWN' is assigned.
  function get_var_type(c_var_address, c_var_type) result(bmi_status) &
    bind(C, name="get_var_type")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_var_type
    ! -- modules
    use ConstantsModule, only: LENMEMTYPE
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    character(kind=c_char), intent(out) :: c_var_type(BMI_LENVARTYPE) !< variable type as a string
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    character(len=LENMEMTYPE) :: mem_type
    logical(LGP) :: valid

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    call get_mem_type(var_name, mem_path, mem_type)
    c_var_type(1:len(trim(mem_type)) + 1) = &
      string_to_char_array(trim(mem_type), len(trim(mem_type)))

    if (mem_type == 'UNKNOWN') then
      write (bmi_last_error, fmt_general_err) 'unknown memory type'
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
    end if

  end function get_var_type

  !> @brief Get the variable rank (non-BMI)
  !!
  !! In order to support multi-dimensional arrays, this function gives
  !! access to the rank of the array.
  !<
  function get_var_rank(c_var_address, c_var_rank) result(bmi_status) &
    bind(C, name="get_var_rank")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_var_rank
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    integer(kind=c_int), intent(out) :: c_var_rank !< variable rank
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    call get_mem_rank(var_name, mem_path, c_var_rank)
    if (c_var_rank == -1) then
      bmi_status = BMI_FAILURE
      return
    end if

  end function get_var_rank

  !> @brief Get the shape of the array for the variable (non-BMI)
  !!
  !! The shape is an integer array with size equal to the rank of
  !! the variable (see get_var_rank()) and values that give the
  !! length of the array in each dimension. The target shape array
  !! @p c_var_shape should be allocated before calling this routine.
  !!
  !! Note that the returned shape representation will has been converted
  !! to C-style.
  !<
  function get_var_shape(c_var_address, c_var_shape) result(bmi_status) &
    bind(C, name="get_var_shape")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_var_shape
    ! -- modules
    use ConstantsModule, only: MAXMEMRANK
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_var_address(*) !< memory address string of the variable
    integer(c_int), intent(inout) :: c_var_shape(*) !< 1D array with the variable's shape
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    integer(I4B), dimension(MAXMEMRANK) :: var_shape
    integer(I4B) :: var_rank
    character(len=LENMEMPATH) :: mem_path
    character(len=LENVARNAME) :: var_name
    logical(LGP) :: valid

    bmi_status = BMI_SUCCESS

    call split_address(c_var_address, mem_path, var_name, valid)
    if (.not. valid) then
      bmi_status = BMI_FAILURE
      return
    end if

    var_shape = 0
    var_rank = 0
    call get_mem_rank(var_name, mem_path, var_rank)
    call get_mem_shape(var_name, mem_path, var_shape)
    if (var_shape(1) == -1 .or. var_rank == -1) then
      bmi_status = BMI_FAILURE
      return
    end if

    ! External calls to this BMI are assumed C style, so if the internal shape
    ! is (100,1) we get (100,1,undef) from the call get_mem_shape
    ! This we need to convert to C-style which should be (1,100).
    ! Hence, we reverse the array and drop undef:
    c_var_shape(1:var_rank) = var_shape(var_rank:1:-1)

  end function get_var_shape

end module mf6bmi
