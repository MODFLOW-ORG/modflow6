module MemoryHelperModule
  use KindModule, only: I4B, LGP
  use ConstantsModule, only: LENMEMPATH, LENMEMSEPARATOR, LENMEMADDRESS, &
                             LENVARNAME, LENCOMPONENTNAME, LENCONTEXTNAME
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg

  implicit none

  character(len=LENMEMSEPARATOR), parameter :: memPathSeparator = '/' !< used to build up the memory address for the stored variables

contains

  !> @brief returns the path to the memory object
  !!
  !! Returns the path to the location in the memory manager where
  !! the variables for this (sub)component are stored, the 'memoryPath'
  !!
  !! NB: no need to trim the input parameters
  !<
  function create_mem_path(component, subcomponent, context) result(memory_path)
    character(len=*), intent(in) :: component !< name of the solution, model, or exchange
    character(len=*), intent(in), optional :: subcomponent !< name of the package (optional)
    character(len=*), intent(in), optional :: context !< name of the context (optional)
    character(len=LENMEMPATH) :: memory_path !< the memory path

    call mem_check_length(component, LENCOMPONENTNAME, "solution/model/exchange")

    if (present(subcomponent)) then
      call mem_check_length(subcomponent, LENCOMPONENTNAME, "package")
    end if

    if (present(context)) then
      call mem_check_length(context, LENCONTEXTNAME, "context")
    end if

    memory_path = trim(component)

    if (present(subcomponent)) then
      memory_path = trim(memory_path)//memPathSeparator//trim(subcomponent)
    end if

    if (present(context)) then
      memory_path = trim(context)//memPathSeparator//trim(memory_path)
    end if

  end function create_mem_path

  !> @brief returns the address string of the memory object
  !!
  !! Returns the memory address, i.e. the full path plus name of the stored variable
  !!
  !! NB: no need to trim the input parameters
  !<
  function create_mem_address(mem_path, var_name) result(mem_address)
    character(len=*), intent(in) :: mem_path !< path to the memory object
    character(len=*), intent(in) :: var_name !< name of the stored variable
    character(len=LENMEMADDRESS) :: mem_address !< full address string to the memory object

    call mem_check_length(mem_path, LENMEMPATH, "memory path")
    call mem_check_length(var_name, LENVARNAME, "variable")

    mem_address = trim(mem_path)//memPathSeparator//trim(var_name)

  end function create_mem_address

  !> @brief Split a memory address string into memory path and variable name
  !<
  subroutine split_mem_address(mem_address, mem_path, var_name, success)
    character(len=*), intent(in) :: mem_address !< the full memory address string
    character(len=LENMEMPATH), intent(out) :: mem_path !< the memory path
    character(len=LENVARNAME), intent(out) :: var_name !< the variable name
    logical(LGP), intent(out) :: success !< true when successful
    ! local
    integer(I4B) :: idx

    idx = index(mem_address, memPathSeparator, back=.true.)

    ! if no separator, or it's at the end of the string,
    ! the memory address is not valid:
    if (idx < 1 .or. idx == len(mem_address)) then
      success = .false.
      mem_path = ''
      var_name = ''
    else
      success = .true.
      mem_path = mem_address(:idx - 1)
      var_name = mem_address(idx + 1:)
    end if

    ! remove context specifier if prepended to mempath
    !if (success) then
    !  idx = index(mem_path, memPathSeparator, back=.true.)
    !  if (idx > 0 .and. mem_path(1:2) == '__') then
    !    mem_path = mem_path(idx + 1:)
    !  end if
    !end if

  end subroutine split_mem_address

  !> @brief Split the memory path into component(s)
  !!
  !! NB: when there is no subcomponent in the path, the
  !! value for @par subcomponent is set to an empty string.
  !<
  subroutine split_mem_path(mem_path, component, subcomponent)
    character(len=*), intent(in) :: mem_path !< path to the memory object
    character(len=LENCOMPONENTNAME), intent(out) :: component !< name of the component (solution, model, exchange)
    character(len=LENCOMPONENTNAME), intent(out) :: subcomponent !< name of the subcomponent (package)
    ! local
    character(len=LENMEMPATH) :: local_mem_path
    integer(I4B) :: idx

    call strip_context_mem_path(mem_path, local_mem_path)

    idx = index(local_mem_path, memPathSeparator, back=.true.)
    ! if the separator is found at the end of the string,
    ! the path is invalid:
    if (idx == len_trim(local_mem_path)) then
      write (errmsg, '(*(G0))') &
        'Fatal error in Memory Manager, cannot split invalid memory path: ', &
        mem_path

      ! -- store error and stop program execution
      call store_error(errmsg, terminate=.TRUE.)
    end if

    if (idx > 0) then
      ! when found:
      component = local_mem_path(:idx - 1)
      subcomponent = local_mem_path(idx + 1:)
    else
      ! when not found, there apparently is no subcomponent:
      component = local_mem_path(:LENCOMPONENTNAME)
      subcomponent = ''
    end if

  end subroutine split_mem_path

  !> @brief Return the context from the memory path
  !!
  !! NB: when there is no context in the memory path, a
  !! empty character string is returned.
  !<
  function get_mem_path_context(mem_path) result(res)
    character(len=*), intent(in) :: mem_path !< path to the memory object
    character(len=LENMEMPATH) :: res !< memory path context
    ! local
    integer(I4B) :: idx

    ! initialize the memory path context
    res = ' '

    if (mem_path(1:2) == '__') then
      idx = index(mem_path, memPathSeparator)
      if (idx > 0) then
        res = mem_path(:idx)
      end if
    end if
  end function get_mem_path_context

  !> @brief Remove the context from the memory path
  !!
  !! NB: when there is no context in the memory path, the
  !! original memory path is returned.
  !<
  subroutine strip_context_mem_path(mem_path, mem_path_no_context)
    character(len=*), intent(in) :: mem_path !< path to the memory object
    character(len=LENMEMPATH), intent(inout) :: mem_path_no_context !< path to the memory object without the context
    ! local
    integer(I4B) :: idx
    character(len=LENMEMPATH) :: context

    ! initialize the local mem_path
    mem_path_no_context = mem_path

    context = get_mem_path_context(mem_path)

    if (len_trim(context) > 0) then
      idx = len_trim(context)
      mem_path_no_context = mem_path(idx + 1:)
    end if

  end subroutine strip_context_mem_path

  !> @brief Generic routine to check the length of (parts of) the memory address
  !!
  !! The string will be trimmed before the measurement.
  !!
  !! @warning{if the length exceeds the maximum, a message is recorded
  !! and the program will be stopped}
  !!
  !! The description should describe the part of the address that is checked
  !! (variable, package, model, solution, exchange name) or the full memory path
  !! itself
  !<
  subroutine mem_check_length(name, max_length, description)
    character(len=*), intent(in) :: name !< string to be checked
    integer(I4B), intent(in) :: max_length !< maximum length
    character(len=*), intent(in) :: description !< a descriptive string

    if (len(trim(name)) > max_length) then
      write (errmsg, '(*(G0))') &
        'Fatal error in Memory Manager, length of ', description, ' must be ', &
        max_length, ' characters or less: ', name, '(len=', len(trim(name)), ')'

      ! -- store error and stop program execution
      call store_error(errmsg, terminate=.TRUE.)
    end if

  end subroutine mem_check_length
end module MemoryHelperModule
