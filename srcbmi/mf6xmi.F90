!> @brief This module contains the eXtended Model Interface
!!
!! In this module we expose functionality in the MODFLOW 6 shared library,
!! that is _beyond_ the basic model interface: https://bmi-spec.readthedocs.io/en/latest/.
!! The main extensions are
!!
!! - Controlling the kernel at a finer granularity, isolating the call to the linear solve
!!   of the system. This way, the interface can be used to set up a non-linear coupling
!!   with, for example, an external unsaturated zone model.
!! - Expose the concept of subcomponents, which in case of MODFLOW 6 are 'Numerical Solution'
!!   objects, each of which represents a separate linear system to solve. An example here
!!   would be a transport model (GWT) coupled to a groundwater model (GWF).
!!
!! The common BMI control flow is
!!
!! ~~~{.py}
!! initialize()
!!
!! while t < t_end:
!!   update()
!!
!! finalize()
!! ~~~
!!
!! With the XMI you can now also use it as:
!!
!! ~~~{.py}
!! initialize()
!!
!! while(t < t_end):
!!
!!   prepare_time_step()
!!
!!   # modify some values here
!!
!!   do_time_step()
!!
!!   # and maybe something here as well
!!
!!   finalize_time_step()
!!
!! finalize()
!! ~~~
!!
!! Or, when you want to isolate the call to the linear solve, a
!! typical application could look like this:
!!
!! ~~~{.py}
!! initialize()
!!
!! while t < t_end:
!!
!!   prepare_time_step()
!!
!!   for i_sol in solutions:
!!
!!     prepare_solve(i_sol)
!!
!!       while k_iter < max_iter:
!!
!!         # exchange coupled variables
!!         exchange_data()
!!
!!         # the MODFLOW linear solve:
!!         solve()
!!
!!         # maybe solve some other, external model here:
!!         solveExternalModel()
!!
!!         # and exchange back
!!         exchange_data()
!!
!!         # check for convergence
!!         convergence_check()
!!
!!     finalize_solve(i_sol)
!!
!!   finalize_time_step()
!!
!! finalize()
!! ~~~
!!
!! Note that the last example can only work when there is a single Solution Group defined.
!! This will typically not be a problem, though, as applications with multiple Solution Groups
!! should be quite rare.
!<
module mf6xmi
  use mf6bmi
  use mf6bmiUtil
  use mf6bmiError
  use Mf6CoreModule
  use KindModule
  use iso_c_binding, only: c_int, c_char, c_double
  implicit none

  integer(I4B), pointer :: iterationCounter => null() !< the counter for the outer iteration loop, initialized in xmi_prepare_iteration()

contains

#if defined(__WITH_MPI__)
  function xmi_initialize_mpi(mpi_comm) result(bmi_status) &
    bind(C, name="initialize_mpi")
    use MpiWorldModule
    use SimVariablesModule, only: simulation_mode
    !DIR$ ATTRIBUTES DLLEXPORT :: xmi_initialize_mpi
    ! -- dummy variables
    integer(kind=c_int) :: mpi_comm !< the Fortran communicator (as an integer)
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    type(MpiWorldType), pointer :: mpi_world => null()

    ! set parallel
    mpi_world => get_mpi_world()
    call mpi_world%set_comm(mpi_comm)
    simulation_mode = 'PARALLEL'

    ! regular initialize
    bmi_status = bmi_initialize()

  end function xmi_initialize_mpi
#endif

  !> @brief Prepare a single time step
  !!
  !! The routine takes the time step \p dt as an argument. However, MODFLOW (currently)
  !! does not allow to alter this value after initialization, so it is ignored
  !! here.
  !<
  function xmi_prepare_time_step(dt) result(bmi_status) &
    bind(C, name="prepare_time_step")
    !DIR$ ATTRIBUTES DLLEXPORT :: xmi_prepare_time_step
    ! -- dummy variables
    real(kind=c_double), intent(in) :: dt !< time step
    integer(kind=c_int) :: bmi_status !< BMI status code

    call Mf6PrepareTimestep()
    bmi_status = BMI_SUCCESS

  end function xmi_prepare_time_step

  !> @brief Perform a single time step
  !!
  !! It does so by looping over all solution groups, and calling
  !! the calculate function on all solutions in there.
  !<
  function xmi_do_time_step() result(bmi_status) bind(C, name="do_time_step")
    !DIR$ ATTRIBUTES DLLEXPORT :: xmi_do_time_step
    ! -- dummy variables
    integer(kind=c_int) :: bmi_status !< BMI status code

    call Mf6DoTimestep()
    bmi_status = BMI_SUCCESS

  end function xmi_do_time_step

  !> @brief Finalize the time step
  !!
  !! This will mostly write output and messages. It is essential
  !! to call this to finish the time step.
  !<
  function xmi_finalize_time_step() result(bmi_status) &
    bind(C, name="finalize_time_step")
    !DIR$ ATTRIBUTES DLLEXPORT :: xmi_finalize_time_step
    ! -- dummy variables
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    logical :: hasConverged

    hasConverged = Mf6FinalizeTimestep()
    if (hasConverged) then
      bmi_status = BMI_SUCCESS
    else
      write (bmi_last_error, fmt_general_err) 'simulation failed to converge'
      call report_bmi_error(bmi_last_error)
      bmi_status = BMI_FAILURE
    end if

  end function xmi_finalize_time_step

  !> @brief This will get the number of Numerical Solutions in the simulation
  !!
  !! For most applications, this number will be equal to 1. Note that this part
  !! of the XMI only works when the simulation is defined with a single
  !! Solution Group. (If you don't know what a Solution Group is, then
  !! you are most likely not using more than one...)
  !<
  function xmi_get_subcomponent_count(count) result(bmi_status) &
    bind(C, name="get_subcomponent_count")
    !DIR$ ATTRIBUTES DLLEXPORT :: xmi_get_subcomponent_count
    ! -- modules
    use ListsModule, only: solutiongrouplist
    use SimVariablesModule, only: istdout
    ! -- dummy variables
    integer(kind=c_int), intent(out) :: count !< number of solutions
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    class(SolutionGroupType), pointer :: sgp

    ! the following is true for all calls at this level (subcomponent)
    if (solutiongrouplist%Count() /= 1) then
      write (istdout, *) &
        'Error: BMI does not support the use of multiple solution groups'
      count = -1
      bmi_status = BMI_FAILURE
      return
    end if

    sgp => GetSolutionGroupFromList(solutiongrouplist, 1)
    count = sgp%nsolutions
    bmi_status = BMI_SUCCESS

  end function xmi_get_subcomponent_count

  !> @brief Prepare for solving the system
  !!
  !! This preparation mostly consists of advancing the solutions, models, and exchanges
  !! in the simulation. The index \p subcomponent_idx runs from 1 to the value returned
  !! by xmi_get_subcomponent_count().
  !<
  function xmi_prepare_solve(subcomponent_idx) result(bmi_status) &
    bind(C, name="prepare_solve")
    !DIR$ ATTRIBUTES DLLEXPORT :: xmi_prepare_solve
    ! -- modules
    use ListsModule, only: solutiongrouplist
    use BaseSolutionModule, only: BaseSolutionType
    use SimVariablesModule, only: istdout
    ! -- dummy variables
    integer(kind=c_int) :: subcomponent_idx !< index of the subcomponent (i.e. Numerical Solution)
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    class(BaseSolutionType), pointer :: bs

    ! people might not call 'xmi_get_subcomponent_count' first, so let's repeat this:
    if (solutiongrouplist%Count() /= 1) then
      write (istdout, *) &
        'Error: BMI does not support the use of multiple solution groups'
      bmi_status = BMI_FAILURE
      return
    end if

    ! get the solution we are running
    bs => getSolution(subcomponent_idx)

    ! *_ad (model, exg, sln)
    call bs%prepareSolve()

    ! reset counter
    allocate (iterationCounter)
    iterationCounter = 0

    bmi_status = BMI_SUCCESS

  end function xmi_prepare_solve

  !> @brief Build and solve the linear system
  !!
  !! The solve is called on the Numerical Solution indicated by the value of \p subcomponent_idx,
  !! which runs from 1 to the value returned by xmi_get_subcomponent_count(). Before calling
  !! this, a matching call to xmi_prepare_solve() should be done.
  !<
  function xmi_solve(subcomponent_idx, has_converged) result(bmi_status) &
    bind(C, name="solve")
    !DIR$ ATTRIBUTES DLLEXPORT :: xmi_solve
    ! -- modules
    use BaseSolutionModule, only: BaseSolutionType
    use NumericalSolutionModule, only: NumericalSolutionType
    use ExplicitSolutionModule, only: ExplicitSolutionType
    ! -- dummy variables
    integer(kind=c_int), intent(in) :: subcomponent_idx !< index of the subcomponent (i.e. Numerical Solution)
    integer(kind=c_int), intent(out) :: has_converged !< equal to 1 for convergence, 0 otherwise
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    class(BaseSolutionType), pointer :: bs

    ! get the numerical solution we are running
    bs => getSolution(subcomponent_idx)

    ! execute the nth iteration
    iterationCounter = iterationCounter + 1
    call bs%solve(iterationCounter)

    ! the following check is equivalent to that in NumericalSolution%sln_ca
    select type (bs)
    class is (NumericalSolutionType)
      if (bs%icnvg == 1) then
        has_converged = 1
      else
        has_converged = 0
      end if
    class is (ExplicitSolutionType)
      has_converged = 1
    end select

    bmi_status = BMI_SUCCESS

  end function xmi_solve

  !> @brief Finalize the solve of the system
  !!
  !! This will determine convergence, reports, calculate flows and budgets, and more...
  !! It should always follow after a call to xmi_prepare_solve() and xmi_solve().
  !!
  !! @param[in]   subcomponent_idx    the index of the subcomponent (Numerical Solution)
  !! @return      bmi_status          the BMI status code
  !<
  function xmi_finalize_solve(subcomponent_idx) result(bmi_status) &
    bind(C, name="finalize_solve")
    !DIR$ ATTRIBUTES DLLEXPORT :: xmi_finalize_solve
    ! -- modules
    use BaseSolutionModule, only: BaseSolutionType
    ! -- dummy variables
    integer(kind=c_int), intent(in) :: subcomponent_idx !< index of the subcomponent (i.e. Numerical Solution)
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    class(BaseSolutionType), pointer :: bs
    integer(I4B) :: hasConverged

    ! get the numerical solution we are running
    bs => getSolution(subcomponent_idx)

    ! hasConverged is equivalent to the isgcnvg variable which is initialized to 1,
    ! see the body of the picard loop in SolutionGroupType%sgp_ca
    hasConverged = 1

    ! finish up
    call bs%finalizeSolve(iterationCounter, hasConverged, 0)

    ! check convergence on solution
    if (.not. hasConverged == 1) then
      write (bmi_last_error, fmt_fail_cvg_sol) subcomponent_idx
      call report_bmi_error(bmi_last_error)
    end if

    ! non-convergence is no reason to crash the API:
    bmi_status = BMI_SUCCESS

    ! clear this for safety
    deallocate (iterationCounter)

  end function xmi_finalize_solve

  !> @brief Get the version string for this component
  !<
  function xmi_get_version(mf_version) result(bmi_status) &
    bind(C, name="get_version")
    !DIR$ ATTRIBUTES DLLEXPORT :: xmi_get_version
    ! -- modules
    use VersionModule, only: VERSIONNUMBER, IDEVELOPMODE
    ! -- dummy variables
    character(kind=c_char), intent(inout) :: mf_version(BMI_LENVERSION)
    integer(kind=c_int) :: bmi_status !< BMI status code

    mf_version = string_to_char_array(VERSIONNUMBER, len_trim(VERSIONNUMBER))
    bmi_status = BMI_SUCCESS

  end function xmi_get_version

  !> @brief Get the full address string for a variable
  !!
  !! This routine constructs the full address string of a variable using the
  !! exact same logic as the internal memory manager. This routine
  !! should always be used when accessing a variable through the BMI
  !! to assure compatibility with future versions of the library.
  !<
  function get_var_address(c_component_name, c_subcomponent_name, &
                           c_var_name, c_var_address) &
    result(bmi_status) bind(C, name="get_var_address")
    !DIR$ ATTRIBUTES DLLEXPORT :: get_var_address
    ! -- modules
    use MemoryHelperModule, only: create_mem_path, create_mem_address
    use ConstantsModule, only: LENCOMPONENTNAME, LENVARNAME, LENMEMPATH, &
                               LENMEMADDRESS
    ! -- dummy variables
    character(kind=c_char), intent(in) :: c_component_name(*) !< name of the component (a Model or Solution)
    character(kind=c_char), intent(in) :: c_subcomponent_name(*) !< name of the subcomponent (Package), or an empty
                                                                            !! string'' when not applicable
    character(kind=c_char), intent(in) :: c_var_name(*) !< name of the variable
    character(kind=c_char), intent(out) :: c_var_address(BMI_LENVARADDRESS) !< full address of the variable
    integer(kind=c_int) :: bmi_status !< BMI status code
    ! -- local variables
    character(len=LENCOMPONENTNAME) :: component_name
    character(len=LENCOMPONENTNAME) :: subcomponent_name
    character(len=LENVARNAME) :: variable_name
    character(len=LENMEMPATH) :: mem_path
    character(len=LENMEMADDRESS) :: mem_address

    ! convert to Fortran strings
    component_name = char_array_to_string(c_component_name, &
                                          strlen(c_component_name, &
                                                 LENCOMPONENTNAME + 1))
    subcomponent_name = char_array_to_string(c_subcomponent_name, &
                                             strlen(c_subcomponent_name, &
                                                    LENCOMPONENTNAME + 1))
    variable_name = char_array_to_string(c_var_name, &
                                         strlen(c_var_name, &
                                                LENVARNAME + 1))

    ! create memory address
    if (subcomponent_name == '') then
      mem_path = create_mem_path(component_name)
    else
      mem_path = create_mem_path(component_name, subcomponent_name)
    end if
    mem_address = create_mem_address(mem_path, variable_name)

    ! convert to c string:
    c_var_address(1:len(trim(mem_address)) + 1) = &
      string_to_char_array(trim(mem_address), len(trim(mem_address)))

    bmi_status = BMI_SUCCESS

  end function get_var_address

end module mf6xmi
