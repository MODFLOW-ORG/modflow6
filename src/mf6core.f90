!> @brief Core MODFLOW 6 module
!!
!! This module contains the core components for MODFLOW 6. This module
!! is used by the stand-alone executable and the share object versions
!! of MODFLOW 6.
!!
!<
module Mf6CoreModule
  use KindModule, only: I4B, LGP
  use ListsModule, only: basesolutionlist, solutiongrouplist, &
                         basemodellist, baseexchangelist, &
                         baseconnectionlist
  use BaseModelModule, only: BaseModelType, GetBaseModelFromList
  use BaseExchangeModule, only: BaseExchangeType, GetBaseExchangeFromList
  use SpatialModelConnectionModule, only: SpatialModelConnectionType, &
                                          get_smc_from_list
  use BaseSolutionModule, only: BaseSolutionType, GetBaseSolutionFromList
  use SolutionGroupModule, only: SolutionGroupType, GetSolutionGroupFromList
  use RunControlModule, only: RunControlType
  use SimStagesModule
  use ProfilerModule
  implicit none

  class(RunControlType), pointer :: run_ctrl => null() !< the run controller for this simulation

contains

  !> @brief Main controller
    !!
    !! This subroutine is the main controller for MODFLOW 6.
    !!
  !<
  subroutine Mf6Run
    ! -- modules
    use CommandArguments, only: GetCommandLineArguments
    use TdisModule, only: endofsimulation
    ! -- local
    logical(LGP) :: hasConverged
    !
    ! -- parse any command line arguments
    call GetCommandLineArguments()
    !
    ! initialize simulation
    call Mf6Initialize()
    !
    ! -- time loop
    do while (.not. endofsimulation)

      ! perform a time step
      hasConverged = Mf6Update()

      ! if not converged, break
      if (.not. hasConverged) exit

    end do
    !
    ! -- finalize simulation
    call Mf6Finalize()

  end subroutine Mf6Run

  !> @brief Initialize a simulation
    !!
    !! This subroutine initializes a MODFLOW 6 simulation. The subroutine:
    !!   - creates the simulation
    !!   - defines
    !!   - allocates and reads static data
    !!
  !<
  subroutine Mf6Initialize()
    ! -- modules
    use RunControlFactoryModule, only: create_run_control
    use SimulationCreateModule, only: simulation_cr
    use SourceLoadModule, only: export_cr

    call g_prof%pre_init()

    ! -- get the run controller for sequential or parallel builds
    run_ctrl => create_run_control()
    call run_ctrl%start()

    ! -- print info and start timer
    call print_info()

    ! -- create mfsim.lst
    call create_lstfile()

    ! -- load input context
    call static_input_load()

    ! init timer and start
    call g_prof%initialize()
    call g_prof%start("Initialize", g_prof%tmr_init)

    ! -- create
    call simulation_cr()

    ! -- define
    call simulation_df()

    ! -- allocate and read
    call simulation_ar()

    ! -- create model exports
    call export_cr()

    ! -- stop the timer
    call g_prof%stop(g_prof%tmr_init)

  end subroutine Mf6Initialize

  !> @brief Run a time step
    !!
    !! This function runs a single time step to completion.
    !!
    !! @return  hasConverged   boolean indicating if convergence was achieved for the time step
    !!
  !<
  function Mf6Update() result(hasConverged)
    logical(LGP) :: hasConverged
    ! start timer
    call g_prof%start("Update", g_prof%tmr_update)
    !
    ! -- prepare timestep
    call Mf6PrepareTimestep()
    !
    ! -- do timestep
    call Mf6DoTimestep()
    !
    ! -- after timestep
    hasConverged = Mf6FinalizeTimestep()

    ! stop timer
    call g_prof%stop(g_prof%tmr_update)

  end function Mf6Update

  !> @brief Finalize the simulation
    !!
    !! This subroutine finalizes a simulation. Steps include:
    !!   - final processing
    !!   - deallocate memory
    !!
  !<
  subroutine Mf6Finalize()
    ! -- modules
    use, intrinsic :: iso_fortran_env, only: output_unit
    use ListsModule, only: lists_da
    use SimulationCreateModule, only: simulation_da
    use TdisModule, only: tdis_da
    use IdmLoadModule, only: idm_da
    use SourceLoadModule, only: export_da
    use SimVariablesModule, only: iout
    ! -- local variables
    integer(I4B) :: im
    integer(I4B) :: ic
    integer(I4B) :: is
    integer(I4B) :: isg
    class(SolutionGroupType), pointer :: sgp => null()
    class(BaseSolutionType), pointer :: sp => null()
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    class(SpatialModelConnectionType), pointer :: mc => null()
    integer(I4B) :: tmr_dealloc

    ! start timer
    call g_prof%start("Finalize", g_prof%tmr_finalize)

    !
    ! -- FINAL PROCESSING (FP)
    ! -- Final processing for each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_fp()
    end do
    !
    ! -- Final processing for each exchange
    do ic = 1, baseexchangelist%Count()
      ep => GetBaseExchangeFromList(baseexchangelist, ic)
      call ep%exg_fp()
    end do
    !
    ! -- Final processing for each solution
    do is = 1, basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%sln_fp()
    end do

    ! start timer for deallocation
    tmr_dealloc = -1
    call g_prof%start("Deallocate", tmr_dealloc)

    !
    ! -- DEALLOCATE (DA)
    ! -- Deallocate tdis
    call tdis_da()
    !
    ! -- Deallocate for each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_da()
      deallocate (mp)
    end do
    !
    ! -- Deallocate for each exchange
    do ic = 1, baseexchangelist%Count()
      ep => GetBaseExchangeFromList(baseexchangelist, ic)
      call ep%exg_da()
      deallocate (ep)
    end do
    !
    ! -- Deallocate for each connection
    do ic = 1, baseconnectionlist%Count()
      mc => get_smc_from_list(baseconnectionlist, ic)
      call mc%exg_da()
      deallocate (mc)
    end do
    !
    ! -- Deallocate for each solution
    do is = 1, basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%sln_da()
      deallocate (sp)
    end do
    !
    ! -- Deallocate solution group and simulation variables
    do isg = 1, solutiongrouplist%Count()
      sgp => GetSolutionGroupFromList(solutiongrouplist, isg)
      call sgp%sgp_da()
      deallocate (sgp)
    end do
    !
    call idm_da(iout)
    call export_da()
    call simulation_da()
    call lists_da()

    ! stop timer
    call g_prof%stop(tmr_dealloc)

    ! finish gently (No calls after this)
    ! timer is stopped inside because this call does not return
    call run_ctrl%finish()

  end subroutine Mf6Finalize

  !> @brief print initial message
  !<
  subroutine print_info()
    use SimModule, only: initial_message
    use TimerModule, only: print_start_time

    ! print initial message
    call initial_message()

    ! get start time
    call print_start_time()

  end subroutine print_info

  !> @brief Set up mfsim list file output logging
    !!
    !! This subroutine creates the mfsim list file
    !! and writes the header.
    !!
  !<
  subroutine create_lstfile()
    use ConstantsModule, only: LINELENGTH
    use SimVariablesModule, only: proc_id, nr_procs, simlstfile, iout
    use InputOutputModule, only: getunit, openfile, append_processor_id
    use MessageModule, only: write_message
    use VersionModule, only: write_listfile_header
    character(len=LINELENGTH) :: line
    !
    ! -- Open simulation list file
    iout = getunit()
    !
    if (nr_procs > 1) then
      call append_processor_id(simlstfile, proc_id)
    end if
    !
    call openfile(iout, 0, simlstfile, 'LIST', filstat_opt='REPLACE')
    !
    ! -- write simlstfile to stdout
    write (line, '(2(1x,A))') 'Writing simulation list file:', &
      trim(adjustl(simlstfile))
    !
    call write_message(line)
    call write_listfile_header(iout)
  end subroutine create_lstfile

  !> @brief Create simulation input context
    !!
    !! This subroutine creates the simulation input context
    !!
  !<
  subroutine static_input_load()
    ! -- modules
    use ConstantsModule, only: LENMEMPATH
    use SimVariablesModule, only: iout
    use IdmLoadModule, only: simnam_load, simtdis_load, &
                             load_models, load_exchanges
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: mem_setptr, mem_allocate
    use SimVariablesModule, only: iparamlog
    !
    ! -- load simnam input context
    call simnam_load(iparamlog)
    !
    ! -- load tdis to input context
    call simtdis_load()
    !
    ! -- load in scope models
    call load_models(iout)
    !
    ! -- load in scope exchanges
    call load_exchanges(iout)
  end subroutine static_input_load

  !> @brief Define the simulation
    !!
    !! This subroutine defined the simulation. Steps include:
    !!   - define each model
    !!   - define each solution
    !!
  !<
  subroutine simulation_df()
    ! -- modules
    use IdmLoadModule, only: idm_df
    ! -- local variables
    integer(I4B) :: im
    integer(I4B) :: ic
    integer(I4B) :: is
    class(BaseSolutionType), pointer :: sp => null()
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    class(SpatialModelConnectionType), pointer :: mc => null()

    ! -- init virtual data environment
    call run_ctrl%at_stage(STG_BFR_MDL_DF)

    ! -- Define each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_df()
    end do
    !
    ! -- synchronize
    call run_ctrl%at_stage(STG_AFT_MDL_DF)
    !
    ! -- Define each exchange
    do ic = 1, baseexchangelist%Count()
      ep => GetBaseExchangeFromList(baseexchangelist, ic)
      call ep%exg_df()
    end do
    !
    ! -- synchronize
    call run_ctrl%at_stage(STG_AFT_EXG_DF)
    !
    ! -- when needed, this is where the interface models are
    !     created and added to the numerical solutions
    call connections_cr()
    !
    ! -- synchronize
    call run_ctrl%at_stage(STG_AFT_CON_CR)
    !
    ! -- synchronize
    call run_ctrl%at_stage(STG_BFR_CON_DF)
    !
    ! -- Define each connection
    do ic = 1, baseconnectionlist%Count()
      mc => get_smc_from_list(baseconnectionlist, ic)
      call mc%exg_df()
    end do
    !
    ! -- synchronize
    call run_ctrl%at_stage(STG_AFT_CON_DF)
    !
    ! -- Define each solution
    do is = 1, basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%sln_df()
    end do

    ! idm df
    call idm_df()

  end subroutine simulation_df

  !> @brief Simulation allocate and read
    !!
    !! This subroutine allocates and reads static data for the simulation.
    !! Steps include:
    !!   - allocate and read for each model
    !!   - allocate and read for each exchange
    !!   - allocate and read for each solution
    !!
  !<
  subroutine simulation_ar()
    use DistVariableModule
    ! -- local variables
    integer(I4B) :: im
    integer(I4B) :: ic
    integer(I4B) :: is
    class(BaseSolutionType), pointer :: sp => null()
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    class(SpatialModelConnectionType), pointer :: mc => null()

    ! -- Allocate and read each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_ar()
    end do
    !
    ! -- Allocate and read each exchange
    do ic = 1, baseexchangelist%Count()
      ep => GetBaseExchangeFromList(baseexchangelist, ic)
      call ep%exg_ar()
    end do
    !
    ! -- Synchronize
    call run_ctrl%at_stage(STG_BFR_CON_AR)
    !
    ! -- Allocate and read all model connections
    do ic = 1, baseconnectionlist%Count()
      mc => get_smc_from_list(baseconnectionlist, ic)
      call mc%exg_ar()
    end do
    !
    ! -- Synchronize
    call run_ctrl%at_stage(STG_AFT_CON_AR)
    !
    ! -- Allocate and read each solution
    do is = 1, basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%sln_ar()
    end do
    !
  end subroutine simulation_ar

  !> @brief Create the model connections from the exchanges
    !!
    !! This will upgrade the numerical exchanges in the solution,
    !! whenever the configuration requires this, to Connection
    !! objects. Currently we anticipate:
    !!
    !!   GWF-GWF => GwfGwfConnection
    !!   GWT-GWT => GwtGwtConecction
  !<
  subroutine connections_cr()
    use ConnectionBuilderModule
    use SimVariablesModule, only: iout
    use VersionModule, only: IDEVELOPMODE
    integer(I4B) :: isol
    type(ConnectionBuilderType) :: connectionBuilder
    class(BaseSolutionType), pointer :: sol => null()
    integer(I4B) :: status
    character(len=16) :: envvar

    write (iout, '(/a)') 'PROCESSING MODEL CONNECTIONS'

    if (baseexchangelist%Count() == 0) then
      ! if this is not a coupled simulation in any way,
      ! then we will not need model connections
      return
    end if

    if (IDEVELOPMODE == 1) then
      call get_environment_variable('DEV_ALWAYS_USE_IFMOD', &
                                    value=envvar, status=status)
      if (status == 0 .and. envvar == '1') then
        connectionBuilder%dev_always_ifmod = .true.
        write (iout, '(/a)') "Development option: forcing interface model"
      end if
    end if

    do isol = 1, basesolutionlist%Count()
      sol => GetBaseSolutionFromList(basesolutionlist, isol)
      call connectionBuilder%processSolution(sol)
    end do

    write (iout, '(a)') 'END OF MODEL CONNECTIONS'
  end subroutine connections_cr

  !> @brief Read and prepare time step
    !!
    !! This subroutine reads and prepares period data for the simulation.
    !! Steps include:
    !!   - read and prepare for each model
    !!   - read and prepare for each exchange
    !!   - reset convergence flag
    !!   - calculate maximum time step for each model
    !!   - calculate maximum time step for each exchange
    !!   - calculate maximum time step for each solution
    !!   - set time discretization timestep using smallest maximum timestep
    !!
  !<
  subroutine Mf6PrepareTimestep()
    ! -- modules
    use KindModule, only: I4B
    use ConstantsModule, only: LINELENGTH, MNORMAL, MVALIDATE
    use TdisModule, only: tdis_set_counters, tdis_set_timestep, &
                          kstp, kper
    use ListsModule, only: basemodellist, baseexchangelist
    use BaseModelModule, only: BaseModelType, GetBaseModelFromList
    use BaseExchangeModule, only: BaseExchangeType, GetBaseExchangeFromList
    use BaseSolutionModule, only: BaseSolutionType, GetBaseSolutionFromList
    use SimModule, only: converge_reset
    use SimVariablesModule, only: isim_mode
    use IdmLoadModule, only: idm_rp
    use SourceLoadModule, only: export_post_prepare
    ! -- local variables
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    class(SpatialModelConnectionType), pointer :: mc => null()
    class(BaseSolutionType), pointer :: sp => null()
    character(len=LINELENGTH) :: line
    character(len=LINELENGTH) :: fmt
    integer(I4B) :: im
    integer(I4B) :: ie
    integer(I4B) :: ic
    integer(I4B) :: is

    ! start timer
    call g_prof%start("Prepare time step", g_prof%tmr_prep_tstp)

    !
    ! -- initialize fmt
    fmt = "(/,a,/)"
    !
    ! -- period update
    call tdis_set_counters()
    !
    ! -- set base line
    write (line, '(a,i0,a,i0,a)') &
      'start timestep kper="', kper, '" kstp="', kstp, '" mode="'
    !
    ! -- evaluate simulation mode
    select case (isim_mode)
    case (MVALIDATE)
      line = trim(line)//'validate"'
    case (MNORMAL)
      line = trim(line)//'normal"'
    end select

    ! -- load dynamic input
    call idm_rp()

    ! -- Read and prepare each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_message(line, fmt=fmt)
      call mp%model_rp()
    end do
    !
    ! -- Synchronize
    call run_ctrl%at_stage(STG_BFR_EXG_RP)
    !
    ! -- Read and prepare each exchange
    do ie = 1, baseexchangelist%Count()
      ep => GetBaseExchangeFromList(baseexchangelist, ie)
      call ep%exg_rp()
    end do
    !
    ! -- Read and prepare each connection
    do ic = 1, baseconnectionlist%Count()
      mc => get_smc_from_list(baseconnectionlist, ic)
      call mc%exg_rp()
    end do
    !
    ! -- Synchronize
    call run_ctrl%at_stage(STG_AFT_CON_RP)
    !
    ! -- reset simulation convergence flag
    call converge_reset()
    !
    ! -- time update for each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_dt()
    end do
    !
    ! -- time update for each exchange
    do ie = 1, baseexchangelist%Count()
      ep => GetBaseExchangeFromList(baseexchangelist, ie)
      call ep%exg_dt()
    end do
    !
    ! -- time update for each connection
    do ic = 1, baseconnectionlist%Count()
      mc => get_smc_from_list(baseconnectionlist, ic)
      call mc%exg_dt()
    end do
    !
    ! -- time update for each solution
    do is = 1, basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%sln_dt()
    end do
    !
    ! -- update exports
    call export_post_prepare()
    !
    ! -- set time step
    call tdis_set_timestep()

    ! stop timer
    call g_prof%stop(g_prof%tmr_prep_tstp)

  end subroutine Mf6PrepareTimestep

  !> @brief Run time step
    !!
    !! This subroutine runs a single time step for the simulation.
    !! Steps include:
    !!   - formulate the system of equations for each model and exchange
    !!   - solve each solution
    !!
  !<
  subroutine Mf6DoTimestep()
    ! --  modules
    use KindModule, only: I4B
    use ListsModule, only: solutiongrouplist
    use SimVariablesModule, only: iFailedStepRetry
    use SolutionGroupModule, only: SolutionGroupType, GetSolutionGroupFromList
    use IdmLoadModule, only: idm_ad
    ! -- local variables
    class(SolutionGroupType), pointer :: sgp => null()
    integer(I4B) :: isg
    logical :: finishedTrying

    ! start timer
    call g_prof%start("Do time step", g_prof%tmr_do_tstp)

    ! -- By default, the solution groups will be solved once, and
    !    may fail.  But if adaptive stepping is active, then
    !    the solution groups may be solved over and over with
    !    progressively smaller time steps to see if convergence
    !    can be obtained.
    iFailedStepRetry = 0
    retryloop: do

      ! -- idm advance
      call idm_ad()

      do isg = 1, solutiongrouplist%Count()
        sgp => GetSolutionGroupFromList(solutiongrouplist, isg)
        call sgp%sgp_ca()
      end do

      call sim_step_retry(finishedTrying)
      if (finishedTrying) exit retryloop
      iFailedStepRetry = iFailedStepRetry + 1

    end do retryloop

    ! stop timer
    call g_prof%stop(g_prof%tmr_do_tstp)

  end subroutine Mf6DoTimestep

  !> @brief Rerun time step
    !!
    !! This subroutine reruns a single time step for the simulation when
    !! the adaptive time step option is used.
    !!
  !<
  subroutine sim_step_retry(finishedTrying)
    ! -- modules
    use KindModule, only: DP
    use SimVariablesModule, only: lastStepFailed
    use SimModule, only: converge_reset
    use TdisModule, only: kstp, kper, delt, tdis_delt_reset
    use AdaptiveTimeStepModule, only: ats_reset_delt
    ! -- dummy variables
    logical, intent(out) :: finishedTrying !< boolean that indicates if no
    ! additional reruns of the time step are required
    !
    ! -- Check with ats to reset delt and keep trying
    finishedTrying = .true.
    call ats_reset_delt(kstp, kper, lastStepFailed, delt, finishedTrying)
    !
    if (.not. finishedTrying) then
      !
      ! -- Reset delt, which requires updating pertim, totim
      !    and end of period and simulation indicators
      call tdis_delt_reset(delt)
      !
      ! -- Reset state of the simulation convergence flag
      call converge_reset()

    end if
  end subroutine sim_step_retry

  !> @brief Finalize time step
    !!
    !! This function finalizes a single time step for the simulation
    !! and writes output for the time step. Steps include:
    !!   - write output for each model
    !!   - write output for each exchange
    !!   - write output for each solutions
    !!   - perform a final convergence check and whether the simulation
    !!     can continue if convergence was not achieved
    !!
    !! @return  hasConverged   boolean indicating if convergence was achieved for the time step
    !!
  !<
  function Mf6FinalizeTimestep() result(hasConverged)
    ! -- modules
    use KindModule, only: I4B
    use ConstantsModule, only: LINELENGTH, MNORMAL, MVALIDATE
    use ListsModule, only: basesolutionlist, basemodellist, baseexchangelist
    use BaseModelModule, only: BaseModelType, GetBaseModelFromList
    use BaseExchangeModule, only: BaseExchangeType, GetBaseExchangeFromList
    use BaseSolutionModule, only: BaseSolutionType, GetBaseSolutionFromList
    use SimModule, only: converge_check
    use SimVariablesModule, only: isim_mode
    use SourceLoadModule, only: export_post_step
    ! -- return variable
    logical(LGP) :: hasConverged
    ! -- local variables
    class(BaseSolutionType), pointer :: sp => null()
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    class(SpatialModelConnectionType), pointer :: mc => null()
    character(len=LINELENGTH) :: line
    character(len=LINELENGTH) :: fmt
    integer(I4B) :: im
    integer(I4B) :: ix
    integer(I4B) :: ic
    integer(I4B) :: is
    !
    ! -- initialize format and line
    fmt = "(/,a,/)"
    line = 'end timestep'

    ! start timer
    call g_prof%start("Finalize time step", g_prof%tmr_final_tstp)

    !
    ! -- evaluate simulation mode
    select case (isim_mode)
    case (MVALIDATE)
      !
      ! -- Write final message for timestep for each model
      do im = 1, basemodellist%Count()
        mp => GetBaseModelFromList(basemodellist, im)
        call mp%model_message(line, fmt=fmt)
      end do
    case (MNORMAL)

      call g_prof%start("Write output", g_prof%tmr_output)
      !
      ! -- Write output and final message for timestep for each model
      do im = 1, basemodellist%Count()
        mp => GetBaseModelFromList(basemodellist, im)
        call mp%model_ot()
        call mp%model_message(line, fmt=fmt)
      end do
      !
      ! -- Write output for each exchange
      do ix = 1, baseexchangelist%Count()
        ep => GetBaseExchangeFromList(baseexchangelist, ix)
        call ep%exg_ot()
      end do
      !
      ! -- Write output for each connection
      do ic = 1, baseconnectionlist%Count()
        mc => get_smc_from_list(baseconnectionlist, ic)
        call mc%exg_ot()
      end do
      !
      ! -- Write output for each solution
      do is = 1, basesolutionlist%Count()
        sp => GetBaseSolutionFromList(basesolutionlist, is)
        call sp%sln_ot()
      end do
      !
      ! -- update exports
      call g_prof%start("NetCDF export", g_prof%tmr_nc_export)
      call export_post_step()
      call g_prof%stop(g_prof%tmr_nc_export)

      call g_prof%stop(g_prof%tmr_output)
    end select
    !
    ! -- Check if we're done
    call converge_check(hasConverged)

    ! stop timer
    call g_prof%stop(g_prof%tmr_final_tstp)

  end function Mf6FinalizeTimestep

end module Mf6CoreModule
