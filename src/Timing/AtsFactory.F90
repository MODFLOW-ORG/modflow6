module AtsFactoryModule
  use SimModule, only: ustop
  use AdaptiveTimeStepModule, only: AtsType
#if defined(__WITH_MPI__)
  use ParallelAtsModule, only: ParallelAtsType
#endif
  implicit none
  private

  public :: create_ats

contains

!> @brief Factory method to create ATS object, to hide
!< the parallel details for the synchronization across ranks
  function create_ats() result(ats)
    use SimVariablesModule, only: simulation_mode
    class(AtsType), pointer :: ats
    ! local
#if defined(__WITH_MPI__)    
    class(ParallelAtsType), pointer :: par_ats
#endif

    if (simulation_mode == 'SEQUENTIAL') then
      allocate (ats)
#if defined(__WITH_MPI__)
    else if (simulation_mode == 'PARALLEL') then
      allocate (par_ats)
      ats => par_ats
#endif
    else
      call ustop('Unsupported simulation mode for creating ATS: '&
                 &//trim(simulation_mode))
    end if

  end function create_ats

end module AtsFactoryModule
