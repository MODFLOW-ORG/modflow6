!> @brief This module contains the ObsUtilityModule module
!!
!! This module contains subroutines for writing simulated values stored
!! in objects of ObserveType to output files.  The subroutines handle
!! continuous observations, and can write values to either formatted or
!! unformatted files.
!!
!<
module ObsUtilityModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENOBSNAME, LENBIGLINE
  use ObserveModule, only: ObserveType
  use ObsOutputListModule, only: ObsOutputListType
  use ObsOutputModule, only: ObsOutputType
  use TdisModule, only: totim

  implicit none

  private
  public :: write_fmtd_obs, write_unfmtd_obs

contains

  !> @ brief Write formatted observation
  !!
  !!  Subroutine to write observation data for the end of a time step to
  !!  a formatted file. If the simulation time has not been written to
  !!  for the current time step, totim is written. The simulated value is
  !!  written in the format specified in the fmtc argument.
  !!
  !<
  subroutine write_fmtd_obs(fmtc, obsrv, obsOutputList, value)
    ! -- dummy
    character(len=*), intent(in) :: fmtc !< observation format
    type(ObserveType), intent(inout) :: obsrv !< observation type
    type(ObsOutputListType), pointer, intent(inout) :: obsOutputList !< observation list
    real(DP), intent(in) :: value !< observation
    ! -- local
    integer(I4B) :: indx
    integer(I4B) :: nunit
    character(len=20) :: ctotim
    character(len=50) :: cval
    type(ObsOutputType), pointer :: ObsOutput => null()
    ! -- output unit
    nunit = obsrv%UnitNumber
    !
    indx = obsrv%indxObsOutput
    ObsOutput => obsOutputList%Get(indx)
    if (obsOutput%empty_line) then
      ObsOutput%empty_line = .FALSE.
      write (ctotim, '(G20.13)') totim
    else
      ctotim = ''
    end if
    ! -- append value to output line
    write (cval, fmtc) value
    write (nunit, '(3a)', advance='NO') &
      trim(adjustl(ctotim)), ',', trim(adjustl(cval))
    !
    ! -- flush the file
    !    Added flush after each non-advancing write to resolve
    !    issue with ifort (IFORT) 19.1.0.166 20191121 for Linux
    !    that occurred on some Linux systems.
    flush (nunit)
  end subroutine write_fmtd_obs

  !> @ brief Write unformatted observation
  !!
  !!  Subroutine to write observation data for the end of a time step to
  !!  a unformatted file. If the simulation time has not been written for
  !!  the current time step, totim is written. The simulated value is
  !!  written using the precision specified in the iprec argument.
  !!
  !!  iprec = 1: real32 specifies 32-bit real = 4 bytes = single precision.
  !!  iprec = 2: real64 specifies 64-bit real = 8 bytes = double precision.
  !!
  !<
  subroutine write_unfmtd_obs(obsrv, iprec, obsOutputList, value)
    use iso_fortran_env, only: real32, real64
    ! -- dummy
    type(ObserveType), intent(inout) :: obsrv !< observation type
    integer(I4B), intent(in) :: iprec !< observation precision
    type(ObsOutputListType), pointer, intent(inout) :: obsOutputList !< observation list
    real(DP), intent(in) :: value !< observation
    ! -- local
    integer(I4B) :: indx, nunit
    real(real32) :: totimsngl, valsngl
    real(real64) :: totimdbl, valdbl
    type(ObsOutputType), pointer :: obsOutput => null()
    !
    ! -- output unit
    nunit = obsrv%UnitNumber
    ! -- continuous observation
    indx = obsrv%indxObsOutput
    obsOutput => obsOutputList%Get(indx)
    if (obsOutput%empty_line) then
      obsOutput%empty_line = .FALSE.
      if (iprec == 1) then
        totimsngl = real(totim, real32)
        write (nunit) totimsngl
      elseif (iprec == 2) then
        totimdbl = totim
        write (nunit) totimdbl
      end if
    end if
    ! -- write value to unformatted output
    if (iprec == 1) then
      valsngl = real(value, real32)
      write (nunit) valsngl
    elseif (iprec == 2) then
      valdbl = value
      write (nunit) valdbl
    end if
  end subroutine write_unfmtd_obs

end module ObsUtilityModule
