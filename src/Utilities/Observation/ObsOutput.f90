!> @brief This module defines the derived type ObsOutputType
!!
!! This module contains information and methods needed for writing
!! a line of simulated values for observations to an output file.  Each
!! block of type continuous in an observation file is
!! associated with an ObsOutputType object. However, the methods are
!! needed only for continuous observations.
!!
!<
module ObsOutputModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENBIGLINE, LENOBSNAME
  use ListModule, only: ListType

  implicit none

  private
  public :: ObsOutputType, ConstructObsOutput, AddObsOutputToList, &
            GetObsOutputFromList

  type :: ObsOutputType
    ! -- Public members
    ! kind specified to ensure consistent binary output
    integer(kind=4), public :: nobs = 0 !< number of observations
    integer(I4B), public :: nunit = 0 !< observation output unit
    character(len=500), public :: filename = '' !< observation output filename
    logical(LGP), public :: empty_line = .TRUE. !< logical indicating if the line for a time step is empty
    character(len=LENOBSNAME), public :: header = '' !< observation header string
    logical, public :: FormattedOutput = .true. !< logical indicating if writing formatted output
  contains
    ! -- Public procedures
    procedure, public :: ResetObsEmptyLine
    procedure, public :: WriteObsLineReturn
  end type ObsOutputType

contains

  ! Procedures bound to ObsOutputType

  !> @ brief Reset empty line logical
  !!
  !!  Subroutine to reset the empty line logical.
  !!
  !<
  subroutine ResetObsEmptyLine(this)
    ! -- dummy
    class(ObsOutputType), intent(inout) :: this
    !
    this%empty_line = .TRUE.
  end subroutine ResetObsEmptyLine

  !> @ brief Write line return for observation
  !!
  !!  Subroutine to write a line return for a time step in an observation
  !!  output file.
  !!
  !<
  subroutine WriteObsLineReturn(this)
    ! -- dummy
    class(ObsOutputType), intent(inout) :: this
    ! -- write a line return to end of observation output line
    !    for this totim
    write (this%nunit, '(a)', advance='YES') ''
  end subroutine WriteObsLineReturn

  ! Non-type-bound procedures

  !> @ brief Cast as ObsOutputType
  !!
  !!  Cast an object as an ObsOutputType.
  !!
  !<
  function CastAsObsOutputType(obj) result(res)
    ! -- dummy
    class(*), pointer, intent(inout) :: obj !< input object
    type(ObsOutputType), pointer :: res !< ObsOutputType
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    type is (ObsOutputType)
      res => obj
    class default
      continue
    end select
  end function CastAsObsOutputType

  !> @ brief Construct and assign ObsOutputType object
  !!
  !!  Subroutine to construct an ObsOutputType object and assign
  !!  the observation output file name and unit number.
  !!
  !<
  subroutine ConstructObsOutput(newObsOutput, fname, nunit)
    ! -- dummy
    type(ObsOutputType), pointer, intent(out) :: newObsOutput
    character(len=*), intent(in) :: fname !< observation output file name
    integer(I4B), intent(in) :: nunit !< observation output unit number
    !
    allocate (newObsOutput)
    newObsOutput%filename = fname
    newObsOutput%nunit = nunit
  end subroutine ConstructObsOutput

  !> @ brief Add observation output to a list
  !!
  !!  Subroutine to add observation output to a observation list.
  !!
  !<
  subroutine AddObsOutputToList(list, obsOutput)
    ! -- dummy
    type(ListType), intent(inout) :: list !< observation list
    type(ObsOutputType), pointer, intent(inout) :: obsOutput !< observation output
    ! -- local
    class(*), pointer :: obj
    !
    obj => obsOutput
    call list%Add(obj)
  end subroutine AddObsOutputToList

  !> @ brief Get observation output from a list
  !!
  !!  Subroutine to get observation output from a observation list.
  !!
  !<
  function GetObsOutputFromList(list, idx) result(res)
    implicit none
    ! -- dummy
    type(ListType), intent(inout) :: list !< observation list
    integer(I4B), intent(in) :: idx !< observation index
    type(ObsOutputType), pointer :: res !< observation output
    ! -- local
    class(*), pointer :: obj
    !
    obj => list%GetItem(idx)
    res => CastAsObsOutputType(obj)
  end function GetObsOutputFromList

end module ObsOutputModule
