!> @brief This module defines the derived type ObsOutputListType
!!
!! This module contains a list of ObsOutputType objects and
!! methods needed for coordinating between an ObsType object and its
!! ObsOutputType objects.  Like ObsOutputType, ObsOutputListType is
!! needed only for processing continuous observations.
!!
!<
module ObsOutputListModule

  use KindModule, only: DP, I4B
  use InputOutputModule, only: same_word
  use ListModule, only: ListType
  use ObsOutputModule, only: ObsOutputType, ConstructObsOutput, &
                             AddObsOutputToList, GetObsOutputFromList

  implicit none

  private
  public :: ObsOutputListType

  type ObsOutputListType
    ! -- Private members
    type(ListType), private :: ObsOutputs
  contains
    ! -- Public procedures
    procedure, public :: Add
    procedure, public :: ResetAllObsEmptyLines
    procedure, public :: ContainsFile
    procedure, public :: Count
    procedure, public :: Get
    procedure, public :: WriteAllObsLineReturns
    procedure, public :: Clear
    procedure, public :: DeallocObsOutputList
  end type ObsOutputListType

contains

  !> @ brief Reset empty line logical for all observations
  !!
  !!  Subroutine to reset the empty line logical for all ObsOutputType
  !!  objects in the list.
  !!
  !<
  subroutine ResetAllObsEmptyLines(this)
    ! -- dummy
    class(ObsOutputListType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i, num
    type(ObsOutputType), pointer :: obsOutput => null()
    !
    num = this%Count()
    do i = 1, num
      obsOutput => this%Get(i)
      call obsOutput%ResetObsEmptyLine()
    end do
  end subroutine ResetAllObsEmptyLines

  !> @ brief Count the number of ObsOutputType objects
  !!
  !!  Subroutine to return the number of ObsOutputType objects in the list.
  !!
  !<
  function Count(this)
    ! -- return
    integer(I4B) :: count !< number of ObsOutputType objects
    ! -- dummy
    class(ObsOutputListType), intent(inout) :: this
    !
    Count = this%ObsOutputs%Count()
  end function Count

  !> @ brief Determine if a file name is in the list of ObsOutputType objects
  !!
  !!  Function to determine if a file name is in the list of
  !!  ObsOutptType objects.
  !!
  !<
  logical function ContainsFile(this, fname)
    ! -- dummy
    class(ObsOutputListType), intent(inout) :: this
    character(len=*), intent(in) :: fname !< observation output file name
    ! -- local
    type(ObsOutputType), pointer :: obsOutput => null()
    integer(I4B) :: i, n
    !
    ContainsFile = .false.
    n = this%Count()
    loop1: do i = 1, n
      obsOutput => this%Get(i)
      if (same_word(obsOutput%filename, fname)) then
        ContainsFile = .true.
        exit loop1
      end if
    end do loop1
  end function ContainsFile

  !> @ brief Add a ObsOutputType object to the list
  !!
  !!  Subroutine to add a new ObsOutputType object to the ObsOutputList and
  !!  assign ObsOutputType members.
  !!
  !<
  subroutine Add(this, fname, nunit)
    ! -- dummy
    class(ObsOutputListType), intent(inout) :: this
    character(len=*), intent(in) :: fname !< observation output file name
    integer(I4B), intent(in) :: nunit !< observation output unit number
    ! -- local
    type(ObsOutputType), pointer :: obsOutput => null()
    !
    call ConstructObsOutput(obsOutput, fname, nunit)
    call AddObsOutputToList(this%ObsOutputs, obsOutput)
  end subroutine Add

  !> @ brief Write line returns for all ObsOutputListType
  !!
  !!  Subroutine to write line returns for a time step for all observation
  !!  output files in a ObsOutputListType.
  !!
  !<
  subroutine WriteAllObsLineReturns(this)
    ! -- dummy
    class(ObsOutputListType), intent(inout) :: this
    ! -- local
    type(ObsOutputType), pointer :: obsOutput => null()
    integer(I4B) :: i, num
    !
    num = this%Count()
    do i = 1, num
      obsOutput => this%Get(i)
      if (obsOutput%FormattedOutput) then
        call obsOutput%WriteObsLineReturn()
      end if
    end do
  end subroutine WriteAllObsLineReturns

  !> @ brief Get an item from a ObsOutputListType
  !!
  !!  Function to get a ObsOutputType from a ObsOutputListType list.
  !!
  !<
  function Get(this, indx) result(obsOutput)
    ! -- dummy
    class(ObsOutputListType), intent(inout) :: this
    integer(I4B), intent(in) :: indx !< index for ObsOutputType object
    ! result
    type(ObsOutputType), pointer :: obsOutput
    !
    obsOutput => GetObsOutputFromList(this%ObsOutputs, indx)
  end function Get

  !> @ brief Clear a ObsOutputListType
  !!
  !!  Subroutine to clear a ObsOutputListType list.
  !!
  !<
  subroutine Clear(this)
    ! -- dummy
    class(ObsOutputListType), intent(inout) :: this
    !
    call this%ObsOutputs%Clear()
  end subroutine Clear

  !> @ brief Deallocate a ObsOutputListType
  !!
  !!  Subroutine to deallocate a ObsOutputListType list.
  !!
  !<
  subroutine DeallocObsOutputList(this)
    ! -- dummy
    class(ObsOutputListType), intent(inout) :: this
    ! -- local
    integer :: i, n
    type(ObsOutputType), pointer :: obsoutput => null()
    !
    n = this%Count()
    do i = 1, n
      obsoutput => GetObsOutputFromList(this%ObsOutputs, i)
      !call obsoutput%DeallocObsOutput()
    end do
    !
    call this%ObsOutputs%Clear(.true.)
  end subroutine DeallocObsOutputList

end module ObsOutputListModule
