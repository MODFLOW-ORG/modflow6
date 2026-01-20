!> @brief This module contains the derived types ObserveType and ObsDataType
!!
!! This module contains the derived types ObserveType and ObsDataType.
!!
!! - ObserveType -- is designed to contain all information and
!!   functionality needed for one observation. ObserveType contains a
!!   pointer to an ObsDataType object.
!!
!! - ObsDataType -- is for storing package ID, observation type, and a
!!   pointer to a subroutine that will be called to process the IDstring
!!   provided in Obs input.  The ProcessIdPtr member of ObsDataType
!!   requires a pointer to an ObserveType object.
!!
!<
module ObserveModule

  use KindModule, only: DP, I4B
  use BaseDisModule, only: DisBaseType
  use ConstantsModule, only: LENBOUNDNAME, LENOBSNAME, LENOBSTYPE, &
                             MAXOBSTYPES, DNODATA, DZERO
  use TableModule, only: TableType
  use InputOutputModule, only: urword
  use ListModule, only: ListType
  use SimModule, only: store_warning, store_error, &
                       store_error_unit
  use TdisModule, only: totim, totalsimtime
  use ArrayHandlersModule, only: ExpandArrayWrapper

  implicit none

  private
  public :: ObserveType, ConstructObservation, ObsDataType, &
            AddObsToList, GetObsFromList

  type :: ObserveType
    ! -- Public members
    !
    ! -- For all observations
    integer(I4B), public :: NodeNumber = 0 !< observation node number
    integer(I4B), public :: UnitNumber = 0 !< observation output unit number
    character(len=LENOBSNAME), public :: Name = '' !< observation name
    character(len=LENOBSTYPE), public :: ObsTypeId = '' !< observation type id
    character(len=200), public :: IDstring = '' !< observation id string
    character(len=LENBOUNDNAME), public :: FeatureName = '' !< observation feature name
    character(len=LENBOUNDNAME), public :: FeatureName2 = '' !< observation feature name 2
    !
    ! -- members specific to NPF intercell-flow observations
    integer(I4B), public :: NodeNumber2 = 0 !< observation second nod number
    integer(I4B), public :: JaIndex = -2 !< observation JA index
    !
    ! -- members that can be used as needed by packages or models
    integer(I4B), public :: intPak1 = 0 !<
    real(DP), public :: Obsdepth = DZERO !<
    real(DP), public :: dblPak1 = DZERO !<
    !
    ! -- indxbnds is intended to hold indices of position(s) in bound
    !    array of boundaries included in the observation.
    integer(I4B), public :: indxbnds_count = 0 !< number of observations indexes when using boundname
    integer(I4B), allocatable, dimension(:), public :: indxbnds !< node numbers for observations when using boundname
    !
    ! -- Set FormattedOutput false if output unit is opened for unformatted i/o
    logical, public :: FormattedOutput = .true. !< logical indicating if observation output is formatted
    logical, public :: BndFound = .false. !< logical indicating if a boundname was found
    real(DP), public :: CurrentTimeStepEndValue = DZERO !< observation value
    real(DP), public :: CurrentTimeStepEndTime = DZERO !< observation time
    !
    ! -- Members specific to continuous observations
    integer(I4B), public :: indxObsOutput = -1 !< index for observation output
    !
    ! -- Private members
    type(ObsDataType), pointer, private :: obsDatum => null() !< observation Datum
  contains
    ! -- Public procedures
    procedure, public :: ResetCurrentValue
    procedure, public :: WriteTo
    procedure, public :: AddObsIndex
    procedure, public :: ResetObsIndex
    procedure, public :: da
  end type ObserveType

  type :: ObsDataType
    ! -- Public members
    character(len=LENOBSTYPE), public :: ObsTypeID = '' !< observation type id
    logical, public :: Cumulative = .false. !< logical indicating if observations should be summed
    procedure(ProcessIdSub), nopass, pointer, public :: ProcessIdPtr => null() !< process id pointer
  end type ObsDataType

  abstract interface

    !> @ brief Process user-provided IDstring
  !!
  !!  Subroutine that processes the user-provided IDstring, which identifies
  !!  the grid location or model feature to be observed.
  !!
    !<
    subroutine ProcessIdSub(obsrv, dis, inunitobs, iout)
      use KindModule, only: DP, I4B
      import :: ObserveType
      import :: DisBaseType
      ! -- dummy
      type(ObserveType), intent(inout) :: obsrv !< observation type
      class(DisBaseType), intent(in) :: dis !< discretization object
      integer(I4B), intent(in) :: inunitobs !< observation input file unit
      integer(I4B), intent(in) :: iout !< model list file unit
    end subroutine ProcessIdSub
  end interface

contains

  ! Procedures bound to ObserveType

  !> @ brief Reset current observation value
  !!
  !!  Subroutine to reset the current observation value.
  !!
  !<
  subroutine ResetCurrentValue(this)
    ! -- dummy
    class(ObserveType), intent(inout) :: this
    !
    ! -- Reset current value to zero.
    this%CurrentTimeStepEndValue = DZERO
  end subroutine ResetCurrentValue

  !> @ brief Write observation input data
  !!
  !!  Subroutine to write observation input data to a table in the model
  !!  list file.
  !!
  !<
  subroutine WriteTo(this, obstab, btagfound, fnamein)
    ! -- dummy
    class(ObserveType), intent(inout) :: this
    type(TableType), intent(inout) :: obstab !< observation table
    character(len=*), intent(in) :: btagfound !< logical indicating if boundname was found
    character(len=*), intent(in) :: fnamein !< observation input file name
    ! -- local
    character(len=12) :: tag
    character(len=80) :: fnameout
    !
    ! -- write btagfound to tag
    if (len_trim(btagfound) > 12) then
      tag = btagfound(1:12)
    else
      write (tag, '(a12)') btagfound
    end if
    !
    ! -- write fnamein to fnameout
    if (len_trim(fnamein) > 80) then
      fnameout = fnamein(1:80)
    else
      write (fnameout, '(a80)') fnamein
    end if
    !
    ! -- write data to observation table
    call obstab%add_term(this%Name)
    call obstab%add_term(tag//trim(this%ObsTypeId))
    call obstab%add_term('ALL TIMES')
    call obstab%add_term('"'//trim(this%IDstring)//'"')
    call obstab%add_term(fnameout)
  end subroutine WriteTo

  !> @ brief Reset a observation index
  !!
  !!  Subroutine to reset the observation index count and array.
  !!
  !<
  subroutine ResetObsIndex(this)
    ! -- dummy
    class(ObserveType), intent(inout) :: this
    !
    ! -- Reset the index count
    this%indxbnds_count = 0
    !
    ! -- Deallocate observation index array, if necessary
    if (allocated(this%indxbnds)) then
      deallocate (this%indxbnds)
    end if
    !
    ! -- Allocate observation index array to size 0
    allocate (this%indxbnds(0))
  end subroutine ResetObsIndex

  !> @ brief Add a observation index
  !!
  !!  Subroutine to add the observation index to the observation index
  !!  array (indxbnds). The observation index count (indxbnds_count) is
  !!  also incremented by one and the observation index array is
  !!  expanded, if necessary.
  !!
  !<
  subroutine AddObsIndex(this, indx)
    ! -- dummy
    class(ObserveType), intent(inout) :: this
    integer(I4B), intent(in) :: indx !< observation index
    !
    ! -- Increment the index count
    this%indxbnds_count = this%indxbnds_count + 1
    !
    ! -- Expand the observation index array, if necessary
    call ExpandArrayWrapper(this%indxbnds_count, this%indxbnds, loginc=.TRUE.)
    !
    ! -- add index to observation index
    this%indxbnds(this%indxbnds_count) = indx
  end subroutine AddObsIndex

  !> @ brief Deallocate a observation
  !!
  !!  Subroutine to deallocated a observation (ObserveType).
  !!
  !<
  subroutine da(this)
    ! -- dummy
    class(ObserveType), intent(inout) :: this
    if (allocated(this%indxbnds)) then
      deallocate (this%indxbnds)
    end if
  end subroutine da

  ! Non-type-bound procedures

  !> @ brief Construct a new ObserveType
  !!
  !!  Subroutine to construct and return an ObserveType object based
  !!  on the contents of defLine.
  !!
  !<
  subroutine ConstructObservation(newObservation, defLine, numunit, &
                                  formatted, indx, obsData, inunit)
    ! -- dummy variables
    type(ObserveType), pointer :: newObservation !< new ObserveType
    character(len=*), intent(in) :: defLine !< string with observation data
    integer(I4B), intent(in) :: numunit !< Output unit number
    logical, intent(in) :: formatted !< logical indicating if formatted output will be written
    integer(I4B), intent(in) :: indx !< Index in ObsOutput array
    type(ObsDataType), dimension(:), pointer, intent(in) :: obsData !< obsData type
    integer(I4B), intent(in) :: inunit !< observation input file unit
    ! -- local
    real(DP) :: r
    integer(I4B) :: i
    integer(I4B) :: icol
    integer(I4B) :: iout
    integer(I4B) :: istart
    integer(I4B) :: istop
    integer(I4B) :: n
    !
    ! -- initialize
    iout = 0
    icol = 1
    !
    ! -- Allocate an ObserveType object.
    allocate (newObservation)
    allocate (newObservation%indxbnds(0))
    !
    ! -- Set indxbnds_count to 0
    newObservation%indxbnds_count = 0
    !
    ! -- Define the contents of the ObservationSingleType object based on the
    !    contents of defLine.
    !
    ! -- Get observation name and store it
    call urword(defLine, icol, istart, istop, 1, n, r, iout, inunit)
    newObservation%Name = defLine(istart:istop)
    !
    ! -- Get observation type, convert it to uppercase, and store it.
    call urword(defLine, icol, istart, istop, 1, n, r, iout, inunit)
    newObservation%ObsTypeId = defLine(istart:istop)
    !
    ! -- Look up package ID for this observation type and store it
    do i = 1, MAXOBSTYPES
      if (obsData(i)%ObsTypeID == newObservation%ObsTypeId) then
        newObservation%obsDatum => obsData(i)
        exit
      elseif (obsData(i)%ObsTypeID == '') then
        exit
      end if
    end do
    !
    ! -- Remaining text is ID [and ID2]; store the remainder of the string
    istart = istop + 1
    istop = len_trim(defLine)
    if (istart > istop) then
      istart = istop
    end if
    newObservation%IDstring = defLine(istart:istop)
    !
    ! Store UnitNumber, FormattedOutput, and IndxObsOutput
    newObservation%UnitNumber = numunit
    newObservation%FormattedOutput = formatted
    newObservation%IndxObsOutput = indx
  end subroutine ConstructObservation

  !> @ brief Cast a object as a ObserveType
  !!
  !!  Function to cast an object as a ObserveType object.
  !!
  !<
  function CastAsObserveType(obj) result(res)
    ! -- dummy
    class(*), pointer, intent(inout) :: obj !< object
    ! -- return
    type(ObserveType), pointer :: res !< returned ObserveType object
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    type is (ObserveType)
      res => obj
    end select
  end function CastAsObserveType

  !> @ brief Add a ObserveType to a list
  !!
  !!  Subroutine to add a ObserveType to a list.
  !!
  !<
  subroutine AddObsToList(list, obs)
    ! -- dummy
    type(ListType), intent(inout) :: list !< ObserveType list
    type(ObserveType), pointer, intent(inout) :: obs !< ObserveType
    ! -- local
    class(*), pointer :: obj
    !
    obj => obs
    call list%Add(obj)
  end subroutine AddObsToList

  !> @ brief Get an ObserveType from a list
  !!
  !!  Function to get an ObserveType from a list.
  !!
  !<
  function GetObsFromList(list, idx) result(res)
    ! -- dummy
    type(ListType), intent(inout) :: list !< ObserveType list
    integer(I4B), intent(in) :: idx !< ObserveType list index
    ! -- return
    type(ObserveType), pointer :: res !< returned ObserveType
    ! -- local
    class(*), pointer :: obj
    !
    obj => list%GetItem(idx)
    res => CastAsObserveType(obj)
  end function GetObsFromList

end module ObserveModule
