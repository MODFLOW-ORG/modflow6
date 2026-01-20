! This module defines two derived types:
!
! ObserveType -- is designed to contain all information and
!   functionality needed for one observation.  The observation can be
!   of either "single" or "continuous" type.  ObserveType contains a 
!   pointer to an ObsDataType object.
!
! ObsDataType -- is for storing package ID, observation type, and a
!   pointer to a subroutine that will be called to process the IDstring
!   provided in Obs input.  The ProcessIdPtr member of ObsDataType
!   requires a pointer to an ObserveType object.
!-----------------------------------------------------------------------
module ObserveModule

  use ConstantsModule,        only: DONE, DZERO, LENOBSNAME,  &
                                    LENOBSTYPE, MAXCHARLEN
  use ConstantsPHMFModule,    only: LENOBSNAMENEW, HUGEDBL, HDRYDEFAULT
  use MathUtilModule,         only: is_close
  use ListModule,             only: ListType
  use SimPHMFModule,          only: store_warning, store_error, &
                                    store_error_unit, ustop

  implicit none

  private
  public :: ObserveType, CastAsObserveType, ConstructObservation, &
            ObsDataType, obs_type, AddObserveToList, GetObserveFromList
            
  character(len=10), dimension(2) :: obs_type = (/'SINGLE    ','CONTINUOUS'/)

  type :: ObserveType
    ! -- Type for all observations, either single or continuous
    integer                    :: isorc = 0 ! 1=single, 2=continuous
    character(len=LENOBSNAME)  :: Name = ''
    character(len=LENOBSTYPE)  :: ObsTypeId = ''
    type(ObsDataType), pointer :: obsDatum => null()
    character(len=200) :: IDstring = ''
    integer            :: NodeNumber = 0
    integer            :: UnitNumber = 0
    ! -- members specific to NPF intercell-flow observations
    integer            :: NodeNumber2 = 0
    integer            :: JaIndex = -2
    ! -- members that can be used as needed by packages or models
    integer            :: intPak1 = 0
    integer            :: intPak2 = 0
    double precision   :: Hdry = hdrydefault
    double precision   :: Obsdepth = 0.0d0
    double precision   :: dblPak1 = 0.0d0
    ! -- indxbnds is intended to hold indices of position(s) in bound
    !    array of boundaries included in the observation.
    integer, allocatable, dimension(:) :: indxbnds
    ! Set FormattedOutput false if output unit is opened for unformatted i/o
    logical            :: FormattedOutput = .true.
    ! -- Members specific to single observations
    double precision   :: ObsTime = 0.0d0   ! Modflow time of observation
    double precision   :: CurrentTimeStepEndValue = 0.0d0
    double precision   :: CurrentTimeStepEndTime = 0.0d0
    double precision   :: PrecedingTimeStepEndValue = 0.0d0
    double precision   :: PrecedingTimeStepEndTime = 0.0d0
    double precision   :: SimulatedValue = 0.0d0
    ! -- Done is set true when time step containing a single observation
    !    has been simulated. Written is set true when the simulated value
    !    for a single observation has been written.
    logical            :: Done = .false.
    logical            :: Written = .false.
    logical            :: Summ = .false.
    ! -- Members specific to continuous observations
    integer            :: indxObsOutput = -1
    ! -- For PostObsMF, to support combining continuous observations
    integer :: ksrc
    character(len=LENOBSNAME), allocatable, dimension(:) :: srcobsnames
    double precision, allocatable, dimension(:) :: weights
    double precision, allocatable, dimension(:,:) :: srcvals ! (ntimes,nsrc)
    ! -- For PreHeadsMF
    ! Ned todo: may not need Xcoord and Ycoord?
    ! Instead, need xoff & yoff and indices for adjacent cells.
    ! Also need DELR and DELC of adjacent cells.
    character(len=LENOBSNAMENEW)  :: NameP = ''
    character(len=LENOBSNAMENEW)  :: NameR = ''
    character(len=LENOBSNAMENEW)  :: NameC = ''
    character(len=LENOBSNAMENEW)  :: NameD = ''
    double precision :: xoff = DZERO
    double precision :: yoff = DZERO
    double precision :: delradj = DZERO
    double precision :: delcadj = DZERO
    ! Contributing factors
    double precision :: cfp = DONE
    double precision :: cfr = DZERO
    double precision :: cfc = DZERO
    double precision :: cfd = DZERO
    ! Grid coordinates of observation point
    double precision :: gridX = DZERO
    double precision :: gridY  = DZERO
    integer :: irow = 0     ! row index, primary cell
    integer :: jcol = 0     ! column index, primary cell
    integer :: jcoladj = 0  ! column injex, adjacent column
    integer :: irowadj = 0  ! row index, adjacent row
    integer :: idiag = 0    ! = 1 if head in diagonally adjacent cell is needed
    integer :: Layer
  contains
    procedure, public :: CalcSimVal
    procedure, public :: WriteTo
    procedure, public :: WriteGroup
  end type ObserveType

  type :: ObsDataType
    character(len=LENOBSTYPE)                :: ObsTypeID = ''
    logical                                  :: Cumulative = .false.
  end type ObsDataType

  interface ConstructObservation
    module procedure ConstructObservationFromLine, &
                     ConstructObservationSimple
  end interface ConstructObservation

contains

  ! Procedures bound to ObserveType

  subroutine WriteTo(this, iout)
! **************************************************************************
! WriteTo -- Write information about this observation to table in list file.
! **************************************************************************
!
!    SPECIFICATIONS:
! --------------------------------------------------------------------------
    implicit none
    ! -- dummy arguments
    class(ObserveType), intent(inout) :: this
    integer, intent(in) :: iout
    ! formats
    10 format(a,2x,a,a,t41,G12.5,      t55,'"',a,'"')
    20 format(a,2x,a,a,t42,'All times',t55,'"',a,'"')
    !
    select case (this%isorc)
    case (1)
      write(iout,10)this%Name, 'Single ', this%ObsTypeId, &
                    this%ObsTime, trim(this%IDstring)
    case (2)
      write(iout,20)this%Name, 'Continuous ', this%ObsTypeId, &
                    trim(this%IDstring)
    end select
    return
  end subroutine WriteTo

  subroutine WriteGroup(this, iout)
    ! Write lines to be read by PostObsMF to define one combined observation.
    ! Depending on the need for contributions from adjacent cells, from 1 to 4
    ! source observations will be needed to define the combined observation.
    ! dummy
    class(ObserveType), intent(inout) :: this
    ! Iout needs to be for the PostObsMF input file being created.
    integer,            intent(in) :: iout
    ! local
    ! formats
    10 format(2x,'NEW',2x,a)
    !
    ! Ned todo: Call routine to calculate weights using bilinear interpolation.
    !
    ! Write NEW line. Outname is the same as source name for primary cell.
    write(iout,10)trim(this%Name)
    !
    ! Always write source line for primary cell
    !
    ! Optionally write source line for cell in adjacent row
    if (this%irowadj > 0) then
    endif
    !
    ! Optionally write source line for cell in adjacent column
    if  (this%jcoladj > 0) then
    endif
    !
    ! Optionally write source line for cell in diagonal position
    if (this%idiag == 1) then
    endif
    !
    return
  end subroutine WriteGroup

  subroutine CalcSimVal(this, itime)
    ! Return combined simulated value, recalculating weights as needed to 
    ! account for source simulated values that equal HDRY.
    ! dummy
    class(ObserveType) :: this
    integer, intent(in) :: itime
    ! local
    integer :: i, k, nsrc
    double precision :: simval
    double precision :: factor, sumweights
    double precision, allocatable, dimension(:) :: weights
    !
    simval = DZERO
    nsrc = size(this%srcobsnames)
    !
    if (this%hdry == hugedbl) then
      ! HDRY has not been specified, so don't bother with recalculating weights
      do i=1,nsrc
        simval = simval + this%srcvals(itime, i) * this%weights(i)
      enddo
      this %SimulatedValue = simval
    else
      ! allocate local array
      allocate(weights(nsrc))
      !
      ! store local copy of weights, and count values that equal HDRY
      sumweights = DZERO
      k = 0
      do i=1,nsrc
        if (is_close(this%srcvals(itime, i), this%hdry)) then
          k = k + 1
          weights(i) = DZERO
        else
          weights(i) = this%weights(i)
          sumweights = sumweights + weights(i)
        endif
      enddo
      ! If SUM has not been specified and any value = HDRY,
      ! recalculate weights. No need to revise weight values if 
      ! SUM has been specified.
      if (k > 0 .and. .not. this%Summ) then
        factor = DONE / sumweights
        do i=1,nsrc
          weights(i) = weights(i) * factor
        enddo
      endif
      ! Calculate simulated value using revised weights
      do i=1,nsrc
        simval = simval + this%srcvals(itime, i) * weights(i)
      enddo
      this %SimulatedValue = simval
      !
      ! deallocate local array
      deallocate(weights)
    endif
    !
    return
  end subroutine CalcSimVal

  ! Non-type-bound procedures

  function CastAsObserveType(obj) result(res)
    implicit none
    class(*), pointer, intent(inout) :: obj
    type(ObserveType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    type is (ObserveType)
      res => obj
    end select
    return
  end function CastAsObserveType

  subroutine ConstructObservationFromLine(newObservation, isorc, defLine, &
                                  numunit, formatted, indx)
! **************************************************************************
! ConstructObservationFromLine -- Construct and return an ObserveType
! object based on the contents of defLine.
! **************************************************************************
!
!    SPECIFICATIONS:
! --------------------------------------------------------------------------
    use InputOutputModule, only: urword
    use SimPHMFModule, only: store_error, ustop
    implicit none
    ! -- dummy variables
    type(ObserveType), pointer :: newObservation
    integer,           intent(in) :: isorc
    character(len=*),  intent(in) :: defLine
    integer,           intent(in) :: numunit   ! Output unit number
    logical,           intent(in) :: formatted ! Formatted output?
    integer,           intent(in) :: indx      ! Index in ObsOutput array
    ! -- local variables
    double precision :: r
    integer :: icol, inunit, iout, istart, istop, ltrim, n
    ! -- formats
    10 format('Time for observation "',a, &
              '" exceeds simulation time.')
    ! --------------------------------------------------------------------------
    !
    ! -- initialize
    inunit = 0
    iout = 0
    icol = 1
    !
    ! -- Allocate an ObserveType object.
    allocate(newObservation)
    allocate(newObservation%indxbnds(0))
    newObservation%isorc = isorc
    newObservation%ksrc = 0
    allocate(newObservation%srcobsnames(0))
    allocate(newObservation%weights(0))
    !
    ! -- Define the contents of the ObservationSingleType object based on the
    !    contents of defLine.
    !
    if (defLine /= '') then
      ! -- Get observation name and store it
      call urword(defLine,icol,istart,istop,1,n,r,iout,inunit)
      newObservation%Name = defLine(istart:istop)
      !
      ! -- Get observation type, convert it to uppercase, and store it.
      call urword(defLine,icol,istart,istop,1,n,r,iout,inunit)
      newObservation%ObsTypeId = defLine(istart:istop)
      !
      ! -- Remaining text is ID [and ID2]; store it
      ltrim = len_trim(defLine)
      call urword(defLine,icol,istart,istop,1,n,r,iout,inunit)
      newObservation%IDstring = (defLine(istart:ltrim))
    endif
    !
    ! Store UnitNumber, FormattedOutput, and IndxObsOutput
    newObservation%UnitNumber = numunit
    newObservation%FormattedOutput = formatted
    newObservation%IndxObsOutput = indx
    !
    return
  end subroutine ConstructObservationFromLine

  subroutine ConstructObservationSimple(newObservation, isorc, &
               time, layer, numunit, formatted)
! **************************************************************************
! ConstructObservationWithCoords -- Construct and return an ObserveType
! object with X and Y coordinates.
! **************************************************************************
!
!    SPECIFICATIONS:
! --------------------------------------------------------------------------
    use InputOutputModule, only: urword
    use SimPHMFModule, only: store_error, ustop
    implicit none
    ! -- dummy variables
    type(ObserveType), pointer :: newObservation
    integer,           intent(in) :: isorc
    double precision,  intent(in) :: time
    integer,           intent(in) :: layer
    integer,           intent(in) :: numunit   ! Output unit number
    logical,           intent(in) :: formatted ! Formatted output?
    ! -- local variables
    integer :: icol, inunit, iout
    ! -- formats
    10 format('Time for observation "',a, &
              '" exceeds simulation time.')
    ! --------------------------------------------------------------------------
    !
    ! -- initialize
    inunit = 0
    iout = 0
    icol = 1
    !
    ! -- Allocate an ObserveType object.
    allocate(newObservation)
    allocate(newObservation%indxbnds(0))
    newObservation%isorc = isorc
    newObservation%ksrc = 0
    allocate(newObservation%srcobsnames(0))
    allocate(newObservation%weights(0))
    !
    newObservation%ObsTime = time
    newObservation%Layer = layer
    !
    ! Store UnitNumber, FormattedOutput, and IndxObsOutput
    newObservation%UnitNumber = numunit
    newObservation%FormattedOutput = formatted
!    newObservation%IndxObsOutput = indx
    !
    return
  end subroutine ConstructObservationSimple

  subroutine AddObserveToList(list, obs)
    ! dummy
    type(ListType), intent(inout) :: list
    type(ObserveType), pointer, intent(inout) :: obs
    ! local
    class(*), pointer :: obj => null()
    !
    obj => obs
    call list%Add(obj)
    !
    return
  end subroutine AddObserveToList

  function GetObserveFromList(list, indx) result(res)
    ! dummy
    type(ListType), pointer, intent(inout) :: list
    integer, intent(in) :: indx
    type(ObserveType), pointer :: res
    ! local
    class(*), pointer :: obj => null()
    !
    obj => list%GetItem(indx)
    select type (obj)
    type is (ObserveType)
      res => obj
    end select
    !
    return
  end function GetObserveFromList

end module ObserveModule
