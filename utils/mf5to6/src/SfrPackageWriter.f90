module SfrPackageWriterModule

  use ConstantsModule, only: DONE, LENMODELNAME, &
                             LENPACKAGENAME, LINELENGTH, MAXCHARLEN
  use ConstantsPHMFModule, only: FCINPUT, LENCTYPE, FCUNKNOWN
  use FileTypeModule, only: FileType
  use GLOBAL, only: NPER, IUNIT, NCOL, NROW, NLAY, IOUT, IBOUND
  use GlobalVariablesModule, only: masteridomain, echo
  use GWFSFRMODULE, only: CONST, DLEAK, ISTCB1, ISTCB2, NSTRM, ISEG, &
                          SEG, ISTRM, STRM, NeedLgrWaterMover, &
                          NeedLakWaterMover, NeedSfrWaterMover, ITMP
  use GwfSfrSubsModule, only: GWF2SFR7AR, GWF2SFR7RP
  use LineListModule, only: LineListType, same_lines
  use ListModule, only: ListType
  use ModelPackageModule, only: ModelPackageType, GetModelPack
  use MoverModule, only: MoverType, ConstructWaterMover, AddMoverToList, &
                         GetMoverFromList
  use MvrPackageWriterModule, only: MvrPackageWriterType
  use PackageWriterModule, only: PackageWriterType
  use SfrDiversionModule, only: SfrDiversionType, AddDiversionToList, &
                                GetDiversionFromList
  use SfrReachModule, only: SfrReachType, AddReachToList, GetReachFromList
  use SfrSegmentModule, only: SfrSegmentType, AddSegmentToList, &
                              GetSegmentFromList
  use SimPHMFModule, only: ustop, store_error, store_note, store_warning, &
                           count_errors
  use UtilitiesModule, only: RemoveElement
  use utl7module, only: U1DREL, U2DREL, urword, URDCOM, ULSTRD

  implicit none

  private
  public :: AddSfrPackageWriter, AllSfrPkgWriters, GetSfrPackageWriter, &
            GetSfrPackageWriterByIgrid, SfrPackageWriterType

  double precision, parameter :: DEFAULTUNITCONV = 1.0d0
  double precision, parameter :: DEFAULTDMAXDEPTHCHANGE = 1.0d-5
  type(ListType) :: AllSfrPkgWriters

  type, extends(PackageWriterType) :: SfrPackageWriterType
    integer :: MaxIteration = 100  ! Default in mf6
    integer :: SeqNum = 0
    double precision :: DmaxDepthChange = DEFAULTDMAXDEPTHCHANGE
    ! UnitConv = 1.0d0 is for meters & seconds; revise according to model
    ! units and SFR input instructions in user guide.
    double precision :: UnitConv = DEFAULTUNITCONV
    character(len=MAXCHARLEN) :: SfrOutFile = ''
    type(ListType), pointer :: Reaches => null()
    type(ListType), pointer :: Segments => null()
    type(ListType), pointer :: Diversions => null()
    type(MvrPackageWriterType), pointer :: MvrWriter => null()
  contains
    ! public procedures
    procedure, public :: BuildModelMovers
    procedure, public :: BuildSimMovers
    procedure, public :: GetSegment
    procedure, public :: MyType
    procedure, public :: ProcessAllocate
    procedure, public :: ProcessStressLoop
    procedure, public :: ReviseReachNetwork
    procedure, public :: WriteOptions
    procedure, public :: WriteDimensions
    procedure, public :: WriteStressPeriodListData
    ! private procedures
    procedure, private :: AddReach
    procedure, private :: AddDiversion
    procedure, private :: AssignReachData
    procedure, private :: ConnectReachesToSegments
    procedure, private :: DefineAllConnections
    procedure, private :: DefineReachDiversions
    procedure, private :: GetReach
    procedure, private :: GetDiversion
    procedure, private :: process_options
    procedure, private :: write_reaches
    procedure, private :: write_connectivity
    procedure, private :: write_diversions
  end type SfrPackageWriterType

contains

  subroutine AddSfrPackageWriter(sfrPkgWriter)
    implicit none
    ! dummy
    type(SfrPackageWriterType), pointer, intent(inout) :: sfrPkgWriter
    ! local
    class(*), pointer :: obj => null()
    !
    obj => sfrPkgWriter
    call AllSfrPkgWriters%Add(obj)
    !
    return
  end subroutine AddSfrPackageWriter

  function GetSfrPackageWriter(idx) result (res)
    implicit none
    ! dummy
    integer, intent(in) :: idx
    type(SfrPackageWriterType), pointer :: res
    ! local
    class(*), pointer :: obj => null()
    !
    res => null()
    obj => AllSfrPkgWriters%GetItem(idx)
    select type (obj)
    type is (SfrPackageWriterType)
      res => obj
    end select
    !
    return
  end function GetSfrPackageWriter

  function GetSfrPackageWriterByIgrid(igrid) result (res)
    implicit none
    ! dummy
    integer, intent(in) :: igrid
    type(SfrPackageWriterType), pointer :: res
    ! local
    integer :: i, n
    type(SfrPackageWriterType), pointer :: spw => null()
    !
    res => null()
    n = AllSfrPkgWriters%Count()
    do i=1,n
      spw => GetSfrPackageWriter(i)
      if (spw%Igrid == igrid) then
        res => spw
        exit
      endif
    enddo
    !
    return
  end function GetSfrPackageWriterByIgrid

  subroutine ProcessAllocate(this, igrid)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    character(len=MAXCHARLEN) :: fname
    character(len=100) :: ctemp
    ! format
    10 format(A,'_',i0)
    !
    this%Igrid = igrid
    call this%AllocatePointers()
    if (associated(this%Reaches)) then
      call this%Reaches%Clear(.true.)
      deallocate(this%Reaches)
    endif
    allocate(this%Reaches)
    this%Reaches%name = 'Reaches'
    if (associated(this%Segments)) then
      call this%Segments%Clear(.true.)
      deallocate(this%Segments)
    endif
    allocate(this%Segments)
    this%Segments%name = 'Segments'
    if (associated(this%Diversions)) then
      call this%Diversions%Clear(.true.)
      deallocate(this%Diversions)
    endif
    allocate(this%Diversions)
    this%Diversions%name = 'Diversions'
    this%Active = .true.
    this%DefaultBudgetText = 'SFR'
    ! Ned todo: not sure that this will work for LGR--need to look at LGR source
    if (iunit(44) > 0) then
      this%IuOrig = IUNIT(44)    ! Specific to SFR package
      this%OrigPkg = 'SFR'
    else
      call store_warning('STR package is not yet supported')
      this%Active = .false.
      return
    endif
    write(this%PackageName,10)trim(this%OrigPkg), this%Igrid
    this%fileobj%PkgName = this%PackageName
    this%Segments%name = this%PackageName
    this%fileobj%FCode = FCINPUT
    !
    this%fileobj%FType = 'SFR6'
    this%PkgType = 'SFR'
    fname = trim(this%ModelBasename) // '.sfr'
    call this%FileWriterType%InitializeFile(fname, this%fileobj%FType, &
                                            this%PackageName)
    this%IuOrig = iunit(44)
    this%IunitBcf = iunit(1)
    this%IunitGwt = iunit(15)
    this%IunitHuf = iunit(37)
    this%IunitLpf = iunit(23)
    this%IunitUzf = iunit(55)
    call GWF2SFR7AR(this%IuOrig, this%IunitBcf, this%IunitLpf, this%IunitHuf, &
                    this%IunitGwt, this%Nsol, this%Iouts, this%IunitUzf, &
                    igrid, this%Reaches, this%Segments, this%Diversions, &
                    this%Active)
    if (NeedLgrWaterMover) this%NeedWaterMover = .true.
    if (.not. this%Active) then
      this%fileobj%FCode = FCUNKNOWN
      return
    endif
    !
    ! NStressDim is number of values after k,i,j required to specify boundary stress
    !this%NStressDim = 3     ! Specific to RIV package; for Stage, Cond, Rbot
    !
    ! NVALP1 is rlist index for element following minimum required data
    !this%nvalp1 = this%NStressDim + 4
    !
    ! Assign package-specific pointers
    this%NBndPeriod => NSTRM
    !
    this%UnitConv = CONST
    this%DmaxDepthChange = DLEAK
    this%Istcb2 = ISTCB2
    this%ICbc => ISTCB1
    !
    ! define format for printing boundary data
    write(ctemp,'(i0)')this%NAux + this%NStressDim
    this%fmat = '(3(2x,i0),' // trim(ctemp) // '(2x,g15.8))'
    !
    return
  end subroutine ProcessAllocate

  subroutine WriteOptions(this)
    ! Overrides PackageWriterType%WriteOptions
    implicit none
    class(SfrPackageWriterType) :: this
    ! local
    integer :: i, iu
    ! formats
    5  format()
    10 format(2x,a,2x,G15.8)
    20 format(2x,a,2x,a)
    50 format(2x,a)
    60 format(a)
    !
    call this%process_options()
    if (.not. this%Active) return
    !
    iu = this%fileobj%IUnit
    ! Write BEGIN Options
    write(iu,5)
    write(iu,60)'BEGIN Options'
    !
    ! Aux variable names
    if (this%NAux > 0) then
      do i=1,this%NAux
        write(iu,20) 'AUXILIARY', trim(this%Aux(i))
      enddo
    endif
    !
!    if (echo) then
!      if (this%ipr==1) then
        write(iu,50)'PRINT_INPUT'
!      endif
!    endif
    !
    if (associated(this%ICbc)) then
      if (this%ICbc > 0) then
        write(iu,50)'SAVE_FLOWS'
      elseif (this%ICbc < 0 .or. this%SfrOutFile /= '') then
        write(iu,50)'PRINT_STAGE'
        write(iu,50)'PRINT_FLOWS'
      endif
    endif
    !
    ! Ned TODO: Write MAXIMUM_DEPTH_CHANGE, UNIT_CONVERSION
    ! SFR_OUTPUT as needed here.
    ! (Allow MAXIMUM_ITERATION to default.)
    !
    if (this%DmaxDepthChange /= DEFAULTDMAXDEPTHCHANGE) then
      write(iu,10)'MAXIMUM_DEPTH_CHANGE',this%DmaxDepthChange
    endif
    !
    if (this%UnitConv /= DEFAULTUNITCONV) then
      write(iu,10)'UNIT_CONVERSION',this%UnitConv
    endif
    !
    if (this%NeedWaterMover .or. NeedSfrWaterMover) then
      write(iu,50)'MOVER'
    endif
    !
    ! Write END Options
    write(iu,60)'END Options'
    !
    return
  end subroutine WriteOptions

  subroutine WriteDimensions(this)
    ! Overrides PackageWriterType%WriteDimensions
    implicit none
    class(SfrPackageWriterType) :: this
    ! local
    integer :: iu
    ! formats
    5  format()
    30 format(2x,a,2x,i0)
    60 format(a)
    !
    if (.not. this%NeedDimensionsBlock) return
    !
    iu = this%fileobj%IUnit
    ! Write BEGIN Dimensions
    write(iu,5)
    write(iu,60)'BEGIN Dimensions'
    !
    write(iu,30)'NREACHES',this%MaxActiveBnd
    !
    ! Write END Dimensions
    write(iu,60)'END Dimensions'
    !
    return
  end subroutine WriteDimensions

  subroutine ProcessStressLoop(this, igrid)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: kper
    logical :: currentA
    ! formats
    10 format(' Processing ',a,' package input...')
    20 format(/,'Processing ',a,' package for stress period ',i0)
    !
    ! Write Options and Dimensions blocks
    if (.not. this%Active) return
!    call this%WriteOptions()
    !
    ! Initially, current block is BlockA; alternate each stress period
    this%CurrentBlock => this%BlockA
    this%PreviousBlock => this%BlockB
    currentA = .true.
    !
    do kper=1,nper
      if (kper==1) write(*,10)trim(this%OrigPkg)
      write(iout,20)trim(this%OrigPkg), kper
      ! Read MF2005 input for current stress period
      !                                                               kstp
      call GWF2SFR7RP(this%IuOrig, this%IunitGwt, this%IunitLak, kper, -9, &
                      this%Nsol, this%Iouts, this%IunitBcf, this%IunitLpf, &
                      this%IunitHuf, this%IunitUzf, igrid, this%Reaches, &
                      this%Segments, this%Diversions, this%Active, &
                      this%SeqNum)
      if (NeedLakWaterMover) this%NeedWaterMover = .true.
      if (.not. this%Active) then
        this%fileobj%FCode = FCUNKNOWN
        return
      endif
      !
      if (itmp >= 0) then
        ! Set up connections and diversions
        call this%DefineAllConnections()
      endif
      !
      if (kper==1) then
        !
        if (associated(masteridomain)) then
          call this%ReviseReachNetwork()
        endif
        !
        call this%WriteOptions()
        call this%WriteDimensions()
        call this%write_reaches()
        call this%write_connectivity()
        call this%write_diversions()
      endif
      !
      ! Populate reaches with data that need
      ! to be written to PERIOD block.
      call this%AssignReachData()
      !
      ! Write MF6 input for current stress period to a LineList
      call this%CurrentBlock%Clear(.true.)
      call this%WriteStressPeriodListData(this%CurrentBlock)
      ! Write block to MF6 input file if needed
      call this%WriteBlockIfNeeded(kper)
      !
      ! switcheroo
      if (currentA) then
        this%CurrentBlock => this%BlockB
        this%PreviousBlock => this%BlockA
        currentA = .false.
      else
        this%CurrentBlock => this%BlockA
        this%PreviousBlock => this%BlockB
        currentA = .true.
      endif
    enddo
    !
    call this%BlockA%Clear(.true.)
    call this%BlockB%Clear(.true.)
    !
    return
  end subroutine ProcessStressLoop

  subroutine AssignReachData(this)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    ! local
    integer :: icalc, iprior, iupseg
    integer :: i, irch, nrch, nseg
    integer :: rchnumpkg, rchnumseg
    real :: area
    real :: q
    type(SfrSegmentType), pointer :: segptr
    type(SfrReachType), pointer :: rch
    !
    ! Iterate through segments
    nseg = this%Segments%Count()
    do i=1,nseg
      segptr => this%GetSegment(i)
      icalc = iseg(1,i)
      iprior = segptr%Iprior
      iupseg = segptr%Iupseg
      ! Iterate through reaches in segment
      nrch = segptr%SegReaches%Count()
      do irch=1,nrch
        rch => segptr%GetReach(irch)
        rch%icalc = icalc
        rchnumpkg = rch%rnopkg
        rchnumseg = rch%rnoseg
        area = rch%rwid * rch%rlen
        ! Get interpolated values from SFR2 arrays.
        if (rchnumseg == 1) then
          ! Get segment-based values
          if (iupseg > 0) then
            ! this reach diverts from segment iupseg
            if (iprior == -2) then
              rch%ustrf = SEG(2,i)   ! fraction of available flow in upstream reach to be diverted
              rch%inflow = 0.0
            else
              rch%ustrf = 0.0
              rch%inflow = 0.0
              !rch%inflow = SEG(2,i)  ! flow diverted to 1st reach of segment
              ! Ned todo: Assign diversion +/or upstream fraction
              !rch%diversion = SEG(2,i)
            endif
          else
            q = SEG(2,i)  ! inflow to 1st reach of segment
            if (q /= rch%inflow) then
              rch%inflow = q
              rch%inflow_iprn = 1
            else
              rch%inflow_iprn = 0
            end if
          endif
        else
          rch%ustrf = 0.0
          rch%inflow = 0.0
        endif
        q = STRM(14,rchnumpkg) / area   ! convert rainfall Q to flux
        if (q /= rch%rainfall) then
          rch%rainfall = q
          rch%rainfall_iprn = 1
        else
          rch%rainfall_iprn = 0
        end if
        q = STRM(13,rchnumpkg) / area       ! convert evaporation Q to flux
        if (q /= rch%evap) then
          rch%evap = q
          rch%evap_iprn = 1
        else
          rch%evap_iprn = 0
        end if
        q = STRM(12,rchnumpkg)
        if (q /= rch%runoff) then
          rch%runoff = q
          rch%runoff_iprn = 1
        else
          rch%runoff_iprn = 0
        end if
      enddo
    enddo
    !
    return
  end subroutine AssignReachData

  subroutine WriteStressPeriodListData(this, lineList)
    ! Overrides PackageWriterType%WriteStressPeriodListData
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    type(LineListType), pointer :: lineList
    ! local
    character(len=MAXCHARLEN) :: line
    integer :: icalc, idiv, irch, ndiv, nrch
    type(SfrDiversionType), pointer :: div
    type(SfrReachType), pointer :: rch
    ! formats
    10 format(i0,1x,a,1x,i0,1x,g22.15)
    20 format(i0,1x,a,1x,a)
    40 format(i0,1x,a,1x,g22.15)
    !
    if (.not. this%Active) return
    !
    ! Ned todo: Write stress period data for MF6.
    ! Write keywords and values for each reach.
    ! Keywords include: inflow, rainfall, evaporation, runoff,
    !                   diversion, upstream_fraction, auxiliary.
    !
    ! Iterate through  reaches
    nrch = this%Reaches%Count()
    reachloop: do irch=1,nrch
      rch => this%GetReach(irch)
      if (.not. rch%KeepReachActive) cycle reachloop
      icalc = rch%icalc
      ! inflow
      if (rch%inflow_iprn /= 0) then
        write(line,40)rch%newReachNum,'inflow',rch%inflow
        call lineList%AddLine(line)
      end if
      ! rainfall
      if (rch%rainfall_iprn /= 0) then
        write(line,40)rch%newReachNum,'rainfall',rch%rainfall
        call lineList%AddLine(line)
      end if
      ! evaporation
      if (rch%evap_iprn /= 0) then
        write(line,40)rch%newReachNum,'evaporation',rch%evap
        call lineList%AddLine(line)
      end if
      ! runoff
      if (rch%runoff_iprn /= 0) then
        write(line,40)rch%newReachNum,'runoff',rch%runoff
        call lineList%AddLine(line)
      end if
      ! iterate through diversions (volumetric or as fraction of upstream flow)
      ndiv = rch%Diversions%Count()
      divloop: do idiv=1,ndiv
        div => rch%GetDiversion(idiv)
        if (.not. div%KeepDiversionActive) cycle divloop
        if (div%Iprior == -2) then
          write(line,10)rch%newReachNum,'diversion',idiv,div%Ustrf
        else
          write(line,10)rch%newReachNum,'diversion',idiv,div%DivFlow
        endif
        call lineList%AddLine(line)
      enddo divloop
      !
      ! Simple status reaches
      if (icalc <= 0) then
        write(line,20)rch%newReachNum, 'status', 'simple'
        call lineList%AddLine(line)
        write(line,40)rch%newReachNum, 'stage', rch%stage
        call lineList%AddLine(line)
      endif
    enddo reachloop
    !
    return
  end subroutine WriteStressPeriodListData

  subroutine process_options(this)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    ! local
    type(FileType), pointer :: filtyp => null()
    character(len=LINELENGTH) :: msg
    !
    ! If MaxActiveBnd has been assigned by ReviseReachNetwork, do not change it.
    if (this%MaxActiveBnd < 0) then
      this%MaxActiveBnd = NSTRM
    endif
    !
    ! Process ISTCB2
    if (this%Istcb2 > 0) then
      ! Calls for formatted output of info for reaches.
      ! Find file associated with ISTCB2
      filtyp => this%Mf6Files%GetFileByUnit(this%Istcb2, .false.)
      if (associated(filtyp)) then
        this%SfrOutFile = filtyp%FName
      else
        this%SfrOutFile = trim(this%ModelBasename) // '_sfr_flows.out'
      endif
    elseif (this%Istcb2 < 0) then
      ! Calls for unformatted output of streamflow for reaches.
      ! SFR6 does not support this option, I think.
      msg = 'SFR6 does not support unformatted output of streamflows as ' // &
            'specified in SFR2 by ISTCB2 < 0.'
      call store_note(msg)
    endif
    !
    return
  end subroutine process_options

  subroutine AddReach(this, reach)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    type(SfrReachType), pointer :: reach
    !
    call AddReachToList(this%Reaches, reach)
    !
    return
  end subroutine AddReach

  function GetReach(this, idx) result (res)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    integer, intent(in) :: idx
    type(SfrReachType), pointer :: res
    !
    res => GetReachFromList(this%Reaches, idx)
    !
    return
  end function GetReach

  function GetSegment(this, idx) result (res)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    integer, intent(in) :: idx
    type(SfrSegmentType), pointer :: res
    !
    res => GetSegmentFromList(this%Segments, idx)
    !
    return
  end function GetSegment

  subroutine AddDiversion(this, diversion)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    type(SfrDiversionType), pointer :: diversion
    !
    call AddDiversionToList(this%Diversions, diversion)
    !
    return
  end subroutine AddDiversion

  function GetDiversion(this, idx) result (res)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    integer, intent(in) :: idx
    type(SfrDiversionType), pointer :: res
    !
    res => GetDiversionFromList(this%Diversions, idx)
    !
    return
  end function GetDiversion

  subroutine ConnectReachesToSegments(this)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    ! local
    integer :: irch, nrch, segnum
    type(SfrReachType), pointer :: rch => null()
    type(SfrSegmentType), pointer :: seg => null()
    !
    nrch = this%Reaches%Count()
    do irch=1,nrch
      ! Connect reach to segment
      rch => this%GetReach(irch)
      segnum = rch%segnum
      seg => this%GetSegment(segnum)
      call seg%AddReach(rch)
    enddo
    !
    return
  end subroutine ConnectReachesToSegments

  subroutine DefineAllConnections(this)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    ! local
    integer :: i, idiv, iseg, n, ndiv, nrch, nseg
    integer :: rchnumDs, rchnumUs, segnumDs, segnumUs
    integer :: outseg
    type(SfrReachType), pointer :: rchDn => null()
    type(SfrReachType), pointer :: rchUp => null()
    type(SfrReachType), pointer :: rchLast => null()
    type(SfrReachType), pointer :: rchFirst => null()
    type(SfrSegmentType), pointer :: seg => null()
    type(SfrSegmentType), pointer :: segDown => null()
    type(SfrSegmentType), pointer :: segUp => null()
    type(SfrDiversionType), pointer :: div => null()
    !
    call this%ConnectReachesToSegments()
    call this%DefineReachDiversions()
    !
    nseg = this%Segments%Count()
    ! Iterate through segments
    do iseg=1,nseg
      seg => this%GetSegment(iseg)
      ! Connect reaches within this segment to each other
      n = seg%SegReaches%Count()
      do i=1,n
        rchDn => seg%GetReach(i)
        if (i>1) then
          ! Make reach-reach connection
          call rchUp%AddConnection(-rchDn%rnopkg)
          call rchDn%AddConnection(rchUp%rnopkg)
        endif
        rchUp => rchDn
      enddo
      !
      ! Connnect last reach in this segment to first reach in segment outseg
      outseg = seg%Outseg
      if (outseg > 0) then
        ! Get last reach in current segment
        nrch = seg%SegReaches%Count()
        rchLast => seg%GetReach(nrch)
        ! Get first reach in downstream segment
        segDown => this%GetSegment(outseg)
        rchFirst => segDown%GetReach(1)
        ! Make reach-reach connection
        call rchLast%AddConnection(-rchFirst%rnopkg)
        call rchFirst%AddConnection(rchLast%rnopkg)
      endif
    enddo
    !
    ndiv = this%Diversions%Count()
    ! Iterate through diversions
    do idiv=1, ndiv
      div => this%GetDiversion(idiv)
      ! Make connections between reach from which flow
      !     is diverted and reach receiving flow.
      !
      ! Find upstream reach number
      segnumUs = div%SegnumUs
      segUp => this%GetSegment(segnumUs)
      nrch = segUp%SegReaches%Count()
      rchUp => segUp%GetReach(nrch)
      rchnumUs = rchUp%rnopkg
      ! Find downstream reach number
      segnumDs = div%SegnumDs
      segDown => this%GetSegment(segnumDs)
      rchDn => segDown%GetReach(1)
      rchnumDs = rchDn%rnopkg
      ! Make connections
      call rchUp%AddConnection(-rchnumDs)
      call rchDn%AddConnection(rchnumUs)
    enddo
    !
    return
  end subroutine DefineAllConnections

  subroutine DefineReachDiversions(this)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    ! local
    !
    integer :: idiv, iprior, ndiv, nrch, &
               rchnumDs, rchnumUs, segnumDs, segnumUs
    type(SfrReachType), pointer :: rchDs => null()
    type(SfrReachType), pointer :: rchUs => null()
    type(SfrDiversionType), pointer :: div => null()
    type(SfrSegmentType), pointer :: segDs => null()
    type(SfrSegmentType), pointer :: segUs => null()
    !
    ! Iterate through all diversions defined for package
    ndiv = this%Diversions%Count()
    do idiv=1,ndiv
      div => this%GetDiversion(idiv)
      ! Get segment from which flow is diverted
      segnumUs = div%SegnumUs
      segUs => this%GetSegment(segnumUs)
      ! Get segment receiving diverted flow
      segnumDs = div%SegnumDs
      segDs => this%GetSegment(segnumDs)
      ! Get reach from which flow is diverted
      nrch = segUs%SegReaches%Count()
      rchUs => segUs%GetReach(nrch)
      rchnumUs = rchUs%rnopkg
      ! Get reach receiving diverted flow
      rchDs => segDs%GetReach(1)
      rchnumDs = rchDs%rnopkg
      ! The reach from which flow is diverted is identified
      ! by "rno" in the reach_diversions block.
      div%Rno = rchnumUs
      ! The reach receiving diverted flow is identified
      ! by "iconr" in the reach_diversions block.
      div%Iconr = rchnumDs
      ! Assign the priority and, if iprior=-2, the flow fraction
      iprior = div%Iprior
      if (iprior == -2) then
        ! Value entered as FLOW is fraction of available flow
        ! in upstream segment.
        segUs%Ustrf = segUs%Flow  ! should be a fraction between 0.0 and 1.0
      else
        ! By default, set ustrf for reach receiving diverted flow to 0.0.
        ! Activation of diversion will be handled in PERIOD blocks.
        rchDs%ustrf = 0.0
        ! Assign diversion flow
        div%DivFlow = segDs%Flow
      endif
      ! Add diversion to list of diversions belonging to upstream reach
      call rchUs%AddDiversion(div)
    enddo
    !
    return
  end subroutine DefineReachDiversions

  subroutine write_reaches(this)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    ! local
    integer :: i, iu, ncon, ndiv, nrch
    integer :: irow, jcol, klay
    type(SfrReachType), pointer :: r
    ! formats
    5  format()
    10 format(2x,a,2x,G14.7)
    20 format(2x,a,2x,a)
    30 format(1x,4(1x,i0),7(2x,g14.7),2x,i0,2x,f6.4,2x,i0)
    35 format(2x,i0,2x,'NONE',7(2x,g14.7),2x,i0,2x,f6.4,2x,i0)
    50 format(2x,a)
    60 format(a)
    !
    iu = this%fileobj%IUnit
    nrch = this%Reaches%Count()
    !
    ! Write BEGIN REACHES
    write(iu,5)
    write(iu,60)'BEGIN PACKAGEDATA'
    write(iu,50)'# rno cellid rlen           rwid            ' //  &
                'rgrd            rbtp            rbth            ' // &
                'rbhk             man       ncon  ustrf  ndv'
    !
    do i=1,nrch
      r => this%GetReach(i)
      if (.not. r%KeepReachActive) cycle
      ncon = size(r%iconn)
      ndiv = r%Diversions%Count()
      klay = r%klay
      irow = r%irow
      jcol = r%jcol
      if (IBOUND(jcol,irow,klay) /= 0) then
        write(iu,30)r%newReachNum, klay, irow, jcol, r%rlen, r%rwid, r%rgrd, &
                    r%rtp, r%rbth, r%rbhk, r%man, ncon, r%ustrf, ndiv
      else
        write(iu,35)r%newReachNum, r%rlen, r%rwid, r%rgrd, &
                    r%rtp, r%rbth, r%rbhk, r%man, ncon, r%ustrf, ndiv
      endif
    enddo
    !
    ! Write END
    write(iu,60)'END PACKAGEDATA'
    !
    return
  end subroutine write_reaches

  subroutine write_connectivity(this)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    ! local
    integer :: i, iu, j, ncon, nrch
    type(SfrReachType), pointer :: reach => null()
    character(len=MAXCHARLEN) :: line
    character(len=10) :: ic
    ! formats
    5  format()
    10 format(2x,a,2x,G15.8)
    20 format(2x,a,2x,a)
    50 format(2x,a)
    60 format(a)
    70 format(2x,i0)
    80 format(i0)
    !
    iu = this%fileobj%IUnit
    ! Write BEGIN REACH_CONNECTIVITY
    write(iu,5)
    write(iu,60)'BEGIN CONNECTIONDATA'
    nrch = this%Reaches%Count()
    do i=1,nrch
      reach => this%GetReach(i)
      if (.not. reach%KeepReachActive) cycle
      ncon = size(reach%iconn)
      if (ncon > 0) then
        write(line,70)reach%NewReachNum
        do j=1,ncon
          write(ic,80)reach%iconn(j)
          line = trim(line) // ' ' // ic
        enddo
        write(iu,60)trim(line)
      endif
    enddo
    !
    ! Write END
    write(iu,60)'END CONNECTIONDATA'
    !
    return
  end subroutine write_connectivity

  subroutine write_diversions(this)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    ! local
    character (len=10) :: cprior
    integer :: i, iu, j, ndiv, nrch
    logical :: need
    type(SfrDiversionType), pointer :: div => null()
    type(SfrReachType), pointer :: reach => null()
    ! formats
    5  format()
    10 format(3(2x,i0),2x,a10)
    60 format(a)
    !
    need = .false.
    nrch = this%Reaches%Count()
    do i=1,nrch
      reach => this%GetReach(i)
      if (.not. reach%KeepReachActive) cycle
      ndiv = reach%Diversions%Count()
      if (ndiv > 0) then
        need = .true.
        exit
      endif
    enddo
    !
    if (need) then
      iu = this%fileobj%IUnit
      ! Write BEGIN REACH_DIVERSIONS
      write(iu,5)
      write(iu,60)'BEGIN DIVERSIONS'
      !
      write(iu,60)'  # rno idv iconr iprior'
      do i=1,nrch
        reach => this%GetReach(i)
        if (.not. reach%KeepReachActive) cycle
        ndiv = reach%Diversions%Count()
        if (ndiv > 0) then
          do j=1,ndiv
            div => reach%GetDiversion(j)
            if (associated(div)) then
              select case (div%Iprior)
                case (0)
                  cprior = 'UPTO'
                case (-1)
                  cprior = 'THRESHOLD'
                case (-2)
                  cprior = 'FRACTION'
                case (-3)
                  cprior = 'EXCESS'
              end select
              write(iu, 10) reach%rnopkg, j, div%Iconr, cprior
            endif
          enddo
        endif
      enddo
      !
      ! Write END
      write(iu,60)'END DIVERSIONS'
    endif
    !
    return
  end subroutine write_diversions

  subroutine BuildModelMovers(this)
    ! dummy
    class(SfrPackageWriterType) :: this
    ! local
    integer :: i, idProv, idRec, nsegs, iupseg, ioutseg
    type(MoverType), pointer :: mover => null()
    type(SfrSegmentType), pointer :: segmnt => null()
    type(SfrReachType), pointer   :: reach => null()
    type(ModelPackageType), pointer :: modelPack
    character(len=LENMODELNAME) :: provModelName, recModelName
    character(len=LENPACKAGENAME) :: provPkgName, recPkgName
    ! format
    10 format('No SFR package found with IGRID = ',i0)
    12 format(a,'_',i0)
    !
    if (.not. this%NeedWaterMover) return
    !
    nsegs = this%Segments%Count()
    segmentloop: do i=1,nsegs
      segmnt => GetSegmentFromList(this%Segments, i)
      iupseg = segmnt%Iupseg
      ioutseg = segmnt%Outseg
      ! First, consider inflow to this segment.
      ! If LgrGrid > 0, water is provided from SFR of another model.
      if (segmnt%LgrGrid > 0) then
        ! Inter-model flow to this segment will be handled
        ! by BuildSimMovers, not here. Do nothing.
      else
        ! Iupseg is the ID number of the provider segment (or lake).
        if (iupseg < 0) then
          ! Inflow for this segment comes from a lake in the same model
          ! as this segment.  Provider is lake number abs(iupseg), and
          ! receiver is first reach in this segment.
          ! [Really, the provider is an outlet ID, but ID number should
          ! be the same as for the lake?]
          idProv = abs(iupseg)
          reach => segmnt%GetReach(1)
          idRec = reach%rnopkg
          ! Assign package name and model name of provider lake package
          write(provPkgName,12)'LAK',this%Igrid
          modelPack => GetModelPack(provPkgName)
          provModelName = modelPack%ModelName
          recPkgName = this%PackageName
          recModelName = provModelName
          ! Construct a LAK -> SFR water mover
          call ConstructWaterMover(mover, 'FACTOR', provModelName, &
                   recModelName, provPkgName, recPkgName, idProv, idRec, &
                   this%Igrid, this%Igrid, 'LAK', 'SFR', DONE)
          call this%AddMover(mover)
        endif
      endif
      ! Next, consider outflow from this segment.
      ! Ioutseg is the ID number of the receiver segment (or lake).
      if (ioutseg < 0) then
        ! Outflow from this segment goes to a lake in the same model as
        ! this segment.  Receiver is lake number abs(ioutseg), and
        ! provider is first reach in this segment.
        idRec = abs(ioutseg)         ! ID number of lake receiving flow
        reach => segmnt%GetReach(1)  ! First reach in this segment
        idProv = reach%rnopkg        ! Reach number in package
        ! Assign package name and model name of receiver lake package
        write(recPkgName,12)'LAK',this%Igrid
        modelPack => GetModelPack(recPkgName)
        recModelName = modelPack%ModelName
        provPkgName = this%PackageName
        provModelName = recModelName
        ! Construct a SFR -> LAK water mover
        call ConstructWaterMover(mover, 'FACTOR', provModelName, &
                 recModelName, provPkgName, recPkgName, idProv, idRec, &
                 this%Igrid, this%Igrid, 'SFR', 'LAK', DONE)
        call this%AddMover(mover)
      endif
    enddo segmentloop
    !
    if (count_errors() > 0) call ustop()
    return
  end subroutine BuildModelMovers

  subroutine BuildSimMovers(this)
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    ! local
    integer :: i, idProv, idRec, nreaches, nsegs
    type(MoverType), pointer :: mover => null()
    type(SfrSegmentType), pointer :: segmnt => null()
    type(SfrSegmentType), pointer :: provSegmnt => null()
    type(SfrReachType), pointer   :: reach => null()
    type(SfrPackageWriterType), pointer :: provSfrPkgWriter => null()
    type(ModelPackageType), pointer :: modelPack
    character(len=LINELENGTH) :: ermsg
    character(len=LENMODELNAME) :: provModelName, recModelName
    character(len=LENPACKAGENAME) :: provPkgName, recPkgName
    ! format
    10 format('No SFR package found with IGRID = ',i0)
    !
    if (.not. this%NeedWaterMover) return
    !
    nsegs = this%Segments%Count()
    do i=1,nsegs
      segmnt => GetSegmentFromList(this%Segments, i)
      ! if LgrGrid > 0, water is provided from SFR of
      ! grid LgrGrid to this segment, which is a receiver.
      if (segmnt%LgrGrid > 0) then
        ! Ned todo: Need to find last reach of provider segment and put
        ! its number (in list of all reaches in SFR6 package) as
        ! idProvider argument.
        ! segmnt%LgrSeg is the segment number of the provider segment.
        provSfrPkgWriter => GetSfrPackageWriterByIgrid(segmnt%LgrGrid)
        if (associated(provSfrPkgWriter)) then
          provSegmnt => provSfrPkgWriter%GetSegment(segmnt%LgrSeg) ! provider segment
          nreaches = provSegmnt%SegReaches%Count() ! # reaches in provider segment
          reach => provSegmnt%GetReach(nreaches) ! Last reach in provider segment
          if (associated(reach)) then
            idProv = reach%newReachNum ! Reach number in package of last reach in provider segment
          else
            call store_error('programmer error in BuildSimMovers')
          endif
          ! Ned todo: Need to get first reach of receiver segment and put
          ! its number (in list of all reaches in SFR6 package) as
          ! idReceiver argument. Segmnt is the receiver segment.
          reach => segmnt%GetReach(1)
          if (associated(reach)) then
            idRec = reach%newReachNum ! Reach number in package of last reach in provider segment
          else
            call store_error('programmer error in BuildSimMovers')
          endif
          provPkgName = provSfrPkgWriter%PackageName
          modelPack => GetModelPack(provPkgName)
          provModelName = modelPack%ModelName
          recPkgName = this%PackageName
          recModelName = this%ModelPack%ModelName
          ! Construct a SFR -> SFR water mover
          call ConstructWaterMover(mover, 'FACTOR', provModelName, &
                                   recModelName, provPkgName, recPkgName, &
                                   idProv, idRec, segmnt%LgrGrid, &
                                   this%Igrid, 'SFR', 'SFR', DONE)
          call this%AddMover(mover)
        else
          write(ermsg,10)segmnt%LgrGrid
          call store_error(ermsg)
        endif
      endif
    enddo
    !
    if (count_errors() > 0) call ustop()
    return
  end subroutine BuildSimMovers

  subroutine ReviseReachNetwork(this)
    ! To be run after masteridomain has been defined, taking into
    ! account areas deactivated in child grid areas.
    implicit none
    ! dummy
    class(SfrPackageWriterType) :: this
    ! local
    integer :: i, icon, icr, idiv, ncon, ndiv, nreaches
    integer :: iconabs, krch
    integer :: klay, irow, jcol
    integer, allocatable, dimension(:) :: newReach
    type(SfrReachType), pointer :: reach => null()
    type(SfrReachType), pointer :: creach => null()
    type(SfrDiversionType), pointer :: diversion => null()
    !
    nreaches = this%Reaches%Count()
    allocate(newReach(nreaches))
    newReach = 0
    krch = 0
    do i=1,nreaches
      reach => GetReachFromList(this%Reaches, i)
      if (associated(reach)) then
        klay = reach%klay
        irow = reach%irow
        jcol = reach%jcol
        if (masteridomain(jcol,irow,klay) == 0) then
          ! Deactivate this reach
          reach%KeepReachActive = .false.
          reach%newReachNum = 0
          ndiv = reach%Diversions%Count()
          ! Deactivate all diversions from this reach
          do idiv=1,ndiv
            diversion => GetDiversionFromList(reach%Diversions, idiv)
            if (associated(diversion)) then
              diversion%KeepDiversionActive = .false.
            endif
          enddo
          ! Remove all connections to this reach from other reaches
          do icr=1,nreaches
            if (icr == i) cycle
            creach => GetReachFromList(this%Reaches, icr)
            if (associated(creach)) then
              if (.not. creach%KeepReachActive) cycle
              ncon = size(creach%iconn)
              connectionsloop: do icon=1,ncon
                ! If creach connects to this reach, remove the
                ! connection by removing the creach%iconn element
                ! that creates the connection.
                iconabs = abs(creach%iconn(icon))
                if (iconabs == i) then
                  ! Remove connection and exit connections loop
                  call RemoveElement(creach%iconn, icon)
                  exit connectionsloop
                endif
              enddo connectionsloop
            endif
          enddo
        else
          krch = krch + 1
          newReach(i) = krch
          reach%newReachNum = krch
        endif
      endif
    enddo
    !
    if (this%MaxActiveBnd /= krch) then
      ! Some reaches have been removed.
      ! Renumber the reaches, diversions, and connections.
      this%MaxActiveBnd = krch
      do i=1,nreaches
        reach => GetReachFromList(this%Reaches, i)
        if (associated(reach)) then
          ! Renumber the reaches in the active diversions
          ndiv = reach%Diversions%Count()
          do idiv=1,ndiv
            diversion => GetDiversionFromList(reach%Diversions, idiv)
            if (associated(diversion)) then
              if (diversion%KeepDiversionActive) then
                diversion%Iconr = newReach(diversion%Iconr)
              endif
            endif
          enddo
          ! Renumber connected reaches
          ncon = size(reach%iconn)
          connectionsloop2: do icon=1,ncon
            iconabs = abs(reach%iconn(icon))
            reach%iconn(icon) = sign(newReach(iconabs), reach%iconn(icon))
          enddo connectionsloop2
        endif
      enddo
    endif
    !
    return
  end subroutine ReviseReachNetwork

  function MyType(this) result (ctype)
    ! dummy
    class(SfrPackageWriterType) :: this
    character(len=LENCTYPE) :: ctype
    !
    ctype = 'SfrPackageWriterType'
    !
    return
  end function MyType

end module SfrPackageWriterModule

module junk
! ---SFR6 input-----------------   ---SFR2 input-----------------------
! MAXIMUM_ITERATION                -- (No equivalent in SFR2 input)
! MAXIMUM_DEPTH_CHANGE             DLEAK
! UNIT_CONVERSION                  CONST
! SFR_OUTPUT                       file associated with ISTCB2 (if > 0)
! NREACHES (# reaches)             NSTRM
!
! REACH VALUES:
! rno (reach number)               --
! cellid                           KRCH / IRCH / JRCH
! rlen (length)                    RCHLEN
! rwid (width)                     WIDTH1 and WIDTH2 (by segment)
! rgrd (gradient)                  SLOPE
! rtp (streambed top)              STRTOP
! rbth (streambed thickness)       STRTHICK
! rhk (bed conductivity)           STRHC1
! man (Manning's roughness)        ROUGHCH (by segment)
! ncon (# connected reaches)       -- (base on stream network)
! ustrf (upstream flow fraction)   -- (always 1?)
! ndv (# downstream diversions)    -- (base on stream network)
! [x,y,z] (aux var. values)        --

!                                                     SfrReachType member
!                                                     --------------------
! STRM(1, ) = RCHLEN                                  rlen
! STRM(2, ) = SLOPE                                   rgrd
! STRM(3, ) = STRTOP                                  rtp
! STRM(4, ) = streambed bottom elevation              --
! STRM(5, ) = interpolated width                      rwid
! STRM(6, ) = STRHC1                                  rbhk
! STRM(7, ) = interpolated depth                      --
! STRM(8, ) = STRTHICK                                rbth
! STRM(9, ) = --
! STRM(10, ) = --
! STRM(11, ) = --
! STRM(12, ) = volumetric runoff                      runoff (Q)
! STRM(13, ) = ET Q                                   [divide by area] -> evap (flux)
! STRM(14, ) = PPT Q                                  [divide by area] -> rainfall (flux)
! STRM(15, ) = water surface elevation (icalc=0)
! STRM(16, ) = streambed conductance
! STRM(17, ) =
! STRM(18, ) =
! STRM(19, ) =
! STRM(20, ) =
! STRM(21, ) = (UNSATURATED FLOW RATE FOR GAGE PACKAGE)
! STRM(22, ) = (UNSATURATED FLOW RATE FOR GAGE PACKAGE)
! STRM(23, ) = (UNSATURATED FLOW RATE FOR GAGE PACKAGE)
! STRM(24, ) =
! STRM(25, ) =
! STRM(26, ) =
! STRM(27, ) =
! STRM(28, ) =
! STRM(29, ) = Totdelstor (for UZF)
! STRM(30, ) =

! ISTRM(1, ) = krch (layer)
! ISTRM(2, ) = irch (row)
! ISTRM(3, ) = jrch (column)
! ISTRM(4, ) = jseg (segment #)
! ISTRM(5, ) = ireach (reach number in segment)
! ISTRM(6, ) = IFACE

! ISEG(1, ) = icalc
! ISEG(2, ) = # points in rating table
! ISEG(3, ) = Flags:
!             5: segment does not receive trib. flow,
!             6: segment diverts flow,
!             7: segment does receive trib. flow
! ISEG(4, ) = # of reaches in segment

!                                                     SfrReachType member
!                                                     --------------------
! SEG(1, ) = segment length
! SEG(2, ) = inflow rate                              inflow (to 1st reach in segment)
!             (or diversion fraction when             ustrf (for 1st reach in segment)
!              idivar(2, ) == -2)
! SEG(3, ) = runoff
! SEG(4, ) = etsw
! SEG(5, ) = pptsw
! SEG(6, ) = HC1FACT
! SEG(7, ) = THICKM1
! SEG(8, ) = ELEVUP
! SEG(9, ) = width1                                   rwid
! SEG(10, ) = depth1
! SEG(11, ) = HC2FACT
! SEG(12, ) = THICKM2
! SEG(13, ) = ELEVDN
! SEG(14, ) = width2
! SEG(15, ) = depth2
! SEG(16, ) = roughch (Manning's n)                   man
! SEG(17, ) = roughbk
! SEG(18, ) = (UZF)
! SEG(19, ) = (UZF)
! SEG(20, ) = (UZF)
! SEG(21, ) = (UZF)
! SEG(22, ) = (UZF)
! SEG(23, ) = (UZF)
! SEG(24, ) = (UZF)
! SEG(25, ) = (UZF)
! SEG(26, ) =

! IDIVAR(1, ) = IUPSEG
! IDIVAR(2, ) = IPRIOR

! XSEC(1-8, ) = XCPT
! XSEC(9-16, ) = ZCPT

! Reach roughness (constant over segment) is given by SEG(16,ISTRM(4,n)) where n is reach #
end module junk
