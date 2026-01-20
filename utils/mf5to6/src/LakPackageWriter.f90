module LakPackageWriterModule

  use ConstantsModule, only: DONE, DZERO, LINELENGTH, &
                             MAXCHARLEN, LENMODELNAME, LENPACKAGENAME
  use ConstantsPHMFModule, only: FCINPUT, FCUNKNOWN, LENCTYPE
  use FileTypeModule, only: FileType
  use GLOBAL, only: NPER, IUNIT, IOUT, ITMUNI, LENUNI
  use GlobalVariablesModule, only: masteridomain, echo
  use GWFLAKMODULE, only: ILKCB, NLAKES, NeedLakWaterMover, STAGES, SURFDEPTH, &
                          IRDTAB, PRCPLK, EVAPLK, RNF, WTHDRW, ITRB, IDIV, &
                          SGWF2LAK7PNT
  use GwfLakSubs, only: GWF2LAK7AR, GWF2LAK7RP
  use GWFSFRMODULE, only: IDIVAR, IOTSG, ISEG, ISTRM, NSS, SEG, STRM, &
                          SGWF2SFR7PNT, SfrPackageName
  use LakeModule, only: LakeType, AddLakeToList, GetLakeFromList, &
                        ConstructLakeType
  use LakeConnectionModule, only: LakeConnectionType, AddLakeConnectionToList, &
                                  GetConnectionFromList
  use LakeOutletModule, only: LakeOutletType, AddLakeOutletToList, &
                              GetOutletFromList, ConstructLakeOutlet
  use LakeTributaryModule, only: LakeTributaryType, AddLakeTributaryToList, &
                                 GetTributaryFromList, ConstructLakeTributary
  use LineListModule, only: LineListType, same_lines
  use ListModule, only: ListType
  use ModelPackageModule, only: GetModelPack, ModelPackageType
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
  public :: AllLakPkgWriters, LakPackageWriterType

  double precision, parameter :: DEFAULTUNITCONV = DONE
  type(ListType) :: AllLakPkgWriters

  type, extends(PackageWriterType) :: LakPackageWriterType
    integer :: MaxIteration = 100  ! Default in mf6
    integer :: SeqNum = 0
    integer :: Noutlets = 0
    integer :: Ntables = 0
    ! LengthConv = 1.0d0 is for meters
    ! TimeConv = 1.0d0 is for seconds
    ! Revise according to model units and
    ! LAK input instructions in user guide.
    double precision :: LengthConv = DEFAULTUNITCONV
    double precision :: TimeConv = DEFAULTUNITCONV
    double precision :: SurfDep = DZERO
    type(ListType), pointer :: Lakes => null()
    type(ListType), pointer :: Outlets => null()
    type(ListType), pointer :: Tribs => null()
    type(MvrPackageWriterType), pointer :: MvrWriter => null()
  contains
    ! public procedures
    procedure, public :: BuildModelMovers
!    procedure, public :: BuildLAKtoSFRmovers
!    procedure, public :: BuildSFRtoLAKmovers
    procedure, public :: MyType
    procedure, public :: ProcessAllocate
    procedure, public :: ProcessStressLoop
    procedure, public :: WriteDimensions
    procedure, public :: WriteOptions
    procedure, public :: WriteStressPeriodListData
    ! private procedures
    procedure, private :: AddLake
    procedure, private :: AddOutlet
    procedure, private :: AddTrib
    procedure, private :: AnyDivs
    procedure, private :: AnyTribs
    procedure, private :: DefineLakes
    procedure, private :: DefineAllConnections
    procedure, private :: DefineAllOutlets
    procedure, private :: DefineAllTributaries
    procedure, private :: FindNDiv
    procedure, private :: FindNTribs
    procedure, private :: FindReach1Index
    procedure, private :: GetLake
    procedure, private :: GetOutlet
    procedure, private :: GetTributary
    procedure, private :: process_options
    procedure, private :: write_connections
    procedure, private :: write_lakes
    procedure, private :: write_outlets
    procedure, private :: write_tables
  end type LakPackageWriterType

contains

  subroutine ProcessAllocate(this, igrid)
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    character(len=MAXCHARLEN) :: ermsg, fname
    character(len=100) :: ctemp
    ! format
    10 format(A,'_',i0)
    !
    this%Active = .true.
    this%Igrid = igrid
    call this%AllocatePointers()
    allocate(this%Lakes)
    allocate(this%Outlets)
    allocate(this%Tribs)
    this%DefaultBudgetText = 'LAKE'
    this%ipr = 1
    this%IuOrig = IUNIT(22)    ! Specific to LAK package
    this%OrigPkg = 'LAK'
    this%PkgType = 'LAK'
    write(this%PackageName,10)trim(this%PkgType),igrid
    this%fileobj%FCode = FCINPUT
    this%fileobj%FType = 'LAK6'
    this%fileobj%PkgName = this%PackageName
    fname = trim(this%ModelBasename) // '.lak'
    call this%FileWriterType%InitializeFile(fname, this%fileobj%FType, &
                                            this%PackageName)
    this%IunitBcf = iunit(1)
    this%IunitGwt = iunit(15)
    this%IunitLpf = iunit(23)
    this%IunitHuf = iunit(37)
    this%IunitSfr = iunit(44)
    this%IunitUzf = iunit(55)
    this%IunitUpw = iunit(62)
    call GWF2LAK7AR(this%IuOrig, this%IunitSfr, this%IunitGwt, this%IunitUzf, &
                    this%Nsol, igrid)
    if (NeedLakWaterMover) this%NeedWaterMover = .true.
    if (.not. this%Active) then
      this%fileobj%FCode = FCUNKNOWN
      return
    endif
    !
    this%NBndPeriod => NLAKES
    this%ICbc => ILKCB
    !
    select case (ITMUNI)
    case (0)
      ! undefined -- issue warning
      ermsg = 'ITMUNI is undefined, so mf5to6 is assuming time unit is seconds for LAK8 package input.'
      call store_warning(ermsg)
    case (1)
      ! seconds
      this%TimeConv = 1.0
    case (2)
      ! minutes
      this%TimeConv = 60.0
    case (3)
      ! hours
      this%TimeConv = 3600.0
    case (4)
      ! days
      this%TimeConv = 86400.0
    case (5)
      ! years
      this%TimeConv = 31557600.0
    end select
    !
    select case (LENUNI)
    case (0)
      ! undefined -- issue warning
      ermsg = 'LENUNI is undefined, so mf5to6 is assuming length unit is meters for LAK8 package input.'
      call store_warning(ermsg)
    case (1)
      ! feet
      this%LengthConv = 3.28081
    case (2)
      ! meters
      this%LengthConv = 1.0
    case (3)
      ! centimeters
      this%LengthConv = 100.0
    end select
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
    class(LakPackageWriterType) :: this
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
    if (echo) then
      if (this%ipr==1) then
        write(iu,50)'PRINT_INPUT'
      endif
    endif
    !
    if (associated(this%ICbc)) then
      if (this%ICbc > 0) then
        write(iu,50)'SAVE_FLOWS'
      elseif (this%ICbc < 0) then
        write(iu,50)'PRINT_FLOWS'
      endif
    endif
    !
    if (NeedLakWaterMover .or. this%NeedWaterMover) then
      write(iu,50)'MOVER'
    endif
    !
    if (this%SurfDep /= DZERO) then
      write(iu,10)'SURFDEP',this%SurfDep
    endif
    !
    if (this%TimeConv /= DEFAULTUNITCONV) then
      write(iu,10)'TIME_CONVERSION',this%TimeConv
    endif
    !
    if (this%LengthConv /= DEFAULTUNITCONV) then
      write(iu,10)'LENGTH_CONVERSION',this%LengthConv
    endif
    !!
    !if (.not. this%Newton) write(iu,50)'NO_NEWTON'
    !
    ! Write END Options
    write(iu,60)'END Options'
    !
    return
  end subroutine WriteOptions

  subroutine WriteDimensions(this)
    implicit none
    class(LakPackageWriterType) :: this
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
    write(iu,30)'NLAKES', this%MaxActiveBnd
    !
    this%Noutlets = this%Outlets%Count()
    if (this%Noutlets > 0) then
      write(iu,30)'NOUTLETS', this%Noutlets
    endif
    !
    if (this%Ntables > 0) then
      write(iu,30)'NTABLES', this%Ntables
    endif
    !
    ! Write END Dimensions
    write(iu,60)'END Dimensions'
    !
    return
  end subroutine WriteDimensions

  subroutine ProcessStressLoop(this, igrid)
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: kper
    logical :: currentA
!    integer :: kc, kp
    ! formats
    10 format(' Processing ',a,' package input...')
    20 format(/,'Processing ',a,' package for stress period ',i0)
    !
    ! Write Options and Dimensions blocks
    if (.not. this%Active) return
    !
    ! Check to see if SFR specifies potential for flow between any
    ! lakes and streams. If so, water mover(s) will be needed.
    call SGWF2LAK7PNT(this%Igrid)
    call SGWF2SFR7PNT(this%Igrid)
    if (this%AnyDivs() .or. this%AnyTribs()) then
      this%NeedWaterMover = .true.
    endif
    !
    ! Write the options block
    call this%WriteOptions()
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
      call GWF2LAK7RP(this%IuOrig, this%IunitBcf, this%IunitGwt, &
                      this%IunitLpf, this%IunitHuf, this%IunitSfr, &
                      this%IunitUzf, this%IunitUpw, kper, this%Nsol, &
                      this%Iouts, igrid)
      if (.not. this%Active) then
        this%fileobj%FCode = FCUNKNOWN
        return
      endif
      !
      if (kper==1) then
        ! Define lakes
        call this%DefineLakes()
        ! Define lake connections
        call this%DefineAllConnections()
        ! Define lake outlets
        call this%DefineAllOutlets()
        ! Define lake tributaries
        call this%DefineAllTributaries()
        !
!        if (associated(masteridomain)) then
!          call this%ReviseReachNetwork()
!        endif
        !
        call this%WriteDimensions()
        call this%write_lakes()
        call this%write_connections()
        call this%write_tables()
        call this%write_outlets()
      endif
      !
      ! Populate reaches with data that need
      ! to be written to PERIOD block.
!      call this%AssignReachData()
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

  subroutine WriteStressPeriodListData(this, lineList)
    ! Overrides PackageWriterType%WriteStressPeriodListData
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    type(LineListType), pointer :: lineList
    ! local
    character(len=MAXCHARLEN) :: line
    integer :: i, iu, nlak
    ! formats
    10 format(2x,a,2x,i0,2x,a)
    20 format(2x,i0,2x,a,2x,g15.8)
    !
    if (.not. this%Active) return
    !
    iu = this%fileobj%IUnit
    nlak = this%Lakes%Count()
    !
    do i=1,nlak
      !lake => this%GetLake(i)
      line = ' '
      write(line,20) i, 'rainfall', PRCPLK(i)
      call lineList%AddLine(line)
      line = ' '
      write(line,20) i, 'evaporation', EVAPLK(i)
      call lineList%AddLine(line)
      line = ' '
      write(line,20) i, 'runoff', RNF(i)
      call lineList%AddLine(line)
      line = ' '
      write(line,20) i, 'withdrawal', WTHDRW(i)
      call lineList%AddLine(line)
    enddo
    !
    return
  end subroutine WriteStressPeriodListData

  subroutine process_options(this)
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    ! local
    ! format
    10 format('In processing LAK input, no file found', &
              ' corresponding to unit number: ',i0)
    !
    this%MaxActiveBnd = NLAKES
    !
    ! If TABLEINPUT has been specified, write tables (for LAK8) for all lakes.
    if (IRDTAB == 1) this%Ntables = NLAKES
    this%SurfDep = SURFDEPTH
    !
    return
  end subroutine process_options

  subroutine AddLake(this, lake)
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    type(LakeType), pointer :: lake
    !
    call AddLakeToList(this%Lakes, lake)
    !
    return
  end subroutine AddLake

  subroutine AddOutlet(this, outlet)
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    type(LakeOutletType), pointer :: outlet
    !
    call AddLakeOutletToList(this%Outlets, outlet)
    !
    return
  end subroutine AddOutlet

  subroutine AddTrib(this, tributary)
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    type(LakeTributaryType), pointer :: tributary
    !
    call AddLakeTributaryToList(this%Tribs, tributary)
    !
    return
  end subroutine AddTrib

  function GetLake(this, idx) result (res)
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    integer, intent(in) :: idx
    type(LakeType), pointer :: res
    !
    res => GetLakeFromList(this%Lakes, idx)
    !
    return
  end function GetLake

  function GetOutlet(this, idx) result (res)
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    integer, intent(in) :: idx
    type(LakeOutletType), pointer :: res
    !
    res => GetOutletFromList(this%Outlets, idx)
    !
    return
  end function GetOutlet

  function GetTributary(this, idx) result (res)
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    integer, intent(in) :: idx
    type(LakeTributaryType), pointer :: res
    !
    res => GetTributaryFromList(this%Tribs, idx)
    !
    return
  end function GetTributary

  subroutine DefineLakes(this)
    ! dummy
    class(LakPackageWriterType) :: this
    ! local
    integer :: i
    type(LakeType), pointer :: newLake => null()
    !
    do i=1,NLAKES
      call ConstructLakeType(newLake, this%Igrid)
      newLake%ModelBasename = this%ModelBasename
      newLake%LakeNum = i
      newLake%Strt = STAGES(i)
      call this%AddLake(newLake)
    enddo
    !
    return
  end subroutine DefineLakes

  subroutine DefineAllConnections(this)
    ! dummy
    class(LakPackageWriterType) :: this
    ! local
    integer :: i, nlak
    type(LakeType), pointer :: lake => null()
    !
    nlak = this%Lakes%Count()
    do i=1,nlak
      lake => this%GetLake(i)
      call lake%DefineConnections()
    enddo
    !
    return
  end subroutine DefineAllConnections

  subroutine DefineAllOutlets(this)
    ! Using IDIV, find diversions from lakes to stream segments.
    ! Define these as lake outlets.
    ! dummy
    class(LakPackageWriterType) :: this
    ! local
    integer :: ilak, irch, isg, j, k, ndiv
    type(LakeOutletType), pointer :: newOutlet => null()
    !
    call SGWF2SFR7PNT(this%Igrid)
    if (.not. associated(NSS)) then
      ! SFR is not used in the IGRID model
      return
    endif
    call SGWF2LAK7PNT(this%Igrid)
    !
    k = 0
    do ilak=1,NLAKES
      ndiv = this%FindNDiv(ilak)
      do j=1,ndiv
        ! Construct new outlet and assign its properties
        call ConstructLakeOutlet(newOutlet)
        newOutlet%rcvrPkgName = SfrPackageName
        newOutlet%IgridRcvr = this%Igrid
        k = k + 1
        newOutlet%iOutletNum = k
        ! Lake number
        newOutlet%LakeIn = ilak  ! Lake number
        ! Assign cOutType, width, roughness, slope, invert elevation
        newOutlet%cOutType = 'MANNING'
        ! Stream segment data
        isg = IDIV(ilak,j)             ! Segment number
        newOutlet%iSegnum = isg
        newOutlet%Width = SEG(9,isg)   ! WIDTH1 in SFR input
        newOutlet%Rough = SEG(16,isg)  ! Manning's roughness (ROUGHCH in SFR)
        newOutlet%dInvert = SEG(8,isg) ! ELEVUP in SFR input
        ! Stream reach data
        irch = this%FindReach1Index(isg)  ! Index for reach 1 of segment isg
        newOutlet%Slope = STRM(2,irch)    ! Slope of reach
        ! Add outlet to list of outlets
        call this%AddOutlet(newOutlet)
      enddo
    enddo
    !
    return
  end subroutine DefineAllOutlets

  subroutine DefineAllTributaries(this)
    ! Using ITRB, find stream-segment tributaries to lakes.
    ! Define these as lake tributaries.
    ! dummy
    class(LakPackageWriterType) :: this
    ! local
    integer :: ilak, isg, j, k, ntribs
    type(LakeTributaryType), pointer :: newTrib => null()
    !
    call SGWF2SFR7PNT(this%Igrid)
    if (.not. associated(NSS)) then
      ! SFR is not used in the IGRID model
      return
    endif
    call SGWF2LAK7PNT(this%Igrid)
    !
    k = 0
    do ilak=1,NLAKES
      ntribs = this%FindNDiv(ilak)
      do j=1,ntribs
        ! Construct new tributary and assign its properties
        call ConstructLakeTributary(newTrib)
        newTrib%provPkgName = SfrPackageName
        newTrib%IgridProv = this%Igrid
        k = k + 1
        newTrib%iTribNum = k
        ! Lake number
        newTrib%LakeOut = ilak  ! Lake number
        ! Stream segment data
        isg = ITRB(ilak,j)             ! Segment number
        newTrib%iSegnum = isg
        ! Add tributary to list of tributaries
        call this%AddTrib(newTrib)
      enddo
    enddo
    !
    return
  end subroutine DefineAllTributaries

  function FindNDiv(this,ilak) result (ndiv)
    ! dummy
    class(LakPackageWriterType) :: this
    integer, intent(in) :: ilak
    integer :: ndiv
    ! local
    integer :: i
    !
    ndiv = 0
    if (associated(NSS)) then
      loop: do i=1,NSS
        if (IDIV(ilak,i) > 0) then
          ndiv = i
        endif
      enddo loop
    endif
    !
    return
  end function FindNDiv

  function FindNTribs(this,ilak) result (ntribs)
    ! dummy
    class(LakPackageWriterType) :: this
    integer, intent(in) :: ilak
    integer :: ntribs
    ! local
    integer :: i
    !
    ntribs = 0
    if (associated(NSS)) then
      loop: do i=1,NSS
        if (ITRB(ilak,i) > 0) then
          ntribs = i
        endif
      enddo loop
    endif
    !
    return
  end function FindNTribs

  logical function AnyDivs(this)
    ! dummy
    class(LakPackageWriterType) :: this
    ! local
    integer :: j
    !
    AnyDivs = .false.
    if (associated(NSS)) then
      loop: do j=1,NSS
        if (IDIVAR(1,j) < 0) then
          AnyDivs = .true.
          exit loop
        endif
      enddo loop
    endif
    !
    return
  end function AnyDivs

  logical function AnyTribs(this)
    ! dummy
    class(LakPackageWriterType) :: this
    ! local
    integer :: j
    !
    AnyTribs = .false.
    if (associated(NSS)) then
      loop: do j=1,NSS
        if (IOTSG(j) < 0) then
          AnyTribs = .true.
          exit loop
        endif
      enddo loop
    endif
    !
    return
  end function AnyTribs

  function FindReach1Index(this,isg) result (idx1)
    ! Find index in ISTRM, STRM for reach 1 of segment isg
    ! dummy
    class(LakPackageWriterType) :: this
    integer, intent(in) :: isg
    integer :: idx1, nrch
    ! local
    integer :: i
    !
    nrch = size(ISTRM,2)
    idx1 = 0
    loop: do i=1,nrch
      if (ISTRM(4,i) == isg) then
        if (ISTRM(5,i) == 1) then
          idx1 = i
          exit loop
        endif
      endif
    enddo loop
    !
    return
  end function FindReach1Index

  subroutine write_lakes(this)
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    ! local
    integer :: i, iu, ncon, nlak
    type(LakeType), pointer :: lake
    double precision :: strt
    ! formats
    5  format()
    10 format(2x,a,2x,G15.8)
    20 format(2x,a,2x,a)
    30 format(2x,i0,2x,g15.8,2x,i0)
    50 format(2x,a)
    60 format(a)
    !
    iu = this%fileobj%IUnit
    nlak = this%Lakes%Count()
    !
    write(iu,5)
    write(iu,60)'BEGIN PACKAGEDATA'
    write(iu,50)'#   lakeno    strt        nlakeconn'
    !
    do i=1,nlak
      lake => this%GetLake(i)
      strt = lake%Strt
      ncon = lake%Connections%Count()
      write(iu,30) i,strt,ncon
    enddo
    !
    write(iu,60)'END PACKAGEDATA'
    !
    return
  end subroutine write_lakes

  subroutine write_tables(this)
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    ! local
    integer :: i, iu, nlak
    type(LakeType), pointer :: lake
    ! formats
    5  format()
    20 format(2x,i0,2x,a,2x,a)
    50 format(2x,a)
    60 format(a)
    !
    if (this%Ntables == 0) return
    !
    iu = this%fileobj%IUnit
    nlak = this%Lakes%Count()
    !
    write(iu,5)
    write(iu,60)'BEGIN TABLES'
    write(iu,50)'#   lakeno              ctabname'
    !
    do i=1,nlak
      lake => this%GetLake(i)
      call lake%WriteBathFile()
      write(iu,20) lake%LakeNum, '  TAB6 FILEIN', trim(lake%TableFile)
    enddo
    !
    write(iu,60)'END TABLES'
    !
    return
  end subroutine write_tables

  subroutine write_connections(this)
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    ! local
    integer :: i, iconn, iu, j, ncon, nlak, lakenum
    integer :: irow, jcol, klay
    character(len=10) :: claktype
    double precision :: bleak, belev, telev, clen, cwid
    type(LakeType), pointer           :: lake => null()
    type(LakeConnectionType), pointer :: conn => null()
    ! formats
    5  format()
    10 format(5(2x,i0),2x,a,5(2x,g14.7))
    20 format(2x,a,2x,a)
    50 format(2x,a)
    60 format(a)
    70 format(2x,i0)
    80 format(i0)
    !
    iu = this%fileobj%IUnit
    write(iu,5)
    write(iu,60)'BEGIN CONNECTIONDATA'
    write(iu,60)'  # lakeno iconn cellid claktype     bedleak         belev' // &
       '           telev           connlen         connwidth'
    !
    nlak = this%Lakes%Count()
    do i=1,nlak
      lake => this%GetLake(i)
      lakenum = lake%LakeNum
      ncon = lake%Connections%Count()
      do j=1,ncon
        conn => lake%GetConnection(j)
        iconn = conn%Iconn
        irow = conn%Irow
        jcol = conn%Jcol
        klay = conn%Klay
        claktype = conn%CLakType
        bleak = conn%BedLeak
        belev = conn%Belev
        telev = conn%Telev
        clen = conn%ConnLen
        cwid = conn%ConnWidth
        write(iu,10) lakenum, iconn, klay, irow, jcol, trim(claktype), &
                    bleak, belev, telev, clen, cwid
      enddo
    enddo
    !
    write(iu,60)'END CONNECTIONDATA'
    !
    return
  end subroutine write_connections

  subroutine write_outlets(this)
    implicit none
    ! dummy
    class(LakPackageWriterType) :: this
    ! local
    integer :: i, iu, lakein, lakeout, nout, outletno
    character(len=10) :: outtype
    double precision :: invert, width, rough, slope
    type(LakeOutletType), pointer :: outlet => null()
    ! formats
    5  format()
    10 format(3(2x,i0),2x,a,4(2x,g15.8))
    20 format(2x,a,2x,a)
    50 format(2x,a)
    60 format(a)
    70 format(2x,i0)
    80 format(i0)
    !
    nout = this%Outlets%Count()
    if (nout == 0) return
    !
    iu = this%fileobj%IUnit
    write(iu,5)
    write(iu,60)'BEGIN OUTLETS'
    write(iu,60)'  #   outletno lakein lakeout couttype invert width rough slope'
    do i=1,nout
      outlet => this%GetOutlet(i)
      outletno = outlet%iOutletNum
      lakein = outlet%LakeIn
      lakeout = outlet%LakeOut
      outtype = outlet%cOutType
      invert = outlet%dInvert
      width = outlet%Width
      rough = outlet%Rough
      slope = outlet%Slope
      write(iu,10) outletno, lakein, lakeout, trim(outtype), &
                  invert, width, rough, slope
    enddo
    !
    write(iu,60)'END OUTLETS'
    !
    return
  end subroutine write_outlets

  subroutine BuildModelMovers(this)
    ! dummy
    class(LakPackageWriterType) :: this
    !
    ! ***Let SFR handle LAK -> SFR water movers***
    !call this%BuildLAKtoSFRmovers()
    ! ***Let SFR handle SFR -> LAK water movers***
    !call this%BuildSFRtoLAKmovers()
    !
    return
  end subroutine BuildModelMovers

!  subroutine BuildLAKtoSFRmovers(this)
!    ! Define a water mover for each outlet in this%Outlets
!    implicit none
!    ! dummy
!    class(LakPackageWriterType) :: this
!    ! local
!    integer :: i, idProv, idRec, noutlets
!    integer :: igridProv, igridRcvr
!    type(MoverType), pointer :: newMover => null()
!    type(LakeOutletType), pointer :: outlet => null()
!    character(len=LINELENGTH) :: ermsg
!    character(len=9) :: mvrType
!    character(len=LENPACKAGENAME) :: provPkgName, recPkgName
!    character(len=LENMODELNAME) :: provModelName, recModelName
!    type(ModelPackageType), pointer :: modelPack => null()
!    type(SfrReachType), pointer :: reach => null()
!    type(SfrSegmentType), pointer :: segmnt =>  null()
!    !
!    mvrType = 'FACTOR'
!    provPkgName = this%PackageName
!    provModelName = this%ModelPack%ModelName
!    igridProv = this%Igrid
!    !
!    ! Cycle through all outlets defined in this LAK package
!    noutlets = this%Outlets%Count()
!    do i=1,noutlets
!      outlet => this%GetOutlet(i)
!      recPkgName = outlet%rcvrPkgName
!      modelPack => GetModelPack(recPkgName)
!      recModelName = modelPack%ModelName
!      idProv = outlet%iOutletNum  ! outlet number is the provider ID
!       = outlet%iSegnum      ! segment number is the receiver ID   ****NO**** -- Reach number is the receiver ID !!!
!       idRec = 
!      igridRcvr = outlet%IgridRcvr
!      ! Construct a LAK -> SFR water mover for this outlet
!      call ConstructWaterMover(newMover, mvrType, provModelName, recModelName, &
!                               provPkgName, recPkgName, idProv, idRec, &
!                               igridProv, igridRcvr, 'LAK', 'SFR', DONE)
!      call this%AddMover(newMover)
!    enddo
!    !
!    return
!  end subroutine BuildLAKtoSFRmovers

!  subroutine BuildSFRtoLAKmovers(this)
!    ! Define a water mover for each tributary in this%Tribs
!    implicit none
!    ! dummy
!    class(LakPackageWriterType) :: this
!    ! local
!    integer :: i, idProv, idRec, j, ntribs
!    integer :: igridProv, igridRcvr
!    type(MoverType), pointer :: newMover => null()
!    type(LakeTributaryType), pointer :: trib => null()
!    character(len=LINELENGTH) :: ermsg
!    character(len=9) :: mvrType
!    character(len=LENMODELNAME) :: provModelName, recModelName
!    character(len=LENPACKAGENAME) :: provPkgName, recPkgName
!    type(ModelPackageType), pointer :: modelPack => null()
!    !
!    mvrType = 'FACTOR'
!    recPkgName = this%PackageName
!    recModelName = this%ModelPack%ModelName
!    igridRcvr = this%Igrid
!    !
!    ntribs = this%Tribs%Count()
!    ! Cycle through all tributaries
!    ! Nah--Let SFR do this
!    tribloop: do i=1,ntribs
!      trib => this%GetTributary(i)
!      idProv = trib%iSegnum  ! segment number is the provider ID
!      idRec = trib%LakeOut   ! lake number is the receiver ID
!      igridProv = trib%IgridProv
!      provPkgName = trib%provPkgName
!      modelPack => GetModelPack(provPkgName)
!      provModelName = modelPack%ModelName
!      !! Construct a SFR -> LAK water mover for this tributary
!      !call ConstructWaterMover(newMover, mvrType, provModelName, recModelName, &
!      !                         provPkgName, recPkgName, idProv, &
!      !                         idRec, igridProv, igridRcvr, 'SFR', &
!      !                         'LAK', DONE)
!      !call this%AddMover(newMover)
!    enddo tribloop
!    !
!    return
!  end subroutine BuildSFRtoLAKmovers

  function MyType(this) result (ctype)
    ! dummy
    class(LakPackageWriterType) :: this
    character(len=LENCTYPE) :: ctype
    !
    ctype = 'LakPackageWriterType'
    !
    return
  end function MyType

end module LakPackageWriterModule

module lakjunk
!
!  subroutine AddDiversion(this, diversion)
!    implicit none
!    ! dummy
!    class(LakPackageWriterType) :: this
!    type(SfrDiversionType), pointer :: diversion
!    !
!    call AddDiversionToList(this%Diversions, diversion)
!    !
!    return
!  end subroutine AddDiversion
!
!  function GetDiversion(this, idx) result (res)
!    implicit none
!    ! dummy
!    class(LakPackageWriterType) :: this
!    integer, intent(in) :: idx
!    type(SfrDiversionType), pointer :: res
!    !
!    res => GetDiversionFromList(this%Diversions, idx)
!    !
!    return
!  end function GetDiversion
!
!  subroutine DefineReachDiversions(this)
!    implicit none
!    ! dummy
!    class(LakPackageWriterType) :: this
!    ! local
!    !
!    integer :: i, idiv, iprior, irch, n, ndiv, nrch, &
!               rchnumDs, rchnumUs, segnumDs, segnumUs
!    type(SfrReachType), pointer :: rchDs => null()
!    type(SfrReachType), pointer :: rchUs => null()
!    type(SfrReachType), pointer :: rchLast => null()
!    type(SfrReachType), pointer :: rchFirst => null()
!    type(SfrDiversionType), pointer :: div => null()
!    type(SfrSegmentType), pointer :: segDs => null()
!    type(SfrSegmentType), pointer :: segUs => null()
!    !
!    ! Iterate through all diversions defined for package
!    ndiv = this%Diversions%Count()
!    do idiv=1,ndiv
!      div => this%GetDiversion(idiv)
!      ! Get segment from which flow is diverted
!      segnumUs = div%SegnumUs
!      segUs => this%GetSegment(segnumUs)
!      ! Get segment receiving diverted flow
!      segnumDs = div%SegnumDs
!      segDs => this%GetSegment(segnumDs)
!      ! Get reach from which flow is diverted
!      nrch = segUs%SegReaches%Count()
!      rchUs => this%GetReach(nrch)
!      rchnumUs = rchUs%rnopkg
!      ! Get reach receiving diverted flow
!      rchDs => segDs%GetReach(1)
!      rchnumDs = rchDs%rnopkg
!      ! The reach from which flow is diverted is identified
!      ! by "rno" in the reach_diversions block.
!      div%Rno = rchnumUs
!      ! The reach receiving diverted flow is identified
!      ! by "iconr" in the reach_diversions block.
!      div%Iconr = rchnumDs
!      ! Assign the priority and, if iprior=-2, the flow fraction
!      iprior = div%Iprior
!      if (iprior == -2) then
!        ! Value entered as FLOW is fraction of available flow
!        ! in upstream segment.
!        segUs%Ustrf = segUs%Flow  ! should be a fraction between 0.0 and 1.0
!      else
!        ! By default, set ustrf for reach receiving diverted flow to 0.0.
!        ! Activation of diversion will be handled in PERIOD blocks.
!        rchDs%ustrf = 0.0
!        ! Assign diversion flow
!        div%DivFlow = segDs%Flow
!      endif
!      ! Add diversion to list of diversions belonging to upstream reach
!      call rchUs%AddDiversion(div)
!    enddo
!    !
!    return
!  end subroutine DefineReachDiversions
!
!  subroutine write_diversions(this)
!    implicit none
!    ! dummy
!    class(LakPackageWriterType) :: this
!    ! local
!    integer :: i, iu, j, ndiv, nrch
!    logical :: need
!    type(SfrDiversionType), pointer :: div => null()
!    type(SfrReachType), pointer :: reach => null()
!    ! formats
!    5  format()
!    10 format(4(2x,i0))
!    60 format(a)
!    !
!    need = .false.
!    nrch = this%Reaches%Count()
!    do i=1,nrch
!      reach => this%GetReach(i)
!      if (.not. reach%KeepReachActive) cycle
!      ndiv = reach%Diversions%Count()
!      if (ndiv > 0) then
!        need = .true.
!        exit
!      endif
!    enddo
!    !
!    if (need) then
!      iu = this%fileobj%IUnit
!      ! Write BEGIN REACH_DIVERSIONS
!      write(iu,5)
!      write(iu,60)'BEGIN REACH_DIVERSIONS'
!      !
!      write(iu,60)'  # rno idv iconr iprior'
!      do i=1,nrch
!        reach => this%GetReach(i)
!        if (.not. reach%KeepReachActive) cycle
!        ndiv = reach%Diversions%Count()
!        if (ndiv > 0) then
!          do j=1,ndiv
!            div => reach%GetDiversion(j)
!            if (associated(div)) then
!              write(iu,10)i, j, div%Iconr, div%Iprior
!            endif
!          enddo
!        endif
!      enddo
!      !
!      ! Write END
!      write(iu,60)'END REACH_DIVERSIONS'
!    endif
!    !
!    return
!  end subroutine write_diversions

!  subroutine ReviseReachNetwork(this)
!    ! To be run after masteridomain has been defined, taking into
!    ! account areas deactivated in child grid areas.
!    implicit none
!    ! dummy
!    class(LakPackageWriterType) :: this
!    ! local
!    integer :: i, icon, icr, idiv, ncon, ndiv, nreaches
!    integer :: iconabs, irch, krch
!    integer :: klay, irow, jcol
!    integer, allocatable, dimension(:) :: newReach
!    type(SfrReachType), pointer :: reach => null()
!    type(SfrReachType), pointer :: creach => null()
!    type(SfrDiversionType), pointer :: diversion => null()
!    !
!    nreaches = this%Reaches%Count()
!    allocate(newReach(nreaches))
!    newReach = 0
!    krch = 0
!    do i=1,nreaches
!      reach => GetReachFromList(this%Reaches, i)
!      if (associated(reach)) then
!        klay = reach%klay
!        irow = reach%irow
!        jcol = reach%jcol
!        if (masteridomain(jcol,irow,klay) == 0) then
!          ! Deactivate this reach
!          reach%KeepReachActive = .false.
!          reach%newReachNum = 0
!          ndiv = reach%Diversions%Count()
!          ! Deactivate all diversions from this reach
!          do idiv=1,ndiv
!            diversion => GetDiversionFromList(reach%Diversions, idiv)
!            if (associated(diversion)) then
!              diversion%KeepDiversionActive = .false.
!            endif
!          enddo
!          ! Remove all connections to this reach from other reaches
!          do icr=1,nreaches
!            if (icr == i) cycle
!            creach => GetReachFromList(this%Reaches, icr)
!            if (associated(creach)) then
!              if (.not. creach%KeepReachActive) cycle
!              ncon = size(creach%iconn)
!              connectionsloop: do icon=1,ncon
!                ! If creach connects to this reach, remove the
!                ! connection by removing the creach%iconn element
!                ! that creates the connection.
!                iconabs = abs(creach%iconn(icon))
!                if (iconabs == i) then
!                  ! Remove connection and exit connections loop
!                  call RemoveElement(creach%iconn, icon)
!                  exit connectionsloop
!                endif
!              enddo connectionsloop
!            endif
!          enddo
!        else
!          krch = krch + 1
!          newReach(i) = krch
!          reach%newReachNum = krch
!        endif
!      endif
!    enddo
!    !
!    if (this%MaxActiveBnd /= krch) then
!      ! Some reaches have been removed.
!      ! Renumber the reaches, diversions, and connections.
!      this%MaxActiveBnd = krch
!      do i=1,nreaches
!        reach => GetReachFromList(this%Reaches, i)
!        if (associated(reach)) then
!          ! Renumber the reaches in the active diversions
!          ndiv = reach%Diversions%Count()
!          do idiv=1,ndiv
!            diversion => GetDiversionFromList(reach%Diversions, idiv)
!            if (associated(diversion)) then
!              if (diversion%KeepDiversionActive) then
!                diversion%Iconr = newReach(diversion%Iconr)
!              endif
!            endif
!          enddo
!          ! Renumber connected reaches
!          ncon = size(reach%iconn)
!          connectionsloop2: do icon=1,ncon
!            iconabs = abs(reach%iconn(icon))
!            reach%iconn(icon) = sign(newReach(iconabs), reach%iconn(icon))
!          enddo connectionsloop2
!        endif
!      enddo
!    endif
!    !
!    return
!  end subroutine ReviseReachNetwork
end module lakjunk
