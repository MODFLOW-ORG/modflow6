module UzfPackageWriterModule

  use ConstantsModule, only: DONE, DZERO, LINELENGTH, &
                             MAXCHARLEN, LENMODELNAME, LENPACKAGENAME
  use ConstantsPHMFModule, only: FCINPUT, FCUNKNOWN, LENCTYPE
  use FileTypeModule, only: FileType
  use GLOBAL, only: NPER, IUNIT, NCOL, NROW, NLAY, IOUT, ISSFLG, IBOUND
  use GWFBASMODULE, only: FindHighestActiveLayer, SGWF2BAS7PNT
  use GWFLAKMODULE, only: SGWF2LAK7PNT, NeedLakWaterMover
  use GWFSFRMODULE, only: SGWF2SFR7PNT, NeedSfrWaterMover
  use GwfUzfModule, only: NeedUzfWaterMover, IUZFBND, SGWF2UZF1PNT, SURFDEP, &
                          IRUNBND, NTRAIL, NSETS, VKS, THTI, THTR, THTS, &
                          EPS, IUZFCB1, IUZFCB2, NUZTOP, PETRATE, ROOTDPTH, &
                          WCWILT, FINF, IRUNFLG
  use GwfUzfSubs, only: GWF2UZF1AR, GWF2UZF1RP
  use LineListModule, only: LineListType
  use ListModule, only: ListType
  use ModelPackageModule, only: ModelPackageType, GetModelPack
  use MoverModule, only: MoverType, ConstructWaterMover, AddMoverToList, &
                         GetMoverFromList
  use MvrPackageWriterModule, only: MvrPackageWriterType
  use PackageWriterModule, only: PackageWriterType
  use SfrPackageWriterModule, only: SfrPackageWriterType, &
                                    GetSfrPackageWriterByIgrid
  use SfrReachModule, only: SfrReachType
  use SfrSegmentModule, only: SfrSegmentType
  use SimPHMFModule, only: store_error, ustop

  implicit none

  private
  public :: AllUzfPkgWriters, UzfPackageWriterType

  type(ListType) :: AllUzfPkgWriters

  type, extends(PackageWriterType) :: UzfPackageWriterType
    type(MvrPackageWriterType), pointer :: MvrWriter => null()
    type(SfrPackageWriterType), pointer :: SfrWriter => null()
    logical :: simulate_et = .true. ! simulate ET in UZ
    logical :: linear_gwet = .true. ! original ET fomulation of mf2005
    logical :: square_gwet = .false.
    logical :: simulate_gwseep = .false.
    logical :: unsat_etwc = .false.
    logical :: unsat_etae = .false.
    logical :: multilayer = .false.
    integer :: nuzfcells = 0
    integer :: nuztop = 0
    integer :: kdim = 0
    double precision :: surfdep = DZERO ! one value per sim in UZF1; one per cell in UZF8
    integer, allocatable, dimension(:,:,:) :: IUzfNum
  contains
    ! Public procedures
    procedure, public :: BuildModelMovers
    procedure, public :: MyType
    procedure, public :: ProcessAllocate
    procedure, public :: ProcessStressLoop
    procedure, public :: WriteDimensions
    procedure, public :: WriteOptions
    procedure, public :: WriteStressPeriodListData
    ! Private procedures
    procedure, private :: CountRunoffLakes
    procedure, private :: CountRunoffSegs
    procedure, private :: PopulateIUzfNum
    procedure, private :: process_options
    procedure, private :: UzfCount
    procedure, private :: WriteDataBlock
  end type UzfPackageWriterType

contains

  subroutine ProcessAllocate(this, igrid)
    implicit none
    ! dummy
    class(UzfPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: krunLak, krunSfr
    character(len=MAXCHARLEN) :: ermsg, fname
    ! format
    12 format(a,'_',i0)
    !
    this%Active = .true.
    this%Igrid = igrid
    call this%AllocatePointers()
    this%DefaultBudgetText = 'UZF CELLS'
    this%ipr = 1
    this%IuOrig = IUNIT(55)    ! Specific to UZF package
    this%OrigPkg = 'UZF'
    this%PkgType = 'UZF'
    write(this%PackageName,12)trim(this%PkgType),igrid
    this%fileobj%FCode = FCINPUT
    this%fileobj%FType = 'UZF6'
    this%fileobj%PkgName = this%PackageName
    fname = trim(this%ModelBasename) // '.uzf'
    call this%FileWriterType%InitializeFile(fname, this%fileobj%FType, &
                                            this%PackageName)
    this%IunitBcf = iunit(1)
    this%IunitGwt = iunit(15)
    this%IunitLpf = iunit(23)
    this%IunitHuf = iunit(37)
    this%IunitSfr = iunit(44)
    this%IunitUzf = iunit(55)
    this%IunitUpw = iunit(62)
    call GWF2UZF1AR(this%IuOrig, this%IunitBcf, this%IunitLpf, &
                    this%IunitUpw, igrid)
    !
    this%nuztop = NUZTOP
    call this%PopulateIUzfNum()
    if (IRUNFLG > 0) then
      krunLak = this%CountRunoffLakes()
      krunSfr = this%CountRunoffSegs()
      if (krunLak > 0) then
        this%NeedWaterMover = .true.
        if (.not. associated(NeedLakWaterMover)) then
          ermsg = 'UZF needs information from LAK package, which has not' // &
              ' been read yet.  Please edit name file so that LAK precedes UZF.'
          call store_error(ermsg)
          call ustop()
        endif
        NeedLakWaterMover = .true.
      endif
      if (krunSfr > 0) then
        this%NeedWaterMover = .true.
        if (.not. associated(NeedSfrWaterMover)) then
          ermsg = 'UZF needs information from SFR package, which has not' // &
              ' been read yet.  Please edit name file so that SFR precedes UZF.'
          call store_error(ermsg)
          call ustop()
        endif
        NeedSfrWaterMover = .true.
      endif
    endif
    !
    if (.not. this%Active) then
      this%fileobj%FCode = FCUNKNOWN
      return
    endif
    !
    if (IUZFCB1 > 0) then
      this%ICbc => IUZFCB1
    elseif (IUZFCB2 > 0) then
      this%ICbc => IUZFCB2
    else
      this%ICbc => IUZFCB1
    endif
    !
    return
  end subroutine ProcessAllocate

  subroutine PopulateIUzfNum(this)
    ! IUzfNum will contain UZF cell ID number where UZF is active, 0 elsewhere
    ! dummy
    class(UzfPackageWriterType) :: this
    ! local
    integer :: i, j, k, kk, ki, kl
    !
    kk = 0
    if (this%nuztop == 1 .or. this%nuztop == 2) then
      this%kdim = 1
    elseif (this%nuztop == 3) then
      this%kdim = NLAY
    endif
    allocate(this%IUzfNum(NCOL,NROW,this%kdim))
    kloop: do kl=1,this%kdim
      rowloop: do i=1,NROW
        colloop: do j=1,NCOL
          this%IUzfNum(j,i,kl) = 0
          if (IUZFBND(j,i) <= 0) cycle colloop
          select case (NUZTOP)
          case (1)
            k = 1
            ki = 1
          case (2)
            k = IUZFBND(j,i)
            ki = 1
          case (3)
            if (IBOUND(j,i,kl) > 0) then
              k = kl
              ki = kl
            else
              k = 0
            endif
            !k = FindHighestActiveLayer(i,j)
          end select
          if (k <= 0) cycle colloop
          if (IBOUND(j,i,k) <= 0) cycle colloop
          kk = kk + 1
          this%IUzfNum(j,i,ki) = kk
        enddo colloop
      enddo rowloop
    enddo kloop
    !
    return
  end subroutine PopulateIUzfNum

  subroutine WriteOptions(this)
    ! Overrides PackageWriterType%WriteOptions
    implicit none
    class(UzfPackageWriterType) :: this
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
    if (this%ipr==1) then
      write(iu,50)'PRINT_INPUT'
    endif
    !
    write(iu,50)'PRINT_FLOWS'
    if (associated(this%ICbc)) then
      if (this%ICbc > 0) then
        write(iu,50)'SAVE_FLOWS'
      endif
    endif
    !
    if (this%simulate_et) then
      write(iu,50)'SIMULATE_ET'
    endif
    !
    if (this%linear_gwet) then
      write(iu,50)'LINEAR_GWET'
    endif
    !
    if (this%simulate_gwseep) then
      write(iu,50)'SIMULATE_GWSEEP'
    endif
    !
    !if (this%multilayer) then
    !  write(iu,50)'MULTILAYER'
    !endif
    !
    if (this%NeedWaterMover) then
      write(iu,50)'MOVER'
    endif
    !
    ! Write END Options
    write(iu,60)'END Options'
    !
    return
  end subroutine WriteOptions

  subroutine process_options(this)
    implicit none
    ! dummy
    class(UzfPackageWriterType) :: this
    ! local
    integer :: irl, irs
    ! format
    !
    this%MaxActiveBnd = this%UzfCount()
    !
    this%SurfDep = SURFDEP
    !
    irl = this%CountRunoffLakes()
    irs = this%CountRunoffSegs()
    if (irl > 0 .or. irs > 0) then
      this%simulate_gwseep = .true.
    endif
    !
    if (NUZTOP == 2 .or. NUZTOP == 3) then
      ! 2 -- Recharge to and discharge from layer specified in IUZFBND
      ! 3 -- Recharge to and discharge from highest active cell
      !      in each vertical column
      !this%multilayer = .true.
    endif
    !
    return
  end subroutine process_options

  subroutine WriteDimensions(this)
    implicit none
    class(UzfPackageWriterType) :: this
    ! local
    integer :: iu
    ! formats
    5  format()
    30 format(2x,a,2x,i0)
    60 format(a)
    !
    if (.not. this%NeedDimensionsBlock) return
    !
    this%MaxActiveBnd = this%UzfCount()
    !
    iu = this%fileobj%IUnit
    ! Write BEGIN Dimensions
    write(iu,5)
    write(iu,60)'BEGIN Dimensions'
    !
    write(iu,30)'NUZFCELLS', this%MaxActiveBnd
    write(iu,30)'NTRAILWAVES', NTRAIL
    write(iu,30)'NWAVESETS', NSETS
    !
    ! Write END Dimensions
    write(iu,60)'END Dimensions'
    !
    return
  end subroutine WriteDimensions

  subroutine ProcessStressLoop(this, igrid)
    implicit none
    ! dummy
    class(UzfPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: kper
    logical :: currentA, forceWrite
    character(len=12) :: sstrOption, lastSstrOption
    ! formats
    10 format(' Processing UZF package input...')
    20 format(/,'Processing UZF package for stress period ',i0)
    !
    ! Write Options and Dimensions blocks
    if (.not. this%Active) return
    !
    ! Check to see if UZF specifies potential for flow between UZF and any
    ! lakes and streams. If so, water mover(s) will be needed.
    call SGWF2LAK7PNT(this%Igrid)
    call SGWF2SFR7PNT(this%Igrid)
    call SGWF2UZF1PNT(this%Igrid)
    !
    ! Write the options block
    call this%WriteOptions()
    !
    ! Initially, current block is BlockA; alternate each stress period
    this%CurrentBlock => this%BlockA
    this%PreviousBlock => this%BlockB
    currentA = .true.
    !
    lastSstrOption = ''
!    sstrOption = 'TRANSIENT'
    forceWrite = .true.
    do kper=1,nper
!! Following Rich's example, use TRANSIENT always -- Not sure about this
!      select case (ISSFLG(kper))
!        case (0)
!          sstrOption = 'TRANSIENT'
!        case (1)
!          sstrOption = 'STEADY-STATE'
!      end select
      sstrOption = ''
!
! UZF requires block every stress period, so this logic is commented out.
! Just let forceWrite = true every stress period.
!      if (sstrOption == lastSstrOption) then
!        forceWrite = .false.
!      else
!        forceWrite = .true.
!      endif
      if (kper==1) write(*,10)
      write(iout,20)kper
      ! Read MF2005 input for current stress period
      call GWF2UZF1RP(this%IuOrig, kper, this%IunitSfr, this%Igrid)
      if (.not. this%Active) then
        this%fileobj%FCode = FCUNKNOWN
        return
      endif
      !
      if (kper==1) then
        ! Write Dimensions block
        call this%WriteDimensions()
        ! Write Data block
        call this%WriteDataBlock()
      endif
      !
      ! Write MF6 input for current stress period to a LineList
      call this%CurrentBlock%Clear(.true.)
      call this%WriteStressPeriodListData(this%CurrentBlock)
      ! Write block to MF6 input file if needed
      call this%WriteBlockIfNeeded(kper, sstrOption, forceWrite)
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
      lastSstrOption = sstrOption
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
    class(UzfPackageWriterType) :: this
    type(LineListType), pointer :: lineList
    ! local
    character(len=MAXCHARLEN) :: line
    integer :: i, iu, j, kk, kl, kha
    ! formats
    5 format('#',4x,'iuzno',6x,'finf',13x,'pet',12x,'extdp',12x,'extwc',12x, &
             ' ha  ',12x,'hroot',12x,'rootact')
    10 format(2x,a,2x,i0,2x,a)
    20 format(4x,i0,7(2x,g15.8))
    !
    if (.not. this%Active) return
    call SGWF2BAS7PNT(this%Igrid)
    call SGWF2UZF1PNT(this%Igrid)
    !
    iu = this%fileobj%IUnit
    !
    ! Write uzfid finf pet extdp extwc ha hroot rootact
    write(line,5)
    call lineList%AddLine(line)
    kloop: do kl=1,this%kdim
      rowloop: do i=1,NROW
        colloop: do j=1,NCOL
          kk = this%IUzfNum(j,i,kl)
          if (kk == 0) cycle colloop
          ! Find highest active layer for this row,column location
          kha = FindHighestActiveLayer(i,j)
          if (this%nuztop == 3 .and. kha /= kl) cycle colloop
          ! DZEROes at end are for ha, hroot, and rootact.
          write(line,20)kk, FINF(j,i), PETRATE(j,i), ROOTDPTH(j,i), &
                        WCWILT(j,i), DZERO, DZERO, DZERO
          call lineList%AddLine(line)
        enddo colloop
      enddo rowloop
    enddo kloop
    !
    return
  end subroutine WriteStressPeriodListData

  subroutine WriteDataBlock(this)
    ! dummy
    class(UzfPackageWriterType) :: this
    ! local
    integer :: iuzno, i, j, iu, k, landflag, vertcon
    integer :: kl, kha
    double precision :: vk, thr, ths, thi, ep
    ! formats
    5  format()
    20 format('# iuzno cellid landflag vertcon surfdep',9x,'vks',14x,'thtr',13x, &
              'thts',13x,'thti',14x,'eps')
    30 format(2x,i0, 5(2x,i0), 6(2x,g15.8))
    60 format(a)
    !
    call SGWF2BAS7PNT(this%igrid)
    call SGWF2UZF1PNT(this%Igrid)
    landflag = 1
    vertcon = 0
    !
    ! Write BEGIN Data
    iu = this%fileobj%IUnit
    write(iu,5)
    write(iu,60)'BEGIN PACKAGEDATA'
    write(iu,20)
    !
    ! For each active UZF cell, need to write:
    ! layer row column landflag vertcon surfdep vks thtr thts thti eps
    iuzno = 1
    kloop: do kl=1,this%kdim
      rowloop: do i=1,NROW
        colloop: do j=1,NCOL
          if (this%IUzfNum(j,i,kl) > 0) then
            ! this is an active UZF cell
            ! Get layer index
            select case (NUZTOP)
            case (1)
              ! Recharge to and discharge from only the top model layer
              k = 1
            case (2)
              ! Recharge to and discharge from layer specified in IUZFBND
              k = IUZFBND(j,i)
            case (3)
              ! Find highest active layer for this row,column location
              kha = FindHighestActiveLayer(i,j)
              ! Recharge to and discharge from highest active cell
              ! in each vertical column
              k = kl
              if (kl == kha) then
                landflag = 1
              else
                landflag = 0
              endif
            end select
            if (k == 0) cycle colloop
            vk = VKS(j,i)
            thr = THTR(j,i)
            ths = THTS(j,i)
            thi = THTI(j,i)
            ep = EPS(j,i)
            write(iu,30) iuzno, k, i, j, landflag, vertcon, SURFDEP, vk, thr, ths, thi, ep
            iuzno = iuzno + 1
          endif
        enddo colloop
      enddo rowloop
    enddo kloop
    !
    write(iu,60)'END PACKAGEDATA'
    !
    return
  end subroutine WriteDataBlock

  function MyType(this) result (ctype)
    ! dummy
    class(UzfPackageWriterType) :: this
    character(len=LENCTYPE) :: ctype
    !
    ctype = 'UzfPackageWriterType'
    !
    return
  end function MyType

  integer function UzfCount(this)
    ! dummy
    class(UzfPackageWriterType) :: this
    ! local
    integer :: i, j, k, kk, kl
    !
    call SGWF2UZF1PNT(this%Igrid)
    call SGWF2BAS7PNT(this%Igrid)
    kk = 0
    kloop: do kl=1,this%kdim
      rowloop: do i=1,NROW
        colloop: do j=1,NCOL
          if (this%IUzfNum(j,i,kl) /= 0) then
            k = FindHighestActiveLayer(i, j)
            if (k == 0) cycle colloop
            kk = kk + 1
          endif
        enddo colloop
      enddo rowloop
    enddo kloop
    !
    UzfCount = kk
    !
    return
  end function UzfCount

  integer function CountRunoffLakes(this)
    ! dummy
    class(UzfPackageWriterType) :: this
    ! local
    integer :: i, j, kk, kl
    !
    call SGWF2UZF1PNT(this%Igrid)
    kk = 0
    rowloop: do i=1,NROW
      colloop: do j=1,NCOL
        kloop: do kl=1,this%kdim
          if (this%IUzfNum(j,i,kl) == 0) cycle kloop
          if (IRUNBND(j,i) < 0) then
            kk = kk + 1
            exit kloop
          endif
        enddo kloop
      enddo colloop
    enddo rowloop
    !
    CountRunoffLakes = kk
    !
    return
  end function CountRunoffLakes

  integer function CountRunoffSegs(this)
    ! dummy
    class(UzfPackageWriterType) :: this
    ! local
    integer :: i, j, kk, kl
    !
    call SGWF2UZF1PNT(this%Igrid)
    kk = 0
    rowloop: do i=1,NROW
      colloop: do j=1,NCOL
        kloop: do kl=1,this%kdim
          if (this%IUzfNum(j,i,kl) == 0) cycle kloop
          if (IRUNBND(j,i) > 0) then
            kk = kk + 1
            exit kloop
          endif
        enddo kloop
      enddo colloop
    enddo rowloop
    !
    CountRunoffSegs = kk
    !
    return
  end function CountRunoffSegs

  subroutine BuildModelMovers(this)
    ! Define a water mover for each outlet in this%Outlets
    implicit none
    ! dummy
    class(UzfPackageWriterType) :: this
    ! local
    integer :: i, j, kl
    integer :: igridProv, igridRcvr, idProv, idRec, isegnum
    type(MoverType), pointer :: newMover => null()
    type(SfrSegmentType), pointer :: seg => null()
    type(SfrReachType), pointer :: reach => null()
    character(len=MAXCHARLEN) :: ermsg
    character(len=9) :: mvrType
    character(len=LENPACKAGENAME) :: provPkgName, recLakPkgName, &
                                     recSfrPkgName
    character(len=LENMODELNAME) :: provModelName, recModelName
    logical :: needSfrWriter
    ! formats
    10 format(a,i0)
    20 format('SFR package for IGRID = ',i0,' not found. Please reorder', &
              ' name file entries so that SFR file precedes UZF file.')
    !
    if (.not. this%NeedWaterMover) return
    !
    mvrType = 'FACTOR'
    provPkgName = this%PackageName
    provModelName = this%ModelPack%ModelName
    recModelName = provModelName
    igridProv = this%Igrid
    igridRcvr = igridProv
    write(recLakPkgName,10)'LAK_',this%Igrid
    write(recSfrPkgName,10)'SFR_',this%Igrid
    !
    ! Need access to SfrPackageWriter for this grid if any IRUNBND > 0
    needSfrWriter = .false.
    outerloop: do i=1,NROW
      innerloop: do j=1,NCOL
        if (IRUNBND(j,i) > 0) then
          needSfrWriter = .true.
          exit outerloop
        endif
      enddo innerloop
    enddo outerloop
    if (needSfrWriter) then
      ! Get SfrPackageWriter for this igrid
      this%SfrWriter => GetSfrPackageWriterByIgrid(this%Igrid)
      if (.not. associated(this%SfrWriter)) then
        write(ermsg,20)this%Igrid
        call store_error(ermsg)
        call ustop()
      endif
    endif
    !
    ! Cycle through grid and build movers as needed
    rowloop: do i=1,NROW
      colloop: do j=1,NCOL
        kloop: do kl=1,this%kdim
          if (this%IUzfNum(j,i,kl) == 0) cycle colloop
          if (IRUNBND(j,i) < 0) then
            ! Construct a UZF -> LAK water mover for this cell
            idProv = this%IUzfNum(j,i,kl) ! UZF cell ID number is the provider ID
            idRec = ABS(IRUNBND(j,i))  ! Lake number is the receiver ID
            call ConstructWaterMover(newMover, mvrType, provModelName, &
                     recModelName, provPkgName, recLakPkgName, idProv, idRec, &
                     igridProv, igridRcvr, 'UZF', 'LAK', DONE)
          elseif (IRUNBND(j,i) > 0) then
            ! Construct a UZF -> SFR water mover for this cell
            idProv = this%IUzfNum(j,i,kl) ! UZF cell ID number is the provider ID
            ! Reach number will be the receiver ID. First find segment number.
            isegnum = IRUNBND(j,i)
            ! Get reach number of first reach in this segment from SFR
            seg => this%SfrWriter%GetSegment(isegnum)  ! Get segment indicated in IRUNBND
            reach => seg%GetReach(1)                   ! Get first reach in this segment
            idRec = reach%rnopkg                       ! Get reach number in package
            call ConstructWaterMover(newMover, mvrType, provModelName, &
                     recModelName, provPkgName, recSfrPkgName, idProv, idRec, &
                     igridProv, igridRcvr, 'UZF', 'SFR', DONE)
          endif
          call this%AddMover(newMover)
        enddo kloop
      enddo colloop
    enddo rowloop
    !
    return
  end subroutine BuildModelMovers

end module UzfPackageWriterModule

module uzfjunk
! Variables and arrays of UZF1 (SZ = saturated zone; LS = land surface)
!
! Read by GWF2UZF1AR:
! ITHTRFLG = 1 when SPECIFYTHTR is used (THTR will be read for 1st transient SP)
! ITHTIFLG = 1 when SPECIFYTHTI is used (THTI will be read for 1st SP)
! Iseepsupress = 1 when NOSURFLEAK is used (surface leakage not simulated)
! NUZTOP     - similar to NRCHOP, signifies where recharge/discharge is simulated.
!        : 1 - RECHARGE IN UZF TO TOP LAYER ONLY
!        : 2 - RECHARGE IN UZF TO SPECIFIED NODE IN EACH VERTICAL COLUMN
!        : 3 - RECHARGE IN UZF TO HIGHEST ACTIVE NODE IN EACH VERTICAL COLUMN
! IUZFOPT: 0 - flow will not be routed through the UZ, goes directly to SZ.
!          1 - vertical K will be specified in UZF1 input using array VKS.
!          2 - vertical K will be specified in either BCF or LPF input.
! IRUNFLG: 0 - GW discharge to LS is removed from model.
!          > 0 - GW discharge goes to LAK3 or SFR2.
! IETFLG:  0 - ET will not be simulated.
!          not 0 - ET will be simulated.
! IUZFCB1 = IUZFB11 - >0 is unit # for writing unformatted cell-by-cell recharge, ET, discharge using UBUDSV
! IUZFCB2 = IUZFB22 - >0 is unit # for writing unformatted cell-by-cell recharge, ET, discharge using UBDSV3
! NTRAIL (default = 1) - # trailing waves
! NSETS (default = 1) - # wave sets
! NUZGAG - # cells for printing detailed information
! SURFDEP - average undulation depth within a finite-difference cell
! IUZFBND(ncol,nrow) - defines cells where recharge and discharge will be simulated
! NUMACTIVE - # cells connected to stream segments
! IRUNBND(ncol,nrow) - stream segment (>0) or lake number (<0) that can receive discharge or excess infiltration.
! VKS(ncol,nrow) - sat'd vertical K of UZ
! EPS(ncol,nrow) - Brooks-Corey epsilon
! THTS(ncol,nrow) - saturated water content of UZ
! THTR(ncol,nrow) - residual water content
! THTI(ncol,nrow) - initial water content
! IUZLIST(1,) - row index for cell where detailed information is printed
!        (2,) - column index for cell where detailed information is printed
!        (3,) - unit # of output file; or, if <0, output is summed over all model cells
!        (4,) - flag for type of expanded listing
!
! Read (or calculated) by GWF2UZF1RP:
! NUZF - Flag for reading infiltration rates; if <0, use rates from previous SP
! FINF(ncol,nrow) - infiltration rates
! EXCESPP(ncol,nrow) - excess infiltration greater than VKS
! RECHSAVE(ncol,nrow) - FINF used for SFR interaction somehow
! NUZF - Flag for reading ET demand rates; if <0, use rates from previous SP
! PETRATE(ncol,nrow) - ET demand rates
! NUZF - Flag for reading ET extinction depths
! ROOTDPTH(ncol,nrow) - ET extinction depths
! NUZF - Flag for reading extinction water contents
! WCWILT - Extinction water contents
!
end module uzfjunk

