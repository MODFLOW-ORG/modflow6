module PackageWriterModule

  use ConstantsModule, only: LENFTYPE, LENPACKAGENAME, LINELENGTH, &
                             MAXCHARLEN
  use ConstantsPHMFModule, only: FCUNKNOWN, LENCTYPE
  use FileListModule, only: FileListType
  use FileTypeModule, only: FileType, ConstructFileType
  use FileWriterModule, only: FileWriterType
  use GLOBAL, only: NPER, NCOL, NROW, NLAY, IUNIT, IOUT, DELR, DELC, IBOUND
  use GlobalVariablesModule, only: echo, masteridomain, verbose
  use LineListModule, only: LineListType, same_lines
  use ListModule, only: ListType
  use ModelPackageModule, only: ModelPackageType, ConstructModelPackageType
  use MoverModule, only: MoverType, AddMoverToList
  use ObsWriterModule, only: ObsWriterType
  use SimPHMFModule, only: store_error, store_note, ustop
  use SimListVariablesModule, only: SimMovers
  use UtilitiesModule, only: ConstantReal2D

  implicit none

  private
  public :: PackageWriterType, CastAsPackageWriterClass

  type, extends(FileWriterType) :: PackageWriterType
    integer, pointer :: ICbc => null()       ! CBC unit number
    integer, pointer :: ipr => null()        ! print flag
    integer, pointer :: NBndPeriod => null() ! number of boundaries used in stress period
    integer, pointer :: IPB => null()
    integer, pointer :: NBNDVL => null()
    integer, pointer :: Nlaynew => null()
    integer, dimension(:), pointer :: Layptr => null()
    integer, dimension(:,:,:), pointer :: Idomain => null()
    integer :: Igrid = 0
    integer :: IuOrig = 0        ! unit number for package input file of MF2005 model
    integer :: MaxActiveBnd = -1
    integer :: NAux = 0          ! number of auxiliary variables
    integer :: nstart = 4        ! default: first element in rlist following k,i,j
    integer :: nstop = 0         ! last element in rlist to be read
    integer :: nvalp1            ! Rlist index for element following minimum required data
    integer :: FieldsPerLine     ! number of fields available per line of array output
    integer :: DataWidth = 14    ! characters available for each array data value
    integer :: NStressDim  ! Number of values after k,i,j required to specify boundary stress
    integer :: Iouts = 0  ! For solute output? Never assigned in mf2005
    integer :: Istcb2 = 0
    integer :: IunitBcf = 0
    integer :: IunitGwt = 0
    integer :: IunitHuf = 0
    integer :: IunitLak = 0
    integer :: IunitLpf = 0
    integer :: IunitSfr = 0
    integer :: IunitUpw = 0
    integer :: IunitUzf = 0
    integer :: iulist = 0
    integer :: Nsol = 0
    character(len=MAXCHARLEN)     :: ModelBasename
    character(len=LENPACKAGENAME) :: PackageName = ''
    character(len=LENPACKAGENAME) :: DefaultBudgetText = ''
    character(len=6) :: source = ''
    character(len=100) :: fmat      ! format for printing boundary data
    character(len=MAXCHARLEN) :: TsFile = ''
    character(len=4) :: OrigPkg = ''
    character(len=16), pointer, dimension(:)  :: Aux => null()
    double precision, pointer, dimension(:,:) :: rlist => null()
    double precision, pointer, dimension(:,:) :: CellArea => null()
    double precision, pointer, dimension(:,:) :: work => null()
    logical :: PrintInput = .true.
    logical :: PrintFlows = .true.
    logical :: SaveFlows = .false.
    logical :: NeedDimensionsBlock = .true.
    logical :: NeedWaterMover = .false.
    logical :: Newton = .false.
    type(FileListType), pointer     :: Mf6Files => null()
    type(ListType), pointer         :: ModelMovers => null()
    type(ModelPackageType), pointer :: ModelPack => null()
    type(LineListType), pointer     :: BlockA => null()
    type(LineListType), pointer     :: BlockB => null()
    type(LineListType), pointer     :: CurrentBlock => null()
    type(LineListType), pointer     :: PreviousBlock => null()
    ! for observations
    logical :: ObsActive = .false.
    class(ObsWriterType), pointer   :: PkgObsWriter => null()
  contains
    procedure :: ProcessAllocate
    procedure :: ProcessStressLoop
    procedure :: MyType
    procedure :: AddMover
    procedure :: AllocatePointers
    procedure :: BuildModelMovers
    procedure :: BuildSimMovers
    procedure :: GetBudgetText
    procedure :: SameType
    procedure :: WriteOptions
    procedure :: WriteDimensions
    procedure :: WriteStressPeriodListData
    procedure :: WriteStressPeriodArrayData
    procedure :: WriteBlockIfNeeded
    procedure :: WriteArrayDataIfNeeded
  end type PackageWriterType

contains

  !=============================================================================
  ! Procedures that need to be overridden by types
  ! that inherit from PackageWriterType.

  subroutine ProcessAllocate(this, igrid)
    class(PackageWriterType) :: this
    integer, intent(in) :: igrid
    stop 'Error: PackageWriterType%ProcessAllocate not overridden in child type.'
  end subroutine ProcessAllocate

  subroutine ProcessStressLoop(this, igrid)
    class(PackageWriterType) :: this
    integer, intent(in) :: igrid
    stop 'Error: PackageWriterType%ProcessStressLoop not overridden in child type.'
  end subroutine ProcessStressLoop

  function MyType(this) result(ctype)
    class(PackageWriterType) :: this
    character(len=30) :: ctype
    ctype = ' '
    stop 'Error: PackageWriterType%MyType not overridden in child type.'
  end function MyType

  !=============================================================================

  subroutine AllocatePointers(this)
    implicit none
    ! dummy
    class(PackageWriterType) :: this
    ! local
    type(FileType), pointer :: newFile
    !
    allocate(this%ICbc)
    this%ICbc = 0
    allocate(this%ipr)
    allocate(this%NBndPeriod)
    allocate(this%IPB)
    allocate(this%NBNDVL)
    allocate(this%Nlaynew)
    !
    allocate(this%BlockA)
    allocate(this%BlockB)
    allocate(this%CurrentBlock)
    allocate(this%PreviousBlock)
    !
    call this%BlockA%InitializeLineList()
    call this%BlockB%InitializeLineList()
    call this%CurrentBlock%InitializeLineList()
    call this%PreviousBlock%InitializeLineList()
    !
    if (.not. associated(this%fileobj)) then
      call ConstructFileType(newFile)
      this%fileobj => newFile
    endif
    !
    if (.not. associated(this%PkgObsWriter)) then
      allocate(this%PkgObsWriter)
    endif
    !
    return
  end subroutine AllocatePointers

  subroutine WriteOptions(this)
    ! Override if needed
    implicit none
    class(PackageWriterType) :: this
    ! local
    integer :: i, iu
    character(len=200) :: ermsg
    character(len=MAXCHARLEN) :: auxline
    ! formats
    5  format()
    20 format(2x,a,2x,a,2x,a)
    50 format(2x,a)
60     format(a)
70 format(a,1x,i0,1x,a)
    !
    iu = this%fileobj%IUnit
    ! Write BEGIN Options
    write(iu,5)
    write(iu,60)'BEGIN Options'
    !
!    ! Aux variable names
!    if (this%NAux > 0) then
!      do i=1,this%NAux
!        write(iu,20) 'AUXILIARY', trim(this%Aux(i))
!      enddo
!    endif
    ! Aux variable names
    if (this%NAux > 0) then
      auxline = 'AUXILIARY'
      do i=1,this%NAux
        auxline = trim(auxline) // ' ' // trim(this%Aux(i))
      enddo
      write(iu,20)trim(auxline)
    endif
    !
    if (echo) then
      if (this%ipr==1) then
        this%PrintInput = .true.
      endif
    endif
    !
    if (associated(this%ICbc)) then
      if (this%ICbc > 0) then
        this%SaveFlows = .true.
      elseif (this%ICbc < 0) then
        this%PrintFlows = .true.
      endif
    endif
    !
    ! Print_Input
    if (echo) then
      if (this%PrintInput) write(iu,50)'PRINT_INPUT'
    endif
    !
    ! Print_Flows
    if (this%PrintFlows) write(iu,50)'PRINT_FLOWS'
    !
    ! Save_Flows
    if (this%SaveFlows) write(iu,50)'SAVE_FLOWS'
    !!
    !! Newton
    !if (.not. this%Newton) write(iu,50)'NO_NEWTON'
    !
    if (this%ObsActive) then
      if (this%PkgObsWriter%fileobj%FName == '') then
        write(ermsg,70)'Programmer error: OBS6 file for package opened on unit', &
                       this%fileobj%IUnit, 'is unnamed.'
        call store_error(ermsg)
        call ustop()
      endif
      write(iu,20)'OBS6', 'FILEIN', trim(this%PkgObsWriter%fileobj%FName)
    endif
    !
    if (this%TsFile /= '') then
      write(iu,20)'TS6 FILEIN',trim(this%TsFile)
    endif
    !
    ! Write END Options
    write(iu,60)'END Options'
    !
    return
  end subroutine WriteOptions

  subroutine WriteDimensions(this)
    implicit none
    class(PackageWriterType) :: this
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
    write(iu,30)'MAXBOUND',this%MaxActiveBnd
    !
    ! Write END Dimensions
    write(iu,60)'END Dimensions'
    !
    return
  end subroutine WriteDimensions

  subroutine WriteBlockIfNeeded(this, kper, option, forceWrite)
    implicit none
    ! dummy
    class(PackageWriterType) :: this
    integer, intent(in) ::  kper
    character(len=*), optional :: option
    logical, optional :: forceWrite
    ! local
    integer :: kcurrent
    logical :: needToWrite, forceWriteLocal
    integer :: i, iu
    character(len=200) :: line
    ! formats
    5  format()
    10 format(a)
    20 format('BEGIN PERIOD ',i0)
    25 format('BEGIN PERIOD ',i0,2x,a)
    30 format('END PERIOD')
    !
    ! Determine if a block is needed for current stress period,
    ! and if it is, write it.
    forceWriteLocal = .false.
    if (present(forceWrite)) then
      forceWriteLocal = forceWrite
    endif
    !
    if (forceWriteLocal) then
      needToWrite = .true.
      kcurrent = this%CurrentBlock%CountLines()
    else
      needToWrite = .false.
      kcurrent = this%CurrentBlock%CountLines()
      if (kper==1) then
        if (kcurrent>0) needToWrite = .true.
      else
        ! kper > 1
        needToWrite = .not. same_lines(this%BlockA, this%BlockB)
      endif
    endif
    !
    if (needToWrite) then
      ! write block of stress period input for MF6
      iu = this%fileobj%IUnit
      write(iu,5)
      if (present(option)) then
        write(iu,25)kper,trim(option)
      else
        write(iu,20)kper
      endif
      do i=1,kcurrent
        call this%CurrentBlock%GetLine(i, line)
        write(iu,10)trim(line)
      enddo
      write(iu,30)
    endif
    !
    return
  end subroutine WriteBlockIfNeeded

  subroutine WriteStressPeriodListData(this, lineList)
    ! Override if needed
    implicit none
    ! dummy
    class(PackageWriterType) :: this
    type(LineListType), pointer :: lineList
    ! local
    integer :: i, ii, j, k, kold, n
    character(len=MAXCHARLEN) :: line, msg
    ! format
    10 format(a,' boundary removed at (',i0,',',i0,',',i0,')')
    !
    ! Write stress period data for MF6
    do ii=1,this%NBndPeriod
      kold = nint(this%rlist(1,ii))
      k = this%Layptr(kold)
      i = nint(this%rlist(2,ii))
      j = nint(this%rlist(3,ii))
      if (associated(masteridomain)) then
        if (masteridomain(j,i,k)==0) then
          if (verbose) then
            msg = ''
            write(msg,10)trim(this%PkgType),k,i,j
!            write(*,*)trim(msg)
            call store_note(msg)
          endif
          cycle
        endif
      else
        if (IBOUND(j,i,kold)==0) then
          if (verbose) then
            msg = ''
            write(msg,10)trim(this%PkgType),k,i,j
!            write(*,*)trim(msg)
            call store_note(msg)
          endif
          cycle
        endif
      endif
      line = ' '
      write(line,this%fmat)k, i, j, (this%rlist(n,ii),n=this%nstart,this%nstop)
      call lineList%AddLine(line)
    enddo
    !
    return
  end subroutine WriteStressPeriodListData

  subroutine WriteStressPeriodArrayData(this, array, lineList, convertToFlux)
    implicit none
    ! dummy
    class(PackageWriterType)    :: this
    double precision, dimension(NCOL, NROW) :: array
    type(LineListType), pointer :: lineList
    logical                     :: convertToFlux
    ! local
    integer :: fpl, i, j, jstart, jend, n, nrem, nrep
    character(len=MAXCHARLEN) :: line
    !
    ! Write array data to lineList
    do i=1,NROW
      do j=1,NCOL
        ! If specified, convert flow rates to flux rates by dividing by cell area
        if (convertToFlux) then
          this%work(j,i) = array(j,i)/this%CellArea(j,i)
        else
          this%work(j,i) = array(j,i)
        endif
      enddo
    enddo
    !
    fpl = this%FieldsPerLine
    nrep = NCOL / fpl
    nrem = NCOL - (nrep * this%FieldsPerLine)
    do i=1,NROW
      jstart = 1
      do n=1,nrep
        line = ''
        jend = jstart + fpl - 1
        write(line,this%fmat)(this%work(j,i),j=jstart,jend)
        call lineList%AddLine(line)
        jstart = jend + 1
      enddo
      if (nrem > 0) then
        jend = jstart + nrem - 1
        line = ''
        write(line,this%fmat)(this%work(j,i),j=jstart,jend)
        call lineList%AddLine(line)
      endif
    enddo
    !
    return
  end subroutine WriteStressPeriodArrayData

  subroutine WriteArrayDataIfNeeded(this, iu, kper, array, currentBlock, previousBlock, &
                                    writing, convertToFlux, label)
    implicit none
    ! dummy
    class(PackageWriterType), intent(inout)    :: this
    integer, intent(in)                        :: iu, kper
    double precision, dimension(NCOL, NROW), intent(inout) :: array
    type(LineListType), pointer, intent(inout) :: currentBlock, previousBlock
    logical, intent(in)                        :: convertToFlux
    logical, intent(inout)                     :: writing
    character(len=*), intent(in)               :: label
    ! local
    integer :: i, iprnr, n
    logical :: constant, needToWrite
    double precision    :: rval
    character(len=MAXCHARLEN) :: line
    ! formats
    10 format(a)
    20 format(2x,a)
    30 format(/,'BEGIN PERIOD ',i0)
    60 format(4x,a,2x,g15.8)
    80 format(4x,a,2x,i0)
    !
    if (echo) then
      iprnr = 12
    else
      iprnr = -12
    endif
    ! Determine if max ET rate array differs from previous stress period,
    ! and write array if needed.
    call this%WriteStressPeriodArrayData(array, currentBlock, convertToFlux)
    needToWrite = .not. same_lines(currentBlock, previousBlock)
    if (needToWrite) then
      call ConstantReal2D(NCOL, NROW, this%work, constant, rval)
      ! If kper=1 and rval is zero, no need to write
      if (kper == 1 .and. constant .and. rval == 0.0) needToWrite = .false.
      if (needToWrite) then
        if (.not. writing) then
          write(iu,30)kper
          writing = .true.
        endif
        write(iu,20) trim(label)
        if (constant) then
          write(iu,60)'CONSTANT', rval
        else
          ! write array control line and array to input file
          write(iu,80)'INTERNAL  FACTOR  1.0  IPRN ',iprnr
          n = currentBlock%CountLines()
          do i=1,n
            call currentBlock%GetLine(i, line)
            write(iu,10)trim(line)
          enddo
        endif
      endif
    endif
    !
    return
  end subroutine WriteArrayDataIfNeeded

  function CastAsPackageWriterClass(obj) result(res)
    ! Cast an unlimited polymorphic object as class(PackageWriterType)
    implicit none
    class(*), pointer, intent(inout) :: obj
    class(PackageWriterType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    class is (PackageWriterType)
      res => obj
    end select
    return
  end function CastAsPackageWriterClass

  subroutine GetBudgetText(this, text)
    ! Return text that will be associated with this
    ! package in MF6 volume budget tables.
    implicit none
    ! dummy
    class(PackageWriterType) :: this
    character(len=*), intent(out) :: text
    !
    if (this%PackageName /= '') then
      text = this%PackageName
    elseif (this%DefaultBudgetText /= '') then
      text = this%DefaultBudgetText
    else
      text = '(undefined)'
    endif
    !
    return
  end subroutine GetBudgetText

  subroutine AddMover(this, mover)
    implicit none
    ! dummy
    class(PackageWriterType) :: this
    type(MoverType), pointer :: mover
    !
    if (mover%ProvModelName == mover%RecModelName) then
      call AddMoverToList(this%ModelMovers, mover)
    else
      call AddMoverToList(SimMovers, mover)
    endif
    !
    return
  end subroutine AddMover

  subroutine BuildModelMovers(this)
    ! Most packages do not support water movers, so this procedure does nothing.
    ! Packages that support water movers must override this procedure.
    implicit none
    ! dummy
    class(PackageWriterType) :: this
    !
    return
  end subroutine BuildModelMovers

  subroutine BuildSimMovers(this)
    ! Most packages do not support water movers, so this procedure does nothing.
    ! Packages that support water movers must override this procedure.
    implicit none
    ! dummy
    class(PackageWriterType) :: this
    !
    return
  end subroutine BuildSimMovers

  logical function SameType(this, pkgwriter)
    ! dummy
    class(PackageWriterType) :: this
    class(PackageWriterType) :: pkgwriter
    !
    character(len=LENCTYPE) :: thistype, othertype
    !
    SameType = .false.
    thistype = this%MyType()
    othertype = pkgwriter%MyType()
    if (thistype == othertype) then
      SameType = .true.
    endif
    !
    return
  end function SameType

end module PackageWriterModule
