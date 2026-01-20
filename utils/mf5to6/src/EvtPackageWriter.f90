module EvtPackageWriterModule

  use ConstantsModule, only: MAXCHARLEN, DZERO
  use ConstantsPHMFModule, only: FCINPUT, LENCTYPE
  use GLOBAL, only: NPER, NCOL, NROW, NLAY, IUNIT, IOUT, DELR, DELC, IBOUND
  use GlobalVariablesModule, only: echo, masteridomain, verbose
  use GWFETSMODULE, only: NETSOP, IETSCB, NPETS, IETSPF, ETSR, ETSX, &
                          ETSS, PXDP, PETM, IETS, NETSEG
  use GwfEtsSubs, only: GWF2ETS7AR, GWF2ETS7RP
  use GWFEVTMODULE, only: NEVTOP, IEVTCB, NPEVT, IEVTPF, EVTR, EXDP, &
                          SURF, IEVT
  use GwfEvtSubs, only: GWF2EVT7AR, GWF2EVT7RP
  use LineListModule, only: LineListType, same_lines
  use PackageWriterModule, only: PackageWriterType
  use SimPHMFModule, only: ustop, store_error, store_warning, store_note
  use UtilitiesModule, only: Write1Drel, Write2Drel, Write3Drel, &
                             ConstantInt2D, ConstantReal2D, &
                             BuildArrayFormat
  use utl7module, only: U1DREL, U2DREL, & !UBDSV1, UBDSV2, UBDSVA, &
                        urword, URDCOM, & !UBDSV4, UBDSVB, &
                        ULSTRD

  type, extends(PackageWriterType) :: EvtPackageWriterType
    logical                          :: HighestActive
    integer, pointer, dimension(:,:) :: iet => null()
    integer, pointer                 :: netop => null()
    integer                          :: nseg = 0
    double precision, pointer, dimension(:,:)    :: rate => null()
    double precision, pointer, dimension(:,:)    :: exdp => null()
    double precision, pointer, dimension(:,:)    :: surf => null()
    character(len=1000) :: evtfmat      ! format for printing EVT list data
    type(LineListType), pointer      :: BlockRateA => null()
    type(LineListType), pointer      :: BlockRateB => null()
    type(LineListType), pointer      :: CurrentRateBlock => null()
    type(LineListType), pointer      :: PreviousRateBlock => null()
    type(LineListType), pointer      :: BlockExdpA => null()
    type(LineListType), pointer      :: BlockExdpB => null()
    type(LineListType), pointer      :: CurrentExdpBlock => null()
    type(LineListType), pointer      :: PreviousExdpBlock => null()
    type(LineListType), pointer      :: BlockSurfA => null()
    type(LineListType), pointer      :: BlockSurfB => null()
    type(LineListType), pointer      :: CurrentSurfBlock => null()
    type(LineListType), pointer      :: PreviousSurfBlock => null()
    type(LineListType), pointer      :: BlockSegA => null()
    type(LineListType), pointer      :: BlockSegB => null()
    type(LineListType), pointer      :: CurrentSegBlock => null()
    type(LineListType), pointer      :: PreviousSegBlock => null()
  contains
    procedure :: MyType
    procedure :: ProcessAllocate
    procedure :: ProcessStressLoop
    procedure :: WriteOptions
    procedure :: WriteDimensions
    procedure :: WriteSegArrayDataIfNeeded
    procedure :: WriteListDataIfNeeded
    procedure :: WriteStressPeriodListData
    procedure :: AllocateEvtPointers
    procedure :: BuildEvtArrayFormat
  end type EvtPackageWriterType

contains

  subroutine ProcessAllocate(this, igrid)
    implicit none
    ! dummy
    class(EvtPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: i, j
    character(len=MAXCHARLEN) :: fname
    ! format
    !
    select case (this%source)
    case ('EVT')
      this%IuOrig = IUNIT(5)        ! Specific to EVT package
    case ('ETS')
      this%IuOrig = IUNIT(39)       ! Specific to ETS package
    case default
      call store_error('EVT source not valid: "'//trim(this%source)//'"')
    end select
    call this%AllocatePointers()
    call this%AllocateEvtPointers()
    this%Active = .true.
    this%DefaultBudgetText = 'EVAPORATION'
    ! Ned todo: not sure that this will work for LGR--need to look at LGR source
    this%fileobj%FCode = FCINPUT
    !
    this%fileobj%FType = 'EVT6'
    this%PkgType = 'EVT'
!    this%PackageName = 'EVT'
    fname = trim(this%ModelBasename) // '.evt'
    call this%FileWriterType%InitializeFile(fname, this%fileobj%FType, &
                                            this%PackageName)
    select case (this%source)
    case ('EVT')
      call GWF2EVT7AR(this%IuOrig, igrid)
      this%netop => NEVTOP
      this%iet => IEVT
      this%nseg = 1
    case ('ETS')
      call GWF2ETS7AR(this%IuOrig, igrid)
      this%netop => NETSOP
      this%iet => IETS
      this%nseg = NETSEG
    end select
    !
    this%HighestActive = .false.
    if (this%netop == 3) this%HighestActive = .true.
    !
    ! NStressDim is number of values after k,i,j required to specify boundary stress
    this%NStressDim = 0     ! EVT uses this as flag to indicate array-type package(?)
    !
    ! NVALP1 is rlist index for element following minimum required data
    this%nvalp1 = 0         ! Does not apply to EVT
    !
    ! Assign package-specific pointers
    select case (this%source)
    case ('EVT')
      this%rate => EVTR
      this%exdp => EXDP
      this%surf => SURF
      this%ipr => IEVTPF
      this%ICbc => IEVTCB
    case ('ETS')
      this%rate => ETSR
      this%exdp => ETSX
      this%surf => ETSS
      this%ipr => IETSPF
      this%ICbc => IETSCB
    end select
    allocate(this%NBndPeriod)
    this%NBndPeriod = NROW     ! EVT uses this to print array(?)
    this%Aux => null()
    this%IPB => null()
    this%MaxActiveBnd = NCOL * NROW
    this%NBNDVL => null()
    this%NAux = 0
    this%nstart = 1
    this%nstop = NCOL          ! EVT uses this to print array(?)
    !
    ! Allocate and populate CellArea array
    allocate(this%CellArea(NCOL, NROW))
    do i=1,NROW
      do j=1,NCOL
        this%CellArea(j,i) = DELC(i)*DELR(j)
      enddo
    enddo
    allocate(this%work(NCOL, NROW))
    !
    ! define format for printing boundary data
    if (this%nseg == 1) then
      call BuildArrayFormat(this%DataWidth, 7, 100, this%fmat, this%FieldsPerLine)
    else
      call this%BuildEvtArrayFormat() !this%DataWidth, 7, 100, this%evtfmat, this%FieldsPerLine)
    endif
    !
    ! Write Options block
    call this%WriteOptions()
    !
    ! Dimensions block is not needed if READASARRAYS is used.
    ! READASARRAYS is used unless ETS function uses more than 1 segment.
    if (this%nseg > 1) then
      call this%WriteDimensions()
    endif
    !
    return
  end subroutine ProcessAllocate

  subroutine BuildEvtArrayFormat(this)
    ! dummy
    class(EvtPackageWriterType) :: this
    ! local
    integer :: i
    !
    this%evtfmat = '(1x,3(1x,i0),3(2x,g12.5)'
    if (this%nseg > 1) then
      do i=1,this%nseg-1
        this%evtfmat = trim(this%evtfmat) // ',2(2x,g12.5)'
      enddo
    endif
    this%evtfmat = trim(this%evtfmat) // ')'
    !
    return
  end subroutine BuildEvtArrayFormat

  subroutine ProcessStressLoop(this, igrid)
    implicit none
    ! dummy
    class(EvtPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    logical :: constant, currentA, writing
    integer :: i, iprni, iprnr, iu, j, kper, nlval
    ! formats
    1 format(/,'Processing EVT package for stress period ',i0)
    20 format(2x,a)
    30 format(/,'BEGIN PERIOD ',i0)
    40 format('END PERIOD')
    50 format(4x,a,2x,i0)
    80 format(20(i3,1x))
    !
    iu = this%fileobj%IUnit
    if (echo) then
      iprni = 4
      iprnr = 12
    else
      iprni = -4
      iprnr = -12
    endif
    ! Initially, current blocks are Block*A; alternate each stress period
    this%CurrentRateBlock => this%BlockRateA
    this%CurrentExdpBlock => this%BlockExdpA
    this%CurrentSurfBlock => this%BlockSurfA
    this%PreviousRateBlock => this%BlockRateB
    this%PreviousExdpBlock => this%BlockExdpB
    this%PreviousSurfBlock => this%BlockSurfB
    this%CurrentSegBlock => this%BlockSegA
    this%PreviousSegBlock => this%BlockSegB
    currentA = .true.
    !
    do kper=1,nper  ! stress period loop
      write(iout,1)kper
      ! Read MF2005 input for current stress period
      select case (this%source)
      case ('EVT')
        if (kper==1) write(*,*)'Processing EVT package input...'
        call GWF2EVT7RP(this%IuOrig, igrid)
      case ('ETS')
        if (kper==1) write(*,*)'Processing ETS package input...'
        call GWF2ETS7RP(this%IuOrig, igrid)
      end select
      !
      ! if kper = 1, write IEVT; if kper > 1, don't.
      writing = .false.
      if (this%nseg == 1) then
        if (kper==1) then
          writing = .true.
          write(iu,30)kper
          select case (this%netop)
          case (1)  ! evapotranspiration is only from layer 1
            constant = .true.
            nlval = 1
          case (2)  ! evapotranspiration is vertically distributed according to IEVT array
            ! determine if IEVT is constant
            call ConstantInt2D(NCOL, NROW, this%iet, constant, nlval)
          case (3)  ! evapotranspiration is from highest active cell in each vertical column
            constant = .true.
            nlval = 1
          end select
          ! write IEVT
          if (this%netop == 2) then
            write(iu,20)'IEVT'
            if (constant) then
              write(iu,50)'CONSTANT', nlval
            else
              write(iu,50)'INTERNAL  FACTOR  1  IPRN ',iprni
              do i=1,nrow
                write(iu,80)(this%iet(j,i),j=1,ncol)
              enddo
            endif
          endif
        endif
      endif
      !
      call this%CurrentRateBlock%Clear(.true.)
      call this%CurrentExdpBlock%Clear(.true.)
      call this%CurrentSurfBlock%Clear(.true.)
      call this%CurrentSegBlock%Clear(.true.)
      !
      if (this%nseg == 1) then
        !
        ! Determine if max ET rate array differs from previous stress period,
        ! and write array if needed.
        call this%WriteArrayDataIfNeeded(iu, kper, this%rate, this%CurrentRateBlock, &
                                         this%PreviousRateBlock, writing, .true., 'RATE')
        !
        ! Deal with EXDP similarly.
        call this%WriteArrayDataIfNeeded(iu, kper, this%exdp, this%CurrentExdpBlock, &
                                         this%PreviousExdpBlock, writing, .false., 'DEPTH')
        ! Deal with SURF similarly.
        call this%WriteArrayDataIfNeeded(iu, kper, this%surf, this%CurrentSurfBlock, &
                                         this%PreviousSurfBlock, writing, .false., 'SURFACE')
        if (this%source == 'ETS') then
          ! Deal with segment data similarly
          call this%WriteSegArrayDataIfNeeded(iu, kper, this%CurrentSegBlock, &
                                            this%PreviousSegBlock, writing)
        endif
        if (writing) write(iu,40)
      else
        call this%WriteListDataIfNeeded(iu, kper, this%CurrentBlock, this%PreviousBlock)
      endif
      !
      ! switcheroo
      if (currentA) then
        this%CurrentRateBlock => this%BlockRateB
        this%PreviousRateBlock => this%BlockRateA
        this%CurrentExdpBlock => this%BlockExdpB
        this%PreviousExdpBlock => this%BlockExdpA
        this%CurrentSurfBlock => this%BlockSurfB
        this%PreviousSurfBlock => this%BlockSurfA
        this%CurrentSegBlock => this%BlockSegB
        this%PreviousSegBlock => this%BlockSegA
        currentA = .false.
      else
        this%CurrentRateBlock => this%BlockRateA
        this%PreviousRateBlock => this%BlockRateB
        this%CurrentExdpBlock => this%BlockExdpA
        this%PreviousExdpBlock => this%BlockExdpB
        this%CurrentSurfBlock => this%BlockSurfA
        this%PreviousSurfBlock => this%BlockSurfB
        this%CurrentSegBlock => this%BlockSegA
        this%PreviousSegBlock => this%BlockSegB
        currentA = .true.
      endif
    enddo  ! stress period loop
    !
    call this%BlockA%Clear(.true.)
    call this%BlockB%Clear(.true.)
    call this%BlockExdpA%Clear(.true.)
    call this%BlockExdpB%Clear(.true.)
    call this%BlockRateA%Clear(.true.)
    call this%BlockRateB%Clear(.true.)
    call this%BlockSegA%Clear(.true.)
    call this%BlockSegB%Clear(.true.)
    call this%BlockSurfA%Clear(.true.)
    call this%BlockSurfB%Clear(.true.)
    !
    return
  end subroutine ProcessStressLoop

  subroutine WriteOptions(this)
    ! Override if needed
    implicit none
    class(EvtPackageWriterType) :: this
    ! local
    integer :: i, iu
    character(len=MAXCHARLEN) :: auxline
    ! formats
    5  format()
20     format(2x,a,2x,a)
30     format(2x,a,1x,i0)
    50 format(2x,a)
    60 format(a)
    !
    iu = this%fileobj%IUnit
    ! Write BEGIN Options
    write(iu,5)
    write(iu,60)'BEGIN Options'
    !
    if (this%nseg == 1) then
      write(iu,50)'READASARRAYS'
    endif
    !
    ! Aux variable names
    if (this%NAux > 0) then
      auxline = 'AUXILIARY'
      do i=1,this%NAux
        auxline = trim(auxline) // ' ' // trim(this%Aux(i))
      enddo
      write(iu,20)trim(auxline)
    endif
    !
!    if (echo) then
!      if (this%ipr==1) then
        write(iu,50)'PRINT_INPUT'
!      endif
!    endif
    !
    if (this%ICbc > 0) then
      write(iu,50)'SAVE_FLOWS'
    elseif (this%ICbc < 0) then
      write(iu,50)'PRINT_FLOWS'
    endif
    !
    if (.not. this%HighestActive) then
      write(iu,50)'FIXED_CELL'
    endif
    !
    ! Write END Options
    write(iu,60)'END Options'
    !
    return
  end subroutine WriteOptions

  subroutine WriteDimensions(this)
    implicit none
    class(EvtPackageWriterType) :: this
    ! local
    integer :: iu
    ! formats
    5  format()
    30 format(2x,a,2x,i0)
    60 format(a)
    !
    iu = this%fileobj%IUnit
    ! Write BEGIN Dimensions
    write(iu,5)
    write(iu,60)'BEGIN Dimensions'
    !
    !write(iu,30)'MAXBOUND',this%MaxActiveBnd
    !
    ! Number of segments
    if (this%nseg > 1) then
      write(iu,30)'MAXBOUND',NCOL*NROW
      write(iu,30)'NSEG',this%nseg
    endif
    !
    ! Write END Dimensions
    write(iu,60)'END Dimensions'
    !
    return
  end subroutine WriteDimensions

  subroutine WriteStressPeriodListData(this, lineList)
    implicit none
    ! dummy
    class(EvtPackageWriterType) :: this
    type(LineListType), pointer :: lineList
    ! local
    integer :: i, j, k, kold, n
    character(len=MAXCHARLEN) :: line, msg
    ! format
    10 format(a,' boundary removed at (',i0,',',i0,',',i0,')')
    !
    ! Write stress period data for MF6
    ! Loop through rows and columns to generate list of EVT input
    rowloop: do i=1,NROW
      colloop: do j=1,NCOL
        if (NETSOP == 1) then
          k = 1
        else
          kold = IETS(j,i)
          k = this%Layptr(kold)
        endif
        if (IBOUND(j,i,k) == 0) cycle colloop
        if (associated(masteridomain)) then
          if (masteridomain(j,i,k)==0) then
            if (verbose) then
              msg = ''
              write(msg,10)trim(this%PkgType),k,i,j
  !            write(*,*)trim(msg)
              call store_note(msg)
            endif
            cycle colloop
          endif
        endif
        line = ' '
        ! list data are: cellid, surface, rate, depth, pxdp(nseg-1), petm(nseg-1)
        if (this%rate(j,i) /= DZERO) then
          write(line,this%evtfmat)k, i, j, this%surf(j,i), &
                this%rate(j,i) / this%CellArea(j,i), &
                this%exdp(j,i), (PXDP(j,i,n),n=1,this%nseg-1), &
                (PETM(j,i,n),n=1,this%nseg-1)
          call lineList%AddLine(line)
        endif
      enddo colloop
    enddo rowloop
    !
    return
  end subroutine WriteStressPeriodListData

  subroutine WriteSegArrayDataIfNeeded(this, iu, kper, currentBlock, &
                                       previousBlock, writing)
    implicit none
    ! dummy
    class(EvtPackageWriterType), intent(inout)    :: this
    integer, intent(in)                        :: iu, kper
    type(LineListType), pointer, intent(inout) :: currentBlock, previousBlock
    logical, intent(inout)                     :: writing
    ! local
    integer :: i, ii, inc, iprni, iprnr, j, n, nlines
    logical :: constantpe, constantpx, needToWrite
    logical :: convertToFlux
    double precision    :: rvalpe, rvalpx
    character(len=10) :: label
    character(len=MAXCHARLEN) :: line
    ! formats
    10 format(a)
    20 format(2x,a)
    30 format(/,'BEGIN PERIOD ',i0)
    60 format(4x,a,2x,g15.8)
    70 format(4x,a,2x,i0)
    !
    if (echo) then
      iprni = 4
      iprnr = 12
    else
      iprni = -4
      iprnr = -12
    endif
    ! Determine if segment data differ from previous stress period,
    ! and write arrays if needed.
    ! Start by writing all data for this stress period to currentBlock
    convertToFlux = .false.
    do j=1,this%nseg-1
      call this%WriteStressPeriodArrayData(PXDP(:,:,j), currentBlock, convertToFlux)
      call this%WriteStressPeriodArrayData(PETM(:,:,j), currentBlock, convertToFlux)
    enddo
    needToWrite = .not. same_lines(currentBlock, previousBlock)
    nlines = currentBlock%CountLines()
    if (nlines > 0) then
      inc = ( nlines / (this%nseg-1) ) / 2
      if (needToWrite) then
        if (.not. writing) then
          write(iu,30)kper
          writing = .true.
        endif
        ! for each segment intersection, write PXDP and PETM
        ii = 1   ! first line of an array in currentBlock
        n = inc  ! last line of an array in currentBlock
        do j=1,this%nseg-1
          ! Write PXDP
          call ConstantReal2D(NCOL, NROW, PXDP(:,:,j), constantpx, rvalpx)
          label = 'PXDP'
          write(iu,20) trim(label)
          if (constantpx) then
            write(iu,60)'CONSTANT', rvalpx
          else
            ! write array control line and array to input file
            write(iu,70)'INTERNAL  FACTOR  1.0  IPRN ',iprnr
            do i=ii,n
              call currentBlock%GetLine(i, line)
              write(iu,10)trim(line)
            enddo
          endif
          ! Write PETM
          ii = n + 1
          n = n + inc
          call ConstantReal2D(NCOL, NROW, PETM(:,:,j), constantpe, rvalpe)
          label = 'PETM'
          write(iu,20) trim(label)
          if (constantpe) then
            write(iu,60)'CONSTANT', rvalpe
          else
            ! write array control line and array to input file
            write(iu,70)'INTERNAL  FACTOR  1.0  IPRN ',iprnr
            do i=ii,n
              call currentBlock%GetLine(i, line)
              write(iu,10)trim(line)
            enddo
          endif
          ii = n + 1
          n = n + inc
        enddo
      endif
    endif
    !
    return
  end subroutine WriteSegArrayDataIfNeeded

  subroutine WriteListDataIfNeeded(this, iu, kper, currentBlock, previousBlock)
    implicit none
    ! dummy
    class(EvtPackageWriterType), intent(inout) :: this
    integer, intent(in)                        :: iu, kper
    type(LineListType), pointer, intent(inout) :: currentBlock, previousBlock
    ! local
    integer :: i, n
    logical :: needToWrite
    character(len=MAXCHARLEN) :: line
    character(len=60) :: label = '# cellid  surface  rate  depth  pxdp(nseg-1)  petm(nseg-1)'
    ! formats
    10 format(a)
    20 format(2x,a)
    30 format(/,'BEGIN PERIOD ',i0)
    40 format('END PERIOD ')
    60 format(4x,a,2x,g15.8)
    80 format(4x,a,2x,i0)
    !
    call this%WriteStressPeriodListData(currentBlock)
    needToWrite = .not. same_lines(currentBlock, previousBlock)
    if (needToWrite) then
      write(iu,30)kper
      write(iu,20) trim(label)
      n = currentBlock%CountLines()
      do i=1,n
        call currentBlock%GetLine(i, line)
        write(iu,10)trim(line)
      enddo
      write(iu,40)
    endif
    !
    return
  end subroutine WriteListDataIfNeeded

  subroutine AllocateEvtPointers(this)
    implicit none
    ! dummy
    class(EvtPackageWriterType) :: this
    !
    ! Allocate pointers to handle RATE
    allocate(this%BlockRateA)
    allocate(this%BlockRateB)
    allocate(this%CurrentRateBlock)
    allocate(this%PreviousRateBlock)
    call this%BlockRateA%InitializeLineList()
    call this%BlockRateB%InitializeLineList()
    call this%CurrentRateBlock%InitializeLineList()
    call this%PreviousRateBlock%InitializeLineList()
    !
    ! Allocate pointers to handle EXDP
    allocate(this%BlockExdpA)
    allocate(this%BlockExdpB)
    allocate(this%CurrentExdpBlock)
    allocate(this%PreviousExdpBlock)
    call this%BlockExdpA%InitializeLineList()
    call this%BlockExdpB%InitializeLineList()
    call this%CurrentExdpBlock%InitializeLineList()
    call this%PreviousExdpBlock%InitializeLineList()
    !
    ! Allocate pointers to handle SURF
    allocate(this%BlockSurfA)
    allocate(this%BlockSurfB)
    allocate(this%CurrentSurfBlock)
    allocate(this%PreviousSurfBlock)
    call this%BlockSurfA%InitializeLineList()
    call this%BlockSurfB%InitializeLineList()
    call this%CurrentSurfBlock%InitializeLineList()
    call this%PreviousSurfBlock%InitializeLineList()
    !
    ! Allocate pointers to handle segment blocks
    allocate (this%BlockSegA)
    allocate (this%BlockSegB)
    allocate (this%CurrentSegBlock)
    allocate (this%PreviousSegBlock)
    call this%BlockSegA%InitializeLineList()
    call this%BlockSegB%InitializeLineList()
    call this%CurrentSegBlock%InitializeLineList()
    call this%PreviousSegBlock%InitializeLineList()
    !
    return
  end subroutine AllocateEvtPointers

  function MyType(this) result (ctype)
    ! dummy
    class(EvtPackageWriterType) :: this
    character(len=LENCTYPE) :: ctype
    !
    ctype = 'EvtPackageWriterType'
    !
    return
  end function MyType

end module EvtPackageWriterModule

