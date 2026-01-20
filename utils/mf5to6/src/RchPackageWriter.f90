module RchPackageWriterModule

  use ConstantsModule, ONLY: MAXCHARLEN
  use ConstantsPHMFModule, only: LENCTYPE
  use GLOBAL, only: NPER, IUNIT, NCOL, NROW, NLAY, IOUT, DELR, DELC
  use GlobalVariablesModule, only: echo
  use GWFRCHMODULE, only: NRCHOP, IRCHCB, NPRCH, IRCHPF, RECH, IRCH
  use GwfRchSubs, only: GWF2RCH7AR, GWF2RCH7RP
  use LineListModule, only: LineListType, same_lines
  use PackageWriterModule, only: PackageWriterType
  use SimPHMFModule, only: ustop, store_error, store_warning
  use UtilitiesModule, only: Write1Drel, Write2Drel, Write3Drel, &
                             ConstantInt2D, ConstantReal2D, &
                             BuildArrayFormat
  use utl7module, only: U1DREL, U2DREL,  &
                        urword, URDCOM,  &
                        ULSTRD

  type, extends(PackageWriterType) :: RchPackageWriterType
    logical :: HighestActive
  contains
    procedure :: MyType
    procedure :: ProcessAllocate
    procedure :: ProcessStressLoop
    procedure :: WriteOptions
!    procedure :: WriteStressPeriodData
  end type RchPackageWriterType

contains

  subroutine ProcessAllocate(this, igrid)
    use ConstantsPHMFModule, only: FCINPUT
    implicit none
    ! dummy
    class(RchPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: i, j, nfield
    character(len=MAXCHARLEN) :: fname
    ! format
    !
    call this%AllocatePointers()
    this%Active = .true.
    this%DefaultBudgetText = 'RECHARGE'
    ! Ned todo: not sure that this will work for LGR--need to look at LGR source
    this%IuOrig = IUNIT(8)        ! Specific to RCH package
    this%fileobj%FCode = FCINPUT
    !
    this%fileobj%FType = 'RCH6'
    this%PkgType = 'RCH'
    fname = trim(this%ModelBasename) // '.rch'
    call this%FileWriterType%InitializeFile(fname, this%fileobj%FType, &
                                            this%PackageName)
    call GWF2RCH7AR(this%IuOrig, igrid)
    !
    this%HighestActive = .false.
    if (NRCHOP == 3) this%HighestActive = .true.
    !
    ! NStressDim is number of values after k,i,j required to specify boundary stress
    this%NStressDim = 0     ! RCH uses this as flag to indicate array-type package(?)
    !
    ! NVALP1 is rlist index for element following minimum required data
    this%nvalp1 = 0         ! Does not apply to RCH
    !
    ! Assign package-specific pointers
    this%rlist => RECH
    this%ipr => IRCHPF
    this%ICbc => IRCHCB
    allocate(this%NBndPeriod)
    this%NBndPeriod = NROW     ! RCH uses this to print array(?)
    this%Aux => null()
    this%IPB => null()
    this%MaxActiveBnd = NCOL * NROW
    this%NBNDVL => null()
    this%NAux = 0
    this%nstart = 1
    this%nstop = NCOL          ! RCH uses this to print array(?)
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
    nfield = MAXCHARLEN / 15
    ! define format for printing boundary data
    call BuildArrayFormat(this%DataWidth, 7, 100, this%fmat, this%FieldsPerLine)
    !
    ! Write Options and Dimensions blocks
    call this%WriteOptions()
    this%NeedDimensionsBlock = .false.
    call this%WriteDimensions()
    !
    return
  end subroutine ProcessAllocate

  subroutine ProcessStressLoop(this, igrid)
    implicit none
    ! dummy
    class(RchPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    double precision :: rechval
    logical :: constant, currentA, needToWrite, writing
    integer :: i, iprni, iprnr, iu, j, kper, n, nlval
    character(len=MAXCHARLEN) :: line
    ! formats
    1 format(/,'Processing RCH package for stress period ',i0)
    10 format(a)
    20 format(2x,a)
    30 format(/,'BEGIN PERIOD ',i0)
    40 format('END PERIOD')
    50 format(4x,a,2x,i0)
    60 format(4x,a,2x,g15.8)
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
    ! Initially, current block is BlockA; alternate each stress period
    this%CurrentBlock => this%BlockA
    this%PreviousBlock => this%BlockB
    currentA = .true.
    !
    do kper=1,nper  ! stress period loop
      if (kper==1) write(*,*)'Processing RCH package input...'
      write(iout,1)kper
      ! Read MF2005 input for current stress period
      call GWF2RCH7RP(this%IuOrig, igrid)
      !
      ! if kper = 1, write IRCH; if kper > 1, don't.
      writing = .false.
      if (kper==1) then
        writing = .true.
        write(iu,30)kper
        select case (NRCHOP)
        case (1)  ! recharge is only to layer 1
          constant = .true.
          nlval = 1
        case (2)  ! recharge is vertically distributed according to IRCH array
          ! determine if IRCH is constant
          write(iu,20)'IRCH'
          call ConstantInt2D(NCOL, NROW, IRCH, constant, nlval)
        case (3)  ! recharge is to highest active cell in each vertical column
          constant = .true.
          nlval = 1
        end select
        ! write IRCH
        if (NRCHOP==2) then
          if (constant) then
            write(iu,50)'CONSTANT', nlval
          else
            write(iu,50)'INTERNAL  FACTOR  1  IPRN ',iprni
            do i=1,nrow
              write(iu,80)(IRCH(j,i),j=1,ncol)
            enddo
          endif
        endif
      endif
      !
      call this%CurrentBlock%Clear(.true.)
      ! Write recharge rate for current stress period to a LineList
      call this%WriteStressPeriodArrayData(this%rlist, this%CurrentBlock, .true.)
!      call this%WriteStressPeriodData(this%CurrentBlock)
      ! determine if recharge array differs from previous stress period
      needToWrite = .not. same_lines(this%CurrentBlock, this%PreviousBlock)
      if (needToWrite) then
        call ConstantReal2D(NCOL, NROW, this%work, constant, rechval)
        ! If kper=1 and rech is zero, no need to write
        if (kper == 1 .and. constant .and. rechval == 0.0) needToWrite = .false.
        if (needToWrite) then
          if (.not. writing) then
            write(iu,30)kper
            writing = .true.
          endif
          write(iu,20)'RECHARGE'
          if (constant) then
            write(iu,60)'CONSTANT', rechval
          else
            ! write array control line and RECH array to input file
            write(iu,50)'INTERNAL  FACTOR  1.0  IPRN ',iprnr
            n = this%CurrentBlock%CountLines()
            do i=1,n
              call this%CurrentBlock%GetLine(i, line)
              write(iu,10)trim(line)
            enddo
          endif
        endif
      endif
      if (writing) write(iu,40)
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
    enddo  ! stress period loop
    !
    call this%BlockA%Clear(.true.)
    call this%BlockB%Clear(.true.)
    !
    return
  end subroutine ProcessStressLoop

  subroutine WriteOptions(this)
    ! Override if needed
    implicit none
    class(RchPackageWriterType) :: this
    ! local
    integer :: i, iu
    ! formats
    5  format()
    20 format(2x,a,2x,a)
    50 format(2x,a)
    60 format(a)
    !
    iu = this%fileobj%IUnit
    ! Write BEGIN Options
    write(iu,5)
    write(iu,60)'BEGIN Options'
    !
    write(iu,50)'READASARRAYS'
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

  function MyType(this) result (ctype)
    ! dummy 
    class(RchPackageWriterType) :: this
    character(len=LENCTYPE) :: ctype
    !
    ctype = 'RchPackageWriterType'
    !
    return
  end function MyType
  
end module RchPackageWriterModule
