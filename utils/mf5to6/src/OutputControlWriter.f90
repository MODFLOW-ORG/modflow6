module OutputControlWriterModule

  use ConstantsModule, only: LINELENGTH, MAXCHARLEN
  use ConstantsPHMFModule, only: FCINPUT, LENCTYPE
  use FileListModule, only: FileListType
  use FileTypeModule, only: FileType
  use GLOBAL, only: IOUT, IUNIT, NPER, NSTP, NLAY, cbcfilename
  use GWFBASMODULE, only: IHDDFL, IBUDFL, ICBCFL, IPEROC, ITSOC, IBDOPT, IOFLG
  use GwfBasOcSubsModule, only: GWF2BAS7OC
  use LineListModule, only: LineListType
  use PackageWriterModule, only: PackageWriterType
  use SimPHMFModule, only: store_error, store_note, ustop

  integer, parameter :: UNASSIGNED = -99
  integer, parameter :: BUDGET   = 1
  integer, parameter :: HEAD     = 2
  integer, parameter :: UNUSED = 3
  integer, parameter :: IBOUND   = 4

  type OcSaverType
    character(len=8) :: DType = ''
    logical :: PrtAllSteps = .false.
    logical :: SavAllSteps = .false.
    logical :: PrtLast = .false.
    logical :: SavLast = .false.
    logical, allocatable, dimension(:) :: Sav
    logical, allocatable, dimension(:) :: Prt
  contains
  end type OcSaverType

  type, extends(PackageWriterType) :: OutputControlWriterType
    ! -- Budget options
    character(len=LINELENGTH) :: BdSvFil = ''
    ! -- Head options
    character(len=LINELENGTH) :: HdSvFil = ''
    character(len=LINELENGTH) :: HdSvFmt = ''
    integer                   :: HdPrnNCol = 0
    integer                   :: HdPrnWid  = 0
    integer                   :: HdPrnNDig = 0
    logical                   :: HdLbl = .false.
    logical                   :: HdWrp = .false.
    character(len=11)         :: HdEfgs = ''
    integer                   :: IHEDFM = UNASSIGNED
    integer                   :: IHEDUN = UNASSIGNED
    ! -- Drawdown options
    character(len=LINELENGTH) :: DdSvFil = ''
    character(len=LINELENGTH) :: DdSvFmt = ''
    integer                   :: DdPrnNCol = 0
    integer                   :: DdPrnWid  = 0
    integer                   :: DdPrnNDig = 0
    logical                   :: DdLbl = .false.
    logical                   :: DdWrp = .false.
    character(len=11)         :: DdEfgs = ''
    integer                   :: IDDNFM = UNASSIGNED
    integer                   :: IDDNUN = UNASSIGNED
    ! -- Ibound options
    character(len=LINELENGTH) :: IbSvFil = ''
    character(len=LINELENGTH) :: IbSvFmt = '(20I4)'
    integer                   :: IbPrnNCol = 0
    integer                   :: IbPrnWid  = 0
    integer                   :: IbPrnNDig = 0
    logical                   :: IbLbl = .false.
    integer                   :: IBOUUN = 0
    ! -- Stress period data
    type(OcSaverType), allocatable, dimension(:) :: Savers
  contains
!    ! Public procedures
    procedure, public :: MyType
    procedure, public :: ProcessAllocate
    procedure, public :: ProcessStressLoop => process_stress_loop
    procedure, public :: WriteFile
    procedure, public :: WriteStressPeriodListData
    ! Private procedures
    procedure, private :: set_saver_defaults
    procedure, private :: process_options
    procedure, private :: write_options
  end type OutputControlWriterType

contains

  subroutine set_saver_defaults(this, nsteps)
    implicit none
    ! dummy
    class(OutputControlWriterType) :: this
    integer, intent(in) :: nsteps
    ! local
    integer :: i
    !
    do i=1,4
      this%Savers(i)%PrtAllSteps = .false.
      this%Savers(i)%SavAllSteps = .false.
      if (allocated(this%Savers(i)%Prt)) then
        deallocate(this%Savers(i)%Prt)
      endif
      if (allocated(this%Savers(i)%Sav)) then
        deallocate(this%Savers(i)%Sav)
      endif
      allocate(this%Savers(i)%Prt(nsteps))
      allocate(this%Savers(i)%Sav(nsteps))
      this%Savers(i)%Prt = .false.
      this%Savers(i)%Sav = .false.
      select case (this%Savers(i)%DType)
      case ('BUDGET', 'HEAD')
        this%Savers(i)%PrtLast = .true.
        this%Savers(i)%Prt(nsteps) = .true.
      case default
        this%Savers(i)%PrtLast = .false.
      end select
    enddo
    !
    return
  end subroutine set_saver_defaults

  subroutine ProcessAllocate(this, igrid)
    ! Implements deferred procedure of PackageWriterType
    implicit none
    ! dummy
    class(OutputControlWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    !
    if (.not. allocated(this%Savers)) then
      allocate(this%Savers(4))
      this%Savers(BUDGET)%DType = 'BUDGET'
      this%Savers(HEAD)%DType = 'HEAD'
!      this%Savers(DRAWDOWN)%DType = 'DRAWDOWN'
      this%Savers(IBOUND)%DType = 'IBOUND'
      !
      call this%AllocatePointers()
    endif
    this%Active = .true.
    !
    this%IuOrig = IUNIT(12)        ! Specific to Output Control
    this%fileobj%FCode = FCINPUT
    !
    this%fileobj%FType = 'OC6'
    !
    return
  end subroutine ProcessAllocate

  subroutine WriteFile(this, Mf6FileList, igrid)
    implicit none
    ! dummy
    class(OutputControlWriterType) :: this
    type(FileListType) :: Mf6FileList
    integer, intent(in) :: igrid
    !
    if (this%Active) then
      write(*,*)'Processing OC package input...'
      ! Process some of the OC options
      call this%process_options(Mf6FileList)
      ! Write the Options block
      call this%write_options()
      ! Generate stress period input
      call this%ProcessStressLoop(igrid)
    else
      ! Ned todo: There probably should be something in an OC file
      ! even if the original model did not have one.
    endif
    !
    return
  end subroutine WriteFile

  subroutine process_options(this, Mf6FileList)
    implicit none
    ! dummy
    class(OutputControlWriterType) :: this
    type(FileListType) :: Mf6FileList
    ! local
    type(FileType), pointer :: filtyp
    character(len=LINELENGTH) :: msg
    ! format
    10 format('In processing output control, no file found', &
              ' corresponding to unit number: ',i0)
    !
    ! Decode IHEDFM
    call decode_ifm(this%IHEDFM, this%HdPrnNCol, this%HdPrnWid, &
                    this%HdPrnNDig, this%HdWrp, this%HdEfgs)
    !
    ! Decode IDDNFM
    call decode_ifm(this%IDDNFM, this%DdPrnNCol, this%DdPrnWid, &
                    this%DdPrnNDig, this%DdWrp, this%DdEfgs)
    !
    ! Find file associated with IHEDUN
    if (this%IHEDUN > 0) then
      filtyp => Mf6FileList%GetFileByUnit(this%IHEDUN)
      if (associated(filtyp)) then
        this%HdSvFil = filtyp%FName
      else
        write(msg,10)this%IHEDUN
        call store_error(msg)
        call ustop()
      endif
    endif
    !
    ! Find file associated with IDDNUN
    if (this%IDDNUN > 0) then
      msg = 'Calculation of drawdown is not supported in MODFLOW 6.'
      call store_note(msg)
!      filtyp => Mf6FileList%GetFileByUnit(this%IDDNUN)
!      if (associated(filtyp)) then
!        this%DdSvFil = filtyp%FName
!      else
!        write(msg,10)this%IDDNUN
!        call store_error(msg)
!        call ustop()
!      endif
    endif
    !
    ! Find CBC file - identified by I**CB option in one or more packages.
    ! E.g. ILPFCB for LPF package.
    this%BdSvFil = cbcfilename(1:LINELENGTH)
    !
!    ! Ensure that HdSvFil and DdSvFil are not the same file.
!    if (this%HdSvFil /= '') then
!      if (this%HdSvFil == this%DdSvFil) then
!        this%DdSvFil = 'ddn_' // this%DdSvFil
!      endif
!    endif
    !
    return
  end subroutine process_options

  subroutine write_options(this)
    implicit none
    ! dummy
    class(OutputControlWriterType) :: this
    ! local
    integer :: iout
    character(len=20) :: options
    ! formats
    1 format()
    10 format(a)
    20 format(2x,a,2x,a)
    30 format(2x,a,2x,a,2x,a)
    40 format(2x,a,1x,i0,2x,a,1x,i0,2x,a,1x,i0,2x,a)
    !
    iout = this%fileobj%IUnit
    write(iout,1)
    write(iout,10)'BEGIN Options'
    !
    ! Budget save file
    if (this%BdSvFil /= '') then
      write(iout,20)'BUDGET FILEOUT', trim(this%BdSvFil)
    endif
    !
    ! Head save file
    if (this%HdSvFil /= '') then
      write(iout,20)'HEAD FILEOUT', trim(this%HdSvFil)
    endif
    !
    ! Head save format
    if (this%HdSvFmt /= '') then
!      if (this%HdLbl) then
!        write(iout,30)'HEAD SAVE FORMAT',this%HdSvFmt,'LABEL'
!      else
!        write(iout,20)'HEAD SAVE FORMAT',this%HdSvFmt
!      endif
    endif
    !
    ! Head print format
    if (this%IHEDFM /= UNASSIGNED) then
      options = this%HdEfgs
!      if (this%HdWrp) then
!        options = trim(options) // ' WRAP'
!      else
!        options = trim(options) // ' STRIP'
!      endif
      write(iout,40)'HEAD PRINT_FORMAT  COLUMNS',this%HdPrnNCol, 'WIDTH', &
                    this%HdPrnWid, 'DIGITS', this%HdPrnNDig, options
    endif
    !
!    ! Drawdown save file
!    if (this%DdSvFil /= '') then
!      write(iout,20)'DRAWDOWN SAVE FILE', trim(this%DdSvFil)
!    endif
!    !
!    ! Drawdown save format
!    if (this%ddSvFmt /= '') then
!      if (this%DdLbl) then
!        write(iout,30)'DRAWDOWN SAVE FORMAT',this%DdSvFmt,'LABEL'
!      else
!        write(iout,20)'DRAWDOWN SAVE FORMAT',this%DdSvFmt
!      endif
!    endif
!    !
!    ! Drawdown print format
!    if (this%IDDNFM /= UNASSIGNED) then
!      options = this%DdEfgs
!      if (this%DdWrp) then
!        options = trim(options) // ' WRAP'
!      else
!        options = trim(options) // ' STRIP'
!      endif
!      write(iout,40)'DRAWDOWN PRINT FORMAT  COLUMNS',this%DdPrnNCol, 'WIDTH', &
!                    this%DdPrnWid, 'DIGITS', this%DdPrnNDig, options
!    endif
    !
    write(iout,10)'END Options'
    !
    return
  end subroutine write_options

  subroutine process_stress_loop(this, igrid)
    ! Implements deferred procedure of PackageWriterType
    implicit none
    ! dummy
    class(OutputControlWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: i, icnvg, inoc, k, kper, kstp, nstps
    integer :: kstpHS, kstpHP, kstpDS, kstpDP, kstpBS, kstpBP
    logical :: currentA, UseTSDefaults
    logical :: HdPrt, HdSav, DdPrt, DdSav, BdPrt, BdSav
    logical :: HdPrtStp, HdSavStp, DdPrtStp, DdSavStp
    ! formats
    10 format(/,'Processing OC package for stress period ',i0)
    !
    inoc = this%fileobj%IUnit
    icnvg = 1  ! behave as if solution converged
    !
    ! Initially, current block is BlockA; alternate each stress period
    this%CurrentBlock => this%BlockA
    this%PreviousBlock => this%BlockB
    currentA = .true.
    !
    periodloop: do kper=1,nper
      write(iout,10)kper
      nstps = NSTP(kper)
      call this%set_saver_defaults(nstps)
      kstpHS = 0
      kstpHP = 0
      kstpDS = 0
      kstpDP = 0
      kstpBS = 0
      kstpBP = 0
      HdPrt = .false.
      HdSav = .false.
      DdPrt = .false.
      DdSav = .false.
      BdPrt = .false.
      BdSav = .false.
      steploop: do kstp = 1,nstps
        ! -- Determine which output is specified for this time step
        !    This assigns IHDDFL, IBUDFL, ICBCFL, IOFLG array
        call GWF2BAS7OC(KSTP, KPER, ICNVG, this%IuOrig, IGRID, UseTSDefaults)
        if (.not. UseTSDefaults) then
          this%Savers(BUDGET)%PrtLast = .false.
          this%Savers(HEAD)%PrtLast = .false.
          this%Savers(UNUSED)%PrtLast = .false.
          ! Cycle through IOFLG by layer; set save (or print) if any layer
          ! head (or drawdown) is to be saved (or printed) this time step.
          HdPrtStp = .false.
          HdSavStp = .false.
          DdPrtStp = .false.
          DdSavStp = .false.
          dataloop: do i=1,4
            layerloop: do k=1,nlay
              if (IOFLG(k,i) /= 0) then
                select case (i)
                case (1) ! Head print
                  if (.not. HdPrt) then
                    ! If only the last time step is to be printed,
                    ! set LAST option.
                    if (kstp == nstps) then
                      this%Savers(HEAD)%PrtLast = .true.
                    endif
                  endif
                  HdPrt = .true.
                  HdPrtStp = .true.
                case (2) ! Drawdown print
!                  if (.not. DdPrt) then
!                    ! If only the last time step is to be printed,
!                    ! set LAST option.
!                    if (kstp == nstps) then
!                      this%Savers(DRAWDOWN)%PrtLast = .true.
!                    endif
!                  endif
!                  DdPrt = .true.
!                  DdPrtStp = .true.
                case (3) ! Head save
                  if (.not. HdSav) then
                    ! If only the last time step is to be saved,
                    ! set LAST option.
                    if (kstp == nstps) then
                      this%Savers(HEAD)%SavLast = .true.
                    endif
                  endif
                  HdSav = .true.
                  HdSavStp = .true.
                case (4) ! Drawdown save
!                  if (.not. DdSav) then
!                    ! If only the last time step is to be saved,
!                    ! set LAST option.
!                    if (kstp == nstps) then
!                      this%Savers(DRAWDOWN)%SavLast = .true.
!                    endif
!                  endif
!                  DdSav = .true.
!                  DdSavStp = .true.
                end select
              endif
            enddo layerloop
          enddo dataloop
          ! Accumulate time steps for each type of output
          if (HdPrtStp) then
            kstpHP = kstpHP + 1
          endif
          if (HdSavStp) then
            kstpHS = kstpHS + 1
          endif
          if (DdPrtStp) then
            kstpDP = kstpDP + 1
          endif
          if (DdSavStp) then
            kstpDS = kstpDS + 1
          endif
          !
          ! Head and drawdown control
          this%Savers(HEAD)%Prt(kstp) = HdPrtStp
          this%Savers(HEAD)%Sav(kstp) = HdSavStp
!          this%Savers(DRAWDOWN)%Prt(kstp) = DdPrtStp
!          this%Savers(DRAWDOWN)%Sav(kstp) = DdSavStp
          !
          ! Budget control - printing
          this%Savers(BUDGET)%Prt(kstp) = (IBUDFL /= 0)
          if (IBUDFL /= 0) then
            kstpBP = kstpBP + 1
            if (.not. BdPrt) then
              if (kstp == nstps) then
                this%Savers(BUDGET)%PrtLast = .true.
              endif
              BdPrt = .true.
            endif
          endif
          !
          ! Budget control - saving
          if (cbcfilename /= '') then
            this%Savers(BUDGET)%Sav(kstp) = (ICBCFL /= 0)
            if  (ICBCFL /= 0) then
              kstpBS = kstpBS + 1
              if (.not. BdSav) then
                if (kstp == nstps) then
                  this%Savers(BUDGET)%SavLast = .true.
                endif
                BdSav = .true.
              endif
            endif
          endif
        endif
      enddo steploop
      !
      ! If all steps are to be printed/saved, set ALL option
      this%Savers(HEAD)%PrtAllSteps = (kstpHP == nstps)
      this%Savers(HEAD)%SavAllSteps = (kstpHS == nstps)
!      this%Savers(DRAWDOWN)%PrtAllSteps = (kstpDP == nstps)
!      this%Savers(DRAWDOWN)%SavAllSteps = (kstpDS == nstps)
      this%Savers(BUDGET)%PrtAllSteps = (kstpBP == nstps)
      this%Savers(BUDGET)%SavAllSteps = (kstpBS == nstps)
      !
      ! If print/save ALL is set, unset print/save LAST
      do i=BUDGET,HEAD
        if (this%Savers(i)%PrtAllSteps) then
          this%Savers(i)%PrtLast = .false.
        endif
        if (this%Savers(i)%SavAllSteps) then
          this%Savers(i)%SavLast = .false.
        endif
      enddo
      !
      ! Write block contents for current stress period
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
    enddo periodloop
    !
    call this%BlockA%Clear(.true.)
    call this%BlockB%Clear(.true.)
    !
    return
  end subroutine process_stress_loop

  subroutine WriteStressPeriodListData(this, lineList)
    ! Write stress period data for MF6 Output Control
    implicit none
    ! dummy
    class(OutputControlWriterType) :: this
    type(LineListType), pointer :: lineList
    ! local
    integer :: i, j, nsteps, nstpprt, nstpsav
    character(len=MAXCHARLEN) :: line
    character(len=10) :: stepnum
    ! formats
10  format(2x,a)
20  format(2x,a,1x,a)
    !
    ! Initialize variables
    nsteps = size(this%Savers(1)%Prt)
    !
    ! Process mf2005 defaults
    do i=BUDGET,IBOUND
      line = ''
      if (this%Savers(i)%PrtLast) then
        write(line, 10) 'PRINT ' // trim(this%Savers(i)%DType) // ' LAST'
      endif
      if (line/='') then
        call lineList%AddLine(line)
      endif
      line = ''
      if (this%Savers(i)%SavLast) then
        write(line, 10) 'SAVE ' // trim(this%Savers(i)%DType) // ' LAST'
      endif
      if (line/='') then
        call lineList%AddLine(line)
      endif
    enddo
    !
    ! Process step-by-step control
    do i=BUDGET,IBOUND
      ! Process SAVE control
      if (.not. this%Savers(i)%SavLast) then
        line = ''
        ! Count steps to save
        nstpsav = 0
        do j=1,nsteps
          if (this%Savers(i)%Sav(j)) nstpsav = nstpsav + 1
        enddo
        ! If all steps are to be saved, just specify that--don't list all steps.
        if (this%Savers(i)%SavAllSteps) then
          nstpsav = 0
          line = '  SAVE ' // trim(this%Savers(i)%DType) // ' ALL'
          call lineList%AddLine(line)
        endif
        ! If only some steps are to be saved, specify that.
        if (nstpsav>0) then
          write(line,20)'  SAVE ' // trim(this%Savers(i)%DType), ' STEPS'
          do j=1,nsteps
            if (this%Savers(i)%Sav(j)) then
              write(stepnum,'(i0)')j
              line = trim(line) // ' ' // trim(stepnum)
            endif
          enddo
          call lineList%AddLine(line)
        endif
      endif
      !
      ! Process PRINT control
      if (.not. this%Savers(i)%PrtLast) then
        line = ''
        ! Count steps to print
        nstpprt = 0
        do j=1,nsteps
          if (this%Savers(i)%Prt(j)) nstpprt = nstpprt + 1
        enddo
        ! If all steps are to be printed, just specify that--don't list all steps.
        if (this%Savers(i)%PrtAllSteps) then
          nstpprt = 0
          line = '  PRINT ' // trim(this%Savers(i)%DType) // ' ALL'
          call lineList%AddLine(line)
        endif
        ! If only some steps are to be save, specify that.
        if (nstpprt>0) then
          write(line,20)'  PRINT ' // trim(this%Savers(i)%DType), ' STEPS'
          do j=1,nsteps
            if (this%Savers(i)%Prt(j)) then
              write(stepnum,'(i0)')j
              line = trim(line) // ' ' // trim(stepnum)
            endif
          enddo
          call lineList%AddLine(line)
        endif
      endif
    enddo
    !
    return
  end subroutine WriteStressPeriodListData

  subroutine decode_ifm(ifm, ncol, nwid, ndig, wrap, efgs)
    implicit none
    ! dummy
    integer,          intent(in)  :: ifm
    integer,          intent(out) :: ncol, nwid, ndig
    logical,          intent(out) :: wrap
    character(len=*), intent(out) :: efgs
    ! local
    integer :: ifmabs
    character(len=200) :: msg
    !
    wrap = (ifm >= 0)
    if (.not. wrap) then
      msg = 'STRIP format for printing is not supported in MODFLOW 6'
      call store_note(msg)
    endif
    ifmabs = abs(ifm)
    select case (ifmabs)
    case (0)
      ncol = 10
      nwid = 11
      ndig = 4
      efgs = 'GENERAL'
    case (1)
      ncol = 11
      nwid = 10
      ndig = 3
      efgs = 'GENERAL'
    case (2)
      ncol = 9
      nwid = 13
      ndig = 6
      efgs = 'GENERAL'
    case (3)
      ncol = 15
      nwid = 7
      ndig = 1
      efgs = 'FIXED'
    case (4)
      ncol = 15
      nwid = 7
      ndig = 2
      efgs = 'FIXED'
    case (5)
      ncol = 15
      nwid = 7
      ndig = 3
      efgs = 'FIXED'
    case (6)
      ncol = 15
      nwid = 7
      ndig = 4
      efgs = 'FIXED'
    case (7)
      ncol = 20
      nwid = 5
      ndig = 0
      efgs = 'FIXED'
    case (8)
      ncol = 20
      nwid = 5
      ndig = 1
      efgs = 'FIXED'
    case (9)
      ncol = 20
      nwid = 5
      ndig = 2
      efgs = 'FIXED'
    case (10)
      ncol = 20
      nwid = 5
      ndig = 3
      efgs = 'FIXED'
    case (11)
      ncol = 20
      nwid = 5
      ndig = 4
      efgs = 'FIXED'
    case (12)
      ncol = 10
      nwid = 11
      ndig = 4
      efgs = 'GENERAL'
    case (13)
      ncol = 10
      nwid = 6
      ndig = 0
      efgs = 'FIXED'
    case (14)
      ncol = 10
      nwid = 6
      ndig = 1
      efgs = 'FIXED'
    case (15)
      ncol = 10
      nwid = 6
      ndig = 2
      efgs = 'FIXED'
    case (16)
      ncol = 10
      nwid = 6
      ndig = 3
      efgs = 'FIXED'
    case (17)
      ncol = 10
      nwid = 6
      ndig = 4
      efgs = 'FIXED'
    case (18)
      ncol = 10
      nwid = 6
      ndig = 5
      efgs = 'FIXED'
    case (19)
      ncol = 5
      nwid = 12
      ndig = 5
      efgs = 'GENERAL'
    case (20)
      ncol = 6
      nwid = 11
      ndig = 4
      efgs = 'GENERAL'
    case (21)
      ncol = 7
      nwid = 9
      ndig = 2
      efgs = 'GENERAL'
    case default
      ncol = 10
      nwid = 11
      ndig = 4
      efgs = 'GENERAL'
    end select
    !
    return
  end subroutine decode_ifm

  function MyType(this) result (ctype)
    ! dummy 
    class(OutputControlWriterType) :: this
    character(len=LENCTYPE) :: ctype
    !
    ctype = 'OutputControlWriterType'
    !
    return
  end function MyType
  
end module OutputControlWriterModule
