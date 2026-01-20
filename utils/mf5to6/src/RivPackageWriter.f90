module RivPackageWriterModule

  use ConstantsModule, only: MAXCHARLEN
  use ConstantsPHMFModule, only: LENCTYPE
  use GLOBAL, only: NPER, IUNIT, NCOL, NROW, NLAY, IOUT
  use GWFRIVMODULE, only: NRIVER, MXRIVR, NRIVVL, IPRRIV, NPRIV, &
                           IRIVPB, NNPRIV, RIVAUX, RIVR, IRIVCB
  use GwfRivSubs, only: GWF2RIV7AR, GWF2RIV7RP
  use LineListModule, only: same_lines
  use PackageWriterModule, only: PackageWriterType
  use SimPHMFModule, only: ustop, store_error, store_warning
  use utl7module, only: U1DREL, U2DREL,  &
                        urword, URDCOM,  &
                        ULSTRD

  type, extends(PackageWriterType) :: RivPackageWriterType
  contains
    procedure :: MyType
    procedure :: ProcessAllocate
    procedure :: ProcessStressLoop
  end type RivPackageWriterType

contains

  subroutine ProcessAllocate(this, igrid)
    use ConstantsPHMFModule, only: FCINPUT
    implicit none
    ! dummy
    class(RivPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    character(len=MAXCHARLEN) :: fname
    character(len=100) :: ctemp
    ! format
    !
    this%PkgType = 'RIV'
    call this%AllocatePointers()
    this%Active = .true.
    this%DefaultBudgetText = 'RIVER'
    ! Ned todo: not sure that this will work for LGR--need to look at LGR source
    this%IuOrig = IUNIT(4)        ! Specific to RIV package
    this%fileobj%FCode = FCINPUT
    !
    this%fileobj%FType = 'RIV6'
    fname = trim(this%ModelBasename) // '.riv6'
    call this%FileWriterType%InitializeFile(fname, this%fileobj%FType, &
                                            this%PackageName)
    call GWF2RIV7AR(this%IuOrig, igrid)
    !
    ! NStressDim is number of values after k,i,j required to specify boundary stress
    this%NStressDim = 3     ! Specific to RIV package; for Stage, Cond, Rbot
    !
    ! NVALP1 is rlist index for element following minimum required data
    this%nvalp1 = this%NStressDim + 4
    !
    ! Assign package-specific pointers
    this%rlist => RIVR
    this%ipr => IPRRIV
    this%ICbc => IRIVCB
    this%NBndPeriod => NRIVER
    this%Aux => RIVAUX
    this%IPB => IRIVPB
    this%MaxActiveBnd = this%IPB - 1
    this%NBNDVL => NRIVVL
    this%NAux = this%NBNDVL - this%nvalp1
    this%nstop = this%nvalp1 - 1 + this%NAux
    !
    ! define format for printing boundary data
    write(ctemp,'(i0)')this%NAux + this%NStressDim
    this%fmat = '(3(2x,i0),' // trim(ctemp) // '(2x,g15.8))'
    !
    return
  end subroutine ProcessAllocate

  subroutine ProcessStressLoop(this, igrid)
    implicit none
    ! dummy
    class(RivPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: kper
    logical :: currentA
!    integer :: kc, kp
    ! formats
    10 format(/,'Processing RIV package for stress period ',i0)
    !
    ! Write Options and Dimensions blocks
    call this%WriteOptions()
    call this%WriteDimensions()
    !
    ! Initially, current block is BlockA; alternate each stress period
    this%CurrentBlock => this%BlockA
    this%PreviousBlock => this%BlockB
    currentA = .true.
    !
    do kper=1,nper
      if (kper==1) write(*,*)'Processing RIV package input...'
      write(iout,10)kper
      ! Read MF2005 input for current stress period
      call GWF2RIV7RP(this%IuOrig, igrid)
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

  function MyType(this) result (ctype)
    ! dummy 
    class(RivPackageWriterType) :: this
    character(len=LENCTYPE) :: ctype
    !
    ctype = 'RivPackageWriterType'
    !
    return
  end function MyType
  
end module RivPackageWriterModule
