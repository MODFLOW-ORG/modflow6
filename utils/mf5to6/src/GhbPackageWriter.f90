module GhbPackageWriterModule

  use ConstantsModule, only: MAXCHARLEN
  use ConstantsPHMFModule, only: FCINPUT, LENCTYPE
  use GLOBAL, only: NPER, IUNIT, NCOL, NROW, NLAY, IOUT
  use GWFGHBMODULE, only: NBOUND, MXBND, NGHBVL, IPRGHB, NPGHB, &
                          IGHBPB, NNPGHB, GHBAUX, BNDS, IGHBCB
  use GwfGhbSubs, only: GWF2GHB7AR, GWF2GHB7RP
  use LineListModule, only: same_lines
  use PackageWriterModule, only: PackageWriterType
  use SimPHMFModule, only: ustop, store_error, store_warning
  use utl7module, only: U1DREL, U2DREL,  &
                        urword, URDCOM,  &
                        ULSTRD

  type, extends(PackageWriterType) :: GhbPackageWriterType
  contains
    procedure :: MyType
    procedure :: ProcessAllocate
    procedure :: ProcessStressLoop
  end type GhbPackageWriterType

contains

  subroutine ProcessAllocate(this, igrid)
    implicit none
    ! dummy
    class(GhbPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    character(len=MAXCHARLEN) :: fname
    character(len=100) :: ctemp
    ! format
    !
    call this%AllocatePointers()
    this%Active = .true.
    this%DefaultBudgetText = 'HEAD DEP BOUNDS'
    ! Ned todo: not sure that this will work for LGR--need to look at LGR source
    this%IuOrig = IUNIT(7)        ! Specific to GHB package
    this%fileobj%FCode = FCINPUT
    !
    this%fileobj%FType = 'GHB6'
    this%PkgType = 'GHB'
    fname = trim(this%ModelBasename) // '.ghb'
    call this%FileWriterType%InitializeFile(fname, this%fileobj%FType, &
                                            this%PackageName)
    call GWF2GHB7AR(this%IuOrig, igrid)
    !
    ! NStressDim is number of values after k,i,j required to specify boundary stress
    this%NStressDim = 2     ! Specific to GHB package; for Bhead, Cond
    !
    ! NVALP1 is rlist index for element following minimum required data
    this%nvalp1 = this%NStressDim + 4
    !
    ! Assign package-specific pointers
    this%rlist => BNDS
    this%ipr => IPRGHB
    this%ICbc => IGHBCB
    this%NBndPeriod => NBOUND
    this%Aux => GHBAUX
    this%IPB => IGHBPB
    this%MaxActiveBnd = this%IPB - 1
    this%NBNDVL => NGHBVL
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
    class(GhbPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: kper
    logical :: currentA
    ! formats
    10 format(/,'Processing GHB package for stress period ',i0)
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
      if (kper==1) write(*,*)'Processing GHB package input...'
      write(iout,10)kper
      ! Read MF2005 input for current stress period
      call GWF2GHB7RP(this%IuOrig, igrid)
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
    class(GhbPackageWriterType) :: this
    character(len=LENCTYPE) :: ctype
    !
    ctype = 'GhbPackageWriterType'
    !
    return
  end function MyType
  
end module GhbPackageWriterModule
