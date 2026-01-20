module WelPackageWriterModule

  use ConstantsModule, only: MAXCHARLEN
  use ConstantsPHMFModule, only: FCINPUT, LENCTYPE
  use GLOBAL, only: NPER, IUNIT, NCOL, NROW, NLAY, IOUT
  use GWFFHBMODULE, only: IFLLOChfb => IFLLOC, NFLWhfb => NFLW
  USE GWFWELMODULE, ONLY: NWELLS, MXWELL, NWELVL, IPRWEL, NPWEL, &
                          IWELPB,NNPWEL,WELAUX,WELL, IWELCB
  use GwfWelNwt, only: GWF2WEL7AR, GWF2WEL7RP
  use LineListModule, only: same_lines
  use PackageWriterModule, only: PackageWriterType
  use SimPHMFModule, only: ustop, store_error, store_warning
  use utl7module, only: U1DREL, U2DREL,  &
                        urword, URDCOM,  &
                        ULSTRD

  type, extends(PackageWriterType) :: WelPackageWriterType
  contains
    procedure :: MyType
    procedure :: ProcessAllocate
    procedure :: ProcessStressLoop
    procedure :: WriteFileUsingTS
  end type WelPackageWriterType

contains

  subroutine ProcessAllocate(this, igrid)
    implicit none
    ! dummy
    class(WelPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    character(len=MAXCHARLEN) :: fname
    character(len=100) :: ctemp
    integer :: iunitnwt
    ! format
    !
    call this%AllocatePointers()
    this%Active = .true.
    this%DefaultBudgetText = 'WELLS'
    ! Ned todo: not sure that this will work for LGR--need to look at LGR source
    this%IuOrig = IUNIT(2)        ! Specific to WEL package
    this%fileobj%FCode = FCINPUT
    !
    this%fileobj%FType = 'WEL6'
    this%PkgType = 'WEL'
    fname = trim(this%ModelBasename) // '.wel'
    call this%FileWriterType%InitializeFile(fname, this%fileobj%FType, &
                                            this%PackageName)
    iunitnwt = 0
    call GWF2WEL7AR(this%IuOrig, IUNITNWT, igrid)
    !
    ! NStressDim is number of values after k,i,j required to specify boundary stress
    this%NStressDim = 1     ! Specific to WEL package; for Q
    !
    ! NVALP1 is rlist index for element following minimum required data
    this%nvalp1 = this%NStressDim + 4
    !
    ! Assign package-specific pointers
    this%rlist => WELL
    this%ipr => IPRWEL
    this%ICbc => IWELCB
    this%NBndPeriod => NWELLS
    this%Aux => WELAUX
    this%IPB => IWELPB
    this%MaxActiveBnd = this%IPB - 1
    this%NBNDVL => NWELVL
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
    class(WelPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: kper
    logical :: currentA
    ! formats
    10 format(/,'Processing WEL package for stress period ',i0)
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
      if (kper==1) write(*,*)'Processing WEL package input...'
      write(iout,10)kper
      ! Read MF2005 input for current stress period
      call GWF2WEL7RP(this%IuOrig, igrid)
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

  subroutine WriteFileUsingTS(this, igrid, ndim, tsnames)
    implicit none
    ! dummy
    class(WelPackageWriterType) :: this
    integer, intent(in) :: igrid
    integer, intent(in) :: ndim
    character(len=*), dimension(ndim) :: tsnames
    ! local
    integer :: i, iu, j, k, m
    ! formats
    1 format()
    5 format(a)
    10 format(2x,3(i0,2x),a)
    !
    call this%WriteOptions()
    call this%WriteDimensions()
    !
    ! Ned todo: write WEL info here (needed for FHB-based wells)
    if (ndim > 0) then
      iu = this%fileobj%IUnit
      write(iu,1)
      write(iu,5)'BEGIN PERIOD 1'
      do m=1, ndim
        k = IFLLOChfb(1,m)   ! layer
        i = IFLLOChfb(2,m)   ! row
        j = IFLLOChfb(3,m)   ! column
        write(iu,10)k, i, j, tsnames(m)
      enddo
      write(iu,5)'END PERIOD'
    endif
    call this%CloseFile()
    !
    return
  end subroutine WriteFileUsingTS

  function MyType(this) result (ctype)
    ! dummy 
    class(WelPackageWriterType) :: this
    character(len=LENCTYPE) :: ctype
    !
    ctype = 'WelPackageWriterType'
    !
    return
  end function MyType
  
end module WelPackageWriterModule
