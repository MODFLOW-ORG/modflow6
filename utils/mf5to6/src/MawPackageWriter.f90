module MawPackageWriterModule

  use ConstantsModule, only: DZERO, MAXCHARLEN
  use ConstantsPHMFModule, only: BIGHKVALUE, FCINPUT, LENCTYPE
  use GLOBAL, only: NPER, IUNIT, NCOL, NROW, NLAY, IOUT, &
                    BOTM, LBOTM
  use GlobalVariablesModule, only: echo
  use GWFMNW2MODULE, only: NMNW2, MNWMAX, NMNWVL, IWL2CB, MNWPRNT, &
                           NODTOT, INTTOT, NTOTNOD, WELLID, MNWAUX, &
                           MNW2, MNWNOD, MNWINT, CapTable
  use GwfMnwSubs, only: GWF2MNW27AR, GWF2MNW27RP
  use LineListModule, only: same_lines, LineListType
  use PackageWriterModule, only: PackageWriterType
  use SimPHMFModule, only: count_errors, store_error, store_warning, ustop
  use utl7module, only: U1DREL, U2DREL, &
                        urword, URDCOM,  &
                        ULSTRD

  private
  public :: MawPackageWriterType

  type, extends(PackageWriterType) :: MawPackageWriterType
    integer :: igwtunit = 0  ! To add support for GWT, assign this as GWT input unit number
  contains
    procedure :: MyType
    procedure :: ProcessAllocate
    procedure :: ProcessStressLoop
    procedure :: WriteOptions
    procedure :: WriteDimensions
    procedure :: WriteStressPeriodListData
    procedure :: WriteWellsBlock
    procedure :: WriteWellConnections
    procedure :: WriteInitialAttributes
  end type MawPackageWriterType

  character(len=10), dimension(0:4) :: LossType
  data LossType/'NONE      ', 'THIEM     ', 'SKIN      ', &
                   'GENERAL   ', 'SPECIFYCWC'/

contains

  subroutine ProcessAllocate(this, igrid)
    implicit none
    ! dummy
    class(MawPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: kperdummy = 1
    character(len=MAXCHARLEN) :: fname
    character(len=100) :: ctemp
    ! formats
    10 format(/,'Processing MAW package for stress period ',i0)
    12 format(a,'_',i0)
    !
    call this%AllocatePointers()
    this%Active = .true.
    this%DefaultBudgetText = 'MULTI-AQ WELLS'
    ! Ned todo: not sure that this will work for LGR--need to look at LGR source
    this%IuOrig = IUNIT(50)        ! Specific to MNW2 package
    this%fileobj%FCode = FCINPUT
    !
    this%fileobj%FType = 'MAW6'
    this%PkgType = 'MAW'
    write(this%PackageName,12)trim(this%PkgType),igrid
    fname = trim(this%ModelBasename) // '.maw'
    call this%FileWriterType%InitializeFile(fname, this%fileobj%FType, &
                                            this%PackageName)
    call GWF2MNW27AR(this%IuOrig, igrid)
    write(iout,10)kperdummy
    call GWF2MNW27RP(this%IuOrig,kperdummy,this%igwtunit,igrid)
    !
    ! NStressDim is number of values after k,i,j required to specify boundary stress
    this%NStressDim = 30                ! Specific to MNW2 package
    !
    ! NVALP1 is rlist index for element following minimum required data
    this%nvalp1 = this%NStressDim  + 4  ! Just a guess for MAW
    !
    ! Assign package-specific pointers
    this%rlist => null()     ! MNW2 is double-precision
    this%ipr => MNWPRNT
    this%ICbc => IWL2CB
    this%NBndPeriod => NMNW2
    this%Aux => MNWAUX
    this%IPB => MNWMAX
    this%MaxActiveBnd = this%IPB
    this%NBNDVL => NMNWVL
    this%NAux = this%NBNDVL - this%nvalp1
    this%nstop = this%nvalp1 - 1 + this%NAux
    !
    ! define format for printing boundary data
    write(ctemp,'(i0)')this%NAux + this%NStressDim
    this%fmat = '(3(2x,i0),' // trim(ctemp) // '(2x,g16.9))'
    !
    ! Write Options and Dimensions blocks
    call this%WriteOptions()
    call this%WriteDimensions()
    !
    ! Write Wells block
    call this%WriteWellsBlock()
    !
    ! Write Well_Connections block
    call this%WriteWellConnections()
    !
    !! Write Well_Initial_Attributes block
    !call this%WriteInitialAttributes()
    !
    return
  end subroutine ProcessAllocate

  subroutine ProcessStressLoop(this, igrid)
    implicit none
    ! dummy
    class(MawPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: kper
    logical :: currentA
!    integer :: kc, kp
    ! formats
    10 format(/,'Processing MAW package for stress period ',i0)
    !
    ! Initially, current block is BlockA; alternate each stress period
    this%CurrentBlock => this%BlockA
    this%PreviousBlock => this%BlockB
    currentA = .true.
    !
    do kper=1,nper
      if (kper==1) write(*,*)'Processing MNW package input...'
      ! Read MF2005 input for current stress period
      if (kper > 1) then
        write(iout,10)kper
        call GWF2MNW27RP(this%IuOrig, kper, this%igwtunit, igrid)
      endif
      ! Write MF6 input for current stress period to a LineList
      call this%CurrentBlock%Clear(.true.)
      call this%WriteStressPeriodListData(this%CurrentBlock)
      ! Write block to MF6 input file if needed
      call this%WriteBlockIfNeeded(kper, 'STEADY-STATE')
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

  subroutine WriteOptions(this)
    implicit none
    class(MawPackageWriterType) :: this
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
    ! Aux variable names
    if (this%NAux > 0) then
      do i=1,this%NAux
        write(iu,20) 'AUXILIARY', trim(this%Aux(i))
      enddo
    endif
    !
    if (this%ipr >= 1) then
      if (echo) then
        write(iu,50)'PRINT_INPUT'
      endif
      if (this%ipr >= 2) then
        write(iu,50)'PRINT_FLOWS'
      endif
    endif
    !
    if (this%ICbc > 0) then
      write(iu,50)'SAVE_FLOWS'
    elseif (this%ICbc < 0 .and. this%ipr < 2) then
      write(iu,50)'PRINT_FLOWS'
    endif
    !!
    !if (.not. this%Newton) then
    !  write(iu,50)'NO_NEWTON'
    !endif
    write(iu,50)'BOUNDNAMES'
    !
    ! Write END Options
    write(iu,60)'END Options'
    !
    return
  end subroutine WriteOptions

  subroutine WriteDimensions(this)
    implicit none
    class(MawPackageWriterType) :: this
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
    write(iu,30)'NMAWWELLS',this%MaxActiveBnd
    !
    ! Write END Dimensions
    write(iu,60)'END Dimensions'
    !
    return
  end subroutine WriteDimensions

  subroutine WriteStressPeriodListData(this, lineList)
    implicit none
    ! dummy
    class(MawPackageWriterType) :: this
    type(LineListType), pointer :: lineList
    ! local
    integer :: i
    double precision :: Qdes
    character(len=100) :: line
    ! formats
    50 format(2x,i0,2x,'RATE',2x,g16.9)
    !
    do i = 1,this%MaxActiveBnd
      Qdes = MNW2(5,i)
      write(line,50)i,Qdes
      call lineList%AddLine(line)
    enddo
    !
    return
  end subroutine WriteStressPeriodListData

  subroutine WriteWellsBlock(this)
    ! For each well, write:
    ! WELL wellno radius strt bottom condeqn ngwfnodes [name]
    implicit none
    ! dummy
    class(MawPackageWriterType) :: this
    ! local
    integer :: i, iloss, iu, nnodes, nnsigned, firstnode, lastnode
    integer :: firstint, lastint
    integer :: kbot, ktop, nodebot, nodetop
    integer :: il, ir, ic, n
    double precision :: radius, strt, bottom
    character(len=10) :: condeqn
    character(len=MAXCHARLEN) :: cwcmsg, genmsg, nonemsg
    ! formats
    1 format()
    10 format(a)
    20 format(2x,i0,3(2x,g15.8),2x,a,2x,i0,2x,a)
    !
    cwcmsg = 'MNW2 LossType "NONE" is not supported. In its place, ' // &
            'the "SPECIFIED" option of MAW8 is used with a large ' // &
            'cell-to-well conductance for any well with LossType ' // &
            '"NONE". If the solution does not ' // &
            'converge, the conductance may be too large, or another ' // &
            'option may need to be used.'
    genmsg = 'MNW2 LossType "GENERAL" has no equivalent in MAW8. ' // &
             'Any MNW2 well with LossType "GENERAL" will need to to be ' // &
             'added by the user to the MAW8 input file with a CONDEQN ' // &
             'option to be determined by the user.'
    nonemsg = 'MNW2 LossType "NONE" has no equivalent in MAW8. ' // &
             'Any MNW2 well with LossType "NONE" will need to to be ' // &
             'added by the user to the MAW8 input file with a CONDEQN ' // &
             'option to be determined by the user.'
    iu = this%fileobj%IUnit
    write(iu,1)
    write(iu,10) 'BEGIN PACKAGEDATA'
    !
    do i=1,this%MaxActiveBnd
      nnsigned = nint(MNW2(2,i))
      nnodes = abs(nnsigned)         ! Number of nodes
      firstnode = nint(MNW2(4,i))           ! Location of first node in MNWNOD array
      lastnode = firstnode + nnodes - 1     ! Location of last node in MNWNOD array
      firstint = nint(MNWNOD(12,firstnode)) ! Location of first interval in MNWINT array
      lastint = nint(MNWNOD(13,lastnode))   ! Location of last interval in MNWINT array
      if (nnsigned > 0) then
        radius = MNWNOD(5,firstnode) ! Rw for first node only
        ! Ztop and Zbotm are not read when nnsigned > 0. Instead, use TOP of
        ! uppermost node and BOTM of lowermost node.
        ! Order of entry of nodes is not prescribed, so need to find uppermost
        ! node and use its TOP, and find lowermost node and use its BOTM.
        ktop = 999999
        kbot = 0
        nodetop = 999999  ! May be able to eliminate
        nodebot = 0       ! May be able to eliminate
        do n=firstnode,lastnode
          il = nint(mnwnod(1,n))
          ir = nint(mnwnod(2,n))
          ic = nint(mnwnod(3,n))
          if (il < ktop) then
            ktop = il
            nodetop = n
          endif
          if (il > kbot) then
            kbot = il
            nodebot = n
          endif
        enddo
        strt = BOTM(ic,ir,LBOTM(ktop)-1)
        bottom = BOTM(ic,ir,LBOTM(kbot))
      else
        radius = MNWINT(5,firstint) ! Rw for first interval
        ! MNW2 does not read starting head. Use Ztop (mnwint(1,*) instead.
        strt = mnwint(1,firstint) ! First interval is always the highest
        bottom = mnwint(2,lastint) ! Last interval is always the lowest
      endif
      !
      ! iloss is the value stored to correspond to LOSSTYPE entered by user.
      iloss = nint(mnw2(3,i))
      select case (iloss)
      case (0)
        ! "NONE" in mnw2 input: single-node well
        ! There are no well corrections and the head in the well is assumed to
        ! equal the head in the cell. This option (hWELL = hn) is only valid
        ! for a single-node well (NNODES = 1). (This is equivalent to using
        ! the original WEL Package of MODFLOW, but specifying the single-node
        ! well within the MNW2 Package enables the use of constraints.)
        ! ERB: Probably just specify a large CWC to approximate no loss.
        condeqn = 'SPECIFIED'
        if (radius <= 0.0d0) radius = 1.0d0
        if (cwcmsg /= '') then
          call store_warning(cwcmsg)
          cwcmsg = ''
        endif
        !if (nonemsg /= '') then
        !  call store_error(nonemsg)
        !  nonemsg = ''
        !endif
        !condeqn = 'INVALID'
        !mnw2(3,i) = -9  ! flags invalid LossType
      case (1)
        ! "THIEM" in mnw2 input
        ! This option allows for only the cell-to-well correction at the well
        ! based on the Thiem (1906) equation; head in the well is determined
        ! from equation 2 as (hWELL = hn + AQn), and the model computes A on
        ! the basis of the user-specified well radius (Rw) and previously
        ! defined values of cell transmissivity and grid spacing. Coefficients
        ! B and C in equation 2 are automatically set = 0.0. User must define
        ! Rw in dataset 2c or 2d.
        ! If LOSSTYPE = THIEM, then specify Rw (the radius of the well).
        condeqn = 'THIEM'
      case (2)
        ! "SKIN" in mnw2 input: user specified characteristics of skin.
        ! This option allows for formation damage or skin corrections at the
        ! well: hWELL = hn + AQn + BQn (from equation 2), where A is determined
        ! by the model from the value of Rw, and B is determined by the model
        ! from Rskin and Kskin. User must define Rw, Rskin, and Kskin in
        ! dataset 2c or 2d.
        ! If LOSSTYPE = SKIN, then specify Rw, Rskin (the radius to the outer
        ! limit of the skin), and Kskin (the hydraulic conductivity of the
        ! skin).
        condeqn = 'MEAN'
      case (3)
        ! "GENERAL" in mnw2 input: user specifies B in eqn. 10.
        ! Head loss is defined with coefficients A, B, and C and power exponent
        ! P (hWELL = hn + AQn + BQn + CQnP). A is determined by the model from
        ! the value of Rw. User must define Rw, B, C, and P in dataset 2c or 2d.
        ! A value of P = 2.0 is suggested if no other data are available (the
        ! model allows 1.0 <= P <= 3.5). Entering a value of C = 0 will result
        ! in a "linear" model in which the value of B is entered directly
        ! (rather than entering properties of the skin, as with the SKIN
        ! option).
        ! ERB: The MEAN option is probably closest--it's the most general.
        ! If LOSSTYPE = GENERAL, then specify Rw, B, C, and P, where the last
        ! three parameters are coefficients in the well-loss equation
        ! (equation 2).
        if (genmsg /= '') then
          call store_error(genmsg)
          genmsg = ''
        endif
        condeqn = 'INVALID'
        mnw2(3,i) = -9  ! flags invalid LossType
      case (4)
        ! "SPECIFYCWC" in mnw2 input.
        ! The user specifies an effective conductance value (equivalent to the
        ! combined effects of the A, B, and C well-loss coefficients expressed
        ! in equation 15) between the well and the cell representing the
        ! aquifer, CWC. User must define CWC in dataset 2c or 2d. If there are
        ! multiple screens within the grid cell or if partial penetration
        ! corrections are to be made, then the effective value of CWC for the
        ! node may be further adjusted automatically by MNW2.
        ! If LOSSTYPE = SPECIFYcwc, then specify CWC (the cell-to-well
        ! conductance; see equation 15).
        condeqn = 'SPECIFIED'
        if (radius <= 0.0d0) radius = 1.0d0
        ! Ned todo: Problem? Ztop and Zbotm are needed, but are they read by MNW2?
      case default
      end select
      !
      write(iu,20)i, radius, bottom, strt, trim(condeqn), nnodes, &
                  trim(WELLID(i))
    enddo
    !
    write(iu,10) 'END PACKAGEDATA'
    !
    if (count_errors() > 0) call ustop()
    !
    return
  end subroutine WriteWellsBlock

  subroutine WriteWellConnections(this)
    ! For each well, write:
    ! WELL wellno iconn layer row column scrn_top scrn_bot hk_skin radius_skin
    implicit none
    class(MawPackageWriterType) :: this
    !
    ! local
    integer :: ic, ii, il, ir, iu, iw, n, nn, nnodes, nnsigned, nodnum
    integer :: firstint, lastint, nintvl, firstnode, lastnode, iloss
    integer :: klay
    double precision :: scrn_top, scrn_bot, hk_skin, radius_skin
    ! formats
    1 format()
    10 format(a)
    30 format(5(2x,i0),2x,4(2x,g15.8))
    !
    iu = this%fileobj%IUnit
    write(iu,1)
    write(iu,10) 'BEGIN CONNECTIONDATA'
    !
    do iw=1,this%MaxActiveBnd
      iloss =  nint(mnw2(3,iw))
      if (iloss == -9) cycle ! indicates invalid LossType
      nnsigned = nint(mnw2(2,iw))
      nnodes = abs(nnsigned)
      nodnum = nint(mnw2(4,iw))
      firstnode = nint(mnw2(4, iw))
      lastnode = firstnode + abs(nint(mnw2(2, iw))) - 1
      if (nnsigned > 0) then
        ! Nodes are specified in MNW2 input
        do n=1,nnodes
          nn = nodnum+n-1
          ! Get layer, row, column
          il = nint(mnwnod(1, nn))
          ir = nint(mnwnod(2, nn))
          ic = nint(mnwnod(3, nn))
          ! Get scrn_top, scrn_bot, hk_skin, radius_skin
          !scrn_top = mnwnod(20, nn)
          !scrn_bot = mnwnod(21, nn)
          ! Set screen top and bottom based on the top and
          ! bottom of the cell
          klay = LBOTM(il)
          scrn_top = BOTM(ic,ir,klay-1)
          scrn_bot = BOTM(ic,ir,klay)
          if (iloss == 0) then
            ! LossType = NONE
            hk_skin = BIGHKVALUE
          elseif (iloss == 2) then
            ! LossType = SKIN
            hk_skin = mnwnod(7, nn)   ! Kskin
          elseif (iloss == 4) then
            ! LossType = SPECIFYCWC
            hk_skin = mnwnod(11, nn)  ! CWC
          else
            ! hk_skin value ignored for LossType THIEM (1) or GENERAL (3)
            hk_skin = DZERO
          endif
          radius_skin =  mnwnod(6, nn)
          write (iu,30) iw, n, il, ir, ic, scrn_top, scrn_bot, hk_skin, &
                        radius_skin
        enddo
      else
        ! Intervals are specified in MNW2 input
        ! Get first and last interval for this well
        firstint = nint(mnwnod(12, firstnode))
        lastint = nint(mnwnod(13, lastnode))
        nintvl = lastint - firstint + 1
        ! loop through intervals
        do ii=firstint,lastint
          ! Get row, column
          ir = nint(mnwint(3, ii))
          ic = nint(mnwint(4, ii))
          ! Get scrn_top, scrn_bot, hk_skin, radius_skin
          scrn_top = mnwint(1, ii)
          scrn_bot = mnwint(2, ii)
          if (iloss == 0) then
            ! LossType = NONE
            hk_skin = BIGHKVALUE
          elseif (iloss == 2) then
            ! LossType = SKIN
            hk_skin = mnwint(7, ii)   ! Kskin
          elseif (iloss == 4) then
            ! LossType = SPECIFYCWC
            hk_skin = mnwint(11, ii)  ! CWC
          else
            ! hk_skin value ignored for LossType THIEM (1) or GENERAL (3)
            hk_skin = DZERO
          endif
          radius_skin = mnwint(6, ii)
          ! loop through nodes intersected by current interval
          do nn = firstnode,lastnode
            ! Get layer
            il = nint(mnwnod(1, nn))
            write(iu,30) iw, nn - firstnode + 1, il, ir, ic, scrn_top, &
                         scrn_bot, hk_skin, radius_skin
          enddo
        enddo
      endif
    enddo
    !
    write(iu,10) 'END CONNECTIONDATA'
    !
    return
  end subroutine WriteWellConnections

  subroutine WriteInitialAttributes(this)
    implicit none
    ! dummy
    class(MawPackageWriterType) :: this
    ! local
    integer :: i, inode, iu, j, nv, n, nnodes, nodnum, numint, nodes
    integer :: pumploc, Qlimit, pumplay, pumprow, pumpcol
    double precision :: pump_elevation, radius, radiustemp, reduction_length
    double precision :: Hlim
    logical :: radius_warning_given
    character(len=1000) :: line
    character(len=100) :: msg
    character(len=20) :: cval
    character(len=10) :: LossTypeLocal
    ! formats
    1 format()
    10 format(a)
    50 format(2x,'WELL',2x,i0,2x,a,2x,g16.9)
    60 format(2x,'WELL',2x,i0,2x,a)
    70 format(2x,g16.9)
    80 format(2x,i0,2x,a,2x,g16.9,2x,g16.9)
    !
    radius_warning_given = .false.
    iu = this%fileobj%IUnit
    write(iu,1)
    write(iu,10) 'BEGIN Well_Initial_Attributes'
    !
    do i=1,this%MaxActiveBnd
      !
      ! One value per well:
      ! IBOUND, default = 1
!      write(iu,30)i,'IBOUND',1
      ! Determine if well discretization is by nodes or by intervals
      nnodes = nint(MNW2(2,i))
      numint = 0
      if (nnodes < 0) then
        nodes = abs(nnodes)
        numint = abs(nint(MNW2(13,i)))
        nv = numint
      else
        nv = abs(nnodes)
      endif
      nodnum = nint(MNW2(4,i))
      radius = -99.9
      LossTypeLocal = LossType(nint(MNW2(3,i)))
      do n=1,nv
        inode = NODNUM + n - 1
        if (n==1) then
          if (numint>0) then
            radius =  MNWINT(5,inode)
          else
            radius = MNWNOD(5,inode)
          endif
        else
          ! Check for varying radius within well
          if (.not. radius_warning_given) then
            if (numint>0) then
              radiustemp =  MNWINT(5,inode)
            else
              radiustemp = MNWNOD(5,inode)
            endif
            if (radiustemp /= radius) then
              msg = 'MAW8 does not support varying radius within a well.'
              call store_warning(msg)
              radius_warning_given = .true.
            endif
          endif
        endif
      enddo
      ! RADIUS, required
      write(iu,50)i,'RADIUS',radius
      ! BOTTOM_ELEVATION, optional
      ! RATE_SCALING pump_elevation reduction_length, optional
      !     (MNW2 similar: Qlimit/=0, items 2f & 4b)
      pumploc = nint(MNW2(11,i))
      ! PUMP_ELEVATION = Hlim - REDUCTION_LENGTH
      if (pumploc > 0) then
        ! cell in which the intake (or outflow) is located will be specified in
        ! dataset 2e as a LAY-ROW-COL grid location.
        pumplay = nint(MNW2(14,i))
        pumprow = nint(MNW2(15,i))
        pumpcol = nint(MNW2(16,i))
        pump_elevation = BOTM(pumpcol, pumprow, LBOTM(pumplay))
        reduction_length = radius
        ! support for rate_scaling.
        write(iu,80)i, 'RATE_SCALING', pump_elevation, reduction_length
      elseif (pumploc < 0) then
        ! vertical position of the pump intake (or outflow) is an elevation in
        ! dataset 2e
        pump_elevation = MNW2(31,i)
        reduction_length = radius
        ! support for rate_scaling.
        write(iu,80)i, 'RATE_SCALING', pump_elevation, reduction_length
      endif
! Currently, nothing is done with Qlimit - erb
      Qlimit = nint(MNW2(6,i))
      if (Qlimit /= 0.0) then
!        msg = 'MAW8 does not support recharge or discharge constraints' &
!              // ' to limit water level in well.'
!        call store_warning(msg)
! If Qlimit /= 0 -- then pumpage will be limited (constrained) by the water level
!                   in the well:
!                   Constant: dataset 2f
!                   Varying:  dataset 4b
      endif
      if (Qlimit > 0) then
        ! pumping will be limited (constrained) by the water level in
        ! the well, and relevant parameters are constant in time. (dataset 2f)
        ! Obtain info from dataset 2f:
        ! MNW2(7,i)=Hlim | the limiting water level (head) in the well, which
        !                  is a minimum for discharging wells and a maximum
        !                  for injection wells
        ! Only use Qlimit to specify RATE_SCALING if pumploc==0
        if (pumploc == 0) then
          Hlim = MNW2(7,i)
          pump_elevation = Hlim
          reduction_length = radius
          write(iu,80)i, 'RATE_SCALING', pump_elevation, reduction_length
        endif
        ! MNW2(8,i)=QCUT
!        Qcut = MNW2(8,i)
!        if (Qcut /= 0.0) then
!          Qfrcmn = MNW2(9,i)
!          Qfrcmx = MNW2(10,i)
!          if (Qcut < 0.0) then
!            ! Qfrcmn and Qfrcmx are fractions of Qdes
!            Qfrcmn = MNW2(9,i) * MNW2(5,i)
!            Qfrcmx = MNW2(10,i) * MNW2(5,i)
!          endif
!        endif
        ! MNW2(9,i)=Qfrcmn | the minimum pumping rate or fraction of original
        !                    pumping rate (a choice that depends on QCUT)
        !                    that a well must exceed to remain active during
        !                    a stress period. The absolute value of Qfrcmn
        !                    must be less than the absolute value of Qfrcmx.
        !                    Only specify if QCUT /= 0.
        ! MNW2(10,i)=Qfrcmx | the minimum pumping rate or fraction of original
        !                     pumping rate that must be exceeded to reactivate
        !                     a well that had been shut off based on Qfrcmn
        !                     during a stress period. The absolute value of
        !                     Qfrcmx must be greater than the absolute value of
        !                     Qfrcmn. Only specify if QCUT /= 0.
      elseif (Qlimit < 0) then
        ! pumping will be limited (constrained) by the water level in
        ! the well, and relevant parameters can vary with time. (dataset 4b)
        ! Obtain info from dataset 4b each stress period.
      endif
      ! STRT, default is top of the smallest GWF node connected to a MAW well
      !   -- MNW2 has no starting head, so leave as default, which is top of first node.
      !
      ! List of values needed for connections:
      write(line,60)i,'HK'
      do n=1,nv
        inode = NODNUM + n - 1
        if (numint > 0) then
          do j=1,nodes
            write(cval,70)MNWINT(7,inode)
            line = trim(line) // trim(cval)
          enddo
        else
          write(cval,70)MNWNOD(7,inode)
          line = trim(line) // trim(cval)
        endif
      enddo
      write(iu,'(a)')trim(line)
      ! If SCREEN TOP and SCREEN BOTTOM are not specified, the Thiem equation
      !     is used to calculate the well conductance for each GWF node connected
      !     to the MAW well.
      select case (LossTypeLocal)
      case ('NONE')
        ! Warn that NONE loss type is not supported.
        msg = 'Loss type NONE is not supported--use THIEM or SKIN.'
        call store_warning(msg)
      case ('THIEM')
        ! No additional attributes are needed.
      case ('SKIN')
        ! RADIUS_SKIN values (do not need to be specified if the Thiem equation is
        !     used to calculate the well conductance)  [MNWNOD(6,__) = Rskin]
        write(line,60)i,'RADIUS_SKIN'
        do n=1,nv
          inode = NODNUM + n - 1
          if (numint > 0) then
            do j=1,nodes
              write(cval,70)MNWINT(6,inode)
              line = trim(line) // trim(cval)
            enddo
          else
            write(cval,70)MNWNOD(6,inode)
            line = trim(line) // trim(cval)
          endif
        enddo
        write(iu,'(a)')trim(line)
        ! SCREEN_TOP values [MNWNOD(20,INODE)=topscreen]
        write(line,60)i,'SCREEN_TOP'
        do n=1,nv
          inode = NODNUM + n - 1
          if (numint > 0) then
            do j=1,nodes
              write(cval,70)MNWINT(1,inode)
              line = trim(line) // trim(cval)
            enddo
          else
            write(cval,70)MNWNOD(20,inode)
            line = trim(line) // trim(cval)
          endif
        enddo
        write(iu,'(a)')trim(line)
        ! SCREEN_BOTTOM values [MNWNOD(21,INODE)=bottomscreen]
        write(line,60)i,'SCREEN_BOTTOM'
        do n=1,nv
          inode = NODNUM + n - 1
          if (numint > 0) then
            do j=1,nodes
              write(cval,70)MNWINT(2,inode)
              line = trim(line) // trim(cval)
            enddo
          else
            write(cval,70)MNWNOD(21,inode)
            line = trim(line) // trim(cval)
          endif
        enddo
        write(iu,'(a)')trim(line)
      case ('GENERAL')
        ! Warn that GENERAL loss type is not supported.
        msg = 'Loss type GENERAL is not supported--use THIEM or SKIN.'
        call store_warning(msg)
      case ('SPECIFYCWC')
        ! Warn that SPECIFYcwc loss type is not supported.
        msg = 'Loss type SPECIFYcwc is not supported--use THIEM or SKIN.'
        call store_warning(msg)
      end select
    enddo
    !
    write(iu,10) 'END Well_Initial_Attributes'
    !
    return
  end subroutine WriteInitialAttributes

  function MyType(this) result (ctype)
    ! dummy 
    class(MawPackageWriterType) :: this
    character(len=LENCTYPE) :: ctype
    !
    ctype = 'MawPackageWriterType'
    !
    return
  end function MyType

end module MawPackageWriterModule
  
subroutine junk()
! From MNW2 source code:
!
! IF (NODTOT.EQ.0) THEN
!   NODTOT=(MNWMAX*NLAY)+(10*NLAY)+25
! END IF
! ALLOCATE (MNWNOD(34,NODTOT))
! ALLOCATE (MNWINT(11,NODTOT))
! ALLOCATE (CapTable(mnwmax,27,2))
! ALLOCATE (WELLID(mnwmax+1))
! ALLOCATE(LIMQ(3,MNWMAX))
!
!     MNW2 array:
! NMNWVL=30+NAUX
! ALLOCATE (MNW2(NMNWVL,MNWMAX))
!     The rows of the MNW2 array store:
!      Row #  = Description
!------------------------------------------------------------------
!         1   = IACTIV (0: inactive; 1: active)
!         2   = NNODES (number of nodes in this well)
!         3   = LOSSTYPE (0: none; 1:THIEM; 2: SKIN; 3: GENERAL; 4:SPEC. COND[=SPECIFYCWC])
!         4   = NODNUM (number, in node list (MNWNOD), of first node of well)
!         5   = QDES (desired flow rate for this stress period)
!         6   = QLIMIT (pumpage constraint flag, QLIMIT>0 turns on constraint)
!         7   = HLIM (limiting water level for pumpage constraint)
!         8   = QCUT (pump cutoff flag: QCUT>0 limit by rate, QCUT<0 limit by
!                     fraction of QDES)
!         9   = Qfrcmn (minimum rate for well to remain active)
!        10   = Qfrcmx (maximum rate to reactivate)
! MNW2(11,MNWID)=PUMPLOC
! MNW2(12,MNWID)=Cprime
! MNW2(13,MNWID)=INTNUM [= interval number]
! MNW2(14,MNWID)=PUMPLAY
! MNW2(15,MNWID)=PUMPROW
! MNW2(16,MNWID)=PUMPCOL
! hwell = mnw2(17,iw)
! mnw2(18,iw)=mnw2(5,iw)   [ = Qdes ]
!     set PPFLAG flag.  PPFLAG>0 means calculate partial penetration effect.
! MNW2(19,MNWID)=PPFLAG
!     set HWflag (horizontal well flag) to 0 as default
!     can be set to 1 (true) if NNODES>0 and any R,C is different than first
!c  if qlimit was hit in a previous FM routine, set qnet to updated Q, then reset flags
!          if(mnw2(20,iw).ne.0.or.MNW2(27,iw).gt.0) then
!            qnet=mnw2(18,iw)
!            mnw2(20,iw)=0
!            mnw2(27,iw)=0
! MNW2(21,MNWID)=0
! MNW2(22,MNWID)=PUMPCAP
! mnw2(23,MNWID)=Hlift
! MNW2(24,MNWID)=CapMult
! mnw2(28,MNWID)=HWtol
! mnw2(30,iw) = Qpot
! mnw2(31,mnwid) = Zpump   !erb: added 3/9/15
!
!     MNWNOD array:
! ALLOCATE (MNWNOD(34,NODTOT))
! MNWNOD(1,NODNUM+INODE-1)=IL       -- regardless of NNODES sign
! MNWNOD(2,NODNUM+INODE-1)=IR       -- regardless of NNODES sign
! MNWNOD(3,NODNUM+INODE-1)=IC       -- regardless of NNODES sign
! MNWNOD(4,NODNUM+INODE-1)=0.0
! MNWNOD(4,firstnode)=Qdes
! MNWNOD(4,firstnode)=Qpot
! MNWNOD(5,NODNUM+INODE-1)=RwNode     -- if NNODES > 0
! MNWNOD(6,NODNUM+INODE-1)=RskinNode  -- if NNODES > 0
! MNWNOD(7,NODNUM+INODE-1)=KskinNode  -- if NNODES > 0
! MNWNOD(8,NODNUM+INODE-1)=BNode      -- if NNODES > 0
! MNWNOD(9,NODNUM+INODE-1)=CNode      -- if NNODES > 0
! MNWNOD(10,NODNUM+INODE-1)=PNode     -- if NNODES > 0
! MNWNOD(11,NODNUM+INODE-1)=CWCNode   -- if NNODES > 0
!     set first and last intervals in this node= interval #1
!     INTNUM is location in int list of 1st interval for this well
!     IINT is loop counter for intervals
! MNWNOD(12,NODNUM)=INTNUM+IINT-1
! MNWNOD(13,NODNUM)=INTNUM+IINT-1
!     get first and last interval intersecting this node
! firstint=MNWNOD(12,INODE)
! lastint=MNWNOD(13,INODE)
! mnwnod(14,firstnode)=cond
! hlim = MNWNOD(15,INODE)   !!  Water level in wellbore, [pos 15 used with seepage face calc.]
! MNWNOD(16,INODE)=Txx
! MNWNOD(17,INODE)=Tyy
!     store partial penetration effect (dhp)
! MNWNOD(18,INODE)=dhp
! MNWNOD(19,NODNUM+INODE-1)=PP    (partial penetration)
! alpha=MNWNOD(19,INODE)
! MNWNOD(20,INODE)=topscreen
! MNWNOD(21,INODE)=bottomscreen
! MNWNOD(22,*) is unused
! MNWNOD(23,INODE+1)=lxf
! MNWNOD(24,INODE)=lxf
! MNWNOD(25,INODE)=lbf
! MNWNOD(25,INODE)=bot1-top2
! MNWNOD(26,INODE)=z1
! MNWNOD(26,INODE+1)=z2
! MNWNOD(27,firstnode)=-Qnet
! MNWNOD(28,INODE)=omega
! MNWNOD(29,INODE)=theta
! MNWNOD(30,INODE)=cond1
! MNWNOD(31,INODE)=cond2
! MNWNOD(32,*) is unused
! MNWNOD(33,INODE)=Kz
! MNWNOD(34,INODE)=SS
!
!
! MNWNOD will access MNWINT by pointing to the first and last interval
!     that intersects the node
!
!
!     MNWINT array:
! ALLOCATE (MNWINT(11,NODTOT))
! MNWINT(1,INTNUM+IINT-1)=Ztop        -- if NNODES < 0
! MNWINT(2,INTNUM+IINT-1)=Zbotm       -- if NNODES < 0
! MNWINT(3,INTNUM+IINT-1)=IR          -- if NNODES < 0
! MNWINT(4,INTNUM+IINT-1)=IC          -- if NNODES < 0
! MNWINT(5,INTNUM+IINT-1)=RwNode      -- if NNODES < 0
! MNWINT(6,INTNUM+IINT-1)=RskinNode   -- if NNODES < 0
! MNWINT(7,INTNUM+IINT-1)=KskinNode   -- if NNODES < 0
! MNWINT(8,INTNUM+IINT-1)=BNode       -- if NNODES < 0
! MNWINT(9,INTNUM+IINT-1)=CNode       -- if NNODES < 0
! MNWINT(10,INTNUM+IINT-1)=PNode      -- if NNODES < 0
! MNWINT(11,INTNUM+IINT-1)=CWCNode    -- if NNODES < 0
end subroutine junk
