module ExchangeModule
  
  use ConnectionModule, only: ConnectionType
  use ConstantsModule, only: DZERO, DHALF, DONE, DTWO, MAXCHARLEN
  use FileWriterModule, only: FileWriterType
  use GLOBAL, only: NLAY, NROW, NCOL, IBOUND, BOTM, DELC, DELR, LBOTM, NBOTM
  use GlobalVariablesModule, only: LgrBilinear
  use GeomUtilModule, only: get_ijk, get_node
  use LGRMODULE, only: IBFLG, NPLBEG, NPRBEG, NPCBEG, NPLEND, &
                       NPREND, NPCEND, NCPP, NCPPL
  use ModelModule, only: ModelType
  use SimPHMFModule, only: store_note, store_warning
  
  type :: ExchangeType
    character(len=9)  :: ExgType = 'GWF6-GWF6'
    character(len=11) :: CellAveraging = 'HARMONIC'
    integer           :: Nexg = 0
    integer           :: Numalphaj = 1
    logical :: PrintInput = .true.
    logical :: PrintFlows = .true.
    logical :: SaveFlows = .true.
    logical :: VariableCV = .false.
    logical :: Dewatered = .false.
    logical :: Newton = .false.
    logical :: GhostNodes = .true.
    logical :: Bilinear = .false.
    logical, dimension(5) :: faceconnect
    type(ModelType), pointer  :: Model1 => null()
    type(ModelType), pointer  :: Model2 => null()
    type(FileWriterType)      :: FileWriter
    type(FileWriterType)      :: GncFileWriter
    type(ConnectionType), pointer, dimension(:) :: Connections
  contains
    ! public procedures
    procedure, public :: InitializeExchange
    procedure, public :: ModifyIdomainForChild
    procedure, public :: WriteExchangeFile
    procedure, public :: WriteGhostNodeFile
    ! private procedures
    procedure, private :: CalcContribFactors
    procedure, private :: CalcNexg
    procedure, private :: DefineFaceconnect
    procedure, private :: DefineConnections
    procedure, private :: DefineGhostNodeCorrection
    procedure, private :: DefineOneConnection
    procedure, private :: FaceProc
    procedure, private :: NextFace
    procedure, private :: WriteDimensionsBlock
    procedure, private :: WriteExchangeDataBlock
    procedure, private :: WriteGhostNodeDimensions
    procedure, private :: WriteGhostNodeOptions
    procedure, private :: WriteGncdata
    procedure, private :: WriteOptionsBlock
  end type ExchangeType
  
contains

  ! Type-bound procedures

  subroutine InitializeExchange(this, Model1, Model2)
    implicit none
    ! dummy
    class(ExchangeType) :: this
    type(ModelType), pointer :: Model1
    type(ModelType), pointer :: Model2
    ! formats
    10 format(a,i0,a,i0,a)
    !
    this%Model1 => Model1
    this%Model2 => Model2
    this%Bilinear = LgrBilinear
    !
    !! Assign filename for exchange input file.
    !write(this%FileWriter%fileobj%FName,10) &
    !    'gwf',Model1%IGrid,'-gwf',Model2%IGrid,'.exg'
    !
    return
  end subroutine InitializeExchange
  
  subroutine WriteExchangeFile(this, mvrFileName)
    implicit none
    ! dummy
    class(ExchangeType), intent(inout) :: this
    character(len=*), intent(in) :: mvrFileName
    ! local
    !
    call this%WriteOptionsBlock(mvrFileName)
    call this%CalcNexg()
    call this%WriteDimensionsBlock()
    call this%DefineConnections()
    call this%WriteExchangeDataBlock()
    !
    return
  end subroutine WriteExchangeFile
  
  subroutine WriteOptionsBlock(this, mvrFileName)
    implicit none
    ! dummy
    class(ExchangeType), intent(inout) :: this
    character(len=*), intent(in) :: mvrFileName
    ! local
    integer :: iu
    ! formats
    1 format()
    10 format(a)
    15 format(2x,a)
    20 format(2x,a,2x,a)
    !
    iu = this%FileWriter%fileobj%IUnit
    !
    write(iu,1)
    ! Write Options block
    write(iu,10)'Begin Options'
    ! Print_Input
    if (this%PrintInput) write(iu,15)'PRINT_INPUT'
    ! Print_Flows
    if (this%PrintFlows) write(iu,15)'PRINT_FLOWS'
    ! Save_Flows
    if (this%SaveFlows) write(iu,15)'SAVE_FLOWS'
    ! Cell_Averaging
    if (this%CellAveraging /= 'HARMONIC') then
      write(iu,20)'ALTERNATIVE_CELL_AVERAGING', trim(this%CellAveraging)
    endif
    ! VariableCV [Dewatered]
    if (this%VariableCV) then
      if (this%Dewatered) then
        write(iu,20)'VARIABLECV','DEWATERED'
      else
        write(iu,15)'VARIABLECV'
      endif
    endif
    ! Newton
    if (this%Newton) then
      write(iu,15)'NEWTON'
    endif
    ! Water movers
    if (mvrFileName /= '') then
      write(iu,20)'MVR6 FILEIN',trim(mvrFileName)
    endif
    ! Ghost_Nodes
    if (this%GhostNodes) then
      write(iu,20)'GNC6 FILEIN', trim(this%GncFileWriter%fileobj%FName)
    endif
    write(iu,10)'End Options'
    !
    return
  end subroutine WriteOptionsBlock
  
  subroutine WriteDimensionsBlock(this)
    implicit none
    ! dummy
    class(ExchangeType), intent(inout) :: this
    ! local
    integer :: iu
    ! formats
    1 format()
    10 format(a)
    20 format(2x,a,2x,i0)
    !
    iu = this%FileWriter%fileobj%IUnit
    !
    write(iu,1)
    write(iu,10)'Begin Dimensions'
    write(iu,20)'NEXG', this%Nexg
    write(iu,10)'End Dimensions'
    !
    return
  end subroutine WriteDimensionsBlock

  subroutine CalcNexg(this)
    ! Count connections required by exchange and
    ! assign number to this%Nexg.
    ! Assumes Model1 is parent and Model2 is child in LGR
    implicit none
    ! dummy
    class(ExchangeType), intent(inout) :: this
    ! local
    integer :: nlayp, nrowp, ncolp  ! parent grid dimensions
    integer :: nlayc, nrowc, ncolc  ! child grid dimensions
    integer :: i, ibflag, j, k, ncon
    integer, pointer, dimension(:,:,:) :: ibp => null()
    integer, pointer, dimension(:,:,:) :: ibc => null()
    integer :: kp
    character(len=MAXCHARLEN) :: msg
    !
    ! Set up grid-specific variables
    ! parent
    call this%Model1%PointToGrid()
    nlayp = this%Model1%Nlaynew
    nrowp = NROW
    ncolp = NCOL
    ibp => IBOUND
    ! child
    call this%Model2%PointToGrid()
    ! From here down, LGRMODULE variables point to child grid.
    nlayc = this%Model2%Nlaynew
    nrowc = NROW
    ncolc = NCOL
    ibc => IBOUND
    this%Model2%Ibflg = IBFLG
    ibflag = IBFLG
    !
    ! Determine if ghost-node correction problem is 1D (default) or 2D
    kp = NPLBEG  ! current parent layer index. In LGR 2.0, NPLBEG = 1.
    do kp=NPLBEG,NPLEND
      if (NCPPL(kp) > 1) then
        if (this%Bilinear) then
          ! For bilinear interpolation, Numalphaj = 3
          this%Numalphaj = 3 ! GNC problem is 2D, use bilinear interpolation.
          msg = 'Using bilinear interpolation for ghost nodes in LGR conversion.'
          call store_note(msg)
        else
          ! For triangular interpolation, Numalphaj = 2
          this%Numalphaj = 2 ! GNC problem is 2D, use triangular interpolation.
        endif
        exit
      endif
    enddo
    !
    ncon = 0
    do k=1,nlayc
      ! Horizontal connections (corners and edges are 
      ! counted two or three times as needed):
      ! Left side (column 1) if child column 1 is not at
      ! parent left side
      if (NPCBEG > 1) then
        do i=1,nrowc
          if (ibc(1,i,k) == ibflag) ncon = ncon + 1
        enddo
      endif
      ! Right side (column ncolc) if child column ncolc
      ! is not at parent right side
      if (NPCEND < ncolp) then
        do i=1,nrowc
          if (ibc(ncolc,i,k) == ibflag) ncon = ncon + 1
        enddo
      endif
      ! Back side (row 1) if child row 1 is not at
      ! parent back side
      if (NPRBEG > 1) then
        do j=1,ncolc
          if (ibc(j,1,k) == ibflag) ncon = ncon + 1
        enddo
      endif
      ! Front side (row nrowc) if child row nrowc
      ! is not at parent front side
      if (NPREND < nrowp) then
        do j=1,ncolc
          if (ibc(j,nrowc,k) == ibflag) ncon = ncon + 1
        enddo
      endif
      ! Vertical connections at bottom of child grid
      ! if child layer nlayc is above parent bottom.
      ! Count if NPLEND < nlayp
      if (k == nlayc .and. NPLEND < nlayp) then
        do i=1,nrowc
          do j=1,ncolc
            if (ibc(j,i,k) == ibflag) ncon = ncon + 1
          enddo
        enddo
      endif
    enddo
    !
    this%Nexg = ncon
    ! 
    return
  end subroutine CalcNexg

  subroutine DefineConnections(this)
    implicit none
    ! dummy
    class(ExchangeType), intent(inout) :: this
    ! local
    integer :: nlayp, nrowp, ncolp  ! parent grid dimensions
    integer :: nlayc, nrowc, ncolc  ! child grid dimensions
    integer :: ic, ibflag, jc, kc, icon  ! for child model
    integer :: ip, jp, kp   ! for parent model
    integer :: mlayc, mrowc, mcolc ! child-cell counters for incrementing parent cell indices
    integer :: ncpplk ! # child layers in current parent layer
    integer :: idir
    integer :: ii, jj, kk ! indices used to define parent cell for connection
    integer, pointer, dimension(:,:,:) :: ibp => null() ! IBOUND of parent
    integer, pointer, dimension(:,:,:) :: ibc => null() ! IBOUND of child
    double precision, pointer, dimension(:) :: delcp => null() ! DELR of parent
    double precision, pointer, dimension(:) :: delrp => null() ! DELC of parent
    double precision, pointer, dimension(:,:,:) :: botmp => null() ! BOTM of parent
    integer, pointer :: nbotmp => null() ! NBOTM of parent
    integer, pointer, dimension(:) :: lbotmp => null() ! LBOTM of parent
    double precision, pointer, dimension(:) :: delcc => null() ! DELR of child
    double precision, pointer, dimension(:) :: delrc => null() ! DELC of child
    double precision, pointer, dimension(:,:,:) :: botmc => null() ! BOTM of child
    integer, pointer :: nbotmc => null() ! NBOTM of child
    integer, pointer, dimension(:) :: lbotmc => null() ! LBOTM of child
    !
    ! Allocate connections
    allocate(this%Connections(this%Nexg))
    !
    ! Set up grid-specific variables
    ! parent
    call this%Model1%PointToGrid()
    nlayp = this%Model1%Nlaynew
    nrowp = NROW
    ncolp = NCOL
    ibp => IBOUND
    delcp => DELC
    delrp => DELR
    botmp => BOTM
    nbotmp => NBOTM
    lbotmp => LBOTM
    ! child
    call this%Model2%PointToGrid()
    nlayc = this%Model2%Nlaynew
    nrowc = NROW
    ncolc = NCOL
    ibc => IBOUND
    delcc => DELC
    delrc => DELR
    botmc => BOTM
    nbotmc => NBOTM
    lbotmc => LBOTM
    this%Model2%Ibflg = IBFLG
    ibflag = IBFLG
    ! From here down, LGRMODULE variables point to child grid.
    !
    call this%DefineFaceconnect()
    !
    icon = 0
    kp = NPLBEG  ! current parent layer index
    mlayc = 1
    do kc=1,nlayc
      ncpplk = NCPPL(kp)
      ! Horizontal connections (corner and edge cells have
      ! two or three connections, as needed):
      !
      ! Left side (column 1) if child column 1 is not at
      ! parent left side
      idir = 1
      mrowc = 1
      mcolc = 1
      if (NPCBEG > 1) then
        ip = NPRBEG   ! current parent row index
        jp = NPCBEG   ! current parent column index
        jj = jp - 1
        mrowc = 1
        jc = 1
        ! Iterate over child rows
        do ic=1,nrowc
          if (ibc(jc,ic,kc) == ibflag) then
            icon = icon + 1
            ! Populate a connection
            call this%DefineOneConnection(ip, jj, kp, ic, jc, kc, idir, &
                                          this%Connections(icon))
            ! Set up ghost node correction
            call this%DefineGhostNodeCorrection(idir, mlayc, mrowc, mcolc, ncpplk, &
                          this%Connections(icon),ncolp, nrowp, nlayp, ibp, &
                          delcp, delrp, nbotmp, botmp, lbotmp, ncolc, nrowc, &
                          nlayc, delcc, delrc, nbotmc, botmc, lbotmc)
          endif
          ! Increment counter of child row within current parent row
          mrowc = mrowc + 1
          if (mrowc > NCPP) then
            ! Reset counter of child row within current parent row
            mrowc = 1
            ! Increment parent row
            ip = ip + 1
          endif
        enddo
      endif
      !
      ! Right side (column ncolc) if child column ncolc
      ! is not at parent right side
      mrowc = 1
      mcolc = 1
      if (NPCEND < ncolp) then
        ip = NPRBEG   ! current parent row index
        jp = NPCEND   ! current parent column index
        jj = jp + 1
        mrowc = 1
        jc = ncolc
        ! Iterate over child rows
        do ic=1,nrowc
          if (ibc(jc,ic,kc) == ibflag) then
            icon = icon + 1
            ! Populate a connection
            call this%DefineOneConnection(ip, jj, kp, ic, jc, kc, idir, &
                                          this%Connections(icon))
            ! Set up ghost node correction
            call this%DefineGhostNodeCorrection(idir, mlayc, mrowc, mcolc, ncpplk, &
                          this%Connections(icon),ncolp, nrowp, nlayp, ibp, &
                          delcp, delrp, nbotmp, botmp, lbotmp, ncolc, nrowc, &
                          nlayc, delcc, delrc, nbotmc, botmc, lbotmc)
          endif
          ! Increment counter of child row within current parent row
          mrowc = mrowc + 1
          if (mrowc > NCPP) then
            ! Reset counter of child row within current parent row
            mrowc = 1
            ! Increment parent row
            ip = ip + 1
          endif
        enddo
      endif
      !
      ! Back side (row 1) if child row 1 is not at
      ! parent back side
      idir = 2
      mrowc = 1
      mcolc = 1
      if (NPRBEG > 1) then
        ip = NPRBEG   ! current parent row index
        jp = NPCBEG   ! current parent column index
        ii = ip - 1
        mcolc = 1
        ic = 1
        ! Iterate over child columns
        do jc=1,ncolc
          if (ibc(jc,ic,kc) == ibflag) then
            icon = icon + 1
            ! Populate a connection
            call this%DefineOneConnection(ii, jp, kp, ic, jc, kc, idir, &
                                          this%Connections(icon))
            ! Set up ghost node correction
            call this%DefineGhostNodeCorrection(idir, mlayc, mrowc, mcolc, ncpplk, &
                          this%Connections(icon),ncolp, nrowp, nlayp, ibp, &
                          delcp, delrp, nbotmp, botmp, lbotmp, ncolc, nrowc, &
                          nlayc, delcc, delrc, nbotmc, botmc, lbotmc)
          endif
          ! Increment counter of child column within current parent column
          mcolc = mcolc + 1
          if (mcolc > NCPP) then
            ! Reset counter of child column within current parent column
            mcolc = 1
            ! Increment parent column
            jp = jp + 1
          endif
        enddo
      endif
      !
      ! Front side (row nrowc) if child row nrowc
      ! is not at parent front side
      mrowc = 1
      mcolc = 1
      if (NPREND < nrowp) then
        ip = NPREND   ! current parent row index
        jp = NPCBEG   ! current parent column index
        ii = ip + 1
        mcolc = 1
        ic = nrowc
        ! Iterate over child columns
        do jc=1,ncolc
          if (ibc(jc,ic,kc) == ibflag) then
            icon = icon + 1
            ! Populate a connection
            call this%DefineOneConnection(ii, jp, kp, ic, jc, kc, idir, &
                                          this%Connections(icon))
            ! Set up ghost node correction
            call this%DefineGhostNodeCorrection(idir, mlayc, mrowc, mcolc, ncpplk, &
                          this%Connections(icon),ncolp, nrowp, nlayp, ibp, &
                          delcp, delrp, nbotmp, botmp, lbotmp, ncolc, nrowc, &
                          nlayc, delcc, delrc, nbotmc, botmc, lbotmc)
          endif
          ! Increment counter of child column within current parent column
          mcolc = mcolc + 1
          if (mcolc > NCPP) then
            ! Reset counter of child column within current parent column
            mcolc = 1
            ! Increment parent column
            jp = jp + 1
          endif
        enddo
      endif
      !
      ! Vertical connections at bottom of child grid
      ! if child layer nlayc is above parent bottom.
      ! Count if NPLEND < nlayp
      idir = 0
      mrowc = 1
      mcolc = 1
      if (kc == nlayc .and. NPLEND < nlayp) then
        ip = NPRBEG   ! current parent row index
        kk = kp + 1
        mrowc = 1
        do ic=1,nrowc
          jp = NPCBEG   ! current parent column index
          mcolc = 1
          do jc=1,ncolc
            if (ibc(jc,ic,kc) == ibflag) then
              icon = icon + 1
              ! Populate a connection
              call this%DefineOneConnection(ip, jp, kk, ic, jc, kc, idir, &
                                            this%Connections(icon))
              ! Set up ghost node correction
              call this%DefineGhostNodeCorrection(idir, mlayc, mrowc, mcolc, ncpplk, &
                          this%Connections(icon),ncolp, nrowp, nlayp, ibp, &
                          delcp, delrp, nbotmp, botmp, lbotmp, ncolc, nrowc, &
                          nlayc, delcc, delrc, nbotmc, botmc, lbotmc)
            endif
            ! Increment counter of child column within current parent column
            mcolc = mcolc + 1
            if (mcolc > NCPP) then
              ! Reset counter of child column within current parent column
              mcolc = 1
              ! Increment parent column
              jp = jp + 1
            endif
          enddo
          ! Increment counter of child row within current parent row
          mrowc = mrowc + 1
          if (mrowc > NCPP) then
            ! Reset counter of child row within current parent row
            mrowc = 1
            ! Increment parent row
            ip = ip + 1
          endif
        enddo
      endif
      !
      ! Increment counter of child layer within current parent layer
      mlayc = mlayc + 1
      if (mlayc > ncpplk) then
        ! Reset counter of child layer within current parent layer
        mlayc = 1
        ! Increment parent layer
        kp = kp + 1
      endif
    enddo
    !
    return
  end subroutine DefineConnections
  
  subroutine DefineFaceconnect(this)
    implicit none
    ! dummy
    class(ExchangeType) :: this
    !
    ! Which faces have connections?
    ! left side
    this%faceconnect(1) = (NPCBEG > 1)
    ! right side
    this%faceconnect(2) = (NPCEND < this%Model1%Ncol)
    ! back
    this%faceconnect(3) = (NPRBEG > 1)
    ! front
    this%faceconnect(4) = (NPREND < this%Model1%Nrow)
    ! bottom
    this%faceconnect(5) = (NPLEND < this%Model1%Nlaynew)
    !
    return
  end subroutine DefineFaceconnect
  
  subroutine DefineOneConnection(this, ip, jp, kp, ic, jc, kc, idir, con)
    implicit none
    ! dummy
    class(ExchangeType), intent(inout) :: this
    integer, intent(in) :: ip, jp, kp, ic, jc, kc, idir
    type(ConnectionType), intent(inout) :: con
    ! local
    integer :: ihc, lbotc, lbotp, nodep, nodec
    double precision :: hwva ! hoizontal width or vertical area
    double precision :: clp, clc  ! distances to shared face
    double precision :: heightc, heightp, widthc, lengthc
    double precision :: botmc, botmp, topc, topp
    !
    ! idir = 0 : Connection is vertical
    ! idir = 1 : Connection is in row direction (left or right boundary)
    ! idir = 2 : Connection is in column direction (front or back boundary)
    !
    ! Find parent node number
    nodep = get_node(kp, ip, jp, this%Model1%Nlaynew, this%Model1%Nrow, &
                     this%Model1%Ncol)
    !
    ! Find child node number
    nodec = get_node(kc, ic, jc, this%Model2%Nlaynew, this%Model2%Nrow, &
                     this%Model2%Ncol)
    !
    ! Calculate child cell dimensions
    lbotc = this%Model2%lbotm(kc)
    topc = this%Model2%botm(jc,ic,lbotc-1)
    botmc = this%Model2%botm(jc,ic,lbotc)
    heightc = topc - botmc
    !
    ! Calculate parent cell dimensions
    lbotp = this%Model1%lbotm(kp)
    topp = this%Model1%botm(jp,ip,lbotp-1)
    botmp = this%Model1%botm(jp,ip,lbotp)
    heightp = topp - botmp
    !
    ! Find distance from child node to shared face, distance from
    ! parent node to shared face, and area of flow.
    select case (idir)
    case (0)
      ! connection is vertical
      ihc = 0
      clc = heightc / DTWO
      clp = heightp / DTWO
      widthc = this%Model2%delr(jc)
      lengthc = this%Model2%delc(ic)
      hwva = widthc * lengthc
    case (1)
      ! connection is in row direction (one column to another)
      if (NCPPL(kp) == 1) then
        ihc = 1
      else
        ihc = 2
      endif
      clc = this%Model2%delr(jc) / DTWO
      clp = this%Model1%delr(jp) / DTWO
      widthc = this%Model2%delc(ic)
      hwva = widthc
    case (2)
      ! connection is in column direction (one row to another)
      if (NCPPL(kp) == 1) then
        ihc = 1
      else
        ihc = 2
      endif
      clc = this%Model2%delc(ic) / DTWO
      clp = this%Model1%delc(ip) / DTWO
      widthc = this%Model2%delr(jc)
      hwva = widthc
    end select
    !
    ! Populate connection
    !   parent cell
    con%noden = nodep
    con%kp = kp
    con%ip = ip
    con%jp = jp
    !   child cell
    con%nodem = nodec
    con%kc = kc
    con%ic = ic
    con%jc = jc
    !
    con%ihc = ihc
    con%cl1 = clp
    con%cl2 = clc
    con%hwva =  hwva
    !
    return
  end subroutine DefineOneConnection
  
  subroutine DefineGhostNodeCorrection(this, idir, mlayc, mrowc, mcolc, mcppl, conn, &
                 ncolp, nrowp, nlayp, iboundp, delcp, delrp, nbotmp, botmp, lbotmp, &
                 ncolc, nrowc, nlayc, delcc, delrc, nbotmc, botmc, lbotmc)
    implicit none
    ! dummy
    class(ExchangeType), intent(inout) :: this
    integer, intent(in) :: ncolp, nrowp, nlayp
    ! idir: 0 -- connection is vertical
    !       1 -- connection is in row direction (one column to another)
    !       2 -- connection is in column direction (one row to another)
    ! mlayc: Layer index within parent cell of child or ghost node 
    !       (ranges from 1 to mcppl in vertical downward direction)
    ! mrowc: Row index within parent cell of child or ghost node 
    !       (ranges from 1 to NCPP in column direction)
    ! mcolc: Column index within parent cell of child or ghost node 
    !       (ranges from 1 to NCPP in row direction)
    ! mcppl: NCPPL(current parent layer) = # child cells per current parent layer
    integer, intent(in) :: idir, mlayc, mrowc, mcolc, mcppl
    type(ConnectionType), intent(inout) :: conn
    integer, dimension(ncolp,nrowp,nlayp), intent(in) :: iboundp
    double precision, dimension(nrowp), intent(in) :: delcp
    double precision, dimension(ncolp), intent(in) :: delrp
    integer, intent(in) :: nbotmp
    double precision, dimension(ncolp,nrowp,0:nbotmp), intent(in) :: botmp
    integer, dimension(nlayp), intent(in) :: lbotmp
    integer, intent(in) :: ncolc, nrowc, nlayc
    double precision, dimension(nrowc), intent(in) :: delcc
    double precision, dimension(ncolc), intent(in) :: delrc
    integer, intent(in) :: nbotmc
    double precision, dimension(ncolc,nrowc,0:nbotmc), intent(in) :: botmc
    integer, dimension(nlayc), intent(in) :: lbotmc
    !
    ! local
    ! row, column, layer indices of cell containing ghost node in parent
    integer :: igp = 0
    integer :: jgp = 0
    integer :: kgp = 0
    ! row, column, layer indices of contributing cells 1 and 2 in parent
    integer :: icp1 = 0
    integer :: jcp1 = 0
    integer :: kcp1 = 0
    integer :: icp2 = 0
    integer :: jcp2 = 0
    integer :: kcp2 = 0
    ! row, column, layer indices of child cell
    integer :: ic = 0
    integer :: jc = 0
    integer :: kc = 0
    ! indices used for off-corner cell
    integer :: i12 = 0
    integer :: j12 = 0
    integer :: k12 = 0
    ! position and # child cells w/r/t contributing cell directions 1 and 2
    integer :: mc1 = 0
    integer :: mc2 = 0
    integer :: ncp1 = 0
    integer :: ncp2 = 0
    !real :: rmc1, rmc2, rnpc1, rnpc2
    ! node locations in horizontal and vertical directions
    !integer :: nodeloch, nodelocv
    !logical :: oddh, oddv
    ! Is offset to contributing cells 1 and 2 forward or backward or at node?
    integer :: ioff1 = 0
    integer :: ioff2 = 0
    integer :: ibound12 = 0 ! Iboundp for off-corner cell
    ! distances
    double precision :: d1c, d1p, d2c, d2p
    !  d1c = Distance from parent node (in cell where ghost node is
    !        needed) to child node in first direction.
    !  d1p = Distance from parent node to contributing node in first
    !        direction.
    !  d2c, d2p = as above, but for second direction.
    character(len=MAXCHARLEN) :: ermsg
    ! formats
    10 format('Warning: In DefineGhostNodeCorrection, node number ', &
              ' for parent cell containing ghost node is invalid: ',i0)
    20 format('Warning: In DefineGhostNodeCorrection, node number ', &
              ' for child cell is invalid: ',i0)
    !
    ! Start by assuming no contributing cells, no ghost node correction.
    conn%nodej1 = 0
    conn%nodej2 = 0
    conn%alphaj1 = DZERO
    conn%alphaj2 = DZERO
    d1c = DZERO
    d2c = DZERO
    d1p = DONE
    d2p = DONE
    !
    ! Determine offsets (direction to contributing cells)
    ! relative to node of cell where ghost node resides.
    select case (idir)
    case (0) ! connection is vertical (layer index fixed)
      ! 1 is direction in which row index increases
      ! 2 is direction in which column index increases
      mc1 = mrowc
      ncp1 = NCPP
      mc2 = mcolc
      ncp2 = NCPP
    case (1) ! connection is in row direction (column index fixed)
      ! 1 is direction in which row index increases
      ! 2 is direction in which layer index increases
      mc1 = mrowc
      ncp1 = NCPP
      mc2 = mlayc
      ncp2 = mcppl
    case (2) ! connection is in column direction (row index fixed)
      ! 1 is direction in which column index increases
      ! 2 is direction in which layer index increases
      mc1 = mcolc
      ncp1 = NCPP
      mc2 = mlayc
      ncp2 = mcppl
    end select
    !
    ! Determine offset directions (perpendicular to
    ! connection direction) to contributing cells 1 and 2.
    ! GetIoffset returns -1, 0, or +1.
    ioff1 = GetIoffset(mc1, ncp1)
    ioff2 = GetIoffset(mc2, ncp2)
    if (ioff1 == 0 .and. ioff2 == 0) then
      ! Ghost node not needed because offsets from parent node are zero.
      return
    endif
    !
    ! Assign indices for parent cell containing ghost node
    call get_ijk(conn%noden, nrowp, ncolp, nlayp, igp, jgp, kgp)
    if (igp < 0) then
      ! Node number for parent cell containing ghost node is invalid
      write(ermsg,10)conn%noden
      call store_warning(ermsg)
      return
    endif
    !
    ! Assign indices for child cell
    call get_ijk(conn%nodem, nrowc, ncolc, nlayc, ic, jc, kc)
    if (ic < 0) then
      write(ermsg,20)conn%nodem
      call store_warning(ermsg)
      return
    endif
    !
    ! Determine indices of contributing cells
    ibound12 = 0
    select case (idir)
    case (0) ! connection is vertical (layer index fixed)
      ! Contributing cells have same layer index as ghost node.
      kcp1 = kgp
      ! 1 is direction in which row index increases
      if (ioff1 /= 0) then
        ! Consider contributing cell offset in direction 1 (direction in which row index increases)
        icp1 = igp + ioff1
        ! Column index for contributing cell is same as for cell where ghost node resides
        jcp1 = jgp
        ! Ensure incremented row index is valid
        if (icp1 >= 1 .and. icp1 <= nrowp) then
          if (iboundp(jcp1,icp1,kcp1) /= 0) then
            ! contributing cell 1 is active
            conn%nodej1 = get_node(kcp1, icp1, jcp1, nlayp, nrowp, ncolp)
            conn%k1 = kcp1
            conn%i1 = icp1
            conn%j1 = jcp1
            ! Define d1c, d1p
            d1c = GetGhostHorizDistance(ic, nrowc, delcc, mc1)
            d1p = GetInternodeHorizDistance(igp, nrowp, delcp, ioff1)
          endif
        endif
      endif
      ! Contributing cells have same layer index as ghost node.
      kcp2 = kgp
      ! 2 is direction in which column index increases
      if (ioff2 /= 0) then
        ! Consider contributing cell offset in direction 2 (direction in which column index increases)
        jcp2 = jgp + ioff2
        ! Row index for contributing cell is same as for cell where ghost node resides
        icp2 = igp
        ! Ensure incremented column index is valid
        if (jcp2 >= 1 .and. jcp2 <= ncolp) then
          if (iboundp(jcp2,icp2,kcp2) /= 0) then
            ! contributing cell 2 is active
            conn%nodej2 = get_node(kcp2, icp2, jcp2, nlayp, nrowp, ncolp)
            conn%k2 = kcp2
            conn%i2 = icp2
            conn%j2 = jcp2
            ! Define d2c, d2p
            d2c = GetGhostHorizDistance(jc, ncolc, delrc, mc2)
            d2p = GetInternodeHorizDistance(jgp, ncolp, delrp, ioff2)
          endif
        endif
      endif
      !
      if (ioff1 /= 0 .and. ioff2 /= 0) then
        ! For bilinear interpolation, need nodej12 and iboundp for off-corner cell
        i12 = igp + ioff1
        j12 = jgp + ioff2
        if (i12 >= 1 .and. i12 <= nrowp .and. j12 >= 1 .and. j12 <= ncolp) then
          conn%nodej12 = get_node(kgp, i12, j12, nlayp, nrowp, ncolp)
          conn%k12 = kgp
          conn%i12 = i12
          conn%j12 = j12
          ibound12 = iboundp(j12,i12,kgp)
        endif
      endif
      !
    case (1) ! connection is in direction in which column index increases (column index fixed)
      ! Contributing cells have same column index as ghost node.
      jcp1 = jgp
      ! 1 is direction in which row index increases
      if (ioff1 /= 0) then
        ! Consider contributing cell offset in direction 1 (direction in which row index increases)
        icp1 = igp + ioff1
        ! Layer index for contributing cell is same as for cell where ghost node resides
        kcp1 = kgp
        ! Ensure incremented row index is valid
        if (icp1 >= 1 .and. icp1 <= nrowp) then
          if (iboundp(jcp1,icp1,kcp1) /= 0) then
            ! contributing cell is active
            conn%nodej1 = get_node(kcp1, icp1, jcp1, nlayp, nrowp, ncolp)
            conn%k1 = kcp1
            conn%i1 = icp1
            conn%j1 = jcp1
            ! Define d1c, d1p
            d1c = GetGhostHorizDistance(ic, nrowc, delcc, mc1)
            d1p = GetInternodeHorizDistance(igp, nrowp, delcp, ioff1) ! just edited
          endif
        endif
      endif
      ! Contributing cells have same column index as ghost node.
      jcp2 = jgp
      ! 2 is direction in which layer index increases
      if (ioff2 /= 0) then
        ! Consider contributing cell offset in direction 2 (direction in which layer index increases)
        kcp2 = kgp + ioff2
        ! Row index for contributing cell is same as for cell where ghost node resides
        icp2 = igp
        ! Ensure incremented layer index is valid
        if (kcp2 >= 1 .and. kcp2 <= nlayp) then
          if (iboundp(jcp2,icp2,kcp2) /= 0) then
            ! contributing cell 2 is active
            conn%nodej2 = get_node(kcp2, icp2, jcp2, nlayp, nrowp, ncolp)
            conn%k2 = kcp2
            conn%i2 = icp2
            conn%j2 = jcp2
            ! Define d2c, d2p
            d2c = GetGhostVertDistance(ic, jc, kc, nbotmc, ncolc, nrowc, &
                                       nlayc, botmc, lbotmc, mc2, ncp2)
            d2p = GetInternodeVertDistance(igp, jgp, kgp, ioff2, nbotmp, &
                                           ncolp, nrowp, nlayp, botmp, lbotmp)
          endif
        endif
      endif
      !
      if (ioff1 /= 0 .and. ioff2 /= 0) then
        ! For bilinear interpolation, need nodej12 and iboundp for off-corner cell
        i12 = igp + ioff1
        k12 = kgp + ioff2
        if (i12 >= 1 .and. i12 <= nrowp .and. k12 >= 1 .and. k12 <= nlayp) then
          conn%nodej12 = get_node(k12, i12, jgp, nlayp, nrowp, ncolp)
          conn%k12 = k12
          conn%i12 = i12
          conn%j12 = jgp
          ibound12 = iboundp(jgp,i12,k12)
        endif
      endif
      !
    case (2) ! connection is in direction in which row index increases (row index fixed)
      ! Contributing cells have same row index as ghost node.
      icp1 = igp
      ! 1 is direction in which column index increases
      if (ioff1 /= 0) then
        ! Consider contributing cell offset in direction 1 (direction in which row index increases)
        jcp1 = jgp + ioff1
        ! Layer index for contributing cell is same as for cell where ghost node resides
        kcp1 = kgp
        ! Ensure incremented column index is valid
        if (jcp1 >= 1 .and. jcp1 <= ncolp) then
          if (iboundp(jcp1,icp1,kcp1) /= 0) then
            ! contributing cell is active
            conn%nodej1= get_node(kcp1, icp1, jcp1, nlayp, nrowp, ncolp)
            conn%k1 = kcp1
            conn%i1 = icp1
            conn%j1 = jcp1
            ! Define d1c, d1p
            d1c = GetGhostHorizDistance(jc, ncolc, delrc, mc1)
            d1p = GetInternodeHorizDistance(jgp, ncolp, delrp, ioff1)
          endif
        endif
      endif
      ! Contributing cells have same row index as ghost node.
      icp2 = igp
      ! 2 is direction in which layer index increases
      if (ioff2 /= 0) then
        ! Consider contributing cell offset in direction 2 (direction in which layer index increases)
        kcp2 = kgp + ioff2
        ! Column index for contributing cell is same as for cell where ghost node resides
        jcp2 = jgp
        ! Ensure incremented layer index is valid
        if (kcp2 >= 1 .and. kcp2 <= nlayp) then
          if (iboundp(jcp2,icp2,kcp2) /= 0) then
            ! contributing cell 2 is active
            conn%nodej2 = get_node(kcp2, icp2, jcp2, nlayp, nrowp, ncolp)
            conn%k2 = kcp2
            conn%i2 = icp2
            conn%j2 = jcp2
            ! Define d2c, d2p
            d2c = GetGhostVertDistance(ic, jc, kc, nbotmc, ncolc, nrowc, &
                                       nlayc, botmc, lbotmc, mc2, ncp2)
            d2p = GetInternodeVertDistance(igp, jgp, kgp, ioff2, nbotmp, &
                                           ncolp, nrowp, nlayp, botmp, lbotmp)
          endif
        endif
      endif
      !
      if (ioff1 /= 0 .and. ioff2 /= 0) then
        ! For bilinear interpolation, need nodej12 and iboundp for off-corner cell
        j12 = jgp + ioff1
        k12 = kgp + ioff2
        if (j12 >= 1 .and. j12 <= ncolp .and. k12 >= 1 .and. k12 <= nlayp) then
          conn%nodej12 = get_node(k12, igp, j12, nlayp, nrowp, ncolp)
          conn%k12 = k12
          conn%i12 = igp
          conn%j12 = j12
          ibound12 = iboundp(j12,igp,k12)
        endif
      endif
      !
    end select
    !
    ! Calculate contributing factor(s).
    call this%CalcContribFactors(d1c, d1p, d2c, d2p, ibound12, conn%alphaj1, &
                                 conn%alphaj2, conn%alphaj12)
    !
    return
  end subroutine DefineGhostNodeCorrection

  subroutine WriteExchangeDataBlock(this)
    implicit none
    ! dummy
    class(ExchangeType), intent(inout) :: this
    ! local
    integer :: i, iu, j, k, m
    integer :: nrowc, ncolc, nlayc
    integer :: iface ! 1=left side, 2=right side, 3=back, 4=front, 5=bottom
    character(len=25), dimension(5) :: comment
    logical :: newface
    data comment /'# left side, child layer ', &
                  '# right side, child layer', &
                  '# back, child layer      ', &
                  '# front, child layer     ', &
                  '# bottom, child layer    '/
    ! formats
    1 format()
    10 format(a)
    20 format(2x,a,2x,i0)
    30 format(2x,a,1x,i0)
    !
    iu = this%FileWriter%fileobj%IUnit
    nrowc = this%Model2%Nrow
    ncolc = this%Model2%Ncol
    nlayc = this%Model2%Nlaynew
    !
    write(iu,1)
    write(iu,10)'Begin ExchangeData'
    i = 1
    j = 1
    k = 1
    !
    ! On which face is first connection located?
    iface = this%NextFace(0)
    !
    newface = .true.
    do m=1,this%Nexg
      if (newface) then
        write(iu,30)trim(comment(iface)),k
        newface = .false.
      endif
      call this%Connections(m)%WriteConnection(iu)
      call this%FaceProc(i, j, k, newface, iface)
    enddo
    write(iu,10)'End ExchangeData'
    !
    return
  end subroutine WriteExchangeDataBlock

  subroutine WriteGhostNodeFile(this)
    implicit none
    ! dummy
    class(ExchangeType) :: this
    ! local
    !
    call this%WriteGhostNodeOptions()
    call this%WriteGhostNodeDimensions()
    call this%WriteGncdata()
    !
    return
  end subroutine WriteGhostNodeFile
  
  subroutine WriteGhostNodeOptions(this)
    implicit none
    ! dummy
    class(ExchangeType), intent(inout) :: this
    ! local
    integer :: iu
    ! formats
    1 format()
    10 format(a)
    20 format(2x,a)
    !
    iu = this%GncFileWriter%fileobj%IUnit
    !
    write(iu,1)
    write(iu,10)'BEGIN OPTIONS'
    write(iu,20)'PRINT_INPUT'
    write(iu,20)'PRINT_FLOWS'
    !write(iu,20)'EXPLICIT'
    write(iu,10)'END OPTIONS'
    !
    return
  end subroutine WriteGhostNodeOptions
  
  subroutine WriteGhostNodeDimensions(this)
    implicit none
    ! dummy
    class(ExchangeType), intent(inout) :: this
    ! local
    integer :: iu
    ! formats
    1 format()
    10 format(a)
    20 format(2x,a,2x,i0)
    !
    iu = this%GncFileWriter%fileobj%IUnit
    !
    write(iu,1)
    write(iu,10)'BEGIN DIMENSIONS'
    write(iu,20)'NUMGNC',this%Nexg
    write(iu,20)'NUMALPHAJ',this%Numalphaj
    write(iu,10)'END DIMENSIONS'
    !
    return
  end subroutine WriteGhostNodeDimensions
  
  subroutine WriteGncdata(this)
    implicit none
    ! dummy
    class(ExchangeType), intent(inout) :: this
    ! local
    integer :: i, iu
    integer :: j, k, m
    integer :: nrowc, ncolc, nlayc
    integer :: iface ! 1=left side, 2=right side, 3=back, 4=front, 5=bottom
    character(len=25), dimension(5) :: comment
    logical :: newface
    data comment /'# left side, child layer ', &
                  '# right side, child layer', &
                  '# back, child layer      ', &
                  '# front, child layer     ', &
                  '# bottom, child layer    '/
    ! formats
    1 format()
    10 format(a)
    20 format(3(2x,i0),g0)
    30 format(2x,a,1x,i0)
    !
    iu = this%GncFileWriter%fileobj%IUnit
    nrowc = this%Model2%Nrow
    ncolc = this%Model2%Ncol
    nlayc = this%Model2%Nlaynew
    !
    write(iu,1)
    write(iu,10)'BEGIN GNCDATA'
    i = 1
    j = 1
    k = 1
    !
    ! On which face is first connection located?
    iface = this%NextFace(0)
    !
    newface = .true.
    do m=1,this%Nexg
      if (newface) then
        write(iu,30)trim(comment(iface)),k
        newface = .false.
      endif
      call this%Connections(m)%WriteGhostNodeCorrection(iu,this%Numalphaj)
      call this%FaceProc(i, j, k, newface, iface)
    enddo
    write(iu,10)'END GNCDATA'
    !
    return
  end subroutine WriteGncdata
  
  subroutine FaceProc(this, i, j, k, newface, iface)
    implicit none
    ! dummy
    class(ExchangeType) :: this
    integer, intent(inout) :: i, j, k
    logical, intent(inout) :: newface
    integer, intent(inout) :: iface
    ! local
    integer :: nrowc, ncolc, nlayc
    !
    nrowc = this%Model2%Nrow
    ncolc = this%Model2%Ncol
    nlayc = this%Model2%Nlaynew
    !
    select case (iface)
    case (1)
      ! left side
      i = i + 1
      if (i > nrowc) then
        i = 1
        if (this%NextFace(iface) == 5) then
          k = k + 1
          if (k > nlayc) then
            k = nlayc
            iface = this%NextFace(iface)
          else
            iface = this%NextFace(0)
          endif
        else
          iface = this%NextFace(iface)
        endif
        newface = .true.
      endif
    case (2)
      ! right side
      i = i + 1
      if (i > nrowc) then
        i = 1
        if (this%NextFace(iface) == 5) then
          k = k + 1
          if (k > nlayc) then
            k = nlayc
            iface = this%NextFace(iface)
          else
            iface = this%NextFace(0)
          endif
        else
          iface = this%NextFace(iface)
        endif
        newface = .true.
      endif
    case (3)
      ! back
      j = j + 1
      if (j > ncolc) then
        j = 1
        if (this%NextFace(iface) == 5) then
          k = k + 1
          if (k > nlayc) then
            k = nlayc
            iface = this%NextFace(iface)
          else
            iface = this%NextFace(0)
          endif
        else
          iface = this%NextFace(iface)
        endif
        newface = .true.
      endif
    case (4)
      ! front
      j = j + 1
      if (j > ncolc) then
        j = 1
        if (this%NextFace(iface) == 5) then
          k = k + 1
          if (k > nlayc) then
            k = nlayc
            iface = this%NextFace(iface)
          else
            iface = this%NextFace(0)
          endif
        else
          iface = this%NextFace(iface)
        endif
        newface = .true.
      endif
    case (5)
      ! bottom
      ! nothing more to do
    end select
    !
    return
  end subroutine FaceProc
  
  function NextFace(this, iface) result(iresult)
    ! Return number of next face that has connections.
    implicit none
    ! dummy
    class(ExchangeType) this
    integer, intent(in) :: iface
    ! local
    integer :: i, iresult, istrt
    ! 
    istrt = iface + 1
    do i=istrt,5
      if (this%faceconnect(i)) then
        iresult = i
        return
      endif
    enddo
    !
    ! If control reaches here, no next face.
    iresult = -1
    !
    return
  end function NextFace
    
  subroutine CalcContribFactors(this, d1c, d1p, d2c, d2p, ibound12, cf1, cf2, cf12)
    ! Calculate contributing factors using either triangular or bilinear
    ! interpolation. 
    ! If this%bilinear is true and off-corner cell is active, bilinear is used.
    ! If this%bilinear is false, triangular is used.
    ! If off-corner cell is inactive, triangular is used.
    implicit none
    ! dummy
    class(ExchangeType) :: this
    double precision, intent(in) :: d1c, d1p, d2c, d2p
    ! Where: d1c = Distance from parent node (in cell where ghost node is
    !              needed) to child node in first dimension.
    !        d1p = Distance from parent node to contributing node in first
    !              dimension.
    !        d2c, d2p = as above, but for second dimension.
    !
    integer, intent(in) :: ibound12   ! Ibound for off-corner cell
    ! CF1 is contributing factor for contributing cell
    !     adjacent in dimension 1.
    double precision, intent(inout) :: cf1   
    ! CF2 is contributing factor for contributing cell
    !     adjacent in dimension 2.
    double precision, intent(inout) :: cf2
    ! CF12 is contributing factor for contributing cell
    !      off corner in dimensions 1 and 2.
    double precision, intent(inout) :: cf12
    !
    if (this%Bilinear .and. ibound12 /= 0) then
      call CalcContribFactorsBilinear(d1c, d1p, d2c, d2p, cf1, cf2, cf12)
    else
      cf12 = DZERO
      call CalcContribFactorsTriangular(d1c, d1p, d2c, d2p, cf1, cf2)
    endif
    !
    return
  end subroutine CalcContribFactors
    
    
  ! Non-type-bound procedures
  
  function GetGhostHorizDistance(i, maxi, del, mc) result(hdist)
    ! Returns distance in row or column direction
    ! between parent node and child node.
    implicit none
    ! dummy
    integer,          intent(in) :: i     ! row or column index in child grid
    integer,          intent(in) :: maxi  ! NROW or NCOL of child grid
    double precision, dimension(maxi), intent(in) :: del   ! DELC or DELR of child grid
    integer,          intent(in) :: mc    ! position of child node within parent cell (ranges from 1 to NCPP)
    double precision :: hdist
    ! local
    double precision :: cellwid, pcellwid, rmc, rncpp, px, cx
    !
    rmc = real(mc,kind(DZERO))
    rncpp = real(NCPP,kind(DZERO))
    cellwid = del(i)             ! child cell width
    pcellwid = cellwid * rncpp   ! parent cell width
    px = pcellwid / DTWO         ! parent node X relative to parent cell side
    cx = cellwid * (rmc - DHALF) ! child node X relative to parent cell side
    hdist = abs(px - cx)         ! distance between parent node and child node
    !
    return
  end function GetGhostHorizDistance
  
  function GetGhostVertDistance(ic, jc, kc, nbotmc, ncolc, nrowc, nlayc, &
                                botmc, lbotmc, mc, mcppl) result(vdist)
    ! Returns vertical distance between parent node and child node.
    implicit none
    ! dummy
    integer, intent(in) :: ic, jc, kc ! row, column, layer indices of child cell
    integer, intent(in) :: nbotmc     ! NBOTM of child grid
    integer, intent(in) :: ncolc, nrowc, nlayc  ! dimensions of child grid
    double precision, dimension(ncolc, nrowc, 0:nbotmc), intent(in) :: botmc ! BOTM of child grid
    integer, dimension(nlayc) :: lbotmc  ! LBOTM of child grid
    integer, intent(in) :: mc     ! position of child node within parent cell (ranges from 1 to mcppl)
    integer, intent(in) :: mcppl  ! NCPPL of current layer in parent grid
    double precision :: vdist
    ! local
    double precision :: cellheight, pcellheight, rmc, rmcppl, pz, cz, ctop, cbot
    !
    rmc = real(mc,kind(DZERO))
    rmcppl = real(mcppl,kind(DZERO))
    ctop = botmc(jc,ic,lbotmc(kc)-1)
    cbot = botmc(jc,ic,lbotmc(kc))
    cellheight = ctop - cbot
    pcellheight = cellheight * rmcppl
    pz = pcellheight / DTWO
    cz = cellheight * (rmc - DHALF)
    vdist = abs(pz - cz)
    !
    return
  end function GetGhostVertDistance
  
  function GetInternodeHorizDistance(i, maxi, del, ioff) result(hdist)
    ! Returns horizontal distance between node with row or column 
    ! index i and adjacent node defined by ioff (= -1 or +1). 
    ! Intended for use with parent grid, but also could be used with child.
    implicit none
    ! dummy
    integer,          intent(in) :: i     ! row or column index in child grid
    integer,          intent(in) :: maxi  ! NROW or NCOL of child grid
    double precision, dimension(maxi), intent(in) :: del   ! DELC or DELR of child grid
    integer, intent(in) :: ioff  ! Row or column index offset, either -1 or +1
    double precision :: hdist
    !
    hdist = DHALF * (del(i) + del(i+ioff))
    !
    return
  end function GetInternodeHorizDistance
  
  function GetInternodeVertDistance(ip, jp, kp, ioff, nbotmp, ncolp, nrowp, &
                                    nlayp, botmp, lbotmp) result(vdist)
    ! Returns vertical distance between specified node and vertically 
    ! adjacent node defined by ioff (= -1 or +1).
    ! Intended for use with parent grid, but also could be used with child.
    implicit none
    ! dummy
    integer, intent(in) :: ip, jp, kp ! row, column, layer indices of cell of interest (in parent grid)
    integer, intent(in) :: ioff       ! Layer index offset, either -1 or +1
    integer, intent(in) :: nbotmp     ! NBOTM of parent grid
    integer, intent(in) :: ncolp, nrowp, nlayp  ! dimensions of parent grid
    double precision, dimension(ncolp, nrowp, 0:nbotmp), intent(in) :: botmp ! BOTM of parent grid
    integer, dimension(nlayp) :: lbotmp  ! LBOTM of parent grid
    double precision :: vdist
    ! local
    integer :: kadj
    double precision :: top, bot, topadj, botadj, z, zadj
    !
    ! Top, bottom, and Z coordinate of node of cell of interest
    top = botmp(jp,ip,lbotmp(kp)-1)
    bot = botmp(jp,ip,lbotmp(kp))
    z = (top + bot) / DTWO
    ! Top, bottom, and Z coordinate of node of adjacent cell
    kadj = kp + ioff
    topadj = botmp(jp,ip,lbotmp(kadj)-1)
    botadj = botmp(jp,ip,lbotmp(kadj))
    zadj = (topadj + botadj) / DTWO
    ! Vertical distance between nodes
    vdist = abs(z - zadj)
    !
    return
  end function GetInternodeVertDistance
    
  function GetIoffset(mc,ncp) result(ioff)
    ! mc is position of child cell within parent cell among ncp children.
    ! If mc is at midpoint, return 0.
    ! If mc < midpoint, return -1.
    ! If mc > midpoint, return +1.
    implicit none
    ! dummy
    integer, intent(in) :: mc, ncp
    integer :: ioff
    ! local
    integer :: mid
    double precision :: rmc, rncp
    logical :: odd
    !
    ! Check if mc is at parent-cell midpoint (coincident with node)
    odd = modulo(ncp,2) == 1
    if (odd) then
      mid = (ncp + 1) / 2
      if (mc == mid) then
        ioff = 0
        return
      endif
    endif
    !
    rmc = real(mc,kind(DZERO))
    rncp = real(ncp,kind(DZERO))
    if ((rmc/rncp) > 0.5) then
      ioff = 1
    else
      ioff = -1
    endif
    !
    return
  end function GetIoffset
  
  subroutine CalcContribFactorsTriangular(d1c, d1p, d2c, d2p, cf1, cf2)
    ! Calculate contributing factors using barycentric interpolation, 2D case.
    ! Method:
    ! https://classes.soe.ucsc.edu/cmps160/Fall10/resources/barycentricInterpolation.pdf
    !
    ! This routine can be used for 1D case by providing:
    !       d2c = 0; results in: a2 = 0, contributing factor returned in cf2 = 0 
    !                (d2p value is arbitrary, non-zero). 
    ! OR:
    !       d1c = 0; results in: a1 = 0, contributing factor returned in cf1 = 0
    !                (d1p value is arbitrary, non-zero).
    implicit none
    ! dummy
    double precision, intent(in) :: d1c, d1p, d2c, d2p
    ! Where: d1c = Distance from parent node (in cell where ghost node is
    !              needed) to child node in first dimension.
    !        d1p = Distance from parent node to contributing node in first
    !              dimension.
    !        d2c, d2p = as above, but for second dimension.
    double precision, intent(inout) :: cf1   ! contributing factor for contributing
    !                                          cell adjacent in first dimension.
    double precision, intent(inout) :: cf2   ! as above, but for second dimension
    ! local
    !
    ! Barycentric calculations are:
    !    a = DHALF * d1p * d2p   ! area of enclosing triangle
    !    a1 = DHALF * d2p * d1c  ! area of subtriangle 1
    !    a2 = DHALF * d1p * d2c  ! area of subtriangle 2
    !    cf1 = a1 / a     ! contributing factor for cell 1
    !    cf2 = a2 / a     ! contributing factor for cell 2
    !
    ! These calculations are simplified to:
    cf1 = d1c / d1p
    cf2 = d2c / d2p
    ! 
    return
  end subroutine CalcContribFactorsTriangular

  subroutine CalcContribFactorsBilinear(d1c, d1p, d2c, d2p, cf1, cf2, cf12)
    ! Calculate contributing factors using bilinear interpolation.
    ! Method: 
    ! http://bmia.bmt.tue.nl/people/BRomeny/Courses/8C080/Interpolation.pdf
    !
    implicit none
    ! dummy
    double precision, intent(in) :: d1c, d1p, d2c, d2p
    ! Where: d1c = Distance from parent node (in cell where ghost node is
    !              needed) to child node in first dimension.
    !        d1p = Distance from parent node to contributing node in first
    !              dimension.
    !        d2c, d2p = as above, but for second dimension.
    !
    ! CF1 is contributing factor for contributing cell
    !     adjacent in dimension 1.
    double precision, intent(inout) :: cf1   
    ! CF2 is contributing factor for contributing cell
    !     adjacent in dimension 2.
    double precision, intent(inout) :: cf2
    ! CF12 is contributing factor for contributing cell
    !      off corner in dimensions 1 and 2.
    double precision, intent(inout) :: cf12
    ! local
    double precision :: ap  ! Area of rectangle formed by parent nodes
    !
    cf1 = DZERO
    cf2 = DZERO
    cf12 = DZERO
    ap = d1p * d2p
    !
    ! Allow d1p or d2p to equal zero if not needed.
    if (d1c > DZERO) then
      if (ap > DZERO) then
        cf1 = d1c * (d2p - d2c) / ap
      else
        ! Use linear interpolation in dimension 1.
        cf1 = d1c / d1p
      endif
    endif
    if (d2c > DZERO) then
      if (ap > DZERO) then
        cf2 = d2c * (d1p - d1c) / ap
      else
        ! Use linear interpolation in dimension 2.
        cf2 = d2c / d2p
      endif
    endif
    if (d1c > DZERO .and. d2c > DZERO) then
      cf12 = d1c * d2c / ap
    endif
    ! 
    return
  end subroutine CalcContribFactorsBilinear

  subroutine ModifyIdomainForChild(this, nc, nr, nl, idomain)
    implicit none
    ! dummy
    class(ExchangeType), intent(inout) :: this
    integer, intent(in) :: nc, nr, nl
    integer, dimension(nc,nr,nl), intent(inout) :: idomain
    ! local
    integer :: i, j, k
    !
    ! Ensure index locations where child goes into parent grid are for the child
    call this%Model2%PointToGrid()
    do k=NPLBEG, NPLEND
      do i=NPRBEG, NPREND
        do j=NPCBEG, NPCEND
          idomain(j,i,k) = 0
        enddo
      enddo
    enddo
    !
    return
  end subroutine ModifyIdomainForChild

end module ExchangeModule
