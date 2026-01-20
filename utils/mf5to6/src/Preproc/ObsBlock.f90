module ObsBlockModule
  
  use BlockParserModule,         only: BlockParserType
  use ConstantsModule,           only: DONE, DZERO, &
                                       LINELENGTH, MAXCHARLEN, LENOBSNAME
  use MathUtilModule,            only: is_close
  use ConstantsPHMFModule,       only: CONTINUOUS, SINGLE, LENOBSNAMENEW
  use DnmDis3dModule,            only: Dis3dType
  use GlobalVariablesPHMFModule, only: verbose
  use InputOutputModule,         only: UPCASE, URWORD
  use ListModule,                only: ListType
  use ObserveModule,             only: ObserveType, AddObserveToList, &
                                       GetObserveFromList, ConstructObservation
  use SimPHMFModule,             only: store_error, store_error_unit, ustop
  use UtilitiesModule,           only: CalcContribFactors
  
  implicit none
  
  private
  public :: ObsBlockType, ConstructObsBlockType, AddObsBlockToList, &
            GetObsBlockFromList
  
  type :: ObsBlockType
    integer                   :: isorc = 0      ! 1=SINGLE, 2=CONTINUOUS
    integer                   :: inunit = 0
    integer                   :: iout = 0
    integer                   :: ioutMFobs = 0
    integer                   :: IoutPostObs = 0
    character(len=10)         :: SorC = ''      ! SINGLE or CONTINUOUS
    character(len=MAXCHARLEN) :: OutputBaseName = ''
    logical                   :: Binary = .false.
    type(ListType), pointer   :: ObsList => null()
    type(Dis3dType), pointer  :: dis3d => null()
  contains
    ! Public procedures
    procedure, public :: AddObs
    procedure, public :: GetObs
    procedure, public :: process_block
    procedure, public :: write_postobs_input
    ! Private procedures
    procedure, private :: calc_contrib_factors
  end type ObsBlockType
  
contains

  ! Procedures bound to ObsBlockType
  
  subroutine AddObs(this, obs)
    ! dummy
    class(ObsBlockType), intent(inout) :: this
    type(ObserveType), pointer, intent(inout) :: obs
    !
    call AddObserveToList(this%ObsList, obs)
    !
    return
  end subroutine AddObs
  
  function GetObs(this, indx) result(res)
    ! dummy
    class(ObsBlockType), intent(inout) :: this
    integer,             intent(in)    :: indx
    type(ObserveType), pointer :: res
    !
    res => GetObserveFromList(this%ObsList, indx)
    !
    return
  end function GetObs

  subroutine process_block(this, insertLine, WriteBeginEnd, parser)
    ! dummy
    class(ObsBlockType) :: this
    logical :: insertLine
    logical, intent(in) :: WriteBeginEnd
    type(BlockParserType), intent(inout) :: parser
    ! local
    integer :: ierr, ioutmf, layer
    integer :: irow, jcol, lloc2, iadjrow, jadjcol, ncol, nrow
    double precision :: time, xcoord, xnode, ycoord, ynode, xoff, yoff
    double precision :: gridX, gridY
    character(len=MAXCHARLEN) :: ermsg
    character(len=LINELENGTH) :: word, outfilename
    logical :: eof, fmtd, ignore_comments, endOfBlock
    type(ObserveType), pointer :: newobs
    type(Dis3dType), pointer :: dis3d => null()
    class(*), pointer :: obj => null()
    ! format
    1 format()
    10 format('Point (',f9.2,', ',f9.2,') is in row ',i0,', column ',i0)
    20 format('For this point, xoff = ',f9.2,',  yoff = ',f9.2)
    30 format(a,1x,a,2x,a,2x,a)
    40 format(a,1x,a,2x,a)
    50 format(a,1x,a)
    70 format(2x,a,2x,'head',2x,g14.7,3(2x,i0))  ! For single obs
    80 format(2x,a,2x,'head',2x,3(2x,i0))        ! For continuous obs
    90 format('BEGIN',1x,a,1x,'FILEOUT',1x,a)
    95 format('END',1x,a)
    !
    dis3d => this%dis3d
    ioutmf  = this%IoutMFobs
    eof = .false.
    ignore_comments = .true.
    ierr = 0
    fmtd = .not. this%Binary
    nrow = dis3d%nrow
    ncol = dis3d%ncol
    !
    if (this%Binary) then
      outfilename = trim(this%OutputBaseName) // '.bsv'
    else
      outfilename = trim(this%OutputBaseName) // '.csv'
    endif
    !
    ! Write BEGIN line
    if (insertLine) write(ioutmf,1)
    if (WriteBeginEnd) then
      write(ioutmf,90)trim(this%SorC),trim(outfilename)
    endif
    !
    ! Read and process each line contained in block
    loop: do
      call parser%GetNextLine(endOfBlock)
      if (endOfBlock) exit
      lloc2 = 1
      if (eof) then
        ! End-of-file encountered. Encountered here, this is an error.
        ermsg = 'End-of-file encountered in ' // trim(this%SorC) // ' block.'
        call store_error(ermsg)
        call parser%StoreErrorUnit()
        call ustop()
      endif
      !
      ! -- get the first word on the line, which is an observation name or "END"
!      call urword(oline, lloc2, istart, istop, 1, idum, rdum, this%iout, iu)
!      word = oline(istart:istop)
      call parser%GetStringCaps(word)
      select case (word)
!      case ('END')
!        call uterminate_block(iu, this%Iout, word, this%SorC, lloc2, &
!                              oline, ierr, iuext)
!        exit loop
      case default
        ! If word is not END, assume it is an observation name. 
        ! Get time, coordinates, and layer.
        if (this%isorc == SINGLE) then
!          call urword(oline, lloc2, istart, istop, 3, idum, time, this%iout, iu)
          time = parser%GetDouble()
        else
          time = -DONE
        endif
!        call urword(oline, lloc2, istart, istop, 3, idum, xcoord, this%iout, iu)
!        call urword(oline, lloc2, istart, istop, 3, idum, ycoord, this%iout, iu)
!        call urword(oline, lloc2, istart, istop, 2, layer, rdum, this%iout, iu)
        xcoord = parser%GetDouble()
        ycoord = parser%GetDouble()
        layer = parser%GetInteger()
        !
        ! Construct an ObserveType object for this observation.
        call ConstructObservation(newobs, this%isorc, time, layer, 0, fmtd)
        newobs%Name = word(1:LENOBSNAME)
        !
        ! Determine row and column indices based on X and Y coordinates.
        call dis3d%get_cell(xcoord, ycoord, irow, jcol, gridX, gridY)
        if (verbose) then
          write (this%Iout,10)xcoord, ycoord, irow, jcol
        endif
        newobs%irow = irow
        newobs%jcol = jcol
        newobs%gridX = gridX
        newobs%gridY = gridY
        !
        ! Get coordinates of node and calculate X and Y offsets from node.
        call dis3d%get_node_coords_idx2(irow, jcol, xnode, ynode)
        xoff = gridX - xnode
        yoff = gridY - ynode
        if (verbose) then
          write(this%Iout,20)xoff, yoff
        endif
        newobs%xoff = xoff
        newobs%yoff = yoff
        !
        ! If xoff or yoff is not zero, may need contribution(s)
        ! from adjacent cell(s). Up to 3 adjacent cells may be needed. 
        ! Determine which adjacent cells are needed, if any.
        iadjrow = 0
        jadjcol = 0
        !
        if (.not. is_close(xoff, DZERO)) then
          if (xoff > DZERO) then
            if (jcol < ncol) then
              if (dis3d%idomain(jcol+1, irow, layer) == 1) then
                jadjcol = jcol + 1
                newobs%jcoladj = jadjcol
              endif
            endif
          else
            if (jcol > 1) then
              if (dis3d%idomain(jcol-1, irow, layer) == 1) then
                jadjcol = jcol - 1
                newobs%jcoladj = jadjcol
              endif
            endif
          endif
        endif
        !
        if (.not. is_close(yoff, DZERO)) then
          if (yoff > DZERO) then
            if (irow > 1) then
              if (dis3d%idomain(jcol, irow-1, layer) == 1) then
                iadjrow = irow - 1
                newobs%irowadj = iadjrow
              endif
            endif
          else
            if (irow < nrow) then
              if (dis3d%idomain(jcol, irow+1, layer) == 1) then
                iadjrow = irow + 1
                newobs%irowadj = iadjrow
              endif
            endif
          endif
        endif
        !
        ! If both jadjcol and iadjrow are not 0, and diagonal cell is active,
        ! assign newobs%idiag flag.
        if (jadjcol /= 0 .and. iadjrow /= 0) then
          if (dis3d%idomain(jadjcol,iadjrow,layer) == 1) then
            newobs%idiag = 1
          endif
        endif
        !
        ! If interpolation is not needed in grid-X or grid-Y direction, set
        ! xoff or yoff back to zero.
        if (jadjcol == 0) newobs%xoff = DZERO
        if (iadjrow == 0) newobs%yoff = DZERO
        !
        ! Get delradj and delcadj as needed.
        if (jadjcol /= 0) then
          newobs%delradj = dis3d%delr(jadjcol)
        endif
        if (iadjrow /= 0) then
          newobs%delcadj = dis3d%delc(iadjrow)
        endif
        !
        ! Add this observation to the list
        obj => newobs
        call this%ObsList%Add(obj)
        !
        ! Write lines to define head observations required for MF6.
        !
        ! Always need obs in cell where X,Y is located; append "-p" to name (primary).
        newobs%NameP = trim(newobs%Name) // '-P'
        newobs%NameR = trim(newobs%Name) // '-R'
        newobs%NameC = trim(newobs%Name) // '-C'
        newobs%NameD = trim(newobs%Name) // '-D'
        select case (this%isorc)
        case (SINGLE)
          write(ioutmf,70)trim(newobs%NameP),time,layer,irow,jcol
        case (CONTINUOUS)
          write(ioutmf,80)trim(newobs%NameP),layer,irow,jcol
        end select
        !
        ! Is obs for cell in adjacent row needed?
        if (newobs%irowadj /= 0) then
          select case (this%isorc)
          case (SINGLE)
            write(ioutmf,70)trim(newobs%NameR),time,layer,newobs%irowadj,jcol
          case (CONTINUOUS)
            write(ioutmf,80)trim(newobs%NameR),layer,newobs%irowadj,jcol
          end select
        endif
        !
        ! Is obs for cell in adjacent column needed?
        if (newobs%jcoladj /= 0) then
          select case (this%isorc)
          case (SINGLE)
            write(ioutmf,70)trim(newobs%NameC),time,layer,irow,newobs%jcoladj
          case (CONTINUOUS)
            write(ioutmf,80)trim(newobs%NameC),layer,irow,newobs%jcoladj
          end select
        endif
        !
        ! Is obs for cell in diagonal position needed?
        if (newobs%idiag == 1) then
          select case (this%isorc)
          case (SINGLE)
            write(ioutmf,70)trim(newobs%NameD),time,layer,newobs%irowadj,newobs%jcoladj
          case (CONTINUOUS)
            write(ioutmf,80)trim(newobs%NameD),layer,newobs%irowadj,newobs%jcoladj
          end select
        endif
      end select
    enddo loop
    !
    ! Write END line
    if (WriteBeginEnd) then
      write(ioutmf,95)trim(this%SorC)
    endif
    !
    return
  end subroutine process_block
  
  subroutine calc_contrib_factors(this)
    ! dummy
    class(ObsBlockType) :: this
    ! local
    integer :: i, ibound12, irow, irowadj, jcol, jcoladj, nobs
    ! Contributing factors for primary cell, adjacent row, 
    ! adjacent column, and diagonal cell.
    double precision :: cfp, cfr, cfc, cfd
    ! Distances, as defined on CalcContribFactors diagram.
    double precision :: d1c, d1p, d2c, d2p
    ! Coordinates of observation point
    double precision :: gridX, gridY
    ! Coordinates of nodes (parent, adjacent row, adjacent column
    double precision :: gridXp, gridXr, gridXc, gridYp, gridYr, gridYc
    type(ObserveType), pointer :: obs => null()
    type(Dis3dType), pointer :: dis3d => null()
    logical :: bilinear = .true.
    !
    dis3d => this%dis3d
    ! 
    nobs = this%ObsList%Count()
    !
    ! Calculate contributing factors for all cells involved in interpolation.
    do i=1,nobs
      d1p = DZERO
      d2p = DZERO
      obs => this%GetObs(i)
      irow = obs%irow
      jcol = obs%jcol
      gridX = obs%gridX
      gridY = obs%gridY
      call dis3d%get_node_coords_idx2(irow, jcol, gridXp, gridYp)
      d1c = abs(gridX - gridXp)
      d2c = abs(gridY - gridYp)
      !
      ! Dimension 1 is along row
      jcoladj = obs%jcoladj
      if (jcoladj /= 0) then
        call dis3d%get_node_coords_idx2(irow, jcoladj, gridXc, gridYc)
        d1p = abs(gridXp - gridXc)
      endif
      !
      ! Dimension 2 is along column
      irowadj = obs%irowadj
      if (irowadj /= 0) then
        call dis3d%get_node_coords_idx2(irowadj, jcol, gridXr, gridYr)
        d2p = abs(gridYp - gridYr)
      endif
      !
      ! Consider diagonal?
      if (irowadj /= 0 .and. jcoladj /= 0) then
        ibound12 = 1
      else
        ibound12 = 0
      endif
      !
      ! Calculate contributing factors
      cfp = DONE   ! for primary cell
      cfr = DZERO  ! for cell in adjacent row
      cfc = DZERO  ! for cell in adjacent column
      cfd = DZERO  ! for diagonally adjacent cell
      call CalcContribFactors(bilinear, d1c, d1p, d2c, d2p, ibound12, cfr, cfc, cfd)
      cfp = DONE - cfr - cfc - cfd
      obs%cfp = cfp
      obs%cfr = cfr
      obs%cfc = cfc
      obs%cfd = cfd
    enddo
    !
    return
  end subroutine calc_contrib_factors

  subroutine write_postobs_input(this, outputfilename)
    ! dummy
    class(ObsBlockType) :: this
    character(len=MAXCHARLEN), intent(out) :: outputfilename
    ! local
    integer :: i, iout, nobs
    type(ObserveType), pointer :: obs => null()
    ! formats
    1 format()
    10 format(a,1x,a,1x,'FILEOUT',1x,a)
    12 format(a,1x,a,1x,'FILEOUT',1x,a,2x,'BINARY')
    15 format(a,1x,a)
    20 format(2x,a,2x,a)
    30 format(4x,a,2x,G14.7)
    !
    ! Calculate contributing factors for components of all observations
    call this%calc_contrib_factors()
    !
    ! Write BEGIN line
    iout = this%IoutPostObs
    write(iout,1)
    if (this%Binary) then
      outputfilename = trim(this%OutputBaseName) // '.post.bsv'
      write(iout,12)'BEGIN', trim(this%SorC), trim(outputfilename)
    else
      outputfilename = trim(this%OutputBaseName) // '.post.csv'
      write(iout,10)'BEGIN', trim(this%SorC), trim(outputfilename)
    endif
    !
    ! Write observation information
    nobs = this%ObsList%Count()
    do i=1,nobs
      obs => GetObserveFromList(this%ObsList, i)
      ! Write NEW line
      write(iout,20)'NEW',trim(obs%Name)
      ! Write source line for primary cell
      write(iout,30)trim(obs%NameP), obs%cfp
      ! Write source line for cell in adjacent column (in row direction)
      if (obs%cfr /= DZERO) then
        write(iout,30)trim(obs%NameC), obs%cfr
      endif
      ! Write source line for cell in adjacent row (in column direction)
      if (obs%cfc /= DZERO) then
        write(iout,30)trim(obs%NameR), obs%cfc
      endif
      ! Write source line for diagonally adjacent cell
      if (obs%cfd /= DZERO) then
        write(iout,30)trim(obs%NameD), obs%cfd
      endif
    enddo
    !
    ! Write END line
    write(iout,15)'END',trim(this%SorC)
    !
    return
  end subroutine write_postobs_input

  ! Non-type-bound procedures

  subroutine ConstructObsBlockType(obsblock, sorc, outputbasename, binary)
    ! dummy
    type(ObsBlockType), pointer, intent(inout) :: obsblock
    character(len=*), intent(inout) :: sorc
    character(len=*), intent(in)    :: outputbasename
    logical, intent(in)             :: binary
    ! local
    character(len=MAXCHARLEN) :: ermsg
    !
    allocate(obsblock)
    call UPCASE(sorc)
    obsblock%SorC = sorc
    select case (sorc)
    case ('SINGLE')
      obsblock%isorc = SINGLE
    case ('CONTINUOUS')
      obsblock%isorc = CONTINUOUS
    case default
      ermsg = 'Expecting SINGLE or CONTINUOUS but found: ' // trim(sorc)
      call store_error(ermsg)
      call ustop()
    end select
    obsblock%OutputBaseName = outputbasename
    obsblock%Binary = binary
    allocate(obsblock%ObsList)
    !
    return
  end subroutine ConstructObsBlockType
  
  subroutine AddObsBlockToList(list, obsblock)
    ! dummy
    type(ListType),              intent(inout) :: list
    type(ObsBlockType), pointer, intent(inout) :: obsblock
    ! local
    class(*), pointer :: obj => null()
    !
    obj => obsblock
    call list%Add(obj)
    !
    return
  end subroutine AddObsBlockToList

  function GetObsBlockFromList(list, indx) result(res)
    ! dummy
    type(ListType), pointer, intent(inout) :: list
    integer, intent(in) :: indx
    type(ObsBlockType), pointer :: res
    ! local
    class(*), pointer :: obj => null()
    !
    obj => list%GetItem(indx)
    select type (obj)
    type is (ObsBlockType)
      res => obj
    end select
    !
    return
  end function GetObsBlockFromList
  
  
end module ObsBlockModule
