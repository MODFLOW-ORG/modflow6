module DnmDisBaseModule
  use ArrayReadersMF5Module, only: ReadArray
  use BlockParserModule,  only: BlockParserType
  use ConstantsModule,    only: LENMODELNAME, LENMEMPATH, LINELENGTH, DONE, DZERO
  use InputOutputModule,  only: URWORD
  use SimPHMFModule,          only: count_errors, store_error, store_error_unit, &
                                ustop
  use GeomUtilModule, only: get_node
  implicit none
  private
  public :: DisBaseType
  public :: dis_cr

  type :: DisBaseType
    character(len=LENMEMPATH)                           :: memoryPath            !< path for allocation in memory manager
    character(len=LENMODELNAME), pointer                :: name_model => null()  !< name of the model
    integer, pointer                                    :: inunit     => null()  !< unit number for input file
    integer, pointer                                    :: iout       => null()  !< unit number for output file
    integer, pointer                                    :: nodes      => null()  !< number of nodes in solution
    integer, pointer                                    :: nodesuser  => null()  !< number of user nodes (same as nodes for unstructured model)
    integer, pointer                                    :: ndim       => null()  !< number of spatial model dimensions (1 for unstructured)
    integer, pointer, dimension(:), contiguous          :: mshape     => null()  !< shape of the model; (nodes) for DisBaseType
    integer, pointer                                    :: nja        => null()  !< number of connections plus number of nodes
    integer, pointer                                    :: njas       => null()  !< (nja-nodes)/2
    integer, pointer                                    :: lenuni     => null()  !< length unit
    integer, pointer                                    :: idsymrd    => null()  !< indicates how symmetric arrays are read
    double precision, dimension(:), pointer, contiguous :: top        => null()  !< (size:nodes) cell top elevation
    double precision, dimension(:), pointer, contiguous :: bot        => null()  !< (size:nodes) cell bottom elevation
    double precision, dimension(:), pointer, contiguous :: area       => null()  !< (size:nodes) cell area, in plan view
    ! for PreHeadsMF
    ! For structured DIS, origin is outside corner of cell at row 1, column 1.
    double precision :: Xorigin = DZERO !< X coordinate of grid origin
    double precision :: Yorigin = DZERO !< Y coordinate of grid origin
    double precision :: Theta = DZERO   !< Rotation of rows CCW from horizontal, in radians
    double precision :: ConvertFactor = DONE
    !
    type(BlockParserType) :: parser
  contains
    procedure :: dis_df
    procedure :: dis_da
    ! -- helper functions
    !
    ! -- get_nodenumber is an overloaded integer function that will always
    !    return the reduced nodenumber.  For all grids, get_nodenumber can
    !    be passed the user nodenumber.  For some other grids, it can also
    !    be passed an index.  For dis3d the index is k, i, j, and for 
    !    disv the index is k, n.
    generic   :: get_nodenumber => get_nodenumber_idx1,                       &
                                   get_nodenumber_idx2,                       &
                                   get_nodenumber_idx3
    procedure :: get_nodenumber_idx1
    procedure :: get_nodenumber_idx2
    procedure :: get_nodenumber_idx3
    procedure :: get_nodeuser
    procedure :: nodeu_to_string
    procedure :: connection_normal
    procedure :: connection_vector
    procedure :: noder_from_string
    procedure :: supports_layers
!    procedure :: get_node_coords
    ! -- private
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: read_options
    procedure :: read_dimensions
    procedure :: read_data
  end type DisBaseType

  contains

  subroutine dis_cr(dis, name_model, inunit, iout)
! ******************************************************************************
! dis_cr -- Create discretization object
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    class(DisBaseType), pointer :: dis
    character(len=*), intent(in) :: name_model
    integer, intent(in) :: inunit
    integer, intent(in) :: iout
    type(DisBaseType), pointer :: disnew
! ------------------------------------------------------------------------------
    !
    ! -- Create a new discretization object
    allocate(disnew)
    dis => disnew
    !
    ! -- Allocate scalars and assign data
    call dis%allocate_scalars(name_model)
    dis%inunit = inunit
    dis%iout = iout
    !
    ! Initialize block parser
    call dis%parser%Initialize(inunit, iout)
    !
    ! -- Return
    return
  end subroutine dis_cr

  subroutine dis_df(this)
! ******************************************************************************
! dis_df -- Read discretization information from file
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(DisBaseType) :: this
! ------------------------------------------------------------------------------
    !
    ! -- Identify
    write(this%iout,1) this%inunit
  1 format(1X,/1X,'DISU -- UNSTRUCTURED GRID DISCRETIZATION PACKAGE,',         &
                  ' VERSION 2 : 3/27/2014 - INPUT READ FROM UNIT ',I0,//)
    !
    call this%read_options()
    call this%read_dimensions()
    call this%allocate_arrays()
    call this%read_data()
    !
    ! -- Create and fill the connections object
!    allocate(this%con)
!    call this%con%read_from_block(this%name_model, this%nodes, this%nja,         &
!                                  this%inunit, this%iout, this%top, this%bot)
!    this%njas = this%con%njas
    !
    ! -- Return
    return
  end subroutine dis_df

  subroutine dis_da(this)
! ******************************************************************************
! dis_da -- Deallocate discretization object
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(DisBaseType) :: this
! ------------------------------------------------------------------------------
    !
    ! -- Strings
    deallocate(this%name_model)
    !
    ! -- Scalars
    call mem_deallocate(this%inunit)
    call mem_deallocate(this%iout)
    call mem_deallocate(this%nodes)
    call mem_deallocate(this%nodesuser)
    call mem_deallocate(this%ndim)
    call mem_deallocate(this%nja)
    call mem_deallocate(this%njas)
    call mem_deallocate(this%lenuni)
    call mem_deallocate(this%idsymrd)
    !
    ! -- Arrays
    call mem_deallocate(this%mshape)
    call mem_deallocate(this%top)
    call mem_deallocate(this%bot)
    call mem_deallocate(this%area)
    !
!    ! -- Connections
!    call this%con%con_da()
!    deallocate(this%con)
    !
    ! -- Return
    return
  end subroutine dis_da

  subroutine nodeu_to_string(this, nodeu, str)
! ******************************************************************************
! noder_to_string -- Convert user node number to a string in the form of
! (nodenumber) or (k,i,j)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(DisBaseType) :: this
    integer, intent(in) :: nodeu
    character(len=*), intent(inout) :: str
    ! -- local
    character(len=10) :: nstr
! ------------------------------------------------------------------------------
    !
    write(nstr, '(i0)') nodeu
    str = '(' // trim(adjustl(nstr)) // ')'
    !
    ! -- return
    return
  end subroutine nodeu_to_string

  integer function get_nodeuser(this, noder) &
    result(nodenumber)
! ******************************************************************************
! get_nodeuser -- Return the user nodenumber from the reduced node number
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    class(DisBaseType) :: this
    integer, intent(in) :: noder
! ------------------------------------------------------------------------------
    !
    nodenumber = noder
    !
    ! -- return
    return
  end function get_nodeuser

  subroutine read_options(this)
! ******************************************************************************
! read_options -- Read discretization options
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use MemoryManagerModule, only: mem_allocate
    use ConstantsModule, only: LINELENGTH
    use InputOutputModule, only: urword
    use ArrayReadersMF5Module, only: ReadArray
    use SimPHMFModule, only: ustop, count_errors, store_error
    implicit none
    class(DisBaseType) :: this
    character(len=LINELENGTH) :: errmsg, keyword
    integer :: ierr, nerr
    logical :: isfound, endOfBlock
! ------------------------------------------------------------------------------
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, supportOpenClose = .true.)
    !
    ! -- set default options
    this%lenuni = 0
    !
    ! -- parse options block if detected
    if (isfound) then
      write(this%iout,'(1x,a)')'PROCESSING DISCRETIZATION OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('LENGTH_UNITS')
          call this%parser%GetStringCaps(keyword)
          select case (keyword)
          case ('FEET')
            this%lenuni = 1
            write(this%iout,'(4x,a)') 'MODEL LENGTH UNIT IS FEET'
          case ('METERS')
            this%lenuni = 2
            write(this%iout,'(4x,a)') 'MODEL LENGTH UNIT IS METERS'
          case ('CENTIMETERS')
            this%lenuni = 3
            write(this%iout,'(4x,a)') 'MODEL LENGTH UNIT IS CENTIMETERS'
          case default
            write(this%iout,'(4x,a)')'UNKNOWN UNIT: ',trim(keyword)
            write(this%iout,'(4x,a)')'SETTING TO: ','UNDEFINED'
          end select
        case default
          write(errmsg,'(4x,a,a)')'****ERROR. UNKNOWN DIS OPTION: ',      &
                                   trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write(this%iout,'(1x,a)')'END OF DISCRETIZATION OPTIONS'
    else
      write(this%iout,'(1x,a)')'NO OPTION BLOCK DETECTED.'
    end if
    if(this%lenuni==0) write(this%iout,'(1x,a)') 'MODEL LENGTH UNIT IS UNDEFINED'
    !
    nerr = count_errors()
    if(nerr > 0) then
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    !
    ! -- Return
    return
  end subroutine read_options
  
  subroutine read_dimensions(this)
! ******************************************************************************
! read_dimensions -- Read discretization information from file
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use MemoryManagerModule, only: mem_allocate
    use ConstantsModule, only: LINELENGTH
    use InputOutputModule, only: urword
    use SimPHMFModule, only: ustop, count_errors, store_error
    implicit none
    class(DisBaseType) :: this
    character(len=LINELENGTH) :: errmsg, keyword
    integer :: ierr
    logical :: isfound, endOfBlock
! ------------------------------------------------------------------------------
    !
    ! -- Initialize dimensions
    this%nodes = -1
    this%nja = -1
    !
    ! -- get options block
    call this%parser%GetBlock('DIMENSIONS', isfound, ierr, supportOpenClose=.true.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write(this%iout,'(1x,a)')'PROCESSING DISCRETIZATION DIMENSIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('NODES')
          this%nodes = this%parser%GetInteger()
          write(this%iout,'(4x,a,i7)')'NODES = ', this%nodes
        case ('NJA')
          this%nja = this%parser%GetInteger()
          write(this%iout,'(4x,a,i7)')'NJA   = ', this%nja
        case default
          write(errmsg,'(4x,a,a)')'****ERROR. UNKNOWN DIS DIMENSION: ',      &
                                    trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
          call ustop()
        end select
      end do
      write(this%iout,'(1x,a)')'END OF DISCRETIZATION OPTIONS'
    else
      call store_error('ERROR.  REQUIRED DIMENSIONS BLOCK NOT FOUND.')
      call this%parser%StoreErrorUnit()
      call ustop()
    end if
    !
    ! -- Set nodesuser to nodes
    this%nodesuser = this%nodes
    !
    ! -- verify dimensions were set
    if(this%nodes < 1) then
      call store_error( &
          'ERROR.  NODES WAS NOT SPECIFIED OR WAS SPECIFIED INCORRECTLY.')
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    if(this%nja < 1) then
      call store_error( &
          'ERROR.  NJA WAS NOT SPECIFIED OR WAS SPECIFIED INCORRECTLY.')
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    !
    ! -- Return
    return
  end subroutine read_dimensions

  subroutine read_data(this)
! ******************************************************************************
! read_data -- Read discretization data
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use ConstantsModule, only: LINELENGTH, DZERO
    use SimPHMFModule, only: ustop, count_errors, store_error
    ! -- dummy
    class(DisBaseType) :: this
    ! -- local
    character(len=LINELENGTH) :: errmsg, keyword
    integer :: n
    integer :: ierr
    logical :: isfound, endOfBlock
    double precision :: dz
    integer, parameter :: nname = 3
    logical,dimension(nname) :: lname
    character(len=24),dimension(nname) :: aname(nname)
    ! -- formats
    character(len=*), parameter :: fmtdz = &
      "('ERROR. CELL ', i0, ' WITH THICKNESS <= 0. TOP, BOT: ', 2(1pg15.6))"
    character(len=*), parameter :: fmtarea = &
      "('ERROR. CELL ', i0, ' WITH AREA <= 0. AREA: ', 1(1pg15.6))"
    ! -- data
    data aname(1) /'                     TOP'/
    data aname(2) /'                     BOT'/
    data aname(3) /'                    AREA'/
! ------------------------------------------------------------------------------
    !
    ! -- get griddata block
    call this%parser%GetBlock('GRIDDATA', isfound, ierr, supportOpenClose=.true.)
    lname(:) = .false.
    if(isfound) then
      write(this%iout,'(1x,a)') 'PROCESSING GRIDDATA'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('TOP')
          call ReadArray(this%top, aname(1), this%nodes, this%inunit, this%iout)
          lname(1) = .true.
        case ('BOT')
          call ReadArray(this%bot, aname(2), this%nodes, this%inunit, this%iout)
          lname(2) = .true.
        case ('AREA')
          call ReadArray(this%area, aname(3), this%nodes, this%inunit, this%iout)
          lname(3) = .true.
        case default
          write(errmsg,'(4x,a,a)')'ERROR. UNKNOWN GRIDDATA TAG: ', &
                                   trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
          call ustop()
        end select
      end do
      write(this%iout,'(1x,a)')'END PROCESSING GRIDDATA'
    else
      call store_error('ERROR.  REQUIRED GRIDDATA BLOCK NOT FOUND.')
      call this%parser%StoreErrorUnit()
      call ustop()
    end if
    !
    ! -- verify all items were read
    do n = 1, nname
      if(.not. lname(n)) then
        write(errmsg,'(1x,a,a)') &
          'ERROR.  REQUIRED INPUT WAS NOT SPECIFIED: ', aname(n)
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      endif
    enddo
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    !
    ! -- Check for zero and negative thickness and zero or negative areas
    do n = 1, this%nodes
      dz = this%top(n) - this%bot(n)
      if (dz <= DZERO) then
        write(errmsg, fmt=fmtdz) n, this%top(n), this%bot(n)
        call store_error(errmsg)
      endif
      if (this%area(n) <= DZERO) then
        write(errmsg, fmt=fmtarea) n, this%area(n)
        call store_error(errmsg)
      endif
    enddo
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    !
    ! -- Return
    return
  end subroutine read_data

  function get_nodenumber_idx1(this, nodeu, icheck) result(nodenumber)
! ******************************************************************************
! get_nodenumber -- Return a nodenumber from the user specified node number
!                   with an option to perform a check.  This subroutine
!                   can be overridden by child classes to perform mapping
!                   to a model node number
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use ConstantsModule, only: LINELENGTH
    use SimPHMFModule, only: store_error
    implicit none
    class(DisBaseType), intent(in) :: this
    integer, intent(in) :: nodeu
    integer, intent(in) :: icheck
    character(len=LINELENGTH) :: errmsg
    integer :: nodenumber
! ------------------------------------------------------------------------------
    !
    if(icheck /= 0) then
      if(nodeu < 1 .or. nodeu > this%nodes) then
        write(errmsg, '(a,i10)') &
          'Nodenumber less than 1 or greater than nodes:', nodeu
        call store_error(errmsg)
      endif
    endif
    !
    ! -- set node number to passed in nodenumber since there is a one to one
    !    mapping for an unstructured grid
    nodenumber = nodeu
    !
    ! -- return
    return
  end function get_nodenumber_idx1

  function get_nodenumber_idx2(this, k, j, icheck) result(nodenumber)
! ******************************************************************************
! get_nodenumber_idx2 -- This function should never be called.  It must be
!   overriden by a child class.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use SimPHMFModule, only: ustop, store_error
    implicit none
    ! Dummy arguments
    class(DisBaseType), intent(in) :: this
    integer, intent(in) :: k, j
    integer, intent(in) :: icheck
    integer :: nodenumber
! ------------------------------------------------------------------------------
    !
    nodenumber = 0
    call store_error('Program error: get_nodenumber_idx2 not implemented.')
    call ustop()
    !
    ! -- Return
    return
  end function get_nodenumber_idx2

  function get_nodenumber_idx3(this, k, i, j, icheck) result(nodenumber)
! ******************************************************************************
! get_nodenumber_idx3 -- This function will not be invoked for an unstructured
! model, but it may be from a Discretization3dType model.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use SimPHMFModule, only: ustop, store_error
    implicit none
    ! Dummy arguments
    class(DisBaseType), intent(in) :: this
    integer, intent(in) :: k, i, j
    integer, intent(in) :: icheck
    integer :: nodenumber
! ------------------------------------------------------------------------------
    !
    nodenumber = 0
    call store_error('Program error: get_nodenumber_idx3 not implemented.')
    call ustop()
    !
    ! -- Return
    return
  end function get_nodenumber_idx3

  subroutine connection_normal(this, noden, nodem, ihc, xcomp, ycomp, zcomp)
! ******************************************************************************
! connection_normal -- calculate the normal vector components for reduced 
!   nodenumber cell (noden) and its shared face with cell nodem.  ihc is the
!   horizontal connection flag.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use SimPHMFModule, only: ustop, store_error
    ! -- dummy
    class(DisBaseType) :: this
    integer, intent(in) :: noden
    integer, intent(in) :: nodem
    integer, intent(in) :: ihc
    double precision, intent(inout) :: xcomp
    double precision, intent(inout) :: ycomp
    double precision, intent(inout) :: zcomp
! ------------------------------------------------------------------------------
    !
    call store_error('Program error: connection_normal not implemented.')
    call ustop()
    !
    ! -- return
    return
  end subroutine connection_normal
    
  subroutine connection_vector(this, noden, nodem, nozee, satn, satm, ihc,   &
                               xcomp, ycomp, zcomp, conlen)                      ! xt3d
! ******************************************************************************
! connection_vector -- calculate the unit vector components from reduced 
!   nodenumber cell (noden) to its neighbor cell (nodem).  The saturation for
!   for these cells are also required so that the vertical position of the cell
!   cell centers can be calculated.  ihc is the horizontal flag.  Also return    ! xt3d
!   the straight-line connection length.                                         ! xt3d
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use SimPHMFModule, only: ustop, store_error
    ! -- dummy
    class(DisBaseType) :: this
    integer, intent(in) :: noden
    integer, intent(in) :: nodem
    logical, intent(in) :: nozee
    double precision, intent(in) :: satn
    double precision, intent(in) :: satm
    integer, intent(in) :: ihc
    double precision, intent(inout) :: xcomp
    double precision, intent(inout) :: ycomp
    double precision, intent(inout) :: zcomp
    double precision, intent(inout) :: conlen
    ! -- local
! ------------------------------------------------------------------------------
    !
    call store_error('Program error: connection_vector not implemented.')
    call ustop()
    !
    ! -- return
    return
  end subroutine connection_vector
    
  subroutine allocate_scalars(this, name_model)
! ******************************************************************************
! allocate_scalars -- Allocate and initialize scalar variables in this class
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use MemoryHelperModule, only: create_mem_path
    ! -- dummy
    class(DisBaseType) :: this
    character(len=*), intent(in) :: name_model
    ! -- local
! ------------------------------------------------------------------------------
    !
    ! -- Assign memory path
    this%memoryPath = create_mem_path(name_model, 'DIS')
    !
    ! -- Allocate
    allocate(this%name_model)
    call mem_allocate(this%inunit, 'INUNIT', this%memoryPath)
    call mem_allocate(this%iout, 'IOUT', this%memoryPath)
    call mem_allocate(this%nodes, 'NODES', this%memoryPath)
    call mem_allocate(this%nodesuser, 'NODESUSER', this%memoryPath)
    call mem_allocate(this%ndim, 'NDIM', this%memoryPath)
    call mem_allocate(this%nja, 'NJA', this%memoryPath)
    call mem_allocate(this%njas, 'NJAS', this%memoryPath)
    call mem_allocate(this%lenuni, 'LENUNI', this%memoryPath)
    call mem_allocate(this%idsymrd, 'IDSYMRD', this%memoryPath)
    !
    ! -- Initialize
    this%name_model = name_model
    this%inunit = 0
    this%iout = 0
    this%nodes = 0
    this%nodesuser = 0
    this%ndim = 1
    this%nja = 0
    this%njas = 0
    this%lenuni = 0
    this%idsymrd = 0
    !
    ! -- Return
    return
  end subroutine allocate_scalars

  subroutine allocate_arrays(this)
! ******************************************************************************
! allocate_arrays -- Read discretization information from file
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(DisBaseType) :: this
! ------------------------------------------------------------------------------
    !
    ! -- Allocate
    call mem_allocate(this%mshape, this%ndim, 'MSHAPE', this%memoryPath)
    call mem_allocate(this%top, this%nodes, 'TOP', this%memoryPath)
    call mem_allocate(this%bot, this%nodes, 'BOT', this%memoryPath)
    call mem_allocate(this%area, this%nodes, 'AREA', this%memoryPath)
    !
    ! -- Initialize
    this%mshape(1) = this%nodes
    !
    ! -- Return
    return
  end subroutine allocate_arrays
  
  function noder_from_string(this, lloc, istart, istop, in, iout, line, &
                             flag_string) result(noder)
! ******************************************************************************
! noder_from_string -- Receive a string and convert the string to a reduced
!   nodenumber.  The model is unstructured; just read user nodenumber.
!   If flag_string argument is present and true, the first token in string 
!   is allowed to be a string (e.g. boundary name). In this case, if a string
!   is encountered, return value as -2.
! ******************************************************************************
    implicit none
    ! dummy
    class(DisBaseType)               :: this
    integer,           intent(inout) :: lloc
    integer,           intent(inout) :: istart
    integer,           intent(inout) :: istop
    integer,           intent(in)    :: in
    integer,           intent(in)    :: iout
    character(len=*),  intent(inout) :: line
    logical, optional, intent(in)    :: flag_string
    integer                          :: noder
    ! local
    integer :: nodes, nodeu
    integer :: lloclocal, n, istat
    double precision :: r
    character(len=LINELENGTH) :: ermsg, nodestr, fname
    !
    if (present(flag_string)) then
      if (flag_string) then
        ! Check to see if first token in line can be read as an integer.
        lloclocal = lloc
        call urword(line, lloclocal, istart, istop, 1, nodeu, r, iout, in)
        read(line(istart:istop),*,iostat=istat)n
        if (istat /= 0) then
          ! First token in line is not an integer; return flag to this effect.
          noder = -2
          return
        endif
      endif
    endif
    !
    nodes = this%nodes
    !
    call urword(line, lloc, istart, istop, 2, nodeu, r, iout, in)
    if(nodeu < 1 .or. nodeu > nodes) then
      write(ermsg, *) ' Node number in list is outside of the grid', nodeu
      call store_error(ermsg)
    end if
    !
    ! -- Check for errors
    if(count_errors() > 0) then
      inquire(unit=in, name=fname)
      call store_error('Error converting in file: ')
      call store_error(trim(adjustl(fname)))
      call store_error('Cell number cannot be determined in line: ')
      call store_error(trim(adjustl(line)))
      call store_error_unit(in)
      call ustop()
    endif
    !
    ! -- Convert user-based nodenumber to reduced node number
    noder = this%get_nodenumber(nodeu, 0)
    if(noder <= 0) then
      call this%nodeu_to_string(nodeu, nodestr)
      write(ermsg, *) &
              ' Cell is outside active grid domain: ' // &
              trim(adjustl(nodestr))
      call store_error(ermsg)
    endif
    !
    ! -- return
    return
    
  end function noder_from_string
  
  logical function supports_layers(this)
    implicit none
    ! dummy
    class(DisBaseType) :: this
    !
    supports_layers = .false.
    return
  end function supports_layers
  
!  subroutine get_node_coords(this, nodeuser, xnode, ynode)
!    ! Find X- and Y-coordinates (in grid space) of a user node.
!    ! Override for classes that extend DisBaseType.
!    ! dummy
!    class(DisBaseType), intent(inout) :: this
!    integer, intent(in) :: nodeuser
!    double precision, intent(out) :: xnode, ynode
!    !
!    xnode = 0.0d0
!    ynode = 0.0d0
!    !
!    return
!  end subroutine get_node_coords
  
end module DnmDisBaseModule
