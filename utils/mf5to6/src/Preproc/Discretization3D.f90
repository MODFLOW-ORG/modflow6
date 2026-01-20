module DnmDis3dModule
  use ArrayReadersMF5Module,        only: ReadArray
  use BlockParserModule,         only: BlockParserType
  use ConstantsModule,           only: DTWO, DZERO, LINELENGTH, MAXCHARLEN
  use ConstantsPHMFModule,       only: PI
  use DnmDisBaseModule,          only: DisBaseType
  use GlobalVariablesPHMFModule, only: verbose
  use GeomUtilModule,            only: get_ijk, get_node
  use InputOutputModule,         only: URWORD
  use SimPHMFModule,                 only: count_errors, store_error, &
                                       store_error_unit, ustop
  implicit none
  private
  public dis3d_cr, Dis3dType, CastAsDis3dType

  type, extends(DisBaseType) :: Dis3dType
    integer, pointer                            :: nlay  => null()        !number of layers
    integer, pointer                            :: nrow   => null()       !number of rows
    integer, pointer                            :: ncol   => null()       !number of columns
    integer, dimension(:), pointer, contiguous  :: nodereduced  => null() !(size:nodesuser) given user nodenumber provide the reduced nodenumber (size 0 if not reduced)
    integer, dimension(:), pointer, contiguous  :: nodeuser  => null()    !(size:nodes) given a reduced nodenumber, provide the user nodenumber (size 0 if not reduced)
    ! for PreHeadsMF
    double precision                            :: gridXmax = DZERO
    double precision                            :: gridYmin = DZERO
    double precision, dimension(:), allocatable :: delr
    double precision, dimension(:), allocatable :: delc
    integer,      dimension(:,:,:), allocatable :: idomain
  contains
    procedure :: dis_df => dis3d_df
    procedure :: dis_da => dis3d_da
    ! -- helper functions
    procedure :: get_nodenumber_idx1
    procedure :: get_nodenumber_idx3
    procedure :: get_nodeuser
    procedure :: nodeu_to_string
    procedure :: noder_from_string
    procedure :: supports_layers
    procedure :: get_cell
    procedure :: get_node_coords_idx2
    ! -- private
    procedure :: read_options
    procedure :: read_dimensions
    procedure :: read_data
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: transform_geo_to_grid
  end type Dis3dType

  contains

  subroutine dis3d_cr(dis, name_model, inunit, iout)
! ******************************************************************************
! dis3d_cr -- Create a new discretization 3d object
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    class(DisBaseType), pointer :: dis
    character(len=*), intent(in) :: name_model
    integer, intent(in) :: inunit
    integer, intent(in) :: iout
    type(Dis3dType), pointer :: disnew
! ------------------------------------------------------------------------------
    allocate(disnew)
    dis => disnew
    call disnew%allocate_scalars(name_model)
    dis%inunit = inunit
    dis%iout = iout
    !
    ! Initialize block parser
    call dis%parser%Initialize(inunit, iout)
    !
    ! -- Return
    return
  end subroutine dis3d_cr

  subroutine dis3d_df(this)
! ******************************************************************************
! read_from_file -- Allocate and read discretization information
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(Dis3dType) :: this
    ! -- locals
! ------------------------------------------------------------------------------
    !
    ! -- Identify
    write(this%iout,1) this%inunit
  1 format(1X,/1X,'DIS -- STRUCTURED GRID DISCRETIZATION PACKAGE,',            &
                  ' VERSION 2 : 3/27/2014 - INPUT READ FROM UNIT ',I0,//)
    !
    ! -- Read options
    call this%read_options()
    !
    ! -- Read dimensions block
    call this%read_dimensions()
    !
    ! -- Read GRIDDATA block
    call this%read_data()
    !
    ! -- Return
    return
  end subroutine dis3d_df

  subroutine dis3d_da(this)
! ******************************************************************************
! dis3d_da -- Deallocate discretization data
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(Dis3dType) :: this
    ! -- locals
! ------------------------------------------------------------------------------
    !
    ! -- DisBaseType deallocate
    call this%DisBaseType%dis_da()
    !
    ! -- Deallocate scalars
    call mem_deallocate(this%nlay)
    call mem_deallocate(this%nrow)
    call mem_deallocate(this%ncol)
    !
    ! -- Deallocate Arrays
    call mem_deallocate(this%nodereduced)
    call mem_deallocate(this%nodeuser)
    !
    ! -- Return
    return
  end subroutine dis3d_da

  subroutine read_options(this)
! ******************************************************************************
! read_options -- Read options
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use ConstantsModule, only: LINELENGTH
    !dummy
    class(Dis3dType) :: this
    !locals
    character(len=LINELENGTH) :: errmsg, keyword
    integer :: lloc, istart, istop
    integer :: ierr, ival
    double precision :: rval
    logical :: isfound, endOfBlock
    integer :: iuext
! ------------------------------------------------------------------------------
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, supportOpenClose=.true.)
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
          case('METERS')
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
          write(errmsg,'(4x,a,a)')'****ERROR. UNKNOWN DIS OPTION: ',       &
                                   trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
          call ustop()
        end select
      end do
      write(this%iout,'(1x,a)')'END OF DISCRETIZATION OPTIONS'
    else
      write(this%iout,'(1x,a)')'NO OPTION BLOCK DETECTED.'
    end if
    if(this%lenuni==0) write(this%iout,'(1x,a)') 'MODEL LENGTH UNIT IS UNDEFINED'
    !
    ! -- Return
    return
  end subroutine read_options

  subroutine read_dimensions(this)
! ******************************************************************************
! read_dimensions -- Read dimensions
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use ConstantsModule,  only: LINELENGTH
    !dummy
    class(Dis3dType) :: this
    !locals
    character(len=LINELENGTH) :: errmsg, keyword
    integer :: ierr
    logical :: isfound, endOfBlock
! ------------------------------------------------------------------------------
    !
    ! -- get dimensions block
    call this%parser%GetBlock('DIMENSIONS', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse dimensions block if detected
    if (isfound) then
      write(this%iout,'(1x,a)')'PROCESSING DISCRETIZATION DIMENSIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('NLAY')
          this%nlay = this%parser%GetInteger()
          write(this%iout,'(4x,a,i7)')'NLAY = ', this%nlay
        case ('NROW')
          this%nrow = this%parser%GetInteger()
          write(this%iout,'(4x,a,i7)')'NROW = ', this%nrow
        case ('NCOL')
          this%ncol = this%parser%GetInteger()
          write(this%iout,'(4x,a,i7)')'NCOL = ', this%ncol
        case default
          write(errmsg,'(4x,a,a)')'****ERROR. UNKNOWN DIS DIMENSION: ',      &
                                    trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
          call ustop()
        end select
      end do
      write(this%iout,'(1x,a)')'END OF DISCRETIZATION DIMENSIONS'
    else
      call store_error('ERROR.  REQUIRED DIMENSIONS BLOCK NOT FOUND.')
      call this%parser%StoreErrorUnit()
      call ustop()
    end if
    !
    ! -- verify dimensions were set
    if(this%nlay < 1) then
      call store_error( &
          'ERROR.  NLAY WAS NOT SPECIFIED OR WAS SPECIFIED INCORRECTLY.')
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    if(this%nrow < 1) then
      call store_error( &
          'ERROR.  NROW WAS NOT SPECIFIED OR WAS SPECIFIED INCORRECTLY.')
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    if(this%ncol < 1) then
      call store_error( &
          'ERROR.  NCOL WAS NOT SPECIFIED OR WAS SPECIFIED INCORRECTLY.')
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    !
    ! -- Return
    return
  end subroutine read_dimensions

  subroutine read_data(this)
! ******************************************************************************
! read_data -- Read data
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use SimPHMFModule, only: ustop, count_errors, store_error
    use ConstantsModule,   only: LINELENGTH, DZERO
    ! -- dummy
    class(Dis3dType) :: this
    ! -- locals
!    double precision, dimension(:),     allocatable  :: delr
!    double precision, dimension(:),     allocatable  :: delc
    double precision, dimension(:,:,:), allocatable  :: botm
!    integer,          dimension(:,:,:), allocatable  :: idomain
    character(len=LINELENGTH) :: keyword
    integer :: lloc, istart, istop, n, node, noder, i, j, k
    integer :: ierr, ival
    integer :: nrsize
    logical :: isfound, endOfBlock
    double precision :: rval, dz
    integer, parameter :: nname = 5
    logical, dimension(nname) :: lname
    character(len=24),dimension(nname) :: aname
    character(len=300) :: ermsg
    integer :: iuext
    ! -- formats
    character(len=*), parameter :: fmtdz = &
      "('ERROR. CELL (',i0,',',i0,',',i0,') THICKNESS <= 0. ', " //             &
      "'TOP, BOT: ',2(1pg15.6))"
    character(len=*), parameter :: fmtnr = &
      "(/1x, 'THE SPECIFIED IDOMAIN RESULTS IN A REDUCED NUMBER OF CELLS.'," // &
      "/1x, 'NUMBER OF USER NODES: ',I7," // &
      "/1X, 'NUMBER OF NODES IN SOLUTION: ', I7, //)"
    ! -- data
    data aname(1) /'                    DELR'/
    data aname(2) /'                    DELC'/
    data aname(3) /'TOP ELEVATION OF LAYER 1'/
    data aname(4) /'  MODEL LAYER BOTTOM EL.'/
    data aname(5) /'                 IDOMAIN'/
! ------------------------------------------------------------------------------
    !
    ! -- Allocate temporary arrays used in this subroutine
    allocate(this%delr(this%ncol))
    allocate(this%delc(this%nrow))
    allocate(this%idomain(this%ncol, this%nrow, this%nlay))
    allocate(botm(this%ncol, this%nrow, 0:this%nlay))
    this%nodesuser = this%nlay * this%nrow * this%ncol
    !
    ! --Read GRIDDATA block
    call this%parser%GetBlock('GRIDDATA', isfound, ierr, supportOpenClose=.true.)
    lname(:) = .false.
    if(isfound) then
      write(this%iout,'(1x,a)')'PROCESSING GRIDDATA'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('DELR')
          call ReadArray(this%delr, aname(1), this%ncol, this%inunit, this%iout)
          lname(1) = .true.
        case ('DELC')
          call ReadArray(this%delc, aname(2), this%nrow, this%inunit, this%iout)
          lname(2) = .true.
        case ('TOP')
          call ReadArray(botm(:,:,0), aname(3), this%nrow, this%ncol, 0,        &
                      this%inunit, this%iout)
          lname(3) = .true.
        case ('BOTM')
          call this%parser%GetStringCaps(keyword)
          if (keyword.EQ.'LAYERED') then
            call ReadArray(this%inunit, this%iout, this%nlay, this%nrow,        &
                           this%ncol, this%nodesuser, botm(:,:,1:this%nlay),    &
                           aname(4))
          else
            call ReadArray(botm(:, :, 1:this%nlay), aname(4), this%nodesuser,   &
                           this%inunit, this%iout)
          end if
          lname(4) = .true.
        case ('IDOMAIN')
          call this%parser%GetStringCaps(keyword)
          if (keyword.EQ.'LAYERED') then
            call ReadArray(this%inunit, this%iout, this%nlay, this%nrow,        &
                           this%ncol, this%nodesuser, this%idomain(:, :, :),         &
                           aname(5))
          else
            call ReadArray(this%idomain(:, :, :), aname(5), this%nodesuser,          &
                        this%inunit, this%iout)
          end if
          lname(5) = .true.
        case default
          write(ermsg,'(4x,a,a)')'ERROR. UNKNOWN GRIDDATA TAG: ',            &
                                   trim(keyword)
          call store_error(ermsg)
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
    ! -- Verify all required items were read (IDOMAIN not required)
    do n = 1, nname - 1
      if(.not. lname(n)) then
        write(ermsg,'(1x,a,a)') &
          'ERROR.  REQUIRED INPUT WAS NOT SPECIFIED: ',aname(n)
        call store_error(ermsg)
      endif
    enddo
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    !
    ! -- If IDOMAIN was not read, then set all values to 1, otherwise
    !    count active cells
    if(.not. lname(5)) then
      do k = 1, this%nlay
        do i = 1, this%nrow
          do j = 1, this%ncol
            this%idomain(j, i, k) = 1
          enddo
        enddo
      enddo
      this%nodes = this%nodesuser
    else
      this%nodes = 0
      do k = 1, this%nlay
        do i = 1, this%nrow
          do j = 1, this%ncol
            if(this%idomain(j, i, k) > 0) this%nodes = this%nodes + 1
          enddo
        enddo
      enddo
    endif
    !
    ! -- Check cell thicknesses
    do k = 1, this%nlay
      do i = 1, this%nrow
        do j = 1, this%ncol
          if (this%idomain(j, i, k) == 0) cycle
          if (this%idomain(j, i, k) > 0) then
            dz = botm(j, i, k - 1) - botm(j, i, k)
            if (dz <= DZERO) then
              write(ermsg, fmt=fmtdz) k, i, j, botm(j, i, k - 1), botm(j, i, k)
              call store_error(ermsg)
            endif
          endif
        enddo
      enddo
    enddo
    if (count_errors() > 0) call ustop()
    !
    ! -- Write message if reduced grid
    if(this%nodes < this%nodesuser) then
      write(this%iout, fmtnr) this%nodesuser, this%nodes
    endif
    !
    ! -- Array size is now known, so allocate
    call this%allocate_arrays()
    !
    ! -- Fill the nodereduced array with the reduced nodenumber, or
    !    a negative number to indicate it is a pass-through cell, or
    !    a zero to indicate that the cell is excluced from the
    !    solution.
    if(this%nodes < this%nodesuser) then
      node = 1
      noder = 1
      do k = 1, this%nlay
        do i = 1, this%nrow
          do j = 1, this%ncol
            if(this%idomain(j, i, k) > 0) then
              this%nodereduced(node) = noder
              noder = noder + 1
            elseif(this%idomain(j, i, k) < 0) then
              this%nodereduced(node) = -1
            else
              this%nodereduced(node) = 0
            endif
            node = node + 1
          enddo
        enddo
      enddo
    endif
    !
    ! -- allocate and fill nodeuser if a reduced grid
    if(this%nodes < this%nodesuser) then
      node = 1
      noder = 1
      do k = 1, this%nlay
        do i = 1, this%nrow
          do j = 1, this%ncol
            if(this%idomain(j, i, k) > 0) then
              this%nodeuser(noder) = node
              noder = noder + 1
            endif
            node = node + 1
          enddo
        enddo
      enddo
    endif
    !
    ! -- Move botm into top and bot, and calculate area
    node = 0
    do k=1,this%nlay
      do i=1,this%nrow
        do j=1,this%ncol
          node = node + 1
          noder = node
          if(this%nodes < this%nodesuser) noder = this%nodereduced(node)
          if(noder <= 0) cycle
          this%top(noder) = botm(j, i, k - 1)
          this%bot(noder) = botm(j, i, k)
          this%area(noder) = this%delr(j) * this%delc(i)
        enddo
      enddo
    enddo
    !
    ! -- create and fill the connections object
    nrsize = 0
    if(this%nodes < this%nodesuser) nrsize = this%nodes
!    allocate(this%con)
!    call this%con%structuredconnections(this%name_model, this%nodes,           &
!                                        this%ncol, this%nrow, this%nlay,       &
!                                        nrsize, delr, delc, this%top,          &
!                                        this%bot, this%nodereduced)
!    this%nja = this%con%nja
!    this%njas = this%con%njas
    !
    ! -- deallocate local arrays
!    deallocate(idomain)
    deallocate(botm)
!    deallocate(delr)
!    deallocate(delc)
    !
    ! Compute gridXmax and gridYmin
    this%gridXmax = DZERO
    do j=1,this%ncol
      this%gridXmax = this%gridXmax + this%delr(j)
    enddo
    ! This is grid-coordinate space, so Y becomes more negative
    ! as row index increases.
    this%gridYmin = DZERO
    do i=1,this%nrow
      this%gridYmin = this%gridYmin - this%delc(i)
    enddo
    write(this%iout,*)
    write(this%iout,*)'GridXmax = ',this%gridXmax
    write(this%iout,*)'GridYmin = ',this%gridYmin
    !
    ! -- Return
    return
  end subroutine read_data

  subroutine nodeu_to_string(this, nodeu, str)
! ******************************************************************************
! nodeu_to_string -- Convert user node number to a string in the form of
! (nodenumber) or (k,i,j)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(Dis3dType) :: this
    integer, intent(in) :: nodeu
    character(len=*), intent(inout) :: str
    !local
    integer :: i, j, k
    character(len=10) :: kstr, istr, jstr
! ------------------------------------------------------------------------------
    !
    call get_ijk(nodeu, this%nrow, this%ncol, this%nlay, i, j, k)
    write(kstr, '(i10)') k
    write(istr, '(i10)') i
    write(jstr, '(i10)') j
    str = '(' // trim(adjustl(kstr)) // ',' // &
                 trim(adjustl(istr)) // ',' // &
                 trim(adjustl(jstr)) // ')'
    !
    ! -- return
    return
  end subroutine nodeu_to_string

  integer function get_nodenumber_idx1(this, nodeu, icheck) result(nodenumber)
! ******************************************************************************
! get_nodenumber -- Return a nodenumber from the user specified node number
!                   with an option to perform a check.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    ! -- dummy
    class(Dis3dType), intent(in) :: this
    integer, intent(in) :: nodeu
    integer, intent(in) :: icheck
    ! -- local
    character(len=LINELENGTH) :: errmsg
! ------------------------------------------------------------------------------
    !
    ! -- check the node number if requested
    if(icheck /= 0) then
      !
      ! -- If within valid range, convert to reduced nodenumber
      if(nodeu < 1 .or. nodeu > this%nodesuser) then
        write(errmsg, '(a,i10)') &
          'Nodenumber less than 1 or greater than nodes:', nodeu
        call store_error(errmsg)
        nodenumber = 0
      else
        nodenumber = nodeu
        if(this%nodes < this%nodesuser) nodenumber = this%nodereduced(nodeu)
      endif
    else
      nodenumber = nodeu
      if(this%nodes < this%nodesuser) nodenumber = this%nodereduced(nodeu)
    endif
    !
    ! -- return
    return
  end function get_nodenumber_idx1

  integer function get_nodenumber_idx3(this, k, i, j, icheck)                  &
    result(nodenumber)
! ******************************************************************************
! get_nodenumber_idx3 -- Return a nodenumber from the user specified layer, row,
!                        and column with an option to perform a check.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use ConstantsModule, only: LINELENGTH
    implicit none
    ! dummy
    class(Dis3dType), intent(in) :: this
    integer, intent(in) :: k, i, j
    integer, intent(in) :: icheck
    ! local
    character(len=LINELENGTH) :: errmsg
    integer :: nodeu
    ! formats
    character(len=*), parameter :: fmterr = &
    "('Error in structured-grid cell indices: layer = ',i0,', row = ',i0,', column = ',i0)"
! ------------------------------------------------------------------------------
    !
    nodeu = get_node(k, i, j, this%nlay, this%nrow, this%ncol)
    if (nodeu < 1) then
      write(errmsg, fmterr) k, i, j
      call store_error(errmsg)
      call ustop()
    endif
    nodenumber = nodeu
    if(this%nodes < this%nodesuser) nodenumber = this%nodereduced(nodeu)
    !
    ! -- check the node number if requested
    if(icheck /= 0) then
      !
      if(k < 1 .or. k > this%nlay) &
        call store_error('Layer less than one or greater than nlay')
      if(i < 1 .or. i > this%nrow) &
        call store_error('Row less than one or greater than nrow')
      if(j < 1 .or. j > this%ncol) &
        call store_error('Column less than one or greater than ncol')
      !
      ! -- Error if outside of range
      if(nodeu < 1 .or. nodeu > this%nodesuser) then
        write(errmsg, '(a,i10)') &
          'Nodenumber less than 1 or greater than nodes:', nodeu
        call store_error(errmsg)
      endif
    endif
    !
    ! -- return
    return
  end function get_nodenumber_idx3

  integer function get_nodeuser(this, noder) &
    result(nodenumber)
! ******************************************************************************
! get_nodeuser -- Return the user nodenumber from the reduced node number
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    class(Dis3dType) :: this
    integer, intent(in) :: noder
! ------------------------------------------------------------------------------
    !
    if(this%nodes < this%nodesuser) then
      nodenumber = this%nodeuser(noder)
    else
      nodenumber = noder
    endif
    !
    ! -- return
    return
  end function get_nodeuser

  subroutine allocate_scalars(this, name_model)
! ******************************************************************************
! allocate_scalars -- Allocate and initialize scalars
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(Dis3dType) :: this
    character(len=*), intent(in) :: name_model
! ------------------------------------------------------------------------------
    !
    ! -- Allocate parent scalars
    call this%DisBaseType%allocate_scalars(name_model)
    !
    ! -- Allocate
    call mem_allocate(this%nlay, 'NLAY', this%memoryPath)
    call mem_allocate(this%nrow, 'NROW', this%memoryPath)
    call mem_allocate(this%ncol, 'NCOL', this%memoryPath)
    !
    ! -- Initialize
    this%nlay = 0
    this%nrow = 0
    this%ncol = 0
    this%ndim = 3
    !
    ! -- Return
    return
  end subroutine allocate_scalars

  subroutine allocate_arrays(this)
! ******************************************************************************
! allocate_arrays -- Allocate arrays
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(Dis3dType) :: this
! ------------------------------------------------------------------------------
    !
    ! -- Allocate arrays in DisBaseType (mshape, top, bot, area)
    call this%DisBaseType%allocate_arrays()
    !
    ! -- Allocate arrays for Dis3dType
    if(this%nodes < this%nodesuser) then
      call mem_allocate(this%nodeuser, this%nodes, 'NODEUSER', this%memoryPath)
      call mem_allocate(this%nodereduced, this%nodesuser, 'NODEREDUCED',       &
                        this%memoryPath)
    else
      call mem_allocate(this%nodeuser, 1, 'NODEUSER', this%memoryPath)
      call mem_allocate(this%nodereduced, 1, 'NODEREDUCED', this%memoryPath)
    endif
    !
    ! -- Initialize
    this%mshape(1) = this%nlay
    this%mshape(2) = this%nrow
    this%mshape(3) = this%ncol
    !
    ! -- Return
    return
  end subroutine allocate_arrays

  function noder_from_string(this, lloc, istart, istop, in, iout, line, &
                             flag_string) result(noder)
! ******************************************************************************
! noder_from_string -- Receive a string and convert the string to a reduced
!   nodenumber.  The model discretization is DIS; read layer, row, and column.
!   If flag_string argument is present and true, the first token in string
!   is allowed to be a string (e.g. boundary name). In this case, if a string
!   is encountered, return value as -2.
! ******************************************************************************
    implicit none
    ! dummy
    class(Dis3dType)                 :: this
    integer,           intent(inout) :: lloc
    integer,           intent(inout) :: istart
    integer,           intent(inout) :: istop
    integer,           intent(in)    :: in
    integer,           intent(in)    :: iout
    character(len=*),  intent(inout) :: line
    logical, optional, intent(in)    :: flag_string
    integer                          :: noder
    ! local
    integer :: nlay, nrow, ncol, nodeu
    integer :: k, i, j, lloclocal, n, istat
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
    nlay = this%mshape(1)
    nrow = this%mshape(2)
    ncol = this%mshape(3)
    !
    call urword(line, lloc, istart, istop, 2, k, r, iout, in)
    call urword(line, lloc, istart, istop, 2, i, r, iout, in)
    call urword(line, lloc, istart, istop, 2, j, r, iout, in)
    if(k < 1 .or. k > nlay) then
        write(ermsg, *) ' Layer number in list is outside of the grid', k
        call store_error(ermsg)
    end if
    if(i < 1 .or. i > nrow) then
        write(ermsg, *) ' Row number in list is outside of the grid', i
        call store_error(ermsg)
    end if
    if(j < 1 .or. j > ncol) then
        write(ermsg, *) ' Column number in list is outside of the grid', j
        call store_error(ermsg)
    end if
    nodeu = get_node(k, i, j, nlay, nrow, ncol)
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
    class(Dis3dType) :: this
    !
    supports_layers = .true.
    return
  end function supports_layers

  subroutine get_cell(this, x, y, irow, jcol, gridX, gridY)
    ! Using DELR, DELC, Xorigin, Yorigin, Angle, X, and Y,
    ! determine row and column indices of cell containing
    ! point (x,y).
    ! dummy
    class(Dis3dType) :: this
    double precision, intent(in)  :: x, y
    integer,          intent(out) :: irow, jcol
    double precision, intent(out) :: gridX, gridY
    ! local
    integer :: i, j
    double precision :: xloc, yloc
    character(len=MAXCHARLEN) :: ermsg
    ! formats
    10 format('Point with X = ',f10.2,' and Y = ',f10.2,' is outside of grid.')
    20 format(/,'(X,Y) = ','(',f10.2,', ',f10.2,')',/, &
              'transforms to (gridX, gridY) = ', &
              '(',f9.2,', ',f9.2,') -- inside grid')
    30 format(/,'(X,Y) = ','(',f10.2,', ',f10.2,')',/, &
              'transforms to (gridX, gridY) = ', &
              '(',f9.2,', ',f9.2,') -- outside grid')
    !
    irow = 0
    jcol = 0
    ! Transform X and Y to gridX and gridY
    call this%transform_geo_to_grid(x, y, gridX, gridY)
    if (gridX < DZERO .or. gridX > this%gridXmax .or. &
        gridY > DZERO .or. gridY < this%gridYmin) then
      write(this%iout,30)x,y,gridX,gridY
      write(*,30)x,y,gridX,gridY
      write(ermsg,10)x,y
      call store_error(ermsg)
      call ustop()
    endif
    if (verbose) then
      write(this%iout,20)x,y,gridX,gridY
!      write(*,20)x,y,gridX,gridY
    endif
    !
    ! Find row and column indices
    xloc = DZERO
    do j=1,this%ncol
      xloc = xloc + this%delr(j)
      if (xloc >= gridX) then
        jcol = j
        exit
      endif
    enddo
    yloc = DZERO
    do i=1,this%nrow
      yloc = yloc - this%delc(i)
      if (yloc <= gridY) then
        irow = i
        exit
      endif
    enddo
    !
    return
  end subroutine get_cell

  subroutine transform_geo_to_grid(this, x, y, gridX, gridY)
    ! Transform geospatial projected coordinates to grid coordinates.
    ! Reference:
    ! http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/node153.html
    ! Grid origin is outside corner of cell at row 1, column 1.
    ! GridX increases as column index increases.
    ! GridY decreases (larger negative numbers) as row index increases.
    ! dummy
    class(Dis3dType) :: this
    double precision, intent(in) :: x, y
    double precision, intent(out) :: gridX, gridY
    ! local
    double precision :: xt, yt, xo, yo
    double precision :: theta ! rotation angle, in radians
    !
    ! Coordinates of grid origin in geospatial coordinates
    xo = this%Xorigin
    yo = this%Yorigin
    !
    ! Translate geospatial coordinates to be relative to grid origin
    xt = x - xo
    yt = y - yo
    !
    ! Rotate point about grid origin through theta
    theta = this%Theta
    gridX = xt * cos(theta) + yt * sin(theta)
    gridY = -xt * sin(theta) + yt * cos(theta)
    !
    ! Convert coordinates if length units differ
    gridX = gridX * this%ConvertFactor
    gridY = gridY * this%ConvertFactor
    !
    return
  end subroutine transform_geo_to_grid

  subroutine get_node_coords_idx2(this, irow, jcol, xnode, ynode)
    ! Find X- and Y-coordinates (in grid space) of a user node.
    ! dummy
    class(Dis3dType), intent(inout) :: this
    integer, intent(in) :: irow, jcol
    double precision, intent(out) :: xnode, ynode
    ! local
    integer :: i, j
    double precision :: halfx, halfy
    !
    ! Compute X coordinate; iterate across columns
    halfx = this%delr(1) / DTWO
    xnode = halfx
    if (jcol > 1) then
      do j=2,jcol
        xnode = xnode + halfx
        halfx = this%delr(j) / DTWO
        xnode = xnode + halfx
      enddo
    endif
    !
    ! Compute Y coordinate; iterate down through rows
    halfy = this%delc(1) / DTWO
    ynode = -halfy
    if (irow > 1) then
      do i=2,irow
        ynode = ynode - halfy
        halfy = this%delc(i) / DTWO
        ynode = ynode - halfy
      enddo
    endif
    !
    return
  end subroutine get_node_coords_idx2

  function CastAsDis3dType(obj) result (res)
    implicit none
    class(*), pointer, intent(inout) :: obj
    type(Dis3dType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    type is (Dis3dType)
      res => obj
    end select
    return
  end function

end module DnmDis3dModule






