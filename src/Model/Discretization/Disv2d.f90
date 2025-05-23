module Disv2dModule

  use ArrayReadersModule, only: ReadArray
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, LENVARNAME, &
                             DZERO, DONE, DHALF
  use BaseDisModule, only: DisBaseType
  use GeomUtilModule, only: get_node, get_ijk, get_jk
  use InputOutputModule, only: URWORD, ulasav, &
                               ulaprufw, ubdsv1, ubdsv06, getunit, openfile
  use SimModule, only: count_errors, store_error, store_error_unit, &
                       store_error_filename
  use SimVariablesModule, only: errmsg, idm_context
  use DisvGeom, only: DisvGeomType, line_unit_vector
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, mem_setptr
  use MemoryManagerExtModule, only: mem_set_value, memorystore_remove
  use TdisModule, only: kstp, kper, pertim, totim, delt

  implicit none
  private
  public disv2d_cr, Disv2dType

  !> @brief Vertex grid discretization
  type, extends(DisBaseType) :: Disv2dType
    integer(I4B), pointer :: nvert => null() !< number of x,y vertices
    real(DP), dimension(:, :), pointer, contiguous :: vertices => null() !< cell vertices stored as 2d array of x and y
    real(DP), dimension(:, :), pointer, contiguous :: cellxy => null() !< cell center stored as 2d array of x and y
    integer(I4B), dimension(:), pointer, contiguous :: iavert => null() !< cell vertex pointer ia array
    integer(I4B), dimension(:), pointer, contiguous :: javert => null() !< cell vertex pointer ja array
    real(DP), dimension(:), pointer, contiguous :: bottom => null() !< bottom elevations for each cell (nodes)
  contains

    procedure :: dis_df => disv2d_df
    procedure :: dis_da => disv2d_da
    procedure :: disv2d_load
    procedure :: get_dis_type => get_dis_type
    procedure :: get_dis_enum => get_dis_enum
    procedure, public :: record_array
    procedure, public :: record_srcdst_list_header
    ! -- helper functions
    procedure :: get_nodenumber_idx1
    procedure :: nodeu_to_string
    procedure :: nodeu_to_array
    procedure :: nodeu_from_string
    procedure :: nodeu_from_cellid
    procedure :: connection_normal
    procedure :: connection_vector
    procedure :: get_polyverts
    ! -- private
    procedure :: source_options
    procedure :: source_dimensions
    procedure :: source_griddata
    procedure :: log_options
    procedure :: log_dimensions
    procedure :: log_griddata
    procedure :: source_vertices
    procedure :: source_cell2d
    procedure :: define_cellverts
    procedure :: grid_finalize
    procedure :: connect
    procedure :: write_grb
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: get_cell2d_area

  end type Disv2dType

  type DisvFoundType
    logical :: length_units = .false.
    logical :: nogrb = .false.
    logical :: xorigin = .false.
    logical :: yorigin = .false.
    logical :: angrot = .false.
    logical :: nodes = .false.
    logical :: nvert = .false.
    logical :: bottom = .false.
    logical :: idomain = .false.
    logical :: iv = .false.
    logical :: xv = .false.
    logical :: yv = .false.
    logical :: icell2d = .false.
    logical :: xc = .false.
    logical :: yc = .false.
    logical :: ncvert = .false.
    logical :: icvert = .false.
  end type DisvFoundType

contains

  !> @brief Create a new discretization by vertices object
  !<
  subroutine disv2d_cr(dis, name_model, input_mempath, inunit, iout)
    ! -- dummy
    class(DisBaseType), pointer :: dis
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    ! -- local
    type(Disv2dType), pointer :: disnew
    ! -- formats
    character(len=*), parameter :: fmtheader = &
      "(1X, /1X, 'DISV -- VERTEX GRID DISCRETIZATION PACKAGE,', &
      &' VERSION 1 : 12/23/2015 - INPUT READ FROM MEMPATH: ', A, //)"
    !
    allocate (disnew)
    dis => disnew
    call disnew%allocate_scalars(name_model, input_mempath)
    dis%inunit = inunit
    dis%iout = iout
    !
    ! -- If disv enabled
    if (inunit > 0) then
      !
      ! -- Identify package
      if (iout > 0) then
        write (iout, fmtheader) dis%input_mempath
      end if
      !
      ! -- load disv
      call disnew%disv2d_load()
    end if
    !
  end subroutine disv2d_cr

  !> @brief Transfer IDM data into this discretization object
  !<
  subroutine disv2d_load(this)
    ! -- dummy
    class(Disv2dType) :: this
    !
    ! -- source input data
    call this%source_options()
    call this%source_dimensions()
    call this%source_griddata()
    call this%source_vertices()
    call this%source_cell2d()
    !
  end subroutine disv2d_load

  !> @brief Define the discretization
  !<
  subroutine disv2d_df(this)
    ! -- dummy
    class(Disv2dType) :: this
    !
    call this%grid_finalize()
    !
  end subroutine disv2d_df

  subroutine disv2d_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    ! -- dummy
    class(Disv2dType) :: this
    ! -- local

    ! -- Deallocate idm memory
    call memorystore_remove(this%name_model, 'DISV2D', idm_context)

    ! -- scalars
    call mem_deallocate(this%nvert)

    ! -- arrays
    call mem_deallocate(this%nodeuser)
    call mem_deallocate(this%nodereduced)
    call mem_deallocate(this%bottom)
    call mem_deallocate(this%idomain)

    ! -- cdl hack for arrays for vertices and cell2d blocks
    call mem_deallocate(this%vertices)
    call mem_deallocate(this%cellxy)
    call mem_deallocate(this%iavert)
    call mem_deallocate(this%javert)
    !
    ! -- DisBaseType deallocate
    call this%DisBaseType%dis_da()
  end subroutine disv2d_da

  ! !> @brief Deallocate variables
  ! !<
  ! subroutine disv2d_da(this)
  !   ! -- dummy
  !   class(Disv2dType) :: this
  !   !
  !   ! -- Deallocate idm memory
  !   call memorystore_remove(this%name_model, 'DISV', idm_context)
  !   call memorystore_remove(component=this%name_model, &
  !                          context=idm_context)
  !   !
  !   ! -- DisBaseType deallocate
  !   call this%DisBaseType%dis_da()
  !   !
  !   ! -- Deallocate scalars
  !   call mem_deallocate(this%nvert)
  !   !
  !   ! -- Deallocate Arrays
  !   call mem_deallocate(this%nodereduced)
  !   call mem_deallocate(this%nodeuser)
  !   call mem_deallocate(this%vertices)
  !   call mem_deallocate(this%cellxy)
  !   call mem_deallocate(this%iavert)
  !   call mem_deallocate(this%javert)
  !   call mem_deallocate(this%bottom)
  !   call mem_deallocate(this%idomain)
  !   !
  ! end subroutine disv2d_da

  !> @brief Copy options from IDM into package
  !<
  subroutine source_options(this)
    ! -- dummy
    class(Disv2dType) :: this
    ! -- locals
    character(len=LENVARNAME), dimension(3) :: lenunits = &
      &[character(len=LENVARNAME) :: 'FEET', 'METERS', 'CENTIMETERS']
    type(DisvFoundType) :: found
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%lenuni, 'LENGTH_UNITS', this%input_mempath, &
                       lenunits, found%length_units)
    call mem_set_value(this%nogrb, 'NOGRB', this%input_mempath, found%nogrb)
    call mem_set_value(this%xorigin, 'XORIGIN', this%input_mempath, found%xorigin)
    call mem_set_value(this%yorigin, 'YORIGIN', this%input_mempath, found%yorigin)
    call mem_set_value(this%angrot, 'ANGROT', this%input_mempath, found%angrot)
    !
    ! -- log values to list file
    if (this%iout > 0) then
      call this%log_options(found)
    end if
    !
  end subroutine source_options

  !> @brief Write user options to list file
  !<
  subroutine log_options(this, found)
    ! -- dummy
    class(Disv2dType) :: this
    type(DisvFoundType), intent(in) :: found
    !
    write (this%iout, '(1x,a)') 'Setting Discretization Options'
    !
    if (found%length_units) then
      write (this%iout, '(4x,a,i0)') 'Model length unit [0=UND, 1=FEET, &
      &2=METERS, 3=CENTIMETERS] set as ', this%lenuni
    end if
    !
    if (found%nogrb) then
      write (this%iout, '(4x,a,i0)') 'Binary grid file [0=GRB, 1=NOGRB] &
        &set as ', this%nogrb
    end if
    !
    if (found%xorigin) then
      write (this%iout, '(4x,a,G0)') 'XORIGIN = ', this%xorigin
    end if
    !
    if (found%yorigin) then
      write (this%iout, '(4x,a,G0)') 'YORIGIN = ', this%yorigin
    end if
    !
    if (found%angrot) then
      write (this%iout, '(4x,a,G0)') 'ANGROT = ', this%angrot
    end if
    !
    write (this%iout, '(1x,a,/)') 'End Setting Discretization Options'
    !
  end subroutine log_options

  !> @brief Copy dimensions from IDM into package
  !<
  subroutine source_dimensions(this)
    ! -- dummy
    class(Disv2dType) :: this
    ! -- locals
    integer(I4B) :: j
    type(DisvFoundType) :: found
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%nodes, 'NODES', this%input_mempath, found%nodes)
    call mem_set_value(this%nvert, 'NVERT', this%input_mempath, found%nvert)
    !
    ! -- log simulation values
    if (this%iout > 0) then
      call this%log_dimensions(found)
    end if
    !
    ! -- verify dimensions were set
    if (this%nodes < 1) then
      call store_error( &
        'NODES was not specified or was specified incorrectly.')
      call store_error_filename(this%input_fname)
    end if
    if (this%nvert < 1) then
      call store_error( &
        'NVERT was not specified or was specified incorrectly.')
      call store_error_filename(this%input_fname)
    end if
    !
    ! -- Calculate nodesuser
    this%nodesuser = this%nodes
    !
    ! -- Allocate non-reduced vectors for disv
    call mem_allocate(this%idomain, this%nodes, 'IDOMAIN', &
                      this%memoryPath)
    call mem_allocate(this%bottom, this%nodes, 'BOTTOM', &
                      this%memoryPath)
    !
    ! -- Allocate vertices array
    call mem_allocate(this%vertices, 2, this%nvert, 'VERTICES', this%memoryPath)
    call mem_allocate(this%cellxy, 2, this%nodesuser, 'CELLXY', this%memoryPath)
    !
    ! -- initialize all cells to be active (idomain = 1)
    do j = 1, this%nodesuser
      this%idomain(j) = 1
    end do
    !
  end subroutine source_dimensions

  !> @brief Write dimensions to list file
  !<
  subroutine log_dimensions(this, found)
    ! -- dummy
    class(Disv2dType) :: this
    type(DisvFoundType), intent(in) :: found
    !
    write (this%iout, '(1x,a)') 'Setting Discretization Dimensions'
    !
    if (found%nodes) then
      write (this%iout, '(4x,a,i0)') 'NODES = ', this%nodesuser
    end if
    !
    if (found%nvert) then
      write (this%iout, '(4x,a,i0)') 'NVERT = ', this%nvert
    end if
    !
    write (this%iout, '(1x,a,/)') 'End Setting Discretization Dimensions'
    !
  end subroutine log_dimensions

  !> @brief Copy grid data from IDM into package
  !<
  subroutine source_griddata(this)
    ! -- dummy
    class(Disv2dType) :: this
    ! -- locals
    type(DisvFoundType) :: found
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%bottom, 'BOTTOM', this%input_mempath, found%bottom)
    call mem_set_value(this%idomain, 'IDOMAIN', this%input_mempath, found%idomain)
    !
    ! -- log simulation values
    if (this%iout > 0) then
      call this%log_griddata(found)
    end if
    !
  end subroutine source_griddata

  !> @brief Write griddata found to list file
  !<
  subroutine log_griddata(this, found)
    ! -- dummy
    class(Disv2dType) :: this
    type(DisvFoundType), intent(in) :: found
    !
    write (this%iout, '(1x,a)') 'Setting Discretization Griddata'
    !
    if (found%bottom) then
      write (this%iout, '(4x,a)') 'BOTTOM set from input file'
    end if
    !
    if (found%idomain) then
      write (this%iout, '(4x,a)') 'IDOMAIN set from input file'
    end if
    !
    write (this%iout, '(1x,a,/)') 'End Setting Discretization Griddata'
    !
  end subroutine log_griddata

  !> @brief Finalize grid (check properties, allocate arrays, compute connections)
  !<
  subroutine grid_finalize(this)
    ! dummy
    class(Disv2dType) :: this
    ! locals
    integer(I4B) :: node, noder, j, ncell_count
    ! formats
    character(len=*), parameter :: fmtnr = &
      "(/1x, 'The specified IDOMAIN results in a reduced number of cells.',&
      &/1x, 'Number of user nodes: ',I0,&
      &/1X, 'Number of nodes in solution: ', I0, //)"

    ! count active cells and set nodes to that number
    ncell_count = 0
    do j = 1, this%nodesuser
      if (this%idomain(j) > 0) ncell_count = ncell_count + 1
    end do
    this%nodes = ncell_count

    ! Check to make sure nodes is a valid number
    if (ncell_count == 0) then
      call store_error('Model does not have any active nodes. &
                       &Ensure IDOMAIN array has some values greater &
                       &than zero.')
      call store_error_filename(this%input_fname)
    end if

    ! Write message if reduced grid
    if (this%nodes < this%nodesuser) then
      write (this%iout, fmtnr) this%nodesuser, this%nodes
    end if

    ! Array size is now known, so allocate
    call this%allocate_arrays()

    ! Fill the nodereduced array with the reduced nodenumber, or
    ! a negative number to indicate it is a pass-through cell, or
    ! a zero to indicate that the cell is excluded from the
    ! solution.
    if (this%nodes < this%nodesuser) then
      node = 1
      noder = 1
      do j = 1, this%nodesuser
        if (this%idomain(j) > 0) then
          this%nodereduced(node) = noder
          noder = noder + 1
        else
          this%nodereduced(node) = 0
        end if
        node = node + 1
      end do
    end if

    ! allocate and fill nodeuser if a reduced grid
    if (this%nodes < this%nodesuser) then
      node = 1
      noder = 1
      do j = 1, this%nodesuser
        if (this%idomain(j) > 0) then
          this%nodeuser(noder) = node
          noder = noder + 1
        end if
        node = node + 1
      end do
    end if

    ! Move bottom into bot
    ! and set x and y center coordinates
    node = 0
    do j = 1, this%nodesuser
      node = node + 1
      noder = node
      if (this%nodes < this%nodesuser) noder = this%nodereduced(node)
      if (noder <= 0) cycle
      this%bot(noder) = this%bottom(j)
      this%xc(noder) = this%cellxy(1, j)
      this%yc(noder) = this%cellxy(2, j)
    end do

    ! Build connections
    call this%connect()

  end subroutine grid_finalize

  !> @brief Load grid vertices from IDM into package
  !<
  subroutine source_vertices(this)
    ! -- dummy
    class(Disv2dType) :: this
    ! -- local
    integer(I4B) :: i
    real(DP), dimension(:), contiguous, pointer :: vert_x => null()
    real(DP), dimension(:), contiguous, pointer :: vert_y => null()
    !
    ! -- set pointers to memory manager input arrays
    call mem_setptr(vert_x, 'XV', this%input_mempath)
    call mem_setptr(vert_y, 'YV', this%input_mempath)
    !
    ! -- set vertices 2d array
    if (associated(vert_x) .and. associated(vert_y)) then
      do i = 1, this%nvert
        this%vertices(1, i) = vert_x(i)
        this%vertices(2, i) = vert_y(i)
      end do
    else
      call store_error('Required Vertex arrays not found.')
    end if
    !
    ! -- log
    if (this%iout > 0) then
      write (this%iout, '(1x,a)') 'Discretization Vertex data loaded'
    end if
    !
  end subroutine source_vertices

  !> @brief Build data structures to hold cell vertex info
  !<
  subroutine define_cellverts(this, icell2d, ncvert, icvert)
    ! -- modules
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(Disv2dType) :: this
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: icell2d
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: ncvert
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: icvert
    ! -- locals
    type(sparsematrix) :: vert_spm
    integer(I4B) :: i, j, ierr
    integer(I4B) :: icv_idx, startvert, maxnnz = 5
    !
    ! -- initialize sparse matrix
    call vert_spm%init(this%nodes, this%nvert, maxnnz)
    !
    ! -- add sparse matrix connections from input memory paths
    icv_idx = 1
    do i = 1, this%nodes
      if (icell2d(i) /= i) call store_error('ICELL2D input sequence violation.')
      do j = 1, ncvert(i)
        call vert_spm%addconnection(i, icvert(icv_idx), 0)
        if (j == 1) then
          startvert = icvert(icv_idx)
        elseif (j == ncvert(i) .and. (icvert(icv_idx) /= startvert)) then
          call vert_spm%addconnection(i, startvert, 0)
        end if
        icv_idx = icv_idx + 1
      end do
    end do
    !
    ! -- allocate and fill iavert and javert
    call mem_allocate(this%iavert, this%nodes + 1, 'IAVERT', this%memoryPath)
    call mem_allocate(this%javert, vert_spm%nnz, 'JAVERT', this%memoryPath)
    call vert_spm%filliaja(this%iavert, this%javert, ierr)
    call vert_spm%destroy()
    !
  end subroutine define_cellverts

  !> @brief Copy cell2d data from IDM into package
  !<
  subroutine source_cell2d(this)
    ! -- dummy
    class(Disv2dType) :: this
    ! -- locals
    integer(I4B), dimension(:), contiguous, pointer :: icell2d => null()
    integer(I4B), dimension(:), contiguous, pointer :: ncvert => null()
    integer(I4B), dimension(:), contiguous, pointer :: icvert => null()
    real(DP), dimension(:), contiguous, pointer :: cell_x => null()
    real(DP), dimension(:), contiguous, pointer :: cell_y => null()
    integer(I4B) :: i
    !
    ! -- set pointers to input path ncvert and icvert
    call mem_setptr(icell2d, 'ICELL2D', this%input_mempath)
    call mem_setptr(ncvert, 'NCVERT', this%input_mempath)
    call mem_setptr(icvert, 'ICVERT', this%input_mempath)
    !
    ! --
    if (associated(icell2d) .and. associated(ncvert) &
        .and. associated(icvert)) then
      call this%define_cellverts(icell2d, ncvert, icvert)
    else
      call store_error('Required cell vertex array(s) [ICELL2D, NCVERT, ICVERT] &
                       &not found.')
    end if
    !
    ! -- copy cell center idm sourced values to local arrays
    call mem_setptr(cell_x, 'XC', this%input_mempath)
    call mem_setptr(cell_y, 'YC', this%input_mempath)
    !
    ! -- set cell centers
    if (associated(cell_x) .and. associated(cell_y)) then
      do i = 1, this%nodesuser
        this%cellxy(1, i) = cell_x(i)
        this%cellxy(2, i) = cell_y(i)
      end do
    else
      call store_error('Required cell center arrays not found.')
    end if
    !
    ! -- log
    if (this%iout > 0) then
      write (this%iout, '(1x,a)') 'Discretization Cell2d data loaded'
    end if
    !
  end subroutine source_cell2d

  !> @brief Build grid connections
  !<
  subroutine connect(this)
    ! -- dummy
    class(Disv2dType) :: this
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: noder, nrsize
    integer(I4B) :: narea_eq_zero
    integer(I4B) :: narea_lt_zero
    real(DP) :: area
    !
    ! -- Initialize
    narea_eq_zero = 0
    narea_lt_zero = 0
    !
    ! -- Assign the cell area
    do j = 1, this%nodesuser
      area = this%get_cell2d_area(j)
      noder = this%get_nodenumber(j, 0)
      if (noder > 0) this%area(noder) = area
      if (area < DZERO) then
        narea_lt_zero = narea_lt_zero + 1
        write (errmsg, '(a,i0,a)') &
          &'Calculated CELL2D area less than zero for cell ', j, '.'
        call store_error(errmsg)
      end if
      if (area == DZERO) then
        narea_eq_zero = narea_eq_zero + 1
        write (errmsg, '(a,i0,a)') &
          'Calculated CELL2D area is zero for cell ', j, '.'
        call store_error(errmsg)
      end if
    end do
    !
    ! -- check for errors
    if (count_errors() > 0) then
      if (narea_lt_zero > 0) then
        write (errmsg, '(i0,a)') narea_lt_zero, &
          ' cell(s) have an area less than zero. Calculated cell &
          &areas must be greater than zero. Negative areas often &
          &mean vertices are not listed in clockwise order.'
        call store_error(errmsg)
      end if
      if (narea_eq_zero > 0) then
        write (errmsg, '(i0,a)') narea_eq_zero, &
          ' cell(s) have an area equal to zero. Calculated cell &
          &areas must be greater than zero. Calculated cell &
          &areas equal to zero indicate that the cell is not defined &
          &by a valid polygon.'
        call store_error(errmsg)
      end if
      call store_error_filename(this%input_fname)
    end if
    !
    ! -- create and fill the connections object
    nrsize = 0
    if (this%nodes < this%nodesuser) nrsize = this%nodes
    allocate (this%con)
    call this%con%disvconnections(this%name_model, this%nodes, &
                                  this%nodesuser, 1, nrsize, &
                                  this%nvert, this%vertices, this%iavert, &
                                  this%javert, this%cellxy, &
                                  this%bot, this%bot, &
                                  this%nodereduced, this%nodeuser)
    this%nja = this%con%nja
    this%njas = this%con%njas
    !
  end subroutine connect

  !> @brief Write a binary grid file
  !<
  subroutine write_grb(this, icelltype)
    ! -- modules
    use OpenSpecModule, only: access, form
    ! -- dummy
    class(Disv2dType) :: this
    integer(I4B), dimension(:), intent(in) :: icelltype
    ! -- local
    integer(I4B) :: iunit, i, ntxt
    integer(I4B), parameter :: lentxt = 100
    character(len=50) :: txthdr
    character(len=lentxt) :: txt
    character(len=LINELENGTH) :: fname
    ! -- formats
    character(len=*), parameter :: fmtgrdsave = &
      "(4X,'BINARY GRID INFORMATION WILL BE WRITTEN TO:', &
       &/,6X,'UNIT NUMBER: ', I0,/,6X, 'FILE NAME: ', A)"
    !
    ! -- Initialize
    ntxt = 18
    !
    ! -- Open the file
    fname = trim(this%output_fname)
    iunit = getunit()
    write (this%iout, fmtgrdsave) iunit, trim(adjustl(fname))
    call openfile(iunit, this%iout, trim(adjustl(fname)), 'DATA(BINARY)', &
                  form, access, 'REPLACE')
    !
    ! -- write header information
    write (txthdr, '(a)') 'GRID DISV2D'
    txthdr(50:50) = new_line('a')
    write (iunit) txthdr
    write (txthdr, '(a)') 'VERSION 1'
    txthdr(50:50) = new_line('a')
    write (iunit) txthdr
    write (txthdr, '(a, i0)') 'NTXT ', ntxt
    txthdr(50:50) = new_line('a')
    write (iunit) txthdr
    write (txthdr, '(a, i0)') 'LENTXT ', lentxt
    txthdr(50:50) = new_line('a')
    write (iunit) txthdr
    !
    ! -- write variable definitions
    write (txt, '(3a, i0)') 'NCELLS ', 'INTEGER ', 'NDIM 0 # ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'NODES ', 'INTEGER ', 'NDIM 0 # ', this%nodes
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'NVERT ', 'INTEGER ', 'NDIM 0 # ', this%nvert
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'NJAVERT ', 'INTEGER ', 'NDIM 0 # ', size(this%javert)
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'NJA ', 'INTEGER ', 'NDIM 0 # ', this%con%nja
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, 1pg25.15e3)') &
      'XORIGIN ', 'DOUBLE ', 'NDIM 0 # ', this%xorigin
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, 1pg25.15e3)') &
      'YORIGIN ', 'DOUBLE ', 'NDIM 0 # ', this%yorigin
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, 1pg25.15e3)') 'ANGROT ', 'DOUBLE ', 'NDIM 0 # ', this%angrot
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'BOTM ', 'DOUBLE ', 'NDIM 1 ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'VERTICES ', 'DOUBLE ', 'NDIM 2 2 ', this%nvert
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'CELLX ', 'DOUBLE ', 'NDIM 1 ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'CELLY ', 'DOUBLE ', 'NDIM 1 ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'IAVERT ', 'INTEGER ', 'NDIM 1 ', this%nodesuser + 1
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'JAVERT ', 'INTEGER ', 'NDIM 1 ', size(this%javert)
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'IA ', 'INTEGER ', 'NDIM 1 ', this%nodesuser + 1
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'JA ', 'INTEGER ', 'NDIM 1 ', size(this%con%jausr)
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'IDOMAIN ', 'INTEGER ', 'NDIM 1 ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'ICELLTYPE ', 'INTEGER ', 'NDIM 1 ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    !
    ! -- write data
    write (iunit) this%nodesuser ! ncells
    write (iunit) this%nodes ! nodes
    write (iunit) this%nvert ! nvert
    write (iunit) size(this%javert) ! njavert
    write (iunit) this%nja ! nja
    write (iunit) this%xorigin ! xorigin
    write (iunit) this%yorigin ! yorigin
    write (iunit) this%angrot ! angrot
    write (iunit) this%bottom ! botm
    write (iunit) this%vertices ! vertices
    write (iunit) (this%cellxy(1, i), i=1, this%nodesuser) ! cellx
    write (iunit) (this%cellxy(2, i), i=1, this%nodesuser) ! celly
    write (iunit) this%iavert ! iavert
    write (iunit) this%javert ! javert
    write (iunit) this%con%iausr ! iausr
    write (iunit) this%con%jausr ! jausr
    write (iunit) this%idomain ! idomain
    write (iunit) icelltype ! icelltype
    !
    ! -- Close the file
    close (iunit)
    !
  end subroutine write_grb

  !> @brief Convert a user nodenumber to a string (nodenumber) or (k,j)
  !<
  subroutine nodeu_to_string(this, nodeu, str)
    ! -- dummy
    class(Disv2dType) :: this
    integer(I4B), intent(in) :: nodeu
    character(len=*), intent(inout) :: str
    ! -- local
    integer(I4B) :: i, j, k
    character(len=10) :: jstr
    !
    call get_ijk(nodeu, 1, this%nodes, 1, i, j, k)
    write (jstr, '(i10)') j
    str = '('//trim(adjustl(jstr))//')'
    !
  end subroutine nodeu_to_string

  !> @brief Convert a user nodenumber to an array (nodenumber) or (k,j)
  !<
  subroutine nodeu_to_array(this, nodeu, arr)
    ! -- dummy
    class(Disv2dType) :: this
    integer(I4B), intent(in) :: nodeu
    integer(I4B), dimension(:), intent(inout) :: arr
    ! -- local
    integer(I4B) :: isize
    integer(I4B) :: i, j, k
    !
    ! -- check the size of arr
    isize = size(arr)
    if (isize /= this%ndim) then
      write (errmsg, '(a,i0,a,i0,a)') &
        'Program error: nodeu_to_array size of array (', isize, &
        ') is not equal to the discretization dimension (', this%ndim, ').'
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- get k, i, j
    call get_ijk(nodeu, 1, this%nodes, 1, i, j, k)
    !
    ! -- fill array
    arr(1) = j
    !
  end subroutine nodeu_to_array

  !> @brief Get reduced node number from user node number
  !<
  function get_nodenumber_idx1(this, nodeu, icheck) result(nodenumber)
    ! return
    integer(I4B) :: nodenumber
    ! dummy
    class(Disv2dType), intent(in) :: this
    integer(I4B), intent(in) :: nodeu
    integer(I4B), intent(in) :: icheck
    ! local

    ! check the node number if requested
    if (icheck /= 0) then

      ! If within valid range, convert to reduced nodenumber
      if (nodeu < 1 .or. nodeu > this%nodesuser) then
        nodenumber = 0
        write (errmsg, '(a,i0,a,i0,a)') &
          'Node number (', nodeu, ') is less than 1 or greater than nodes (', &
          this%nodesuser, ').'
        call store_error(errmsg)
      else
        nodenumber = nodeu
        if (this%nodes < this%nodesuser) nodenumber = this%nodereduced(nodeu)
      end if
    else
      nodenumber = nodeu
      if (this%nodes < this%nodesuser) nodenumber = this%nodereduced(nodeu)
    end if

  end function get_nodenumber_idx1

  !> @brief Get normal vector components between the cell and a given neighbor
  !!
  !! The normal points outward from the shared face between noden and nodem.
  !<
  subroutine connection_normal(this, noden, nodem, ihc, xcomp, ycomp, zcomp, &
                               ipos)
    ! -- dummy
    class(Disv2dType) :: this
    integer(I4B), intent(in) :: noden !< cell (reduced nn)
    integer(I4B), intent(in) :: nodem !< neighbor (reduced nn)
    integer(I4B), intent(in) :: ihc !< horizontal connection flag
    real(DP), intent(inout) :: xcomp
    real(DP), intent(inout) :: ycomp
    real(DP), intent(inout) :: zcomp
    integer(I4B), intent(in) :: ipos
    ! -- local
    real(DP) :: angle, dmult
    !
    ! -- Set vector components based on ihc
    if (ihc == 0) then
      xcomp = DZERO
      ycomp = DZERO
      if (nodem < noden) then
        !
        ! -- nodem must be above noden, so upward connection
        zcomp = DONE
      else
        !
        ! -- nodem must be below noden, so downward connection
        zcomp = -DONE
      end if
    else
      ! -- find from anglex, since anglex is symmetric, need to flip vector
      !    for lower triangle (nodem < noden)
      !ipos = this%con%getjaindex(noden, nodem)
      angle = this%con%anglex(this%con%jas(ipos))
      dmult = DONE
      if (nodem < noden) dmult = -DONE
      xcomp = cos(angle) * dmult
      ycomp = sin(angle) * dmult
      zcomp = DZERO
    end if
    !
  end subroutine connection_normal

  !> @brief Get unit vector components between the cell and a given neighbor
  !!
  !! Saturation must be provided to compute cell center vertical coordinates.
  !! Also return the straight-line connection length.
  !<
  subroutine connection_vector(this, noden, nodem, nozee, satn, satm, ihc, &
                               xcomp, ycomp, zcomp, conlen)
    ! -- dummy
    class(Disv2dType) :: this
    integer(I4B), intent(in) :: noden !< cell (reduced nn)
    integer(I4B), intent(in) :: nodem !< neighbor (reduced nn)
    logical, intent(in) :: nozee !< do not use z in calculations
    real(DP), intent(in) :: satn !< not used for disv1d
    real(DP), intent(in) :: satm !< not used for disv1d
    integer(I4B), intent(in) :: ihc !< horizontal connection flag
    real(DP), intent(inout) :: xcomp !< x component of connection vector
    real(DP), intent(inout) :: ycomp !< y component of connection vector
    real(DP), intent(inout) :: zcomp !< z component of connection vector
    real(DP), intent(inout) :: conlen !< calculated straight-line distance between cell centers
    ! -- local
    integer(I4B) :: nodeun, nodeum
    real(DP) :: xn, xm, yn, ym, zn, zm

    ! horizontal connection, with possible z component due to cell offsets
    ! and/or water table conditions
    if (nozee) then
      zn = DZERO
      zm = DZERO
    else
      zn = this%bot(noden)
      zm = this%bot(nodem)
    end if
    nodeun = this%get_nodeuser(noden)
    nodeum = this%get_nodeuser(nodem)
    xn = this%cellxy(1, nodeun)
    yn = this%cellxy(2, nodeun)
    xm = this%cellxy(1, nodeum)
    ym = this%cellxy(2, nodeum)
    call line_unit_vector(xn, yn, zn, xm, ym, zm, xcomp, ycomp, zcomp, &
                          conlen)

  end subroutine connection_vector

  !> @brief Get the discretization type
  !<
  subroutine get_dis_type(this, dis_type)
    ! -- dummy
    class(Disv2dType), intent(in) :: this
    character(len=*), intent(out) :: dis_type
    !
    dis_type = "DISV2D"
    !
  end subroutine get_dis_type

  !> @brief Get the discretization type enumeration
  function get_dis_enum(this) result(dis_enum)
    use ConstantsModule, only: DISV2D
    class(Disv2dType), intent(in) :: this
    integer(I4B) :: dis_enum
    dis_enum = DISV2D
  end function get_dis_enum

  !> @brief Allocate and initialize scalars
  !<
  subroutine allocate_scalars(this, name_model, input_mempath)
    ! -- dummy
    class(Disv2dType) :: this
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    !
    ! -- Allocate parent scalars
    call this%DisBaseType%allocate_scalars(name_model, input_mempath)
    !
    ! -- Allocate
    call mem_allocate(this%nvert, 'NVERT', this%memoryPath)
    !
    ! -- Initialize
    this%nvert = 0
    this%ndim = 1
    !
  end subroutine allocate_scalars

  !> @brief Allocate and initialize arrays
  !<
  subroutine allocate_arrays(this)
    ! dummy
    class(Disv2dType) :: this

    ! Allocate arrays in DisBaseType (mshape, top, bot, area)
    call this%DisBaseType%allocate_arrays()
    !
    ! Allocate arrays for DisvType
    if (this%nodes < this%nodesuser) then
      call mem_allocate(this%nodeuser, this%nodes, 'NODEUSER', this%memoryPath)
      call mem_allocate(this%nodereduced, this%nodesuser, 'NODEREDUCED', &
                        this%memoryPath)
    else
      call mem_allocate(this%nodeuser, 1, 'NODEUSER', this%memoryPath)
      call mem_allocate(this%nodereduced, 1, 'NODEREDUCED', this%memoryPath)
    end if

    ! Initialize
    this%mshape(1) = this%nodesuser

  end subroutine allocate_arrays

  !> @brief Get the signed area of the cell
  !!
  !! A negative result means points are in counter-clockwise orientation.
  !! Area is computed from the formula:
  !! a = 1/2 *[(x1*y2 + x2*y3 + x3*y4 + ... + xn*y1) -
  !!           (x2*y1 + x3*y2 + x4*y3 + ... + x1*yn)]
  !<
  function get_cell2d_area(this, icell2d) result(area)
    ! -- dummy
    class(Disv2dType) :: this
    integer(I4B), intent(in) :: icell2d
    ! -- return
    real(DP) :: area
    ! -- local
    integer(I4B) :: ivert
    integer(I4B) :: nvert
    integer(I4B) :: icount
    integer(I4B) :: iv1
    real(DP) :: x
    real(DP) :: y
    real(DP) :: x1
    real(DP) :: y1
    !
    area = DZERO
    nvert = this%iavert(icell2d + 1) - this%iavert(icell2d)
    icount = 1
    iv1 = this%javert(this%iavert(icell2d))
    x1 = this%vertices(1, iv1)
    y1 = this%vertices(2, iv1)
    do ivert = this%iavert(icell2d), this%iavert(icell2d + 1) - 1
      x = this%vertices(1, this%javert(ivert))
      if (icount < nvert) then
        y = this%vertices(2, this%javert(ivert + 1))
      else
        y = this%vertices(2, this%javert(this%iavert(icell2d)))
      end if
      area = area + (x - x1) * (y - y1)
      icount = icount + 1
    end do
    !
    icount = 1
    do ivert = this%iavert(icell2d), this%iavert(icell2d + 1) - 1
      y = this%vertices(2, this%javert(ivert))
      if (icount < nvert) then
        x = this%vertices(1, this%javert(ivert + 1))
      else
        x = this%vertices(1, this%javert(this%iavert(icell2d)))
      end if
      area = area - (x - x1) * (y - y1)
      icount = icount + 1
    end do
    !
    area = -DONE * area * DHALF
    !
  end function get_cell2d_area

  !> @brief Convert a string to a user nodenumber
  !!
  !! Parse layer and within-layer cell number and return user nodenumber.
  !! If flag_string is present and true, the first token may be
  !! non-numeric (e.g. boundary name). In this case, return -2.
  !<
  function nodeu_from_string(this, lloc, istart, istop, in, iout, line, &
                             flag_string, allow_zero) result(nodeu)
    ! -- dummy
    class(Disv2dType) :: this
    integer(I4B), intent(inout) :: lloc
    integer(I4B), intent(inout) :: istart
    integer(I4B), intent(inout) :: istop
    integer(I4B), intent(in) :: in
    integer(I4B), intent(in) :: iout
    character(len=*), intent(inout) :: line
    logical, optional, intent(in) :: flag_string
    logical, optional, intent(in) :: allow_zero
    integer(I4B) :: nodeu
    ! -- local
    integer(I4B) :: j, nodes
    integer(I4B) :: lloclocal, ndum, istat, n
    real(DP) :: r
    !
    if (present(flag_string)) then
      if (flag_string) then
        ! Check to see if first token in line can be read as an integer.
        lloclocal = lloc
        call urword(line, lloclocal, istart, istop, 1, ndum, r, iout, in)
        read (line(istart:istop), *, iostat=istat) n
        if (istat /= 0) then
          ! First token in line is not an integer; return flag to this effect.
          nodeu = -2
          return
        end if
      end if
    end if
    !
    nodes = this%mshape(1)
    !
    call urword(line, lloc, istart, istop, 2, j, r, iout, in)
    !
    if (j == 0) then
      if (present(allow_zero)) then
        if (allow_zero) then
          nodeu = 0
          return
        end if
      end if
    end if
    !
    errmsg = ''
    !
    if (j < 1 .or. j > nodes) then
      write (errmsg, '(a,1x,a,i0,a)') &
        trim(adjustl(errmsg)), 'Cell number in list (', j, &
        ') is outside of the grid.'
    end if
    !
    nodeu = get_node(1, 1, j, 1, 1, nodes)
    !
    if (nodeu < 1 .or. nodeu > this%nodesuser) then
      write (errmsg, '(a,1x,a,i0,a)') &
        trim(adjustl(errmsg)), &
        "Node number in list (", nodeu, ") is outside of the grid. "// &
        "Cell number cannot be determined in line '"// &
        trim(adjustl(line))//"'."
    end if
    !
    if (len_trim(adjustl(errmsg)) > 0) then
      call store_error(errmsg)
      call store_error_unit(in)
    end if
    !
  end function nodeu_from_string

  !> @brief Convert a cellid string to a user nodenumber
  !!
  !! If flag_string is present and true, the first token may be
  !! non-numeric (e.g. boundary name). In this case, return -2.
  !!
  !! If allow_zero is present and true, and all indices are zero, the
  !! result can be zero. If allow_zero is false, a zero in any index is an error.
  !<
  function nodeu_from_cellid(this, cellid, inunit, iout, flag_string, &
                             allow_zero) result(nodeu)
    ! -- return
    integer(I4B) :: nodeu
    ! -- dummy
    class(Disv2dType) :: this
    character(len=*), intent(inout) :: cellid
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    logical, optional, intent(in) :: flag_string
    logical, optional, intent(in) :: allow_zero
    ! -- local
    integer(I4B) :: j, nodes
    integer(I4B) :: lloclocal, ndum, istat, n
    integer(I4B) :: istart, istop
    real(DP) :: r
    !
    if (present(flag_string)) then
      if (flag_string) then
        ! Check to see if first token in cellid can be read as an integer.
        lloclocal = 1
        call urword(cellid, lloclocal, istart, istop, 1, ndum, r, iout, inunit)
        read (cellid(istart:istop), *, iostat=istat) n
        if (istat /= 0) then
          ! First token in cellid is not an integer; return flag to this effect.
          nodeu = -2
          return
        end if
      end if
    end if
    !
    nodes = this%mshape(1)
    !
    lloclocal = 1
    call urword(cellid, lloclocal, istart, istop, 2, j, r, iout, inunit)
    !
    if (j == 0) then
      if (present(allow_zero)) then
        if (allow_zero) then
          nodeu = 0
          return
        end if
      end if
    end if
    !
    errmsg = ''
    !
    if (j < 1 .or. j > nodes) then
      write (errmsg, '(a,1x,a,i0,a)') &
        trim(adjustl(errmsg)), 'Cell2d number in list (', j, &
        ') is outside of the grid.'
    end if
    !
    nodeu = get_node(1, 1, j, 1, 1, nodes)
    !
    if (nodeu < 1 .or. nodeu > this%nodesuser) then
      write (errmsg, '(a,1x,a,i0,a)') &
        trim(adjustl(errmsg)), &
        "Cell number cannot be determined for cellid ("// &
        trim(adjustl(cellid))//") and results in a user "// &
        "node number (", nodeu, ") that is outside of the grid."
    end if
    !
    if (len_trim(adjustl(errmsg)) > 0) then
      call store_error(errmsg)
      call store_error_unit(inunit)
    end if
    !
  end function nodeu_from_cellid

  !> @brief Get a 2D array of polygon vertices, listed in clockwise order
  !! beginning with the lower left corner
  !<
  subroutine get_polyverts(this, ic, polyverts, closed)
    ! -- dummy
    class(Disv2dType), intent(inout) :: this
    integer(I4B), intent(in) :: ic !< cell number (reduced)
    real(DP), allocatable, intent(out) :: polyverts(:, :) !< polygon vertices (column-major indexing)
    logical(LGP), intent(in), optional :: closed !< whether to close the polygon, duplicating a vertex (default false)
    ! -- local
    integer(I4B) :: icu, icu2d, iavert, nverts, m, j
    logical(LGP) :: lclosed
    !
    ! count vertices
    icu = this%get_nodeuser(ic)
    icu2d = icu - ((icu - 1) / this%nodes) * this%nodes
    nverts = this%iavert(icu2d + 1) - this%iavert(icu2d) - 1
    if (nverts .le. 0) nverts = nverts + size(this%javert)
    !
    ! check closed option
    if (.not. (present(closed))) then
      lclosed = .false.
    else
      lclosed = closed
    end if
    !
    ! allocate vertices array
    if (lclosed) then
      allocate (polyverts(2, nverts + 1))
    else
      allocate (polyverts(2, nverts))
    end if
    !
    ! set vertices
    iavert = this%iavert(icu2d)
    do m = 1, nverts
      j = this%javert(iavert - 1 + m)
      polyverts(:, m) = (/this%vertices(1, j), this%vertices(2, j)/)
    end do
    !
    ! close if enabled
    if (lclosed) &
      polyverts(:, nverts + 1) = polyverts(:, 1)
    !
  end subroutine

  !> @brief Record a double precision array
  !!
  !! The array is written to a formatted or unformatted external file depending
  !! on the arguments.
  !<
  subroutine record_array(this, darray, iout, iprint, idataun, aname, &
                          cdatafmp, nvaluesp, nwidthp, editdesc, dinact)
    ! -- dummy
    class(Disv2dType), intent(inout) :: this
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: darray !< double precision array to record
    integer(I4B), intent(in) :: iout !< ascii output unit number
    integer(I4B), intent(in) :: iprint !< whether to print the array
    integer(I4B), intent(in) :: idataun !< binary output unit number, if negative don't write by layers, write entire array
    character(len=*), intent(in) :: aname !< text descriptor
    character(len=*), intent(in) :: cdatafmp !< write format
    integer(I4B), intent(in) :: nvaluesp !< values per line
    integer(I4B), intent(in) :: nwidthp !< number width
    character(len=*), intent(in) :: editdesc !< format type (I, G, F, S, E)
    real(DP), intent(in) :: dinact !< double precision value for cells excluded from model domain
    ! -- local
    integer(I4B) :: k, ifirst
    integer(I4B) :: nlay
    integer(I4B) :: nrow
    integer(I4B) :: ncol
    integer(I4B) :: nval
    integer(I4B) :: nodeu, noder
    integer(I4B) :: istart, istop
    real(DP), dimension(:), pointer, contiguous :: dtemp
    ! -- formats
    character(len=*), parameter :: fmthsv = &
      "(1X,/1X,a,' WILL BE SAVED ON UNIT ',I4, &
       &' AT END OF TIME STEP',I5,', STRESS PERIOD ',I4)"
    !
    ! -- set variables
    nlay = 1
    nrow = 1
    ncol = this%mshape(1)
    !
    ! -- If this is a reduced model, then copy the values from darray into
    !    dtemp.
    if (this%nodes < this%nodesuser) then
      nval = this%nodes
      dtemp => this%dbuff
      do nodeu = 1, this%nodesuser
        noder = this%get_nodenumber(nodeu, 0)
        if (noder <= 0) then
          dtemp(nodeu) = dinact
          cycle
        end if
        dtemp(nodeu) = darray(noder)
      end do
    else
      nval = this%nodes
      dtemp => darray
    end if
    !
    ! -- Print to iout if iprint /= 0
    if (iprint /= 0) then
      istart = 1
      do k = 1, nlay
        istop = istart + nrow * ncol - 1
        call ulaprufw(ncol, nrow, kstp, kper, k, iout, dtemp(istart:istop), &
                      aname, cdatafmp, nvaluesp, nwidthp, editdesc)
        istart = istop + 1
      end do
    end if
    !
    ! -- Save array to an external file.
    if (idataun > 0) then
      ! -- write to binary file by layer
      ifirst = 1
      istart = 1
      do k = 1, nlay
        istop = istart + nrow * ncol - 1
        if (ifirst == 1) write (iout, fmthsv) &
          trim(adjustl(aname)), idataun, &
          kstp, kper
        ifirst = 0
        call ulasav(dtemp(istart:istop), aname, kstp, kper, &
                    pertim, totim, ncol, nrow, k, idataun)
        istart = istop + 1
      end do
    elseif (idataun < 0) then
      !
      ! -- write entire array as one record
      call ubdsv1(kstp, kper, aname, -idataun, dtemp, ncol, nrow, nlay, &
                  iout, delt, pertim, totim)
    end if
    !
  end subroutine record_array

  !> @brief Record list header for imeth=6
  !<
  subroutine record_srcdst_list_header(this, text, textmodel, textpackage, &
                                       dstmodel, dstpackage, naux, auxtxt, &
                                       ibdchn, nlist, iout)
    ! -- dummy
    class(Disv2dType) :: this
    character(len=16), intent(in) :: text
    character(len=16), intent(in) :: textmodel
    character(len=16), intent(in) :: textpackage
    character(len=16), intent(in) :: dstmodel
    character(len=16), intent(in) :: dstpackage
    integer(I4B), intent(in) :: naux
    character(len=16), dimension(:), intent(in) :: auxtxt
    integer(I4B), intent(in) :: ibdchn
    integer(I4B), intent(in) :: nlist
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B) :: nlay, nrow, ncol
    !
    nlay = 1
    nrow = 1
    ncol = this%mshape(1)
    !
    ! -- Use ubdsv06 to write list header
    call ubdsv06(kstp, kper, text, textmodel, textpackage, dstmodel, dstpackage, &
                 ibdchn, naux, auxtxt, ncol, nrow, nlay, &
                 nlist, iout, delt, pertim, totim)
    !
  end subroutine record_srcdst_list_header

end module Disv2dModule
