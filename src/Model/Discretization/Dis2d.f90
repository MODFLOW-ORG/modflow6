module Dis2dModule

  use ArrayReadersModule, only: ReadArray
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, DHALF, DONE, DZERO, &
                             LENMEMPATH, LENVARNAME
  use BaseDisModule, only: DisBaseType
  use GeomUtilModule, only: get_node, get_ijk
  use InputOutputModule, only: URWORD, ulasav, ulaprufw, &
                               ubdsv1, ubdsv06, urword, getunit, openfile
  use SimModule, only: count_errors, store_error, store_error_unit, &
                       store_error_filename
  use SimVariablesModule, only: errmsg, idm_context
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use MemoryManagerExtModule, only: mem_set_value, memorystore_remove
  use TdisModule, only: kstp, kper, pertim, totim, delt

  implicit none
  private
  public dis2d_cr, Dis2dType

  !> @brief Structured grid discretization
  type, extends(DisBaseType) :: Dis2dType
    integer(I4B), pointer :: nrow => null() !< number of rows
    integer(I4B), pointer :: ncol => null() !< number of columns
    real(DP), dimension(:), pointer, contiguous :: delr => null() !< spacing along a row
    real(DP), dimension(:), pointer, contiguous :: delc => null() !< spacing along a column
    real(DP), dimension(:, :), pointer, contiguous :: bottom => null() !< bottom elevations for each cell (ncol, nrow)
    real(DP), dimension(:), pointer, contiguous :: cellx => null() !< cell center x coordinate for column j
    real(DP), dimension(:), pointer, contiguous :: celly => null() !< cell center y coordinate for row i

  contains

    procedure :: dis_df => dis3d_df
    procedure :: dis_da => dis3d_da
    procedure :: get_dis_type => get_dis_type
    procedure :: get_dis_enum => get_dis_enum
    procedure, public :: record_array
    procedure, public :: read_layer_array
    procedure, public :: record_srcdst_list_header
    procedure, public :: nlarray_to_nodelist
    ! -- helper functions
    procedure :: get_nodenumber_idx1
    procedure :: get_nodenumber_idx2
    procedure :: nodeu_to_string
    procedure :: nodeu_to_array
    procedure :: nodeu_from_string
    procedure :: nodeu_from_cellid
    procedure :: supports_layers
    procedure :: get_ncpl
    procedure :: get_polyverts
    procedure :: connection_vector
    procedure :: connection_normal
    ! -- private
    procedure :: source_options
    procedure :: source_dimensions
    procedure :: source_griddata
    procedure :: log_options
    procedure :: log_dimensions
    procedure :: log_griddata
    procedure :: grid_finalize
    procedure :: write_grb
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    !
    ! -- Read a node-sized model array (reduced or not)
    procedure :: read_int_array
    procedure :: read_dbl_array
  end type Dis2dType

  !> @brief Simplifies tracking parameters sourced from the input context.
  type DisFoundType
    logical :: length_units = .false.
    logical :: nogrb = .false.
    logical :: xorigin = .false.
    logical :: yorigin = .false.
    logical :: angrot = .false.
    logical :: nrow = .false.
    logical :: ncol = .false.
    logical :: delr = .false.
    logical :: delc = .false.
    logical :: bottom = .false.
    logical :: idomain = .false.
  end type DisFoundtype

contains

  !> @brief Create a new structured discretization object
  !<
  subroutine dis2d_cr(dis, name_model, input_mempath, inunit, iout)
    ! -- dummy
    class(DisBaseType), pointer :: dis
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    ! -- locals
    type(Dis2dType), pointer :: disnew
    ! -- formats
    character(len=*), parameter :: fmtheader = &
      "(1X, /1X, 'DIS -- STRUCTURED GRID DISCRETIZATION PACKAGE,', &
      &' VERSION 2 : 3/27/2014 - INPUT READ FROM MEMPATH: ', A, /)"
    !
    allocate (disnew)
    dis => disnew
    call disnew%allocate_scalars(name_model, input_mempath)
    dis%inunit = inunit
    dis%iout = iout
    !
    ! -- If dis enabled
    if (inunit > 0) then
      !
      ! -- Identify package
      if (iout > 0) then
        write (iout, fmtheader) dis%input_mempath
      end if
    end if
    !
  end subroutine dis2d_cr

  !> @brief Define the discretization
  !<
  subroutine dis3d_df(this)
    ! -- dummy
    class(Dis2dType) :: this
    !
    ! -- Transfer the data from the memory manager into this package object
    if (this%inunit /= 0) then
      !
      ! -- source input options
      call this%source_options()
      !
      ! -- source input dimensions
      call this%source_dimensions()
      !
      ! -- source input griddata
      call this%source_griddata()
    end if
    !
    ! -- Final grid initialization
    call this%grid_finalize()
    !
  end subroutine dis3d_df

  !> @brief Deallocate variables
  !<
  subroutine dis3d_da(this)
    ! -- dummy
    class(Dis2dType) :: this
    !
    ! -- Deallocate idm memory
    call memorystore_remove(this%name_model, 'DIS2D', idm_context)
    !
    ! -- DisBaseType deallocate
    call this%DisBaseType%dis_da()
    !
    ! -- Deallocate scalars
    call mem_deallocate(this%nrow)
    call mem_deallocate(this%ncol)
    call mem_deallocate(this%delr)
    call mem_deallocate(this%delc)
    call mem_deallocate(this%cellx)
    call mem_deallocate(this%celly)
    !
    ! -- Deallocate Arrays
    call mem_deallocate(this%nodereduced)
    call mem_deallocate(this%nodeuser)
    call mem_deallocate(this%bottom)
    call mem_deallocate(this%idomain)
    !
  end subroutine dis3d_da

  !> @brief Copy options from IDM into package
  !<
  subroutine source_options(this)
    ! -- dummy
    class(Dis2dType) :: this
    ! -- locals
    character(len=LENVARNAME), dimension(3) :: lenunits = &
      &[character(len=LENVARNAME) :: 'FEET', 'METERS', 'CENTIMETERS']
    type(DisFoundType) :: found
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
    class(Dis2dType) :: this
    type(DisFoundType), intent(in) :: found
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
    class(Dis2dType) :: this
    ! -- locals
    type(DisFoundType) :: found
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%nrow, 'NROW', this%input_mempath, found%nrow)
    call mem_set_value(this%ncol, 'NCOL', this%input_mempath, found%ncol)
    !
    ! -- log simulation values
    if (this%iout > 0) then
      call this%log_dimensions(found)
    end if
    !
    ! -- verify dimensions were set
    if (this%nrow < 1) then
      call store_error( &
        'NROW was not specified or was specified incorrectly.')
      call store_error_filename(this%input_fname)
    end if
    if (this%ncol < 1) then
      call store_error( &
        'NCOL was not specified or was specified incorrectly.')
      call store_error_filename(this%input_fname)
    end if
    !
    ! -- calculate nodesuser
    this%nodesuser = this%nrow * this%ncol
    !
    ! -- Allocate delr, delc, and non-reduced vectors for dis
    call mem_allocate(this%delr, this%ncol, 'DELR', this%memoryPath)
    call mem_allocate(this%delc, this%nrow, 'DELC', this%memoryPath)
    call mem_allocate(this%idomain, this%nodesuser, 'IDOMAIN', this%memoryPath)
    call mem_allocate(this%bottom, this%ncol, this%nrow, 'BOTTOM', &
                      this%memoryPath)
    call mem_allocate(this%cellx, this%ncol, 'CELLX', this%memoryPath)
    call mem_allocate(this%celly, this%nrow, 'CELLY', this%memoryPath)
    !
  end subroutine source_dimensions

  !> @brief Write dimensions to list file
  !<
  subroutine log_dimensions(this, found)
    ! -- dummy
    class(Dis2dType) :: this
    type(DisFoundType), intent(in) :: found
    !
    write (this%iout, '(1x,a)') 'Setting Discretization Dimensions'
    !
    if (found%nrow) then
      write (this%iout, '(4x,a,i0)') 'NROW = ', this%nrow
    end if
    !
    if (found%ncol) then
      write (this%iout, '(4x,a,i0)') 'NCOL = ', this%ncol
    end if
    !
    write (this%iout, '(1x,a,/)') 'End Setting Discretization Dimensions'
    !
  end subroutine log_dimensions

  !> @brief Copy grid data from IDM into package
  !<
  subroutine source_griddata(this)
    ! -- dummy
    class(Dis2dType) :: this
    ! -- local
    type(DisFoundType) :: found
    integer(I4B) :: i, j
    integer(I4B), contiguous, pointer :: idomain(:, :)
    !
    ! -- initialize all cells to be active (idomain = 1)
    allocate (idomain(this%ncol, this%nrow))
    do i = 1, this%nrow
      do j = 1, this%ncol
        idomain(j, i) = 1
      end do
    end do
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%delr, 'DELR', this%input_mempath, found%delr)
    call mem_set_value(this%delc, 'DELC', this%input_mempath, found%delc)
    call mem_set_value(this%bottom, 'BOTTOM', this%input_mempath, found%bottom)
    call mem_set_value(idomain, 'IDOMAIN', this%input_mempath, found%idomain)
    this%idomain = reshape(idomain, [this%nodesuser])
    deallocate (idomain)
    !
    ! -- log simulation values
    if (this%iout > 0) then
      call this%log_griddata(found)
    end if
    !
  end subroutine source_griddata

  !> @brief Write dimensions to list file
  !<
  subroutine log_griddata(this, found)
    ! -- dummy
    class(Dis2dType) :: this
    type(DisFoundType), intent(in) :: found
    !
    write (this%iout, '(1x,a)') 'Setting Discretization Griddata'
    !
    if (found%delr) then
      write (this%iout, '(4x,a)') 'DELR set from input file'
    end if
    !
    if (found%delc) then
      write (this%iout, '(4x,a)') 'DELC set from input file'
    end if
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
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(Dis2dType) :: this
    ! -- locals
    integer(I4B) :: i, j
    integer(I4B) :: node
    integer(I4B) :: noder
    integer(I4B) :: nrsize
    ! -- formats
    character(len=*), parameter :: fmtdz = &
      "('CELL (',i0,',',i0,',',i0,') THICKNESS <= 0. ', &
      &'TOP, BOT: ',2(1pg24.15))"
    character(len=*), parameter :: fmtnr = &
      "(/1x, 'The specified IDOMAIN results in a reduced number of cells.',&
      &/1x, 'Number of user nodes: ',I0,&
      &/1X, 'Number of nodes in solution: ', I0, //)"
    !
    ! -- count active cells
    this%nodes = 0
    do i = 1, this%nodesuser
      if (this%idomain(i) > 0) this%nodes = this%nodes + 1
    end do
    !
    ! -- Check to make sure nodes is a valid number
    if (this%nodes == 0) then
      call store_error('Model does not have any active nodes. &
                       &Ensure IDOMAIN array has some values greater &
                       &than zero.')
      call store_error_filename(this%input_fname)
    end if
    !
    ! -- Write message if reduced grid
    if (this%nodes < this%nodesuser) then
      write (this%iout, fmtnr) this%nodesuser, this%nodes
    end if
    !
    ! -- Array size is now known, so allocate
    call this%allocate_arrays()
    !
    ! -- Fill the nodereduced array with the reduced nodenumber, or
    !    a negative number to indicate it is a pass-through cell, or
    !    a zero to indicate that the cell is excluded from the
    !    solution.
    if (this%nodes < this%nodesuser) then
      node = 1
      noder = 1
      do i = 1, this%nrow
        do j = 1, this%ncol
          if (this%idomain(node) > 0) then
            this%nodereduced(node) = noder
            noder = noder + 1
          elseif (this%idomain(node) < 0) then
            this%nodereduced(node) = -1
          else
            this%nodereduced(node) = 0
          end if
          node = node + 1
        end do
      end do
    end if
    !
    ! -- allocate and fill nodeuser if a reduced grid
    if (this%nodes < this%nodesuser) then
      node = 1
      noder = 1
      do i = 1, this%nrow
        do j = 1, this%ncol
          if (this%idomain(node) > 0) then
            this%nodeuser(noder) = node
            noder = noder + 1
          end if
          node = node + 1
        end do
      end do
    end if
    !
    ! -- fill x,y coordinate arrays
    this%cellx(1) = DHALF * this%delr(1)
    this%celly(this%nrow) = DHALF * this%delc(this%nrow)
    do j = 2, this%ncol
      this%cellx(j) = this%cellx(j - 1) + DHALF * this%delr(j - 1) + &
                      DHALF * this%delr(j)
    end do
    ! -- row number increases in negative y direction:
    do i = this%nrow - 1, 1, -1
      this%celly(i) = this%celly(i + 1) + DHALF * this%delc(i + 1) + &
                      DHALF * this%delc(i)
    end do
    !
    ! -- Move bottom into bot, and calculate area
    node = 0
    do i = 1, this%nrow
      do j = 1, this%ncol
        node = node + 1
        noder = node
        if (this%nodes < this%nodesuser) noder = this%nodereduced(node)
        if (noder <= 0) cycle
        this%bot(noder) = this%bottom(j, i)
        this%area(noder) = this%delr(j) * this%delc(i)
        this%xc(noder) = this%cellx(j)
        this%yc(noder) = this%celly(i)
      end do
    end do
    !
    ! -- create and fill the connections object
    nrsize = 0
    if (this%nodes < this%nodesuser) nrsize = this%nodes
    allocate (this%con)
    call this%con%disconnections(this%name_model, this%nodes, &
                                 this%ncol, this%nrow, 1, &
                                 nrsize, this%delr, this%delc, &
                                 this%top, this%bot, this%nodereduced, &
                                 this%nodeuser)
    this%nja = this%con%nja
    this%njas = this%con%njas
    !
  end subroutine grid_finalize

  !> @brief Write a binary grid file
  !<
  subroutine write_grb(this, icelltype)
    ! -- modules
    use OpenSpecModule, only: access, form
    ! -- dummy
    class(Dis2dType) :: this
    integer(I4B), dimension(:), intent(in) :: icelltype
    ! -- local
    integer(I4B) :: iunit, ntxt
    integer(I4B), parameter :: lentxt = 100
    character(len=50) :: txthdr
    character(len=lentxt) :: txt
    character(len=LINELENGTH) :: fname
    character(len=*), parameter :: fmtgrdsave = &
      "(4X,'BINARY GRID INFORMATION WILL BE WRITTEN TO:', &
       &/,6X,'UNIT NUMBER: ', I0,/,6X, 'FILE NAME: ', A)"
    !
    ! -- Initialize
    ntxt = 14
    !
    ! -- Open the file
    fname = trim(this%output_fname)
    iunit = getunit()
    write (this%iout, fmtgrdsave) iunit, trim(adjustl(fname))
    call openfile(iunit, this%iout, trim(adjustl(fname)), 'DATA(BINARY)', &
                  form, access, 'REPLACE')
    !
    ! -- write header information
    write (txthdr, '(a)') 'GRID DIS2D'
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
    write (txt, '(3a, i0)') 'NROW ', 'INTEGER ', 'NDIM 0 # ', this%nrow
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'NCOL ', 'INTEGER ', 'NDIM 0 # ', this%ncol
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'NJA ', 'INTEGER ', 'NDIM 0 # ', this%nja
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, 1pg24.15)') 'XORIGIN ', 'DOUBLE ', 'NDIM 0 # ', this%xorigin
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, 1pg24.15)') 'YORIGIN ', 'DOUBLE ', 'NDIM 0 # ', this%yorigin
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, 1pg24.15)') 'ANGROT ', 'DOUBLE ', 'NDIM 0 # ', this%angrot
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'DELR ', 'DOUBLE ', 'NDIM 1 ', this%ncol
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'DELC ', 'DOUBLE ', 'NDIM 1 ', this%nrow
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'BOTM ', 'DOUBLE ', 'NDIM 1 ', this%nodesuser
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
    write (iunit) this%nrow ! nrow
    write (iunit) this%ncol ! ncol
    write (iunit) this%nja ! nja
    write (iunit) this%xorigin ! xorigin
    write (iunit) this%yorigin ! yorigin
    write (iunit) this%angrot ! angrot
    write (iunit) this%delr ! delr
    write (iunit) this%delc ! delc
    write (iunit) this%bottom ! bottom
    write (iunit) this%con%iausr ! iausr
    write (iunit) this%con%jausr ! jausr
    write (iunit) this%idomain ! idomain
    write (iunit) icelltype ! icelltype
    !
    ! -- Close the file
    close (iunit)
    !
  end subroutine write_grb

  !> @brief Convert a user nodenumber to a string (nodenumber) or (k,i,j)
  !<
  subroutine nodeu_to_string(this, nodeu, str)
    ! -- dummy
    class(Dis2dType) :: this
    integer(I4B), intent(in) :: nodeu
    character(len=*), intent(inout) :: str
    ! -- local
    integer(I4B) :: i, j, k
    character(len=10) :: istr, jstr
    !
    call get_ijk(nodeu, this%nrow, this%ncol, 1, i, j, k)
    write (istr, '(i10)') i
    write (jstr, '(i10)') j
    str = '('//trim(adjustl(istr))//','// &
          trim(adjustl(jstr))//')'
    !
  end subroutine nodeu_to_string

  !> @brief Convert a user nodenumber to an array (nodenumber) or (i,j)
  !<
  subroutine nodeu_to_array(this, nodeu, arr)
    ! -- dummy
    class(Dis2dType) :: this
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
        ') is not equal to the discretization dimension (', this%ndim, ')'
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- get k, i, j
    call get_ijk(nodeu, this%nrow, this%ncol, 1, i, j, k)
    !
    ! -- fill array
    arr(1) = i
    arr(2) = j
    !
  end subroutine nodeu_to_array

  !> @brief Get reduced node number from user node number
  !<
  function get_nodenumber_idx1(this, nodeu, icheck) result(nodenumber)
    ! -- return
    integer(I4B) :: nodenumber
    ! -- dummy
    class(Dis2dType), intent(in) :: this
    integer(I4B), intent(in) :: nodeu
    integer(I4B), intent(in) :: icheck
    !
    ! -- check the node number if requested
    if (icheck /= 0) then
      !
      ! -- If within valid range, convert to reduced nodenumber
      if (nodeu < 1 .or. nodeu > this%nodesuser) then
        write (errmsg, '(a,i0,a)') &
          'Node number (', nodeu, &
          ') less than 1 or greater than the number of nodes.'
        call store_error(errmsg)
        nodenumber = 0
      else
        nodenumber = nodeu
        if (this%nodes < this%nodesuser) nodenumber = this%nodereduced(nodeu)
      end if
    else
      nodenumber = nodeu
      if (this%nodes < this%nodesuser) nodenumber = this%nodereduced(nodeu)
    end if
    !
  end function get_nodenumber_idx1

  !> @brief Get reduced node number from layer, row and column indices
  !<
  function get_nodenumber_idx2(this, k, j, icheck) result(nodenumber)
    ! -- return
    integer(I4B) :: nodenumber
    ! -- dummy
    class(Dis2dType), intent(in) :: this
    integer(I4B), intent(in) :: k, j
    integer(I4B), intent(in) :: icheck
    ! -- local
    integer(I4B) :: nodeu, i
    ! formats
    character(len=*), parameter :: fmterr = &
           "('Error in structured-grid cell indices: row = ',i0,&
           &', column = ',i0)"
    !
    i = k
    nodeu = get_node(1, i, j, 1, this%nrow, this%ncol)
    if (nodeu < 1) then
      write (errmsg, fmterr) i, j
      call store_error(errmsg, terminate=.TRUE.)
    end if
    nodenumber = nodeu
    if (this%nodes < this%nodesuser) nodenumber = this%nodereduced(nodeu)
    !
    ! -- check the node number if requested
    if (icheck /= 0) then
      !
      if (i < 1 .or. i > this%nrow) &
        call store_error('Row less than one or greater than nrow')
      if (j < 1 .or. j > this%ncol) &
        call store_error('Column less than one or greater than ncol')
      !
      ! -- Error if outside of range
      if (nodeu < 1 .or. nodeu > this%nodesuser) then
        write (errmsg, '(a,i0,a)') &
          'Node number (', nodeu, ')less than 1 or greater than nodes.'
        call store_error(errmsg)
      end if
    end if
    !
  end function get_nodenumber_idx2

  !> @brief Allocate and initialize scalar variables
  !<
  subroutine allocate_scalars(this, name_model, input_mempath)
    ! -- dummy
    class(Dis2dType) :: this
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    !
    ! -- Allocate parent scalars
    call this%DisBaseType%allocate_scalars(name_model, input_mempath)
    !
    ! -- Allocate
    call mem_allocate(this%nrow, 'NROW', this%memoryPath)
    call mem_allocate(this%ncol, 'NCOL', this%memoryPath)
    !
    ! -- Initialize
    this%nrow = 0
    this%ncol = 0
    this%ndim = 2
    !
  end subroutine allocate_scalars

  !> @brief Allocate and initialize arrays
  !<
  subroutine allocate_arrays(this)
    ! -- dummy
    class(Dis2dType) :: this
    !
    ! -- Allocate arrays in DisBaseType (mshape, top, bot, area)
    call this%DisBaseType%allocate_arrays()
    !
    ! -- Allocate arrays for DisType
    if (this%nodes < this%nodesuser) then
      call mem_allocate(this%nodeuser, this%nodes, 'NODEUSER', this%memoryPath)
      call mem_allocate(this%nodereduced, this%nodesuser, 'NODEREDUCED', &
                        this%memoryPath)
    else
      call mem_allocate(this%nodeuser, 1, 'NODEUSER', this%memoryPath)
      call mem_allocate(this%nodereduced, 1, 'NODEREDUCED', this%memoryPath)
    end if
    !
    ! -- Initialize
    this%mshape(1) = this%nrow
    this%mshape(2) = this%ncol
    !
  end subroutine allocate_arrays

  !> @brief Convert a string to a user nodenumber
  !!
  !! Parse layer, row and column and return user nodenumber.
  !! If flag_string is present and true, the first token may be
  !! non-numeric (e.g. boundary name). In this case, return -2.
  !<
  function nodeu_from_string(this, lloc, istart, istop, in, iout, line, &
                             flag_string, allow_zero) result(nodeu)
    ! -- dummy
    class(Dis2dType) :: this
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
    integer(I4B) :: i, j, nrow, ncol
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
    nrow = this%mshape(1)
    ncol = this%mshape(2)
    !
    call urword(line, lloc, istart, istop, 2, i, r, iout, in)
    call urword(line, lloc, istart, istop, 2, j, r, iout, in)
    !
    if (i == 0 .and. j == 0) then
      if (present(allow_zero)) then
        if (allow_zero) then
          nodeu = 0
          return
        end if
      end if
    end if
    !
    errmsg = ""
    !
    if (i < 1 .or. i > nrow) then
      write (errmsg, '(a,1x,a,i0,a)') &
        trim(adjustl(errmsg)), 'Row number in list (', i, &
        ') is outside of the grid.'
    end if
    if (j < 1 .or. j > ncol) then
      write (errmsg, '(a,1x,a,i0,a)') &
        trim(adjustl(errmsg)), 'Column number in list (', j, &
        ') is outside of the grid.'
    end if
    !
    nodeu = get_node(1, i, j, 1, nrow, ncol)
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
    class(Dis2dType) :: this
    character(len=*), intent(inout) :: cellid
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    logical, optional, intent(in) :: flag_string
    logical, optional, intent(in) :: allow_zero
    ! -- local
    integer(I4B) :: lloclocal, istart, istop, ndum, n
    integer(I4B) :: i, j, nrow, ncol
    integer(I4B) :: istat
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
    nrow = this%mshape(1)
    ncol = this%mshape(2)
    !
    lloclocal = 1
    call urword(cellid, lloclocal, istart, istop, 2, i, r, iout, inunit)
    call urword(cellid, lloclocal, istart, istop, 2, j, r, iout, inunit)
    !
    if (i == 0 .and. j == 0) then
      if (present(allow_zero)) then
        if (allow_zero) then
          nodeu = 0
          return
        end if
      end if
    end if
    !
    errmsg = ""
    !
    if (i < 1 .or. i > nrow) then
      write (errmsg, '(a,1x,a,i0,a)') &
        trim(adjustl(errmsg)), 'Row number in list (', i, &
        ') is outside of the grid.'
    end if
    if (j < 1 .or. j > ncol) then
      write (errmsg, '(a,1x,a,i0,a)') &
        trim(adjustl(errmsg)), 'Column number in list (', j, &
        ') is outside of the grid.'
    end if
    !
    nodeu = get_node(1, i, j, 1, nrow, ncol)
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

  !> @brief Indicates whether the grid discretization supports layers
  !<
  logical function supports_layers(this)
    ! -- dummy
    class(Dis2dType) :: this
    !
    supports_layers = .true.
    !
  end function supports_layers

  !> @brief Return number of cells per layer (nrow * ncol)
  !<
  function get_ncpl(this)
    integer(I4B) :: get_ncpl
    class(Dis2dType) :: this
    get_ncpl = this%nrow * this%ncol
  end function get_ncpl

  !> @brief Get normal vector components between the cell and a given neighbor
  !!
  !! The normal points outward from the shared face between noden and nodem.
  !<
  subroutine connection_normal(this, noden, nodem, ihc, xcomp, ycomp, zcomp, &
                               ipos)
    ! -- dummy
    class(Dis2dType) :: this
    integer(I4B), intent(in) :: noden !< cell (reduced nn)
    integer(I4B), intent(in) :: nodem !< neighbor (reduced nn)
    integer(I4B), intent(in) :: ihc !< horizontal connection flag
    real(DP), intent(inout) :: xcomp
    real(DP), intent(inout) :: ycomp
    real(DP), intent(inout) :: zcomp
    integer(I4B), intent(in) :: ipos
    ! -- local
    integer(I4B) :: nodeu1, i1, j1, k1
    integer(I4B) :: nodeu2, i2, j2, k2
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
      xcomp = DZERO
      ycomp = DZERO
      zcomp = DZERO
      nodeu1 = this%get_nodeuser(noden)
      nodeu2 = this%get_nodeuser(nodem)
      call get_ijk(nodeu1, this%nrow, this%ncol, 1, i1, j1, k1)
      call get_ijk(nodeu2, this%nrow, this%ncol, 1, i2, j2, k2)
      if (i2 < i1) then ! back
        ycomp = DONE
      elseif (j2 < j1) then ! left
        xcomp = -DONE
      elseif (j2 > j1) then ! right
        xcomp = DONE
      else ! front
        ycomp = -DONE
      end if
      !
    end if
    !
  end subroutine connection_normal

  !> @brief Get unit vector components between the cell and a given neighbor
  !!
  !< For dis2d there is no z component
  subroutine connection_vector(this, noden, nodem, nozee, satn, satm, ihc, &
                               xcomp, ycomp, zcomp, conlen)
    ! modules
    use DisvGeom, only: line_unit_vector
    ! dummy
    class(Dis2dType) :: this
    integer(I4B), intent(in) :: noden !< cell (reduced nn)
    integer(I4B), intent(in) :: nodem !< neighbor (reduced nn)
    logical, intent(in) :: nozee !< not used for dis2d
    real(DP), intent(in) :: satn !< not used for dis2d
    real(DP), intent(in) :: satm !< not used for dis2d
    integer(I4B), intent(in) :: ihc !< not used for dis2d (always horizontal)
    real(DP), intent(inout) :: xcomp !< x component of the connection vector
    real(DP), intent(inout) :: ycomp !< y component of the connection vector
    real(DP), intent(inout) :: zcomp !< z component, which is always zero
    real(DP), intent(inout) :: conlen !< calculated connection length
    ! local
    real(DP) :: z1, z2
    real(DP) :: x1, y1, x2, y2
    real(DP) :: ds
    integer(I4B) :: i1, i2, j1, j2, k1, k2
    integer(I4B) :: nodeu1, nodeu2, ipos

    ! Calculate vector components
    z1 = DZERO
    z2 = DZERO
    ipos = this%con%getjaindex(noden, nodem)
    ds = this%con%cl1(this%con%jas(ipos)) + this%con%cl2(this%con%jas(ipos))
    nodeu1 = this%get_nodeuser(noden)
    nodeu2 = this%get_nodeuser(nodem)
    call get_ijk(nodeu1, this%nrow, this%ncol, 1, i1, j1, k1)
    call get_ijk(nodeu2, this%nrow, this%ncol, 1, i2, j2, k2)
    x1 = DZERO
    x2 = DZERO
    y1 = DZERO
    y2 = DZERO
    if (i2 < i1) then ! back
      y2 = ds
    elseif (j2 < j1) then ! left
      x2 = -ds
    elseif (j2 > j1) then ! right
      x2 = ds
    else ! front
      y2 = -ds
    end if
    call line_unit_vector(x1, y1, z1, x2, y2, z2, xcomp, ycomp, zcomp, conlen)
  end subroutine connection_vector

  !> @brief Get the discretization type
  !<
  subroutine get_dis_type(this, dis_type)
    ! -- dummy
    class(Dis2dType), intent(in) :: this
    character(len=*), intent(out) :: dis_type
    !
    dis_type = "DIS2D"
    !
  end subroutine get_dis_type

  !> @brief Get the discretization type enumeration
  function get_dis_enum(this) result(dis_enum)
    use ConstantsModule, only: DIS2D
    class(Dis2dType), intent(in) :: this
    integer(I4B) :: dis_enum
    dis_enum = DIS2D
  end function get_dis_enum

  !> @brief Get a 2D array of polygon vertices, listed in
  !!
  !! clockwise order beginning with the lower left corner
  !<
  subroutine get_polyverts(this, ic, polyverts, closed)
    ! -- dummy
    class(Dis2dType), intent(inout) :: this
    integer(I4B), intent(in) :: ic !< cell number (reduced)
    real(DP), allocatable, intent(out) :: polyverts(:, :) !< polygon vertices (column-major indexing)
    logical(LGP), intent(in), optional :: closed !< whether to close the polygon, duplicating a vertex
    ! -- local
    integer(I4B) :: icu, nverts, irow, jcol, klay
    real(DP) :: cellx, celly, dxhalf, dyhalf
    logical(LGP) :: lclosed
    !
    nverts = 4
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
    icu = this%get_nodeuser(ic)
    call get_ijk(icu, this%nrow, this%ncol, 1, irow, jcol, klay)
    cellx = this%cellx(jcol)
    celly = this%celly(irow)
    dxhalf = DHALF * this%delr(jcol)
    dyhalf = DHALF * this%delc(irow)
    polyverts(:, 1) = (/cellx - dxhalf, celly - dyhalf/) ! SW
    polyverts(:, 2) = (/cellx - dxhalf, celly + dyhalf/) ! NW
    polyverts(:, 3) = (/cellx + dxhalf, celly + dyhalf/) ! NE
    polyverts(:, 4) = (/cellx + dxhalf, celly - dyhalf/) ! SE
    !
    ! close if enabled
    if (lclosed) &
      polyverts(:, nverts + 1) = polyverts(:, 1)
    !
  end subroutine

  !> @brief Read an integer array
  !< TODO: REMOVE?
  subroutine read_int_array(this, line, lloc, istart, istop, iout, in, &
                            iarray, aname)
    ! -- dummy
    class(Dis2dType), intent(inout) :: this
    character(len=*), intent(inout) :: line
    integer(I4B), intent(inout) :: lloc
    integer(I4B), intent(inout) :: istart
    integer(I4B), intent(inout) :: istop
    integer(I4B), intent(in) :: in
    integer(I4B), intent(in) :: iout
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: iarray
    character(len=*), intent(in) :: aname
    ! -- local

    !
  end subroutine read_int_array

  !> @brief Read a double precision array
  !< TODO: REMOVE?
  subroutine read_dbl_array(this, line, lloc, istart, istop, iout, in, &
                            darray, aname)
    ! -- dummy
    class(Dis2dType), intent(inout) :: this
    character(len=*), intent(inout) :: line
    integer(I4B), intent(inout) :: lloc
    integer(I4B), intent(inout) :: istart
    integer(I4B), intent(inout) :: istop
    integer(I4B), intent(in) :: in
    integer(I4B), intent(in) :: iout
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: darray
    character(len=*), intent(in) :: aname
    !
  end subroutine read_dbl_array

  !> @brief Read a 2d double array into col icolbnd of darray
  !!
  !! For cells that are outside of the active domain,
  !! do not copy the array value into darray.
  !< TODO: REMOVE?
  subroutine read_layer_array(this, nodelist, darray, ncolbnd, maxbnd, &
                              icolbnd, aname, inunit, iout)
    ! -- dummy
    class(Dis2dType) :: this
    integer(I4B), intent(in) :: maxbnd
    integer(I4B), dimension(maxbnd) :: nodelist
    integer(I4B), intent(in) :: ncolbnd
    real(DP), dimension(ncolbnd, maxbnd), intent(inout) :: darray
    integer(I4B), intent(in) :: icolbnd
    character(len=*), intent(in) :: aname
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    !
  end subroutine read_layer_array

  !> @brief Record a double precision array.
  !!
  !! The array is written to a formatted or unformatted external file
  !! depending on the arguments.
  !<
  subroutine record_array(this, darray, iout, iprint, idataun, aname, &
                          cdatafmp, nvaluesp, nwidthp, editdesc, dinact)
    ! -- dummy
    class(Dis2dType), intent(inout) :: this
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
    nrow = this%mshape(1)
    ncol = this%mshape(2)
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
    class(Dis2dType) :: this
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
    nrow = this%mshape(1)
    ncol = this%mshape(2)
    !
    ! -- Use ubdsv06 to write list header
    call ubdsv06(kstp, kper, text, textmodel, textpackage, dstmodel, dstpackage, &
                 ibdchn, naux, auxtxt, ncol, nrow, nlay, &
                 nlist, iout, delt, pertim, totim)
    !
  end subroutine record_srcdst_list_header

  !> @brief Convert an integer array (layer numbers) to nodelist
  !<
  subroutine nlarray_to_nodelist(this, darray, nodelist, maxbnd, nbound, aname)
    ! -- dummy
    class(Dis2dType) :: this
    integer(I4B), intent(in) :: maxbnd
    integer(I4B), dimension(:), pointer, contiguous :: darray
    integer(I4B), dimension(maxbnd), intent(inout) :: nodelist
    integer(I4B), intent(inout) :: nbound
    character(len=*), intent(in) :: aname
    ! -- local
    integer(I4B) :: il, ir, ic, ncol, nrow, nlay, nval, nodeu, noder, ipos, ierr
    !
    ! -- set variables
    nlay = 1
    nrow = this%mshape(1)
    ncol = this%mshape(2)
    !
    if (this%ndim > 1) then
      !
      nval = ncol * nrow
      !
      ! -- Copy array into nodelist
      ipos = 1
      ierr = 0
      do ir = 1, nrow
        do ic = 1, ncol
          nodeu = get_node(1, ir, ic, nlay, nrow, ncol)
          il = darray(nodeu)
          if (il < 1 .or. il > nlay) then
            write (errmsg, '(a,1x,i0)') 'Invalid layer number:', il
            call store_error(errmsg, terminate=.TRUE.)
          end if
          nodeu = get_node(il, ir, ic, nlay, nrow, ncol)
          noder = this%get_nodenumber(nodeu, 0)
          if (ipos > maxbnd) then
            ierr = ipos
          else
            nodelist(ipos) = noder
          end if
          ipos = ipos + 1
        end do
      end do
      !
      ! -- Check for errors
      nbound = ipos - 1
      if (ierr > 0) then
        write (errmsg, '(a,1x,i0)') &
          'MAXBOUND dimension is too small.'// &
          'INCREASE MAXBOUND TO:', ierr
        call store_error(errmsg, terminate=.TRUE.)
      end if
      !
      ! -- If nbound < maxbnd, then initialize nodelist to zero in this range
      if (nbound < maxbnd) then
        do ipos = nbound + 1, maxbnd
          nodelist(ipos) = 0
        end do
      end if
      !
    else
      !
      ! -- For unstructured, read nodelist directly, then check node numbers
      nodelist = darray
      do noder = 1, maxbnd
        if (noder < 1 .or. noder > this%nodes) then
          write (errmsg, '(a,1x,i0)') 'Invalid node number:', noder
          call store_error(errmsg, terminate=.TRUE.)
        end if
      end do
      nbound = maxbnd
      !
    end if
    !
  end subroutine nlarray_to_nodelist

end module Dis2dModule
