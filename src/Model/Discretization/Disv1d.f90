module Disv1dModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENMEMPATH, LENVARNAME, DZERO, DONE, LINELENGTH
  use SimVariablesModule, only: errmsg, warnmsg
  use MemoryHelperModule, only: create_mem_path
  use MemoryManagerModule, only: mem_allocate
  use SimModule, only: count_errors, store_error, store_warning, &
                       store_error_filename
  use InputOutputModule, only: urword
  use BaseDisModule, only: DisBaseType
  use Disv1dGeom, only: calcdist
  use DisvGeom, only: line_unit_vector

  implicit none

  private
  public :: disv1d_cr
  public :: Disv1dType

  type, extends(DisBaseType) :: Disv1dType
    integer(I4B), pointer :: nvert => null() !< number of x,y vertices
    real(DP), dimension(:), pointer, contiguous :: length => null() !< length of each reach (of size nodesuser)
    real(DP), dimension(:), pointer, contiguous :: width => null() !< reach width
    real(DP), dimension(:), pointer, contiguous :: bottom => null() !< reach bottom elevation
    real(DP), dimension(:, :), pointer, contiguous :: vertices => null() !< cell vertices stored as 2d array with columns of x, y
    real(DP), dimension(:, :), pointer, contiguous :: cellxy => null() !< reach midpoints stored as 2d array with columns of x, y
    real(DP), dimension(:), pointer, contiguous :: fdc => null() !< fdc stored as array
    integer(I4B), dimension(:), pointer, contiguous :: iavert => null() !< cell vertex pointer ia array
    integer(I4B), dimension(:), pointer, contiguous :: javert => null() !< cell vertex pointer ja array
  contains
    procedure :: disv1d_load => disv1d_load
    procedure :: dis_df => disv1d_df
    procedure :: dis_da => disv1d_da
    procedure :: get_dis_type => get_dis_type
    procedure :: get_dis_enum => get_dis_enum
    procedure :: get_flow_width => get_flow_width
    procedure, public :: record_array
    procedure, public :: record_srcdst_list_header
    ! -- private
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: source_options
    procedure :: source_dimensions
    procedure :: source_griddata
    procedure :: source_vertices
    procedure :: source_cell1d
    procedure :: log_options
    procedure :: log_dimensions
    procedure :: log_griddata
    procedure :: define_cellverts
    procedure :: grid_finalize
    !procedure :: connect
    procedure :: create_connections
    procedure :: write_grb
    procedure :: get_nodenumber_idx1
    procedure :: nodeu_to_string
    procedure :: nodeu_from_string
    procedure :: connection_normal
    procedure :: connection_vector

  end type Disv1dType

  !> @brief Simplifies tracking parameters sourced from the input context.
  type DisFoundType
    logical :: length_units = .false.
    logical :: nogrb = .false.
    logical :: xorigin = .false.
    logical :: yorigin = .false.
    logical :: angrot = .false.
    logical :: nodes = .false.
    logical :: nvert = .false.
    logical :: width = .false.
    logical :: bottom = .false.
    logical :: idomain = .false.
    logical :: iv = .false.
    logical :: xv = .false.
    logical :: yv = .false.
    logical :: icell1d = .false.
    logical :: fdc = .false.
    logical :: ncvert = .false.
    logical :: icvert = .false.
  end type DisFoundType

contains

  subroutine disv1d_cr(dis, name_model, input_mempath, inunit, iout)
    ! -- modules
    use KindModule, only: LGP
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy
    class(DisBaseType), pointer :: dis
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    ! -- locals
    type(Disv1dType), pointer :: disnew
    logical(LGP) :: found_fname
    character(len=*), parameter :: fmtheader = &
      "(1X, /1X, 'DISV1D -- DISCRETIZATION BY VERTICES IN 1D PACKAGE,', &
      &' VERSION 1 : 4/2/2024 - INPUT READ FROM MEMPATH: ', A, /)"
    allocate (disnew)
    dis => disnew
    call disnew%allocate_scalars(name_model, input_mempath)
    dis%input_mempath = input_mempath
    dis%inunit = inunit
    dis%iout = iout
    !
    ! -- set name of input file
    call mem_set_value(dis%input_fname, 'INPUT_FNAME', dis%input_mempath, &
                       found_fname)
    !
    ! -- If dis enabled
    if (inunit > 0) then

      ! -- Identify package
      if (iout > 0) then
        write (iout, fmtheader) dis%input_mempath
      end if

    end if
  end subroutine disv1d_cr

  !> @brief Define the discretization
  !<
  subroutine disv1d_df(this)
    ! -- dummy
    class(Disv1dType) :: this
    !
    ! -- Transfer the data from the memory manager into this package object
    if (this%inunit /= 0) then
      call this%disv1d_load()
    end if

    ! finalize the grid
    call this%grid_finalize()

  end subroutine disv1d_df

  !> @brief Get normal vector components between the cell and a given neighbor
  !!
  !! The normal points outward from the shared face between noden and nodem.
  !<
  subroutine connection_normal(this, noden, nodem, ihc, xcomp, ycomp, zcomp, &
                               ipos)
    ! -- dummy
    class(Disv1dType) :: this
    integer(I4B), intent(in) :: noden !< cell (reduced nn)
    integer(I4B), intent(in) :: nodem !< neighbor (reduced nn)
    integer(I4B), intent(in) :: ihc !< horizontal connection flag (not used)
    real(DP), intent(inout) :: xcomp !< x component of outward normal vector
    real(DP), intent(inout) :: ycomp !< y component of outward normal vector
    real(DP), intent(inout) :: zcomp !< z component of outward normal vector
    integer(I4B), intent(in) :: ipos !< connection position
    ! -- local
    real(DP) :: angle, dmult

    ! find from anglex, since anglex is symmetric, need to flip vector
    ! for lower triangle (nodem < noden)
    angle = this%con%anglex(this%con%jas(ipos))
    dmult = DONE
    if (nodem < noden) dmult = -DONE
    xcomp = cos(angle) * dmult
    ycomp = sin(angle) * dmult
    zcomp = DZERO
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
    class(Disv1dType) :: this
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

  !> @brief Get the discretization type (DIS, DIS2D, DISV, DISV1D, DISU)
  subroutine get_dis_type(this, dis_type)
    class(Disv1dType), intent(in) :: this
    character(len=*), intent(out) :: dis_type
    dis_type = "DISV1D"
  end subroutine get_dis_type

  !> @brief Get the discretization type enumeration
  function get_dis_enum(this) result(dis_enum)
    use ConstantsModule, only: DISV1D
    class(Disv1dType), intent(in) :: this
    integer(I4B) :: dis_enum
    dis_enum = DISV1D
  end function get_dis_enum

  !> @brief Allocate scalar variables
  !<
  subroutine allocate_scalars(this, name_model, input_mempath)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use ConstantsModule, only: DONE
    ! -- dummy
    class(Disv1dType) :: this
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
  end subroutine allocate_scalars

  subroutine disv1d_load(this)
    ! -- dummy
    class(Disv1dType) :: this
    ! -- locals
    !
    ! -- source input data
    call this%source_options()
    call this%source_dimensions()
    call this%source_griddata()

    ! If vertices provided by user, read and store vertices
    if (this%nvert > 0) then
      call this%source_vertices()
      call this%source_cell1d()
    end if

  end subroutine disv1d_load

  !> @brief Copy options from IDM into package
  !<
  subroutine source_options(this)
    ! -- modules
    use KindModule, only: LGP
    use MemoryManagerExtModule, only: mem_set_value
    use SimVariablesModule, only: idm_context
    ! -- dummy
    class(Disv1dType) :: this
    ! -- locals
    character(len=LENVARNAME), dimension(3) :: lenunits = &
      &[character(len=LENVARNAME) :: 'FEET', 'METERS', 'CENTIMETERS']
    character(len=LENMEMPATH) :: idmMemoryPath
    type(DisFoundType) :: found
    !
    ! -- set memory path
    idmMemoryPath = create_mem_path(this%name_model, 'DISV1D', idm_context)
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%lenuni, 'LENGTH_UNITS', &
                       idmMemoryPath, lenunits, found%length_units)
    call mem_set_value(this%nogrb, 'NOGRB', &
                       idmMemoryPath, found%nogrb)
    call mem_set_value(this%xorigin, 'XORIGIN', &
                       idmMemoryPath, found%xorigin)
    call mem_set_value(this%yorigin, 'YORIGIN', &
                       idmMemoryPath, found%yorigin)
    call mem_set_value(this%angrot, 'ANGROT', &
                       idmMemoryPath, found%angrot)
    !
    ! -- log values to list file
    if (this%iout > 0) then
      call this%log_options(found)
    end if
  end subroutine source_options

  !> @brief Write user options to list file
  !<
  subroutine log_options(this, found)
    class(Disv1dType) :: this
    type(DisFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'Setting Discretization Options'

    if (found%length_units) then
      write (this%iout, '(4x,a,i0)') 'Model length unit [0=UND, 1=FEET, &
      &2=METERS, 3=CENTIMETERS] set as ', this%lenuni
    end if

    if (found%nogrb) then
      write (this%iout, '(4x,a,i0)') 'Binary grid file [0=GRB, 1=NOGRB] &
        &set as ', this%nogrb
    end if

    if (found%xorigin) then
      write (this%iout, '(4x,a,G0)') 'XORIGIN = ', this%xorigin
    end if

    if (found%yorigin) then
      write (this%iout, '(4x,a,G0)') 'YORIGIN = ', this%yorigin
    end if

    if (found%angrot) then
      write (this%iout, '(4x,a,G0)') 'ANGROT = ', this%angrot
    end if

    write (this%iout, '(1x,a,/)') 'End Setting Discretization Options'

  end subroutine log_options

  !> @brief Copy dimensions from IDM into package
  !<
  subroutine source_dimensions(this)
    use KindModule, only: LGP
    use MemoryManagerExtModule, only: mem_set_value
    use SimVariablesModule, only: idm_context
    ! -- dummy
    class(Disv1dType) :: this
    ! -- locals
    character(len=LENMEMPATH) :: idmMemoryPath
    integer(I4B) :: n
    type(DisFoundType) :: found
    !
    ! -- set memory path
    idmMemoryPath = create_mem_path(this%name_model, 'DISV1D', idm_context)
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%nodes, 'NODES', idmMemoryPath, found%nodes)
    call mem_set_value(this%nvert, 'NVERT', idmMemoryPath, found%nvert)
    !
    ! -- for now assume nodes = nodesuser
    this%nodesuser = this%nodes
    !
    ! -- log simulation values
    if (this%iout > 0) then
      call this%log_dimensions(found)
    end if
    !
    ! -- verify dimensions were set
    if (this%nodesuser < 1) then
      call store_error( &
        'NODES was not specified or was specified incorrectly.')
      call store_error_filename(this%input_fname)
    end if
    if (this%nvert < 1) then
      call store_warning( &
        'NVERT was not specified or was specified as zero.  The &
        &VERTICES and CELL1D blocks will not be read for the DISV1D6 &
        &Package in model '//trim(this%memoryPath)//'.')
    end if
    !
    ! -- Allocate non-reduced vectors for disv1d
    call mem_allocate(this%length, this%nodesuser, &
                      'LENGTH', this%memoryPath)
    call mem_allocate(this%width, this%nodesuser, &
                      'WIDTH', this%memoryPath)
    call mem_allocate(this%bottom, this%nodesuser, &
                      'BOTTOM', this%memoryPath)
    call mem_allocate(this%idomain, this%nodesuser, &
                      'IDOMAIN', this%memoryPath)
    !
    ! -- Allocate vertices array
    if (this%nvert > 0) then
      call mem_allocate(this%vertices, 2, this%nvert, &
                        'VERTICES', this%memoryPath)
      call mem_allocate(this%fdc, this%nodesuser, &
                        'FDC', this%memoryPath)
      call mem_allocate(this%cellxy, 2, this%nodesuser, &
                        'CELLXY', this%memoryPath)
    end if
    !
    ! -- initialize all cells to be active (idomain = 1)
    do n = 1, this%nodesuser
      this%length(n) = DZERO
      this%width(n) = DZERO
      this%bottom(n) = DZERO
      this%idomain(n) = 1
    end do
  end subroutine source_dimensions

  !> @brief Write dimensions to list file
  !<
  subroutine log_dimensions(this, found)
    class(Disv1dType) :: this
    type(DisFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'Setting Discretization Dimensions'

    if (found%nodes) then
      write (this%iout, '(4x,a,i0)') 'NODES = ', this%nodesuser
    end if

    if (found%nvert) then
      write (this%iout, '(4x,a,i0)') 'NVERT = ', this%nvert
    end if

    write (this%iout, '(1x,a,/)') 'End Setting Discretization Dimensions'

  end subroutine log_dimensions

  subroutine source_griddata(this)
    ! modules
    use MemoryManagerExtModule, only: mem_set_value
    use SimVariablesModule, only: idm_context
    ! dummy
    class(Disv1dType) :: this
    ! locals
    character(len=LENMEMPATH) :: idmMemoryPath
    type(DisFoundType) :: found
    ! formats

    ! set memory path
    idmMemoryPath = create_mem_path(this%name_model, 'DISV1D', idm_context)

    call mem_set_value(this%width, 'WIDTH', idmMemoryPath, &
                       found%width)
    call mem_set_value(this%bottom, 'BOTTOM', idmMemoryPath, &
                       found%bottom)
    call mem_set_value(this%idomain, 'IDOMAIN', idmMemoryPath, found%idomain)

    if (.not. found%width) then
      write (errmsg, '(a)') 'Error in GRIDDATA block: WIDTH not found.'
      call store_error(errmsg)
    end if

    if (.not. found%bottom) then
      write (errmsg, '(a)') 'Error in GRIDDATA block: BOTTOM not found.'
      call store_error(errmsg)
    end if

    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if

    ! log simulation values
    if (this%iout > 0) then
      call this%log_griddata(found)
    end if

  end subroutine source_griddata

  !> @brief Write griddata found to list file
  !<
  subroutine log_griddata(this, found)
    class(Disv1dType) :: this
    type(DisFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'Setting Discretization Griddata'

    if (found%width) then
      write (this%iout, '(4x,a)') 'WIDTH set from input file'
    end if

    if (found%bottom) then
      write (this%iout, '(4x,a)') 'BOTTOM set from input file'
    end if

    if (found%idomain) then
      write (this%iout, '(4x,a)') 'IDOMAIN set from input file'
    end if

    write (this%iout, '(1x,a,/)') 'End Setting Discretization Griddata'

  end subroutine log_griddata

  !> @brief Copy vertex information from input data context
  !! to model context
  !<
  subroutine source_vertices(this)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    use MemoryManagerExtModule, only: mem_set_value
    use SimVariablesModule, only: idm_context
    ! -- dummy
    class(Disv1dType) :: this
    ! -- local
    integer(I4B) :: i
    character(len=LENMEMPATH) :: idmMemoryPath
    real(DP), dimension(:), contiguous, pointer :: vert_x => null()
    real(DP), dimension(:), contiguous, pointer :: vert_y => null()
    ! -- formats
    !
    ! -- set memory path
    idmMemoryPath = create_mem_path(this%name_model, 'DISV1D', idm_context)
    !
    ! -- set pointers to memory manager input arrays
    call mem_setptr(vert_x, 'XV', idmMemoryPath)
    call mem_setptr(vert_y, 'YV', idmMemoryPath)
    !
    ! -- set vertices 3d array
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
      write (this%iout, '(1x,a)') 'Setting Discretization Vertices'
      write (this%iout, '(1x,a,/)') 'End setting discretization vertices'
    end if
  end subroutine source_vertices

  !> @brief Copy cell1d information from input data context
  !! to model context
  !<
  subroutine source_cell1d(this)
    ! -- modules
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: mem_setptr
    use MemoryManagerExtModule, only: mem_set_value
    use SimVariablesModule, only: idm_context
    ! -- dummy
    class(Disv1dType) :: this
    ! -- locals
    character(len=LENMEMPATH) :: idmMemoryPath
    integer(I4B), dimension(:), contiguous, pointer :: icell1d => null()
    integer(I4B), dimension(:), contiguous, pointer :: ncvert => null()
    integer(I4B), dimension(:), contiguous, pointer :: icvert => null()
    real(DP), dimension(:), contiguous, pointer :: fdc => null()
    integer(I4B) :: i
    ! -- formats
    !
    ! -- set memory path
    idmMemoryPath = create_mem_path(this%name_model, 'DISV1D', idm_context)
    !
    ! -- set pointers to input path ncvert and icvert
    call mem_setptr(icell1d, 'ICELL1D', idmMemoryPath)
    call mem_setptr(ncvert, 'NCVERT', idmMemoryPath)
    call mem_setptr(icvert, 'ICVERT', idmMemoryPath)
    !
    ! --
    if (associated(icell1d) .and. associated(ncvert) &
        .and. associated(icvert)) then
      call this%define_cellverts(icell1d, ncvert, icvert)
    else
      call store_error('Required cell vertex arrays not found.')
    end if
    !
    ! -- set pointers to cell center arrays
    call mem_setptr(fdc, 'FDC', idmMemoryPath)
    !
    ! -- set fractional distance to cell center
    if (associated(fdc)) then
      do i = 1, this%nodesuser
        this%fdc(i) = fdc(i)
      end do
    else
      call store_error('Required fdc array not found.')
    end if

    ! calculate length from vertices
    call calculate_cell_length(this%vertices, this%iavert, this%javert, &
                               this%length)

    ! calculate cellxy from vertices and fdc
    call calculate_cellxy(this%vertices, this%fdc, this%iavert, &
                          this%javert, this%length, this%cellxy)

    ! log
    if (this%iout > 0) then
      write (this%iout, '(1x,a)') 'Setting Discretization CELL1D'
      write (this%iout, '(1x,a,/)') 'End Setting Discretization CELL1D'
    end if
  end subroutine source_cell1d

  !> @brief Construct the iavert and javert integer vectors which
  !! are compressed sparse row index arrays that relate the vertices
  !! to reaches
  !<
  subroutine define_cellverts(this, icell1d, ncvert, icvert)
    ! -- modules
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(Disv1dType) :: this
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: icell1d
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: ncvert
    integer(I4B), dimension(:), contiguous, pointer, intent(in) :: icvert
    ! -- locals
    type(sparsematrix) :: vert_spm
    integer(I4B) :: i, j, ierr
    integer(I4B) :: icv_idx, startvert, maxnnz = 2
    !
    ! -- initialize sparse matrix
    call vert_spm%init(this%nodesuser, this%nvert, maxnnz)
    !
    ! -- add sparse matrix connections from input memory paths
    icv_idx = 1
    do i = 1, this%nodesuser
      if (icell1d(i) /= i) call store_error('ICELL1D input sequence violation.')
      do j = 1, ncvert(i)
        call vert_spm%addconnection(i, icvert(icv_idx), 0)
        if (j == 1) then
          startvert = icvert(icv_idx)
        end if
        icv_idx = icv_idx + 1
      end do
    end do
    !
    ! -- allocate and fill iavert and javert
    call mem_allocate(this%iavert, this%nodesuser + 1, 'IAVERT', this%memoryPath)
    call mem_allocate(this%javert, vert_spm%nnz, 'JAVERT', this%memoryPath)
    call vert_spm%filliaja(this%iavert, this%javert, ierr)
    call vert_spm%destroy()
  end subroutine define_cellverts

  !> @brief Calculate x, y, coordinates of reach midpoint
  !<
  subroutine calculate_cellxy(vertices, fdc, iavert, javert, length, cellxy)
    ! -- dummy
    real(DP), dimension(:, :), intent(in) :: vertices !< 2d array of vertices with x, y as columns
    real(DP), dimension(:), intent(in) :: fdc !< fractional distance to reach midpoint (normally 0.5)
    integer(I4B), dimension(:), intent(in) :: iavert !< csr mapping of vertices to cell reaches
    integer(I4B), dimension(:), intent(in) :: javert !< csr mapping of vertices to cell reaches
    real(DP), dimension(:), intent(in) :: length !< vector of cell lengths
    real(DP), dimension(:, :), intent(inout) :: cellxy !< 2d array of reach midpoint with x, y as columns
    ! -- local
    integer(I4B) :: nodes !< number of nodes
    integer(I4B) :: n !< node index
    integer(I4B) :: j !< vertex index
    integer(I4B) :: iv0 !< index for line reach start
    integer(I4B) :: iv1 !< index for linen reach end
    integer(I4B) :: ixy !< x, y column index
    real(DP) :: fd0 !< fractional distance to start of this line reach
    real(DP) :: fd1 !< fractional distance to end of this line reach
    real(DP) :: fd !< fractional distance where midpoint (defined by fdc) is located
    real(DP) :: d !< distance

    nodes = size(iavert) - 1
    do n = 1, nodes

      ! find vertices that span midpoint
      iv0 = 0
      iv1 = 0
      fd0 = DZERO
      do j = iavert(n), iavert(n + 1) - 2
        d = calcdist(vertices, javert(j), javert(j + 1))
        fd1 = fd0 + d / length(n)

        ! if true, then we know the midpoint is some fractional distance (fd)
        ! from vertex j to vertex j + 1
        if (fd1 >= fdc(n)) then
          iv0 = javert(j)
          iv1 = javert(j + 1)
          fd = (fdc(n) - fd0) / (fd1 - fd0)
          exit
        end if
        fd0 = fd1
      end do

      ! find x, y position of point on line
      do ixy = 1, 2
        cellxy(ixy, n) = (DONE - fd) * vertices(ixy, iv0) + &
                         fd * vertices(ixy, iv1)
      end do

    end do
  end subroutine calculate_cellxy

  !> @brief Calculate x, y, coordinates of reach midpoint
  !<
  subroutine calculate_cell_length(vertices, iavert, javert, length)
    ! -- dummy
    real(DP), dimension(:, :), intent(in) :: vertices !< 2d array of vertices with x, y as columns
    integer(I4B), dimension(:), intent(in) :: iavert !< csr mapping of vertices to cell reaches
    integer(I4B), dimension(:), intent(in) :: javert !< csr mapping of vertices to cell reaches
    real(DP), dimension(:), intent(inout) :: length !< 2d array of reach midpoint with x, y as columns
    ! -- local
    integer(I4B) :: nodes !< number of nodes
    integer(I4B) :: n !< node index
    integer(I4B) :: j !< vertex index
    real(DP) :: dlen !< length

    nodes = size(iavert) - 1
    do n = 1, nodes

      ! calculate length of this reach
      dlen = DZERO
      do j = iavert(n), iavert(n + 1) - 2
        dlen = dlen + calcdist(vertices, javert(j), javert(j + 1))
      end do
      length(n) = dlen

    end do
  end subroutine calculate_cell_length

  !> @brief Finalize grid construction
  !<
  subroutine grid_finalize(this)
    ! modules
    use SimModule, only: ustop, count_errors, store_error
    use ConstantsModule, only: LINELENGTH, DZERO, DONE
    ! dummy
    class(Disv1dType) :: this
    ! local
    integer(I4B) :: node, noder, k
    ! format
    character(len=*), parameter :: fmtnr = &
      "(/1x, 'The specified IDOMAIN results in a reduced number of cells.',&
      &/1x, 'Number of user nodes: ',I0,&
      &/1X, 'Number of nodes in solution: ', I0, //)"

    ! count active cells
    this%nodes = 0
    do k = 1, this%nodesuser
      if (this%idomain(k) > 0) this%nodes = this%nodes + 1
    end do
    !
    ! Check to make sure nodes is a valid number
    if (this%nodes == 0) then
      call store_error('Model does not have any active nodes.  Make sure &
                       &IDOMAIN has some values greater than zero.')
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
      do k = 1, this%nodesuser
        if (this%idomain(k) > 0) then
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
      do k = 1, this%nodesuser
        if (this%idomain(k) > 0) then
          this%nodeuser(noder) = node
          noder = noder + 1
        end if
        node = node + 1
      end do
    end if

    ! Move bottom into bot and put length into disbase%area
    ! and set x and y center coordinates
    do node = 1, this%nodesuser
      noder = node
      if (this%nodes < this%nodesuser) noder = this%nodereduced(node)
      if (noder <= 0) cycle
      this%bot(noder) = this%bottom(node)
      this%area(noder) = this%length(node)
    end do

    ! create connectivity using vertices and cell1d
    call this%create_connections()
  end subroutine grid_finalize

  subroutine allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(Disv1dType) :: this
    !
    ! -- Allocate arrays in DisBaseType (mshape, top, bot, area)
    ! todo: disbasetype will have memory allocated for unneeded arrays
    call this%DisBaseType%allocate_arrays()
    !
    ! -- Allocate arrays
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
    this%mshape(1) = this%nodesuser
  end subroutine allocate_arrays

  subroutine create_connections(this)
    ! modules
    ! dummy
    class(Disv1dType) :: this
    ! local
    integer(I4B) :: nrsize

    ! create and fill the connections object
    nrsize = 0
    if (this%nodes < this%nodesuser) nrsize = this%nodes

    ! Allocate connections object
    allocate (this%con)

    ! Build connectivity based on vertices
    call this%con%disv1dconnections_verts(this%name_model, this%nodes, &
                                          this%nodesuser, nrsize, this%nvert, &
                                          this%vertices, this%iavert, &
                                          this%javert, this%cellxy, this%fdc, &
                                          this%nodereduced, this%nodeuser, &
                                          this%length)

    this%nja = this%con%nja
    this%njas = this%con%njas

  end subroutine create_connections

  !> @brief Write binary grid file
  !<
  subroutine write_grb(this, icelltype)
    ! -- modules
    use InputOutputModule, only: getunit, openfile
    use OpenSpecModule, only: access, form
    ! -- dummy
    class(Disv1dType) :: this
    integer(I4B), dimension(:), intent(in) :: icelltype
    ! -- local
    integer(I4B) :: i, iunit, ntxt
    integer(I4B), parameter :: lentxt = 100
    character(len=50) :: txthdr
    character(len=lentxt) :: txt
    character(len=LINELENGTH) :: fname
    character(len=*), parameter :: fmtgrdsave = &
      "(4X,'BINARY GRID INFORMATION WILL BE WRITTEN TO:', &
       &/,6X,'UNIT NUMBER: ', I0,/,6X, 'FILE NAME: ', A)"
    !
    ! -- Initialize
    ntxt = 10
    if (this%nvert > 0) ntxt = ntxt + 5
    !
    ! -- Open the file
    fname = trim(this%output_fname)
    iunit = getunit()
    write (this%iout, fmtgrdsave) iunit, trim(adjustl(fname))
    call openfile(iunit, this%iout, trim(adjustl(fname)), 'DATA(BINARY)', &
                  form, access, 'REPLACE')
    !
    ! -- write header information
    write (txthdr, '(a)') 'GRID DISV1D'
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
    write (txt, '(3a, i0)') 'NJA ', 'INTEGER ', 'NDIM 0 # ', this%con%nja
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
    write (txt, '(3a, i0)') 'BOTM ', 'DOUBLE ', 'NDIM 1 ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'IA ', 'INTEGER ', 'NDIM 1 ', this%nodesuser + 1
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'JA ', 'INTEGER ', 'NDIM 1 ', this%con%nja
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'ICELLTYPE ', 'INTEGER ', 'NDIM 1 ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    write (txt, '(3a, i0)') 'IDOMAIN   ', 'INTEGER ', 'NDIM 1 ', this%nodesuser
    txt(lentxt:lentxt) = new_line('a')
    write (iunit) txt
    !
    ! -- if vertices have been read then write additional header information
    if (this%nvert > 0) then
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
    end if
    !
    ! -- write data
    write (iunit) this%nodesuser ! nodes
    write (iunit) this%nja ! nja
    write (iunit) this%xorigin ! xorigin
    write (iunit) this%yorigin ! yorigin
    write (iunit) this%angrot ! angrot
    write (iunit) this%bottom ! botm
    write (iunit) this%con%iausr ! ia
    write (iunit) this%con%jausr ! ja
    write (iunit) icelltype ! icelltype
    write (iunit) this%idomain ! idomain
    !
    ! -- if vertices have been read then write additional data
    if (this%nvert > 0) then
      write (iunit) this%vertices ! vertices
      write (iunit) (this%cellxy(1, i), i=1, this%nodesuser) ! cellx
      write (iunit) (this%cellxy(2, i), i=1, this%nodesuser) ! celly
      write (iunit) this%iavert ! iavert
      write (iunit) this%javert ! javert
    end if
    !
    ! -- Close the file
    close (iunit)
  end subroutine write_grb

  !>
  !! Return a nodenumber from the user specified node number with an
  !! option to perform a check.  This subroutine can be overridden by
  !! child classes to perform mapping to a model node number
  !<
  function get_nodenumber_idx1(this, nodeu, icheck) result(nodenumber)
    class(Disv1dType), intent(in) :: this
    integer(I4B), intent(in) :: nodeu
    integer(I4B), intent(in) :: icheck
    integer(I4B) :: nodenumber
    !
    if (icheck /= 0) then
      if (nodeu < 1 .or. nodeu > this%nodes) then
        write (errmsg, '(a,i10)') &
          'Nodenumber less than 1 or greater than nodes:', nodeu
        call store_error(errmsg)
      end if
    end if
    !
    ! -- set node number based on whether it is reduced or not
    if (this%nodes == this%nodesuser) then
      nodenumber = nodeu
    else
      nodenumber = this%nodereduced(nodeu)
    end if
  end function get_nodenumber_idx1

  subroutine nodeu_to_string(this, nodeu, str)
    ! -- dummy
    class(Disv1dType) :: this
    integer(I4B), intent(in) :: nodeu
    character(len=*), intent(inout) :: str
    ! -- local
    character(len=10) :: nstr
    !
    write (nstr, '(i0)') nodeu
    str = '('//trim(adjustl(nstr))//')'
  end subroutine nodeu_to_string

  !>
  !! nodeu_from_string -- Receive a string and convert the string to a user
  !!   nodenumber.  The model is unstructured; just read user nodenumber.
  !!   If flag_string argument is present and true, the first token in string
  !!   is allowed to be a string (e.g. boundary name). In this case, if a string
  !!   is encountered, return value as -2.
  !<
  function nodeu_from_string(this, lloc, istart, istop, in, iout, line, &
                             flag_string, allow_zero) result(nodeu)
    ! -- dummy
    class(Disv1dType) :: this
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
    call urword(line, lloc, istart, istop, 2, nodeu, r, iout, in)
    !
    if (nodeu == 0) then
      if (present(allow_zero)) then
        if (allow_zero) then
          return
        end if
      end if
    end if
    !
    if (nodeu < 1 .or. nodeu > this%nodesuser) then
      write (errmsg, '(a,i0,a)') &
        "Node number in list (", nodeu, ") is outside of the grid. "// &
        "Cell number cannot be determined in line '"// &
        trim(adjustl(line))//"'."
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if
  end function nodeu_from_string

  subroutine disv1d_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    ! -- dummy
    class(Disv1dType) :: this
    ! -- local
    logical(LGP) :: deallocate_vertices
    !
    ! -- Deallocate idm memory
    call memorystore_remove(this%name_model, 'DISV1D', idm_context)
    !
    ! -- scalars
    deallocate_vertices = (this%nvert > 0)
    call mem_deallocate(this%nvert)
    !
    ! -- arrays
    call mem_deallocate(this%nodeuser)
    call mem_deallocate(this%nodereduced)
    call mem_deallocate(this%length)
    call mem_deallocate(this%width)
    call mem_deallocate(this%bottom)
    call mem_deallocate(this%idomain)
    !
    ! -- cdl hack for arrays for vertices and cell1d blocks
    if (deallocate_vertices) then
      call mem_deallocate(this%vertices)
      call mem_deallocate(this%fdc)
      call mem_deallocate(this%cellxy)
      call mem_deallocate(this%iavert)
      call mem_deallocate(this%javert)
    end if
    !
    ! -- DisBaseType deallocate
    call this%DisBaseType%dis_da()
  end subroutine disv1d_da

  !> @brief Record a double precision array
  !!
  !!   Record a double precision array.  The array will be
  !!   printed to an external file and/or written to an unformatted external file
  !!   depending on the argument specifications.
  !<
  subroutine record_array(this, darray, iout, iprint, idataun, aname, &
                          cdatafmp, nvaluesp, nwidthp, editdesc, dinact)
    ! -- modules
    use TdisModule, only: kstp, kper, pertim, totim, delt
    use InputOutputModule, only: ulasav, ulaprufw, ubdsv1
    ! -- dummy
    class(Disv1dType), intent(inout) :: this
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: darray !< double precision array to record
    integer(I4B), intent(in) :: iout !< unit number for ascii output
    integer(I4B), intent(in) :: iprint !< flag indicating whether or not to print the array
    integer(I4B), intent(in) :: idataun !< unit number to which the array will be written in binary
    character(len=*), intent(in) :: aname !< text descriptor of the array
    character(len=*), intent(in) :: cdatafmp ! fortran format for writing the array
    integer(I4B), intent(in) :: nvaluesp !< number of values per line for printing
    integer(I4B), intent(in) :: nwidthp !< width of the number for printing
    character(len=*), intent(in) :: editdesc !< format type (I, G, F, S, E)
    real(DP), intent(in) :: dinact !< double precision value to use for cells that are excluded from model domain
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
  end subroutine record_array

  !> @brief Record list header using ubdsv06
  !<
  subroutine record_srcdst_list_header(this, text, textmodel, textpackage, &
                                       dstmodel, dstpackage, naux, auxtxt, &
                                       ibdchn, nlist, iout)
    ! -- module
    use TdisModule, only: kstp, kper, pertim, totim, delt
    use InputOutputModule, only: ubdsv06
    ! -- dummy
    class(Disv1dType) :: this
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
  end subroutine record_srcdst_list_header

  !> @ brief Calculate the flow width between two cells
  !!
  !! This should only be called for connections with IHC > 0.
  !! Routine is needed, so it can be overridden by the linear
  !! network discretization, which allows for a separate flow
  !< width for each cell.
  subroutine get_flow_width(this, n, m, idx_conn, width_n, width_m)
    ! dummy
    class(Disv1dType) :: this
    integer(I4B), intent(in) :: n !< cell node number
    integer(I4B), intent(in) :: m !< cell node number
    integer(I4B), intent(in) :: idx_conn !< connection index
    real(DP), intent(out) :: width_n !< flow width for cell n
    real(DP), intent(out) :: width_m !< flow width for cell m

    ! For disv1d case, width_n and width_m can be different
    width_n = this%width(n)
    width_m = this%width(m)

  end subroutine get_flow_width

end module Disv1dModule
