!> @brief This module contains the CHF junction (JNC) package
!!
!! Explicit channel-flow junctions.  A junction is a DISV1D vertex
!! that is touched by two or more reaches.  Each such vertex is promoted to an
!! explicit node carrying its own stage unknown and a continuity equation
!! (only inflows/outflows, no storage).  Junctions are auto-detected from the
!! DISV1D connectivity and wired into the CHF model (chf_df appends the junction
!! equations; chf_ac/mc/fc/cq assemble them), following the MAW/advanced-package
!! convention for a package that adds its own matrix rows (ioffset + per-junction
!! row index; connections added to the solution sparse matrix, not to dis%con).
!!
!! The package is enabled by listing JNC6 in the CHF name file: present means
!! junctions are detected and activated; absent means none of the junction logic
!! runs.  The input file carries only an OPTIONS block (currently SAVE_FLOWS);
!! junctions themselves are always auto-detected from the grid, never listed.
!<
module ChfJncModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DHALF
  use SimModule, only: store_error, store_error_filename
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use MatrixBaseModule, only: MatrixBaseType
  use Disv1dModule, only: Disv1dType
  use SwfDfwModule, only: SwfDfwType

  implicit none
  private
  public :: ChfJncType, chf_jnc_cr

  !> @brief CHF junction (JNC) package type
  !!
  !! Junction data are stored struct-of-arrays, indexed by junction k.  The
  !! ragged junction->reach connectivity uses a CSR layout: iajunc(k)..iajunc(k+1)
  !! -1 index into the flat jareach / reach_nodes arrays.
  !<
  type, extends(NumericalPackageType) :: ChfJncType
    integer(I4B), pointer :: njunctions => null() !< total number of junctions
    integer(I4B), pointer :: nconn => null() !< total junction-reach connections (size of flat arrays)
    integer(I4B), pointer :: ioffset => null() !< offset of junction rows in the model (row = dis%nodes + ioffset + k)
    type(Disv1dType), pointer :: disv1d => null() !< concrete DISV1D grid (narrowed in set_pointers)
    type(SwfDfwType), pointer :: dfw => null() !< DFW package, for reach-junction conductance
    integer(I4B), dimension(:), pointer, contiguous :: ivert => null() !< DISV1D vertex number, per junction
    integer(I4B), dimension(:), pointer, contiguous :: nreaches => null() !< number of connected reaches, per junction
    real(DP), dimension(:), pointer, contiguous :: stage => null() !< current stage, per junction
    integer(I4B), dimension(:), pointer, contiguous :: iajunc => null() !< CSR row pointer (size njunctions + 1) into jareach / reach_nodes
    integer(I4B), dimension(:), pointer, contiguous :: jareach => null() !< flat connected reach (cell) numbers (size nconn)
    integer(I4B), dimension(:), pointer, contiguous :: reach_nodes => null() !< flat connected reach reduced node numbers (size nconn)
    ! matrix position caches (filled in jnc_mc, MAW convention)
    integer(I4B), dimension(:), pointer, contiguous :: idxjdglo => null() !< global position of each junction-row diagonal (size njunctions)
    integer(I4B), dimension(:), pointer, contiguous :: idxjoffdglo => null() !< global position of each junction-row off-diagonal to its reach (size nconn)
    integer(I4B), dimension(:), pointer, contiguous :: idxrdglo => null() !< global position of each reach-row diagonal touched by a junction (size nconn)
    integer(I4B), dimension(:), pointer, contiguous :: idxroffdglo => null() !< global position of each reach-row off-diagonal to its junction (size nconn)
  contains
    procedure :: jnc_df
    procedure :: jnc_ac
    procedure :: jnc_mc
    procedure :: jnc_fc
    procedure :: jnc_cq
    procedure :: jnc_ot
    procedure :: jnc_da
    procedure :: allocate_scalars
    procedure :: set_pointers
    procedure :: detect_junctions
    procedure :: mask_reach_connections
    procedure, private :: source_options
    procedure, private :: log_options
    procedure, private :: qcalc_rj
  end type ChfJncType

contains

  !> @brief Create a new JNC package object
  !!
  !! Called unconditionally by the CHF model; inunit > 0 only when JNC6 is listed
  !! in the name file.  When inunit <= 0 the package stays dormant: the model
  !! skips detect_junctions / the junction rows, so no junction logic runs.  When
  !! active, the OPTIONS block is sourced from the input context at mempath.
  !<
  subroutine chf_jnc_cr(jncobj, name_model, mempath, inunit, iout)
    ! dummy
    type(ChfJncType), pointer :: jncobj !< object to create
    character(len=*), intent(in) :: name_model !< name of the CHF model
    character(len=*), intent(in) :: mempath !< input context mem path
    integer(I4B), intent(in) :: inunit !< package input file unit (0 if JNC6 absent)
    integer(I4B), intent(in) :: iout !< unit number for model listing output

    ! create the object
    allocate (jncobj)

    ! create name and memory path
    call jncobj%set_names(1, name_model, 'JNC', 'JNC', mempath)

    ! allocate scalars
    call jncobj%allocate_scalars()

    ! set variables
    jncobj%inunit = inunit
    jncobj%iout = iout

    ! source the OPTIONS block when the package is active
    if (inunit > 0) then
      call jncobj%source_options()
    end if

  end subroutine chf_jnc_cr

  !> @brief Source input options for the package
  !!
  !! Reads the JNC OPTIONS block from the input context.  Only SAVE_FLOWS is
  !! supported; it maps to the standard ipakcb flag on NumericalPackageType.
  !<
  subroutine source_options(this)
    ! modules
    use MemoryManagerExtModule, only: mem_set_value
    use ChfJncInputModule, only: ChfJncParamFoundType
    ! dummy
    class(ChfJncType) :: this !< this instance
    ! local
    type(ChfJncParamFoundType) :: found

    call mem_set_value(this%ipakcb, 'IPAKCB', this%input_mempath, found%ipakcb)

    if (found%ipakcb) then
      this%ipakcb = -1
    end if

    call this%log_options(found)

  end subroutine source_options

  !> @brief Log the options found in the input
  !<
  subroutine log_options(this, found)
    ! modules
    use ChfJncInputModule, only: ChfJncParamFoundType
    ! dummy
    class(ChfJncType) :: this !< this instance
    type(ChfJncParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'PROCESSING JNC OPTIONS'

    if (found%ipakcb) then
      write (this%iout, '(4x,a)') &
        'JUNCTION FLOWS WILL BE SAVED TO BUDGET FILE SPECIFIED IN OUTPUT CONTROL'
    end if

    write (this%iout, '(1x,a)') 'END OF JNC OPTIONS'

  end subroutine log_options

  !> @brief Allocate scalar members
  !<
  subroutine allocate_scalars(this)
    ! modules
    use MemoryManagerModule, only: mem_allocate
    ! dummy
    class(ChfJncType) :: this

    ! allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()

    ! allocate scalars
    call mem_allocate(this%njunctions, 'NJUNCTIONS', this%memoryPath)
    call mem_allocate(this%nconn, 'NCONN', this%memoryPath)
    call mem_allocate(this%ioffset, 'IOFFSET', this%memoryPath)

    ! allocate junction arrays at zero size; detect_junctions sizes them
    call mem_allocate(this%ivert, 0, 'IVERT', this%memoryPath)
    call mem_allocate(this%nreaches, 0, 'NREACHES', this%memoryPath)
    call mem_allocate(this%stage, 0, 'STAGE', this%memoryPath)
    call mem_allocate(this%iajunc, 0, 'IAJUNC', this%memoryPath)
    call mem_allocate(this%jareach, 0, 'JAREACH', this%memoryPath)
    call mem_allocate(this%reach_nodes, 0, 'REACH_NODES', this%memoryPath)

    ! allocate matrix position caches at zero size; jnc_mc sizes them
    call mem_allocate(this%idxjdglo, 0, 'IDXJDGLO', this%memoryPath)
    call mem_allocate(this%idxjoffdglo, 0, 'IDXJOFFDGLO', this%memoryPath)
    call mem_allocate(this%idxrdglo, 0, 'IDXRDGLO', this%memoryPath)
    call mem_allocate(this%idxroffdglo, 0, 'IDXROFFDGLO', this%memoryPath)

    ! initialize scalars
    this%njunctions = 0
    this%nconn = 0
    this%ioffset = 0

  end subroutine allocate_scalars

  !> @brief Store pointers to objects needed by the package
  !!
  !! Narrow the model's polymorphic dis to a concrete DISV1D pointer once and
  !! cache it, so later routines need no type guard (as in the SWF ZDG package).
  !! Also cache the DFW package, used for DFW-consistent reach-junction
  !! conductance in jnc_fc / jnc_cq.
  !<
  subroutine set_pointers(this, dis, dfw)
    ! dummy
    class(ChfJncType) :: this !< this instance
    class(DisBaseType), pointer, intent(in) :: dis !< model discretization (expected DISV1D)
    type(SwfDfwType), pointer, intent(in) :: dfw !< the DFW package

    select type (dis)
    type is (Disv1dType)
      this%disv1d => dis
    class default
      call store_error('JNC package requires a DISV1D discretization.')
      call store_error_filename(this%input_fname)
    end select

    this%dfw => dfw

  end subroutine set_pointers

  !> @brief Detect junctions from DISV1D connectivity
  !!
  !! Promotes every vertex touched by two or more reaches to a junction (a vertex
  !! touched by one reach is a boundary endpoint).  Fills the junction arrays and
  !! njunctions / nconn; the solution-row mapping is assigned later during wiring.
  !! Requires set_pointers first.
  !!
  !! A vertex->cell map is needed but the grid does not persist one (its builder,
  !! disv1dconnections_verts, computes an equivalent map only as a transient
  !! local), so this routine rebuilds it locally from iavert/javert.
  !!
  !! LIMITATION: reduced DISV1D grids (dis%nodes < dis%nodesuser) are not yet
  !! supported; detection runs in user-node space while the model solves in
  !! reduced space, so this routine guards against them rather than implementing
  !! unverified mapping.
  !! TODO(jnc-future): support reduced DISV1D grids - translate connected reaches
  !! through nodereduced, exclude inactive reaches, base Nr on dis%nodes, and add
  !! a dedicated reduced-grid test; then remove the guard below.
  !<
  subroutine detect_junctions(this)
    ! modules
    use MemoryManagerModule, only: mem_reallocate
    ! dummy
    class(ChfJncType) :: this !< this instance
    ! local
    integer(I4B), dimension(:), pointer, contiguous :: iavert
    integer(I4B), dimension(:), pointer, contiguous :: javert
    integer(I4B) :: nodesuser !< number of user reaches
    integer(I4B) :: nvert !< number of vertices
    integer(I4B), allocatable :: vertcount(:) !< reaches touching each vertex
    integer(I4B), allocatable :: vfill(:) !< per-vertex fill cursor
    integer(I4B), allocatable :: iavertcells(:) !< CSR index: vertex -> cells
    integer(I4B), allocatable :: javertcells(:) !< CSR data: cells touching each vertex
    integer(I4B) :: n !< reach (cell) index
    integer(I4B) :: iv !< vertex index
    integer(I4B) :: j !< javert position
    integer(I4B) :: ipos !< javertcells position
    integer(I4B) :: nvc !< running total of vertex->cell connections
    integer(I4B) :: k !< junction index
    integer(I4B) :: njunc !< number of detected junctions
    integer(I4B) :: njconn !< number of junction-reach connections
    integer(I4B) :: ja !< running fill cursor into flat junction arrays

    ! the concrete DISV1D grid must have been cached by set_pointers
    if (.not. associated(this%disv1d)) then
      call store_error('JNC detect_junctions called before set_pointers; &
                       &no DISV1D grid is associated.')
      call store_error_filename(this%input_fname)
      return
    end if

    ! guard: reduced DISV1D grids are not yet supported (see routine header)
    if (this%disv1d%nodes < this%disv1d%nodesuser) then
      call store_error('JNC junctions are not yet supported on a reduced &
                       &DISV1D grid (IDOMAIN excludes one or more reaches). &
                       &Remove IDOMAIN exclusions to use channel junctions.')
      call store_error_filename(this%input_fname)
      return
    end if

    iavert => this%disv1d%iavert
    javert => this%disv1d%javert
    nodesuser = this%disv1d%nodesuser
    nvert = this%disv1d%nvert

    ! first pass: count how many reaches touch each vertex
    allocate (vertcount(nvert))
    do iv = 1, nvert
      vertcount(iv) = 0
    end do
    do n = 1, nodesuser
      do j = iavert(n), iavert(n + 1) - 1
        iv = javert(j)
        vertcount(iv) = vertcount(iv) + 1
      end do
    end do

    ! build CSR index array iavertcells (vertex -> cells)
    allocate (iavertcells(nvert + 1))
    nvc = 0
    iavertcells(1) = 1
    do iv = 1, nvert
      nvc = nvc + vertcount(iv)
      iavertcells(iv + 1) = iavertcells(iv) + vertcount(iv)
    end do

    ! second pass: fill the javertcells data array
    allocate (javertcells(nvc))
    allocate (vfill(nvert))
    do iv = 1, nvert
      vfill(iv) = iavertcells(iv)
    end do
    do n = 1, nodesuser
      do j = iavert(n), iavert(n + 1) - 1
        iv = javert(j)
        javertcells(vfill(iv)) = n
        vfill(iv) = vfill(iv) + 1
      end do
    end do

    ! count junctions (vertices touched by >= 2 reaches) and total connections
    njunc = 0
    njconn = 0
    do iv = 1, nvert
      if (vertcount(iv) >= 2) then
        njunc = njunc + 1
        njconn = njconn + vertcount(iv)
      end if
    end do

    ! size the package arrays now that the counts are known
    this%njunctions = njunc
    this%nconn = njconn
    call mem_reallocate(this%ivert, njunc, 'IVERT', this%memoryPath)
    call mem_reallocate(this%nreaches, njunc, 'NREACHES', this%memoryPath)
    call mem_reallocate(this%stage, njunc, 'STAGE', this%memoryPath)
    call mem_reallocate(this%iajunc, njunc + 1, 'IAJUNC', this%memoryPath)
    call mem_reallocate(this%jareach, njconn, 'JAREACH', this%memoryPath)
    call mem_reallocate(this%reach_nodes, njconn, 'REACH_NODES', this%memoryPath)

    ! fill the junction arrays and the CSR connectivity (iajunc / jareach)
    k = 0
    ja = 1
    this%iajunc(1) = 1
    do iv = 1, nvert
      if (vertcount(iv) < 2) cycle
      k = k + 1
      this%ivert(k) = iv
      this%nreaches(k) = vertcount(iv)
      this%stage(k) = DZERO
      do ipos = iavertcells(iv), iavertcells(iv + 1) - 1
        this%jareach(ja) = javertcells(ipos)
        ! reduced node == user reach number here (reduced grids are guarded out)
        this%reach_nodes(ja) = javertcells(ipos)
        ja = ja + 1
      end do
      this%iajunc(k + 1) = ja
    end do

    ! report detected junctions to the listing file
    if (this%iout > 0) then
      write (this%iout, '(1x,a,i0,a)') &
        'CHF JNC: detected ', this%njunctions, ' junction(s).'
      do k = 1, this%njunctions
        write (this%iout, '(4x,a,i0,a,i0,a)') &
          'Junction at vertex ', this%ivert(k), &
          ' connects ', this%nreaches(k), ' reaches.'
      end do
    end if

    ! cleanup local maps
    deallocate (vertcount)
    deallocate (vfill)
    deallocate (iavertcells)
    deallocate (javertcells)

  end subroutine detect_junctions

  !> @brief Mask the direct reach-reach connections at each junction
  !!
  !! Two reaches that share a junction vertex are directly connected in the
  !! DISV1D connectivity (dis%con), so DFW would compute a pairwise flow between
  !! them.  With explicit junctions that pairwise coupling is replaced by routing
  !! through the junction node, so the direct reach-reach connections must be
  !! masked out (both directions) to avoid double-counting the exchange.  DFW
  !! honors con%mask in both its formulate and flow-calculation loops.
  !<
  subroutine mask_reach_connections(this)
    ! dummy
    class(ChfJncType) :: this !< this instance
    ! local
    integer(I4B) :: k !< junction index
    integer(I4B) :: ia1, ia2 !< flat-array cursors for the two reaches
    integer(I4B) :: ireach, jreach !< the two reach node numbers
    integer(I4B) :: ipos !< connection position in dis%con

    do k = 1, this%njunctions
      ! for every unordered pair of reaches meeting at this junction
      do ia1 = this%iajunc(k), this%iajunc(k + 1) - 1
        ireach = this%reach_nodes(ia1)
        do ia2 = ia1 + 1, this%iajunc(k + 1) - 1
          jreach = this%reach_nodes(ia2)

          ! mask ireach -> jreach
          ipos = this%disv1d%con%getjaindex(ireach, jreach)
          if (ipos > 0) call this%disv1d%con%set_mask(ipos, 0)

          ! mask jreach -> ireach
          ipos = this%disv1d%con%getjaindex(jreach, ireach)
          if (ipos > 0) call this%disv1d%con%set_mask(ipos, 0)
        end do
      end do
    end do

  end subroutine mask_reach_connections

  !> @brief Define the junction package
  !<
  subroutine jnc_df(this)
    ! dummy
    class(ChfJncType) :: this !< this instance
    ! Options are sourced in chf_jnc_cr and junctions are auto-detected in
    ! chf_df, so there is nothing to define here.
  end subroutine jnc_df

  !> @brief Add junction rows and reach-junction connections to the sparse matrix
  !!
  !! Each junction k owns global row (moffset + nreach + ioffset + k).  For every
  !! connected reach it adds the symmetric pair of off-diagonals coupling the
  !! junction row and the reach row, plus the junction-row diagonal.  Follows the
  !! MAW convention; the model's own ia/ja (dis%con) is left untouched.
  !<
  subroutine jnc_ac(this, moffset, nreach, sparse)
    use SparseModule, only: sparsematrix
    ! dummy
    class(ChfJncType) :: this !< this instance
    integer(I4B), intent(in) :: moffset !< model offset in the solution
    integer(I4B), intent(in) :: nreach !< number of reach equations (dis%nodes)
    type(sparsematrix), intent(inout) :: sparse !< sparse matrix structure
    ! local
    integer(I4B) :: k !< junction index
    integer(I4B) :: ipos !< position in flat connection arrays
    integer(I4B) :: jglo !< global junction row number
    integer(I4B) :: rglo !< global reach row number

    do k = 1, this%njunctions
      jglo = moffset + nreach + this%ioffset + k
      call sparse%addconnection(jglo, jglo, 1)
      do ipos = this%iajunc(k), this%iajunc(k + 1) - 1
        rglo = moffset + this%reach_nodes(ipos)
        call sparse%addconnection(jglo, rglo, 1)
        call sparse%addconnection(rglo, jglo, 1)
      end do
    end do

  end subroutine jnc_ac

  !> @brief Cache solution-matrix positions for the junction connections
  !!
  !! Mirrors MAW: find and store the global A positions this package writes to,
  !! so jnc_fc can fill coefficients by position.
  !<
  subroutine jnc_mc(this, moffset, nreach, matrix_sln)
    use MemoryManagerModule, only: mem_reallocate
    ! dummy
    class(ChfJncType) :: this !< this instance
    integer(I4B), intent(in) :: moffset !< model offset in the solution
    integer(I4B), intent(in) :: nreach !< number of reach equations (dis%nodes)
    class(MatrixBaseType), pointer :: matrix_sln !< solution matrix
    ! local
    integer(I4B) :: k !< junction index
    integer(I4B) :: ipos !< position in flat connection arrays
    integer(I4B) :: jglo !< global junction row number
    integer(I4B) :: rglo !< global reach row number

    ! size the position caches
    call mem_reallocate(this%idxjdglo, this%njunctions, 'IDXJDGLO', &
                        this%memoryPath)
    call mem_reallocate(this%idxjoffdglo, this%nconn, 'IDXJOFFDGLO', &
                        this%memoryPath)
    call mem_reallocate(this%idxrdglo, this%nconn, 'IDXRDGLO', this%memoryPath)
    call mem_reallocate(this%idxroffdglo, this%nconn, 'IDXROFFDGLO', &
                        this%memoryPath)

    do k = 1, this%njunctions
      jglo = moffset + nreach + this%ioffset + k
      this%idxjdglo(k) = matrix_sln%get_position_diag(jglo)
      do ipos = this%iajunc(k), this%iajunc(k + 1) - 1
        rglo = moffset + this%reach_nodes(ipos)
        ! junction-row entries: diagonal (cached above) and off-diagonal to reach
        this%idxjoffdglo(ipos) = matrix_sln%get_position(jglo, rglo)
        ! reach-row entries: its own diagonal and off-diagonal back to junction
        this%idxrdglo(ipos) = matrix_sln%get_position_diag(rglo)
        this%idxroffdglo(ipos) = matrix_sln%get_position(rglo, jglo)
      end do
    end do

  end subroutine jnc_mc

  !> @brief Flow from reach r into junction k
  !!
  !! q = C_rk (h_r - h_k), where C_rk is the half-cell conductance between the
  !! reach center and the junction.  Positive q is flow from the reach into the
  !! junction.
  !!
  !! The junction sits at the reach endpoint and is a zero-length point, so the
  !! only half-cell distance is the reach's own half-length (dx) and there is no
  !! second half-cell to harmonically average with.  This routine prepares the
  !! half-cell inputs (upstream-weighted depth, friction gradient over the half-
  !! length, smoothed depth, flow width) the same way the DFW reach-reach path
  !! does, then calls the shared DFW half-cell conductance kernel get_cond_n so
  !! the Manning conductance physics lives in one place.  The reach-junction
  !! gradient uses the half-length dx (abs(stage_r - stage_j) / dx), not a
  !! center-to-center distance, because the junction is at the face.
  !<
  function qcalc_rj(this, ireach, stage_r, stage_j) result(q)
    ! modules
    use ConstantsModule, only: DPREC
    use SmoothingModule, only: sQuadratic
    ! dummy
    class(ChfJncType) :: this !< this instance
    integer(I4B), intent(in) :: ireach !< reduced reach node number
    real(DP), intent(in) :: stage_r !< stage in the reach
    real(DP), intent(in) :: stage_j !< stage at the junction
    ! return
    real(DP) :: q
    ! local
    real(DP) :: dx
    real(DP) :: depth
    real(DP) :: dhds
    real(DP) :: width
    real(DP) :: width_dummy
    real(DP) :: range = 1.d-6
    real(DP) :: dydx
    real(DP) :: smooth_factor
    real(DP) :: cond

    ! reach half-length to the junction endpoint
    dx = DHALF * this%disv1d%length(ireach)

    cond = DZERO
    if (dx > DPREC) then

      ! TODO(jnc-depth): the half-cell conductance between a reach and a
      ! junction always computes depth using the reach's own bottom,
      ! regardless of which side is upstream.  But when the junction is upstream
      ! (stage_j > stage_r) the junction should lend its own depth, as the
      ! upstream reach would in the reach-reach case (SwfDfwType%get_cond,
      ! icentral==0).  The junction currently has no bottom of its own.
      ! To match the reach-reach case, a junction connecting two reaches,
      ! A and B, must be assigned the bottom of reach B when calculating the
      ! half-cell conductance for reach A (and vice versa).  Generalizing to a
      ! multi-reach junction, the junction must be assigned the bottom of one of
      ! the reaches other than A.  Note that this implies that a static junction
      ! bottom (e.g., user-specified or min/max over all connected reaches)
      ! cannot work in general.
      !   One possibility is a "most-upstream" rule: when the junction is
      ! upstream of reach A, use the bottom of whichever other reach connected
      ! to this junction currently has the highest stage.  (For an unconverged
      ! solution, a junction might be upstream of all of its reaches, in which
      ! case this could more properly be called a "least-downstream" rule.) This
      ! selection is always well-defined (no empty-set/fallback case needed)
      ! and reduces exactly to the selection in the reach-reach case.
      !   This selection makes q_rk jump discontinuously when the most-
      ! upstream reach changes, if bottoms differ.  This is not a new class of
      ! issue, since get_cond's icentral==0 is subject to the same kind of
      ! discontinuity.  But a multi-reach junction may flip more often, so
      ! Newton convergence at junctions might be worse in practice.  If so,
      ! consider a "smooth" assignment (e.g. weighted by relative stage).

      ! depth with upstream weighting (use the higher of reach/junction stage)
      if (stage_r >= stage_j) then
        depth = stage_r - this%disv1d%bot(ireach)
      else
        depth = stage_j - this%disv1d%bot(ireach)
      end if

      ! friction gradient between reach center and junction over the half-length
      dhds = abs(stage_r - stage_j) / dx

      ! smoothed depth that goes to zero over the specified range
      call sQuadratic(depth, range, dydx, smooth_factor)
      depth = depth * smooth_factor

      ! reach flow width (the junction endpoint uses the reach width)
      call this%disv1d%get_flow_width(ireach, ireach, 1, width, width_dummy)

      ! half-cell conductance from the reach center to the junction, using the
      ! shared DFW Manning conductance kernel
      cond = this%dfw%get_cond_n(ireach, depth, dx, width, dhds)

    end if

    q = cond * (stage_r - stage_j)

  end function qcalc_rj

  !> @brief Formulate junction continuity rows (Newton)
  !!
  !! For each junction k, assembles continuity sum_r q_rk = 0 where q_rk is the
  !! flow from reach r into the junction.  The conductance is stage-dependent, so
  !! derivatives are formed by numerical perturbation, matching dfw_qnm_fc_nr.
  !! The reach equations receive the opposite contribution (-q_rk).
  !<
  subroutine jnc_fc(this, matrix_sln, rhs, stage, nreach)
    ! modules
    use MathUtilModule, only: get_perturbation
    ! dummy
    class(ChfJncType) :: this !< this instance
    class(MatrixBaseType), pointer :: matrix_sln !< solution matrix
    real(DP), intent(inout), dimension(:) :: rhs !< solution right-hand side (model-local)
    real(DP), intent(inout), dimension(:) :: stage !< solution dependent-variable (model-local)
    integer(I4B), intent(in) :: nreach !< number of reach equations (dis%nodes)
    ! local
    integer(I4B) :: k !< junction index
    integer(I4B) :: ipos !< position in flat connection arrays
    integer(I4B) :: jeq !< model-local junction equation number
    integer(I4B) :: req !< model-local reach equation number
    real(DP) :: hj !< junction stage
    real(DP) :: hr !< reach stage
    real(DP) :: q !< flow from reach into junction
    real(DP) :: qeps
    real(DP) :: eps
    real(DP) :: dqdhr !< d(q)/d(reach stage)
    real(DP) :: dqdhj !< d(q)/d(junction stage)

    do k = 1, this%njunctions
      jeq = nreach + this%ioffset + k
      hj = stage(jeq)
      do ipos = this%iajunc(k), this%iajunc(k + 1) - 1
        req = this%reach_nodes(ipos)
        hr = stage(req)

        ! flow from reach into junction and its stage derivatives
        q = this%qcalc_rj(req, hr, hj)
        eps = get_perturbation(hr)
        qeps = this%qcalc_rj(req, hr + eps, hj)
        dqdhr = (qeps - q) / eps
        eps = get_perturbation(hj)
        qeps = this%qcalc_rj(req, hr, hj + eps)
        dqdhj = (qeps - q) / eps

        ! junction continuity row: + q_rk (inflow positive), Newton-linearized
        rhs(jeq) = rhs(jeq) - q + dqdhr * hr + dqdhj * hj
        call matrix_sln%add_value_pos(this%idxjoffdglo(ipos), dqdhr)
        call matrix_sln%add_value_pos(this%idxjdglo(k), dqdhj)

        ! reach row receives the opposite contribution: - q_rk
        rhs(req) = rhs(req) + q - dqdhr * hr - dqdhj * hj
        call matrix_sln%add_value_pos(this%idxrdglo(ipos), -dqdhr)
        call matrix_sln%add_value_pos(this%idxroffdglo(ipos), -dqdhj)
      end do
    end do

  end subroutine jnc_fc

  !> @brief Capture junction stage and add reach-junction flows to flowja
  !!
  !! Each reach-junction connection is off dis%con, so csr_diagsum never sees it.
  !! Following the MAW convention (bnd_cq_simrate), add each connection's flow
  !! onto the reach cell's flowja diagonal so the reach residual closes and the
  !! model budget balances.  Positive q_rk is flow from the reach into the
  !! junction, i.e. flow leaving the reach, so -q_rk is added to the reach
  !! diagonal (flowja is positive into a cell).  The junction node itself has
  !! zero storage and net-zero through-flow, so it needs no diagonal term.
  !<
  subroutine jnc_cq(this, stage, flowja, nreach)
    ! dummy
    class(ChfJncType) :: this !< this instance
    real(DP), intent(in), dimension(:) :: stage !< solution dependent-variable (model-local)
    real(DP), intent(inout), dimension(:) :: flowja !< model connection flows (CSR)
    integer(I4B), intent(in) :: nreach !< number of reach equations (dis%nodes)
    ! local
    integer(I4B) :: k !< junction index
    integer(I4B) :: ipos !< position in flat connection arrays
    integer(I4B) :: jeq !< model-local junction equation number
    integer(I4B) :: req !< model-local reach equation number
    integer(I4B) :: idiag !< reach diagonal position in flowja
    real(DP) :: hj !< junction stage
    real(DP) :: q !< flow from reach into junction

    do k = 1, this%njunctions
      jeq = nreach + this%ioffset + k
      hj = stage(jeq)
      this%stage(k) = hj
      do ipos = this%iajunc(k), this%iajunc(k + 1) - 1
        req = this%reach_nodes(ipos)
        q = this%qcalc_rj(req, stage(req), hj)
        idiag = this%disv1d%con%ia(req)
        flowja(idiag) = flowja(idiag) - q
      end do
    end do

  end subroutine jnc_cq

  !> @brief Output junction results
  !<
  subroutine jnc_ot(this)
    ! dummy
    class(ChfJncType) :: this !< this instance
    ! TODO(jnc-future): output junction observations
  end subroutine jnc_ot

  !> @brief Deallocate junction storage
  !<
  subroutine jnc_da(this)
    ! modules
    use MemoryManagerModule, only: mem_deallocate
    ! dummy
    class(ChfJncType) :: this !< this instance

    ! deallocate junction arrays
    call mem_deallocate(this%ivert)
    call mem_deallocate(this%nreaches)
    call mem_deallocate(this%stage)
    call mem_deallocate(this%iajunc)
    call mem_deallocate(this%jareach)
    call mem_deallocate(this%reach_nodes)

    ! deallocate matrix position caches
    call mem_deallocate(this%idxjdglo)
    call mem_deallocate(this%idxjoffdglo)
    call mem_deallocate(this%idxrdglo)
    call mem_deallocate(this%idxroffdglo)

    ! nullify borrowed pointers (not owned by this package)
    nullify (this%disv1d)
    nullify (this%dfw)

    ! deallocate scalars
    call mem_deallocate(this%njunctions)
    call mem_deallocate(this%nconn)
    call mem_deallocate(this%ioffset)

    ! deallocate parent
    call this%NumericalPackageType%da()

  end subroutine jnc_da

end module ChfJncModule
