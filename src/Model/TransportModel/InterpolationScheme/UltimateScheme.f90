module UltimateSchemeModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE, DTWO, DTHREE, DSIX, DHALF, DPREC, &
                             LINELENGTH
  use InterpolationSchemeInterfaceModule, only: InterpolationSchemeInterface, &
                                                CoefficientsType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use GeomUtilModule, only: get_ijk
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg

  implicit none
  private

  public :: UltimateSchemeType

  !> Slack allowed on the Courant number before the time step is rejected.  The
  !! stability limit is one, but ATS arrives at the time step by dividing the
  !! pore volume by the flow rate and the check here multiplies them back
  !! together, so a time step that is exactly at the limit can land just above
  !! it.  The slack is far smaller than any time step length that would matter
  !! for stability.
  real(DP), parameter :: COURANT_SLACK = 1.0e-8_DP

  !> @brief Third-order TVD (ULTIMATE) interpolation scheme.
  !!
  !! Explicit, third-order accurate advection scheme for structured (DIS) grids,
  !! ported from the ULTIMATE scheme of MT3DMS.  The face concentration is the
  !! QUICKEST reconstruction, built from the two cells straddling the face plus
  !! the second cell upstream, and then clamped by Leonard's universal limiter to
  !! enforce the TVD condition.
  !!
  !! Because the reconstruction traces characteristics backward over the time
  !! step, the face value depends on the time step length and on the
  !! concentration at the beginning of the time step.  The scheme is therefore
  !! explicit: it contributes nothing to the coefficient matrix and everything to
  !! the right-hand side.  Face values are computed once per time step by
  !! prepare() and reused for every outer iteration.
  !!
  !! Being explicit, the scheme is stable only for a Courant number of one or
  !! less, which is why the ADV Package requires the adaptive time stepping (ATS)
  !! Package and the ATS_PERCEL option when this scheme is selected.  Those
  !! requirements do not fully guarantee the constraint is met, so prepare() also
  !! verifies the Courant number using the time step actually taken.
  !!
  !! The reconstruction is the full multidimensional one: the face value is the
  !! average of the field over the volume swept through the face during the time
  !! step, which carries transverse advection, transverse curvature and twist
  !! (mixed derivative) terms in addition to the terms normal to the face.
  !<
  type, extends(InterpolationSchemeInterface) :: UltimateSchemeType
    private
    class(DisBaseType), pointer :: dis => null() !< pointer to model discretization
    type(TspFmiType), pointer :: fmi => null() !< pointer to fmi object
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    real(DP), dimension(:), pointer, contiguous :: porosity => null() !< pointer to mobile domain porosity
    real(DP), dimension(:), pointer, contiguous :: retardation => null() !< pointer to retardation factor
    real(DP), dimension(:), pointer, contiguous :: phiold => null() !< field at the start of the time step
    integer(I4B), dimension(:), allocatable :: nodefar !< (nja) neighbor of n opposite m, or 0
    real(DP), dimension(:), allocatable :: distfar !< (nja) center distance from n to nodefar
    real(DP), dimension(:), allocatable :: phiface !< (njas) explicit face value
    integer(I4B), dimension(:), allocatable :: faceaxis !< (njas) grid axis the face is normal to
    integer(I4B), dimension(:, :), allocatable :: nbr !< (6, nodes) neighbor in each axis direction, or 0
    integer(I4B), dimension(:, :), allocatable :: nbrpos !< (6, nodes) connection position of nbr
    real(DP), dimension(:, :), allocatable :: distnbr !< (6, nodes) center distance to nbr
    real(DP), dimension(:, :), allocatable :: span !< (3, nodes) distance across both neighbors
    real(DP), dimension(:, :), allocatable :: width !< (3, nodes) cell width along each axis
    real(DP), dimension(:, :), allocatable :: courant_axis !< (3, nodes) signed cell Courant number
    real(DP), dimension(:, :), allocatable :: grad_axis !< (3, nodes) centered first derivative
    real(DP), dimension(:, :), allocatable :: curv_axis !< (3, nodes) centered second derivative
    real(DP), dimension(:, :), allocatable :: cross_axis !< (3, nodes) mixed derivative in the normal plane
  contains
    procedure :: compute
    procedure :: set_field
    procedure :: prepare

    procedure, private :: build_stencil
    procedure, private :: check_courant
    procedure, private :: calc_derivatives
    procedure, private :: diagonal_node
    procedure, private :: face_value
  end type UltimateSchemeType

  interface UltimateSchemeType
    module procedure constructor
  end interface UltimateSchemeType

contains

  !> @brief Create a new ULTIMATE interpolation scheme
  !<
  function constructor(dis, fmi, ibound, porosity, retardation) &
    result(interpolation_scheme)
    ! -- return
    type(UltimateSchemeType) :: interpolation_scheme
    ! -- dummy
    class(DisBaseType), pointer, intent(in) :: dis !< model discretization
    type(TspFmiType), pointer, intent(in) :: fmi !< flow model interface
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: ibound !< model ibound
    real(DP), dimension(:), pointer, contiguous, intent(in) :: porosity !< mobile domain porosity
    real(DP), dimension(:), pointer, contiguous, intent(in) :: retardation !< retardation factor

    interpolation_scheme%dis => dis
    interpolation_scheme%fmi => fmi
    interpolation_scheme%ibound => ibound
    interpolation_scheme%porosity => porosity
    interpolation_scheme%retardation => retardation

    allocate (interpolation_scheme%nodefar(dis%con%nja))
    allocate (interpolation_scheme%distfar(dis%con%nja))
    allocate (interpolation_scheme%phiface(dis%njas))
    allocate (interpolation_scheme%faceaxis(dis%njas))
    allocate (interpolation_scheme%nbr(6, dis%nodes))
    allocate (interpolation_scheme%nbrpos(6, dis%nodes))
    allocate (interpolation_scheme%distnbr(6, dis%nodes))
    allocate (interpolation_scheme%span(3, dis%nodes))
    allocate (interpolation_scheme%width(3, dis%nodes))
    allocate (interpolation_scheme%courant_axis(3, dis%nodes))
    allocate (interpolation_scheme%grad_axis(3, dis%nodes))
    allocate (interpolation_scheme%curv_axis(3, dis%nodes))
    allocate (interpolation_scheme%cross_axis(3, dis%nodes))

    call interpolation_scheme%build_stencil()
  end function constructor

  !> @brief Map the grid onto the axis-aligned stencils the scheme needs
  !!
  !! Two mappings are built, both of which depend only on the grid and are
  !! therefore built once.
  !!
  !! The first resolves each connection onto a grid axis and records, for every
  !! cell, its neighbor in each of the six axis directions along with the
  !! distances involved.  The transverse terms differentiate sideways from a
  !! face, which the connection lists alone cannot express.
  !!
  !! The second is the second upstream cell.  For each connection from n to m,
  !! it is the neighbor of n in the direction opposite m, which is the cell on
  !! the far side of n along the same grid axis.  That cell is the second
  !! upstream cell whenever flow is from n to m.  Taking it from the axis map
  !! rather than from arithmetic on cell indices is what makes vertical
  !! pass-through cells work: a pass-through is skipped by the connection list,
  !! so the cells either side of it are more than one layer apart and reflecting
  !! one through the other in index space lands on the wrong cell.  A value of
  !! zero means no such cell exists, either because the connection is at the
  !! edge of the grid or because the cell was removed by idomain, in which case
  !! the scheme reverts to first-order upwinding.
  !<
  subroutine build_stencil(this)
    ! -- dummy
    class(UltimateSchemeType) :: this
    ! -- local
    integer(I4B) :: n, m, nfar
    integer(I4B) :: ipos, isympos
    integer(I4B) :: nlay, nrow, ncol
    integer(I4B) :: nodeun, nodeum
    integer(I4B) :: irown, icoln, ilayn
    integer(I4B) :: irowm, icolm, ilaym
    integer(I4B) :: iaxis, idir, idirfar, idirm, idirp
    real(DP) :: cl_own

    this%nodefar(:) = 0
    this%distfar(:) = DZERO
    this%faceaxis(:) = 0
    this%nbr(:, :) = 0
    this%nbrpos(:, :) = 0
    this%distnbr(:, :) = DZERO
    this%span(:, :) = DZERO
    this%width(:, :) = DZERO

    nlay = this%dis%mshape(1)
    nrow = this%dis%mshape(2)
    ncol = this%dis%mshape(3)

    do n = 1, this%dis%nodes
      nodeun = this%dis%get_nodeuser(n)
      call get_ijk(nodeun, nrow, ncol, nlay, irown, icoln, ilayn)
      !
      ! -- Resolve every connection of n onto a grid axis and a direction along
      !    it, numbered so that odd directions decrease the cell index and even
      !    directions increase it.  Axis 1 is the column direction, axis 2 the
      !    row direction and axis 3 the layer direction.
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        nodeum = this%dis%get_nodeuser(m)
        call get_ijk(nodeum, nrow, ncol, nlay, irowm, icolm, ilaym)
        isympos = this%dis%con%jas(ipos)
        !
        if (icolm /= icoln) then
          iaxis = 1
          idir = 2 * iaxis - 1
          if (icolm > icoln) idir = 2 * iaxis
        else if (irowm /= irown) then
          iaxis = 2
          idir = 2 * iaxis - 1
          if (irowm > irown) idir = 2 * iaxis
        else if (ilaym /= ilayn) then
          iaxis = 3
          idir = 2 * iaxis - 1
          if (ilaym > ilayn) idir = 2 * iaxis
        else
          cycle
        end if
        !
        if (n < m) then
          cl_own = this%dis%con%cl1(isympos)
        else
          cl_own = this%dis%con%cl2(isympos)
        end if
        this%faceaxis(isympos) = iaxis
        this%nbr(idir, n) = m
        this%nbrpos(idir, n) = ipos
        this%distnbr(idir, n) = this%dis%con%cl1(isympos) + &
                                this%dis%con%cl2(isympos)
        !
        ! -- On a structured grid the connection length is half the cell width,
        !    so either connection along an axis gives the width along that axis
        this%width(iaxis, n) = DTWO * cl_own
      end do
      !
      ! -- The second upstream cell for a connection is the neighbor of n
      !    opposite the cell the connection leads to.  The axis map for n is
      !    complete at this point, so a vertical pass-through cell, which leaves
      !    the cells either side of it more than one layer apart, needs no
      !    special case here.
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        idir = 0
        do idirfar = 1, 6
          if (this%nbr(idirfar, n) == m) then
            idir = idirfar
            exit
          end if
        end do
        if (idir == 0) cycle ! not an axis-aligned connection
        !
        ! -- step to the opposite direction along the same axis
        if (mod(idir, 2) == 1) then
          idirfar = idir + 1
        else
          idirfar = idir - 1
        end if
        !
        nfar = this%nbr(idirfar, n)
        if (nfar == 0) cycle ! no cell on the far side of n
        this%nodefar(ipos) = nfar
        this%distfar(ipos) = this%distnbr(idirfar, n)
      end do
      !
      ! -- Distance across both neighbors, the interval a centered derivative
      !    at n is taken over.  Zero where the cell lacks a neighbor on either
      !    side, which switches off the derivatives along that axis.
      do iaxis = 1, 3
        idirm = 2 * iaxis - 1
        idirp = 2 * iaxis
        if (this%nbr(idirm, n) == 0 .or. this%nbr(idirp, n) == 0) cycle
        this%span(iaxis, n) = this%distnbr(idirm, n) + this%distnbr(idirp, n)
      end do
    end do
  end subroutine build_stencil

  !> @brief Set the scalar field for which interpolation will be computed
  !!
  !! Not used by this scheme.  The ULTIMATE face values are explicit, so they do
  !! not depend on the field being iterated on.  They are calculated by prepare()
  !! at the start of the time step from the field at the start of the time step.
  !<
  subroutine set_field(this, phi)
    ! -- dummy
    class(UltimateSchemeType), target :: this
    real(DP), intent(in), dimension(:), pointer :: phi
  end subroutine set_field

  !> @brief Calculate the explicit face values for this time step
  !!
  !! Verify that the time step satisfies the Courant constraint and then
  !! calculate and store a face value for every active connection.  Face values
  !! are stored by symmetric connection position so that both cells sharing a
  !! face see an identical value.
  !<
  subroutine prepare(this, phiold)
    ! -- dummy
    class(UltimateSchemeType) :: this
    real(DP), dimension(:), pointer, contiguous, intent(in) :: phiold !< field at start of time step
    ! -- local
    integer(I4B) :: n, m, ipos, isympos

    this%phiold => phiold
    call this%check_courant()
    call this%calc_derivatives()

    this%phiface(:) = DZERO
    do n = 1, this%dis%nodes
      if (this%ibound(n) == 0) cycle
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        if (m < n) cycle ! each face is visited once
        if (this%ibound(m) == 0) cycle
        isympos = this%dis%con%jas(ipos)
        this%phiface(isympos) = this%face_value(n, m, ipos)
      end do
    end do
  end subroutine prepare

  !> @brief Verify the Courant stability constraint
  !!
  !! The scheme is explicit, so it is stable only if no cell exchanges more than
  !! its own mobile pore volume during the time step.  The ADV Package requires
  !! ATS and ATS_PERCEL when this scheme is used, but that is not by itself
  !! enough: ATS does not apply the submitted stability constraint on the first
  !! time step of a stress period, and it will not reduce the time step below
  !! DTMIN.  The constraint is therefore verified here using the time step
  !! actually taken.
  !<
  subroutine check_courant(this)
    ! -- modules
    use TdisModule, only: delt, kstp, kper
    ! -- dummy
    class(UltimateSchemeType) :: this
    ! -- local
    integer(I4B) :: n, m, ipos, nrmax
    integer(I4B) :: iaxis, idir, nface
    real(DP) :: flownm, flowsumin, flowsumout, flowmax, flowaxis
    real(DP) :: capacity, courant, courantmax
    character(len=LINELENGTH) :: cellstr

    courantmax = DZERO
    nrmax = 0
    this%courant_axis(:, :) = DZERO

    do n = 1, this%dis%nodes
      if (this%ibound(n) == 0) cycle
      flowsumin = DZERO
      flowsumout = DZERO
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        if (this%dis%con%mask(ipos) == 0) cycle
        m = this%dis%con%ja(ipos)
        if (this%ibound(m) == 0) cycle
        flownm = this%fmi%gwfflowja(ipos)
        if (flownm < DZERO) then
          flowsumout = flowsumout - flownm
        else
          flowsumin = flowsumin + flownm
        end if
      end do
      flowmax = max(flowsumin, flowsumout)
      if (flowmax < DPREC) cycle
      capacity = this%dis%get_cell_volume(n, this%dis%top(n)) * &
                 this%fmi%gwfsat(n) * this%porosity(n) * this%retardation(n)
      if (capacity < DPREC) cycle
      courant = flowmax * delt / capacity
      !
      ! -- Resolve the flow through the cell onto each grid axis by averaging
      !    the flow through the two faces normal to it, signed positive in the
      !    direction of increasing cell index.  The transverse terms need a
      !    velocity along axes the face itself is not normal to.
      do iaxis = 1, 3
        flowaxis = DZERO
        nface = 0
        do idir = 2 * iaxis - 1, 2 * iaxis
          ipos = this%nbrpos(idir, n)
          if (ipos == 0) cycle
          if (this%ibound(this%nbr(idir, n)) == 0) cycle
          !
          ! -- Flow is positive into n, so flow arriving from the neighbor at
          !    the lower index travels along the axis and flow arriving from
          !    the neighbor at the higher index travels against it
          if (mod(idir, 2) == 1) then
            flowaxis = flowaxis + this%fmi%gwfflowja(ipos)
          else
            flowaxis = flowaxis - this%fmi%gwfflowja(ipos)
          end if
          nface = nface + 1
        end do
        if (nface > 0) then
          this%courant_axis(iaxis, n) = flowaxis * delt / (nface * capacity)
        end if
      end do
      if (courant > courantmax) then
        courantmax = courant
        nrmax = n
      end if
    end do

    if (courantmax > DONE + COURANT_SLACK) then
      call this%dis%noder_to_string(nrmax, cellstr)
      call store_error('The ULTIMATE advection scheme is explicit and requires &
                       &a Courant number of 1.0 or less.')
      write (errmsg, '(a,g0,a,a,a,i0,a,i0,a)') &
        'Courant number of ', courantmax, ' was calculated in cell ', &
        trim(adjustl(cellstr)), ' for time step ', kstp, &
        ' of stress period ', kper, '.'
      call store_error(errmsg)
      write (errmsg, '(a,g0,a)') &
        'Time step length must not exceed ', delt / courantmax, &
        '.  Check the DT0 and DTMIN settings in the ATS Package and the &
        &ATS_PERCEL setting in the ADV Package.'
      call store_error(errmsg, terminate=.TRUE.)
    end if
  end subroutine check_courant

  !> @brief Calculate the derivatives of the field at every cell center
  !!
  !! The transverse terms of the reconstruction need first and second
  !! derivatives along the two axes a face is not normal to, and the mixed
  !! derivative in the plane of the face.  These are properties of a cell rather
  !! than of a face, so they are calculated once per time step and reused by
  !! every face of the cell.
  !!
  !! Each derivative is centered and is set to zero unless the whole stencil it
  !! needs is present and active.  A missing neighbor therefore switches off the
  !! terms that would have used it and leaves the rest of the reconstruction
  !! intact.
  !<
  subroutine calc_derivatives(this)
    ! -- dummy
    class(UltimateSchemeType) :: this
    ! -- local
    integer(I4B) :: n, iaxis, jaxis, kaxis
    integer(I4B) :: idirm, idirp, nminus, nplus
    integer(I4B) :: node_mm, node_mp, node_pm, node_pp
    real(DP) :: phic, phiminus, phiplus
    real(DP) :: distminus, distplus, area

    this%grad_axis(:, :) = DZERO
    this%curv_axis(:, :) = DZERO
    this%cross_axis(:, :) = DZERO

    do n = 1, this%dis%nodes
      if (this%ibound(n) == 0) cycle
      phic = this%phiold(n)
      !
      ! -- First and second derivative along each axis
      do iaxis = 1, 3
        if (this%span(iaxis, n) < DPREC) cycle
        idirm = 2 * iaxis - 1
        idirp = 2 * iaxis
        nminus = this%nbr(idirm, n)
        nplus = this%nbr(idirp, n)
        if (this%ibound(nminus) == 0 .or. this%ibound(nplus) == 0) cycle
        phiminus = this%phiold(nminus)
        phiplus = this%phiold(nplus)
        distminus = this%distnbr(idirm, n)
        distplus = this%distnbr(idirp, n)
        this%grad_axis(iaxis, n) = (phiplus - phiminus) / this%span(iaxis, n)
        this%curv_axis(iaxis, n) = ((phiplus - phic) / distplus - &
                                    (phic - phiminus) / distminus) / &
                                   this%width(iaxis, n)
      end do
      !
      ! -- Mixed derivative in the plane normal to each axis, taken over the
      !    four cells diagonally adjacent to n in that plane
      do iaxis = 1, 3
        jaxis = mod(iaxis, 3) + 1
        kaxis = mod(iaxis + 1, 3) + 1
        area = this%span(jaxis, n) * this%span(kaxis, n)
        if (area < DPREC) cycle
        node_mm = this%diagonal_node(n, 2 * jaxis - 1, 2 * kaxis - 1)
        node_mp = this%diagonal_node(n, 2 * jaxis - 1, 2 * kaxis)
        node_pm = this%diagonal_node(n, 2 * jaxis, 2 * kaxis - 1)
        node_pp = this%diagonal_node(n, 2 * jaxis, 2 * kaxis)
        if (node_mm == 0 .or. node_mp == 0) cycle
        if (node_pm == 0 .or. node_pp == 0) cycle
        this%cross_axis(iaxis, n) = (this%phiold(node_pp) - &
                                     this%phiold(node_pm) - &
                                     this%phiold(node_mp) + &
                                     this%phiold(node_mm)) / area
      end do
    end do
  end subroutine calc_derivatives

  !> @brief Find the cell diagonally adjacent to n in two axis directions
  !!
  !! Returns zero if either step leaves the grid or lands on an inactive cell.
  !<
  function diagonal_node(this, n, idir1, idir2) result(ndiag)
    ! -- return
    integer(I4B) :: ndiag !< diagonally adjacent cell, or zero
    ! -- dummy
    class(UltimateSchemeType) :: this
    integer(I4B), intent(in) :: n !< cell to step from
    integer(I4B), intent(in) :: idir1 !< first axis direction
    integer(I4B), intent(in) :: idir2 !< second axis direction
    ! -- local
    integer(I4B) :: nmid

    ndiag = 0
    nmid = this%nbr(idir1, n)
    if (nmid == 0) return
    if (this%ibound(nmid) == 0) return
    ndiag = this%nbr(idir2, nmid)
    if (ndiag == 0) return
    if (this%ibound(ndiag) == 0) ndiag = 0
  end function diagonal_node

  !> @brief Calculate the explicit face value for one connection
  !!
  !! Build the QUICKEST reconstruction of the face concentration and clamp it
  !! with the universal limiter.  Falls back to first-order upwinding where the
  !! stencil cannot be formed.  Node n must be smaller than node m so that cl1
  !! and cl2 refer to n and m respectively.
  !<
  function face_value(this, n, m, iposnm) result(phi_face)
    ! -- modules
    use TdisModule, only: delt
    ! -- return
    real(DP) :: phi_face !< explicit concentration at the face
    ! -- dummy
    class(UltimateSchemeType) :: this
    integer(I4B), intent(in) :: n !< lower numbered cell
    integer(I4B), intent(in) :: m !< higher numbered cell
    integer(I4B), intent(in) :: iposnm !< position of the n-m connection
    ! -- local
    integer(I4B) :: iup, idn, i2up, isympos, iposup
    integer(I4B) :: iaxis, jaxis, kaxis
    real(DP) :: qnm, cl_up, cl_dn
    real(DP) :: phi_up, phi_dn, phi_2up
    real(DP) :: dist_updn, dist_2upup, width_up
    real(DP) :: weight_up, weight_dn, courant, capacity
    real(DP) :: gradient_updn, gradient_2upup, curvature
    real(DP) :: displacement_normal, gradient_trans, curvature_cross
    real(DP), dimension(3) :: displacement

    isympos = this%dis%con%jas(iposnm)
    qnm = this%fmi%gwfflowja(iposnm)
    !
    ! -- Find the upstream cell.  A positive flow is into n, so m is upstream.
    if (qnm > DZERO) then
      iup = m
      idn = n
      cl_up = this%dis%con%cl2(isympos)
      cl_dn = this%dis%con%cl1(isympos)
      iposup = this%dis%con%isym(iposnm)
    else
      iup = n
      idn = m
      cl_up = this%dis%con%cl1(isympos)
      cl_dn = this%dis%con%cl2(isympos)
      iposup = iposnm
    end if
    phi_up = this%phiold(iup)
    phi_dn = this%phiold(idn)
    !
    ! -- Default to the first-order upwind value
    phi_face = phi_up
    !
    ! -- Nothing to interpolate if there is no flow across the face
    if (abs(qnm) < DPREC) return
    !
    ! -- Revert to upwinding if there is no second upstream cell
    i2up = this%nodefar(iposup)
    if (i2up == 0) return
    if (this%ibound(i2up) == 0) return
    phi_2up = this%phiold(i2up)
    !
    ! -- Courant number for this face, based on the storage capacity of the
    !    upstream cell.  The capacity includes the retardation factor because
    !    the reconstruction traces the front backward over the time step and the
    !    front advances at the retarded velocity, not the water velocity.  This
    !    is the same measure of advective travel distance that the ADV Package
    !    submits to the ATS Package for this scheme.
    capacity = this%dis%get_cell_volume(iup, this%dis%top(iup)) * &
               this%fmi%gwfsat(iup) * this%porosity(iup) * this%retardation(iup)
    if (capacity < DPREC) return
    courant = abs(qnm) * delt / capacity
    !
    ! -- Distances.  For a structured grid the distance from the cell center to
    !    the face is half the cell width, so the width of the upstream cell
    !    along the connection is twice the connection length.
    dist_updn = cl_up + cl_dn
    dist_2upup = this%distfar(iposup)
    width_up = DTWO * cl_up
    weight_up = cl_dn / dist_updn
    weight_dn = DONE - weight_up
    !
    ! -- QUICKEST reconstruction: linear interpolation to the face, corrected
    !    for the distance traveled over the time step and for the curvature of
    !    the concentration profile on the upstream side.  On a uniform grid this
    !    reduces to
    !      (phi_up + phi_dn)/2 - courant/2*(phi_dn - phi_up)
    !        - (1 - courant^2)/6*(phi_dn - 2*phi_up + phi_2up)
    gradient_updn = (phi_dn - phi_up) / dist_updn
    gradient_2upup = (phi_up - phi_2up) / dist_2upup
    curvature = (gradient_updn - gradient_2upup) / width_up
    displacement_normal = courant * dist_updn
    phi_face = weight_up * phi_up + weight_dn * phi_dn &
               - DHALF * displacement_normal * gradient_updn &
               - (dist_updn * dist_updn - &
                  displacement_normal * displacement_normal) * curvature / DSIX
    !
    ! -- Transverse terms.  The face value is the average of the field over the
    !    volume swept through the face during the time step.  That volume is a
    !    slab of thickness displacement_normal, sheared sideways by the flow
    !    along the other two axes, so the average also picks up the variation of
    !    the field along and across those axes.  Distances along an axis are
    !    measured as the Courant number for that axis times the cell width,
    !    which is the distance the front travels over the time step.
    !
    !    Each transverse term is invariant under reversing the direction of any
    !    axis, so the normal direction is oriented downstream while the two
    !    transverse directions keep the orientation of the cell indexing.
    iaxis = this%faceaxis(isympos)
    displacement(:) = DZERO
    do jaxis = 1, 3
      if (jaxis == iaxis) cycle
      displacement(jaxis) = &
        weight_up * this%courant_axis(jaxis, iup) * this%width(jaxis, iup) + &
        weight_dn * this%courant_axis(jaxis, idn) * this%width(jaxis, idn)
      !
      ! -- Transverse advection uses the gradient interpolated to the face, and
      !    transverse curvature the second derivative in the upstream cell.  The
      !    mixed derivative is the change in the transverse gradient across the
      !    face, oriented downstream to match the normal displacement.
      gradient_trans = weight_up * this%grad_axis(jaxis, iup) + &
                       weight_dn * this%grad_axis(jaxis, idn)
      curvature_cross = (this%grad_axis(jaxis, idn) - &
                         this%grad_axis(jaxis, iup)) / dist_updn
      phi_face = phi_face &
                 - DHALF * displacement(jaxis) * gradient_trans &
                 + displacement(jaxis) * displacement(jaxis) * &
                 this%curv_axis(jaxis, iup) / DSIX &
                 + displacement_normal * displacement(jaxis) * &
                 curvature_cross / DTHREE
    end do
    !
    ! -- Twist in the plane of the face, which is nonzero only where flow has a
    !    component along both of the transverse axes
    jaxis = mod(iaxis, 3) + 1
    kaxis = mod(iaxis + 1, 3) + 1
    phi_face = phi_face + displacement(jaxis) * displacement(kaxis) * &
               this%cross_axis(iaxis, iup) / DTHREE
    !
    ! -- Enforce the TVD condition along the flow direction
    phi_face = universal_limiter(phi_face, phi_2up, phi_up, phi_dn, courant)
  end function face_value

  !> @brief Apply Leonard's universal limiter to a face value
  !!
  !! Clamp the high-order face value to the range that keeps the explicit update
  !! monotonic, expressed in terms of the three cells along the flow direction
  !! and the Courant number.  Where the upstream profile has a local extremum
  !! there is no monotonic range to preserve and the face value reverts to
  !! first-order upwinding.
  !<
  function universal_limiter(phi_face, phi_2up, phi_up, phi_dn, courant) &
    result(phi_limited)
    ! -- return
    real(DP) :: phi_limited !< limited face value
    ! -- dummy
    real(DP), intent(in) :: phi_face !< unlimited face value
    real(DP), intent(in) :: phi_2up !< second upstream cell
    real(DP), intent(in) :: phi_up !< upstream cell
    real(DP), intent(in) :: phi_dn !< downstream cell
    real(DP), intent(in) :: courant !< Courant number for the face
    ! -- local
    real(DP) :: bound

    if (phi_2up >= phi_up .and. phi_up >= phi_dn) then
      !
      ! -- Profile decreases monotonically along the flow direction
      bound = phi_dn
      if (courant > DPREC) then
        bound = max(phi_dn, phi_2up + (phi_up - phi_2up) / courant)
      end if
      phi_limited = max(min(phi_face, phi_up), bound)
    else if (phi_2up <= phi_up .and. phi_up <= phi_dn) then
      !
      ! -- Profile increases monotonically along the flow direction
      bound = phi_dn
      if (courant > DPREC) then
        bound = min(phi_dn, phi_2up + (phi_up - phi_2up) / courant)
      end if
      phi_limited = min(max(phi_face, phi_up), bound)
    else
      !
      ! -- Local extremum; revert to first-order upwinding
      phi_limited = phi_up
    end if
  end function universal_limiter

  !> @brief Return the interpolation coefficients for a face
  !!
  !! The face value is explicit, so it carries no dependence on the cell values
  !! being solved for and is returned entirely through the right-hand side term.
  !<
  function compute(this, n, m, iposnm) result(phi_face)
    ! -- return
    type(CoefficientsType) :: phi_face !< coefficients for the n-m face
    ! -- dummy
    class(UltimateSchemeType), target :: this
    integer(I4B), intent(in) :: n !< index of the first cell
    integer(I4B), intent(in) :: m !< index of the second cell
    integer(I4B), intent(in) :: iposnm !< position of the n-m connection

    phi_face%c_n = DZERO
    phi_face%c_m = DZERO
    phi_face%rhs = -this%phiface(this%dis%con%jas(iposnm))
  end function compute

end module UltimateSchemeModule
