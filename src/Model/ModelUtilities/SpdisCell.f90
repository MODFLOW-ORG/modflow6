module SpdisCellModule
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: DZERO, DONE, C3D_VERTICAL, &
                             DPI, DTWOPI, DNODATA, DEM10
  use BaseDisModule, only: DisBaseType
  use DisUtilsModule, only: number_connected_faces
  use SimModule, only: store_error
  implicit none
  private

  public :: SpdisCellType

  ! TODO_MJR: how accurately can we expect angles to be defined,
  ! considering exchanges with conversion from degrees to radians,
  ! vertex information parsed from file, etc.
  real(DP), parameter :: TINY_ANGLE = DEM10

  !> @brief Container for spdis calculation in cell
  !!
  !! face index: 1 = top face, 2 = bottom face, and
  !! {3, .., nr_faces} vertical faces in clockwise order
  !! starting from lower left corner.
  !<
  type SpdisCellType
    integer(I4B) :: nr_faces
    real(DP), dimension(:), allocatable :: anglex !< face angle w.r.t. x-axis
    real(DP), dimension(:), allocatable :: dist !< distance
    real(DP), dimension(:, :), allocatable :: norm !< face normal
    integer(I4B), dimension(:), allocatable :: matched !< =1 when there is an internal or external flow associated, 0 otherwise

    class(DisBaseType), pointer :: dis => null() !< discretization
    real(DP), dimension(:), pointer :: flowja => null() !< internal flows
    ! exchanges, TODO_MJR: this should probably be grouped in a type
    integer(I4B), dimension(:), pointer, contiguous :: ihc_edge => null() !< face type (horizontal or vertical)
    real(DP), dimension(:, :), pointer, contiguous :: props_edge => null() !< face properties (Q, area, nx, ny, distance)
    integer(I4B), dimension(:), pointer, contiguous :: iedge_ptr => null() !< csr pointer into edge index array
    integer(I4B), dimension(:), pointer, contiguous :: edge_idxs => null() !< sorted edge indexes for faster lookup
  contains
    procedure :: create
    procedure :: load_cell_boundaries
    procedure :: destroy
  end type SpdisCellType

contains

  !> @brief Create the worker data structure
  !<
  subroutine create(this, dis, flowja, ihc, props, iedge_ptr, edge_idxs)
    class(SpdisCellType) :: this
    class(DisBaseType), pointer :: dis !< discretization
    real(DP), dimension(:), target :: flowja !< the internal flows
    integer(I4B), dimension(:), pointer, contiguous :: ihc !< ihc flag for faces
    real(DP), dimension(:, :), pointer, contiguous :: props !< properties for faces
    integer(I4B), dimension(:), pointer, contiguous :: iedge_ptr !< csr pointer into edge index array
    integer(I4B), dimension(:), pointer, contiguous :: edge_idxs !< sorted edge indexes for faster lookup
    ! local
    integer(I4B) :: max_nr_faces

    this%dis => dis
    this%flowja => flowja
    this%ihc_edge => ihc
    this%props_edge => props
    this%iedge_ptr => iedge_ptr
    this%edge_idxs => edge_idxs

    ! max. nr. of faces is equal to max. nr. of vertices, plus top and bottom face
    max_nr_faces = dis%get_max_npolyverts(closed=.false.) + 2

    allocate (this%anglex(max_nr_faces))
    allocate (this%dist(max_nr_faces))
    allocate (this%norm(3, max_nr_faces))
    allocate (this%matched(max_nr_faces))

  end subroutine create

  !> @brief Load boundary faces for cell
  !<
  subroutine load_cell_boundaries(this, n, bnd_faces)
    class(SpdisCellType) :: this
    integer(I4B), intent(in) :: n !< reduced node number
    real(DP), dimension(:, :), allocatable, intent(inout) :: bnd_faces !< the boundary faces for the cell: dist, nx, ny, nz
    ! local
    integer(I4B) :: ivert, ipos, isym, m, iface, iedge, ibnd
    integer(I4B) :: npoly, nr_boundaries, icand
    logical(LGP) :: have_match
    real(DP), dimension(:, :), allocatable :: pverts
    real(DP), dimension(2) :: v
    real(DP) :: length, nx, ny, nz, alpha, dangle
    real(DP) :: fcx, fcy

    ! reset
    this%matched(:) = 0

    ! determine nr of faces
    npoly = this%dis%get_npolyverts(n, closed=.false.)
    this%nr_faces = npoly + 2

    ! horizontal faces
    this%norm(:, 1) = (/DZERO, DZERO, DONE/)
    this%anglex(1) = DNODATA
    this%dist(1) = (this%dis%top(n) - this%dis%bot(n)) / 2.0_DP
    this%norm(:, 2) = (/DZERO, DZERO, -DONE/)
    this%anglex(2) = DNODATA
    this%dist(2) = (this%dis%top(n) - this%dis%bot(n)) / 2.0_DP

    ! vertical faces
    call this%dis%get_polyverts(n, pverts)
    do ivert = 1, npoly
      iface = 2 + ivert ! start after horizontal faces

      v = pverts(:, modulo(ivert, npoly) + 1) - pverts(:, ivert)
      length = sqrt(v(1) * v(1) + v(2) * v(2))

      this%norm(:, iface) = (/-v(2) / length, v(1) / length, DZERO/)
      this%anglex(iface) = acos(this%norm(1, iface)) ! between 0 and pi
      if (this%norm(2, iface) < DZERO) then ! if y is negative, take complement
        this%anglex(iface) = DTWOPI - this%anglex(iface)
      end if

      fcx = 0.5_DP * (pverts(1, modulo(ivert, npoly) + 1) + pverts(1, ivert))
      fcy = 0.5_DP * (pverts(2, modulo(ivert, npoly) + 1) + pverts(2, ivert))
      this%dist(iface) = sqrt((fcx - this%dis%xc(n)) * (fcx - this%dis%xc(n)) + &
                              (fcy - this%dis%yc(n)) * (fcy - this%dis%yc(n)))
    end do

    ! match internal connections
    have_match = .true.
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      have_match = .false.
      m = this%dis%con%ja(ipos)
      isym = this%dis%con%jas(ipos)
      if (this%dis%con%ihc(isym) == C3D_VERTICAL) then
        if (m > n) then ! bottom
          this%matched(2) = this%matched(2) + 1
          have_match = .true.
        else
          this%matched(1) = this%matched(1) + 1
          have_match = .true.
        end if
      else
        ! match angle for horizontal connections
        dangle = huge(dangle)
        icand = 0
        do iface = 3, this%nr_faces
          if (n > m) then ! flip the normal out of its symmetric storage
            alpha = modulo(DPI + this%dis%con%anglex(isym), DTWOPI)
          else
            alpha = this%dis%con%anglex(isym)
          end if
          ! find the best match
          if (abs(alpha - this%anglex(iface)) < dangle) then
            dangle = abs(alpha - this%anglex(iface))
            icand = iface
          end if
        end do

        if (icand > 0 .and. dangle < TINY_ANGLE) then
          this%matched(icand) = this%matched(icand) + 1
          have_match = .true.
        end if

      end if

      if (.not. have_match) then
        call store_error( &
          "Invalid connection encountered in specific discharge calculation", &
          terminate=.true.)
      end if
    end do

    ! match exchange flows
    have_match = .true.
    if (size(this%edge_idxs) > 0) then
      do ipos = this%iedge_ptr(n), this%iedge_ptr(n + 1) - 1
        have_match = .false.
        iedge = this%edge_idxs(ipos)
        if (this%ihc_edge(iedge) == C3D_VERTICAL) then
          nz = this%props_edge(5, iedge)
          if (nz > DZERO) then ! top face
            this%matched(1) = this%matched(1) + 1
            have_match = .true.
          else
            this%matched(2) = this%matched(2) + 1
            have_match = .true.
          end if
        else
          ! get angle
          nx = this%props_edge(3, iedge)
          ny = this%props_edge(4, iedge)
          alpha = acos(nx) ! between 0 and pi
          if (ny < DZERO .and. alpha > DZERO) then ! if y is negative, take complement
            alpha = DTWOPI - alpha
          end if
          ! try match to face
          do iface = 3, this%nr_faces
            if (abs(alpha - this%anglex(iface)) < TINY_ANGLE) then
              ! TODO_MJR: we should match uniquely, need more geometry unfortunately...
              this%matched(iface) = this%matched(iface) + 1
              have_match = .true.
              exit
            end if
          end do
        end if

        if (.not. have_match) then
          call store_error( &
            "Invalid exchange flow in specific discharge calculation", &
            terminate=.true.)
        end if

      end do

    end if

    ! what's left are boundary faces
    nr_boundaries = this%nr_faces
    do iface = 1, this%nr_faces
      if (this%matched(iface) > 0) then
        nr_boundaries = nr_boundaries - 1
      end if
    end do

    allocate (bnd_faces(4, nr_boundaries))
    ibnd = 0
    do iface = 1, this%nr_faces
      if (this%matched(iface) == 0) then
        ibnd = ibnd + 1
        bnd_faces(:, ibnd) = &
          [this%dist(iface), this%norm(1, iface), this%norm(2, iface), this%norm(3, iface)]
      end if
    end do

  end subroutine load_cell_boundaries

  subroutine destroy(this)
    class(SpdisCellType) :: this

    deallocate (this%anglex)
    deallocate (this%dist)
    deallocate (this%norm)
    deallocate (this%matched)

  end subroutine destroy

end module SpdisCellModule
