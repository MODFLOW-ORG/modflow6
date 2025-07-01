module MethodCellTernaryModule

  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use MethodModule, only: MethodType
  use MethodCellModule, only: MethodCellType
  use MethodSubcellPoolModule
  use CellPolyModule
  use CellDefnModule
  use SubcellTriModule, only: SubcellTriType, create_subcell_tri
  use ParticleModule
  use GeomUtilModule, only: area, transform
  use ConstantsModule, only: DZERO, DONE, DTWO, DSIX
  use GeomUtilModule, only: point_in_polygon
  implicit none

  private
  public :: MethodCellTernaryType
  public :: create_method_cell_ternary

  type, extends(MethodCellType) :: MethodCellTernaryType
    private
    integer(I4B) :: nverts !< number of vertices
    real(DP), allocatable, dimension(:) :: xvert
    real(DP), allocatable, dimension(:) :: yvert !< cell vertex coordinates
    real(DP), allocatable, dimension(:) :: vne !< cell edge normal velocities
    real(DP), allocatable, dimension(:) :: vv0x
    real(DP), allocatable, dimension(:) :: vv0y
    real(DP), allocatable, dimension(:) :: vv1x
    real(DP), allocatable, dimension(:) :: vv1y !< cell vertex velocities
    real(DP) :: xctr
    real(DP) :: yctr !< cell center coordinates
    real(DP) :: vctrx
    real(DP) :: vctry !< cell center velocities
    real(DP) :: ztop
    real(DP) :: zbot !< cell top and bottom elevations
    real(DP) :: dz !< cell thickness
    real(DP) :: vztop
    real(DP) :: vzbot !< cell top and bottom velocities
    integer(I4B), allocatable, dimension(:) :: iprev !< array of shifted indices
    real(DP), allocatable, dimension(:) :: xvertnext
    real(DP), allocatable, dimension(:) :: yvertnext !< arrays of shifted cell vertex coordinates
    integer(I4B), public, pointer :: zeromethod
  contains
    procedure, public :: apply => apply_mct
    procedure, public :: deallocate => destroy_mct
    procedure, public :: load => load_mct
    procedure, public :: load_subcell
    procedure, public :: pass => pass_mct
    procedure :: vertvelo
    procedure :: calc_thru_hcsum
  end type MethodCellTernaryType

contains

  !> @brief Create a tracking method
  subroutine create_method_cell_ternary(method)
    ! dummy
    type(MethodCellTernaryType), pointer :: method
    ! local
    type(CellPolyType), pointer :: cell
    type(SubcellTriType), pointer :: subcell

    allocate (method)
    call create_cell_poly(cell)
    method%cell => cell
    method%name => method%cell%type
    method%delegates = .true.
    call create_subcell_tri(subcell)
    method%subcell => subcell
  end subroutine create_method_cell_ternary

  !> @brief Destroy the tracking method
  subroutine destroy_mct(this)
    class(MethodCellTernaryType), intent(inout) :: this
    deallocate (this%name)
  end subroutine destroy_mct

  !> @brief Load subcell into tracking method
  subroutine load_mct(this, particle, next_level, submethod)
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer(I4B), intent(in) :: next_level
    class(MethodType), pointer, intent(inout) :: submethod

    select type (subcell => this%subcell)
    type is (SubcellTriType)
      call this%load_subcell(particle, subcell)
    end select

    call method_subcell_tern%init( &
      fmi=this%fmi, &
      cell=this%cell, &
      subcell=this%subcell, &
      trackctl=this%trackctl, &
      tracktimes=this%tracktimes)
    submethod => method_subcell_tern
  end subroutine load_mct

  !> @brief Pass particle to next subcell if there is one, or to the cell face
  subroutine pass_mct(this, particle)
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    ! local
    integer(I4B) :: isc
    integer(I4B) :: exitFace
    integer(I4B) :: inface

    exitFace = particle%iboundary(3)
    isc = particle%idomain(3)

    select case (exitFace)
    case (0)
      ! Subcell interior (cell interior)
      inface = -1
    case (1)
      ! Subcell face 1 (cell face)
      inface = isc
      if (inface .eq. 0) inface = this%nverts
    case (2)
      ! Subcell face --> next subcell in "cycle" (cell interior)
      isc = isc + 1
      if (isc .gt. this%nverts) isc = 1
      particle%idomain(3) = isc
      particle%iboundary(3) = 3
      inface = 0
    case (3)
      ! Subcell face --> preceding subcell in "cycle" (cell interior)
      isc = isc - 1
      if (isc .lt. 1) isc = this%nverts
      particle%idomain(3) = isc
      particle%iboundary(3) = 2
      inface = 0
    case (4)
      ! Subcell bottom (cell bottom)
      inface = this%nverts + 2
    case (5)
      ! Subcell top (cell top)
      inface = this%nverts + 3
    end select

    if (inface .eq. -1) then
      particle%iboundary(2) = 0
    else if (inface .eq. 0) then
      particle%iboundary(2) = 0
    else
      particle%iboundary(2) = inface
    end if
  end subroutine pass_mct

  !> @brief Apply the ternary method to a polygonal cell
  subroutine apply_mct(this, particle, tmax)
    use ConstantsModule, only: DZERO, DONE, DHALF
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local
    integer(I4B) :: i
    real(DP) :: x, y, z, xO, yO
    real(DP), allocatable :: xs(:), ys(:)

    ! (Re)allocate type-bound arrays
    select type (cell => this%cell)
    type is (CellPolyType)
      ! Check termination/reporting conditions
      call this%check(particle, this%cell%defn, tmax)
      if (.not. particle%advancing) return

      ! Number of vertices
      this%nverts = cell%defn%npolyverts

      ! (Re)allocate type-bound arrays
      if (allocated(this%xvert)) then
        deallocate (this%xvert)
        deallocate (this%yvert)
        deallocate (this%vne)
        deallocate (this%vv0x)
        deallocate (this%vv0y)
        deallocate (this%vv1x)
        deallocate (this%vv1y)
        deallocate (this%iprev)
        deallocate (this%xvertnext)
        deallocate (this%yvertnext)
      end if

      allocate (this%xvert(this%nverts))
      allocate (this%yvert(this%nverts))
      allocate (this%vne(this%nverts))
      allocate (this%vv0x(this%nverts))
      allocate (this%vv0y(this%nverts))
      allocate (this%vv1x(this%nverts))
      allocate (this%vv1y(this%nverts))
      allocate (this%iprev(this%nverts))
      allocate (this%xvertnext(this%nverts))
      allocate (this%yvertnext(this%nverts))

      ! New origin is the corner of the cell's
      ! bounding box closest to the old origin
      allocate (xs(this%nverts))
      allocate (ys(this%nverts))
      xs = cell%defn%polyvert(1, :)
      ys = cell%defn%polyvert(2, :)
      xO = xs(minloc(abs(xs), dim=1))
      yO = ys(minloc(abs(ys), dim=1))
      deallocate (xs)
      deallocate (ys)

      ! Cell vertices
      do i = 1, this%nverts
        x = cell%defn%polyvert(1, i)
        y = cell%defn%polyvert(2, i)
        call transform(x, y, DZERO, x, y, z, xO, yO)
        this%xvert(i) = x
        this%yvert(i) = y
      end do

      ! Top, bottom, and thickness
      this%ztop = cell%defn%top
      this%zbot = cell%defn%bot
      this%dz = this%ztop - this%zbot

      ! Shifted arrays
      do i = 1, this%nverts
        this%iprev(i) = i
      end do
      this%iprev = cshift(this%iprev, -1)
      this%xvertnext = cshift(this%xvert, 1)
      this%yvertnext = cshift(this%yvert, 1)

      ! Calculate vertex velocities.
      call this%vertvelo()

      ! Transform particle coordinates
      call particle%transform(xO, yO)

      ! Track the particle across the cell.
      call this%track(particle, 2, tmax)

      ! Transform particle coordinates back
      call particle%transform(xO, yO, invert=.true.)
      call particle%reset_transform()

    end select

  end subroutine apply_mct

  !> @brief Loads a triangular subcell from the polygonal cell
  subroutine load_subcell(this, particle, subcell)
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(SubcellTriType), intent(inout) :: subcell
    ! local
    integer(I4B) :: ic
    integer(I4B) :: isc
    integer(I4B) :: iv0
    real(DP) :: x0
    real(DP) :: y0
    real(DP) :: x1
    real(DP) :: y1
    real(DP) :: x2
    real(DP) :: y2
    real(DP) :: x1rel
    real(DP) :: y1rel
    real(DP) :: x2rel
    real(DP) :: y2rel
    real(DP) :: xi
    real(DP) :: yi
    real(DP) :: di2
    real(DP) :: d02
    real(DP) :: d12
    real(DP) :: di1
    real(DP) :: d01
    real(DP) :: alphai
    real(DP) :: betai

    select type (cell => this%cell)
    type is (CellPolyType)
      ic = cell%defn%icell
      subcell%icell = ic
      isc = particle%idomain(3)
      if (isc .le. 0) then
        xi = particle%x
        yi = particle%y
        do iv0 = 1, this%nverts
          x0 = this%xvert(iv0)
          y0 = this%yvert(iv0)
          x1 = this%xvertnext(iv0)
          y1 = this%yvertnext(iv0)
          x2 = this%xctr
          y2 = this%yctr
          x1rel = x1 - x0
          y1rel = y1 - y0
          x2rel = x2 - x0
          y2rel = y2 - y0
          di2 = xi * y2rel - yi * x2rel
          d02 = x0 * y2rel - y0 * x2rel
          d12 = x1rel * y2rel - y1rel * x2rel
          di1 = xi * y1rel - yi * x1rel
          d01 = x0 * y1rel - y0 * x1rel
          alphai = (di2 - d02) / d12
          betai = -(di1 - d01) / d12
          ! assumes particle is in the cell, so no check needed for beta
          if ((alphai .ge. DZERO) .and. &
              (alphai + betai .le. DONE)) then
            isc = iv0
            exit
          end if
        end do
        if (isc .le. 0) then
          print *, "error -- initial triangle not found in cell ", ic, &
            " for particle at ", particle%x, particle%y, particle%z

          call pstop(1)
        else
          particle%idomain(3) = isc
        end if
      end if
      subcell%isubcell = isc

      ! Set coordinates and velocities at vertices of triangular subcell
      iv0 = isc
      subcell%x0 = this%xvert(iv0)
      subcell%y0 = this%yvert(iv0)
      subcell%x1 = this%xvertnext(iv0)
      subcell%y1 = this%yvertnext(iv0)
      subcell%x2 = this%xctr
      subcell%y2 = this%yctr
      subcell%v0x = this%vv0x(iv0)
      subcell%v0y = this%vv0y(iv0)
      subcell%v1x = this%vv1x(iv0) ! the indices here actually refer to subcells, not vertices
      subcell%v1y = this%vv1y(iv0)
      subcell%v2x = this%vctrx
      subcell%v2y = this%vctry
      subcell%ztop = this%ztop
      subcell%zbot = this%zbot
      subcell%dz = this%dz
      subcell%vzbot = this%vzbot
      subcell%vztop = this%vztop
    end select
  end subroutine load_subcell

  !> @brief Calculate vertex velocities
  subroutine vertvelo(this)
    use ConstantsModule, only: DZERO, DONE, DHALF
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    ! local
    real(DP) :: term
    integer(I4B) :: i
    real(DP) :: perturb
    real(DP), allocatable, dimension(:) :: xvals
    real(DP), allocatable, dimension(:) :: yvals
    real(DP) :: sixa
    real(DP) :: vm0i0
    real(DP) :: vm0ival
    real(DP) :: hcsum0
    real(DP) :: hcsum
    real(DP) :: jac
    real(DP), allocatable, dimension(:) :: wk1
    real(DP), allocatable, dimension(:) :: wk2
    real(DP), allocatable, dimension(:) :: unextxnext
    real(DP), allocatable, dimension(:) :: unextynext
    real(DP), allocatable, dimension(:) :: le
    real(DP), allocatable, dimension(:) :: unextx
    real(DP), allocatable, dimension(:) :: unexty
    real(DP) :: areacell
    real(DP), allocatable, dimension(:) :: areasub
    real(DP) :: divcell
    real(DP), allocatable, dimension(:) :: li
    real(DP), allocatable, dimension(:) :: unintx
    real(DP), allocatable, dimension(:) :: uninty
    real(DP), allocatable, dimension(:) :: xmid
    real(DP), allocatable, dimension(:) :: ymid
    real(DP), allocatable, dimension(:) :: lm
    real(DP), allocatable, dimension(:) :: umx
    real(DP), allocatable, dimension(:) :: umy
    real(DP), allocatable, dimension(:) :: kappax
    real(DP), allocatable, dimension(:) :: kappay
    real(DP), allocatable, dimension(:) :: vm0x
    real(DP), allocatable, dimension(:) :: vm0y
    real(DP), allocatable, dimension(:) :: vm1x
    real(DP), allocatable, dimension(:) :: vm1y
    ! real(DP), allocatable, dimension(:, :) :: poly
    integer(I4B) :: nvert

    select type (cell => this%cell)
    type is (CellPolyType)

      ! Allocate local arrays
      allocate (le(this%nverts)) ! lengths of exterior (cell) edges
      allocate (unextx(this%nverts)) ! x components of unit normals to exterior edges
      allocate (unexty(this%nverts)) ! y components of unit normals to exterior edges
      allocate (areasub(this%nverts)) ! subcell areas
      allocate (li(this%nverts)) ! lengths of interior edges ("spokes")
      allocate (unintx(this%nverts)) ! x components of unit normals to interior edges
      allocate (uninty(this%nverts)) ! y components of unit normals to interior edges
      allocate (xmid(this%nverts)) ! x coordinates of midpoints
      allocate (ymid(this%nverts)) ! y coordinates of midpoints
      allocate (lm(this%nverts)) ! lengths of midpoint connectors
      allocate (umx(this%nverts)) ! x components of midpoint-connector (cw) unit vectors
      allocate (umy(this%nverts)) ! y components of midpoint-connector (cw) unit vectors
      allocate (kappax(this%nverts)) ! x components of kappa vectors
      allocate (kappay(this%nverts)) ! y components of kappa vectors
      allocate (vm0x(this%nverts)) ! x component of vm0
      allocate (vm0y(this%nverts)) ! y component of vm0
      allocate (vm1x(this%nverts)) ! x component of vm1
      allocate (vm1y(this%nverts)) ! y component of vm1
      allocate (unextxnext(this%nverts)) ! vector of "next" interior unit-normal x coordinates defined for convenience
      allocate (unextynext(this%nverts)) ! vector of "next" interior unit-normal y coordinates defined for convenience
      allocate (wk1(this%nverts))
      allocate (wk2(this%nverts))
      allocate (xvals(3))
      allocate (yvals(3))

      ! Exterior edge unit normals (outward) and lengths
      wk1 = this%xvertnext - this%xvert
      wk2 = this%yvertnext - this%yvert
      le = dsqrt(wk1 * wk1 + wk2 * wk2)
      unextx = -wk2 / le
      unexty = wk1 / le

      ! Cell centroid (in general, this is NOT the average of the vertex coordinates)
      areacell = area(this%xvert, this%yvert)
      sixa = areacell * DSIX
      wk1 = -(this%xvert * this%yvertnext - this%xvertnext * this%yvert)
      nvert = size(this%xvert)
      this%xctr = sum((this%xvert + this%xvertnext) * wk1) / sixa
      this%yctr = sum((this%yvert + this%yvertnext) * wk1) / sixa

      ! TODO: can we use some of the terms of the centroid calculation
      ! to do a cheap point in polygon check?

      ! Subcell areas
      do i = 1, this%nverts
        xvals(1) = this%xvert(i)
        xvals(2) = this%xvertnext(i)
        xvals(3) = this%xctr
        yvals(1) = this%yvert(i)
        yvals(2) = this%yvertnext(i)
        yvals(3) = this%yctr
        areasub(i) = area(xvals, yvals)
      end do

      ! Cell-edge normal velocities (outward)
      term = DONE / (cell%defn%porosity * cell%defn%retfactor * this%dz)
      do i = 1, this%nverts
        this%vne(i) = -cell%defn%faceflow(i) * term / le(i)
      end do

      ! Cell divergence (2D)
      divcell = sum(le * this%vne) / areacell

      ! Interior edge (cw) unit normals and lengths
      wk1 = this%xvert - this%xctr
      wk2 = this%yvert - this%yctr
      li = dsqrt(wk1 * wk1 + wk2 * wk2)
      unintx = wk2 / li
      uninty = -wk1 / li
      ! Shifted arrays for convenience
      unextxnext = cshift(unintx, 1)
      unextynext = cshift(uninty, 1)

      ! Midpoints of interior edges
      xmid = 5.d-1 * (this%xvert + this%xctr)
      ymid = 5.d-1 * (this%yvert + this%yctr)

      ! Unit midpoint-connector (cw) vectors and lengths
      wk1 = cshift(xmid, 1) - xmid
      wk2 = cshift(ymid, 1) - ymid
      lm = dsqrt(wk1 * wk1 + wk2 * wk2)
      umx = wk1 / lm
      umy = wk2 / lm

      ! Kappa vectors (K tensor times unit midpoint-connector vectors) do not
      ! account for anisotropy, which is consistent with the way internal face
      ! flow calculations are done in MP7. The isotropic value of K does not
      ! matter in this case because it cancels out of the calculations, so
      ! K = 1 is assumed for simplicity.
      kappax = umx
      kappay = umy

      ! Use linearity to find vm0i[0] such that curl of the head gradient
      ! is zero
      perturb = 1.d-2
      ! Calculations at base value
      vm0i0 = 0.d0
      call this%calc_thru_hcsum(vm0i0, divcell, le, li, lm, areasub, areacell, &
                                unintx, uninty, unextx, unexty, &
                                unextxnext, unextynext, &
                                kappax, kappay, vm0x, vm0y, vm1x, vm1y, hcsum0)
      ! Calculations at perturbed value
      vm0ival = vm0i0 + perturb
      call this%calc_thru_hcsum(vm0ival, divcell, le, li, lm, areasub, areacell, &
                                unintx, uninty, unextx, unexty, &
                                unextxnext, unextynext, &
                                kappax, kappay, vm0x, vm0y, vm1x, vm1y, hcsum)
      ! Calculations at root value
      jac = (hcsum - hcsum0) / perturb
      vm0ival = vm0i0 - hcsum0 / jac
      call this%calc_thru_hcsum(vm0ival, divcell, le, li, lm, areasub, areacell, &
                                unintx, uninty, unextx, unexty, &
                                unextxnext, unextynext, &
                                kappax, kappay, vm0x, vm0y, vm1x, vm1y, hcsum)

      ! Project linearly to get corner (vertex) velocities. Note that velocity
      ! vv1 is at the next vertex cw from vv0, so vv0(i) and vv1(i) are the
      ! two vertex velocities used by triangular subcell i.
      this%vv0x = 2.d0 * vm0x - this%vctrx
      this%vv0y = 2.d0 * vm0y - this%vctry
      this%vv1x = 2.d0 * vm1x - this%vctrx
      this%vv1y = 2.d0 * vm1y - this%vctry

      ! Set top and bottom velocities
      term = DONE / (cell%defn%retfactor * cell%defn%porosity * areacell)
      this%vzbot = cell%defn%faceflow(this%nverts + 2) * term
      this%vztop = -cell%defn%faceflow(this%nverts + 3) * term

      ! Deallocate local arrays
      deallocate (le)
      deallocate (unextx)
      deallocate (unexty)
      deallocate (areasub)
      deallocate (li)
      deallocate (unintx)
      deallocate (uninty)
      deallocate (xmid)
      deallocate (ymid)
      deallocate (lm)
      deallocate (umx)
      deallocate (umy)
      deallocate (kappax)
      deallocate (kappay)
      deallocate (vm0x)
      deallocate (vm0y)
      deallocate (vm1x)
      deallocate (vm1y)
      deallocate (unextxnext)
      deallocate (unextynext)
      deallocate (wk1)
      deallocate (wk2)
      deallocate (xvals)
      deallocate (yvals)

    end select
  end subroutine vertvelo

  subroutine calc_thru_hcsum(this, vm0ival, divcell, &
                             le, li, lm, areasub, areacell, &
                             unintx, uninty, unextx, unexty, &
                             unintxnext, unintynext, &
                             kappax, kappay, vm0x, vm0y, vm1x, vm1y, hcsum)
    ! dummy
    class(MethodCellTernaryType), intent(inout) :: this
    real(DP) :: vm0ival
    real(DP) :: divcell
    real(DP) :: hcsum
    real(DP), dimension(:) :: le
    real(DP), dimension(:) :: li
    real(DP), dimension(:) :: lm
    real(DP), dimension(:) :: areasub
    real(DP) :: areacell
    real(DP), dimension(:) :: unintx
    real(DP), dimension(:) :: uninty
    real(DP), dimension(:) :: unextx
    real(DP), dimension(:) :: unexty
    real(DP), dimension(:) :: unintxnext
    real(DP), dimension(:) :: unintynext
    real(DP), dimension(:) :: kappax
    real(DP), dimension(:) :: kappay
    real(DP), dimension(:) :: vm0x
    real(DP), dimension(:) :: vm0y
    real(DP), dimension(:) :: vm1x
    real(DP), dimension(:) :: vm1y
    ! local
    real(DP), allocatable, dimension(:) :: vm0i
    real(DP), allocatable, dimension(:) :: vm0e
    real(DP), allocatable, dimension(:) :: vm1i
    real(DP), allocatable, dimension(:) :: vm1e
    real(DP), allocatable, dimension(:) :: uprod
    real(DP), allocatable, dimension(:) :: det
    real(DP), allocatable, dimension(:) :: wt
    real(DP), allocatable, dimension(:) :: bi0x
    real(DP), allocatable, dimension(:) :: be0x
    real(DP), allocatable, dimension(:) :: bi0y
    real(DP), allocatable, dimension(:) :: be0y
    real(DP), allocatable, dimension(:) :: bi1x
    real(DP), allocatable, dimension(:) :: be1x
    real(DP), allocatable, dimension(:) :: bi1y
    real(DP), allocatable, dimension(:) :: be1y
    real(DP), allocatable, dimension(:) :: be01x
    real(DP), allocatable, dimension(:) :: be01y
    real(DP) :: emxx
    real(DP) :: emxy
    real(DP) :: emyx
    real(DP) :: emyy
    real(DP) :: rx
    real(DP) :: ry
    real(DP) :: emdet
    integer(I4B) :: i
    integer(I4B) :: ip

    ! Allocate local arrays
    allocate (vm0i(this%nverts))
    allocate (vm0e(this%nverts))
    allocate (vm1i(this%nverts))
    allocate (vm1e(this%nverts))
    allocate (uprod(this%nverts))
    allocate (det(this%nverts))
    allocate (wt(this%nverts))
    allocate (bi0x(this%nverts))
    allocate (be0x(this%nverts))
    allocate (bi0y(this%nverts))
    allocate (be0y(this%nverts))
    allocate (bi1x(this%nverts))
    allocate (be1x(this%nverts))
    allocate (bi1y(this%nverts))
    allocate (be1y(this%nverts))
    allocate (be01x(this%nverts))
    allocate (be01y(this%nverts))

    ! Set vm0i(1)
    vm0i(1) = vm0ival

    ! Get remaining vm0i's sequentially using divergence conditions
    do i = 2, this%nverts
      ip = this%iprev(i)
      vm0i(i) = (li(ip) * vm0i(ip) - le(ip) * this%vne(ip) &
                 + areasub(ip) * divcell) / li(i)
    end do

    ! Get vm1i's from vm0i's using continuity conditions
    vm1i = cshift(vm0i, 1)

    ! Get centroid velocity by setting up and solving 2x2 linear system
    uprod = unintx * unextx + uninty * unexty
    det = DONE - uprod * uprod
    bi0x = (unintx - unextx * uprod) / det
    be0x = (unextx - unintx * uprod) / det
    bi0y = (uninty - unexty * uprod) / det
    be0y = (unexty - uninty * uprod) / det
    uprod = unintxnext * unextx + unintynext * unexty
    det = DONE - uprod * uprod
    bi1x = (unintxnext - unextx * uprod) / det
    be1x = (unextx - unintxnext * uprod) / det
    bi1y = (unintynext - unexty * uprod) / det
    be1y = (unexty - unintynext * uprod) / det
    be01x = 5.d-1 * (be0x + be1x)
    be01y = 5.d-1 * (be0y + be1y)
    wt = areasub / areacell
    emxx = DTWO - sum(wt * be01x * unextx)
    emxy = -sum(wt * be01x * unexty)
    emyx = -sum(wt * be01y * unextx)
    emyy = DTWO - sum(wt * be01y * unexty)
    rx = sum(wt * (bi0x * vm0i + bi1x * vm1i + be01x * this%vne))
    ry = sum(wt * (bi0y * vm0i + bi1y * vm1i + be01y * this%vne))
    emdet = emxx * emyy - emxy * emyx
    this%vctrx = (emyy * rx - emxy * ry) / emdet
    this%vctry = (emxx * ry - emyx * rx) / emdet

    ! Get vm0e's using "known" conditions
    vm0e = 5.d-1 * (this%vne + unextx * this%vctrx + unexty * this%vctry)

    ! Get vm1e's from uniformity along exterior edges
    vm1e = vm0e

    ! Transform vm0 and vm1 to (x, y) coordinates
    vm0x = bi0x * vm0i + be0x * vm0e
    vm0y = bi0y * vm0i + be0y * vm0e
    vm1x = bi1x * vm1i + be1x * vm0e
    vm1y = bi1y * vm1i + be1y * vm0e

    ! Calculate head-cycle summation (which is proportional to
    ! the curl of the head gradient)
    hcsum = sum(lm * (kappax * (vm0x + vm1x) + kappay * (vm0y + vm1y)))

    ! Deallocate local arrays
    deallocate (vm0i)
    deallocate (vm0e)
    deallocate (vm1i)
    deallocate (vm1e)
    deallocate (uprod)
    deallocate (det)
    deallocate (wt)
    deallocate (bi0x)
    deallocate (be0x)
    deallocate (bi0y)
    deallocate (be0y)
    deallocate (bi1x)
    deallocate (be1x)
    deallocate (bi1y)
    deallocate (be1y)
    deallocate (be01x)
    deallocate (be01y)

  end subroutine calc_thru_hcsum

end module MethodCellTernaryModule

