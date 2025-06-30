module MethodCellPollockQuadModule

  use KindModule, only: DP, I4B
  use ErrorUtilModule, only: pstop
  use ConstantsModule, only: DONE, DZERO
  use MethodModule, only: MethodType
  use MethodCellModule, only: MethodCellType
  use MethodSubcellPoolModule, only: method_subcell_plck
  use CellRectQuadModule, only: CellRectQuadType, create_cell_rect_quad
  use CellDefnModule, only: CellDefnType
  use SubcellRectModule, only: SubcellRectType, create_subcell_rect
  use ParticleModule, only: ParticleType
  implicit none

  private
  public :: MethodCellPollockQuadType
  public :: create_method_cell_quad

  type, extends(MethodCellType) :: MethodCellPollockQuadType
  contains
    procedure, public :: apply => apply_mcpq
    procedure, public :: deallocate
    procedure, public :: load => load_mcpq
    procedure, public :: load_subcell
    procedure, public :: pass => pass_mcpq
  end type MethodCellPollockQuadType

contains

  !> @brief Create a new Pollock quad-refined cell method
  subroutine create_method_cell_quad(method)
    ! dummy
    type(MethodCellPollockQuadType), pointer :: method
    ! local
    type(CellRectQuadType), pointer :: cell
    type(SubcellRectType), pointer :: subcell

    allocate (method)
    call create_cell_rect_quad(cell)
    method%cell => cell
    method%name => method%cell%type
    method%delegates = .true.
    call create_subcell_rect(subcell)
    method%subcell => subcell
  end subroutine create_method_cell_quad

  !> @brief Deallocate the Pollock quad-refined cell method
  subroutine deallocate (this)
    class(MethodCellPollockQuadType), intent(inout) :: this
    deallocate (this%name)
  end subroutine deallocate

  !> @brief Load subcell into tracking method
  subroutine load_mcpq(this, particle, next_level, submethod)
    class(MethodCellPollockQuadType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: next_level
    class(MethodType), pointer, intent(inout) :: submethod

    select type (subcell => this%subcell)
    type is (SubcellRectType)
      call this%load_subcell(particle, subcell)
    end select
    call method_subcell_plck%init( &
      fmi=this%fmi, &
      cell=this%cell, &
      subcell=this%subcell, &
      trackctl=this%trackctl, &
      tracktimes=this%tracktimes)
    submethod => method_subcell_plck
  end subroutine load_mcpq

  !> @brief Pass particle to next subcell if there is one, or to the cell face
  subroutine pass_mcpq(this, particle)
    ! dummy
    class(MethodCellPollockQuadType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    ! local
    integer(I4B) :: isc, exitFace, npolyverts, inface, infaceoff

    select type (cell => this%cell)
    type is (CellRectQuadType)
      exitFace = particle%iboundary(3)
      isc = particle%idomain(3)
      npolyverts = cell%defn%npolyverts

      ! exitFace uses MODPATH 7 iface convention here
      select case (exitFace)
      case (0)
        ! Subcell interior (cell interior)
        inface = -1
      case (1)
        select case (isc)
        case (1)
          ! W face, subcell 1 --> E face, subcell 4  (cell interior)
          particle%idomain(3) = 4
          particle%iboundary(3) = 2
          inface = 0 ! want Domain(2) unchanged; Boundary(2) = 0
        case (2)
          ! W face, subcell 2 --> E face, subcell 3 (cell interior)
          particle%idomain(3) = 3
          particle%iboundary(3) = 2
          inface = 0 ! want Domain(2) unchanged; Boundary(2) = 0
        case (3)
          ! W face, subcell 3 (cell face)
          inface = 1 ! want Domain(2) = -Domain(2); Boundary(2) = inface
          infaceoff = 0
        case (4)
          ! W face, subcell 4 (cell face)
          inface = 2 ! want Domain(2) = -Domain(2); Boundary(2) = inface
          infaceoff = -1
        end select
      case (2)
        select case (isc)
        case (1)
          ! E face, subcell 1 (cell face)
          inface = 3 ! want Domain(2) = -Domain(2); Boundary(2) = inface
          infaceoff = 0
        case (2)
          ! E face, subcell 2 (cell face)
          inface = 4 ! want Domain(2) = -Domain(2); Boundary(2) = inface
          infaceoff = -1
        case (3)
          ! E face, subcell 3 --> W face, subcell 2 (cell interior)
          particle%idomain(3) = 2
          particle%iboundary(3) = 1
          inface = 0 ! want Domain(2) unchanged; Boundary(2) = 0
        case (4)
          ! E face, subcell 4 --> W face subcell 1 (cell interior)
          particle%idomain(3) = 1
          particle%iboundary(3) = 1
          inface = 0 ! want Domain(2) unchanged; Boundary(2) = 0
        end select
      case (3)
        select case (isc)
        case (1)
          ! S face, subcell 1 --> N face, subcell 2 (cell interior)
          particle%idomain(3) = 2
          particle%iboundary(3) = 4
          inface = 0 ! want Domain(2) unchanged; Boundary(2) = 0
        case (2)
          ! S face, subcell 2 (cell face)
          inface = 4 ! want Domain(2) = -Domain(2); Boundary(2) = inface
          infaceoff = 0
        case (3)
          ! S face, subcell 3 (cell face)
          inface = 1 ! want Domain(2) = -Domain(2); Boundary(2) = inface
          infaceoff = -1
        case (4)
          ! S face, subcell 4 --> N face, subcell 3 (cell interior)
          particle%idomain(3) = 3
          particle%iboundary(3) = 4
          inface = 0 ! want Domain(2) unchanged; Boundary(2) = 0
        end select
      case (4)
        select case (isc)
        case (1)
          ! N face, subcell 1 (cell face)
          inface = 3 ! want Domain(2) = -Domain(2); Boundary(2) = inface
          infaceoff = -1
        case (2)
          ! N face, subcell 2 --> S face, subcell 1 (cell interior)
          particle%idomain(3) = 1
          particle%iboundary(3) = 3
          inface = 0 ! want Domain(2) unchanged; Boundary(2) = 0
        case (3)
          ! N face, subcell 3 --> S face, subcell 4 (cell interior)
          particle%idomain(3) = 4
          particle%iboundary(3) = 3
          inface = 0 ! want Domain(2) unchanged; Boundary(2) = 0
        case (4)
          ! N face, subcell 4 (cell face)
          inface = 2 ! want Domain(2) = -Domain(2); Boundary(2) = inface
          infaceoff = 0
        end select
      case (5)
        ! Subcell bottom (cell bottom)
        inface = npolyverts + 2 ! want Domain(2) = -Domain(2); Boundary(2) = inface
      case (6)
        ! Subcell top (cell top)
        inface = npolyverts + 3 ! want Domain(2) = -Domain(2); Boundary(2) = inface
      end select

      if (inface .eq. -1) then
        particle%iboundary(2) = 0
      else if (inface .eq. 0) then
        particle%iboundary(2) = 0
      else
        if ((inface .ge. 1) .and. (inface .le. 4)) then
          ! Account for local cell rotation
          inface = inface + cell%irvOrigin - 1
          if (inface .gt. 4) inface = inface - 4
          inface = cell%irectvert(inface) + infaceoff
          if (inface .lt. 1) inface = inface + npolyverts
        end if
        particle%iboundary(2) = inface
      end if
    end select
  end subroutine pass_mcpq

  !> @brief Solve the quad-rectangular cell via Pollock's method
  subroutine apply_mcpq(this, particle, tmax)
    ! dummy
    class(MethodCellPollockQuadType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local
    double precision :: xOrigin, yOrigin, zOrigin, sinrot, cosrot

    select type (cell => this%cell)
    type is (CellRectQuadType)
      ! Check termination/reporting conditions
      call this%check(particle, cell%defn, tmax)
      if (.not. particle%advancing) return

      ! Transform model coordinates to local cell coordinates
      ! (translated/rotated but not scaled relative to model)
      xOrigin = cell%xOrigin
      yOrigin = cell%yOrigin
      zOrigin = cell%zOrigin
      sinrot = cell%sinrot
      cosrot = cell%cosrot
      call particle%transform(xOrigin, yOrigin, zOrigin, &
                              sinrot, cosrot)

      ! Track the particle over the cell
      call this%track(particle, 2, tmax)

      ! Transform cell coordinates back to model coordinates
      call particle%transform(xOrigin, yOrigin, zOrigin, &
                              sinrot, cosrot, invert=.true.)
      call particle%reset_transform()
    end select
  end subroutine apply_mcpq

  !> @brief Load the rectangular subcell from the rectangular cell
  subroutine load_subcell(this, particle, subcell)
    ! dummy
    class(MethodCellPollockQuadType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(SubcellRectType), intent(inout) :: subcell
    ! local
    real(DP) :: dx, dy, dz, areax, areay, areaz
    real(DP) :: dxprel, dyprel
    integer(I4B) :: isc, npolyverts, m1, m2
    real(DP) :: qextl1, qextl2, qintl1, qintl2
    real(DP) :: factor, term

    select type (cell => this%cell)
    type is (CellRectQuadType)
      factor = DONE / cell%defn%retfactor
      factor = factor / cell%defn%porosity
      npolyverts = cell%defn%npolyverts

      isc = particle%idomain(3)
      ! Subcells 1, 2, 3, and 4 are Pollock's subcells A, B, C, and D,
      ! respectively

      dx = cell%dx
      dy = cell%dy
      ! If not already known, determine subcell number
      if (isc .le. 0) then
        dxprel = particle%x / dx
        dyprel = particle%y / dy

        if (dyprel .ge. 5d-1) then
          if (dxprel .le. 5d-1) then
            isc = 4
          else
            isc = 1
          end if
        else
          if (dxprel .le. 5d-1) then
            isc = 3
          else
            isc = 2
          end if
        end if

        subcell%isubcell = isc
        particle%idomain(3) = isc
      end if
      dx = 5d-1 * dx
      dy = 5d-1 * dy
      dz = cell%defn%top - &
           cell%defn%bot
      areax = dy * dz
      areay = dx * dz
      areaz = dx * dy
      qintl1 = cell%qintl(isc)
      ! qintl list wraps around, so isc+1=5 is ok
      qintl2 = cell%qintl(isc + 1)
      qextl1 = cell%qextl1(isc)
      qextl2 = cell%qextl2(isc)

      subcell%dx = dx
      subcell%dy = dy
      subcell%dz = dz
      subcell%sinrot = DZERO
      subcell%cosrot = DONE
      subcell%zOrigin = DZERO
      select case (isc)
      case (1)
        subcell%xOrigin = dx
        subcell%yOrigin = dy
        term = factor / areax
        subcell%vx1 = qintl1 * term
        subcell%vx2 = -qextl2 * term
        term = factor / areay
        subcell%vy1 = -qintl2 * term
        subcell%vy2 = -qextl1 * term
      case (2)
        subcell%xOrigin = dx
        subcell%yOrigin = DZERO
        term = factor / areax
        subcell%vx1 = -qintl2 * term
        subcell%vx2 = -qextl1 * term
        term = factor / areay
        subcell%vy1 = qextl2 * term
        subcell%vy2 = -qintl1 * term
      case (3)
        subcell%xOrigin = DZERO
        subcell%yOrigin = DZERO
        term = factor / areax
        subcell%vx1 = qextl2 * term
        subcell%vx2 = -qintl1 * term
        term = factor / areay
        subcell%vy1 = qextl1 * term
        subcell%vy2 = qintl2 * term
      case (4)
        subcell%xOrigin = DZERO
        subcell%yOrigin = dy
        term = factor / areax
        subcell%vx1 = qextl1 * term
        subcell%vx2 = qintl2 * term
        term = factor / areay
        subcell%vy1 = qintl1 * term
        subcell%vy2 = -qextl2 * term
      end select
      m1 = npolyverts + 2
      m2 = m1 + 1
      term = factor / areaz
      subcell%vz1 = 2.5d-1 * cell%defn%faceflow(m1) * term
      subcell%vz2 = -2.5d-1 * cell%defn%faceflow(m2) * term
    end select
  end subroutine load_subcell

end module MethodCellPollockQuadModule
