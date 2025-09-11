module SpdisCellModule
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE, C3D_VERTICAL, DTWOPI, DSAME
  use BaseDisModule, only: DisBaseType
  use DisUtilsModule, only: number_connected_faces
  implicit none
  private

  !> @brief Container for spdis calculation in cell
  !!
  !! face index: 1 = top face, 2 = bottom face, and
  !< {3, .., nr_faces} vertical faces in clockwise order
  type SpdisCellType
    integer(I4B) :: nr_faces
    real(DP), dimension(:), allocatable :: hwva
    real(DP), dimension(:), allocatable :: anglex
    real(DP), dimension(:,:), allocatable :: n
    integer(I4B), dimension(:), allocatable :: matched
    ! flow data
    real(DP), dimension(:), pointer :: flowja => null() !< internal flows
    integer(I4B), dimension(:), pointer, contiguous :: ihc_face => null() !< face type (horizontal or vertical)
    real(DP), dimension(:, :), pointer, contiguous :: props_face => null() !< face properties (Q, area, nx, ny, distance)
  contains
    procedure :: create
    procedure :: load_cell
    procedure :: destroy
  end type SpdisCellType

contains

  subroutine create(this, dis, flowja, ihc, props)
    class(SpdisCellType) :: this
    class(DisBaseType), intent(inout) :: dis !< discretization
    real(DP), dimension(:), pointer :: flowja
    integer(I4B), dimension(:), pointer, contiguous :: ihc !< ihc flag for faces
    real(DP), dimension(:, :), pointer, contiguous :: props !< properties for faces
    ! local
    integer(I4B) :: max_nr_faces

    this%flowja => flowja
    this%ihc_face => ihc
    this%props_face => props

    ! max. nr. of faces is equal to max. nr. of vertices, plus top and bottom face
    max_nr_faces = dis%get_max_npolyverts(closed = .false.) + 2

    allocate(this%hwva(max_nr_faces))
    allocate(this%anglex(max_nr_faces))
    allocate(this%n(3, max_nr_faces))
    allocate(this%matched(max_nr_faces))

  end subroutine create

  subroutine load_cell(this, dis, n)
    class(SpdisCellType) :: this
    class(DisBaseType), intent(inout) :: dis
    integer(I4B), intent(in) :: n !< reduced node number
    ! local
    integer(I4B) :: ivert, ipos, m, iface
    integer(I4B) :: npoly
    real(DP), dimension(:,:), allocatable :: pverts
    real(DP), dimension(2) :: v
    real(DP) :: length

    ! reset
    this%matched = 0

    ! determine nr of faces
    npoly = dis%get_npolyverts(n, closed=.false.) 
    this%nr_faces = npoly + 2

    ! horizontal faces
    this%n(:, 1) = (/ DZERO, DZERO, DONE /) ! TODO_MJR: + or - ...
    this%hwva(1) = dis%area(n)
    this%n(:, 2) = (/ DZERO, DZERO, -DONE /)
    this%hwva(2) = dis%area(n)

    ! vertical faces
    call dis%get_polyverts(n, pverts)
    do ivert = 1, npoly
      iface = 2 + ivert ! start after horizontal faces

      v = pverts(:, modulo(ivert, npoly) + 1) - pverts(:,ivert)
      length = sqrt(v(1)*v(1) + v(2)*v(2))

      this%n(:, iface) = (/ v(2)/length, -v(1)/length, DZERO /)
      this%hwva(iface) = length
      this%anglex(iface) = acos(v(2)/length) ! between 0 and pi
      if (this%n(2, iface) < DZERO) then ! if y is negative, take complement
        this%anglex(iface) = DTWOPI - this%anglex(iface) ! TODO_MJR: MODFLOW takes the largest angle, no?
      end if
    end do

    ! match internal connections
    do ipos = dis%con%ia(n), dis%con%ia(n + 1) - 1
      m = dis%con%ja(ipos)
      if (dis%con%ihc(ipos) == C3D_VERTICAL) then
        if (m > n) then ! bottom
          this%matched(2) = 1
        else
          this%matched(1) = 1
        end if
      else
        ! match angle + hwva for horizontal connections
        do iface = 3, this%nr_faces
          if (this%matched(iface) == 1) cycle
          if (abs(dis%con%anglex - this%anglex(iface)) < DSAME) then
            if (abs(this%hwva(iface) - dis%con%hwva(iface)) < DSAME) then
              this%matched(iface) = 1
            end if
          end if
        end do
      end if
    end do

    ! match exchange flows

    ! what's left are boundary faces

  end subroutine load_cell

  subroutine destroy(this)
    class(SpdisCellType) :: this

    deallocate(this%hwva)
    deallocate(this%anglex)
    deallocate(this%n)
    deallocate(this%matched)


  end subroutine destroy

end module SpdisCellModule