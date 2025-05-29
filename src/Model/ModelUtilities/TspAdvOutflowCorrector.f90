module TspAdvOutflowCorrectorModule
  use TspFmiModule, only: TspFmiType
  use IGradient, only: IGradientType
  use BaseDisModule, only: DisBaseType
  use KindModule, only: DP, I4B
  use SimModule, only: store_warning
  use BoundaryFacesModule, only: BoundaryFacesType
  use ExtendedLeastSquaredGradientModule, only: ExtendedLeastSquaredGradientType
  use ConstantsModule, only: DONE, DZERO
  use PackageBudgetModule, only: PackageBudgetType

  implicit none
  private

  public :: TspAdvOutflowCorrectorType

  !> @brief High-order outflow corrector for advection transport.
  !!
  !! This type computes the high-order correction term for sink (outflow) terms in
  !! in the transport model. The total outflow
  !! can be written as the sum of the original lower-order (upstream) treatment
  !! and a new higher-order contribution. This class is responsible for computing
  !! only the higher-order contribution, which is then added to the standard sink treatment.
  !!
  !! The regular sink treatment uses upstream information and assumes the cell
  !! concentration is constant. With high-order cell value reconstruction, the
  !! outflow can be determined more accurately by accounting for local gradients
  !! within the cell. This module implements such a correction, improving the
  !! accuracy of sink (outflow) terms in the transport model.
  !!
  !! The corrector is only applied when a high-order interpolation scheme is in use
  !! (i.e., when `iadvwt >= 2`). In addition, the specific discharge must be available
  !! in the flow model. At runtime, the corrector checks these requirements and will
  !! issue a warning if high-order correction cannot be applied.
  !<
  type TspAdvOutflowCorrectorType
    private
    class(DisBaseType), pointer :: dis !< Pointer to discretization base type object
    type(TspFmiType), pointer :: fmi => null() !< Pointer to FMI (Flow Model Interface) object
    integer(I4B), pointer :: iadvwt => null() !< Pointer to integer array for advection weighting
    real(DP), pointer :: eqnsclfac => null() !< Pointer to equation scaling factor

    type(BoundaryFacesType), allocatable :: boundary_faces !< Contains information about boundary faces
    class(IGradientType), allocatable :: gradient
    logical :: is_initialized = .false. !< flag to indicate if the corrector is initialized
    logical :: correction_enabled = .false. !< flag to indicate if the corrector can compute the outflow correction
  contains
    procedure :: fc !< method to calculate the outflow correction

    procedure, private :: apply_correction_for_package
    procedure, private :: initialize_if_needed
    procedure, private :: get_cell_position

  end type TspAdvOutflowCorrectorType

  interface TspAdvOutflowCorrectorType
    module procedure Constructor
  end interface TspAdvOutflowCorrectorType

contains

  function Constructor(dis, fmi, iadvwt, eqnsclfac) result(corrector)
    ! -- return
    type(TspAdvOutflowCorrectorType) :: corrector
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    type(TspFmiType), pointer, intent(in) :: fmi
    integer(I4B), pointer, intent(in) :: iadvwt
    real(DP), pointer, intent(in) :: eqnsclfac

    corrector%dis => dis
    corrector%fmi => fmi
    corrector%iadvwt => iadvwt
    corrector%eqnsclfac => eqnsclfac

  end function Constructor

  subroutine fc(this, cnew, rhs)
    ! -- dummy
    class(TspAdvOutflowCorrectorType) :: this
    real(DP), intent(in), dimension(:) :: cnew
    real(DP), dimension(:), intent(inout) :: rhs
    ! -- local
    integer(I4B) :: ip
    type(PackageBudgetType), pointer :: package

    call this%initialize_if_needed()
    if (.not. this%correction_enabled) return
    !
    ! Compute the higher order contribution for each package
    do ip = 1, this%fmi%nflowpack
      package => this%fmi%gwfpackages(ip)
      if (this%fmi%iatp(ip) /= 0) cycle
      call this%apply_correction_for_package(package, cnew, rhs)
    end do
  end subroutine fc

  subroutine initialize_if_needed(this)
    ! -- dummy
    class(TspAdvOutflowCorrectorType) :: this

    if (this%is_initialized) return

    ! The gwfspdis is initialized after the models are initialized. Therefore we have to
    ! do a runtime check to see if the gwfspdis is associated before we can use it.
    if (this%iadvwt >= 2 .and. associated(this%fmi%gwfspdis)) then
      this%gradient = ExtendedLeastSquaredGradientType(this%dis, this%fmi)
      this%boundary_faces = BoundaryFacesType(this%dis)
      this%correction_enabled = .true.
    else
      call store_warning( &
        'High order sinks cannot be applied due to missing specific '// &
        'discharge. To use it enable the save_specific_discharge '// &
        'flag in the NPF package.' &
        )
      this%correction_enabled = .false.
    end if

    this%is_initialized = .true.

  end subroutine initialize_if_needed

  !!> @brief Apply high-order outflow correction for a single flow package.
  !!
  !! This subroutine computes and applies the high-order correction term for sink (outflow)
  !! terms associated with a specific flow package. It loops over all boundary cells in the
  !! package, calculates the outflow correction based on local gradients and face geometry,
  !! and updates the right-hand side vector accordingly.
  !!
  !! The correction is computed by assuming that the flow is leaving through the boundary
  !! of the domain, allowing a higher-order outflow to be constructed. In reality, the flow
  !! is not leaving through the physical boundary but through a sink; therefore, the computed
  !! "boundary flux" is added to the sink term to improve accuracy.
  !!
  !<
  subroutine apply_correction_for_package(this, package, cnew, rhs)
    ! -- dummy
    class(TspAdvOutflowCorrectorType) :: this
    type(PackageBudgetType), pointer :: package
    real(DP), intent(in) :: cnew(:)
    real(DP), intent(inout) :: rhs(:)
    ! -- local
    integer(I4B) :: i, n, ipos
    real(DP), dimension(3) :: normal !< Normal vector of the boundary face
    real(DP), dimension(3) :: q !< Specific discharge vector at cell center
    real(DP), dimension(3) :: grad_c !< Concentration/temperature gradient at cell center
    real(DP), dimension(3) :: dnm !< Vector from cell center to face center
    real(DP), dimension(3) :: xn !< Position of cell center
    real(DP), dimension(3) :: xf !< Position of face center
    real(DP) :: area !< Area of the boundary face
    real(DP) :: qn !< Normal component of specific discharge at the face
    real(DP) :: flow !< Sink (outflow) rate for the boundary cell
    real(DP) :: qnm !< Outflow through the face (scaled by area and eqnsclfac)

    if (package%naux > 0) return

    do i = 1, package%nbound
      n = package%nodelist(i)
      flow = package%flow(i)
      if (n <= 0) cycle
      if (flow >= DZERO) cycle

      xn = this%get_cell_position(n)
      !
      ! -- Loop over all boundary faces of the cell
      do ipos = this%boundary_faces%ia(n), this%boundary_faces%ia(n + 1) - 1
        area = this%boundary_faces%get_area(ipos)
        normal = this%boundary_faces%get_normal(ipos)
        !
        ! -- Compute the boundary flowja
        q = this%fmi%gwfspdis(:, n)
        qn = dot_product(q, normal)
        qnm = qn * area * this%eqnsclfac
        !
        ! -- Only apply correction if there is flow leaving the cell
        if (qnm <= DZERO) cycle
        !
        ! -- Get the face position and compute the correction term
        xf = this%boundary_faces%get_xf(ipos)
        grad_c = this%gradient%get(n, cnew)
        dnm = xf - xn
        !
        ! -- The low order term is already included in the rhs
        ! -- in the ssm package, so we only need to add the
        ! -- high order term here.
        rhs(n) = rhs(n) + qnm * dot_product(grad_c, dnm)
      end do
    end do
  end subroutine apply_correction_for_package

  function get_cell_position(this, n) result(xn)
    ! -- return
    real(DP), dimension(3) :: xn
    ! -- dummy
    class(TspAdvOutflowCorrectorType) :: this
    integer(I4B), intent(in) :: n
    ! -- local

    xn(1) = this%dis%xc(n)
    xn(2) = this%dis%yc(n)
    xn(3) = (this%dis%top(n) + this%dis%bot(n)) / 2.0_dp

  end function get_cell_position

end module TspAdvOutflowCorrectorModule
