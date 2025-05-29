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

  type TspAdvOutflowCorrectorType
    private
    class(DisBaseType), pointer :: dis
    type(TspFmiType), pointer :: fmi => null() !< pointer to fmi object
    integer(I4B), pointer :: iadvwt => null()
    real(DP), pointer :: eqnsclfac => null()

    type(BoundaryFacesType), allocatable :: boundary_faces
    class(IGradientType), allocatable :: gradient2
    logical :: is_initialized = .false. !< flag to indicate if the corrector is initialized
    logical :: can_compute = .false. !< flag to indicate if the corrector can compute the outflow correction
  contains
    procedure :: fc !< method to calculate the outflow correction

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
    use TdisModule, only: delt
    ! -- dummy
    class(TspAdvOutflowCorrectorType) :: this
    real(DP), intent(in), dimension(:) :: cnew
    real(DP), dimension(:), intent(inout) :: rhs
    ! -- local
    integer(I4B) :: ip, n, i, ipos
    real(DP), dimension(3) :: normal, q, grad_c, dnm, xn, xf
    real(DP) :: area, qn, flow, qnm
    type(PackageBudgetType), pointer :: package

    ! The gwfspdis in initialized after the models are initialized. Therefore we have to
    ! do a runtime check to see if the fmi is associated before we can use it.
    if (.not. this%is_initialized) then
      if (this%iadvwt >= 2 .and. associated(this%fmi%gwfspdis)) then
        this%gradient2 = ExtendedLeastSquaredGradientType(this%dis, this%fmi)
        this%boundary_faces = BoundaryFacesType(this%dis)

        this%can_compute = .true.
      else
        call store_warning( &
          'High order sinks cannot be applied due to missing specific '// &
          'discharge. To use it enable the save_specific_discharge '// &
          'flag in the NPF package.' &
          )

        this%can_compute = .false.
      end if

      this%is_initialized = .true.
    end if

    ! Nothing to do. Return
    if (.not. this%can_compute) return

    ! Compute the higher order terms
    do ip = 1, this%fmi%nflowpack
      package => this%fmi%gwfpackages(ip)
      if (this%fmi%iatp(ip) /= 0) cycle
      ! If there is an auxiliary variable, then the concentration may be fixed in the cell,
      ! so we cannot use the high order correction
      if (package%naux > 0) cycle

      do i = 1, package%nbound
        n = package%nodelist(i)
        flow = package%flow(i)
        if (n <= 0) cycle
        if (flow >= DZERO) cycle

        xn = this%get_cell_position(n)
        do ipos = this%boundary_faces%ia(n), this%boundary_faces%ia(n + 1) - 1
          area = this%boundary_faces%get_area(ipos)
          normal = this%boundary_faces%get_normal(ipos)
          xf = this%boundary_faces%get_xf(ipos)

          q = this%fmi%gwfspdis(:, n)
          qn = dot_product(q, normal)
          qnm = qn * area * this%eqnsclfac

          if (qnm <= DZERO) cycle

          ! High order flux
          grad_c = this%gradient2%get(n, cnew)
          dnm = xf - xn - 0.5_dp * normal * qn * delt

          rhs(n) = rhs(n) + qnm * dot_product(grad_c, dnm)
        end do
        !
      end do
    end do
  end subroutine fc

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