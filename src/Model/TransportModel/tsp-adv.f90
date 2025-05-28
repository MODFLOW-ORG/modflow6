module TspAdvOutflowCorrectorModule
  use TspFmiModule, only: TspFmiType
  use IGradient, only: IGradientType
  use BaseDisModule, only: DisBaseType
  use KindModule, only: DP, I4B
  use SimModule, only: store_warning
  use BoundaryFacesModule, only: BoundaryFacesType
  use ExtendedLeastSquaredGradientModule, only: ExtendedLeastSquaredGradientType
  use ConstantsModule, only: DONE, DZERO

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
    ! -- dummy
    class(TspAdvOutflowCorrectorType) :: this
    real(DP), intent(in), dimension(:) :: cnew
    real(DP), dimension(:), intent(inout) :: rhs
    ! -- local
    integer(I4B) :: ip, n, i, ipos
    real(DP), dimension(3) :: normal, q, grad_c, dnm, xn, xf
    real(DP) :: area, qn, flow, qnm

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
      if (this%fmi%iatp(ip) /= 0) cycle
      do i = 1, this%fmi%gwfpackages(ip)%nbound
        n = this%fmi%gwfpackages(ip)%nodelist(i)
        flow = this%fmi%gwfpackages(ip)%flow(i)
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
          dnm = xf - xn

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

module TspAdvModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, DNODATA, DPREC, LINELENGTH
  use SimModule, only: store_warning
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use TspAdvOptionsModule, only: TspAdvOptionsType
  use MatrixBaseModule
  use IGradient, only: IGradientType
  use LeastSquaredGradientBoundaryModule, only: LeastSquaredGradientBoundaryType
  use IInterpolationSchemeModule, only: IInterpolationSchemeType, CoefficientsType
  use UpwindSchemeModule, only: UpwindSchemeType
  use CentralDifferenceSchemeModule, only: CentralDifferenceSchemeType
  use TVDSchemeModule, only: TVDSchemeType
  use BarthJespersenSchemeModule, only: BarthJespersenSchemeType
  use TspAdvOutflowCorrectorModule, only: TspAdvOutflowCorrectorType

  implicit none
  private
  public :: TspAdvType
  public :: adv_cr

  type, extends(NumericalPackageType) :: TspAdvType

    integer(I4B), pointer :: iadvwt => null() !< advection scheme (0 up, 1 central, 2 tvd)
    logical :: use_high_order_sinks = .false. !< Add a correction to the sink term
    real(DP), pointer :: ats_percel => null() !< user-specified fractional number of cells advection can move a particle during one time step
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    type(TspFmiType), pointer :: fmi => null() !< pointer to fmi object
    real(DP), pointer :: eqnsclfac => null() !< governing equation scale factor; =1. for solute; =rhow*cpw for energy

    class(IInterpolationSchemeType), allocatable :: face_interpolation
    class(IGradientType), allocatable :: gradient
    type(TspAdvOutflowCorrectorType), allocatable :: outflow_corrector
  contains

    procedure :: adv_df
    procedure :: adv_ar
    procedure :: adv_dt
    procedure :: adv_fc
    procedure :: adv_cq
    procedure :: adv_da

    procedure :: allocate_scalars
    procedure, private :: read_options

  end type TspAdvType

contains

  !> @ brief Create a new ADV object
  !!
  !!  Create a new ADV package
  !<
  subroutine adv_cr(advobj, name_model, inunit, iout, fmi, eqnsclfac)
    ! -- dummy
    type(TspAdvType), pointer :: advobj
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    type(TspFmiType), intent(in), target :: fmi
    real(DP), intent(in), pointer :: eqnsclfac !< governing equation scale factor
    !
    ! -- Create the object
    allocate (advobj)
    !
    ! -- create name and memory path
    call advobj%set_names(1, name_model, 'ADV', 'ADV')
    !
    ! -- Allocate scalars
    call advobj%allocate_scalars()
    !
    ! -- Set variables
    advobj%inunit = inunit
    advobj%iout = iout
    advobj%fmi => fmi
    advobj%eqnsclfac => eqnsclfac
  end subroutine adv_cr

  !> @brief Define ADV object
  !!
  !! Define the ADV package
  !<
  subroutine adv_df(this, adv_options)
    ! -- dummy
    class(TspAdvType) :: this
    type(TspAdvOptionsType), optional, intent(in) :: adv_options !< the optional options, for when not constructing from file
    ! -- local
    character(len=*), parameter :: fmtadv = &
      "(1x,/1x,'ADV-- ADVECTION PACKAGE, VERSION 1, 8/25/2017', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! -- Read or set advection options
    if (.not. present(adv_options)) then
      !
      ! -- Initialize block parser (adv has no define, so it's
      ! not done until here)
      call this%parser%Initialize(this%inunit, this%iout)
      !
      ! --print a message identifying the advection package.
      write (this%iout, fmtadv) this%inunit
      !
      ! --read options from file
      call this%read_options()
    else
      !
      ! --set options from input arg
      this%iadvwt = adv_options%iAdvScheme
    end if
  end subroutine adv_df

  !> @brief Allocate and read method for package
  !!
  !!  Method to allocate and read static data for the ADV package.
  !<
  subroutine adv_ar(this, dis, ibound)
    ! -- modules
    ! -- dummy
    class(TspAdvType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: ibound
    ! -- local
    integer(I4B) :: nodes
    ! -- formats
    !
    ! -- adv pointers to arguments that were passed in
    this%dis => dis
    this%ibound => ibound

    ! -- Compute the gradient operator
    nodes = dis%nodes

    this%gradient = LeastSquaredGradientBoundaryType(this%dis, this%fmi)

    ! -- Select interpolation scheme
    if (this%iadvwt == 0) then
      this%face_interpolation = &
        UpwindSchemeType(this%dis, this%fmi, this%gradient)
    elseif (this%iadvwt == 1) then
      this%face_interpolation = &
        CentralDifferenceSchemeType(this%dis, this%fmi, this%gradient)
    elseif (this%iadvwt == 2) then
      this%face_interpolation = &
        TVDSchemeType(this%dis, this%fmi, this%gradient)
    elseif (this%iadvwt == 3) then
      this%face_interpolation = &
        BarthJespersenSchemeType(this%dis, this%fmi, this%gradient)
    end if

    ! -- Determine if the high order sink term should be used
    this%outflow_corrector = &
      TspAdvOutflowCorrectorType(this%dis, this%fmi, this%iadvwt, this%eqnsclfac)

  end subroutine adv_ar

  !> @brief  Calculate maximum time step length
  !!
  !!  Return the largest time step that meets stability constraints
  !<
  subroutine adv_dt(this, dtmax, msg, thetam)
    ! dummy
    class(TspAdvType) :: this !< this instance
    real(DP), intent(out) :: dtmax !< maximum allowable dt subject to stability constraint
    character(len=*), intent(inout) :: msg !< package/cell dt constraint message
    real(DP), dimension(:), intent(in) :: thetam !< porosity
    ! local
    integer(I4B) :: n
    integer(I4B) :: m
    integer(I4B) :: ipos
    integer(I4B) :: nrmax
    character(len=LINELENGTH) :: cellstr
    real(DP) :: dt
    real(DP) :: flowmax
    real(DP) :: flowsumpos
    real(DP) :: flowsumneg
    real(DP) :: flownm
    real(DP) :: cell_volume
    dtmax = DNODATA
    nrmax = 0
    msg = ''

    ! If ats_percel not specified by user, then return without making
    ! the courant time step calculation
    if (this%ats_percel == DNODATA) then
      return
    end if

    ! Calculate time step lengths based on stability constraint for each cell
    ! and store the smallest one
    do n = 1, this%dis%nodes
      if (this%ibound(n) == 0) cycle
      flowsumneg = DZERO
      flowsumpos = DZERO
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        if (this%dis%con%mask(ipos) == 0) cycle
        m = this%dis%con%ja(ipos)
        if (this%ibound(m) == 0) cycle
        flownm = this%fmi%gwfflowja(ipos)
        if (flownm < DZERO) then
          flowsumneg = flowsumneg - flownm
        else
          flowsumpos = flowsumpos + flownm
        end if
      end do
      flowmax = max(flowsumneg, flowsumpos)
      if (flowmax < DPREC) cycle
      cell_volume = this%dis%get_cell_volume(n, this%dis%top(n))
      dt = cell_volume * this%fmi%gwfsat(n) * thetam(n) / flowmax
      dt = dt * this%ats_percel
      if (dt < dtmax) then
        dtmax = dt
        nrmax = n
      end if
    end do
    if (nrmax > 0) then
      call this%dis%noder_to_string(nrmax, cellstr)
      write (msg, *) adjustl(trim(this%memoryPath))//'-'//trim(cellstr)
    end if
  end subroutine adv_dt

  !> @brief  Fill coefficient method for ADV package
  !!
  !!  Method to calculate coefficients and fill amat and rhs.
  !<
  subroutine adv_fc(this, nodes, matrix_sln, idxglo, cnew, rhs)
    ! -- modules
    ! -- dummy
    class(TspAdvType) :: this
    integer, intent(in) :: nodes
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(in), dimension(:) :: cnew
    real(DP), dimension(:), intent(inout) :: rhs
    ! -- local
    integer(I4B) :: n, m, idiag, ipos
    real(DP) :: qnm
    type(CoefficientsType) :: coefficients

    do n = 1, nodes
      if (this%ibound(n) == 0) cycle
      idiag = this%dis%con%ia(n)
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        if (this%ibound(m) == 0) cycle
        if (this%dis%con%mask(ipos) == 0) cycle

        qnm = this%fmi%gwfflowja(ipos) * this%eqnsclfac
        coefficients = this%face_interpolation%compute(n, m, ipos, cnew)

        call matrix_sln%add_value_pos(idxglo(idiag), qnm * coefficients%c_n)
        call matrix_sln%add_value_pos(idxglo(ipos), qnm * coefficients%c_m)
        rhs(n) = rhs(n) + qnm * coefficients%rhs
      end do
    end do

    ! -- Calculate and add high order flux contribution for sinks
    call this%outflow_corrector%fc(cnew, rhs)

  end subroutine adv_fc

  !> @brief Calculate advection contribution to flowja
  !<
  subroutine adv_cq(this, cnew, flowja)
    ! -- modules
    ! -- dummy
    class(TspAdvType) :: this
    real(DP), intent(in), dimension(:) :: cnew
    real(DP), intent(inout), dimension(:) :: flowja
    ! -- local
    integer(I4B) :: nodes
    integer(I4B) :: n, m, ipos
    real(DP) :: qnm
    type(CoefficientsType) :: coefficients
    !
    ! -- Calculate advection and add to flowja. qnm is the volumetric flow
    !    rate and has dimensions of L^/T.
    nodes = this%dis%nodes

    do n = 1, nodes
      if (this%ibound(n) == 0) cycle
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        if (this%ibound(m) == 0) cycle
        qnm = this%fmi%gwfflowja(ipos) * this%eqnsclfac

        coefficients = this%face_interpolation%compute(n, m, ipos, cnew)
        flowja(ipos) = flowja(ipos) &
                       + qnm * coefficients%c_n * cnew(n) &
                       + qnm * coefficients%c_m * cnew(m) &
                       - qnm * coefficients%rhs
      end do
    end do

  end subroutine adv_cq

  !> @brief Deallocate memory
  !<
  subroutine adv_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(TspAdvType) :: this
    !
    ! -- Deallocate arrays if package was active
    if (this%inunit > 0) then
    end if
    !
    ! -- nullify pointers
    this%ibound => null()
    !
    ! -- Scalars
    call mem_deallocate(this%iadvwt)
    call mem_deallocate(this%ats_percel)
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()
  end subroutine adv_da

  !> @brief Allocate scalars specific to the streamflow energy transport (SFE)
  !! package.
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy
    class(TspAdvType) :: this
    ! -- local
    !
    ! -- allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%iadvwt, 'IADVWT', this%memoryPath)
    call mem_allocate(this%ats_percel, 'ATS_PERCEL', this%memoryPath)
    !
    ! -- Initialize
    this%iadvwt = 0
    this%ats_percel = DNODATA
    !
    ! -- Advection creates an asymmetric coefficient matrix
    this%iasym = 1
  end subroutine allocate_scalars

  !> @brief Read options
  !!
  !! Read the options block
  !<
  subroutine read_options(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error
    ! -- dummy
    class(TspAdvType) :: this
    ! -- local
    character(len=LINELENGTH) :: errmsg, keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    ! -- formats
    character(len=*), parameter :: fmtiadvwt = &
      &"(4x,'ADVECTION WEIGHTING SCHEME HAS BEEN SET TO: ', a)"
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, blockRequired=.false., &
                              supportOpenClose=.true.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING ADVECTION OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('SCHEME')
          call this%parser%GetStringCaps(keyword)
          select case (keyword)
          case ('UPSTREAM')
            this%iadvwt = 0
            write (this%iout, fmtiadvwt) 'UPSTREAM'
          case ('CENTRAL')
            this%iadvwt = 1
            write (this%iout, fmtiadvwt) 'CENTRAL'
          case ('TVD')
            this%iadvwt = 2
            write (this%iout, fmtiadvwt) 'TVD'
          case ('BARTH')
            this%iadvwt = 3
            write (this%iout, fmtiadvwt) 'BARTH'
          case default
            write (errmsg, '(a, a)') &
              'Unknown scheme: ', trim(keyword)
            call store_error(errmsg)
            write (errmsg, '(a, a)') &
              'Scheme must be "UPSTREAM", "CENTRAL" or "TVD"'
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end select
        case ('ATS_PERCEL')
          this%ats_percel = this%parser%GetDouble()
          if (this%ats_percel == DZERO) this%ats_percel = DNODATA
          write (this%iout, '(4x,a,1pg15.6)') &
            'User-specified fractional cell distance for adaptive time &
            &steps: ', this%ats_percel
        case default
          write (errmsg, '(a,a)') 'Unknown ADVECTION option: ', &
            trim(keyword)
          call store_error(errmsg, terminate=.TRUE.)
        end select
      end do
      write (this%iout, '(1x,a)') 'END OF ADVECTION OPTIONS'
    end if
  end subroutine read_options

end module TspAdvModule
