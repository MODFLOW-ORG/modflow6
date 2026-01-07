module TspAdvModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, DNODATA, DPREC, LINELENGTH
  use FeatureFlagsModule, only: developmode
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use TspAdvOptionsModule, only: TspAdvOptionsType
  use MatrixBaseModule, only: MatrixBaseType
  ! -- Gradient schemes
  use IGradient, only: IGradientType
  use LeastSquaresGradientModule, only: LeastSquaresGradientType
  use CachedGradientModule, only: CachedGradientType
  ! -- Interpolation schemes
  use InterpolationSchemeInterfaceModule, only: InterpolationSchemeInterface, &
                                                CoefficientsType
  use AdvSchemeEnumModule
  use UpstreamSchemeModule, only: UpstreamSchemeType
  use CentralDifferenceSchemeModule, only: CentralDifferenceSchemeType
  use TVDSchemeModule, only: TVDSchemeType
  use UTVDSchemeModule, only: UTVDSchemeType

  implicit none
  private
  public :: TspAdvType
  public :: adv_cr

  type, extends(NumericalPackageType) :: TspAdvType
    integer(I4B), pointer :: iadvwt => null() !< advection scheme. See ADV_SCHEME_* constants
    real(DP), pointer :: ats_percel => null() !< user-specified fractional number of cells advection can move a particle during one time step
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    type(TspFmiType), pointer :: fmi => null() !< pointer to fmi object
    real(DP), pointer :: eqnsclfac => null() !< governing equation scale factor; =1. for solute; =rhow*cpw for energy

    class(InterpolationSchemeInterface), allocatable :: face_interpolation !< interpolation scheme for face values
    class(IGradientType), allocatable :: gradient !< cell centered gradient
  contains

    procedure :: adv_df
    procedure :: adv_ar
    procedure :: adv_dt
    procedure :: adv_fc
    procedure :: adv_cq
    procedure :: adv_da

    procedure :: allocate_scalars
    procedure, private :: source_options

  end type TspAdvType

contains

  !> @ brief Create a new ADV object
  !!
  !!  Create a new ADV package
  !<
  subroutine adv_cr(advobj, name_model, input_mempath, inunit, iout, fmi, &
                    eqnsclfac)
    ! -- dummy
    type(TspAdvType), pointer :: advobj
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    type(TspFmiType), intent(in), target :: fmi
    real(DP), intent(in), pointer :: eqnsclfac !< governing equation scale factor
    !
    ! -- Create the object
    allocate (advobj)
    !
    ! -- create name and memory path
    call advobj%set_names(1, name_model, 'ADV', 'ADV', input_mempath)
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
      "(1x,/1x,'ADV -- ADVECTION PACKAGE, VERSION 1, 8/25/2017', &
      &' INPUT READ FROM MEMPATH: ', A, //)"
    !
    ! -- Read or set advection options
    if (.not. present(adv_options)) then
      !
      ! --print a message identifying the advection package.
      write (this%iout, fmtadv) this%input_mempath
      !
      ! --read options from file
      call this%source_options()
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
    use SimModule, only: store_error
    ! -- dummy
    class(TspAdvType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: ibound
    ! -- local
    integer(I4B) :: iadvwt_value
    class(IGradientType), allocatable :: gradient
    ! -- adv pointers to arguments that were passed in
    this%dis => dis
    this%ibound => ibound
    !
    ! -- Create interpolation scheme
    iadvwt_value = this%iadvwt ! Dereference iadvwt to work with case statement
    select case (iadvwt_value)
    case (ADV_SCHEME_UPSTREAM)
      this%face_interpolation = &
        UpstreamSchemeType(this%dis, this%fmi)
    case (ADV_SCHEME_CENTRAL)
      this%face_interpolation = &
        CentralDifferenceSchemeType(this%dis, this%fmi)
    case (ADV_SCHEME_TVD)
      this%face_interpolation = &
        TVDSchemeType(this%dis, this%fmi, this%ibound)
    case (ADV_SCHEME_UTVD)
      gradient = LeastSquaresGradientType(this%dis)
      this%gradient = CachedGradientType(gradient, this%dis)
      this%face_interpolation = &
        UTVDSchemeType(this%dis, this%fmi, this%gradient)
    case default
      call store_error("Unknown advection scheme", terminate=.TRUE.)
    end select
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
    class(TspAdvType) :: this !< this instance
    integer(I4B), intent(in) :: nodes !< number of nodes
    class(MatrixBaseType), pointer :: matrix_sln !< pointer to solution matrix
    integer(I4B), intent(in), dimension(:) :: idxglo !< global indices for matrix
    real(DP), intent(in), dimension(:), target :: cnew !< new concentration/temperature values
    real(DP), dimension(:), intent(inout) :: rhs !< right-hand side vector
    ! -- local
    integer(I4B) :: n, m, idiag, ipos
    real(DP) :: qnm !< volumetric flow rate
    type(CoefficientsType) :: coefficients

    ! Calculate internal domain fluxes and add to matrix_sln and rhs.
    call this%face_interpolation%set_field(cnew)
    do n = 1, nodes
      if (this%ibound(n) == 0) cycle ! skip inactive nodes
      idiag = this%dis%con%ia(n)
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        if (this%ibound(m) == 0) cycle ! skip inactive nodes
        if (this%dis%con%mask(ipos) == 0) cycle ! skip masked connections

        qnm = this%fmi%gwfflowja(ipos) * this%eqnsclfac
        coefficients = this%face_interpolation%compute(n, m, ipos)

        call matrix_sln%add_value_pos(idxglo(idiag), qnm * coefficients%c_n)
        call matrix_sln%add_value_pos(idxglo(ipos), qnm * coefficients%c_m)
        rhs(n) = rhs(n) + qnm * coefficients%rhs
      end do
    end do
  end subroutine adv_fc

  !> @brief Calculate advection contribution to flowja
  !<
  subroutine adv_cq(this, cnew, flowja)
    ! -- modules
    ! -- dummy
    class(TspAdvType) :: this
    real(DP), intent(in), dimension(:), target :: cnew
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

    call this%face_interpolation%set_field(cnew)
    do n = 1, nodes
      if (this%ibound(n) == 0) cycle
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        if (this%ibound(m) == 0) cycle
        qnm = this%fmi%gwfflowja(ipos) * this%eqnsclfac

        coefficients = this%face_interpolation%compute(n, m, ipos)
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
    this%iadvwt = ADV_SCHEME_UPSTREAM
    this%ats_percel = DNODATA
    !
    ! -- Advection creates an asymmetric coefficient matrix
    this%iasym = 1
  end subroutine allocate_scalars

  !> @brief Source input options
  !<
  subroutine source_options(this)
    ! -- modules
    use KindModule, only: LGP
    use ConstantsModule, only: LENVARNAME
    use SimVariablesModule, only: errmsg
    use SimModule, only: store_error, store_error_filename
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy
    class(TspAdvType) :: this
    ! -- locals
    character(len=LENVARNAME), dimension(4) :: scheme = &
      &[character(len=LENVARNAME) :: 'UPSTREAM', 'CENTRAL', 'TVD', 'UTVD']
    logical(LGP) :: found_scheme, found_atspercel
    ! -- formats
    character(len=*), parameter :: fmtiadvwt = &
      &"(4x,'ADVECTION WEIGHTING SCHEME HAS BEEN SET TO: ', a)"

    ! update defaults with input sourced values
    call mem_set_value(this%iadvwt, 'SCHEME', this%input_mempath, &
                       scheme, found_scheme)
    call mem_set_value(this%ats_percel, 'ATS_PERCEL', this%input_mempath, &
                       found_atspercel)

    if (found_scheme) then
      ! should currently be set to index of scheme names
      if (this%iadvwt == 0) then
        write (errmsg, '(a, a)') &
          'Unknown scheme, must be "UPSTREAM", "CENTRAL", "TVD" or "UTVD"'
        call store_error(errmsg)
        call store_error_filename(this%input_fname)
      else
        ! scheme parameters are 0 based
        this%iadvwt = this%iadvwt - 1
      end if
    end if

    if (this%iadvwt == ADV_SCHEME_UTVD) then
      call developmode('UTVD is still under development, install the &
          &nightly build or compile from source with IDEVELOPMODE = 1.')
    end if

    if (found_atspercel) then
      if (this%ats_percel == DZERO) this%ats_percel = DNODATA
    end if

    ! log options
    write (this%iout, '(1x,a)') 'PROCESSING ADVECTION OPTIONS'
    if (found_scheme) then
      write (this%iout, fmtiadvwt) trim(scheme(this%iadvwt + 1))
    end if
    if (found_atspercel) then
      write (this%iout, '(4x,a,1pg15.6)') &
        'User-specified fractional cell distance for adaptive time &
        &steps: ', this%ats_percel
    end if
    write (this%iout, '(1x,a)') 'END OF ADVECTION OPTIONS'
  end subroutine source_options

end module TspAdvModule
