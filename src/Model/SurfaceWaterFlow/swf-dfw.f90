!> @brief Stream Network Flow (SWF) Diffusive Wave (DFW) Module
!!
!! This module solves one-dimensional flow routing using a diffusive
!< wave approach.
module SwfDfwModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENMEMPATH, LINELENGTH, &
                             DZERO, DHALF, DONE, DTWO, &
                             DTWOTHIRDS, DP9, DONETHIRD, &
                             DPREC, DEM10
  use MemoryHelperModule, only: create_mem_path
  use MemoryManagerModule, only: mem_allocate, mem_setptr, get_isize, &
                                 mem_reallocate
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, store_error_unit, &
                       store_error_filename
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use SwfCxsModule, only: SwfCxsType
  use ObsModule, only: ObsType, obs_cr
  use ObserveModule, only: ObserveType
  use MatrixBaseModule

  implicit none
  private
  public :: SwfDfwType, dfw_cr

  type, extends(NumericalPackageType) :: SwfDfwType

    ! user-provided input
    integer(I4B), pointer :: is2d => null() !< flag to indicate this model is 2D overland flow and not 1d channel flow
    integer(I4B), pointer :: icentral => null() !< flag to use central in space weighting (default is upstream weighting)
    integer(I4B), pointer :: iswrcond => null() !< flag to activate the dev SWR conductance formulation
    real(DP), pointer :: unitconv !< conversion factor used in mannings equation; calculated from timeconv and lengthconv
    real(DP), pointer :: timeconv !< conversion factor from model length units to meters (1.0 if model uses meters for length)
    real(DP), pointer :: lengthconv !< conversion factor from model time units to seconds (1.0 if model uses seconds for time)
    real(DP), dimension(:), pointer, contiguous :: hnew => null() !< pointer to model xnew
    real(DP), dimension(:), pointer, contiguous :: manningsn => null() !< mannings roughness for each reach
    integer(I4B), dimension(:), pointer, contiguous :: idcxs => null() !< cross section id for each reach
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    integer(I4B), dimension(:), pointer, contiguous :: icelltype => null() !< set to 1 and is accessed by chd for checking

    ! velocity
    integer(I4B), pointer :: icalcvelocity => null() !< flag to indicate velocity will be calculated (always on)
    integer(I4B), pointer :: isavvelocity => null() !< flag to indicate velocity will be saved
    real(DP), dimension(:, :), pointer, contiguous :: vcomp => null() !< velocity components: vx, vy, vz (nodes, 3)
    real(DP), dimension(:), pointer, contiguous :: vmag => null() !< velocity magnitude (of size nodes)
    integer(I4B), pointer :: nedges => null() !< number of cell edges
    integer(I4B), pointer :: lastedge => null() !< last edge number
    integer(I4B), dimension(:), pointer, contiguous :: nodedge => null() !< array of node numbers that have edges
    integer(I4B), dimension(:), pointer, contiguous :: ihcedge => null() !< edge type (horizontal or vertical)
    real(DP), dimension(:, :), pointer, contiguous :: propsedge => null() !< edge properties (Q, area, nx, ny, distance)
    real(DP), dimension(:), pointer, contiguous :: grad_dhds_mag => null() !< magnitude of the gradient (of size nodes)
    real(DP), dimension(:), pointer, contiguous :: dhdsja => null() !< gradient for each connection (of size njas)

    ! observation data
    integer(I4B), pointer :: inobspkg => null() !< unit number for obs package
    type(ObsType), pointer :: obs => null() !< observation package

    ! pointer to cross section data
    type(SwfCxsType), pointer :: cxs

  contains

    procedure :: dfw_df
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: dfw_load
    procedure :: source_options
    procedure :: log_options
    procedure :: source_griddata
    procedure :: log_griddata
    procedure :: dfw_ar
    procedure :: dfw_rp
    procedure :: dfw_ad
    procedure :: dfw_fc
    procedure :: dfw_qnm_fc_nr
    !procedure :: dfw_qnm_fc
    procedure :: dfw_fn
    procedure :: dfw_nur
    procedure :: dfw_cq
    procedure :: dfw_bd
    procedure :: dfw_save_model_flows
    procedure :: dfw_print_model_flows
    procedure :: dfw_da
    procedure :: dfw_df_obs
    procedure :: dfw_rp_obs
    procedure :: dfw_bd_obs
    procedure :: qcalc
    procedure :: get_cond
    procedure :: get_cond_swr
    procedure :: get_cond_n
    procedure :: get_flow_area_nm
    procedure :: calc_velocity
    procedure :: sav_velocity
    procedure, public :: increase_edge_count
    procedure, public :: set_edge_properties
    procedure :: calc_dhds
    procedure :: write_cxs_tables

  end type SwfDfwType

contains

  !> @brief create package
  !<
  subroutine dfw_cr(dfwobj, name_model, input_mempath, inunit, iout, &
                    cxs)
    ! modules
    use MemoryManagerExtModule, only: mem_set_value
    ! dummy
    type(SwfDfwType), pointer :: dfwobj !< object to create
    character(len=*), intent(in) :: name_model !< name of the SWF model
    character(len=*), intent(in) :: input_mempath !< memory path
    integer(I4B), intent(in) :: inunit !< flag to indicate if package is active
    integer(I4B), intent(in) :: iout !< unit number for output
    type(SwfCxsType), pointer, intent(in) :: cxs !< the pointer to the cxs package
    ! locals
    logical(LGP) :: found_fname
    ! formats
    character(len=*), parameter :: fmtheader = &
      "(1x, /1x, 'DFW --  DIFFUSIVE WAVE (DFW) PACKAGE, VERSION 1, 9/25/2023', &
       &' INPUT READ FROM MEMPATH: ', A, /)"
    !
    ! Create the object
    allocate (dfwobj)

    ! create name and memory path
    call dfwobj%set_names(1, name_model, 'DFW', 'DFW')

    ! Allocate scalars
    call dfwobj%allocate_scalars()

    ! Set variables
    dfwobj%input_mempath = input_mempath
    dfwobj%inunit = inunit
    dfwobj%iout = iout

    ! set name of input file
    call mem_set_value(dfwobj%input_fname, 'INPUT_FNAME', dfwobj%input_mempath, &
                       found_fname)

    ! Set a pointers to passed in objects
    dfwobj%cxs => cxs

    ! create obs package
    call obs_cr(dfwobj%obs, dfwobj%inobspkg)

    ! check if dfw is enabled
    if (inunit > 0) then

      ! Print a message identifying the package.
      write (iout, fmtheader) input_mempath

    end if

  end subroutine dfw_cr

  !> @brief load data from IDM to package
  !<
  subroutine dfw_df(this, dis)
    ! dummy
    class(SwfDfwType) :: this !< this instance
    class(DisBaseType), pointer, intent(inout) :: dis !< the pointer to the discretization

    ! Set a pointers to passed in objects
    this%dis => dis

    ! Set the distype (either DISV1D or DIS2D)
    if (this%dis%is_2d()) then
      this%is2d = 1
    end if

    ! check if dfw is enabled
    ! this will need to become if (.not. present(dfw_options)) then
    !if (inunit > 0) then

    ! allocate arrays
    call this%allocate_arrays()

    ! load dfw
    call this%dfw_load()

    !end if

  end subroutine dfw_df

  !> @ brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the package. The base model
  !! allocate scalars method is also called.
  !!
  !<
  subroutine allocate_scalars(this)
    ! modules
    ! dummy
    class(SwfDfwtype) :: this !< this instance
    !
    ! allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! Allocate scalars
    call mem_allocate(this%is2d, 'IS2D', this%memoryPath)
    call mem_allocate(this%icentral, 'ICENTRAL', this%memoryPath)
    call mem_allocate(this%iswrcond, 'ISWRCOND', this%memoryPath)
    call mem_allocate(this%unitconv, 'UNITCONV', this%memoryPath)
    call mem_allocate(this%lengthconv, 'LENGTHCONV', this%memoryPath)
    call mem_allocate(this%timeconv, 'TIMECONV', this%memoryPath)
    call mem_allocate(this%inobspkg, 'INOBSPKG', this%memoryPath)
    call mem_allocate(this%icalcvelocity, 'ICALCVELOCITY', this%memoryPath)
    call mem_allocate(this%isavvelocity, 'ISAVVELOCITY', this%memoryPath)
    call mem_allocate(this%nedges, 'NEDGES', this%memoryPath)
    call mem_allocate(this%lastedge, 'LASTEDGE', this%memoryPath)

    this%is2d = 0
    this%icentral = 0
    this%iswrcond = 0
    this%unitconv = DONE
    this%lengthconv = DONE
    this%timeconv = DONE
    this%inobspkg = 0
    this%icalcvelocity = 0
    this%isavvelocity = 0
    this%nedges = 0
    this%lastedge = 0

  end subroutine allocate_scalars

  !> @brief allocate memory for arrays
  !<
  subroutine allocate_arrays(this)
    ! dummy
    class(SwfDfwType) :: this !< this instance
    ! locals
    integer(I4B) :: n
    !
    ! user-provided input
    call mem_allocate(this%manningsn, this%dis%nodes, &
                      'MANNINGSN', this%memoryPath)
    call mem_allocate(this%idcxs, this%dis%nodes, &
                      'IDCXS', this%memoryPath)
    call mem_allocate(this%icelltype, this%dis%nodes, &
                      'ICELLTYPE', this%memoryPath)

    ! optional arrays
    call mem_allocate(this%nodedge, 0, 'NODEDGE', this%memoryPath)
    call mem_allocate(this%ihcedge, 0, 'IHCEDGE', this%memoryPath)
    call mem_allocate(this%propsedge, 0, 0, 'PROPSEDGE', this%memoryPath)

    ! Specific discharge is (re-)allocated when nedges is known
    call mem_allocate(this%vcomp, 3, 0, 'VCOMP', this%memoryPath)
    call mem_allocate(this%vmag, 0, 'VMAG', this%memoryPath)

    do n = 1, this%dis%nodes
      this%manningsn(n) = DZERO
      this%idcxs(n) = 0
      this%icelltype(n) = 1
    end do

    ! for 2d models, need to calculate and store dhds magnitude
    if (this%is2d == 1) then
      call mem_allocate(this%grad_dhds_mag, this%dis%nodes, &
                        'GRAD_DHDS_MAG', this%memoryPath)
      call mem_allocate(this%dhdsja, this%dis%njas, &
                        'DHDSJA', this%memoryPath)
      do n = 1, this%dis%nodes
        this%grad_dhds_mag(n) = DZERO
      end do
      do n = 1, this%dis%njas
        this%dhdsja(n) = DZERO
      end do
    end if

  end subroutine allocate_arrays

  !> @brief load data from IDM to package
  !<
  subroutine dfw_load(this)
    ! dummy
    class(SwfDfwType) :: this !< this instance

    ! source input data
    call this%source_options()
    call this%source_griddata()

  end subroutine dfw_load

  !> @brief Copy options from IDM into package
  !<
  subroutine source_options(this)
    ! modules
    use KindModule, only: LGP
    use InputOutputModule, only: getunit, openfile
    use MemoryManagerExtModule, only: mem_set_value
    use CharacterStringModule, only: CharacterStringType
    use SwfDfwInputModule, only: SwfDfwParamFoundType
    ! dummy
    class(SwfDfwType) :: this !< this instance
    ! locals
    integer(I4B) :: isize
    type(SwfDfwParamFoundType) :: found
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: obs6_fnames

    ! update defaults with idm sourced values
    call mem_set_value(this%icentral, 'ICENTRAL', &
                       this%input_mempath, found%icentral)
    call mem_set_value(this%iswrcond, 'ISWRCOND', &
                       this%input_mempath, found%iswrcond)
    call mem_set_value(this%lengthconv, 'LENGTHCONV', &
                       this%input_mempath, found%lengthconv)
    call mem_set_value(this%timeconv, 'TIMECONV', &
                       this%input_mempath, found%timeconv)
    call mem_set_value(this%iprflow, 'IPRFLOW', &
                       this%input_mempath, found%iprflow)
    call mem_set_value(this%ipakcb, 'IPAKCB', &
                       this%input_mempath, found%ipakcb)
    call mem_set_value(this%isavvelocity, 'ISAVVELOCITY', &
                       this%input_mempath, found%isavvelocity)

    ! save flows option active
    if (found%icentral) this%icentral = 1
    if (found%ipakcb) this%ipakcb = -1

    ! calculate unit conversion
    this%unitconv = this%lengthconv**DONETHIRD
    this%unitconv = this%unitconv * this%timeconv

    ! save velocity active
    if (found%isavvelocity) this%icalcvelocity = this%isavvelocity

    ! check for obs6_filename
    call get_isize('OBS6_FILENAME', this%input_mempath, isize)
    if (isize > 0) then
      !
      if (isize /= 1) then
        errmsg = 'Multiple OBS6 keywords detected in OPTIONS block.'// &
                 ' Only one OBS6 entry allowed.'
        call store_error(errmsg)
        call store_error_filename(this%input_fname)
      end if

      call mem_setptr(obs6_fnames, 'OBS6_FILENAME', this%input_mempath)

      found%obs6_filename = .true.
      this%obs%inputFilename = obs6_fnames(1)
      this%obs%active = .true.
      this%inobspkg = GetUnit()
      this%obs%inUnitObs = this%inobspkg
      call openfile(this%inobspkg, this%iout, this%obs%inputFilename, 'OBS')
      call this%obs%obs_df(this%iout, this%packName, this%filtyp, this%dis)
      call this%dfw_df_obs()
    end if

    ! log values to list file
    if (this%iout > 0) then
      call this%log_options(found)
    end if

  end subroutine source_options

  !> @brief Write user options to list file
  !<
  subroutine log_options(this, found)
    use SwfDfwInputModule, only: SwfDfwParamFoundType
    class(SwfDfwType) :: this !< this instance
    type(SwfDfwParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'Setting DFW Options'

    if (found%lengthconv) then
      write (this%iout, '(4x,a, G0)') 'Mannings length conversion value &
                                  &specified as ', this%lengthconv
    end if

    if (found%timeconv) then
      write (this%iout, '(4x,a, G0)') 'Mannings time conversion value &
                                  &specified as ', this%timeconv
    end if

    if (found%lengthconv .or. found%timeconv) then
      write (this%iout, '(4x,a, G0)') 'Mannings conversion value calculated &
                                  &from user-provided length_conversion and &
                                  &time_conversion is ', this%unitconv
    end if

    if (found%iprflow) then
      write (this%iout, '(4x,a)') 'Cell-by-cell flow information will be printed &
                                  &to listing file whenever ICBCFL is not zero.'
    end if

    if (found%ipakcb) then
      write (this%iout, '(4x,a)') 'Cell-by-cell flow information will be printed &
                                  &to listing file whenever ICBCFL is not zero.'
    end if

    if (found%obs6_filename) then
      write (this%iout, '(4x,a)') 'Observation package is active.'
    end if

    if (found%isavvelocity) &
      write (this%iout, '(4x,a)') 'Velocity will be calculated at cell &
                                  &centers and written to DATA-VCOMP in budget &
                                  &file when requested.'

    if (found%iswrcond) then
      write (this%iout, '(4x,a, G0)') 'Conductance will be calculated using &
                                       &the SWR development option.'
    end if

    write (this%iout, '(1x,a,/)') 'End Setting DFW Options'

  end subroutine log_options

  !> @brief copy griddata from IDM to package
  !<
  subroutine source_griddata(this)
    ! modules
    use SimModule, only: count_errors, store_error
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerExtModule, only: mem_set_value
    use SimVariablesModule, only: idm_context
    use SwfDfwInputModule, only: SwfDfwParamFoundType
    ! dummy
    class(SwfDfwType) :: this !< this instance
    ! locals
    character(len=LENMEMPATH) :: idmMemoryPath
    type(SwfDfwParamFoundType) :: found
    integer(I4B), dimension(:), pointer, contiguous :: map

    ! set memory path
    idmMemoryPath = create_mem_path(this%name_model, 'DFW', idm_context)

    ! set map to convert user input data into reduced data
    map => null()
    if (this%dis%nodes < this%dis%nodesuser) map => this%dis%nodeuser

    ! update defaults with idm sourced values
    call mem_set_value(this%manningsn, 'MANNINGSN', &
                       idmMemoryPath, map, found%manningsn)
    call mem_set_value(this%idcxs, 'IDCXS', idmMemoryPath, map, found%idcxs)

    ! ensure MANNINGSN was found
    if (.not. found%manningsn) then
      write (errmsg, '(a)') 'Error in GRIDDATA block: MANNINGSN not found.'
      call store_error(errmsg)
    end if

    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if

    ! log griddata
    if (this%iout > 0) then
      call this%log_griddata(found)
    end if

  end subroutine source_griddata

  !> @brief log griddata to list file
  !<
  subroutine log_griddata(this, found)
    use SwfDfwInputModule, only: SwfDfwParamFoundType
    class(SwfDfwType) :: this !< this instance
    type(SwfDfwParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'Setting DFW Griddata'

    if (found%manningsn) then
      write (this%iout, '(4x,a)') 'MANNINGSN set from input file'
    end if

    if (found%idcxs) then
      write (this%iout, '(4x,a)') 'IDCXS set from input file'
    end if

    call this%write_cxs_tables()

    write (this%iout, '(1x,a,/)') 'End Setting DFW Griddata'

  end subroutine log_griddata

  subroutine write_cxs_tables(this)
    ! modules
    ! dummy
    class(SwfDfwType) :: this !< this instance
    ! local
    ! integer(I4B) :: idcxs
    ! integer(I4B) :: n

    !-- TODO: write cross section tables
    ! do n = 1, this%dis%nodes
    !   idcxs = this%idcxs(n)
    !   if (idcxs > 0) then
    !     call this%cxs%write_cxs_table(idcxs, this%width(n), this%slope(n), &
    !                                   this%manningsn(n), this%unitconv)
    !   end if
    ! end do
  end subroutine write_cxs_tables

  !> @brief allocate memory
  !<
  subroutine dfw_ar(this, ibound, hnew)
    ! modules
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B), dimension(:), pointer, contiguous :: ibound !< model ibound array
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: hnew !< pointer to model head array
    ! local
    integer(I4B) :: n

    ! store pointer to ibound
    this%ibound => ibound
    this%hnew => hnew

    if (this%icalcvelocity == 1) then
      call mem_reallocate(this%vcomp, 3, this%dis%nodes, 'VCOMP', this%memoryPath)
      call mem_reallocate(this%vmag, this%dis%nodes, 'VMAG', this%memoryPath)
      call mem_reallocate(this%nodedge, this%nedges, 'NODEDGE', this%memoryPath)
      call mem_reallocate(this%ihcedge, this%nedges, 'IHCEDGE', this%memoryPath)
      call mem_reallocate(this%propsedge, 5, this%nedges, 'PROPSEDGE', &
                          this%memoryPath)
      do n = 1, this%dis%nodes
        this%vcomp(:, n) = DZERO
        this%vmag(n) = DZERO
      end do
    end if

    ! observation data
    call this%obs%obs_ar()

  end subroutine dfw_ar

  !> @brief allocate memory
  !<
  subroutine dfw_rp(this)
    ! modules
    ! dummy
    class(SwfDfwType) :: this !< this instance

    ! read observations
    call this%dfw_rp_obs()

  end subroutine dfw_rp

  !> @brief advance
  !<
  subroutine dfw_ad(this, irestore)
    class(SwfDfwType) :: this !< this instance
    integer(I4B), intent(in) :: irestore !< ATS flag for retrying time step (1) or advancing (0)

    ! Push simulated values to preceding time/subtime step
    call this%obs%obs_ad()

  end subroutine dfw_ad

  !> @brief fill coefficients
  !!
  !! The DFW Package is entirely Newton based.  All matrix and rhs terms
  !! are added from thish routine.
  !!
  !<
  subroutine dfw_fc(this, kiter, matrix_sln, idxglo, rhs, stage, stage_old)
    ! modules
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    real(DP), intent(inout), dimension(:) :: stage
    real(DP), intent(inout), dimension(:) :: stage_old
    ! local

    ! calculate dhds at cell center for 2d case
    if (this%is2d == 1) then
      call this%calc_dhds()
    end if

    ! add qnm contributions to matrix equations
    call this%dfw_qnm_fc_nr(kiter, matrix_sln, idxglo, rhs, stage, stage_old)

  end subroutine dfw_fc

  !> @brief fill coefficients
  !!
  !< Add qnm contributions to matrix equations
  subroutine dfw_qnm_fc_nr(this, kiter, matrix_sln, idxglo, rhs, stage, stage_old)
    ! modules
    use MathUtilModule, only: get_perturbation
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    real(DP), intent(inout), dimension(:) :: stage
    real(DP), intent(inout), dimension(:) :: stage_old
    ! local
    integer(I4B) :: n, m, ii, idiag
    real(DP) :: qnm
    real(DP) :: qeps
    real(DP) :: eps
    real(DP) :: derv

    ! Calculate conductance and put into amat
    do n = 1, this%dis%nodes

      ! Find diagonal position for row n
      idiag = this%dis%con%ia(n)

      ! Loop through connections adding matrix terms
      do ii = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1

        ! skip for masked cells
        if (this%dis%con%mask(ii) == 0) cycle

        ! connection variables
        m = this%dis%con%ja(ii)

        ! Fill the qnm term on the right-hand side
        qnm = this%qcalc(n, m, stage(n), stage(m), ii)
        rhs(n) = rhs(n) - qnm

        ! Derivative calculation and fill of n terms
        eps = get_perturbation(stage(n))
        qeps = this%qcalc(n, m, stage(n) + eps, stage(m), ii)
        derv = (qeps - qnm) / eps
        call matrix_sln%add_value_pos(idxglo(idiag), derv)
        rhs(n) = rhs(n) + derv * stage(n)

        ! Derivative calculation and fill of m terms
        eps = get_perturbation(stage(m))
        qeps = this%qcalc(n, m, stage(n), stage(m) + eps, ii)
        derv = (qeps - qnm) / eps
        call matrix_sln%add_value_pos(idxglo(ii), derv)
        rhs(n) = rhs(n) + derv * stage(m)

      end do
    end do

  end subroutine dfw_qnm_fc_nr

  !> @brief fill newton
  !<
  subroutine dfw_fn(this, kiter, matrix_sln, idxglo, rhs, stage)
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    real(DP), intent(inout), dimension(:) :: stage
    ! local

    ! add newton terms to solution matrix
    ! todo: add newton terms here instead?
    ! this routine is probably not necessary as method is fully newton

  end subroutine dfw_fn

  !> @brief calculate flow between cells n and m
  !<
  function qcalc(this, n, m, stage_n, stage_m, ipos) result(qnm)
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B), intent(in) :: n !< number for cell n
    integer(I4B), intent(in) :: m !< number for cell m
    real(DP), intent(in) :: stage_n !< stage in reach n
    real(DP), intent(in) :: stage_m !< stage in reach m
    integer(I4B), intent(in) :: ipos !< connection number
    ! local
    integer(I4B) :: isympos
    real(DP) :: qnm
    real(DP) :: cond
    real(DP) :: cl1
    real(DP) :: cl2

    ! Set connection lengths
    isympos = this%dis%con%jas(ipos)
    if (n < m) then
      cl1 = this%dis%con%cl1(isympos)
      cl2 = this%dis%con%cl2(isympos)
    else
      cl1 = this%dis%con%cl2(isympos)
      cl2 = this%dis%con%cl1(isympos)
    end if

    ! Calculate conductance
    if (this%iswrcond == 0) then
      cond = this%get_cond(n, m, ipos, stage_n, stage_m, cl1, cl2)
    else if (this%iswrcond == 1) then
      cond = this%get_cond_swr(n, m, ipos, stage_n, stage_m, cl1, cl2)
    end if

    ! calculate flow between n and m
    qnm = cond * (stage_m - stage_n)

  end function qcalc

  !> @brief calculate effective conductance between cells n and m
  !!
  !! Calculate half-cell conductances for cell n and cell m and then use
  !! harmonic averaging to calculate the effective conductance between the
  !< two cells.
  function get_cond(this, n, m, ipos, stage_n, stage_m, cln, clm) result(cond)
    ! modules
    use SmoothingModule, only: sQuadratic
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B), intent(in) :: n !< number for cell n
    integer(I4B), intent(in) :: m !< number for cell m
    integer(I4B), intent(in) :: ipos !< connection number
    real(DP), intent(in) :: stage_n !< stage in reach n
    real(DP), intent(in) :: stage_m !< stage in reach m
    real(DP), intent(in) :: cln !< distance from cell n to shared face with m
    real(DP), intent(in) :: clm !< distance from cell m to shared face with n
    ! local
    real(DP) :: depth_n
    real(DP) :: depth_m
    real(DP) :: dhds_n
    real(DP) :: dhds_m
    real(DP) :: width_n
    real(DP) :: width_m
    real(DP) :: range = 1.d-6
    real(DP) :: dydx
    real(DP) :: smooth_factor
    real(DP) :: length_nm
    real(DP) :: cond
    real(DP) :: cn
    real(DP) :: cm

    ! we are using a harmonic conductance approach here; however
    ! the SWR Process for MODFLOW-2005/NWT uses length-weighted
    ! average areas and hydraulic radius instead.
    length_nm = cln + clm
    cond = DZERO
    if (length_nm > DPREC) then

      ! Calculate depth in each reach
      depth_n = stage_n - this%dis%bot(n)
      depth_m = stage_m - this%dis%bot(m)

      ! assign gradients
      if (this%is2d == 0) then
        dhds_n = abs(stage_m - stage_n) / (cln + clm)
        dhds_m = dhds_n
      else
        dhds_n = this%grad_dhds_mag(n)
        dhds_m = this%grad_dhds_mag(m)
      end if

      ! Assign upstream depth, if not central
      if (this%icentral == 0) then
        ! use upstream weighting
        if (stage_n > stage_m) then
          depth_m = depth_n
        else
          depth_n = depth_m
        end if
      end if

      ! Calculate a smoothed depth that goes to zero over
      ! the specified range
      call sQuadratic(depth_n, range, dydx, smooth_factor)
      depth_n = depth_n * smooth_factor
      call sQuadratic(depth_m, range, dydx, smooth_factor)
      depth_m = depth_m * smooth_factor

      ! Get the flow widths for n and m from dis package
      call this%dis%get_flow_width(n, m, ipos, width_n, width_m)

      ! Calculate half-cell conductance for reach
      ! n and m
      cn = this%get_cond_n(n, depth_n, cln, width_n, dhds_n)
      cm = this%get_cond_n(m, depth_m, clm, width_m, dhds_m)

      ! Use harmonic mean to calculate weighted
      ! conductance between the centers of reaches
      ! n and m
      if ((cn + cm) > DPREC) then
        cond = cn * cm / (cn + cm)
      else
        cond = DZERO
      end if

    end if

  end function get_cond

  !> @brief Calculate half cell conductance
  !!
  !! Calculate half-cell conductance for cell n
  !< using conveyance and Manning's equation
  function get_cond_n(this, n, depth, dx, width, dhds) result(c)
    ! modules
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B), intent(in) :: n !< reach number
    real(DP), intent(in) :: depth !< simulated depth (stage - elevation) in reach n for this iteration
    real(DP), intent(in) :: dx !< half-cell distance
    real(DP), intent(in) :: width !< width of the reach perpendicular to flow
    real(DP), intent(in) :: dhds !< gradient
    ! return
    real(DP) :: c
    ! local
    real(DP) :: rough
    real(DP) :: dhds_sqr
    real(DP) :: conveyance

    ! Calculate conveyance, which is a * r**DTWOTHIRDS / roughc
    rough = this%manningsn(n)
    conveyance = this%cxs%get_conveyance(this%idcxs(n), width, depth, rough)
    dhds_sqr = dhds**DHALF
    if (dhds_sqr < DEM10) then
      dhds_sqr = DEM10
    end if

    ! Multiply by unitconv and divide conveyance by sqrt of friction slope and dx
    c = this%unitconv * conveyance / dx / dhds_sqr

  end function get_cond_n

  !> @brief Calculate effective conductance for cells n and m using SWR method
  !!
  !! The SWR Process for MODFLOW uses average cell parameters from cell n and
  !! m to calculate an effective conductance.  This is different from the
  !! default approach used in SWF, which uses harmonic averaging on two half-
  !< cell conductances.
  function get_cond_swr(this, n, m, ipos, stage_n, stage_m, cln, clm) result(cond)
    ! modules
    use SmoothingModule, only: sQuadratic
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B), intent(in) :: n !< number for cell n
    integer(I4B), intent(in) :: m !< number for cell m
    integer(I4B), intent(in) :: ipos !< connection number
    real(DP), intent(in) :: stage_n !< stage in reach n
    real(DP), intent(in) :: stage_m !< stage in reach m
    real(DP), intent(in) :: cln !< distance from cell n to shared face with m
    real(DP), intent(in) :: clm !< distance from cell m to shared face with n
    ! local
    real(DP) :: depth_n
    real(DP) :: depth_m
    real(DP) :: dhds_n
    real(DP) :: dhds_m
    real(DP) :: dhds_nm
    real(DP) :: dhds_sqr
    real(DP) :: width_n
    real(DP) :: width_m
    real(DP) :: range = 1.d-6
    real(DP) :: dydx
    real(DP) :: smooth_factor
    real(DP) :: length_nm
    real(DP) :: cond
    real(DP) :: ravg
    real(DP) :: rinv_avg
    real(DP) :: area_n, area_m, area_avg
    real(DP) :: rhn, rhm, rhavg
    real(DP) :: weight_n
    real(DP) :: weight_m
    real(DP) :: rough_n
    real(DP) :: rough_m

    ! Use harmonic weighting for 1/manningsn, but using length-weighted
    ! averaging for other terms
    length_nm = cln + clm
    cond = DZERO
    if (length_nm > DPREC) then

      ! Calculate depth in each reach
      depth_n = stage_n - this%dis%bot(n)
      depth_m = stage_m - this%dis%bot(m)

      ! Assign upstream depth, if not central
      if (this%icentral == 0) then
        ! use upstream weighting
        if (stage_n > stage_m) then
          depth_m = depth_n
        else
          depth_n = depth_m
        end if
      end if

      ! Calculate a smoothed depth that goes to zero over
      !    the specified range
      call sQuadratic(depth_n, range, dydx, smooth_factor)
      depth_n = depth_n * smooth_factor
      call sQuadratic(depth_m, range, dydx, smooth_factor)
      depth_m = depth_m * smooth_factor

      ! Get the flow widths for n and m from dis package
      call this%dis%get_flow_width(n, m, ipos, width_n, width_m)

      ! linear weight toward node closer to shared face
      weight_n = clm / length_nm
      weight_m = DONE - weight_n

      ! average cross sectional flow area
      area_n = this%cxs%get_area(this%idcxs(n), width_n, depth_n)
      area_m = this%cxs%get_area(this%idcxs(m), width_m, depth_m)
      area_avg = weight_n * area_n + weight_m * area_m

      ! average hydraulic radius
      if (this%is2d == 0) then
        rhn = this%cxs%get_hydraulic_radius(this%idcxs(n), width_n, &
                                            depth_n, area_n)
        rhm = this%cxs%get_hydraulic_radius(this%idcxs(m), width_m, &
                                            depth_m, area_m)
        rhavg = weight_n * rhn + weight_m * rhm
      else
        rhavg = area_avg / width_n
      end if
      rhavg = rhavg**DTWOTHIRDS

      ! average gradient
      if (this%is2d == 0) then
        dhds_nm = abs(stage_m - stage_n) / (length_nm)
      else
        dhds_n = this%grad_dhds_mag(n)
        dhds_m = this%grad_dhds_mag(m)
        dhds_nm = weight_n * dhds_n + weight_m * dhds_m
      end if
      dhds_sqr = dhds_nm**DHALF
      if (dhds_sqr < DEM10) then
        dhds_sqr = DEM10
      end if

      ! weighted harmonic mean for inverse mannings value
      weight_n = cln / length_nm
      weight_m = DONE - weight_n
      rough_n = this%cxs%get_roughness(this%idcxs(n), width_n, depth_n, &
                                       this%manningsn(n), dhds_nm)
      rough_m = this%cxs%get_roughness(this%idcxs(m), width_m, depth_m, &
                                       this%manningsn(m), dhds_nm)
      ravg = (weight_n + weight_m) / &
             (weight_n / rough_n + weight_m / rough_m)
      rinv_avg = DONE / ravg

      ! calculate conductance using averaged values
      cond = this%unitconv * rinv_avg * area_avg * rhavg / dhds_sqr / length_nm

    end if

  end function get_cond_swr

  !> @brief Calculate flow area between cell n and m
  !!
  !! Calculate an average flow area between cell n and m.
  !! First calculate a flow area for cell n and then for
  !! cell m and linearly weight the areas using the connection
  !< distances.
  function get_flow_area_nm(this, n, m, stage_n, stage_m, cln, clm, &
                            ipos) result(area_avg)
    ! module
    use SmoothingModule, only: sQuadratic
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    real(DP), intent(in) :: stage_n
    real(DP), intent(in) :: stage_m
    real(DP), intent(in) :: cln
    real(DP), intent(in) :: clm
    integer(I4B), intent(in) :: ipos
    ! local
    real(DP) :: depth_n
    real(DP) :: depth_m
    real(DP) :: width_n
    real(DP) :: width_m
    real(DP) :: area_n
    real(DP) :: area_m
    real(DP) :: weight_n
    real(DP) :: weight_m
    real(DP) :: length_nm
    real(DP) :: range = 1.d-6
    real(DP) :: dydx
    real(DP) :: smooth_factor
    ! return
    real(DP) :: area_avg

    ! depths
    depth_n = stage_n - this%dis%bot(n)
    depth_m = stage_m - this%dis%bot(m)

    ! Assign upstream depth, if not central
    if (this%icentral == 0) then
      ! use upstream weighting
      if (stage_n > stage_m) then
        depth_m = depth_n
      else
        depth_n = depth_m
      end if
    end if

    ! Calculate a smoothed depth that goes to zero over
    ! the specified range
    call sQuadratic(depth_n, range, dydx, smooth_factor)
    depth_n = depth_n * smooth_factor
    call sQuadratic(depth_m, range, dydx, smooth_factor)
    depth_m = depth_m * smooth_factor

    ! Get the flow widths for n and m from dis package
    call this%dis%get_flow_width(n, m, ipos, width_n, width_m)

    ! linear weight toward node closer to shared face
    length_nm = cln + clm
    weight_n = clm / length_nm
    weight_m = DONE - weight_n

    ! average cross sectional flow area
    area_n = this%cxs%get_area(this%idcxs(n), width_n, depth_n)
    area_m = this%cxs%get_area(this%idcxs(m), width_m, depth_m)
    area_avg = weight_n * area_n + weight_m * area_m

  end function get_flow_area_nm

  !> @brief Calculate average hydraulic gradient magnitude for each cell
  !!
  !! Go through each cell and calculate the average hydraulic gradient using
  !! an xt3d-style gradient interpolation.  This is used for 2D grids in the
  !! calculation of an effective conductance.  For 1D grids, gradients are
  !< calculated between cell centers.
  subroutine calc_dhds(this)
    ! modules
    use VectorInterpolationModule, only: vector_interpolation_2d
    ! dummy
    class(SwfDfwType) :: this !< this instance
    ! local
    integer(I4B) :: n
    integer(I4B) :: m
    integer(I4B) :: ipos
    integer(I4B) :: isympos
    real(DP) :: cl1
    real(DP) :: cl2

    do n = 1, this%dis%nodes
      this%grad_dhds_mag(n) = DZERO
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        isympos = this%dis%con%jas(ipos)

        ! determine cl1 and cl2
        if (n < m) then
          cl1 = this%dis%con%cl1(isympos)
          cl2 = this%dis%con%cl2(isympos)
        else
          cl1 = this%dis%con%cl2(isympos)
          cl2 = this%dis%con%cl1(isympos)
        end if

        ! store for n < m in upper right triangular part of symmetric dhdsja array
        if (n < m) then
          if (cl1 + cl2 > DPREC) then
            this%dhdsja(isympos) = (this%hnew(m) - this%hnew(n)) / (cl1 + cl2)
          else
            this%dhdsja(isympos) = DZERO
          end if
        end if
      end do
    end do

    ! pass dhdsja into the vector interpolation to get the components
    ! of the gradient at the cell center
    call vector_interpolation_2d(this%dis, this%dhdsja, vmag=this%grad_dhds_mag)

  end subroutine calc_dhds

  !> @ brief Newton under relaxation
  !!
  !<  If stage is below the bottom, then pull it up a bit
  subroutine dfw_nur(this, neqmod, x, xtemp, dx, inewtonur, dxmax, locmax)
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B), intent(in) :: neqmod !< number of equations
    real(DP), dimension(neqmod), intent(inout) :: x !< dependent variable
    real(DP), dimension(neqmod), intent(in) :: xtemp !< temporary dependent variable
    real(DP), dimension(neqmod), intent(inout) :: dx !< change in dependent variable
    integer(I4B), intent(inout) :: inewtonur !< flag to indication relaxation was applied
    real(DP), intent(inout) :: dxmax !< max change in x
    integer(I4B), intent(inout) :: locmax !< location of max change
    ! local
    integer(I4B) :: n
    real(DP) :: botm
    real(DP) :: xx
    real(DP) :: dxx

    ! Newton-Raphson under-relaxation
    do n = 1, this%dis%nodes
      if (this%ibound(n) < 1) cycle
      if (this%icelltype(n) > 0) then
        botm = this%dis%bot(n)
        ! only apply Newton-Raphson under-relaxation if
        ! solution head is below the bottom of the model
        if (x(n) < botm) then
          inewtonur = 1
          xx = xtemp(n) * (DONE - DP9) + botm * DP9
          dxx = x(n) - xx
          if (abs(dxx) > abs(dxmax)) then
            locmax = n
            dxmax = dxx
          end if
          x(n) = xx
          dx(n) = DZERO
        end if
      end if
    end do

  end subroutine dfw_nur

  !> @ brief Calculate flow for each connection and store in flowja
  !<
  subroutine dfw_cq(this, stage, stage_old, flowja)
    ! dummy
    class(SwfDfwType) :: this !< this instance
    real(DP), intent(inout), dimension(:) :: stage !< calculated head
    real(DP), intent(inout), dimension(:) :: stage_old !< calculated head from previous time step
    real(DP), intent(inout), dimension(:) :: flowja !< vector of flows in CSR format
    ! local
    integer(I4B) :: n, ipos, m
    real(DP) :: qnm

    do n = 1, this%dis%nodes
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        if (m < n) cycle
        qnm = this%qcalc(n, m, stage(n), stage(m), ipos)
        flowja(ipos) = qnm
        flowja(this%dis%con%isym(ipos)) = -qnm
      end do
    end do

  end subroutine dfw_cq

  !> @ brief Model budget calculation for package
  !<
  subroutine dfw_bd(this, isuppress_output, model_budget)
    ! modules
    use BudgetModule, only: BudgetType
    ! dummy variables
    class(SwfDfwType) :: this !< this instance
    integer(I4B), intent(in) :: isuppress_output !< flag to suppress model output
    type(BudgetType), intent(inout) :: model_budget !< model budget object

    ! Add any DFW budget terms
    ! none

  end subroutine dfw_bd

  !> @ brief save flows for package
  !<
  subroutine dfw_save_model_flows(this, flowja, icbcfl, icbcun)
    ! dummy
    class(SwfDfwType) :: this !< this instance
    real(DP), dimension(:), intent(in) :: flowja !< vector of flows in CSR format
    integer(I4B), intent(in) :: icbcfl !< flag to indicate if flows should be saved
    integer(I4B), intent(in) :: icbcun !< unit number for flow output
    ! local
    integer(I4B) :: ibinun

    ! Set unit number for binary output
    if (this%ipakcb < 0) then
      ibinun = icbcun
    elseif (this%ipakcb == 0) then
      ibinun = 0
    else
      ibinun = this%ipakcb
    end if
    if (icbcfl == 0) ibinun = 0

    ! Write the face flows if requested
    if (ibinun /= 0) then
      ! flowja
      call this%dis%record_connection_array(flowja, ibinun, this%iout)
    end if

    ! Calculate velocities at cell centers and write, if requested
    if (this%isavvelocity /= 0) then
      if (ibinun /= 0) call this%sav_velocity(ibinun)
    end if

  end subroutine dfw_save_model_flows

  !> @ brief print flows for package
  !<
  subroutine dfw_print_model_flows(this, ibudfl, flowja)
    ! modules
    use TdisModule, only: kper, kstp
    use ConstantsModule, only: LENBIGLINE
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B), intent(in) :: ibudfl !< print flag
    real(DP), intent(inout), dimension(:) :: flowja !< vector of flows in CSR format
    ! local
    character(len=LENBIGLINE) :: line
    character(len=30) :: tempstr
    integer(I4B) :: n, ipos, m
    real(DP) :: qnm
    ! formats
    character(len=*), parameter :: fmtiprflow = &
      &"(/,4x,'CALCULATED INTERCELL FLOW FOR PERIOD ', i0, ' STEP ', i0)"

    ! Write flowja to list file if requested
    if (ibudfl /= 0 .and. this%iprflow > 0) then
      write (this%iout, fmtiprflow) kper, kstp
      do n = 1, this%dis%nodes
        line = ''
        call this%dis%noder_to_string(n, tempstr)
        line = trim(tempstr)//':'
        do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
          m = this%dis%con%ja(ipos)
          call this%dis%noder_to_string(m, tempstr)
          line = trim(line)//' '//trim(tempstr)
          qnm = flowja(ipos)
          write (tempstr, '(1pg15.6)') qnm
          line = trim(line)//' '//trim(adjustl(tempstr))
        end do
        write (this%iout, '(a)') trim(line)
      end do
    end if

  end subroutine dfw_print_model_flows

  !> @brief deallocate memory
  !<
  subroutine dfw_da(this)
    ! modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    ! dummy
    class(SwfDfwType) :: this !< this instance

    ! Deallocate input memory
    call memorystore_remove(this%name_model, 'DFW', idm_context)

    ! Deallocate arrays
    call mem_deallocate(this%manningsn)
    call mem_deallocate(this%idcxs)
    call mem_deallocate(this%icelltype)
    call mem_deallocate(this%nodedge)
    call mem_deallocate(this%ihcedge)
    call mem_deallocate(this%propsedge)
    call mem_deallocate(this%vcomp)
    call mem_deallocate(this%vmag)
    if (this%is2d == 1) then
      call mem_deallocate(this%grad_dhds_mag)
      call mem_deallocate(this%dhdsja)
    end if

    ! Scalars
    call mem_deallocate(this%is2d)
    call mem_deallocate(this%icentral)
    call mem_deallocate(this%iswrcond)
    call mem_deallocate(this%unitconv)
    call mem_deallocate(this%lengthconv)
    call mem_deallocate(this%timeconv)
    call mem_deallocate(this%isavvelocity)
    call mem_deallocate(this%icalcvelocity)
    call mem_deallocate(this%nedges)
    call mem_deallocate(this%lastedge)

    ! obs package
    call mem_deallocate(this%inobspkg)
    call this%obs%obs_da()
    deallocate (this%obs)
    nullify (this%obs)
    nullify (this%cxs)

    ! deallocate parent
    call this%NumericalPackageType%da()

    ! pointers
    this%hnew => null()

  end subroutine dfw_da

  !> @brief Calculate the 3 components of velocity at the cell center
  !!
  !! todo: duplicated from NPF; should consolidate
  !<
  subroutine calc_velocity(this, flowja)
    ! modules
    ! dummy
    class(SwfDfwType) :: this !< this instance
    real(DP), intent(in), dimension(:) :: flowja !< vector of flows in CSR format
    ! local
    integer(I4B) :: n
    integer(I4B) :: m
    integer(I4B) :: ipos
    integer(I4B) :: isympos
    integer(I4B) :: ihc
    integer(I4B) :: ic
    integer(I4B) :: iz
    integer(I4B) :: nc
    integer(I4B) :: ncz
    real(DP) :: vx
    real(DP) :: vy
    real(DP) :: vz
    real(DP) :: xn
    real(DP) :: yn
    real(DP) :: zn
    real(DP) :: xc
    real(DP) :: yc
    real(DP) :: zc
    real(DP) :: cl1
    real(DP) :: cl2
    real(DP) :: dltot
    real(DP) :: ooclsum
    real(DP) :: dsumx
    real(DP) :: dsumy
    real(DP) :: dsumz
    real(DP) :: denom
    real(DP) :: area
    real(DP) :: axy
    real(DP) :: ayx
    real(DP), allocatable, dimension(:) :: vi
    real(DP), allocatable, dimension(:) :: di
    real(DP), allocatable, dimension(:) :: viz
    real(DP), allocatable, dimension(:) :: diz
    real(DP), allocatable, dimension(:) :: nix
    real(DP), allocatable, dimension(:) :: niy
    real(DP), allocatable, dimension(:) :: wix
    real(DP), allocatable, dimension(:) :: wiy
    real(DP), allocatable, dimension(:) :: wiz
    real(DP), allocatable, dimension(:) :: bix
    real(DP), allocatable, dimension(:) :: biy
    logical :: nozee = .true.

    ! Ensure dis has necessary information
    ! todo: do we need this for SWF?
    if (this%icalcvelocity /= 0 .and. this%dis%con%ianglex == 0) then
      call store_error('Error.  ANGLDEGX not provided in '// &
                       'discretization file.  ANGLDEGX required for '// &
                       'calculation of velocity.', terminate=.TRUE.)
    end if

    ! Find max number of connections and allocate weight arrays
    nc = 0
    do n = 1, this%dis%nodes

      ! Count internal model connections
      ic = this%dis%con%ia(n + 1) - this%dis%con%ia(n) - 1

      ! Count edge connections
      do m = 1, this%nedges
        if (this%nodedge(m) == n) then
          ic = ic + 1
        end if
      end do

      ! Set max number of connections for any cell
      if (ic > nc) nc = ic
    end do

    ! Allocate storage arrays needed for cell-centered calculation
    allocate (vi(nc))
    allocate (di(nc))
    allocate (viz(nc))
    allocate (diz(nc))
    allocate (nix(nc))
    allocate (niy(nc))
    allocate (wix(nc))
    allocate (wiy(nc))
    allocate (wiz(nc))
    allocate (bix(nc))
    allocate (biy(nc))

    ! Go through each cell and calculate specific discharge
    do n = 1, this%dis%nodes

      ! first calculate geometric properties for x and y directions and
      !    the specific discharge at a face (vi)
      ic = 0
      iz = 0
      vi(:) = DZERO
      di(:) = DZERO
      viz(:) = DZERO
      diz(:) = DZERO
      nix(:) = DZERO
      niy(:) = DZERO
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        isympos = this%dis%con%jas(ipos)
        ihc = this%dis%con%ihc(isympos)
        ic = ic + 1
        call this%dis%connection_normal(n, m, ihc, xn, yn, zn, ipos)
        call this%dis%connection_vector(n, m, nozee, DONE, DONE, &
                                        ihc, xc, yc, zc, dltot)
        cl1 = this%dis%con%cl1(isympos)
        cl2 = this%dis%con%cl2(isympos)
        if (m < n) then
          cl1 = this%dis%con%cl2(isympos)
          cl2 = this%dis%con%cl1(isympos)
        end if
        ooclsum = DONE / (cl1 + cl2)
        nix(ic) = -xn
        niy(ic) = -yn
        di(ic) = dltot * cl1 * ooclsum
        area = this%get_flow_area_nm(n, m, this%hnew(n), this%hnew(m), &
                                     cl1, cl2, ipos)
        if (area > DZERO) then
          vi(ic) = flowja(ipos) / area
        else
          vi(ic) = DZERO
        end if

      end do

      ! Look through edge flows that may have been provided by an exchange
      ! and incorporate them into the averaging arrays
      do m = 1, this%nedges
        if (this%nodedge(m) == n) then

          ! propsedge: (Q, area, nx, ny, distance)
          ihc = this%ihcedge(m)
          area = this%propsedge(2, m)

          ic = ic + 1
          nix(ic) = -this%propsedge(3, m)
          niy(ic) = -this%propsedge(4, m)
          di(ic) = this%propsedge(5, m)
          if (area > DZERO) then
            vi(ic) = this%propsedge(1, m) / area
          else
            vi(ic) = DZERO
          end if

        end if
      end do

      ! Assign number of vertical and horizontal connections
      ncz = iz
      nc = ic

      ! calculate z weight (wiz) and z velocity
      if (ncz == 1) then
        wiz(1) = DONE
      else
        dsumz = DZERO
        do iz = 1, ncz
          dsumz = dsumz + diz(iz)
        end do
        denom = (ncz - DONE)
        if (denom < DZERO) denom = DZERO
        dsumz = dsumz + DEM10 * dsumz
        do iz = 1, ncz
          if (dsumz > DZERO) wiz(iz) = DONE - diz(iz) / dsumz
          if (denom > 0) then
            wiz(iz) = wiz(iz) / denom
          else
            wiz(iz) = DZERO
          end if
        end do
      end if
      vz = DZERO
      do iz = 1, ncz
        vz = vz + wiz(iz) * viz(iz)
      end do

      ! distance-based weighting
      nc = ic
      dsumx = DZERO
      dsumy = DZERO
      dsumz = DZERO
      do ic = 1, nc
        wix(ic) = di(ic) * abs(nix(ic))
        wiy(ic) = di(ic) * abs(niy(ic))
        dsumx = dsumx + wix(ic)
        dsumy = dsumy + wiy(ic)
      end do

      ! Finish computing omega weights.  Add a tiny bit
      ! to dsum so that the normalized omega weight later
      ! evaluates to (essentially) 1 in the case of a single
      ! relevant connection, avoiding 0/0.
      dsumx = dsumx + DEM10 * dsumx
      dsumy = dsumy + DEM10 * dsumy
      do ic = 1, nc
        wix(ic) = (dsumx - wix(ic)) * abs(nix(ic))
        wiy(ic) = (dsumy - wiy(ic)) * abs(niy(ic))
      end do

      ! compute B weights
      dsumx = DZERO
      dsumy = DZERO
      do ic = 1, nc
        bix(ic) = wix(ic) * sign(DONE, nix(ic))
        biy(ic) = wiy(ic) * sign(DONE, niy(ic))
        dsumx = dsumx + wix(ic) * abs(nix(ic))
        dsumy = dsumy + wiy(ic) * abs(niy(ic))
      end do
      if (dsumx > DZERO) dsumx = DONE / dsumx
      if (dsumy > DZERO) dsumy = DONE / dsumy
      axy = DZERO
      ayx = DZERO
      do ic = 1, nc
        bix(ic) = bix(ic) * dsumx
        biy(ic) = biy(ic) * dsumy
        axy = axy + bix(ic) * niy(ic)
        ayx = ayx + biy(ic) * nix(ic)
      end do

      ! Calculate specific discharge.  The divide by zero checking below
      ! is problematic for cells with only one flow, such as can happen
      ! with triangular cells in corners.  In this case, the resulting
      ! cell velocity will be calculated as zero.  The method should be
      ! improved so that edge flows of zero are included in these
      ! calculations.  But this needs to be done with consideration for LGR
      ! cases in which flows are submitted from an exchange.
      vx = DZERO
      vy = DZERO
      do ic = 1, nc
        vx = vx + (bix(ic) - axy * biy(ic)) * vi(ic)
        vy = vy + (biy(ic) - ayx * bix(ic)) * vi(ic)
      end do
      denom = DONE - axy * ayx
      if (denom /= DZERO) then
        vx = vx / denom
        vy = vy / denom
      end if

      this%vcomp(1, n) = vx
      this%vcomp(2, n) = vy
      this%vcomp(3, n) = vz
      this%vmag(n) = sqrt(vx**2 + vy**2 + vz**2)

    end do

    ! cleanup
    deallocate (vi)
    deallocate (di)
    deallocate (nix)
    deallocate (niy)
    deallocate (wix)
    deallocate (wiy)
    deallocate (wiz)
    deallocate (bix)
    deallocate (biy)

  end subroutine calc_velocity

  !> @brief Reserve space for nedges cells that have an edge on them.
  !!
  !! todo: duplicated from NPF; should consolidate
  !! This must be called before the swf%allocate_arrays routine, which is
  !< called from swf%ar.
  subroutine increase_edge_count(this, nedges)
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B), intent(in) :: nedges

    this%nedges = this%nedges + nedges

  end subroutine increase_edge_count

  !> @brief Provide the swf package with edge properties
  !!
  !! todo: duplicated from NPF; should consolidate
  !<
  subroutine set_edge_properties(this, nodedge, ihcedge, q, area, nx, ny, &
                                 distance)
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B), intent(in) :: nodedge
    integer(I4B), intent(in) :: ihcedge
    real(DP), intent(in) :: q
    real(DP), intent(in) :: area
    real(DP), intent(in) :: nx
    real(DP), intent(in) :: ny
    real(DP), intent(in) :: distance
    ! local
    integer(I4B) :: lastedge

    this%lastedge = this%lastedge + 1
    lastedge = this%lastedge
    this%nodedge(lastedge) = nodedge
    this%ihcedge(lastedge) = ihcedge
    this%propsedge(1, lastedge) = q
    this%propsedge(2, lastedge) = area
    this%propsedge(3, lastedge) = nx
    this%propsedge(4, lastedge) = ny
    this%propsedge(5, lastedge) = distance

    ! If this is the last edge, then the next call must be starting a new
    ! edge properties assignment loop, so need to reset lastedge to 0
    if (this%lastedge == this%nedges) this%lastedge = 0

  end subroutine set_edge_properties

  !> @brief Save specific discharge in binary format to ibinun
  !!
  !! todo: should write 2d velocity; what about for 1D channel?
  !<
  subroutine sav_velocity(this, ibinun)
    ! dummy
    class(SwfDfwType) :: this !< this instance
    integer(I4B), intent(in) :: ibinun
    ! local
    character(len=16) :: text
    character(len=16), dimension(3) :: auxtxt
    integer(I4B) :: n
    integer(I4B) :: naux

    ! Write the header
    text = '      DATA-VCOMP'
    naux = 3
    auxtxt(:) = ['              vx', '              vy', '              vz']
    call this%dis%record_srcdst_list_header(text, this%name_model, &
                                            this%packName, this%name_model, &
                                            this%packName, naux, auxtxt, ibinun, &
                                            this%dis%nodes, this%iout)

    ! Write a zero for Q, and then write qx, qy, qz as aux variables
    do n = 1, this%dis%nodes
      call this%dis%record_mf6_list_entry(ibinun, n, n, DZERO, naux, &
                                          this%vcomp(:, n))
    end do

  end subroutine sav_velocity

  !> @brief Define the observation types available in the package
  !!
  !< Method to define the observation types available in the package.
  subroutine dfw_df_obs(this)
    ! dummy variables
    class(SwfDfwType) :: this !< this instance
    ! local variables
    integer(I4B) :: indx

    ! Store obs type and assign procedure pointer
    !    for ext-outflow observation type.
    call this%obs%StoreObsType('ext-outflow', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => dfwobsidprocessor

  end subroutine dfw_df_obs

  subroutine dfwobsidprocessor(obsrv, dis, inunitobs, iout)
    ! dummy
    type(ObserveType), intent(inout) :: obsrv
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: inunitobs
    integer(I4B), intent(in) :: iout
    ! local
    integer(I4B) :: n
    character(len=LINELENGTH) :: string

    ! Initialize variables
    string = obsrv%IDstring
    read (string, *) n

    if (n > 0) then
      obsrv%NodeNumber = n
    else
      errmsg = 'Error reading data from ID string'
      call store_error(errmsg)
      call store_error_unit(inunitobs)
    end if

  end subroutine dfwobsidprocessor

  !> @brief Save observations for the package
  !!
  !< Method to save simulated values for the package.
  subroutine dfw_bd_obs(this)
    ! dummy variables
    class(SwfDfwType) :: this !< this instance
    ! local variables
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: n
    real(DP) :: v
    character(len=100) :: msg
    type(ObserveType), pointer :: obsrv => null()

    ! Write simulated values for all observations
    if (this%obs%npakobs > 0) then
      call this%obs%obs_bd_clear()
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv
        do j = 1, obsrv%indxbnds_count
          n = obsrv%indxbnds(j)
          v = DZERO
          select case (obsrv%ObsTypeId)
          case default
            msg = 'Unrecognized observation type: '//trim(obsrv%ObsTypeId)
            call store_error(msg)
          end select
          call this%obs%SaveOneSimval(obsrv, v)
        end do
      end do

      ! write summary of package error messages
      if (count_errors() > 0) then
        call store_error_filename(this%input_fname)
      end if
    end if

  end subroutine dfw_bd_obs

  !> @brief Read and prepare observations for a package
  !!
  !< Method to read and prepare observations for a package.
  subroutine dfw_rp_obs(this)
    ! modules
    use TdisModule, only: kper
    ! dummy
    class(SwfDfwType), intent(inout) :: this !< this instance
    ! local
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: nn1
    class(ObserveType), pointer :: obsrv => null()

    ! process each package observation
    ! only done the first stress period since boundaries are fixed
    ! for the simulation
    if (kper == 1) then
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv

        ! get node number 1
        nn1 = obsrv%NodeNumber
        if (nn1 < 1 .or. nn1 > this%dis%nodes) then
          write (errmsg, '(a,1x,a,1x,i0,1x,a,1x,i0,a)') &
            trim(adjustl(obsrv%ObsTypeId)), &
            'reach must be greater than 0 and less than or equal to', &
            this%dis%nodes, '(specified value is ', nn1, ')'
          call store_error(errmsg)
        else
          if (obsrv%indxbnds_count == 0) then
            call obsrv%AddObsIndex(nn1)
          else
            errmsg = 'Programming error in dfw_rp_obs'
            call store_error(errmsg)
          end if
        end if

        ! check that node number 1 is valid; call store_error if not
        do j = 1, obsrv%indxbnds_count
          nn1 = obsrv%indxbnds(j)
          if (nn1 < 1 .or. nn1 > this%dis%nodes) then
            write (errmsg, '(a,1x,a,1x,i0,1x,a,1x,i0,a)') &
              trim(adjustl(obsrv%ObsTypeId)), &
              'reach must be greater than 0 and less than or equal to', &
              this%dis%nodes, '(specified value is ', nn1, ')'
            call store_error(errmsg)
          end if
        end do
      end do

      ! evaluate if there are any observation errors
      if (count_errors() > 0) then
        call store_error_filename(this%input_fname)
      end if

    end if

  end subroutine dfw_rp_obs

end module SwfDfwModule
