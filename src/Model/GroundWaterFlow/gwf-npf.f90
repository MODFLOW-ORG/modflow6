module GwfNpfModule
  use KindModule, only: DP, I4B
  use SimVariablesModule, only: errmsg
  use ConstantsModule, only: DZERO, DEM9, DEM8, DEM7, DEM6, DEM2, &
                             DHALF, DP9, DONE, DTWO, &
                             DHNOFLO, DHDRY, DEM10, &
                             LENMEMPATH, LENVARNAME, LINELENGTH, &
                             C3D_VERTICAL
  use SmoothingModule, only: sQuadraticSaturation, &
                             sQuadraticSaturationDerivative
  use NumericalPackageModule, only: NumericalPackageType
  use GwfNpfOptionsModule, only: GwfNpfOptionsType
  use BaseDisModule, only: DisBaseType
  use GwfIcModule, only: GwfIcType
  use GwfVscModule, only: GwfVscType
  use Xt3dModule, only: Xt3dType
  use SpdisWorkArrayModule, only: SpdisWorkArrayType
  use InputOutputModule, only: GetUnit, openfile
  use TvkModule, only: TvkType, tvk_cr
  use MemoryManagerModule, only: mem_allocate, mem_reallocate, &
                                 mem_deallocate, mem_setptr, &
                                 mem_reassignptr
  use MatrixBaseModule
  use HGeoUtilModule, only: hyeff
  use GwfConductanceUtilsModule, only: hcond, vcond, &
                                       condmean, thksatnm, &
                                       CCOND_HMEAN

  implicit none

  private
  public :: GwfNpfType
  public :: npf_cr

  type, extends(NumericalPackageType) :: GwfNpfType

    type(GwfIcType), pointer :: ic => null() !< initial conditions object
    type(GwfVscType), pointer :: vsc => null() !< viscosity object
    type(Xt3dType), pointer :: xt3d => null() !< xt3d pointer
    integer(I4B), pointer :: iname => null() !< length of variable names
    character(len=24), dimension(:), pointer :: aname => null() !< variable names
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    real(DP), dimension(:), pointer, contiguous :: hnew => null() !< pointer to model xnew
    integer(I4B), pointer :: ixt3d => null() !< xt3d flag (0 is off, 1 is lhs, 2 is rhs)
    integer(I4B), pointer :: ixt3drhs => null() !< xt3d rhs flag, xt3d rhs is set active if 1
    integer(I4B), pointer :: iperched => null() !< vertical flow corrections if 1
    integer(I4B), pointer :: ivarcv => null() !< CV is function of water table
    integer(I4B), pointer :: idewatcv => null() !< CV may be a discontinuous function of water table
    integer(I4B), pointer :: ithickstrt => null() !< thickstrt option flag
    integer(I4B), pointer :: igwfnewtonur => null() !< newton head dampening using node bottom option flag
    integer(I4B), pointer :: icalcspdis => null() !< Calculate specific discharge at cell centers
    integer(I4B), pointer :: isavspdis => null() !< Save specific discharge at cell centers
    integer(I4B), pointer :: isavsat => null() !< Save sat to budget file
    real(DP), pointer :: hnoflo => null() !< default is 1.e30
    real(DP), pointer :: satomega => null() !< newton-raphson saturation omega
    integer(I4B), pointer :: irewet => null() !< rewetting (0:off, 1:on)
    integer(I4B), pointer :: iwetit => null() !< wetting interval (default is 1)
    integer(I4B), pointer :: ihdwet => null() !< (0 or not 0)
    integer(I4B), pointer :: icellavg => null() !< harmonic(0), logarithmic(1), or arithmetic thick-log K (2)
    real(DP), pointer :: wetfct => null() !< wetting factor
    real(DP), pointer :: hdry => null() !< default is -1.d30
    integer(I4B), dimension(:), pointer, contiguous :: icelltype => null() !< confined (0) or convertible (1)
    integer(I4B), dimension(:), pointer, contiguous :: ithickstartflag => null() !< array of flags for handling the thickstrt option
    !
    ! K properties
    real(DP), dimension(:), pointer, contiguous :: k11 => null() !< hydraulic conductivity; if anisotropic, then this is Kx prior to rotation
    real(DP), dimension(:), pointer, contiguous :: k22 => null() !< hydraulic conductivity; if specified then this is Ky prior to rotation
    real(DP), dimension(:), pointer, contiguous :: k33 => null() !< hydraulic conductivity; if specified then this is Kz prior to rotation
    real(DP), dimension(:), pointer, contiguous :: k11input => null() !< hydraulic conductivity originally specified by user prior to TVK or VSC modification
    real(DP), dimension(:), pointer, contiguous :: k22input => null() !< hydraulic conductivity originally specified by user prior to TVK or VSC modification
    real(DP), dimension(:), pointer, contiguous :: k33input => null() !< hydraulic conductivity originally specified by user prior to TVK or VSC modification
    integer(I4B), pointer :: iavgkeff => null() !< effective conductivity averaging (0: harmonic, 1: arithmetic)
    integer(I4B), pointer :: ik22 => null() !< flag that k22 is specified
    integer(I4B), pointer :: ik33 => null() !< flag that k33 is specified
    integer(I4B), pointer :: ik22overk => null() !< flag that k22 is specified as anisotropy ratio
    integer(I4B), pointer :: ik33overk => null() !< flag that k33 is specified as anisotropy ratio
    integer(I4B), pointer :: iangle1 => null() !< flag to indicate angle1 was read
    integer(I4B), pointer :: iangle2 => null() !< flag to indicate angle2 was read
    integer(I4B), pointer :: iangle3 => null() !< flag to indicate angle3 was read
    real(DP), dimension(:), pointer, contiguous :: angle1 => null() !< k ellipse rotation in xy plane around z axis (yaw)
    real(DP), dimension(:), pointer, contiguous :: angle2 => null() !< k ellipse rotation up from xy plane around y axis (pitch)
    real(DP), dimension(:), pointer, contiguous :: angle3 => null() !< k tensor rotation around x axis (roll)
    !
    integer(I4B), pointer :: iwetdry => null() !< flag to indicate angle1 was read
    real(DP), dimension(:), pointer, contiguous :: wetdry => null() !< wetdry array
    real(DP), dimension(:), pointer, contiguous :: sat => null() !< saturation (0. to 1.) for each cell
    real(DP), dimension(:), pointer, contiguous :: condsat => null() !< saturated conductance (symmetric array)
    integer(I4B), dimension(:), pointer, contiguous :: ibotnode => null() !< bottom node used if igwfnewtonur /= 0
    ! spdis machinery:
    real(DP), dimension(:, :), pointer, contiguous :: spdis => null() !< specific discharge : qx, qy, qz (nodes, 3)
    integer(I4B), pointer :: nedges => null() !< number of cell edges
    integer(I4B), pointer :: lastedge => null() !< last edge number
    integer(I4B), dimension(:), pointer, contiguous :: nodedge => null() !< array of node numbers that have edges
    integer(I4B), dimension(:), pointer, contiguous :: ihcedge => null() !< edge type (horizontal or vertical)
    real(DP), dimension(:, :), pointer, contiguous :: propsedge => null() !< edge properties (Q, area, nx, ny, distance)
    integer(I4B), dimension(:), pointer, contiguous :: iedge_ptr => null() !< csr pointer into edge index array
    integer(I4B), dimension(:), pointer, contiguous :: edge_idxs => null() !< sorted edge indexes for faster lookup
    type(SpdisWorkArrayType), pointer :: spdis_wa => null() !< work arrays for spdis calculation
    !
    integer(I4B), pointer :: intvk => null() ! TVK (time-varying K) unit number (0 if unused)
    integer(I4B), pointer :: invsc => null() ! VSC (viscosity) unit number (0 if unused); viscosity leads to time-varying K's
    type(TvkType), pointer :: tvk => null() ! TVK object
    integer(I4B), pointer :: kchangeper => null() ! last stress period in which any node K (or K22, or K33) values were changed (0 if unchanged from start of simulation)
    integer(I4B), pointer :: kchangestp => null() ! last time step in which any node K (or K22, or K33) values were changed (0 if unchanged from start of simulation)
    integer(I4B), dimension(:), pointer, contiguous :: nodekchange => null() ! grid array of flags indicating for each node whether its K (or K22, or K33) value changed (1) at (kchangeper, kchangestp) or not (0)

  contains

    procedure :: npf_df
    procedure :: npf_ac
    procedure :: npf_mc
    procedure :: npf_ar
    procedure :: npf_rp
    procedure :: npf_ad
    procedure :: npf_cf
    procedure :: npf_fc
    procedure :: npf_fn
    procedure :: npf_cq
    procedure :: npf_save_model_flows
    procedure :: npf_nur
    procedure :: npf_print_model_flows
    procedure :: npf_da
    procedure, private :: thksat => sgwf_npf_thksat
    procedure, private :: qcalc => sgwf_npf_qcalc
    procedure, private :: wd => sgwf_npf_wetdry
    procedure, private :: wdmsg => sgwf_npf_wdmsg
    procedure :: allocate_scalars
    procedure, private :: store_original_k_arrays
    procedure, private :: allocate_arrays
    procedure, private :: source_options
    procedure, private :: source_griddata
    procedure, private :: log_options
    procedure, private :: log_griddata
    procedure, private :: set_options
    procedure, private :: check_options
    procedure, private :: prepcheck
    procedure, private :: preprocess_input
    procedure, private :: calc_condsat
    procedure, private :: calc_initial_sat
    procedure, public :: rewet_check
    procedure, public :: hy_eff
    procedure, public :: calc_spdis
    procedure, public :: sav_spdis
    procedure, public :: sav_sat
    procedure, public :: increase_edge_count
    procedure, public :: set_edge_properties
    procedure, public :: calcSatThickness
    procedure, private :: calc_max_conns
    procedure, private :: prepare_edge_lookup
  end type

contains

  !> @brief Create a new NPF object. Pass a inunit value of 0 if npf data will
  !! initialized from memory
  !<
  subroutine npf_cr(npfobj, name_model, input_mempath, inunit, iout)
    ! -- modules
    use KindModule, only: LGP
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy
    type(GwfNpfType), pointer :: npfobj
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    ! -- formats
    character(len=*), parameter :: fmtheader = &
      "(1x, /1x, 'NPF -- NODE PROPERTY FLOW PACKAGE, VERSION 1, 3/30/2015', &
       &' INPUT READ FROM MEMPATH: ', A, /)"
    !
    ! -- Create the object
    allocate (npfobj)
    !
    ! -- create name and memory path
    call npfobj%set_names(1, name_model, 'NPF', 'NPF', input_mempath)
    !
    ! -- Allocate scalars
    call npfobj%allocate_scalars()
    !
    ! -- Set variables
    npfobj%inunit = inunit
    npfobj%iout = iout
    !
    ! -- check if npf is enabled
    if (inunit > 0) then
      !
      ! -- Print a message identifying the node property flow package.
      write (iout, fmtheader) input_mempath
    end if

    ! allocate spdis structure
    allocate (npfobj%spdis_wa)

  end subroutine npf_cr

  !> @brief Define the NPF package instance
  !!
  !! This is a hybrid routine: it either reads the options for this package
  !! from the input file, or the optional argument @param npf_options
  !! should be passed. A consistency check is performed, and finally
  !! xt3d_df is called, when enabled.
  !<
  subroutine npf_df(this, dis, xt3d, ingnc, invsc, npf_options)
    ! -- modules
    use SimModule, only: store_error
    use Xt3dModule, only: xt3d_cr
    ! -- dummy
    class(GwfNpftype) :: this !< instance of the NPF package
    class(DisBaseType), pointer, intent(inout) :: dis !< the pointer to the discretization
    type(Xt3dType), pointer :: xt3d !< the pointer to the XT3D 'package'
    integer(I4B), intent(in) :: ingnc !< ghostnodes enabled? (>0 means yes)
    integer(I4B), intent(in) :: invsc !< viscosity enabled? (>0 means yes)
    type(GwfNpfOptionsType), optional, intent(in) :: npf_options !< the optional options, for when not constructing from file
    !
    ! -- Set a pointer to dis
    this%dis => dis
    !
    ! -- Set flag signifying whether vsc is active
    if (invsc > 0) this%invsc = invsc
    !
    if (.not. present(npf_options)) then
      !
      ! -- source options
      call this%source_options()
      !
      ! -- allocate arrays
      call this%allocate_arrays(this%dis%nodes, this%dis%njas)
      !
      ! -- source griddata, set, and convert/check the input
      call this%source_griddata()
      call this%prepcheck()
    else
      call this%set_options(npf_options)
      !
      ! -- allocate arrays
      call this%allocate_arrays(this%dis%nodes, this%dis%njas)
    end if
    !
    call this%check_options()
    !
    ! -- Save pointer to xt3d object
    this%xt3d => xt3d
    if (this%ixt3d /= 0) xt3d%ixt3d = this%ixt3d
    call this%xt3d%xt3d_df(dis)
    !
    ! -- Ensure GNC and XT3D are not both on at the same time
    if (this%ixt3d /= 0 .and. ingnc > 0) then
      call store_error('Error in model '//trim(this%name_model)// &
                       '.  The XT3D option cannot be used with the GNC &
                       &Package.', terminate=.TRUE.)
    end if
  end subroutine npf_df

  !> @brief Add connections for extended neighbors to the sparse matrix
  !<
  subroutine npf_ac(this, moffset, sparse)
    ! -- modules
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(GwfNpftype) :: this
    integer(I4B), intent(in) :: moffset
    type(sparsematrix), intent(inout) :: sparse
    !
    ! -- Add extended neighbors (neighbors of neighbors)
    if (this%ixt3d /= 0) call this%xt3d%xt3d_ac(moffset, sparse)
  end subroutine npf_ac

  !> @brief Map connections and construct iax, jax, and idxglox
  !<
  subroutine npf_mc(this, moffset, matrix_sln)
    ! -- dummy
    class(GwfNpftype) :: this
    integer(I4B), intent(in) :: moffset
    class(MatrixBaseType), pointer :: matrix_sln
    !
    if (this%ixt3d /= 0) call this%xt3d%xt3d_mc(moffset, matrix_sln)
  end subroutine npf_mc

  !> @brief Allocate and read this NPF instance
  !!
  !! Allocate remaining package arrays, preprocess the input data and
  !! call *_ar on xt3d, when active.
  !<
  subroutine npf_ar(this, ic, vsc, ibound, hnew)
    ! -- modules
    use MemoryManagerModule, only: mem_reallocate
    ! -- dummy
    class(GwfNpftype) :: this !< instance of the NPF package
    type(GwfIcType), pointer, intent(in) :: ic !< initial conditions
    type(GwfVscType), pointer, intent(in) :: vsc !< viscosity package
    integer(I4B), dimension(:), pointer, contiguous, intent(inout) :: ibound !< model ibound array
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: hnew !< pointer to model head array
    ! -- local
    integer(I4B) :: n
    !
    ! -- Store pointers to arguments that were passed in
    this%ic => ic
    this%ibound => ibound
    this%hnew => hnew
    !
    if (this%icalcspdis == 1) then
      call mem_reallocate(this%spdis, 3, this%dis%nodes, 'SPDIS', this%memoryPath)
      call mem_reallocate(this%nodedge, this%nedges, 'NODEDGE', this%memoryPath)
      call mem_reallocate(this%ihcedge, this%nedges, 'IHCEDGE', this%memoryPath)
      call mem_reallocate(this%propsedge, 5, this%nedges, 'PROPSEDGE', &
                          this%memoryPath)
      call mem_reallocate(this%iedge_ptr, this%dis%nodes + 1, &
                          'NREDGESNODE', this%memoryPath)
      call mem_reallocate(this%edge_idxs, this%nedges, &
                          'EDGEIDXS', this%memoryPath)

      do n = 1, this%nedges
        this%edge_idxs(n) = 0
      end do
      do n = 1, this%dis%nodes
        this%iedge_ptr(n) = 0
        this%spdis(:, n) = DZERO
      end do
    end if
    !
    ! -- Store pointer to VSC if active
    if (this%invsc /= 0) then
      this%vsc => vsc
    end if
    !
    ! -- allocate arrays to store original user input in case TVK/VSC modify them
    if (this%invsc > 0) then
      !
      ! -- Reallocate arrays that store user-input values.
      call mem_reallocate(this%k11input, this%dis%nodes, 'K11INPUT', &
                          this%memoryPath)
      call mem_reallocate(this%k22input, this%dis%nodes, 'K22INPUT', &
                          this%memoryPath)
      call mem_reallocate(this%k33input, this%dis%nodes, 'K33INPUT', &
                          this%memoryPath)
      ! Allocate arrays that will store the original K values.  When VSC active,
      ! the current Kxx arrays carry the viscosity-adjusted K values.
      ! This approach leverages existing functionality that makes use of K.
      call this%store_original_k_arrays(this%dis%nodes, this%dis%njas)
    end if
    !
    ! -- preprocess data
    call this%preprocess_input()
    !
    ! -- xt3d
    if (this%ixt3d /= 0) then
      call this%xt3d%xt3d_ar(ibound, this%k11, this%ik33, this%k33, &
                             this%sat, this%ik22, this%k22, &
                             this%iangle1, this%iangle2, this%iangle3, &
                             this%angle1, this%angle2, this%angle3, &
                             this%inewton, this%icelltype)
    end if
    !
    ! -- TVK
    if (this%intvk /= 0) then
      call this%tvk%ar(this%dis)
    end if
  end subroutine npf_ar

  !> @brief Read and prepare method for package
  !!
  !! Read and prepare NPF stress period data.
  !<
  subroutine npf_rp(this)
    implicit none
    ! -- dummy
    class(GwfNpfType) :: this
    !
    ! -- TVK
    if (this%intvk /= 0) then
      call this%tvk%rp()
    end if
  end subroutine npf_rp

  !> @brief Advance
  !!
  !! Sets hold (head old) to bot whenever a wettable cell is dry
  !<
  subroutine npf_ad(this, nodes, hold, hnew, irestore)
    ! -- modules
    use TdisModule, only: kper, kstp
    !
    implicit none
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: nodes
    real(DP), dimension(nodes), intent(inout) :: hold
    real(DP), dimension(nodes), intent(inout) :: hnew
    integer(I4B), intent(in) :: irestore
    ! -- local
    integer(I4B) :: n
    !
    ! -- loop through all cells and set hold=bot if wettable cell is dry
    if (this%irewet > 0) then
      do n = 1, this%dis%nodes
        if (this%wetdry(n) == DZERO) cycle
        if (this%ibound(n) /= 0) cycle
        hold(n) = this%dis%bot(n)
      end do
      !
      ! -- if restore state, then set hnew to DRY if it is a dry wettable cell
      do n = 1, this%dis%nodes
        if (this%wetdry(n) == DZERO) cycle
        if (this%ibound(n) /= 0) cycle
        hnew(n) = DHDRY
      end do
    end if
    !
    ! -- TVK
    if (this%intvk /= 0) then
      call this%tvk%ad()
    end if
    !
    ! -- VSC
    ! -- Hit the TVK-updated K's with VSC correction before calling/updating condsat
    if (this%invsc /= 0) then
      call this%vsc%update_k_with_vsc()
    end if
    !
    ! -- If any K values have changed, we need to update CONDSAT or XT3D arrays
    if (this%kchangeper == kper .and. this%kchangestp == kstp) then
      if (this%ixt3d == 0) then
        !
        ! -- Update the saturated conductance for all connections
        ! -- of the affected nodes
        do n = 1, this%dis%nodes
          if (this%nodekchange(n) == 1) then
            call this%calc_condsat(n, .false.)
          end if
        end do
      else
        !
        ! -- Recompute XT3D coefficients for permanently confined connections
        if (this%xt3d%lamatsaved .and. .not. this%xt3d%ldispersion) then
          call this%xt3d%xt3d_fcpc(this%dis%nodes, .true.)
        end if
      end if
    end if
  end subroutine npf_ad

  !> @brief Routines associated fill coefficients
  !<
  subroutine npf_cf(this, kiter, nodes, hnew)
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B) :: kiter
    integer(I4B), intent(in) :: nodes
    real(DP), intent(inout), dimension(nodes) :: hnew
    ! -- local
    integer(I4B) :: n
    real(DP) :: satn
    !
    ! -- Perform wetting and drying
    if (this%inewton /= 1) then
      call this%wd(kiter, hnew)
    end if
    !
    ! -- Calculate saturation for convertible cells
    do n = 1, this%dis%nodes
      if (this%icelltype(n) /= 0) then
        if (this%ibound(n) == 0) then
          satn = DZERO
        else
          call this%thksat(n, hnew(n), satn)
        end if
        this%sat(n) = satn
      end if
    end do
  end subroutine npf_cf

  !> @brief Formulate coefficients
  !<
  subroutine npf_fc(this, kiter, matrix_sln, idxglo, rhs, hnew)
    ! -- modules
    use ConstantsModule, only: DONE
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    real(DP), intent(inout), dimension(:) :: hnew
    ! -- local
    integer(I4B) :: n, m, ii, idiag, ihc
    integer(I4B) :: isymcon, idiagm
    real(DP) :: hyn, hym
    real(DP) :: cond
    !
    ! -- Calculate conductance and put into amat
    !
    if (this%ixt3d /= 0) then
      call this%xt3d%xt3d_fc(kiter, matrix_sln, idxglo, rhs, hnew)
    else
      do n = 1, this%dis%nodes
        do ii = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
          if (this%dis%con%mask(ii) == 0) cycle

          m = this%dis%con%ja(ii)
          !
          ! -- Calculate conductance only for upper triangle but insert into
          !    upper and lower parts of amat.
          if (m < n) cycle
          ihc = this%dis%con%ihc(this%dis%con%jas(ii))
          hyn = this%hy_eff(n, m, ihc, ipos=ii)
          hym = this%hy_eff(m, n, ihc, ipos=ii)
          !
          ! -- Vertical connection
          if (ihc == C3D_VERTICAL) then
            !
            ! -- Calculate vertical conductance
            cond = vcond(this%ibound(n), this%ibound(m), &
                         this%icelltype(n), this%icelltype(m), this%inewton, &
                         this%ivarcv, this%idewatcv, &
                         this%condsat(this%dis%con%jas(ii)), hnew(n), hnew(m), &
                         hyn, hym, &
                         this%sat(n), this%sat(m), &
                         this%dis%top(n), this%dis%top(m), &
                         this%dis%bot(n), this%dis%bot(m), &
                         this%dis%con%hwva(this%dis%con%jas(ii)))
            !
            ! -- Vertical flow for perched conditions
            if (this%iperched /= 0) then
              if (this%icelltype(m) /= 0) then
                if (hnew(m) < this%dis%top(m)) then
                  !
                  ! -- Fill row n
                  idiag = this%dis%con%ia(n)
                  rhs(n) = rhs(n) - cond * this%dis%bot(n)
                  call matrix_sln%add_value_pos(idxglo(idiag), -cond)
                  !
                  ! -- Fill row m
                  isymcon = this%dis%con%isym(ii)
                  call matrix_sln%add_value_pos(idxglo(isymcon), cond)
                  rhs(m) = rhs(m) + cond * this%dis%bot(n)
                  !
                  ! -- cycle the connection loop
                  cycle
                end if
              end if
            end if
            !
          else
            !
            ! -- Horizontal conductance
            cond = hcond(this%ibound(n), this%ibound(m), &
                         this%icelltype(n), this%icelltype(m), &
                         this%inewton, &
                         this%dis%con%ihc(this%dis%con%jas(ii)), &
                         this%icellavg, &
                         this%condsat(this%dis%con%jas(ii)), &
                         hnew(n), hnew(m), this%sat(n), this%sat(m), hyn, hym, &
                         this%dis%top(n), this%dis%top(m), &
                         this%dis%bot(n), this%dis%bot(m), &
                         this%dis%con%cl1(this%dis%con%jas(ii)), &
                         this%dis%con%cl2(this%dis%con%jas(ii)), &
                         this%dis%con%hwva(this%dis%con%jas(ii)))
          end if
          !
          ! -- Fill row n
          idiag = this%dis%con%ia(n)
          call matrix_sln%add_value_pos(idxglo(ii), cond)
          call matrix_sln%add_value_pos(idxglo(idiag), -cond)
          !
          ! -- Fill row m
          isymcon = this%dis%con%isym(ii)
          idiagm = this%dis%con%ia(m)
          call matrix_sln%add_value_pos(idxglo(isymcon), cond)
          call matrix_sln%add_value_pos(idxglo(idiagm), -cond)
        end do
      end do
      !
    end if
  end subroutine npf_fc

  !> @brief Fill newton terms
  !<
  subroutine npf_fn(this, kiter, matrix_sln, idxglo, rhs, hnew)
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    real(DP), intent(inout), dimension(:) :: hnew
    ! -- local
    integer(I4B) :: nodes, nja
    integer(I4B) :: n, m, ii, idiag
    integer(I4B) :: isymcon, idiagm
    integer(I4B) :: iups
    integer(I4B) :: idn
    real(DP) :: cond
    real(DP) :: consterm
    real(DP) :: filledterm
    real(DP) :: derv
    real(DP) :: hds
    real(DP) :: term
    real(DP) :: topup
    real(DP) :: botup
    !
    ! -- add newton terms to solution matrix
    nodes = this%dis%nodes
    nja = this%dis%con%nja
    if (this%ixt3d /= 0) then
      call this%xt3d%xt3d_fn(kiter, nodes, nja, matrix_sln, idxglo, rhs, hnew)
    else
      !
      do n = 1, nodes
        idiag = this%dis%con%ia(n)
        do ii = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
          if (this%dis%con%mask(ii) == 0) cycle

          m = this%dis%con%ja(ii)
          isymcon = this%dis%con%isym(ii)
          ! work on upper triangle
          if (m < n) cycle
          if (this%dis%con%ihc(this%dis%con%jas(ii)) == 0 .and. &
              this%ivarcv == 0) then
            !call this%vcond(n,m,hnew(n),hnew(m),ii,cond)
            ! do nothing
          else
            ! determine upstream node
            iups = m
            if (hnew(m) < hnew(n)) iups = n
            idn = n
            if (iups == n) idn = m
            !
            ! -- no newton terms if upstream cell is confined
            if (this%icelltype(iups) == 0) cycle
            !
            ! -- Set the upstream top and bot, and then recalculate for a
            !    vertically staggered horizontal connection
            topup = this%dis%top(iups)
            botup = this%dis%bot(iups)
            if (this%dis%con%ihc(this%dis%con%jas(ii)) == 2) then
              topup = min(this%dis%top(n), this%dis%top(m))
              botup = max(this%dis%bot(n), this%dis%bot(m))
            end if
            !
            ! get saturated conductivity for derivative
            cond = this%condsat(this%dis%con%jas(ii))
            !
            ! compute additional term
            consterm = -cond * (hnew(iups) - hnew(idn)) !needs to use hwadi instead of hnew(idn)
            !filledterm = cond
            filledterm = matrix_sln%get_value_pos(idxglo(ii))
            derv = sQuadraticSaturationDerivative(topup, botup, hnew(iups), &
                                                  this%satomega)
            idiagm = this%dis%con%ia(m)
            ! fill jacobian for n being the upstream node
            if (iups == n) then
              hds = hnew(m)
              !isymcon =  this%dis%con%isym(ii)
              term = consterm * derv
              rhs(n) = rhs(n) + term * hnew(n) !+ amat(idxglo(isymcon)) * (dwadi * hds - hds) !need to add dwadi
              rhs(m) = rhs(m) - term * hnew(n) !- amat(idxglo(isymcon)) * (dwadi * hds - hds) !need to add dwadi
              ! fill in row of n
              call matrix_sln%add_value_pos(idxglo(idiag), term)
              ! fill newton term in off diagonal if active cell
              if (this%ibound(n) > 0) then
                filledterm = matrix_sln%get_value_pos(idxglo(ii))
                call matrix_sln%set_value_pos(idxglo(ii), filledterm) !* dwadi !need to add dwadi
              end if
              !fill row of m
              filledterm = matrix_sln%get_value_pos(idxglo(idiagm))
              call matrix_sln%set_value_pos(idxglo(idiagm), filledterm) !- filledterm * (dwadi - DONE) !need to add dwadi
              ! fill newton term in off diagonal if active cell
              if (this%ibound(m) > 0) then
                call matrix_sln%add_value_pos(idxglo(isymcon), -term)
              end if
              ! fill jacobian for m being the upstream node
            else
              hds = hnew(n)
              term = -consterm * derv
              rhs(n) = rhs(n) + term * hnew(m) !+ amat(idxglo(ii)) * (dwadi * hds - hds) !need to add dwadi
              rhs(m) = rhs(m) - term * hnew(m) !- amat(idxglo(ii)) * (dwadi * hds - hds) !need to add dwadi
              ! fill in row of n
              filledterm = matrix_sln%get_value_pos(idxglo(idiag))
              call matrix_sln%set_value_pos(idxglo(idiag), filledterm) !- filledterm * (dwadi - DONE) !need to add dwadi
              ! fill newton term in off diagonal if active cell
              if (this%ibound(n) > 0) then
                call matrix_sln%add_value_pos(idxglo(ii), term)
              end if
              !fill row of m
              call matrix_sln%add_value_pos(idxglo(idiagm), -term)
              ! fill newton term in off diagonal if active cell
              if (this%ibound(m) > 0) then
                filledterm = matrix_sln%get_value_pos(idxglo(isymcon))
                call matrix_sln%set_value_pos(idxglo(isymcon), filledterm) !* dwadi  !need to add dwadi
              end if
            end if
          end if

        end do
      end do
      !
    end if
  end subroutine npf_fn

  !> @brief Under-relaxation
  !!
  !! Under-relaxation of Groundwater Flow Model Heads for current outer
  !! iteration using the cell bottoms at the bottom of the model
  !<
  subroutine npf_nur(this, neqmod, x, xtemp, dx, inewtonur, dxmax, locmax)
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: neqmod
    real(DP), dimension(neqmod), intent(inout) :: x
    real(DP), dimension(neqmod), intent(in) :: xtemp
    real(DP), dimension(neqmod), intent(inout) :: dx
    integer(I4B), intent(inout) :: inewtonur
    real(DP), intent(inout) :: dxmax
    integer(I4B), intent(inout) :: locmax
    ! -- local
    integer(I4B) :: n
    real(DP) :: botm
    real(DP) :: xx
    real(DP) :: dxx
    !
    ! -- Newton-Raphson under-relaxation
    do n = 1, this%dis%nodes
      if (this%ibound(n) < 1) cycle
      if (this%icelltype(n) > 0) then
        botm = this%dis%bot(this%ibotnode(n))
        ! -- only apply Newton-Raphson under-relaxation if
        !    solution head is below the bottom of the model
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
  end subroutine npf_nur

  !> @brief Calculate flowja
  !<
  subroutine npf_cq(this, hnew, flowja)
    ! -- dummy
    class(GwfNpfType) :: this
    real(DP), intent(inout), dimension(:) :: hnew
    real(DP), intent(inout), dimension(:) :: flowja
    ! -- local
    integer(I4B) :: n, ipos, m
    real(DP) :: qnm
    !
    ! -- Calculate the flow across each cell face and store in flowja
    !
    if (this%ixt3d /= 0) then
      call this%xt3d%xt3d_flowja(hnew, flowja)
    else
      !
      do n = 1, this%dis%nodes
        do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
          m = this%dis%con%ja(ipos)
          if (m < n) cycle
          call this%qcalc(n, m, hnew(n), hnew(m), ipos, qnm)
          flowja(ipos) = qnm
          flowja(this%dis%con%isym(ipos)) = -qnm
        end do
      end do
      !
    end if
  end subroutine npf_cq

  !> @brief Fractional cell saturation
  !<
  subroutine sgwf_npf_thksat(this, n, hn, thksat)
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: n
    real(DP), intent(in) :: hn
    real(DP), intent(inout) :: thksat
    !
    ! -- Standard Formulation
    if (hn >= this%dis%top(n)) then
      thksat = DONE
    else
      thksat = (hn - this%dis%bot(n)) / (this%dis%top(n) - this%dis%bot(n))
    end if
    !
    ! -- Newton-Raphson Formulation
    if (this%inewton /= 0) then
      thksat = sQuadraticSaturation(this%dis%top(n), this%dis%bot(n), hn, &
                                    this%satomega)
    end if
  end subroutine sgwf_npf_thksat

  !> @brief Flow between two cells
  !<
  subroutine sgwf_npf_qcalc(this, n, m, hn, hm, icon, qnm)
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    real(DP), intent(in) :: hn
    real(DP), intent(in) :: hm
    integer(I4B), intent(in) :: icon
    real(DP), intent(inout) :: qnm
    ! -- local
    real(DP) :: hyn, hym
    real(DP) :: condnm
    real(DP) :: hntemp, hmtemp
    integer(I4B) :: ihc
    !
    ! -- Initialize
    ihc = this%dis%con%ihc(this%dis%con%jas(icon))
    hyn = this%hy_eff(n, m, ihc, ipos=icon)
    hym = this%hy_eff(m, n, ihc, ipos=icon)
    !
    ! -- Calculate conductance
    if (ihc == C3D_VERTICAL) then
      condnm = vcond(this%ibound(n), this%ibound(m), &
                     this%icelltype(n), this%icelltype(m), this%inewton, &
                     this%ivarcv, this%idewatcv, &
                     this%condsat(this%dis%con%jas(icon)), hn, hm, &
                     hyn, hym, &
                     this%sat(n), this%sat(m), &
                     this%dis%top(n), this%dis%top(m), &
                     this%dis%bot(n), this%dis%bot(m), &
                     this%dis%con%hwva(this%dis%con%jas(icon)))
    else
      condnm = hcond(this%ibound(n), this%ibound(m), &
                     this%icelltype(n), this%icelltype(m), &
                     this%inewton, &
                     this%dis%con%ihc(this%dis%con%jas(icon)), &
                     this%icellavg, &
                     this%condsat(this%dis%con%jas(icon)), &
                     hn, hm, this%sat(n), this%sat(m), hyn, hym, &
                     this%dis%top(n), this%dis%top(m), &
                     this%dis%bot(n), this%dis%bot(m), &
                     this%dis%con%cl1(this%dis%con%jas(icon)), &
                     this%dis%con%cl2(this%dis%con%jas(icon)), &
                     this%dis%con%hwva(this%dis%con%jas(icon)))
    end if
    !
    ! -- Initialize hntemp and hmtemp
    hntemp = hn
    hmtemp = hm
    !
    ! -- Check and adjust for dewatered conditions
    if (this%iperched /= 0) then
      if (this%dis%con%ihc(this%dis%con%jas(icon)) == 0) then
        if (n > m) then
          if (this%icelltype(n) /= 0) then
            if (hn < this%dis%top(n)) hntemp = this%dis%bot(m)
          end if
        else
          if (this%icelltype(m) /= 0) then
            if (hm < this%dis%top(m)) hmtemp = this%dis%bot(n)
          end if
        end if
      end if
    end if
    !
    ! -- Calculate flow positive into cell n
    qnm = condnm * (hmtemp - hntemp)
  end subroutine sgwf_npf_qcalc

  !> @brief Record flowja and calculate specific discharge if requested
  !<
  subroutine npf_save_model_flows(this, flowja, icbcfl, icbcun)
    ! -- dummy
    class(GwfNpfType) :: this
    real(DP), dimension(:), intent(in) :: flowja
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: icbcun
    ! -- local
    integer(I4B) :: ibinun
    !
    ! -- Set unit number for binary output
    if (this%ipakcb < 0) then
      ibinun = icbcun
    elseif (this%ipakcb == 0) then
      ibinun = 0
    else
      ibinun = this%ipakcb
    end if
    if (icbcfl == 0) ibinun = 0
    !
    ! -- Write the face flows if requested
    if (ibinun /= 0) then
      call this%dis%record_connection_array(flowja, ibinun, this%iout)
    end if
    !
    ! -- Calculate specific discharge at cell centers and write, if requested
    if (this%isavspdis /= 0) then
      if (ibinun /= 0) call this%sav_spdis(ibinun)
    end if
    !
    ! -- Save saturation, if requested
    if (this%isavsat /= 0) then
      if (ibinun /= 0) call this%sav_sat(ibinun)
    end if
  end subroutine npf_save_model_flows

  !> @brief Print budget
  !<
  subroutine npf_print_model_flows(this, ibudfl, flowja)
    ! -- modules
    use TdisModule, only: kper, kstp
    use ConstantsModule, only: LENBIGLINE
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: ibudfl
    real(DP), intent(inout), dimension(:) :: flowja
    ! -- local
    character(len=LENBIGLINE) :: line
    character(len=30) :: tempstr
    integer(I4B) :: n, ipos, m
    real(DP) :: qnm
    ! -- formats
    character(len=*), parameter :: fmtiprflow = &
      &"(/,4x,'CALCULATED INTERCELL FLOW FOR PERIOD ', i0, ' STEP ', i0)"
    !
    ! -- Write flowja to list file if requested
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
  end subroutine npf_print_model_flows

  !> @brief Deallocate variables
  !<
  subroutine npf_da(this)
    ! -- modules
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    ! -- dummy
    class(GwfNpftype) :: this

    ! free spdis work structure
    if (this%icalcspdis == 1 .and. this%spdis_wa%is_created()) &
      call this%spdis_wa%destroy()
    deallocate (this%spdis_wa)
    !
    ! -- Deallocate input memory
    call memorystore_remove(this%name_model, 'NPF', idm_context)
    !
    ! -- TVK
    if (this%intvk /= 0) then
      call this%tvk%da()
      deallocate (this%tvk)
    end if
    !
    ! -- VSC
    if (this%invsc /= 0) then
      nullify (this%vsc)
    end if
    !
    ! -- Strings
    !
    ! -- Scalars
    call mem_deallocate(this%iname)
    call mem_deallocate(this%ixt3d)
    call mem_deallocate(this%ixt3drhs)
    call mem_deallocate(this%satomega)
    call mem_deallocate(this%hnoflo)
    call mem_deallocate(this%hdry)
    call mem_deallocate(this%icellavg)
    call mem_deallocate(this%iavgkeff)
    call mem_deallocate(this%ik22)
    call mem_deallocate(this%ik33)
    call mem_deallocate(this%iperched)
    call mem_deallocate(this%ivarcv)
    call mem_deallocate(this%idewatcv)
    call mem_deallocate(this%ithickstrt)
    call mem_deallocate(this%isavspdis)
    call mem_deallocate(this%isavsat)
    call mem_deallocate(this%icalcspdis)
    call mem_deallocate(this%irewet)
    call mem_deallocate(this%wetfct)
    call mem_deallocate(this%iwetit)
    call mem_deallocate(this%ihdwet)
    call mem_deallocate(this%ibotnode)
    call mem_deallocate(this%iwetdry)
    call mem_deallocate(this%iangle1)
    call mem_deallocate(this%iangle2)
    call mem_deallocate(this%iangle3)
    call mem_deallocate(this%nedges)
    call mem_deallocate(this%lastedge)
    call mem_deallocate(this%ik22overk)
    call mem_deallocate(this%ik33overk)
    call mem_deallocate(this%intvk)
    call mem_deallocate(this%invsc)
    call mem_deallocate(this%kchangeper)
    call mem_deallocate(this%kchangestp)
    !
    ! -- Deallocate arrays
    deallocate (this%aname)
    call mem_deallocate(this%ithickstartflag)
    call mem_deallocate(this%icelltype)
    call mem_deallocate(this%k11)
    call mem_deallocate(this%k22)
    call mem_deallocate(this%k33)
    call mem_deallocate(this%k11input)
    call mem_deallocate(this%k22input)
    call mem_deallocate(this%k33input)
    call mem_deallocate(this%sat, 'SAT', this%memoryPath)
    call mem_deallocate(this%condsat)
    call mem_deallocate(this%wetdry)
    call mem_deallocate(this%angle1)
    call mem_deallocate(this%angle2)
    call mem_deallocate(this%angle3)
    call mem_deallocate(this%nodedge)
    call mem_deallocate(this%ihcedge)
    call mem_deallocate(this%propsedge)
    call mem_deallocate(this%iedge_ptr)
    call mem_deallocate(this%edge_idxs)
    call mem_deallocate(this%spdis, 'SPDIS', this%memoryPath)
    call mem_deallocate(this%nodekchange)
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()

    ! pointers
    this%hnew => null()

  end subroutine npf_da

  !> @ brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the VSC package. The base model
  !! allocate scalars method is also called.
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryHelperModule, only: create_mem_path
    ! -- dummy
    class(GwfNpftype) :: this
    !
    ! -- allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate scalars
    call mem_allocate(this%iname, 'INAME', this%memoryPath)
    call mem_allocate(this%ixt3d, 'IXT3D', this%memoryPath)
    call mem_allocate(this%ixt3drhs, 'IXT3DRHS', this%memoryPath)
    call mem_allocate(this%satomega, 'SATOMEGA', this%memoryPath)
    call mem_allocate(this%hnoflo, 'HNOFLO', this%memoryPath)
    call mem_allocate(this%hdry, 'HDRY', this%memoryPath)
    call mem_allocate(this%icellavg, 'ICELLAVG', this%memoryPath)
    call mem_allocate(this%iavgkeff, 'IAVGKEFF', this%memoryPath)
    call mem_allocate(this%ik22, 'IK22', this%memoryPath)
    call mem_allocate(this%ik33, 'IK33', this%memoryPath)
    call mem_allocate(this%ik22overk, 'IK22OVERK', this%memoryPath)
    call mem_allocate(this%ik33overk, 'IK33OVERK', this%memoryPath)
    call mem_allocate(this%iperched, 'IPERCHED', this%memoryPath)
    call mem_allocate(this%ivarcv, 'IVARCV', this%memoryPath)
    call mem_allocate(this%idewatcv, 'IDEWATCV', this%memoryPath)
    call mem_allocate(this%ithickstrt, 'ITHICKSTRT', this%memoryPath)
    call mem_allocate(this%icalcspdis, 'ICALCSPDIS', this%memoryPath)
    call mem_allocate(this%isavspdis, 'ISAVSPDIS', this%memoryPath)
    call mem_allocate(this%isavsat, 'ISAVSAT', this%memoryPath)
    call mem_allocate(this%irewet, 'IREWET', this%memoryPath)
    call mem_allocate(this%wetfct, 'WETFCT', this%memoryPath)
    call mem_allocate(this%iwetit, 'IWETIT', this%memoryPath)
    call mem_allocate(this%ihdwet, 'IHDWET', this%memoryPath)
    call mem_allocate(this%iangle1, 'IANGLE1', this%memoryPath)
    call mem_allocate(this%iangle2, 'IANGLE2', this%memoryPath)
    call mem_allocate(this%iangle3, 'IANGLE3', this%memoryPath)
    call mem_allocate(this%iwetdry, 'IWETDRY', this%memoryPath)
    call mem_allocate(this%nedges, 'NEDGES', this%memoryPath)
    call mem_allocate(this%lastedge, 'LASTEDGE', this%memoryPath)
    call mem_allocate(this%intvk, 'INTVK', this%memoryPath)
    call mem_allocate(this%invsc, 'INVSC', this%memoryPath)
    call mem_allocate(this%kchangeper, 'KCHANGEPER', this%memoryPath)
    call mem_allocate(this%kchangestp, 'KCHANGESTP', this%memoryPath)
    !
    ! -- set pointer to inewtonur
    call mem_setptr(this%igwfnewtonur, 'INEWTONUR', &
                    create_mem_path(this%name_model))
    !
    ! -- Initialize value
    this%iname = 8
    this%ixt3d = 0
    this%ixt3drhs = 0
    this%satomega = DZERO
    this%hnoflo = DHNOFLO !1.d30
    this%hdry = DHDRY !-1.d30
    this%icellavg = CCOND_HMEAN
    this%iavgkeff = 0
    this%ik22 = 0
    this%ik33 = 0
    this%ik22overk = 0
    this%ik33overk = 0
    this%iperched = 0
    this%ivarcv = 0
    this%idewatcv = 0
    this%ithickstrt = 0
    this%icalcspdis = 0
    this%isavspdis = 0
    this%isavsat = 0
    this%irewet = 0
    this%wetfct = DONE
    this%iwetit = 1
    this%ihdwet = 0
    this%iangle1 = 0
    this%iangle2 = 0
    this%iangle3 = 0
    this%iwetdry = 0
    this%nedges = 0
    this%lastedge = 0
    this%intvk = 0
    this%invsc = 0
    this%kchangeper = 0
    this%kchangestp = 0
    !
    ! -- If newton is on, then NPF creates asymmetric matrix
    this%iasym = this%inewton
  end subroutine allocate_scalars

  !> @ brief Store backup copy of hydraulic conductivity when the VSC
  !!         package is activate
  !!
  !! The K arrays (K11, etc.) get multiplied by the viscosity ratio so that
  !! subsequent uses of K already take into account the effect of viscosity.
  !! Thus the original user-specified K array values are lost unless they are
  !! backed up in k11input, for example.  In a new stress period/time step,
  !! the values in k11input are multiplied by the viscosity ratio, not k11
  !! since it contains viscosity-adjusted hydraulic conductivity values.
  !<
  subroutine store_original_k_arrays(this, ncells, njas)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwfNpftype) :: this
    integer(I4B), intent(in) :: ncells
    integer(I4B), intent(in) :: njas
    ! -- local
    integer(I4B) :: n
    !
    ! -- Retain copy of user-specified K arrays
    do n = 1, ncells
      this%k11input(n) = this%k11(n)
      this%k22input(n) = this%k22(n)
      this%k33input(n) = this%k33(n)
    end do
  end subroutine store_original_k_arrays

  !> @brief Allocate npf arrays
  !<
  subroutine allocate_arrays(this, ncells, njas)
    ! -- dummy
    class(GwfNpftype) :: this
    integer(I4B), intent(in) :: ncells
    integer(I4B), intent(in) :: njas
    ! -- local
    integer(I4B) :: n
    !
    call mem_allocate(this%ithickstartflag, ncells, 'ITHICKSTARTFLAG', &
                      this%memoryPath)
    call mem_allocate(this%icelltype, ncells, 'ICELLTYPE', this%memoryPath)
    call mem_allocate(this%k11, ncells, 'K11', this%memoryPath)
    call mem_allocate(this%sat, ncells, 'SAT', this%memoryPath)
    call mem_allocate(this%condsat, njas, 'CONDSAT', this%memoryPath)
    !
    ! -- Optional arrays dimensioned to full size initially
    call mem_allocate(this%k22, ncells, 'K22', this%memoryPath)
    call mem_allocate(this%k33, ncells, 'K33', this%memoryPath)
    call mem_allocate(this%wetdry, ncells, 'WETDRY', this%memoryPath)
    call mem_allocate(this%angle1, ncells, 'ANGLE1', this%memoryPath)
    call mem_allocate(this%angle2, ncells, 'ANGLE2', this%memoryPath)
    call mem_allocate(this%angle3, ncells, 'ANGLE3', this%memoryPath)
    !
    ! -- Optional arrays
    call mem_allocate(this%ibotnode, 0, 'IBOTNODE', this%memoryPath)
    call mem_allocate(this%nodedge, 0, 'NODEDGE', this%memoryPath)
    call mem_allocate(this%ihcedge, 0, 'IHCEDGE', this%memoryPath)
    call mem_allocate(this%propsedge, 0, 0, 'PROPSEDGE', this%memoryPath)
    call mem_allocate(this%iedge_ptr, 0, 'NREDGESNODE', this%memoryPath)
    call mem_allocate(this%edge_idxs, 0, 'EDGEIDXS', this%memoryPath)
    !
    ! -- Optional arrays only needed when vsc package is active
    call mem_allocate(this%k11input, 0, 'K11INPUT', this%memoryPath)
    call mem_allocate(this%k22input, 0, 'K22INPUT', this%memoryPath)
    call mem_allocate(this%k33input, 0, 'K33INPUT', this%memoryPath)
    !
    ! -- Specific discharge is (re-)allocated when nedges is known
    call mem_allocate(this%spdis, 3, 0, 'SPDIS', this%memoryPath)
    !
    ! -- Time-varying property flag arrays
    call mem_allocate(this%nodekchange, ncells, 'NODEKCHANGE', this%memoryPath)
    !
    ! -- initialize iangle1, iangle2, iangle3, and wetdry
    do n = 1, ncells
      this%angle1(n) = DZERO
      this%angle2(n) = DZERO
      this%angle3(n) = DZERO
      this%wetdry(n) = DZERO
      this%nodekchange(n) = DZERO
    end do
    !
    ! -- allocate variable names
    allocate (this%aname(this%iname))
    this%aname = ['               ICELLTYPE', '                       K', &
                  '                     K33', '                     K22', &
                  '                  WETDRY', '                  ANGLE1', &
                  '                  ANGLE2', '                  ANGLE3']
  end subroutine allocate_arrays

  !> @brief Log npf options sourced from the input mempath
  !<
  subroutine log_options(this, found)
    ! -- modules
    use KindModule, only: LGP
    use GwfNpfInputModule, only: GwfNpfParamFoundType
    ! -- dummy
    class(GwfNpftype) :: this
    ! -- locals
    type(GwfNpfParamFoundType), intent(in) :: found
    !
    write (this%iout, '(1x,a)') 'Setting NPF Options'
    if (found%iprflow) &
      write (this%iout, '(4x,a)') 'Cell-by-cell flow information will be printed &
                                  &to listing file whenever ICBCFL is not zero.'
    if (found%ipakcb) &
      write (this%iout, '(4x,a)') 'Cell-by-cell flow information will be saved &
                                  &to binary file whenever ICBCFL is not zero.'
    if (found%cellavg) &
      write (this%iout, '(4x,a,i0)') 'Alternative cell averaging [1=logarithmic, &
                                     &2=AMT-LMK, 3=AMT-HMK] set to: ', &
                                     this%icellavg
    if (found%ithickstrt) &
      write (this%iout, '(4x,a)') 'THICKSTRT option has been activated.'
    if (found%iperched) &
      write (this%iout, '(4x,a)') 'Vertical flow will be adjusted for perched &
                                  &conditions.'
    if (found%ivarcv) &
      write (this%iout, '(4x,a)') 'Vertical conductance varies with water table.'
    if (found%idewatcv) &
      write (this%iout, '(4x,a)') 'Vertical conductance accounts for dewatered &
                                  &portion of an underlying cell.'
    if (found%ixt3d) write (this%iout, '(4x,a)') 'XT3D formulation is selected.'
    if (found%ixt3drhs) &
      write (this%iout, '(4x,a)') 'XT3D RHS formulation is selected.'
    if (found%isavspdis) &
      write (this%iout, '(4x,a)') 'Specific discharge will be calculated at cell &
                                  &centers and written to DATA-SPDIS in budget &
                                  &file when requested.'
    if (found%isavsat) &
      write (this%iout, '(4x,a)') 'Saturation will be written to DATA-SAT in &
                                  &budget file when requested.'
    if (found%ik22overk) &
      write (this%iout, '(4x,a)') 'Values specified for K22 are anisotropy &
                                  &ratios and will be multiplied by K before &
                                  &being used in calculations.'
    if (found%ik33overk) &
      write (this%iout, '(4x,a)') 'Values specified for K33 are anisotropy &
                                  &ratios and will be multiplied by K before &
                                  &being used in calculations.'
    if (found%inewton) &
      write (this%iout, '(4x,a)') 'NEWTON-RAPHSON method disabled for unconfined &
                                  &cells'
    if (found%satomega) &
      write (this%iout, '(4x,a,1pg15.6)') 'Saturation omega: ', this%satomega
    if (found%irewet) &
      write (this%iout, '(4x,a)') 'Rewetting is active.'
    if (found%wetfct) &
      write (this%iout, '(4x,a,1pg15.6)') &
      'Wetting factor (WETFCT) has been set to: ', this%wetfct
    if (found%iwetit) &
      write (this%iout, '(4x,a,i5)') &
      'Wetting iteration interval (IWETIT) has been set to: ', this%iwetit
    if (found%ihdwet) &
      write (this%iout, '(4x,a,i5)') &
      'Head rewet equation (IHDWET) has been set to: ', this%ihdwet
    write (this%iout, '(1x,a,/)') 'End Setting NPF Options'
  end subroutine log_options

  !> @brief Update simulation options from input mempath
  !<
  subroutine source_options(this)
    ! -- modules
    use SimModule, only: store_error, store_error_filename
    use MemoryManagerModule, only: mem_setptr, get_isize
    use MemoryManagerExtModule, only: mem_set_value
    use CharacterStringModule, only: CharacterStringType
    use GwfNpfInputModule, only: GwfNpfParamFoundType
    use SourceCommonModule, only: filein_fname
    ! -- dummy
    class(GwfNpftype) :: this
    ! -- locals
    character(len=LENVARNAME), dimension(3) :: cellavg_method = &
      &[character(len=LENVARNAME) :: 'LOGARITHMIC', 'AMT-LMK', 'AMT-HMK']
    type(GwfNpfParamFoundType) :: found
    character(len=LINELENGTH) :: tvk6_filename
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%iprflow, 'IPRFLOW', this%input_mempath, found%iprflow)
    call mem_set_value(this%ipakcb, 'IPAKCB', this%input_mempath, found%ipakcb)
    call mem_set_value(this%icellavg, 'CELLAVG', this%input_mempath, &
                       cellavg_method, found%cellavg)
    call mem_set_value(this%ithickstrt, 'ITHICKSTRT', this%input_mempath, &
                       found%ithickstrt)
    call mem_set_value(this%iperched, 'IPERCHED', this%input_mempath, &
                       found%iperched)
    call mem_set_value(this%ivarcv, 'IVARCV', this%input_mempath, found%ivarcv)
    call mem_set_value(this%idewatcv, 'IDEWATCV', this%input_mempath, &
                       found%idewatcv)
    call mem_set_value(this%ixt3d, 'IXT3D', this%input_mempath, found%ixt3d)
    call mem_set_value(this%ixt3drhs, 'IXT3DRHS', this%input_mempath, &
                       found%ixt3drhs)
    call mem_set_value(this%isavspdis, 'ISAVSPDIS', this%input_mempath, &
                       found%isavspdis)
    call mem_set_value(this%isavsat, 'ISAVSAT', this%input_mempath, found%isavsat)
    call mem_set_value(this%ik22overk, 'IK22OVERK', this%input_mempath, &
                       found%ik22overk)
    call mem_set_value(this%ik33overk, 'IK33OVERK', this%input_mempath, &
                       found%ik33overk)
    call mem_set_value(this%inewton, 'INEWTON', this%input_mempath, found%inewton)
    call mem_set_value(this%satomega, 'SATOMEGA', this%input_mempath, &
                       found%satomega)
    call mem_set_value(this%irewet, 'IREWET', this%input_mempath, found%irewet)
    call mem_set_value(this%wetfct, 'WETFCT', this%input_mempath, found%wetfct)
    call mem_set_value(this%iwetit, 'IWETIT', this%input_mempath, found%iwetit)
    call mem_set_value(this%ihdwet, 'IHDWET', this%input_mempath, found%ihdwet)
    !
    ! -- save flows option active
    if (found%ipakcb) this%ipakcb = -1
    !
    ! -- xt3d active with rhs
    if (found%ixt3d .and. found%ixt3drhs) this%ixt3d = 2
    !
    ! -- save specific discharge active
    if (found%isavspdis) this%icalcspdis = this%isavspdis
    !
    ! -- no newton specified
    if (found%inewton) then
      this%inewton = 0
      this%iasym = 0
    end if
    !
    ! -- enforce 0 or 1 TVK6_FILENAME entries in option block
    if (filein_fname(tvk6_filename, 'TVK6_FILENAME', this%input_mempath, &
                     this%input_fname)) then
      call openfile(this%intvk, this%iout, tvk6_filename, 'TVK')
      call tvk_cr(this%tvk, this%name_model, this%intvk, this%iout)
    end if
    !
    ! -- log options
    if (this%iout > 0) then
      call this%log_options(found)
    end if
  end subroutine source_options

  !> @brief Set options in the NPF object
  !<
  subroutine set_options(this, options)
    ! -- dummy
    class(GwfNpftype) :: this
    type(GwfNpfOptionsType), intent(in) :: options
    !
    this%icellavg = options%icellavg
    this%ithickstrt = options%ithickstrt
    this%iperched = options%iperched
    this%ivarcv = options%ivarcv
    this%idewatcv = options%idewatcv
    this%irewet = options%irewet
    this%wetfct = options%wetfct
    this%iwetit = options%iwetit
    this%ihdwet = options%ihdwet
  end subroutine set_options

  !> @brief Check for conflicting NPF options
  !<
  subroutine check_options(this)
    ! -- modules
    use SimModule, only: store_error, count_errors, store_error_filename
    use ConstantsModule, only: LINELENGTH
    ! -- dummy
    class(GwfNpftype) :: this
    !
    ! -- set omega value used for saturation calculations
    if (this%inewton > 0) then
      this%satomega = DEM6
    end if
    !
    if (this%inewton > 0) then
      if (this%iperched > 0) then
        write (errmsg, '(a)') 'ERROR IN NPF OPTIONS. NEWTON OPTION CANNOT '// &
          'BE USED WITH PERCHED OPTION.'
        call store_error(errmsg)
      end if
      if (this%ivarcv > 0) then
        write (errmsg, '(a)') 'ERROR IN NPF OPTIONS. NEWTON OPTION CANNOT '// &
          'BE USED WITH VARIABLECV OPTION.'
        call store_error(errmsg)
      end if
      if (this%irewet > 0) then
        write (errmsg, '(a)') 'ERROR IN NPF OPTIONS. NEWTON OPTION CANNOT '// &
          'BE USED WITH REWET OPTION.'
        call store_error(errmsg)
      end if
    end if
    !
    if (this%ixt3d /= 0) then
      if (this%icellavg > 0) then
        write (errmsg, '(a)') 'ERROR IN NPF OPTIONS. '// &
          'ALTERNATIVE_CELL_AVERAGING OPTION '// &
          'CANNOT BE USED WITH XT3D OPTION.'
        call store_error(errmsg)
      end if
      if (this%ithickstrt > 0) then
        write (errmsg, '(a)') 'ERROR IN NPF OPTIONS. THICKSTRT OPTION '// &
          'CANNOT BE USED WITH XT3D OPTION.'
        call store_error(errmsg)
      end if
      if (this%iperched > 0) then
        write (errmsg, '(a)') 'ERROR IN NPF OPTIONS. PERCHED OPTION '// &
          'CANNOT BE USED WITH XT3D OPTION.'
        call store_error(errmsg)
      end if
      if (this%ivarcv > 0) then
        write (errmsg, '(a)') 'ERROR IN NPF OPTIONS. VARIABLECV OPTION '// &
          'CANNOT BE USED WITH XT3D OPTION.'
        call store_error(errmsg)
      end if
    end if
    !
    ! -- Terminate if errors
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine check_options

  !> @brief Write dimensions to list file
  !<
  subroutine log_griddata(this, found)
    ! -- modules
    use GwfNpfInputModule, only: GwfNpfParamFoundType
    ! -- dummy
    class(GwfNpfType) :: this
    type(GwfNpfParamFoundType), intent(in) :: found
    !
    write (this%iout, '(1x,a)') 'Setting NPF Griddata'
    !
    if (found%icelltype) then
      write (this%iout, '(4x,a)') 'ICELLTYPE set from input file'
    end if
    !
    if (found%k) then
      write (this%iout, '(4x,a)') 'K set from input file'
    end if
    !
    if (found%k33) then
      write (this%iout, '(4x,a)') 'K33 set from input file'
    else
      write (this%iout, '(4x,a)') 'K33 not provided.  Setting K33 = K.'
    end if
    !
    if (found%k22) then
      write (this%iout, '(4x,a)') 'K22 set from input file'
    else
      write (this%iout, '(4x,a)') 'K22 not provided.  Setting K22 = K.'
    end if
    !
    if (found%wetdry) then
      write (this%iout, '(4x,a)') 'WETDRY set from input file'
    end if
    !
    if (found%angle1) then
      write (this%iout, '(4x,a)') 'ANGLE1 set from input file'
    end if
    !
    if (found%angle2) then
      write (this%iout, '(4x,a)') 'ANGLE2 set from input file'
    end if
    !
    if (found%angle3) then
      write (this%iout, '(4x,a)') 'ANGLE3 set from input file'
    end if
    !
    write (this%iout, '(1x,a,/)') 'End Setting NPF Griddata'
  end subroutine log_griddata

  !> @brief Update simulation griddata from input mempath
  !<
  subroutine source_griddata(this)
    ! -- modules
    use SimModule, only: count_errors, store_error
    use MemoryManagerModule, only: mem_reallocate
    use MemoryManagerExtModule, only: mem_set_value
    use GwfNpfInputModule, only: GwfNpfParamFoundType
    ! -- dummy
    class(GwfNpftype) :: this
    ! -- locals
    character(len=LINELENGTH) :: errmsg
    type(GwfNpfParamFoundType) :: found
    logical, dimension(2) :: afound
    integer(I4B), dimension(:), pointer, contiguous :: map
    !
    ! -- set map to convert user input data into reduced data
    map => null()
    if (this%dis%nodes < this%dis%nodesuser) map => this%dis%nodeuser
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%icelltype, 'ICELLTYPE', this%input_mempath, map, &
                       found%icelltype)
    call mem_set_value(this%k11, 'K', this%input_mempath, map, found%k)
    call mem_set_value(this%k33, 'K33', this%input_mempath, map, found%k33)
    call mem_set_value(this%k22, 'K22', this%input_mempath, map, found%k22)
    call mem_set_value(this%wetdry, 'WETDRY', this%input_mempath, map, &
                       found%wetdry)
    call mem_set_value(this%angle1, 'ANGLE1', this%input_mempath, map, &
                       found%angle1)
    call mem_set_value(this%angle2, 'ANGLE2', this%input_mempath, map, &
                       found%angle2)
    call mem_set_value(this%angle3, 'ANGLE3', this%input_mempath, map, &
                       found%angle3)
    !
    ! -- ensure ICELLTYPE was found
    if (.not. found%icelltype) then
      write (errmsg, '(a)') 'Error in GRIDDATA block: ICELLTYPE not found.'
      call store_error(errmsg)
    end if
    !
    ! -- ensure K was found
    if (.not. found%k) then
      write (errmsg, '(a)') 'Error in GRIDDATA block: K not found.'
      call store_error(errmsg)
    end if
    !
    ! -- set error if ik33overk set with no k33
    if (.not. found%k33 .and. this%ik33overk /= 0) then
      write (errmsg, '(a)') 'K33OVERK option specified but K33 not specified.'
      call store_error(errmsg)
    end if
    !
    ! -- set error if ik22overk set with no k22
    if (.not. found%k22 .and. this%ik22overk /= 0) then
      write (errmsg, '(a)') 'K22OVERK option specified but K22 not specified.'
      call store_error(errmsg)
    end if
    !
    ! -- handle found side effects
    if (found%k33) this%ik33 = 1
    if (found%k22) this%ik22 = 1
    if (found%wetdry) this%iwetdry = 1
    if (found%angle1) this%iangle1 = 1
    if (found%angle2) this%iangle2 = 1
    if (found%angle3) this%iangle3 = 1
    !
    ! -- handle not found side effects
    if (.not. found%k33) then
      call mem_set_value(this%k33, 'K', this%input_mempath, map, afound(1))
    end if
    if (.not. found%k22) then
      call mem_set_value(this%k22, 'K', this%input_mempath, map, afound(2))
    end if
    if (.not. found%wetdry) call mem_reallocate(this%wetdry, 1, 'WETDRY', &
                                                trim(this%memoryPath))
    if (.not. found%angle1 .and. this%ixt3d == 0) &
      call mem_reallocate(this%angle1, 0, 'ANGLE1', trim(this%memoryPath))
    if (.not. found%angle2 .and. this%ixt3d == 0) &
      call mem_reallocate(this%angle2, 0, 'ANGLE2', trim(this%memoryPath))
    if (.not. found%angle3 .and. this%ixt3d == 0) &
      call mem_reallocate(this%angle3, 0, 'ANGLE3', trim(this%memoryPath))
    !
    ! -- log griddata
    if (this%iout > 0) then
      call this%log_griddata(found)
    end if
  end subroutine source_griddata

  !> @brief Initialize and check NPF data
  !<
  subroutine prepcheck(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH, DPIO180
    use SimModule, only: store_error, count_errors, store_error_filename
    ! -- dummy
    class(GwfNpfType) :: this
    ! -- local
    character(len=24), dimension(:), pointer :: aname
    character(len=LINELENGTH) :: cellstr, errmsg
    integer(I4B) :: nerr, n
    ! -- format
    character(len=*), parameter :: fmtkerr = &
      &"(1x, 'Hydraulic property ',a,' is <= 0 for cell ',a, ' ', 1pg15.6)"
    character(len=*), parameter :: fmtkerr2 = &
      &"(1x, '... ', i0,' additional errors not shown for ',a)"
    !
    ! -- initialize
    aname => this%aname
    !
    ! -- check k11
    nerr = 0
    do n = 1, size(this%k11)
      if (this%k11(n) <= DZERO) then
        nerr = nerr + 1
        if (nerr <= 20) then
          call this%dis%noder_to_string(n, cellstr)
          write (errmsg, fmtkerr) trim(adjustl(aname(2))), trim(cellstr), &
            this%k11(n)
          call store_error(errmsg)
        end if
      end if
    end do
    if (nerr > 20) then
      write (errmsg, fmtkerr2) nerr, trim(adjustl(aname(2)))
      call store_error(errmsg)
    end if
    !
    ! -- check k33 because it was read
    if (this%ik33 /= 0) then
      !
      ! -- Check to make sure values are greater than or equal to zero
      nerr = 0
      do n = 1, size(this%k33)
        if (this%ik33overk /= 0) this%k33(n) = this%k33(n) * this%k11(n)
        if (this%k33(n) <= DZERO) then
          nerr = nerr + 1
          if (nerr <= 20) then
            call this%dis%noder_to_string(n, cellstr)
            write (errmsg, fmtkerr) trim(adjustl(aname(3))), trim(cellstr), &
              this%k33(n)
            call store_error(errmsg)
          end if
        end if
      end do
      if (nerr > 20) then
        write (errmsg, fmtkerr2) nerr, trim(adjustl(aname(3)))
        call store_error(errmsg)
      end if
    end if
    !
    ! -- check k22 because it was read
    if (this%ik22 /= 0) then
      !
      ! -- Check to make sure that angles are available
      if (this%dis%con%ianglex == 0) then
        write (errmsg, '(a)') 'Error.  ANGLDEGX not provided in '// &
          'discretization file, but K22 was specified. '
        call store_error(errmsg)
      end if
      !
      ! -- Check to make sure values are greater than or equal to zero
      nerr = 0
      do n = 1, size(this%k22)
        if (this%ik22overk /= 0) this%k22(n) = this%k22(n) * this%k11(n)
        if (this%k22(n) <= DZERO) then
          nerr = nerr + 1
          if (nerr <= 20) then
            call this%dis%noder_to_string(n, cellstr)
            write (errmsg, fmtkerr) trim(adjustl(aname(4))), trim(cellstr), &
              this%k22(n)
            call store_error(errmsg)
          end if
        end if
      end do
      if (nerr > 20) then
        write (errmsg, fmtkerr2) nerr, trim(adjustl(aname(4)))
        call store_error(errmsg)
      end if
    end if
    !
    ! -- check for wetdry conflicts
    if (this%irewet == 1) then
      if (this%iwetdry == 0) then
        write (errmsg, '(a, a, a)') 'Error in GRIDDATA block: ', &
          trim(adjustl(aname(5))), ' not found.'
        call store_error(errmsg)
      end if
    end if
    !
    ! -- Check for angle conflicts
    if (this%iangle1 /= 0) then
      do n = 1, size(this%angle1)
        this%angle1(n) = this%angle1(n) * DPIO180
      end do
    else
      if (this%ixt3d /= 0) then
        this%iangle1 = 1
        write (this%iout, '(a)') 'XT3D IN USE, BUT ANGLE1 NOT SPECIFIED. '// &
          'SETTING ANGLE1 TO ZERO.'
        do n = 1, size(this%angle1)
          this%angle1(n) = DZERO
        end do
      end if
    end if
    if (this%iangle2 /= 0) then
      if (this%iangle1 == 0) then
        write (errmsg, '(a)') 'ANGLE2 SPECIFIED BUT NOT ANGLE1. '// &
          'ANGLE2 REQUIRES ANGLE1. '
        call store_error(errmsg)
      end if
      if (this%iangle3 == 0) then
        write (errmsg, '(a)') 'ANGLE2 SPECIFIED BUT NOT ANGLE3. '// &
          'SPECIFY BOTH OR NEITHER ONE. '
        call store_error(errmsg)
      end if
      do n = 1, size(this%angle2)
        this%angle2(n) = this%angle2(n) * DPIO180
      end do
    end if
    if (this%iangle3 /= 0) then
      if (this%iangle1 == 0) then
        write (errmsg, '(a)') 'ANGLE3 SPECIFIED BUT NOT ANGLE1. '// &
          'ANGLE3 REQUIRES ANGLE1. '
        call store_error(errmsg)
      end if
      if (this%iangle2 == 0) then
        write (errmsg, '(a)') 'ANGLE3 SPECIFIED BUT NOT ANGLE2. '// &
          'SPECIFY BOTH OR NEITHER ONE. '
        call store_error(errmsg)
      end if
      do n = 1, size(this%angle3)
        this%angle3(n) = this%angle3(n) * DPIO180
      end do
    end if
    !
    ! -- terminate if data errors
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine prepcheck

  !> @brief preprocess the NPF input data
  !!
  !! This routine consists of the following steps:
  !!
  !! 1. convert cells to noflow when all transmissive parameters equal zero
  !! 2. perform initial wetting and drying
  !! 3. initialize cell saturation
  !! 4. calculate saturated conductance (when not xt3d)
  !! 5. If NEWTON under-relaxation, determine lower most node
  !<
  subroutine preprocess_input(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error, count_errors, store_error_filename
    ! -- dummy
    class(GwfNpfType) :: this !< the instance of the NPF package
    ! -- local
    integer(I4B) :: n, m, ii, nn
    real(DP) :: hyn, hym
    real(DP) :: satn, topn, botn
    integer(I4B) :: nextn
    real(DP) :: minbot, botm
    logical :: finished
    character(len=LINELENGTH) :: cellstr, errmsg
    ! -- format
    character(len=*), parameter :: fmtcnv = &
      "(1X,'CELL ', A, &
      &' ELIMINATED BECAUSE ALL HYDRAULIC CONDUCTIVITIES TO NODE ARE 0.')"
    character(len=*), parameter :: fmtnct = &
      &"(1X,'Negative cell thickness at cell ', A)"
    character(len=*), parameter :: fmtihbe = &
      &"(1X,'Initial head, bottom elevation:',1P,2G13.5)"
    character(len=*), parameter :: fmttebe = &
      &"(1X,'Top elevation, bottom elevation:',1P,2G13.5)"
    !
    do n = 1, this%dis%nodes
      this%ithickstartflag(n) = 0
    end do
    !
    ! -- Insure that each cell has at least one non-zero transmissive parameter
    !    Note that a cell can be deactivated even if it has a valid connection
    !    to another model.
    nodeloop: do n = 1, this%dis%nodes
      !
      ! -- Skip if already inactive
      if (this%ibound(n) == 0) then
        if (this%irewet /= 0) then
          if (this%wetdry(n) == DZERO) cycle nodeloop
        else
          cycle nodeloop
        end if
      end if
      !
      ! -- Cycle if k11 is not zero
      if (this%k11(n) /= DZERO) cycle nodeloop
      !
      ! -- Cycle if at least one vertical connection has non-zero k33
      !    for n and m
      do ii = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ii)
        if (this%dis%con%ihc(this%dis%con%jas(ii)) == 0) then
          hyn = this%k11(n)
          if (this%ik33 /= 0) hyn = this%k33(n)
          if (hyn /= DZERO) then
            hym = this%k11(m)
            if (this%ik33 /= 0) hym = this%k33(m)
            if (hym /= DZERO) cycle
          end if
        end if
      end do
      !
      ! -- If this part of the loop is reached, then all connections have
      !    zero transmissivity, so convert to noflow.
      this%ibound(n) = 0
      this%hnew(n) = this%hnoflo
      if (this%irewet /= 0) this%wetdry(n) = DZERO
      call this%dis%noder_to_string(n, cellstr)
      write (this%iout, fmtcnv) trim(adjustl(cellstr))
      !
    end do nodeloop
    !
    ! -- Preprocess cell status and heads based on initial conditions
    if (this%inewton == 0) then
      !
      ! -- For standard formulation (non-Newton) call wetdry routine
      call this%wd(0, this%hnew)
    else
      !
      ! -- Newton formulation, so adjust heads to be above bottom
      !    (Not used in present formulation because variable cv
      !    cannot be used with Newton)
      if (this%ivarcv == 1) then
        do n = 1, this%dis%nodes
          if (this%hnew(n) < this%dis%bot(n)) then
            this%hnew(n) = this%dis%bot(n) + DEM6
          end if
        end do
      end if
    end if
    !
    ! -- If THCKSTRT is not active, then loop through icelltype and replace
    !    any negative values with 1.
    if (this%ithickstrt == 0) then
      do n = 1, this%dis%nodes
        if (this%icelltype(n) < 0) then
          this%icelltype(n) = 1
        end if
      end do
    end if
    !
    ! -- Initialize sat to zero for ibound=0 cells, unless the cell can
    !    rewet.  Initialize sat to the saturated fraction based on strt
    !    if icelltype is negative and the THCKSTRT option is in effect.
    !    Initialize sat to 1.0 for all other cells in order to calculate
    !    condsat in next section.
    do n = 1, this%dis%nodes
      if (this%ibound(n) == 0) then
        this%sat(n) = DONE
        if (this%icelltype(n) < 0 .and. this%ithickstrt /= 0) then
          this%ithickstartflag(n) = 1
          this%icelltype(n) = 0
        end if
      else
        topn = this%dis%top(n)
        botn = this%dis%bot(n)
        if (this%icelltype(n) < 0 .and. this%ithickstrt /= 0) then
          call this%thksat(n, this%ic%strt(n), satn)
          if (botn > this%ic%strt(n)) then
            call this%dis%noder_to_string(n, cellstr)
            write (errmsg, fmtnct) trim(adjustl(cellstr))
            call store_error(errmsg)
            write (errmsg, fmtihbe) this%ic%strt(n), botn
            call store_error(errmsg)
          end if
          this%ithickstartflag(n) = 1
          this%icelltype(n) = 0
        else
          satn = DONE
          if (botn > topn) then
            call this%dis%noder_to_string(n, cellstr)
            write (errmsg, fmtnct) trim(adjustl(cellstr))
            call store_error(errmsg)
            write (errmsg, fmttebe) topn, botn
            call store_error(errmsg)
          end if
        end if
        this%sat(n) = satn
      end if
    end do
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
    !
    ! -- Calculate condsat, but only if xt3d is not active.  If xt3d is
    !    active, then condsat is allocated to size of zero.
    if (this%ixt3d == 0) then
      !
      ! -- Calculate the saturated conductance for all connections assuming
      !    that saturation is 1 (except for case where icelltype was entered
      !    as a negative value and THCKSTRT option in effect)
      do n = 1, this%dis%nodes
        call this%calc_condsat(n, .true.)
      end do
      !
    end if
    !
    ! -- Determine the lower most node
    if (this%igwfnewtonur /= 0) then
      call mem_reallocate(this%ibotnode, this%dis%nodes, 'IBOTNODE', &
                          trim(this%memoryPath))
      do n = 1, this%dis%nodes
        !
        minbot = this%dis%bot(n)
        nn = n
        finished = .false.
        do while (.not. finished)
          nextn = 0
          !
          ! -- Go through the connecting cells
          do ii = this%dis%con%ia(nn) + 1, this%dis%con%ia(nn + 1) - 1
            !
            ! -- Set the m cell number
            m = this%dis%con%ja(ii)
            botm = this%dis%bot(m)
            !
            ! -- select vertical connections: ihc == 0
            if (this%dis%con%ihc(this%dis%con%jas(ii)) == 0) then
              if (m > nn .and. botm < minbot) then
                nextn = m
                minbot = botm
              end if
            end if
          end do
          if (nextn > 0) then
            nn = nextn
          else
            finished = .true.
          end if
        end do
        this%ibotnode(n) = nn
      end do
    end if
    !
    ! -- nullify unneeded gwf pointers
    this%igwfnewtonur => null()
  end subroutine preprocess_input

  !> @brief Calculate CONDSAT array entries for the given node
  !!
  !! Calculate saturated conductances for all connections of the given node,
  !! or optionally for the upper portion of the matrix only.
  !<
  subroutine calc_condsat(this, node, upperOnly)
    ! -- dummy variables
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: node
    logical, intent(in) :: upperOnly
    ! -- local variables
    integer(I4B) :: ii, m, n, ihc, jj
    real(DP) :: topm, topn, topnode, botm, botn, botnode, satm, satn, satnode
    real(DP) :: hyn, hym, hn, hm, fawidth, csat
    !
    satnode = this%calc_initial_sat(node)
    !
    topnode = this%dis%top(node)
    botnode = this%dis%bot(node)
    !
    ! -- Go through the connecting cells
    do ii = this%dis%con%ia(node) + 1, this%dis%con%ia(node + 1) - 1
      !
      ! -- Set the m cell number and cycle if lower triangle connection and
      ! -- we're not updating both upper and lower matrix parts for this node
      m = this%dis%con%ja(ii)
      jj = this%dis%con%jas(ii)
      if (m < node) then
        if (upperOnly) cycle
        ! m => node, n => neighbour
        n = m
        m = node
        topm = topnode
        botm = botnode
        satm = satnode
        topn = this%dis%top(n)
        botn = this%dis%bot(n)
        satn = this%calc_initial_sat(n)
      else
        ! n => node, m => neighbour
        n = node
        topn = topnode
        botn = botnode
        satn = satnode
        topm = this%dis%top(m)
        botm = this%dis%bot(m)
        satm = this%calc_initial_sat(m)
      end if
      !
      ihc = this%dis%con%ihc(jj)
      hyn = this%hy_eff(n, m, ihc, ipos=ii)
      hym = this%hy_eff(m, n, ihc, ipos=ii)
      if (this%ithickstartflag(n) == 0) then
        hn = topn
      else
        hn = this%ic%strt(n)
      end if
      if (this%ithickstartflag(m) == 0) then
        hm = topm
      else
        hm = this%ic%strt(m)
      end if
      !
      ! -- Calculate conductance depending on whether connection is
      !    vertical (0), horizontal (1), or staggered horizontal (2)
      if (ihc == C3D_VERTICAL) then
        !
        ! -- Vertical conductance for fully saturated conditions
        csat = vcond(1, 1, 1, 1, 0, 1, 1, DONE, &
                     botn, botm, &
                     hyn, hym, &
                     satn, satm, &
                     topn, topm, &
                     botn, botm, &
                     this%dis%con%hwva(jj))
      else
        !
        ! -- Horizontal conductance for fully saturated conditions
        fawidth = this%dis%con%hwva(jj)
        csat = hcond(1, 1, 1, 1, 0, &
                     ihc, &
                     this%icellavg, &
                     DONE, &
                     hn, hm, satn, satm, hyn, hym, &
                     topn, topm, &
                     botn, botm, &
                     this%dis%con%cl1(jj), &
                     this%dis%con%cl2(jj), &
                     fawidth)
      end if
      this%condsat(jj) = csat
    end do
  end subroutine calc_condsat

  !> @brief Calculate initial saturation for the given node
  !!
  !! Calculate saturation as a fraction of thickness for the given node, used
  !! for saturated conductance calculations: full thickness by default (1.0) or
  !! saturation based on initial conditions if the THICKSTRT option is used.
  !!
  !<
  function calc_initial_sat(this, n) result(satn)
    ! -- dummy variables
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: n
    ! -- Return
    real(DP) :: satn
    !
    satn = DONE
    if (this%ibound(n) /= 0 .and. this%ithickstartflag(n) /= 0) then
      call this%thksat(n, this%ic%strt(n), satn)
    end if
  end function calc_initial_sat

  !> @brief Perform wetting and drying
  !<
  subroutine sgwf_npf_wetdry(this, kiter, hnew)
    ! -- modules
    use TdisModule, only: kstp, kper
    use SimModule, only: store_error, store_error_filename
    use ConstantsModule, only: LINELENGTH
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: kiter
    real(DP), intent(inout), dimension(:) :: hnew
    ! -- local
    integer(I4B) :: n, m, ii, ihc
    real(DP) :: ttop, bbot, thick
    integer(I4B) :: ncnvrt, ihdcnv
    character(len=30), dimension(5) :: nodcnvrt
    character(len=30) :: nodestr
    character(len=3), dimension(5) :: acnvrt
    character(len=LINELENGTH) :: errmsg
    integer(I4B) :: irewet
    ! -- formats
    character(len=*), parameter :: fmtnct = &
      "(1X,/1X,'Negative cell thickness at (layer,row,col)', &
      &I4,',',I5,',',I5)"
    character(len=*), parameter :: fmttopbot = &
      &"(1X,'Top elevation, bottom elevation:',1P,2G13.5)"
    character(len=*), parameter :: fmttopbotthk = &
      &"(1X,'Top elevation, bottom elevation, thickness:',1P,3G13.5)"
    character(len=*), parameter :: fmtdrychd = &
      &"(1X,/1X,'CONSTANT-HEAD CELL WENT DRY -- SIMULATION ABORTED')"
    character(len=*), parameter :: fmtni = &
      &"(1X,'CELLID=',a,' ITERATION=',I0,' TIME STEP=',I0,' STRESS PERIOD=',I0)"
    !
    ! -- Initialize
    ncnvrt = 0
    ihdcnv = 0
    !
    ! -- Convert dry cells to wet
    do n = 1, this%dis%nodes
      do ii = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ii)
        ihc = this%dis%con%ihc(this%dis%con%jas(ii))
        call this%rewet_check(kiter, n, hnew(m), this%ibound(m), ihc, hnew, &
                              irewet)
        if (irewet == 1) then
          call this%wdmsg(2, ncnvrt, nodcnvrt, acnvrt, ihdcnv, kiter, n)
        end if
      end do
    end do
    !
    ! -- Perform drying
    do n = 1, this%dis%nodes
      !
      ! -- cycle if inactive or confined
      if (this%ibound(n) == 0) cycle
      if (this%icelltype(n) == 0) cycle
      !
      ! -- check for negative cell thickness
      bbot = this%dis%bot(n)
      ttop = this%dis%top(n)
      if (bbot > ttop) then
        write (errmsg, fmtnct) n
        call store_error(errmsg)
        write (errmsg, fmttopbot) ttop, bbot
        call store_error(errmsg)
        call store_error_filename(this%input_fname)
      end if
      !
      ! -- Calculate saturated thickness
      if (this%icelltype(n) /= 0) then
        if (hnew(n) < ttop) ttop = hnew(n)
      end if
      thick = ttop - bbot
      !
      ! -- If thick<0 print message, set hnew, and ibound
      if (thick <= DZERO) then
        call this%wdmsg(1, ncnvrt, nodcnvrt, acnvrt, ihdcnv, kiter, n)
        hnew(n) = this%hdry
        if (this%ibound(n) < 0) then
          errmsg = 'CONSTANT-HEAD CELL WENT DRY -- SIMULATION ABORTED'
          call store_error(errmsg)
          write (errmsg, fmttopbotthk) ttop, bbot, thick
          call store_error(errmsg)
          call this%dis%noder_to_string(n, nodestr)
          write (errmsg, fmtni) trim(adjustl(nodestr)), kiter, kstp, kper
          call store_error(errmsg)
          call store_error_filename(this%input_fname)
        end if
        this%ibound(n) = 0
      end if
    end do
    !
    ! -- Print remaining cell conversions
    call this%wdmsg(0, ncnvrt, nodcnvrt, acnvrt, ihdcnv, kiter, n)
    !
    ! -- Change ibound from 30000 to 1
    do n = 1, this%dis%nodes
      if (this%ibound(n) == 30000) this%ibound(n) = 1
    end do
  end subroutine sgwf_npf_wetdry

  !> @brief Determine if a cell should rewet
  !!
  !! This method can be called from any external object that has a head that
  !! can be used to rewet the GWF cell node.  The ihc value is used to
  !! determine if it is a vertical or horizontal connection, which can operate
  !! differently depending on user settings.
  !<
  subroutine rewet_check(this, kiter, node, hm, ibdm, ihc, hnew, irewet)
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: node
    real(DP), intent(in) :: hm
    integer(I4B), intent(in) :: ibdm
    integer(I4B), intent(in) :: ihc
    real(DP), intent(inout), dimension(:) :: hnew
    integer(I4B), intent(out) :: irewet
    ! -- local
    integer(I4B) :: itflg
    real(DP) :: wd, awd, turnon, bbot
    !
    irewet = 0
    !
    ! -- Convert a dry cell to wet if it meets the criteria
    if (this%irewet > 0) then
      itflg = mod(kiter, this%iwetit)
      if (itflg == 0) then
        if (this%ibound(node) == 0 .and. this%wetdry(node) /= DZERO) then
          !
          ! -- Calculate wetting elevation
          bbot = this%dis%bot(node)
          wd = this%wetdry(node)
          awd = wd
          if (wd < 0) awd = -wd
          turnon = bbot + awd
          !
          ! -- Check head in adjacent cells to see if wetting elevation has
          !    been reached
          if (ihc == C3D_VERTICAL) then
            !
            ! -- check cell below
            if (ibdm > 0 .and. hm >= turnon) irewet = 1
          else
            if (wd > DZERO) then
              !
              ! -- check horizontally adjacent cells
              if (ibdm > 0 .and. hm >= turnon) irewet = 1
            end if
          end if
          !
          if (irewet == 1) then
            ! -- rewet cell; use equation 3a if ihdwet=0; use equation 3b if
            !    ihdwet is not 0.
            if (this%ihdwet == 0) then
              hnew(node) = bbot + this%wetfct * (hm - bbot)
            else
              hnew(node) = bbot + this%wetfct * awd !(hm - bbot)
            end if
            this%ibound(node) = 30000
          end if
        end if
      end if
    end if
  end subroutine rewet_check

  !> @brief Print wet/dry message
  !<
  subroutine sgwf_npf_wdmsg(this, icode, ncnvrt, nodcnvrt, acnvrt, ihdcnv, &
                            kiter, n)
    ! -- modules
    use TdisModule, only: kstp, kper
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: icode
    integer(I4B), intent(inout) :: ncnvrt
    character(len=30), dimension(5), intent(inout) :: nodcnvrt
    character(len=3), dimension(5), intent(inout) :: acnvrt
    integer(I4B), intent(inout) :: ihdcnv
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: n
    ! -- local
    integer(I4B) :: l
    ! -- formats
    character(len=*), parameter :: fmtcnvtn = &
      "(1X,/1X,'CELL CONVERSIONS FOR ITER.=',I0, &
       &'  STEP=',I0,'  PERIOD=',I0,'   (NODE or LRC)')"
    character(len=*), parameter :: fmtnode = "(1X,3X,5(A4, A20))"
    !
    ! -- Keep track of cell conversions
    if (icode > 0) then
      ncnvrt = ncnvrt + 1
      call this%dis%noder_to_string(n, nodcnvrt(ncnvrt))
      if (icode == 1) then
        acnvrt(ncnvrt) = 'DRY'
      else
        acnvrt(ncnvrt) = 'WET'
      end if
    end if
    !
    ! -- Print a line if 5 conversions have occurred or if icode indicates that a
    !    partial line should be printed
    if (ncnvrt == 5 .or. (icode == 0 .and. ncnvrt > 0)) then
      if (ihdcnv == 0) write (this%iout, fmtcnvtn) kiter, kstp, kper
      ihdcnv = 1
      write (this%iout, fmtnode) &
        (acnvrt(l), trim(adjustl(nodcnvrt(l))), l=1, ncnvrt)
      ncnvrt = 0
    end if
  end subroutine sgwf_npf_wdmsg

  !> @brief Calculate the effective hydraulic conductivity for the n-m connection
  !!
  !! n is primary node node number
  !! m is connected node (not used if vg is provided)
  !! ihc is horizontal indicator (0 vertical, 1 horizontal, 2 vertically
  !!   staggered)
  !! ipos_opt is position of connection in ja array
  !! vg is the global unit vector that expresses the direction from which to
  !!   calculate an effective hydraulic conductivity.
  !<
  function hy_eff(this, n, m, ihc, ipos, vg) result(hy)
    ! -- return
    real(DP) :: hy
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: ihc
    integer(I4B), intent(in), optional :: ipos
    real(DP), dimension(3), intent(in), optional :: vg
    ! -- local
    integer(I4B) :: iipos
    real(DP) :: hy11, hy22, hy33
    real(DP) :: ang1, ang2, ang3
    real(DP) :: vg1, vg2, vg3
    !
    ! -- Initialize
    iipos = 0
    if (present(ipos)) iipos = ipos
    hy11 = this%k11(n)
    hy22 = this%k11(n)
    hy33 = this%k11(n)
    hy22 = this%k22(n)
    hy33 = this%k33(n)
    !
    ! -- Calculate effective K based on whether connection is vertical
    !    or horizontal
    if (ihc == C3D_VERTICAL) then
      !
      ! -- Handle rotated anisotropy case that would affect the effective
      !    vertical hydraulic conductivity
      hy = hy33
      if (this%iangle2 > 0) then
        if (present(vg)) then
          vg1 = vg(1)
          vg2 = vg(2)
          vg3 = vg(3)
        else
          call this%dis%connection_normal(n, m, ihc, vg1, vg2, vg3, iipos)
        end if
        ang1 = this%angle1(n)
        ang2 = this%angle2(n)
        ang3 = DZERO
        if (this%iangle3 > 0) ang3 = this%angle3(n)
        hy = hyeff(hy11, hy22, hy33, ang1, ang2, ang3, vg1, vg2, vg3, &
                   this%iavgkeff)
      end if
      !
    else
      !
      ! -- Handle horizontal case
      hy = hy11
      if (this%ik22 > 0) then
        if (present(vg)) then
          vg1 = vg(1)
          vg2 = vg(2)
          vg3 = vg(3)
        else
          call this%dis%connection_normal(n, m, ihc, vg1, vg2, vg3, iipos)
        end if
        ang1 = DZERO
        ang2 = DZERO
        ang3 = DZERO
        if (this%iangle1 > 0) then
          ang1 = this%angle1(n)
          if (this%iangle2 > 0) then
            ang2 = this%angle2(n)
            if (this%iangle3 > 0) ang3 = this%angle3(n)
          end if
        end if
        hy = hyeff(hy11, hy22, hy33, ang1, ang2, ang3, vg1, vg2, vg3, &
                   this%iavgkeff)
      end if
      !
    end if
  end function hy_eff

  !> @brief Calculate the 3 components of specific discharge at the cell center
  !<
  subroutine calc_spdis(this, flowja)
    ! -- modules
    use SimModule, only: store_error
    ! -- dummy
    class(GwfNpfType) :: this
    real(DP), intent(in), dimension(:) :: flowja
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: m
    integer(I4B) :: ipos
    integer(I4B) :: iedge
    integer(I4B) :: isympos
    integer(I4B) :: ihc
    integer(I4B) :: ic
    integer(I4B) :: iz
    integer(I4B) :: nc
    integer(I4B) :: ncz
    real(DP) :: qz
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
    real(DP) :: dz
    real(DP) :: axy
    real(DP) :: ayx
    logical :: nozee = .true.
    type(SpdisWorkArrayType), pointer :: swa => null() !< pointer to spdis work arrays structure
    !
    ! -- Ensure dis has necessary information
    if (this%icalcspdis /= 0 .and. this%dis%con%ianglex == 0) then
      call store_error('Error.  ANGLDEGX not provided in '// &
                       'discretization file.  ANGLDEGX required for '// &
                       'calculation of specific discharge.', terminate=.TRUE.)
    end if

    swa => this%spdis_wa
    if (.not. swa%is_created()) then
      ! prepare work arrays
      call this%spdis_wa%create(this%calc_max_conns())

      ! prepare lookup table
      if (this%nedges > 0) call this%prepare_edge_lookup()
    end if
    !
    ! -- Go through each cell and calculate specific discharge
    do n = 1, this%dis%nodes
      !
      ! -- first calculate geometric properties for x and y directions and
      !    the specific discharge at a face (vi)
      ic = 0
      iz = 0

      ! reset work arrays
      call swa%reset()

      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        isympos = this%dis%con%jas(ipos)
        ihc = this%dis%con%ihc(isympos)
        area = this%dis%con%hwva(isympos)
        if (ihc == C3D_VERTICAL) then
          !
          ! -- vertical connection
          iz = iz + 1
          !call this%dis%connection_normal(n, m, ihc, xn, yn, zn, ipos)
          call this%dis%connection_vector(n, m, nozee, this%sat(n), this%sat(m), &
                                          ihc, xc, yc, zc, dltot)
          cl1 = this%dis%con%cl1(isympos)
          cl2 = this%dis%con%cl2(isympos)
          if (m < n) then
            cl1 = this%dis%con%cl2(isympos)
            cl2 = this%dis%con%cl1(isympos)
          end if
          ooclsum = DONE / (cl1 + cl2)
          swa%diz(iz) = dltot * cl1 * ooclsum
          qz = flowja(ipos)
          if (n > m) qz = -qz
          swa%viz(iz) = qz / area
        else
          !
          ! -- horizontal connection
          ic = ic + 1
          dz = thksatnm(this%ibound(n), this%ibound(m), &
                        this%icelltype(n), this%icelltype(m), &
                        this%inewton, ihc, &
                        this%hnew(n), this%hnew(m), this%sat(n), this%sat(m), &
                        this%dis%top(n), this%dis%top(m), this%dis%bot(n), &
                        this%dis%bot(m))
          area = area * dz
          call this%dis%connection_normal(n, m, ihc, xn, yn, zn, ipos)
          call this%dis%connection_vector(n, m, nozee, this%sat(n), this%sat(m), &
                                          ihc, xc, yc, zc, dltot)
          cl1 = this%dis%con%cl1(isympos)
          cl2 = this%dis%con%cl2(isympos)
          if (m < n) then
            cl1 = this%dis%con%cl2(isympos)
            cl2 = this%dis%con%cl1(isympos)
          end if
          ooclsum = DONE / (cl1 + cl2)
          swa%nix(ic) = -xn
          swa%niy(ic) = -yn
          swa%di(ic) = dltot * cl1 * ooclsum
          if (area > DZERO) then
            swa%vi(ic) = flowja(ipos) / area
          else
            swa%vi(ic) = DZERO
          end if
        end if
      end do

      ! add contribution from edge flows (i.e. from exchanges)
      if (this%nedges > 0) then
        do ipos = this%iedge_ptr(n), this%iedge_ptr(n + 1) - 1
          iedge = this%edge_idxs(ipos)

          ! propsedge: (Q, area, nx, ny, distance)
          ihc = this%ihcedge(iedge)
          area = this%propsedge(2, iedge)
          if (ihc == C3D_VERTICAL) then
            iz = iz + 1
            swa%viz(iz) = this%propsedge(1, iedge) / area
            swa%diz(iz) = this%propsedge(5, iedge)
          else
            ic = ic + 1
            swa%nix(ic) = -this%propsedge(3, iedge)
            swa%niy(ic) = -this%propsedge(4, iedge)
            swa%di(ic) = this%propsedge(5, iedge)
            if (area > DZERO) then
              swa%vi(ic) = this%propsedge(1, iedge) / area
            else
              swa%vi(ic) = DZERO
            end if
          end if
        end do
      end if
      !
      ! -- Assign number of vertical and horizontal connections
      ncz = iz
      nc = ic
      !
      ! -- calculate z weight (wiz) and z velocity
      if (ncz == 1) then
        swa%wiz(1) = DONE
      else
        dsumz = DZERO
        do iz = 1, ncz
          dsumz = dsumz + swa%diz(iz)
        end do
        denom = (ncz - DONE)
        if (denom < DZERO) denom = DZERO
        dsumz = dsumz + DEM10 * dsumz
        do iz = 1, ncz
          if (dsumz > DZERO) swa%wiz(iz) = DONE - swa%diz(iz) / dsumz
          if (denom > 0) then
            swa%wiz(iz) = swa%wiz(iz) / denom
          else
            swa%wiz(iz) = DZERO
          end if
        end do
      end if
      vz = DZERO
      do iz = 1, ncz
        vz = vz + swa%wiz(iz) * swa%viz(iz)
      end do
      !
      ! -- distance-based weighting
      nc = ic
      dsumx = DZERO
      dsumy = DZERO
      dsumz = DZERO
      do ic = 1, nc
        swa%wix(ic) = swa%di(ic) * abs(swa%nix(ic))
        swa%wiy(ic) = swa%di(ic) * abs(swa%niy(ic))
        dsumx = dsumx + swa%wix(ic)
        dsumy = dsumy + swa%wiy(ic)
      end do
      !
      ! -- Finish computing omega weights.  Add a tiny bit
      !    to dsum so that the normalized omega weight later
      !    evaluates to (essentially) 1 in the case of a single
      !    relevant connection, avoiding 0/0.
      dsumx = dsumx + DEM10 * dsumx
      dsumy = dsumy + DEM10 * dsumy
      do ic = 1, nc
        swa%wix(ic) = (dsumx - swa%wix(ic)) * abs(swa%nix(ic))
        swa%wiy(ic) = (dsumy - swa%wiy(ic)) * abs(swa%niy(ic))
      end do
      !
      ! -- compute B weights
      dsumx = DZERO
      dsumy = DZERO
      do ic = 1, nc
        swa%bix(ic) = swa%wix(ic) * sign(DONE, swa%nix(ic))
        swa%biy(ic) = swa%wiy(ic) * sign(DONE, swa%niy(ic))
        dsumx = dsumx + swa%wix(ic) * abs(swa%nix(ic))
        dsumy = dsumy + swa%wiy(ic) * abs(swa%niy(ic))
      end do
      if (dsumx > DZERO) dsumx = DONE / dsumx
      if (dsumy > DZERO) dsumy = DONE / dsumy
      axy = DZERO
      ayx = DZERO
      do ic = 1, nc
        swa%bix(ic) = swa%bix(ic) * dsumx
        swa%biy(ic) = swa%biy(ic) * dsumy
        axy = axy + swa%bix(ic) * swa%niy(ic)
        ayx = ayx + swa%biy(ic) * swa%nix(ic)
      end do
      !
      ! -- Calculate specific discharge.  The divide by zero checking below
      !    is problematic for cells with only one flow, such as can happen
      !    with triangular cells in corners.  In this case, the resulting
      !    cell velocity will be calculated as zero.  The method should be
      !    improved so that edge flows of zero are included in these
      !    calculations.  But this needs to be done with consideration for LGR
      !    cases in which flows are submitted from an exchange.
      vx = DZERO
      vy = DZERO
      do ic = 1, nc
        vx = vx + (swa%bix(ic) - axy * swa%biy(ic)) * swa%vi(ic)
        vy = vy + (swa%biy(ic) - ayx * swa%bix(ic)) * swa%vi(ic)
      end do
      denom = DONE - axy * ayx
      if (denom /= DZERO) then
        vx = vx / denom
        vy = vy / denom
      end if
      !
      this%spdis(1, n) = vx
      this%spdis(2, n) = vy
      this%spdis(3, n) = vz
      !
    end do

  end subroutine calc_spdis

  !> @brief Save specific discharge in binary format to ibinun
  !<
  subroutine sav_spdis(this, ibinun)
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: ibinun
    ! -- local
    character(len=16) :: text
    character(len=16), dimension(3) :: auxtxt
    integer(I4B) :: n
    integer(I4B) :: naux
    !
    ! -- Write the header
    text = '      DATA-SPDIS'
    naux = 3
    auxtxt(:) = ['              qx', '              qy', '              qz']
    call this%dis%record_srcdst_list_header(text, this%name_model, &
                                            this%packName, this%name_model, &
                                            this%packName, naux, auxtxt, ibinun, &
                                            this%dis%nodes, this%iout)
    !
    ! -- Write a zero for Q, and then write qx, qy, qz as aux variables
    do n = 1, this%dis%nodes
      call this%dis%record_mf6_list_entry(ibinun, n, n, DZERO, naux, &
                                          this%spdis(:, n))
    end do
  end subroutine sav_spdis

  !> @brief Save saturation in binary format to ibinun
  !<
  subroutine sav_sat(this, ibinun)
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: ibinun
    ! -- local
    character(len=16) :: text
    character(len=16), dimension(1) :: auxtxt
    real(DP), dimension(1) :: a
    integer(I4B) :: n
    integer(I4B) :: naux
    !
    ! -- Write the header
    text = '        DATA-SAT'
    naux = 1
    auxtxt(:) = ['             sat']
    call this%dis%record_srcdst_list_header(text, this%name_model, &
                                            this%packName, this%name_model, &
                                            this%packName, naux, auxtxt, ibinun, &
                                            this%dis%nodes, this%iout)
    !
    ! -- Write a zero for Q, and then write saturation as an aux variables
    do n = 1, this%dis%nodes
      a(1) = this%sat(n)
      call this%dis%record_mf6_list_entry(ibinun, n, n, DZERO, naux, a)
    end do
  end subroutine sav_sat

  !> @brief Reserve space for nedges cells that have an edge on them.
  !!
  !! This must be called before the npf%allocate_arrays routine, which is
  !! called from npf%ar.
  !<
  subroutine increase_edge_count(this, nedges)
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: nedges
    !
    this%nedges = this%nedges + nedges
  end subroutine increase_edge_count

  !> @brief Calculate the maximum number of connections for any cell
  !<
  function calc_max_conns(this) result(max_conns)
    class(GwfNpfType) :: this
    integer(I4B) :: max_conns
    ! local
    integer(I4B) :: n, m, ic

    max_conns = 0
    do n = 1, this%dis%nodes

      ! Count internal model connections
      ic = this%dis%con%ia(n + 1) - this%dis%con%ia(n) - 1

      ! Add edge connections
      do m = 1, this%nedges
        if (this%nodedge(m) == n) then
          ic = ic + 1
        end if
      end do

      ! Set max number of connections for any cell
      if (ic > max_conns) max_conns = ic
    end do

  end function calc_max_conns

  !> @brief Provide the npf package with edge properties
  !<
  subroutine set_edge_properties(this, nodedge, ihcedge, q, area, nx, ny, &
                                 distance)
    ! -- dummy
    class(GwfNpfType) :: this
    integer(I4B), intent(in) :: nodedge
    integer(I4B), intent(in) :: ihcedge
    real(DP), intent(in) :: q
    real(DP), intent(in) :: area
    real(DP), intent(in) :: nx
    real(DP), intent(in) :: ny
    real(DP), intent(in) :: distance
    ! -- local
    integer(I4B) :: lastedge
    !
    this%lastedge = this%lastedge + 1
    lastedge = this%lastedge
    this%nodedge(lastedge) = nodedge
    this%ihcedge(lastedge) = ihcedge
    this%propsedge(1, lastedge) = q
    this%propsedge(2, lastedge) = area
    this%propsedge(3, lastedge) = nx
    this%propsedge(4, lastedge) = ny
    this%propsedge(5, lastedge) = distance
    !
    ! -- If this is the last edge, then the next call must be starting a new
    !    edge properties assignment loop, so need to reset lastedge to 0
    if (this%lastedge == this%nedges) this%lastedge = 0
  end subroutine set_edge_properties

  subroutine prepare_edge_lookup(this)
    class(GwfNpfType) :: this
    ! local
    integer(I4B) :: i, inode, iedge
    integer(I4B) :: n, start, end
    integer(I4B) :: prev_cnt, strt_idx, ipos

    do i = 1, size(this%iedge_ptr)
      this%iedge_ptr(i) = 0
    end do
    do i = 1, size(this%edge_idxs)
      this%edge_idxs(i) = 0
    end do

    ! count
    do iedge = 1, this%nedges
      n = this%nodedge(iedge)
      this%iedge_ptr(n) = this%iedge_ptr(n) + 1
    end do

    ! determine start indexes
    prev_cnt = this%iedge_ptr(1)
    this%iedge_ptr(1) = 1
    do inode = 2, this%dis%nodes + 1
      strt_idx = this%iedge_ptr(inode - 1) + prev_cnt
      prev_cnt = this%iedge_ptr(inode)
      this%iedge_ptr(inode) = strt_idx
    end do

    ! loop over edges to fill lookup table
    do iedge = 1, this%nedges
      n = this%nodedge(iedge)
      start = this%iedge_ptr(n)
      end = this%iedge_ptr(n + 1) - 1
      do ipos = start, end
        if (this%edge_idxs(ipos) > 0) cycle ! go to next
        this%edge_idxs(ipos) = iedge
        exit
      end do
    end do

  end subroutine prepare_edge_lookup

  !> Calculate saturated thickness between cell n and m
  !<
  function calcSatThickness(this, n, m, ihc) result(satThickness)
    ! -- dummy
    class(GwfNpfType) :: this !< this NPF instance
    integer(I4B) :: n !< node n
    integer(I4B) :: m !< node m
    integer(I4B) :: ihc !< 1 = horizontal connection, 0 for vertical
    ! -- return
    real(DP) :: satThickness !< saturated thickness
    !
    satThickness = thksatnm(this%ibound(n), &
                            this%ibound(m), &
                            this%icelltype(n), &
                            this%icelltype(m), &
                            this%inewton, &
                            ihc, &
                            this%hnew(n), &
                            this%hnew(m), &
                            this%sat(n), &
                            this%sat(m), &
                            this%dis%top(n), &
                            this%dis%top(m), &
                            this%dis%bot(n), &
                            this%dis%bot(m))
  end function calcSatThickness

end module GwfNpfModule
