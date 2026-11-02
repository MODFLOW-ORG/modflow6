module GwfSwiModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, DONE, DZERO, DHALF, LENBUDTXT, &
                             MNORMAL, DHNOFLO, C3D_VERTICAL, LENMODELNAME, &
                             DEM6, LENMEMPATH
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use GwfNpfModule, only: GwfNpfType
  use GwfStoModule, only: GwfStoType
  use SmoothingModule, only: sQuadraticSaturation, &
                             sQuadraticSaturationDerivative, sSCurve, &
                             sQSaturation, sQSaturationDerivative
  use MemoryManagerModule, only: mem_setptr, mem_allocate
  use MemoryHelperModule, only: create_mem_path
  use MatrixBaseModule
  use GwfConductanceUtilsModule, only: hcond, vcond
  use UtlTvgModule, only: UtlTvgType
  use GwfNpfFormulationModule, only: GwfNpfFormulationType, SWI_FLOW
  use GwfStoExtModule, only: GwfStoFormulationType, SWI_STORAGE
  use FlowReductionInterfaceModule, only: FlowReductionType

  implicit none
  private
  public :: GwfSwiType
  public :: swi_cr

  character(len=LENBUDTXT), dimension(1) :: budtxt = & !< text labels for budget terms
    &['         STORAGE']
  integer(I4B), parameter :: FRESHWATER_ONE_FLUID = 0 !< freshwater model without swi-swi exchange
  integer(I4B), parameter :: FRESHWATER_TWO_FLUID = 1 !< freshwater model with swi-swi exchange
  integer(I4B), parameter :: SALTWATER_ONE_FLUID = 2 !< saltwater model without swi-swi exchange, not supported
  integer(I4B), parameter :: SALTWATER_TWO_FLUID = 3 !< saltwater model with swi-swi exchange
  real(DP), parameter :: DEFAULT_RELDENSDIFF = 0.025D0 !< default relative density difference, (rhos - rhof) / rhof, for 1000. and 1025.
  real(DP), parameter :: ZONE_MINIMUM_THICKNESS = 1.d-3 !< minimum thickness of saltwater or freshwater zone to avoid numerical issues

  type, extends(NumericalPackageType) :: GwfSwiType

    character(len=LENMODELNAME) :: freshwater_model_name = ''
    character(len=LENMODELNAME) :: saltwater_model_name = ''

    integer(I4B), pointer :: iconfiguration => null() !< 0 is freshwater, 1 is freshwater two fluid, 2 is saltwater, 3 is saltwater two fluid
    integer(I4B), pointer :: isaltwater => null() !< 0 is freshwater, 1 is saltwater
    integer(I4B), pointer :: izetaout => null() !< unit number for binary zeta output file
    integer(I4B), pointer :: intvg => null() !< flag: a TVG6 file supplies hsalt
    integer(I4B), pointer :: inoupstream => null() !< flag: average (do not upstream-weight) the fluid-zone conductance; dev option
    character(len=LENMEMPATH) :: tvg6_mempath = '' !< input memory path of the TVG subpackage
    real(DP), dimension(:), pointer, contiguous :: zeta => null() !< interface elevation from the converged heads
    real(DP), dimension(:), pointer, contiguous :: storage => null() !< calculated swi storage

    ! information needed for full implementation
    real(DP), pointer :: reldensdiff => null() !< relative density difference, (rhos - rhof) / rhof; user input, default 0.025
    real(DP), pointer :: alphaf => null() !< rhof / (rhos - rhof), derived from reldensdiff, default is 40
    real(DP), pointer :: alphas => null() !< rhos / (rhos - rhof), default is 41
    real(DP), pointer :: hsalt_user => null() !< user specified head of saltwater (zero by default)
    real(DP), dimension(:), pointer, contiguous :: hfresh => null() !< head of freshwater
    real(DP), dimension(:), pointer, contiguous :: hsalt => null() !< head of saltwater
    real(DP), dimension(:), pointer, contiguous :: hfreshold => null() !< old head of freshwater
    real(DP), dimension(:), pointer, contiguous :: hsaltold => null() !< old head of saltwater
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound

    ! pointers for transient simulations
    integer(I4B), pointer :: insto => null() !< pointer to check of storage package is on
    integer(I4B), pointer :: iss => null() !< pointer to gwf steady state flag
    real(DP), dimension(:), pointer, contiguous :: sy => null() !< pointer to storage package specific yield
    real(DP), dimension(:), pointer, contiguous :: oldsy => null() !< pointer to storage package previous specific yield
    real(DP), dimension(:), pointer, contiguous :: oldss => null() !< pointer to storage package previous specific storage

    ! objects
    type(UtlTvgType) :: tvg_reader

    ! SWI flow formulation (overrides the default NPF conductance)
    class(GwfNpfFormulationType), pointer :: npf_form => null()
    ! SWI storage formulation (overrides the default STO terms; adds bd/save_flows)
    class(GwfStoFormulationType), pointer :: sto_form => null()
    ! SWI flow-reduction provider (answers stress packages that taper a
    ! withdrawal, such as WEL AUTO_FLOW_REDUCE, with the fluid-zone factor)
    class(FlowReductionType), pointer :: flowred => null()

  contains

    procedure :: swi_df
    procedure :: swi_ar
    procedure, private :: check_compatibility
    procedure :: swi_rp
    procedure :: swi_ad
    procedure :: swi_cq
    procedure :: swi_ot_dv
    procedure :: swi_da
    procedure, private :: swi_load
    procedure, private :: source_options
    procedure, private :: log_options
    procedure, private :: set_alphas
    procedure :: allocate_scalars
    procedure, private :: allocate_arrays
    procedure, private :: update_zeta

    procedure :: set_configuration
    procedure :: set_head_pointers
    procedure :: get_zetanew
    procedure :: get_zetaold

  end type GwfSwiType

  !> @brief NPF flow formulation for SWI. Overrides the default NPF connection
  !! conductance with the fluid-zone conductance (freshwater S^f = S^w - S^s, or
  !! saltwater S^s).
  !<
  type, extends(GwfNpfFormulationType) :: SwiNpfFormulationType
    class(GwfSwiType), pointer :: swi => null() !< owning SWI package
    type(GwfNpfType), pointer :: npf => null() !< the model NPF package
  contains
    procedure :: fc => swinpf_fc
    procedure :: fn => swinpf_fn
    procedure :: cq => swinpf_cq
  end type SwiNpfFormulationType

  !> @brief Flow-reduction provider for SWI. Answers a stress package that
  !! tapers a withdrawal (WEL AUTO_FLOW_REDUCE) with the factor for the fluid
  !! zone this model simulates, so the withdrawal tapers as the zone thins
  !! rather than as the cell empties, in confined cells too. The interface
  !! moves with the head, and the derivative carries that dependence, so the
  !! well term is consistent with the Newton formulation and needs no lag.
  !<
  type, extends(FlowReductionType) :: SwiFlowReductionType
    class(GwfSwiType), pointer :: swi => null() !< owning SWI package
    type(GwfNpfType), pointer :: npf => null() !< the model NPF package (cell type)
  contains
    procedure :: reduce => swi_reduce
  end type SwiFlowReductionType

  !> @brief STO storage formulation for SWI. In a freshwater model it assembles
  !! the freshwater storage as the water column (S^w, via the STO default) minus
  !! the salt column (S^s); in a saltwater model it supplies the salt-column
  !! storage. fc fills the residual-consistent terms, fn adds the smoothed
  !! interface tangent, and cq/bd/save_flows report the storage and the SWI
  !! storage budget term. Dispatched from sto_fc/fn/cq for every active cell;
  !! swi_ar flags the cells SWI_STORAGE only in a saltwater model, so that the
  !! STO default skips them.
  !<
  type, extends(GwfStoFormulationType) :: SwiStoFormulationType
    class(GwfSwiType), pointer :: swi => null() !< owning SWI package
    type(GwfStoType), pointer :: sto => null() !< the model STO package
  contains
    procedure :: is_active => swisto_is_active
    procedure :: fc => swisto_fc
    procedure :: fn => swisto_fn
    procedure :: cq => swisto_cq
    procedure :: bd => swisto_bd
    procedure :: save_flows => swisto_save_flows
  end type SwiStoFormulationType

contains

  !> @brief Create a new swi package object
  !<
  subroutine swi_cr(swi, name_model, input_mempath, inunit, iout, dis)
    ! dummy
    type(GwfSwiType), pointer :: swi
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    class(DisBaseType), pointer, intent(in) :: dis
    ! formats
    character(len=*), parameter :: fmtswi = &
      "(1x, /1x, 'SWI -- Seawater Intrusion Package, Version 8, 1/8/2024', &
      &' input read from mempath: ', A, //)"
    !
    ! create SWI object
    allocate (swi)
    !
    ! create name and memory path
    call swi%set_names(1, name_model, 'SWI', 'SWI', input_mempath)
    !
    ! allocate scalars
    call swi%allocate_scalars()
    !
    ! set variables
    swi%inunit = inunit
    swi%iout = iout
    !
    ! set pointers
    swi%dis => dis
    !
    ! check if pkg is enabled,
    if (inunit > 0) then
      ! print message identifying pkg
      write (swi%iout, fmtswi) input_mempath
    end if
    !
    ! return
    return
  end subroutine swi_cr

  !> @brief Allocate arrays, load from IDM, and assign head
  !<
  subroutine swi_df(this)
    ! dummy
    class(GwfSwiType) :: this
    ! local

    ! allocate arrays
    call this%allocate_arrays(this%dis%nodes)

    ! load from IDM
    call this%swi_load()

    ! setup the tvg reader
    if (this%intvg > 0) then
      call this%tvg_reader%initialize(this%dis, this%tvg6_mempath, this%iout, &
                                      this%name_model, 'SWI')
      if (.not. this%tvg_reader%has('SALTWATER_HEAD')) then
        write (errmsg, '(a)') 'TVG Input file does not have SALTWATER_HEAD &
          &variable.  Specify SALTWATER_HEAD as a time varying grid array in &
          &the TVG input file.'
        call store_error(errmsg, terminate=.false.)
        call store_error_filename(this%input_fname)
      end if
    end if

  end subroutine swi_df

  !> @brief Setup pointers
  !<
  subroutine swi_ar(this, ibound, npf, sto)
    ! dummy
    class(GwfSwiType), target :: this
    integer(I4B), dimension(:), pointer, contiguous :: ibound !< model ibound
    type(GwfNpfType), pointer, intent(inout) :: npf !< model NPF package
    type(GwfStoType), pointer, intent(inout) :: sto !< model STO package
    ! local
    integer(I4B) :: n
    integer(I4B) :: ipos
    type(SwiNpfFormulationType), pointer :: npf_form
    type(SwiStoFormulationType), pointer :: sto_form
    type(SwiFlowReductionType), pointer :: flowred

    ! set pointer to ibound
    this%ibound => ibound

    ! stop if the model uses a package the SWI formulation does not support
    call this%check_compatibility(npf)

    ! set pointer to gwf steady state flag
    call mem_setptr(this%insto, 'INSTO', &
                    create_mem_path(this%name_model))
    call mem_setptr(this%iss, 'ISS', &
                    create_mem_path(this%name_model))
    if (this%insto > 0) then
      call mem_setptr(this%sy, 'SY', &
                      create_mem_path(this%name_model, 'STO'))
      ! The SWI Package uses specific yield (SY) as the drainable porosity for
      ! the interface-movement storage. A moving interface changes the freshwater
      ! and saltwater storage even where the aquifer is confined, so SY is
      ! required whenever the STO Package is active with SWI. STO zero-fills SY
      ! when it is not provided, so an all-zero SY means it was not specified.
      if (maxval(this%sy) <= DZERO) then
        write (errmsg, '(a)') &
          'The SWI Package requires that specific yield (SY) be specified in '// &
          'the STO Package. A moving interface produces a drainable '// &
          '(specific-yield) storage change in both the freshwater and '// &
          'saltwater zones, even where the aquifer is confined.'
        call store_error(errmsg, terminate=.false.)
        call store_error_filename(this%input_fname)
      end if
    end if

    ! If hfresh is not associated, then this swi package must not
    ! be part of a model that has a swi-swi exchange with another
    ! model.  In this case, isaltwater must be 0 because we assume
    ! that if a model is not connected with a swi-swi exchange then
    ! it must be a freshwater model.
    if (.not. associated(this%hfresh)) then
      call mem_setptr(this%hfresh, 'X', &
                      create_mem_path(this%name_model))
    end if
    if (.not. associated(this%hfreshold)) then
      call mem_setptr(this%hfreshold, 'XOLD', &
                      create_mem_path(this%name_model))
    end if

    ! the saltwater head is simulated in a two-fluid model, so it cannot also
    ! be supplied through TVG
    if (this%intvg > 0 .and. this%iconfiguration /= FRESHWATER_ONE_FLUID) then
      call store_error('A TVG file supplies SALTWATER_HEAD, but the model is &
                       &part of a two-fluid (SWI-SWI) simulation, in which the &
                       &saltwater head is simulated. Remove the TVG file.', &
                       terminate=.false.)
      call store_error_filename(this%input_fname)
    end if

    ! If this is a single-fluid freshwater model
    ! then memory needs to be allocated for hsalt
    if (this%iconfiguration == FRESHWATER_ONE_FLUID) then

      ! allocate saltwater head array
      call mem_allocate(this%hsalt, size(this%hfresh), 'HSALT', this%memoryPath)
      do n = 1, size(this%hsalt)
        this%hsalt(n) = this%hsalt_user
      end do

      ! allocate old saltwater head array
      call mem_allocate(this%hsaltold, size(this%hfresh), 'HSALTOLD', &
                        this%memoryPath)
      do n = 1, size(this%hsaltold)
        this%hsaltold(n) = this%hsalt_user
      end do
    end if

    ! register the SWI NPF flow formulation and flag ALL of its connections
    ! (horizontal and vertical) so NPF dispatches them to the SWI override; the
    ! formulation applies the fluid-zone conductance for horizontal connections
    ! and the buoyancy-restricted conductance for vertical connections
    allocate (npf_form)
    npf_form%swi => this
    npf_form%npf => npf
    this%npf_form => npf_form
    call npf%add_flow_formulation(this%npf_form, SWI_FLOW)
    do n = 1, npf%dis%nodes
      do ipos = npf%dis%con%ia(n) + 1, npf%dis%con%ia(n + 1) - 1
        npf%iformulation(ipos) = SWI_FLOW
      end do
    end do

    ! the flow-reduction provider; the model hands it to the stress packages
    allocate (flowred)
    flowred%swi => this
    flowred%npf => npf
    this%flowred => flowred

    ! register the SWI STO storage formulation. Storage formulations compose
    ! additively, so the claim mask is only set where SWI must replace the
    ! standard storage entirely. A freshwater model assembles its storage as
    ! the water column (S^w) minus the salt column (S^s): the default storage
    ! formulation supplies the water column and SWI leaves its cells unclaimed
    ! so that it runs, adding only the salt-column correction. A saltwater
    ! model has no water-column term of its own, so it claims its cells and the
    ! default skips them. The SWI storage budget term is the
    ! interface-movement (fresh<->salt conversion) volume.
    if (this%insto > 0) then
      allocate (sto_form)
      sto_form%swi => this
      sto_form%sto => sto
      this%sto_form => sto_form
      call sto%add_sto_formulation(this%sto_form, SWI_STORAGE)
      if (this%isaltwater == 1) then
        do n = 1, sto%dis%nodes
          sto%iformulation(n) = SWI_STORAGE
        end do
      end if
    end if

  end subroutine swi_ar

  !> @brief Store an error for each model package the SWI formulation
  !! does not support.
  !!
  !! XT3D replaces the NPF conductance terms, so the SWI fluid-zone
  !! conductance would be dropped. GNC and HFB adjust the standard NPF
  !! conductance, and BUY and VSC change the density and viscosity terms
  !! that SWI already represents with its own formulation.
  !<
  subroutine check_compatibility(this, npf)
    ! dummy
    class(GwfSwiType) :: this
    type(GwfNpfType), pointer, intent(in) :: npf
    ! local
    integer(I4B), pointer :: inpkg => null()
    character(len=4), dimension(4) :: pkgs = ['GNC ', 'HFB ', 'BUY ', 'VSC ']
    integer(I4B) :: i
    logical(LGP) :: found
    !
    found = .false.
    if (npf%ixt3d /= 0) then
      found = .true.
      call store_error('The SWI Package cannot be used with XT3D in the '// &
                       'NPF Package.', terminate=.false.)
    end if
    do i = 1, size(pkgs)
      call mem_setptr(inpkg, 'IN'//trim(pkgs(i)), &
                      create_mem_path(this%name_model))
      if (inpkg > 0) then
        found = .true.
        call store_error('The SWI Package cannot be used with the '// &
                         trim(pkgs(i))//' Package.', terminate=.false.)
      end if
    end do
    if (found) call store_error_filename(this%input_fname)
  end subroutine check_compatibility

  !> @brief Read and prepare
  !<
  subroutine swi_rp(this)
    ! dummy
    class(GwfSwiType) :: this
    ! local

  end subroutine swi_rp

  !> @brief Advance
  !<
  subroutine swi_ad(this, irestore)
    use TdisModule, only: kper, kstp
    ! dummy
    class(GwfSwiType) :: this
    integer(I4B), intent(in) :: irestore
    ! local
    integer(I4B) :: n

    if (this%iconfiguration == FRESHWATER_ONE_FLUID) then
      if (irestore == 0) then
        ! advance hsaltold
        do n = 1, size(this%hsalt)
          this%hsaltold(n) = this%hsalt(n)
        end do
      else
        ! redoing time step so reset hsalt to hsaltold
        do n = 1, size(this%hsalt)
          this%hsalt(n) = this%hsaltold(n)
        end do
      end if
    end if

    ! the Input Data Model has loaded this period's arrays; copy the
    ! saltwater head into the package's own array
    if (this%intvg > 0) then
      if (this%iconfiguration == FRESHWATER_ONE_FLUID) then
        call this%tvg_reader%fill(this%hsalt, 'SALTWATER_HEAD')
        ! the first saltwater head is the initial condition, so the
        ! starting interface is computed from it rather than from the
        ! scalar SALTWATER_HEAD
        if (kper == 1 .and. kstp == 1) then
          do n = 1, size(this%hsalt)
            this%hsaltold(n) = this%hsalt(n)
          end do
        end if
      end if
    end if

  end subroutine swi_ad

  !> @ brief Calculate flows and storages for SWI package and adjust flowja
  !!
  !<
  subroutine swi_cq(this)
    ! dummy variables
    class(GwfSwiType) :: this !< GwfSwiType object
    ! local
    integer(I4B) :: n
    !
    ! Refresh the stored zeta array from the converged heads. zeta is not used by
    ! the solve (it is computed on the fly via get_zetanew wherever needed); the
    ! stored array is a derived quantity kept current here -- once per time step,
    ! post-convergence, before the budget/output phases -- so it can be written to
    ! the zeta output file (swi_ot_dv) and inspected by API/BMI users.
    call this%update_zeta()
    !
    ! The STO Package skips the storage formulations in a steady-state period,
    ! so clear the SWI storage term rather than report the last transient one.
    if (this%insto > 0) then
      if (this%iss == 1) then
        do n = 1, size(this%storage)
          this%storage(n) = DZERO
        end do
      end if
    end if
    !
    ! Intercell freshwater/saltwater flows (FLOW-JA-FACE) are assembled by the SWI
    ! NPF formulation (swinpf_cq, dispatched from npf_cq), consistent with the
    ! matrix fill in swinpf_fc; storage flows by the SWI STO formulation
    ! (swisto_cq). Both horizontal and vertical connections are handled: swinpf_cq
    ! shares swinpf_condf with swinpf_fc, which routes vertical connections to the
    ! buoyancy-restricted conductance (test_gwf_swi_vcol2 exercises an all-vertical
    ! saltwater column and confirms the budget balances).
  end subroutine swi_cq

  !> @brief Save zeta array to binary file
  !<
  subroutine swi_ot_dv(this, idvfl)
    ! dummy
    class(GwfSwiType) :: this
    integer(I4B), intent(in) :: idvfl
    ! local
    character(len=1) :: cdatafmp = ' ', editdesc = ' '
    integer(I4B) :: ibinun
    integer(I4B) :: iprint
    integer(I4B) :: nvaluesp
    integer(I4B) :: nwidthp
    real(DP) :: dinact
    !
    ! Set unit number for zeta output
    if (this%izetaout /= 0) then
      ibinun = 1
    else
      ibinun = 0
    end if
    if (idvfl == 0) ibinun = 0
    !
    ! save zeta array
    if (ibinun /= 0) then
      iprint = 0
      dinact = DHNOFLO
      !
      ! write zeta to binary file
      if (this%izetaout /= 0) then
        ibinun = this%izetaout
        call this%dis%record_array(this%zeta, this%iout, iprint, ibinun, &
                                   '            ZETA', cdatafmp, nvaluesp, &
                                   nwidthp, editdesc, dinact)
      end if
    end if
    !
    ! Return
    return
  end subroutine swi_ot_dv

  !> @brief Deallocate
  !<
  subroutine swi_da(this)
    ! modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    ! dummy
    class(GwfSwiType) :: this

    ! deallocate IDM memory
    call memorystore_remove(this%name_model, 'SWI', idm_context)

    if (this%intvg > 0) then
      call this%tvg_reader%tvg_da()
    end if

    ! deallocate arrays
    call mem_deallocate(this%zeta)
    call mem_deallocate(this%storage)

    if (this%iconfiguration == FRESHWATER_ONE_FLUID) then
      call mem_deallocate(this%hsalt)
      call mem_deallocate(this%hsaltold)
    end if

    ! deallocate scalars
    call mem_deallocate(this%iconfiguration)
    call mem_deallocate(this%isaltwater)
    call mem_deallocate(this%izetaout)
    call mem_deallocate(this%intvg)
    call mem_deallocate(this%inoupstream)
    call mem_deallocate(this%reldensdiff)
    call mem_deallocate(this%alphaf)
    call mem_deallocate(this%alphas)
    call mem_deallocate(this%hsalt_user)

    ! deallocate parent
    call this%NumericalPackageType%da()

    ! the flow and storage formulation objects; NPF and STO keep pointers to
    ! them only, so they are freed here
    if (associated(this%npf_form)) then
      deallocate (this%npf_form)
      nullify (this%npf_form)
    end if
    if (associated(this%sto_form)) then
      deallocate (this%sto_form)
      nullify (this%sto_form)
    end if

    ! the flow-reduction provider
    if (associated(this%flowred)) then
      deallocate (this%flowred)
      nullify (this%flowred)
    end if

    ! nullify pointers
    nullify (this%hfresh)
    nullify (this%hsalt)
    nullify (this%hfreshold)
    nullify (this%hsaltold)
    nullify (this%oldss)
    nullify (this%oldsy)

  end subroutine swi_da

  !> @brief Load data from IDM into package
  !<
  subroutine swi_load(this)
    ! modules
    ! dummy
    class(GwfSwiType) :: this

    call this%source_options()

  end subroutine swi_load

  !> @ brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the SWI package. The base model
  !! allocate scalars method is also called.
  !<
  subroutine allocate_scalars(this)
    ! modules
    ! dummy
    class(GwfSwiType) :: this

    ! allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()

    ! Allocate scalars
    call mem_allocate(this%iconfiguration, 'ICONFIGURATION', this%memoryPath)
    call mem_allocate(this%isaltwater, 'ISALTWATER', this%memoryPath)
    call mem_allocate(this%izetaout, 'IZETAOUT', this%memoryPath)
    call mem_allocate(this%intvg, 'INTVG', this%memoryPath)
    call mem_allocate(this%inoupstream, 'INOUPSTREAM', this%memoryPath)
    call mem_allocate(this%alphaf, 'ALPHAF', this%memoryPath)
    call mem_allocate(this%alphas, 'ALPHAS', this%memoryPath)
    call mem_allocate(this%reldensdiff, 'RELDENSDIFF', this%memoryPath)
    call mem_allocate(this%hsalt_user, 'HSALT_USER', this%memoryPath)

    ! Initialize value
    this%izetaout = 0
    this%intvg = 0
    this%inoupstream = 0
    this%reldensdiff = DEFAULT_RELDENSDIFF
    call this%set_alphas()
    this%hsalt_user = DZERO

    ! initialize iconfiguration and isaltwater
    call this%set_configuration("FRESHWATER_ONE_FLUID")

  end subroutine allocate_scalars

  !> @brief Allocate arrays
  !<
  subroutine allocate_arrays(this, nodes)
    ! modules
    use MemoryManagerModule, only: mem_allocate
    ! dummy
    class(GwfSwiType) :: this
    integer(I4B), intent(in) :: nodes
    ! local
    integer(I4B) :: n

    ! Allocate
    call mem_allocate(this%zeta, nodes, 'ZETA', this%memoryPath)
    call mem_allocate(this%storage, nodes, 'STORAGE', this%memoryPath)

    ! initialize
    do n = 1, nodes
      this%zeta(n) = DZERO
      this%storage(n) = DZERO
    end do

  end subroutine allocate_arrays

  !> @brief Update simulation options from input mempath
  !<
  subroutine source_options(this)
    ! modules
    use SourceCommonModule, only: filein_fname
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: getunit, openfile
    use MemoryManagerModule, only: mem_setptr
    use MemoryManagerExtModule, only: mem_set_value
    use CharacterStringModule, only: CharacterStringType
    use GwfSwiInputModule, only: GwfSwiParamFoundType
    ! dummy
    class(GwfSwiType) :: this
    character(len=LINELENGTH) :: zeta_fname
    character(len=LINELENGTH) :: tvg_fname
    type(CharacterStringType), dimension(:), pointer, contiguous :: tvg6_mempaths
    ! locals
    type(GwfSwiParamFoundType) :: found

    ! update defaults with idm sourced values
    zeta_fname = ''
    call mem_set_value(zeta_fname, 'ZETAFILE', this%input_mempath, &
                       found%zetafile)
    call mem_set_value(this%hsalt_user, 'HSALT_USER', this%input_mempath, &
                       found%hsalt_user)
    call mem_set_value(this%reldensdiff, 'RELDENSDIFF', this%input_mempath, &
                       found%reldensdiff)
    call mem_set_value(this%inoupstream, 'INOUPSTREAM', this%input_mempath, &
                       found%inoupstream)

    ! the interface coefficients follow from the relative density difference
    if (found%reldensdiff) then
      if (this%reldensdiff <= DZERO) then
        write (errmsg, '(a,G0,a)') &
          'RELATIVE_DENSITY_DIFFERENCE must be greater than zero. Value &
          &specified is ', this%reldensdiff, '.'
        call store_error(errmsg, terminate=.false.)
        call store_error_filename(this%input_fname)
      end if
      call this%set_alphas()
    end if

    ! open zeta file
    if (zeta_fname /= '') then
      this%izetaout = getunit()
      call openfile(this%izetaout, this%iout, zeta_fname, 'DATA(BINARY)', &
                    form, access, 'REPLACE', MNORMAL)
      write (this%iout, '(4x,a)') &
        'ZETA information will be written to ', trim(zeta_fname)
    end if

    ! open hsalt file
    ! TVG6 subpackage. The Input Data Model loads the file and records its
    ! memory path here; the reader is attached to that path in swi_df.
    tvg_fname = ''
    if (filein_fname(tvg_fname, 'TVG6_FILENAME', this%input_mempath, &
                     this%input_fname)) then
      call mem_setptr(tvg6_mempaths, 'TVG6_MEMPATH', this%input_mempath)
      this%tvg6_mempath = tvg6_mempaths(1)
      this%intvg = 1
      write (this%iout, '(4x,a,a)') &
        'Time varying grid array information will be read from ', trim(tvg_fname)
    end if

    ! log options
    if (this%iout > 0) then
      call this%log_options(found)
    end if

  end subroutine source_options

  !> @brief Log options sourced from the input mempath
  !<
  subroutine log_options(this, found)
    ! modules
    use KindModule, only: LGP
    use GwfSwiInputModule, only: GwfSwiParamFoundType
    ! dummy
    class(GwfSwiType) :: this
    ! locals
    type(GwfSwiParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'Setting SWI Options'
    if (found%zetafile) then
      write (this%iout, '(4x,a)') &
        'Zeta will be written to a binary output file.'
    end if
    if (found%hsalt_user) then
      write (this%iout, '(4x,a,G0,a)') &
        'Saltwater head was set to ', this%hsalt_user, '.'
    end if
    if (found%reldensdiff) then
      write (this%iout, '(4x,a,G0,a)') &
        'Relative density difference was set to ', this%reldensdiff, '.'
    end if
    if (found%tvg6_filename) then
      write (this%iout, '(4x,a)') &
        'Saltwater head will be read from file.'
    end if
    if (found%inoupstream) then
      write (this%iout, '(4x,a)') &
        'Upstream weighting of the fluid-zone conductance is DISABLED &
        &(dev option). The NPF cell-averaging method is used instead and &
        &the flow-term Newton derivatives are omitted.'
    end if
    write (this%iout, '(1x,a,/)') 'End Setting SWI Options'

  end subroutine log_options

  !> @brief Set the interface coefficients from the relative density difference
  !!
  !! alphaf = rhof / (rhos - rhof) = 1 / reldensdiff, and alphas = alphaf + 1.
  !! The two differ by exactly one, so the formulation has a single degree of
  !! freedom and the relative density difference is the only density input
  !! needed. The guard keeps a non-positive user value, which is reported as an
  !! input error, from dividing by zero before the run terminates.
  !<
  subroutine set_alphas(this)
    ! dummy
    class(GwfSwiType) :: this

    if (this%reldensdiff > DZERO) then
      this%alphaf = DONE / this%reldensdiff
      this%alphas = this%alphaf + DONE
    end if

  end subroutine set_alphas

  !> @brief Set this model configuration.
  !! This can be called from an exchange, such as the
  !! SwiSwi exchange, which knows which model/swi package
  !! is the saltwater model.
  !<
  subroutine set_configuration(this, config)
    ! dummy
    class(GwfSwiType) :: this
    character(len=*), intent(in) :: config

    select case (config)
    case ("FRESHWATER_ONE_FLUID")
      this%iconfiguration = FRESHWATER_ONE_FLUID
      this%isaltwater = 0
    case ("SALTWATER_ONE_FLUID")
      this%iconfiguration = SALTWATER_ONE_FLUID
      this%isaltwater = 1
    case ("FRESHWATER_TWO_FLUID")
      this%iconfiguration = FRESHWATER_TWO_FLUID
      this%isaltwater = 0
    case ("SALTWATER_TWO_FLUID")
      this%iconfiguration = SALTWATER_TWO_FLUID
      this%isaltwater = 1
    end select
  end subroutine set_configuration

  !> @brief Set head pointers
  !! This will likely be called from the SWI-SWI exchange
  !! which has access to both fresh and salt models.
  !<
  subroutine set_head_pointers(this, hfresh, hsalt, hfreshold, hsaltold)
    class(GwfSwiType) :: this
    real(DP), dimension(:), target, intent(in) :: hfresh
    real(DP), dimension(:), target, intent(in) :: hsalt
    real(DP), dimension(:), target, intent(in) :: hfreshold
    real(DP), dimension(:), target, intent(in) :: hsaltold
    this%hfresh => hfresh
    this%hsalt => hsalt
    this%hfreshold => hfreshold
    this%hsaltold => hsaltold
  end subroutine set_head_pointers

  !> @brief Get new value for zeta
  !!
  !! Include hsalt contribution if available and use perturbations if
  !! provided.
  !!
  !<
  function get_zetanew(this, n, eps_fresh, eps_salt) result(zetanew)
    ! dummy
    class(GwfSwiType) :: this !< this instance
    integer(I4B) :: n !< cell number
    real(DP), optional, intent(in) :: eps_fresh !< perturbation for fresh head
    real(DP), optional, intent(in) :: eps_salt !< perturbation for salt head
    ! return
    real(DP) :: zetanew
    ! locals
    real(DP) :: eps_f
    real(DP) :: eps_s
    ! The head perturbations are used only by the two-fluid exchange when it forms
    ! the numerical cross-derivatives; in single-fluid mode they are never passed
    ! (eps_s stays zero and hsalt holds the prescribed saltwater head), so the
    ! same expression serves both configurations. zeta is the unclamped
    ! Ghyben-Herzberg map -- a vanished fluid zone is handled by the smoothed
    ! saturation and the minimum-zone conductance floor, not by clamping zeta.
    eps_f = DZERO
    eps_s = DZERO
    if (present(eps_fresh)) eps_f = eps_fresh
    if (present(eps_salt)) eps_s = eps_salt
    zetanew = calc_zeta(this%alphaf, this%hfresh(n) + eps_f, &
                        this%alphas, this%hsalt(n) + eps_s)
  end function get_zetanew

  !> @brief Get zeta from the previous time-step heads
  !<
  function get_zetaold(this, n) result(zetaold)
    ! dummy
    class(GwfSwiType) :: this !< this instance
    integer(I4B) :: n !< cell number
    ! return
    real(DP) :: zetaold
    zetaold = calc_zeta(this%alphaf, this%hfreshold(n), &
                        this%alphas, this%hsaltold(n))
  end function get_zetaold

  !> @brief Calculate zeta surface
  !<
  subroutine update_zeta(this)
    ! modules
    ! dummy
    class(GwfSwiType) :: this
    ! locals
    integer(I4B) :: n

    ! Loop through each node and calculate zeta
    do n = 1, this%dis%nodes
      ! inactive cells are written as no-flow cells
      if (this%ibound(n) == 0) then
        this%zeta(n) = DHNOFLO
        cycle
      end if

      ! Calculate zeta
      this%zeta(n) = this%get_zetanew(n)

      ! todo: Do we want to constrain zeta to top and
      !    bot of cell?
      if (this%zeta(n) > this%dis%top(n)) then
        this%zeta(n) = this%dis%top(n)
      end if
      if (this%zeta(n) < this%dis%bot(n)) then
        this%zeta(n) = this%dis%bot(n)
      end if

    end do

  end subroutine update_zeta

  function calc_zeta(alphaf, hf, alphas, hs) result(zeta)
    real(DP), intent(in) :: alphaf
    real(DP), intent(in) :: hf
    real(DP), intent(in) :: alphas
    real(DP), intent(in) :: hs
    real(DP) :: zeta
    zeta = -alphaf * hf + alphas * hs
  end function calc_zeta

  !> @brief Fraction of the cell below the interface (saltwater fraction)
  !<
  subroutine swi_thksat(top, bot, zeta, thksat, inewton)
    ! dummy
    real(DP), intent(in) :: top
    real(DP), intent(in) :: bot
    real(DP), intent(in) :: zeta
    real(DP), intent(inout) :: thksat
    integer(I4B), intent(in) :: inewton
    !
    ! Standard Formulation
    if (zeta >= top) then
      thksat = DONE
    else if (zeta <= bot) then
      thksat = DZERO
    else
      thksat = (zeta - bot) / (top - bot)
    end if
    !
    ! Smoothed thickness
    if (inewton /= 0) then
      thksat = sQuadraticSaturation(top, bot, zeta)
    end if
    !
    ! Return
    return
  end subroutine swi_thksat

  !> @brief SWI flow-reduction provider: factor for a withdrawal from cell n,
  !! and its derivative with respect to this model's head.
  !!
  !! Freshwater model: the freshwater zone runs from the interface, or the cell
  !! bottom where the interface is below it, up to the water table or the cell
  !! top. The factor tapers up from the zone bottom over the taper length and
  !! is evaluated at the zone top, so a confined cell whose head is far above
  !! the cell top is still reduced as the interface fills the cell. The
  !! interface moves with the head (dzeta/dh^f = -alphaf), so where it is the
  !! zone bottom the derivative carries 1 + alphaf; where the water table is
  !! above the cell top the zone top is fixed and the head term drops out.
  !!
  !! Saltwater model: the saltwater zone runs from the cell bottom up to the
  !! interface or the cell top. The factor tapers up from the cell bottom and
  !! is evaluated at the zone top; dzeta/dh^s = alphas.
  !<
  subroutine swi_reduce(this, n, h, taper, f, dfdh)
    class(SwiFlowReductionType), intent(in) :: this
    integer(I4B), intent(in) :: n
    real(DP), intent(in) :: h
    real(DP), intent(in) :: taper
    real(DP), intent(out) :: f
    real(DP), intent(out) :: dfdh
    ! local
    real(DP) :: tp, bt, zeta, x, dxdh, bz, dbzdh
    f = DONE
    dfdh = DZERO
    if (taper <= DZERO) return
    tp = this%npf%dis%top(n)
    bt = this%npf%dis%bot(n)
    zeta = this%swi%get_zetanew(n)
    if (this%swi%isaltwater == 0) then
      ! freshwater zone: bottom at the interface (or the cell bottom), top at
      ! the water table (or the cell top)
      if (zeta > bt) then
        bz = zeta
        dbzdh = -this%swi%alphaf
      else
        bz = bt
        dbzdh = DZERO
      end if
      if (this%npf%icelltype(n) /= 0 .and. h < tp) then
        x = h
        dxdh = DONE
      else
        x = tp
        dxdh = DZERO
      end if
    else
      ! saltwater zone: bottom at the cell bottom, top at the interface (or
      ! the cell top)
      bz = bt
      dbzdh = DZERO
      if (zeta < tp) then
        x = zeta
        dxdh = this%swi%alphas
      else
        x = tp
        dxdh = DZERO
      end if
    end if
    f = sQSaturation(bz + taper, bz, x)
    dfdh = sQSaturationDerivative(bz + taper, bz, x) * (dxdh - dbzdh)
  end subroutine swi_reduce

  !> @brief SWI NPF formulation: vertical conductance for a connection between an
  !! upper and a lower cell.
  !!
  !! Single-fluid (freshwater-only): the interface is a diagnostic derived from a
  !! continuous freshwater potential, so the standard (unrestricted) vertical
  !! conductance is returned. SWI still owns the connection -- so the treatment
  !! can be changed later -- but does not adjust it. Restricting the vertical flow
  !! here would decouple the deep cells below the interface and leave the head
  !! field underdetermined.
  !!
  !! Two-fluid (buoyancy-restricted): each fluid moves freely in its buoyant
  !! direction -- freshwater up, saltwater down -- and is restricted in the other
  !! by the fluid present at the shared face: freshwater flows down only if the
  !! interface is below the face (fresh at the face); saltwater flows up only if
  !! the interface is above the face (salt at the face). Both switches -- the
  !! interface position and the flow direction -- are smoothed with sSCurve and
  !! combined with a smooth OR, so the conductance stays differentiable as the
  !! interface crosses the face. The smoothing band is a fraction of the thickness
  !! of the cell whose interface is being tested (the upper cell for freshwater,
  !! the lower cell for saltwater), so it is meaningful when the two cells differ
  !! in thickness. The restricted case is reduced to a residual (1e-6)
  !! conductance rather than zero, to avoid a singular row.
  !<
  function swinpf_vcondf(this, n, m, ipos, hnew) result(condf)
    class(SwiNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: ipos
    real(DP), dimension(:), intent(in) :: hnew
    real(DP) :: condf
    ! local
    type(GwfNpfType), pointer :: npf
    integer(I4B) :: iu, il, jas, ihc
    real(DP) :: zface, band, hband, dydx, factor, sface, sup, hyn, hym
    !
    npf => this%npf
    jas = npf%dis%con%jas(ipos)
    !
    ! -- single-fluid: standard (unrestricted) vertical conductance
    if (this%swi%iconfiguration == FRESHWATER_ONE_FLUID) then
      ihc = npf%dis%con%ihc(jas)
      hyn = npf%hy_eff(n, m, ihc, ipos=ipos)
      hym = npf%hy_eff(m, n, ihc, ipos=ipos)
      condf = vcond(npf%ibound(n), npf%ibound(m), &
                    npf%icelltype(n), npf%icelltype(m), npf%inewton, &
                    npf%ivarcv, npf%idewatcv, npf%condsat(jas), &
                    hnew(n), hnew(m), hyn, hym, npf%sat(n), npf%sat(m), &
                    npf%dis%top(n), npf%dis%top(m), &
                    npf%dis%bot(n), npf%dis%bot(m), npf%dis%con%hwva(jas))
      return
    end if
    ! -- no flow to or from an inactive cell, as in vcond
    condf = DZERO
    if (npf%ibound(n) == 0 .or. npf%ibound(m) == 0) return
    ! -- identify the upper (iu) and lower (il) cell; the shared face is at
    !    bot(iu) = top(il)
    if (npf%dis%bot(n) > npf%dis%bot(m)) then
      iu = n
      il = m
    else
      iu = m
      il = n
    end if
    zface = npf%dis%bot(iu)
    !
    ! -- each fluid passes if it is present at the face OR flows in its buoyant
    !    direction (smooth "OR"): freshwater passes if fresh-at-face OR upflow;
    !    saltwater passes if salt-at-face OR downflow (= 1 - upflow). Smoothing
    !    scales: the interface switch smooths over a fraction of the thickness of
    !    the cell that holds the tested interface; the flow-direction switch
    !    smooths over the equivalent head band (a zeta change of `band` is a head
    !    change of band/alpha).  TODO: expose the 0.1 fraction as input.
    if (this%swi%isaltwater == 0) then
      ! freshwater: the tested interface zeta(iu) is in the UPPER cell
      band = 1.0d-1 * (npf%dis%top(iu) - npf%dis%bot(iu))
      hband = band / this%swi%alphaf
      ! smoothed upflow indicator (1 for upflow, h_lower > h_upper; 0 for down)
      call sSCurve((hnew(il) - hnew(iu)) + DHALF * hband, hband, dydx, sup)
      ! sface = smoothed "fresh at the face" (interface in upper cell below face)
      call sSCurve((zface - this%swi%get_zetanew(iu)) + DHALF * band, band, &
                   dydx, sface)
      factor = DONE - (DONE - sface) * (DONE - sup)
    else
      ! saltwater: the tested interface zeta(il) is in the LOWER cell
      band = 1.0d-1 * (npf%dis%top(il) - npf%dis%bot(il))
      hband = band / this%swi%alphas
      ! smoothed upflow indicator (1 for upflow, h_lower > h_upper; 0 for down)
      call sSCurve((hnew(il) - hnew(iu)) + DHALF * hband, hband, dydx, sup)
      ! sface = smoothed "salt at the face" (interface in lower cell above face)
      call sSCurve((this%swi%get_zetanew(il) - zface) + DHALF * band, band, &
                   dydx, sface)
      factor = DONE - (DONE - sface) * sup
    end if
    ! -- reduce (not zero) to avoid a singular row where a fluid is blocked
    condf = npf%condsat(jas) * max(factor, DEM6)
  end function swinpf_vcondf

  !> @brief SWI NPF formulation: fluid-zone conductance for connection ipos.
  !! Freshwater zone (S^w - S^s, between zeta and the water table) or saltwater
  !! zone (S^s, between the bottom and zeta). The conductance follows the model's
  !! Newton setting (via hcond): the NPF cell averaging for Picard, upstream
  !! weighting for Newton. Floored at the conductance of a minimum-thickness zone
  !! to guard a vanished fluid zone. Vertical connections are passed to
  !! swinpf_vcondf. Shared by fc and cq so the matrix fill and the cell-by-cell
  !! flow use the same conductance.
  !<
  function swinpf_condf(this, n, m, ipos, hnew) result(condf)
    class(SwiNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: ipos
    real(DP), dimension(:), intent(in) :: hnew
    real(DP) :: condf
    ! local
    type(GwfNpfType), pointer :: npf
    integer(I4B) :: ihc, jas, ictn, ictm
    real(DP) :: hyn, hym, cond, condsw
    real(DP) :: satn_s, satm_s, tn, tm, bn, bm, tthkn, tthkm
    real(DP) :: satn_f, satm_f, swn, swm
    !
    npf => this%npf
    condf = DZERO
    jas = npf%dis%con%jas(ipos)
    ihc = npf%dis%con%ihc(jas)
    if (ihc == C3D_VERTICAL) then
      condf = swinpf_vcondf(this, n, m, ipos, hnew)
      return
    end if
    tn = npf%dis%top(n)
    tm = npf%dis%top(m)
    bn = npf%dis%bot(n)
    bm = npf%dis%bot(m)
    tthkn = tn - bn
    tthkm = tm - bm
    hyn = npf%hy_eff(n, m, ihc, ipos=ipos)
    hym = npf%hy_eff(m, n, ihc, ipos=ipos)
    !
    ! -- DEV_NO_UPSTREAM_WEIGHTING: average the fluid-zone fraction of the two
    !    cells with the NPF cell-averaging method (ICELLAVG) instead of taking
    !    the upstream cell's fraction. The fraction is averaged directly rather
    !    than as the water-column minus salt-column difference used below,
    !    because that difference equals the fluid-zone conductance only for
    !    upstream (or arithmetic-mean) weighting; a harmonic mean of the two
    !    columns is not the harmonic mean of their difference. Cell type is
    !    forced convertible so the fraction is honored for confined cells too.
    !    Smoothing follows the Newton setting, as elsewhere.
    if (this%swi%inoupstream /= 0) then
      call swi_thksat(tn, bn, this%swi%get_zetanew(n), satn_s, npf%inewton)
      call swi_thksat(tm, bm, this%swi%get_zetanew(m), satm_s, npf%inewton)
      if (this%swi%isaltwater == 1) then
        satn_f = satn_s
        satm_f = satm_s
      else
        swn = DONE
        swm = DONE
        if (npf%icelltype(n) /= 0) swn = npf%sat(n)
        if (npf%icelltype(m) /= 0) swm = npf%sat(m)
        satn_f = max(swn - satn_s, DZERO)
        satm_f = max(swm - satm_s, DZERO)
      end if
      condf = hcond(npf%ibound(n), npf%ibound(m), 1, 1, 0, ihc, &
                    npf%icellavg, npf%condsat(jas), hnew(n), hnew(m), &
                    satn_f, satm_f, hyn, hym, tn, tm, bn, bm, &
                    npf%dis%con%cl1(jas), npf%dis%con%cl2(jas), &
                    npf%dis%con%hwva(jas))
      condf = max(condf, npf%condsat(jas) * ZONE_MINIMUM_THICKNESS / &
                  min(tthkn, tthkm))
      return
    end if
    !
    ! -- full water-column (S^w) conductance, as in fc_default_flow
    cond = hcond(npf%ibound(n), npf%ibound(m), &
                 npf%icelltype(n), npf%icelltype(m), npf%inewton, ihc, &
                 npf%icellavg, npf%condsat(jas), hnew(n), hnew(m), &
                 npf%sat(n), npf%sat(m), hyn, hym, tn, tm, bn, bm, &
                 npf%dis%con%cl1(jas), npf%dis%con%cl2(jas), &
                 npf%dis%con%hwva(jas))
    !
    ! -- saltwater-column (S^s) conductance below the interface (zeta); celltype
    !    forced to 1 so zeta acts as the top surface. zeta is recomputed on the
    !    fly from the current heads (get_zetanew) so the conductance is consistent
    !    with the storage terms and not lagged.
    ictn = 1
    ictm = 1
    call swi_thksat(tn, bn, this%swi%get_zetanew(n), satn_s, npf%inewton)
    call swi_thksat(tm, bm, this%swi%get_zetanew(m), satm_s, npf%inewton)
    condsw = hcond(npf%ibound(n), npf%ibound(m), ictn, ictm, npf%inewton, ihc, &
                   npf%icellavg, npf%condsat(jas), hnew(n), hnew(m), &
                   satn_s, satm_s, hyn, hym, tn, tm, bn, bm, &
                   npf%dis%con%cl1(jas), npf%dis%con%cl2(jas), &
                   npf%dis%con%hwva(jas))
    !
    if (this%swi%isaltwater == 1) then
      condf = condsw
    else
      condf = cond - condsw
    end if
    ! -- guard: floor the model's own fluid-zone conductance at the conductance of
    !    a minimum-thickness zone (condsat scaled by ZONE_MINIMUM_THICKNESS/tthk).
    !    Symmetric for both models, it keeps a fully drained fluid zone from
    !    producing a zero/negative row -- a singular system, or an unphysical
    !    negative conductance where S^s > S^w -- without perturbing normally
    !    saturated cells. TODO: replacing this fixed threshold with proper
    !    dry-zone handling would let ZONE_MINIMUM_THICKNESS be removed entirely.
    condf = max(condf, npf%condsat(jas) * ZONE_MINIMUM_THICKNESS / &
                min(tthkn, tthkm))
  end function swinpf_condf

  subroutine swinpf_fc_conn(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)
    class(SwiNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: ipos
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: hnew
    ! local
    type(GwfNpfType), pointer :: npf
    integer(I4B) :: idiag, isymcon, idiagm
    real(DP) :: condf
    !
    npf => this%npf
    condf = swinpf_condf(this, n, m, ipos, hnew)
    idiag = npf%dis%con%ia(n)
    isymcon = npf%dis%con%isym(ipos)
    idiagm = npf%dis%con%ia(m)
    call matrix_sln%add_value_pos(idxglo(ipos), condf)
    call matrix_sln%add_value_pos(idxglo(idiag), -condf)
    call matrix_sln%add_value_pos(idxglo(isymcon), condf)
    call matrix_sln%add_value_pos(idxglo(idiagm), -condf)
  end subroutine swinpf_fc_conn

  !> @brief SWI NPF formulation: Newton terms for the freshwater flow override.
  !! Mirrors fn_default_flow but uses the freshwater saturation derivative
  !! dS^f/dh = S'(h) + alphaf*S'(zeta) of the upstream cell. Horizontal only.
  !<
  subroutine swinpf_fn_conn(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)
    use SmoothingModule, only: sQuadraticSaturationDerivative
    class(SwiNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: ipos
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: hnew
    ! local
    type(GwfNpfType), pointer :: npf
    integer(I4B) :: isymcon, idiag, idiagm, iups, idn, jas
    real(DP) :: cond, consterm, derv, term, topup, botup, zetaup
    !
    npf => this%npf
    jas = npf%dis%con%jas(ipos)
    ! Vertical connections use no Newton term (Picard only). A Newton term for
    ! the two-fluid buoyancy-restricted vertical conductance was implemented and
    ! tested (the restriction factor depends on both connected heads, via the
    ! smoothed flow-direction switch and, through zeta, the interface-at-face
    ! switch). It did not materially improve convergence -- outer-iteration counts
    ! were within a few percent of Picard (and slightly higher) across interface
    ! speeds of 0.1-0.4 m/step, with and without under-relaxation -- because the
    ! vertical restriction is not the bottleneck (interface storage is). It
    ! was therefore not included. Single-fluid vertical flow is unrestricted.
    if (npf%dis%con%ihc(jas) == C3D_VERTICAL) return
    ! DEV_NO_UPSTREAM_WEIGHTING: the derivative below is the upstream cell's
    ! fraction derivative. An averaged conductance depends on both cells'
    ! fractions; that derivative is not implemented, so the flow term is
    ! Picard and only the (weighting-independent) storage Newton terms remain.
    if (this%swi%inoupstream /= 0) return
    idiag = npf%dis%con%ia(n)
    isymcon = npf%dis%con%isym(ipos)
    ! determine upstream node (by freshwater head)
    iups = m
    if (hnew(m) < hnew(n)) iups = n
    idn = n
    if (iups == n) idn = m
    topup = npf%dis%top(iups)
    botup = npf%dis%bot(iups)
    if (npf%dis%con%ihc(jas) == 2) then
      topup = min(npf%dis%top(n), npf%dis%top(m))
      botup = max(npf%dis%bot(n), npf%dis%bot(m))
    end if
    ! saturated conductance and the fluid saturation derivative of the upstream
    ! cell: freshwater dS^f/dh^f = S'(h) + alphaf*S'(zeta); saltwater
    ! dS^s/dh^s = alphas*S'(zeta)
    cond = npf%condsat(jas)
    consterm = -cond * (hnew(iups) - hnew(idn))
    zetaup = this%swi%get_zetanew(iups)
    if (this%swi%isaltwater == 1) then
      derv = this%swi%alphas * &
             sQuadraticSaturationDerivative(topup, botup, zetaup, npf%satomega)
    else
      derv = sQuadraticSaturationDerivative(topup, botup, hnew(iups), &
                                            npf%satomega) + &
             this%swi%alphaf * &
             sQuadraticSaturationDerivative(topup, botup, zetaup, npf%satomega)
    end if
    idiagm = npf%dis%con%ia(m)
    if (iups == n) then
      term = consterm * derv
      rhs(n) = rhs(n) + term * hnew(n)
      rhs(m) = rhs(m) - term * hnew(n)
      call matrix_sln%add_value_pos(idxglo(idiag), term)
      if (npf%ibound(m) > 0) then
        call matrix_sln%add_value_pos(idxglo(isymcon), -term)
      end if
    else
      term = -consterm * derv
      rhs(n) = rhs(n) + term * hnew(m)
      rhs(m) = rhs(m) - term * hnew(m)
      if (npf%ibound(n) > 0) then
        call matrix_sln%add_value_pos(idxglo(ipos), term)
      end if
      call matrix_sln%add_value_pos(idxglo(idiagm), -term)
    end if
  end subroutine swinpf_fn_conn

  !> @brief SWI NPF formulation: intercell (FLOW-JA-FACE) flow for connection
  !! ipos, using the same fluid-zone conductance as swinpf_fc so the cell-by-cell
  !! flows are consistent with the matrix fill (and hence the CHD/GHB budgets,
  !! which are recovered from the flowja residual). Handles both horizontal and
  !! vertical connections: swinpf_condf routes vertical connections to the
  !! buoyancy-restricted conductance (swinpf_vcondf).
  !<
  subroutine swinpf_cq_conn(this, n, m, ipos, flowja, h_new)
    class(SwiNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: ipos
    real(DP), dimension(:), intent(inout) :: flowja
    real(DP), dimension(:), intent(in) :: h_new
    ! local
    type(GwfNpfType), pointer :: npf
    real(DP) :: condf, qnm
    !
    npf => this%npf
    condf = swinpf_condf(this, n, m, ipos, h_new)
    qnm = condf * (h_new(m) - h_new(n))
    flowja(ipos) = qnm
    flowja(npf%dis%con%isym(ipos)) = -qnm
  end subroutine swinpf_cq_conn

  !> @brief SWI NPF formulation: fill coefficients over all SWI connections.
  !!
  !! The formulation owns its traversal: it visits every connection and
  !! contributes only on the faces SWI claims through npf%iformulation, which
  !! the default conductance formulation skips.
  !<
  subroutine swinpf_fc(this, kiter, matrix_sln, idxglo, rhs, hnew)
    class(SwiNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(inout) :: rhs
    real(DP), dimension(:), intent(inout) :: hnew
    ! local
    type(GwfNpfType), pointer :: npf
    integer(I4B) :: n, m, ipos
    !
    npf => this%npf
    do n = 1, npf%dis%nodes
      do ipos = npf%dis%con%ia(n) + 1, npf%dis%con%ia(n + 1) - 1
        if (npf%dis%con%mask(ipos) == 0) cycle
        m = npf%dis%con%ja(ipos)
        if (m < n) cycle
        if (npf%iformulation(ipos) /= SWI_FLOW) cycle
        call swinpf_fc_conn(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)
      end do
    end do
  end subroutine swinpf_fc

  !> @brief SWI NPF formulation: Newton terms over all SWI connections.
  !<
  subroutine swinpf_fn(this, kiter, matrix_sln, idxglo, rhs, hnew)
    class(SwiNpfFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(inout) :: rhs
    real(DP), dimension(:), intent(inout) :: hnew
    ! local
    type(GwfNpfType), pointer :: npf
    integer(I4B) :: n, m, ipos
    !
    npf => this%npf
    do n = 1, npf%dis%nodes
      do ipos = npf%dis%con%ia(n) + 1, npf%dis%con%ia(n + 1) - 1
        if (npf%dis%con%mask(ipos) == 0) cycle
        m = npf%dis%con%ja(ipos)
        if (m < n) cycle
        if (npf%iformulation(ipos) /= SWI_FLOW) cycle
        call swinpf_fn_conn(this, n, m, ipos, matrix_sln, rhs, idxglo, hnew)
      end do
    end do
  end subroutine swinpf_fn

  !> @brief SWI NPF formulation: flows over all SWI connections.
  !!
  !! Masked connections are visited here, matching the default flow
  !! formulation, which does not apply the connection mask in cq.
  !<
  subroutine swinpf_cq(this, hnew, flowja)
    class(SwiNpfFormulationType), intent(inout) :: this
    real(DP), dimension(:), intent(inout) :: hnew
    real(DP), dimension(:), intent(inout) :: flowja
    ! local
    type(GwfNpfType), pointer :: npf
    integer(I4B) :: n, m, ipos
    !
    npf => this%npf
    do n = 1, npf%dis%nodes
      do ipos = npf%dis%con%ia(n) + 1, npf%dis%con%ia(n + 1) - 1
        m = npf%dis%con%ja(ipos)
        if (m < n) cycle
        if (npf%iformulation(ipos) /= SWI_FLOW) cycle
        call swinpf_cq_conn(this, n, m, ipos, flowja, hnew)
      end do
    end do
  end subroutine swinpf_cq

  !> @brief SWI STO formulation: is this cell handled by SWI storage?
  !<
  function swisto_is_active(this, n) result(is_active)
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    logical(LGP) :: is_active
    is_active = .true.
  end function swisto_is_active

  !> @brief SWI STO formulation: freshwater storage fill for cell n, in override
  !! form. The freshwater storage is the difference of two bottom-referenced
  !! columns, S^f = S^w - S^s: the water column (S^w) is filled by the STO
  !! default, and the salt column (S^s) is subtracted here. Each contribution is
  !! residual-consistent (its aterm*h - rhsterm equals the smoothed storage rate),
  !! so it converges under Picard on its own; the smoothed-tangent upgrade for the
  !! interface-movement (drainable) term is added in swisto_fn (Newton only).
  !<
  subroutine swisto_fc_node(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
    use TdisModule, only: delt
    use GwfStorageUtilsModule, only: SsCapacity, SyCapacity, SsTerms
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: h_old
    real(DP), dimension(:), intent(in) :: h_new
    ! local
    type(GwfStoType), pointer :: sto
    integer(I4B) :: idiag
    real(DP) :: tled, tp, bt, tthk
    real(DP) :: sc1, rho1, rho1old, sc2, rho2, rho2old
    real(DP) :: snew, sold, ssnew, ssold, aterm, rhsterm, rate_s
    !
    sto => this%sto
    if (this%swi%ibound(n) < 1) return
    ! -- saltwater model: storage over the salt column S^s (driven by h^s)
    if (this%swi%isaltwater == 1) then
      call swisto_fc_salt(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
      return
    end if
    tled = DONE / delt
    idiag = sto%dis%con%ia(n)
    tp = sto%dis%top(n)
    bt = sto%dis%bot(n)
    tthk = tp - bt
    !
    ! -- water column (S^w) is assembled by the default storage formulation,
    !    which runs first because a freshwater model does not claim its cells;
    !    SWI adds only the salt-column correction below
    !
    ! -- salt column (S^s) smoothed saturations: S^s = sat(zeta)
    ssnew = sQuadraticSaturation(tp, bt, this%swi%get_zetanew(n), sto%satomega)
    ssold = sQuadraticSaturation(tp, bt, this%swi%get_zetaold(n), sto%satomega)
    !
    ! -- compressible (Ss): subtract the saltwater fraction of the S^w term
    if (sto%iconvert(n) == 0) then
      sold = DONE
      snew = DONE
    else
      sold = sQuadraticSaturation(tp, bt, h_old(n), sto%satomega)
      snew = sQuadraticSaturation(tp, bt, h_new(n), sto%satomega)
    end if
    sc1 = SsCapacity(sto%istor_coef, tp, bt, sto%dis%area(n), sto%ss(n))
    rho1 = sc1 * tled
    if (sto%integratechanges /= 0) then
      rho1old = SsCapacity(sto%istor_coef, tp, bt, sto%dis%area(n), &
                           this%swi%oldss(n)) * tled
    else
      rho1old = rho1
    end if
    call SsTerms(sto%iconvert(n), sto%iorig_ss, sto%iconf_ss, tp, bt, &
                 rho1, rho1old, snew, sold, h_new(n), h_old(n), aterm, rhsterm)
    aterm = aterm * ssnew
    rhsterm = rhsterm * ssnew
    call matrix_sln%add_value_pos(idxglo(idiag), -aterm)
    rhs(n) = rhs(n) - rhsterm
    !
    ! -- drainable (Sy) salt column, residual-consistent chord. The residual
    !    contribution is -rate_s (freshwater = water - salt), applied always;
    !    the chord diagonal d(-rate_s)/dh ~ -alphaf*rho2 (interior) is applied
    !    only while the interface is within the cell. swisto_fn upgrades this
    !    chord to the smoothed interface tangent under Newton.
    sc2 = SyCapacity(sto%dis%area(n), this%swi%sy(n))
    rho2 = sc2 * tled
    if (sto%integratechanges /= 0) then
      rho2old = SyCapacity(sto%dis%area(n), this%swi%oldsy(n)) * tled
    else
      rho2old = rho2
    end if
    rate_s = rho2old * tthk * ssold - rho2 * tthk * ssnew
    if (ssnew > DZERO .and. ssnew < DONE) then
      aterm = -rho2 * this%swi%alphaf
    else
      aterm = DZERO
    end if
    rhsterm = aterm * h_new(n) + rate_s
    call matrix_sln%add_value_pos(idxglo(idiag), aterm)
    rhs(n) = rhs(n) + rhsterm
  end subroutine swisto_fc_node

  !> @brief SWI STO formulation: storage Newton terms for cell n. Upgrades the
  !! residual-consistent chord diagonals filled in swisto_fc to the exact smoothed
  !! tangents: the water column via the STO default, and the salt column via the
  !! interface-movement delta. The delta cancels at convergence (the residual is
  !! set in swisto_fc), so it only sharpens the Jacobian. Called only when Newton
  !! is active (sto_fn dispatch).
  !<
  subroutine swisto_fn_node(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
    use TdisModule, only: delt
    use GwfStorageUtilsModule, only: SyCapacity
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: h_old
    real(DP), dimension(:), intent(in) :: h_new
    ! local
    type(GwfStoType), pointer :: sto
    integer(I4B) :: idiag
    real(DP) :: tled, tp, bt, tthk
    real(DP) :: sc2, rho2, ssnew, sderiv, a_fn
    !
    sto => this%sto
    if (this%swi%ibound(n) < 1) return
    ! -- saltwater model: salt-column interface tangent (driven by h^s)
    if (this%swi%isaltwater == 1) then
      call swisto_fn_salt(this, n, matrix_sln, rhs, idxglo, h_new)
      return
    end if
    !
    ! -- water column (S^w) Newton tangent comes from the default storage
    !    formulation; SWI adds only the salt-column tangent below
    !
    ! -- salt column (S^s): smoothed interface tangent delta over the chord in
    !    swisto_fc, applied only while the interface is within the cell
    tled = DONE / delt
    idiag = sto%dis%con%ia(n)
    tp = sto%dis%top(n)
    bt = sto%dis%bot(n)
    tthk = tp - bt
    ssnew = sQuadraticSaturation(tp, bt, this%swi%get_zetanew(n), sto%satomega)
    if (ssnew > DZERO .and. ssnew < DONE) then
      sc2 = SyCapacity(sto%dis%area(n), this%swi%sy(n))
      rho2 = sc2 * tled
      sderiv = sQuadraticSaturationDerivative(tp, bt, this%swi%get_zetanew(n), &
                                              sto%satomega)
      ! net salt tangent d(-rate_s)/dh = -alphaf*rho2*tthk*sderiv; the fc chord
      ! was -alphaf*rho2, so the fn delta is alphaf*rho2*(1 - tthk*sderiv)
      a_fn = this%swi%alphaf * rho2 * (DONE - tthk * sderiv)
      call matrix_sln%add_value_pos(idxglo(idiag), a_fn)
      rhs(n) = rhs(n) + a_fn * h_new(n)
    end if
  end subroutine swisto_fn_node

  !> @brief SWI STO formulation: freshwater storage rate for cell n, in override
  !! form. The water column (S^w) rate comes from the STO default (fills
  !! strgss/strgsy and flowja); the salt-column (S^s) correction is subtracted
  !! here. The interface-movement (drainable) rate is reported separately as the
  !! SWI storage term (swi%storage) -- the fresh<->salt conversion volume -- and
  !! also added to flowja.
  !<
  subroutine swisto_cq_node(this, n, flowja, h_new, h_old)
    use TdisModule, only: delt
    use GwfStorageUtilsModule, only: SsCapacity, SyCapacity, SsTerms, SyTerms
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(inout) :: flowja
    real(DP), dimension(:), intent(in) :: h_new
    real(DP), dimension(:), intent(in) :: h_old
    ! local
    type(GwfStoType), pointer :: sto
    integer(I4B) :: idiag
    real(DP) :: tled, tp, bt, tthk
    real(DP) :: sc1, rho1, rho1old, sc2, rho2, rho2old
    real(DP) :: snew, sold, ssnew, ssold, aterm, rhsterm, rate, rate_s
    !
    sto => this%sto
    this%swi%storage(n) = DZERO
    ! -- saltwater model: storage over the salt column S^s (driven by h^s)
    if (this%swi%isaltwater == 1) then
      call swisto_cq_salt(this, n, flowja, h_new, h_old)
      return
    end if
    ! -- water column (S^w) rate (strgss/strgsy) is filled by the default
    !    storage formulation; SWI accumulates its salt-column correction
    if (sto%iss == 1) return
    if (this%swi%ibound(n) < 1) return
    tled = DONE / delt
    idiag = sto%dis%con%ia(n)
    tp = sto%dis%top(n)
    bt = sto%dis%bot(n)
    tthk = tp - bt
    ssnew = sQuadraticSaturation(tp, bt, this%swi%get_zetanew(n), sto%satomega)
    ssold = sQuadraticSaturation(tp, bt, this%swi%get_zetaold(n), sto%satomega)
    !
    ! -- compressible (Ss): subtract the saltwater fraction from strgss and flowja
    if (sto%iconvert(n) == 0) then
      sold = DONE
      snew = DONE
    else
      sold = sQuadraticSaturation(tp, bt, h_old(n), sto%satomega)
      snew = sQuadraticSaturation(tp, bt, h_new(n), sto%satomega)
    end if
    sc1 = SsCapacity(sto%istor_coef, tp, bt, sto%dis%area(n), sto%ss(n))
    rho1 = sc1 * tled
    if (sto%integratechanges /= 0) then
      rho1old = SsCapacity(sto%istor_coef, tp, bt, sto%dis%area(n), &
                           this%swi%oldss(n)) * tled
    else
      rho1old = rho1
    end if
    call SsTerms(sto%iconvert(n), sto%iorig_ss, sto%iconf_ss, tp, bt, &
                 rho1, rho1old, snew, sold, h_new(n), h_old(n), &
                 aterm, rhsterm, rate)
    rate = ssnew * rate
    sto%strgss(n) = sto%strgss(n) - rate
    flowja(idiag) = flowja(idiag) - rate
    !
    ! -- drainable (Sy): interface-movement storage. The SWI storage budget term
    !    is the salt-column rate (the fresh<->salt conversion volume), which is
    !    also added to flowja to complete the freshwater storage.
    sc2 = SyCapacity(sto%dis%area(n), this%swi%sy(n))
    rho2 = sc2 * tled
    if (sto%integratechanges /= 0) then
      rho2old = SyCapacity(sto%dis%area(n), this%swi%oldsy(n)) * tled
    else
      rho2old = rho2
    end if
    rate_s = rho2old * tthk * ssold - rho2 * tthk * ssnew
    this%swi%storage(n) = -rate_s
    flowja(idiag) = flowja(idiag) + this%swi%storage(n)
  end subroutine swisto_cq_node

  !> @brief SWI STO formulation (saltwater model): storage fill for cell n. The
  !! saltwater storage is over the salt column S^s = sat(zeta), which is
  !! bottom-referenced and driven by h^s (dzeta/dh^s = alphas). Residual-consistent
  !! chord in fc; the smoothed interface tangent is added in swisto_fn_salt.
  !<
  subroutine swisto_fc_salt(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
    use TdisModule, only: delt
    use GwfStorageUtilsModule, only: SsCapacity, SyCapacity, SsTerms
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: h_old
    real(DP), dimension(:), intent(in) :: h_new
    ! local
    type(GwfStoType), pointer :: sto
    integer(I4B) :: idiag
    real(DP) :: tled, tp, bt, tthk
    real(DP) :: sc1, rho1, rho1old, sc2, rho2, rho2old
    real(DP) :: ssnew, ssold, aterm, rhsterm, rate_s
    !
    sto => this%sto
    tled = DONE / delt
    idiag = sto%dis%con%ia(n)
    tp = sto%dis%top(n)
    bt = sto%dis%bot(n)
    tthk = tp - bt
    ssnew = sQuadraticSaturation(tp, bt, this%swi%get_zetanew(n), sto%satomega)
    ssold = sQuadraticSaturation(tp, bt, this%swi%get_zetaold(n), sto%satomega)
    !
    ! -- compressible (Ss) over the salt column, driven by h^s
    sc1 = SsCapacity(sto%istor_coef, tp, bt, sto%dis%area(n), sto%ss(n))
    rho1 = sc1 * tled
    if (sto%integratechanges /= 0) then
      rho1old = SsCapacity(sto%istor_coef, tp, bt, sto%dis%area(n), &
                           this%swi%oldss(n)) * tled
    else
      rho1old = rho1
    end if
    call SsTerms(sto%iconvert(n), sto%iorig_ss, sto%iconf_ss, tp, bt, &
                 rho1, rho1old, ssnew, ssold, h_new(n), h_old(n), aterm, rhsterm)
    call matrix_sln%add_value_pos(idxglo(idiag), aterm)
    rhs(n) = rhs(n) + rhsterm
    !
    ! -- drainable (Sy) over the salt column, residual-consistent chord. The
    !    residual contribution is +rate_s (the saltwater's own storage); the
    !    chord diagonal d(rate_s)/dh^s ~ -alphas*rho2 (interior) is applied only
    !    while the interface is within the cell.
    sc2 = SyCapacity(sto%dis%area(n), this%swi%sy(n))
    rho2 = sc2 * tled
    if (sto%integratechanges /= 0) then
      rho2old = SyCapacity(sto%dis%area(n), this%swi%oldsy(n)) * tled
    else
      rho2old = rho2
    end if
    rate_s = rho2old * tthk * ssold - rho2 * tthk * ssnew
    if (ssnew > DZERO .and. ssnew < DONE) then
      aterm = -rho2 * this%swi%alphas
    else
      aterm = DZERO
    end if
    rhsterm = aterm * h_new(n) - rate_s
    call matrix_sln%add_value_pos(idxglo(idiag), aterm)
    rhs(n) = rhs(n) + rhsterm
  end subroutine swisto_fc_salt

  !> @brief SWI STO formulation (saltwater model): storage Newton tangent for
  !! cell n. Upgrades the salt-column chord diagonal in swisto_fc_salt to the
  !! smoothed interface tangent (alphas chain rule). Newton only.
  !<
  subroutine swisto_fn_salt(this, n, matrix_sln, rhs, idxglo, h_new)
    use TdisModule, only: delt
    use GwfStorageUtilsModule, only: SyCapacity
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: h_new
    ! local
    type(GwfStoType), pointer :: sto
    integer(I4B) :: idiag
    real(DP) :: tled, tp, bt, tthk
    real(DP) :: sc2, rho2, ssnew, sderiv, a_fn
    !
    sto => this%sto
    tled = DONE / delt
    idiag = sto%dis%con%ia(n)
    tp = sto%dis%top(n)
    bt = sto%dis%bot(n)
    tthk = tp - bt
    ssnew = sQuadraticSaturation(tp, bt, this%swi%get_zetanew(n), sto%satomega)
    if (ssnew > DZERO .and. ssnew < DONE) then
      sc2 = SyCapacity(sto%dis%area(n), this%swi%sy(n))
      rho2 = sc2 * tled
      sderiv = sQuadraticSaturationDerivative(tp, bt, this%swi%get_zetanew(n), &
                                              sto%satomega)
      ! net salt tangent d(rate_s)/dh^s = -alphas*rho2*tthk*sderiv; the fc chord
      ! was -alphas*rho2, so the fn delta is alphas*rho2*(1 - tthk*sderiv)
      a_fn = this%swi%alphas * rho2 * (DONE - tthk * sderiv)
      call matrix_sln%add_value_pos(idxglo(idiag), a_fn)
      rhs(n) = rhs(n) + a_fn * h_new(n)
    end if
  end subroutine swisto_fn_salt

  !> @brief SWI STO formulation (saltwater model): storage rate for cell n. The
  !! salt-column storage is the saltwater model's own storage. It is reported
  !! through the SWI storage term (swi%storage) and flowja; the STO storage
  !! arrays (strgss/strgsy) are set to zero for the cell.
  !<
  subroutine swisto_cq_salt(this, n, flowja, h_new, h_old)
    use TdisModule, only: delt
    use GwfStorageUtilsModule, only: SsCapacity, SyCapacity, SsTerms
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(inout) :: flowja
    real(DP), dimension(:), intent(in) :: h_new
    real(DP), dimension(:), intent(in) :: h_old
    ! local
    type(GwfStoType), pointer :: sto
    integer(I4B) :: idiag
    real(DP) :: tled, tp, bt, tthk
    real(DP) :: sc1, rho1, rho1old, sc2, rho2, rho2old
    real(DP) :: ssnew, ssold, aterm, rhsterm, ratess, ratesy
    !
    sto => this%sto
    sto%strgss(n) = DZERO
    sto%strgsy(n) = DZERO
    if (sto%iss == 1) return
    if (this%swi%ibound(n) < 1) return
    tled = DONE / delt
    idiag = sto%dis%con%ia(n)
    tp = sto%dis%top(n)
    bt = sto%dis%bot(n)
    tthk = tp - bt
    ssnew = sQuadraticSaturation(tp, bt, this%swi%get_zetanew(n), sto%satomega)
    ssold = sQuadraticSaturation(tp, bt, this%swi%get_zetaold(n), sto%satomega)
    !
    ! -- compressible (Ss) rate
    sc1 = SsCapacity(sto%istor_coef, tp, bt, sto%dis%area(n), sto%ss(n))
    rho1 = sc1 * tled
    if (sto%integratechanges /= 0) then
      rho1old = SsCapacity(sto%istor_coef, tp, bt, sto%dis%area(n), &
                           this%swi%oldss(n)) * tled
    else
      rho1old = rho1
    end if
    call SsTerms(sto%iconvert(n), sto%iorig_ss, sto%iconf_ss, tp, bt, &
                 rho1, rho1old, ssnew, ssold, h_new(n), h_old(n), &
                 aterm, rhsterm, ratess)
    !
    ! -- drainable (Sy) rate over the salt column
    sc2 = SyCapacity(sto%dis%area(n), this%swi%sy(n))
    rho2 = sc2 * tled
    if (sto%integratechanges /= 0) then
      rho2old = SyCapacity(sto%dis%area(n), this%swi%oldsy(n)) * tled
    else
      rho2old = rho2
    end if
    ratesy = rho2old * tthk * ssold - rho2 * tthk * ssnew
    !
    ! -- the saltwater storage is reported through the always-on SWI budget term
    !    (matches how the freshwater interface storage is reported), so it does
    !    not depend on the STO iusesy flag
    this%swi%storage(n) = ratess + ratesy
    flowja(idiag) = flowja(idiag) + this%swi%storage(n)
  end subroutine swisto_cq_salt

  !> @brief SWI STO formulation: storage coefficients over all SWI cells.
  !!
  !! The formulation owns its traversal. A freshwater model leaves its cells
  !! unclaimed so the default storage formulation assembles the water column
  !! first and SWI adds the salt-column correction; a saltwater model claims
  !! its cells, so the default skips them and SWI supplies the salt column.
  !<
  subroutine swisto_fc(this, kiter, matrix_sln, rhs, idxglo, h_old, h_new)
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: h_old
    real(DP), dimension(:), intent(in) :: h_new
    ! local
    integer(I4B) :: n
    !
    call swisto_set_old_ss_sy(this)
    do n = 1, this%sto%dis%nodes
      if (.not. this%is_active(n)) cycle
      call swisto_fc_node(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
    end do
  end subroutine swisto_fc

  !> @brief SWI STO formulation: storage Newton terms over all SWI cells.
  !<
  subroutine swisto_fn(this, kiter, matrix_sln, rhs, idxglo, h_old, h_new)
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer, intent(inout) :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: idxglo
    real(DP), dimension(:), intent(in) :: h_old
    real(DP), dimension(:), intent(in) :: h_new
    ! local
    integer(I4B) :: n
    !
    call swisto_set_old_ss_sy(this)
    do n = 1, this%sto%dis%nodes
      if (.not. this%is_active(n)) cycle
      call swisto_fn_node(this, n, matrix_sln, rhs, idxglo, h_old, h_new)
    end do
  end subroutine swisto_fn

  !> @brief SWI STO formulation: storage flows over all SWI cells.
  !<
  subroutine swisto_cq(this, flowja, h_new, h_old)
    class(SwiStoFormulationType), intent(inout) :: this
    real(DP), dimension(:), intent(inout) :: flowja
    real(DP), dimension(:), intent(in) :: h_new
    real(DP), dimension(:), intent(in) :: h_old
    ! local
    integer(I4B) :: n
    !
    call swisto_set_old_ss_sy(this)
    do n = 1, this%sto%dis%nodes
      if (.not. this%is_active(n)) cycle
      call swisto_cq_node(this, n, flowja, h_new, h_old)
    end do
  end subroutine swisto_cq

  !> @brief SWI STO formulation: point at the STO Package's previous SS and SY
  !!
  !! When storage changes are integrated (TVS), STO keeps the previous time
  !! step's SS and SY in arrays it allocates on first use, during sto_rp, so
  !! the pointers are set here rather than in swi_ar.
  !<
  subroutine swisto_set_old_ss_sy(this)
    class(SwiStoFormulationType), intent(inout) :: this
    if (this%sto%integratechanges == 0) return
    if (.not. associated(this%swi%oldss)) then
      call mem_setptr(this%swi%oldss, 'OLDSS', this%sto%memoryPath)
    end if
    if (.not. associated(this%swi%oldsy)) then
      ! STO keeps a previous SY only when some cell is convertible; otherwise
      ! SY is used by SWI alone and does not change between time steps
      if (this%sto%iusesy == 1) then
        call mem_setptr(this%swi%oldsy, 'OLDSY', this%sto%memoryPath)
      else
        this%swi%oldsy => this%swi%sy
      end if
    end if
  end subroutine swisto_set_old_ss_sy

  !> @brief SWI STO formulation: add SWI storage to the model budget
  !<
  subroutine swisto_bd(this, isuppress_output, model_budget)
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType, rate_accumulator
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: isuppress_output
    type(BudgetType), intent(inout) :: model_budget
    real(DP) :: rin, rout
    call rate_accumulator(this%swi%storage, rin, rout)
    call model_budget%addentry(rin, rout, delt, budtxt(1), &
                               isuppress_output, '     SWI')
  end subroutine swisto_bd

  !> @brief SWI STO formulation: save SWI storage cell-by-cell flows
  !<
  subroutine swisto_save_flows(this, iprint, ibinun)
    class(SwiStoFormulationType), intent(inout) :: this
    integer(I4B), intent(in) :: iprint
    integer(I4B), intent(in) :: ibinun
    integer(I4B) :: nvaluesp, nwidthp
    character(len=1) :: cdatafmp = ' ', editdesc = ' '
    real(DP) :: dinact
    if (ibinun /= 0) then
      dinact = DZERO
      call this%swi%dis%record_array(this%swi%storage, this%swi%iout, iprint, &
                                     -ibinun, budtxt(1), cdatafmp, nvaluesp, &
                                     nwidthp, editdesc, dinact)
    end if
  end subroutine swisto_save_flows

end module GwfSwiModule
