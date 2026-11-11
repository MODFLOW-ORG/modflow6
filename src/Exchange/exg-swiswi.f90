!> @brief This module contains the SwiSwiExchangeModule Module
!!
!! This module contains the code for connecting freshwater and
!! saltwater GWF models.
!!
!<
module SwiSwiExchangeModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, LINELENGTH, DONE, DTWO, C3D_VERTICAL
  use SimModule, only: store_error, store_error_filename
  use SimVariablesModule, only: errmsg, model_loc_idx
  use MemoryHelperModule, only: create_mem_path
  use BaseModelModule, only: BaseModelType, GetBaseModelFromList
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use ListsModule, only: basemodellist, baseexchangelist
  use NumericalExchangeModule, only: NumericalExchangeType
  use GwfModule, only: GwfModelType
  use MatrixBaseModule
  use SmoothingModule, only: sQuadraticSaturation

  implicit none

  private
  public :: swiswi_cr
  public :: SwiSwiExchangeType

  !> @brief Head perturbation used for the central-difference derivatives of
  !! the saturated cell fraction in get_dsfdhs and get_dssdhf
  !<
  real(DP), parameter :: DEREPS = 1.D-10

  !> @brief Derived type for SwiSwiExchangeType
  !!
  !! This derived type contains information and methods for
  !! connecting freshwater and saltwater GWF models.
  !<
  type, extends(NumericalExchangeType) :: SwiSwiExchangeType

    ! model pointers
    class(GwfModelType), pointer :: gwf_fresh => null() !< freshwater GWF model
    class(GwfModelType), pointer :: gwf_salt => null() !< saltwater GWF model

    ! options
    character(len=LINELENGTH), pointer :: filename => null() !< name of the input file
    integer(I4B), pointer :: inocrossstorage => null() !< dev flag: disable cross-fluid storage newton terms
    integer(I4B), pointer :: inocrossflow => null() !< dev flag: disable cross-fluid flow newton terms

    ! number of connections
    integer(I4B), pointer :: nexg => null() !< number of connections (number of cells in fresh model or salt model)

    ! matrix position index arrays
    integer(I4B), dimension(:), pointer, contiguous :: idxglo => null() !< mapping to global (solution) amat
    integer(I4B), dimension(:), pointer, contiguous :: idxsymglo => null() !< mapping to global (solution) symmetric amat
    integer(I4B), dimension(:), pointer, contiguous :: idxjasalt => null() !< mapping for freshwater nodes to saltwater ja positions
    integer(I4B), dimension(:), pointer, contiguous :: idxjafresh => null() !< mapping for saltwater nodes to freshwater ja positions

  contains

    procedure :: exg_df => swi_swi_df
    procedure :: exg_ac => swi_swi_ac
    procedure :: exg_mc => swi_swi_mc
    procedure :: exg_ar => swi_swi_ar
    procedure :: exg_fc => swi_swi_fc
    procedure :: exg_fn => swi_swi_fn
    procedure :: exg_da => swi_swi_da
    procedure :: swi_fn_cross_storage
    procedure :: swi_fn_cross_flow
    procedure :: swi_cross_storage
    procedure, private :: same_grid
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: source_options
    procedure :: get_dsfdhs
    procedure :: get_dssdhf
    procedure :: connects_model => swi_swi_connects_model

  end type SwiSwiExchangeType

contains

  !> @brief Create SWI-SWI Exchange
  !!
  !! Create a new SWI-SWI Exchange object, add it to the list of exchanges,
  !! and set pointers to the freshwater and saltwater GWF models.
  !<
  subroutine swiswi_cr(filename, name, id, m1_id, m2_id, input_mempath)
    ! dummy
    character(len=*), intent(in) :: filename !< filename for reading
    character(len=*) :: name !< exchange name
    integer(I4B), intent(in) :: id !< id for the exchange
    integer(I4B), intent(in) :: m1_id !< id for model 1
    integer(I4B), intent(in) :: m2_id !< id for model 2
    character(len=*), intent(in) :: input_mempath
    ! local
    type(SwiSwiExchangeType), pointer :: exchange
    class(BaseModelType), pointer :: mb
    class(BaseExchangeType), pointer :: baseexchange
    integer(I4B) :: m1_index, m2_index

    ! Create a new exchange and add it to the baseexchangelist container
    allocate (exchange)
    baseexchange => exchange
    call AddBaseExchangeToList(baseexchangelist, baseexchange)

    ! Assign id and name
    exchange%id = id
    exchange%name = name
    exchange%memoryPath = create_mem_path(exchange%name)
    exchange%input_mempath = input_mempath

    ! allocate scalars and set defaults
    call exchange%allocate_scalars()
    exchange%filename = filename
    exchange%typename = 'SWI-SWI'

    ! -- set gwf_fresh
    m1_index = model_loc_idx(m1_id)
    if (m1_index > 0) then
      mb => GetBaseModelFromList(basemodellist, m1_index)
      select type (mb)
      type is (GwfModelType)
        exchange%gwf_fresh => mb
      end select
    end if
    !
    ! -- set gwf_salt
    m2_index = model_loc_idx(m2_id)
    if (m2_index > 0) then
      mb => GetBaseModelFromList(basemodellist, m2_index)
      select type (mb)
      type is (GwfModelType)
        exchange%gwf_salt => mb
      end select
    end if
    !
    ! -- Verify that gwf model1 is of the correct type
    if (.not. associated(exchange%gwf_fresh) .and. m1_index > 0) then
      write (errmsg, '(3a)') 'Problem with SWI-SWI Exchange ', &
        trim(exchange%name), &
        '.  First specified GWF Model does not appear to be of the &
        &correct type.'
      call store_error(errmsg, terminate=.true.)
    end if
    !
    ! -- Verify that gwf model2 is of the correct type
    if (.not. associated(exchange%gwf_salt) .and. m2_index > 0) then
      write (errmsg, '(3a)') 'Problem with SWI-SWI Exchange ', &
        trim(exchange%name), &
        '.  Second specified GWF Model does not appear to be of the &
        &correct type.'
      call store_error(errmsg, terminate=.true.)
    end if

  end subroutine swiswi_cr

  !> @ brief Define SWI-SWI Exchange
  !!
  !! Define SWI-SWI Exchange object.
  !<
  subroutine swi_swi_df(this)
    ! modules
    use SimVariablesModule, only: iout
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    ! local

    ! log the exchange
    write (iout, '(/a,a)') ' Creating exchange: ', this%name

    ! -- Ensure models are in same solution
    if (associated(this%gwf_fresh) .and. associated(this%gwf_salt)) then
      if (this%gwf_fresh%idsoln /= this%gwf_salt%idsoln) then
        call store_error('Two models are connected in a SWI-SWI Exchange '// &
                         'but they are in different solutions. The '// &
                         'freshwater and saltwater GWF models must be in '// &
                         'the same solution: '// &
                         trim(this%gwf_fresh%name)//' '// &
                         trim(this%gwf_salt%name))
        call store_error_filename(this%filename)
      end if
    end if

    ! Ensure fresh model has active SWI Package
    if (this%gwf_fresh%inswi == 0) then
      call store_error( &
        'A SWI-SWI Exchange is active, but the freshwater GWF model &
        &does not have an active Seawater Intrusion (SWI) Package.  Activate &
        &the SWI Package for the freshwater GWF model.')
      call store_error_filename(this%filename)
    end if

    ! Ensure salt model has active SWI Package
    if (this%gwf_salt%inswi == 0) then
      call store_error( &
        'A SWI-SWI Exchange is active, but the saltwater GWF model &
        &does not have an active Seawater Intrusion (SWI) Package.  Activate &
        &the SWI Package for the saltwater GWF model.')
      call store_error_filename(this%filename)
    end if

    ! Set the configuration for each SWI package
    call this%gwf_fresh%swi%set_configuration("FRESHWATER_TWO_FLUID")
    call this%gwf_salt%swi%set_configuration("SALTWATER_TWO_FLUID")

    ! Check to make sure fresh and salt models have same number of nodes
    if (this%gwf_fresh%dis%nodes /= this%gwf_salt%dis%nodes) then
      call store_error( &
        'A SWI-SWI Exchange is active, but the freshwater GWF model &
        &does not have the same number of nodes as the saltwater &
        &GWF model.  Both models must have the same discretization &
        &properties.')
      call store_error_filename(this%filename)
    end if

    ! The exchange pairs cell n and connection ipos of one model with the same
    ! cell and connection of the other, so the two grids must be identical
    if (this%gwf_fresh%dis%nodes == this%gwf_salt%dis%nodes) then
      if (.not. this%same_grid()) then
        call store_error( &
          'A SWI-SWI Exchange is active, but the freshwater and saltwater &
          &GWF models do not have the same grid. Both models must have the &
          &same cell connections, cell tops, cell bottoms, and cell areas.')
        call store_error_filename(this%filename)
      end if
    end if

    ! Set the number of exchanges equal to the number of equations in
    ! the fresh or saltwater models
    this%nexg = this%gwf_fresh%dis%nodes

    ! source options
    call this%source_options(iout)

    ! the cross-flow Newton terms differentiate the upstream cell's fraction;
    ! without upstream weighting (SWI dev option) the flow term is Picard in
    ! both models, so drop the cross-flow terms and their sparsity
    if (this%gwf_fresh%swi%inoupstream /= 0 .or. &
        this%gwf_salt%swi%inoupstream /= 0) then
      if (this%inocrossflow == 0) then
        this%inocrossflow = 1
        write (iout, '(4x,a)') &
          'Cross-fluid flow Newton terms are DISABLED because a SWI Package &
          &runs without upstream weighting (dev option).'
      end if
    end if

    ! allocate arrays
    call this%allocate_arrays()

    return
  end subroutine swi_swi_df

  !> @ brief Add connections
  !!
  !! Add the freshwater-saltwater connections to the solution sparse matrix.
  !<
  subroutine swi_swi_ac(this, sparse)
    ! modules
    use SparseModule, only: sparsematrix
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    type(sparsematrix), intent(inout) :: sparse
    ! local
    integer(I4B) :: n, iglo, jglo
    integer(I4B) :: ipos, m

    ! add same-cell fresh<->salt connections (used by cross-storage and by
    ! cross-flow's self-upstream term); skip only if both are disabled
    if (this%inocrossstorage == 0 .or. this%inocrossflow == 0) then
      do n = 1, this%nexg
        iglo = n + this%gwf_fresh%moffset
        jglo = n + this%gwf_salt%moffset
        call sparse%addconnection(iglo, jglo, 1)
        call sparse%addconnection(jglo, iglo, 1)
      end do
    end if

    ! add cross-flow neighbor connections (fresh<->salt) only if cross-flow is
    ! active, so a disabled cross-flow does not expand the matrix sparsity
    if (this%inocrossflow == 0) then
      ! fresh node to salt nodes
      do n = 1, this%gwf_fresh%dis%nodes
        iglo = n + this%gwf_fresh%moffset
        do ipos = this%gwf_salt%dis%con%ia(n), this%gwf_salt%dis%con%ia(n + 1) - 1
          m = this%gwf_salt%dis%con%ja(ipos)
          jglo = m + this%gwf_salt%moffset
          call sparse%addconnection(iglo, jglo, 1)
        end do
      end do
      ! salt node to fresh nodes
      do n = 1, this%gwf_salt%dis%nodes
        iglo = n + this%gwf_salt%moffset
        do ipos = this%gwf_fresh%dis%con%ia(n), &
          this%gwf_fresh%dis%con%ia(n + 1) - 1
          m = this%gwf_fresh%dis%con%ja(ipos)
          jglo = m + this%gwf_fresh%moffset
          call sparse%addconnection(iglo, jglo, 1)
        end do
      end do
    end if

  end subroutine swi_swi_ac

  !> @ brief Map connections
  !!
  !! Map the connections in the global matrix
  !<
  subroutine swi_swi_mc(this, matrix_sln)
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    class(MatrixBaseType), pointer :: matrix_sln !< the system matrix
    ! local
    integer(I4B) :: n, iglo, jglo
    integer(I4B) :: ipos, m

    ! map same-cell fresh<->salt connections (must match swi_swi_ac inclusion)
    if (this%inocrossstorage == 0 .or. this%inocrossflow == 0) then
      do n = 1, this%nexg
        iglo = n + this%gwf_fresh%moffset
        jglo = n + this%gwf_salt%moffset
        this%idxglo(n) = matrix_sln%get_position(iglo, jglo)
        this%idxsymglo(n) = matrix_sln%get_position(jglo, iglo)
      end do
    end if

    ! map cross-flow neighbor connections (must match swi_swi_ac inclusion)
    if (this%inocrossflow == 0) then
      ! fresh node to salt nodes, indexed by connection position
      do n = 1, this%gwf_fresh%dis%nodes
        iglo = n + this%gwf_fresh%moffset
        do ipos = this%gwf_salt%dis%con%ia(n), this%gwf_salt%dis%con%ia(n + 1) - 1
          m = this%gwf_salt%dis%con%ja(ipos)
          jglo = m + this%gwf_salt%moffset
          this%idxjasalt(ipos) = matrix_sln%get_position(iglo, jglo)
        end do
      end do
      ! salt node to fresh nodes, indexed by connection position
      do n = 1, this%gwf_salt%dis%nodes
        iglo = n + this%gwf_salt%moffset
        do ipos = this%gwf_fresh%dis%con%ia(n), &
          this%gwf_fresh%dis%con%ia(n + 1) - 1
          m = this%gwf_fresh%dis%con%ja(ipos)
          jglo = m + this%gwf_fresh%moffset
          this%idxjafresh(ipos) = matrix_sln%get_position(iglo, jglo)
        end do
      end do
    end if

  end subroutine swi_swi_mc

  !> @ brief Allocate and read
  !<
  subroutine swi_swi_ar(this)
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType

    ! The interface elevation is computed from both heads using coefficients
    ! derived from the relative density difference. The two models must
    ! therefore agree on it, or each would place the interface differently and
    ! the coupling would be inconsistent without any outward sign of trouble.
    if (this%gwf_fresh%swi%reldensdiff /= this%gwf_salt%swi%reldensdiff) then
      write (errmsg, '(a,a,a,G0,a,a,a,G0,a)') &
        'RELATIVE_DENSITY_DIFFERENCE must be the same in both models of a &
        &two-fluid simulation. Model ', trim(this%gwf_fresh%name), &
        ' specified ', this%gwf_fresh%swi%reldensdiff, &
        ' and model ', trim(this%gwf_salt%name), &
        ' specified ', this%gwf_salt%swi%reldensdiff, '.'
      call store_error(errmsg, terminate=.false.)
      call store_error_filename(this%filename)
    end if

    ! The cross-storage terms use the freshwater model's specific yield for
    ! both fluids, so the two models must agree on it
    if (this%gwf_fresh%insto > 0 .and. this%gwf_salt%insto > 0) then
      if (any(this%gwf_fresh%sto%sy /= this%gwf_salt%sto%sy)) then
        write (errmsg, '(a,a,a,a,a)') &
          'Specific yield (SY) must be the same in both models of a &
          &two-fluid simulation. Models ', trim(this%gwf_fresh%name), &
          ' and ', trim(this%gwf_salt%name), ' differ.'
        call store_error(errmsg, terminate=.false.)
        call store_error_filename(this%filename)
      end if
    end if

    ! Set head pointers for fresh and salt model swi packages.
    call this%gwf_fresh%swi%set_head_pointers(this%gwf_fresh%x, &
                                              this%gwf_salt%x, &
                                              this%gwf_fresh%xold, &
                                              this%gwf_salt%xold)
    call this%gwf_salt%swi%set_head_pointers(this%gwf_fresh%x, &
                                             this%gwf_salt%x, &
                                             this%gwf_fresh%xold, &
                                             this%gwf_salt%xold)
  end subroutine swi_swi_ar

  !> @ brief Fill coefficients
  !!
  !! The SWI-SWI Exchange contributes only Newton (Jacobian) cross terms; the
  !! Picard cross-fluid coupling enters each model through the interface, which
  !! each model computes in its own storage/flow fill from the other model's
  !! head at the current iterate. Dispatch the Newton terms to swi_swi_fn.
  !<
  subroutine swi_swi_fc(this, kiter, matrix_sln, rhs_sln, inwtflag)
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs_sln
    integer(I4B), optional, intent(in) :: inwtflag
    ! local
    integer(I4B) :: inwt

    ! -- set inwt to exchange newton, but shut off if requested by caller
    inwt = this%gwf_fresh%inewton
    if (present(inwtflag)) then
      if (inwtflag == 0) inwt = 0
    end if
    if (inwt /= 0) then
      call this%exg_fn(kiter, matrix_sln)
    end if
  end subroutine swi_swi_fc

  !> @ brief Fill Newton (Jacobian) cross terms for the SWI-SWI Exchange
  !!
  !! Freshwater flow/storage depend on the saltwater head (and vice versa)
  !! through the interface zeta = -alphaf*hf + alphas*hs. These are the
  !! derivative terms coupling each fluid's equation to the other fluid's
  !! head; they cancel at convergence (the Picard residual in the model's own
  !! fc uses the other fluid's head at the current iterate). Horizontal
  !! connections only; vertical flow is Picard.
  !<
  subroutine swi_swi_fn(this, kiter, matrix_sln)
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln

    ! -- cross-storage Jacobian is transient only; skip in a steady-state stress
    !    period. The cross-flow Jacobian below is not storage and applies
    !    always. Either can be turned off with a dev option
    !    (inocrossstorage/inocrossflow).
    if (this%inocrossstorage == 0 .and. this%gwf_fresh%iss == 0) then
      call this%swi_fn_cross_storage(matrix_sln)
    end if

    if (this%inocrossflow == 0) then
      call this%swi_fn_cross_flow(matrix_sln)
    end if

  end subroutine swi_swi_fn

  !> @ brief Deallocate
  !!
  !! Deallocate memory associated with this object
  !<
  subroutine swi_swi_da(this)
    ! modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    ! dummy
    class(SwiSwiExchangeType) :: this !< SwiSwiExchangeType

    ! deallocate IDM memory
    call memorystore_remove(this%name, '', idm_context)

    ! arrays
    call mem_deallocate(this%idxglo)
    call mem_deallocate(this%idxsymglo)
    call mem_deallocate(this%idxjasalt)
    call mem_deallocate(this%idxjafresh)

    ! scalars
    deallocate (this%filename)
    call mem_deallocate(this%inocrossstorage)
    call mem_deallocate(this%inocrossflow)
    call mem_deallocate(this%nexg)

  end subroutine swi_swi_da

  !> @ brief Cross-fluid STORAGE Newton terms
  !!
  !! d(freshwater storage)/d(hs) and d(saltwater storage)/d(hf) from interface
  !! movement, on the same-cell fresh<->salt coupling positions.
  !<
  subroutine swi_fn_cross_storage(this, matrix_sln)
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    class(MatrixBaseType), pointer :: matrix_sln
    ! local
    integer(I4B) :: n
    real(DP) :: termf, terms, rtermf, rterms

    do n = 1, this%nexg
      termf = DZERO
      rtermf = DZERO
      terms = DZERO
      rterms = DZERO

      ! calculate cross storage terms
      call this%swi_cross_storage(n, termf, rtermf, terms, rterms)

      ! fill off-diagonal matrix coefficients
      call matrix_sln%add_value_pos(this%idxglo(n), termf)
      call matrix_sln%add_value_pos(this%idxsymglo(n), terms)

      ! update rhs for the fresh and salt model nodes
      this%gwf_fresh%rhs(n) = this%gwf_fresh%rhs(n) + rtermf
      this%gwf_salt%rhs(n) = this%gwf_salt%rhs(n) + rterms
    end do

  end subroutine swi_fn_cross_storage

  !> @ brief Cross-fluid FLOW Newton terms
  !!
  !! d(freshwater flow)/d(hs) and d(saltwater flow)/d(hf) from the
  !! interface-dependent zone conductance. Horizontal connections only (vertical
  !! flow is Picard). For each connection the derivative is evaluated at the
  !! upstream cell and lands on that cell's fresh<->salt coupling position.
  !<
  subroutine swi_fn_cross_flow(this, matrix_sln)
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    class(MatrixBaseType), pointer :: matrix_sln
    ! local
    integer(I4B) :: n, ipos, m, ihc
    real(DP) :: termf, terms, csat, q

    ! The two models share the same grid (checked in swi_swi_df), so ipos
    ! indexes the same connection in both, and idxjasalt and idxjafresh are
    ! indexed by ipos. The diagonal position ia(n) is skipped.

    ! FRESHWATER EQUATIONS: the saltwater head in a connected cell affects the
    ! freshwater zone conductance, so it enters the freshwater equation.
    associate (con => this%gwf_fresh%dis%con, ibound => this%gwf_fresh%ibound, &
               hf => this%gwf_fresh%x, hs => this%gwf_salt%x)
      do n = 1, this%gwf_fresh%dis%nodes
        if (ibound(n) == 0) cycle
        do ipos = con%ia(n) + 1, con%ia(n + 1) - 1
          m = con%ja(ipos)
          if (ibound(m) == 0) cycle
          ihc = con%ihc(con%jas(ipos))
          if (ihc == C3D_VERTICAL) cycle
          csat = this%gwf_fresh%npf%condsat(con%jas(ipos))
          q = csat * (hf(m) - hf(n))
          if (hf(m) > hf(n)) then
            ! m is upstream: hs(m) affects the conductance
            termf = this%get_dsfdhs(m) * q
            call matrix_sln%add_value_pos(this%idxjasalt(ipos), termf)
            this%gwf_fresh%rhs(n) = this%gwf_fresh%rhs(n) + termf * hs(m)
          else
            ! n is upstream: hs(n) affects the conductance
            termf = this%get_dsfdhs(n) * q
            call matrix_sln%add_value_pos(this%idxglo(n), termf)
            this%gwf_fresh%rhs(n) = this%gwf_fresh%rhs(n) + termf * hs(n)
          end if
        end do
      end do
    end associate

    ! SALTWATER EQUATIONS: the freshwater head in a connected cell affects the
    ! saltwater zone conductance, so it enters the saltwater equation.
    associate (con => this%gwf_salt%dis%con, ibound => this%gwf_salt%ibound, &
               hf => this%gwf_fresh%x, hs => this%gwf_salt%x)
      do n = 1, this%gwf_salt%dis%nodes
        if (ibound(n) == 0) cycle
        do ipos = con%ia(n) + 1, con%ia(n + 1) - 1
          m = con%ja(ipos)
          if (ibound(m) == 0) cycle
          ihc = con%ihc(con%jas(ipos))
          if (ihc == C3D_VERTICAL) cycle
          csat = this%gwf_salt%npf%condsat(con%jas(ipos))
          q = csat * (hs(m) - hs(n))
          if (hs(m) > hs(n)) then
            ! m is upstream: hf(m) affects the conductance
            terms = this%get_dssdhf(m) * q
            call matrix_sln%add_value_pos(this%idxjafresh(ipos), terms)
            this%gwf_salt%rhs(n) = this%gwf_salt%rhs(n) + terms * hf(m)
          else
            ! n is upstream: hf(n) affects the conductance
            terms = this%get_dssdhf(n) * q
            call matrix_sln%add_value_pos(this%idxsymglo(n), terms)
            this%gwf_salt%rhs(n) = this%gwf_salt%rhs(n) + terms * hf(n)
          end if
        end do
      end do
    end associate

  end subroutine swi_fn_cross_flow

  !> @brief Calculate cross-fluid storage coefficients for a cell
  !!
  !! Calculate the matrix coefficients and right-hand-side contributions that
  !! couple freshwater storage in cell n to the saltwater head, and saltwater
  !! storage in cell n to the freshwater head, through interface movement.
  !<
  subroutine swi_cross_storage(this, n, termf, rtermf, terms, rterms)
    ! modules
    use TdisModule, only: delt
    use GwfStorageUtilsModule, only: SyCapacity
    !
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    ! dummy variables
    integer(I4B), intent(in) :: n
    real(DP), intent(out) :: termf
    real(DP), intent(out) :: rtermf
    real(DP), intent(out) :: terms
    real(DP), intent(out) :: rterms
    ! local variables
    real(DP) :: tled
    real(DP) :: sc2
    real(DP) :: rho2
    real(DP) :: tp
    real(DP) :: bt
    real(DP) :: tthk
    real(DP) :: hf, hs
    real(DP) :: dssdhf
    real(DP) :: dsfdhs

    termf = DZERO
    rtermf = DZERO
    terms = DZERO
    rterms = DZERO

    ! set variables
    tled = DONE / delt

    if (this%gwf_fresh%ibound(n) <= 0) return

    ! aquifer elevations and thickness
    tp = this%gwf_fresh%dis%top(n)
    bt = this%gwf_fresh%dis%bot(n)
    tthk = tp - bt
    hf = this%gwf_fresh%x(n)
    hs = this%gwf_salt%x(n)

    ! storage coefficients
    sc2 = SyCapacity(this%gwf_fresh%dis%area(n), this%gwf_fresh%sto%sy(n))
    rho2 = sc2 * tled

    dsfdhs = this%get_dsfdhs(n)
    termf = -rho2 * tthk * dsfdhs
    rtermf = termf * hs

    dssdhf = this%get_dssdhf(n)
    terms = -rho2 * tthk * dssdhf
    rterms = terms * hf

  end subroutine swi_cross_storage

  !> @brief Check that the fresh and salt models share the same grid
  !<
  function same_grid(this) result(same)
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    logical(LGP) :: same
    ! local
    integer(I4B) :: n, ipos
    !
    same = .false.
    associate (dfr => this%gwf_fresh%dis, dsa => this%gwf_salt%dis)
      if (dfr%con%nja /= dsa%con%nja) return
      do n = 1, dfr%nodes + 1
        if (dfr%con%ia(n) /= dsa%con%ia(n)) return
      end do
      do ipos = 1, dfr%con%nja
        if (dfr%con%ja(ipos) /= dsa%con%ja(ipos)) return
      end do
      do n = 1, dfr%nodes
        if (dfr%top(n) /= dsa%top(n)) return
        if (dfr%bot(n) /= dsa%bot(n)) return
        if (dfr%area(n) /= dsa%area(n)) return
      end do
    end associate
    same = .true.
  end function same_grid

  !> @ brief Allocate scalars
  !!
  !! Allocate scalar variables
  !<
  subroutine allocate_scalars(this)
    ! modules
    use MemoryManagerModule, only: mem_allocate
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType

    allocate (this%filename)
    this%filename = ''

    call mem_allocate(this%inocrossstorage, 'INOCROSSSTORAGE', this%memoryPath)
    call mem_allocate(this%inocrossflow, 'INOCROSSFLOW', this%memoryPath)
    call mem_allocate(this%nexg, 'NEXG', this%memoryPath)

    this%inocrossstorage = 0
    this%inocrossflow = 0
    this%nexg = 0

  end subroutine allocate_scalars

  !> @ brief Allocate arrays
  !!
  !! Allocate arrays
  !<
  subroutine allocate_arrays(this)
    ! modules
    use MemoryManagerModule, only: mem_allocate
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    ! local
    integer(I4B) :: i
    integer(I4B) :: idxsize, nglo

    ! -- idxglo/idxsymglo (same-cell fresh<->salt) are used by cross-storage AND
    !    by cross-flow's self-upstream term, so size them if either is active;
    !    idxjasalt/idxjafresh (neighbor fresh<->salt) are cross-flow only. A
    !    disabled term is sized 0 so it expands neither the matrix nor memory.
    nglo = 0
    if (this%inocrossstorage == 0 .or. this%inocrossflow == 0) nglo = this%nexg
    idxsize = 0
    if (this%inocrossflow == 0) idxsize = this%gwf_fresh%dis%con%nja

    call mem_allocate(this%idxglo, nglo, 'IDXGLO', this%memoryPath)
    call mem_allocate(this%idxsymglo, nglo, 'IDXSYMGLO', this%memoryPath)
    call mem_allocate(this%idxjasalt, idxsize, 'IDXJASALT', this%memoryPath)
    call mem_allocate(this%idxjafresh, idxsize, 'IDXJAFRESH', this%memoryPath)

    ! Initialize
    do i = 1, nglo
      this%idxglo(i) = 0
      this%idxsymglo(i) = 0
    end do
    do i = 1, idxsize
      this%idxjasalt(i) = 0
      this%idxjafresh(i) = 0
    end do
    !
    ! -- Return
    return
  end subroutine allocate_arrays

  !> @ brief Source options
  !!
  !! Source the options block
  !<
  subroutine source_options(this, iout)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use ExgSwiswiInputModule, only: ExgSwiswiParamFoundType
    ! -- dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    integer(I4B), intent(in) :: iout
    ! -- local
    type(ExgSwiswiParamFoundType) :: found
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%inocrossstorage, 'INOCROSSSTORAGE', &
                       this%input_mempath, found%inocrossstorage)
    call mem_set_value(this%inocrossflow, 'INOCROSSFLOW', &
                       this%input_mempath, found%inocrossflow)
    !
    write (iout, '(1x,a)') 'Processing SWI-SWI Exchange options'
    !
    if (found%inocrossstorage) then
      write (iout, '(4x,a)') &
        'Cross-fluid storage Newton terms are DISABLED (dev option).'
    end if
    !
    if (found%inocrossflow) then
      write (iout, '(4x,a)') &
        'Cross-fluid flow Newton terms are DISABLED (dev option).'
    end if
    !
    write (iout, '(1x,a)') 'End processing SWI-SWI Exchange options'
    !
    ! -- Return
    return
  end subroutine source_options

  !> @ brief Calculate dSf/dhs for the cross-fluid Newton terms
  !!
  !! Calculate the derivative of the freshwater saturated cell fraction with
  !! respect to the saltwater head in cell n. The derivative is computed as a
  !! central difference of the saturated cell fraction for saltwater head
  !! perturbations of plus and minus DEREPS. A central difference is used
  !! instead of the analytical derivative of the quadratic saturated cell
  !! fraction function because it gives more robust convergence of the
  !! two-fluid Newton iterations.
  !<
  function get_dsfdhs(this, n) result(dsfdhs)
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    integer(I4B), intent(in) :: n !< node number
    real(DP) :: dsfdhs !< derivative of Sf with respect to hs
    ! local
    real(DP) :: tp !< top elevation
    real(DP) :: bt !< bottom elevation
    real(DP) :: z !< zeta with negative perturbation
    real(DP) :: zp !< zeta with positive perturbation
    real(DP) :: ss !< saltwater saturated cell fraction at z
    real(DP) :: ssp !< saltwater saturated cell fraction at zp

    tp = this%gwf_fresh%dis%top(n)
    bt = this%gwf_fresh%dis%bot(n)
    z = this%gwf_fresh%swi%get_zetanew(n, eps_salt=-DEREPS)
    zp = this%gwf_fresh%swi%get_zetanew(n, eps_salt=DEREPS)
    ss = sQuadraticSaturation(tp, bt, z, this%gwf_fresh%sto%satomega)
    ssp = sQuadraticSaturation(tp, bt, zp, this%gwf_fresh%sto%satomega)
    ! sf = 1 - ss, sfp = 1 - ssp; sfp - sf = ss - ssp
    dsfdhs = (ss - ssp) / (DTWO * DEREPS)

  end function get_dsfdhs

  !> @ brief Calculate dSs/dhf for the cross-fluid Newton terms
  !!
  !! Calculate the derivative of the saltwater saturated cell fraction with
  !! respect to the freshwater head in cell n. The derivative is computed as a
  !! central difference of the saturated cell fraction for freshwater head
  !! perturbations of plus and minus DEREPS. A central difference is used for
  !! the same reason as in get_dsfdhs.
  !<
  function get_dssdhf(this, n) result(dssdhf)
    ! dummy
    class(SwiSwiExchangeType) :: this !<  SwiSwiExchangeType
    integer(I4B), intent(in) :: n !< node number
    real(DP) :: dssdhf !< derivative of Ss with respect to hf
    ! local
    real(DP) :: tp !< top elevation
    real(DP) :: bt !< bottom elevation
    real(DP) :: z !< zeta with negative perturbation
    real(DP) :: zp !< zeta with positive perturbation
    real(DP) :: ss !< saltwater saturated cell fraction at z
    real(DP) :: ssp !< saltwater saturated cell fraction at zp

    tp = this%gwf_fresh%dis%top(n)
    bt = this%gwf_fresh%dis%bot(n)
    z = this%gwf_fresh%swi%get_zetanew(n, eps_fresh=-DEREPS)
    zp = this%gwf_fresh%swi%get_zetanew(n, eps_fresh=DEREPS)
    ss = sQuadraticSaturation(tp, bt, z, this%gwf_fresh%sto%satomega)
    ssp = sQuadraticSaturation(tp, bt, zp, this%gwf_fresh%sto%satomega)
    dssdhf = (ssp - ss) / (DTWO * DEREPS)

  end function get_dssdhf

  !> @brief Should return true when the exchange should be added to the
  !! solution where the model resides
  !<
  function swi_swi_connects_model(this, model) result(is_connected)
    ! dummy
    class(SwiSwiExchangeType) :: this !< the instance of the exchange
    class(BaseModelType), pointer, intent(in) :: model !< the model to which the exchange might hold a connection
    ! return
    logical(LGP) :: is_connected !< true, when connected

    is_connected = .false.
    select type (model)
    class is (GwfModelType)
      if (associated(this%gwf_fresh, model)) then
        is_connected = .true.
      else if (associated(this%gwf_salt, model)) then
        is_connected = .true.
      end if
    end select

  end function swi_swi_connects_model

end module SwiSwiExchangeModule
