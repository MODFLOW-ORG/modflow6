module GhostNodeModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH
  use NumericalModelModule, only: NumericalModelType
  use NumericalPackageModule, only: NumericalPackageType
  use BlockParserModule, only: BlockParserType
  use MatrixBaseModule

  implicit none

  private
  public :: GhostNodeType
  public :: gnc_cr

  type, extends(NumericalPackageType) :: GhostNodeType
    logical, pointer :: smgnc => null() ! single model gnc
    logical, pointer :: implicit => null() ! lhs or rhs
    logical, pointer :: i2kn => null() ! not used
    integer(I4B), pointer :: nexg => null() ! number of gncs
    integer(I4B), pointer :: numjs => null() ! number of connecting nodes
    class(NumericalModelType), pointer :: m1 => null() ! pointer to model 1
    class(NumericalModelType), pointer :: m2 => null() ! pointer to model 2
    integer(I4B), dimension(:), pointer, contiguous :: nodem1 => null() ! array of nodes in model 1
    integer(I4B), dimension(:), pointer, contiguous :: nodem2 => null() ! array of nodes in model 2
    integer(I4B), dimension(:, :), pointer, contiguous :: nodesj => null() ! array of interpolation nodes
    real(DP), dimension(:), pointer, contiguous :: cond => null() ! array of conductance
    integer(I4B), dimension(:), pointer, contiguous :: idxglo => null() ! connection position in amat
    integer(I4B), dimension(:), pointer, contiguous :: idxsymglo => null() ! symmetric position in amat
    real(DP), dimension(:, :), pointer, contiguous :: alphasj => null() ! interpolation factors
    integer(I4B), dimension(:), pointer, contiguous :: idiagn => null() ! amat diagonal position of n
    integer(I4B), dimension(:), pointer, contiguous :: idiagm => null() ! amat diagonal position of m
    integer(I4B), dimension(:, :), pointer, contiguous :: jposinrown => null() ! amat j position in row n
    integer(I4B), dimension(:, :), pointer, contiguous :: jposinrowm => null() ! amat j position in row m

  contains

    procedure :: gnc_df
    procedure :: gnc_ac
    procedure :: gnc_mc
    procedure, private :: gnc_fmsav
    procedure :: gnc_fc
    procedure :: gnc_fn
    procedure :: gnc_cq
    procedure :: gnc_ot
    procedure :: gnc_da
    procedure :: deltaQgnc
    procedure :: allocate_scalars
    procedure, private :: allocate_arrays
    procedure, private :: source_options
    procedure, private :: source_dimensions
    procedure, private :: source_data
  end type GhostNodeType

contains

  !> @brief Create new GNC exchange object
  !<
  subroutine gnc_cr(gncobj, name_parent, input_mempath, inunit, iout)
    ! -- dummy
    type(GhostNodeType), pointer, intent(inout) :: gncobj
    character(len=*), intent(in) :: name_parent
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    !
    ! -- Allocate the gnc exchange object
    allocate (gncobj)
    !
    ! -- create name and memory path. name_parent will either be model name or the
    !    exchange name.
    call gncobj%set_names(1, name_parent, 'GNC', 'GNC', input_mempath)
    !
    ! -- allocate scalars
    call gncobj%allocate_scalars()
    !
    ! -- Set variables
    gncobj%inunit = inunit
    gncobj%iout = iout
  end subroutine gnc_cr

  !> @brief Initialize a gnc object
  !<
  subroutine gnc_df(this, m1, m2)
    ! -- modules
    use NumericalModelModule, only: NumericalModelType
    use SimModule, only: store_error, store_error_unit
    use SimVariablesModule, only: errmsg
    ! -- dummy
    class(GhostNodeType) :: this
    class(NumericalModelType), target :: m1
    class(NumericalModelType), target, optional :: m2
    !
    ! -- Point or set attributes
    this%m1 => m1
    this%m2 => m1
    !
    ! -- If m2 is present, then GNC spans two models
    if (present(m2)) then
      this%m2 => m2
      this%smgnc = .false.
    end if
    !
    ! -- Initialize block parser
    call this%parser%Initialize(this%inunit, this%iout)
    !
    ! -- read gnc options
    call this%source_options()
    !
    ! -- read gnc dimensions
    call this%source_dimensions()
    !
    ! -- allocate arrays
    call this%allocate_arrays()
    !
    ! -- Allocate and read the gnc entries
    call this%source_data()
    !
    ! -- Trap for implicit gnc but models are in different solutions
    if (this%m1%idsoln /= this%m2%idsoln) then
      if (this%implicit) then
        write (errmsg, '(a)') 'GNC is implicit but models are in '// &
          'different solutions.'
        call store_error(errmsg)
        call store_error_unit(this%inunit)
      end if
    end if
  end subroutine gnc_df

  !> @brief Single or Two-Model GNC Add Connections
  !!
  !! For implicit GNC, expand the sparse solution matrix
  !<
  subroutine gnc_ac(this, sparse)
    ! -- modules
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(GhostNodeType) :: this
    type(sparsematrix), intent(inout) :: sparse
    ! -- local
    integer(I4B) :: ignc, jidx, noden, nodem, nodej
    !
    ! -- Expand the sparse matrix for ghost node connections.  No need to add
    !    connection between n and m as they must be connected some other way
    !    that will calculate the conductance.
    if (this%implicit) then
      do ignc = 1, this%nexg
        noden = this%nodem1(ignc) + this%m1%moffset
        nodem = this%nodem2(ignc) + this%m2%moffset
        jloop: do jidx = 1, this%numjs
          nodej = this%nodesj(jidx, ignc)
          if (nodej == 0) cycle
          nodej = nodej + this%m1%moffset
          call sparse%addconnection(nodem, nodej, 1)
          call sparse%addconnection(nodej, nodem, 1)
          call sparse%addconnection(noden, nodej, 1)
          call sparse%addconnection(nodej, noden, 1)
        end do jloop
      end do
    end if
  end subroutine gnc_ac

  !> @brief Single or Two-Model GNC Map Connections
  !!
  !! Fill the following mapping arrays:
  !!  this%idiagn, this%idiagm (diagonal positions in solution amat)
  !!  this%idxglo (nm connection in solution amat)
  !!  this%idxsymglo (mn connection in solution amat)
  !!  this%jposinrown (position of j in row n of solution amat)
  !!  this%jposinrowm (position of j in row m of solution amat)
  !<
  subroutine gnc_mc(this, matrix_sln)
    ! -- modules
    use SimModule, only: store_error, store_error_unit, count_errors
    use SimVariablesModule, only: errmsg
    ! -- dummy
    class(GhostNodeType) :: this
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: noden, nodem, ipos, ignc, jidx, nodej
    ! -- formats
    character(len=*), parameter :: fmterr = &
      "('GHOST NODE ERROR.  Cell ', i0, ' in model ', a, &
        &' is not connected to cell ', i0, ' in model ', a)"
    !
    ! -- Find the location of Cnm in the global solution and store it in
    !    this%idxglo
    do ignc = 1, this%nexg
      noden = this%nodem1(ignc) + this%m1%moffset
      nodem = this%nodem2(ignc) + this%m2%moffset
      !
      ! -- store diagonal positions in idiagn and idiagm
      this%idiagn(ignc) = matrix_sln%get_position_diag(noden)
      this%idiagm(ignc) = matrix_sln%get_position_diag(nodem)
      !if(this%implicit) then
      !  this%idiagn(ignc) = iasln(noden)
      !  this%idiagm(ignc) = iasln(nodem)
      !endif
      !
      ! -- find location of m in row n of global solution, and v.v.
      this%idxglo(ignc) = matrix_sln%get_position(noden, nodem)
      this%idxsymglo(ignc) = matrix_sln%get_position(nodem, noden)
      !
      ! -- Check to make sure idxglo is set
      if (this%idxglo(ignc) == -1) then
        write (errmsg, fmterr) this%nodem1(ignc), trim(this%m1%name), &
          this%nodem2(ignc), trim(this%m2%name)
        call store_error(errmsg)
      end if
      !
    end do
    !
    ! -- Stop if errors
    if (count_errors() > 0) then
      call store_error_unit(this%inunit)
    end if
    !
    ! -- find locations of j in rows n and row m of global solution
    if (this%implicit) then
      do ignc = 1, this%nexg
        noden = this%nodem1(ignc) + this%m1%moffset
        nodem = this%nodem2(ignc) + this%m2%moffset
        !
        do jidx = 1, this%numjs
          nodej = this%nodesj(jidx, ignc)
          if (nodej > 0) nodej = nodej + this%m1%moffset
          !
          ! -- search for nodej in row n, unless it is 0
          if (nodej == 0) then
            ipos = 0
            this%jposinrown(jidx, ignc) = ipos
          else
            this%jposinrown(jidx, ignc) = matrix_sln%get_position(noden, nodej)
          end if
          !
          ! -- search for nodej in row m
          if (nodej == 0) then
            ipos = 0
            this%jposinrowm(jidx, ignc) = ipos
          else
            this%jposinrowm(jidx, ignc) = matrix_sln%get_position(nodem, nodej)
          end if
        end do
      end do
    end if
  end subroutine gnc_mc

  !> @brief Store the n-m Picard conductance in cond prior to the Newton terms
  !! terms being added
  !<
  subroutine gnc_fmsav(this, kiter, matrix)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy
    class(GhostNodeType) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix
    ! -- local
    integer(I4B) :: ignc, ipos
    real(DP) :: cond
    !
    ! -- An ipos value of zero indicates that noden is not connected to
    !    nodem, and therefore the conductance is zero.
    gncloop: do ignc = 1, this%nexg
      ipos = this%idxglo(ignc)
      if (ipos > 0) then
        cond = matrix%get_value_pos(ipos)
      else
        cond = DZERO
      end if
      this%cond(ignc) = cond
    end do gncloop
  end subroutine gnc_fmsav

  !> @brief Fill matrix terms
  !!
  !! Add the GNC terms to the solution matrix or model rhs depending on whether
  !! whether GNC is implicit or explicit
  !<
  subroutine gnc_fc(this, kiter, matrix)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy
    class(GhostNodeType) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix
    ! -- local
    integer(I4B) :: ignc, j, noden, nodem, ipos, jidx, iposjn, iposjm
    real(DP) :: cond, alpha, aterm, rterm
    !
    ! -- If this is a single model gnc (not an exchange across models), then
    !    pull conductances out of the system matrix and store them in this%cond
    if (this%smgnc) call this%gnc_fmsav(kiter, matrix)
    !
    ! -- Add gnc terms to rhs or to the matrix depending on whether gnc is implicit
    !    or explicit
    gncloop: do ignc = 1, this%nexg
      noden = this%nodem1(ignc)
      nodem = this%nodem2(ignc)
      if (this%m1%ibound(noden) == 0 .or. &
          this%m2%ibound(nodem) == 0) cycle gncloop
      ipos = this%idxglo(ignc)
      cond = this%cond(ignc)
      jloop: do jidx = 1, this%numjs
        j = this%nodesj(jidx, ignc)
        if (j == 0) cycle
        alpha = this%alphasj(jidx, ignc)
        if (alpha == DZERO) cycle
        aterm = alpha * cond
        if (this%implicit) then
          iposjn = this%jposinrown(jidx, ignc)
          iposjm = this%jposinrowm(jidx, ignc)
          call matrix%add_value_pos(this%idiagn(ignc), aterm)
          call matrix%add_value_pos(iposjn, -aterm)
          call matrix%add_value_pos(this%idxsymglo(ignc), -aterm)
          call matrix%add_value_pos(iposjm, aterm)
        else
          rterm = aterm * (this%m1%x(noden) - this%m1%x(j))
          this%m1%rhs(noden) = this%m1%rhs(noden) - rterm
          this%m2%rhs(nodem) = this%m2%rhs(nodem) + rterm
        end if
      end do jloop
    end do gncloop
  end subroutine gnc_fc

  !> @brief Fill GNC Newton terms
  !!
  !! Required arguments:
  !!  kiter : outer iteration number
  !!  matrix_sln: the solution matrix
  !!  condsat is of size(njas) if single model, otherwise nexg
  !!
  !! Optional arguments:
  !!   ihc_opt : an optional vector of size(nexg), which contains a horizontal
  !!     connection code (0=vertical, 1=horizontal, 2=vertically staggered)
  !!   ivarcv_opt : variable vertical conductance flag (default is 0)
  !!   ictm1_opt : icelltype for model 1 integer vector (default is 1)
  !!   ictm2_opt : icelltype for model 2 integer vector (default is 1)
  !<
  subroutine gnc_fn(this, kiter, matrix_sln, condsat, ihc_opt, &
                    ivarcv_opt, ictm1_opt, ictm2_opt)
    ! -- modules
    use ConstantsModule, only: DZERO
    use SmoothingModule, only: sQuadraticSaturationDerivative
    ! -- dummy
    class(GhostNodeType) :: this
    integer(I4B) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    real(DP), dimension(:), intent(in) :: condsat
    integer(I4B), dimension(:), optional :: ihc_opt
    integer(I4B), optional :: ivarcv_opt
    integer(I4B), dimension(:), optional :: ictm1_opt
    integer(I4B), dimension(:), optional :: ictm2_opt
    ! -- local
    integer(I4B) :: ignc, jidx, ipos, isympos, ihc, ivarcv
    integer(I4B) :: nodej, noden, nodem
    integer(I4B) :: iups, ictup
    real(DP) :: csat, alpha, consterm, term, derv
    real(DP) :: xup, topup, botup
    !
    ! -- Set the ivarcv to indicate whether or not the vertical conductance
    !    is a function of water table
    ivarcv = 0
    if (present(ivarcv_opt)) ivarcv = ivarcv_opt
    !
    gncloop: do ignc = 1, this%nexg
      noden = this%nodem1(ignc)
      nodem = this%nodem2(ignc)
      if (this%m1%ibound(noden) == 0 .or. &
          this%m2%ibound(nodem) == 0) cycle gncloop
      !
      ! -- Assign variables depending on whether single model gnc or exchange
      !    gnc
      if (this%smgnc) then
        ipos = this%m1%dis%con%getjaindex(noden, nodem)
        isympos = this%m1%dis%con%jas(ipos)
        ihc = this%m1%dis%con%ihc(isympos)
        csat = condsat(isympos)
      else
        ihc = ihc_opt(ignc)
        csat = condsat(ignc)
      end if
      !
      ! If vertical connection and not variable cv, then cycle
      if (ihc == 0 .and. ivarcv == 0) cycle
      !
      ! determine upstream node (0 is noden, 1 is nodem)
      iups = 0
      if (this%m2%x(nodem) > this%m1%x(noden)) iups = 1
      !
      ! -- Set the upstream top and bot, and then recalculate for a
      !    vertically staggered horizontal connection
      if (iups == 0) then
        topup = this%m1%dis%top(noden)
        botup = this%m1%dis%bot(noden)
        ictup = 1
        if (present(ictm1_opt)) ictup = ictm1_opt(noden)
        xup = this%m1%x(noden)
      else
        topup = this%m2%dis%top(nodem)
        botup = this%m2%dis%bot(nodem)
        ictup = 1
        if (present(ictm2_opt)) ictup = ictm2_opt(nodem)
        xup = this%m2%x(nodem)
      end if
      !
      ! -- No newton terms if upstream cell is confined
      if (ictup == 0) cycle
      !
      ! -- Handle vertically staggered horizontal connection
      if (ihc == 2) then
        topup = min(this%m1%dis%top(noden), this%m2%dis%top(nodem))
        botup = max(this%m1%dis%bot(noden), this%m2%dis%bot(nodem))
      end if
      !
      ! -- Process each contributing node
      jloop: do jidx = 1, this%numjs
        nodej = this%nodesj(jidx, ignc)
        if (nodej == 0) cycle
        if (this%m1%ibound(nodej) == 0) cycle
        alpha = this%alphasj(jidx, ignc)
        if (alpha == DZERO) cycle
        consterm = csat * alpha * (this%m1%x(noden) - this%m1%x(nodej))
        derv = sQuadraticSaturationDerivative(topup, botup, xup)
        term = consterm * derv
        if (iups == 0) then
          call matrix_sln%add_value_pos(this%idiagn(ignc), term)
          if (this%m2%ibound(nodem) > 0) then
            call matrix_sln%add_value_pos(this%idxsymglo(ignc), -term)
          end if
          this%m1%rhs(noden) = this%m1%rhs(noden) + term * this%m1%x(noden)
          this%m2%rhs(nodem) = this%m2%rhs(nodem) - term * this%m1%x(noden)
        else
          call matrix_sln%add_value_pos(this%idiagm(ignc), -term)
          if (this%m1%ibound(noden) > 0) then
            call matrix_sln%add_value_pos(this%idxglo(ignc), term)
          end if
          this%m1%rhs(noden) = this%m1%rhs(noden) + term * this%m2%x(nodem)
          this%m2%rhs(nodem) = this%m2%rhs(nodem) - term * this%m2%x(nodem)
        end if
      end do jloop
    end do gncloop
  end subroutine gnc_fn

  !> @brief Single Model GNC Output
  !!
  !! Output GNC deltaQgnc values
  !<
  subroutine gnc_ot(this, ibudfl)
    ! -- dummy
    class(GhostNodeType) :: this
    integer(I4B), intent(in) :: ibudfl
    ! -- local
    integer(I4B) :: ignc
    real(DP) :: deltaQgnc
    character(len=LINELENGTH) :: nodenstr, nodemstr
    ! -- format
    character(len=*), parameter :: fmtgnc = "(i10, 2a10, 2(1pg15.6))"
    !
    ! -- Process each gnc and output deltaQgnc
    if (ibudfl /= 0 .and. this%iprflow /= 0) then
      write (this%iout, '(//, a)') 'GHOST NODE CORRECTION RESULTS'
      write (this%iout, '(3a10, 2a15)') 'GNC NUM', 'NODEN', 'NODEM', &
        'DELTAQGNC', 'CONDNM'
      do ignc = 1, this%nexg
        deltaQgnc = this%deltaQgnc(ignc)
        call this%m1%dis%noder_to_string(this%nodem1(ignc), nodenstr)
        call this%m2%dis%noder_to_string(this%nodem2(ignc), nodemstr)
        write (this%iout, fmtgnc) ignc, trim(adjustl(nodenstr)), &
          trim(adjustl(nodemstr)), &
          deltaQgnc, this%cond(ignc)
      end do
    end if
  end subroutine gnc_ot

  !> @brief Add GNC to flowja
  !<
  subroutine gnc_cq(this, flowja)
    ! -- dummy
    class(GhostNodeType) :: this
    real(DP), dimension(:), intent(inout) :: flowja
    ! -- local
    integer(I4B) :: ignc, n1, n2, ipos, isympos
    real(DP) :: deltaQgnc
    !
    ! -- go through each gnc and add deltagnc to flowja
    do ignc = 1, this%nexg
      !
      ! -- calculate correction term between n1 and n2 connection
      n1 = this%nodem1(ignc)
      n2 = this%nodem2(ignc)
      deltaQgnc = this%deltaQgnc(ignc)
      !
      ! -- find the positions of this connection in the csr array
      ipos = this%m1%dis%con%getjaindex(n1, n2)
      isympos = this%m1%dis%con%isym(ipos)
      !
      ! -- add/subtract the corrections
      flowja(ipos) = flowja(ipos) + deltaQgnc
      flowja(isympos) = flowja(isympos) - deltaQgnc
      !
    end do
  end subroutine gnc_cq

  !> @brief Single Model deltaQgnc (ghost node correction flux)
  !!
  !! Calculate the deltaQgnc value for any GNC in the GNC list
  !<
  function deltaQgnc(this, ignc)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- Return
    real(DP) :: deltaQgnc
    ! -- dummy
    class(GhostNodeType) :: this
    integer(I4B), intent(in) :: ignc
    ! -- local
    integer(I4B) :: noden, nodem, nodej, jidx
    real(DP) :: sigalj, alpha, hd, aterm, cond
    !
    ! -- initialize values
    deltaQgnc = DZERO
    sigalj = DZERO
    hd = DZERO
    noden = this%nodem1(ignc)
    nodem = this%nodem2(ignc)
    !
    ! -- calculate deltaQgnc
    if (this%m1%ibound(noden) /= 0 .and. this%m2%ibound(nodem) /= 0) then
      jloop: do jidx = 1, this%numjs
        nodej = this%nodesj(jidx, ignc)
        if (nodej == 0) cycle jloop
        if (this%m1%ibound(nodej) == 0) cycle jloop
        alpha = this%alphasj(jidx, ignc)
        sigalj = sigalj + alpha
        hd = hd + alpha * this%m1%x(nodej)
      end do jloop
      aterm = sigalj * this%m1%x(noden) - hd
      cond = this%cond(ignc)
      deltaQgnc = aterm * cond
    end if
  end function deltaQgnc

  !> @brief Allocate gnc scalar variables
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GhostNodeType) :: this
    !
    ! -- allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    call mem_allocate(this%smgnc, 'SMGNC', this%memoryPath)
    call mem_allocate(this%implicit, 'IMPLICIT', this%memoryPath)
    call mem_allocate(this%i2kn, 'I2KN', this%memoryPath)
    call mem_allocate(this%nexg, 'NEXG', this%memoryPath)
    call mem_allocate(this%numjs, 'NUMJS', this%memoryPath)
    !
    ! -- Initialize values
    this%smgnc = .true.
    this%implicit = .true.
    this%i2kn = .false.
    this%nexg = 0
    this%numjs = 0
  end subroutine allocate_scalars

  !> @brief Allocate gnc scalar variables
  !<
  subroutine allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GhostNodeType) :: this
    !
    ! -- allocate memory for arrays
    call mem_allocate(this%nodem1, this%nexg, 'NODEM1', this%memoryPath)
    call mem_allocate(this%nodem2, this%nexg, 'NODEM2', this%memoryPath)
    call mem_allocate(this%nodesj, this%numjs, this%nexg, 'NODESJ', &
                      this%memoryPath)
    call mem_allocate(this%alphasj, this%numjs, this%nexg, 'ALPHASJ', &
                      this%memoryPath)
    call mem_allocate(this%cond, this%nexg, 'COND', this%memoryPath)
    call mem_allocate(this%idxglo, this%nexg, 'IDXGLO', this%memoryPath)
    call mem_allocate(this%idiagn, this%nexg, 'IDIAGN', this%memoryPath)
    call mem_allocate(this%idiagm, this%nexg, 'IDIAGM', this%memoryPath)
    call mem_allocate(this%idxsymglo, this%nexg, 'IDXSYMGLO', this%memoryPath)
    if (this%implicit) then
      call mem_allocate(this%jposinrown, this%numjs, this%nexg, 'JPOSINROWN', &
                        this%memoryPath)
      call mem_allocate(this%jposinrowm, this%numjs, this%nexg, 'JPOSINROWM', &
                        this%memoryPath)
    else
      call mem_allocate(this%jposinrown, 0, 0, 'JPOSINROWN', this%memoryPath)
      call mem_allocate(this%jposinrowm, 0, 0, 'JPOSINROWM', this%memoryPath)
    end if
  end subroutine allocate_arrays

  !> @brief Deallocate memory
  !<
  subroutine gnc_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GhostNodeType) :: this
    !
    call mem_deallocate(this%smgnc)
    call mem_deallocate(this%implicit)
    call mem_deallocate(this%i2kn)
    call mem_deallocate(this%nexg)
    call mem_deallocate(this%numjs)
    !
    ! -- Arrays
    if (this%inunit > 0) then
      call mem_deallocate(this%nodem1)
      call mem_deallocate(this%nodem2)
      call mem_deallocate(this%nodesj)
      call mem_deallocate(this%alphasj)
      call mem_deallocate(this%cond)
      call mem_deallocate(this%idxglo)
      call mem_deallocate(this%idiagn)
      call mem_deallocate(this%idiagm)
      call mem_deallocate(this%idxsymglo)
      call mem_deallocate(this%jposinrown)
      call mem_deallocate(this%jposinrowm)
    end if
    !
    ! -- deallocate NumericalPackageType
    call this%NumericalPackageType%da()
  end subroutine gnc_da

  !> @ brief Source options for gnc package
  !!
  !!  Source options for GNC package.
  !!
  !<
  subroutine source_options(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use GwfGncInputModule, only: GwfGncParamFoundType
    ! -- dummy
    class(GhostNodeType) :: this
    ! -- local
    type(GwfGncParamFoundType) :: found

    ! -- update defaults from input context
    call mem_set_value(this%iprpak, 'PRINT_INPUT', this%input_mempath, &
                       found%print_input)
    call mem_set_value(this%iprflow, 'PRINT_FLOWS', this%input_mempath, &
                       found%print_flows)
    call mem_set_value(this%implicit, 'EXPLICIT', this%input_mempath, &
                       found%explicit)

    ! -- update internal state
    if (found%explicit) this%implicit = .false.

    ! -- log options
    write (this%iout, '(1x,a)') 'PROCESSING GNC OPTIONS'
    if (found%print_input) &
      write (this%iout, '(4x,a)') &
      'THE LIST OF GHOST-NODE CORRECTIONS WILL BE PRINTED.'
    if (found%print_flows) &
      write (this%iout, '(4x,a)') &
      'DELTAQGNC VALUES WILL BE PRINTED TO THE LIST FILE.'
    if (found%explicit) &
      write (this%iout, '(4x,a)') 'GHOST NODE CORRECTION IS EXPLICIT.'
    write (this%iout, '(1x,a)') 'END OF GNC OPTIONS'
  end subroutine source_options

  !> @ brief Source dimensions for gnc package
  !!
  !!  Source dimensions for GNC package.
  !!
  !<
  subroutine source_dimensions(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use GwfGncInputModule, only: GwfGncParamFoundType
    ! -- dummy
    class(GhostNodeType) :: this
    ! -- local
    type(GwfGncParamFoundType) :: found

    ! -- update defaults from input context
    call mem_set_value(this%nexg, 'NUMGNC', this%input_mempath, &
                       found%numgnc)
    call mem_set_value(this%numjs, 'NUMALPHAJ', this%input_mempath, &
                       found%numalphaj)

    ! -- log dimensions
    write (this%iout, '(/1x,a)') 'PROCESSING GNC DIMENSIONS'
    write (this%iout, '(4x,a,i7)') 'NUMGNC = ', this%nexg
    write (this%iout, '(4x,a,i7)') 'NUMAPHAJ = ', this%numjs
    write (this%iout, '(1x,a)') 'END OF GNC DIMENSIONS'
  end subroutine source_dimensions

  !> @ brief Source GNCDATA block
  !<
  subroutine source_data(this)
    ! -- modules
    use SimModule, only: store_error, count_errors, store_error_filename
    use MemoryManagerModule, only: mem_setptr
    use GeomUtilModule, only: get_node
    ! -- dummy variables
    class(GhostNodeType), intent(inout) :: this
    ! -- locals
    integer(I4B), dimension(:, :), pointer, &
      contiguous :: cellidsm, cellidsn, cellidsj
    real(DP), dimension(:, :), pointer, contiguous :: alphasj
    integer(I4B), dimension(:), pointer :: cellidm, cellidn, cellidj
    integer(I4B), dimension(:), allocatable :: nodesuj
    integer(I4B) :: n, m, nodeun, nodeum, nodeuj, nodern, noderm, noderj, istart
    character(len=LINELENGTH) :: fmtgnc

    ! Construct the fmtgnc format
    write (fmtgnc, '("(2i10,",i0,"i10,",i0, "(1pg15.6))")') this%numjs, &
      this%numjs

    write (this%iout, '(/1x,a)') 'PROCESSING GNC DATA'

    ! allocate temporary nodeu array
    allocate (nodesuj(this%numjs))

    ! set input context pointers
    call mem_setptr(cellidsn, 'CELLIDN', this%input_mempath)
    call mem_setptr(cellidsm, 'CELLIDM', this%input_mempath)
    call mem_setptr(cellidsj, 'CELLIDSJ', this%input_mempath)
    call mem_setptr(alphasj, 'ALPHASJ', this%input_mempath)

    ! update state
    do n = 1, this%nexg

      ! set cellid
      cellidn => cellidsn(:, n)
      cellidm => cellidsm(:, n)

      ! set node user
      if (this%m1%dis%ndim == 1) then
        nodeun = cellidn(1)
      elseif (this%m1%dis%ndim == 2) then
        nodeun = get_node(cellidn(1), 1, cellidn(2), &
                          this%m1%dis%mshape(1), 1, &
                          this%m1%dis%mshape(2))
      else
        nodeun = get_node(cellidn(1), cellidn(2), cellidn(3), &
                          this%m1%dis%mshape(1), &
                          this%m1%dis%mshape(2), &
                          this%m1%dis%mshape(3))
      end if

      ! set node user
      if (this%m2%dis%ndim == 1) then
        nodeum = cellidm(1)
      elseif (this%m2%dis%ndim == 2) then
        nodeum = get_node(cellidm(1), 1, cellidm(2), &
                          this%m2%dis%mshape(1), 1, &
                          this%m2%dis%mshape(2))
      else
        nodeum = get_node(cellidm(1), cellidm(2), cellidm(3), &
                          this%m2%dis%mshape(1), &
                          this%m2%dis%mshape(2), &
                          this%m2%dis%mshape(3))
      end if

      ! set nodes
      nodern = this%m1%dis%get_nodenumber(nodeun, 1)
      noderm = this%m2%dis%get_nodenumber(nodeum, 1)
      if (nodern <= 0 .or. &
          noderm <= 0) then
        cycle
      else
        this%nodem1(n) = nodern
        this%nodem2(n) = noderm
      end if

      ! set nodesj and alphajs
      cellidj => cellidsj(:, n)
      do m = 1, this%numjs
        ! set cellid start index
        istart = (m - 1) * this%m1%dis%ndim + 1

        ! set node user
        if (this%m1%dis%ndim == 1) then
          nodeuj = cellidj(istart)
        elseif (this%m1%dis%ndim == 2) then
          nodeuj = get_node(cellidj(istart), 1, cellidj(istart + 1), &
                            this%m1%dis%mshape(1), 1, &
                            this%m1%dis%mshape(2))
        else
          nodeuj = get_node(cellidj(istart), &
                            cellidj(istart + 1), &
                            cellidj(istart + 2), &
                            this%m1%dis%mshape(1), &
                            this%m1%dis%mshape(2), &
                            this%m1%dis%mshape(3))
        end if

        ! set nodesj
        if (nodeuj > 0) then
          nodesuj(m) = nodeuj
          noderj = this%m1%dis%get_nodenumber(nodeuj, 1)
          if (noderj <= 0) then
            cycle
          else
            this%nodesj(m, n) = noderj
          end if
        else
          nodesuj(m) = 0
          this%nodesj(m, n) = 0
        end if

        ! alphaj
        this%alphasj(m, n) = alphasj(m, n)
      end do

      if (this%iprpak /= 0) then
        write (this%iout, fmtgnc) nodeun, nodeum, &
          (nodesuj(m), m=1, this%numjs), &
          (this%alphasj(m, n), m=1, this%numjs)
      end if
    end do

    write (this%iout, '(1x,a)') 'END OF GNC DATA'

    ! check errors
    if (count_errors() > 0) then
      call store_error('Errors encountered in GNC input file.')
      call store_error_filename(this%input_fname)
    end if

    ! cleanup
    deallocate (nodesuj)
  end subroutine source_data

end module GhostNodeModule
