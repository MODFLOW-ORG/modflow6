!> @brief This module contains the SwfGwfExchangeModule Module
!!
!! This module contains the code for connecting a SWF model with
!! a GWF model.  The SwfGwfExchangeType class is a parent to the CHF and
!! OLF models and should not be used directly.
!!
!<
module SwfGwfExchangeModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, DZERO, LENFTYPE
  use SimVariablesModule, only: errmsg, iout, model_loc_idx
  use SimModule, only: count_errors, store_error, store_error_unit, &
                       store_error_filename
  use InputOutputModule, only: getunit, openfile
  use MemoryManagerModule, only: mem_allocate
  use BaseModelModule, only: BaseModelType, GetBaseModelFromList
  use NumericalModelModule, only: NumericalModelType
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use ListsModule, only: basemodellist, baseexchangelist
  use ObsModule, only: obs_cr, ObsType
  use MemoryHelperModule, only: create_mem_path
  use NumericalExchangeModule, only: NumericalExchangeType
  use GwfModule, only: GwfModelType
  use SwfModule, only: SwfModelType
  use MatrixBaseModule
  use TableModule, only: TableType, table_cr

  implicit none
  private
  public :: SwfGwfExchangeType

  type, extends(NumericalExchangeType) :: SwfGwfExchangeType

    class(NumericalModelType), pointer :: model1 => null() !< model 1
    class(NumericalModelType), pointer :: model2 => null() !< model 2
    class(SwfModelType), pointer :: swfmodel => null() !< pointer to SWF Model
    class(GwfModelType), pointer :: gwfmodel => null() !< pointer to GWF Model

    character(len=LENFTYPE), pointer :: swf_ftype => null() !< type of swf model (CHF or OLF)
    character(len=LINELENGTH), pointer :: filename => null() !< name of the input file
    integer(I4B), pointer :: ipr_input => null() !< flag to print input
    integer(I4B), pointer :: ipr_flow => null() !< print flag for cell by cell flows

    integer(I4B), pointer :: ifixedcond => null() !< conductance is fixed as product of bedleak and cfact

    integer(I4B), pointer :: nexg => null() !< number of exchanges
    integer(I4B), dimension(:), pointer, contiguous :: nodeswf => null() !< node numbers in swf model
    integer(I4B), dimension(:), pointer, contiguous :: nodegwf => null() !< node numbers in gwf model
    real(DP), dimension(:), pointer, contiguous :: bedleak => null() !< bed leakance, size: nexg
    real(DP), dimension(:), pointer, contiguous :: cfact => null() !< factor used in conductance calculation, size: nexg
    integer(I4B), dimension(:), pointer, contiguous :: idxglo => null() !< mapping to global (solution) amat
    integer(I4B), dimension(:), pointer, contiguous :: idxsymglo => null() !< mapping to global (solution) symmetric amat
    real(DP), dimension(:), pointer, contiguous :: simvals => null() !< simulated flow rate for each exchange

    integer(I4B), pointer :: inobs => null() !< unit number for GWF-GWF observations
    type(ObsType), pointer :: obs => null() !< observation object

    ! -- table objects
    type(TableType), pointer :: outputtab1 => null()
    type(TableType), pointer :: outputtab2 => null()

  contains

    ! procedure :: exg_df (delegated to CHF and OLF df routines)
    procedure :: exg_ac => swf_gwf_ac
    procedure :: exg_mc => swf_gwf_mc
    procedure :: exg_fc => swf_gwf_fc
    procedure :: exg_cq => swf_gwf_cq
    procedure :: exg_bd => swf_gwf_bd
    procedure :: exg_ot => swf_gwf_ot
    procedure :: exg_da => swf_gwf_da
    procedure :: initialize
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: swf_gwf_calc_simvals
    procedure :: swf_gwf_save_simvals
    procedure :: qcalc
    procedure :: get_cond
    procedure :: get_wetted_perimeter
    procedure :: swf_gwf_add_to_flowja
    procedure, private :: swf_gwf_chd_bd
    !todo: procedure :: swf_gwf_bdsav_model
    procedure :: swf_gwf_bdsav
    procedure :: connects_model => swf_gwf_connects_model

    procedure, pass(this) :: noder
    procedure, pass(this) :: cellstr

  end type SwfGwfExchangeType

contains

  !> @ brief Initialize SWF GWF exchange
  !<
  subroutine initialize(this, filename, name, swf_ftype, id, m1_id, m2_id, &
                        input_mempath)
    ! dummy
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    character(len=*), intent(in) :: filename !< filename for reading
    character(len=*) :: name !< exchange name
    character(len=*) :: swf_ftype !< type of swf model, CHF or OLF
    integer(I4B), intent(in) :: id !< id for the exchange
    integer(I4B), intent(in) :: m1_id !< id for model 1
    integer(I4B), intent(in) :: m2_id !< id for model 2
    character(len=*), intent(in) :: input_mempath
    ! -- local
    class(BaseModelType), pointer :: mb
    integer(I4B) :: m1_index, m2_index

    ! Assign id and name
    this%id = id
    this%name = name
    this%memoryPath = create_mem_path(this%name)
    this%input_mempath = input_mempath

    ! allocate scalars and set defaults
    call this%allocate_scalars()
    this%filename = filename
    this%swf_ftype = swf_ftype
    this%typename = trim(swf_ftype)//'-GWF'

    ! set swfmodel
    m1_index = model_loc_idx(m1_id)
    if (m1_index > 0) then
      mb => GetBaseModelFromList(basemodellist, m1_index)
      select type (mb)
      class is (SwfModelType)
        this%model1 => mb
        this%swfmodel => mb
      end select
    end if
    ! this%v_model1 => get_virtual_model(m1_id)
    ! this%is_datacopy = .not. this%v_model1%is_local

    ! set gwfmodel
    m2_index = model_loc_idx(m2_id)
    if (m2_index > 0) then
      mb => GetBaseModelFromList(basemodellist, m2_index)
      select type (mb)
      type is (GwfModelType)
        this%model2 => mb
        this%gwfmodel => mb
      end select
    end if
    ! this%v_model2 => get_virtual_model(m2_id)

    ! Verify that the surface water model is of the correct type
    if (.not. associated(this%swfmodel) .and. m1_index > 0) then
      write (errmsg, '(7a)') &
        'Problem with ', &
        trim(this%typename), &
        ' exchange ', &
        trim(this%name), &
        '.  Specified ', &
        trim(this%swf_ftype), &
        ' model does not appear to be of the correct type.'
      call store_error(errmsg, terminate=.true.)
    end if

    ! Verify that gwf model is of the correct type
    if (.not. associated(this%gwfmodel) .and. m2_index > 0) then
      write (errmsg, '(3a)') 'Problem with SWF-GWF exchange ', &
        trim(this%name), &
        '.  Specified GWF model does not appear to be of the correct type.'
      call store_error(errmsg, terminate=.true.)
    end if

    ! Create the obs package
    call obs_cr(this%obs, this%inobs)

  end subroutine initialize

  !> @ brief Add connections
  !!
  !! Override parent exg_ac so that gnc can add connections here.
  !<
  subroutine swf_gwf_ac(this, sparse)
    ! -- modules
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    type(sparsematrix), intent(inout) :: sparse
    ! -- local
    integer(I4B) :: n, iglo, jglo
    !
    ! -- add exchange connections
    do n = 1, this%nexg
      iglo = this%nodeswf(n) + this%swfmodel%moffset
      jglo = this%nodegwf(n) + this%gwfmodel%moffset
      call sparse%addconnection(iglo, jglo, 1)
      call sparse%addconnection(jglo, iglo, 1)
    end do
  end subroutine swf_gwf_ac

  !> @ brief Map connections
  !!
  !! Map the connections in the global matrix
  !<
  subroutine swf_gwf_mc(this, matrix_sln)
    ! -- modules
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    class(MatrixBaseType), pointer :: matrix_sln !< the system matrix
    ! -- local
    integer(I4B) :: n, iglo, jglo
    !
    ! -- map exchange connections
    do n = 1, this%nexg
      iglo = this%nodeswf(n) + this%swfmodel%moffset
      jglo = this%nodegwf(n) + this%gwfmodel%moffset
      this%idxglo(n) = matrix_sln%get_position(iglo, jglo)
      this%idxsymglo(n) = matrix_sln%get_position(jglo, iglo)
    end do
  end subroutine swf_gwf_mc

  !> @ brief Fill coefficients
  !!
  !! Fill conductance into coefficient matrix.  For now assume
  !! all connections are vertical and no newton correction is
  !! needed.
  !<
  subroutine swf_gwf_fc(this, kiter, matrix_sln, rhs_sln, inwtflag)
    ! modules
    use MathUtilModule, only: get_perturbation
    ! dummy
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    real(DP), dimension(:), intent(inout) :: rhs_sln
    integer(I4B), optional, intent(in) :: inwtflag
    ! -- local
    integer(I4B) :: iexg
    integer(I4B) :: nodeswf
    integer(I4B) :: nodegwf
    integer(I4B) :: nodeswf_sln
    integer(I4B) :: nodegwf_sln
    integer(I4B) :: ibdn1
    integer(I4B) :: ibdn2
    real(DP) :: hswf
    real(DP) :: hgwf
    real(DP) :: qnm
    real(DP) :: qeps
    real(DP) :: eps
    real(DP) :: derv

    ! Fill terms into solution matrix and rhs vector
    do iexg = 1, this%nexg

      nodeswf = this%nodeswf(iexg)
      nodegwf = this%nodegwf(iexg)
      nodeswf_sln = this%nodeswf(iexg) + this%swfmodel%moffset
      nodegwf_sln = this%nodegwf(iexg) + this%gwfmodel%moffset
      ibdn1 = this%swfmodel%ibound(nodeswf)
      ibdn2 = this%gwfmodel%ibound(nodegwf)
      hswf = this%swfmodel%x(nodeswf)
      hgwf = this%gwfmodel%x(nodegwf)

      ! First add these terms to the row for the surface water model

      ! Fill the qnm term on the right-hand side
      qnm = this%qcalc(iexg, hswf, hgwf)
      rhs_sln(nodeswf_sln) = rhs_sln(nodeswf_sln) - qnm

      ! Derivative calculation and fill of n terms
      eps = get_perturbation(hswf)
      qeps = this%qcalc(iexg, hswf + eps, hgwf)
      derv = (qeps - qnm) / eps
      call matrix_sln%add_diag_value(nodeswf_sln, derv)
      rhs_sln(nodeswf_sln) = rhs_sln(nodeswf_sln) + derv * hswf

      ! Derivative calculation and fill of m terms
      eps = get_perturbation(hgwf)
      qeps = this%qcalc(iexg, hswf, hgwf + eps)
      derv = (qeps - qnm) / eps
      call matrix_sln%add_value_pos(this%idxglo(iexg), derv)
      rhs_sln(nodeswf_sln) = rhs_sln(nodeswf_sln) + derv * hgwf

      ! now add these terms to the row for the groundwater model

      ! Fill the qnm term on the right-hand side
      qnm = -this%qcalc(iexg, hswf, hgwf)
      rhs_sln(nodegwf_sln) = rhs_sln(nodegwf_sln) - qnm

      ! Derivative calculation and fill of n terms
      eps = get_perturbation(hgwf)
      qeps = -this%qcalc(iexg, hswf, hgwf + eps)
      derv = (qeps - qnm) / eps
      call matrix_sln%add_diag_value(nodegwf_sln, derv)
      rhs_sln(nodegwf_sln) = rhs_sln(nodegwf_sln) + derv * hgwf

      ! Derivative calculation and fill of m terms
      eps = get_perturbation(hswf)
      qeps = -this%qcalc(iexg, hswf + eps, hgwf)
      derv = (qeps - qnm) / eps
      call matrix_sln%add_value_pos(this%idxsymglo(iexg), derv)
      rhs_sln(nodegwf_sln) = rhs_sln(nodegwf_sln) + derv * hswf

    end do

  end subroutine swf_gwf_fc

  !> @ brief Calculate flow
  !!
  !! Calculate flow between two cells and store in simvals, also set
  !! information needed for specific discharge calculation
  !<
  subroutine swf_gwf_cq(this, icnvg, isuppress_output, isolnid)
    ! -- dummy
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    integer(I4B), intent(inout) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    integer(I4B), intent(in) :: isolnid
    !
    ! -- calculate flow and store in simvals
    call this%swf_gwf_calc_simvals()
    !
    ! -- set flows to model edges in NPF
    ! todo: do we add these flows for specific discharge calculation?
    !call this%swf_gwf_set_flow_to_npf()
    !
    ! -- add exchange flows to model's flowja diagonal
    call this%swf_gwf_add_to_flowja()
  end subroutine swf_gwf_cq

  !> @ brief Deallocate
  !!
  !! Deallocate memory associated with this object
  !<
  subroutine swf_gwf_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    !
    ! -- objects
    call this%obs%obs_da()
    deallocate (this%obs)
    !
    ! -- arrays
    call mem_deallocate(this%nodeswf)
    call mem_deallocate(this%nodegwf)
    call mem_deallocate(this%bedleak)
    call mem_deallocate(this%cfact)
    call mem_deallocate(this%idxglo)
    call mem_deallocate(this%idxsymglo)
    call mem_deallocate(this%simvals)
    !
    ! -- scalars
    deallocate (this%swf_ftype)
    deallocate (this%filename)
    call mem_deallocate(this%ipr_input)
    call mem_deallocate(this%ipr_flow)
    call mem_deallocate(this%ifixedcond)
    call mem_deallocate(this%nexg)
    call mem_deallocate(this%inobs)
  end subroutine swf_gwf_da

  !> @ brief Allocate scalars
  !!
  !! Allocate scalar variables
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    ! -- dummy
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    !
    allocate (this%swf_ftype)
    this%swf_ftype = ''
    allocate (this%filename)
    this%filename = ''
    !
    call mem_allocate(this%ipr_input, 'IPR_INPUT', this%memoryPath)
    call mem_allocate(this%ipr_flow, 'IPR_FLOW', this%memoryPath)
    call mem_allocate(this%ifixedcond, 'IFIXEDCOND', this%memoryPath)
    call mem_allocate(this%nexg, 'NEXG', this%memoryPath)
    call mem_allocate(this%inobs, 'INOBS', this%memoryPath)
    !
    this%ipr_input = 0
    this%ipr_flow = 0
    this%ifixedcond = 0
    this%nexg = 0
    this%inobs = 0
  end subroutine allocate_scalars

  !> @brief Allocate array data, using the number of
  !! connected nodes @param nexg
  !<
  subroutine allocate_arrays(this)
    ! -- dummy
    class(SwfGwfExchangeType) :: this !< instance of exchange object
    !
    call mem_allocate(this%nodeswf, this%nexg, 'NODEM1', this%memoryPath)
    call mem_allocate(this%nodegwf, this%nexg, 'NODEM2', this%memoryPath)
    call mem_allocate(this%bedleak, this%nexg, 'BEDLEAK', this%memoryPath)
    call mem_allocate(this%cfact, this%nexg, 'CFACT', this%memoryPath)
    call mem_allocate(this%idxglo, this%nexg, 'IDXGLO', this%memoryPath)
    call mem_allocate(this%idxsymglo, this%nexg, 'IDXSYMGLO', this%memoryPath)
    call mem_allocate(this%simvals, this%nexg, 'SIMVALS', this%memoryPath)
  end subroutine allocate_arrays

  !> @brief
  !<
  function noder(this, model, cellid, iout)
    ! -- modules
    use GeomUtilModule, only: get_node
    ! -- dummy
    class(SwfGwfExchangeType) :: this !< instance of exchange object
    class(NumericalModelType), pointer, intent(in) :: model
    integer(I4B), dimension(:), intent(in) :: cellid
    integer(I4B), intent(in) :: iout !< the output file unit
    integer(I4B) :: noder, node
    !
    if (model%dis%ndim == 1) then
      node = cellid(1)
    elseif (model%dis%ndim == 2) then
      node = get_node(cellid(1), 1, cellid(2), &
                      model%dis%mshape(1), 1, &
                      model%dis%mshape(2))
    else
      node = get_node(cellid(1), cellid(2), cellid(3), &
                      model%dis%mshape(1), &
                      model%dis%mshape(2), &
                      model%dis%mshape(3))
    end if
    noder = model%dis%get_nodenumber(node, 0)
  end function noder

  !> @brief
  !<
  function cellstr(this, model, cellid, iout)
    ! -- modules
    ! -- dummy
    class(SwfGwfExchangeType) :: this !< instance of exchange object
    class(NumericalModelType), pointer, intent(in) :: model
    integer(I4B), dimension(:), intent(in) :: cellid
    integer(I4B), intent(in) :: iout !< the output file unit
    character(len=20) :: cellstr
    character(len=*), parameter :: fmtndim1 = &
                                   "('(',i0,')')"
    character(len=*), parameter :: fmtndim2 = &
                                   "('(',i0,',',i0,')')"
    character(len=*), parameter :: fmtndim3 = &
                                   "('(',i0,',',i0,',',i0,')')"
    !
    cellstr = ''
    !
    select case (model%dis%ndim)
    case (1)
      write (cellstr, fmtndim1) cellid(1)
    case (2)
      write (cellstr, fmtndim2) cellid(1), cellid(2)
    case (3)
      write (cellstr, fmtndim3) cellid(1), cellid(2), cellid(3)
    case default
    end select
  end function cellstr

  !> @brief Calculate flow rates for the exchanges and store them in a member
  !! array
  !<
  subroutine swf_gwf_calc_simvals(this)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    ! -- local
    integer(I4B) :: iexg
    integer(I4B) :: nodeswf, nodegwf
    integer(I4B) :: ibdn1, ibdn2
    real(DP) :: hswf
    real(DP) :: hgwf
    real(DP) :: rrate
    !
    do iexg = 1, this%nexg
      rrate = DZERO
      nodeswf = this%nodeswf(iexg)
      nodegwf = this%nodegwf(iexg)
      ibdn1 = this%swfmodel%ibound(nodeswf)
      ibdn2 = this%gwfmodel%ibound(nodegwf)
      hswf = this%swfmodel%x(nodeswf)
      hgwf = this%gwfmodel%x(nodegwf)
      if (ibdn1 /= 0 .and. ibdn2 /= 0) then
        rrate = this%qcalc(iexg, hswf, hgwf)
      end if
      this%simvals(iexg) = rrate
    end do
  end subroutine swf_gwf_calc_simvals

  !> @ brief Calculate flow
  !!
  !! Calculate the flow for the specified exchange and node numbers.
  !! Flow is positive into the surface water model
  !<
  function qcalc(this, iexg, hswf, hgwf)
    ! -- return
    real(DP) :: qcalc
    ! -- dummy
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    integer(I4B), intent(in) :: iexg
    real(DP), intent(in) :: hswf
    real(DP), intent(in) :: hgwf
    ! -- local
    real(DP) :: cond

    ! Calculate flow between swf and gwf models; positive into swf
    cond = this%get_cond(iexg, hswf, hgwf)
    qcalc = cond * (hgwf - hswf)

  end function qcalc

  !> @ brief Calculate conductance
  !!
  !! Calculate the conductance between the surface water cell
  !! and the underlying groundwater cell.
  !<
  function get_cond(this, iexg, hswf, hgwf)
    ! module
    use SmoothingModule, only: sQuadratic
    ! return
    real(DP) :: get_cond
    ! dummy
    class(SwfGwfExchangeType) :: this !< SwfGwfExchangeType
    integer(I4B), intent(in) :: iexg !< exchange number
    real(DP), intent(in) :: hswf !< surface water model head
    real(DP), intent(in) :: hgwf !< groundwater model head
    ! local
    integer(I4B) :: nodeswf
    real(DP) :: range = 1.d-6
    real(DP) :: depth_ups
    real(DP) :: dydx
    real(DP) :: smooth_factor
    real(DP) :: area
    real(DP) :: perimeter

    ! -- Calculate or return conductance
    area = this%cfact(iexg)
    if (this%ifixedcond == 1) then
      get_cond = this%bedleak(iexg) * area
      return
    end if

    ! Calculate smooth factor between zero, when the upstream-weighted
    ! depth is zero, and 1.0, when the upstream weighted depth is
    ! greater than or equal to the smoothening depth
    nodeswf = this%nodeswf(iexg)
    depth_ups = max(hswf, hgwf) - this%swfmodel%dis%bot(nodeswf)
    call sQuadratic(depth_ups, range, dydx, smooth_factor)

    ! For channel model calculate the interaction area as product
    ! of cfact and upstream-wetted perimeter
    if (this%swfmodel%dfw%is2d == 0) then
      perimeter = this%get_wetted_perimeter(nodeswf, depth_ups)
      area = area * perimeter
    end if

    ! Calculate conductance
    get_cond = smooth_factor * this%bedleak(iexg) * area
  end function get_cond

  !> @ brief Get wetted perimeter for swf channel model
  !<
  function get_wetted_perimeter(this, nodeswf, depth) result(wp)
    ! return
    real(DP) :: wp
    ! dummy
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    integer(I4B), intent(in) :: nodeswf !< node number for surface water model cell
    real(DP), intent(in) :: depth !< water depth in surface water model cell
    ! local
    integer(I4B) :: idcxs
    real(DP) :: width
    real(DP) :: dummy

    idcxs = this%swfmodel%dfw%idcxs(nodeswf)
    call this%swfmodel%dis%get_flow_width(nodeswf, nodeswf, 0, width, dummy)
    wp = this%swfmodel%cxs%get_wetted_perimeter(idcxs, width, depth)

  end function get_wetted_perimeter

  !> @brief Add exchange flow to each model flowja diagonal position so that
  !! residual is calculated correctly.
  !<
  subroutine swf_gwf_add_to_flowja(this)
    ! -- modules
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: n
    integer(I4B) :: idiag
    real(DP) :: flow
    !
    do i = 1, this%nexg
      !
      if (associated(this%swfmodel)) then
        n = this%nodeswf(i)
        if (this%swfmodel%ibound(n) > 0) then
          flow = this%simvals(i)
          idiag = this%swfmodel%ia(n)
          this%swfmodel%flowja(idiag) = this%swfmodel%flowja(idiag) + flow
        end if
      end if
      !
      if (associated(this%gwfmodel)) then
        n = this%nodegwf(i)
        if (this%gwfmodel%ibound(n) > 0) then
          flow = -this%simvals(i)
          idiag = this%gwfmodel%ia(n)
          this%gwfmodel%flowja(idiag) = this%gwfmodel%flowja(idiag) + flow
        end if
      end if
      !
    end do
  end subroutine swf_gwf_add_to_flowja

  !> @ brief Budget
  !!
  !! Accumulate budget terms
  !<
  subroutine swf_gwf_bd(this, icnvg, isuppress_output, isolnid)
    ! -- modules
    use ConstantsModule, only: DZERO, LENBUDTXT, LENPACKAGENAME
    use BudgetModule, only: rate_accumulator
    ! -- dummy
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    integer(I4B), intent(inout) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    integer(I4B), intent(in) :: isolnid
    ! -- local
    character(len=LENBUDTXT), dimension(1) :: budtxt
    real(DP), dimension(2, 1) :: budterm
    real(DP) :: ratin, ratout
    !
    ! -- initialize
    budtxt(1) = '    FLOW-JA-FACE'
    !
    ! -- Calculate ratin/ratout and pass to model budgets
    call rate_accumulator(this%simvals, ratin, ratout)
    !
    ! -- Add the budget terms to model 1
    if (associated(this%swfmodel)) then
      budterm(1, 1) = ratin
      budterm(2, 1) = ratout
      call this%swfmodel%model_bdentry(budterm, budtxt, this%name)
    end if
    !
    ! -- Add the budget terms to model 2
    if (associated(this%gwfmodel)) then
      budterm(1, 1) = ratout
      budterm(2, 1) = ratin
      call this%gwfmodel%model_bdentry(budterm, budtxt, this%name)
    end if
    !
    ! -- Add any flows from one model into a constant head in another model
    !    as a separate budget term called FLOW-JA-FACE-CHD
    call this%swf_gwf_chd_bd()
  end subroutine swf_gwf_bd

  !> @ brief swf-gwf-chd-bd
  !!
  !! Account for flow from an external model into a chd cell
  !<
  subroutine swf_gwf_chd_bd(this)
    ! -- modules
    use ConstantsModule, only: DZERO, LENBUDTXT, LENPACKAGENAME
    ! -- dummy
    class(SwfGwfExchangeType) :: this !<  GwfExchangeType
    ! -- local
    character(len=LENBUDTXT), dimension(1) :: budtxt
    integer(I4B) :: n
    integer(I4B) :: i
    real(DP), dimension(2, 1) :: budterm
    real(DP) :: ratin, ratout
    real(DP) :: q
    !
    ! -- initialize
    budtxt(1) = 'FLOW-JA-FACE-CHD'
    !
    ! -- Add the constant-head budget terms for flow from model 2 into model 1
    if (associated(this%swfmodel)) then
      ratin = DZERO
      ratout = DZERO
      do i = 1, this%nexg
        n = this%nodeswf(i)
        if (this%swfmodel%ibound(n) < 0) then
          q = this%simvals(i)
          if (q > DZERO) then
            ratout = ratout + q
          else
            ratin = ratin - q
          end if
        end if
      end do
      budterm(1, 1) = ratin
      budterm(2, 1) = ratout
      call this%swfmodel%model_bdentry(budterm, budtxt, this%name)
    end if
    !
    ! -- Add the constant-head budget terms for flow from model 1 into model 2
    if (associated(this%gwfmodel)) then
      ratin = DZERO
      ratout = DZERO
      do i = 1, this%nexg
        n = this%nodegwf(i)
        if (this%gwfmodel%ibound(n) < 0) then
          ! -- flip flow sign as flow is relative to model 1
          q = -this%simvals(i)
          if (q > DZERO) then
            ratout = ratout + q
          else
            ratin = ratin - q
          end if
        end if
      end do
      budterm(1, 1) = ratin
      budterm(2, 1) = ratout
      call this%gwfmodel%model_bdentry(budterm, budtxt, this%name)
    end if
  end subroutine swf_gwf_chd_bd

  !> @ brief Budget save
  !!
  !! Output individual flows to listing file and binary budget files
  !<
  subroutine swf_gwf_bdsav(this)
    ! -- dummy
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    ! -- local
    integer(I4B) :: icbcfl, ibudfl
    ! !
    ! ! -- budget for model1
    ! if (associated(this%swfmodel1)) then
    !   call this%swf_gwf_bdsav_model(this%swfmodel1, this%gwfmodel2%name)
    ! end if
    ! !
    ! ! -- budget for model2
    ! if (associated(this%gwfmodel2)) then
    !   call this%swf_gwf_bdsav_model(this%gwfmodel2, this%swfmodel1%name)
    ! end if
    !
    ! -- Set icbcfl, ibudfl to zero so that flows will be printed and
    !    saved, if the options were set in the MVR package
    icbcfl = 1
    ibudfl = 1
    !
    ! -- Calculate and write simulated values for observations
    if (this%inobs /= 0) then
      call this%swf_gwf_save_simvals()
    end if
  end subroutine swf_gwf_bdsav

  ! subroutine swf_gwf_bdsav_model(this, model, neighbor_name)
  !   ! -- modules
  !   use ConstantsModule, only: DZERO, LENBUDTXT, LENPACKAGENAME
  !   use TdisModule, only: kstp, kper
  !   ! -- dummy
  !   class(SwfGwfExchangeType) :: this !< this exchange
  !   class(NumericalModelType), pointer :: model !< the model to save budget for
  !   character(len=*), intent(in) :: neighbor_name !< name of the connected neighbor model
  !   ! -- local
  !   character(len=LENPACKAGENAME + 4) :: packname
  !   character(len=LENBUDTXT), dimension(1) :: budtxt
  !   type(TableType), pointer :: output_tab
  !   character(len=20) :: nodestr
  !   character(len=LENBOUNDNAME) :: bname
  !   integer(I4B) :: ntabrows
  !   integer(I4B) :: nodeu
  !   integer(I4B) :: i, n1, n2, n1u, n2u
  !   integer(I4B) :: ibinun
  !   real(DP) :: ratin, ratout, rrate
  !   logical(LGP) :: is_for_model1
  !   real(DP), dimension(this%naux) :: auxrow
  !   !
  !   budtxt(1) = '    FLOW-JA-FACE'
  !   packname = 'EXG '//this%name
  !   packname = adjustr(packname)
  !   if (associated(model, this%swfmodel1)) then
  !     output_tab => this%outputtab1
  !     is_for_model1 = .true.
  !   else
  !     output_tab => this%outputtab2
  !     is_for_model1 = .false.
  !   end if
  !   !
  !   ! -- update output tables
  !   if (this%ipr_flow /= 0) then
  !     !
  !     ! -- update titles
  !     if (model%oc%oc_save('BUDGET')) then
  !       call output_tab%set_title(packname)
  !     end if
  !     !
  !     ! -- set table kstp and kper
  !     call output_tab%set_kstpkper(kstp, kper)
  !     !
  !     ! -- update maxbound of tables
  !     ntabrows = 0
  !     do i = 1, this%nexg
  !       n1 = this%nodeswf(i)
  !       n2 = this%nodegwf(i)
  !       !
  !       ! -- If both cells are active then calculate flow rate
  !       if (this%swfmodel1%ibound(n1) /= 0 .and. &
  !           this%gwfmodel2%ibound(n2) /= 0) then
  !         ntabrows = ntabrows + 1
  !       end if
  !     end do
  !     if (ntabrows > 0) then
  !       call output_tab%set_maxbound(ntabrows)
  !     end if
  !   end if
  !   !
  !   ! -- Print and write budget terms
  !   !
  !   ! -- Set binary unit numbers for saving flows
  !   if (this%ipakcb /= 0) then
  !     ibinun = model%oc%oc_save_unit('BUDGET')
  !   else
  !     ibinun = 0
  !   end if
  !   !
  !   ! -- If save budget flag is zero for this stress period, then
  !   !    shut off saving
  !   if (.not. model%oc%oc_save('BUDGET')) ibinun = 0
  !   !
  !   ! -- If cell-by-cell flows will be saved as a list, write header.
  !   if (ibinun /= 0) then
  !     call model%dis%record_srcdst_list_header(budtxt(1), &
  !                                              model%name, &
  !                                              this%name, &
  !                                              neighbor_name, &
  !                                              this%name, &
  !                                              this%naux, this%auxname, &
  !                                              ibinun, this%nexg, &
  !                                              model%iout)
  !   end if
  !   !
  !   ! Initialize accumulators
  !   ratin = DZERO
  !   ratout = DZERO
  !   !
  !   ! -- Loop through all exchanges
  !   do i = 1, this%nexg
  !     !
  !     ! -- Assign boundary name
  !     if (this%inamedbound > 0) then
  !       bname = this%boundname(i)
  !     else
  !       bname = ''
  !     end if
  !     !
  !     ! -- Calculate the flow rate between n1 and n2
  !     rrate = DZERO
  !     n1 = this%nodeswf(i)
  !     n2 = this%nodegwf(i)
  !     !
  !     ! -- If both cells are active then calculate flow rate
  !     if (this%v_model1%ibound%get(n1) /= 0 .and. &
  !         this%v_model2%ibound%get(n2) /= 0) then
  !       rrate = this%simvals(i)
  !       !
  !       ! -- Print the individual rates to model list files if requested
  !       if (this%ipr_flow /= 0) then
  !         if (model%oc%oc_save('BUDGET')) then
  !           !
  !           ! -- set nodestr and write outputtab table
  !           if (is_for_model1) then
  !             nodeu = model%dis%get_nodeuser(n1)
  !             call model%dis%nodeu_to_string(nodeu, nodestr)
  !             call output_tab%print_list_entry(i, trim(adjustl(nodestr)), &
  !                                              rrate, bname)
  !           else
  !             nodeu = model%dis%get_nodeuser(n2)
  !             call model%dis%nodeu_to_string(nodeu, nodestr)
  !             call output_tab%print_list_entry(i, trim(adjustl(nodestr)), &
  !                                              -rrate, bname)
  !           end if
  !         end if
  !       end if
  !       if (rrate < DZERO) then
  !         ratout = ratout - rrate
  !       else
  !         ratin = ratin + rrate
  !       end if
  !     end if
  !     !
  !     ! -- If saving cell-by-cell flows in list, write flow
  !     n1u = this%v_model1%dis_get_nodeuser(n1)
  !     n2u = this%v_model2%dis_get_nodeuser(n2)
  !     if (ibinun /= 0) then
  !       if (is_for_model1) then
  ! if (size(auxrow) > 0) then
  !   auxrow(:) = this%auxvar(:, i)
  ! end if
!         call model%dis%record_mf6_list_entry(ibinun, n1u, n2u, rrate, &
  !                                              this%naux, auxrow, &
  !                                              .false., .false.)
  !       else
  !         call model%dis%record_mf6_list_entry(ibinun, n2u, n1u, -rrate, &
  !                                              this%naux, auxrow, &
  !                                              .false., .false.)
  !       end if
  !     end if
  !     !
  !   end do
  !   !
  !   ! -- Return
  !   return
  ! end subroutine swf_gwf_bdsav_model

  !> @ brief Output
  !!
  !! Write output
  !<
  subroutine swf_gwf_ot(this)
    ! -- modules
    use SimVariablesModule, only: iout
    use ConstantsModule, only: DZERO, LINELENGTH
    ! -- dummy
    class(SwfGwfExchangeType) :: this !<  SwfGwfExchangeType
    ! -- local
    integer(I4B) :: iexg, n1, n2
    real(DP) :: flow
    character(len=LINELENGTH) :: node1str, node2str
    ! -- format
    character(len=*), parameter :: fmtheader2 = &
     "(/1x, 'SUMMARY OF EXCHANGE RATES FOR EXCHANGE ', a, ' WITH ID ', i0, /, &
       &2a16, 4a16, /, 96('-'))"
    character(len=*), parameter :: fmtdata = &
                                   "(2a16, 5(1pg16.6))"
    !
    ! -- Call bdsave
    call this%swf_gwf_bdsav()
    !
    ! -- Write a table of exchanges
    if (this%ipr_flow /= 0) then
      write (iout, fmtheader2) trim(adjustl(this%name)), this%id, 'NODEM1', &
        'NODEM2', 'COND', 'X_M1', 'X_M2', 'FLOW'
      do iexg = 1, this%nexg
        n1 = this%nodeswf(iexg)
        n2 = this%nodegwf(iexg)
        flow = this%simvals(iexg)
        call this%swfmodel%dis%noder_to_string(n1, node1str)
        call this%gwfmodel%dis%noder_to_string(n2, node2str)
        write (iout, fmtdata) trim(adjustl(node1str)), &
          trim(adjustl(node2str)), &
          this%bedleak(iexg), this%swfmodel%x(n1), &
          this%gwfmodel%x(n2), flow
      end do
    end if
    !
    ! -- OBS output
    call this%obs%obs_ot()
  end subroutine swf_gwf_ot

  !> @ brief Save simulated flow observations
  !!
  !! Save the simulated flows for each exchange
  !<
  subroutine swf_gwf_save_simvals(this)
    ! -- modules
    use SimModule, only: store_error, store_error_unit
    use SimVariablesModule, only: errmsg
    use ConstantsModule, only: DZERO
    use ObserveModule, only: ObserveType
    ! -- dummy
    class(SwfGwfExchangeType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: n1
    integer(I4B) :: n2
    integer(I4B) :: iexg
    real(DP) :: v
    type(ObserveType), pointer :: obsrv => null()
    !
    ! -- Write simulated values for all gwf-gwf observations
    if (this%obs%npakobs > 0) then
      call this%obs%obs_bd_clear()
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv
        do j = 1, obsrv%indxbnds_count
          iexg = obsrv%indxbnds(j)
          v = DZERO
          select case (obsrv%ObsTypeId)
          case ('FLOW-JA-FACE')
            n1 = this%nodeswf(iexg)
            n2 = this%nodegwf(iexg)
            v = this%simvals(iexg)
          case default
            errmsg = 'Unrecognized observation type: '// &
                     trim(obsrv%ObsTypeId)
            call store_error(errmsg)
            call store_error_unit(this%inobs)
          end select
          call this%obs%SaveOneSimval(obsrv, v)
        end do
      end do
    end if
  end subroutine swf_gwf_save_simvals

  !> @brief Should return true when the exchange should be added to the
  !! solution where the model resides
  !<
  function swf_gwf_connects_model(this, model) result(is_connected)
    ! -- dummy
    class(SwfGwfExchangeType) :: this !< the instance of the exchange
    class(BaseModelType), pointer, intent(in) :: model !< the model to which the exchange might hold a connection
    ! -- return
    logical(LGP) :: is_connected !< true, when connected
    !
    is_connected = .false.
    select type (model)
    class is (GwfModelType)
      if (associated(this%gwfmodel, model)) then
        is_connected = .true.
      end if
    class is (SwfModelType)
      if (associated(this%swfmodel, model)) then
        is_connected = .true.
      end if
    end select
  end function

end module SwfGwfExchangeModule
