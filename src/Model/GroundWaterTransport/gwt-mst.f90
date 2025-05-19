!> -- @ brief Mobile Storage and Transfer (MST) Module
!!
!!    The GwtMstModule contains the GwtMstType, which is the
!!    derived type responsible for adding the effects of
!!      1. Changes in dissolved solute mass
!!      2. Decay of dissolved solute mass
!!      3. Sorption
!!      4. Decay of sorbed solute mass
!<
module GwtMstModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, IZERO, DTWO, DHALF, LENBUDTXT, &
                             MAXCHARLEN, MNORMAL, LINELENGTH, DHNOFLO
  use SimVariablesModule, only: errmsg, warnmsg
  use SimModule, only: store_error, count_errors, &
                       store_warning
  use MatrixBaseModule
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType

  implicit none
  public :: GwtMstType
  public :: mst_cr
  !
  integer(I4B), parameter :: NBDITEMS = 4
  character(len=LENBUDTXT), dimension(NBDITEMS) :: budtxt
  data budtxt/' STORAGE-AQUEOUS', '   DECAY-AQUEOUS', &
    '  STORAGE-SORBED', '    DECAY-SORBED'/

  !> @brief Enumerator that defines the decay options
  !<
  ENUM, BIND(C)
    ENUMERATOR :: DECAY_OFF = 0 !< Decay (or production) of mass inactive (default)
    ENUMERATOR :: DECAY_FIRST_ORDER = 1 !< First-order decay
    ENUMERATOR :: DECAY_ZERO_ORDER = 2 !< Zeroth-order decay
    ENUMERATOR :: SORPTION_OFF = 0 !< Sorption is inactive (default)
    ENUMERATOR :: SORPTION_LINEAR = 1 !< Linear sorption between aqueous and solid phases
    ENUMERATOR :: SORPTION_FREUND = 2 !< Freundlich sorption between aqueous and solid phases
    ENUMERATOR :: SORPTION_LANG = 3 !< Langmuir sorption between aqueous and solid phases
  END ENUM

  !> @ brief Mobile storage and transfer
  !!
  !!  Data and methods for handling changes in solute storage,
  !!  decay of dissolved solute mass, sorption, and decay of
  !!  sorbed mass.
  !<
  type, extends(NumericalPackageType) :: GwtMstType
    !
    ! -- storage
    real(DP), dimension(:), pointer, contiguous :: porosity => null() !< mobile porosity defined as volume mobile voids per volume of mobile domain
    real(DP), dimension(:), pointer, contiguous :: thetam => null() !< mobile porosity defined as volume mobile voids per volume of aquifer
    real(DP), dimension(:), pointer, contiguous :: volfracim => null() !< sum of all immobile domain volume fractions
    real(DP), dimension(:), pointer, contiguous :: ratesto => null() !< rate of mobile storage
    !
    ! -- decay
    integer(I4B), pointer :: idcy => null() !< order of decay rate (0:none, 1:first, 2:zero)
    real(DP), dimension(:), pointer, contiguous :: decay => null() !< first or zero order decay rate (aqueous)
    real(DP), dimension(:), pointer, contiguous :: decay_sorbed => null() !< first or zero order decay rate (sorbed)
    real(DP), dimension(:), pointer, contiguous :: ratedcy => null() !< rate of decay
    real(DP), dimension(:), pointer, contiguous :: decaylast => null() !< decay rate used for last iteration (needed for zero order decay)
    real(DP), dimension(:), pointer, contiguous :: decayslast => null() !< sorbed decay rate used for last iteration (needed for zero order decay)
    !
    ! -- sorption
    integer(I4B), pointer :: isrb => null() !< sorption active flag (0:off, 1:linear, 2:freundlich, 3:langmuir)
    integer(I4B), pointer :: ioutsorbate => null() !< unit number for sorbate concentration output
    real(DP), dimension(:), pointer, contiguous :: bulk_density => null() !< bulk density of mobile domain; mass of mobile domain solid per aquifer volume
    real(DP), dimension(:), pointer, contiguous :: distcoef => null() !< kd distribution coefficient
    real(DP), dimension(:), pointer, contiguous :: sp2 => null() !< second sorption parameter
    real(DP), dimension(:), pointer, contiguous :: ratesrb => null() !< rate of sorption
    real(DP), dimension(:), pointer, contiguous :: ratedcys => null() !< rate of sorbed mass decay
    real(DP), dimension(:), pointer, contiguous :: csrb => null() !< sorbate concentration
    !
    ! -- misc
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    type(TspFmiType), pointer :: fmi => null() !< pointer to fmi object

  contains

    procedure :: mst_ar
    procedure :: mst_fc
    procedure :: mst_fc_sto
    procedure :: mst_fc_dcy
    procedure :: mst_fc_srb
    procedure :: mst_fc_dcy_srb
    procedure :: mst_cq
    procedure :: mst_cq_sto
    procedure :: mst_cq_dcy
    procedure :: mst_cq_srb
    procedure :: mst_cq_dcy_srb
    procedure :: mst_calc_csrb
    procedure :: mst_bd
    procedure :: mst_ot_flow
    procedure :: mst_ot_dv
    procedure :: mst_da
    procedure :: allocate_scalars
    procedure :: addto_volfracim
    procedure :: get_volfracm
    procedure, private :: allocate_arrays
    procedure, private :: read_options
    procedure, private :: read_data

  end type GwtMstType

contains

  !> @ brief Create a new package object
  !!
  !!  Create a new MST object
  !<
  subroutine mst_cr(mstobj, name_model, inunit, iout, fmi)
    ! -- dummy
    type(GwtMstType), pointer :: mstobj !< unallocated new mst object to create
    character(len=*), intent(in) :: name_model !< name of the model
    integer(I4B), intent(in) :: inunit !< unit number of WEL package input file
    integer(I4B), intent(in) :: iout !< unit number of model listing file
    type(TspFmiType), intent(in), target :: fmi !< fmi package for this GWT model
    !
    ! -- Create the object
    allocate (mstobj)
    !
    ! -- create name and memory path
    call mstobj%set_names(1, name_model, 'MST', 'MST')
    !
    ! -- Allocate scalars
    call mstobj%allocate_scalars()
    !
    ! -- Set variables
    mstobj%inunit = inunit
    mstobj%iout = iout
    mstobj%fmi => fmi
    !
    ! -- Initialize block parser
    call mstobj%parser%Initialize(mstobj%inunit, mstobj%iout)
  end subroutine mst_cr

  !> @ brief Allocate and read method for package
  !!
  !!  Method to allocate and read static data for the package.
  !<
  subroutine mst_ar(this, dis, ibound)
    ! -- modules
    ! -- dummy
    class(GwtMstType), intent(inout) :: this !< GwtMstType object
    class(DisBaseType), pointer, intent(in) :: dis !< pointer to dis package
    integer(I4B), dimension(:), pointer, contiguous :: ibound !< pointer to GWT ibound array
    ! -- local
    ! -- formats
    character(len=*), parameter :: fmtmst = &
      "(1x,/1x,'MST -- MOBILE STORAGE AND TRANSFER PACKAGE, VERSION 1, &
      &7/29/2020 INPUT READ FROM UNIT ', i0, //)"
    !
    ! --print a message identifying the mobile storage and transfer package.
    write (this%iout, fmtmst) this%inunit
    !
    ! -- Read options
    call this%read_options()
    !
    ! -- store pointers to arguments that were passed in
    this%dis => dis
    this%ibound => ibound
    !
    ! -- Allocate arrays
    call this%allocate_arrays(dis%nodes)
    !
    ! -- read the data block
    call this%read_data()
  end subroutine mst_ar

  !> @ brief Fill coefficient method for package
  !!
  !!  Method to calculate and fill coefficients for the package.
  !<
  subroutine mst_fc(this, nodes, cold, nja, matrix_sln, idxglo, cnew, &
                    rhs, kiter)
    ! -- modules
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer, intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cold !< concentration at end of last time step
    integer(I4B), intent(in) :: nja !< number of GWT connections
    class(MatrixBaseType), pointer :: matrix_sln !< solution matrix
    integer(I4B), intent(in), dimension(nja) :: idxglo !< mapping vector for model (local) to solution (global)
    real(DP), intent(inout), dimension(nodes) :: rhs !< right-hand side vector for model
    real(DP), intent(in), dimension(nodes) :: cnew !< concentration at end of this time step
    integer(I4B), intent(in) :: kiter !< solution outer iteration number
    !
    ! -- storage contribution
    call this%mst_fc_sto(nodes, cold, nja, matrix_sln, idxglo, rhs)
    !
    ! -- decay contribution
    if (this%idcy /= DECAY_OFF) then
      call this%mst_fc_dcy(nodes, cold, cnew, nja, matrix_sln, idxglo, &
                           rhs, kiter)
    end if
    !
    ! -- sorption contribution
    if (this%isrb /= SORPTION_OFF) then
      call this%mst_fc_srb(nodes, cold, nja, matrix_sln, idxglo, rhs, cnew)
    end if
    !
    ! -- decay sorbed contribution
    if (this%isrb /= SORPTION_OFF .and. this%idcy /= DECAY_OFF) then
      call this%mst_fc_dcy_srb(nodes, cold, nja, matrix_sln, idxglo, rhs, &
                               cnew, kiter)
    end if
  end subroutine mst_fc

  !> @ brief Fill storage coefficient method for package
  !!
  !!  Method to calculate and fill storage coefficients for the package.
  !<
  subroutine mst_fc_sto(this, nodes, cold, nja, matrix_sln, idxglo, rhs)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer, intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cold !< concentration at end of last time step
    integer(I4B), intent(in) :: nja !< number of GWT connections
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    integer(I4B), intent(in), dimension(nja) :: idxglo !< mapping vector for model (local) to solution (global)
    real(DP), intent(inout), dimension(nodes) :: rhs !< right-hand side vector for model
    ! -- local
    integer(I4B) :: n, idiag
    real(DP) :: tled
    real(DP) :: hhcof, rrhs
    real(DP) :: vnew, vold
    !
    ! -- set variables
    tled = DONE / delt
    !
    ! -- loop through and calculate storage contribution to hcof and rhs
    do n = 1, this%dis%nodes
      !
      ! -- skip if transport inactive
      if (this%ibound(n) <= 0) cycle
      !
      ! -- calculate new and old water volumes
      vnew = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n)) * &
             this%fmi%gwfsat(n) * this%thetam(n)
      vold = vnew
      if (this%fmi%igwfstrgss /= 0) vold = vold + this%fmi%gwfstrgss(n) * delt
      if (this%fmi%igwfstrgsy /= 0) vold = vold + this%fmi%gwfstrgsy(n) * delt
      !
      ! -- add terms to diagonal and rhs accumulators
      hhcof = -vnew * tled
      rrhs = -vold * tled * cold(n)
      idiag = this%dis%con%ia(n)
      call matrix_sln%add_value_pos(idxglo(idiag), hhcof)
      rhs(n) = rhs(n) + rrhs
    end do
  end subroutine mst_fc_sto

  !> @ brief Fill decay coefficient method for package
  !!
  !!  Method to calculate and fill decay coefficients for the package.
  !<
  subroutine mst_fc_dcy(this, nodes, cold, cnew, nja, matrix_sln, &
                        idxglo, rhs, kiter)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer, intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cold !< concentration at end of last time step
    real(DP), intent(in), dimension(nodes) :: cnew !< concentration at end of this time step
    integer(I4B), intent(in) :: nja !< number of GWT connections
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    integer(I4B), intent(in), dimension(nja) :: idxglo !< mapping vector for model (local) to solution (global)
    real(DP), intent(inout), dimension(nodes) :: rhs !< right-hand side vector for model
    integer(I4B), intent(in) :: kiter !< solution outer iteration number
    ! -- local
    integer(I4B) :: n, idiag
    real(DP) :: hhcof, rrhs
    real(DP) :: swtpdt
    real(DP) :: vcell
    real(DP) :: decay_rate
    !
    ! -- loop through and calculate decay contribution to hcof and rhs
    do n = 1, this%dis%nodes
      !
      ! -- skip if transport inactive
      if (this%ibound(n) <= 0) cycle
      !
      ! -- calculate new and old water volumes
      vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
      swtpdt = this%fmi%gwfsat(n)
      !
      ! -- add decay rate terms to accumulators
      idiag = this%dis%con%ia(n)
      select case (this%idcy)
      case (DECAY_FIRST_ORDER)
        !
        ! -- first order decay rate is a function of concentration, so add
        !    to left hand side
        if (cnew(n) > DZERO) then
          hhcof = -this%decay(n) * vcell * swtpdt * this%thetam(n)
          call matrix_sln%add_value_pos(idxglo(idiag), hhcof)
        end if
      case (DECAY_ZERO_ORDER)
        !
        ! -- Call function to get zero-order decay rate, which may be changed
        !    from the user-specified rate to prevent negative concentrations
        decay_rate = get_zero_order_decay(this%decay(n), this%decaylast(n), &
                                          kiter, cold(n), cnew(n), delt)
        this%decaylast(n) = decay_rate
        rrhs = decay_rate * vcell * swtpdt * this%thetam(n)
        rhs(n) = rhs(n) + rrhs
      end select
      !
    end do
  end subroutine mst_fc_dcy

  !> @ brief Fill sorption coefficient method for package
  !!
  !!  Method to calculate and fill sorption coefficients for the package.
  !<
  subroutine mst_fc_srb(this, nodes, cold, nja, matrix_sln, idxglo, rhs, &
                        cnew)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer, intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cold !< concentration at end of last time step
    integer(I4B), intent(in) :: nja !< number of GWT connections
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    integer(I4B), intent(in), dimension(nja) :: idxglo !< mapping vector for model (local) to solution (global)
    real(DP), intent(inout), dimension(nodes) :: rhs !< right-hand side vector for model
    real(DP), intent(in), dimension(nodes) :: cnew !< concentration at end of this time step
    ! -- local
    integer(I4B) :: n, idiag
    real(DP) :: tled
    real(DP) :: hhcof, rrhs
    real(DP) :: swt, swtpdt
    real(DP) :: vcell
    real(DP) :: const1
    real(DP) :: const2
    real(DP) :: volfracm
    real(DP) :: rhobm
    !
    ! -- set variables
    tled = DONE / delt
    !
    ! -- loop through and calculate sorption contribution to hcof and rhs
    do n = 1, this%dis%nodes
      !
      ! -- skip if transport inactive
      if (this%ibound(n) <= 0) cycle
      !
      ! -- assign variables
      vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
      swtpdt = this%fmi%gwfsat(n)
      swt = this%fmi%gwfsatold(n, delt)
      idiag = this%dis%con%ia(n)
      const1 = this%distcoef(n)
      const2 = 0.
      if (this%isrb == SORPTION_FREUND .or. this%isrb == SORPTION_LANG) then
        const2 = this%sp2(n)
      end if
      volfracm = this%get_volfracm(n)
      rhobm = this%bulk_density(n)
      call mst_srb_term(this%isrb, volfracm, rhobm, vcell, tled, cnew(n), &
                        cold(n), swtpdt, swt, const1, const2, &
                        hcofval=hhcof, rhsval=rrhs)
      !
      ! -- Add hhcof to diagonal and rrhs to right-hand side
      call matrix_sln%add_value_pos(idxglo(idiag), hhcof)
      rhs(n) = rhs(n) + rrhs
      !
    end do
  end subroutine mst_fc_srb

  !> @ brief Calculate sorption terms
  !!
  !!  Subroutine to calculate sorption terms
  !<
  subroutine mst_srb_term(isrb, volfracm, rhobm, vcell, tled, cnew, cold, &
                          swnew, swold, const1, const2, rate, hcofval, rhsval)
    ! -- dummy
    integer(I4B), intent(in) :: isrb !< sorption flag 1, 2, 3 are linear, freundlich, and langmuir
    real(DP), intent(in) :: volfracm !< volume fraction of mobile domain (fhat_m)
    real(DP), intent(in) :: rhobm !< bulk density of mobile domain (rhob_m)
    real(DP), intent(in) :: vcell !< volume of cell
    real(DP), intent(in) :: tled !< one over time step length
    real(DP), intent(in) :: cnew !< concentration at end of this time step
    real(DP), intent(in) :: cold !< concentration at end of last time step
    real(DP), intent(in) :: swnew !< cell saturation at end of this time step
    real(DP), intent(in) :: swold !< cell saturation at end of last time step
    real(DP), intent(in) :: const1 !< distribution coefficient or freundlich or langmuir constant
    real(DP), intent(in) :: const2 !< zero, freundlich exponent, or langmuir sorption sites
    real(DP), intent(out), optional :: rate !< calculated sorption rate
    real(DP), intent(out), optional :: hcofval !< diagonal contribution to solution coefficient matrix
    real(DP), intent(out), optional :: rhsval !< contribution to solution right-hand-side
    ! -- local
    real(DP) :: term
    real(DP) :: derv
    real(DP) :: cbarnew
    real(DP) :: cbarold
    real(DP) :: cavg
    real(DP) :: cbaravg
    real(DP) :: swavg
    !
    ! -- Calculate based on type of sorption
    if (isrb == SORPTION_LINEAR) then
      ! -- linear
      term = -volfracm * rhobm * vcell * tled * const1
      if (present(hcofval)) hcofval = term * swnew
      if (present(rhsval)) rhsval = term * swold * cold
      if (present(rate)) rate = term * swnew * cnew - term * swold * cold
    else
      !
      ! -- calculate average aqueous concentration
      cavg = DHALF * (cold + cnew)
      !
      ! -- set values based on isotherm
      select case (isrb)
      case (SORPTION_FREUND)
        ! -- freundlich
        cbarnew = get_freundlich_conc(cnew, const1, const2)
        cbarold = get_freundlich_conc(cold, const1, const2)
        derv = get_freundlich_derivative(cavg, const1, const2)
      case (SORPTION_LANG)
        ! -- langmuir
        cbarnew = get_langmuir_conc(cnew, const1, const2)
        cbarold = get_langmuir_conc(cold, const1, const2)
        derv = get_langmuir_derivative(cavg, const1, const2)
      end select
      !
      ! -- calculate hcof, rhs, and rate for freundlich and langmuir
      term = -volfracm * rhobm * vcell * tled
      cbaravg = (cbarold + cbarnew) * DHALF
      swavg = (swnew + swold) * DHALF
      if (present(hcofval)) then
        hcofval = term * derv * swavg
      end if
      if (present(rhsval)) then
        rhsval = term * derv * swavg * cold - term * cbaravg * (swnew - swold)
      end if
      if (present(rate)) then
        rate = term * derv * swavg * (cnew - cold) &
               + term * cbaravg * (swnew - swold)
      end if
    end if
  end subroutine mst_srb_term

  !> @ brief Fill sorption-decay coefficient method for package
  !!
  !!  Method to calculate and fill sorption-decay coefficients for the package.
  !<
  subroutine mst_fc_dcy_srb(this, nodes, cold, nja, matrix_sln, idxglo, &
                            rhs, cnew, kiter)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer, intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cold !< concentration at end of last time step
    integer(I4B), intent(in) :: nja !< number of GWT connections
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    integer(I4B), intent(in), dimension(nja) :: idxglo !< mapping vector for model (local) to solution (global)
    real(DP), intent(inout), dimension(nodes) :: rhs !< right-hand side vector for model
    real(DP), intent(in), dimension(nodes) :: cnew !< concentration at end of this time step
    integer(I4B), intent(in) :: kiter !< solution outer iteration number
    ! -- local
    integer(I4B) :: n, idiag
    real(DP) :: hhcof, rrhs
    real(DP) :: vcell
    real(DP) :: swnew
    real(DP) :: distcoef
    real(DP) :: volfracm
    real(DP) :: rhobm
    real(DP) :: term
    real(DP) :: csrb
    real(DP) :: decay_rate
    real(DP) :: csrbold
    real(DP) :: csrbnew
    !
    ! -- loop through and calculate sorption contribution to hcof and rhs
    do n = 1, this%dis%nodes
      !
      ! -- skip if transport inactive
      if (this%ibound(n) <= 0) cycle
      !
      ! -- set variables
      hhcof = DZERO
      rrhs = DZERO
      vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
      swnew = this%fmi%gwfsat(n)
      distcoef = this%distcoef(n)
      idiag = this%dis%con%ia(n)
      volfracm = this%get_volfracm(n)
      rhobm = this%bulk_density(n)
      term = this%decay_sorbed(n) * volfracm * rhobm * swnew * vcell
      !
      ! -- add sorbed mass decay rate terms to accumulators
      select case (this%idcy)
      case (DECAY_FIRST_ORDER)
        select case (this%isrb)
        case (SORPTION_LINEAR)
          !
          ! -- first order decay rate is a function of concentration, so add
          !    to left hand side
          if (cnew(n) > DZERO) then
            hhcof = -term * distcoef
          end if
        case (SORPTION_FREUND)
          !
          ! -- nonlinear Freundlich sorption, so add to RHS
          csrb = get_freundlich_conc(cnew(n), distcoef, this%sp2(n))
          rrhs = term * csrb
        case (SORPTION_LANG)
          !
          ! -- nonlinear Lanmuir sorption, so add to RHS
          csrb = get_langmuir_conc(cnew(n), distcoef, this%sp2(n))
          rrhs = term * csrb
        end select
      case (DECAY_ZERO_ORDER)
        !
        ! -- call function to get zero-order decay rate, which may be changed
        !    from the user-specified rate to prevent negative concentrations
        if (distcoef > DZERO) then
          select case (this%isrb)
          case (SORPTION_LINEAR)
            csrbold = cold(n) * distcoef
            csrbnew = cnew(n) * distcoef
          case (SORPTION_FREUND)
            csrbold = get_freundlich_conc(cold(n), distcoef, this%sp2(n))
            csrbnew = get_freundlich_conc(cnew(n), distcoef, this%sp2(n))
          case (SORPTION_LANG)
            csrbold = get_langmuir_conc(cold(n), distcoef, this%sp2(n))
            csrbnew = get_langmuir_conc(cnew(n), distcoef, this%sp2(n))
          end select
          !
          decay_rate = get_zero_order_decay(this%decay_sorbed(n), &
                                            this%decayslast(n), &
                                            kiter, csrbold, csrbnew, delt)
          this%decayslast(n) = decay_rate
          rrhs = decay_rate * volfracm * rhobm * swnew * vcell
        end if
      end select
      !
      ! -- Add hhcof to diagonal and rrhs to right-hand side
      call matrix_sln%add_value_pos(idxglo(idiag), hhcof)
      rhs(n) = rhs(n) + rrhs
      !
    end do
  end subroutine mst_fc_dcy_srb

  !> @ brief Calculate flows for package
  !!
  !!  Method to calculate flows for the package.
  !<
  subroutine mst_cq(this, nodes, cnew, cold, flowja)
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer(I4B), intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cnew !< concentration at end of this time step
    real(DP), intent(in), dimension(nodes) :: cold !< concentration at end of last time step
    real(DP), dimension(:), contiguous, intent(inout) :: flowja !< flow between two connected control volumes
    !
    ! - storage
    call this%mst_cq_sto(nodes, cnew, cold, flowja)
    !
    ! -- decay
    if (this%idcy /= DECAY_OFF) then
      call this%mst_cq_dcy(nodes, cnew, cold, flowja)
    end if
    !
    ! -- sorption
    if (this%isrb /= SORPTION_OFF) then
      call this%mst_cq_srb(nodes, cnew, cold, flowja)
    end if
    !
    ! -- decay sorbed
    if (this%isrb /= SORPTION_OFF .and. this%idcy /= DECAY_OFF) then
      call this%mst_cq_dcy_srb(nodes, cnew, cold, flowja)
    end if
    !
    ! -- calculate csrb
    if (this%isrb /= SORPTION_OFF) then
      call this%mst_calc_csrb(cnew)
    end if
  end subroutine mst_cq

  !> @ brief Calculate storage terms for package
  !!
  !!  Method to calculate storage terms for the package.
  !<
  subroutine mst_cq_sto(this, nodes, cnew, cold, flowja)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer(I4B), intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cnew !< concentration at end of this time step
    real(DP), intent(in), dimension(nodes) :: cold !< concentration at end of last time step
    real(DP), dimension(:), contiguous, intent(inout) :: flowja !< flow between two connected control volumes
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: idiag
    real(DP) :: rate
    real(DP) :: tled
    real(DP) :: vnew, vold
    real(DP) :: hhcof, rrhs
    !
    ! -- initialize
    tled = DONE / delt
    !
    ! -- Calculate storage change
    do n = 1, nodes
      this%ratesto(n) = DZERO
      !
      ! -- skip if transport inactive
      if (this%ibound(n) <= 0) cycle
      !
      ! -- calculate new and old water volumes
      vnew = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n)) * &
             this%fmi%gwfsat(n) * this%thetam(n)
      vold = vnew
      if (this%fmi%igwfstrgss /= 0) vold = vold + this%fmi%gwfstrgss(n) * delt
      if (this%fmi%igwfstrgsy /= 0) vold = vold + this%fmi%gwfstrgsy(n) * delt
      !
      ! -- calculate rate
      hhcof = -vnew * tled
      rrhs = -vold * tled * cold(n)
      rate = hhcof * cnew(n) - rrhs
      this%ratesto(n) = rate
      idiag = this%dis%con%ia(n)
      flowja(idiag) = flowja(idiag) + rate
    end do
  end subroutine mst_cq_sto

  !> @ brief Calculate decay terms for package
  !!
  !!  Method to calculate decay terms for the package.
  !<
  subroutine mst_cq_dcy(this, nodes, cnew, cold, flowja)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer(I4B), intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cnew !< concentration at end of this time step
    real(DP), intent(in), dimension(nodes) :: cold !< concentration at end of last time step
    real(DP), dimension(:), contiguous, intent(inout) :: flowja !< flow between two connected control volumes
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: idiag
    real(DP) :: rate
    real(DP) :: swtpdt
    real(DP) :: hhcof, rrhs
    real(DP) :: vcell
    real(DP) :: decay_rate
    !
    ! -- initialize
    !
    ! -- Calculate decay change
    do n = 1, nodes
      !
      ! -- skip if transport inactive
      this%ratedcy(n) = DZERO
      if (this%ibound(n) <= 0) cycle
      !
      ! -- calculate new and old water volumes
      vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
      swtpdt = this%fmi%gwfsat(n)
      !
      ! -- calculate decay gains and losses
      rate = DZERO
      hhcof = DZERO
      rrhs = DZERO
      if (this%idcy == DECAY_FIRST_ORDER) then
        if (cnew(n) > DZERO) then
          hhcof = -this%decay(n) * vcell * swtpdt * this%thetam(n)
        end if
      elseif (this%idcy == DECAY_ZERO_ORDER) then
        decay_rate = get_zero_order_decay(this%decay(n), this%decaylast(n), &
                                          0, cold(n), cnew(n), delt)
        rrhs = decay_rate * vcell * swtpdt * this%thetam(n)
      end if
      rate = hhcof * cnew(n) - rrhs
      this%ratedcy(n) = rate
      idiag = this%dis%con%ia(n)
      flowja(idiag) = flowja(idiag) + rate
      !
    end do
  end subroutine mst_cq_dcy

  !> @ brief Calculate sorption terms for package
  !!
  !!  Method to calculate sorption terms for the package.
  !<
  subroutine mst_cq_srb(this, nodes, cnew, cold, flowja)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer(I4B), intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cnew !< concentration at end of this time step
    real(DP), intent(in), dimension(nodes) :: cold !< concentration at end of last time step
    real(DP), dimension(:), contiguous, intent(inout) :: flowja !< flow between two connected control volumes
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: idiag
    real(DP) :: rate
    real(DP) :: tled
    real(DP) :: swt, swtpdt
    real(DP) :: vcell
    real(DP) :: volfracm
    real(DP) :: rhobm
    real(DP) :: const1
    real(DP) :: const2
    !
    ! -- initialize
    tled = DONE / delt
    !
    ! -- Calculate sorption change
    do n = 1, nodes
      !
      ! -- initialize rates
      this%ratesrb(n) = DZERO
      !
      ! -- skip if transport inactive
      if (this%ibound(n) <= 0) cycle
      !
      ! -- assign variables
      vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
      swtpdt = this%fmi%gwfsat(n)
      swt = this%fmi%gwfsatold(n, delt)
      volfracm = this%get_volfracm(n)
      rhobm = this%bulk_density(n)
      const1 = this%distcoef(n)
      const2 = 0.
      if (this%isrb == SORPTION_FREUND .or. this%isrb == SORPTION_LANG) then
        const2 = this%sp2(n)
      end if
      call mst_srb_term(this%isrb, volfracm, rhobm, vcell, tled, cnew(n), &
                        cold(n), swtpdt, swt, const1, const2, &
                        rate=rate)
      this%ratesrb(n) = rate
      idiag = this%dis%con%ia(n)
      flowja(idiag) = flowja(idiag) + rate
      !
    end do
  end subroutine mst_cq_srb

  !> @ brief Calculate decay-sorption terms for package
  !!
  !!  Method to calculate decay-sorption terms for the package.
  !<
  subroutine mst_cq_dcy_srb(this, nodes, cnew, cold, flowja)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer(I4B), intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cnew !< concentration at end of this time step
    real(DP), intent(in), dimension(nodes) :: cold !< concentration at end of last time step
    real(DP), dimension(:), contiguous, intent(inout) :: flowja !< flow between two connected control volumes
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: idiag
    real(DP) :: rate
    real(DP) :: hhcof, rrhs
    real(DP) :: vcell
    real(DP) :: swnew
    real(DP) :: distcoef
    real(DP) :: volfracm
    real(DP) :: rhobm
    real(DP) :: term
    real(DP) :: csrb
    real(DP) :: csrbnew
    real(DP) :: csrbold
    real(DP) :: decay_rate
    !
    ! -- Calculate sorbed decay change
    !    This routine will only be called if sorption and decay are active
    do n = 1, nodes
      !
      ! -- initialize rates
      this%ratedcys(n) = DZERO
      !
      ! -- skip if transport inactive
      if (this%ibound(n) <= 0) cycle
      !
      ! -- set variables
      hhcof = DZERO
      rrhs = DZERO
      vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
      swnew = this%fmi%gwfsat(n)
      distcoef = this%distcoef(n)
      volfracm = this%get_volfracm(n)
      rhobm = this%bulk_density(n)
      term = this%decay_sorbed(n) * volfracm * rhobm * swnew * vcell
      !
      ! -- add sorbed mass decay rate terms to accumulators
      select case (this%idcy)
      case (DECAY_FIRST_ORDER)

        !
        select case (this%isrb)
        case (SORPTION_LINEAR)
          !
          ! -- first order decay rate is a function of concentration, so add
          !    to left hand side
          if (cnew(n) > DZERO) then
            hhcof = -term * distcoef
          end if
        case (SORPTION_FREUND)
          !
          ! -- nonlinear Freundlich sorption, so add to RHS
          csrb = get_freundlich_conc(cnew(n), distcoef, this%sp2(n))
          rrhs = term * csrb
        case (SORPTION_LANG)
          !
          ! -- nonlinear Lanmuir sorption, so add to RHS
          csrb = get_langmuir_conc(cnew(n), distcoef, this%sp2(n))
          rrhs = term * csrb
        end select
      case (DECAY_ZERO_ORDER)
        !
        ! -- Call function to get zero-order decay rate, which may be changed
        !    from the user-specified rate to prevent negative concentrations
        if (distcoef > DZERO) then
          select case (this%isrb)
          case (SORPTION_LINEAR)
            csrbold = cold(n) * distcoef
            csrbnew = cnew(n) * distcoef
          case (SORPTION_FREUND)
            csrbold = get_freundlich_conc(cold(n), distcoef, this%sp2(n))
            csrbnew = get_freundlich_conc(cnew(n), distcoef, this%sp2(n))
          case (SORPTION_LANG)
            csrbold = get_langmuir_conc(cold(n), distcoef, this%sp2(n))
            csrbnew = get_langmuir_conc(cnew(n), distcoef, this%sp2(n))
          end select
          decay_rate = get_zero_order_decay(this%decay_sorbed(n), &
                                            this%decayslast(n), &
                                            0, csrbold, csrbnew, delt)
          rrhs = decay_rate * volfracm * rhobm * swnew * vcell
        end if
      end select
      !
      ! -- calculate rate
      rate = hhcof * cnew(n) - rrhs
      this%ratedcys(n) = rate
      idiag = this%dis%con%ia(n)
      flowja(idiag) = flowja(idiag) + rate
      !
    end do
  end subroutine mst_cq_dcy_srb

  !> @ brief Calculate sorbed concentration
  !<
  subroutine mst_calc_csrb(this, cnew)
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    real(DP), intent(in), dimension(:) :: cnew !< concentration at end of this time step
    ! -- local
    integer(I4B) :: n
    real(DP) :: distcoef
    real(DP) :: csrb

    ! Calculate sorbed concentration
    do n = 1, size(cnew)
      csrb = DZERO
      if (this%ibound(n) > 0) then
        distcoef = this%distcoef(n)
        select case (this%isrb)
        case (SORPTION_LINEAR)
          csrb = cnew(n) * distcoef
        case (SORPTION_FREUND)
          csrb = get_freundlich_conc(cnew(n), distcoef, this%sp2(n))
        case (SORPTION_LANG)
          csrb = get_langmuir_conc(cnew(n), distcoef, this%sp2(n))
        end select
      end if
      this%csrb(n) = csrb
    end do

  end subroutine mst_calc_csrb

  !> @ brief Calculate budget terms for package
  !<
  subroutine mst_bd(this, isuppress_output, model_budget)
    ! -- modules
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType, rate_accumulator
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer(I4B), intent(in) :: isuppress_output !< flag to suppress output
    type(BudgetType), intent(inout) :: model_budget !< model budget object
    ! -- local
    real(DP) :: rin
    real(DP) :: rout
    !
    ! -- sto
    call rate_accumulator(this%ratesto, rin, rout)
    call model_budget%addentry(rin, rout, delt, budtxt(1), &
                               isuppress_output, rowlabel=this%packName)
    !
    ! -- dcy
    if (this%idcy /= DECAY_OFF) then
      call rate_accumulator(this%ratedcy, rin, rout)
      call model_budget%addentry(rin, rout, delt, budtxt(2), &
                                 isuppress_output, rowlabel=this%packName)
    end if
    !
    ! -- srb
    if (this%isrb /= SORPTION_OFF) then
      call rate_accumulator(this%ratesrb, rin, rout)
      call model_budget%addentry(rin, rout, delt, budtxt(3), &
                                 isuppress_output, rowlabel=this%packName)
    end if
    !
    ! -- srb dcy
    if (this%isrb /= SORPTION_OFF .and. this%idcy /= DECAY_OFF) then
      call rate_accumulator(this%ratedcys, rin, rout)
      call model_budget%addentry(rin, rout, delt, budtxt(4), &
                                 isuppress_output, rowlabel=this%packName)
    end if
  end subroutine mst_bd

  !> @ brief Output flow terms for package
  !!
  !!  Method to output terms for the package.
  !<
  subroutine mst_ot_flow(this, icbcfl, icbcun)
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer(I4B), intent(in) :: icbcfl !< flag and unit number for cell-by-cell output
    integer(I4B), intent(in) :: icbcun !< flag indication if cell-by-cell data should be saved
    ! -- local
    integer(I4B) :: ibinun
    integer(I4B) :: iprint, nvaluesp, nwidthp
    character(len=1) :: cdatafmp = ' ', editdesc = ' '
    real(DP) :: dinact
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
    ! -- Record the storage rate if requested
    if (ibinun /= 0) then
      iprint = 0
      dinact = DZERO
      !
      ! -- sto
      call this%dis%record_array(this%ratesto, this%iout, iprint, -ibinun, &
                                 budtxt(1), cdatafmp, nvaluesp, &
                                 nwidthp, editdesc, dinact)
      !
      ! -- dcy
      if (this%idcy /= DECAY_OFF) &
        call this%dis%record_array(this%ratedcy, this%iout, iprint, -ibinun, &
                                   budtxt(2), cdatafmp, nvaluesp, &
                                   nwidthp, editdesc, dinact)
      !
      ! -- srb
      if (this%isrb /= SORPTION_OFF) &
        call this%dis%record_array(this%ratesrb, this%iout, iprint, -ibinun, &
                                   budtxt(3), cdatafmp, nvaluesp, &
                                   nwidthp, editdesc, dinact)
      !
      ! -- dcy srb
      if (this%isrb /= SORPTION_OFF .and. this%idcy /= DECAY_OFF) &
        call this%dis%record_array(this%ratedcys, this%iout, iprint, -ibinun, &
                                   budtxt(4), cdatafmp, nvaluesp, &
                                   nwidthp, editdesc, dinact)
    end if
  end subroutine mst_ot_flow

  !> @brief Save sorbate concentration array to binary file
  !<
  subroutine mst_ot_dv(this, idvsave)
    ! -- dummy
    class(GwtMstType) :: this
    integer(I4B), intent(in) :: idvsave
    ! -- local
    character(len=1) :: cdatafmp = ' ', editdesc = ' '
    integer(I4B) :: ibinun
    integer(I4B) :: iprint
    integer(I4B) :: nvaluesp
    integer(I4B) :: nwidthp
    real(DP) :: dinact

    ! Set unit number for sorbate output
    if (this%ioutsorbate /= 0) then
      ibinun = 1
    else
      ibinun = 0
    end if
    if (idvsave == 0) ibinun = 0

    ! save sorbate concentration array
    if (ibinun /= 0) then
      iprint = 0
      dinact = DHNOFLO
      if (this%ioutsorbate /= 0) then
        ibinun = this%ioutsorbate
        call this%dis%record_array(this%csrb, this%iout, iprint, ibinun, &
                                   '         SORBATE', cdatafmp, nvaluesp, &
                                   nwidthp, editdesc, dinact)
      end if
    end if

  end subroutine mst_ot_dv

  !> @ brief Deallocate
  !!
  !!  Method to deallocate memory for the package.
  !<
  subroutine mst_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    !
    ! -- Deallocate arrays if package was active
    if (this%inunit > 0) then
      call mem_deallocate(this%porosity)
      call mem_deallocate(this%thetam)
      call mem_deallocate(this%volfracim)
      call mem_deallocate(this%ratesto)
      call mem_deallocate(this%idcy)
      call mem_deallocate(this%decay)
      call mem_deallocate(this%decay_sorbed)
      call mem_deallocate(this%ratedcy)
      call mem_deallocate(this%decaylast)
      call mem_deallocate(this%decayslast)
      call mem_deallocate(this%isrb)
      call mem_deallocate(this%ioutsorbate)
      call mem_deallocate(this%bulk_density)
      call mem_deallocate(this%distcoef)
      call mem_deallocate(this%sp2)
      call mem_deallocate(this%ratesrb)
      call mem_deallocate(this%csrb)
      call mem_deallocate(this%ratedcys)
      this%ibound => null()
      this%fmi => null()
    end if
    !
    ! -- Scalars
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()
  end subroutine mst_da

  !> @ brief Allocate scalar variables for package
  !!
  !!  Method to allocate scalar variables for the package.
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    !
    ! -- Allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%isrb, 'ISRB', this%memoryPath)
    call mem_allocate(this%ioutsorbate, 'IOUTSORBATE', this%memoryPath)
    call mem_allocate(this%idcy, 'IDCY', this%memoryPath)
    !
    ! -- Initialize
    this%isrb = IZERO
    this%ioutsorbate = 0
    this%idcy = IZERO
  end subroutine allocate_scalars

  !> @ brief Allocate arrays for package
  !!
  !!  Method to allocate arrays for the package.
  !<
  subroutine allocate_arrays(this, nodes)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use ConstantsModule, only: DZERO
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer(I4B), intent(in) :: nodes !< number of nodes
    ! -- local
    integer(I4B) :: n
    !
    ! -- Allocate
    ! -- sto
    call mem_allocate(this%porosity, nodes, 'POROSITY', this%memoryPath)
    call mem_allocate(this%thetam, nodes, 'THETAM', this%memoryPath)
    call mem_allocate(this%volfracim, nodes, 'VOLFRACIM', this%memoryPath)
    call mem_allocate(this%ratesto, nodes, 'RATESTO', this%memoryPath)
    !
    ! -- dcy
    if (this%idcy == DECAY_OFF) then
      call mem_allocate(this%ratedcy, 1, 'RATEDCY', this%memoryPath)
      call mem_allocate(this%decay, 1, 'DECAY', this%memoryPath)
      call mem_allocate(this%decaylast, 1, 'DECAYLAST', this%memoryPath)
    else
      call mem_allocate(this%ratedcy, this%dis%nodes, 'RATEDCY', this%memoryPath)
      call mem_allocate(this%decay, nodes, 'DECAY', this%memoryPath)
      call mem_allocate(this%decaylast, nodes, 'DECAYLAST', this%memoryPath)
    end if
    if (this%idcy /= DECAY_OFF .and. this%isrb /= SORPTION_OFF) then
      call mem_allocate(this%ratedcys, this%dis%nodes, 'RATEDCYS', &
                        this%memoryPath)
      call mem_allocate(this%decayslast, this%dis%nodes, 'DECAYSLAST', &
                        this%memoryPath)
    else
      call mem_allocate(this%ratedcys, 1, 'RATEDCYS', this%memoryPath)
      call mem_allocate(this%decayslast, 1, 'DECAYSLAST', this%memoryPath)
    end if
    call mem_allocate(this%decay_sorbed, 1, 'DECAY_SORBED', &
                      this%memoryPath)
    !
    ! -- srb
    if (this%isrb == SORPTION_OFF) then
      call mem_allocate(this%bulk_density, 1, 'BULK_DENSITY', this%memoryPath)
      call mem_allocate(this%sp2, 1, 'SP2', this%memoryPath)
      call mem_allocate(this%distcoef, 1, 'DISTCOEF', this%memoryPath)
      call mem_allocate(this%ratesrb, 1, 'RATESRB', this%memoryPath)
      call mem_allocate(this%csrb, 1, 'CSRB', this%memoryPath)
    else
      call mem_allocate(this%bulk_density, nodes, 'BULK_DENSITY', this%memoryPath)
      call mem_allocate(this%distcoef, nodes, 'DISTCOEF', this%memoryPath)
      call mem_allocate(this%ratesrb, nodes, 'RATESRB', this%memoryPath)
      call mem_allocate(this%csrb, nodes, 'CSRB', this%memoryPath)
      if (this%isrb == SORPTION_LINEAR) then
        call mem_allocate(this%sp2, 1, 'SP2', this%memoryPath)
      else
        call mem_allocate(this%sp2, nodes, 'SP2', this%memoryPath)
      end if
    end if
    !
    ! -- Initialize
    do n = 1, nodes
      this%porosity(n) = DZERO
      this%thetam(n) = DZERO
      this%volfracim(n) = DZERO
      this%ratesto(n) = DZERO
    end do
    do n = 1, size(this%decay)
      this%decay(n) = DZERO
      this%ratedcy(n) = DZERO
      this%decaylast(n) = DZERO
    end do
    do n = 1, size(this%bulk_density)
      this%bulk_density(n) = DZERO
      this%distcoef(n) = DZERO
      this%ratesrb(n) = DZERO
      this%csrb(n) = DZERO
    end do
    do n = 1, size(this%sp2)
      this%sp2(n) = DZERO
    end do
    do n = 1, size(this%ratedcys)
      this%ratedcys(n) = DZERO
      this%decayslast(n) = DZERO
    end do
  end subroutine allocate_arrays

  !> @ brief Read options for package
  !!
  !!  Method to read options for the package.
  !<
  subroutine read_options(this)
    ! -- modules
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: getunit, openfile
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    ! -- local
    character(len=LINELENGTH) :: keyword
    character(len=LINELENGTH) :: sorption
    character(len=MAXCHARLEN) :: fname
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    ! -- formats
    character(len=*), parameter :: fmtisvflow = &
      "(4x,'CELL-BY-CELL FLOW INFORMATION WILL BE SAVED TO BINARY FILE &
      &WHENEVER ICBCFL IS NOT ZERO.')"
    character(len=*), parameter :: fmtlinear = &
      &"(4x,'LINEAR SORPTION IS ACTIVE. ')"
    character(len=*), parameter :: fmtfreundlich = &
      &"(4x,'FREUNDLICH SORPTION IS ACTIVE. ')"
    character(len=*), parameter :: fmtlangmuir = &
      &"(4x,'LANGMUIR SORPTION IS ACTIVE. ')"
    character(len=*), parameter :: fmtidcy1 = &
      &"(4x,'FIRST-ORDER DECAY IS ACTIVE. ')"
    character(len=*), parameter :: fmtidcy2 = &
      &"(4x,'ZERO-ORDER DECAY IS ACTIVE. ')"
    character(len=*), parameter :: fmtfileout = &
      "(4x,'MST ',1x,a,1x,' WILL BE SAVED TO FILE: ',a,/4x,&
      &'OPENED ON UNIT: ',I7)"
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, blockRequired=.false., &
                              supportOpenClose=.true.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING MOBILE STORAGE AND TRANSFER OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('SAVE_FLOWS')
          this%ipakcb = -1
          write (this%iout, fmtisvflow)
        case ('SORBTION')
          call store_error('SORBTION is not a valid option.  Use &
                           &SORPTION instead.')
          call this%parser%StoreErrorUnit()
        case ('SORPTION')
          call this%parser%GetStringCaps(sorption)
          select case (sorption)
          case ('LINEAR', '')
            this%isrb = SORPTION_LINEAR
            write (this%iout, fmtlinear)
          case ('FREUNDLICH')
            this%isrb = SORPTION_FREUND
            write (this%iout, fmtfreundlich)
          case ('LANGMUIR')
            this%isrb = SORPTION_LANG
            write (this%iout, fmtlangmuir)
          case default
            call store_error('Unknown sorption type was specified ' &
                             & //'('//trim(adjustl(sorption))//').'// &
                             &' Sorption must be specified as LINEAR, &
                             &FREUNDLICH, or LANGMUIR.')
            call this%parser%StoreErrorUnit()
          end select
        case ('FIRST_ORDER_DECAY')
          this%idcy = DECAY_FIRST_ORDER
          write (this%iout, fmtidcy1)
        case ('ZERO_ORDER_DECAY')
          this%idcy = DECAY_ZERO_ORDER
          write (this%iout, fmtidcy2)
        case ('SORBATE')
          call this%parser%GetStringCaps(keyword)
          if (keyword == 'FILEOUT') then
            call this%parser%GetString(fname)
            this%ioutsorbate = getunit()
            call openfile(this%ioutsorbate, this%iout, fname, 'DATA(BINARY)', &
                          form, access, 'REPLACE', mode_opt=MNORMAL)
            write (this%iout, fmtfileout) &
              'SORBATE', fname, this%ioutsorbate
          else
            errmsg = 'Optional SORBATE keyword must be '// &
                     'followed by FILEOUT.'
            call store_error(errmsg)
          end if
        case default
          write (errmsg, '(a,a)') 'Unknown MST option: ', trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (this%iout, '(1x,a)') 'END OF MOBILE STORAGE AND TRANSFER OPTIONS'
    end if
  end subroutine read_options

  !> @ brief Read data for package
  !!
  !!  Method to read data for the package.
  !<
  subroutine read_data(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use MemoryManagerModule, only: mem_reallocate, mem_reassignptr
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    ! -- local
    character(len=LINELENGTH) :: keyword
    character(len=:), allocatable :: line
    integer(I4B) :: istart, istop, lloc, ierr, n
    logical :: isfound, endOfBlock
    logical, dimension(6) :: lname
    character(len=24), dimension(6) :: aname
    ! -- data
    data aname(1)/'  MOBILE DOMAIN POROSITY'/
    data aname(2)/'            BULK DENSITY'/
    data aname(3)/'DISTRIBUTION COEFFICIENT'/
    data aname(4)/'              DECAY RATE'/
    data aname(5)/'       DECAY SORBED RATE'/
    data aname(6)/'   SECOND SORPTION PARAM'/
    !
    ! -- initialize
    isfound = .false.
    lname(:) = .false.
    !
    ! -- get griddata block
    call this%parser%GetBlock('GRIDDATA', isfound, ierr)
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING GRIDDATA'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        call this%parser%GetRemainingLine(line)
        lloc = 1
        select case (keyword)
        case ('POROSITY')
          call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                        this%parser%iuactive, this%porosity, &
                                        aname(1))
          lname(1) = .true.
        case ('BULK_DENSITY')
          if (this%isrb == SORPTION_OFF) &
            call mem_reallocate(this%bulk_density, this%dis%nodes, &
                                'BULK_DENSITY', trim(this%memoryPath))
          call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                        this%parser%iuactive, &
                                        this%bulk_density, aname(2))
          lname(2) = .true.
        case ('DISTCOEF')
          if (this%isrb == SORPTION_OFF) &
            call mem_reallocate(this%distcoef, this%dis%nodes, 'DISTCOEF', &
                                trim(this%memoryPath))
          call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                        this%parser%iuactive, this%distcoef, &
                                        aname(3))
          lname(3) = .true.
        case ('DECAY')
          if (this%idcy == DECAY_OFF) &
            call mem_reallocate(this%decay, this%dis%nodes, 'DECAY', &
                                trim(this%memoryPath))
          call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                        this%parser%iuactive, this%decay, &
                                        aname(4))
          lname(4) = .true.
        case ('DECAY_SORBED')
          call mem_reallocate(this%decay_sorbed, this%dis%nodes, &
                              'DECAY_SORBED', trim(this%memoryPath))
          call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                        this%parser%iuactive, &
                                        this%decay_sorbed, aname(5))
          lname(5) = .true.
        case ('SP2')
          if (this%isrb == SORPTION_OFF .or. this%isrb == SORPTION_LINEAR) &
            call mem_reallocate(this%sp2, this%dis%nodes, 'SP2', &
                                trim(this%memoryPath))
          call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                        this%parser%iuactive, this%sp2, &
                                        aname(6))
          lname(6) = .true.
        case default
          write (errmsg, '(a,a)') 'Unknown GRIDDATA tag: ', trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (this%iout, '(1x,a)') 'END PROCESSING GRIDDATA'
    else
      write (errmsg, '(a)') 'Required GRIDDATA block not found.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- Check for required porosity
    if (.not. lname(1)) then
      write (errmsg, '(a)') 'POROSITY not specified in GRIDDATA block.'
      call store_error(errmsg)
    end if
    !
    ! -- Check for required sorption variables
    if (this%isrb == SORPTION_LINEAR .or. this%isrb == SORPTION_FREUND .or. &
        this%isrb == SORPTION_LANG) then
      if (.not. lname(2)) then
        write (errmsg, '(a)') 'Sorption is active but BULK_DENSITY &
          &not specified.  BULK_DENSITY must be specified in GRIDDATA block.'
        call store_error(errmsg)
      end if
      if (.not. lname(3)) then
        write (errmsg, '(a)') 'Sorption is active but distribution &
          &coefficient not specified.  DISTCOEF must be specified in &
          &GRIDDATA block.'
        call store_error(errmsg)
      end if
      if (this%isrb == SORPTION_FREUND .or. this%isrb == SORPTION_LANG) then
        if (.not. lname(6)) then
          write (errmsg, '(a)') 'Freundlich or langmuir sorption is active &
            &but SP2 not specified.  SP2 must be specified in &
            &GRIDDATA block.'
          call store_error(errmsg)
        end if
      end if
    else
      if (lname(2)) then
        write (warnmsg, '(a)') 'Sorption is not active but &
          &BULK_DENSITY was specified.  BULK_DENSITY will have no affect on &
          &simulation results.'
        call store_warning(warnmsg)
        write (this%iout, '(1x,a)') 'WARNING.  '//warnmsg
      end if
      if (lname(3)) then
        write (warnmsg, '(a)') 'Sorption is not active but &
          &distribution coefficient was specified.  DISTCOEF will have &
          &no affect on simulation results.'
        call store_warning(warnmsg)
        write (this%iout, '(1x,a)') 'WARNING.  '//warnmsg
      end if
      if (lname(6)) then
        write (warnmsg, '(a)') 'Sorption is not active but &
          &SP2 was specified.  SP2 will have &
          &no affect on simulation results.'
        call store_warning(warnmsg)
        write (this%iout, '(1x,a)') 'WARNING.  '//warnmsg
      end if
    end if
    !
    ! -- Check for required decay/production rate coefficients
    if (this%idcy /= DECAY_OFF) then
      if (.not. lname(4)) then
        write (errmsg, '(a)') 'First or zero order decay is &
          &active but the first rate coefficient is not specified.  DECAY &
          &must be specified in GRIDDATA block.'
        call store_error(errmsg)
      end if
      if (.not. lname(5)) then
        !
        ! -- If DECAY_SORBED not specified and sorption is active, then
        !    terminate with an error
        if (this%isrb == SORPTION_LINEAR .or. this%isrb == SORPTION_FREUND .or. &
            this%isrb == SORPTION_LANG) then
          write (errmsg, '(a)') 'DECAY_SORBED not provided in GRIDDATA &
            &block but decay and sorption are active.  Specify DECAY_SORBED &
            &in GRIDDATA block.'
          call store_error(errmsg)
        end if
      end if
    else
      if (lname(4)) then
        write (warnmsg, '(a)') 'First- or zero-order decay &
          &is not active but decay was specified.  DECAY will &
          &have no affect on simulation results.'
        call store_warning(warnmsg)
        write (this%iout, '(1x,a)') 'WARNING.  '//warnmsg
      end if
      if (lname(5)) then
        write (warnmsg, '(a)') 'First- or zero-order decay &
          &is not active but DECAY_SORBED was specified.  &
          &DECAY_SORBED will have no affect on simulation results.'
        call store_warning(warnmsg)
        write (this%iout, '(1x,a)') 'WARNING.  '//warnmsg
      end if
    end if
    !
    ! -- terminate if errors
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- initialize thetam from porosity
    do n = 1, size(this%porosity)
      this%thetam(n) = this%porosity(n)
    end do
  end subroutine read_data

  !> @ brief Add volfrac values to volfracim
  !!
  !!  Method to add immobile domain volume fracions, which are stored as a
  !!  cumulative value in volfracim.
  !<
  subroutine addto_volfracim(this, volfracim)
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    real(DP), dimension(:), intent(in) :: volfracim !< immobile domain volume fraction that contributes to total immobile volume fraction
    ! -- local
    integer(I4B) :: n
    !
    ! -- Add to volfracim
    do n = 1, this%dis%nodes
      this%volfracim(n) = this%volfracim(n) + volfracim(n)
    end do
    !
    ! -- An immobile domain is adding a volume fraction, so update thetam
    !    accordingly.
    do n = 1, this%dis%nodes
      this%thetam(n) = this%get_volfracm(n) * this%porosity(n)
    end do
  end subroutine addto_volfracim

  !> @ brief Return mobile domain volume fraction
  !!
  !!  Calculate and return the volume fraction of the aquifer that is mobile
  !<
  function get_volfracm(this, node) result(volfracm)
    ! -- dummy
    class(GwtMstType) :: this !< GwtMstType object
    integer(I4B), intent(in) :: node !< node number
    ! -- return
    real(DP) :: volfracm
    !
    volfracm = DONE - this%volfracim(node)
  end function get_volfracm

  !> @ brief Calculate sorption concentration using Freundlich
  !!
  !!  Function to calculate sorption concentration using Freundlich
  !<
  function get_freundlich_conc(conc, kf, a) result(cbar)
    ! -- dummy
    real(DP), intent(in) :: conc !< solute concentration
    real(DP), intent(in) :: kf !< freundlich constant
    real(DP), intent(in) :: a !< freundlich exponent
    ! -- return
    real(DP) :: cbar
    !
    if (conc > DZERO) then
      cbar = kf * conc**a
    else
      cbar = DZERO
    end if
  end function

  !> @ brief Calculate sorption concentration using Langmuir
  !!
  !!  Function to calculate sorption concentration using Langmuir
  !<
  function get_langmuir_conc(conc, kl, sbar) result(cbar)
    ! -- dummy
    real(DP), intent(in) :: conc !< solute concentration
    real(DP), intent(in) :: kl !< langmuir constant
    real(DP), intent(in) :: sbar !< langmuir sorption sites
    ! -- return
    real(DP) :: cbar
    !
    if (conc > DZERO) then
      cbar = (kl * sbar * conc) / (DONE + kl * conc)
    else
      cbar = DZERO
    end if
  end function

  !> @ brief Calculate sorption derivative using Freundlich
  !!
  !!  Function to calculate sorption derivative using Freundlich
  !<
  function get_freundlich_derivative(conc, kf, a) result(derv)
    ! -- dummy
    real(DP), intent(in) :: conc !< solute concentration
    real(DP), intent(in) :: kf !< freundlich constant
    real(DP), intent(in) :: a !< freundlich exponent
    ! -- return
    real(DP) :: derv
    !
    if (conc > DZERO) then
      derv = kf * a * conc**(a - DONE)
    else
      derv = DZERO
    end if
  end function

  !> @ brief Calculate sorption derivative using Langmuir
  !!
  !!  Function to calculate sorption derivative using Langmuir
  !<
  function get_langmuir_derivative(conc, kl, sbar) result(derv)
    ! -- dummy
    real(DP), intent(in) :: conc !< solute concentration
    real(DP), intent(in) :: kl !< langmuir constant
    real(DP), intent(in) :: sbar !< langmuir sorption sites
    ! -- return
    real(DP) :: derv
    !
    if (conc > DZERO) then
      derv = (kl * sbar) / (DONE + kl * conc)**DTWO
    else
      derv = DZERO
    end if
  end function

  !> @ brief Get effective Freundlich distribution coefficient
  !<
  function get_freundlich_kd(conc, kf, a) result(kd)
    ! -- dummy
    real(DP), intent(in) :: conc !< solute concentration
    real(DP), intent(in) :: kf !< freundlich constant
    real(DP), intent(in) :: a !< freundlich exponent
    ! -- return
    real(DP) :: kd !< effective distribution coefficient
    !
    if (conc > DZERO) then
      kd = kf * conc**(a - DONE)
    else
      kd = DZERO
    end if
  end function get_freundlich_kd

  !> @ brief Get effective Langmuir distribution coefficient
  !<
  function get_langmuir_kd(conc, kl, sbar) result(kd)
    ! -- dummy
    real(DP), intent(in) :: conc !< solute concentration
    real(DP), intent(in) :: kl !< langmuir constant
    real(DP), intent(in) :: sbar !< langmuir sorption sites
    ! -- return
    real(DP) :: kd !< effective distribution coefficient
    !
    if (conc > DZERO) then
      kd = (kl * sbar) / (DONE + kl * conc)
    else
      kd = DZERO
    end if
  end function get_langmuir_kd

  !> @ brief Calculate zero-order decay rate and constrain if necessary
  !!
  !!  Function to calculate the zero-order decay rate from the user specified
  !!  decay rate.  If the decay rate is positive, then the decay rate must
  !!  be constrained so that more mass is not removed than is available.
  !!  Without this constraint, negative concentrations could result from
  !!  zero-order decay.
  !<
  function get_zero_order_decay(decay_rate_usr, decay_rate_last, kiter, &
                                cold, cnew, delt) result(decay_rate)
    ! -- dummy
    real(DP), intent(in) :: decay_rate_usr !< user-entered decay rate
    real(DP), intent(in) :: decay_rate_last !< decay rate used for last iteration
    integer(I4B), intent(in) :: kiter !< Picard iteration counter
    real(DP), intent(in) :: cold !< concentration at end of last time step
    real(DP), intent(in) :: cnew !< concentration at end of this time step
    real(DP), intent(in) :: delt !< length of time step
    ! -- return
    real(DP) :: decay_rate !< returned value for decay rate
    !
    ! -- Return user rate if production, otherwise constrain, if necessary
    if (decay_rate_usr < DZERO) then
      !
      ! -- Production, no need to limit rate
      decay_rate = decay_rate_usr
    else
      !
      ! -- Need to ensure decay does not result in negative
      !    concentration, so reduce the rate if it would result in
      !    removing more mass than is in the cell.
      if (kiter == 1) then
        decay_rate = min(decay_rate_usr, cold / delt)
      else
        decay_rate = decay_rate_last
        if (cnew < DZERO) then
          decay_rate = decay_rate_last + cnew / delt
        else if (cnew > cold) then
          decay_rate = decay_rate_last + cnew / delt
        end if
        decay_rate = min(decay_rate_usr, decay_rate)
      end if
      decay_rate = max(decay_rate, DZERO)
    end if
  end function get_zero_order_decay

end module GwtMstModule
