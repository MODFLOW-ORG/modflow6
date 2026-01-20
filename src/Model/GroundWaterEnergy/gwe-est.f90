!>  @ brief Energy Storage and Transfer (EST) Module
!!
!!  The GweEstModule contains the GweEstType, which is related
!!  to GwtEstModule; however, there are some important differences
!!  owing to the fact that a sorbed phase is not considered.
!!  Instead, a single temperature is simulated for each grid
!!  cell and is representative of both the aqueous and solid
!!  phases (i.e., instantaneous thermal equilibrium is
!!  assumed).  Also, "thermal bleeding" is accommodated, where
!!  conductive processes can transport into, through, or
!!  out of dry cells that are part of the active domain.
!<
module GweEstModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, IZERO, DZERO, DTWO, DHALF, LENBUDTXT, DEP3
  use SimVariablesModule, only: errmsg, warnmsg
  use SimModule, only: store_error, count_errors, &
                       store_warning, store_error_filename
  use MatrixBaseModule
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use GweInputDataModule, only: GweInputDataType

  implicit none
  public :: GweEstType
  public :: est_cr
  !
  integer(I4B), parameter :: NBDITEMS = 3
  character(len=LENBUDTXT), dimension(NBDITEMS) :: budtxt
  data budtxt/' STORAGE-CELLBLK', '   DECAY-AQUEOUS', '     DECAY-SOLID'/

  !> @brief Enumerator that defines the decay options
  !<
  ENUM, BIND(C)
    ENUMERATOR :: DECAY_OFF = 0 !< Decay (or production) of thermal energy inactive (default)
    ENUMERATOR :: DECAY_ZERO_ORDER = 2 !< Zeroth-order decay
    ENUMERATOR :: DECAY_WATER = 1 !< Zeroth-order decay in water only
    ENUMERATOR :: DECAY_SOLID = 2 !< Zeroth-order decay in solid only
    ENUMERATOR :: DECAY_BOTH = 3 !< Zeroth-order decay in water and solid
  END ENUM

  !> @ brief Energy storage and transfer
  !!
  !!  Data and methods for handling changes in temperature
  !<
  type, extends(NumericalPackageType) :: GweEstType
    !
    ! -- storage
    real(DP), pointer :: cpw => null() !< heat capacity of water
    real(DP), pointer :: rhow => null() !< density of water
    real(DP), pointer :: latheatvap => null() !< latent heat of vaporization
    real(DP), dimension(:), pointer, contiguous :: cps => null() !< heat capacity of solid
    real(DP), dimension(:), pointer, contiguous :: rhos => null() !< density of solid
    real(DP), dimension(:), pointer, contiguous :: porosity => null() !< porosity
    real(DP), dimension(:), pointer, contiguous :: ratesto => null() !< rate of energy storage
    !
    ! -- decay
    integer(I4B), pointer :: idcy => null() !< order of decay rate (0:none, 1:first, 2:zero (aqueous and/or solid))
    integer(I4B), pointer :: idcysrc => null() !< decay source (or sink) (1: aqueous only, 2: solid only, 3: both phases
    real(DP), dimension(:), pointer, contiguous :: decay_water => null() !< first or zero order decay rate (aqueous)
    real(DP), dimension(:), pointer, contiguous :: decay_solid => null() !< first or zero order decay rate (solid)
    real(DP), dimension(:), pointer, contiguous :: ratedcyw => null() !< rate of decay in aqueous phase
    real(DP), dimension(:), pointer, contiguous :: ratedcys => null() !< rate of decay in solid phase
    real(DP), dimension(:), pointer, contiguous :: decaylastw => null() !< aqueous phase decay rate used for last iteration (needed for zero order decay)
    real(DP), dimension(:), pointer, contiguous :: decaylasts => null() !< solid phase decay rate used for last iteration (needed for zero order decay)
    !
    ! -- misc
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    type(TspFmiType), pointer :: fmi => null() !< pointer to fmi object
    type(GweInputDataType), pointer :: gwecommon => null() !< pointer to shared gwe data used by multiple packages but set in est
    real(DP), pointer :: eqnsclfac => null() !< governing equation scale factor; =rhow*cpw for energy

  contains

    procedure :: est_ar
    procedure :: est_fc
    procedure :: est_fc_sto
    procedure :: est_fc_dcy_water
    procedure :: est_fc_dcy_solid
    procedure :: est_cq
    procedure :: est_cq_sto
    procedure :: est_cq_dcy
    procedure :: est_cq_dcy_solid
    procedure :: est_bd
    procedure :: est_ot_flow
    procedure :: est_da
    procedure :: allocate_scalars
    procedure, private :: allocate_arrays
    procedure, private :: source_options
    procedure, private :: source_data
    procedure, private :: log_options

  end type GweEstType

contains

  !> @ brief Create a new EST package object
  !!
  !!  Create a new EST package
  !<
  subroutine est_cr(estobj, name_model, input_mempath, inunit, iout, fmi, &
                    eqnsclfac, gwecommon)
    ! -- dummy
    type(GweEstType), pointer :: estobj !< unallocated new est object to create
    character(len=*), intent(in) :: name_model !< name of the model
    character(len=*), intent(in) :: input_mempath !< input mempath of package
    integer(I4B), intent(in) :: inunit !< unit number of WEL package input file
    integer(I4B), intent(in) :: iout !< unit number of model listing file
    type(TspFmiType), intent(in), target :: fmi !< fmi package for this GWE model
    real(DP), intent(in), pointer :: eqnsclfac !< governing equation scale factor
    type(GweInputDataType), intent(in), target :: gwecommon !< shared data container for use by multiple GWE packages
    !
    ! -- Create the object
    allocate (estobj)
    !
    ! -- create name and memory path
    call estobj%set_names(1, name_model, 'EST', 'EST', input_mempath)
    !
    ! -- Allocate scalars
    call estobj%allocate_scalars()
    !
    ! -- Set variables
    estobj%inunit = inunit
    estobj%iout = iout
    estobj%fmi => fmi
    estobj%eqnsclfac => eqnsclfac
    estobj%gwecommon => gwecommon
  end subroutine est_cr

  !> @ brief Allocate and read method for package
  !!
  !!  Method to allocate and read static data for the package.
  !<
  subroutine est_ar(this, dis, ibound)
    ! -- modules
    use GweInputDataModule, only: set_gwe_dat_ptrs
    ! -- dummy
    class(GweEstType), intent(inout) :: this !< GweEstType object
    class(DisBaseType), pointer, intent(in) :: dis !< pointer to dis package
    integer(I4B), dimension(:), pointer, contiguous :: ibound !< pointer to GWE ibound array
    ! -- formats
    character(len=*), parameter :: fmtest = &
      "(1x,/1x,'EST -- ENERGY STORAGE AND TRANSFER PACKAGE, VERSION 1, &
      &7/29/2020 INPUT READ FROM MEMPATH: ', a, /)"
    !
    ! --print a message identifying the energy storage and transfer package.
    write (this%iout, fmtest) this%input_mempath
    !
    ! -- Read options
    call this%source_options()
    !
    ! -- store pointers to arguments that were passed in
    this%dis => dis
    this%ibound => ibound
    !
    ! -- Allocate arrays
    call this%allocate_arrays(dis%nodes)
    !
    ! -- read the gridded data
    call this%source_data()
    !
    ! -- set data required by other packages
    call this%gwecommon%set_gwe_dat_ptrs(this%rhow, this%cpw, this%latheatvap, &
                                         this%rhos, this%cps)
  end subroutine est_ar

  !> @ brief Fill coefficient method for package
  !!
  !!  Method to calculate and fill coefficients for the package.
  !<
  subroutine est_fc(this, nodes, cold, nja, matrix_sln, idxglo, cnew, &
                    rhs, kiter)
    ! -- dummy
    class(GweEstType) :: this !< GweEstType object
    integer, intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cold !< temperature at end of last time step
    integer(I4B), intent(in) :: nja !< number of GWE connections
    class(MatrixBaseType), pointer :: matrix_sln !< solution matrix
    integer(I4B), intent(in), dimension(nja) :: idxglo !< mapping vector for model (local) to solution (global)
    real(DP), intent(inout), dimension(nodes) :: rhs !< right-hand side vector for model
    real(DP), intent(in), dimension(nodes) :: cnew !< temperature at end of this time step
    integer(I4B), intent(in) :: kiter !< solution outer iteration number
    !
    ! -- storage contribution
    call this%est_fc_sto(nodes, cold, nja, matrix_sln, idxglo, rhs)
    !
    ! -- decay contribution
    if (this%idcy == DECAY_ZERO_ORDER) then
      call this%est_fc_dcy_water(nodes, cold, cnew, nja, matrix_sln, idxglo, &
                                 rhs, kiter)
      call this%est_fc_dcy_solid(nodes, cold, nja, matrix_sln, idxglo, rhs, &
                                 cnew, kiter)
    end if
  end subroutine est_fc

  !> @ brief Fill storage coefficient method for package
  !!
  !!  Method to calculate and fill storage coefficients for the package.
  !<
  subroutine est_fc_sto(this, nodes, cold, nja, matrix_sln, idxglo, rhs)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(GweEstType) :: this !< GweEstType object
    integer, intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cold !< temperature at end of last time step
    integer(I4B), intent(in) :: nja !< number of GWE connections
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    integer(I4B), intent(in), dimension(nja) :: idxglo !< mapping vector for model (local) to solution (global)
    real(DP), intent(inout), dimension(nodes) :: rhs !< right-hand side vector for model
    ! -- local
    integer(I4B) :: n, idiag
    real(DP) :: tled
    real(DP) :: hhcof, rrhs
    real(DP) :: vnew, vold, vcell, vsolid, term
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
      ! -- calculate new and old water volumes and solid volume
      vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
      vnew = vcell * this%fmi%gwfsat(n) * this%porosity(n)
      vold = vnew
      if (this%fmi%igwfstrgss /= 0) vold = vold + this%fmi%gwfstrgss(n) * delt
      if (this%fmi%igwfstrgsy /= 0) vold = vold + this%fmi%gwfstrgsy(n) * delt
      vsolid = vcell * (DONE - this%porosity(n))
      !
      ! -- add terms to diagonal and rhs accumulators
      term = (this%rhos(n) * this%cps(n)) * vsolid
      hhcof = -(this%eqnsclfac * vnew + term) * tled
      rrhs = -(this%eqnsclfac * vold + term) * tled * cold(n)
      idiag = this%dis%con%ia(n)
      call matrix_sln%add_value_pos(idxglo(idiag), hhcof)
      rhs(n) = rhs(n) + rrhs
    end do
  end subroutine est_fc_sto

  !> @ brief Fill decay coefficient method for package
  !!
  !!  Method to calculate and fill decay coefficients for the package.
  !<
  subroutine est_fc_dcy_water(this, nodes, cold, cnew, nja, matrix_sln, &
                              idxglo, rhs, kiter)
    ! -- dummy
    class(GweEstType) :: this !< GweEstType object
    integer, intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cold !< temperature at end of last time step
    real(DP), intent(in), dimension(nodes) :: cnew !< temperature at end of this time step
    integer(I4B), intent(in) :: nja !< number of GWE connections
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    integer(I4B), intent(in), dimension(nja) :: idxglo !< mapping vector for model (local) to solution (global)
    real(DP), intent(inout), dimension(nodes) :: rhs !< right-hand side vector for model
    integer(I4B), intent(in) :: kiter !< solution outer iteration number
    ! -- local
    integer(I4B) :: n
    real(DP) :: rrhs
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
      ! -- add zero-order decay rate terms to accumulators
      if (this%idcy == DECAY_ZERO_ORDER .and. (this%idcysrc == DECAY_WATER .or. &
                                               this%idcysrc == DECAY_BOTH)) then
        !
        decay_rate = this%decay_water(n)
        ! -- This term does get divided by eqnsclfac for fc purposes because it
        !    should start out being a rate of energy
        this%decaylastw(n) = decay_rate
        rrhs = decay_rate * vcell * swtpdt * this%porosity(n)
        rhs(n) = rhs(n) + rrhs
      end if
      !
    end do
  end subroutine est_fc_dcy_water

  !> @ brief Fill solid decay coefficient method for package
  !!
  !!  Method to calculate and fill energy decay coefficients for the solid phase.
  !<
  subroutine est_fc_dcy_solid(this, nodes, cold, nja, matrix_sln, idxglo, &
                              rhs, cnew, kiter)
    ! -- dummy
    class(GweEstType) :: this !< GwtMstType object
    integer, intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cold !< temperature at end of last time step
    integer(I4B), intent(in) :: nja !< number of GWE connections
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    integer(I4B), intent(in), dimension(nja) :: idxglo !< mapping vector for model (local) to solution (global)
    real(DP), intent(inout), dimension(nodes) :: rhs !< right-hand side vector for model
    real(DP), intent(in), dimension(nodes) :: cnew !< temperature at end of this time step
    integer(I4B), intent(in) :: kiter !< solution outer iteration number
    ! -- local
    integer(I4B) :: n
    real(DP) :: rrhs
    real(DP) :: vcell
    real(DP) :: decay_rate
    !
    ! -- loop through and calculate sorption contribution to hcof and rhs
    do n = 1, this%dis%nodes
      !
      ! -- skip if transport inactive
      if (this%ibound(n) <= 0) cycle
      !
      ! -- set variables
      rrhs = DZERO
      vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
      !
      ! -- account for zero-order decay rate terms in rhs
      if (this%idcy == DECAY_ZERO_ORDER .and. (this%idcysrc == DECAY_SOLID .or. &
                                               this%idcysrc == DECAY_BOTH)) then
        !
        ! -- negative temps are currently not checked for or prevented since a
        !    user can define a temperature scale of their own choosing.  if
        !    negative temps result from the specified zero-order decay value,
        !    it is up to the user to decide if the calculated temperatures are
        !    acceptable
        decay_rate = this%decay_solid(n)
        this%decaylasts(n) = decay_rate
        rrhs = decay_rate * vcell * (1 - this%porosity(n)) * this%rhos(n)
        rhs(n) = rhs(n) + rrhs
      end if
    end do
  end subroutine est_fc_dcy_solid

  !> @ brief Calculate flows for package
  !!
  !!  Method to calculate flows for the package.
  !<
  subroutine est_cq(this, nodes, cnew, cold, flowja)
    ! -- dummy
    class(GweEstType) :: this !< GweEstType object
    integer(I4B), intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cnew !< temperature at end of this time step
    real(DP), intent(in), dimension(nodes) :: cold !< temperature at end of last time step
    real(DP), dimension(:), contiguous, intent(inout) :: flowja !< flow between two connected control volumes
    ! -- local
    !
    ! - storage
    call this%est_cq_sto(nodes, cnew, cold, flowja)
    !
    ! -- decay
    if (this%idcy == DECAY_ZERO_ORDER) then
      if (this%idcysrc == DECAY_WATER .or. this%idcysrc == DECAY_BOTH) then
        call this%est_cq_dcy(nodes, cnew, cold, flowja)
      end if
      if (this%idcysrc == DECAY_SOLID .or. this%idcysrc == DECAY_BOTH) then
        call this%est_cq_dcy_solid(nodes, cnew, cold, flowja)
      end if
    end if
  end subroutine est_cq

  !> @ brief Calculate storage terms for package
  !!
  !!  Method to calculate storage terms for the package.
  !<
  subroutine est_cq_sto(this, nodes, cnew, cold, flowja)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(GweEstType) :: this !< GweEstType object
    integer(I4B), intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cnew !< temperature at end of this time step
    real(DP), intent(in), dimension(nodes) :: cold !< temperature at end of last time step
    real(DP), dimension(:), contiguous, intent(inout) :: flowja !< flow between two connected control volumes
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: idiag
    real(DP) :: rate
    real(DP) :: tled
    real(DP) :: vwatnew, vwatold, vcell, vsolid, term
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
      ! -- calculate new and old water volumes and solid volume
      vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
      vwatnew = vcell * this%fmi%gwfsat(n) * this%porosity(n)
      vwatold = vwatnew
      if (this%fmi%igwfstrgss /= 0) vwatold = vwatold + this%fmi%gwfstrgss(n) &
                                              * delt
      if (this%fmi%igwfstrgsy /= 0) vwatold = vwatold + this%fmi%gwfstrgsy(n) &
                                              * delt
      vsolid = vcell * (DONE - this%porosity(n))
      !
      ! -- calculate rate
      term = (this%rhos(n) * this%cps(n)) * vsolid
      hhcof = -(this%eqnsclfac * vwatnew + term) * tled
      rrhs = -(this%eqnsclfac * vwatold + term) * tled * cold(n)
      rate = hhcof * cnew(n) - rrhs
      this%ratesto(n) = rate
      idiag = this%dis%con%ia(n)
      flowja(idiag) = flowja(idiag) + rate
    end do
  end subroutine est_cq_sto

  !> @ brief Calculate decay terms for aqueous phase
  !!
  !!  Method to calculate decay terms for the aqueous phase.
  !<
  subroutine est_cq_dcy(this, nodes, cnew, cold, flowja)
    ! -- dummy
    class(GweEstType) :: this !< GweEstType object
    integer(I4B), intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cnew !< temperature at end of this time step
    real(DP), intent(in), dimension(nodes) :: cold !< temperature at end of last time step
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
    ! -- Calculate decay change
    do n = 1, nodes
      !
      ! -- skip if transport inactive
      this%ratedcyw(n) = DZERO
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
      ! -- zero order decay aqueous phase
      if (this%idcy == DECAY_ZERO_ORDER .and. &
          (this%idcysrc == DECAY_WATER .or. this%idcysrc == DECAY_BOTH)) then
        decay_rate = this%decay_water(n)
        ! -- this term does NOT get multiplied by eqnsclfac for cq purposes
        !    because it should already be a rate of energy
        rrhs = decay_rate * vcell * swtpdt * this%porosity(n)
      end if
      rate = hhcof * cnew(n) - rrhs
      this%ratedcyw(n) = rate
      idiag = this%dis%con%ia(n)
      flowja(idiag) = flowja(idiag) + rate
      !
    end do
  end subroutine est_cq_dcy

  !> @ brief Calculate decay terms for solid phase
  !!
  !!  Method to calculate decay terms for the solid phase.
  !<
  subroutine est_cq_dcy_solid(this, nodes, cnew, cold, flowja)
    ! -- dummy
    class(GweEstType) :: this !< GweEstType object
    integer(I4B), intent(in) :: nodes !< number of nodes
    real(DP), intent(in), dimension(nodes) :: cnew !< temperature at end of this time step
    real(DP), intent(in), dimension(nodes) :: cold !< temperature at end of last time step
    real(DP), dimension(:), contiguous, intent(inout) :: flowja !< flow between two connected control volumes
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: idiag
    real(DP) :: rate
    real(DP) :: hhcof, rrhs
    real(DP) :: vcell
    real(DP) :: decay_rate
    !
    ! -- calculate decay change
    do n = 1, nodes
      !
      ! -- skip if transport inactive
      this%ratedcys(n) = DZERO
      if (this%ibound(n) <= 0) cycle
      !
      ! -- calculate new and old water volumes
      vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
      !
      ! -- calculate decay gains and losses
      rate = DZERO
      hhcof = DZERO
      rrhs = DZERO
      ! -- first-order decay (idcy=1) is not supported for temperature modeling
      if (this%idcy == DECAY_ZERO_ORDER .and. &
          (this%idcysrc == DECAY_SOLID .or. this%idcysrc == DECAY_BOTH)) then ! zero order decay in the solid phase
        decay_rate = this%decay_solid(n)
        ! -- this term does NOT get multiplied by eqnsclfac for cq purposes
        !    because it should already be a rate of energy
        rrhs = decay_rate * vcell * (1 - this%porosity(n)) * this%rhos(n)
      end if
      rate = hhcof * cnew(n) - rrhs
      this%ratedcys(n) = rate
      idiag = this%dis%con%ia(n)
      flowja(idiag) = flowja(idiag) + rate
    end do
  end subroutine est_cq_dcy_solid

  !> @ brief Calculate budget terms for package
  !!
  !!  Method to calculate budget terms for the package.
  !<
  subroutine est_bd(this, isuppress_output, model_budget)
    ! -- modules
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType, rate_accumulator
    ! -- dummy
    class(GweEstType) :: this !< GweEstType object
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
    if (this%idcy == DECAY_ZERO_ORDER) then
      if (this%idcysrc == DECAY_WATER .or. this%idcysrc == DECAY_BOTH) then
        ! -- aqueous phase
        call rate_accumulator(this%ratedcyw, rin, rout)
        call model_budget%addentry(rin, rout, delt, budtxt(2), &
                                   isuppress_output, rowlabel=this%packName)
      end if
      if (this%idcysrc == DECAY_SOLID .or. this%idcysrc == DECAY_BOTH) then
        ! -- solid phase
        call rate_accumulator(this%ratedcys, rin, rout)
        call model_budget%addentry(rin, rout, delt, budtxt(3), &
                                   isuppress_output, rowlabel=this%packName)
      end if
    end if
  end subroutine est_bd

  !> @ brief Output flow terms for package
  !!
  !!  Method to output terms for the package.
  !<
  subroutine est_ot_flow(this, icbcfl, icbcun)
    ! -- dummy
    class(GweEstType) :: this !< GweEstType object
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
      if (this%idcy == DECAY_ZERO_ORDER) then
        if (this%idcysrc == DECAY_WATER .or. this%idcysrc == DECAY_BOTH) then
          ! -- aqueous phase
          call this%dis%record_array(this%ratedcyw, this%iout, iprint, &
                                     -ibinun, budtxt(2), cdatafmp, nvaluesp, &
                                     nwidthp, editdesc, dinact)
        end if
        if (this%idcysrc == DECAY_SOLID .or. this%idcysrc == DECAY_BOTH) then
          ! -- solid phase
          call this%dis%record_array(this%ratedcys, this%iout, iprint, &
                                     -ibinun, budtxt(3), cdatafmp, nvaluesp, &
                                     nwidthp, editdesc, dinact)
        end if
      end if
    end if
  end subroutine est_ot_flow

  !> @brief Deallocate memory
  !!
  !!  Method to deallocate memory for the package.
  !<
  subroutine est_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GweEstType) :: this !< GweEstType object
    !
    ! -- Deallocate arrays if package was active
    if (this%inunit > 0) then
      call mem_deallocate(this%porosity)
      call mem_deallocate(this%ratesto)
      call mem_deallocate(this%idcy)
      call mem_deallocate(this%idcysrc)
      call mem_deallocate(this%decay_water)
      call mem_deallocate(this%decay_solid)
      call mem_deallocate(this%ratedcyw)
      call mem_deallocate(this%ratedcys)
      call mem_deallocate(this%decaylastw)
      call mem_deallocate(this%decaylasts)
      call mem_deallocate(this%cpw)
      call mem_deallocate(this%cps)
      call mem_deallocate(this%rhow)
      call mem_deallocate(this%rhos)
      call mem_deallocate(this%latheatvap)
      this%ibound => null()
      this%fmi => null()
    end if
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()
  end subroutine est_da

  !> @ brief Allocate scalar variables for package
  !!
  !!  Method to allocate scalar variables for the package.
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy
    class(GweEstType) :: this !< GweEstType object
    !
    ! -- Allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%cpw, 'CPW', this%memoryPath)
    call mem_allocate(this%rhow, 'RHOW', this%memoryPath)
    call mem_allocate(this%latheatvap, 'LATHEATVAP', this%memoryPath)
    call mem_allocate(this%idcy, 'IDCY', this%memoryPath)
    call mem_allocate(this%idcysrc, 'IDCYSRC', this%memoryPath)
    !
    ! -- Initialize
    this%cpw = DZERO
    this%rhow = DZERO
    this%latheatvap = DZERO
    this%idcy = IZERO
    this%idcysrc = IZERO
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
    class(GweEstType) :: this !< GweEstType object
    integer(I4B), intent(in) :: nodes !< number of nodes
    ! -- local
    integer(I4B) :: n
    !
    ! -- Allocate
    ! -- sto
    call mem_allocate(this%porosity, nodes, 'POROSITY', this%memoryPath)
    call mem_allocate(this%ratesto, nodes, 'RATESTO', this%memoryPath)
    call mem_allocate(this%cps, nodes, 'CPS', this%memoryPath)
    call mem_allocate(this%rhos, nodes, 'RHOS', this%memoryPath)
    !
    ! -- dcy
    if (this%idcy == DECAY_OFF) then
      call mem_allocate(this%ratedcyw, 1, 'RATEDCYW', this%memoryPath)
      call mem_allocate(this%ratedcys, 1, 'RATEDCYS', this%memoryPath)
      call mem_allocate(this%decay_water, 1, 'DECAY_WATER', this%memoryPath)
      call mem_allocate(this%decay_solid, 1, 'DECAY_SOLID', this%memoryPath)
      call mem_allocate(this%decaylastw, 1, 'DECAYLASTW', this%memoryPath)
      call mem_allocate(this%decaylasts, 1, 'DECAYLAST', this%memoryPath)
    else
      call mem_allocate(this%ratedcyw, this%dis%nodes, 'RATEDCYW', &
                        this%memoryPath)
      call mem_allocate(this%ratedcys, this%dis%nodes, 'RATEDCYS', &
                        this%memoryPath)
      call mem_allocate(this%decay_water, nodes, 'DECAY_WATER', this%memoryPath)
      call mem_allocate(this%decay_solid, nodes, 'DECAY_SOLID', this%memoryPath)
      call mem_allocate(this%decaylastw, nodes, 'DECAYLASTW', this%memoryPath)
      call mem_allocate(this%decaylasts, nodes, 'DECAYLASTS', this%memoryPath)
    end if
    !
    ! -- Initialize
    do n = 1, nodes
      this%porosity(n) = DZERO
      this%ratesto(n) = DZERO
      this%cps(n) = DZERO
      this%rhos(n) = DZERO
    end do
    do n = 1, size(this%decay_water)
      this%decay_water(n) = DZERO
      this%decay_solid(n) = DZERO
      this%ratedcyw(n) = DZERO
      this%ratedcys(n) = DZERO
      this%decaylastw(n) = DZERO
      this%decaylasts(n) = DZERO
    end do
  end subroutine allocate_arrays

  !> @brief Update simulation mempath options
  !<
  subroutine source_options(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use GweEstInputModule, only: GweEstParamFoundType
    ! -- dummy
    class(GweEstType) :: this
    ! -- locals
    type(GweEstParamFoundType) :: found
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%ipakcb, 'SAVE_FLOWS', this%input_mempath, &
                       found%save_flows)
    call mem_set_value(this%idcy, 'ORD0_DECAY_WATER', this%input_mempath, &
                       found%ord0_decay_water)
    call mem_set_value(this%idcy, 'ORD0_DECAY_SOLID', this%input_mempath, &
                       found%ord0_decay_solid)
    call mem_set_value(this%cpw, 'CPW', this%input_mempath, &
                       found%cpw)
    call mem_set_value(this%rhow, 'RHOW', this%input_mempath, &
                       found%rhow)
    call mem_set_value(this%latheatvap, 'LATHEATVAP', this%input_mempath, &
                       found%latheatvap)

    ! -- update internal state
    if (found%save_flows) this%ipakcb = -1
    if (found%ord0_decay_water .and. &
        found%ord0_decay_solid) then
      this%idcy = DECAY_ZERO_ORDER
      this%idcysrc = DECAY_BOTH
    else if (found%ord0_decay_water) then
      this%idcy = DECAY_ZERO_ORDER
      this%idcysrc = DECAY_WATER
    else if (found%ord0_decay_solid) then
      this%idcy = DECAY_ZERO_ORDER
      this%idcysrc = DECAY_SOLID
    end if
    if (found%cpw) then
      if (this%cpw <= 0.0) then
        write (errmsg, '(a)') 'Specified value for the heat capacity of &
          &water must be greater than 0.0.'
        call store_error(errmsg)
        call store_error_filename(this%input_fname)
      end if
    end if
    if (found%rhow) then
      if (this%rhow <= 0.0) then
        write (errmsg, '(a)') 'Specified value for the density of &
          &water must be greater than 0.0.'
        call store_error(errmsg)
        call store_error_filename(this%input_fname)
      end if
    end if

    ! -- log options
    call this%log_options(found)
  end subroutine source_options

  !> @brief Write user options to list file
  !<
  subroutine log_options(this, found)
    ! -- modules
    use GweEstInputModule, only: GweEstParamFoundType
    ! -- dummy
    class(GweEstType) :: this
    ! -- local
    type(GweEstParamFoundType), intent(in) :: found
    ! -- formats
    character(len=*), parameter :: fmtisvflow = &
            &"(4x,'CELL-BY-CELL FLOW INFORMATION WILL BE SAVED TO BINARY "// &
            &"FILE WHENEVER ICBCFL IS NOT ZERO.')"
    character(len=*), parameter :: fmtidcy2 = &
                                   "(4x,'ZERO-ORDER DECAY IN THE AQUEOUS "// &
                                   &"PHASE IS ACTIVE. ')"
    character(len=*), parameter :: fmtidcy3 = &
                                   "(4x,'ZERO-ORDER DECAY IN THE SOLID "// &
                                   &"PHASE IS ACTIVE. ')"

    write (this%iout, '(1x,a)') 'PROCESSING ENERGY STORAGE AND TRANSFER OPTIONS'
    if (found%save_flows) write (this%iout, fmtisvflow)
    if (found%ord0_decay_water) then
      write (this%iout, fmtidcy2)
    else if (found%ord0_decay_solid) then
      write (this%iout, fmtidcy3)
    end if
    if (found%cpw) then
      write (this%iout, '(4x,a,1pg15.6)') &
        'Heat capacity of the water has been set to: ', &
        this%cpw
    end if
    if (found%rhow) then
      write (this%iout, '(4x,a,1pg15.6)') &
        'Density of the water has been set to: ', &
        this%rhow
    end if
    if (found%latheatvap) then
      write (this%iout, '(4x,a,1pg15.6)') &
        'Latent heat of vaporization of the water has been set to: ', &
        this%latheatvap
    end if
    write (this%iout, '(1x,a,/)') &
      'END PROCESSING ENERGY STORAGE AND TRANSFER OPTIONS'
  end subroutine log_options

  !> @brief Source EST griddata from input mempath
  !<
  subroutine source_data(this)
    ! -- modules
    use SimModule, only: count_errors, store_error
    use MemoryManagerModule, only: mem_reallocate, get_isize
    use MemoryManagerExtModule, only: mem_set_value
    use ConstantsModule, only: LENMEMPATH, LINELENGTH
    use GweEstInputModule, only: GweEstParamFoundType
    ! -- dummy
    class(GweEstType) :: this
    ! -- locals
    character(len=LINELENGTH) :: errmsg
    type(GweEstParamFoundType) :: found
    integer(I4B), dimension(:), pointer, contiguous :: map
    integer(I4B) :: asize
    ! -- formats

    ! -- set map
    map => null()
    if (this%dis%nodes < this%dis%nodesuser) map => this%dis%nodeuser

    ! -- reallocate
    if (this%idcy == DECAY_OFF) then
      call get_isize('DECAY_WATER', this%input_mempath, asize)
      if (asize > 0) &
        call mem_reallocate(this%decay_water, this%dis%nodes, 'DECAY_WATER', &
                            trim(this%memoryPath))
      call get_isize('DECAY_SOLID', this%input_mempath, asize)
      if (asize > 0) &
        call mem_reallocate(this%decay_solid, this%dis%nodes, 'DECAY_SOLID', &
                            trim(this%memoryPath))
    end if

    ! -- update defaults with idm sourced values
    call mem_set_value(this%porosity, 'POROSITY', this%input_mempath, map, &
                       found%porosity)
    call mem_set_value(this%decay_water, 'DECAY_WATER', this%input_mempath, map, &
                       found%decay_water)
    call mem_set_value(this%decay_solid, 'DECAY_SOLID', this%input_mempath, map, &
                       found%decay_solid)
    call mem_set_value(this%cps, 'CPS', this%input_mempath, map, found%cps)
    call mem_set_value(this%rhos, 'RHOS', this%input_mempath, map, found%rhos)

    ! -- Check for required params
    if (.not. found%porosity) then
      write (errmsg, '(a)') 'Porosity not specified in griddata block.'
      call store_error(errmsg)
    end if
    if (.not. found%cps) then
      write (errmsg, '(a)') 'HEAT_CAPACITY_SOLID not specified in griddata block.'
      call store_error(errmsg)
    end if
    if (.not. found%rhos) then
      write (errmsg, '(a)') 'DENSITY_SOLID not specified in griddata block.'
      call store_error(errmsg)
    end if

    ! -- log griddata
    write (this%iout, '(1x,a)') 'PROCESSING ENERGY STORAGE AND TRANSFER GRIDDATA'
    if (found%porosity) &
      write (this%iout, '(4x,a)') 'POROSITY set from input file'
    if (found%decay_water) &
      write (this%iout, '(4x,a)') 'DECAY_WATER set from input file'
    if (found%decay_solid) &
      write (this%iout, '(4x,a)') 'DECAY_SOLID set from input file'
    if (found%cps) &
      write (this%iout, '(4x,a)') 'HEAT_CAPACITY_SOLID set from input file'
    if (found%rhos) &
      write (this%iout, '(4x,a)') 'DENSITY_SOLID set from input file'
    write (this%iout, '(1x,a)') &
      'END PROCESSING ENERGY STORAGE AND TRANSFER GRIDDATA'

    ! -- Check for required decay/production rate coefficients
    if (this%idcy == DECAY_ZERO_ORDER) then
      if (.not. (found%decay_water .or. found%decay_solid)) then
        write (errmsg, '(a)') 'Zero order decay in either the aqueous &
          &or solid phase is active but the corresponding zero-order &
          &rate coefficient is not specified. Either DECAY_WATER or &
          &DECAY_SOLID must be specified in the griddata block.'
        call store_error(errmsg)
      end if
    else
      if (found%decay_water) then
        write (warnmsg, '(a)') 'Zero order decay in the aqueous phase has &
          &not been activated but DECAY_WATER has been specified. Zero &
          &order decay in the aqueous phase will have no affect on &
          &simulation results.'
        call store_warning(warnmsg)
        write (this%iout, '(/1x,a)') 'WARNING: '//trim(warnmsg)
      else if (found%decay_solid) then
        write (warnmsg, '(a)') 'Zero order decay in the solid phase has not &
          &been activated but DECAY_SOLID has been specified.  Zero order &
          &decay in the solid phase will have no affect on simulation &
          &results.'
        call store_warning(warnmsg)
        write (this%iout, '(/1x,a)') 'WARNING: '//trim(warnmsg)
      end if
    end if

    ! -- terminate if errors
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine source_data

end module GweEstModule
