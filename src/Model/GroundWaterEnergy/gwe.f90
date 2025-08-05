! Groundwater Energy Transport (GWE) Model

module GweModule

  use KindModule, only: DP, I4B
  use InputOutputModule, only: ParseLine, upcase
  use ConstantsModule, only: LENFTYPE, LENMEMPATH, DZERO, DNODATA, &
                             LENPAKLOC, LENVARNAME, LENPACKAGETYPE, &
                             LINELENGTH
  use NumericalModelModule, only: NumericalModelType
  use BaseModelModule, only: BaseModelType
  use BndModule, only: BndType, AddBndToList, GetBndFromList
  use GweCndModule, only: GweCndType
  use GweEstModule, only: GweEstType
  use BudgetModule, only: BudgetType
  use GweInputDataModule, only: GweInputDataType
  use TransportModelModule
  use MatrixBaseModule

  implicit none

  private
  public :: gwe_cr
  public :: GweModelType
  public :: CastAsGweModel

  public :: GWE_NBASEPKG, GWE_NMULTIPKG
  public :: GWE_BASEPKG, GWE_MULTIPKG
  character(len=LENVARNAME), parameter :: dvt = 'TEMPERATURE     ' !< dependent variable type, varies based on model type
  character(len=LENVARNAME), parameter :: dvu = 'ENERGY          ' !< dependent variable unit of measure, either "mass" or "energy"
  character(len=LENVARNAME), parameter :: dvua = 'E               ' !< abbreviation of the dependent variable unit of measure, either "M" or "E"

  type, extends(TransportModelType) :: GweModelType

    type(GweInputDataType), pointer :: gwecommon => null() !< container for data shared with multiple packages
    type(GweEstType), pointer :: est => null() !< mass storage and transfer package
    type(GweCndType), pointer :: cnd => null() !< dispersion package
    integer(I4B), pointer :: inest => null() ! unit number EST
    integer(I4B), pointer :: incnd => null() ! unit number CND

  contains

    procedure :: model_df => gwe_df
    procedure :: model_ac => gwe_ac
    procedure :: model_mc => gwe_mc
    procedure :: model_ar => gwe_ar
    procedure :: model_rp => gwe_rp
    procedure :: model_dt => gwe_dt
    procedure :: model_ad => gwe_ad
    procedure :: model_cf => gwe_cf
    procedure :: model_fc => gwe_fc
    procedure :: model_cc => gwe_cc
    procedure :: model_cq => gwe_cq
    procedure :: model_bd => gwe_bd
    procedure :: tsp_ot_flow => gwe_ot_flow
    procedure :: model_da => gwe_da
    procedure :: model_bdentry => gwe_bdentry
    procedure :: allocate_scalars
    procedure :: get_iasym => gwe_get_iasym
    procedure :: create_packages => create_gwe_packages
    procedure, private :: create_bndpkgs
    procedure, private :: package_create

  end type GweModelType

  !> @brief GWE base package array descriptors
  !!
  !! GWE6 model base package types.  Only listed packages are candidates
  !! for input and these will be loaded in the order specified.
  !<
  integer(I4B), parameter :: GWE_NBASEPKG = 50
  character(len=LENPACKAGETYPE), dimension(GWE_NBASEPKG) :: GWE_BASEPKG
  data GWE_BASEPKG/'DIS6 ', 'DISV6', 'DISU6', '     ', '     ', & !  5
                  &'IC6  ', 'FMI6 ', 'EST6 ', 'ADV6 ', '     ', & ! 10
                  &'CND6 ', 'SSM6 ', 'MVE6 ', 'OC6  ', '     ', & ! 15
                  &'OBS6 ', '     ', '     ', '     ', '     ', & ! 20
                  &30*'     '/ ! 50

  !> @brief GWE multi package array descriptors
  !!
  !! GWE6 model multi-instance package types.  Only listed packages are
  !! candidates for input and these will be loaded in the order specified.
  !<
  integer(I4B), parameter :: GWE_NMULTIPKG = 50
  character(len=LENPACKAGETYPE), dimension(GWE_NMULTIPKG) :: GWE_MULTIPKG
  data GWE_MULTIPKG/'CTP6 ', 'ESL6 ', 'LKE6 ', 'SFE6 ', '     ', & !  5
                   &'MWE6 ', 'UZE6 ', 'API6 ', '     ', '     ', & ! 10
                   &40*'     '/ ! 50

  ! -- Size of supported model package arrays
  integer(I4B), parameter :: NIUNIT_GWE = GWE_NBASEPKG + GWE_NMULTIPKG

contains

  !> @brief Create a new groundwater energy transport model object
  !<
  subroutine gwe_cr(filename, id, modelname)
    ! -- modules
    use ListsModule, only: basemodellist
    use BaseModelModule, only: AddBaseModelToList
    use ConstantsModule, only: LENPACKAGENAME
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerExtModule, only: mem_set_value
    use GweNamInputModule, only: GweNamParamFoundType
    use BudgetModule, only: budget_cr
    use GweInputDataModule, only: gweshared_dat_cr
    ! -- dummy
    character(len=*), intent(in) :: filename !< input file
    integer(I4B), intent(in) :: id !< consecutive model number listed in mfsim.nam
    character(len=*), intent(in) :: modelname !< name of the model
    ! -- local
    integer(I4B) :: indis
    type(GweModelType), pointer :: this
    class(BaseModelType), pointer :: model
    !
    ! -- Allocate a new GWE Model (this) and add it to basemodellist
    allocate (this)
    !
    ! -- Set memory path before allocation in memory manager can be done
    this%memoryPath = create_mem_path(modelname)
    !
    ! -- Allocate scalars and add model to basemodellist
    call this%allocate_scalars(modelname)
    !
    ! -- Set labels for transport model - needed by create_packages() below
    call this%set_tsp_labels(this%macronym, dvt, dvu, dvua)
    !
    model => this
    call AddBaseModelToList(basemodellist, model)
    !
    ! -- Instantiate shared data container
    call gweshared_dat_cr(this%gwecommon)
    !
    ! -- Call parent class routine
    call this%tsp_cr(filename, id, modelname, 'GWE', indis)
    !
    ! -- Create model packages
    call this%create_packages(indis)
  end subroutine gwe_cr

  !> @brief Define packages of the GWE model
  !!
  !! This subroutine defines a gwe model type. Steps include:
  !!   - call df routines for each package
  !!   - set variables and pointers
  !<
  subroutine gwe_df(this)
    ! -- modules
    use SimModule, only: store_error
    ! -- dummy
    class(GweModelType) :: this
    ! -- local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    !
    ! -- Define packages and utility objects
    call this%dis%dis_df()
    call this%fmi%fmi_df(this%dis, 0)
    if (this%inmvt > 0) call this%mvt%mvt_df(this%dis)
    if (this%inadv > 0) call this%adv%adv_df()
    if (this%incnd > 0) call this%cnd%cnd_df(this%dis)
    if (this%inssm > 0) call this%ssm%ssm_df()
    call this%oc%oc_df()
    call this%budget%budget_df(NIUNIT_GWE, this%depvarunit, &
                               this%depvarunitabbrev)
    !
    ! -- Check for SSM package
    if (this%inssm == 0) then
      if (this%fmi%nflowpack > 0) then
        call store_error('Flow model has boundary packages, but there &
          &is no SSM package.  The SSM package must be activated.', &
          terminate=.TRUE.)
      end if
    end if
    !
    ! -- Assign or point model members to dis members
    this%neq = this%dis%nodes
    this%nja = this%dis%nja
    this%ia => this%dis%con%ia
    this%ja => this%dis%con%ja
    !
    ! -- Allocate model arrays, now that neq and nja are assigned
    call this%allocate_arrays()
    !
    ! -- Define packages and assign iout for time series managers
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_df(this%neq, this%dis)
      packobj%TsManager%iout = this%iout
      packobj%TasManager%iout = this%iout
    end do
    !
    ! -- Store information needed for observations
    call this%obs%obs_df(this%iout, this%name, 'GWE', this%dis)
  end subroutine gwe_df

  !> @brief Add the internal connections of this model to the sparse matrix
  !<
  subroutine gwe_ac(this, sparse)
    ! -- modules
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(GweModelType) :: this
    type(sparsematrix), intent(inout) :: sparse
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    !
    ! -- Add the internal connections of this model to sparse
    call this%dis%dis_ac(this%moffset, sparse)
    if (this%incnd > 0) &
      call this%cnd%cnd_ac(this%moffset, sparse)
    !
    ! -- Add any package connections
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ac(this%moffset, sparse)
    end do
  end subroutine gwe_ac

  !> @brief Map the positions of the GWE model connections in the numerical
  !! solution coefficient matrix.
  !<
  subroutine gwe_mc(this, matrix_sln)
    ! -- dummy
    class(GweModelType) :: this
    class(MatrixBaseType), pointer :: matrix_sln !< global system matrix
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    !
    ! -- Find the position of each connection in the global ia, ja structure
    !    and store them in idxglo.
    call this%dis%dis_mc(this%moffset, this%idxglo, matrix_sln)
    !
    if (this%incnd > 0) call this%cnd%cnd_mc(this%moffset, matrix_sln)
    !
    ! -- Map any package connections
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_mc(this%moffset, matrix_sln)
    end do
  end subroutine gwe_mc

  !> @brief GWE Model Allocate and Read
  !!
  !! This subroutine:
  !!   - allocates and reads packages that are part of this model,
  !!   - allocates memory for arrays used by this model object
  !<
  subroutine gwe_ar(this)
    ! -- modules
    use ConstantsModule, only: DHNOFLO
    ! -- dummy
    class(GweModelType) :: this
    ! -- locals
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    !
    ! -- Allocate and read modules attached to model
    call this%fmi%fmi_ar(this%ibound)
    if (this%inmvt > 0) call this%mvt%mvt_ar()
    if (this%inic > 0) call this%ic%ic_ar(this%x)
    if (this%inest > 0) call this%est%est_ar(this%dis, this%ibound)
    if (this%inadv > 0) call this%adv%adv_ar(this%dis, this%ibound)
    if (this%incnd > 0) call this%cnd%cnd_ar(this%ibound, this%est%porosity)
    if (this%inssm > 0) call this%ssm%ssm_ar(this%dis, this%ibound, this%x)
    if (this%inobs > 0) call this%obs%tsp_obs_ar(this%ic, this%x, this%flowja)
    !
    ! -- Set governing equation scale factor
    this%eqnsclfac = this%gwecommon%gwerhow * this%gwecommon%gwecpw
    !
    ! -- Call dis_ar to write binary grid file
    !call this%dis%dis_ar(this%npf%icelltype)
    !
    ! -- Set up output control
    call this%oc%oc_ar(this%x, this%dis, DHNOFLO, this%depvartype)
    call this%budget%set_ibudcsv(this%oc%ibudcsv)
    !
    ! -- Package input files now open, so allocate and read
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%set_pointers(this%dis%nodes, this%ibound, this%x, &
                                this%xold, this%flowja)
      ! -- Read and allocate package
      call packobj%bnd_ar()
    end do
  end subroutine gwe_ar

  !> @brief GWE Model Read and Prepare
  !!
  !! This subroutine calls the attached packages' read and prepare routines
  !<
  subroutine gwe_rp(this)
    ! -- modules
    use TdisModule, only: readnewdata
    ! -- dummy
    class(GweModelType) :: this
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    !
    ! -- In fmi, check for mvt and mvrbudobj consistency
    call this%fmi%fmi_rp(this%inmvt)
    if (this%inmvt > 0) call this%mvt%mvt_rp()
    !
    ! -- Check with TDIS on whether or not it is time to RP
    if (.not. readnewdata) return
    !
    ! -- Read and prepare
    if (this%inoc > 0) call this%oc%oc_rp()
    if (this%inssm > 0) call this%ssm%ssm_rp()
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_rp()
      call packobj%bnd_rp_obs()
    end do
  end subroutine gwe_rp

  !> @brief GWT Model time step size
  !!
  !! Calculate the maximum allowable time step size subject to time-step
  !! constraints.  If adaptive time steps are used, then the time step used
  !! will be no larger than dtmax calculated here.
  !<
  subroutine gwe_dt(this)
    use TdisModule, only: kstp, kper
    use AdaptiveTimeStepModule, only: ats_submit_delt
    ! dummy
    class(GweModelType) :: this
    ! local
    real(DP) :: dtmax
    character(len=LINELENGTH) :: msg
    dtmax = DNODATA

    ! advection package courant stability
    call this%adv%adv_dt(dtmax, msg, this%est%porosity)
    if (msg /= '') then
      call ats_submit_delt(kstp, kper, dtmax, msg)
    end if
  end subroutine gwe_dt

  !> @brief GWE Model Time Step Advance
  !!
  !! This subroutine calls the attached packages' advance subroutines
  !<
  subroutine gwe_ad(this)
    ! -- modules
    use SimVariablesModule, only: isimcheck, iFailedStepRetry
    ! -- dummy
    class(GweModelType) :: this
    class(BndType), pointer :: packobj
    ! -- local
    integer(I4B) :: irestore
    integer(I4B) :: ip, n
    !
    ! -- Reset state variable
    irestore = 0
    if (iFailedStepRetry > 0) irestore = 1
    if (irestore == 0) then
      !
      ! -- Copy x into xold
      do n = 1, this%dis%nodes
        if (this%ibound(n) == 0) then
          this%xold(n) = DZERO
        else
          this%xold(n) = this%x(n)
        end if
      end do
    else
      !
      ! -- Copy xold into x if this time step is a redo
      do n = 1, this%dis%nodes
        this%x(n) = this%xold(n)
      end do
    end if
    !
    ! -- Advance fmi
    call this%fmi%fmi_ad(this%x)
    !
    ! -- Advance
    if (this%incnd > 0) call this%cnd%cnd_ad()
    if (this%inssm > 0) call this%ssm%ssm_ad()
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ad()
      if (isimcheck > 0) then
        call packobj%bnd_ck()
      end if
    end do
    !
    ! -- Push simulated values to preceding time/subtime step
    call this%obs%obs_ad()
  end subroutine gwe_ad

  !> @brief GWE Model calculate coefficients
  !!
  !! This subroutine calls the attached packages' calculate coefficients
  !! subroutines
  !<
  subroutine gwe_cf(this, kiter)
    ! -- dummy
    class(GweModelType) :: this
    integer(I4B), intent(in) :: kiter
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    !
    ! -- Call package cf routines
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_cf()
    end do
  end subroutine gwe_cf

  !> @brief GWE Model fill coefficients
  !!
  !! This subroutine calls the attached packages' fill coefficients
  !! subroutines
  !<
  subroutine gwe_fc(this, kiter, matrix_sln, inwtflag)
    ! -- dummy
    class(GweModelType) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in) :: inwtflag
    ! -- local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    !
    ! -- Call fc routines
    call this%fmi%fmi_fc(this%dis%nodes, this%xold, this%nja, matrix_sln, &
                         this%idxglo, this%rhs)
    if (this%inmvt > 0) then
      call this%mvt%mvt_fc(this%x, this%x)
    end if
    if (this%inest > 0) then
      call this%est%est_fc(this%dis%nodes, this%xold, this%nja, matrix_sln, &
                           this%idxglo, this%x, this%rhs, kiter)
    end if
    if (this%inadv > 0) then
      call this%adv%adv_fc(this%dis%nodes, matrix_sln, this%idxglo, this%x, &
                           this%rhs)
    end if
    if (this%incnd > 0) then
      call this%cnd%cnd_fc(kiter, this%dis%nodes, this%nja, matrix_sln, &
                           this%idxglo, this%rhs, this%x)
    end if
    if (this%inssm > 0) then
      call this%ssm%ssm_fc(matrix_sln, this%idxglo, this%rhs)
    end if
    !
    ! -- Packages
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_fc(this%rhs, this%ia, this%idxglo, matrix_sln)
    end do
  end subroutine gwe_fc

  !> @brief GWE Model Final Convergence Check
  !!
  !! If MVR/MVT is active, this subroutine calls the MVR convergence check
  !! subroutines.
  !<
  subroutine gwe_cc(this, innertot, kiter, iend, icnvgmod, cpak, ipak, dpak)
    ! -- dummy
    class(GweModelType) :: this
    integer(I4B), intent(in) :: innertot
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: iend
    integer(I4B), intent(in) :: icnvgmod
    character(len=LENPAKLOC), intent(inout) :: cpak
    integer(I4B), intent(inout) :: ipak
    real(DP), intent(inout) :: dpak
    !
    ! -- If mover is on, then at least 2 outers required
    if (this%inmvt > 0) call this%mvt%mvt_cc(kiter, iend, icnvgmod, cpak, dpak)
  end subroutine gwe_cc

  !> @brief GWE Model calculate flow
  !!
  !! This subroutine calls the attached packages' intercell flows (flow ja)
  !<
  subroutine gwe_cq(this, icnvg, isuppress_output)
    ! -- modules
    use SparseModule, only: csr_diagsum
    ! -- dummy
    class(GweModelType) :: this
    integer(I4B), intent(in) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    !
    ! -- Construct the flowja array.  Flowja is calculated each time, even if
    !    output is suppressed.  (flowja is positive into a cell.)  The diagonal
    !    position of the flowja array will contain the flow residual after
    !    these routines are called, so each package is responsible for adding
    !    its flow to this diagonal position.
    do i = 1, this%nja
      this%flowja(i) = DZERO
    end do
    if (this%inadv > 0) call this%adv%adv_cq(this%x, this%flowja)
    if (this%incnd > 0) call this%cnd%cnd_cq(this%x, this%flowja)
    if (this%inest > 0) call this%est%est_cq(this%dis%nodes, this%x, this%xold, &
                                             this%flowja)
    if (this%inssm > 0) call this%ssm%ssm_cq(this%flowja)
    if (this%infmi > 0) call this%fmi%fmi_cq(this%x, this%flowja)
    !
    ! -- Go through packages and call cq routines.  cf() routines are called
    !    first to regenerate non-linear terms to be consistent with the final
    !    conc solution.
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_cf()
      call packobj%bnd_cq(this%x, this%flowja)
    end do
    !
    ! -- Finalize calculation of flowja by adding face flows to the diagonal.
    !    This results in the flow residual being stored in the diagonal
    !    position for each cell.
    call csr_diagsum(this%dis%con%ia, this%flowja)
  end subroutine gwe_cq

  !> @brief GWE Model Budget
  !!
  !! This subroutine:
  !!   - calculates intercell flows (flowja)
  !!   - calculates package contributions to the model budget
  !<
  subroutine gwe_bd(this, icnvg, isuppress_output)
    ! -- dummy
    class(GweModelType) :: this
    integer(I4B), intent(in) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    ! -- local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    !
    ! -- Save the solution convergence flag
    this%icnvg = icnvg
    !
    ! -- Budget routines (start by resetting).  Sole purpose of this section
    !    is to add in and outs to model budget.  All ins and out for a model
    !    should be added here to this%budget.  In a subsequent exchange call,
    !    exchange flows might also be added.
    call this%budget%reset()
    if (this%inest > 0) call this%est%est_bd(isuppress_output, this%budget)
    if (this%inssm > 0) call this%ssm%ssm_bd(isuppress_output, this%budget)
    if (this%infmi > 0) call this%fmi%fmi_bd(isuppress_output, this%budget)
    if (this%inmvt > 0) call this%mvt%mvt_bd(this%x, this%x)
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_bd(this%budget)
    end do
  end subroutine gwe_bd

  !> @brief GWE model output routine
  !!
  !! Save and print flows
  !<
  subroutine gwe_ot_flow(this, icbcfl, ibudfl, icbcun)
    ! dummy
    class(GweModelType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    integer(I4B), intent(in) :: icbcun
    ! local

    if (this%inest > 0) call this%est%est_ot_flow(icbcfl, icbcun)
    call this%TransportModelType%tsp_ot_flow(icbcfl, ibudfl, icbcun)

  end subroutine gwe_ot_flow

  !> @brief Deallocate
  !!
  !! Deallocate memory at conclusion of model run
  !<
  subroutine gwe_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    ! -- dummy
    class(GweModelType) :: this
    ! -- local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    !
    ! -- Deallocate idm memory
    call memorystore_remove(this%name, 'NAM', idm_context)
    call memorystore_remove(component=this%name, context=idm_context)
    !
    ! -- Internal flow packages deallocate
    call this%dis%dis_da()
    call this%ic%ic_da()
    call this%fmi%fmi_da()
    call this%adv%adv_da()
    call this%cnd%cnd_da()
    call this%ssm%ssm_da()
    call this%est%est_da()
    call this%mvt%mvt_da()
    call this%budget%budget_da()
    call this%oc%oc_da()
    call this%obs%obs_da()
    call this%gwecommon%gweshared_dat_da()
    !
    ! -- Internal package objects
    deallocate (this%dis)
    deallocate (this%ic)
    deallocate (this%fmi)
    deallocate (this%adv)
    deallocate (this%cnd)
    deallocate (this%ssm)
    deallocate (this%est)
    deallocate (this%mvt)
    deallocate (this%budget)
    deallocate (this%oc)
    deallocate (this%obs)
    nullify (this%gwecommon)
    !
    ! -- Boundary packages
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_da()
      deallocate (packobj)
    end do
    !
    ! -- Scalars
    call mem_deallocate(this%inest)
    call mem_deallocate(this%incnd)
    !
    ! -- Parent class members
    call this%TransportModelType%tsp_da()
    !
    ! -- NumericalModelType
    call this%NumericalModelType%model_da()
  end subroutine gwe_da

  !> @brief GroundWater Energy Transport Model Budget Entry
  !!
  !! This subroutine adds a budget entry to the flow budget.  It was added as
  !! a method for the gwe model object so that the exchange object could add its
  !! contributions.
  !<
  subroutine gwe_bdentry(this, budterm, budtxt, rowlabel)
    ! -- modules
    use ConstantsModule, only: LENBUDTXT
    use TdisModule, only: delt
    ! -- dummy
    class(GweModelType) :: this
    real(DP), dimension(:, :), intent(in) :: budterm
    character(len=LENBUDTXT), dimension(:), intent(in) :: budtxt
    character(len=*), intent(in) :: rowlabel
    !
    call this%budget%addentry(budterm, delt, budtxt, rowlabel=rowlabel)
  end subroutine gwe_bdentry

  !> @brief return 1 if any package causes the matrix to be asymmetric.
  !! Otherwise return 0.
  !<
  function gwe_get_iasym(this) result(iasym)
    class(GweModelType) :: this
    ! -- local
    integer(I4B) :: iasym
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    !
    ! -- Start by setting iasym to zero
    iasym = 0
    !
    ! -- ADV
    if (this%inadv > 0) then
      if (this%adv%iasym /= 0) iasym = 1
    end if
    !
    ! -- CND
    if (this%incnd > 0) then
      if (this%cnd%ixt3d /= 0) iasym = 1
    end if
    !
    ! -- Check for any packages that introduce matrix asymmetry
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      if (packobj%iasym /= 0) iasym = 1
    end do
  end function gwe_get_iasym

  !> Allocate memory for non-allocatable members
  !!
  !! A subroutine for allocating the scalars specific to the GWE model type.
  !! Additional scalars used by the parent class are allocated by the parent
  !! class.
  !<
  subroutine allocate_scalars(this, modelname)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GweModelType) :: this
    character(len=*), intent(in) :: modelname
    !
    ! -- Allocate parent class scalars
    call this%allocate_tsp_scalars(modelname)
    !
    ! -- Allocate members that are part of model class
    call mem_allocate(this%inest, 'INEST', this%memoryPath)
    call mem_allocate(this%incnd, 'INCND', this%memoryPath)
    !
    this%inest = 0
    this%incnd = 0
  end subroutine allocate_scalars

  !> @brief Create boundary condition packages for this model
  !!
  !! This subroutine calls the package create routines for packages activated
  !! by the user.
  !<
  subroutine package_create(this, filtyp, ipakid, ipaknum, pakname, mempath, &
                            inunit, iout)
    ! -- modules
    use SimModule, only: store_error
    use GweCtpModule, only: ctp_create
    use GweEslModule, only: esl_create
    use GweLkeModule, only: lke_create
    use GweSfeModule, only: sfe_create
    use GweMweModule, only: mwe_create
    use GweUzeModule, only: uze_create
    use ApiModule, only: api_create
    ! -- dummy
    class(GweModelType) :: this
    character(len=*), intent(in) :: filtyp
    character(len=LINELENGTH) :: errmsg
    integer(I4B), intent(in) :: ipakid
    integer(I4B), intent(in) :: ipaknum
    character(len=*), intent(in) :: pakname
    character(len=*), intent(in) :: mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    ! -- local
    class(BndType), pointer :: packobj
    class(BndType), pointer :: packobj2
    integer(I4B) :: ip
    !
    ! -- This part creates the package object
    select case (filtyp)
    case ('CTP6')
      call ctp_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, this%depvartype, mempath)
    case ('ESL6')
      call esl_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, this%gwecommon)
    case ('LKE6')
      call lke_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, this%fmi, this%eqnsclfac, this%gwecommon, &
                      this%depvartype, this%depvarunit, this%depvarunitabbrev)
    case ('SFE6')
      call sfe_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, this%fmi, this%eqnsclfac, this%gwecommon, &
                      this%depvartype, this%depvarunit, this%depvarunitabbrev)
    case ('MWE6')
      call mwe_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, this%fmi, this%eqnsclfac, this%gwecommon, &
                      this%depvartype, this%depvarunit, this%depvarunitabbrev)
    case ('UZE6')
      call uze_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, this%fmi, this%eqnsclfac, this%gwecommon, &
                      this%depvartype, this%depvarunit, this%depvarunitabbrev)
      !case('API6')
      !  call api_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
      !                  pakname)
    case default
      write (errmsg, *) 'Invalid package type: ', filtyp
      call store_error(errmsg, terminate=.TRUE.)
    end select
    !
    ! -- Packages is the bndlist that is associated with the parent model
    ! -- The following statement puts a pointer to this package in the ipakid
    ! -- position of packages.
    do ip = 1, this%bndlist%Count()
      packobj2 => GetBndFromList(this%bndlist, ip)
      if (packobj2%packName == pakname) then
        write (errmsg, '(a,a)') 'Cannot create package.  Package name  '// &
          'already exists: ', trim(pakname)
        call store_error(errmsg, terminate=.TRUE.)
      end if
    end do
    call AddBndToList(this%bndlist, packobj)
  end subroutine package_create

  !> @brief Cast to GweModelType
  !<
  function CastAsGweModel(model) result(gwemodel)
    ! -- dummy
    class(*), pointer :: model !< The object to be cast
    ! -- return
    class(GweModelType), pointer :: gwemodel !< The GWE model
    !
    gwemodel => null()
    if (.not. associated(model)) return
    select type (model)
    type is (GweModelType)
      gwemodel => model
    end select
  end function CastAsGweModel

  !> @brief Source package info and begin to process
  !<
  subroutine create_bndpkgs(this, bndpkgs, pkgtypes, pkgnames, &
                            mempaths, inunits)
    ! -- modules
    use ConstantsModule, only: LENPACKAGENAME
    use CharacterStringModule, only: CharacterStringType
    ! -- dummy
    class(GweModelType) :: this
    integer(I4B), dimension(:), allocatable, intent(inout) :: bndpkgs
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(inout) :: pkgtypes
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(inout) :: pkgnames
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(inout) :: mempaths
    integer(I4B), dimension(:), contiguous, &
      pointer, intent(inout) :: inunits
    ! -- local
    integer(I4B) :: ipakid, ipaknum
    character(len=LENFTYPE) :: pkgtype, bndptype
    character(len=LENPACKAGENAME) :: pkgname
    character(len=LENMEMPATH) :: mempath
    integer(I4B), pointer :: inunit
    integer(I4B) :: n
    !
    if (allocated(bndpkgs)) then
      !
      ! -- Create stress packages
      ipakid = 1
      bndptype = ''
      do n = 1, size(bndpkgs)
        !
        pkgtype = pkgtypes(bndpkgs(n))
        pkgname = pkgnames(bndpkgs(n))
        mempath = mempaths(bndpkgs(n))
        inunit => inunits(bndpkgs(n))
        !
        if (bndptype /= pkgtype) then
          ipaknum = 1
          bndptype = pkgtype
        end if
        !
        call this%package_create(pkgtype, ipakid, ipaknum, pkgname, mempath, &
                                 inunit, this%iout)
        ipakid = ipakid + 1
        ipaknum = ipaknum + 1
      end do
      !
      ! -- Cleanup
      deallocate (bndpkgs)
    end if
  end subroutine create_bndpkgs

  !> @brief Source package info and begin to process
  !<
  subroutine create_gwe_packages(this, indis)
    ! -- modules
    use ConstantsModule, only: LENPACKAGENAME
    use CharacterStringModule, only: CharacterStringType
    use ArrayHandlersModule, only: expandarray
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path
    use SimVariablesModule, only: idm_context
    use GweEstModule, only: est_cr
    use GweCndModule, only: cnd_cr
    ! -- dummy
    class(GweModelType) :: this
    integer(I4B), intent(in) :: indis
    ! -- local
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgtypes => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgnames => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mempaths => null()
    integer(I4B), dimension(:), contiguous, &
      pointer :: inunits => null()
    character(len=LENMEMPATH) :: model_mempath
    character(len=LENFTYPE) :: pkgtype
    character(len=LENPACKAGENAME) :: pkgname
    character(len=LENMEMPATH) :: mempath
    integer(I4B), pointer :: inunit
    integer(I4B), dimension(:), allocatable :: bndpkgs
    integer(I4B) :: n
    character(len=LENMEMPATH) :: mempathcnd = ''
    !
    ! -- Set input memory paths, input/model and input/model/namfile
    model_mempath = create_mem_path(component=this%name, context=idm_context)
    !
    ! -- Set pointers to model path package info
    call mem_setptr(pkgtypes, 'PKGTYPES', model_mempath)
    call mem_setptr(pkgnames, 'PKGNAMES', model_mempath)
    call mem_setptr(mempaths, 'MEMPATHS', model_mempath)
    call mem_setptr(inunits, 'INUNITS', model_mempath)
    !
    do n = 1, size(pkgtypes)
      !
      ! -- Attributes for this input package
      pkgtype = pkgtypes(n)
      pkgname = pkgnames(n)
      mempath = mempaths(n)
      inunit => inunits(n)
      !
      ! -- Create dis package as it is a prerequisite for other packages
      select case (pkgtype)
      case ('EST6')
        this%inest = inunit
      case ('CND6')
        this%incnd = 1
        mempathcnd = mempath
      case ('CTP6', 'ESL6', 'LKE6', 'SFE6', &
            'MWE6', 'UZE6', 'API6')
        call expandarray(bndpkgs)
        bndpkgs(size(bndpkgs)) = n
      case default
        ! TODO
      end select
    end do
    !
    ! -- Create packages that are tied directly to model
    call est_cr(this%est, this%name, this%inest, this%iout, this%fmi, &
                this%eqnsclfac, this%gwecommon)
    call cnd_cr(this%cnd, this%name, mempathcnd, this%incnd, this%iout, &
                this%fmi, this%eqnsclfac, this%gwecommon)
    !
    ! -- Check to make sure that required ftype's have been specified
    call this%ftype_check(indis, this%inest)
    !
    call this%create_bndpkgs(bndpkgs, pkgtypes, pkgnames, mempaths, inunits)
  end subroutine create_gwe_packages

end module GweModule
