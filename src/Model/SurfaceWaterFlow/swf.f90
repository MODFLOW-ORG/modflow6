!> @brief Surface Water Flow (SWF) Module
!<
module SwfModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, LENFTYPE, DNODATA, &
                             LENMEMPATH, LENPACKAGETYPE
  use SimModule, only: count_errors, store_error, store_error_filename
  use SimVariablesModule, only: errmsg
  use MemoryManagerModule, only: mem_allocate
  use MemoryHelperModule, only: create_mem_path
  use NumericalModelModule, only: NumericalModelType
  use BndModule, only: BndType, AddBndToList, GetBndFromList
  use SwfIcModule, only: SwfIcType
  use SwfDfwModule, only: SwfDfwType
  use SwfCxsModule, only: SwfCxsType
  use SwfStoModule, only: SwfStoType
  use SwfObsModule, only: SwfObsType, swf_obs_cr
  use SwfOcModule, only: SwfOcType
  use BudgetModule, only: BudgetType
  use MatrixBaseModule

  implicit none

  private
  public :: SwfModelType

  type, abstract, extends(NumericalModelType) :: SwfModelType
    type(SwfIcType), pointer :: ic => null() ! initial conditions package
    type(SwfDfwType), pointer :: dfw => null() !< diffusive wave package
    type(SwfCxsType), pointer :: cxs => null() !< cross section package
    type(SwfStoType), pointer :: sto => null() !< storage package
    type(SwfObsType), pointer :: obs => null() ! observation package
    type(SwfOcType), pointer :: oc => null() !< output control package
    type(BudgetType), pointer :: budget => null() !< budget object
    integer(I4B), pointer :: inic => null() ! unit number IC
    integer(I4B), pointer :: indfw => null() !< unit number DFW
    integer(I4B), pointer :: incxs => null() !< unit number CXS
    integer(I4B), pointer :: insto => null() !< STO enabled flag
    integer(I4B), pointer :: inobs => null() ! unit number OBS
    integer(I4B), pointer :: inoc => null() !< unit number OC
    integer(I4B), pointer :: iss => null() ! steady state flag
    integer(I4B), pointer :: inewtonur => null() ! newton under relaxation flag
  contains
    procedure :: initialize
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: model_df => swf_df
    procedure :: model_ac => swf_ac
    procedure :: model_mc => swf_mc
    procedure :: model_ar => swf_ar
    procedure :: model_rp => swf_rp
    procedure :: model_ad => swf_ad
    procedure :: model_nur => swf_nur
    procedure :: model_cf => swf_cf
    procedure :: model_fc => swf_fc
    procedure :: model_cq => swf_cq
    procedure :: model_bd => swf_bd
    procedure :: model_ot => swf_ot
    procedure :: model_da => swf_da
    procedure :: model_bdentry => swf_bdentry
    procedure :: swf_ot_obs
    procedure :: swf_ot_flow
    procedure :: swf_ot_dv
    procedure :: swf_ot_bdsummary
    procedure :: package_create
    procedure :: ftype_check
    procedure :: get_iasym => swf_get_iasym
    procedure :: create_packages
    procedure, private :: create_bndpkgs
    !procedure :: log_namfile_options
    procedure, private :: steady_period_check
  end type SwfModelType

  ! todo: these should be removed and entirely delegated to CHF and OLF
  integer(I4B), parameter :: SWF_NBASEPKG = 9
  integer(I4B), parameter :: SWF_NMULTIPKG = 50
  integer(I4B), parameter :: NIUNIT_SWF = SWF_NBASEPKG + SWF_NMULTIPKG

contains

  !> @brief Initialize common swf members
  !<
  subroutine initialize(this, modelftype, filename, id, modelname)
    ! modules
    ! dummy
    class(SwfModelType) :: this
    character(len=*), intent(in) :: modelftype !< abbreviation for model type (CHF or OLF)
    character(len=*), intent(in) :: filename !< input file
    integer(I4B), intent(in) :: id !< consecutive model number listed in mfsim.nam
    character(len=*), intent(in) :: modelname !< name of the model
    ! local

    ! Set memory path before allocation in memory manager can be done
    this%memoryPath = create_mem_path(modelname)

    ! allocate scalars
    call this%allocate_scalars(modelname)

    ! Assign values
    this%filename = filename
    this%name = modelname
    this%macronym = trim(modelftype)
    this%id = id

  end subroutine initialize

  !> @brief Allocate memory for scalar members
  !<
  subroutine allocate_scalars(this, modelname)
    ! modules
    ! dummy
    class(SwfModelType) :: this
    character(len=*), intent(in) :: modelname

    ! allocate members from parent class
    call this%NumericalModelType%allocate_scalars(modelname)

    ! allocate members that are part of model class
    call mem_allocate(this%inic, 'INIC', this%memoryPath)
    call mem_allocate(this%indfw, 'INDFW', this%memoryPath)
    call mem_allocate(this%incxs, 'INCXS', this%memoryPath)
    call mem_allocate(this%insto, 'INSTO', this%memoryPath)
    call mem_allocate(this%inobs, 'INOBS', this%memoryPath)
    call mem_allocate(this%inoc, 'INOC', this%memoryPath)
    call mem_allocate(this%iss, 'ISS', this%memoryPath)
    call mem_allocate(this%inewtonur, 'INEWTONUR', this%memoryPath)

    ! initialize
    this%inic = 0
    this%indfw = 0
    this%incxs = 0
    this%insto = 0
    this%inobs = 0
    this%inoc = 0
    this%iss = 1 !default is steady-state (i.e., no STO package)
    this%inewtonur = 0

  end subroutine allocate_scalars

  !> @brief Allocate memory for scalar members
  !<
  subroutine allocate_arrays(this)
    ! modules
    ! dummy
    class(SwfModelType) :: this
    integer(I4B) :: i

    ! allocate members from parent class
    call this%NumericalModelType%allocate_arrays()

    if (this%indfw == 0) then
      ! explicit model, so these must be manually allocated (not used)
      call mem_allocate(this%x, this%dis%nodes, 'X', this%memoryPath)
      call mem_allocate(this%rhs, this%dis%nodes, 'RHS', this%memoryPath)
      call mem_allocate(this%ibound, this%dis%nodes, 'IBOUND', this%memoryPath)
      do i = 1, this%dis%nodes
        this%x(i) = DZERO
        this%rhs(i) = DZERO
        this%ibound(i) = 1
      end do
    end if

  end subroutine allocate_arrays

  !> @brief Define packages of the model
  !<
  subroutine swf_df(this)
    ! modules
    ! dummy
    class(SwfModelType) :: this
    ! local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj

    ! call package df routines
    call this%dis%dis_df()
    call this%dfw%dfw_df(this%dis)
    call this%oc%oc_df()
    call this%budget%budget_df(NIUNIT_SWF, 'VOLUME', 'L**3')

    ! set model sizes
    this%neq = this%dis%nodes
    this%nja = this%dis%nja
    this%ia => this%dis%con%ia
    this%ja => this%dis%con%ja

    ! Allocate model arrays, now that neq and nja are known
    call this%allocate_arrays()

    ! Define packages and assign iout for time series managers
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_df(this%dis%nodes, this%dis)
    end do

    ! Store information needed for observations
    call this%obs%obs_df(this%iout, this%name, 'SWF', this%dis)

  end subroutine swf_df

  !> @brief Add the internal connections of this model to the sparse matrix
  !<
  subroutine swf_ac(this, sparse)
    ! modules
    use SparseModule, only: sparsematrix
    ! dummy
    class(SwfModelType) :: this
    type(sparsematrix), intent(inout) :: sparse
    ! local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! Add the primary grid connections of this model to sparse
    call this%dis%dis_ac(this%moffset, sparse)

    ! Add any additional connections
    !    none

    ! Add any package connections
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ac(this%moffset, sparse)
    end do

  end subroutine swf_ac

  !> @brief Map the positions of this models connections in the
  !< numerical solution coefficient matrix.
  subroutine swf_mc(this, matrix_sln)
    ! dummy
    class(SwfModelType) :: this
    class(MatrixBaseType), pointer :: matrix_sln
    ! local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! Find the position of each connection in the global ia, ja structure
    ! and store them in idxglo.
    call this%dis%dis_mc(this%moffset, this%idxglo, matrix_sln)

    ! Map any additional connections
    ! none

    ! Map any package connections
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_mc(this%moffset, matrix_sln)
    end do

  end subroutine swf_mc

  !> @brief SWF Allocate and Read
  !<
  subroutine swf_ar(this)
    ! dummy
    class(SwfModelType) :: this
    ! locals
    integer(I4B), dimension(:), allocatable :: itemp
    integer(I4B) :: ip
    class(BndType), pointer :: packobj

    ! Allocate and read modules attached to model
    if (this%inic > 0) call this%ic%ic_ar(this%x)

    ! need temporary integer variable to pass to dis_ar
    ! TODO: this should be generalized so dis_ar doesn't have to have it
    allocate (itemp(this%dis%nodes))

    ! Call dis_ar to write binary grid file
    call this%dis%dis_ar(itemp)
    if (this%indfw > 0) call this%dfw%dfw_ar(this%ibound, this%x)
    if (this%insto > 0) call this%sto%sto_ar(this%dis, this%ibound)
    if (this%inobs > 0) call this%obs%swf_obs_ar(this%ic, this%x, this%flowja)
    deallocate (itemp)

    ! set up output control
    if (this%indfw > 0) then
      call this%oc%oc_ar(this%x, this%dis, DNODATA)
    end if
    call this%budget%set_ibudcsv(this%oc%ibudcsv)

    ! Package input files now open, so allocate and read
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%set_pointers(this%dis%nodes, this%ibound, this%x, &
                                this%xold, this%flowja)
      ! Read and allocate package
      call packobj%bnd_ar()
    end do

  end subroutine swf_ar

  !> @brief Surface Water Flow Model Read and Prepare
  !<
  subroutine swf_rp(this)
    ! modules
    use TdisModule, only: readnewdata
    ! dummy
    class(SwfModelType) :: this
    ! local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! Check with TDIS on whether or not it is time to RP
    if (.not. readnewdata) return

    ! Read and prepare
    if (this%indfw > 0) call this%dfw%dfw_rp()
    if (this%inoc > 0) call this%oc%oc_rp()
    if (this%insto > 0) call this%sto%sto_rp()
    !if (this%inmvr > 0) call this%mvr%mvr_rp()
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_rp()
      call packobj%bnd_rp_obs()
    end do

    ! Check for steady state period
    call this%steady_period_check()

  end subroutine swf_rp

  !> @brief Surface Water Flow Model Time Step Advance
  !<
  subroutine swf_ad(this)
    ! modules
    use SimVariablesModule, only: isimcheck, iFailedStepRetry
    ! dummy
    class(SwfModelType) :: this
    class(BndType), pointer :: packobj
    ! local
    integer(I4B) :: irestore
    integer(I4B) :: ip, n

    ! Reset state variable
    irestore = 0
    if (iFailedStepRetry > 0) irestore = 1
    if (irestore == 0) then

      ! copy x into xold
      do n = 1, this%dis%nodes
        this%xold(n) = this%x(n)
      end do
    else

      ! copy xold into x if this time step is a redo
      do n = 1, this%dis%nodes
        this%x(n) = this%xold(n)
      end do
    end if

    ! Advance
    if (this%indfw > 0) call this%dfw%dfw_ad(irestore)
    if (this%insto > 0) call this%sto%sto_ad()
    !if (this%inmvr > 0) call this%mvr%mvr_ad()
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ad()
      if (isimcheck > 0) then
        call packobj%bnd_ck()
      end if
    end do

    ! Push simulated values to preceding time/subtime step
    call this%obs%obs_ad()

  end subroutine swf_ad

  !> @brief Calculate coefficients
  !<
  subroutine swf_cf(this, kiter)
    ! dummy
    class(SwfModelType) :: this
    integer(I4B), intent(in) :: kiter
    ! local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! Call package cf routines
    !if (this%indfw > 0) call this%dfw%dfw_cf(kiter, this%dis%nodes, this%x)
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_cf()
    end do

  end subroutine swf_cf

  !> @brief Fill coefficients
  !<
  subroutine swf_fc(this, kiter, matrix_sln, inwtflag)
    ! dummy
    class(SwfModelType) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in) :: inwtflag
    ! local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
    integer(I4B) :: inwt, inwtpak

    ! newton flags
    inwt = inwtflag
    ! if (inwtflag == 1) inwt = this%dfw%inewton
    ! inwtsto = inwtflag
    ! if (this%insto > 0) then
    !   if (inwtflag == 1) inwtsto = this%sto%inewton
    ! end if
    ! inwtcsub = inwtflag
    ! if (this%incsub > 0) then
    !   if (inwtflag == 1) inwtcsub = this%csub%inewton
    ! end if

    ! Fill standard conductance terms
    if (this%indfw > 0) call this%dfw%dfw_fc(kiter, matrix_sln, this%idxglo, &
                                             this%rhs, this%x, this%xold)
    ! storage
    if (this%insto > 0) then
      call this%sto%sto_fc(kiter, this%xold, this%x, matrix_sln, &
                           this%idxglo, this%rhs)
    end if
    ! if (this%inmvr > 0) call this%mvr%mvr_fc()
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_fc(this%rhs, this%ia, this%idxglo, matrix_sln)
    end do

    !--Fill newton terms
    if (this%indfw > 0) then
      if (inwt /= 0) then
        call this%dfw%dfw_fn(kiter, matrix_sln, this%idxglo, this%rhs, this%x)
      end if
    end if

    ! Fill newton terms for storage
    ! if (this%insto > 0) then
    !   if (inwtsto /= 0) then
    !     call this%sto%sto_fn(kiter, this%xold, this%x, matrix_sln, &
    !                          this%idxglo, this%rhs)
    !   end if
    ! end if

    ! Fill Newton terms for packages
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      inwtpak = inwtflag
      if (inwtflag == 1) inwtpak = packobj%inewton
      if (inwtpak /= 0) then
        call packobj%bnd_fn(this%rhs, this%ia, this%idxglo, matrix_sln)
      end if
    end do

  end subroutine swf_fc

  !> @brief under-relaxation
  !<
  subroutine swf_nur(this, neqmod, x, xtemp, dx, inewtonur, dxmax, locmax)
    ! modules
    use ConstantsModule, only: DONE, DP9
    ! dummy
    class(SwfModelType) :: this
    integer(I4B), intent(in) :: neqmod
    real(DP), dimension(neqmod), intent(inout) :: x
    real(DP), dimension(neqmod), intent(in) :: xtemp
    real(DP), dimension(neqmod), intent(inout) :: dx
    integer(I4B), intent(inout) :: inewtonur
    real(DP), intent(inout) :: dxmax
    integer(I4B), intent(inout) :: locmax
    ! local
    integer(I4B) :: i0
    integer(I4B) :: i1
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! apply Newton-Raphson under-relaxation if model is using
    ! the Newton-Raphson formulation and this Newton-Raphson
    ! under-relaxation is turned on.
    if (this%inewton /= 0 .and. this%inewtonur /= 0) then
      if (this%indfw > 0) then
        call this%dfw%dfw_nur(neqmod, x, xtemp, dx, inewtonur, dxmax, locmax)
      end if

      ! Call package nur routines
      i0 = this%dis%nodes + 1
      do ip = 1, this%bndlist%Count()
        packobj => GetBndFromList(this%bndlist, ip)
        if (packobj%npakeq > 0) then
          i1 = i0 + packobj%npakeq - 1
          call packobj%bnd_nur(packobj%npakeq, x(i0:i1), xtemp(i0:i1), &
                               dx(i0:i1), inewtonur, dxmax, locmax)
          i0 = i1 + 1
        end if
      end do
    end if

  end subroutine swf_nur

  !> @brief Calculate flow
  !<
  subroutine swf_cq(this, icnvg, isuppress_output)
    ! modules
    ! dummy
    class(SwfModelType) :: this
    integer(I4B), intent(in) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    ! local
    integer(I4B) :: i
    integer(I4B) :: ip
    class(BndType), pointer :: packobj

    ! Construct the flowja array.  Flowja is calculated each time, even if
    ! output is suppressed.  (flowja is positive into a cell.)  The diagonal
    ! position of the flowja array will contain the flow residual after
    ! these routines are called, so each package is responsible for adding
    ! its flow to this diagonal position.
    do i = 1, this%dis%nja
      this%flowja(i) = DZERO
    end do
    if (this%indfw > 0) call this%dfw%dfw_cq(this%x, this%xold, this%flowja)
    if (this%insto > 0) call this%sto%sto_cq(this%flowja, this%x, this%xold)

    ! Go through packages and call cq routines.  cf() routines are called
    ! first to regenerate non-linear terms to be consistent with the final
    ! head solution.
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_cf()
      call packobj%bnd_cq(this%x, this%flowja)
    end do

  end subroutine swf_cq

  !> @brief Model Budget
  !<
  subroutine swf_bd(this, icnvg, isuppress_output)
    ! modules
    use SparseModule, only: csr_diagsum
    ! dummy
    class(SwfModelType) :: this
    integer(I4B), intent(in) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    ! local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj

    ! Finalize calculation of flowja by adding face flows to the diagonal.
    ! This results in the flow residual being stored in the diagonal
    ! position for each cell.
    call csr_diagsum(this%dis%con%ia, this%flowja)

    ! Budget routines (start by resetting).  Sole purpose of this section
    ! is to add in and outs to model budget.  All ins and out for a model
    ! should be added here to this%budget.  In a subsequent exchange call,
    ! exchange flows might also be added.
    call this%budget%reset()
    if (this%insto > 0) call this%sto%sto_bd(isuppress_output, this%budget)
    if (this%indfw > 0) call this%dfw%dfw_bd(isuppress_output, this%budget)
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_bd(this%budget)
    end do

    ! dfw velocities have to be calculated here, after swf-swf exchanges
    ! have passed in their contributions from exg_cq()
    if (this%indfw > 0) then
      if (this%dfw%icalcvelocity /= 0) then
        call this%dfw%calc_velocity(this%flowja)
      end if
    end if

  end subroutine swf_bd

  !> @brief Surface Water Flow Model Output
  !<
  subroutine swf_ot(this)
    ! modules
    use TdisModule, only: tdis_ot, endofperiod
    ! dummy
    class(SwfModelType) :: this
    ! local
    integer(I4B) :: idvsave
    integer(I4B) :: idvprint
    integer(I4B) :: icbcfl
    integer(I4B) :: icbcun
    integer(I4B) :: ibudfl
    integer(I4B) :: ipflag
    integer(I4B) :: icnvg = 1

    ! Set write and print flags
    idvsave = 0
    idvprint = 0
    icbcfl = 0
    ibudfl = 0
    if (this%oc%oc_save('QOUTFLOW')) idvsave = 1
    if (this%oc%oc_print('QOUTFLOW')) idvprint = 1
    if (this%oc%oc_save('BUDGET')) icbcfl = 1
    if (this%oc%oc_print('BUDGET')) ibudfl = 1
    icbcun = this%oc%oc_save_unit('BUDGET')

    ! Override ibudfl and idvprint flags for nonconvergence
    ! and end of period
    ibudfl = this%oc%set_print_flag('BUDGET', icnvg, endofperiod)
    idvprint = this%oc%set_print_flag('QOUTFLOW', icnvg, endofperiod)

    ! Calculate and save observations
    call this%swf_ot_obs()

    ! Save and print flows
    call this%swf_ot_flow(icbcfl, ibudfl, icbcun)

    ! Save and print dependent variables
    call this%swf_ot_dv(idvsave, idvprint, ipflag)

    ! Print budget summaries
    call this%swf_ot_bdsummary(ibudfl, ipflag)

    ! Timing Output; if any dependent variables or budgets
    ! are printed, then ipflag is set to 1.
    if (ipflag == 1) call tdis_ot(this%iout)

  end subroutine swf_ot

  !> @brief Surface Water Flow Model output observations
  !<
  subroutine swf_ot_obs(this)
    class(SwfModelType) :: this
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! Calculate and save SWF observations
    call this%obs%obs_bd()
    call this%obs%obs_ot()

    ! Calculate and save dfw observations
    if (this%indfw > 0) then
      call this%dfw%dfw_bd_obs()
      call this%dfw%obs%obs_ot()
    end if

    ! Calculate and save package observations
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_bd_obs()
      call packobj%bnd_ot_obs()
    end do

  end subroutine swf_ot_obs

  !> @brief Surface Water Flow Model output flows
  !<
  subroutine swf_ot_flow(this, icbcfl, ibudfl, icbcun)
    class(SwfModelType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    integer(I4B), intent(in) :: icbcun
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! Save SWF flows
    if (this%insto > 0) then
      call this%sto%sto_save_model_flows(icbcfl, icbcun)
    end if
    if (this%indfw > 0) then
      call this%dfw%dfw_save_model_flows(this%flowja, icbcfl, icbcun)
    end if
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_model_flows(icbcfl=icbcfl, ibudfl=0, icbcun=icbcun)
    end do

    ! Save advanced package flows
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_package_flows(icbcfl=icbcfl, ibudfl=0)
    end do
    ! if (this%inmvr > 0) then
    !   call this%mvr%mvr_ot_saveflow(icbcfl, ibudfl)
    ! end if

    ! Print SWF flows
    if (this%indfw > 0) then
      call this%dfw%dfw_print_model_flows(ibudfl, this%flowja)
    end if
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_model_flows(icbcfl=icbcfl, ibudfl=ibudfl, icbcun=0)
    end do

    ! Print advanced package flows
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_package_flows(icbcfl=0, ibudfl=ibudfl)
    end do
    ! if (this%inmvr > 0) then
    !   call this%mvr%mvr_ot_printflow(icbcfl, ibudfl)
    ! end if

  end subroutine swf_ot_flow

  !> @brief Surface Water Flow Model output dependent variable
  !<
  subroutine swf_ot_dv(this, idvsave, idvprint, ipflag)
    class(SwfModelType) :: this
    integer(I4B), intent(in) :: idvsave
    integer(I4B), intent(in) :: idvprint
    integer(I4B), intent(inout) :: ipflag
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! Print advanced package dependent variables
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_dv(idvsave, idvprint)
    end do
    !
    ! save stage and print stage (if implemented)
    call this%oc%oc_ot(ipflag)

  end subroutine swf_ot_dv

  !> @brief Surface Water Flow Model output budget summary
  !<
  subroutine swf_ot_bdsummary(this, ibudfl, ipflag)
    use TdisModule, only: kstp, kper, totim, delt
    class(SwfModelType) :: this
    integer(I4B), intent(in) :: ibudfl
    integer(I4B), intent(inout) :: ipflag
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! Package budget summary
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_bdsummary(kstp, kper, this%iout, ibudfl)
    end do

    ! mover budget summary
    ! if (this%inmvr > 0) then
    !   call this%mvr%mvr_ot_bdsummary(ibudfl)
    ! end if

    ! model budget summary
    call this%budget%finalize_step(delt)
    if (ibudfl /= 0) then
      ipflag = 1
      call this%budget%budget_ot(kstp, kper, this%iout)
    end if

    ! Write to budget csv every time step
    call this%budget%writecsv(totim)

  end subroutine swf_ot_bdsummary

  !> @brief Deallocate
  !<
  subroutine swf_da(this)
    ! modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    ! dummy
    class(SwfModelType) :: this
    ! local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj

    ! Deallocate idm memory
    call memorystore_remove(this%name, 'NAM', idm_context)
    call memorystore_remove(component=this%name, context=idm_context)

    ! Internal flow packages deallocate
    call this%dis%dis_da()
    if (this%insto > 0) call this%sto%sto_da()
    if (this%inic > 0) call this%ic%ic_da()
    if (this%indfw > 0) call this%dfw%dfw_da()
    call this%cxs%cxs_da()
    call this%obs%obs_da()
    call this%oc%oc_da()
    call this%budget%budget_da()

    ! Internal package objects
    deallocate (this%dis)
    deallocate (this%budget)
    deallocate (this%obs)
    deallocate (this%oc)

    ! Boundary packages
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_da()
      deallocate (packobj)
    end do

    ! Scalars
    call mem_deallocate(this%inic)
    call mem_deallocate(this%indfw)
    call mem_deallocate(this%incxs)
    call mem_deallocate(this%insto)
    call mem_deallocate(this%inobs)
    call mem_deallocate(this%inoc)
    call mem_deallocate(this%iss)
    call mem_deallocate(this%inewtonur)

    ! Arrays

    ! NumericalModelType
    call this%NumericalModelType%model_da()

  end subroutine swf_da

  !> @brief Surface Flow Model Budget Entry
  !!
  !! This subroutine adds a budget entry to the flow budget.  It was added as
  !! a method for the swf model object so that the exchange object could add its
  !< contributions.
  subroutine swf_bdentry(this, budterm, budtxt, rowlabel)
    ! modules
    use ConstantsModule, only: LENBUDTXT
    use TdisModule, only: delt
    ! dummy
    class(SwfModelType) :: this
    real(DP), dimension(:, :), intent(in) :: budterm
    character(len=LENBUDTXT), dimension(:), intent(in) :: budtxt
    character(len=*), intent(in) :: rowlabel

    call this%budget%addentry(budterm, delt, budtxt, rowlabel=rowlabel)

  end subroutine swf_bdentry

  !> @brief Create boundary condition packages for this model
  !<
  subroutine package_create(this, filtyp, ipakid, ipaknum, pakname, mempath, &
                            inunit, iout)
    ! modules
    use SwfFlwModule, only: flw_create
    use ChdModule, only: chd_create
    use SwfCdbModule, only: cdb_create
    use SwfZdgModule, only: zdg_create
    use SwfPcpModule, only: pcp_create
    use SwfEvpModule, only: evp_create
    ! dummy
    class(SwfModelType) :: this
    character(len=*), intent(in) :: filtyp
    integer(I4B), intent(in) :: ipakid
    integer(I4B), intent(in) :: ipaknum
    character(len=*), intent(in) :: pakname
    character(len=*), intent(in) :: mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    ! local
    class(BndType), pointer :: packobj
    class(BndType), pointer :: packobj2
    integer(I4B) :: ip

    ! This part creates the package object
    select case (filtyp)
    case ('FLW6')
      call flw_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, mempath)
      packobj%ictMemPath = ''
    case ('CHD6')
      call chd_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, mempath)
      packobj%ictMemPath = create_mem_path(this%name, 'DFW')
    case ('CDB6')
      call cdb_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, mempath, this%dis, this%cxs, &
                      this%dfw%lengthconv, this%dfw%timeconv)
    case ('ZDG6')
      call zdg_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, mempath, this%dis, this%cxs, this%dfw%unitconv)
    case ('PCP6')
      call pcp_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, mempath, this%dis, this%dfw, this%cxs)
    case ('EVP6')
      call evp_create(packobj, ipakid, ipaknum, inunit, iout, this%name, &
                      pakname, mempath, this%dis, this%dfw, this%cxs)
    case default
      write (errmsg, *) 'Invalid package type: ', filtyp
      call store_error(errmsg)
      call store_error_filename(this%filename)
    end select

    ! Check to make sure that the package name is unique, then store a
    ! pointer to the package in the model bndlist
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

  !> @brief Check to make sure required input files have been specified
  !<
  subroutine ftype_check(this, indis)
    ! modules
    ! dummy
    class(SwfModelType) :: this
    integer(I4B), intent(in) :: indis
    ! local

    ! Check for required packages. Stop if not present.
    if (indis == 0) then
      write (errmsg, '(a)') &
        'Discretization Package (DISV1D6 or DIS2D6) not specified.'
      call store_error(errmsg)
    end if
    if (this%inic == 0 .and. this%indfw /= 0) then
      write (errmsg, '(a)') &
        'Initial Conditions (IC6) must be specified if the Diffusive &
        &Wave (DFW) package is used.'
      call store_error(errmsg)
    end if
    if (this%indfw == 0) then
      write (errmsg, '(1x,a)') &
        'DFW6 Package must be specified.'
      call store_error(errmsg)
    end if
    if (this%incxs > 0) then
      if (this%dis%ndim /= 1) then
        write (errmsg, '(1x,a)') &
          'CXS6 Package can only be used for one-dimensional discretization &
          &package (DISV1D).  Remove CXS6 Package from model or convert &
          &discretization package to DISV1D.'
        call store_error(errmsg)
      end if
    end if
    if (count_errors() > 0) then
      write (errmsg, '(a)') 'One or more required package(s) not specified &
        &and/or there are package incompatibilities.'
      call store_error(errmsg)
      call store_error_filename(this%filename)
    end if

  end subroutine ftype_check

  !> @brief Source package info and begin to process
  !<
  subroutine create_bndpkgs(this, bndpkgs, pkgtypes, pkgnames, &
                            mempaths, inunits)
    ! modules
    use ConstantsModule, only: LENPACKAGENAME
    use CharacterStringModule, only: CharacterStringType
    ! dummy
    class(SwfModelType) :: this
    integer(I4B), dimension(:), allocatable, intent(inout) :: bndpkgs
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(inout) :: pkgtypes
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(inout) :: pkgnames
    type(CharacterStringType), dimension(:), contiguous, &
      pointer, intent(inout) :: mempaths
    integer(I4B), dimension(:), contiguous, &
      pointer, intent(inout) :: inunits
    ! local
    integer(I4B) :: ipakid, ipaknum
    character(len=LENFTYPE) :: pkgtype, bndptype
    character(len=LENPACKAGENAME) :: pkgname
    character(len=LENMEMPATH) :: mempath
    integer(I4B), pointer :: inunit
    integer(I4B) :: n

    if (allocated(bndpkgs)) then

      ! create stress packages
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

      ! cleanup
      deallocate (bndpkgs)
    end if

  end subroutine create_bndpkgs

  !> @brief Source package info and begin to process
  !<
  subroutine create_packages(this)
    ! modules
    use ConstantsModule, only: LINELENGTH, LENPACKAGENAME
    use CharacterStringModule, only: CharacterStringType
    use ArrayHandlersModule, only: expandarray
    use MemoryManagerModule, only: mem_setptr
    use SimVariablesModule, only: idm_context
    use Disv1dModule, only: disv1d_cr
    use Dis2dModule, only: dis2d_cr
    use Disv2dModule, only: disv2d_cr
    use SwfDfWModule, only: dfw_cr
    use SwfCxsModule, only: cxs_cr
    use SwfStoModule, only: sto_cr
    use SwfIcModule, only: ic_cr
    use SwfOcModule, only: oc_cr
    ! dummy
    class(SwfModelType) :: this
    ! local
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgtypes => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: pkgnames => null()
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: mempaths => null()
    integer(I4B), dimension(:), contiguous, &
      pointer :: inunits => null()
    character(len=LENMEMPATH) :: model_mempath
    character(len=LENPACKAGETYPE) :: pkgtype
    character(len=LENPACKAGENAME) :: pkgname
    character(len=LENMEMPATH) :: mempath
    integer(I4B), pointer :: inunit
    integer(I4B), dimension(:), allocatable :: bndpkgs
    integer(I4B) :: n
    integer(I4B) :: indis = 0 ! DIS enabled flag
    character(len=LENMEMPATH) :: mempathic = ''
    character(len=LENMEMPATH) :: mempathdfw = ''
    character(len=LENMEMPATH) :: mempathcxs = ''
    character(len=LENMEMPATH) :: mempathoc = ''
    character(len=LENMEMPATH) :: mempathsto = ''

    ! set input model memory path
    model_mempath = create_mem_path(component=this%name, context=idm_context)

    ! set pointers to model path package info
    call mem_setptr(pkgtypes, 'PKGTYPES', model_mempath)
    call mem_setptr(pkgnames, 'PKGNAMES', model_mempath)
    call mem_setptr(mempaths, 'MEMPATHS', model_mempath)
    call mem_setptr(inunits, 'INUNITS', model_mempath)

    do n = 1, size(pkgtypes)

      ! attributes for this input package
      pkgtype = pkgtypes(n)
      pkgname = pkgnames(n)
      mempath = mempaths(n)
      inunit => inunits(n)

      ! create dis package as it is a prerequisite for other packages
      select case (pkgtype)
      case ('DISV1D6')
        indis = 1
        call disv1d_cr(this%dis, this%name, mempath, indis, this%iout)
      case ('DIS2D6')
        indis = 1
        call dis2d_cr(this%dis, this%name, mempath, indis, this%iout)
      case ('DISV2D6')
        indis = 1
        call disv2d_cr(this%dis, this%name, mempath, indis, this%iout)
      case ('DFW6')
        this%indfw = 1
        mempathdfw = mempath
      case ('CXS6')
        this%incxs = 1
        mempathcxs = mempath
      case ('STO6')
        this%insto = 1
        mempathsto = mempath
      case ('IC6')
        this%inic = 1
        mempathic = mempath
      case ('OC6')
        !this%inoc = inunit
        this%inoc = 1
        mempathoc = mempath
      case ('OBS6')
        this%inobs = inunit
      case ('CHD6', 'FLW6', 'CDB6', 'ZDG6', 'PCP6', 'EVP6')
        call expandarray(bndpkgs)
        bndpkgs(size(bndpkgs)) = n
      case default
        ! TODO
      end select
    end do

    ! Create packages that are tied directly to model
    if (this%inic > 0) then
      call ic_cr(this%ic, this%name, mempathic, this%inic, this%iout, &
                 this%dis)
    end if
    call cxs_cr(this%cxs, this%name, mempathcxs, this%incxs, this%iout, &
                this%dis)
    if (this%indfw > 0) then
      call dfw_cr(this%dfw, this%name, mempathdfw, this%indfw, this%iout, &
                  this%cxs)
    end if
    if (this%insto > 0) then
      call sto_cr(this%sto, this%name, mempathsto, this%insto, this%iout, &
                  this%cxs)
    end if
    call oc_cr(this%oc, this%name, mempathoc, this%inoc, this%iout)
    call swf_obs_cr(this%obs, this%inobs)

    ! Check to make sure that required ftype's have been specified
    call this%ftype_check(indis)

    call this%create_bndpkgs(bndpkgs, pkgtypes, pkgnames, mempaths, inunits)

  end subroutine create_packages

  !> @brief Check for steady state period
  !!
  !! Write warning message if steady state
  !! period and adaptive time stepping is
  !< active for the period
  subroutine steady_period_check(this)
    ! modules
    use TdisModule, only: kper
    use AdaptiveTimeStepModule, only: isAdaptivePeriod
    use SimVariablesModule, only: warnmsg
    use SimModule, only: store_warning
    ! dummy
    class(SwfModelType) :: this

    if (this%iss == 1) then
      if (isAdaptivePeriod(kper)) then
        write (warnmsg, '(a,a,a,i0,a)') &
          'SWF Model (', trim(this%name), ') is steady state for period ', &
          kper, ' and adaptive time stepping is active.  Adaptive time &
          &stepping may not work properly for steady-state conditions.'
        call store_warning(warnmsg)
      end if
    end if

  end subroutine steady_period_check

  !> @brief return 1 if any package causes the matrix to be asymmetric.
  !< Otherwise return 0.
  function swf_get_iasym(this) result(iasym)
    class(SwfModelType) :: this
    ! local
    integer(I4B) :: iasym
    integer(I4B) :: ip
    class(BndType), pointer :: packobj

    ! Start by setting iasym to zero
    iasym = 0

    ! DFW
    if (this%indfw > 0) then
      iasym = 1
    end if

    ! Check for any packages that introduce matrix asymmetry
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      if (packobj%iasym /= 0) iasym = 1
    end do

  end function swf_get_iasym

end module SwfModule
