module PrtModule
  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use InputOutputModule, only: ParseLine, upcase, lowcase
  use ConstantsModule, only: LENFTYPE, LENMEMPATH, DZERO, DONE, &
                             LENPAKLOC, LENPACKAGETYPE, LENBUDTXT, MNORMAL, &
                             LINELENGTH
  use VersionModule, only: write_listfile_header
  use NumericalModelModule, only: NumericalModelType
  use BaseModelModule, only: BaseModelType
  use BndModule, only: BndType, AddBndToList, GetBndFromList
  use DisModule, only: DisType, dis_cr
  use DisvModule, only: DisvType, disv_cr
  use DisuModule, only: DisuType, disu_cr
  use PrtPrpModule, only: PrtPrpType
  use PrtFmiModule, only: PrtFmiType
  use PrtMipModule, only: PrtMipType
  use PrtOcModule, only: PrtOcType
  use BudgetModule, only: BudgetType
  use ListModule, only: ListType
  use ParticleModule, only: ParticleType, create_particle
  use TrackFileModule, only: TrackFileType
  use TrackControlModule, only: TrackControlType
  use SimModule, only: count_errors, store_error, store_error_filename
  use MemoryManagerModule, only: mem_allocate
  use MethodModule, only: MethodType

  implicit none

  private
  public :: prt_cr
  public :: PrtModelType
  public :: PRT_NBASEPKG, PRT_NMULTIPKG
  public :: PRT_BASEPKG, PRT_MULTIPKG

  integer(I4B), parameter :: NBDITEMS = 1
  character(len=LENBUDTXT), dimension(NBDITEMS) :: budtxt
  data budtxt/'         STORAGE'/

  !> @brief Particle tracking (PRT) model
  type, extends(NumericalModelType) :: PrtModelType
    type(PrtFmiType), pointer :: fmi => null() ! flow model interface
    type(PrtMipType), pointer :: mip => null() ! model input package
    type(PrtOcType), pointer :: oc => null() ! output control package
    type(BudgetType), pointer :: budget => null() ! budget object
    class(MethodType), pointer :: method => null() ! tracking method
    type(TrackControlType), pointer :: trackctl ! track control
    integer(I4B), pointer :: infmi => null() ! unit number FMI
    integer(I4B), pointer :: inmip => null() ! unit number MIP
    integer(I4B), pointer :: inmvt => null() ! unit number MVT
    integer(I4B), pointer :: inmst => null() ! unit number MST
    integer(I4B), pointer :: inadv => null() ! unit number ADV
    integer(I4B), pointer :: indsp => null() ! unit number DSP
    integer(I4B), pointer :: inssm => null() ! unit number SSM
    integer(I4B), pointer :: inoc => null() ! unit number OC
    integer(I4B), pointer :: nprp => null() ! number of PRP packages in the model
    real(DP), dimension(:), pointer, contiguous :: masssto => null() !< particle mass storage in cells, new value
    real(DP), dimension(:), pointer, contiguous :: massstoold => null() !< particle mass storage in cells, old value
    real(DP), dimension(:), pointer, contiguous :: ratesto => null() !< particle mass storage rate in cells
  contains
    ! Override BaseModelType procs
    procedure :: model_df => prt_df
    procedure :: model_ar => prt_ar
    procedure :: model_rp => prt_rp
    procedure :: model_ad => prt_ad
    procedure :: model_cq => prt_cq
    procedure :: model_bd => prt_bd
    procedure :: model_ot => prt_ot
    procedure :: model_da => prt_da
    procedure :: model_solve => prt_solve

    ! Private utilities
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure, private :: package_create
    procedure, private :: ftype_check
    procedure, private :: prt_ot_flow
    procedure, private :: prt_ot_saveflow
    procedure, private :: prt_ot_printflow
    procedure, private :: prt_ot_dv
    procedure, private :: prt_ot_bdsummary
    procedure, private :: prt_cq_sto
    procedure, private :: create_packages
    procedure, private :: create_bndpkgs
    procedure, private :: log_namfile_options

  end type PrtModelType

  !> @brief PRT base package array descriptors
  !!
  !! PRT6 model base package types.  Only listed packages are candidates
  !! for input and these will be loaded in the order specified.
  !<
  integer(I4B), parameter :: PRT_NBASEPKG = 50
  character(len=LENPACKAGETYPE), dimension(PRT_NBASEPKG) :: PRT_BASEPKG
  data PRT_BASEPKG/'DIS6 ', 'DISV6', 'DISU6', 'IC6  ', 'MST6 ', & !  5
                  &'ADV6 ', 'DSP6 ', 'SSM6 ', 'MIP6 ', 'CNC6 ', & ! 10
                  &'OC6  ', '     ', 'FMI6 ', '     ', 'IST6 ', & ! 15
                  &'LKT6 ', 'SFT6 ', 'MWT6 ', 'UZT6 ', 'MVT6 ', & ! 20
                  &'API6 ', '     ', '     ', '     ', '     ', & ! 25
                  25*'     '/ ! 50

  !> @brief PRT multi package array descriptors
  !!
  !! PRT6 model multi-instance package types.  Only listed packages are
  !! candidates for input and these will be loaded in the order specified.
  !<
  integer(I4B), parameter :: PRT_NMULTIPKG = 50
  character(len=LENPACKAGETYPE), dimension(PRT_NMULTIPKG) :: PRT_MULTIPKG
  data PRT_MULTIPKG/'PRP6 ', '     ', '     ', '     ', '     ', & !  5
                   &45*'     '/ ! 50

  ! size of supported model package arrays
  integer(I4B), parameter :: NIUNIT_PRT = PRT_NBASEPKG + PRT_NMULTIPKG

contains

  !> @brief Create a new particle tracking model object
  subroutine prt_cr(filename, id, modelname)
    ! modules
    use ListsModule, only: basemodellist
    use BaseModelModule, only: AddBaseModelToList
    use ConstantsModule, only: LINELENGTH, LENPACKAGENAME
    use CompilerVersion
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerExtModule, only: mem_set_value
    use SimVariablesModule, only: idm_context
    use GwfNamInputModule, only: GwfNamParamFoundType
    ! dummy
    character(len=*), intent(in) :: filename
    integer(I4B), intent(in) :: id
    character(len=*), intent(in) :: modelname
    ! local
    type(PrtModelType), pointer :: this
    class(BaseModelType), pointer :: model
    character(len=LENMEMPATH) :: input_mempath
    character(len=LINELENGTH) :: lst_fname
    type(GwfNamParamFoundType) :: found

    ! Allocate a new PRT Model (this)
    allocate (this)

    ! Set this before any allocs in the memory manager can be done
    this%memoryPath = create_mem_path(modelname)

    ! Allocate track control object
    allocate (this%trackctl)

    ! Allocate scalars and add model to basemodellist
    call this%allocate_scalars(modelname)
    model => this
    call AddBaseModelToList(basemodellist, model)

    ! Assign variables
    this%filename = filename
    this%name = modelname
    this%macronym = 'PRT'
    this%id = id

    ! Set input model namfile memory path
    input_mempath = create_mem_path(modelname, 'NAM', idm_context)

    ! Copy options from input context
    call mem_set_value(this%iprpak, 'PRINT_INPUT', input_mempath, &
                       found%print_input)
    call mem_set_value(this%iprflow, 'PRINT_FLOWS', input_mempath, &
                       found%print_flows)
    call mem_set_value(this%ipakcb, 'SAVE_FLOWS', input_mempath, &
                       found%save_flows)

    ! Create the list file
    call this%create_lstfile(lst_fname, filename, found%list, &
                             'PARTICLE TRACKING MODEL (PRT)')

    ! Activate save_flows if found
    if (found%save_flows) then
      this%ipakcb = -1
    end if

    ! Log options
    if (this%iout > 0) then
      call this%log_namfile_options(found)
    end if

    ! Create model packages
    call this%create_packages()
  end subroutine prt_cr

  !> @brief Define packages
  !!
  !! (1) call df routines for each package
  !! (2) set variables and pointers
  !<
  subroutine prt_df(this)
    ! modules
    use PrtPrpModule, only: PrtPrpType
    ! dummy
    class(PrtModelType) :: this
    ! local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj

    ! Define packages and utility objects
    call this%dis%dis_df()
    call this%fmi%fmi_df(this%dis, 1)
    call this%oc%oc_df()
    call this%budget%budget_df(NIUNIT_PRT, 'MASS', 'M')

    ! Define packages and assign iout for time series managers
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_df(this%dis%nodes, this%dis)
      packobj%TsManager%iout = this%iout
      packobj%TasManager%iout = this%iout
    end do

    ! Allocate model arrays
    call this%allocate_arrays()

  end subroutine prt_df

  !> @brief Allocate and read
  !!
  !! (1) allocates and reads packages part of this model,
  !! (2) allocates memory for arrays part of this model object
  !<
  subroutine prt_ar(this)
    ! modules
    use ConstantsModule, only: DHNOFLO
    use PrtPrpModule, only: PrtPrpType
    use PrtMipModule, only: PrtMipType
    use MethodPoolModule, only: method_dis, method_disv
    ! dummy
    class(PrtModelType) :: this
    ! locals
    integer(I4B) :: ip
    class(BndType), pointer :: packobj

    ! Allocate and read modules attached to model
    call this%fmi%fmi_ar(this%ibound)
    if (this%inmip > 0) call this%mip%mip_ar()

    ! set up output control
    call this%oc%oc_ar(this%dis, DHNOFLO)
    call this%budget%set_ibudcsv(this%oc%ibudcsv)

    ! Package input files now open, so allocate and read
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      select type (packobj)
      type is (PrtPrpType)
        call packobj%prp_set_pointers(this%ibound, this%mip%izone, &
                                      this%trackctl)
      end select
      ! Read and allocate package
      call packobj%bnd_ar()
    end do

    ! Initialize tracking method
    select type (dis => this%dis)
    type is (DisType)
      call method_dis%init( &
        fmi=this%fmi, &
        trackctl=this%trackctl, &
        izone=this%mip%izone, &
        flowja=this%flowja, &
        porosity=this%mip%porosity, &
        retfactor=this%mip%retfactor, &
        tracktimes=this%oc%tracktimes)
      this%method => method_dis
    type is (DisvType)
      call method_disv%init( &
        fmi=this%fmi, &
        trackctl=this%trackctl, &
        izone=this%mip%izone, &
        flowja=this%flowja, &
        porosity=this%mip%porosity, &
        retfactor=this%mip%retfactor, &
        tracktimes=this%oc%tracktimes)
      this%method => method_disv
    end select

    ! Initialize track output files and reporting options
    if (this%oc%itrkout > 0) &
      call this%trackctl%init_track_file(this%oc%itrkout)
    if (this%oc%itrkcsv > 0) &
      call this%trackctl%init_track_file(this%oc%itrkcsv, csv=.true.)
    call this%trackctl%set_track_events( &
      this%oc%trackrelease, &
      this%oc%trackexit, &
      this%oc%tracktimestep, &
      this%oc%trackterminate, &
      this%oc%trackweaksink, &
      this%oc%trackusertime)
  end subroutine prt_ar

  !> @brief Read and prepare (calls package read and prepare routines)
  subroutine prt_rp(this)
    use TdisModule, only: readnewdata
    ! dummy
    class(PrtModelType) :: this
    ! local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! Check with TDIS on whether or not it is time to RP
    if (.not. readnewdata) return

    ! Read and prepare
    if (this%inoc > 0) call this%oc%oc_rp()
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_rp()
    end do
  end subroutine prt_rp

  !> @brief Time step advance (calls package advance subroutines)
  subroutine prt_ad(this)
    ! modules
    use SimVariablesModule, only: isimcheck, iFailedStepRetry
    ! dummy
    class(PrtModelType) :: this
    class(BndType), pointer :: packobj
    ! local
    integer(I4B) :: irestore
    integer(I4B) :: ip, n, i

    ! Reset state variable
    irestore = 0
    if (iFailedStepRetry > 0) irestore = 1

    ! Copy masssto into massstoold
    do n = 1, this%dis%nodes
      this%massstoold(n) = this%masssto(n)
    end do

    ! Advance fmi
    call this%fmi%fmi_ad()

    ! Advance
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ad()
      if (isimcheck > 0) then
        call packobj%bnd_ck()
      end if
    end do
    !
    ! Initialize the flowja array.  Flowja is calculated each time,
    !    even if output is suppressed.  (Flowja represents flow of particle
    !    mass and is positive into a cell.  Currently, each particle is assigned
    !    unit mass.)  Flowja is updated continually as particles are tracked
    !    over the time step and at the end of the time step.  The diagonal
    !    position of the flowja array will contain the flow residual.
    do i = 1, this%dis%nja
      this%flowja(i) = DZERO
    end do
  end subroutine prt_ad

  !> @brief Calculate intercell flow (flowja)
  subroutine prt_cq(this, icnvg, isuppress_output)
    ! modules
    use SparseModule, only: csr_diagsum
    use TdisModule, only: delt
    use PrtPrpModule, only: PrtPrpType
    ! dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    ! local
    integer(I4B) :: i
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    real(DP) :: tled

    ! Flowja is calculated each time, even if output is suppressed.
    !    Flowja represents flow of particle mass and is positive into a cell.
    !    Currently, each particle is assigned unit mass.
    !
    ! Reciprocal of time step size.
    tled = DONE / delt
    !
    ! Flowja was updated continually as particles were tracked over the
    !    time step.  At this point, flowja contains the net particle mass
    !    exchanged between cells during the time step.  To convert these to
    !    flow rates (particle mass per time), divide by the time step size.
    do i = 1, this%dis%nja
      this%flowja(i) = this%flowja(i) * tled
    end do

    ! Particle mass storage
    call this%prt_cq_sto()

    ! Go through packages and call cq routines. Just a formality.
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_cq(this%masssto, this%flowja)
    end do

    ! Finalize calculation of flowja by adding face flows to the diagonal.
    !    This results in the flow residual being stored in the diagonal
    !    position for each cell.
    call csr_diagsum(this%dis%con%ia, this%flowja)
  end subroutine prt_cq

  !> @brief Calculate particle mass storage
  subroutine prt_cq_sto(this)
    use ParticleModule, only: LVL_CELL
    ! modules
    use TdisModule, only: delt
    use PrtPrpModule, only: PrtPrpType
    ! dummy
    class(PrtModelType) :: this
    ! local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    integer(I4B) :: n
    integer(I4B) :: np
    integer(I4B) :: idiag
    integer(I4B) :: istatus
    real(DP) :: tled
    real(DP) :: rate

    ! Reciprocal of time step size.
    tled = DONE / delt

    ! Particle mass storage rate
    do n = 1, this%dis%nodes
      this%masssto(n) = DZERO
      this%ratesto(n) = DZERO
    end do
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      select type (packobj)
      type is (PrtPrpType)
        do np = 1, packobj%nparticles
          istatus = packobj%particles%istatus(np)
          ! this may need to change if istatus flags change
          if ((istatus > 0) .and. (istatus /= 8)) then
            n = packobj%particles%itrdomain(np, LVL_CELL)
            ! Each particle currently assigned unit mass
            this%masssto(n) = this%masssto(n) + DONE
          end if
        end do
      end select
    end do
    do n = 1, this%dis%nodes
      rate = -(this%masssto(n) - this%massstoold(n)) * tled
      this%ratesto(n) = rate
      idiag = this%dis%con%ia(n)
      this%flowja(idiag) = this%flowja(idiag) + rate
    end do
  end subroutine prt_cq_sto

  !> @brief Calculate flows and budget
  !!
  !! (1) Calculate intercell flows (flowja)
  !! (2) Calculate package contributions to model budget
  !!
  !<
  subroutine prt_bd(this, icnvg, isuppress_output)
    ! modules
    use TdisModule, only: delt
    use BudgetModule, only: rate_accumulator
    ! dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    ! local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj
    real(DP) :: rin
    real(DP) :: rout

    ! Budget routines (start by resetting).  Sole purpose of this section
    !    is to add in and outs to model budget.  All ins and out for a model
    !    should be added here to this%budget.  In a subsequent exchange call,
    !    exchange flows might also be added.
    call this%budget%reset()
    call rate_accumulator(this%ratesto, rin, rout)
    call this%budget%addentry(rin, rout, delt, budtxt(1), &
                              isuppress_output, '             PRT')
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_bd(this%budget)
    end do
  end subroutine prt_bd

  !> @brief Print and/or save model output
  subroutine prt_ot(this)
    use TdisModule, only: tdis_ot, endofperiod
    ! dummy
    class(PrtModelType) :: this
    ! local
    integer(I4B) :: idvsave
    integer(I4B) :: idvprint
    integer(I4B) :: icbcfl
    integer(I4B) :: icbcun
    integer(I4B) :: ibudfl
    integer(I4B) :: ipflag

    ! Note: particle tracking output is handled elsewhere

    ! Set write and print flags
    idvsave = 0
    idvprint = 0
    icbcfl = 0
    ibudfl = 0
    if (this%oc%oc_save('CONCENTRATION')) idvsave = 1
    if (this%oc%oc_print('CONCENTRATION')) idvprint = 1
    if (this%oc%oc_save('BUDGET')) icbcfl = 1
    if (this%oc%oc_print('BUDGET')) ibudfl = 1
    icbcun = this%oc%oc_save_unit('BUDGET')

    ! Override ibudfl and idvprint flags for nonconvergence
    !    and end of period
    ibudfl = this%oc%set_print_flag('BUDGET', 1, endofperiod)
    idvprint = this%oc%set_print_flag('CONCENTRATION', 1, endofperiod)

    ! Save and print flows
    call this%prt_ot_flow(icbcfl, ibudfl, icbcun)

    ! Save and print dependent variables
    call this%prt_ot_dv(idvsave, idvprint, ipflag)

    ! Print budget summaries
    call this%prt_ot_bdsummary(ibudfl, ipflag)

    ! Timing Output; if any dependent variables or budgets
    !    are printed, then ipflag is set to 1.
    if (ipflag == 1) call tdis_ot(this%iout)
  end subroutine prt_ot

  !> @brief Save flows
  subroutine prt_ot_flow(this, icbcfl, ibudfl, icbcun)
    use PrtPrpModule, only: PrtPrpType
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    integer(I4B), intent(in) :: icbcun
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! Save PRT flows
    call this%prt_ot_saveflow(this%dis%nja, this%flowja, icbcfl, icbcun)
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_model_flows(icbcfl=icbcfl, ibudfl=0, icbcun=icbcun)
    end do

    ! Save advanced package flows
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_package_flows(icbcfl=icbcfl, ibudfl=0)
    end do

    ! Print PRT flows
    call this%prt_ot_printflow(ibudfl, this%flowja)
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_model_flows(icbcfl=icbcfl, ibudfl=ibudfl, icbcun=0)
    end do

    ! Print advanced package flows
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_package_flows(icbcfl=0, ibudfl=ibudfl)
    end do
  end subroutine prt_ot_flow

  !> @brief Save intercell flows
  subroutine prt_ot_saveflow(this, nja, flowja, icbcfl, icbcun)
    ! dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: nja
    real(DP), dimension(nja), intent(in) :: flowja
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: icbcun
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
      call this%dis%record_connection_array(flowja, ibinun, this%iout)
    end if
  end subroutine prt_ot_saveflow

  !> @brief Print intercell flows
  subroutine prt_ot_printflow(this, ibudfl, flowja)
    ! modules
    use TdisModule, only: kper, kstp
    use ConstantsModule, only: LENBIGLINE
    ! dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: ibudfl
    real(DP), intent(inout), dimension(:) :: flowja
    ! local
    character(len=LENBIGLINE) :: line
    character(len=30) :: tempstr
    integer(I4B) :: n, ipos, m
    real(DP) :: qnm
    ! formats
    character(len=*), parameter :: fmtiprflow = &
                "(/,4x,'CALCULATED INTERCELL FLOW &
                &FOR PERIOD ', i0, ' STEP ', i0)"

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
  end subroutine prt_ot_printflow

  !> @brief Print dependent variables
  subroutine prt_ot_dv(this, idvsave, idvprint, ipflag)
    ! dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: idvsave
    integer(I4B), intent(in) :: idvprint
    integer(I4B), intent(inout) :: ipflag
    ! local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! Print advanced package dependent variables
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_dv(idvsave, idvprint)
    end do

    ! save head and print head
    call this%oc%oc_ot(ipflag)
  end subroutine prt_ot_dv

  !> @brief Print budget summary
  subroutine prt_ot_bdsummary(this, ibudfl, ipflag)
    ! modules
    use TdisModule, only: kstp, kper, totim, delt
    ! dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: ibudfl
    integer(I4B), intent(inout) :: ipflag
    ! local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    ! Package budget summary
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_ot_bdsummary(kstp, kper, this%iout, ibudfl)
    end do

    ! model budget summary
    call this%budget%finalize_step(delt)
    if (ibudfl /= 0) then
      ipflag = 1
      ! model budget summary
      call this%budget%budget_ot(kstp, kper, this%iout)
    end if

    ! Write to budget csv
    call this%budget%writecsv(totim)
  end subroutine prt_ot_bdsummary

  !> @brief Deallocate
  subroutine prt_da(this)
    ! modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    use MethodPoolModule, only: destroy_method_pool
    use MethodCellPoolModule, only: destroy_method_cell_pool
    use MethodSubcellPoolModule, only: destroy_method_subcell_pool
    ! dummy
    class(PrtModelType) :: this
    ! local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj

    ! Deallocate idm memory
    call memorystore_remove(this%name, 'NAM', idm_context)
    call memorystore_remove(component=this%name, context=idm_context)

    ! Internal packages
    call this%dis%dis_da()
    call this%fmi%fmi_da()
    call this%mip%mip_da()
    call this%budget%budget_da()
    call this%oc%oc_da()
    deallocate (this%dis)
    deallocate (this%fmi)
    deallocate (this%mip)
    deallocate (this%budget)
    deallocate (this%oc)

    ! Method objects
    call destroy_method_subcell_pool()
    call destroy_method_cell_pool()
    call destroy_method_pool()

    ! Boundary packages
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_da()
      deallocate (packobj)
    end do

    ! Scalars
    call mem_deallocate(this%infmi)
    call mem_deallocate(this%inmip)
    call mem_deallocate(this%inadv)
    call mem_deallocate(this%indsp)
    call mem_deallocate(this%inssm)
    call mem_deallocate(this%inmst)
    call mem_deallocate(this%inmvt)
    call mem_deallocate(this%inoc)

    ! Arrays
    call mem_deallocate(this%masssto)
    call mem_deallocate(this%massstoold)
    call mem_deallocate(this%ratesto)

    deallocate (this%trackctl)

    call this%NumericalModelType%model_da()
  end subroutine prt_da

  !> @brief Allocate memory for scalars
  subroutine allocate_scalars(this, modelname)
    ! dummy
    class(PrtModelType) :: this
    character(len=*), intent(in) :: modelname

    ! allocate members from parent class
    call this%NumericalModelType%allocate_scalars(modelname)

    ! allocate members that are part of model class
    call mem_allocate(this%infmi, 'INFMI', this%memoryPath)
    call mem_allocate(this%inmip, 'INMIP', this%memoryPath)
    call mem_allocate(this%inmvt, 'INMVT', this%memoryPath)
    call mem_allocate(this%inmst, 'INMST', this%memoryPath)
    call mem_allocate(this%inadv, 'INADV', this%memoryPath)
    call mem_allocate(this%indsp, 'INDSP', this%memoryPath)
    call mem_allocate(this%inssm, 'INSSM', this%memoryPath)
    call mem_allocate(this%inoc, 'INOC ', this%memoryPath)

    this%infmi = 0
    this%inmip = 0
    this%inmvt = 0
    this%inmst = 0
    this%inadv = 0
    this%indsp = 0
    this%inssm = 0
    this%inoc = 0
  end subroutine allocate_scalars

  !> @brief Allocate arrays
  subroutine allocate_arrays(this)
    use MemoryManagerModule, only: mem_allocate
    class(PrtModelType) :: this
    integer(I4B) :: n

    ! Allocate arrays in parent type
    this%nja = this%dis%nja
    call this%NumericalModelType%allocate_arrays()

    ! Allocate and initialize arrays
    call mem_allocate(this%masssto, this%dis%nodes, &
                      'MASSSTO', this%memoryPath)
    call mem_allocate(this%massstoold, this%dis%nodes, &
                      'MASSSTOOLD', this%memoryPath)
    call mem_allocate(this%ratesto, this%dis%nodes, &
                      'RATESTO', this%memoryPath)
    ! explicit model, so these must be manually allocated
    call mem_allocate(this%x, this%dis%nodes, 'X', this%memoryPath)
    call mem_allocate(this%rhs, this%dis%nodes, 'RHS', this%memoryPath)
    call mem_allocate(this%ibound, this%dis%nodes, 'IBOUND', this%memoryPath)
    do n = 1, this%dis%nodes
      this%masssto(n) = DZERO
      this%massstoold(n) = DZERO
      this%ratesto(n) = DZERO
      this%x(n) = DZERO
      this%rhs(n) = DZERO
      this%ibound(n) = 1
    end do
  end subroutine allocate_arrays

  !> @brief Create boundary condition packages for this model
  subroutine package_create(this, filtyp, ipakid, ipaknum, pakname, mempath, &
                            inunit, iout)
    ! modules
    use ConstantsModule, only: LINELENGTH
    use PrtPrpModule, only: prp_create
    use ApiModule, only: api_create
    ! dummy
    class(PrtModelType) :: this
    character(len=*), intent(in) :: filtyp
    character(len=LINELENGTH) :: errmsg
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
    case ('PRP6')
      call prp_create(packobj, ipakid, ipaknum, inunit, iout, &
                      this%name, pakname, this%fmi)
    case ('API6')
      call api_create(packobj, ipakid, ipaknum, inunit, iout, &
                      this%name, pakname)
    case default
      write (errmsg, *) 'Invalid package type: ', filtyp
      call store_error(errmsg, terminate=.TRUE.)
    end select

    ! Packages is the bndlist that is associated with the parent model
    ! The following statement puts a pointer to this package in the ipakid
    ! position of packages.
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
  subroutine ftype_check(this, indis)
    ! dummy
    class(PrtModelType) :: this
    integer(I4B), intent(in) :: indis
    ! local
    character(len=LINELENGTH) :: errmsg

    ! Check for DIS(u) and MIP. Stop if not present.
    if (indis == 0) then
      write (errmsg, '(1x,a)') &
        'Discretization (DIS6, DISV6, or DISU6) package not specified.'
      call store_error(errmsg)
    end if
    if (this%inmip == 0) then
      write (errmsg, '(1x,a)') &
        'Model input (MIP6) package not specified.'
      call store_error(errmsg)
    end if

    if (count_errors() > 0) then
      write (errmsg, '(1x,a)') 'One or more required package(s) not specified.'
      call store_error(errmsg)
      call store_error_filename(this%filename)
    end if
  end subroutine ftype_check

  !> @brief Solve the model
  subroutine prt_solve(this)
    ! modules
    use TdisModule, only: kper, kstp, totimc, delt, endofsimulation
    use PrtPrpModule, only: PrtPrpType
    use ParticleModule, only: ACTIVE, TERM_UNRELEASED, TERM_TIMEOUT
    ! dummy variables
    class(PrtModelType) :: this
    ! local variables
    integer(I4B) :: np, ip
    class(BndType), pointer :: packobj
    type(ParticleType), pointer :: particle
    real(DP) :: tmax
    integer(I4B) :: iprp

    call create_particle(particle)

    ! Apply tracking solution to PRP packages
    iprp = 0
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      select type (packobj)
      type is (PrtPrpType)
        ! Update PRP index
        iprp = iprp + 1

        ! Initialize PRP track files
        if (kper == 1 .and. kstp == 1) then
          if (packobj%itrkout > 0) then
            call this%trackctl%init_track_file( &
              packobj%itrkout, &
              iprp=iprp)
          end if
          if (packobj%itrkcsv > 0) then
            call this%trackctl%init_track_file( &
              packobj%itrkcsv, &
              csv=.true., &
              iprp=iprp)
          end if
        end if

        ! Track particles
        do np = 1, packobj%nparticles
          ! Load particle from storage
          call packobj%particles%get(particle, this%id, iprp, np)

          ! If particle is permanently unreleased, cycle.
          ! Save a termination record if we haven't yet.
          ! TODO: when we have generic dynamic vectors,
          ! consider terminating permanently unreleased
          ! in PRP instead of here. For now, status -8
          ! indicates the permanently unreleased event
          ! is not yet recorded, status 8 it has been.
          if (particle%istatus == (-1 * TERM_UNRELEASED)) then
            particle%istatus = TERM_UNRELEASED
            call this%method%save(particle, reason=3)
            call packobj%particles%put(particle, np)
          end if

          ! Skip terminated particles
          if (particle%istatus > ACTIVE) cycle

          ! If particle was released this time step, record release
          particle%istatus = ACTIVE
          if (particle%trelease >= totimc) &
            call this%method%save(particle, reason=0)

          ! Maximum time is the end of the time step or the particle
          ! stop time, whichever comes first, unless it's the final
          ! time step and the extend option is on, in which case
          ! it's just the particle stop time.
          if (endofsimulation .and. particle%iextend > 0) then
            tmax = particle%tstop
          else
            tmax = min(totimc + delt, particle%tstop)
          end if

          ! Get and apply the tracking method
          call this%method%apply(particle, tmax)

          ! Reset previous cell and zone numbers
          particle%icp = 0
          particle%izp = 0

          ! If the particle timed out, terminate it.
          ! "Timeout" means it remains active, but
          !   - it reached its stop time, or
          !   - the simulation is over and no extending.
          ! We can't detect timeout within the tracking
          ! method because the method just receives the
          ! maximum time with no context on what it is.
          ! TODO maybe think about changing that?
          if (particle%istatus <= ACTIVE .and. &
              (particle%ttrack == particle%tstop .or. &
               (endofsimulation .and. particle%iextend == 0))) then
            particle%istatus = TERM_TIMEOUT
            call this%method%save(particle, reason=3)
          end if

          ! Update particle storage
          call packobj%particles%put(particle, np)
        end do
      end select
    end do

    deallocate (particle)
  end subroutine prt_solve

  !> @brief Source package info and begin to process
  subroutine create_bndpkgs(this, bndpkgs, pkgtypes, pkgnames, &
                            mempaths, inunits)
    ! modules
    use ConstantsModule, only: LINELENGTH, LENPACKAGENAME
    use CharacterStringModule, only: CharacterStringType
    ! dummy
    class(PrtModelType) :: this
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
      !
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
      !
      ! cleanup
      deallocate (bndpkgs)
    end if

  end subroutine create_bndpkgs

  !> @brief Source package info and begin to process
  subroutine create_packages(this)
    ! modules
    use ConstantsModule, only: LINELENGTH, LENPACKAGENAME
    use CharacterStringModule, only: CharacterStringType
    use ArrayHandlersModule, only: expandarray
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path
    use SimVariablesModule, only: idm_context
    use BudgetModule, only: budget_cr
    use MethodPoolModule, only: create_method_pool
    use MethodCellPoolModule, only: create_method_cell_pool
    use MethodSubcellPoolModule, only: create_method_subcell_pool
    use PrtMipModule, only: mip_cr
    use PrtFmiModule, only: fmi_cr
    use PrtOcModule, only: oc_cr
    ! dummy
    class(PrtModelType) :: this
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
    character(len=LENFTYPE) :: pkgtype
    character(len=LENPACKAGENAME) :: pkgname
    character(len=LENMEMPATH) :: mempath
    integer(I4B), pointer :: inunit
    integer(I4B), dimension(:), allocatable :: bndpkgs
    integer(I4B) :: n
    integer(I4B) :: indis = 0 ! DIS enabled flag
    character(len=LENMEMPATH) :: mempathmip = ''

    ! set input memory paths, input/model and input/model/namfile
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

      ! create dis package first as it is a prerequisite for other packages
      select case (pkgtype)
      case ('DIS6')
        indis = 1
        call dis_cr(this%dis, this%name, mempath, indis, this%iout)
      case ('DISV6')
        indis = 1
        call disv_cr(this%dis, this%name, mempath, indis, this%iout)
      case ('DISU6')
        indis = 1
        call disu_cr(this%dis, this%name, mempath, indis, this%iout)
      case ('MIP6')
        this%inmip = 1
        mempathmip = mempath
      case ('FMI6')
        this%infmi = inunit
      case ('OC6')
        this%inoc = inunit
      case ('PRP6')
        call expandarray(bndpkgs)
        bndpkgs(size(bndpkgs)) = n
      case default
        call pstop(1, "Unrecognized package type: "//pkgtype)
      end select
    end do

    ! Create budget manager
    call budget_cr(this%budget, this%name)

    ! Create tracking method pools
    call create_method_pool()
    call create_method_cell_pool()
    call create_method_subcell_pool()

    ! Create packages that are tied directly to model
    call mip_cr(this%mip, this%name, mempathmip, this%inmip, this%iout, this%dis)
    call fmi_cr(this%fmi, this%name, this%infmi, this%iout)
    call oc_cr(this%oc, this%name, this%inoc, this%iout)

    ! Check to make sure that required ftype's have been specified
    call this%ftype_check(indis)

    ! Create boundary packages
    call this%create_bndpkgs(bndpkgs, pkgtypes, pkgnames, mempaths, inunits)
  end subroutine create_packages

  !> @brief Write model namfile options to list file
  subroutine log_namfile_options(this, found)
    use GwfNamInputModule, only: GwfNamParamFoundType
    class(PrtModelType) :: this
    type(GwfNamParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'NAMEFILE OPTIONS:'

    if (found%newton) then
      write (this%iout, '(4x,a)') &
        'NEWTON-RAPHSON method enabled for the model.'
      if (found%under_relaxation) then
        write (this%iout, '(4x,a,a)') &
          'NEWTON-RAPHSON UNDER-RELAXATION based on the bottom ', &
          'elevation of the model will be applied to the model.'
      end if
    end if

    if (found%print_input) then
      write (this%iout, '(4x,a)') 'STRESS PACKAGE INPUT WILL BE PRINTED '// &
        'FOR ALL MODEL STRESS PACKAGES'
    end if

    if (found%print_flows) then
      write (this%iout, '(4x,a)') 'PACKAGE FLOWS WILL BE PRINTED '// &
        'FOR ALL MODEL PACKAGES'
    end if

    if (found%save_flows) then
      write (this%iout, '(4x,a)') &
        'FLOWS WILL BE SAVED TO BUDGET FILE SPECIFIED IN OUTPUT CONTROL'
    end if

    write (this%iout, '(1x,a)') 'END NAMEFILE OPTIONS:'
  end subroutine log_namfile_options

end module PrtModule
