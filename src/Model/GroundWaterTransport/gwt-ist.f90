!> -- @ brief Immobile Storage and Transfer (IST) Module
!!
!!    The GwtIstModule is contains the GwtIstType, which is the
!!    derived type responsible for adding the effects of an
!!    immobile domain.  In addition to representing transfer
!!    of mass between the mobile and immobile domain, the IST
!!    Package also represents the following processes within
!!    the immobile domain
!!      1. Changes in dissolved solute mass
!!      2. Decay of dissolved solute mass
!!      3. Sorption
!!      4. Decay of sorbed solute mass
!<
module GwtIstModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, DHALF, LENFTYPE, &
                             LENPACKAGENAME, &
                             LENBUDTXT, DHNOFLO
  use BndModule, only: BndType
  use BudgetModule, only: BudgetType
  use TspFmiModule, only: TspFmiType
  use GwtMstModule, only: GwtMstType, get_zero_order_decay, &
                          get_freundlich_conc, get_freundlich_derivative, &
                          get_langmuir_conc, get_langmuir_derivative, &
                          get_freundlich_kd, get_langmuir_kd
  use OutputControlDataModule, only: OutputControlDataType
  use MatrixBaseModule
  !
  implicit none
  !
  private
  public :: ist_create
  !
  character(len=LENFTYPE) :: ftype = 'IST'
  character(len=LENPACKAGENAME) :: text = ' IMMOBILE DOMAIN'
  integer(I4B), parameter :: NBDITEMS = 5
  character(len=LENBUDTXT), dimension(NBDITEMS) :: budtxt
  data budtxt/' STORAGE-AQUEOUS', '  STORAGE-SORBED', &
    '   DECAY-AQUEOUS', '    DECAY-SORBED', &
    '   MOBILE-DOMAIN'/

  !> @ brief Immobile storage and transfer
  !!
  !!  Data and methods for handling the effects of an
  !!  immobile domain.  Note that there can be as many of these
  !!  domains as necessary.  Each immobile domain represents
  !!  changes in immobile solute storage, decay of dissolved
  !!  immobile solute mass, sorption within the immobile domain,
  !!  and decay of immobile domain sorbed mass.  The immobile
  !!  domain also includes exchange with the mobile domain.
  !<
  type, extends(BndType) :: GwtIstType

    type(TspFmiType), pointer :: fmi => null() !< flow model interface
    type(GwtMstType), pointer :: mst => null() !< mobile storage and transfer
    type(BudgetType), pointer :: budget => null() !< budget
    type(OutputControlDataType), pointer :: ocd => null() !< output control data

    integer(I4B), pointer :: icimout => null() !< unit number for binary cim output
    integer(I4B), pointer :: ibudgetout => null() !< binary budget output file
    integer(I4B), pointer :: ibudcsv => null() !< unit number for csv budget output file
    integer(I4B), pointer :: ioutsorbate => null() !< unit number for sorbate concentration output
    integer(I4B), pointer :: idcy => null() !< order of decay rate (0:none, 1:first, 2:zero)
    integer(I4B), pointer :: isrb => null() !< sorption active flag (0:off, 1:on); only linear is supported in ist
    integer(I4B), pointer :: kiter => null() !< picard iteration counter
    real(DP), pointer, contiguous :: cim(:) => null() !< concentration for immobile domain
    real(DP), pointer, contiguous :: cimnew(:) => null() !< immobile concentration at end of current time step
    real(DP), pointer, contiguous :: cimold(:) => null() !< immobile concentration at end of last time step
    real(DP), pointer, contiguous :: cimsrb(:) => null() !< sorbate concentration in immobile domain
    real(DP), pointer, contiguous :: zetaim(:) => null() !< mass transfer rate to immobile domain
    real(DP), pointer, contiguous :: porosity(:) => null() !< immobile domain porosity defined as volume of immobile voids per volume of immobile domain
    real(DP), pointer, contiguous :: volfrac(:) => null() !< volume fraction of the immobile domain defined as volume of immobile domain per aquifer volume
    real(DP), pointer, contiguous :: bulk_density(:) => null() !< bulk density of immobile domain defined as mass of solids in immobile domain per volume of immobile domain
    real(DP), pointer, contiguous :: distcoef(:) => null() !< distribution coefficient
    real(DP), pointer, contiguous :: sp2(:) => null() !< second sorption parameter
    real(DP), pointer, contiguous :: decay(:) => null() !< first or zero order rate constant for liquid
    real(DP), pointer, contiguous :: decaylast(:) => null() !< decay rate used for last iteration (needed for zero order decay)
    real(DP), pointer, contiguous :: decayslast(:) => null() !< sorbed decay rate used for last iteration (needed for zero order decay)
    real(DP), pointer, contiguous :: decay_sorbed(:) => null() !< first or zero order rate constant for sorbed mass
    real(DP), pointer, contiguous :: strg(:) => null() !< mass transfer rate
    real(DP) :: budterm(2, NBDITEMS) !< immmobile domain mass summaries

  contains

    procedure :: bnd_ar => ist_ar
    procedure :: bnd_rp => ist_rp
    procedure :: bnd_ad => ist_ad
    procedure :: bnd_fc => ist_fc
    procedure :: bnd_cq => ist_cq
    procedure :: bnd_bd => ist_bd
    procedure :: bnd_ot_model_flows => ist_ot_model_flows
    procedure :: bnd_ot_dv => ist_ot_dv
    procedure :: output_immobile_concentration
    procedure :: output_immobile_sorbate_concentration
    procedure :: bnd_ot_bdsummary => ist_ot_bdsummary
    procedure :: bnd_da => ist_da
    procedure :: allocate_scalars
    procedure :: read_dimensions => ist_read_dimensions
    procedure :: source_data
    procedure :: log_data
    procedure :: read_options
    procedure :: source_options
    procedure :: log_options
    procedure :: get_thetaim
    procedure :: ist_calc_csrb
    procedure, private :: ist_allocate_arrays
    procedure, private :: read_data

  end type GwtIstType

contains

  !> @ brief Create a new package object
  !!
  !!  Create a new IST object
  !!
  !<
  subroutine ist_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        mempath, fmi, mst)
    ! -- dummy
    class(BndType), pointer :: packobj !< BndType pointer that will point to new IST Package
    integer(I4B), intent(in) :: id !< name of the model
    integer(I4B), intent(in) :: ibcnum !< consecutive package number
    integer(I4B), intent(in) :: inunit !< unit number of package input file
    integer(I4B), intent(in) :: iout !< unit number of model listing file
    character(len=*), intent(in) :: namemodel !< name of the model
    character(len=*), intent(in) :: pakname !< name of the package
    character(len=*), intent(in) :: mempath
    type(TspFmiType), pointer, intent(in) :: fmi
    type(GwtMstType), pointer, intent(in) :: mst
    ! -- local
    type(GwtIstType), pointer :: istobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (istobj)
    packobj => istobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text
    !
    ! -- allocate scalars
    call packobj%allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()
    !
    ! -- store values
    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 1
    packobj%iscloc = 1
    !
    ! -- Point IST specific variables
    istobj%fmi => fmi
    istobj%mst => mst
  end subroutine ist_create

  !> @ brief Allocate and read method for package
  !!
  !!  Method to allocate and read static data for the package.
  !!
  !<
  subroutine ist_ar(this)
    ! -- modules
    use SimModule, only: store_error, count_errors
    use BudgetModule, only: budget_cr
    ! -- dummy
    class(GwtIstType), intent(inout) :: this !< GwtIstType object
    ! -- local
    integer(I4B) :: n
    !
    ! -- Allocate arrays
    call this%ist_allocate_arrays()
    !
    ! -- Now that arrays are allocated, check in the cimnew array to
    !    the output control manager for subsequent printing/saving
    call this%ocd%init_dbl('CIM', this%cimnew, this%dis, 'PRINT LAST ', &
                           'COLUMNS 10 WIDTH 11 DIGITS 4 GENERAL ', &
                           this%iout, DHNOFLO)
    !
    ! -- read the data block
    call this%read_data()
    !
    ! -- set cimnew to the cim start values read from input
    do n = 1, this%dis%nodes
      this%cimnew(n) = this%cim(n)
    end do
    !
    ! -- add volfrac to the volfracim accumulator in mst package
    call this%mst%addto_volfracim(this%volfrac)
    !
    ! -- setup the immobile domain budget
    call budget_cr(this%budget, this%memoryPath)
    call this%budget%budget_df(NBDITEMS, 'MASS', 'M', bdzone=this%packName)
    call this%budget%set_ibudcsv(this%ibudcsv)
    !
    ! -- Perform a check to ensure that sorption and decay are set
    !    consistently between the MST and IST packages.
    if (this%idcy /= this%mst%idcy) then
      call store_error('DECAY must be activated consistently between the &
        &MST and IST Packages.  Activate or deactivate DECAY for both &
        &Packages.')
    end if
    if (this%isrb /= this%mst%isrb) then
      call store_error('SORPTION must be activated consistently between the &
        &MST and IST Packages.  Activate or deactivate SORPTION for both &
        &Packages.  If activated, the same type of sorption (LINEAR, &
        &FREUNDLICH, or LANGMUIR) must be specified in the options block of &
        &both the MST and IST Packages.')
    end if
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine ist_ar

  !> @ brief Read and prepare method for package
  !!
  !!  Method to read and prepare package data
  !!
  !<
  subroutine ist_rp(this)
    ! -- dummy
    class(GwtIstType), intent(inout) :: this !< GwtIstType object
    ! -- local
    ! -- format
  end subroutine ist_rp

  !> @ brief Advance the ist package
  !!
  !!  Advance the IST Package and handle the adaptive time stepping
  !!  feature by copying from new to old or old to new accordingly
  !!
  !<
  subroutine ist_ad(this)
    ! -- modules
    use SimVariablesModule, only: iFailedStepRetry
    ! -- dummy variables
    class(GwtIstType) :: this !< GwtIstType object
    ! -- local variables
    integer(I4B) :: n
    !
    ! -- Call parent advance
    call this%BndType%bnd_ad()
    !
    ! -- set independent kiter counter to zero
    this%kiter = 0
    !
    ! -- copy cimnew into cimold or vice versa if this is a repeat of
    !    a failed time step
    if (iFailedStepRetry == 0) then
      do n = 1, this%dis%nodes
        this%cimold(n) = this%cimnew(n)
      end do
    else
      do n = 1, this%dis%nodes
        this%cimnew(n) = this%cimold(n)
      end do
    end if
  end subroutine ist_ad

  !> @ brief Fill coefficient method for package
  !<
  subroutine ist_fc(this, rhs, ia, idxglo, matrix_sln)
    ! modules
    use TdisModule, only: delt
    ! dummy
    class(GwtIstType) :: this !< GwtIstType object
    real(DP), dimension(:), intent(inout) :: rhs !< right-hand side vector for model
    integer(I4B), dimension(:), intent(in) :: ia !< solution CRS row pointers
    integer(I4B), dimension(:), intent(in) :: idxglo !< mapping vector for model (local) to solution (global)
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    ! local
    integer(I4B) :: n, idiag
    real(DP) :: tled
    real(DP) :: hhcof, rrhs
    real(DP) :: swt, swtpdt
    real(DP) :: vcell
    real(DP) :: thetaim
    real(DP) :: zetaim
    real(DP) :: kdnew
    real(DP) :: kdold
    real(DP) :: volfracim
    real(DP) :: rhobim
    real(DP) :: lambda1im
    real(DP) :: lambda2im
    real(DP) :: gamma1im
    real(DP) :: gamma2im
    real(DP) :: cimold
    real(DP) :: f
    real(DP) :: cimsrbold
    real(DP) :: cimsrbnew
    real(DP), dimension(9) :: ddterm

    ! set variables
    tled = DONE / delt
    this%kiter = this%kiter + 1

    ! loop through each node and calculate immobile domain contribution
    ! to hcof and rhs
    do n = 1, this%dis%nodes

      ! skip if transport inactive
      if (this%ibound(n) <= 0) cycle

      ! calculate new and old water volumes
      vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
      swtpdt = this%fmi%gwfsat(n)
      swt = this%fmi%gwfsatold(n, delt)
      thetaim = this%get_thetaim(n)
      idiag = ia(n)

      ! set exchange coefficient
      zetaim = this%zetaim(n)

      ! Add dual domain mass transfer contributions to rhs and hcof
      ! dcimsrbdc = DZERO
      kdnew = DZERO
      kdold = DZERO
      volfracim = DZERO
      rhobim = DZERO
      lambda1im = DZERO
      lambda2im = DZERO
      gamma1im = DZERO
      gamma2im = DZERO

      ! set variables for decay of aqueous solute
      if (this%idcy == 1) lambda1im = this%decay(n)
      if (this%idcy == 2) then
        gamma1im = get_zero_order_decay(this%decay(n), this%decaylast(n), &
                                        this%kiter, this%cimold(n), &
                                        this%cimnew(n), delt)
        this%decaylast(n) = gamma1im
      end if

      ! setup sorption variables
      if (this%isrb > 0) then

        ! initialize sorption variables
        volfracim = this%volfrac(n)
        rhobim = this%bulk_density(n)

        ! set isotherm dependent sorption variables
        select case (this%isrb)
        case (1)
          ! linear
          kdnew = this%distcoef(n)
          kdold = this%distcoef(n)
          cimsrbnew = this%cimnew(n) * kdnew
          cimsrbold = this%cimold(n) * kdold
        case (2)
          ! freundlich
          kdnew = get_freundlich_kd(this%cimnew(n), this%distcoef(n), &
                                    this%sp2(n))
          kdold = get_freundlich_kd(this%cimold(n), this%distcoef(n), &
                                    this%sp2(n))
          cimsrbnew = get_freundlich_conc(this%cimnew(n), this%distcoef(n), &
                                          this%sp2(n))
          cimsrbold = get_freundlich_conc(this%cimold(n), this%distcoef(n), &
                                          this%sp2(n))
        case (3)
          ! langmuir
          kdnew = get_langmuir_kd(this%cimnew(n), this%distcoef(n), &
                                  this%sp2(n))
          kdold = get_langmuir_kd(this%cimold(n), this%distcoef(n), &
                                  this%sp2(n))
          cimsrbnew = get_langmuir_conc(this%cimnew(n), this%distcoef(n), &
                                        this%sp2(n))
          cimsrbold = get_langmuir_conc(this%cimold(n), this%distcoef(n), &
                                        this%sp2(n))
        end select

        ! set decay of sorbed solute parameters
        if (this%idcy == 1) then
          lambda2im = this%decay_sorbed(n)
        else if (this%idcy == 2) then
          gamma2im = get_zero_order_decay(this%decay_sorbed(n), &
                                          this%decayslast(n), &
                                          this%kiter, cimsrbold, &
                                          cimsrbnew, delt)
          this%decayslast(n) = gamma2im
        end if
      end if

      ! calculate dual domain terms and get the hcof and rhs contributions
      call get_ddterm(thetaim, vcell, delt, swtpdt, &
                      volfracim, rhobim, kdnew, kdold, lambda1im, lambda2im, &
                      gamma1im, gamma2im, zetaim, ddterm, f)
      cimold = this%cimold(n)
      call get_hcofrhs(ddterm, f, cimold, hhcof, rrhs)

      ! update solution accumulators
      call matrix_sln%add_value_pos(idxglo(idiag), hhcof)
      rhs(n) = rhs(n) + rrhs

    end do

  end subroutine ist_fc

  !> @ brief Calculate package flows.
  !<
  subroutine ist_cq(this, x, flowja, iadv)
    ! modules
    use TdisModule, only: delt
    use ConstantsModule, only: DZERO
    ! dummy
    class(GwtIstType), intent(inout) :: this !< GwtIstType object
    real(DP), dimension(:), intent(in) :: x !< current dependent-variable value
    real(DP), dimension(:), contiguous, intent(inout) :: flowja !< flow between two connected control volumes
    integer(I4B), optional, intent(in) :: iadv !< flag that indicates if this is an advance package
    ! local
    integer(I4B) :: idiag
    integer(I4B) :: n
    real(DP) :: rate
    real(DP) :: swt, swtpdt
    real(DP) :: hhcof, rrhs
    real(DP) :: vcell
    real(DP) :: thetaim
    real(DP) :: zetaim
    real(DP) :: kdnew
    real(DP) :: kdold
    real(DP) :: volfracim
    real(DP) :: rhobim
    real(DP) :: lambda1im
    real(DP) :: lambda2im
    real(DP) :: gamma1im
    real(DP) :: gamma2im
    real(DP) :: cimnew
    real(DP) :: cimold
    real(DP) :: f
    real(DP) :: cimsrbold
    real(DP) :: cimsrbnew
    real(DP), dimension(9) :: ddterm
    ! formats

    ! initialize
    this%budterm(:, :) = DZERO

    ! Calculate immobile domain transfer rate
    do n = 1, this%dis%nodes

      ! skip if transport inactive
      rate = DZERO
      cimnew = DZERO
      if (this%ibound(n) > 0) then

        ! calculate new and old water volumes
        vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
        swtpdt = this%fmi%gwfsat(n)
        swt = this%fmi%gwfsatold(n, delt)
        thetaim = this%get_thetaim(n)

        ! set exchange coefficient
        zetaim = this%zetaim(n)

        ! Calculate exchange with immobile domain
        rate = DZERO
        hhcof = DZERO
        rrhs = DZERO
        kdnew = DZERO
        kdold = DZERO
        volfracim = DZERO
        rhobim = DZERO
        lambda1im = DZERO
        lambda2im = DZERO
        gamma1im = DZERO
        gamma2im = DZERO
        if (this%idcy == 1) lambda1im = this%decay(n)
        if (this%idcy == 2) then
          gamma1im = get_zero_order_decay(this%decay(n), this%decaylast(n), 0, &
                                          this%cimold(n), this%cimnew(n), delt)
        end if

        ! setup sorption variables
        if (this%isrb > 0) then

          ! initialize sorption variables
          volfracim = this%volfrac(n)
          rhobim = this%bulk_density(n)

          ! set isotherm dependent sorption variables
          select case (this%isrb)
          case (1)
            ! linear
            kdnew = this%distcoef(n)
            kdold = this%distcoef(n)
            cimsrbnew = this%cimnew(n) * kdnew
            cimsrbold = this%cimold(n) * kdold
          case (2)
            ! freundlich
            kdnew = get_freundlich_kd(this%cimnew(n), this%distcoef(n), &
                                      this%sp2(n))
            kdold = get_freundlich_kd(this%cimold(n), this%distcoef(n), &
                                      this%sp2(n))
            cimsrbnew = get_freundlich_conc(this%cimnew(n), this%distcoef(n), &
                                            this%sp2(n))
            cimsrbold = get_freundlich_conc(this%cimold(n), this%distcoef(n), &
                                            this%sp2(n))
          case (3)
            ! langmuir
            kdnew = get_langmuir_kd(this%cimnew(n), this%distcoef(n), &
                                    this%sp2(n))
            kdold = get_langmuir_kd(this%cimold(n), this%distcoef(n), &
                                    this%sp2(n))
            cimsrbnew = get_langmuir_conc(this%cimnew(n), this%distcoef(n), &
                                          this%sp2(n))
            cimsrbold = get_langmuir_conc(this%cimold(n), this%distcoef(n), &
                                          this%sp2(n))
          end select

          ! set decay of sorbed solute parameters
          if (this%idcy == 1) then
            lambda2im = this%decay_sorbed(n)
          else if (this%idcy == 2) then
            gamma2im = get_zero_order_decay(this%decay_sorbed(n), &
                                            this%decayslast(n), &
                                            0, cimsrbold, &
                                            cimsrbnew, delt)
          end if
        end if

        ! calculate the terms and then get the hcof and rhs contributions
        call get_ddterm(thetaim, vcell, delt, swtpdt, &
                        volfracim, rhobim, kdnew, kdold, lambda1im, lambda2im, &
                        gamma1im, gamma2im, zetaim, ddterm, f)
        cimold = this%cimold(n)
        call get_hcofrhs(ddterm, f, cimold, hhcof, rrhs)

        ! calculate rate from hcof and rhs
        rate = hhcof * x(n) - rrhs

        ! calculate immobile domain concentration
        cimnew = get_ddconc(ddterm, f, cimold, x(n))

        ! accumulate the budget terms
        call accumulate_budterm(this%budterm, ddterm, cimnew, cimold, x(n), &
                                this%idcy)
      end if

      ! store rate and add to flowja
      this%strg(n) = rate
      idiag = this%dis%con%ia(n)
      flowja(idiag) = flowja(idiag) + rate

      ! store immobile domain concentration
      this%cimnew(n) = cimnew

    end do

    ! calculate csrb
    if (this%isrb /= 0) then
      call this%ist_calc_csrb(this%cimnew)
    end if

  end subroutine ist_cq

  !> @ brief Calculate immobile sorbed concentration
  !<
  subroutine ist_calc_csrb(this, cim)
    ! -- dummy
    class(GwtIstType) :: this !< GwtMstType object
    real(DP), intent(in), dimension(:) :: cim !< immobile domain aqueous concentration at end of this time step
    ! -- local
    integer(I4B) :: n
    real(DP) :: distcoef
    real(DP) :: csrb

    ! Calculate sorbed concentration
    do n = 1, size(cim)
      csrb = DZERO
      if (this%ibound(n) > 0) then
        distcoef = this%distcoef(n)
        if (this%isrb == 1) then
          csrb = cim(n) * distcoef
        else if (this%isrb == 2) then
          csrb = get_freundlich_conc(cim(n), distcoef, this%sp2(n))
        else if (this%isrb == 3) then
          csrb = get_langmuir_conc(cim(n), distcoef, this%sp2(n))
        end if
      end if
      this%cimsrb(n) = csrb
    end do

  end subroutine ist_calc_csrb

  !> @ brief Add package flows to model budget.
  !!
  !!  Add the flow between IST package and the model (ratin and ratout) to the
  !!  model budget.
  !!
  !<
  subroutine ist_bd(this, model_budget)
    ! -- modules
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType, rate_accumulator
    ! -- dummy
    class(GwtIstType) :: this !< GwtIstType object
    type(BudgetType), intent(inout) :: model_budget !< model budget object
    ! -- local
    real(DP) :: ratin
    real(DP) :: ratout
    integer(I4B) :: isuppress_output
    isuppress_output = 0
    call rate_accumulator(this%strg(:), ratin, ratout)
    call model_budget%addentry(ratin, ratout, delt, this%text, &
                               isuppress_output, this%packName)
  end subroutine ist_bd

  !> @ brief Output model flow terms.
  !!
  !!  Output flow terms between the IST package and model to a binary file and/or
  !!  print flows to the model listing file.
  !!
  !<
  subroutine ist_ot_model_flows(this, icbcfl, ibudfl, icbcun, imap)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy
    class(GwtIstType) :: this !< GwtIstType object
    integer(I4B), intent(in) :: icbcfl !< flag for cell-by-cell output
    integer(I4B), intent(in) :: ibudfl !< flag indication if cell-by-cell data should be saved
    integer(I4B), intent(in) :: icbcun !< unit number for cell-by-cell output
    integer(I4B), dimension(:), optional, intent(in) :: imap !< mapping vector
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: ibinun
    integer(I4B) :: nbound
    integer(I4B) :: naux
    real(DP) :: rate
    real(DP), dimension(0) :: auxrow
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
    !
    ! -- If cell-by-cell flows will be saved as a list, write header.
    if (ibinun /= 0) then
      nbound = this%dis%nodes
      naux = 0
      call this%dis%record_srcdst_list_header(this%text, this%name_model, &
                                              this%name_model, this%name_model, &
                                              this%packName, naux, this%auxname, &
                                              ibinun, nbound, this%iout)
    end if
    !
    ! -- Calculate immobile domain rhs and hcof
    do n = 1, this%dis%nodes
      !
      ! -- skip if transport inactive
      rate = DZERO
      if (this%ibound(n) > 0) then
        !
        ! -- set rate from this%strg
        rate = this%strg(n)
      end if
      !
      ! -- If saving cell-by-cell flows in list, write flow
      if (ibinun /= 0) then
        call this%dis%record_mf6_list_entry(ibinun, n, n, rate, &
                                            naux, auxrow, &
                                            olconv=.TRUE., &
                                            olconv2=.TRUE.)
      end if
      !
    end do
  end subroutine ist_ot_model_flows

  !> @ brief Output dependent variables.
  !<
  subroutine ist_ot_dv(this, idvsave, idvprint)
    ! dummy variables
    class(GwtIstType) :: this !< BndType object
    integer(I4B), intent(in) :: idvsave !< flag and unit number for dependent-variable output
    integer(I4B), intent(in) :: idvprint !< flag indicating if dependent-variable should be written to the model listing file

    call this%output_immobile_concentration(idvsave, idvprint)
    call this%output_immobile_sorbate_concentration(idvsave, idvprint)

  end subroutine ist_ot_dv

  !> @ brief Output immobile domain aqueous concentration.
  !<
  subroutine output_immobile_concentration(this, idvsave, idvprint)
    ! modules
    use TdisModule, only: kstp, endofperiod
    ! dummy variables
    class(GwtIstType) :: this !< BndType object
    integer(I4B), intent(in) :: idvsave !< flag and unit number for dependent-variable output
    integer(I4B), intent(in) :: idvprint !< flag indicating if dependent-variable should be written to the model listing file
    ! local
    integer(I4B) :: ipflg
    integer(I4B) :: ibinun
    !
    ! Save cim to a binary file. ibinun is a flag where 1 indicates that
    ! cim should be written to a binary file if a binary file is open for it.
    ipflg = 0
    ibinun = 1
    if (idvsave == 0) ibinun = 0
    if (ibinun /= 0) then
      call this%ocd%ocd_ot(ipflg, kstp, endofperiod, this%iout, &
                           iprint_opt=0, isav_opt=ibinun)
    end if
    !
    ! Print immobile domain concentrations to listing file
    if (idvprint /= 0) then
      call this%ocd%ocd_ot(ipflg, kstp, endofperiod, this%iout, &
                           iprint_opt=idvprint, isav_opt=0)
    end if

  end subroutine output_immobile_concentration

  !> @ brief Output immobile domain sorbate concentration.
  !<
  subroutine output_immobile_sorbate_concentration(this, idvsave, idvprint)
    ! modules
    ! dummy
    class(GwtIstType) :: this !< BndType object
    integer(I4B), intent(in) :: idvsave !< flag and unit number for dependent-variable output
    integer(I4B), intent(in) :: idvprint !< flag indicating if dependent-variable should be written to the model listing file
    ! local
    character(len=1) :: cdatafmp = ' ', editdesc = ' '
    integer(I4B) :: ibinun
    integer(I4B) :: iprint, nvaluesp, nwidthp
    real(DP) :: dinact

    ! Save cimsrb to a binary file. ibinun is a flag where 1 indicates that
    ! cim should be written to a binary file if a binary file is open for it.
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
        call this%dis%record_array(this%cimsrb, this%iout, iprint, ibinun, &
                                   '         SORBATE', cdatafmp, nvaluesp, &
                                   nwidthp, editdesc, dinact)
      end if
    end if

  end subroutine output_immobile_sorbate_concentration

  !> @ brief Output IST package budget summary.
  !!
  !!  Output advanced boundary package budget summary. This method only needs
  !!  to be overridden for advanced packages that save budget summaries
  !!  to the model listing file.
  !!
  !<
  subroutine ist_ot_bdsummary(this, kstp, kper, iout, ibudfl)
    ! -- modules
    use TdisModule, only: delt, totim
    ! -- dummy variables
    class(GwtIstType) :: this !< GwtIstType object
    integer(I4B), intent(in) :: kstp !< time step number
    integer(I4B), intent(in) :: kper !< period number
    integer(I4B), intent(in) :: iout !< flag and unit number for the model listing file
    integer(I4B), intent(in) :: ibudfl !< flag indicating budget should be written
    ! -- local
    integer(I4B) :: isuppress_output = 0
    !
    ! -- Fill budget terms
    call this%budget%reset()
    call this%budget%addentry(this%budterm, delt, budtxt, isuppress_output)
    !
    ! -- Write budget to list file
    call this%budget%finalize_step(delt)
    if (ibudfl /= 0) then
      call this%budget%budget_ot(kstp, kper, iout)
    end if
    !
    ! -- Write budget csv
    call this%budget%writecsv(totim)
  end subroutine ist_ot_bdsummary

  !> @ brief Deallocate package memory
  !!
  !!  Deallocate package scalars and arrays.
  !!
  !<
  subroutine ist_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GwtIstType) :: this !< GwtIstType object
    !
    ! -- Deallocate arrays if package was active
    if (this%inunit > 0) then
      call mem_deallocate(this%icimout)
      call mem_deallocate(this%ibudgetout)
      call mem_deallocate(this%ibudcsv)
      call mem_deallocate(this%ioutsorbate)
      call mem_deallocate(this%idcy)
      call mem_deallocate(this%isrb)
      call mem_deallocate(this%kiter)
      call mem_deallocate(this%cim)
      call mem_deallocate(this%cimnew)
      call mem_deallocate(this%cimold)
      call mem_deallocate(this%cimsrb)
      call mem_deallocate(this%zetaim)
      call mem_deallocate(this%porosity)
      call mem_deallocate(this%volfrac)
      call mem_deallocate(this%bulk_density)
      call mem_deallocate(this%distcoef)
      call mem_deallocate(this%sp2)
      call mem_deallocate(this%decay)
      call mem_deallocate(this%decaylast)
      call mem_deallocate(this%decayslast)
      call mem_deallocate(this%decay_sorbed)
      call mem_deallocate(this%strg)
      this%fmi => null()
      this%mst => null()
    end if
    !
    ! -- Scalars
    !
    ! -- Objects
    call this%budget%budget_da()
    deallocate (this%budget)
    call this%ocd%ocd_da()
    deallocate (this%ocd)
    !
    ! -- deallocate parent
    call this%BndType%bnd_da()
  end subroutine ist_da

  !> @ brief Allocate package scalars
  !!
  !!  Allocate and initialize package scalars.
  !!
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    use OutputControlDataModule, only: ocd_cr
    ! -- dummy
    class(GwtIstType) :: this !< GwtIstType object
    ! -- local
    !
    ! -- call standard BndType allocate scalars
    call this%BndType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%icimout, 'ICIMOUT', this%memoryPath)
    call mem_allocate(this%ibudgetout, 'IBUDGETOUT', this%memoryPath)
    call mem_allocate(this%ibudcsv, 'IBUDCSV', this%memoryPath)
    call mem_allocate(this%ioutsorbate, 'IOUTSORBATE', this%memoryPath)
    call mem_allocate(this%isrb, 'ISRB', this%memoryPath)
    call mem_allocate(this%idcy, 'IDCY', this%memoryPath)
    call mem_allocate(this%kiter, 'KITER', this%memoryPath)
    !
    ! -- Initialize
    this%icimout = 0
    this%ibudgetout = 0
    this%ibudcsv = 0
    this%ioutsorbate = 0
    this%isrb = 0
    this%idcy = 0
    this%kiter = 0
    !
    ! -- Create the ocd object, which is used to manage printing and saving
    !    of the immobile domain concentrations
    call ocd_cr(this%ocd)
  end subroutine allocate_scalars

  !> @ brief Allocate package arrays
  !!
  !!  Allocate and initialize package arrays.
  !!
  !<
  subroutine ist_allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwtIstType), intent(inout) :: this !< GwtIstType object
    ! -- local
    integer(I4B) :: n
    !
    ! -- call standard BndType allocate scalars
    !    nbound and maxbound are 0 in order to keep memory footprint low
    call this%BndType%allocate_arrays()
    !
    ! -- allocate ist arrays of size nodes
    call mem_allocate(this%strg, this%dis%nodes, 'STRG', this%memoryPath)
    call mem_allocate(this%cim, this%dis%nodes, 'CIM', this%memoryPath)
    call mem_allocate(this%cimnew, this%dis%nodes, 'CIMNEW', this%memoryPath)
    call mem_allocate(this%cimold, this%dis%nodes, 'CIMOLD', this%memoryPath)
    call mem_allocate(this%porosity, this%dis%nodes, 'POROSITY', this%memoryPath)
    call mem_allocate(this%zetaim, this%dis%nodes, 'ZETAIM', this%memoryPath)
    call mem_allocate(this%volfrac, this%dis%nodes, 'VOLFRAC', this%memoryPath)
    if (this%isrb == 0) then
      call mem_allocate(this%bulk_density, 1, 'BULK_DENSITY', this%memoryPath)
      call mem_allocate(this%distcoef, 1, 'DISTCOEF', this%memoryPath)
      call mem_allocate(this%sp2, 1, 'SP2', this%memoryPath)
      call mem_allocate(this%cimsrb, 1, 'SORBATE', this%memoryPath)
    else
      call mem_allocate(this%bulk_density, this%dis%nodes, 'BULK_DENSITY', &
                        this%memoryPath)
      call mem_allocate(this%distcoef, this%dis%nodes, 'DISTCOEF', &
                        this%memoryPath)
      call mem_allocate(this%cimsrb, this%dis%nodes, 'SORBATE', &
                        this%memoryPath)
      if (this%isrb == 1) then
        call mem_allocate(this%sp2, 1, 'SP2', this%memoryPath)
      else
        call mem_allocate(this%sp2, this%dis%nodes, 'SP2', this%memoryPath)
      end if
    end if
    if (this%idcy == 0) then
      call mem_allocate(this%decay, 1, 'DECAY', this%memoryPath)
      call mem_allocate(this%decaylast, 1, 'DECAYLAST', this%memoryPath)
    else
      call mem_allocate(this%decay, this%dis%nodes, 'DECAY', this%memoryPath)
      call mem_allocate(this%decaylast, this%dis%nodes, 'DECAYLAST', &
                        this%memoryPath)
    end if
    if (this%isrb == 0 .and. this%idcy == 0) then
      call mem_allocate(this%decayslast, 1, 'DECAYSLAST', this%memoryPath)
    else
      call mem_allocate(this%decayslast, this%dis%nodes, 'DECAYSLAST', &
                        this%memoryPath)
    end if
    call mem_allocate(this%decay_sorbed, 1, 'DECAY_SORBED', this%memoryPath)
    !
    ! -- initialize
    do n = 1, this%dis%nodes
      this%strg(n) = DZERO
      this%cim(n) = DZERO
      this%cimnew(n) = DZERO
      this%cimold(n) = DZERO
      this%porosity(n) = DZERO
      this%zetaim(n) = DZERO
      this%volfrac(n) = DZERO
    end do
    do n = 1, size(this%bulk_density)
      this%bulk_density(n) = DZERO
      this%distcoef(n) = DZERO
      this%cimsrb(n) = DZERO
    end do
    do n = 1, size(this%sp2)
      this%sp2(n) = DZERO
    end do
    do n = 1, size(this%decay)
      this%decay(n) = DZERO
      this%decaylast(n) = DZERO
    end do
    do n = 1, size(this%decayslast)
      this%decayslast(n) = DZERO
    end do
    !
    ! -- Set pointers
    this%ocd%dis => this%dis
  end subroutine ist_allocate_arrays

  !> @ brief Source options for package
  !!
  !!  Method to source options for the package.
  !<
  subroutine source_options(this)
    ! -- modules
    use ConstantsModule, only: LENVARNAME, LINELENGTH, MNORMAL, LENBIGLINE
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: getunit, assign_iounit, openfile
    use MemoryManagerExtModule, only: mem_set_value
    use GwtIstInputModule, only: GwtIstParamFoundType
    ! -- dummy
    class(GwtIstType), intent(inout) :: this
    ! -- locals
    type(GwtIstParamFoundType) :: found
    character(len=LENVARNAME), dimension(3) :: sorption_method = &
      &[character(len=LENVARNAME) :: 'LINEAR', 'FREUNDLICH', 'LANGMUIR']
    character(len=LINELENGTH) :: sorbate_fname, cim6_fname, budget_fname, &
                                 budgetcsv_fname
    character(len=LENBIGLINE) :: print_format
    !
    ! -- update defaults with memory sourced values
    call mem_set_value(this%ipakcb, 'SAVE_FLOWS', this%input_mempath, &
                       found%save_flows)
    call mem_set_value(budget_fname, 'BUDGETFILE', this%input_mempath, &
                       found%budgetfile)
    call mem_set_value(budgetcsv_fname, 'BUDGETCSVFILE', this%input_mempath, &
                       found%budgetcsvfile)
    call mem_set_value(this%isrb, 'SORPTION', this%input_mempath, &
                       sorption_method, found%sorption)
    call mem_set_value(this%idcy, 'ORD1_DECAY', this%input_mempath, &
                       found%ord1_decay)
    call mem_set_value(this%idcy, 'ZERO_ORDER_DECAY', this%input_mempath, &
                       found%zero_order_decay)
    call mem_set_value(cim6_fname, 'CIMFILE', this%input_mempath, &
                       found%cimfile)
    call mem_set_value(print_format, 'PRINT_FORMAT', this%input_mempath, &
                       found%print_format)
    call mem_set_value(sorbate_fname, 'SORBATEFILE', this%input_mempath, &
                       found%sorbatefile)

    ! -- found side effects
    if (found%save_flows) this%ipakcb = -1
    if (found%budgetfile) then
      call assign_iounit(this%ibudgetout, this%inunit, "BUDGET fileout")
      call openfile(this%ibudgetout, this%iout, budget_fname, 'DATA(BINARY)', &
                    form, access, 'REPLACE', mode_opt=MNORMAL)
    end if
    if (found%budgetcsvfile) then
      call assign_iounit(this%ibudcsv, this%inunit, "BUDGETCSV fileout")
      call openfile(this%ibudcsv, this%iout, budgetcsv_fname, 'CSV', &
                    filstat_opt='REPLACE')
    end if
    !if (found%sorption) then
    !  if (this%isrb == 0) then
    !    this%isrb = 1
    !  end if
    !end if
    if (found%ord1_decay) this%idcy = 1
    if (found%zero_order_decay) this%idcy = 2
    if (found%cimfile) then
      call this%ocd%set_option('FILEOUT '//trim(cim6_fname)//" ", this%inunit, &
                               this%iout)
    end if
    if (found%print_format) then
      call this%ocd%set_option('PRINT_FORMAT '//trim(print_format)//" ", &
                               this%inunit, this%iout)
    end if
    if (found%sorbatefile) then
      this%ioutsorbate = getunit()
      call openfile(this%ioutsorbate, this%iout, sorbate_fname, &
                    'DATA(BINARY)', form, access, 'REPLACE', mode_opt=MNORMAL)
    end if
    !
    ! -- log options
    if (this%iout > 0) then
      call this%log_options(found, cim6_fname, budget_fname, &
                            budgetcsv_fname, sorbate_fname)
    end if
  end subroutine source_options

  !> @brief Write user options to list file
  !<
  subroutine log_options(this, found, cim6_fname, budget_fname, &
                         budgetcsv_fname, sorbate_fname)
    use GwtIstInputModule, only: GwtIstParamFoundType
    class(GwTIstType), intent(inout) :: this
    type(GwtIstParamFoundType), intent(in) :: found
    character(len=*), intent(in) :: cim6_fname
    character(len=*), intent(in) :: budget_fname
    character(len=*), intent(in) :: budgetcsv_fname
    character(len=*), intent(in) :: sorbate_fname
    ! -- formats
    character(len=*), parameter :: fmtisvflow = &
      "(4x,'CELL-BY-CELL FLOW INFORMATION WILL BE SAVED TO BINARY FILE &
      &WHENEVER ICBCFL IS NOT ZERO.')"
    character(len=*), parameter :: fmtlinear = &
      &"(4x,'LINEAR SORPTION IS SELECTED. ')"
    character(len=*), parameter :: fmtfreundlich = &
      &"(4x,'FREUNDLICH SORPTION IS ACTIVE. ')"
    character(len=*), parameter :: fmtlangmuir = &
      &"(4x,'LANGMUIR SORPTION IS ACTIVE. ')"
    character(len=*), parameter :: fmtidcy1 = &
      &"(4x,'FIRST-ORDER DECAY IS ACTIVE. ')"
    character(len=*), parameter :: fmtidcy2 = &
      &"(4x,'ZERO-ORDER DECAY IS ACTIVE. ')"
    character(len=*), parameter :: fmtistbin = &
      "(4x, 'IST ', 1x, a, 1x, ' WILL BE SAVED TO FILE: ', a, &
      &/4x, 'OPENED ON UNIT: ', I0)"

    write (this%iout, '(1x,a)') 'PROCESSING IMMOBILE STORAGE AND TRANSFER &
                                &OPTIONS'
    if (found%save_flows) then
      write (this%iout, fmtisvflow)
    end if
    if (found%budgetfile) then
      write (this%iout, fmtistbin) 'BUDGET', trim(adjustl(budget_fname)), &
        this%ibudgetout
    end if
    if (found%budgetcsvfile) then
      write (this%iout, fmtistbin) 'BUDGET CSV', trim(adjustl(budgetcsv_fname)), &
        this%ibudcsv
    end if
    if (found%sorption) then
      select case (this%isrb)
      case (1)
        write (this%iout, fmtlinear)
      case (2)
        write (this%iout, fmtfreundlich)
      case (3)
        write (this%iout, fmtlangmuir)
      end select
    end if
    if (found%ord1_decay) then
      write (this%iout, fmtidcy1)
    end if
    if (found%zero_order_decay) then
      write (this%iout, fmtidcy2)
    end if
    if (found%sorbatefile) then
      write (this%iout, fmtistbin) &
        'SORBATE', sorbate_fname, this%ioutsorbate
    end if
    write (this%iout, '(1x,a)') 'END OF IMMOBILE STORAGE AND TRANSFER &
                                &OPTIONS'
  end subroutine log_options

  !> @ brief Read options for package
  !!
  !!  Read options for boundary packages.
  !!
  !<
  subroutine read_options(this)
    ! -- modules
    ! -- dummy
    class(GwtIstType), intent(inout) :: this !< GwtIstType object

    ! -- source options
    call this%source_options()
  end subroutine read_options

  !> @ brief Source data for package
  !!
  !!  Method to source data for the package.
  !<
  subroutine source_data(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use SimVariablesModule, only: errmsg, warnmsg
    use SimModule, only: count_errors, store_error, store_warning
    use MemoryManagerModule, only: get_isize, mem_reallocate
    use MemoryManagerExtModule, only: mem_set_value
    use GwtIstInputModule, only: GwtIstParamFoundType
    ! -- dummy
    class(GwtIsttype) :: this
    ! -- locals
    !character(len=LINELENGTH) :: errmsg
    type(GwtIstParamFoundType) :: found
    integer(I4B) :: asize
    integer(I4B), dimension(:), pointer, contiguous :: map
    !
    ! -- set map to convert user input data into reduced data
    map => null()
    if (this%dis%nodes < this%dis%nodesuser) map => this%dis%nodeuser
    !
    ! -- reallocate
    if (this%isrb == 0) then
      call get_isize('BULK_DENSITY', this%input_mempath, asize)
      if (asize >= 0) &
        call mem_reallocate(this%bulk_density, this%dis%nodes, &
                            'BULK_DENSITY', trim(this%memoryPath))
      call get_isize('DISTCOEF', this%input_mempath, asize)
      if (asize >= 0) &
        call mem_reallocate(this%distcoef, this%dis%nodes, 'DISTCOEF', &
                            trim(this%memoryPath))
    end if
    if (this%idcy == 0) then
      call get_isize('DECAY', this%input_mempath, asize)
      if (asize >= 0) &
        call mem_reallocate(this%decay, this%dis%nodes, 'DECAY', &
                            trim(this%memoryPath))
    end if
    call get_isize('DECAY_SORBED', this%input_mempath, asize)
    if (asize >= 0) then
      call mem_reallocate(this%decay_sorbed, this%dis%nodes, &
                          'DECAY_SORBED', trim(this%memoryPath))
    end if
    if (this%isrb < 2) then
      call get_isize('SP2', this%input_mempath, asize)
      if (asize >= 0) &
        call mem_reallocate(this%sp2, this%dis%nodes, 'SP2', &
                            trim(this%memoryPath))
    end if
    !
    ! -- update defaults with memory sourced values
    call mem_set_value(this%porosity, 'POROSITY', this%input_mempath, map, &
                       found%porosity)
    call mem_set_value(this%volfrac, 'VOLFRAC', this%input_mempath, map, &
                       found%volfrac)
    call mem_set_value(this%zetaim, 'ZETAIM', this%input_mempath, map, &
                       found%zetaim)
    call mem_set_value(this%cim, 'CIM', this%input_mempath, map, &
                       found%cim)
    call mem_set_value(this%decay, 'DECAY', this%input_mempath, map, &
                       found%decay)
    call mem_set_value(this%decay_sorbed, 'DECAY_SORBED', this%input_mempath, &
                       map, found%decay_sorbed)
    call mem_set_value(this%bulk_density, 'BULK_DENSITY', this%input_mempath, &
                       map, found%bulk_density)
    call mem_set_value(this%distcoef, 'DISTCOEF', this%input_mempath, map, &
                       found%distcoef)
    call mem_set_value(this%sp2, 'SP2', this%input_mempath, map, &
                       found%sp2)

    ! -- log options
    if (this%iout > 0) then
      call this%log_data(found)
    end if

    ! -- Check for required sorption variables
    if (this%isrb > 0) then
      if (.not. found%bulk_density) then
        write (errmsg, '(a)') 'Sorption is active but BULK_DENSITY &
          &not specified.  BULK_DENSITY must be specified in GRIDDATA block.'
        call store_error(errmsg)
      end if
      if (.not. found%distcoef) then
        write (errmsg, '(a)') 'Sorption is active but distribution &
          &coefficient not specified.  DISTCOEF must be specified in &
          &GRIDDATA block.'
        call store_error(errmsg)
      end if
      if (this%isrb > 1) then
        if (.not. found%sp2) then
          write (errmsg, '(a)') 'Freundlich or langmuir sorption is active &
            &but SP2 not specified.  SP2 must be specified in &
            &GRIDDATA block.'
          call store_error(errmsg)
        end if
      end if
    else
      if (found%bulk_density) then
        write (warnmsg, '(a)') 'Sorption is not active but &
          &BULK_DENSITY was specified.  BULK_DENSITY will have no affect on &
          &simulation results.'
        call store_warning(warnmsg)
        write (this%iout, '(1x,a)') 'WARNING.  '//warnmsg
      end if
      if (found%distcoef) then
        write (warnmsg, '(a)') 'Sorption is not active but &
          &distribution coefficient was specified.  DISTCOEF will have &
          &no affect on simulation results.'
        call store_warning(warnmsg)
        write (this%iout, '(1x,a)') 'WARNING.  '//warnmsg
      end if
      if (found%sp2) then
        write (warnmsg, '(a)') 'Sorption is not active but &
          &SP2 was specified.  SP2 will have &
          &no affect on simulation results.'
        call store_warning(warnmsg)
        write (this%iout, '(1x,a)') 'WARNING.  '//warnmsg
      end if
    end if

    ! -- Check for required decay/production rate coefficients
    if (this%idcy > 0) then
      if (.not. found%decay) then
        write (errmsg, '(a)') 'First or zero order decay is &
          &active but the first rate coefficient is not specified.  DECAY &
          &must be specified in GRIDDATA block.'
        call store_error(errmsg)
      end if
      if (.not. found%decay_sorbed) then
        !
        ! -- If DECAY_SORBED not specified and sorption is active, then
        !    terminate with an error
        write (errmsg, '(a)') 'DECAY_SORBED not provided in GRIDDATA &
          &block but decay and sorption are active.  Specify DECAY_SORBED &
          &in GRIDDATA block.'
        call store_error(errmsg)
      end if
    else
      if (found%decay) then
        write (warnmsg, '(a)') 'First- or zero-order decay &
          &is not active but decay was specified.  DECAY will &
          &have no affect on simulation results.'
        call store_warning(warnmsg)
        write (this%iout, '(1x,a)') 'WARNING.  '//warnmsg
      end if
      if (found%decay_sorbed) then
        write (warnmsg, '(a)') 'First- or zero-order decay &
          &is not active but DECAY_SORBED was specified.  &
          &DECAY_SORBED will have no affect on simulation results.'
        call store_warning(warnmsg)
        write (this%iout, '(1x,a)') 'WARNING.  '//warnmsg
      end if
    end if

    ! -- Check for required dual domain arrays or warn if they are specified
    !    but won't be used.
    if (.not. found%cim) then
      write (this%iout, '(1x,a)') 'Warning.  Dual domain is active but &
        &initial immobile domain concentration was not specified.  &
        &Setting CIM to zero.'
    end if
    if (.not. found%zetaim) then
      write (errmsg, '(a)') 'Dual domain is active but dual &
        &domain mass transfer rate (ZETAIM) was not specified.  ZETAIM &
        &must be specified in GRIDDATA block.'
      call store_error(errmsg)
    end if
    if (.not. found%porosity) then
      write (errmsg, '(a)') 'Dual domain is active but &
        &immobile domain POROSITY was not specified.  POROSITY &
        &must be specified in GRIDDATA block.'
      call store_error(errmsg)
    end if
    if (.not. found%volfrac) then
      write (errmsg, '(a)') 'Dual domain is active but &
        &immobile domain VOLFRAC was not specified.  VOLFRAC &
        &must be specified in GRIDDATA block. This is a new &
        &requirement for MODFLOW versions later than version &
        &6.4.1.'
      call store_error(errmsg)
    end if

    ! -- terminate if errors
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine source_data

  !> @brief Log user data to list file
  !<
  subroutine log_data(this, found)
    use GwtIstInputModule, only: GwtIstParamFoundType
    class(GwTIstType), intent(inout) :: this
    type(GwtIstParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'PROCESSING GRIDDATA'
    if (found%porosity) then
      write (this%iout, '(4x,a)') 'MOBILE DOMAIN POROSITY set from input file'
    end if
    if (found%volfrac) then
      write (this%iout, '(4x,a)') 'VOLFRAC set from input file'
    end if
    if (found%zetaim) then
      write (this%iout, '(4x,a)') 'ZETAIM set from input file'
    end if
    if (found%cim) then
      write (this%iout, '(4x,a)') 'CIM set from input file'
    end if
    if (found%decay) then
      write (this%iout, '(4x,a)') 'DECAY RATE set from input file'
    end if
    if (found%decay_sorbed) then
      write (this%iout, '(4x,a)') 'DECAY SORBED RATE set from input file'
    end if
    if (found%bulk_density) then
      write (this%iout, '(4x,a)') 'BULK DENSITY set from input file'
    end if
    if (found%distcoef) then
      write (this%iout, '(4x,a)') 'DISTRIBUTION COEFFICIENT set from input file'
    end if
    if (found%sp2) then
      write (this%iout, '(4x,a)') 'SECOND SORPTION PARAM set from input file'
    end if
    write (this%iout, '(1x,a)') 'END PROCESSING GRIDDATA'
  end subroutine log_data

  !> @ brief Read dimensions for package
  !!
  !!  Read dimensions for package.
  !!
  !<
  subroutine ist_read_dimensions(this)
    ! -- dummy
    class(GwtIstType), intent(inout) :: this !< GwtIstType object
    ! -- local
    ! -- format
  end subroutine ist_read_dimensions

  !> @ brief Read data for package
  !!
  !!  Read data for package.
  !!
  !<
  subroutine read_data(this)
    ! -- modules
    ! -- dummy
    class(GwtIstType) :: this !< GwtIstType object

    ! -- source data
    call this%source_data
  end subroutine read_data

  !> @ brief Return thetaim
  !!
  !!  Calculate and return thetaim, volume of immobile voids per aquifer volume
  !!
  !<
  function get_thetaim(this, node) result(thetaim)
    ! -- modules
    ! -- dummy
    class(GwtIstType) :: this !< GwtIstType object
    integer(I4B), intent(in) :: node !< node number
    ! -- return
    real(DP) :: thetaim
    !
    thetaim = this%volfrac(node) * this%porosity(node)
  end function get_thetaim

  !> @ brief Calculate immobile domain equation terms
  !!
  !!  This subroutine calculates the immobile domain (or dual domain) terms.
  !!  The resulting ddterm and f terms are described in the GWT model report.
  !!  The terms are concentration coefficients used in the balance equation
  !!  for the immobile domain.
  !!
  !<
  subroutine get_ddterm(thetaim, vcell, delt, swtpdt, &
                        volfracim, rhobim, kdnew, kdold, &
                        lambda1im, lambda2im, &
                        gamma1im, gamma2im, zetaim, ddterm, f)
    ! -- dummy
    real(DP), intent(in) :: thetaim !< immobile domain porosity
    real(DP), intent(in) :: vcell !< volume of cell
    real(DP), intent(in) :: delt !< length of time step
    real(DP), intent(in) :: swtpdt !< cell saturation at end of time step
    real(DP), intent(in) :: volfracim !< volume fraction of immobile domain
    real(DP), intent(in) :: rhobim !< bulk density for the immobile domain (fim * rhob)
    real(DP), intent(in) :: kdnew !< effective distribution coefficient for new time
    real(DP), intent(in) :: kdold !< effective distribution coefficient for old time
    real(DP), intent(in) :: lambda1im !< first-order decay rate in aqueous phase
    real(DP), intent(in) :: lambda2im !< first-order decay rate in sorbed phase
    real(DP), intent(in) :: gamma1im !< zero-order decay rate in aqueous phase
    real(DP), intent(in) :: gamma2im !< zero-order decay rate in sorbed phase
    real(DP), intent(in) :: zetaim !< transfer coefficient between mobile and immobile domains
    real(DP), dimension(:), intent(inout) :: ddterm !< nine terms comprising the balance equation of the immobile domain
    real(DP), intent(inout) :: f !< the f term used to calculate the immobile domain concentration
    ! -- local
    real(DP) :: tled
    !
    ! -- initialize
    tled = DONE / delt
    !
    ! -- Calculate terms.  These terms correspond to the concentration
    !    coefficients in equation 7-4 of the GWT model report.  However,
    !    an updated equation is presented as 9-9 in the supplemental technical
    !    information guide (mf6suptechinfo.pdf)
    ddterm(1) = thetaim * vcell * tled
    ddterm(2) = thetaim * vcell * tled
    ddterm(3) = volfracim * rhobim * vcell * kdnew * tled
    ddterm(4) = volfracim * rhobim * vcell * kdold * tled
    ddterm(5) = thetaim * lambda1im * vcell
    ddterm(6) = lambda2im * volfracim * rhobim * kdnew * vcell
    ddterm(7) = thetaim * gamma1im * vcell
    ddterm(8) = gamma2im * volfracim * rhobim * vcell
    ddterm(9) = vcell * swtpdt * zetaim
    !
    ! -- calculate denominator term, f
    f = ddterm(1) + ddterm(3) + ddterm(5) + ddterm(6) + ddterm(9)
  end subroutine get_ddterm

  !> @ brief Calculate the hcof and rhs terms for immobile domain
  !!
  !!  This subroutine calculates the hcof and rhs terms that must
  !!  be added to the solution system of equations
  !!
  !<
  subroutine get_hcofrhs(ddterm, f, cimold, hcof, rhs)
    ! -- dummy
    real(DP), dimension(:), intent(in) :: ddterm !< terms comprising the balance equation of the immobile domain
    real(DP), intent(in) :: f !< the f term used to calculate the immobile domain concentration
    real(DP), intent(in) :: cimold !< immobile domain concentration at end of last time step
    real(DP), intent(inout) :: hcof !< calculated contribution for the a-matrix diagonal position
    real(DP), intent(inout) :: rhs !< calculated contribution for the solution right-hand side
    !
    ! -- calculate hcof
    hcof = ddterm(9)**2 / f - ddterm(9)
    !
    ! -- calculate rhs, and switch the sign because this term needs to
    !    be moved to the left hand side
    rhs = (ddterm(2) + ddterm(4)) * cimold - ddterm(7) - ddterm(8)
    rhs = rhs * ddterm(9) / f
    rhs = -rhs
  end subroutine get_hcofrhs

  !> @ brief Calculate the immobile domain concentration
  !!
  !!  This function calculates the concentration of the immobile domain.
  !!
  !<
  function get_ddconc(ddterm, f, cimold, cnew) result(cimnew)
    ! -- dummy
    real(DP), dimension(:), intent(in) :: ddterm !< terms comprising the balance equation of the immobile domain
    real(DP), intent(in) :: f !< the f term used to calculate the immobile domain concentration
    real(DP), intent(in) :: cimold !< immobile domain concentration at end of last time step
    real(DP), intent(in) :: cnew !< concentration of the mobile domain at the end of the time step
    ! -- result
    real(DP) :: cimnew !< calculated concentration of the immobile domain
    ! -- local
    !
    ! -- calculate ddconc
    cimnew = (ddterm(2) + ddterm(4)) * cimold + ddterm(9) * cnew - ddterm(7) &
             - ddterm(8)
    cimnew = cimnew / f
  end function get_ddconc

  !> @ brief Calculate the immobile domain budget terms
  !!
  !!  This subroutine calculates and accumulates the immobile domain
  !!  budget terms into the budterm accumulator
  !!
  !<
  subroutine accumulate_budterm(budterm, ddterm, cimnew, cimold, cnew, idcy)
    ! -- modules
    ! -- dummy
    real(DP), dimension(:, :), intent(inout) :: budterm !<
    real(DP), dimension(:), intent(in) :: ddterm !< terms comprising the balance equation of the immobile domain
    real(DP), intent(in) :: cimnew !< immobile domain concenration at the end of this time step
    real(DP), intent(in) :: cimold !< immobile domain concentration at end of last time step
    real(DP), intent(in) :: cnew !< mobile domain concentration at the end of this time step
    integer(I4B), intent(in) :: idcy !< order of decay rate (0:none, 1:first, 2:zero)
    ! -- local
    real(DP) :: rate
    integer(I4B) :: i
    !
    ! -- calculate STORAGE-AQUEOUS
    i = 1
    rate = -ddterm(1) * cimnew + ddterm(2) * cimold
    if (rate > DZERO) then
      budterm(1, i) = budterm(1, i) + rate
    else
      budterm(2, i) = budterm(2, i) - rate
    end if
    !
    ! -- calculate STORAGE-SORBED
    i = 2
    rate = -ddterm(3) * cimnew + ddterm(4) * cimold
    if (rate > DZERO) then
      budterm(1, i) = budterm(1, i) + rate
    else
      budterm(2, i) = budterm(2, i) - rate
    end if
    !
    ! -- calculate DECAY-AQUEOUS
    i = 3
    rate = DZERO
    if (idcy == 1) then
      rate = -ddterm(5) * cimnew
    else if (idcy == 2) then
      rate = -ddterm(7)
    else
      rate = DZERO
    end if
    if (rate > DZERO) then
      budterm(1, i) = budterm(1, i) + rate
    else
      budterm(2, i) = budterm(2, i) - rate
    end if
    !
    ! -- calculate DECAY-SORBED
    i = 4
    if (idcy == 1) then
      rate = -ddterm(6) * cimnew
    else if (idcy == 2) then
      rate = -ddterm(8)
    else
      rate = DZERO
    end if
    if (rate > DZERO) then
      budterm(1, i) = budterm(1, i) + rate
    else
      budterm(2, i) = budterm(2, i) - rate
    end if
    !
    ! -- calculate MOBILE-DOMAIN
    i = 5
    rate = ddterm(9) * cnew - ddterm(9) * cimnew
    if (rate > DZERO) then
      budterm(1, i) = budterm(1, i) + rate
    else
      budterm(2, i) = budterm(2, i) - rate
    end if
    !
  end subroutine accumulate_budterm

end module GwtIstModule
