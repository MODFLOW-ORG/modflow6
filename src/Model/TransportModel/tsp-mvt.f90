! -- Groundwater Transport Mover Module
! -- This module is responsible for sending mass from providers into
! -- receiver qmfrommvr arrays and writing a mover transport budget

module TspMvtModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH, MAXCHARLEN, DZERO, LENPAKLOC, &
                             DNODATA, LENPACKAGENAME, TABCENTER, &
                             LENMODELNAME, LENVARNAME

  use SimModule, only: store_error
  use BaseDisModule, only: DisBaseType
  use NumericalPackageModule, only: NumericalPackageType
  use TspFmiModule, only: TspFmiType
  use BudgetModule, only: BudgetType, budget_cr
  use BudgetObjectModule, only: BudgetObjectType, budgetobject_cr
  use TableModule, only: TableType, table_cr

  implicit none

  private
  public :: TspMvtType
  public :: mvt_cr

  type, extends(NumericalPackageType) :: TspMvtType

    character(len=LENMODELNAME) :: gwfmodelname1 = '' !< name of model 1
    character(len=LENMODELNAME) :: gwfmodelname2 = '' !< name of model 2 (set to modelname 1 for single model MVT)
    integer(I4B), pointer :: maxpackages !< max number of packages
    integer(I4B), pointer :: ibudgetout => null() !< unit number for budget output file
    integer(I4B), pointer :: ibudcsv => null() !< unit number for csv budget output file
    real(DP), pointer :: eqnsclfac => null() !< governing equation scale factor; =1. for solute; =rhow*cpw for energy
    type(TspFmiType), pointer :: fmi1 => null() !< pointer to fmi object for model 1
    type(TspFmiType), pointer :: fmi2 => null() !< pointer to fmi object for model 2 (set to fmi1 for single model)
    type(BudgetType), pointer :: budget => null() !< mover transport budget object (used to write balance table)
    type(BudgetObjectType), pointer :: budobj => null() !< budget container (used to write binary file)
    type(BudgetObjectType), pointer :: mvrbudobj => null() !< pointer to the water mover budget object
    character(len=LENPACKAGENAME), &
      dimension(:), pointer, contiguous :: paknames => null() !< array of package names
    character(len=LENVARNAME) :: depvartype = ''
    !
    ! -- table objects
    type(TableType), pointer :: outputtab => null()

  contains

    procedure :: mvt_df
    procedure :: mvt_ar
    procedure :: mvt_rp
    procedure :: mvt_fc
    procedure :: mvt_cc
    procedure :: mvt_bd
    procedure :: mvt_ot_saveflow
    procedure :: mvt_ot_printflow
    procedure :: mvt_ot_bdsummary
    procedure :: mvt_da
    procedure :: allocate_scalars
    procedure :: read_options
    procedure :: mvt_setup_budobj
    procedure :: mvt_fill_budobj
    procedure :: mvt_scan_mvrbudobj
    procedure :: set_pointer_mvrbudobj
    procedure :: set_fmi_pr_rc
    procedure, private :: mvt_setup_outputtab
    procedure, private :: mvt_print_outputtab
  end type TspMvtType

contains

  !> @brief Create a new mover transport object
  !<
  subroutine mvt_cr(mvt, name_model, inunit, iout, fmi1, eqnsclfac, &
                    depvartype, gwfmodelname1, gwfmodelname2, fmi2)
    ! -- dummy
    type(TspMvtType), pointer :: mvt
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    type(TspFmiType), intent(in), target :: fmi1
    real(DP), intent(in), pointer :: eqnsclfac !< governing equation scale factor
    character(len=LENVARNAME), intent(in) :: depvartype !< dependent variable type ('concentration' or 'temperature')
    character(len=*), intent(in), optional :: gwfmodelname1
    character(len=*), intent(in), optional :: gwfmodelname2
    type(TspFmiType), intent(in), target, optional :: fmi2
    !
    ! -- Create the object
    allocate (mvt)
    !
    ! -- Create name and memory path
    call mvt%set_names(1, name_model, 'MVT', 'MVT')
    !
    ! -- Allocate scalars
    call mvt%allocate_scalars()
    !
    mvt%inunit = inunit
    mvt%iout = iout
    !
    ! -- Assume that this MVT is owned by a GWT Model
    mvt%fmi1 => fmi1
    mvt%fmi2 => fmi1
    !
    ! -- Set pointers
    if (present(fmi2)) then
      mvt%fmi2 => fmi2
    end if
    !
    ! -- Set model names
    if (present(gwfmodelname1)) then
      mvt%gwfmodelname1 = gwfmodelname1
    end if
    if (present(gwfmodelname2)) then
      mvt%gwfmodelname2 = gwfmodelname2
    end if
    !
    ! -- Create the budget object
    call budgetobject_cr(mvt%budobj, 'TRANSPORT MOVER')
    !
    ! -- Store pointer to governing equation scale factor
    mvt%eqnsclfac => eqnsclfac
    !
    ! -- Store pointer to labels associated with the current model so that the
    !    package has access to the corresponding dependent variable type
    mvt%depvartype = depvartype
  end subroutine mvt_cr

  !> @brief Define mover transport object
  !<
  subroutine mvt_df(this, dis)
    ! -- dummy
    class(TspMvtType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    ! -- formats
    character(len=*), parameter :: fmtmvt = &
      "(1x,/1x,'MVT -- MOVER TRANSPORT PACKAGE, VERSION 1, 4/15/2020', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! -- Set pointer to dis
    this%dis => dis
    !
    ! -- Print a message identifying the MVT package.
    write (this%iout, fmtmvt) this%inunit
    !
    ! -- Initialize block parser
    call this%parser%Initialize(this%inunit, this%iout)
    !
    ! -- Initialize the budget table writer
    call budget_cr(this%budget, this%memoryPath)
    !
    ! -- Read mvt options
    call this%read_options()
  end subroutine mvt_df

  !> @ brief Set pointer to mvrbudobj
  !!
  !! Store a pointer to mvrbudobj, which contains the simulated water
  !! mover flows from either a gwf model MVR package or from a gwf-gwf
  !! exchange MVR package.
  !<
  subroutine set_pointer_mvrbudobj(this, mvrbudobj)
    class(TspMvtType) :: this
    type(BudgetObjectType), intent(in), target :: mvrbudobj
    this%mvrbudobj => mvrbudobj
  end subroutine set_pointer_mvrbudobj

  !> @brief Allocate and read mover-for-transport information
  !<
  subroutine mvt_ar(this)
    ! -- dummy
    class(TspMvtType) :: this
    !
    ! -- Setup the output table
    call this%mvt_setup_outputtab()
  end subroutine mvt_ar

  !> @brief Read and prepare mover transport object
  !<
  subroutine mvt_rp(this)
    ! -- modules
    use TdisModule, only: kper, kstp
    ! -- dummy
    class(TspMvtType) :: this
    !
    ! -- At this point, the mvrbudobj is available to set up the mvt budobj
    if (kper * kstp == 1) then
      !
      ! -- If mvt is for a single model then point to fmi1
      !cdl todo: this needs to be called from GwtGwtExg somehow for the 2 model case
      if (associated(this%fmi1, this%fmi2)) then
        call this%set_pointer_mvrbudobj(this%fmi1%mvrbudobj)
      end if
      !
      ! -- Set up the mvt budobject
      call this%mvt_scan_mvrbudobj()
      call this%mvt_setup_budobj()
      !
      ! -- Define the budget object to be the size of maxpackages
      call this%budget%budget_df(this%maxpackages, 'TRANSPORT MOVER', bddim='M')
      call this%budget%set_ibudcsv(this%ibudcsv)
    end if
  end subroutine mvt_rp

  !> @brief Calculate coefficients and fill amat and rhs
  !!
  !! The mvt package adds the mass flow rate to the provider qmfrommvr array.
  !! The advanced packages know enough to subtract any mass that is leaving, so
  !! the mvt just adds mass coming in from elsewhere.  Because the movers
  !! change by stress period, their solute effects must be added to the right-
  !! hand side of the transport matrix equations.
  !<
  subroutine mvt_fc(this, cnew1, cnew2)
    ! -- dummy
    class(TspMvtType) :: this
    real(DP), intent(in), dimension(:), contiguous, target :: cnew1
    real(DP), intent(in), dimension(:), contiguous, target :: cnew2
    ! -- local
    integer(I4B) :: i, n
    integer(I4B) :: id1, id2, nlist
    integer(I4B) :: ipr, irc
    integer(I4B) :: igwtnode
    integer(I4B) :: nbudterm
    real(DP) :: q, cp
    real(DP), dimension(:), pointer :: concpak
    real(DP), dimension(:), contiguous, pointer :: cnew
    type(TspFmiType), pointer :: fmi_pr !< pointer to provider model fmi package
    type(TspFmiType), pointer :: fmi_rc !< pointer to receiver model fmi package
    !
    ! -- Add mover QC terms to the receiver packages
    nbudterm = this%mvrbudobj%nbudterm
    do i = 1, nbudterm
      nlist = this%mvrbudobj%budterm(i)%nlist
      if (nlist > 0) then
        !
        ! -- Set pointers to the fmi packages for the provider and the receiver
        call this%set_fmi_pr_rc(i, fmi_pr, fmi_rc)
        !
        ! -- Set a pointer to the GWT model concentration (or temperature)
        !    associated with the provider
        cnew => cnew1
        if (associated(fmi_pr, this%fmi2)) then
          cnew => cnew2
        end if
        !
        !-- Get the package index for the provider
        call fmi_pr%get_package_index(this%mvrbudobj%budterm(i)%text2id1, ipr)
        !
        ! -- Get the package index for the receiver
        call fmi_rc%get_package_index(this%mvrbudobj%budterm(i)%text2id2, irc)
        !
        ! -- If provider is an advanced package, then set a pointer to its simulated concentration
        if (fmi_pr%iatp(ipr) /= 0) then
          concpak => fmi_pr%datp(ipr)%concpack
        end if
        !
        ! -- Process flows for each entry in the list and add mass to receivers
        do n = 1, nlist
          !
          ! -- lak/sfr/maw/uzf id1 (provider) and id2 (receiver)
          id1 = this%mvrbudobj%budterm(i)%id1(n)
          id2 = this%mvrbudobj%budterm(i)%id2(n)
          !
          ! -- Obtain mover flow rate from the mover flow budget object
          q = this%mvrbudobj%budterm(i)%flow(n)
          !
          ! -- Assign concentration of the provider
          cp = DZERO
          if (fmi_pr%iatp(ipr) /= 0) then
            !
            ! -- Provider package is being represented by an APT (SFT, LKT, MWT, UZT,
            !    SFE, LKE, MWE, UZE) so set the dependent variable (concentration or
            !    temperature) to the simulated dependent variable of the APT.
            cp = concpak(id1)
          else
            !
            ! -- Provider is a regular stress package (WEL, DRN, RIV, etc.) or the
            !    provider is an advanced stress package but is not represented with
            !    SFT, LKT, MWT, or UZT, so use the GWT cell concentration
            igwtnode = fmi_pr%gwfpackages(ipr)%nodelist(id1)
            cp = cnew(igwtnode)
            !
          end if
          !
          ! -- Add the mover rate times the provider concentration into the receiver
          !    make sure these are accumulated since multiple providers can move
          !    water into the same receiver
          if (fmi_rc%iatp(irc) /= 0) then
            fmi_rc%datp(irc)%qmfrommvr(id2) = fmi_rc%datp(irc)%qmfrommvr(id2) - &
                                              q * cp * this%eqnsclfac
          end if
        end do
      end if
    end do
  end subroutine mvt_fc

  !> @ brief Set the fmi_pr and fmi_rc pointers
  !!
  !! The fmi_pr and fmi_rc arguments are pointers to the provider and receiver
  !! FMI Packages.  If this MVT Package is owned by a single GWT model, then
  !! these pointers are both set to the FMI Package of this GWT model's FMI
  !! package.  If this MVT package is owned by a GWTGWT exchange, then the
  !! fmi_pr and fmi_rc pointers may be assigned to FMI Packages in different
  !! models.
  !<
  subroutine set_fmi_pr_rc(this, ibudterm, fmi_pr, fmi_rc)
    ! -- dummy
    class(TspMvtType) :: this
    integer(I4B), intent(in) :: ibudterm
    type(TspFmiType), pointer :: fmi_pr
    type(TspFmiType), pointer :: fmi_rc
    !
    fmi_pr => null()
    fmi_rc => null()
    if (this%gwfmodelname1 == '' .and. this%gwfmodelname2 == '') then
      fmi_pr => this%fmi1
      fmi_rc => this%fmi1
    else
      ! -- Modelname for provider is this%mvrbudobj%budterm(i)%text1id1
      if (this%mvrbudobj%budterm(ibudterm)%text1id1 == this%gwfmodelname1) then
        ! -- Model 1 is the provider
        fmi_pr => this%fmi1
      else if (this%mvrbudobj%budterm(ibudterm)%text1id1 == &
               this%gwfmodelname2) then
        ! -- Model 2 is the provider
        fmi_pr => this%fmi2
      else
        ! -- Must be an error
        !cdl todo: programming error
        print *, this%mvrbudobj%budterm(ibudterm)%text1id1
        print *, this%gwfmodelname1
        print *, this%gwfmodelname2
        stop "error in set_fmi_pr_rc"
      end if
      !
      ! -- Modelname for receiver is this%mvrbudobj%budterm(i)%text1id2
      if (this%mvrbudobj%budterm(ibudterm)%text1id2 == this%gwfmodelname1) then
        ! -- Model 1 is the receiver
        fmi_rc => this%fmi1
      else if (this%mvrbudobj%budterm(ibudterm)%text1id2 == &
               this%gwfmodelname2) then
        ! -- Model 2 is the receiver
        fmi_rc => this%fmi2
      else
        ! -- Must be an error
        !cdl todo: programming error
        print *, this%mvrbudobj%budterm(ibudterm)%text1id2
        print *, this%gwfmodelname1
        print *, this%gwfmodelname2
        stop "error in set_fmi_pr_rc"
      end if
    end if
    !
    if (.not. associated(fmi_pr)) then
      print *, 'Could not find FMI Package...'
      stop "error in set_fmi_pr_rc"
    end if
    if (.not. associated(fmi_rc)) then
      print *, 'Could not find FMI Package...'
      stop "error in set_fmi_pr_rc"
    end if
  end subroutine set_fmi_pr_rc

  !> @brief Extra convergence check for mover
  !<
  subroutine mvt_cc(this, kiter, iend, icnvgmod, cpak, dpak)
    ! -- dummy
    class(TspMvtType) :: this
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: iend
    integer(I4B), intent(in) :: icnvgmod
    character(len=LENPAKLOC), intent(inout) :: cpak
    real(DP), intent(inout) :: dpak
    ! -- formats
    character(len=*), parameter :: fmtmvrcnvg = &
      "(/,1x,'MOVER PACKAGE REQUIRES AT LEAST TWO OUTER ITERATIONS. CONVERGE &
      &FLAG HAS BEEN RESET TO FALSE.')"
    !
    ! -- If there are active movers, then at least 2 outers required
    if (associated(this%mvrbudobj)) then
      if (icnvgmod == 1 .and. kiter == 1) then
        dpak = DNODATA
        cpak = trim(this%packName)
        write (this%iout, fmtmvrcnvg)
      end if
    end if
  end subroutine mvt_cc

  !> @brief Write mover terms to listing file
  !<
  subroutine mvt_bd(this, cnew1, cnew2)
    ! -- dummy
    class(TspMvtType) :: this
    real(DP), dimension(:), contiguous, intent(in) :: cnew1
    real(DP), dimension(:), contiguous, intent(in) :: cnew2
    !
    ! -- Fill the budget object
    call this%mvt_fill_budobj(cnew1, cnew2)
  end subroutine mvt_bd

  !> @brief Write mover budget terms
  !<
  subroutine mvt_ot_saveflow(this, icbcfl, ibudfl)
    ! -- modules
    use TdisModule, only: kstp, kper, delt, pertim, totim
    ! -- dummy
    class(TspMvttype) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    ! -- locals
    integer(I4B) :: ibinun
    !
    ! -- Save the mover flows from the budobj to a mover binary file
    ibinun = 0
    if (this%ibudgetout /= 0) then
      ibinun = this%ibudgetout
    end if
    if (icbcfl == 0) ibinun = 0
    if (ibinun > 0) then
      call this%budobj%save_flows(this%dis, ibinun, kstp, kper, delt, &
                                  pertim, totim, this%iout)
    end if
  end subroutine mvt_ot_saveflow

  !> @brief Print mover flow table
  !<
  subroutine mvt_ot_printflow(this, icbcfl, ibudfl)
    ! -- dummy
    class(TspMvtType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    !
    ! -- Print the mover flow table
    if (ibudfl /= 0 .and. this%iprflow /= 0) then
      call this%mvt_print_outputtab()
    end if
  end subroutine mvt_ot_printflow

  !> @brief Write mover budget to listing file
  !<
  subroutine mvt_ot_bdsummary(this, ibudfl)
    ! -- modules
    use TdisModule, only: kstp, kper, delt, totim
    ! -- dummy
    class(TspMvtType) :: this
    integer(I4B), intent(in) :: ibudfl
    ! -- locals
    integer(I4B) :: i, j, n
    real(DP), allocatable, dimension(:) :: ratin, ratout
    !
    ! -- Allocate and initialize ratin/ratout
    allocate (ratin(this%maxpackages), ratout(this%maxpackages))
    do j = 1, this%maxpackages
      ratin(j) = DZERO
      ratout(j) = DZERO
    end do
    !
    ! -- Accumulate the rates
    do i = 1, this%maxpackages
      do j = 1, this%budobj%nbudterm
        do n = 1, this%budobj%budterm(j)%nlist
          !
          ! -- Provider is inflow to mover
          if (this%paknames(i) == this%budobj%budterm(j)%text2id1) then
            ratin(i) = ratin(i) + this%budobj%budterm(j)%flow(n)
          end if
          !
          ! -- Receiver is outflow from mover
          if (this%paknames(i) == this%budobj%budterm(j)%text2id2) then
            ratout(i) = ratout(i) + this%budobj%budterm(j)%flow(n)
          end if
        end do
      end do
    end do
    !
    ! -- Send rates to budget object
    call this%budget%reset()
    do j = 1, this%maxpackages
      call this%budget%addentry(ratin(j), ratout(j), delt, this%paknames(j))
    end do
    !
    ! -- Write the budget
    call this%budget%finalize_step(delt)
    if (ibudfl /= 0) then
      call this%budget%budget_ot(kstp, kper, this%iout)
    end if
    !
    ! -- Write budget csv
    call this%budget%writecsv(totim)
    !
    ! -- Deallocate
    deallocate (ratin, ratout)
    !
    ! -- Output mvr budget
    !    Not using budobj write_table here because it would result
    !    in a table that has one entry.  A custom table looks
    !    better here with a row for each package.
    !call this%budobj%write_budtable(kstp, kper, this%iout)
  end subroutine mvt_ot_bdsummary

  !> @ brief Deallocate memory
  !!
  !!  Method to deallocate memory for the package.
  !<
  subroutine mvt_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(TspMvtType) :: this
    !
    ! -- Deallocate arrays if package was active
    if (this%inunit > 0) then
      !
      ! -- Character array
      deallocate (this%paknames)
      !
      ! -- Budget object
      call this%budget%budget_da()
      deallocate (this%budget)
      !
      ! -- Budobj
      call this%budobj%budgetobject_da()
      deallocate (this%budobj)
      nullify (this%budobj)
      !
      ! -- Output table object
      if (associated(this%outputtab)) then
        call this%outputtab%table_da()
        deallocate (this%outputtab)
        nullify (this%outputtab)
      end if
    end if
    !
    ! -- Scalars
    this%fmi1 => null()
    this%fmi1 => null()
    this%mvrbudobj => null()
    call mem_deallocate(this%maxpackages)
    call mem_deallocate(this%ibudgetout)
    call mem_deallocate(this%ibudcsv)
    !
    ! -- Deallocate scalars in NumericalPackageType
    call this%NumericalPackageType%da()
  end subroutine mvt_da

  !> @ brief Allocate scalar variables for package
  !!
  !!  Method to allocate scalar variables for the MVT package.
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy
    class(TspMvtType) :: this
    !
    ! -- Allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%maxpackages, 'MAXPACKAGES', this%memoryPath)
    call mem_allocate(this%ibudgetout, 'IBUDGETOUT', this%memoryPath)
    call mem_allocate(this%ibudcsv, 'IBUDCSV', this%memoryPath)
    !
    ! -- Initialize
    this%maxpackages = 0
    this%ibudgetout = 0
    this%ibudcsv = 0
  end subroutine allocate_scalars

  !> @brief Read mover-for-transport options block
  !<
  subroutine read_options(this)
    ! -- modules
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: assign_iounit, openfile
    ! -- dummy
    class(TspMvtType) :: this
    ! -- local
    character(len=LINELENGTH) :: errmsg, keyword
    character(len=MAXCHARLEN) :: fname
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    ! -- formats
    character(len=*), parameter :: fmtflow = &
      "(4x, a, 1x, a, 1x, ' WILL BE SAVED TO FILE: ', a, &
      &/4x, 'OPENED ON UNIT: ', I0)"
    character(len=*), parameter :: fmtflow2 = &
      &"(4x, 'FLOWS WILL BE SAVED TO BUDGET FILE')"
    !
    ! -- Get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, blockRequired=.false., &
                              supportOpenClose=.true.)
    !
    ! -- Parse options block if detected
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING MVT OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('SAVE_FLOWS')
          this%ipakcb = -1
          write (this%iout, fmtflow2)
        case ('PRINT_INPUT')
          this%iprpak = 1
          write (this%iout, '(4x,a)') 'MVT INPUT WILL BE PRINTED.'
        case ('PRINT_FLOWS')
          this%iprflow = 1
          write (this%iout, '(4x,a)') &
            'MVT FLOWS WILL BE PRINTED TO LISTING FILE.'
        case ('BUDGET')
          call this%parser%GetStringCaps(keyword)
          if (keyword == 'FILEOUT') then
            call this%parser%GetString(fname)
            call assign_iounit(this%ibudgetout, this%inunit, "BUDGET fileout")
            call openfile(this%ibudgetout, this%iout, fname, 'DATA(BINARY)', &
                          form, access, 'REPLACE')
            write (this%iout, fmtflow) 'MVT', 'BUDGET', trim(adjustl(fname)), &
              this%ibudgetout
          else
            call store_error('Optional BUDGET keyword must &
                             &be followed by FILEOUT')
          end if
        case ('BUDGETCSV')
          call this%parser%GetStringCaps(keyword)
          if (keyword == 'FILEOUT') then
            call this%parser%GetString(fname)
            call assign_iounit(this%ibudcsv, this%inunit, "BUDGETCSV fileout")
            call openfile(this%ibudcsv, this%iout, fname, 'CSV', &
                          filstat_opt='REPLACE')
            write (this%iout, fmtflow) 'MVT', 'BUDGET CSV', &
              trim(adjustl(fname)), this%ibudcsv
          else
            call store_error('Optional BUDGETCSV keyword must be followed by &
              &FILEOUT')
          end if
        case default
          write (errmsg, '(a,a)') 'Unknown MVT option: ', &
            trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (this%iout, '(1x,a)') 'END OF MVT OPTIONS'
    end if
  end subroutine read_options

  !> @brief Set up the budget object that stores all the mvr flows
  !<
  subroutine mvt_setup_budobj(this)
    ! -- modules
    use ConstantsModule, only: LENBUDTXT
    ! -- dummy
    class(TspMvtType) :: this
    ! -- local
    integer(I4B) :: nbudterm
    integer(I4B) :: ncv
    integer(I4B) :: maxlist
    integer(I4B) :: i
    integer(I4B) :: naux
    character(len=LENMODELNAME) :: modelname1, modelname2
    character(len=LENPACKAGENAME) :: packagename1, packagename2
    character(len=LENBUDTXT) :: text
    !
    ! -- Assign terms to set up the mover budget object
    nbudterm = this%mvrbudobj%nbudterm
    ncv = 0
    naux = 0
    if (this%depvartype == 'CONCENTRATION') then
      text = '        MVT-FLOW'
    else
      text = '        MVE-FLOW'
    end if
    !
    ! -- Set up budobj
    call this%budobj%budgetobject_df(ncv, nbudterm, 0, 0, bddim_opt='M')
    !
    ! -- Go through the water mover budget terms and set up the transport
    !    mover budget terms
    do i = 1, nbudterm
      modelname1 = this%mvrbudobj%budterm(i)%text1id1
      packagename1 = this%mvrbudobj%budterm(i)%text2id1
      modelname2 = this%mvrbudobj%budterm(i)%text1id2
      packagename2 = this%mvrbudobj%budterm(i)%text2id2
      maxlist = this%mvrbudobj%budterm(i)%maxlist
      call this%budobj%budterm(i)%initialize(text, &
                                             modelname1, &
                                             packagename1, &
                                             modelname2, &
                                             packagename2, &
                                             maxlist, .false., .false., &
                                             naux)
    end do
  end subroutine mvt_setup_budobj

  !> @brief Copy mover-for-transport flow terms into this%budobj
  !<
  subroutine mvt_fill_budobj(this, cnew1, cnew2)
    ! -- dummy
    class(TspMvtType) :: this
    real(DP), intent(in), dimension(:), contiguous, target :: cnew1
    real(DP), intent(in), dimension(:), contiguous, target :: cnew2
    ! -- local
    type(TspFmiType), pointer :: fmi_pr
    type(TspFmiType), pointer :: fmi_rc
    real(DP), dimension(:), contiguous, pointer :: cnew
    integer(I4B) :: nbudterm
    integer(I4B) :: nlist
    integer(I4B) :: ipr
    integer(I4B) :: irc
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: n1, n2
    integer(I4B) :: igwtnode
    real(DP) :: cp
    real(DP) :: q
    real(DP) :: rate
    !
    ! -- Go through the water mover budget terms and set up the transport
    !    mover budget terms
    nbudterm = this%mvrbudobj%nbudterm
    do i = 1, nbudterm
      nlist = this%mvrbudobj%budterm(i)%nlist
      call this%set_fmi_pr_rc(i, fmi_pr, fmi_rc)
      cnew => cnew1
      if (associated(fmi_pr, this%fmi2)) then
        cnew => cnew2
      end if
      call fmi_pr%get_package_index(this%mvrbudobj%budterm(i)%text2id1, ipr)
      call fmi_rc%get_package_index(this%mvrbudobj%budterm(i)%text2id2, irc)
      call this%budobj%budterm(i)%reset(nlist)
      do j = 1, nlist
        n1 = this%mvrbudobj%budterm(i)%id1(j)
        n2 = this%mvrbudobj%budterm(i)%id2(j)
        q = this%mvrbudobj%budterm(i)%flow(j)
        cp = DZERO
        if (fmi_pr%iatp(ipr) /= 0) then
          cp = fmi_pr%datp(ipr)%concpack(n1)
        else
          ! -- Must be a regular stress package
          igwtnode = fmi_pr%gwfpackages(ipr)%nodelist(n1)
          !cdl todo: need to set cnew to model 1; right now it is coming in as argument
          cp = cnew(igwtnode)
        end if
        !
        ! -- Calculate solute mover rate
        rate = DZERO
        if (fmi_rc%iatp(irc) /= 0) then
          rate = -q * cp * this%eqnsclfac
        end if
        !
        ! -- Add the rate to the budterm
        call this%budobj%budterm(i)%update_term(n1, n2, rate)
      end do
    end do
    !
    ! --Terms are filled, now accumulate them for this time step
    call this%budobj%accumulate_terms()
  end subroutine mvt_fill_budobj

  !> @brief Determine max number of packages in use
  !!
  !! Scan through the gwf water mover budget object and determine the maximum
  !! number of packages and unique package names
  !<
  subroutine mvt_scan_mvrbudobj(this)
    class(TspMvtType) :: this
    integer(I4B) :: nbudterm
    integer(I4B) :: maxpackages
    integer(I4B) :: i, j
    integer(I4B) :: ipos
    logical :: found
    !
    ! -- Calculate maxpackages, which is the the square of nbudterm
    nbudterm = this%mvrbudobj%nbudterm
    do i = 1, nbudterm
      if (i * i == nbudterm) then
        maxpackages = i
        exit
      end if
    end do
    this%maxpackages = maxpackages
    !
    ! -- Allocate paknames
    allocate (this%paknames(this%maxpackages))
    do i = 1, this%maxpackages
      this%paknames(i) = ''
    end do
    !
    ! -- Scan through mvrbudobj and create unique paknames
    ipos = 1
    do i = 1, nbudterm
      found = .false.
      do j = 1, ipos
        if (this%mvrbudobj%budterm(i)%text2id1 == this%paknames(j)) then
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        this%paknames(ipos) = this%mvrbudobj%budterm(i)%text2id1
        ipos = ipos + 1
      end if
    end do
  end subroutine mvt_scan_mvrbudobj

  !> @brief Set up the mover-for-transport output table
  !<
  subroutine mvt_setup_outputtab(this)
    ! -- dummy
    class(TspMvtType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    integer(I4B) :: ntabcol
    integer(I4B) :: maxrow
    integer(I4B) :: ilen
    !
    ! -- Allocate and initialize the output table
    if (this%iprflow /= 0) then
      !
      ! -- Dimension table
      ntabcol = 7
      maxrow = 0
      !
      ! -- Initialize the output table object
      title = 'TRANSPORT MOVER PACKAGE ('//trim(this%packName)// &
              ') FLOW RATES'
      call table_cr(this%outputtab, this%packName, title)
      call this%outputtab%table_df(maxrow, ntabcol, this%iout, &
                                   transient=.TRUE.)
      text = 'NUMBER'
      call this%outputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'PROVIDER LOCATION'
      ilen = LENMODELNAME + LENPACKAGENAME + 1
      call this%outputtab%initialize_column(text, ilen)
      text = 'PROVIDER ID'
      call this%outputtab%initialize_column(text, 10)
      text = 'PROVIDER FLOW RATE'
      call this%outputtab%initialize_column(text, 10)
      text = 'PROVIDER TRANSPORT RATE'
      call this%outputtab%initialize_column(text, 10)
      text = 'RECEIVER LOCATION'
      ilen = LENMODELNAME + LENPACKAGENAME + 1
      call this%outputtab%initialize_column(text, ilen)
      text = 'RECEIVER ID'
      call this%outputtab%initialize_column(text, 10)
      !
    end if
  end subroutine mvt_setup_outputtab

  !> @brief Set up mover-for-transport output table
  !<
  subroutine mvt_print_outputtab(this)
    ! -- module
    use TdisModule, only: kstp, kper
    ! -- dummy
    class(TspMvttype), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LENMODELNAME + LENPACKAGENAME + 1) :: cloc1, cloc2
    integer(I4B) :: i
    integer(I4B) :: n
    integer(I4B) :: inum
    integer(I4B) :: ntabrows
    integer(I4B) :: nlist
    !
    ! -- Determine number of table rows
    ntabrows = 0
    do i = 1, this%budobj%nbudterm
      nlist = this%budobj%budterm(i)%nlist
      ntabrows = ntabrows + nlist
    end do
    !
    ! -- Set table kstp and kper
    call this%outputtab%set_kstpkper(kstp, kper)
    !
    ! -- Add terms and print the table
    title = 'TRANSPORT MOVER PACKAGE ('//trim(this%packName)// &
            ') FLOW RATES'
    call this%outputtab%set_title(title)
    call this%outputtab%set_maxbound(ntabrows)
    !
    ! -- Process each table row
    inum = 1
    do i = 1, this%budobj%nbudterm
      nlist = this%budobj%budterm(i)%nlist
      do n = 1, nlist
        cloc1 = trim(adjustl(this%budobj%budterm(i)%text1id1))//' '// &
                trim(adjustl(this%budobj%budterm(i)%text2id1))
        cloc2 = trim(adjustl(this%budobj%budterm(i)%text1id2))//' '// &
                trim(adjustl(this%budobj%budterm(i)%text2id2))
        call this%outputtab%add_term(inum)
        call this%outputtab%add_term(cloc1)
        call this%outputtab%add_term(this%budobj%budterm(i)%id1(n))
        call this%outputtab%add_term(-this%mvrbudobj%budterm(i)%flow(n))
        call this%outputtab%add_term(this%budobj%budterm(i)%flow(n))
        call this%outputtab%add_term(cloc2)
        call this%outputtab%add_term(this%budobj%budterm(i)%id2(n))
        inum = inum + 1
      end do
    end do
  end subroutine mvt_print_outputtab

end module TspMvtModule

