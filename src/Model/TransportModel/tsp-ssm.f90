!> @brief This module contains the TspSsm Module
!!
!! This module contains the code for handling sources and sinks
!! associated with groundwater flow model stress packages.
!!
!! todo: need observations for SSM terms
!<
module TspSsmModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DONE, DZERO, LENAUXNAME, LENFTYPE, &
                             LENPACKAGENAME, LINELENGTH, &
                             TABLEFT, TABCENTER, LENBUDROWLABEL, LENVARNAME
  use SimModule, only: store_error, count_errors, store_error_unit
  use SimVariablesModule, only: errmsg
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use GweInputDataModule, only: GweInputDataType
  use TableModule, only: TableType, table_cr
  use TspSpcModule, only: TspSpcType
  use MatrixBaseModule

  implicit none
  public :: TspSsmType
  public :: ssm_cr

  character(len=LENFTYPE) :: ftype = 'SSM'
  character(len=LENPACKAGENAME) :: text = ' SOURCE-SINK MIX'

  !> @brief Derived type for the SSM Package
  !!
  !! This derived type corresponds to the SSM Package, which adds
  !! the effects of groundwater sources and sinks to the solute transport
  !! equation.
  !<
  type, extends(NumericalPackageType) :: TspSsmType

    integer(I4B), pointer :: nbound !< total number of flow boundaries in this time step
    integer(I4B), dimension(:), pointer, contiguous :: isrctype => null() !< source type 0 is unspecified, 1 is aux, 2 is auxmixed, 3 is ssmi, 4 is ssmimixed
    integer(I4B), dimension(:), pointer, contiguous :: iauxpak => null() !< aux col for concentration/temperature
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    real(DP), dimension(:), pointer, contiguous :: cnew => null() !< pointer to gwt%x
    type(TspFmiType), pointer :: fmi => null() !< pointer to fmi object
    type(TableType), pointer :: outputtab => null() !< output table object
    type(TspSpcType), dimension(:), pointer :: ssmivec => null() !< array of stress package concentration or temperature objects
    real(DP), pointer :: eqnsclfac => null() !< governing equation scale factor; =1. for solute; =rhow*cpw for energy
    character(len=LENVARNAME) :: depvartype = ''

  contains

    procedure :: ssm_df
    procedure :: ssm_ar
    procedure :: ssm_rp
    procedure :: ssm_ad
    procedure :: ssm_fc
    procedure :: ssm_cq
    procedure :: ssm_bd
    procedure :: ssm_ot_flow
    procedure :: ssm_da
    procedure :: allocate_scalars
    procedure, private :: ssm_term
    procedure, private :: allocate_arrays
    procedure, private :: read_options
    procedure, private :: read_data
    procedure, private :: read_sources_aux
    procedure, private :: read_sources_fileinput
    procedure, private :: pak_setup_outputtab
    procedure, private :: set_iauxpak
    procedure, private :: set_ssmivec
    procedure, private :: get_ssm_conc

  end type TspSsmType

contains

  !> @ brief Create a new SSM package
  !!
  !!  Create a new SSM package by defining names, allocating scalars
  !!  and initializing the parser.
  !<
  subroutine ssm_cr(ssmobj, name_model, inunit, iout, fmi, eqnsclfac, &
                    depvartype)
    ! -- dummy
    type(TspSsmType), pointer :: ssmobj !< TspSsmType object
    character(len=*), intent(in) :: name_model !< name of the model
    integer(I4B), intent(in) :: inunit !< fortran unit for input
    integer(I4B), intent(in) :: iout !< fortran unit for output
    type(TspFmiType), intent(in), target :: fmi !< Transport FMI package
    real(DP), intent(in), pointer :: eqnsclfac !< governing equation scale factor
    character(len=LENVARNAME), intent(in) :: depvartype !< dependent variable type ('concentration' or 'temperature')
    !
    ! -- Create the object
    allocate (ssmobj)
    !
    ! -- create name and memory path
    call ssmobj%set_names(1, name_model, 'SSM', 'SSM')
    !
    ! -- Allocate scalars
    call ssmobj%allocate_scalars()
    !
    ! -- Set variables
    ssmobj%inunit = inunit
    ssmobj%iout = iout
    ssmobj%fmi => fmi
    ssmobj%eqnsclfac => eqnsclfac
    !
    ! -- Initialize block parser
    call ssmobj%parser%Initialize(ssmobj%inunit, ssmobj%iout)
    !
    ! -- Store pointer to labels associated with the current model so that the
    !    package has access to the corresponding dependent variable type
    ssmobj%depvartype = depvartype
  end subroutine ssm_cr

  !> @ brief Define SSM Package
  !!
  !! This routine is called from gwt_df(), but does not do anything because
  !! df is typically used to set up dimensions.  For the ssm package, the
  !! total number of ssm entries is defined by the flow model.
  !<
  subroutine ssm_df(this)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
  end subroutine ssm_df

  !> @ brief Allocate and read SSM Package
  !!
  !! This routine is called from gwt_ar().  It allocates arrays, reads
  !! options and data, and sets up the output table.
  !<
  subroutine ssm_ar(this, dis, ibound, cnew)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
    class(DisBaseType), pointer, intent(in) :: dis !< discretization package
    integer(I4B), dimension(:), pointer, contiguous :: ibound !< GWT model ibound
    real(DP), dimension(:), pointer, contiguous :: cnew !< GWT model dependent variable
    ! -- formats
    character(len=*), parameter :: fmtssm = &
      "(1x,/1x,'SSM -- SOURCE-SINK MIXING PACKAGE, VERSION 1, 8/25/2017', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! -- print a message identifying the storage package.
    write (this%iout, fmtssm) this%inunit
    !
    ! -- store pointers to arguments that were passed in
    this%dis => dis
    this%ibound => ibound
    this%cnew => cnew
    !
    ! -- Check to make sure that there are flow packages
    if (this%fmi%nflowpack == 0) then
      write (errmsg, '(a)') 'SSM package does not detect any boundary flows &
                            &that require SSM terms.  Activate GWF-GWT (or &
                            &GWF-GWE, as appropriate) exchange or activate &
                            &FMI package and provide a budget file that &
                            &contains boundary flows.  If no boundary flows &
                            &are present in corresponding GWF model then this &
                            &SSM package should be removed.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- Allocate arrays
    call this%allocate_arrays()
    !
    ! -- Read ssm options
    call this%read_options()
    !
    ! -- read the data block
    call this%read_data()
    !
    ! -- setup the output table
    call this%pak_setup_outputtab()
  end subroutine ssm_ar

  !> @ brief Read and prepare this SSM Package
  !!
  !! This routine is called from gwt_rp().  It is called at the beginning of
  !! each stress period.  If any SPC input files are used to provide source
  !! and sink concentrations (or temperatures), then period blocks for the
  !! current stress period are read.
  !<
  subroutine ssm_rp(this)
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
    ! -- local
    integer(I4B) :: ip
    type(TspSpcType), pointer :: ssmiptr
    ! -- formats
    !
    ! -- Call the rp method on any ssm input files
    do ip = 1, this%fmi%nflowpack
      if (this%fmi%iatp(ip) /= 0) cycle
      if (this%isrctype(ip) == 3 .or. this%isrctype(ip) == 4) then
        ssmiptr => this%ssmivec(ip)
        call ssmiptr%spc_rp()
      end if
    end do
  end subroutine ssm_rp

  !> @ brief Advance the SSM Package
  !!
  !! This routine is called from gwt_ad().  It is called at the beginning of
  !! each time step.  The total number of flow boundaries is counted and stored
  !! in this%nbound.  Also, if any SPC input files are used to provide source
  !! and sink concentrations (or temperatures) and time series are referenced
  !! in those files, then ssm concenrations must be interpolated for the time
  !! step.
  !<
  subroutine ssm_ad(this)
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
    ! -- local
    integer(I4B) :: ip
    type(TspSpcType), pointer :: ssmiptr
    integer(I4B) :: i
    integer(I4B) :: node
    !
    ! -- Calculate total number of existing flow boundaries. It is possible
    !    that a node may equal zero.  In this case, the bound should be
    !    skipped and not written to ssm output.
    this%nbound = 0
    do ip = 1, this%fmi%nflowpack
      if (this%fmi%iatp(ip) /= 0) cycle
      do i = 1, this%fmi%gwfpackages(ip)%nbound
        node = this%fmi%gwfpackages(ip)%nodelist(i)
        if (node > 0) then
          this%nbound = this%nbound + 1
        end if
      end do
    end do
    !
    ! -- Call the ad method on any ssm input files so that values for
    !    time-series are interpolated
    do ip = 1, this%fmi%nflowpack
      if (this%fmi%iatp(ip) /= 0) cycle
      if (this%isrctype(ip) == 3 .or. this%isrctype(ip) == 4) then
        ssmiptr => this%ssmivec(ip)
        call ssmiptr%spc_ad(this%fmi%gwfpackages(ip)%nbound, &
                            this%fmi%gwfpackages(ip)%budtxt)
      end if
    end do
  end subroutine ssm_ad

  !> @ brief Calculate the SSM mass flow rate and hcof and rhs values
  !!
  !! This is the primary SSM routine that calculates the matrix coefficient
  !! and right-hand-side value for any package and package entry.  It returns
  !! several different optional variables that are used throughout this
  !! package to update matrix terms, budget calculations, and output tables.
  !<
  subroutine ssm_term(this, ipackage, ientry, rrate, rhsval, hcofval, &
                      cssm, qssm)
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType
    integer(I4B), intent(in) :: ipackage !< package number
    integer(I4B), intent(in) :: ientry !< bound number
    real(DP), intent(out), optional :: rrate !< calculated mass flow rate
    real(DP), intent(out), optional :: rhsval !< calculated rhs value
    real(DP), intent(out), optional :: hcofval !< calculated hcof value
    real(DP), intent(out), optional :: cssm !< calculated source concentration/temperature depending on flow direction
    real(DP), intent(out), optional :: qssm !< water flow rate into model cell from boundary package
    ! -- local
    logical(LGP) :: lauxmixed
    integer(I4B) :: n
    integer(I4B) :: nbound_flow
    real(DP) :: qbnd
    real(DP) :: ctmp
    real(DP) :: omega
    real(DP) :: hcoftmp
    real(DP) :: rhstmp
    !
    ! -- retrieve node number, qbnd and iauxpos
    hcoftmp = DZERO
    rhstmp = DZERO
    ctmp = DZERO
    qbnd = DZERO
    nbound_flow = this%fmi%gwfpackages(ipackage)%nbound
    n = this%fmi%gwfpackages(ipackage)%nodelist(ientry)
    !
    ! -- If cell is active (ibound > 0) then calculate values
    if (this%ibound(n) > 0) then
      !
      ! -- retrieve qbnd and iauxpos
      qbnd = this%fmi%gwfpackages(ipackage)%get_flow(ientry)
      call this%get_ssm_conc(ipackage, ientry, nbound_flow, ctmp, lauxmixed)
      !
      ! -- assign values for hcoftmp, rhstmp, and ctmp for subsequent assignment
      !    of hcof, rhs, and rate
      if (.not. lauxmixed) then
        !
        ! -- If qbnd is positive, then concentration (or temperature) represents
        !    the inflow concentration (or temperature).  If qbnd is negative,
        !    then the outflow concentration (or temperature) is set equal to the
        !    simulated cell's concentration (or temperature).
        if (qbnd >= DZERO) then
          omega = DZERO ! rhs
        else
          ctmp = this%cnew(n)
          omega = DONE ! lhs
          if (ctmp < DZERO) then
            omega = DZERO ! concentration/temperature is negative, so set mass flux to zero
          end if
        end if
      else
        !
        ! -- lauxmixed value indicates that this is a mixed sink type where
        !    the concentration (or temperature) value represents the injected
        !    concentration (or temperature) if qbnd is positive. If qbnd is
        !    negative, then the withdrawn water is equal to the minimum of the
        !    aux concentration (or temperature) and the cell concentration
        !    (or temperature).
        if (qbnd >= DZERO) then
          omega = DZERO ! rhs (ctmp is aux value)
        else
          if (ctmp < this%cnew(n)) then
            omega = DZERO ! rhs (ctmp is aux value)
          else
            omega = DONE ! lhs (ctmp is cell concentration)
            ctmp = this%cnew(n)
          end if
        end if
      end if
      !
      ! -- Add terms based on qbnd sign
      if (qbnd <= DZERO) then
        hcoftmp = qbnd * omega * this%eqnsclfac
      else
        rhstmp = -qbnd * ctmp * (DONE - omega) * this%eqnsclfac
      end if
      !
      ! -- end of active ibound
    end if
    !
    ! -- set requested values
    if (present(hcofval)) hcofval = hcoftmp
    if (present(rhsval)) rhsval = rhstmp
    if (present(rrate)) rrate = hcoftmp * ctmp - rhstmp
    if (present(cssm)) cssm = ctmp
    if (present(qssm)) qssm = qbnd
  end subroutine ssm_term

  !> @ brief Provide bound concentration (or temperature) and mixed flag
  !!
  !! SSM concentrations and temperatures can be provided in auxiliary variables
  !! or through separate SPC files.  If not provided, the default
  !! concentration (or temperature) is zero.  This single routine provides
  !! the SSM bound concentration (or temperature) based on these different
  !! approaches. The mixed flag indicates whether or not the boundary as a
  !! mixed type.
  !<
  subroutine get_ssm_conc(this, ipackage, ientry, nbound_flow, conc, &
                          lauxmixed)
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType
    integer(I4B), intent(in) :: ipackage !< package number
    integer(I4B), intent(in) :: ientry !< bound number
    integer(I4B), intent(in) :: nbound_flow !< size of flow package bound list
    real(DP), intent(out) :: conc !< user-specified concentration/temperature for this bound
    logical(LGP), intent(out) :: lauxmixed !< user-specified flag for marking this as a mixed boundary
    ! -- local
    integer(I4B) :: isrctype
    integer(I4B) :: iauxpos
    !
    conc = DZERO
    lauxmixed = .false.
    isrctype = this%isrctype(ipackage)
    !
    select case (isrctype)
    case (1, 2)
      iauxpos = this%iauxpak(ipackage)
      conc = this%fmi%gwfpackages(ipackage)%auxvar(iauxpos, ientry)
      if (isrctype == 2) lauxmixed = .true.
    case (3, 4)
      conc = this%ssmivec(ipackage)%get_value(ientry, nbound_flow)
      if (isrctype == 4) lauxmixed = .true.
    end select
  end subroutine get_ssm_conc

  !> @ brief Fill coefficients
  !!
  !! This routine adds the effects of the SSM to the matrix equations by
  !! updating the a matrix and right-hand side vector.
  !<
  subroutine ssm_fc(this, matrix_sln, idxglo, rhs)
    ! -- dummy
    class(TspSsmType) :: this
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    ! -- local
    integer(I4B) :: ip
    integer(I4B) :: i
    integer(I4B) :: n
    integer(I4B) :: idiag
    integer(I4B) :: nflowpack
    integer(I4B) :: nbound
    real(DP) :: hcofval
    real(DP) :: rhsval
    !
    ! -- do for each flow package
    nflowpack = this%fmi%nflowpack
    do ip = 1, nflowpack
      if (this%fmi%iatp(ip) /= 0) cycle
      !
      ! -- do for each entry in package (ip)
      nbound = this%fmi%gwfpackages(ip)%nbound
      do i = 1, nbound
        n = this%fmi%gwfpackages(ip)%nodelist(i)
        if (n <= 0) cycle
        call this%ssm_term(ip, i, rhsval=rhsval, hcofval=hcofval)
        idiag = idxglo(this%dis%con%ia(n))
        call matrix_sln%add_value_pos(idiag, hcofval)
        rhs(n) = rhs(n) + rhsval
        !
      end do
      !
    end do
  end subroutine ssm_fc

  !> @ brief Calculate flow
  !!
  !! Calculate the resulting mass flow between the boundary and the connected
  !! GWT/GWE model cell.  Update the diagonal position of the flowja array so
  !! that it ultimately contains the solute balance residual.
  !<
  subroutine ssm_cq(this, flowja)
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
    real(DP), dimension(:), contiguous, intent(inout) :: flowja !< flow across each face in the model grid
    ! -- local
    integer(I4B) :: ip
    integer(I4B) :: i
    integer(I4B) :: n
    integer(I4B) :: idiag
    real(DP) :: rate
    !
    ! -- do for each flow package
    do ip = 1, this%fmi%nflowpack
      !
      ! -- cycle if package is being managed as an advanced package
      if (this%fmi%iatp(ip) /= 0) cycle
      !
      ! -- do for each boundary
      do i = 1, this%fmi%gwfpackages(ip)%nbound
        n = this%fmi%gwfpackages(ip)%nodelist(i)
        if (n <= 0) cycle
        call this%ssm_term(ip, i, rrate=rate)
        idiag = this%dis%con%ia(n)
        flowja(idiag) = flowja(idiag) + rate
        !
      end do
      !
    end do
  end subroutine ssm_cq

  !> @ brief Calculate the global SSM budget terms
  !!
  !! Calculate the global SSM budget terms using separate in and out entries
  !! for each flow package.
  !<
  subroutine ssm_bd(this, isuppress_output, model_budget)
    ! -- modules
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
    integer(I4B), intent(in) :: isuppress_output !< flag to suppress output
    type(BudgetType), intent(inout) :: model_budget !< budget object for the GWT model
    ! -- local
    character(len=LENBUDROWLABEL) :: rowlabel
    integer(I4B) :: ip
    integer(I4B) :: i
    integer(I4B) :: n
    real(DP) :: rate
    real(DP) :: rin
    real(DP) :: rout
    !
    ! -- do for each flow package, unless it is being handled by an advanced
    !    transport package
    do ip = 1, this%fmi%nflowpack
      !
      ! -- cycle if package is being managed as an advanced package
      if (this%fmi%iatp(ip) /= 0) cycle
      !
      ! -- Initialize the rate accumulators
      rin = DZERO
      rout = DZERO
      !
      ! -- do for each boundary
      do i = 1, this%fmi%gwfpackages(ip)%nbound
        n = this%fmi%gwfpackages(ip)%nodelist(i)
        if (n <= 0) cycle
        call this%ssm_term(ip, i, rrate=rate)
        if (rate < DZERO) then
          rout = rout - rate
        else
          rin = rin + rate
        end if
        !
      end do
      !
      rowlabel = 'SSM_'//adjustl(trim(this%fmi%flowpacknamearray(ip)))
      call model_budget%addentry(rin, rout, delt, &
                                 this%fmi%gwfpackages(ip)%budtxt, &
                                 isuppress_output, rowlabel=rowlabel)
    end do
  end subroutine ssm_bd

  !> @ brief Output flows
  !!
  !! Based on user-specified controls, print SSM mass flow rates to the GWT
  !! listing file and/or write the SSM mass flow rates to the GWT binary
  !! budget file.
  !<
  subroutine ssm_ot_flow(this, icbcfl, ibudfl, icbcun)
    ! -- modules
    use TdisModule, only: kstp, kper
    use ConstantsModule, only: LENPACKAGENAME, LENBOUNDNAME, LENAUXNAME, DZERO
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
    integer(I4B), intent(in) :: icbcfl !< flag for writing binary budget terms
    integer(I4B), intent(in) :: ibudfl !< flag for printing budget terms to list file
    integer(I4B), intent(in) :: icbcun !< fortran unit number for binary budget file
    ! -- local
    character(len=LINELENGTH) :: title
    integer(I4B) :: node, nodeu
    character(len=20) :: nodestr
    integer(I4B) :: maxrows
    integer(I4B) :: ip
    integer(I4B) :: i, n2, ibinun
    real(DP) :: rrate
    real(DP) :: qssm
    real(DP) :: cssm
    integer(I4B) :: naux
    real(DP), dimension(0) :: auxrow
    character(len=LENAUXNAME), dimension(0) :: auxname
    ! -- for observations
    character(len=LENBOUNDNAME) :: bname
    ! -- formats
    character(len=*), parameter :: fmttkk = &
      &"(1X,/1X,A,'   PERIOD ',I0,'   STEP ',I0)"
    !
    ! -- initialize
    naux = 0
    maxrows = 0
    if (ibudfl /= 0 .and. this%iprflow /= 0) then
      call this%outputtab%set_kstpkper(kstp, kper)
      do ip = 1, this%fmi%nflowpack
        if (this%fmi%iatp(ip) /= 0) cycle
        !
        ! -- do for each boundary
        do i = 1, this%fmi%gwfpackages(ip)%nbound
          node = this%fmi%gwfpackages(ip)%nodelist(i)
          if (node > 0) then
            maxrows = maxrows + 1
          end if
        end do
      end do
      if (maxrows > 0) then
        call this%outputtab%set_maxbound(maxrows)
      end if
      title = 'SSM PACKAGE ('//trim(this%packName)// &
              ') FLOW RATES'
      call this%outputtab%set_title(title)
    end if
    !
    ! -- Set unit number for binary output
    if (this%ipakcb < 0) then
      ibinun = icbcun
    else if (this%ipakcb == 0) then
      ibinun = 0
    else
      ibinun = this%ipakcb
    end if
    if (icbcfl == 0) ibinun = 0
    !
    ! -- If cell-by-cell flows will be saved as a list, write header.
    if (ibinun /= 0) then
      call this%dis%record_srcdst_list_header(text, this%name_model, &
                                              this%name_model, this%name_model, &
                                              this%packName, naux, auxname, &
                                              ibinun, this%nbound, this%iout)
    end if
    !
    ! -- If no boundaries, skip flow calculations.
    if (this%nbound > 0) then
      !
      ! -- Loop through each boundary calculating flow.
      do ip = 1, this%fmi%nflowpack
        if (this%fmi%iatp(ip) /= 0) cycle
        !
        ! -- do for each boundary
        do i = 1, this%fmi%gwfpackages(ip)%nbound
          !
          ! -- Calculate rate for this entry
          node = this%fmi%gwfpackages(ip)%nodelist(i)
          if (node <= 0) cycle
          call this%ssm_term(ip, i, rrate=rrate, qssm=qssm, cssm=cssm)
          !
          ! -- Print the individual rates if the budget is being printed
          !    and PRINT_FLOWS was specified (this%iprflow<0)
          if (ibudfl /= 0) then
            if (this%iprflow /= 0) then
              !
              ! -- set nodestr and write outputtab table
              nodeu = this%dis%get_nodeuser(node)
              call this%dis%nodeu_to_string(nodeu, nodestr)
              bname = this%fmi%gwfpackages(ip)%name
              call this%outputtab%add_term(i)
              call this%outputtab%add_term(trim(adjustl(nodestr)))
              call this%outputtab%add_term(qssm)
              call this%outputtab%add_term(cssm)
              call this%outputtab%add_term(rrate)
              call this%outputtab%add_term(bname)
            end if
          end if
          !
          ! -- If saving cell-by-cell flows in list, write flow
          if (ibinun /= 0) then
            n2 = i
            call this%dis%record_mf6_list_entry(ibinun, node, n2, rrate, &
                                                naux, auxrow, &
                                                olconv2=.FALSE.)
          end if
          !
        end do
        !
      end do
    end if
    if (ibudfl /= 0) then
      if (this%iprflow /= 0) then
        write (this%iout, '(1x)')
      end if
    end if
  end subroutine ssm_ot_flow

  !> @ brief Deallocate
  !!
  !! Deallocate the memory associated with this derived type
  !<
  subroutine ssm_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
    ! -- local
    integer(I4B) :: ip
    type(TspSpcType), pointer :: ssmiptr
    !
    ! -- Deallocate the ssmi objects if package was active
    if (this%inunit > 0) then
      do ip = 1, size(this%ssmivec)
        if (this%isrctype(ip) == 3 .or. this%isrctype(ip) == 4) then
          ssmiptr => this%ssmivec(ip)
          call ssmiptr%spc_da()
        end if
      end do
      deallocate (this%ssmivec)
    end if
    !
    ! -- Deallocate arrays if package was active
    if (this%inunit > 0) then
      call mem_deallocate(this%iauxpak)
      call mem_deallocate(this%isrctype)
      this%ibound => null()
      this%fmi => null()
    end if
    !
    ! -- output table object
    if (associated(this%outputtab)) then
      call this%outputtab%table_da()
      deallocate (this%outputtab)
      nullify (this%outputtab)
    end if
    !
    ! -- Scalars
    call mem_deallocate(this%nbound)
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()
  end subroutine ssm_da

  !> @ brief Allocate scalars
  !!
  !! Allocate scalar variables for this derived type
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
    !
    ! -- allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%nbound, 'NBOUND', this%memoryPath)
    !
    ! -- Initialize
    this%nbound = 0
  end subroutine allocate_scalars

  !> @ brief Allocate arrays
  !!
  !! Allocate array variables for this derived type
  !<
  subroutine allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
    ! -- local
    integer(I4B) :: nflowpack
    integer(I4B) :: i
    !
    ! -- Allocate
    nflowpack = this%fmi%nflowpack
    call mem_allocate(this%iauxpak, nflowpack, 'IAUXPAK', this%memoryPath)
    call mem_allocate(this%isrctype, nflowpack, 'ISRCTYPE', this%memoryPath)
    !
    ! -- Initialize
    do i = 1, nflowpack
      this%iauxpak(i) = 0
      this%isrctype(i) = 0
    end do
    !
    ! -- Allocate the ssmivec array
    allocate (this%ssmivec(nflowpack))
  end subroutine allocate_arrays

  !> @ brief Read package options
  !!
  !! Read and set the SSM Package options
  !<
  subroutine read_options(this)
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
    ! -- local
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    ! -- formats
    character(len=*), parameter :: fmtiprflow = &
      "(4x,'SSM FLOW INFORMATION WILL BE PRINTED TO LISTING FILE &
      &WHENEVER ICBCFL IS NOT ZERO.')"
    character(len=*), parameter :: fmtisvflow = &
      "(4x,'CELL-BY-CELL FLOW INFORMATION WILL BE SAVED TO BINARY FILE &
      &WHENEVER ICBCFL IS NOT ZERO.')"
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, blockRequired=.false., &
                              supportOpenClose=.true.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING SSM OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('PRINT_FLOWS')
          this%iprflow = 1
          write (this%iout, fmtiprflow)
        case ('SAVE_FLOWS')
          this%ipakcb = -1
          write (this%iout, fmtisvflow)
        case default
          write (errmsg, '(a,a)') 'Unknown SSM option: ', trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (this%iout, '(1x,a)') 'END OF SSM OPTIONS'
    end if
  end subroutine read_options

  !> @ brief Read package data
  !!
  !! Read and set the SSM Package data
  !<
  subroutine read_data(this)
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
    !
    ! -- read and process required SOURCES block
    call this%read_sources_aux()
    !
    ! -- read and process optional FILEINPUT block
    call this%read_sources_fileinput()
  end subroutine read_data

  !> @ brief Read SOURCES block
  !!
  !! Read SOURCES block and look for auxiliary columns in
  !! corresponding flow data.
  !<
  subroutine read_sources_aux(this)
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
    ! -- local
    character(len=LINELENGTH) :: keyword
    character(len=20) :: srctype
    integer(I4B) :: ierr
    integer(I4B) :: ip
    integer(I4B) :: nflowpack
    integer(I4B) :: isrctype
    logical :: isfound, endOfBlock
    logical :: pakfound
    logical :: lauxmixed
    !
    ! -- initialize
    isfound = .false.
    lauxmixed = .false.
    nflowpack = this%fmi%nflowpack
    !
    ! -- get sources block
    call this%parser%GetBlock('SOURCES', isfound, ierr, &
                              supportOpenClose=.true., &
                              blockrequired=.true.)
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING SOURCES'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        !
        ! -- read package name and make sure it can be found
        call this%parser%GetStringCaps(keyword)
        pakfound = .false.
        do ip = 1, nflowpack
          if (trim(adjustl(this%fmi%gwfpackages(ip)%name)) == keyword) then
            pakfound = .true.
            exit
          end if
        end do
        if (.not. pakfound) then
          write (errmsg, '(a,a)') 'Flow package cannot be found: ', &
            trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
        !
        ! -- Ensure package was not specified more than once in SOURCES block
        if (this%isrctype(ip) /= 0) then
          write (errmsg, '(a, a)') &
            'A package cannot be specified more than once in the SSM SOURCES &
            &block.  The following package was specified more than once: ', &
            trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
        !
        ! -- read the source type
        call this%parser%GetStringCaps(srctype)
        select case (srctype)
        case ('AUX')
          write (this%iout, '(1x,a)') 'AUX SOURCE DETECTED.'
          isrctype = 1
        case ('AUXMIXED')
          write (this%iout, '(1x,a)') 'AUXMIXED SOURCE DETECTED.'
          lauxmixed = .true.
          isrctype = 2
        case default
          write (errmsg, '(a, a)') &
            'SRCTYPE must be AUX or AUXMIXED.  Found: ', trim(srctype)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
        !
        ! -- Store the source type (1 or 2)
        this%isrctype(ip) = isrctype
        !
        ! -- Find and store the auxiliary column
        call this%set_iauxpak(ip, trim(keyword))

      end do
      write (this%iout, '(1x,a)') 'END PROCESSING SOURCES'
    else
      write (errmsg, '(a)') 'Required SOURCES block not found.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- terminate if errors
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine read_sources_aux

  !> @ brief Read FILEINPUT block
  !!
  !! Read optional FILEINPUT block and initialize an
  !! SPC input file reader for each entry.
  !<
  subroutine read_sources_fileinput(this)
    ! -- dummy
    class(TspSsmType) :: this !< TspSsmType object
    ! -- local
    character(len=LINELENGTH) :: keyword
    character(len=LINELENGTH) :: keyword2
    character(len=20) :: srctype
    integer(I4B) :: ierr
    integer(I4B) :: ip
    integer(I4B) :: nflowpack
    integer(I4B) :: isrctype
    logical :: isfound, endOfBlock
    logical :: pakfound
    logical :: lauxmixed
    !
    ! -- initialize
    isfound = .false.
    lauxmixed = .false.
    nflowpack = this%fmi%nflowpack
    !
    ! -- get sources_file block
    call this%parser%GetBlock('FILEINPUT', isfound, ierr, &
                              supportOpenClose=.true., &
                              blockrequired=.false.)
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING FILEINPUT'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        !
        ! -- read package name and make sure it can be found
        call this%parser%GetStringCaps(keyword)
        pakfound = .false.
        do ip = 1, nflowpack
          if (trim(adjustl(this%fmi%gwfpackages(ip)%name)) == keyword) then
            pakfound = .true.
            exit
          end if
        end do
        if (.not. pakfound) then
          write (errmsg, '(a,a)') 'Flow package cannot be found: ', &
            trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
        !
        ! -- Ensure package was not specified more than once in SOURCES block
        if (this%isrctype(ip) /= 0) then
          write (errmsg, '(a, a)') &
            'A package cannot be specified more than once in the SSM SOURCES &
            &and SOURCES_FILES blocks.  The following package was specified &
            &more than once: ', &
            trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
        !
        ! -- read the source type
        call this%parser%GetStringCaps(srctype)
        select case (srctype)
        case ('SPC6')
          write (this%iout, '(1x,a)') 'SPC6 SOURCE DETECTED.'
          isrctype = 3
          !
          ! verify filein is next
          call this%parser%GetStringCaps(keyword2)
          if (trim(adjustl(keyword2)) /= 'FILEIN') then
            errmsg = 'SPC6 keyword must be followed by "FILEIN" '// &
                     'then by filename and optionally by <MIXED>.'
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end if
          !
          ! -- Use set_ssmivec to read file name and set up
          !    ssmi file object
          call this%set_ssmivec(ip, trim(keyword))
          !
          ! -- check for optional MIXED keyword and set isrctype to 4 if found
          call this%parser%GetStringCaps(keyword2)
          if (trim(keyword2) == 'MIXED') then
            isrctype = 4
            write (this%iout, '(1x,a,a)') 'ASSIGNED MIXED SSM TYPE TO PACKAGE ', &
              trim(keyword)
          end if
        case default
          write (errmsg, '(a,a)') &
            'SRCTYPE must be SPC6.  Found: ', trim(srctype)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
        !
        ! -- Store the source type (3 or 4)
        this%isrctype(ip) = isrctype
        !
      end do
      write (this%iout, '(1x,a)') 'END PROCESSING FILEINPUT'
    else
      write (this%iout, '(1x,a)') &
        'OPTIONAL FILEINPUT BLOCK NOT FOUND.  CONTINUING.'
    end if
    !
    ! -- terminate if errors
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine read_sources_fileinput

  !> @ brief Set iauxpak array value for package ip
  !!
  !!  The next call to parser will return the auxiliary name for
  !!  package ip in the SSM SOURCES block.  The routine searches
  !!  through the auxiliary names in package ip and sets iauxpak
  !!  to the column number corresponding to the correct auxiliary
  !!  column.
  !<
  subroutine set_iauxpak(this, ip, packname)
    ! -- dummy
    class(TspSsmType), intent(inout) :: this !< TspSsmType
    integer(I4B), intent(in) :: ip !< package number
    character(len=*), intent(in) :: packname !< name of package
    ! -- local
    character(len=LENAUXNAME) :: auxname
    logical :: auxfound
    integer(I4B) :: iaux
    !
    ! -- read name of auxiliary column
    call this%parser%GetStringCaps(auxname)
    auxfound = .false.
    do iaux = 1, this%fmi%gwfpackages(ip)%naux
      if (trim(this%fmi%gwfpackages(ip)%auxname(iaux)) == &
          trim(auxname)) then
        auxfound = .true.
        exit
      end if
    end do
    if (.not. auxfound) then
      write (errmsg, '(a, a)') &
        'Auxiliary name cannot be found: ', trim(auxname)
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- set iauxpak and write message
    this%iauxpak(ip) = iaux
    write (this%iout, '(4x, a, i0, a, a)') 'USING AUX COLUMN ', &
      iaux, ' IN PACKAGE ', trim(packname)
  end subroutine set_iauxpak

  !> @ brief Set ssmivec array value for package ip
  !!
  !!  The next call to parser will return the input file name for
  !!  package ip in the SSM SOURCES block.  The routine then
  !!  initializes the SPC input file.
  !<
  subroutine set_ssmivec(this, ip, packname)
    ! -- module
    use InputOutputModule, only: openfile, getunit
    ! -- dummy
    class(TspSsmType), intent(inout) :: this !< TspSsmType
    integer(I4B), intent(in) :: ip !< package number
    character(len=*), intent(in) :: packname !< name of package
    ! -- local
    character(len=LINELENGTH) :: filename
    type(TspSpcType), pointer :: ssmiptr
    integer(I4B) :: inunit
    !
    ! -- read file name
    call this%parser%GetString(filename)
    inunit = getunit()
    call openfile(inunit, this%iout, filename, 'SPC', filstat_opt='OLD')
    !
    ! -- Create the SPC file object
    ssmiptr => this%ssmivec(ip)
    call ssmiptr%initialize(this%dis, ip, inunit, this%iout, this%name_model, &
                            trim(packname), this%depvartype)
    !
    write (this%iout, '(4x, a, a, a, a, a)') 'USING SPC INPUT FILE ', &
      trim(filename), ' TO SET ', trim(this%depvartype), &
      'S FOR PACKAGE ', trim(packname)
  end subroutine set_ssmivec

  !> @ brief Setup the output table
  !!
  !! Setup the output table by creating the column headers.
  !<
  subroutine pak_setup_outputtab(this)
    ! -- dummy
    class(TspSsmType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    integer(I4B) :: ntabcol
    !
    ! -- allocate and initialize the output table
    if (this%iprflow /= 0) then
      !
      ! -- dimension table
      ntabcol = 6
      !if (this%inamedbound > 0) then
      !  ntabcol = ntabcol + 1
      !end if
      !
      ! -- initialize the output table object
      title = 'SSM PACKAGE ('//trim(this%packName)// &
              ') FLOW RATES'
      call table_cr(this%outputtab, this%packName, title)
      call this%outputtab%table_df(1, ntabcol, this%iout, transient=.TRUE.)
      text = 'NUMBER'
      call this%outputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'CELLID'
      call this%outputtab%initialize_column(text, 20, alignment=TABLEFT)
      text = 'BOUND Q'
      call this%outputtab%initialize_column(text, 15, alignment=TABCENTER)
      text = 'SSM CONC'
      call this%outputtab%initialize_column(text, 15, alignment=TABCENTER)
      text = 'RATE'
      call this%outputtab%initialize_column(text, 15, alignment=TABCENTER)
      text = 'PACKAGE NAME'
      call this%outputtab%initialize_column(text, 16, alignment=TABCENTER)
      !if (this%inamedbound > 0) then
      !  text = 'NAME'
      !  call this%outputtab%initialize_column(text, 20, alignment=TABLEFT)
      !end if
    end if
  end subroutine pak_setup_outputtab

end module TspSsmModule
