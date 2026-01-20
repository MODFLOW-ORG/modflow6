module TspFmiModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DONE, DZERO, DHALF, LINELENGTH, LENBUDTXT, &
                             LENPACKAGENAME, LENVARNAME
  use SimModule, only: store_error, count_errors, store_error_filename
  use SimVariablesModule, only: errmsg
  use FlowModelInterfaceModule, only: FlowModelInterfaceType
  use BaseDisModule, only: DisBaseType
  use ListModule, only: ListType
  use BudgetFileReaderModule, only: BudgetFileReaderType
  use HeadFileReaderModule, only: HeadFileReaderType
  use PackageBudgetModule, only: PackageBudgetType
  use BudgetObjectModule, only: BudgetObjectType, budgetobject_cr_bfr
  use MatrixBaseModule

  implicit none
  private
  public :: TspFmiType
  public :: fmi_cr

  character(len=LENPACKAGENAME) :: text = '    GWTFMI'

  integer(I4B), parameter :: NBDITEMS = 2
  character(len=LENBUDTXT), dimension(NBDITEMS) :: budtxt
  data budtxt/'      FLOW-ERROR', ' FLOW-CORRECTION'/

  type :: DataAdvancedPackageType
    real(DP), dimension(:), contiguous, pointer :: concpack => null()
    real(DP), dimension(:), contiguous, pointer :: qmfrommvr => null()
  end type

  type :: BudObjPtrArray
    type(BudgetObjectType), pointer :: ptr
  end type BudObjPtrArray

  type, extends(FlowModelInterfaceType) :: TspFmiType

    integer(I4B), dimension(:), pointer, contiguous :: iatp => null() !< advanced transport package applied to gwfpackages
    integer(I4B), pointer :: iflowerr => null() !< add the flow error correction
    real(DP), dimension(:), pointer, contiguous :: flowcorrect => null() !< mass flow correction
    real(DP), pointer :: eqnsclfac => null() !< governing equation scale factor; =1. for solute; =rhow*cpw for energy
    type(DataAdvancedPackageType), &
      dimension(:), pointer, contiguous :: datp => null()
    type(BudObjPtrArray), dimension(:), allocatable :: aptbudobj !< flow budget objects for the advanced packages

  contains

    procedure :: allocate_arrays => gwtfmi_allocate_arrays
    procedure :: allocate_gwfpackages => gwtfmi_allocate_gwfpackages
    procedure :: allocate_scalars => gwtfmi_allocate_scalars
    procedure :: deallocate_gwfpackages => gwtfmi_deallocate_gwfpackages
    procedure :: fmi_rp
    procedure :: fmi_ad
    procedure :: fmi_fc
    procedure :: fmi_cq
    procedure :: fmi_bd
    procedure :: fmi_ot_flow
    procedure :: fmi_da => gwtfmi_da
    procedure :: gwfsatold
    procedure :: initialize_gwfterms_from_bfr
    procedure :: initialize_gwfterms_from_gwfbndlist
    procedure :: source_options => gwtfmi_source_options
    procedure :: set_aptbudobj_pointer
    procedure :: source_packagedata => gwtfmi_source_packagedata
    procedure :: set_active_status

  end type TspFmiType

contains

  !> @brief Create a new FMI object
  !<
  subroutine fmi_cr(fmiobj, name_model, input_mempath, inunit, iout, eqnsclfac, &
                    depvartype)
    ! -- dummy
    type(TspFmiType), pointer :: fmiobj
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    real(DP), intent(in), pointer :: eqnsclfac !< governing equation scale factor
    character(len=LENVARNAME), intent(in) :: depvartype
    !
    ! -- Create the object
    allocate (fmiobj)
    !
    ! -- create name and memory path
    call fmiobj%set_names(1, name_model, 'FMI', 'FMI', input_mempath)
    fmiobj%text = text
    !
    ! -- Allocate scalars
    call fmiobj%allocate_scalars()
    !
    ! -- Set variables
    fmiobj%inunit = inunit
    fmiobj%iout = iout
    !
    ! -- Assign label based on dependent variable
    fmiobj%depvartype = depvartype
    !
    ! -- Store pointer to governing equation scale factor
    fmiobj%eqnsclfac => eqnsclfac
  end subroutine fmi_cr

  !> @brief Read and prepare
  !<
  subroutine fmi_rp(this, inmvr)
    ! -- modules
    use TdisModule, only: kper, kstp
    ! -- dummy
    class(TspFmiType) :: this
    integer(I4B), intent(in) :: inmvr
    ! -- local
    ! -- formats
    !
    ! --Check to make sure MVT Package is active if mvr flows are available.
    !   This cannot be checked until RP because exchange doesn't set a pointer
    !   to mvrbudobj until exg_ar().
    if (kper * kstp == 1) then
      if (associated(this%mvrbudobj) .and. inmvr == 0) then
        write (errmsg, '(a)') 'GWF water mover is active but the GWT MVT &
          &package has not been specified.  activate GWT MVT package.'
        call store_error(errmsg, terminate=.TRUE.)
      end if
      if (.not. associated(this%mvrbudobj) .and. inmvr > 0) then
        write (errmsg, '(a)') 'GWF water mover terms are not available &
          &but the GWT MVT package has been activated.  Activate GWF-GWT &
          &exchange or specify GWFMOVER in FMI PACKAGEDATA.'
        call store_error(errmsg, terminate=.TRUE.)
      end if
    end if
  end subroutine fmi_rp

  !> @brief Advance routine for FMI object
  !<
  subroutine fmi_ad(this, cnew)
    ! -- modules
    use ConstantsModule, only: DHDRY
    ! -- dummy
    class(TspFmiType) :: this
    real(DP), intent(inout), dimension(:) :: cnew
    ! -- local
    integer(I4B) :: n
    character(len=*), parameter :: fmtdry = &
     &"(/1X,'WARNING: DRY CELL ENCOUNTERED AT ',a,';  RESET AS INACTIVE &
     &WITH DRY CONCENTRATION = ', G13.5)"
    character(len=*), parameter :: fmtrewet = &
     &"(/1X,'DRY CELL REACTIVATED AT ', a,&
     &' WITH STARTING CONCENTRATION =',G13.5)"
    !
    ! -- Set flag to indicated that flows are being updated.  For the case where
    !    flows may be reused (only when flows are read from a file) then set
    !    the flag to zero to indicate that flows were not updated
    this%iflowsupdated = 1
    !
    ! -- If reading flows from a budget file, read the next set of records
    if (this%iubud /= 0) then
      call this%advance_bfr()
    end if
    !
    ! -- If reading heads from a head file, read the next set of records
    if (this%iuhds /= 0) then
      call this%advance_hfr()
    end if
    !
    ! -- If mover flows are being read from file, read the next set of records
    if (this%iumvr /= 0) then
      call this%mvrbudobj%bfr_advance(this%dis, this%iout)
    end if
    !
    ! -- If advanced package flows are being read from file, read the next set of records
    if (this%flows_from_file .and. this%inunit /= 0) then
      do n = 1, size(this%aptbudobj)
        call this%aptbudobj(n)%ptr%bfr_advance(this%dis, this%iout)
      end do
    end if
    !
    ! -- set inactive transport cell status
    if (this%idryinactive /= 0) then
      call this%set_active_status(cnew)
    end if
  end subroutine fmi_ad

  !> @brief Calculate coefficients and fill matrix and rhs terms associated
  !! with FMI object
  !<
  subroutine fmi_fc(this, nodes, cold, nja, matrix_sln, idxglo, rhs)
    ! -- dummy
    class(TspFmiType) :: this
    integer, intent(in) :: nodes
    real(DP), intent(in), dimension(nodes) :: cold
    integer(I4B), intent(in) :: nja
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(nja) :: idxglo
    real(DP), intent(inout), dimension(nodes) :: rhs
    ! -- local
    integer(I4B) :: n, idiag, idiag_sln
    real(DP) :: qcorr
    !
    ! -- Calculate the flow imbalance error and make a correction for it
    if (this%iflowerr /= 0) then
      !
      ! -- Correct the transport solution for the flow imbalance by adding
      !    the flow residual to the diagonal
      do n = 1, nodes
        idiag = this%dis%con%ia(n)
        idiag_sln = idxglo(idiag)
        !call matrix_sln%add_value_pos(idiag_sln, -this%gwfflowja(idiag))
        qcorr = -this%gwfflowja(idiag) * this%eqnsclfac
        call matrix_sln%add_value_pos(idiag_sln, qcorr)
      end do
    end if
  end subroutine fmi_fc

  !> @brief Calculate flow correction
  !!
  !! Where there is a flow imbalance for a given cell, a correction may be
  !! applied if selected
  !<
  subroutine fmi_cq(this, cnew, flowja)
    ! -- modules
    ! -- dummy
    class(TspFmiType) :: this
    real(DP), intent(in), dimension(:) :: cnew
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: idiag
    real(DP) :: rate
    !
    ! -- If not adding flow error correction, return
    if (this%iflowerr /= 0) then
      !
      ! -- Accumulate the flow correction term
      do n = 1, this%dis%nodes
        rate = DZERO
        idiag = this%dis%con%ia(n)
        if (this%ibound(n) > 0) then
          rate = -this%gwfflowja(idiag) * cnew(n) * this%eqnsclfac
        end if
        this%flowcorrect(n) = rate
        flowja(idiag) = flowja(idiag) + rate
      end do
    end if
  end subroutine fmi_cq

  !> @brief Calculate budget terms associated with FMI object
  !<
  subroutine fmi_bd(this, isuppress_output, model_budget)
    ! -- modules
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType, rate_accumulator
    ! -- dummy
    class(TspFmiType) :: this
    integer(I4B), intent(in) :: isuppress_output
    type(BudgetType), intent(inout) :: model_budget
    ! -- local
    real(DP) :: rin
    real(DP) :: rout
    !
    ! -- flow correction
    if (this%iflowerr /= 0) then
      call rate_accumulator(this%flowcorrect, rin, rout)
      call model_budget%addentry(rin, rout, delt, budtxt(2), isuppress_output)
    end if
  end subroutine fmi_bd

  !> @brief Save budget terms associated with FMI object
  !<
  subroutine fmi_ot_flow(this, icbcfl, icbcun)
    ! -- dummy
    class(TspFmiType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: icbcun
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
    ! -- Do not save flow corrections if not active
    if (this%iflowerr == 0) ibinun = 0
    !
    ! -- Record the storage rates if requested
    if (ibinun /= 0) then
      iprint = 0
      dinact = DZERO
      !
      ! -- flow correction
      call this%dis%record_array(this%flowcorrect, this%iout, iprint, -ibinun, &
                                 budtxt(2), cdatafmp, nvaluesp, &
                                 nwidthp, editdesc, dinact)
    end if
  end subroutine fmi_ot_flow

  !> @brief Deallocate variables
  !!
  !! Deallocate memory associated with FMI object
  !<
  subroutine gwtfmi_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(TspFmiType) :: this
    ! -- todo: finalize hfr and bfr either here or in a finalize routine
    !
    ! -- deallocate any memory stored with gwfpackages
    call this%deallocate_gwfpackages()
    !
    ! -- deallocate fmi arrays
    if (associated(this%datp)) then
      deallocate (this%datp)
      deallocate (this%gwfpackages)
      deallocate (this%flowpacknamearray)
      call mem_deallocate(this%iatp)
      call mem_deallocate(this%igwfmvrterm)
    end if

    deallocate (this%aptbudobj)
    call mem_deallocate(this%flowcorrect)
    call mem_deallocate(this%ibdgwfsat0)
    if (this%flows_from_file) then
      call mem_deallocate(this%gwfstrgss)
      call mem_deallocate(this%gwfstrgsy)
    end if
    !
    ! -- special treatment, these could be from mem_checkin
    call mem_deallocate(this%gwfhead, 'GWFHEAD', this%memoryPath)
    call mem_deallocate(this%gwfsat, 'GWFSAT', this%memoryPath)
    call mem_deallocate(this%gwfspdis, 'GWFSPDIS', this%memoryPath)
    call mem_deallocate(this%gwfflowja, 'GWFFLOWJA', this%memoryPath)
    !
    ! -- deallocate scalars
    call mem_deallocate(this%flows_from_file)
    call mem_deallocate(this%iflowsupdated)
    call mem_deallocate(this%iflowerr)
    call mem_deallocate(this%igwfstrgss)
    call mem_deallocate(this%igwfstrgsy)
    call mem_deallocate(this%iubud)
    call mem_deallocate(this%iuhds)
    call mem_deallocate(this%iumvr)
    call mem_deallocate(this%iugrb)
    call mem_deallocate(this%nflowpack)
    call mem_deallocate(this%idryinactive)
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()
  end subroutine gwtfmi_da

  !> @ brief Allocate scalars
  !!
  !! Allocate scalar variables for an FMI object
  !<
  subroutine gwtfmi_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy
    class(TspFmiType) :: this
    ! -- local
    !
    ! -- allocate scalars in parent
    call this%FlowModelInterfaceType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%iflowerr, 'IFLOWERR', this%memoryPath)
    !
    ! -- Although not a scalar, allocate the advanced package transport
    !    budget object to zero so that it can be dynamically resized later
    allocate (this%aptbudobj(0))
    !
    ! -- Initialize
    this%iflowerr = 0
  end subroutine gwtfmi_allocate_scalars

  !> @ brief Allocate arrays for FMI object
  !!
  !!  Method to allocate arrays for the FMI package.
  !<
  subroutine gwtfmi_allocate_arrays(this, nodes)
    use MemoryManagerModule, only: mem_allocate
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy
    class(TspFmiType) :: this
    integer(I4B), intent(in) :: nodes
    ! -- local
    integer(I4B) :: n
    !
    ! -- allocate parent arrays
    call this%FlowModelInterfaceType%allocate_arrays(nodes)
    !
    ! -- Allocate variables needed for all cases
    if (this%iflowerr == 0) then
      call mem_allocate(this%flowcorrect, 1, 'FLOWCORRECT', this%memoryPath)
    else
      call mem_allocate(this%flowcorrect, nodes, 'FLOWCORRECT', this%memoryPath)
    end if
    do n = 1, size(this%flowcorrect)
      this%flowcorrect(n) = DZERO
    end do
  end subroutine gwtfmi_allocate_arrays

  !> @brief Set gwt transport cell status
  !!
  !! Dry GWF cells are treated differently by GWT and GWE.  Transport does not
  !! occur in deactivated GWF cells; however, GWE still simulates conduction
  !! through dry cells.
  !<
  subroutine set_active_status(this, cnew)
    ! -- modules
    use ConstantsModule, only: DHDRY
    ! -- dummy
    class(TspFmiType) :: this
    real(DP), intent(inout), dimension(:) :: cnew
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: m
    integer(I4B) :: ipos
    real(DP) :: crewet, tflow, flownm
    character(len=15) :: nodestr
    ! -- formats
    character(len=*), parameter :: fmtoutmsg1 = &
      "(1x,'WARNING: DRY CELL ENCOUNTERED AT ', a,'; RESET AS INACTIVE WITH &
      &DRY ', a, '=', G13.5)"
    character(len=*), parameter :: fmtoutmsg2 = &
      &"(1x,'DRY CELL REACTIVATED AT', a, 'WITH STARTING', a, '=', G13.5)"
    !
    do n = 1, this%dis%nodes
      ! -- Calculate the ibound-like array that has 0 if saturation
      !    is zero and 1 otherwise
      if (this%gwfsat(n) > DZERO) then
        this%ibdgwfsat0(n) = 1
      else
        this%ibdgwfsat0(n) = 0
      end if
      !
      ! -- Check if active transport cell is inactive for flow
      if (this%ibound(n) > 0) then
        if (this%gwfhead(n) == DHDRY) then
          ! -- transport cell should be made inactive
          this%ibound(n) = 0
          cnew(n) = DHDRY
          call this%dis%noder_to_string(n, nodestr)
          write (this%iout, fmtoutmsg1) &
            trim(nodestr), trim(adjustl(this%depvartype)), DHDRY
        end if
      end if
    end do
    !
    ! -- if flow cell is dry, then set gwt%ibound = 0 and conc to dry
    do n = 1, this%dis%nodes
      !
      ! -- Convert dry transport cell to active if flow has rewet
      if (cnew(n) == DHDRY) then
        if (this%gwfhead(n) /= DHDRY) then
          !
          ! -- obtain weighted concentration/temperature
          crewet = DZERO
          tflow = DZERO
          do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
            m = this%dis%con%ja(ipos)
            flownm = this%gwfflowja(ipos)
            if (flownm > 0) then
              if (this%ibound(m) /= 0) then
                crewet = crewet + cnew(m) * flownm ! kluge note: apparently no need to multiply flows by eqnsclfac
                tflow = tflow + this%gwfflowja(ipos) !             since it will divide out below anyway
              end if
            end if
          end do
          if (tflow > DZERO) then
            crewet = crewet / tflow
          else
            crewet = DZERO
          end if
          !
          ! -- cell is now wet
          this%ibound(n) = 1
          cnew(n) = crewet
          call this%dis%noder_to_string(n, nodestr)
          write (this%iout, fmtoutmsg2) &
            trim(nodestr), trim(adjustl(this%depvartype)), crewet
        end if
      end if
    end do
  end subroutine set_active_status

  !> @brief Calculate the previous saturation level
  !!
  !! Calculate the groundwater cell head saturation for the end of
  !! the last time step
  !<
  function gwfsatold(this, n, delt) result(satold)
    ! -- modules
    ! -- dummy
    class(TspFmiType) :: this
    integer(I4B), intent(in) :: n
    real(DP), intent(in) :: delt
    ! -- result
    real(DP) :: satold
    ! -- local
    real(DP) :: vcell
    real(DP) :: vnew
    real(DP) :: vold
    !
    ! -- calculate the value
    vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
    vnew = vcell * this%gwfsat(n)
    vold = vnew
    if (this%igwfstrgss /= 0) vold = vold + this%gwfstrgss(n) * delt
    if (this%igwfstrgsy /= 0) vold = vold + this%gwfstrgsy(n) * delt
    satold = vold / vcell
  end function gwfsatold

  !> @ brief Source input options for package
  !<
  subroutine gwtfmi_source_options(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy
    class(TspFmiType) :: this
    ! -- local
    logical(LGP) :: found_ipakcb, found_flowerr
    character(len=*), parameter :: fmtisvflow = &
      "(4x,'CELL-BY-CELL FLOW INFORMATION WILL BE SAVED TO BINARY FILE &
      &WHENEVER ICBCFL IS NOT ZERO AND FLOW IMBALANCE CORRECTION ACTIVE.')"
    character(len=*), parameter :: fmtifc = &
      &"(4x,'MASS WILL BE ADDED OR REMOVED TO COMPENSATE FOR FLOW IMBALANCE.')"

    write (this%iout, '(1x,a)') 'PROCESSING FMI OPTIONS'

    call mem_set_value(this%ipakcb, 'SAVE_FLOWS', this%input_mempath, &
                       found_ipakcb)
    call mem_set_value(this%iflowerr, 'IMBALANCECORRECT', this%input_mempath, &
                       found_flowerr)

    if (found_ipakcb) then
      this%ipakcb = -1
      write (this%iout, fmtisvflow)
    end if
    if (found_flowerr) write (this%iout, fmtifc)

    write (this%iout, '(1x,a)') 'END OF FMI OPTIONS'
  end subroutine gwtfmi_source_options

  !> @ brief Source input options for package
  !<
  subroutine gwtfmi_source_packagedata(this)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    use MemoryManagerExtModule, only: mem_set_value
    use CharacterStringModule, only: CharacterStringType
    use OpenSpecModule, only: ACCESS, FORM
    use ConstantsModule, only: LINELENGTH, DEM6, LENPACKAGENAME
    use InputOutputModule, only: getunit, openfile, urdaux
    ! -- dummy
    class(TspFmiType) :: this
    ! -- local
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: flowtypes
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: fileops
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: fnames
    type(BudgetObjectType), pointer :: budobjptr
    type(BudObjPtrArray), dimension(:), allocatable :: tmpbudobj
    character(len=LINELENGTH) :: flowtype, fileop, fname
    integer(I4B) :: iapt, inunit, n, i
    logical(LGP) :: exist

    iapt = 0

    call mem_setptr(flowtypes, 'FLOWTYPE', this%input_mempath)
    call mem_setptr(fileops, 'FILEIN', this%input_mempath)
    call mem_setptr(fnames, 'FNAME', this%input_mempath)

    do n = 1, size(flowtypes)
      flowtype = flowtypes(n)
      fileop = fileops(n)
      fname = fnames(n)

      inquire (file=trim(fname), exist=exist)
      if (.not. exist) then
        call store_error('Could not find file '//trim(fname))
        cycle
      end if

      if (fileop /= 'FILEIN') then
        call store_error('Unexpected packagedata input keyword read: "' &
                         //trim(fileop)//'".')
        cycle
      end if

      select case (flowtype)
      case ('GWFBUDGET')
        inunit = getunit()
        call openfile(inunit, this%iout, fname, 'DATA(BINARY)', FORM, &
                      ACCESS, 'UNKNOWN')
        this%iubud = inunit
        call this%initialize_bfr()
      case ('GWFHEAD')
        inunit = getunit()
        call openfile(inunit, this%iout, fname, 'DATA(BINARY)', FORM, &
                      ACCESS, 'UNKNOWN')
        this%iuhds = inunit
        call this%initialize_hfr()
      case ('GWFMOVER')
        inunit = getunit()
        call openfile(inunit, this%iout, fname, 'DATA(BINARY)', FORM, &
                      ACCESS, 'UNKNOWN')
        this%iumvr = inunit
        call budgetobject_cr_bfr(this%mvrbudobj, 'MVT', this%iumvr, &
                                 this%iout)
        call this%mvrbudobj%fill_from_bfr(this%dis, this%iout)
      case default
        !
        ! --expand the size of aptbudobj, which stores a pointer to the budobj
        allocate (tmpbudobj(iapt))
        do i = 1, size(this%aptbudobj)
          tmpbudobj(i)%ptr => this%aptbudobj(i)%ptr
        end do
        deallocate (this%aptbudobj)
        allocate (this%aptbudobj(iapt + 1))
        do i = 1, size(tmpbudobj)
          this%aptbudobj(i)%ptr => tmpbudobj(i)%ptr
        end do
        deallocate (tmpbudobj)
        !
        ! -- Open the budget file and start filling it
        iapt = iapt + 1
        inunit = getunit()
        call openfile(inunit, this%iout, fname, 'DATA(BINARY)', FORM, &
                      ACCESS, 'UNKNOWN')
        call budgetobject_cr_bfr(budobjptr, flowtype, inunit, &
                                 this%iout, colconv2=['GWF             '])
        call budobjptr%fill_from_bfr(this%dis, this%iout)
        this%aptbudobj(iapt)%ptr => budobjptr
      end select
    end do

    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine gwtfmi_source_packagedata

  !> @brief Set the pointer to a budget object
  !!
  !! An advanced transport can pass in a name and a
  !! pointer budget object, and this routine will look through the budget
  !! objects managed by FMI and point to the one with the same name, such as
  !! LAK-1, SFR-1, etc.
  !<
  subroutine set_aptbudobj_pointer(this, name, budobjptr)
    ! -- modules
    class(TspFmiType) :: this
    ! -- dumm
    character(len=*), intent(in) :: name
    type(BudgetObjectType), pointer :: budobjptr
    ! -- local
    integer(I4B) :: i
    !
    ! -- find and set the pointer
    do i = 1, size(this%aptbudobj)
      if (this%aptbudobj(i)%ptr%name == name) then
        budobjptr => this%aptbudobj(i)%ptr
        exit
      end if
    end do
  end subroutine set_aptbudobj_pointer

  !> @brief Initialize the groundwater flow terms based on the budget file
  !! reader
  !!
  !! Initialize terms and figure out how many different terms and packages
  !! are contained within the file
  !<
  subroutine initialize_gwfterms_from_bfr(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(TspFmiType) :: this
    ! -- local
    integer(I4B) :: nflowpack
    integer(I4B) :: i, ip
    integer(I4B) :: naux
    logical :: found_flowja
    logical :: found_dataspdis
    logical :: found_datasat
    logical :: found_stoss
    logical :: found_stosy
    integer(I4B), dimension(:), allocatable :: imap
    !
    ! -- Calculate the number of gwf flow packages
    allocate (imap(this%bfr%nbudterms))
    imap(:) = 0
    nflowpack = 0
    found_flowja = .false.
    found_dataspdis = .false.
    found_datasat = .false.
    found_stoss = .false.
    found_stosy = .false.
    do i = 1, this%bfr%nbudterms
      select case (trim(adjustl(this%bfr%budtxtarray(i))))
      case ('FLOW-JA-FACE')
        found_flowja = .true.
      case ('DATA-SPDIS')
        found_dataspdis = .true.
      case ('DATA-SAT')
        found_datasat = .true.
      case ('STO-SS')
        found_stoss = .true.
        this%igwfstrgss = 1
      case ('STO-SY')
        found_stosy = .true.
        this%igwfstrgsy = 1
      case default
        nflowpack = nflowpack + 1
        imap(i) = 1
      end select
    end do
    !
    ! -- allocate gwfpackage arrays (gwfpackages, iatp, datp, ...)
    call this%allocate_gwfpackages(nflowpack)
    !
    ! -- Copy the package name and aux names from budget file reader
    !    to the gwfpackages derived-type variable
    ip = 1
    do i = 1, this%bfr%nbudterms
      if (imap(i) == 0) cycle
      call this%gwfpackages(ip)%set_name(this%bfr%dstpackagenamearray(i), &
                                         this%bfr%budtxtarray(i))
      naux = this%bfr%nauxarray(i)
      call this%gwfpackages(ip)%set_auxname(naux, &
                                            this%bfr%auxtxtarray(1:naux, i))
      ip = ip + 1
    end do
    !
    ! -- Copy just the package names for the boundary packages into
    !    the flowpacknamearray
    ip = 1
    do i = 1, size(imap)
      if (imap(i) == 1) then
        this%flowpacknamearray(ip) = this%bfr%dstpackagenamearray(i)
        ip = ip + 1
      end if
    end do
    !
    ! -- Error if specific discharge, saturation or flowja not found
    if (.not. found_dataspdis) then
      write (errmsg, '(a)') 'Specific discharge not found in &
                            &budget file. SAVE_SPECIFIC_DISCHARGE and &
                            &SAVE_FLOWS must be activated in the NPF package.'
      call store_error(errmsg)
    end if
    if (.not. found_datasat) then
      write (errmsg, '(a)') 'Saturation not found in &
                            &budget file. SAVE_SATURATION and &
                            &SAVE_FLOWS must be activated in the NPF package.'
      call store_error(errmsg)
    end if
    if (.not. found_flowja) then
      write (errmsg, '(a)') 'FLOWJA not found in &
                            &budget file. SAVE_FLOWS must &
                            &be activated in the NPF package.'
      call store_error(errmsg)
    end if
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine initialize_gwfterms_from_bfr

  !> @brief Initialize groundwater flow terms from the groundwater budget
  !!
  !! Flows are coming from a gwf-gwt exchange object
  !<
  subroutine initialize_gwfterms_from_gwfbndlist(this)
    ! -- modules
    use BndModule, only: BndType, GetBndFromList
    ! -- dummy
    class(TspFmiType) :: this
    ! -- local
    integer(I4B) :: ngwfpack
    integer(I4B) :: ngwfterms
    integer(I4B) :: ip
    integer(I4B) :: imover
    integer(I4B) :: ntomvr
    integer(I4B) :: iterm
    character(len=LENPACKAGENAME) :: budtxt
    class(BndType), pointer :: packobj => null()
    !
    ! -- determine size of gwf terms
    ngwfpack = this%gwfbndlist%Count()
    !
    ! -- Count number of to-mvr terms, but do not include advanced packages
    !    as those mover terms are not losses from the cell, but rather flows
    !    within the advanced package
    ntomvr = 0
    do ip = 1, ngwfpack
      packobj => GetBndFromList(this%gwfbndlist, ip)
      imover = packobj%imover
      if (packobj%isadvpak /= 0) imover = 0
      if (imover /= 0) then
        ntomvr = ntomvr + 1
      end if
    end do
    !
    ! -- Allocate arrays in fmi of size ngwfterms, which is the number of
    !    packages plus the number of packages with mover terms.
    ngwfterms = ngwfpack + ntomvr
    call this%allocate_gwfpackages(ngwfterms)
    !
    ! -- Assign values in the fmi package
    iterm = 1
    do ip = 1, ngwfpack
      !
      ! -- set and store names
      packobj => GetBndFromList(this%gwfbndlist, ip)
      budtxt = adjustl(packobj%text)
      call this%gwfpackages(iterm)%set_name(packobj%packName, budtxt)
      this%flowpacknamearray(iterm) = packobj%packName
      call this%gwfpackages(iterm)%set_auxname(packobj%naux, &
                                               packobj%auxname)
      iterm = iterm + 1
      !
      ! -- if this package has a mover associated with it, then add another
      !    term that corresponds to the mover flows
      imover = packobj%imover
      if (packobj%isadvpak /= 0) imover = 0
      if (imover /= 0) then
        budtxt = trim(adjustl(packobj%text))//'-TO-MVR'
        call this%gwfpackages(iterm)%set_name(packobj%packName, budtxt)
        this%flowpacknamearray(iterm) = packobj%packName
        call this%gwfpackages(iterm)%set_auxname(packobj%naux, &
                                                 packobj%auxname)
        this%igwfmvrterm(iterm) = 1
        iterm = iterm + 1
      end if
    end do
  end subroutine initialize_gwfterms_from_gwfbndlist

  !> @brief Initialize an array for storing PackageBudget objects.
  !!
  !! This routine allocates gwfpackages (an array of PackageBudget
  !! objects) to the proper size and initializes member variables.
  !<
  subroutine gwtfmi_allocate_gwfpackages(this, ngwfterms)
    ! -- modules
    use ConstantsModule, only: LENMEMPATH
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(TspFmiType) :: this
    integer(I4B), intent(in) :: ngwfterms
    ! -- local
    integer(I4B) :: n
    character(len=LENMEMPATH) :: memPath
    !
    ! -- direct allocate
    allocate (this%gwfpackages(ngwfterms))
    allocate (this%flowpacknamearray(ngwfterms))
    allocate (this%datp(ngwfterms))
    !
    ! -- mem_allocate
    call mem_allocate(this%iatp, ngwfterms, 'IATP', this%memoryPath)
    call mem_allocate(this%igwfmvrterm, ngwfterms, 'IGWFMVRTERM', this%memoryPath)
    !
    ! -- initialize
    this%nflowpack = ngwfterms
    do n = 1, this%nflowpack
      this%iatp(n) = 0
      this%igwfmvrterm(n) = 0
      this%flowpacknamearray(n) = ''
      !
      ! -- Create a mempath for each individual flow package data set
      !    of the form, MODELNAME/FMI-FTn
      write (memPath, '(a, i0)') trim(this%memoryPath)//'-FT', n
      call this%gwfpackages(n)%initialize(memPath)
    end do
  end subroutine gwtfmi_allocate_gwfpackages

  !> @brief Deallocate memory
  !!
  !! Deallocate memory that stores the gwfpackages array
  !<
  subroutine gwtfmi_deallocate_gwfpackages(this)
    ! -- modules
    ! -- dummy
    class(TspFmiType) :: this
    ! -- local
    integer(I4B) :: n
    !
    ! -- initialize
    do n = 1, this%nflowpack
      call this%gwfpackages(n)%da()
    end do
  end subroutine gwtfmi_deallocate_gwfpackages

end module TspFmiModule
