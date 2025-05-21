! Comprehensive budget object that stores all of the
! intercell flows, and the inflows and the outflows for
! an advanced package.
module BudgetObjectModule

  use KindModule, only: I4B, DP
  use ConstantsModule, only: LENBUDTXT, LINELENGTH, &
                             TABLEFT, TABCENTER, TABRIGHT, &
                             TABSTRING, TABUCSTRING, TABINTEGER, TABREAL, &
                             DZERO, DHALF, DHUNDRED
  use BudgetModule, only: BudgetType, budget_cr
  use BudgetTermModule, only: BudgetTermType
  use TableModule, only: TableType, table_cr
  use BaseDisModule, only: DisBaseType
  use BudgetFileReaderModule, only: BudgetFileReaderType

  implicit none

  public :: BudgetObjectType
  public :: budgetobject_cr
  public :: budgetobject_cr_bfr

  type :: BudgetObjectType
    !
    ! -- name, number of control volumes, and number of budget terms
    character(len=LENBUDTXT) :: name
    integer(I4B) :: ncv
    integer(I4B) :: nbudterm
    !
    ! -- state variables
    real(DP), dimension(:), pointer :: xnew => null()
    real(DP), dimension(:), pointer :: xold => null()
    !
    ! -- csr intercell flows
    integer(I4B) :: iflowja
    real(DP), dimension(:), pointer :: flowja => null()
    !
    ! -- storage
    integer(I4B) :: nsto
    real(DP), dimension(:, :), pointer :: qsto => null()
    !
    ! -- array of budget terms, with one separate entry for each term
    !    such as rainfall, et, leakage, etc.
    integer(I4B) :: iterm
    type(BudgetTermType), dimension(:), allocatable :: budterm
    !
    ! -- budget table object, for writing the typical MODFLOW budget
    type(BudgetType), pointer :: budtable => null()
    !
    ! -- flow table object, for writing the flow budget for
    !    each control volume
    logical, pointer :: add_cellids => null()
    integer(I4B), pointer :: icellid => null()
    integer(I4B), pointer :: nflowterms => null()
    integer(I4B), dimension(:), pointer :: istart => null()
    integer(I4B), dimension(:), pointer :: iflowterms => null()
    type(TableType), pointer :: flowtab => null()
    !
    ! -- budget file reader, for reading flows from a binary file
    type(BudgetFileReaderType), pointer :: bfr => null()

  contains

    procedure :: budgetobject_df
    procedure :: flowtable_df
    procedure :: accumulate_terms
    procedure :: write_budtable
    procedure :: write_flowtable
    procedure :: save_flows
    procedure :: read_flows
    procedure :: budgetobject_da
    procedure :: bfr_init
    procedure :: bfr_advance
    procedure :: fill_from_bfr

  end type BudgetObjectType

contains

  !> @brief Create a new budget object
  !<
  subroutine budgetobject_cr(this, name)
    ! -- dummy
    type(BudgetObjectType), pointer :: this
    character(len=*), intent(in) :: name
    !
    ! -- Create the object
    allocate (this)
    !
    ! -- Initialize variables
    this%name = name
    this%ncv = 0
    this%nbudterm = 0
    this%iflowja = 0
    this%nsto = 0
    this%iterm = 0
    !
    ! -- Initialize budget table
    call budget_cr(this%budtable, name)
  end subroutine budgetobject_cr

  !> @brief Define the new budget object
  !<
  subroutine budgetobject_df(this, ncv, nbudterm, iflowja, nsto, &
                             bddim_opt, labeltitle_opt, bdzone_opt, &
                             ibudcsv)
    ! -- dummy
    class(BudgetObjectType) :: this
    integer(I4B), intent(in) :: ncv
    integer(I4B), intent(in) :: nbudterm
    integer(I4B), intent(in) :: iflowja
    integer(I4B), intent(in) :: nsto
    character(len=*), optional :: bddim_opt
    character(len=*), optional :: labeltitle_opt
    character(len=*), optional :: bdzone_opt
    integer(I4B), intent(in), optional :: ibudcsv
    ! -- local
    character(len=20) :: bdtype
    character(len=5) :: bddim
    character(len=16) :: labeltitle
    character(len=20) :: bdzone
    !
    ! -- Set values
    this%ncv = ncv
    this%nbudterm = nbudterm
    this%iflowja = iflowja
    this%nsto = nsto
    !
    ! -- Allocate space for budterm
    allocate (this%budterm(nbudterm))
    !
    ! -- Set the budget type to name
    bdtype = this%name
    !
    ! -- Set the budget dimension
    if (present(bddim_opt)) then
      bddim = bddim_opt
    else
      bddim = 'L**3'
    end if
    !
    ! -- Set the budget zone
    if (present(bdzone_opt)) then
      bdzone = bdzone_opt
    else
      bdzone = 'ENTIRE MODEL'
    end if
    !
    ! -- Set the label title
    if (present(labeltitle_opt)) then
      labeltitle = labeltitle_opt
    else
      labeltitle = 'PACKAGE NAME'
    end if
    !
    ! -- Setup the budget table object
    call this%budtable%budget_df(nbudterm, bdtype, bddim, labeltitle, bdzone)
    !
    ! -- Trigger csv output
    if (present(ibudcsv)) then
      call this%budtable%set_ibudcsv(ibudcsv)
    end if
  end subroutine budgetobject_df

  !> @brief Define the new flow table object
  !<
  subroutine flowtable_df(this, iout, cellids)
    ! -- dummy
    class(BudgetObjectType) :: this
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in), optional :: cellids
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    character(len=LENBUDTXT) :: flowtype
    character(len=LENBUDTXT) :: tag
    character(len=LENBUDTXT) :: coupletype
    logical :: lfound
    logical :: add_cellids
    integer(I4B) :: maxcol
    integer(I4B) :: idx
    integer(I4B) :: ipos
    integer(I4B) :: i
    !
    ! -- Process optional variables
    if (present(cellids)) then
      add_cellids = .TRUE.
      coupletype = cellids
    else
      add_cellids = .FALSE.
    end if
    !
    ! -- Allocate scalars
    allocate (this%add_cellids)
    allocate (this%icellid)
    allocate (this%nflowterms)
    !
    ! -- Initialize scalars
    this%add_cellids = add_cellids
    this%nflowterms = 0
    this%icellid = 0
    !
    ! -- Determine the number of columns in the table
    maxcol = 3
    if (add_cellids) then
      maxcol = maxcol + 1
    end if
    do i = 1, this%nbudterm
      lfound = .FALSE.
      flowtype = this%budterm(i)%get_flowtype()
      if (trim(adjustl(flowtype)) == 'FLOW-JA-FACE') then
        lfound = .TRUE.
        maxcol = maxcol + 2
      else if (trim(adjustl(flowtype)) /= 'AUXILIARY') then
        lfound = .TRUE.
        maxcol = maxcol + 1
      end if
      if (lfound) then
        this%nflowterms = this%nflowterms + 1
        if (add_cellids) then
          if (trim(adjustl(flowtype)) == trim(adjustl(coupletype))) then
            this%icellid = i
          end if
        end if
      end if
    end do
    !
    ! -- Allocate arrays
    allocate (this%istart(this%nflowterms))
    allocate (this%iflowterms(this%nflowterms))
    !
    ! -- Set up flow tableobj
    title = trim(this%name)//' PACKAGE - SUMMARY OF FLOWS FOR '// &
            'EACH CONTROL VOLUME'
    call table_cr(this%flowtab, this%name, title)
    call this%flowtab%table_df(this%ncv, maxcol, iout, transient=.TRUE.)
    !
    ! -- Go through and set up flow table budget terms
    text = 'NUMBER'
    call this%flowtab%initialize_column(text, 10, alignment=TABCENTER)
    if (add_cellids) then
      text = 'CELLID'
      call this%flowtab%initialize_column(text, 20, alignment=TABLEFT)
    end if
    idx = 1
    do i = 1, this%nbudterm
      lfound = .FALSE.
      flowtype = this%budterm(i)%get_flowtype()
      tag = trim(adjustl(flowtype))
      ipos = index(tag, '-')
      if (ipos > 0) then
        tag(ipos:ipos) = ' '
      end if
      if (trim(adjustl(flowtype)) == 'FLOW-JA-FACE') then
        lfound = .TRUE.
        text = 'INFLOW'
        call this%flowtab%initialize_column(text, 12, alignment=TABCENTER)
        text = 'OUTFLOW'
        call this%flowtab%initialize_column(text, 12, alignment=TABCENTER)
      else if (trim(adjustl(flowtype)) /= 'AUXILIARY') then
        lfound = .TRUE.
        call this%flowtab%initialize_column(tag, 12, alignment=TABCENTER)
      end if
      if (lfound) then
        this%iflowterms(idx) = i
        idx = idx + 1
      end if
    end do
    text = 'IN - OUT'
    call this%flowtab%initialize_column(text, 12, alignment=TABCENTER)
    text = 'PERCENT DIFFERENCE'
    call this%flowtab%initialize_column(text, 12, alignment=TABCENTER)
  end subroutine flowtable_df

  !> @brief Add up accumulators and submit to budget table
  !<
  subroutine accumulate_terms(this)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(BudgetObjectType) :: this
    ! -- local
    character(len=LENBUDTXT) :: flowtype
    integer(I4B) :: i
    real(DP) :: ratin, ratout
    !
    ! -- Reset the budget table
    call this%budtable%reset()
    !
    ! -- Calculate the budget table terms
    do i = 1, this%nbudterm
      !
      ! -- Accumulate positive and negative flows for each budget term
      flowtype = this%budterm(i)%flowtype
      select case (trim(adjustl(flowtype)))
      case ('FLOW-JA-FACE')
        ! -- Skip
      case default
        !
        ! -- Calculate sum of positive and negative flows
        call this%budterm(i)%accumulate_flow(ratin, ratout)
        !
        ! -- Pass accumulators into the budget table
        call this%budtable%addentry(ratin, ratout, delt, flowtype)
      end select
    end do
  end subroutine accumulate_terms

  !> @brief Write the flow table for each advanced package control volume
  !<
  subroutine write_flowtable(this, dis, kstp, kper, cellidstr)
    ! -- dummy
    class(BudgetObjectType) :: this
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: kstp
    integer(I4B), intent(in) :: kper
    character(len=20), dimension(:), optional :: cellidstr
    ! -- local
    character(len=LENBUDTXT) :: flowtype
    character(len=20) :: cellid
    integer(I4B) :: nlist
    integer(I4B) :: id1
    integer(I4B) :: id2
    integer(I4B) :: icv
    integer(I4B) :: idx
    integer(I4B) :: i
    integer(I4B) :: j
    real(DP) :: v
    real(DP) :: qin
    real(DP) :: qout
    real(DP) :: q
    real(DP) :: qinflow
    real(DP) :: qoutflow
    real(DP) :: qerr
    real(DP) :: qavg
    real(DP) :: qpd
    !
    ! -- Reset starting position
    do j = 1, this%nflowterms
      this%istart(j) = 1
    end do
    !
    ! -- Set table kstp and kper
    call this%flowtab%set_kstpkper(kstp, kper)
    !
    ! -- Write the table
    do icv = 1, this%ncv
      call this%flowtab%add_term(icv)
      !
      ! -- Initialize flow terms for the control volume
      qin = DZERO
      qout = DZERO
      !
      ! -- Add cellid if required
      if (this%add_cellids) then
        if (present(cellidstr)) then
          !
          ! -- If there are not maxbound entries for this%budterm(idx),
          !    which can happen for sfr, for example, if 'none' connections
          !    are specified, then cellidstr should be passed in if the flow
          !    table needs a cellid label.
          cellid = cellidstr(icv)
        else
          !
          ! -- Determine the cellid for this entry.  The cellid, such as
          !    (1, 10, 10), is assumed to be in the id2 column of this budterm.
          j = this%icellid
          idx = this%iflowterms(j)
          i = this%istart(j)
          id2 = this%budterm(idx)%get_id2(i)
          if (id2 > 0) then
            call dis%noder_to_string(id2, cellid)
          else
            cellid = 'NONE'
          end if
        end if
        call this%flowtab%add_term(cellid)
      end if
      !
      ! -- Iterate over the flow terms
      do j = 1, this%nflowterms
        !
        ! -- Initialize flow terms for the row
        q = DZERO
        qinflow = DZERO
        qoutflow = DZERO
        !
        ! -- Determine the index, flowtype and length of
        !    the flowterm
        idx = this%iflowterms(j)
        flowtype = this%budterm(idx)%get_flowtype()
        nlist = this%budterm(idx)%get_nlist()
        !
        ! -- Iterate over the entries in the flowtype.  If id1 is not ordered
        !    then need to look through the entire list each time
        colterm: do i = this%istart(j), nlist
          id1 = this%budterm(idx)%get_id1(i)
          if (this%budterm(idx)%ordered_id1) then
            if (id1 > icv) then
              this%istart(j) = i
              exit colterm
            end if
          else
            if (id1 /= icv) then
              cycle colterm
            end if
          end if
          if (id1 /= icv) then
            v = DZERO
          else
            v = this%budterm(idx)%get_flow(i)
          end if
          if (trim(adjustl(flowtype)) == 'FLOW-JA-FACE') then
            if (v < DZERO) then
              qoutflow = qoutflow + v
            else
              qinflow = qinflow + v
            end if
          end if
          !
          ! -- Accumulators
          q = q + v
          if (v < DZERO) then
            qout = qout + v
          else
            qin = qin + v
          end if
        end do colterm
        !
        ! -- Add entry to table
        if (trim(adjustl(flowtype)) == 'FLOW-JA-FACE') then
          call this%flowtab%add_term(qinflow)
          call this%flowtab%add_term(qoutflow)
        else
          call this%flowtab%add_term(q)
        end if
      end do
      !
      ! -- Calculate in-out and percent difference
      qerr = qin + qout
      qavg = DHALF * (qin - qout)
      qpd = DZERO
      if (qavg > DZERO) then
        qpd = DHUNDRED * qerr / qavg
      end if
      call this%flowtab%add_term(qerr)
      call this%flowtab%add_term(qpd)
    end do
  end subroutine write_flowtable

  !> @brief Write the budget table
  !<
  subroutine write_budtable(this, kstp, kper, iout, ibudfl, totim, delt)
    ! -- dummy
    class(BudgetObjectType) :: this
    integer(I4B), intent(in) :: kstp
    integer(I4B), intent(in) :: kper
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(in) :: ibudfl
    real(DP), intent(in) :: totim
    real(DP), intent(in) :: delt
    !
    ! -- Write the table
    call this%budtable%finalize_step(delt)
    if (ibudfl /= 0) then
      call this%budtable%budget_ot(kstp, kper, iout)
    end if
    call this%budtable%writecsv(totim)
  end subroutine write_budtable

  !> @brief Write the budget table
  !<
  subroutine save_flows(this, dis, ibinun, kstp, kper, delt, &
                        pertim, totim, iout)
    ! -- dummy
    class(BudgetObjectType) :: this
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: ibinun
    integer(I4B), intent(in) :: kstp
    integer(I4B), intent(in) :: kper
    real(DP), intent(in) :: delt
    real(DP), intent(in) :: pertim
    real(DP), intent(in) :: totim
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B) :: i
    !
    ! -- Save flows for each budget term
    do i = 1, this%nbudterm
      call this%budterm(i)%save_flows(dis, ibinun, kstp, kper, delt, &
                                      pertim, totim, iout)
    end do
  end subroutine save_flows

  !> @brief Read from a binary file into this BudgetObjectType
  !<
  subroutine read_flows(this, dis, ibinun)
    ! -- dummy
    class(BudgetObjectType) :: this
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: ibinun
    ! -- local
    integer(I4B) :: kstp
    integer(I4B) :: kper
    real(DP) :: delt
    real(DP) :: pertim
    real(DP) :: totim
    integer(I4B) :: i
    !
    ! -- Read flows for each budget term
    do i = 1, this%nbudterm
      call this%budterm(i)%read_flows(dis, ibinun, kstp, kper, delt, &
                                      pertim, totim)
    end do
  end subroutine read_flows

  !> @brief Deallocate
  !<
  subroutine budgetobject_da(this)
    ! -- dummy
    class(BudgetObjectType) :: this
    ! -- local
    integer(I4B) :: i
    !
    ! -- Save flows for each budget term
    do i = 1, this%nbudterm
      call this%budterm(i)%deallocate_arrays()
    end do
    !
    ! -- Destroy the flow table
    if (associated(this%flowtab)) then
      deallocate (this%add_cellids)
      deallocate (this%icellid)
      deallocate (this%nflowterms)
      deallocate (this%istart)
      deallocate (this%iflowterms)
      call this%flowtab%table_da()
      deallocate (this%flowtab)
      nullify (this%flowtab)
    end if
    !
    ! -- Destroy the budget object table
    if (associated(this%budtable)) then
      call this%budtable%budget_da()
      deallocate (this%budtable)
      nullify (this%budtable)
    end if
  end subroutine budgetobject_da

  !> @brief Create a new budget object from a binary flow file
  !<
  subroutine budgetobject_cr_bfr(this, name, ibinun, iout, colconv1, colconv2)
    ! -- dummy
    type(BudgetObjectType), pointer :: this
    character(len=*), intent(in) :: name
    integer(I4B), intent(in) :: ibinun
    integer(I4B), intent(in) :: iout
    character(len=16), dimension(:), optional :: colconv1
    character(len=16), dimension(:), optional :: colconv2
    ! -- local
    integer(I4B) :: ncv, nbudterm
    integer(I4B) :: iflowja, nsto
    integer(I4B) :: i, j
    !
    ! -- Create the object
    call budgetobject_cr(this, name)
    !
    ! -- Initialize the budget file reader
    call this%bfr_init(ibinun, ncv, nbudterm, iout)
    !
    ! -- Define this budget object using number of control volumes and number
    !    of budget terms read from ibinun
    iflowja = 0
    nsto = 0
    call this%budgetobject_df(ncv, nbudterm, iflowja, nsto)
    !
    ! -- Set the conversion flags, which cause id1 or id2 to be converted from
    !    user node numbers to reduced node numbers
    if (present(colconv1)) then
      do i = 1, nbudterm
        do j = 1, size(colconv1)
          if (colconv1(j) == adjustl(this%bfr%budtxtarray(i))) then
            this%budterm(i)%olconv1 = .true.
            exit
          end if
        end do
      end do
    end if
    if (present(colconv2)) then
      do i = 1, nbudterm
        do j = 1, size(colconv2)
          if (colconv2(j) == adjustl(this%bfr%budtxtarray(i))) then
            this%budterm(i)%olconv2 = .true.
            exit
          end if
        end do
      end do
    end if
  end subroutine budgetobject_cr_bfr

  !> @brief Initialize the budget file reader
  !<
  subroutine bfr_init(this, ibinun, ncv, nbudterm, iout)
    ! -- dummy
    class(BudgetObjectType) :: this
    integer(I4B), intent(in) :: ibinun
    integer(I4B), intent(inout) :: ncv
    integer(I4B), intent(inout) :: nbudterm
    integer(I4B), intent(in) :: iout
    !
    ! -- Initialize budget file reader
    allocate (this%bfr)
    call this%bfr%initialize(ibinun, iout, ncv)
    nbudterm = this%bfr%nbudterms
  end subroutine bfr_init

  !> @brief Advance the binary file readers for setting the budget terms of
  !! the next time step
  !<
  subroutine bfr_advance(this, dis, iout)
    ! -- modules
    use TdisModule, only: kstp, kper
    ! -- dummy
    class(BudgetObjectType) :: this
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: iout
    ! -- local
    logical :: readnext
    ! -- formats
    character(len=*), parameter :: fmtkstpkper = &
      &"(1x,/1x, a, ' READING BUDGET TERMS FOR KSTP ', i0, ' KPER ', i0)"
    character(len=*), parameter :: fmtbudkstpkper = &
      "(1x,/1x, a, ' SETTING BUDGET TERMS FOR KSTP ', i0, ' AND KPER ', &
      &i0, ' TO BUDGET FILE TERMS FROM KSTP ', i0, ' AND KPER ', i0)"
    !
    ! -- Do not read the budget if the budget is at end of file or if the next
    !    record in the budget file is the first timestep of the next stress
    !    period.  Also do not read if it is the very first time step because
    !    the first chunk of data is read as part of the initialization
    readnext = .true.
    if (kstp * kper == 1) then
      readnext = .false.
    else if (kstp * kper > 1) then
      if (this%bfr%endoffile) then
        readnext = .false.
      else
        if (this%bfr%kpernext == kper + 1 .and. this%bfr%kstpnext == 1) &
          readnext = .false.
      end if
    end if
    !
    ! -- Read the next record
    if (readnext) then
      !
      ! -- Write the current time step and stress period
      if (iout > 0) &
        write (iout, fmtkstpkper) this%name, kstp, kper
      !
      ! -- Read flows from the binary file and copy them into this%budterm(:)
      call this%fill_from_bfr(dis, iout)
    else
      if (iout > 0) &
        write (iout, fmtbudkstpkper) trim(this%name), kstp, kper, &
        this%bfr%kstp, this%bfr%kper
    end if
  end subroutine bfr_advance

  !> @brief Copy the information from the binary file into budterms
  !<
  subroutine fill_from_bfr(this, dis, iout)
    ! -- dummy
    class(BudgetObjectType) :: this
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B) :: i
    logical :: success
    !
    ! -- Read flows from the binary file and copy them into this%budterm(:)
    do i = 1, this%nbudterm
      call this%bfr%read_record(success, iout)
      call this%budterm(i)%fill_from_bfr(this%bfr, dis)
    end do
  end subroutine fill_from_bfr

end module BudgetObjectModule
