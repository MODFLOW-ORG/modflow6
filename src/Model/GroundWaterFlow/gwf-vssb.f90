!> @brief Variably Saturated Surface Boundary (VSSB) package
!!
!! The land-surface boundary for variably saturated (Richards) flow models.
!! Precipitation is applied to the top face of a fixed set of surface cells,
!! infiltration is limited by the soil, water that cannot infiltrate is stored
!! in a pond of limited depth, groundwater can exfiltrate into the pond, and
!! water that cannot be stored leaves the package as runoff, which can be
!! routed with the water mover.  Infiltration, exfiltration, pond depth, and
!! runoff are results of the solution; none of them is prescribed.
!!
!! The exchange between the pond and the top cell is calculated by
!! VsfSurfaceExchangeModule.
!<
module VssbModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE, DHALF, DNODATA, DHNOFLO, &
                             LINELENGTH, LENFTYPE, LENPACKAGENAME, &
                             LENBOUNDNAME, LENBUDTXT, LENPAKLOC, LENMEMPATH, &
                             NAMEDBOUNDFLAG, IWETLAKE, TABLEFT, TABCENTER
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, mem_setptr, &
                                 get_isize
  use MemoryManagerExtModule, only: mem_set_value, memorystore_release
  use MemoryHelperModule, only: create_mem_path
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use CharacterStringModule, only: CharacterStringType
  use BudgetObjectModule, only: BudgetObjectType, budgetobject_cr
  use ObserveModule, only: ObserveType
  use TableModule, only: TableType, table_cr
  use InputOutputModule, only: urword, getunit, openfile
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, store_error_filename, &
                       store_error_unit
  use BaseDisModule, only: DisBaseType
  use VsfSurfaceExchangeModule, only: SurfaceExchangeResultType, &
                                      surface_exchange
  use MatrixBaseModule

  implicit none
  private

  public :: vssb_create
  public :: VssbType

  character(len=LENFTYPE) :: ftype = 'VSSB'
  character(len=LENPACKAGENAME) :: text = '            VSSB'

  !> Width over which the minimum of the exchange terms is smoothed, as a
  !! fraction of the thickness of the cell, expressed as a head.
  real(DP), parameter :: SMOOTH_HEAD_FRACTION = 1.0e-3_DP

  type, extends(BndExtType) :: VssbType
    ! -- scalars
    integer(I4B), pointer :: nodes => null() !< number of surface cells
    integer(I4B), pointer :: ioutpond => null() !< unit for binary pond depth output (0 when inactive)
    ! -- surface cell data
    integer(I4B), dimension(:), pointer, contiguous :: igwfnode => null() !< reduced node number of each surface cell
    real(DP), dimension(:), pointer, contiguous :: maxpond => null() !< maximum pond depth
    real(DP), dimension(:), pointer, contiguous :: area => null() !< horizontal area of the surface cell
    character(len=LENBOUNDNAME), dimension(:), pointer, contiguous :: &
      vssbname => null() !< boundname of each surface cell
    ! -- state and results
    real(DP), dimension(:), pointer, contiguous :: pond0 => null() !< pond depth at the start of the time step
    real(DP), dimension(:), pointer, contiguous :: pond => null() !< pond depth at the end of the time step
    real(DP), dimension(:), pointer, contiguous :: pondcf => null() !< pond depth used to formulate the current iteration
    real(DP), dimension(:), pointer, contiguous :: qp => null() !< precipitation rate (volume per time)
    real(DP), dimension(:), pointer, contiguous :: qinf => null() !< exchange with the aquifer, positive into the aquifer
    real(DP), dimension(:), pointer, contiguous :: qro => null() !< runoff rate, before the mover
    real(DP), dimension(:), pointer, contiguous :: qrotomvr => null() !< runoff rate taken by the mover
    real(DP), dimension(:), pointer, contiguous :: qsto => null() !< pond storage change rate
    real(DP), dimension(:), pointer, contiguous :: qevap_pond => null() !< pond evaporation rate (>= 0)
    real(DP), dimension(:), pointer, contiguous :: qevap_soil => null() !< bare-soil evaporation rate (>= 0)
    real(DP), dimension(:), pointer, contiguous :: deriv => null() !< Newton term
    real(DP), dimension(:), pointer, contiguous :: pondout => null() !< pond depth of each model cell, for output
    ! -- input context pointers (time-series aware)
    real(DP), dimension(:), pointer, contiguous :: precip => null() !< precipitation rate
    real(DP), dimension(:), pointer, contiguous :: evap => null() !< potential evaporation rate
    real(DP), dimension(:, :), pointer, contiguous :: uauxvar => null() !< auxiliary variables of each surface cell
    ! -- pointers to model variables (not owned)
    integer(I4B), pointer :: gwfiss => null() !< steady-state flag of the model
    real(DP), dimension(:), pointer, contiguous :: gwfk11 => null() !< NPF horizontal hydraulic conductivity
    real(DP), dimension(:), pointer, contiguous :: gwfk33 => null() !< NPF vertical hydraulic conductivity
    integer(I4B), pointer :: gwfik33 => null() !< NPF flag, nonzero if K33 was specified
    real(DP), dimension(:), pointer, contiguous :: gwfkrel => null() !< NPF relative permeability
    real(DP), dimension(:), pointer, contiguous :: gwfdkrdh => null() !< NPF head derivative of the relative permeability
    ! -- budget
    type(BudgetObjectType), pointer :: budobj => null() !< package budget
  contains
    procedure :: vssb_allocate_vars
    procedure :: vssb_source_packagedata
    procedure :: vssb_print_input
    procedure :: vssb_check_cell_above
    procedure :: vssb_evaluate
    procedure :: vssb_setup_budobj
    procedure :: vssb_fill_budobj
    procedure :: allocate_scalars => vssb_allocate_scalars
    procedure :: source_options => vssb_source_options
    procedure :: source_dimensions => vssb_source_dimensions
    procedure :: bnd_ar => vssb_ar
    procedure :: bnd_rp => vssb_rp
    procedure :: bnd_ad => vssb_ad
    procedure :: bnd_cf => vssb_cf
    procedure :: bnd_fc => vssb_fc
    procedure :: bnd_fn => vssb_fn
    procedure :: bnd_cc => vssb_cc
    procedure :: bnd_cq => vssb_cq
    procedure :: bnd_ot_package_flows => vssb_ot_package_flows
    procedure :: bnd_ot_dv => vssb_ot_dv
    procedure :: bnd_ot_bdsummary => vssb_ot_bdsummary
    procedure :: bnd_da => vssb_da
    procedure :: define_listlabel
    ! -- observations
    procedure, public :: bnd_obs_supported => vssb_obs_supported
    procedure, public :: bnd_df_obs => vssb_df_obs
    procedure, public :: bnd_rp_obs => vssb_rp_obs
    procedure, public :: bnd_bd_obs => vssb_bd_obs
  end type VssbType

contains

  !> @brief Create a new VSSB package
  !<
  subroutine vssb_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                         mempath)
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    character(len=*), intent(in) :: mempath
    ! local
    type(VssbType), pointer :: vssbobj

    allocate (vssbobj)
    packobj => vssbobj

    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text

    call vssbobj%allocate_scalars()
    call packobj%pack_initialize()

    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 1
    packobj%iscloc = 0 ! not supported
    packobj%isadvpak = 1
    packobj%ictMemPath = create_mem_path(namemodel, 'NPF')

  end subroutine vssb_create

  !> @brief Allocate scalars
  !<
  subroutine vssb_allocate_scalars(this)
    class(VssbType) :: this

    call this%BndExtType%allocate_scalars()

    call mem_allocate(this%nodes, 'NODES', this%memoryPath)
    call mem_allocate(this%ioutpond, 'IOUTPOND', this%memoryPath)

    this%nodes = 0
    this%ioutpond = 0

  end subroutine vssb_allocate_scalars

  !> @brief Source the OPTIONS block from the input context
  !<
  subroutine vssb_source_options(this)
    use OpenSpecModule, only: access, form
    use GwfVssbInputModule, only: GwfVssbParamFoundType
    class(VssbType), intent(inout) :: this
    ! local
    type(GwfVssbParamFoundType) :: found
    character(len=LINELENGTH) :: fname

    ! AUXILIARY, BOUNDNAMES, PRINT_INPUT, PRINT_FLOWS, SAVE_FLOWS, OBS6,
    ! and AUXMULTNAME
    call this%BndExtType%source_options()

    write (this%iout, '(1x,a)') &
      'PROCESSING '//trim(adjustl(this%text))//' OPTIONS'

    call mem_set_value(this%imover, 'MOVER', this%input_mempath, found%mover)
    if (found%mover) then
      write (this%iout, '(4x,a)') 'WATER MOVER IS ACTIVE FOR RUNOFF'
    end if

    call mem_set_value(fname, 'PONDDEPTHFILE', this%input_mempath, &
                       found%ponddepthfile)
    if (found%ponddepthfile) then
      this%ioutpond = getunit()
      call openfile(this%ioutpond, this%iout, trim(fname), 'DATA(BINARY)', &
                    form, access, 'REPLACE')
      write (this%iout, '(4x,3a)') 'POND DEPTH WILL BE WRITTEN TO ', &
        trim(fname), ' BASED ON HEAD SAVE RECORD IN OUTPUT CONTROL'
    end if

    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%text))//' OPTIONS'

  end subroutine vssb_source_options

  !> @brief Source the DIMENSIONS and PACKAGEDATA blocks from the input context
  !<
  subroutine vssb_source_dimensions(this)
    use GwfVssbInputModule, only: GwfVssbParamFoundType
    class(VssbType), intent(inout) :: this
    ! local
    type(GwfVssbParamFoundType) :: found

    this%nodes = -1

    write (this%iout, '(/1x,a)') &
      'PROCESSING '//trim(adjustl(this%text))//' DIMENSIONS'
    call mem_set_value(this%nodes, 'MAXBOUND', this%input_mempath, &
                       found%maxbound)
    if (found%maxbound) then
      write (this%iout, '(4x,a,i0)') 'NSURFCELLS = ', this%nodes
    end if
    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%text))//' DIMENSIONS'

    this%maxbound = this%maxbound + this%nodes
    this%nbound = this%maxbound

    if (this%nodes <= 0) then
      write (errmsg, '(a)') &
        'NSURFCELLS (MAXBOUND) was not specified or was specified incorrectly.'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if

    call this%define_listlabel()
    call this%vssb_allocate_vars()
    call this%vssb_source_packagedata()
    if (this%iprpak /= 0) call this%vssb_print_input()
    call this%vssb_setup_budobj()

  end subroutine vssb_source_dimensions

  !> @brief Allocate the surface cell arrays and set the input context pointers
  !<
  subroutine vssb_allocate_vars(this)
    class(VssbType), intent(inout) :: this
    ! local
    integer(I4B) :: i, j

    call mem_allocate(this%igwfnode, this%nodes, 'IGWFNODE', this%memoryPath)
    call mem_allocate(this%maxpond, this%nodes, 'MAXPOND', this%memoryPath)
    call mem_allocate(this%area, this%nodes, 'AREA', this%memoryPath)
    call mem_allocate(this%pond0, this%nodes, 'POND0', this%memoryPath)
    call mem_allocate(this%pond, this%nodes, 'POND', this%memoryPath)
    call mem_allocate(this%pondcf, this%nodes, 'PONDCF', this%memoryPath)
    call mem_allocate(this%qp, this%nodes, 'QP', this%memoryPath)
    call mem_allocate(this%qinf, this%nodes, 'QINF', this%memoryPath)
    call mem_allocate(this%qro, this%nodes, 'QRO', this%memoryPath)
    call mem_allocate(this%qrotomvr, this%nodes, 'QROTOMVR', this%memoryPath)
    call mem_allocate(this%qsto, this%nodes, 'QSTO', this%memoryPath)
    call mem_allocate(this%qevap_pond, this%nodes, 'QEVAP_POND', this%memoryPath)
    call mem_allocate(this%qevap_soil, this%nodes, 'QEVAP_SOIL', this%memoryPath)
    call mem_allocate(this%deriv, this%nodes, 'DERIV', this%memoryPath)
    call mem_allocate(this%pondout, this%dis%nodes, 'PONDOUT', this%memoryPath)
    allocate (this%vssbname(this%nodes))

    call mem_setptr(this%precip, 'PRECIPITATION', this%input_mempath)
    call mem_setptr(this%evap, 'EVAPORATION', this%input_mempath)
    call mem_setptr(this%uauxvar, 'AUXVAR', this%input_mempath)

    do i = 1, this%nodes
      this%igwfnode(i) = 0
      this%maxpond(i) = DZERO
      this%area(i) = DZERO
      this%pond0(i) = DZERO
      this%pond(i) = DZERO
      this%pondcf(i) = DZERO
      this%qp(i) = DZERO
      this%qinf(i) = DZERO
      this%qro(i) = DZERO
      this%qrotomvr(i) = DZERO
      this%qsto(i) = DZERO
      this%qevap_pond(i) = DZERO
      this%qevap_soil(i) = DZERO
      this%deriv(i) = DZERO
      this%vssbname(i) = ''
      do j = 1, this%naux
        if (this%iauxmultcol > 0 .and. j == this%iauxmultcol) then
          this%uauxvar(j, i) = DONE
        else
          this%uauxvar(j, i) = DZERO
        end if
      end do
    end do
    this%pondout(:) = DZERO

  end subroutine vssb_allocate_vars

  !> @brief Source the PACKAGEDATA block from the input context
  !<
  subroutine vssb_source_packagedata(this)
    class(VssbType), intent(inout) :: this
    ! local
    integer(I4B), dimension(:), pointer, contiguous :: ifno => null()
    integer(I4B), dimension(:, :), pointer, contiguous :: cellid_idm => null()
    real(DP), dimension(:), pointer, contiguous :: max_pond_depth => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      bndname => null()
    character(len=LINELENGTH) :: cellid_str
    character(len=30) :: nodestr
    integer(I4B), dimension(:), allocatable :: nboundchk
    integer(I4B), dimension(:), allocatable :: nodecount
    integer(I4B) :: n, i, node

    write (this%iout, '(/1x,3a)') 'PROCESSING ', trim(adjustl(this%text)), &
      ' PACKAGEDATA'

    call mem_setptr(ifno, 'PACKAGEDATA_IFNO', this%input_mempath)
    call mem_setptr(cellid_idm, 'CELLID', this%input_mempath)
    call mem_setptr(max_pond_depth, 'MAX_POND_DEPTH', this%input_mempath)
    if (this%inamedbound /= 0) then
      call mem_setptr(bndname, 'BOUNDNAME', this%input_mempath)
    end if

    allocate (nboundchk(this%nodes))
    allocate (nodecount(this%dis%nodes))
    nboundchk(:) = 0
    nodecount(:) = 0

    do n = 1, size(ifno)
      i = this%validate_ifno(ifno(n), this%nodes, nboundchk, &
                             'surface cell', 'PACKAGEDATA')
      if (i == 0) cycle

      write (cellid_str, '(10(i0,1x))') cellid_idm(1:this%dis%ndim, n)
      node = this%dis%noder_from_cellid(cellid_str, this%inunit, this%iout)
      this%igwfnode(i) = node

      ! the cell set is fixed, so a duplicate cell is reported once
      if (node > 0) then
        nodecount(node) = nodecount(node) + 1
        if (nodecount(node) == 2) then
          call this%dis%noder_to_string(node, nodestr)
          write (errmsg, '(a,1x,a,1x,a)') &
            'Cell', trim(adjustl(nodestr)), &
            'is specified more than once in the VSSB PACKAGEDATA block.'
          call store_error(errmsg)
        end if
      end if

      if (max_pond_depth(n) < DZERO) then
        write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
          'MAX_POND_DEPTH for surface cell', i, &
          'must be greater than or equal to 0 (specified value is', &
          max_pond_depth(n), ').'
        call store_error(errmsg)
      end if
      this%maxpond(i) = max_pond_depth(n)

      if (this%inamedbound == 1) then
        this%vssbname(i) = bndname(n)
      end if
    end do

    write (this%iout, '(1x,3a)') &
      'END OF ', trim(adjustl(this%text)), ' PACKAGEDATA'

    call this%report_ifno_coverage(nboundchk, this%nodes, 'surface cell', &
                                   'PACKAGEDATA')
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if

    deallocate (nboundchk)
    deallocate (nodecount)

    call memorystore_release('PACKAGEDATA_IFNO', this%input_mempath)
    call memorystore_release('CELLID', this%input_mempath)
    call memorystore_release('MAX_POND_DEPTH', this%input_mempath)
    call memorystore_release('BOUNDNAME', this%input_mempath)

  end subroutine vssb_source_packagedata

  !> @brief Write the PACKAGEDATA to the listing file
  !<
  subroutine vssb_print_input(this)
    class(VssbType), intent(inout) :: this
    ! local
    character(len=20) :: cellid
    character(len=LINELENGTH) :: title, tag
    integer(I4B) :: ntabcols, i

    ntabcols = 3
    if (this%inamedbound == 1) ntabcols = ntabcols + 1
    title = trim(adjustl(this%text))//' PACKAGE ('// &
            trim(adjustl(this%packName))//') STATIC SURFACE CELL DATA'
    call table_cr(this%inputtab, this%packName, title)
    call this%inputtab%table_df(this%nodes, ntabcols, this%iout)
    tag = 'NUMBER'
    call this%inputtab%initialize_column(tag, 10)
    tag = 'CELLID'
    call this%inputtab%initialize_column(tag, 20, alignment=TABLEFT)
    tag = 'MAX_POND_DEPTH'
    call this%inputtab%initialize_column(tag, 16)
    if (this%inamedbound == 1) then
      tag = 'BOUNDNAME'
      call this%inputtab%initialize_column(tag, LENBOUNDNAME, alignment=TABLEFT)
    end if
    do i = 1, this%nodes
      if (this%igwfnode(i) > 0) then
        call this%dis%noder_to_string(this%igwfnode(i), cellid)
      else
        cellid = 'none'
      end if
      call this%inputtab%add_term(i)
      call this%inputtab%add_term(cellid)
      call this%inputtab%add_term(this%maxpond(i))
      if (this%inamedbound == 1) then
        call this%inputtab%add_term(this%vssbname(i))
      end if
    end do
    call this%inputtab%table_da()
    deallocate (this%inputtab)
    nullify (this%inputtab)

  end subroutine vssb_print_input

  !> @brief Allocate and read
  !<
  subroutine vssb_ar(this)
    class(VssbType), intent(inout) :: this
    ! local
    character(len=LENMEMPATH) :: npfmempath, vsfmempath
    integer(I4B) :: i, node, isize

    call this%obs%obs_ar()

    call this%BndExtType%allocate_arrays()

    call this%copy_boundname()
    if (this%inamedbound /= 0) then
      do i = 1, this%nodes
        this%boundname(i) = this%vssbname(i)
      end do
    end if

    ! VSSB applies to Richards cells and reads the relative permeability that
    ! the VSF package populates
    vsfmempath = create_mem_path(this%name_model, 'VSF')
    call get_isize('IVARSAT', vsfmempath, isize)
    if (isize < 0) then
      write (errmsg, '(a)') 'The VSSB package requires the VSF package to be '// &
        'active in the same GWF model.'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if

    npfmempath = create_mem_path(this%name_model, 'NPF')
    call mem_setptr(this%gwfk11, 'K11', npfmempath)
    call mem_setptr(this%gwfk33, 'K33', npfmempath)
    call mem_setptr(this%gwfik33, 'IK33', npfmempath)
    call mem_setptr(this%gwfkrel, 'KREL', npfmempath)
    call mem_setptr(this%gwfdkrdh, 'DKRDH', npfmempath)
    call mem_setptr(this%gwfiss, 'ISS', create_mem_path(this%name_model))

    do i = 1, this%nodes
      node = this%igwfnode(i)
      this%nodelist(i) = node
      this%area(i) = this%dis%get_area(node)
    end do

    call this%vssb_check_cell_above()

    if (this%imover /= 0) then
      allocate (this%pakmvrobj)
      call this%pakmvrobj%ar(this%maxbound, 0, this%memoryPath)
    end if

  end subroutine vssb_ar

  !> @brief Check that no active cell is directly above a surface cell
  !!
  !! A surface cell must be at the top of the active domain.  The check is
  !! skipped for grids without vertical connections.
  !<
  subroutine vssb_check_cell_above(this)
    class(VssbType), intent(inout) :: this
    ! local
    character(len=30) :: nodestr, nodestr2
    integer(I4B) :: i, node, ipos, m

    do i = 1, this%nodes
      node = this%igwfnode(i)
      if (node < 1) cycle
      do ipos = this%dis%con%ia(node) + 1, this%dis%con%ia(node + 1) - 1
        if (this%dis%con%ihc(this%dis%con%jas(ipos)) /= 0) cycle
        m = this%dis%con%ja(ipos)
        if (this%ibound(m) == 0) cycle
        if (this%dis%bot(m) >= this%dis%top(node)) then
          call this%dis%noder_to_string(node, nodestr)
          call this%dis%noder_to_string(m, nodestr2)
          write (errmsg, '(a,1x,i0,1x,a,1x,a,a,1x,a,1x,a)') &
            'VSSB surface cell', i, 'in cell', trim(adjustl(nodestr)), &
            ' has an active cell directly above it (cell', &
            trim(adjustl(nodestr2)), &
            '). Surface cells must be at the top of the active domain.'
          call store_error(errmsg)
        end if
      end do
    end do
    if (count_errors() > 0) call store_error_filename(this%input_fname)

  end subroutine vssb_check_cell_above

  !> @brief Read and prepare period data
  !<
  subroutine vssb_rp(this)
    use TdisModule, only: kper
    class(VssbType), intent(inout) :: this
    ! local
    integer(I4B), pointer :: nlist => null()
    integer(I4B), dimension(:), pointer, contiguous :: ifno => null()
    character(len=20) :: cellid
    character(len=LINELENGTH) :: title, tag
    integer(I4B) :: n, i, ntabcols
    ! formats
    character(len=*), parameter :: fmtlsp = &
      &"(1X,/1X,'REUSING ',A,'S FROM LAST STRESS PERIOD')"

    this%nbound = this%maxbound

    if (this%iper /= kper) then
      write (this%iout, fmtlsp) trim(this%filtyp)
      return
    end if

    call mem_setptr(nlist, 'NBOUND', this%input_mempath)
    call mem_setptr(ifno, 'IFNO', this%input_mempath)

    if (this%iprpak /= 0) then
      ntabcols = 4
      if (this%inamedbound == 1) ntabcols = ntabcols + 1
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(adjustl(this%packName))//') DATA FOR PERIOD'
      write (title, '(a,1x,i6)') trim(adjustl(title)), kper
      call table_cr(this%inputtab, this%packName, title)
      call this%inputtab%table_df(1, ntabcols, this%iout, finalize=.false.)
      tag = 'NUMBER'
      call this%inputtab%initialize_column(tag, 10)
      tag = 'CELLID'
      call this%inputtab%initialize_column(tag, 20, alignment=TABLEFT)
      tag = 'PRECIPITATION'
      call this%inputtab%initialize_column(tag, 15)
      tag = 'EVAPORATION'
      call this%inputtab%initialize_column(tag, 15)
      if (this%inamedbound == 1) then
        tag = 'BOUNDNAME'
        call this%inputtab%initialize_column(tag, LENBOUNDNAME, &
                                             alignment=TABLEFT)
      end if
    end if

    do n = 1, nlist
      i = ifno(n)
      if (i < 1 .or. i > this%nodes) then
        write (errmsg, '(a,1x,i0,1x,a,1x,i0,a)') &
          'VSSB PERIOD data IFNO', i, &
          'must be greater than 0 and less than or equal to', this%nodes, '.'
        call store_error(errmsg)
        cycle
      end if
      if (this%iprpak /= 0) then
        if (this%igwfnode(i) > 0) then
          call this%dis%noder_to_string(this%igwfnode(i), cellid)
        else
          cellid = 'none'
        end if
        call this%inputtab%add_term(i)
        call this%inputtab%add_term(cellid)
        call this%inputtab%add_term(this%precip(i))
        call this%inputtab%add_term(this%evap(i))
        if (this%inamedbound == 1) then
          call this%inputtab%add_term(this%vssbname(i))
        end if
      end if
    end do

    if (this%iprpak /= 0) then
      call this%inputtab%finalize_table()
      call this%inputtab%table_da()
      deallocate (this%inputtab)
      nullify (this%inputtab)
    end if

    if (count_errors() > 0) call store_error_filename(this%input_fname)

  end subroutine vssb_rp

  !> @brief Advance the package one time step
  !<
  subroutine vssb_ad(this)
    use SimVariablesModule, only: iFailedStepRetry
    class(VssbType) :: this
    ! local
    integer(I4B) :: i, iaux

    call this%BndExtType%bnd_ad()

    ! the pond state at the end of a step becomes the state at the start of the
    ! next step, unless the step is repeated
    if (iFailedStepRetry == 0) then
      do i = 1, this%nodes
        this%pond0(i) = this%pond(i)
      end do
    end if

    ! precipitation, including time-series values, must not be negative
    do i = 1, this%nodes
      if (this%precip(i) < DZERO) then
        write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
          'PRECIPITATION for surface cell', i, &
          'must be greater than or equal to 0 (specified value is', &
          this%precip(i), ').'
        call store_error(errmsg)
      end if
    end do

    ! potential evaporation, including time-series values, must not be negative
    do i = 1, this%nodes
      if (this%evap(i) < DZERO) then
        write (errmsg, '(a,1x,i0,1x,a,1x,g0,a)') &
          'EVAPORATION for surface cell', i, &
          'must be greater than or equal to 0 (specified value is', &
          this%evap(i), ').'
        call store_error(errmsg)
      end if
    end do
    if (count_errors() > 0) call store_error_filename(this%input_fname)

    ! auxiliary variables are written to the budget file
    if (this%naux > 0) then
      do i = 1, this%nodes
        do iaux = 1, this%naux
          if (this%noupdateauxvar(iaux) /= 0) cycle
          this%auxvar(iaux, i) = this%uauxvar(iaux, i)
        end do
      end do
    end if

    if (this%imover == 1) call this%pakmvrobj%ad()

  end subroutine vssb_ad

  !> @brief Evaluate the surface exchange of a surface cell for a given head
  !!
  !! The precipitation rate, the pond state, and the runoff of the cell are
  !! updated by the caller from the results.
  !<
  subroutine vssb_evaluate(this, i, h, qp, res)
    use TdisModule, only: delt
    class(VssbType) :: this
    integer(I4B), intent(in) :: i !< surface cell
    real(DP), intent(in) :: h !< head of the cell
    real(DP), intent(out) :: qp !< precipitation rate (volume per time)
    type(SurfaceExchangeResultType), intent(out) :: res !< results
    ! local
    integer(I4B) :: node
    logical :: steady
    real(DP) :: dz, kv, cbase, epot

    node = this%igwfnode(i)
    steady = (this%gwfiss /= 0)

    qp = this%precip(i) * this%area(i)
    if (this%iauxmultcol > 0) qp = qp * this%uauxvar(this%iauxmultcol, i)
    epot = this%evap(i) * this%area(i)

    if (this%ibound(node) <= 0 .or. this%ibound(node) == IWETLAKE) then
      ! the cell cannot exchange water, so all water leaves as runoff and no
      ! evaporation is simulated (the result defaults leave the evaporation
      ! terms at zero)
      res%qinf = DZERO
      res%pond = DZERO
      res%hcof = DZERO
      res%rhs = DZERO
      res%deriv = DZERO
      res%qro = qp
      if (.not. steady) res%qro = res%qro + this%pond0(i) * this%area(i) / delt
      return
    end if

    dz = this%dis%top(node) - this%dis%bot(node)
    if (this%gwfik33 == 0) then
      kv = this%gwfk11(node)
    else
      kv = this%gwfk33(node)
    end if
    cbase = kv * this%area(i) / (DHALF * dz)

    call surface_exchange(h, this%dis%top(node), this%maxpond(i), &
                          this%area(i), cbase, this%gwfkrel(node), &
                          this%gwfdkrdh(node), qp, epot, this%pond0(i), delt, &
                          SMOOTH_HEAD_FRACTION * dz, steady, res)

  end subroutine vssb_evaluate

  !> @brief Formulate the HCOF and RHS terms
  !<
  subroutine vssb_cf(this)
    class(VssbType) :: this
    ! local
    type(SurfaceExchangeResultType) :: res
    integer(I4B) :: i
    real(DP) :: qp

    do i = 1, this%nodes
      call this%vssb_evaluate(i, this%xnew(this%igwfnode(i)), qp, res)
      this%qp(i) = qp
      this%qinf(i) = res%qinf
      this%pond(i) = res%pond
      this%pondcf(i) = res%pond
      this%qro(i) = res%qro
      this%qevap_pond(i) = res%qevap_pond
      this%qevap_soil(i) = res%qevap_soil
      this%deriv(i) = res%deriv
      this%hcof(i) = res%hcof
      this%rhs(i) = res%rhs
    end do

  end subroutine vssb_cf

  !> @brief Copy rhs and hcof into the solution rhs and amat
  !<
  subroutine vssb_fc(this, rhs, ia, idxglo, matrix_sln)
    class(VssbType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! local
    integer(I4B) :: i, n, ipos

    if (this%imover == 1) call this%pakmvrobj%fc()

    do i = 1, this%nodes
      n = this%nodelist(i)
      rhs(n) = rhs(n) + this%rhs(i)
      ipos = ia(n)
      call matrix_sln%add_value_pos(idxglo(ipos), this%hcof(i))

      ! the runoff is the water that is available to the mover
      if (this%imover == 1) then
        call this%pakmvrobj%accumulate_qformvr(i, this%qro(i))
      end if
    end do

  end subroutine vssb_fc

  !> @brief Fill the Newton terms
  !!
  !! The terms account for the dependence of the relative permeability on
  !! head.  The pond depth is frozen when linearizing with respect to head.
  !<
  subroutine vssb_fn(this, rhs, ia, idxglo, matrix_sln)
    class(VssbType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! local
    integer(I4B) :: i, n, ipos

    do i = 1, this%nodes
      n = this%nodelist(i)
      ipos = ia(n)
      call matrix_sln%add_value_pos(idxglo(ipos), this%deriv(i))
      rhs(n) = rhs(n) + this%deriv(i) * this%xnew(n)
    end do

  end subroutine vssb_fn

  !> @brief Final convergence check of the pond depth
  !!
  !! The pond depth that follows from the new heads is compared with the pond
  !! depth that was used to formulate the iteration.
  !<
  subroutine vssb_cc(this, innertot, kiter, iend, icnvgmod, cpak, ipak, dpak)
    class(VssbType), intent(inout) :: this
    integer(I4B), intent(in) :: innertot
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: iend
    integer(I4B), intent(in) :: icnvgmod
    character(len=LENPAKLOC), intent(inout) :: cpak
    integer(I4B), intent(inout) :: ipak
    real(DP), intent(inout) :: dpak
    ! local
    type(SurfaceExchangeResultType) :: res
    character(len=LENPAKLOC) :: cloc
    integer(I4B) :: i
    real(DP) :: qp, dpond

    ! check only if the model is otherwise considered converged
    if (icnvgmod == 0) return

    do i = 1, this%nodes
      call this%vssb_evaluate(i, this%xnew(this%igwfnode(i)), qp, res)
      dpond = res%pond - this%pondcf(i)
      if (abs(dpond) > abs(dpak)) then
        ipak = i
        dpak = dpond
        write (cloc, "(a,'-',a)") trim(this%packName), 'pond'
        cpak = trim(cloc)
      end if
    end do

  end subroutine vssb_cc

  !> @brief Calculate the final flows from the converged heads
  !!
  !! The exchange, the pond state, and the runoff are calculated from the final
  !! heads so that the water balance of the package closes exactly.
  !<
  subroutine vssb_cq(this, x, flowja, iadv)
    use TdisModule, only: delt
    class(VssbType), intent(inout) :: this
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    integer(I4B), optional, intent(in) :: iadv
    ! local
    type(SurfaceExchangeResultType) :: res
    integer(I4B) :: i
    real(DP) :: qp

    do i = 1, this%nodes
      call this%vssb_evaluate(i, x(this%igwfnode(i)), qp, res)
      this%qp(i) = qp
      this%qinf(i) = res%qinf
      this%pond(i) = res%pond
      this%qro(i) = res%qro
      this%qevap_pond(i) = res%qevap_pond
      this%qevap_soil(i) = res%qevap_soil
      this%deriv(i) = res%deriv
      this%hcof(i) = res%hcof
      this%rhs(i) = res%rhs
    end do

    ! the exchange with the aquifer is flowja and simvals; advanced package
    ! option, as the mover flows are handled below
    call this%BndType%bnd_cq(x, flowja, iadv=1)

    do i = 1, this%nodes
      this%qrotomvr(i) = DZERO
      if (this%imover == 1) then
        this%qrotomvr(i) = min(this%pakmvrobj%get_qtomvr(i), this%qro(i))
      end if
      if (this%gwfiss /= 0) then
        this%qsto(i) = DZERO
      else
        this%qsto(i) = (this%pond(i) - this%pond0(i)) * this%area(i) / delt
      end if
    end do

    call this%vssb_fill_budobj()

  end subroutine vssb_cq

  !> @brief Output the package flow terms
  !<
  subroutine vssb_ot_package_flows(this, icbcfl, ibudfl)
    use TdisModule, only: kstp, kper
    class(VssbType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl

    if (ibudfl /= 0 .and. this%iprflow /= 0) then
      call this%budobj%write_flowtable(this%dis, kstp, kper)
    end if

  end subroutine vssb_ot_package_flows

  !> @brief Save the pond depth to a binary file
  !<
  subroutine vssb_ot_dv(this, idvsave, idvprint)
    class(VssbType) :: this
    integer(I4B), intent(in) :: idvsave
    integer(I4B), intent(in) :: idvprint
    ! local
    character(len=1) :: cdatafmp = ' ', editdesc = ' '
    integer(I4B) :: i, nvaluesp, nwidthp

    if (idvsave == 0 .or. this%ioutpond == 0) return

    this%pondout(:) = DZERO
    do i = 1, this%nodes
      this%pondout(this%igwfnode(i)) = this%pond(i)
    end do
    call this%dis%record_array(this%pondout, this%iout, 0, this%ioutpond, &
                               '      POND DEPTH', cdatafmp, nvaluesp, &
                               nwidthp, editdesc, DHNOFLO)

  end subroutine vssb_ot_dv

  !> @brief Write the package budget to the listing file
  !<
  subroutine vssb_ot_bdsummary(this, kstp, kper, iout, ibudfl)
    use TdisModule, only: totim, delt
    class(VssbType) :: this
    integer(I4B), intent(in) :: kstp
    integer(I4B), intent(in) :: kper
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(in) :: ibudfl

    call this%budobj%write_budtable(kstp, kper, iout, ibudfl, totim, delt)

  end subroutine vssb_ot_bdsummary

  !> @brief Set up the budget object that stores the package flows
  !!
  !! The terms must correspond in number and order to the terms in
  !! vssb_fill_budobj.
  !<
  subroutine vssb_setup_budobj(this)
    class(VssbType) :: this
    ! local
    integer(I4B) :: nbudterm, idx
    character(len=LENBUDTXT) :: text
    character(len=LENBUDTXT), dimension(1) :: auxtxt

    nbudterm = 5
    if (this%imover == 1) nbudterm = nbudterm + 1
    if (this%naux > 0) nbudterm = nbudterm + 1

    call budgetobject_cr(this%budobj, this%packName)
    call this%budobj%budgetobject_df(this%maxbound, nbudterm, 0, 0)
    idx = 0

    text = '   PRECIPITATION'
    idx = idx + 1
    call this%budobj%budterm(idx)%initialize(text, this%name_model, &
                                             this%packName, this%name_model, &
                                             this%packName, this%nodes, &
                                             .false., .false., 0)

    text = '             GWF'
    idx = idx + 1
    call this%budobj%budterm(idx)%initialize(text, this%name_model, &
                                             this%packName, this%name_model, &
                                             this%name_model, this%nodes, &
                                             .false., .true., 0)

    text = '          RUNOFF'
    idx = idx + 1
    call this%budobj%budterm(idx)%initialize(text, this%name_model, &
                                             this%packName, this%name_model, &
                                             this%packName, this%nodes, &
                                             .false., .false., 0)

    text = '         STORAGE'
    idx = idx + 1
    auxtxt(1) = '          VOLUME'
    call this%budobj%budterm(idx)%initialize(text, this%name_model, &
                                             this%packName, this%name_model, &
                                             this%packName, this%nodes, &
                                             .false., .false., 1, auxtxt)

    text = '     EVAPORATION'
    idx = idx + 1
    call this%budobj%budterm(idx)%initialize(text, this%name_model, &
                                             this%packName, this%name_model, &
                                             this%packName, this%nodes, &
                                             .false., .false., 0)

    if (this%imover == 1) then
      text = '   RUNOFF-TO-MVR'
      idx = idx + 1
      call this%budobj%budterm(idx)%initialize(text, this%name_model, &
                                               this%packName, this%name_model, &
                                               this%packName, this%nodes, &
                                               .false., .false., 0)
    end if

    if (this%naux > 0) then
      text = '       AUXILIARY'
      idx = idx + 1
      call this%budobj%budterm(idx)%initialize(text, this%name_model, &
                                               this%packName, this%name_model, &
                                               this%packName, this%nodes, &
                                               .false., .false., this%naux, &
                                               this%auxname)
    end if

    if (this%iprflow /= 0) then
      call this%budobj%flowtable_df(this%iout, cellids='GWF')
    end if

  end subroutine vssb_setup_budobj

  !> @brief Copy the flow terms into the budget object
  !!
  !! The package budget is closed for each time step:
  !! PRECIPITATION = GWF + RUNOFF + RUNOFF-TO-MVR + STORAGE + EVAPORATION
  !! and flows out of the package are negative.  The GWF term is the net
  !! exchange with the aquifer, which includes the bare-soil evaporation that is
  !! extracted from the top cell.
  !<
  subroutine vssb_fill_budobj(this)
    class(VssbType) :: this
    ! local
    integer(I4B) :: idx, i
    real(DP) :: q
    real(DP), dimension(1) :: vol

    idx = 0

    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nodes)
    do i = 1, this%nodes
      call this%budobj%budterm(idx)%update_term(i, i, this%qp(i))
    end do

    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nodes)
    do i = 1, this%nodes
      call this%budobj%budterm(idx)%update_term(i, this%igwfnode(i), &
                                                -(this%qinf(i) - &
                                                  this%qevap_soil(i)))
    end do

    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nodes)
    do i = 1, this%nodes
      q = -(this%qro(i) - this%qrotomvr(i))
      call this%budobj%budterm(idx)%update_term(i, i, q)
    end do

    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nodes)
    do i = 1, this%nodes
      vol(1) = this%pond(i) * this%area(i)
      call this%budobj%budterm(idx)%update_term(i, i, -this%qsto(i), vol)
    end do

    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%nodes)
    do i = 1, this%nodes
      call this%budobj%budterm(idx)%update_term(i, i, &
                                                -(this%qevap_pond(i) + &
                                                  this%qevap_soil(i)))
    end do

    if (this%imover == 1) then
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%nodes)
      do i = 1, this%nodes
        call this%budobj%budterm(idx)%update_term(i, i, -this%qrotomvr(i))
      end do
    end if

    if (this%naux > 0) then
      idx = idx + 1
      call this%budobj%budterm(idx)%reset(this%nodes)
      do i = 1, this%nodes
        call this%budobj%budterm(idx)%update_term(i, i, DZERO, &
                                                  this%auxvar(:, i))
      end do
    end if

    call this%budobj%accumulate_terms()

  end subroutine vssb_fill_budobj

  !> @brief Define the list heading that is written when PRINT_INPUT is used
  !<
  subroutine define_listlabel(this)
    class(VssbType), intent(inout) :: this

    this%listlabel = trim(this%filtyp)//' NO.'
    if (this%dis%ndim == 3) then
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'LAYER'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'ROW'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'COL'
    elseif (this%dis%ndim == 2) then
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'LAYER'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'CELL2D'
    else
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'NODE'
    end if
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'MAX_POND_DEPTH'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDNAME'
    end if

  end subroutine define_listlabel

  !> @brief Deallocate package memory
  !<
  subroutine vssb_da(this)
    class(VssbType) :: this

    if (associated(this%budobj)) then
      call this%budobj%budgetobject_da()
      deallocate (this%budobj)
      nullify (this%budobj)
    end if

    deallocate (this%vssbname)

    call mem_deallocate(this%nodes)
    call mem_deallocate(this%ioutpond)
    call mem_deallocate(this%igwfnode)
    call mem_deallocate(this%maxpond)
    call mem_deallocate(this%area)
    call mem_deallocate(this%pond0)
    call mem_deallocate(this%pond)
    call mem_deallocate(this%pondcf)
    call mem_deallocate(this%qp)
    call mem_deallocate(this%qinf)
    call mem_deallocate(this%qro)
    call mem_deallocate(this%qrotomvr)
    call mem_deallocate(this%qsto)
    call mem_deallocate(this%qevap_pond)
    call mem_deallocate(this%qevap_soil)
    call mem_deallocate(this%deriv)
    call mem_deallocate(this%pondout)

    call this%BndExtType%bnd_da()

  end subroutine vssb_da

  ! -- Procedures related to observations

  !> @brief Return true because the VSSB package supports observations
  !<
  logical function vssb_obs_supported(this)
    class(VssbType) :: this

    vssb_obs_supported = .true.

  end function vssb_obs_supported

  !> @brief Store the observation types supported by the VSSB package
  !<
  subroutine vssb_df_obs(this)
    class(VssbType) :: this
    ! local
    integer(I4B) :: indx

    call this%obs%StoreObsType('precipitation', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => vssb_process_obsID

    call this%obs%StoreObsType('infiltration', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => vssb_process_obsID

    call this%obs%StoreObsType('runoff', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => vssb_process_obsID

    call this%obs%StoreObsType('runoff-to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => vssb_process_obsID

    call this%obs%StoreObsType('pond-storage', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => vssb_process_obsID

    call this%obs%StoreObsType('evaporation', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => vssb_process_obsID

    call this%obs%StoreObsType('pond-depth', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => vssb_process_obsID

  end subroutine vssb_df_obs

  !> @brief Calculate the observations of this time step
  !!
  !! Flows into the package are positive and flows out of the package are
  !! negative, except for the infiltration, which is positive into the aquifer.
  !<
  subroutine vssb_bd_obs(this)
    class(VssbType) :: this
    ! local
    class(ObserveType), pointer :: obsrv => null()
    integer(I4B) :: i, ii, n
    real(DP) :: v

    if (this%obs%npakobs > 0) then
      call this%obs%obs_bd_clear()
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv
        do ii = 1, obsrv%indxbnds_count
          n = obsrv%indxbnds(ii)
          v = DNODATA
          select case (obsrv%ObsTypeId)
          case ('PRECIPITATION')
            v = this%qp(n)
          case ('INFILTRATION')
            v = this%qinf(n)
          case ('RUNOFF')
            v = -(this%qro(n) - this%qrotomvr(n))
          case ('RUNOFF-TO-MVR')
            if (this%imover == 1) v = -this%qrotomvr(n)
          case ('POND-STORAGE')
            v = -this%qsto(n)
          case ('EVAPORATION')
            v = -(this%qevap_pond(n) + this%qevap_soil(n))
          case ('POND-DEPTH')
            v = this%pond(n)
          case default
            errmsg = 'Unrecognized observation type: '//trim(obsrv%ObsTypeId)
            call store_error(errmsg)
          end select
          call this%obs%SaveOneSimval(obsrv, v)
        end do
      end do
      if (count_errors() > 0) call store_error_filename(this%input_fname)
    end if

  end subroutine vssb_bd_obs

  !> @brief Process each observation
  !!
  !! Only done in the first stress period since the surface cells are fixed for
  !! the simulation.
  !<
  subroutine vssb_rp_obs(this)
    use TdisModule, only: kper
    class(VssbType), intent(inout) :: this
    ! local
    class(ObserveType), pointer :: obsrv => null()
    character(len=LENBOUNDNAME) :: bname
    integer(I4B) :: i, j, nn

    if (kper /= 1) return

    do i = 1, this%obs%npakobs
      obsrv => this%obs%pakobs(i)%obsrv
      nn = obsrv%NodeNumber
      if (nn == NAMEDBOUNDFLAG) then
        bname = obsrv%FeatureName
        do j = 1, this%nodes
          if (this%vssbname(j) == bname) then
            obsrv%BndFound = .true.
            obsrv%CurrentTimeStepEndValue = DZERO
            call obsrv%AddObsIndex(j)
          end if
        end do
        ! a state cannot be summed over several cells
        if (obsrv%ObsTypeId == 'POND-DEPTH' .and. obsrv%indxbnds_count > 1) then
          write (errmsg, '(a,3(1x,a))') &
            trim(adjustl(obsrv%ObsTypeId)), 'for observation', &
            trim(adjustl(obsrv%Name)), &
            'must be assigned to a surface cell with a unique boundname.'
          call store_error(errmsg)
        end if
      else
        if (nn < 1 .or. nn > this%nodes) then
          write (errmsg, '(a,1x,a,1x,i0,1x,a,1x,i0,a)') &
            trim(adjustl(obsrv%ObsTypeId)), &
            'surface cell number must be greater than 0 and less than or '// &
            'equal to', this%nodes, '(specified value is', nn, ').'
          call store_error(errmsg)
        else
          obsrv%BndFound = .true.
          obsrv%CurrentTimeStepEndValue = DZERO
          call obsrv%AddObsIndex(nn)
        end if
      end if
    end do

    if (count_errors() > 0) call store_error_unit(this%inunit)

  end subroutine vssb_rp_obs

  !> @brief Process the ID string of an observation definition
  !!
  !! Pointed to by ObsDataType%ProcessIdPtr.  The ID is a surface cell number
  !! or a boundname.
  !<
  subroutine vssb_process_obsID(obsrv, dis, inunitobs, iout)
    type(ObserveType), intent(inout) :: obsrv
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: inunitobs
    integer(I4B), intent(in) :: iout
    ! local
    character(len=LINELENGTH) :: string
    integer(I4B) :: n, nn, icol, istart, istop, istat
    real(DP) :: r
    ! formats
30  format(i10)

    string = obsrv%IDstring
    icol = 1
    call urword(string, icol, istart, istop, 1, n, r, iout, inunitobs)
    read (string(istart:istop), 30, iostat=istat) nn
    if (istat == 0) then
      obsrv%NodeNumber = nn
    else
      ! not an integer, so the ID is a boundname (already upper case)
      obsrv%FeatureName = string(istart:istop)
      obsrv%NodeNumber = NAMEDBOUNDFLAG
    end if

  end subroutine vssb_process_obsID

end module VssbModule
