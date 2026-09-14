!> @brief Period block keystring-based input loader
!!
!! Each keystring item maps to a typed column in a StructArrayType.
!! A dispatch keyword on each input row selects the target column.
!!
!!   Simple dispatch: keyword matches a DOUBLE/STRING/INTEGER column;
!!   one value token is read into that column.
!!
!!   Compound dispatch: keyword matches a KEYWORD-type column (e.g.
!!   FLOWING_WELL).  The keyword token is stored directly; subsequent
!!   non-KEYWORD body columns are read in order.
!!
!<
module Mf6FileKeystringModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENVARNAME, DZERO, DNODATA
  use InputDefinitionModule, only: InputParamDefinitionType
  use ModflowInputModule, only: ModflowInputType
  use CharacterStringModule, only: CharacterStringType
  use MemoryManagerModule, only: mem_setptr, get_isize
  use TimeSeriesManagerModule, only: TimeSeriesManagerType, tsmanager_cr, &
                                     read_value_or_time_series_adv, &
                                     remove_existing_link
  use StructArrayModule, only: StructArrayType, constructStructArray, &
                               idm_input_varname, is_auxval, &
                               find_auxname_index, destructStructArray
  use AsciiInputLoadTypeModule, only: AsciiDynamicPkgLoadBaseType
  use LoadContextModule, only: LoadContextType, is_advanced, &
                               KeystringItemType, ADDR_FEATURE, ADDR_NODE, &
                               ADDR_INDEXED_BODY
  use LoadMf6FileModule, only: LoadMf6FileType
  use BlockParserModule, only: BlockParserType

  implicit none
  private
  public :: KeystringLoadType

  !> @brief Wraps an allocatable int array for a per-column cache table.
  !<
  type :: IntArrayType
    integer(I4B), dimension(:), allocatable :: vals
  end type IntArrayType

  !> @brief Keystring period block loader
  !!
  !! Leading fixed columns (e.g. CELLID) followed by a dispatch
  !! keyword that routes each input row to a typed item column.
  !!
  !<
  type, extends(AsciiDynamicPkgLoadBaseType) :: KeystringLoadType
    type(TimeSeriesManagerType), pointer :: tsmanager => null()
    type(StructArrayType), pointer :: structarray => null()
    type(LoadContextType) :: ctx !< input load context
    type(LoadMf6FileType) :: static_loader !< persistent static loader
    logical(LGP) :: ts_active !< .true. if TS files are loaded
    integer(I4B) :: nleading !< number of leading (pre-keystring) columns
    ! cached once in allocate_indexed_body (df()-time), since
    ! ctx%indexed_body_dependency's named dimension may be released by rp()-time
    integer(I4B), dimension(:), allocatable :: indexed_body_nfeatures !< per-column feature count (0 = n/a)
    type(IntArrayType), dimension(:), allocatable :: indexed_body_offsets !< per-column offset table
    integer(I4B), dimension(:), allocatable :: indexed_body_index_icol !< per-column sibling index SA column
    integer(I4B), dimension(:), allocatable :: indexed_body_head_icol !< per-column indexed-body head SA column
    character(len=LENVARNAME) :: indexed_body_id_varname = '' !< leading column's mf6varname
  contains
    procedure :: ainit
    procedure :: df
    procedure :: ts_advance
    procedure :: rp
    procedure :: allocate_indexed_body
    procedure :: allocate_permanent_array
    procedure :: allocate_items
    procedure, private :: valid_address
    procedure :: apply_auxiliary
    procedure :: apply_indexed_body
    procedure :: apply_items
    procedure :: resolve_row_addr
    procedure :: resolve_indexed_body_offsets
    procedure :: resolve_indexed_body_address
    procedure :: reset
    procedure :: destroy
    procedure :: create_structarray
  end type KeystringLoadType

contains

  subroutine ainit(this, mf6_input, component_name, component_input_name, &
                   input_name, iperblock, parser, iout)
    use InputOutputModule, only: getunit
    use MemoryManagerModule, only: get_isize, mem_setptr
    use CharacterStringModule, only: CharacterStringType
    use LoadMf6FileModule, only: LoadMf6FileType
    class(KeystringLoadType), intent(inout) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), intent(in) :: component_name
    character(len=*), intent(in) :: component_input_name
    character(len=*), intent(in) :: input_name
    integer(I4B), intent(in) :: iperblock
    type(BlockParserType), pointer, intent(inout) :: parser
    integer(I4B), intent(in) :: iout
    type(CharacterStringType), dimension(:), pointer, contiguous :: ts_fnames
    character(len=LINELENGTH) :: fname
    character(len=LENVARNAME) :: named_bound
    logical(LGP) :: has_named_bound
    integer(I4B) :: n, isize

    call this%DynamicPkgLoadType%init(mf6_input, component_name, &
                                      component_input_name, input_name, &
                                      iperblock, iout)
    this%ts_active = .false.
    this%nleading = 0

    allocate (this%tsmanager)
    call tsmanager_cr(this%tsmanager, iout)

    ! load static input (TS6_FILENAME tag sets static_loader%ts_active)
    call this%static_loader%load(parser, mf6_input, this%nc_vars, &
                                 this%input_name, iout)

    ! add declared TS files to tsmanager
    if (this%static_loader%ts_active) then
      this%ts_active = .true.
      call get_isize('TS6_FILENAME', mf6_input%mempath, isize)
      if (isize > 0) then
        call mem_setptr(ts_fnames, 'TS6_FILENAME', mf6_input%mempath)
        do n = 1, size(ts_fnames)
          fname = ts_fnames(n)
          call this%tsmanager%add_tsfile(fname, getunit())
        end do
      end if
    end if

    ! find a DIMENSIONS parameter to alias as maxbound (skipped for
    ! advanced packages, which get their count from PACKAGEDATA)
    has_named_bound = .false.
    if (.not. is_advanced(mf6_input)) then
      do n = 1, size(mf6_input%param_dfns)
        if (mf6_input%param_dfns(n)%blockname == 'DIMENSIONS') then
          named_bound = trim(mf6_input%param_dfns(n)%mf6varname)
          has_named_bound = .true.
          exit
        end if
      end do
    end if

    ! init load context
    if (has_named_bound) then
      call this%ctx%init(mf6_input, named_bound=named_bound)
    else
      call this%ctx%init(mf6_input)
    end if

    ! params is fully elaborated: leading cols + item names
    this%param_names = this%ctx%params
    this%nparam = size(this%ctx%params)
    this%nleading = this%ctx%nleading
    call this%ctx%check_developmode(this%input_name)

    ! finalize context setup (allocates NBOUND, NODEULIST, etc.)
    call this%ctx%allocate_arrays()

    ! pre-allocate structarray; reused across all periods
    call this%create_structarray()

  end subroutine ainit

  subroutine df(this)
    use StructArrayModule, only: StructArrayType
    use MemoryManagerModule, only: mem_setptr, get_isize
    use CharacterStringModule, only: CharacterStringType
    class(KeystringLoadType), intent(inout) :: this
    type(StructArrayType), pointer :: sa
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      auxnames => null()
    integer(I4B), dimension(:), pointer, contiguous :: pkg_ifno => null()
    integer(I4B) :: n, naux
    ! init tsmanager (TDIS now available)
    call this%tsmanager%tsmanager_df()
    ! resolve aux names for PACKAGEDATA AUX TS registration
    call get_isize('AUXILIARY', this%mf6_input%mempath, naux)
    if (naux > 0) call mem_setptr(auxnames, 'AUXILIARY', this%mf6_input%mempath)
    ! advanced packages: address AUX TS links by feature number (not
    ! PACKAGEDATA row position), so a later PERIOD override finds it
    if (this%ctx%is_advanced) then
      call mem_setptr(pkg_ifno, 'PACKAGEDATA_IFNO', this%mf6_input%mempath)
    end if
    ! link static TS strlocs; preserve for re-registration after reset()
    do n = 1, this%static_loader%ts_sa_count()
      sa => this%static_loader%get_ts_sa(n)
      if (associated(sa)) then
        if (associated(pkg_ifno)) then
          call sa%ts_update(this%tsmanager, &
                            this%mf6_input%subcomponent_name, &
                            this%ctx%iprpak, this%input_name, &
                            clear_strlocs=.false., auxname_cst=auxnames, &
                            ifno_map=pkg_ifno)
        else
          call sa%ts_update(this%tsmanager, &
                            this%mf6_input%subcomponent_name, &
                            this%ctx%iprpak, this%input_name, &
                            clear_strlocs=.false., auxname_cst=auxnames)
        end if
      end if
    end do
    ! allocate feature-addressed (DZERO) and node-addressed (DNODATA) items
    call this%allocate_items()
    ! indexed-body df-time dimension fields + permanent arrays
    if (this%ctx%is_advanced) call this%allocate_indexed_body()
  end subroutine df

  subroutine ts_advance(this)
    class(KeystringLoadType), intent(inout) :: this
    call this%tsmanager%ad()
  end subroutine ts_advance

  subroutine rp(this, parser)
    use IdmLoggerModule, only: idm_log_header, idm_log_close
    class(KeystringLoadType), intent(inout) :: this
    type(BlockParserType), pointer, intent(inout) :: parser

    call this%reset()

    call idm_log_header(this%mf6_input%component_name, &
                        this%mf6_input%subcomponent_name, this%iout)

    this%ctx%nbound = &
      this%structarray%read_from_parser_keystring(parser, this%ts_active, &
                                                  this%nleading, this%iout, &
                                                  this%input_name)

    if (this%ctx%is_advanced) call this%apply_auxiliary()
    if (this%ctx%is_advanced) call this%apply_indexed_body()
    ! apply feature- and node-addressed items
    call this%apply_items()

    if (this%ts_active) then
      call this%structarray%ts_update(this%tsmanager, &
                                      this%mf6_input%subcomponent_name, &
                                      this%ctx%iprpak, this%input_name)
    end if

    call idm_log_close(this%mf6_input%component_name, &
                       this%mf6_input%subcomponent_name, this%iout)
  end subroutine rp

  !> @brief Allocate ctx%indexed_body_dependency targets (indexed-body
  !! params excluded from allocate_items), and cache
  !! everything apply_indexed_body needs every period.
  !<
  subroutine allocate_indexed_body(this)
    use DefinitionSelectModule, only: get_param_definition_type
    class(KeystringLoadType), intent(inout) :: this
    type(InputParamDefinitionType), pointer :: idt, id_idt
    character(len=LENVARNAME) :: dimname, index_tagname
    integer(I4B) :: icol, sa_icol, padj, nfeatures, n
    logical(LGP) :: found

    if (.not. allocated(this%indexed_body_nfeatures)) then
      allocate (this%indexed_body_nfeatures(this%structarray%count()))
      allocate (this%indexed_body_offsets(this%structarray%count()))
      allocate (this%indexed_body_index_icol(this%structarray%count()))
      allocate (this%indexed_body_head_icol(this%structarray%count()))
      this%indexed_body_nfeatures = 0
      this%indexed_body_index_icol = 0
      this%indexed_body_head_icol = 0
    end if

    padj = 0
    if (this%ctx%has_setting_dispatch) padj = 1

    id_idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                        this%mf6_input%component_type, &
                                        this%mf6_input%subcomponent_type, &
                                        'PERIOD', this%param_names(1), &
                                        this%input_name)
    this%indexed_body_id_varname = trim(id_idt%mf6varname)

    do icol = this%nleading + 1, this%nparam
      if (.not. this%ctx%items(icol - this%nleading)%is_body) cycle
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       'PERIOD', this%param_names(icol), &
                                       this%input_name)
      found = this%ctx%indexed_body_dependency(idt%tagname, dimname, &
                                               index_tagname)
      if (.not. found) cycle

      sa_icol = icol + padj
      nfeatures = this%ctx%resolve_item_nfeatures(idt%tagname, dimname, 0)
      if (nfeatures < 1) cycle
      this%indexed_body_nfeatures(sa_icol) = nfeatures
      this%indexed_body_offsets(sa_icol)%vals = &
        this%resolve_indexed_body_offsets(dimname)
      call this%allocate_permanent_array( &
        trim(idt%mf6varname), nfeatures, DZERO)

      ! resolve the sibling index column (by ctx-given tag) and this
      ! target's indexed-body head, both structurally (body_start/head_nbody)
      do n = 1, this%structarray%count()
        if (trim(this%structarray%struct_vectors(n)%idt%tagname) == &
            trim(index_tagname)) this%indexed_body_index_icol(sa_icol) = n
        if (this%structarray%struct_vectors(n)%head_nbody > 0) then
          if (sa_icol >= this%structarray%struct_vectors(n)%body_start .and. &
              sa_icol < this%structarray%struct_vectors(n)%body_start + &
              this%structarray%struct_vectors(n)%head_nbody) &
            this%indexed_body_head_icol(sa_icol) = n
        end if
      end do
      ! neither sibling found (misconfigured ctx entry): not a target
      if (this%indexed_body_index_icol(sa_icol) == 0 .or. &
          this%indexed_body_head_icol(sa_icol) == 0) &
        this%indexed_body_nfeatures(sa_icol) = 0
    end do
  end subroutine allocate_indexed_body

  !> @brief Allocate a permanent per-feature array named `name` with
  !! init_value, unless already allocated.
  !<
  subroutine allocate_permanent_array(this, name, nfeatures, init_value)
    use MemoryManagerModule, only: mem_allocate
    class(KeystringLoadType), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer(I4B), intent(in) :: nfeatures
    real(DP), intent(in) :: init_value
    real(DP), dimension(:), pointer, contiguous :: featarr => null()
    integer(I4B) :: isize

    call get_isize(name, this%mf6_input%mempath, isize)
    if (isize > 0) return ! already allocated (shouldn't happen; df() runs once)
    call mem_allocate(featarr, nfeatures, name, this%mf6_input%mempath)
    featarr = init_value
  end subroutine allocate_permanent_array

  !> @brief Allocate a permanent array for every IDM-managed item, using
  !! the descriptor's per-item feature count and init value.
  !<
  subroutine allocate_items(this)
    use SimModule, only: count_errors, store_error_filename
    class(KeystringLoadType), intent(inout) :: this
    integer(I4B) :: k

    if (.not. allocated(this%ctx%items)) return

    do k = 1, size(this%ctx%items)
      if (.not. this%ctx%items(k)%idm_managed) cycle
      ! indexed-body arrays are allocated in allocate_indexed_body (df-time
      ! dimension); skip here
      if (this%ctx%items(k)%addr_mode == ADDR_INDEXED_BODY) cycle
      if (this%ctx%items(k)%nfeatures < 1) cycle
      call this%allocate_permanent_array(this%ctx%items(k)%mf6varname, &
                                         this%ctx%items(k)%nfeatures, &
                                         this%ctx%items(k)%init_value)
    end do
    if (count_errors() > 0) then
      call store_error_filename(this%input_name)
    end if
  end subroutine allocate_items

  !> @brief Validate ifno against nfeatures, storing an error keyed by
  !! ifno_tagname (the leading column's public tag) if out of range.
  !<
  function valid_address(this, ifno, nfeatures, ifno_tagname, row) result(valid)
    use SimModule, only: store_error
    use SimVariablesModule, only: errmsg
    class(KeystringLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: ifno
    integer(I4B), intent(in) :: nfeatures
    character(len=*), intent(in) :: ifno_tagname
    integer(I4B), intent(in) :: row
    logical(LGP) :: valid

    valid = (ifno >= 1 .and. ifno <= nfeatures)
    if (.not. valid) then
      write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a,1x,i0,a)') &
        trim(ifno_tagname), ifno, 'on row', row, &
        'must be greater than 0 and less than or equal to', nfeatures, '.'
      call store_error(errmsg)
    end if
  end function valid_address

  !> @brief Apply PERIOD AUXILIARY settings to the permanent AUX array.
  !! TS-linked rows resolve via the struct array's own ts_strlocs;
  !! literal rows resolve AUXNAME here and clear any stale TS link.
  !<
  subroutine apply_auxiliary(this)
    use DefinitionSelectModule, only: get_param_definition_type
    use StructVectorModule, only: TSStringLocType
    use SimModule, only: count_errors, store_error_filename
    class(KeystringLoadType), intent(inout) :: this
    integer(I4B), pointer :: nbound => null()
    integer(I4B), dimension(:), pointer, contiguous :: period_ifno => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      period_setting => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      period_auxname => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      auxnames => null()
    real(DP), dimension(:, :), pointer, contiguous :: aux => null()
    real(DP), pointer :: bndElem
    type(InputParamDefinitionType), pointer :: idt
    type(TSStringLocType), pointer :: ts_strloc
    integer(I4B) :: i, n, ifno, jj, isize, naux, nfeatures, sa_icol, k, nts
    logical(LGP) :: found
    logical(LGP), dimension(:), allocatable :: handled
    character(len=LINELENGTH) :: setting, auxname, thisauxname
    character(len=LENVARNAME) :: ifno_tagname

    call get_isize('AUXILIARY', this%mf6_input%mempath, naux)
    if (naux <= 0) return

    call get_isize('NBOUND', this%mf6_input%mempath, isize)
    if (isize < 1) return
    call mem_setptr(nbound, 'NBOUND', this%mf6_input%mempath)
    if (nbound <= 0) return

    call get_isize('AUXNAME', this%mf6_input%mempath, isize)
    if (isize < 1) return

    ! leading column's public tag (e.g. MAWNO for MWE), for the error
    ! message below -- its memory-manager key is always IFNO (MF6INTERNAL)
    idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                     this%mf6_input%component_type, &
                                     this%mf6_input%subcomponent_type, &
                                     'PERIOD', this%param_names(1), &
                                     this%input_name)
    ifno_tagname = trim(idt%tagname)

    ! AUXILIARY dispatch keyword's own mf6varname (e.g. LAK's
    ! PERIOD_AUXILIARY), since SETTING stores mf6varname, not tagname
    idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                     this%mf6_input%component_type, &
                                     this%mf6_input%subcomponent_type, &
                                     'PERIOD', 'AUXILIARY', this%input_name)

    call mem_setptr(period_ifno, 'IFNO', this%mf6_input%mempath)
    call mem_setptr(period_setting, 'SETTING', this%mf6_input%mempath)
    call mem_setptr(period_auxname, 'AUXNAME', this%mf6_input%mempath)
    call mem_setptr(auxnames, 'AUXILIARY', this%mf6_input%mempath)
    call mem_setptr(aux, 'AUX', this%mf6_input%mempath)
    nfeatures = size(aux, 2)

    sa_icol = 0
    do n = 1, this%structarray%count()
      if (is_auxval(this%structarray%struct_vectors(n)%idt)) then
        sa_icol = n
        exit
      end if
    end do
    if (sa_icol == 0) return

    allocate (handled(nbound))
    handled = .false.

    ! TS-linked rows: resolve AUXNAME here too, same as literal rows
    nts = this%structarray%struct_vectors(sa_icol)%ts_strlocs%count()
    do k = 1, nts
      ts_strloc => this%structarray%struct_vectors(sa_icol)%get_ts_strloc(k)
      i = ts_strloc%row
      ifno = period_ifno(i)
      if (.not. this%valid_address(ifno, nfeatures, ifno_tagname, i)) cycle
      auxname = period_auxname(i)
      jj = find_auxname_index(auxname, auxnames, naux)
      if (jj < 1) cycle
      thisauxname = auxnames(jj)
      bndElem => aux(jj, ifno)
      call read_value_or_time_series_adv(ts_strloc%token, ifno, jj, bndElem, &
                                         this%mf6_input%subcomponent_name, &
                                         'AUX', this%tsmanager, &
                                         this%ctx%iprpak, trim(thisauxname))
      handled(i) = .true.
    end do

    ! literal rows: resolve AUXNAME here, clear any stale link, assign
    do i = 1, nbound
      if (handled(i)) cycle
      setting = period_setting(i)
      if (trim(setting) /= trim(idt%mf6varname)) cycle
      ifno = period_ifno(i)
      if (.not. this%valid_address(ifno, nfeatures, ifno_tagname, i)) cycle
      auxname = period_auxname(i)
      jj = find_auxname_index(auxname, auxnames, naux)
      if (jj < 1) cycle
      thisauxname = auxnames(jj)
      found = remove_existing_link(this%tsmanager, ifno, jj, &
                                   this%mf6_input%subcomponent_name, &
                                   'AUX', trim(thisauxname))
      aux(jj, ifno) = this%structarray%struct_vectors(sa_icol)%dbl1d(i)
    end do
    if (count_errors() > 0) then
      call store_error_filename(this%input_name)
    end if
    deallocate (handled)
    call this%structarray%struct_vectors(sa_icol)%clear()
  end subroutine apply_auxiliary

  !> @brief Apply PERIOD settings for ctx%indexed_body_dependency's targets,
  !! via ts_update_indexed -- same mechanism as BEDK/MANNING, but each
  !! row's target index comes from a named sibling field through the
  !! cached offset table rather than the leading id column alone.
  !<
  subroutine apply_indexed_body(this)
    class(KeystringLoadType), intent(inout) :: this
    type(InputParamDefinitionType), pointer :: idt, head_idt
    real(DP), dimension(:), pointer, contiguous :: featarr => null()
    integer(I4B), pointer :: nbound => null()
    integer(I4B), dimension(:), allocatable :: row_addr
    integer(I4B) :: icol, sa_icol, padj, isize

    if (.not. allocated(this%indexed_body_nfeatures)) return

    call get_isize('NBOUND', this%mf6_input%mempath, isize)
    if (isize < 1) return
    call mem_setptr(nbound, 'NBOUND', this%mf6_input%mempath)
    if (nbound <= 0) return

    padj = 0
    if (this%ctx%has_setting_dispatch) padj = 1

    do icol = this%nleading + 1, this%nparam
      sa_icol = icol + padj
      if (sa_icol > size(this%indexed_body_nfeatures)) cycle
      if (this%indexed_body_nfeatures(sa_icol) < 1) cycle

      idt => this%structarray%struct_vectors(sa_icol)%idt
      head_idt => &
        this%structarray%struct_vectors(this%indexed_body_head_icol(sa_icol))%idt
      call mem_setptr(featarr, trim(idt%mf6varname), this%mf6_input%mempath)
      row_addr = &
        this%resolve_indexed_body_address( &
        this%indexed_body_id_varname, trim(head_idt%mf6varname), &
        this%indexed_body_index_icol(sa_icol), &
        this%indexed_body_offsets(sa_icol)%vals, &
        this%indexed_body_nfeatures(sa_icol), nbound)
      call this%structarray%ts_update_indexed( &
        sa_icol, this%tsmanager, this%mf6_input%subcomponent_name, &
        this%ctx%iprpak, nbound, row_addr, trim(idt%tagname), featarr)
    end do
  end subroutine apply_indexed_body

  !> @brief Apply IDM-managed feature- and node-addressed settings to their
  !! permanent arrays; each row's address comes from resolve_row_addr.
  !! Indexed-body and AUX are handled in their own routines.
  !<
  subroutine apply_items(this)
    use SimModule, only: count_errors, store_error_filename
    use DefinitionSelectModule, only: get_param_definition_type
    class(KeystringLoadType), intent(inout) :: this
    integer(I4B), pointer :: nbound => null()
    integer(I4B), dimension(:, :), pointer, contiguous :: cellid => null()
    integer(I4B), dimension(:), pointer, contiguous :: period_ifno => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      period_setting => null()
    real(DP), dimension(:), pointer, contiguous :: featarr => null()
    integer(I4B), dimension(:), allocatable :: row_addr
    integer(I4B) :: i, k, isize, ndim
    character(len=LINELENGTH) :: setting
    character(len=LENVARNAME) :: dispatch_key
    character(len=LENVARNAME) :: leading_tag
    type(InputParamDefinitionType), pointer :: lead_idt

    if (.not. allocated(this%ctx%items)) return

    call get_isize('NBOUND', this%mf6_input%mempath, isize)
    if (isize < 1) return
    call mem_setptr(nbound, 'NBOUND', this%mf6_input%mempath)
    if (nbound <= 0) return

    call mem_setptr(period_setting, 'SETTING', this%mf6_input%mempath)

    ! addressing sources, resolved once
    ndim = 0
    leading_tag = ''
    if (this%ctx%keystring_by_node) then
      if (.not. associated(this%ctx%mshape)) return
      ndim = size(this%ctx%mshape)
      call mem_setptr(cellid, 'CELLID', this%mf6_input%mempath)
    else
      call mem_setptr(period_ifno, this%ctx%feature_id_varname, &
                      this%mf6_input%mempath)
      ! leading id column's public tag (e.g. MAWNO), for the out-of-range
      ! message -- the invalid value is the row's leading id, not the
      ! setting. Feature-addressed only; node loads validate via CELLID.
      lead_idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                            this%mf6_input%component_type, &
                                            this%mf6_input%subcomponent_type, &
                                            'PERIOD', this%param_names(1), &
                                            this%input_name)
      leading_tag = trim(lead_idt%tagname)
    end if

    do k = 1, size(this%ctx%items)
      if (.not. this%ctx%items(k)%idm_managed) cycle
      if (this%ctx%items(k)%addr_mode == ADDR_INDEXED_BODY) cycle ! separate path
      if (this%ctx%items(k)%nfeatures < 1) cycle
      call mem_setptr(featarr, this%ctx%items(k)%mf6varname, &
                      this%mf6_input%mempath)

      ! dispatch key: a RECORD body matches its owning head's SETTING
      ! keyword (the token on the row); a standalone item matches its own
      ! mf6varname. resolve_row_addr maps the matched row to a feature.
      if (this%ctx%items(k)%is_body) then
        dispatch_key = this%ctx%items(k)%head_setting_varname
      else
        dispatch_key = this%ctx%items(k)%mf6varname
      end if

      allocate (row_addr(nbound))
      do i = 1, nbound
        row_addr(i) = 0
        setting = period_setting(i)
        if (trim(setting) /= trim(dispatch_key)) cycle
        row_addr(i) = this%resolve_row_addr(this%ctx%items(k), i, period_ifno, &
                                            cellid, ndim, leading_tag)
      end do
      call this%structarray%ts_update_indexed( &
        this%ctx%items(k)%sa_icol, this%tsmanager, &
        this%mf6_input%subcomponent_name, this%ctx%iprpak, nbound, row_addr, &
        this%ctx%items(k)%tag, featarr)
      deallocate (row_addr)
    end do
    if (count_errors() > 0) then
      call store_error_filename(this%input_name)
    end if
  end subroutine apply_items

  !> @brief Resolve one row's target index for an item, dispatching on
  !! addr_mode. Returns 0 (skip) if out of range.
  !<
  function resolve_row_addr(this, item, irow, period_ifno, cellid, ndim, &
                            leading_tag) &
    result(addr)
    use GeomUtilModule, only: get_node
    class(KeystringLoadType), intent(inout) :: this
    type(KeystringItemType), intent(in) :: item
    integer(I4B), intent(in) :: irow
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: period_ifno
    integer(I4B), dimension(:, :), pointer, contiguous, intent(in) :: cellid
    integer(I4B), intent(in) :: ndim
    character(len=*), intent(in) :: leading_tag !< leading id column's public tag (for the out-of-range message)
    integer(I4B) :: addr
    integer(I4B) :: ifno, nodeu

    addr = 0
    select case (item%addr_mode)
    case (ADDR_FEATURE)
      ifno = period_ifno(irow)
      if (.not. this%valid_address(ifno, item%nfeatures, leading_tag, irow)) &
        return
      addr = ifno
    case (ADDR_NODE)
      if (ndim == 1) then
        nodeu = cellid(1, irow)
      else if (ndim == 2) then
        nodeu = get_node(cellid(1, irow), 1, cellid(2, irow), &
                         this%ctx%mshape(1), 1, this%ctx%mshape(2))
      else
        nodeu = get_node(cellid(1, irow), cellid(2, irow), cellid(3, irow), &
                         this%ctx%mshape(1), this%ctx%mshape(2), &
                         this%ctx%mshape(3))
      end if
      if (nodeu < 1 .or. nodeu > item%nfeatures) return
      addr = nodeu
    end select
  end function resolve_row_addr

  !> @brief Per-feature cumulative offset table from an array-valued
  !! dimension, permuted into feature-index order first. Indexed by
  !! feature (not the dimension's own summed total).
  !<
  function resolve_indexed_body_offsets(this, dimname) result(offsets)
    class(KeystringLoadType), intent(inout) :: this
    character(len=*), intent(in) :: dimname
    integer(I4B), dimension(:), allocatable :: offsets
    integer(I4B), dimension(:), pointer, contiguous :: counts_raw => null()
    integer(I4B), dimension(:), pointer, contiguous :: pkg_ifno => null()
    integer(I4B), dimension(:), allocatable :: counts
    integer(I4B) :: i, n, running, ndomain

    call mem_setptr(counts_raw, trim(dimname), this%mf6_input%mempath)
    call mem_setptr(pkg_ifno, 'PACKAGEDATA_IFNO', this%mf6_input%mempath)
    ndomain = size(counts_raw)
    allocate (counts(ndomain))
    counts = 0
    do i = 1, size(counts_raw)
      n = pkg_ifno(i)
      if (n < 1 .or. n > ndomain) cycle
      counts(n) = counts_raw(i)
    end do

    allocate (offsets(ndomain))
    running = 1
    do n = 1, ndomain
      offsets(n) = running
      running = running + counts(n)
    end do
  end function resolve_indexed_body_offsets

  !> @brief Row -> resolved feature index for an item whose local
  !! position comes from a sibling index field, via the offset table.
  !<
  function resolve_indexed_body_address(this, id_mf6varname, head_mf6varname, &
                                        index_icol, offsets, nfeatures, nbound) &
    result(row_addr)
    use SimModule, only: store_error
    use SimVariablesModule, only: errmsg
    class(KeystringLoadType), intent(inout) :: this
    character(len=*), intent(in) :: id_mf6varname
    character(len=*), intent(in) :: head_mf6varname
    integer(I4B), intent(in) :: index_icol !< SA column holding the sibling local index (e.g. IDV)
    integer(I4B), dimension(:), intent(in) :: offsets
    integer(I4B), intent(in) :: nfeatures
    integer(I4B), intent(in) :: nbound
    integer(I4B), dimension(:), allocatable :: row_addr
    integer(I4B), dimension(:), pointer, contiguous :: period_ifno => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      period_setting => null()
    integer(I4B) :: i, ifno, idx_local, nlocal
    character(len=LINELENGTH) :: setting
    character(len=LENVARNAME) :: id_tag, idx_tag

    call mem_setptr(period_ifno, id_mf6varname, this%mf6_input%mempath)
    call mem_setptr(period_setting, 'SETTING', this%mf6_input%mempath)
    id_tag = trim(this%structarray%struct_vectors(1)%idt%tagname)
    idx_tag = trim(this%structarray%struct_vectors(index_icol)%idt%tagname)
    allocate (row_addr(nbound))
    do i = 1, nbound
      row_addr(i) = 0
      setting = period_setting(i)
      if (trim(setting) /= trim(head_mf6varname)) cycle
      ifno = period_ifno(i)
      if (ifno < 1 .or. ifno > size(offsets)) cycle
      ! per-id local count from the cumulative offset table
      if (ifno < size(offsets)) then
        nlocal = offsets(ifno + 1) - offsets(ifno)
      else
        nlocal = nfeatures - offsets(ifno) + 1
      end if
      idx_local = this%structarray%struct_vectors(index_icol)%int1d(i)
      if (idx_local < 1 .or. idx_local > nlocal) then
        write (errmsg, '(a,1x,a,1x,i0,1x,a,1x,i0,1x,a,1x,a,1x,i0,a)') &
          'index', trim(idx_tag), idx_local, 'must be between 1 and', nlocal, &
          'for', trim(id_tag), ifno, '.'
        call store_error(errmsg)
        cycle
      end if
      row_addr(i) = offsets(ifno) + idx_local - 1
    end do
  end function resolve_indexed_body_address

  subroutine reset(this)
    use StructArrayModule, only: StructArrayType
    use MemoryManagerModule, only: mem_setptr, get_isize
    use CharacterStringModule, only: CharacterStringType
    class(KeystringLoadType), intent(inout) :: this
    type(StructArrayType), pointer :: sa
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      auxnames => null()
    integer(I4B) :: n, naux
    ! every KEYSTRING subtype with SETTING dispatch: PERIOD settings
    ! persist across periods unless reissued, so TS links never reset
    if (this%ctx%has_setting_dispatch) return
    ! clear TS links
    call this%tsmanager%reset(this%mf6_input%subcomponent_name)
    ! re-register static TS links (strlocs preserved in df)
    if (this%ts_active) then
      call get_isize('AUXILIARY', this%mf6_input%mempath, naux)
      if (naux > 0) call mem_setptr(auxnames, 'AUXILIARY', this%mf6_input%mempath)
      do n = 1, this%static_loader%ts_sa_count()
        sa => this%static_loader%get_ts_sa(n)
        if (associated(sa)) then
          call sa%ts_update(this%tsmanager, &
                            this%mf6_input%subcomponent_name, &
                            this%ctx%iprpak, this%input_name, &
                            clear_strlocs=.false., auxname_cst=auxnames)
        end if
      end do
    end if
  end subroutine reset

  subroutine destroy(this)
    class(KeystringLoadType), intent(inout) :: this

    call this%static_loader%cleanup()

    call this%tsmanager%da()
    deallocate (this%tsmanager)
    nullify (this%tsmanager)

    if (associated(this%structarray)) then
      call destructStructArray(this%structarray)
    end if

    call this%ctx%destroy()
    call this%DynamicPkgLoadType%destroy()
  end subroutine destroy

  subroutine create_structarray(this)
    use DefinitionSelectModule, only: get_param_definition_type
    class(KeystringLoadType), intent(inout) :: this
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: icol, sa_icol, nrow_prealloc, nsub, padj
    logical(LGP) :: has_setting

    has_setting = this%ctx%has_setting_dispatch

    ! use pre-allocated managed memory (maxbound = features * nitems);
    ! fall back to deferred shape (-1) if maxbound is unavailable
    if (associated(this%ctx%maxbound) .and. this%ctx%maxbound > 0) then
      nrow_prealloc = this%ctx%maxbound
    else
      nrow_prealloc = -1
    end if

    ! SETTING column inserted at nleading+1 when has_setting
    padj = 0
    if (has_setting) padj = 1

    if (has_setting .and. nrow_prealloc < 0) then
      ! fallback for a genuinely unresolvable count (e.g. empty
      ! PACKAGEDATA), when PACKAGEDATA/DIMENSIONS can't supply one
      this%structarray => &
        constructStructArray(this%mf6_input, this%nparam + padj, &
                             nrow_prealloc, 0, this%mf6_input%mempath, &
                             this%mf6_input%component_mempath, size_init=64)
    else
      this%structarray => &
        constructStructArray(this%mf6_input, this%nparam + padj, &
                             nrow_prealloc, 0, this%mf6_input%mempath, &
                             this%mf6_input%component_mempath)
    end if

    ! create leading (pre-keystring) columns unchanged
    do icol = 1, this%nleading
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       'PERIOD', &
                                       this%param_names(icol), this%input_name)
      call this%structarray%mem_create_vector(icol, idt)
    end do

    ! create SETTING column (ctx owns setting_idt)
    if (has_setting) then
      sa_icol = this%nleading + 1
      call this%structarray%mem_create_vector(sa_icol, this%ctx%setting_idt, &
                                              charlen=LENVARNAME)
    end if

    ! create item columns
    do icol = this%nleading + 1, this%nparam
      sa_icol = icol + padj
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       'PERIOD', &
                                       this%param_names(icol), this%input_name)
      ! nsub from descriptor: 0 = direct dispatch, N = KEYWORD compound with N body members
      nsub = this%ctx%items(icol - this%nleading)%head_nbody
      if (nsub > 0) then
        ! metadata vector: no data allocated; body_start points to next SA col
        call this%structarray%mem_create_metadata_vector(sa_icol, idt, &
                                                         sa_icol + 1, nsub)
      else if (trim(idt%datatype) == 'STRING') then
        ! string value columns (e.g. STATUS) stored at LENVARNAME
        call this%structarray%mem_create_vector(sa_icol, idt, &
                                                charlen=LENVARNAME)
      else if (idt%datatype == 'DOUBLE' .and. &
               (this%ctx%items(icol - this%nleading)%idm_managed .or. &
                this%ctx%items(icol - this%nleading)%addr_mode == &
                ADDR_INDEXED_BODY)) then
        ! managed/indexed-body double: raw read array uses the suffixed input
        ! name, leaving mf6varname for the permanent array the loader populates
        call this%structarray%mem_create_vector(sa_icol, idt, &
                                                varname=idm_input_varname(idt))
      else
        call this%structarray%mem_create_vector(sa_icol, idt)
      end if
    end do
  end subroutine create_structarray

end module Mf6FileKeystringModule
