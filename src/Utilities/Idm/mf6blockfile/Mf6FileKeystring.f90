!> @brief Period block keystring-based input loader
!!
!! Each keystring member maps to a typed column in a StructArrayType.
!! A dispatch keyword on each input row selects the target column.
!!
!!   Simple dispatch: keyword matches a DOUBLE/STRING/INTEGER column;
!!   one value token is read into that column.
!!
!!   Compound dispatch: keyword matches a KEYWORD-type column (e.g.
!!   FLOWING_WELL).  The keyword token is stored directly; subsequent
!!   non-KEYWORD sub-member columns are read in order.
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
  use LoadContextModule, only: LoadContextType, is_advanced
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
  !! keyword that routes each input row to a typed member column.
  !!
  !<
  type, extends(AsciiDynamicPkgLoadBaseType) :: KeystringLoadType
    type(TimeSeriesManagerType), pointer :: tsmanager => null()
    type(StructArrayType), pointer :: structarray => null()
    type(LoadContextType) :: ctx !< input load context
    type(LoadMf6FileType) :: static_loader !< persistent static loader
    logical(LGP) :: ts_active !< .true. if TS files are loaded
    integer(I4B) :: nleading !< number of leading (pre-keystring) columns
    ! cached once in allocate_record_targets (df()-time), since
    ! ctx%record_dependency's named dimension may be released by rp()-time
    integer(I4B), dimension(:), allocatable :: record_nfeatures !< per-column feature count (0 = n/a)
    type(IntArrayType), dimension(:), allocatable :: record_offsets !< per-column offset table
    integer(I4B), dimension(:), allocatable :: record_index_icol !< per-column sibling index SA column
    integer(I4B), dimension(:), allocatable :: record_head_icol !< per-column record head SA column
    character(len=LENVARNAME) :: record_id_varname = '' !< leading column's mf6varname
  contains
    procedure :: ainit
    procedure :: df
    procedure :: ts_advance
    procedure :: rp
    procedure :: allocate_record_targets
    procedure :: allocate_permanent_array
    procedure :: allocate_settings
    procedure :: allocate_node_settings
    procedure, private :: valid_ifno
    procedure :: apply_auxiliary
    procedure :: apply_record_targets
    procedure :: apply_settings
    procedure :: resolve_nfeatures
    procedure :: resolve_member_nfeatures
    procedure, private :: shape_param_is_array
    procedure :: resolve_member_offsets
    procedure :: resolve_member_row_addr
    procedure :: resolve_in_scope_setting
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

    ! params is fully elaborated: leading cols + member names
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
    ! identifier-addressed (advanced + DIMENSIONS-scoped, e.g. SPC):
    ! allocate permanent storage for sticky, TS-continuable persistence
    if (this%ctx%keystring_by_id) call this%allocate_settings()
    ! record-follower targets ctx%record_dependency identifies:
    ! excluded from the generic loop above, allocated here
    if (this%ctx%is_advanced) call this%allocate_record_targets()
    ! CELLID-addressed packages (TVK/TVS): same persistence goal, but
    ! node-indexed and package-resolved -- see allocate_node_settings
    if (this%ctx%keystring_by_node) call this%allocate_node_settings()
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
    if (this%ctx%is_advanced) call this%apply_record_targets()
    if (this%ctx%keystring_by_id) &
      call this%apply_settings(node_addressed=.false.)
    if (this%ctx%keystring_by_node) &
      call this%apply_settings(node_addressed=.true.)

    if (this%ts_active) then
      call this%structarray%ts_update(this%tsmanager, &
                                      this%mf6_input%subcomponent_name, &
                                      this%ctx%iprpak, this%input_name)
    end if

    call idm_log_close(this%mf6_input%component_name, &
                       this%mf6_input%subcomponent_name, this%iout)
  end subroutine rp

  !> @brief Allocate ctx%record_dependency's targets (record-follower
  !! params excluded from allocate_settings), and cache
  !! everything apply_record_targets needs every period.
  !<
  subroutine allocate_record_targets(this)
    use DefinitionSelectModule, only: get_param_definition_type
    class(KeystringLoadType), intent(inout) :: this
    type(InputParamDefinitionType), pointer :: idt, id_idt
    character(len=LENVARNAME) :: dimname, index_tagname
    integer(I4B) :: icol, sa_icol, padj, nfeatures, n
    logical(LGP) :: found

    if (.not. allocated(this%record_nfeatures)) then
      allocate (this%record_nfeatures(this%structarray%count()))
      allocate (this%record_offsets(this%structarray%count()))
      allocate (this%record_index_icol(this%structarray%count()))
      allocate (this%record_head_icol(this%structarray%count()))
      this%record_nfeatures = 0
      this%record_index_icol = 0
      this%record_head_icol = 0
    end if

    padj = 0
    if (this%ctx%has_setting_dispatch) padj = 1

    id_idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                        this%mf6_input%component_type, &
                                        this%mf6_input%subcomponent_type, &
                                        'PERIOD', this%param_names(1), &
                                        this%input_name)
    this%record_id_varname = trim(id_idt%mf6varname)

    do icol = this%nleading + 1, this%nparam
      if (.not. this%ctx%member_is_follower(icol - this%nleading)) cycle
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       'PERIOD', this%param_names(icol), &
                                       this%input_name)
      found = this%ctx%record_dependency(idt%tagname, dimname, index_tagname)
      if (.not. found) cycle

      sa_icol = icol + padj
      nfeatures = this%resolve_member_nfeatures(idt%tagname, dimname, 0)
      if (nfeatures < 1) cycle
      this%record_nfeatures(sa_icol) = nfeatures
      this%record_offsets(sa_icol)%vals = this%resolve_member_offsets(dimname)
      call this%allocate_permanent_array( &
        idt, nfeatures, DZERO, varname=trim(idt%mf6varname)//'_RESOLVED')

      ! resolve the sibling index column (by ctx-given tag) and this
      ! target's record head, both structurally (isubmember/nsubmembers)
      do n = 1, this%structarray%count()
        if (trim(this%structarray%struct_vectors(n)%idt%tagname) == &
            trim(index_tagname)) this%record_index_icol(sa_icol) = n
        if (this%structarray%struct_vectors(n)%nsubmembers > 0) then
          if (sa_icol >= this%structarray%struct_vectors(n)%isubmember .and. &
              sa_icol < this%structarray%struct_vectors(n)%isubmember + &
              this%structarray%struct_vectors(n)%nsubmembers) &
            this%record_head_icol(sa_icol) = n
        end if
      end do
      ! neither sibling found (misconfigured ctx entry): not a target
      if (this%record_index_icol(sa_icol) == 0 .or. &
          this%record_head_icol(sa_icol) == 0) this%record_nfeatures(sa_icol) = 0
    end do
  end subroutine allocate_record_targets

  !> @brief Allocate idt's permanent array with init_value, unless already
  !! allocated.
  !<
  subroutine allocate_permanent_array(this, idt, nfeatures, init_value, &
                                      varname)
    use MemoryManagerModule, only: mem_allocate
    class(KeystringLoadType), intent(inout) :: this
    type(InputParamDefinitionType), intent(in) :: idt
    integer(I4B), intent(in) :: nfeatures
    real(DP), intent(in) :: init_value
    character(len=*), intent(in), optional :: varname !< overrides idt%mf6varname (e.g. when it collides with an existing raw column)
    real(DP), dimension(:), pointer, contiguous :: featarr => null()
    integer(I4B) :: isize
    character(len=LENVARNAME) :: name

    name = trim(idt%mf6varname)
    if (present(varname)) name = trim(varname)
    call get_isize(name, this%mf6_input%mempath, isize)
    if (isize > 0) return ! already allocated (shouldn't happen; df() runs once)
    call mem_allocate(featarr, nfeatures, name, this%mf6_input%mempath)
    featarr = init_value
  end subroutine allocate_permanent_array

  !> @brief Allocate permanent, feature-indexed storage for every in-scope
  !! PERIOD setting, keyed by the field's public tag (e.g. RATE).
  !<
  subroutine allocate_settings(this)
    use SimModule, only: count_errors, store_error_filename
    class(KeystringLoadType), intent(inout) :: this
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: icol, nfeatures, member_nfeatures

    nfeatures = this%resolve_nfeatures()
    if (nfeatures < 1) return

    do icol = this%nleading + 1, this%nparam
      idt => this%resolve_in_scope_setting(icol)
      if (.not. associated(idt)) cycle
      member_nfeatures = &
        this%resolve_member_nfeatures(idt%tagname, idt%shape, nfeatures)
      if (member_nfeatures < 1) cycle
      call this%allocate_permanent_array(idt, member_nfeatures, DZERO)
    end do
    if (count_errors() > 0) then
      call store_error_filename(this%input_name)
    end if
  end subroutine allocate_settings

  !> @brief Allocate permanent, node-indexed storage for every PERIOD
  !! setting in scope, for CELLID-addressed packages (TVK/TVS).
  !! Sized by ctx%nodes; DNODATA marks a node as never set, telling the
  !! package which to copy at its own reduced node number.
  !<
  subroutine allocate_node_settings(this)
    class(KeystringLoadType), intent(inout) :: this
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: icol, nfeatures

    if (.not. associated(this%ctx%nodes)) return
    nfeatures = this%ctx%nodes
    if (nfeatures < 1) return

    do icol = this%nleading + 1, this%nparam
      idt => this%resolve_in_scope_setting(icol)
      if (.not. associated(idt)) cycle
      call this%allocate_permanent_array(idt, nfeatures, DNODATA)
    end do
  end subroutine allocate_node_settings

  !> @brief Validate ifno against nfeatures, storing an error keyed by
  !! ifno_tagname (the leading column's public tag) if out of range.
  !<
  function valid_ifno(this, ifno, nfeatures, ifno_tagname, row) result(valid)
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
  end function valid_ifno

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
      if (.not. this%valid_ifno(ifno, nfeatures, ifno_tagname, i)) cycle
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
      if (.not. this%valid_ifno(ifno, nfeatures, ifno_tagname, i)) cycle
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

  !> @brief Apply PERIOD settings for ctx%record_dependency's targets,
  !! via ts_update_indexed -- same mechanism as BEDK/MANNING, but each
  !! row's target index comes from a named sibling field through the
  !! cached offset table rather than the leading id column alone.
  !<
  subroutine apply_record_targets(this)
    class(KeystringLoadType), intent(inout) :: this
    type(InputParamDefinitionType), pointer :: idt, head_idt
    real(DP), dimension(:), pointer, contiguous :: featarr => null()
    integer(I4B), pointer :: nbound => null()
    integer(I4B), dimension(:), allocatable :: row_addr
    integer(I4B) :: icol, sa_icol, padj, isize

    if (.not. allocated(this%record_nfeatures)) return

    call get_isize('NBOUND', this%mf6_input%mempath, isize)
    if (isize < 1) return
    call mem_setptr(nbound, 'NBOUND', this%mf6_input%mempath)
    if (nbound <= 0) return

    padj = 0
    if (this%ctx%has_setting_dispatch) padj = 1

    do icol = this%nleading + 1, this%nparam
      sa_icol = icol + padj
      if (sa_icol > size(this%record_nfeatures)) cycle
      if (this%record_nfeatures(sa_icol) < 1) cycle

      idt => this%structarray%struct_vectors(sa_icol)%idt
      head_idt => &
        this%structarray%struct_vectors(this%record_head_icol(sa_icol))%idt
      call mem_setptr(featarr, trim(idt%mf6varname)//'_RESOLVED', &
                      this%mf6_input%mempath)
      row_addr = this%resolve_member_row_addr(this%record_id_varname, &
                                              trim(head_idt%mf6varname), &
                                              this%record_index_icol(sa_icol), &
                                              this%record_offsets(sa_icol)%vals, &
                                              nbound)
      call this%structarray%ts_update_indexed( &
        sa_icol, this%tsmanager, this%mf6_input%subcomponent_name, &
        this%ctx%iprpak, nbound, row_addr, trim(idt%tagname), featarr)
    end do
  end subroutine apply_record_targets

  !> @brief Apply PERIOD settings (non-AUX) to their permanent arrays,
  !! so an unrepeated setting in a later period keeps its prior value.
  !! Feature-addressed (IFNO/BNDNO) unless node_addressed (TVK/TVS).
  !<
  subroutine apply_settings(this, node_addressed)
    use DefinitionSelectModule, only: get_param_definition_type
    use GeomUtilModule, only: get_node
    use SimModule, only: count_errors, store_error_filename
    class(KeystringLoadType), intent(inout) :: this
    logical(LGP), intent(in) :: node_addressed
    integer(I4B), pointer :: nbound => null()
    integer(I4B), dimension(:), pointer, contiguous :: period_ifno => null()
    integer(I4B), dimension(:, :), pointer, contiguous :: cellid => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      period_setting => null()
    type(InputParamDefinitionType), pointer :: idt
    real(DP), dimension(:), pointer, contiguous :: featarr => null()
    integer(I4B) :: i, icol, ifno, nodeu, isize, nfeatures, &
                    member_nfeatures, ndim, padj
    integer(I4B), dimension(:), allocatable :: row_addr
    character(len=LINELENGTH) :: setting
    character(len=LENVARNAME) :: ifno_tagname

    call get_isize('NBOUND', this%mf6_input%mempath, isize)
    if (isize < 1) return
    call mem_setptr(nbound, 'NBOUND', this%mf6_input%mempath)
    if (nbound <= 0) return

    if (node_addressed) then
      if (.not. associated(this%ctx%nodes)) return
      nfeatures = this%ctx%nodes
      if (nfeatures < 1) return
      if (.not. associated(this%ctx%mshape)) return
      ndim = size(this%ctx%mshape)
      call mem_setptr(cellid, 'CELLID', this%mf6_input%mempath)
    else
      nfeatures = this%resolve_nfeatures()
      if (nfeatures < 1) return
      ! the leading column is the permanent feature address (IFNO for
      ! advanced packages, BNDNO for SPC), resolved via its own idt
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       'PERIOD', this%param_names(1), &
                                       this%input_name)
      ifno_tagname = trim(idt%tagname)
      call mem_setptr(period_ifno, trim(idt%mf6varname), &
                      this%mf6_input%mempath)
    end if

    ! member columns are offset past the SETTING column in the struct
    ! array (see create_structarray) when a setting dispatch is present
    padj = 0
    if (this%ctx%has_setting_dispatch) padj = 1

    call mem_setptr(period_setting, 'SETTING', this%mf6_input%mempath)

    do icol = this%nleading + 1, this%nparam
      idt => this%resolve_in_scope_setting(icol)
      if (.not. associated(idt)) cycle
      if (.not. node_addressed) then
        member_nfeatures = &
          this%resolve_member_nfeatures(idt%tagname, idt%shape, nfeatures)
        if (member_nfeatures < 1) cycle
      end if
      call mem_setptr(featarr, trim(idt%mf6varname), this%mf6_input%mempath)

      ! resolve each row's address (0 if this row doesn't set idt), then
      ! apply the whole column via ts_update_indexed
      allocate (row_addr(nbound))
      do i = 1, nbound
        row_addr(i) = 0
        setting = period_setting(i)
        if (trim(setting) /= trim(idt%mf6varname)) cycle
        if (node_addressed) then
          if (ndim == 1) then
            nodeu = cellid(1, i)
          else if (ndim == 2) then
            nodeu = get_node(cellid(1, i), 1, cellid(2, i), &
                             this%ctx%mshape(1), 1, this%ctx%mshape(2))
          else
            nodeu = get_node(cellid(1, i), cellid(2, i), cellid(3, i), &
                             this%ctx%mshape(1), this%ctx%mshape(2), &
                             this%ctx%mshape(3))
          end if
          if (nodeu < 1 .or. nodeu > nfeatures) cycle
          row_addr(i) = nodeu
        else
          ifno = period_ifno(i)
          if (.not. this%valid_ifno(ifno, member_nfeatures, ifno_tagname, &
                                    i)) cycle
          row_addr(i) = ifno
        end if
      end do
      call this%structarray%ts_update_indexed( &
        icol + padj, this%tsmanager, this%mf6_input%subcomponent_name, &
        this%ctx%iprpak, nbound, row_addr, trim(idt%tagname), featarr)
      deallocate (row_addr)
    end do
    if (count_errors() > 0) then
      call store_error_filename(this%input_name)
    end if
  end subroutine apply_settings

  !> @brief Resolve the permanent array's feature count.
  !!
  !! Advanced packages: PACKAGEDATA's row count (PACKAGEDATA_IFNO).
  !! DIMENSIONS-scoped packages: ctx%maxbound divided by the keystring
  !! member count (member_nsubs), reversing scale_keystring_maxbound.
  !<
  function resolve_nfeatures(this) result(nfeatures)
    class(KeystringLoadType), intent(inout) :: this
    integer(I4B) :: nfeatures
    integer(I4B) :: isize, nmembers

    call get_isize('PACKAGEDATA_IFNO', this%mf6_input%mempath, isize)
    if (isize > 0) then
      nfeatures = isize
      return
    end if

    nfeatures = 0
    nmembers = 0
    if (allocated(this%ctx%member_nsubs)) nmembers = size(this%ctx%member_nsubs)
    if (nmembers > 0 .and. associated(this%ctx%maxbound)) then
      if (this%ctx%maxbound > 0) nfeatures = this%ctx%maxbound / nmembers
    end if
  end function resolve_nfeatures

  !> @brief Feature count from a named dimension, falling back to
  !! default_nfeatures if unset. Populated-then-released is an error.
  !<
  function resolve_member_nfeatures(this, member_tagname, dimname, &
                                    default_nfeatures) result(nfeatures)
    use SimModule, only: store_error
    class(KeystringLoadType), intent(inout) :: this
    character(len=*), intent(in) :: member_tagname
    character(len=*), intent(in) :: dimname
    integer(I4B), intent(in) :: default_nfeatures
    integer(I4B) :: nfeatures
    integer(I4B), pointer :: shape_val => null()
    integer(I4B), dimension(:), pointer, contiguous :: shape_arr => null()
    integer(I4B) :: isize
    character(len=LINELENGTH) :: errmsg

    nfeatures = default_nfeatures
    if (dimname == '') return
    call get_isize(trim(dimname), this%mf6_input%mempath, isize)
    if (isize < 0) then
      ! -- never populated -- e.g. an omitted, legitimately zero dimension
      nfeatures = 0
      return
    else if (isize == 0) then
      ! -- was populated, then released -- a declared dependency that's
      !    no longer available is an error, not a silent fallback
      write (errmsg, '(a,1x,a,1x,a)') &
        'member', trim(member_tagname)//': DIMENSION', &
        trim(dimname)//' is not defined.'
      call store_error(errmsg)
      nfeatures = 0
      return
    else if (.not. this%shape_param_is_array(trim(dimname))) then
      call mem_setptr(shape_val, trim(dimname), this%mf6_input%mempath)
      nfeatures = shape_val
    else
      ! -- an array-valued dimension (e.g. one count per reach) sums to
      !    the total feature count
      call mem_setptr(shape_arr, trim(dimname), this%mf6_input%mempath)
      nfeatures = sum(shape_arr)
    end if
  end function resolve_member_nfeatures

  !> @brief Is shape_varname (matched by mf6varname) a per-feature array
  !! (PACKAGEDATA) rather than a package-wide scalar (DIMENSIONS)?
  !<
  function shape_param_is_array(this, shape_varname) result(is_array)
    class(KeystringLoadType), intent(inout) :: this
    character(len=*), intent(in) :: shape_varname
    logical(LGP) :: is_array
    integer(I4B) :: i

    is_array = .false.
    do i = 1, size(this%mf6_input%param_dfns)
      if (this%mf6_input%param_dfns(i)%component_type == &
          this%mf6_input%component_type .and. &
          this%mf6_input%param_dfns(i)%subcomponent_type == &
          this%mf6_input%subcomponent_type .and. &
          trim(this%mf6_input%param_dfns(i)%mf6varname) == &
          trim(shape_varname)) then
        is_array = &
          (trim(this%mf6_input%param_dfns(i)%blockname) == 'PACKAGEDATA')
        exit
      end if
    end do
  end function shape_param_is_array

  !> @brief Per-feature cumulative offset table from an array-valued
  !! dimension, permuted into feature-index order first. Indexed by
  !! feature (not the dimension's own summed total).
  !<
  function resolve_member_offsets(this, dimname) result(offsets)
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
  end function resolve_member_offsets

  !> @brief Row -> resolved feature index for a member whose local
  !! position comes from a sibling index field, via the offset table.
  !<
  function resolve_member_row_addr(this, id_mf6varname, head_mf6varname, &
                                   index_icol, offsets, nbound) &
    result(row_addr)
    class(KeystringLoadType), intent(inout) :: this
    character(len=*), intent(in) :: id_mf6varname
    character(len=*), intent(in) :: head_mf6varname
    integer(I4B), intent(in) :: index_icol !< SA column holding the sibling local index (e.g. IDV)
    integer(I4B), dimension(:), intent(in) :: offsets
    integer(I4B), intent(in) :: nbound
    integer(I4B), dimension(:), allocatable :: row_addr
    integer(I4B), dimension(:), pointer, contiguous :: period_ifno => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      period_setting => null()
    integer(I4B) :: i, ifno, idx_local
    character(len=LINELENGTH) :: setting

    call mem_setptr(period_ifno, id_mf6varname, this%mf6_input%mempath)
    call mem_setptr(period_setting, 'SETTING', this%mf6_input%mempath)
    allocate (row_addr(nbound))
    do i = 1, nbound
      row_addr(i) = 0
      setting = period_setting(i)
      if (trim(setting) /= trim(head_mf6varname)) cycle
      ifno = period_ifno(i)
      if (ifno < 1 .or. ifno > size(offsets)) cycle
      idx_local = this%structarray%struct_vectors(index_icol)%int1d(i)
      if (idx_local < 1) cycle
      row_addr(i) = offsets(ifno) + idx_local - 1
    end do
  end function resolve_member_row_addr

  !> @brief Return idt for param_names(icol) if it's an in-scope PERIOD
  !! setting (DOUBLE type with TIME_SERIES TRUE), else disassociated.
  !!
  !! Never matches a RECORD-compound sub-member.
  !<
  function resolve_in_scope_setting(this, icol) result(idt)
    use DefinitionSelectModule, only: get_param_definition_type
    class(KeystringLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: icol
    type(InputParamDefinitionType), pointer :: idt

    idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                     this%mf6_input%component_type, &
                                     this%mf6_input%subcomponent_type, &
                                     'PERIOD', this%param_names(icol), &
                                     this%input_name)
    if (idt%datatype /= 'DOUBLE' .or. .not. idt%timeseries) idt => null()
    ! a RECORD sub-member isn't feature-indexed by its own mf6varname --
    ! SETTING only ever names its record's head
    if (associated(idt)) then
      if (this%ctx%member_is_follower(icol - this%nleading)) idt => null()
    end if
  end function resolve_in_scope_setting

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

    ! use pre-allocated managed memory (maxbound = features * nmembers);
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

    ! create member columns
    do icol = this%nleading + 1, this%nparam
      sa_icol = icol + padj
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       'PERIOD', &
                                       this%param_names(icol), this%input_name)
      ! nsub from context: 0 = direct dispatch, N = KEYWORD compound with N sub-members
      nsub = this%ctx%member_nsubs(icol - this%nleading)
      if (nsub > 0) then
        ! metadata vector: no data allocated; isubmember points to next SA col
        call this%structarray%mem_create_metadata_vector(sa_icol, idt, &
                                                         sa_icol + 1, nsub)
      else if (trim(idt%datatype) == 'STRING') then
        ! string value columns (e.g. STATUS) stored at LENVARNAME
        call this%structarray%mem_create_vector(sa_icol, idt, &
                                                charlen=LENVARNAME)
      else if (idt%datatype == 'DOUBLE' .and. idt%timeseries .and. &
               .not. this%ctx%member_is_follower(icol - this%nleading)) then
        ! TS-capable setting: raw array uses IDM's derived input name;
        ! RECORD sub-members excluded (not feature-indexed by mf6varname)
        call this%structarray%mem_create_vector(sa_icol, idt, &
                                                varname=idm_input_varname(idt))
      else
        call this%structarray%mem_create_vector(sa_icol, idt)
      end if
    end do
  end subroutine create_structarray

end module Mf6FileKeystringModule
