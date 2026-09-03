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
                                     read_value_or_time_series_adv
  use StructArrayModule, only: StructArrayType, constructStructArray, &
                               destructStructArray
  use AsciiInputLoadTypeModule, only: AsciiDynamicPkgLoadBaseType
  use LoadContextModule, only: LoadContextType, is_advanced_package
  use LoadMf6FileModule, only: LoadMf6FileType
  use BlockParserModule, only: BlockParserType

  implicit none
  private
  public :: KeystringLoadType

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
  contains
    procedure :: ainit
    procedure :: df
    procedure :: ts_advance
    procedure :: rp
    procedure :: reset
    procedure :: destroy
    procedure :: create_structarray
    procedure :: apply_period_auxiliary
    procedure :: allocate_period_settings
    procedure :: apply_period_settings
    procedure :: resolve_nfeatures
    procedure :: resolve_member_nfeatures
    procedure :: resolve_in_scope_setting
    procedure :: allocate_permanent_array
    procedure :: apply_setting_value
    procedure :: allocate_period_node_settings
    procedure :: apply_period_node_settings
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

    ! find a DIMENSIONS block parameter name to alias as maxbound for
    ! LoadContext; skipped for advanced packages, whose feature count
    ! comes from PACKAGEDATA's own row count instead. Uses the first
    ! DIMENSIONS parameter found; other parameters (e.g. LAK's
    ! NOUTLETS/NTABLES) are unrelated to feature count and must not be
    ! summed in.
    has_named_bound = .false.
    if (.not. is_advanced_package(mf6_input)) then
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
    ! advanced packages: address AUX TS links by feature number (PACKAGEDATA
    ! row position does not necessarily equal IFNO) so a later PERIOD
    ! override can find and replace this registration
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
    ! feature-indexed packages (advanced packages and, separately,
    ! DIMENSIONS-scoped packages like SPC): allocate permanent,
    ! feature-indexed storage for every PERIOD setting in scope for
    ! sticky, TS-continuable persistence
    if (this%ctx%is_feature_indexed) call this%allocate_period_settings()
    ! CELLID-addressed packages (TVK/TVS): same persistence goal, but
    ! node-indexed and package-resolved -- see allocate_period_node_settings
    if (this%ctx%is_cellid_scoped) call this%allocate_period_node_settings()
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

    if (this%ctx%is_advanced) call this%apply_period_auxiliary()
    if (this%ctx%is_feature_indexed) call this%apply_period_settings()
    if (this%ctx%is_cellid_scoped) call this%apply_period_node_settings()

    if (this%ts_active) then
      call this%structarray%ts_update(this%tsmanager, &
                                      this%mf6_input%subcomponent_name, &
                                      this%ctx%iprpak, this%input_name)
    end if

    call idm_log_close(this%mf6_input%component_name, &
                       this%mf6_input%subcomponent_name, this%iout)
  end subroutine rp

  !> @brief Apply PERIOD AUXILIARY settings to the permanent AUX array
  !!
  !! AUXVAL is read as a string so literal values and TS names resolve
  !! here, directly against AUX(jj,ifno).
  !<
  subroutine apply_period_auxiliary(this)
    use DefinitionSelectModule, only: get_param_definition_type
    use SimModule, only: store_error
    use SimVariablesModule, only: errmsg
    class(KeystringLoadType), intent(inout) :: this
    integer(I4B), pointer :: nbound => null()
    integer(I4B), dimension(:), pointer, contiguous :: period_ifno => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      period_setting => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      period_auxname => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      period_auxval => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      auxnames => null()
    real(DP), dimension(:, :), pointer, contiguous :: aux => null()
    real(DP), pointer :: bndElem
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: i, ifno, jj, isize, naux, nfeatures
    character(len=LINELENGTH) :: token, setting, auxname, thisauxname
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

    call mem_setptr(period_ifno, 'IFNO', this%mf6_input%mempath)
    call mem_setptr(period_setting, 'SETTING', this%mf6_input%mempath)
    call mem_setptr(period_auxname, 'AUXNAME', this%mf6_input%mempath)
    call mem_setptr(period_auxval, 'AUXVAL', this%mf6_input%mempath)
    call mem_setptr(auxnames, 'AUXILIARY', this%mf6_input%mempath)
    call mem_setptr(aux, 'AUX', this%mf6_input%mempath)
    nfeatures = size(aux, 2)

    do i = 1, nbound
      setting = period_setting(i)
      if (trim(setting) /= 'AUXILIARY') cycle
      ifno = period_ifno(i)
      if (ifno < 1 .or. ifno > nfeatures) then
        write (errmsg, '(a,1x,2(a,1x),i0,a)') &
          trim(ifno_tagname), 'must be greater than 0 and', &
          'less than or equal to', nfeatures, '.'
        call store_error(errmsg)
        cycle
      end if
      auxname = period_auxname(i)
      do jj = 1, naux
        thisauxname = auxnames(jj)
        if (trim(auxname) /= trim(thisauxname)) cycle
        token = period_auxval(i)
        if (len_trim(token) == 0) exit
        bndElem => aux(jj, ifno)
        call read_value_or_time_series_adv(token, ifno, jj, bndElem, &
                                           this%mf6_input%subcomponent_name, &
                                           'AUX', this%tsmanager, &
                                           this%ctx%iprpak, trim(thisauxname))
        exit
      end do
    end do
  end subroutine apply_period_auxiliary

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

  !> @brief Feature count for one PERIOD setting member.
  !!
  !! Falls back to default_nfeatures unless the member's own SHAPE names
  !! a DIMENSIONS-block scalar, letting one keystring mix settings from
  !! more than one feature space. A scalar never populated (an omitted,
  !! legitimately zero dimension) is also a no-op; one populated and then
  !! released is an error, since a non-blank SHAPE is a declared
  !! dependency.
  !<
  function resolve_member_nfeatures(this, idt, default_nfeatures) &
    result(nfeatures)
    use SimModule, only: store_error
    class(KeystringLoadType), intent(inout) :: this
    type(InputParamDefinitionType), intent(in) :: idt
    integer(I4B), intent(in) :: default_nfeatures
    integer(I4B) :: nfeatures
    integer(I4B), pointer :: shape_val => null()
    integer(I4B) :: isize
    character(len=LINELENGTH) :: errmsg

    nfeatures = default_nfeatures
    if (idt%shape == '') return
    call get_isize(trim(idt%shape), this%mf6_input%mempath, isize)
    if (isize < 0) then
      ! -- never populated -- e.g. an omitted, legitimately zero dimension
      nfeatures = 0
      return
    else if (isize == 0) then
      ! -- was populated, then released -- a declared SHAPE dependency
      !    that's no longer available is an error, not a silent fallback
      write (errmsg, '(a,1x,a)') &
        'DIMENSION', trim(idt%shape)//' is not defined.'
      call store_error(errmsg)
      nfeatures = 0
      return
    end if
    call mem_setptr(shape_val, trim(idt%shape), this%mf6_input%mempath)
    nfeatures = shape_val
  end function resolve_member_nfeatures

  !> @brief Return idt for param_names(icol) if it's an in-scope PERIOD
  !! setting (STRING type with TIME_SERIES TRUE), else a disassociated
  !! pointer.
  !!
  !! Never matches a RECORD-compound sub-member -- SETTING records the
  !! dispatch keyword's own tagname, not a sub-member's.
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
    if (idt%datatype /= 'STRING' .or. .not. idt%timeseries) idt => null()
  end function resolve_in_scope_setting

  !> @brief Allocate idt's permanent array with init_value, unless already
  !! allocated.
  !<
  subroutine allocate_permanent_array(this, idt, nfeatures, init_value)
    use MemoryManagerModule, only: mem_allocate
    class(KeystringLoadType), intent(inout) :: this
    type(InputParamDefinitionType), intent(in) :: idt
    integer(I4B), intent(in) :: nfeatures
    real(DP), intent(in) :: init_value
    real(DP), dimension(:), pointer, contiguous :: featarr => null()
    integer(I4B) :: isize

    call get_isize(trim(idt%tagname), this%mf6_input%mempath, isize)
    if (isize > 0) return ! already allocated (shouldn't happen; df() runs once)
    call mem_allocate(featarr, nfeatures, trim(idt%tagname), &
                      this%mf6_input%mempath)
    featarr = init_value
  end subroutine allocate_permanent_array

  !> @brief Resolve token (row i of period_val) as a literal or TS name
  !! directly against featarr(address).
  !<
  subroutine apply_setting_value(this, idt, period_val, i, address, featarr)
    class(KeystringLoadType), intent(inout) :: this
    type(InputParamDefinitionType), intent(in) :: idt
    type(CharacterStringType), dimension(:), pointer, contiguous, &
      intent(in) :: period_val
    integer(I4B), intent(in) :: i
    integer(I4B), intent(in) :: address
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: featarr
    real(DP), pointer :: bndElem
    character(len=LINELENGTH) :: token

    token = period_val(i)
    if (len_trim(token) == 0) return
    bndElem => featarr(address)
    call read_value_or_time_series_adv(token, address, 0, bndElem, &
                                       this%mf6_input%subcomponent_name, &
                                       'BND', this%tsmanager, &
                                       this%ctx%iprpak, trim(idt%tagname))
  end subroutine apply_setting_value

  !> @brief Allocate permanent, feature-indexed storage for every in-scope
  !! PERIOD setting, keyed by the field's public tag (e.g. RATE).
  !<
  subroutine allocate_period_settings(this)
    class(KeystringLoadType), intent(inout) :: this
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B) :: icol, nfeatures, member_nfeatures

    nfeatures = this%resolve_nfeatures()
    if (nfeatures < 1) return

    do icol = this%nleading + 1, this%nparam
      idt => this%resolve_in_scope_setting(icol)
      if (.not. associated(idt)) cycle
      member_nfeatures = this%resolve_member_nfeatures(idt, nfeatures)
      if (member_nfeatures < 1) cycle
      call this%allocate_permanent_array(idt, member_nfeatures, DZERO)
    end do
  end subroutine allocate_period_settings

  !> @brief Apply PERIOD settings (non-AUX) to their permanent,
  !! feature-indexed arrays
  !!
  !! Resolves each token (literal or TS name) directly against the
  !! permanent, ifno-addressed array, so a later period that doesn't
  !! repeat the setting leaves the prior value untouched.
  !<
  subroutine apply_period_settings(this)
    use DefinitionSelectModule, only: get_param_definition_type
    use SimModule, only: store_error
    use SimVariablesModule, only: errmsg
    class(KeystringLoadType), intent(inout) :: this
    integer(I4B), pointer :: nbound => null()
    integer(I4B), dimension(:), pointer, contiguous :: period_ifno => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      period_setting => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      period_val => null()
    type(InputParamDefinitionType), pointer :: idt
    real(DP), dimension(:), pointer, contiguous :: featarr => null()
    integer(I4B) :: i, icol, ifno, isize, nfeatures, member_nfeatures
    character(len=LINELENGTH) :: setting
    character(len=LENVARNAME) :: ifno_tagname

    call get_isize('NBOUND', this%mf6_input%mempath, isize)
    if (isize < 1) return
    call mem_setptr(nbound, 'NBOUND', this%mf6_input%mempath)
    if (nbound <= 0) return

    nfeatures = this%resolve_nfeatures()
    if (nfeatures < 1) return

    ! the sole leading column is the permanent, stably-numbered feature
    ! address for every package this routine is gated for -- IFNO for
    ! advanced packages, BNDNO for SPC -- resolved via its own idt, since
    ! the memory-manager key (mf6varname) isn't always the doc-facing tag
    ! (e.g. MWE's leading PERIOD column is tagged MAWNO, MF6INTERNAL IFNO)
    idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                     this%mf6_input%component_type, &
                                     this%mf6_input%subcomponent_type, &
                                     'PERIOD', this%param_names(1), &
                                     this%input_name)
    ifno_tagname = trim(idt%tagname)
    call mem_setptr(period_ifno, trim(idt%mf6varname), this%mf6_input%mempath)
    call mem_setptr(period_setting, 'SETTING', this%mf6_input%mempath)

    do icol = this%nleading + 1, this%nparam
      idt => this%resolve_in_scope_setting(icol)
      if (.not. associated(idt)) cycle
      member_nfeatures = this%resolve_member_nfeatures(idt, nfeatures)
      if (member_nfeatures < 1) cycle
      call mem_setptr(featarr, trim(idt%tagname), this%mf6_input%mempath)
      call mem_setptr(period_val, trim(idt%mf6varname), this%mf6_input%mempath)

      do i = 1, nbound
        setting = period_setting(i)
        if (trim(setting) /= trim(idt%tagname)) cycle
        ifno = period_ifno(i)
        if (ifno < 1 .or. ifno > member_nfeatures) then
          write (errmsg, '(a,1x,2(a,1x),i0,a)') &
            trim(ifno_tagname), 'must be greater than 0 and', &
            'less than or equal to', member_nfeatures, '.'
          call store_error(errmsg)
          cycle
        end if
        call this%apply_setting_value(idt, period_val, i, ifno, featarr)
      end do
    end do
  end subroutine apply_period_settings

  !> @brief Allocate permanent, node-indexed storage for every PERIOD
  !! setting in scope, for CELLID-addressed packages (TVK/TVS)
  !!
  !! Sized by ctx%nodes (product(MSHAPE), the unreduced node count).
  !! DNODATA marks "never set" so the package knows which nodes to
  !! copy into its own target array (e.g. NPF's K11) at the reduced
  !! node number.
  !<
  subroutine allocate_period_node_settings(this)
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
  end subroutine allocate_period_node_settings

  !> @brief Apply PERIOD settings to their permanent, node-indexed arrays,
  !! for CELLID-addressed packages (TVK/TVS)
  !!
  !! Mirrors apply_period_settings, but addresses by the unreduced node
  !! number (nodeu) computed from CELLID via MSHAPE. The reduced-node
  !! lookup stays package-side, where the real dis object is available.
  !<
  subroutine apply_period_node_settings(this)
    use GeomUtilModule, only: get_node
    class(KeystringLoadType), intent(inout) :: this
    integer(I4B), pointer :: nbound => null()
    integer(I4B), dimension(:, :), pointer, contiguous :: cellid => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      period_setting => null()
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      period_val => null()
    type(InputParamDefinitionType), pointer :: idt
    real(DP), dimension(:), pointer, contiguous :: featarr => null()
    integer(I4B) :: i, icol, nodeu, isize, nfeatures, ndim
    character(len=LINELENGTH) :: setting

    call get_isize('NBOUND', this%mf6_input%mempath, isize)
    if (isize < 1) return
    call mem_setptr(nbound, 'NBOUND', this%mf6_input%mempath)
    if (nbound <= 0) return

    if (.not. associated(this%ctx%nodes)) return
    nfeatures = this%ctx%nodes
    if (nfeatures < 1) return
    if (.not. associated(this%ctx%mshape)) return
    ndim = size(this%ctx%mshape)

    call mem_setptr(cellid, 'CELLID', this%mf6_input%mempath)
    call mem_setptr(period_setting, 'SETTING', this%mf6_input%mempath)

    do icol = this%nleading + 1, this%nparam
      idt => this%resolve_in_scope_setting(icol)
      if (.not. associated(idt)) cycle
      call mem_setptr(featarr, trim(idt%tagname), this%mf6_input%mempath)
      call mem_setptr(period_val, trim(idt%mf6varname), this%mf6_input%mempath)

      do i = 1, nbound
        setting = period_setting(i)
        if (trim(setting) /= trim(idt%tagname)) cycle
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
        call this%apply_setting_value(idt, period_val, i, nodeu, featarr)
      end do
    end do
  end subroutine apply_period_node_settings

  subroutine reset(this)
    use StructArrayModule, only: StructArrayType
    use MemoryManagerModule, only: mem_setptr, get_isize
    use CharacterStringModule, only: CharacterStringType
    class(KeystringLoadType), intent(inout) :: this
    type(StructArrayType), pointer :: sa
    type(CharacterStringType), dimension(:), pointer, contiguous :: &
      auxnames => null()
    integer(I4B) :: n, naux
    ! every KEYSTRING subtype with a SETTING dispatch column: PERIOD
    ! settings persist across periods unless reissued, so TS links are
    ! never reset
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
      ! these packages normally resolve maxbound from PACKAGEDATA's or
      ! DIMENSIONS' own count; this remains as a fallback for a genuinely
      ! unresolvable count, e.g. an empty PACKAGEDATA block
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
      else
        call this%structarray%mem_create_vector(sa_icol, idt)
      end if
    end do
  end subroutine create_structarray

end module Mf6FileKeystringModule
