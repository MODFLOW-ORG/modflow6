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

  use KindModule, only: I4B, LGP
  use ConstantsModule, only: LINELENGTH, LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType
  use ModflowInputModule, only: ModflowInputType
  use CharacterStringModule, only: CharacterStringType
  use MemoryManagerModule, only: mem_setptr, get_isize
  use TimeSeriesManagerModule, only: TimeSeriesManagerType, tsmanager_cr
  use StructArrayModule, only: StructArrayType, constructStructArray, &
                               destructStructArray
  use AsciiInputLoadTypeModule, only: AsciiDynamicPkgLoadBaseType
  use LoadContextModule, only: LoadContextType
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
    procedure :: ad
    procedure :: rp
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
    use ArrayHandlersModule, only: expandarray
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
    character(len=LENVARNAME), allocatable :: named_bound(:)
    character(len=LINELENGTH), dimension(:), allocatable :: member_names
    integer(I4B) :: n, nmembers, isize

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

    ! collect DIMENSIONS block parameter names for LoadContext;
    ! absent for TVK/TVS — LoadContext falls back to nodes * nmembers
    do n = 1, size(mf6_input%param_dfns)
      if (mf6_input%param_dfns(n)%blockname == 'DIMENSIONS') then
        call expandarray(named_bound)
        named_bound(size(named_bound)) = trim(mf6_input%param_dfns(n)%mf6varname)
      end if
    end do

    ! init load context
    if (allocated(named_bound)) then
      call this%ctx%init(mf6_input, named_bound=named_bound)
    else
      call this%ctx%init(mf6_input)
    end if

    call this%ctx%tags(this%param_names, this%nparam, this%input_name)
    this%nleading = this%ctx%nleading

    ! append keystring member column names
    call this%ctx%keystring_member_names(member_names, nmembers)
    do n = 1, nmembers
      this%nparam = this%nparam + 1
      call expandarray(this%param_names)
      this%param_names(this%nparam) = trim(member_names(n))
    end do

    ! finalize context setup (allocates NBOUND, NODEULIST, etc.)
    call this%ctx%allocate_arrays()

    ! pre-allocate structarray; reused across all periods
    call this%create_structarray()
  end subroutine ainit

  subroutine df(this)
    use StructArrayModule, only: StructArrayType
    class(KeystringLoadType), intent(inout) :: this
    type(StructArrayType), pointer :: sa
    integer(I4B) :: n
    ! init tsmanager (TDIS now available)
    call this%tsmanager%tsmanager_df()
    ! link static TS strlocs; preserve for re-registration after reset()
    do n = 1, this%static_loader%ts_sa_count()
      sa => this%static_loader%get_ts_sa(n)
      if (associated(sa)) then
        call sa%ts_update(this%tsmanager, &
                          this%mf6_input%subcomponent_name, &
                          this%ctx%iprpak, this%input_name, &
                          clear_strlocs=.false.)
      end if
    end do
  end subroutine df

  subroutine ad(this)
    class(KeystringLoadType), intent(inout) :: this
    call this%tsmanager%ad()
  end subroutine ad

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

    if (this%ts_active) then
      call this%structarray%ts_update(this%tsmanager, &
                                      this%mf6_input%subcomponent_name, &
                                      this%ctx%iprpak, this%input_name)
    end if

    call idm_log_close(this%mf6_input%component_name, &
                       this%mf6_input%subcomponent_name, this%iout)
  end subroutine rp

  subroutine reset(this)
    use StructArrayModule, only: StructArrayType
    class(KeystringLoadType), intent(inout) :: this
    type(StructArrayType), pointer :: sa
    integer(I4B) :: n
    ! clear TS links
    call this%tsmanager%reset(this%mf6_input%subcomponent_name)
    ! re-register static TS links (strlocs preserved in df)
    if (this%ts_active) then
      do n = 1, this%static_loader%ts_sa_count()
        sa => this%static_loader%get_ts_sa(n)
        if (associated(sa)) then
          call sa%ts_update(this%tsmanager, &
                            this%mf6_input%subcomponent_name, &
                            this%ctx%iprpak, this%input_name, &
                            clear_strlocs=.false.)
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
    use DefinitionSelectModule, only: get_param_definition_type, &
                                      idt_parse_rectype
    use InputOutputModule, only: upcase
    class(KeystringLoadType), intent(inout) :: this
    type(InputParamDefinitionType), pointer :: idt, pidt
    character(len=LINELENGTH), allocatable :: rec_cols(:)
    character(len=LINELENGTH) :: kwname, first_col
    integer(I4B) :: icol, nrow_prealloc, jparam, nrec_col, nsub

    ! use pre-allocated managed memory (maxbound = features * nmembers);
    ! fall back to deferred shape (-1) if maxbound is unavailable
    if (associated(this%ctx%maxbound) .and. this%ctx%maxbound > 0) then
      nrow_prealloc = this%ctx%maxbound
    else
      nrow_prealloc = -1
    end if

    this%structarray => constructStructArray(this%mf6_input, this%nparam, &
                                             nrow_prealloc, 0, &
                                             this%mf6_input%mempath, &
                                             this%mf6_input%component_mempath)
    do icol = 1, this%nparam
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       'PERIOD', &
                                       this%param_names(icol), this%input_name)
      call this%structarray%mem_create_vector(icol, idt)

      ! For KEYWORD member columns, store the compound sub-member count
      ! from the RECORD definition so read_from_parser_keystring reads
      ! exactly the right number of sub-values.
      if (trim(idt%datatype) == 'KEYWORD' .and. icol > this%nleading) then
        kwname = trim(idt%tagname)
        call upcase(kwname)
        nsub = 0
        do jparam = 1, size(this%mf6_input%param_dfns)
          pidt => this%mf6_input%param_dfns(jparam)
          if (pidt%blockname /= 'PERIOD') cycle
          if (pidt%datatype(1:6) /= 'RECORD') cycle
          call idt_parse_rectype(pidt, rec_cols, nrec_col)
          if (nrec_col >= 1) then
            first_col = trim(rec_cols(1))
            call upcase(first_col)
            if (trim(first_col) == trim(kwname)) then
              nsub = nrec_col - 1
              if (allocated(rec_cols)) deallocate (rec_cols)
              exit
            end if
          end if
          if (allocated(rec_cols)) deallocate (rec_cols)
        end do
        this%structarray%struct_vectors(icol)%nsubmembers = nsub
      end if
    end do
  end subroutine create_structarray

end module Mf6FileKeystringModule
