!> @brief This module contains the ListLoadModule
!!
!! This module contains the routines for reading period block
!! list based input.
!!
!<
module ListLoadModule

  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LINELENGTH, DNODATA, DZERO
  use InputDefinitionModule, only: InputParamDefinitionType
  use MemoryManagerModule, only: mem_setptr, mem_allocate, get_isize
  use CharacterStringModule, only: CharacterStringType
  use ModflowInputModule, only: ModflowInputType
  use TimeSeriesManagerModule, only: TimeSeriesManagerType, tsmanager_cr, &
                                     read_value_or_time_series_adv
  use StructArrayModule, only: StructArrayType, constructStructArray, &
                               destructStructArray, idm_input_varname
  use AsciiInputLoadTypeModule, only: AsciiDynamicPkgLoadBaseType
  use LoadContextModule, only: LoadContextType
  use LoadMf6FileModule, only: LoadMf6FileType

  implicit none
  private
  public :: ListLoadType

  !> @brief list input loader for dynamic packages.
  !!
  !! Create and update input context for list based period blocks.
  !!
  !<
  type, extends(AsciiDynamicPkgLoadBaseType) :: ListLoadType
    type(TimeSeriesManagerType), pointer :: tsmanager => null()
    type(StructArrayType), pointer :: structarray => null()
    type(LoadContextType) :: ctx
    type(LoadMf6FileType) :: static_loader ! persistent static loader
    logical(LGP) :: ts_active !< .true. if TS files are loaded
  contains
    procedure :: ainit
    procedure :: df
    procedure :: ts_advance
    procedure :: reset
    procedure :: rp
    procedure :: destroy
    procedure :: create_structarray
    procedure :: apply_persistent_settings
  end type ListLoadType

contains

  subroutine ainit(this, mf6_input, component_name, component_input_name, &
                   input_name, iperblock, parser, iout)
    use InputOutputModule, only: getunit
    use MemoryManagerModule, only: get_isize
    use CharacterStringModule, only: CharacterStringType
    use BlockParserModule, only: BlockParserType
    class(ListLoadType), intent(inout) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), intent(in) :: component_name
    character(len=*), intent(in) :: component_input_name
    character(len=*), intent(in) :: input_name
    integer(I4B), intent(in) :: iperblock
    type(BlockParserType), pointer, intent(inout) :: parser
    integer(I4B), intent(in) :: iout
    type(CharacterStringType), dimension(:), pointer, contiguous :: ts_fnames
    character(len=LINELENGTH) :: fname
    integer(I4B) :: ts6_size, n

    ! init loader
    call this%DynamicPkgLoadType%init(mf6_input, component_name, &
                                      component_input_name, input_name, &
                                      iperblock, iout)
    ! initialize scalars
    this%ts_active = .false.

    ! create tsmanager
    allocate (this%tsmanager)
    call tsmanager_cr(this%tsmanager, iout)

    ! load static input (TS6_FILENAME tag sets static_loader%ts_active)
    call this%static_loader%load(parser, mf6_input, this%nc_vars, &
                                 this%input_name, iout)

    ! if TS files were declared, add them to our tsmanager now
    if (this%static_loader%ts_active) then
      this%ts_active = .true.
      call get_isize('TS6_FILENAME', mf6_input%mempath, ts6_size)
      if (ts6_size > 0) then
        call mem_setptr(ts_fnames, 'TS6_FILENAME', mf6_input%mempath)
        do n = 1, size(ts_fnames)
          fname = ts_fnames(n)
          call this%tsmanager%add_tsfile(fname, getunit())
        end do
      end if
    end if

    ! initialize package input context
    call this%ctx%init(mf6_input)

    ! set in-scope param names directly from context
    this%param_names = this%ctx%params
    this%nparam = size(this%ctx%params)
    call this%ctx%check_developmode(this%input_name)

    ! construct and set up the struct array object
    call this%create_structarray()

    ! finalize input context setup
    call this%ctx%allocate_arrays()
  end subroutine ainit

  subroutine df(this)
    use StructArrayModule, only: StructArrayType
    class(ListLoadType), intent(inout) :: this
    type(StructArrayType), pointer :: sa
    integer(I4B) :: n
    ! define tsmanager (TDIS is now available)
    call this%tsmanager%tsmanager_df()
    ! link static TS strlocs; preserve for re-registration after reset()
    do n = 1, this%static_loader%ts_sa_count()
      sa => this%static_loader%get_ts_sa(n)
      if (associated(sa)) then
        call sa%ts_update(this%tsmanager, &
                          this%mf6_input%subcomponent_name, &
                          this%ctx%iprpak, this%input_name, &
                          this%ctx%auxname_cst, &
                          clear_strlocs=.false.)
      end if
    end do
  end subroutine df

  subroutine ts_advance(this)
    class(ListLoadType), intent(inout) :: this
    ! advance timeseries
    call this%tsmanager%ad()
  end subroutine ts_advance

  subroutine reset(this)
    use StructArrayModule, only: StructArrayType
    class(ListLoadType), intent(inout) :: this
    type(StructArrayType), pointer :: sa
    integer(I4B) :: n
    ! clear TS links, unless persistent (unmentioned rows keep their link)
    if (.not. this%ctx%is_advanced) then
      call this%tsmanager%reset(this%mf6_input%subcomponent_name)
    end if
    ! re-register static TS links (strlocs preserved in df)
    if (this%ts_active) then
      do n = 1, this%static_loader%ts_sa_count()
        sa => this%static_loader%get_ts_sa(n)
        if (associated(sa)) then
          call sa%ts_update(this%tsmanager, &
                            this%mf6_input%subcomponent_name, &
                            this%ctx%iprpak, this%input_name, &
                            this%ctx%auxname_cst, &
                            clear_strlocs=.false.)
        end if
      end do
    end if
  end subroutine reset

  subroutine rp(this, parser)
    use BlockParserModule, only: BlockParserType
    use LoadMf6FileModule, only: read_control_record
    use StructVectorModule, only: StructVectorType
    use IdmLoggerModule, only: idm_log_header, idm_log_close
    class(ListLoadType), intent(inout) :: this
    type(BlockParserType), pointer, intent(inout) :: parser
    integer(I4B) :: ibinary
    integer(I4B) :: oc_inunit

    call this%reset()
    ibinary = read_control_record(parser, oc_inunit, this%iout)

    ! log lst file header
    call idm_log_header(this%mf6_input%component_name, &
                        this%mf6_input%subcomponent_name, this%iout)

    if (ibinary == 1) then
      this%ctx%nbound = &
        this%structarray%read_from_binary(oc_inunit, this%iout)
      call parser%terminateblock()
      close (oc_inunit)
    else
      this%ctx%nbound = &
        this%structarray%read_from_parser(parser, this%ts_active, this%iout, &
                                          this%input_name)
    end if

    ! must run before ts_update below, so AUX's ts_strlocs are claimed
    ! by ts_update_adv first, not consumed into the transient array
    call this%apply_persistent_settings()

    ! update ts links for all other columns
    if (this%ts_active) then
      call this%structarray%ts_update(this%tsmanager, &
                                      this%mf6_input%subcomponent_name, &
                                      this%ctx%iprpak, this%input_name, &
                                      this%ctx%auxname_cst)
    end if

    ! close logging statement
    call idm_log_close(this%mf6_input%component_name, &
                       this%mf6_input%subcomponent_name, this%iout)
  end subroutine rp

  subroutine destroy(this)
    class(ListLoadType), intent(inout) :: this
    !
    ! clean up saved static structarrays
    call this%static_loader%cleanup()
    !
    ! deallocate tsmanager
    call this%tsmanager%da()
    deallocate (this%tsmanager)
    nullify (this%tsmanager)
    !
    ! deallocate StructArray
    call destructStructArray(this%structarray)
    call this%ctx%destroy()
  end subroutine destroy

  subroutine create_structarray(this)
    use InputDefinitionModule, only: InputParamDefinitionType
    use DefinitionSelectModule, only: get_param_definition_type
    class(ListLoadType), intent(inout) :: this
    type(InputParamDefinitionType), pointer :: idt
    real(DP), dimension(:), pointer, contiguous :: featarr
    real(DP), dimension(:, :), pointer, contiguous :: featarr2d
    integer(I4B), pointer :: naux
    integer(I4B) :: icol, isize

    ! construct and set up the struct array object
    this%structarray => constructStructArray(this%mf6_input, this%nparam, &
                                             this%ctx%maxbound, 0, &
                                             this%mf6_input%mempath, &
                                             this%mf6_input%component_mempath)
    ! set up struct array
    do icol = 1, this%nparam
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       'PERIOD', &
                                       this%param_names(icol), this%input_name)
      ! persistent array allocated under mf6varname; raw column below
      ! uses IDM's derived input name instead
      if (this%ctx%is_advanced .and. idt%datatype == 'DOUBLE' .and. &
          idt%timeseries) then
        call get_isize(trim(idt%mf6varname), this%mf6_input%mempath, isize)
        if (isize <= 0) then
          call mem_allocate(featarr, this%ctx%maxbound, trim(idt%mf6varname), &
                            this%mf6_input%mempath)
          featarr = DZERO
        end if
      else if (this%ctx%is_advanced .and. idt%datatype == 'DOUBLE1D' .and. &
               idt%timeseries) then
        ! AUX is NAUX-gated, so its bare tag can't be reclaimed the
        ! same way; it keeps its own synthetic tag
        call get_isize(trim(idt%tagname)//'VAR', this%mf6_input%mempath, &
                       isize)
        if (isize <= 0) then
          call mem_setptr(naux, trim(idt%shape), this%mf6_input%mempath)
          call mem_allocate(featarr2d, naux, this%ctx%maxbound, &
                            trim(idt%tagname)//'VAR', this%mf6_input%mempath)
          featarr2d = DZERO
        end if
      end if
      ! allocate variable in memory manager
      if (this%ctx%is_advanced .and. idt%datatype == 'DOUBLE' .and. &
          idt%timeseries) then
        call this%structarray%mem_create_vector(icol, idt, &
                                                varname=idm_input_varname(idt))
      else
        call this%structarray%mem_create_vector(icol, idt)
      end if
    end do
  end subroutine create_structarray

  !> @brief Resolve this period's rows for an advanced package's TS-capable
  !! fields into their permanent, feature-indexed backing arrays.
  !!
  !! Every row's IFNO is validated against maxbound before use.
  !<
  subroutine apply_persistent_settings(this)
    use DefinitionSelectModule, only: get_param_definition_type
    use SimModule, only: store_error, count_errors, store_error_filename
    use SimVariablesModule, only: errmsg
    class(ListLoadType), intent(inout) :: this
    type(InputParamDefinitionType), pointer :: idt
    integer(I4B), dimension(:), pointer, contiguous :: ifno
    integer(I4B), dimension(:), allocatable :: row_ifno
    real(DP), dimension(:), pointer, contiguous :: featarr
    real(DP), dimension(:, :), pointer, contiguous :: featarr2d
    integer(I4B), pointer :: naux
    integer(I4B) :: icol, n, i, j, nfeatures

    if (.not. this%ctx%is_advanced) return

    ! leading column's own name/tag (e.g. IFNO), resolved from the
    ! recarray definition rather than hardcoded
    idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                     this%mf6_input%component_type, &
                                     this%mf6_input%subcomponent_type, &
                                     'PERIOD', this%param_names(1), &
                                     this%input_name)
    call mem_setptr(ifno, trim(idt%mf6varname), this%mf6_input%mempath)

    nfeatures = this%ctx%maxbound
    allocate (row_ifno(this%ctx%nbound))
    do n = 1, this%ctx%nbound
      if (ifno(n) >= 1 .and. ifno(n) <= nfeatures) then
        row_ifno(n) = ifno(n)
      else
        write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a,1x,i0,a)') &
          trim(idt%tagname), ifno(n), 'on row', n, &
          'must be greater than 0 and less than or equal to', nfeatures, '.'
        call store_error(errmsg)
        row_ifno(n) = 0
      end if
    end do

    if (count_errors() > 0) then
      call store_error_filename(this%input_name)
    end if

    do icol = 1, this%nparam
      idt => get_param_definition_type(this%mf6_input%param_dfns, &
                                       this%mf6_input%component_type, &
                                       this%mf6_input%subcomponent_type, &
                                       'PERIOD', &
                                       this%param_names(icol), this%input_name)
      if (idt%datatype == 'DOUBLE' .and. idt%timeseries) then
        call mem_setptr(featarr, trim(idt%mf6varname), this%mf6_input%mempath)
        if (this%ts_active) then
          call this%structarray%ts_update_indexed( &
            icol, this%tsmanager, this%mf6_input%subcomponent_name, &
            this%ctx%iprpak, this%ctx%nbound, row_ifno, &
            varname=trim(idt%tagname), featarr=featarr)
        else
          do n = 1, this%ctx%nbound
            i = row_ifno(n)
            if (i < 1) cycle
            if (this%structarray%struct_vectors(icol)%dbl1d(n) == DNODATA) cycle
            featarr(i) = this%structarray%struct_vectors(icol)%dbl1d(n)
          end do
        end if
      else if (idt%datatype == 'DOUBLE1D' .and. idt%timeseries) then
        call mem_setptr(featarr2d, trim(idt%tagname)//'VAR', &
                        this%mf6_input%mempath)
        if (this%ts_active) then
          call this%structarray%ts_update_adv( &
            icol, this%tsmanager, this%mf6_input%subcomponent_name, &
            this%ctx%iprpak, this%ctx%nbound, row_ifno, &
            auxname_cst=this%ctx%auxname_cst, featarr2d=featarr2d)
        else
          call mem_setptr(naux, trim(idt%shape), this%mf6_input%mempath)
          do n = 1, this%ctx%nbound
            i = row_ifno(n)
            if (i < 1) cycle
            do j = 1, naux
              if (this%structarray%struct_vectors(icol)%dbl2d(j, n) == DNODATA) &
                cycle
              featarr2d(j, i) = &
                this%structarray%struct_vectors(icol)%dbl2d(j, n)
            end do
          end do
        end if
      end if
    end do
    deallocate (row_ifno)
  end subroutine apply_persistent_settings

end module ListLoadModule
