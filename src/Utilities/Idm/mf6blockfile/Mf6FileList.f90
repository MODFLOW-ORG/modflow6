!> @brief This module contains the ListLoadModule
!!
!! This module contains the routines for reading period block
!! list based input.
!!
!<
module ListLoadModule

  use KindModule, only: I4B, LGP
  use ConstantsModule, only: LINELENGTH
  use InputDefinitionModule, only: InputParamDefinitionType
  use MemoryManagerModule, only: mem_setptr
  use CharacterStringModule, only: CharacterStringType
  use ModflowInputModule, only: ModflowInputType
  use TimeSeriesManagerModule, only: TimeSeriesManagerType, tsmanager_cr
  use StructArrayModule, only: StructArrayType, constructStructArray, &
                               destructStructArray
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

    ! store in scope SA cols for list input
    call this%ctx%tags(this%param_names, this%nparam, this%input_name)

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

    ! update ts links
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
    integer(I4B) :: icol

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
      ! allocate variable in memory manager
      call this%structarray%mem_create_vector(icol, idt)
    end do
  end subroutine create_structarray

end module ListLoadModule
