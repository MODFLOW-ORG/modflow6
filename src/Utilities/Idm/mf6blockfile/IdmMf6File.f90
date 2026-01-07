!> @brief This module contains the IdmMf6FileModule
!!
!! This module contains high-level routines for loading
!! MODFLOW 6 ASCII source input. This module implements the
!! loader types that the IdmLoadModule creates and invokes.
!! It also creates and manages dynamic ASCII input loaders
!! for all supported types of MODFLOW 6 ASCII dynamic input.
!!
!<
module IdmMf6FileModule

  use KindModule, only: DP, I4B, LGP
  use SimVariablesModule, only: errmsg
  use ConstantsModule, only: LINELENGTH
  use SimModule, only: store_error, store_error_filename
  use BlockParserModule, only: BlockParserType
  use ModflowInputModule, only: ModflowInputType
  use InputLoadTypeModule, only: StaticPkgLoadBaseType, DynamicPkgLoadBaseType
  use AsciiInputLoadTypeModule, only: AsciiDynamicPkgLoadBaseType
  use NCFileVarsModule, only: NCPackageVarsType

  implicit none
  private
  public :: input_load
  public :: Mf6FileStaticPkgLoadType, Mf6FileDynamicPkgLoadType
  public :: open_mf6file

  !> @brief MF6File static loader type
  !<
  type, extends(StaticPkgLoadBaseType) :: Mf6FileStaticPkgLoadType
  contains
    procedure :: init => static_init
    procedure :: load => static_load
    procedure :: destroy => static_destroy
  end type Mf6FileStaticPkgLoadType

  !> @brief MF6File dynamic loader type
  !<
  type, extends(DynamicPkgLoadBaseType) :: Mf6FileDynamicPkgLoadType
    type(BlockParserType), pointer :: parser !< parser for MF6File period blocks
    integer(I4B), pointer :: iper => null() !< memory managed variable, loader iper
    integer(I4B), pointer :: ionper => null() !< memory managed variable, next load period
    class(AsciiDynamicPkgLoadBaseType), pointer :: rp_loader => null()
  contains
    procedure :: init => dynamic_init
    procedure :: df => dynamic_df
    procedure :: ad => dynamic_ad
    procedure :: rp => dynamic_rp
    procedure :: read_ionper => dynamic_read_ionper
    procedure :: create_loader => dynamic_create_loader
    procedure :: destroy => dynamic_destroy
  end type Mf6FileDynamicPkgLoadType

contains

  !> @brief input load for traditional mf6 simulation static input file
  !<
  subroutine input_load(filename, mf6_input, component_filename, iout, nc_vars)
    use LoadMf6FileModule, only: LoadMf6FileType
    character(len=*), intent(in) :: filename
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), intent(in) :: component_filename !< component (e.g. model) filename
    integer(I4B), intent(in) :: iout !< unit number for output
    type(NCPackageVarsType), pointer, optional, intent(in) :: nc_vars
    type(BlockParserType), allocatable, target :: parser !< block parser
    type(NCPackageVarsType), pointer :: netcdf_vars
    type(LoadMf6FileType) :: loader
    integer(I4B) :: inunit

    if (present(nc_vars)) then
      netcdf_vars => nc_vars
    else
      nullify (netcdf_vars)
    end if

    ! open input file
    inunit = open_mf6file(mf6_input%pkgtype, filename, component_filename, iout)

    ! allocate and initialize parser
    allocate (parser)
    call parser%Initialize(inunit, iout)

    ! invoke the load routine
    call loader%load(parser, mf6_input, netcdf_vars, filename, iout)

    ! clear parser file handles
    call parser%clear()

    ! cleanup
    deallocate (parser)
  end subroutine input_load

  !> @brief static loader init
  !<
  subroutine static_init(this, mf6_input, component_name, component_input_name, &
                         input_name)
    class(Mf6FileStaticPkgLoadType), intent(inout) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), intent(in) :: component_name
    character(len=*), intent(in) :: component_input_name
    character(len=*), intent(in) :: input_name
    ! initialize base type
    call this%StaticPkgLoadType%init(mf6_input, component_name, &
                                     component_input_name, input_name)
  end subroutine static_init

  !> @brief load routine for static loader
  !<
  function static_load(this, iout) result(rp_loader)
    class(Mf6FileStaticPkgLoadType), intent(inout) :: this
    integer(I4B), intent(in) :: iout
    class(DynamicPkgLoadBaseType), pointer :: rp_loader
    class(Mf6FileDynamicPkgLoadType), pointer :: mf6_loader

    ! initialize return pointer
    nullify (rp_loader)

    ! load model package to input context
    if (this%iperblock > 0) then
      ! allocate dynamic loader
      allocate (mf6_loader)
      ! point to nc_vars structure
      mf6_loader%nc_vars => this%nc_vars
      ! nullify nc_vars pointer so it isn't deallocated
      nullify (this%nc_vars)
      ! initialize dynamic loader
      call mf6_loader%init(this%mf6_input, this%component_name, &
                           this%component_input_name, this%input_name, &
                           this%iperblock, iout)
      ! set return pointer to base dynamic loader
      rp_loader => mf6_loader
    else
      ! load static input
      call input_load(this%input_name, this%mf6_input, &
                      this%component_input_name, iout, this%nc_vars)
    end if
  end function static_load

  !> @brief static loader destroy
  !<
  subroutine static_destroy(this)
    class(Mf6FileStaticPkgLoadType), intent(inout) :: this
    ! deallocate base type
    call this%StaticPkgLoadType%destroy()
  end subroutine static_destroy

  !> @brief dynamic loader init
  !<
  subroutine dynamic_init(this, mf6_input, component_name, component_input_name, &
                          input_name, iperblock, iout)
    use MemoryManagerModule, only: mem_allocate
    use InputDefinitionModule, only: InputParamDefinitionType
    class(Mf6FileDynamicPkgLoadType), intent(inout) :: this
    type(ModflowInputType), intent(in) :: mf6_input
    character(len=*), intent(in) :: component_name
    character(len=*), intent(in) :: component_input_name
    character(len=*), intent(in) :: input_name
    integer(I4B), intent(in) :: iperblock
    integer(I4B), intent(in) :: iout
    integer(I4B) :: inunit

    ! initialize base loader
    call this%DynamicPkgLoadType%init(mf6_input, component_name, &
                                      component_input_name, input_name, &
                                      iperblock, iout)
    ! allocate scalars
    call mem_allocate(this%iper, 'IPER', mf6_input%mempath)
    call mem_allocate(this%ionper, 'IONPER', mf6_input%mempath)

    ! initialize
    this%iper = 0
    this%ionper = 0

    ! open input file
    inunit = open_mf6file(mf6_input%pkgtype, input_name, &
                          component_input_name, iout)

    ! allocate and initialize parser
    allocate (this%parser)
    call this%parser%Initialize(inunit, iout)

    ! allocate and initialize loader
    call this%create_loader()
  end subroutine dynamic_init

  !> @brief define routine for dynamic loader
  !<
  subroutine dynamic_df(this)
    class(Mf6FileDynamicPkgLoadType), intent(inout) :: this
    ! invoke loader define
    call this%rp_loader%df()
    ! read first ionper
    call this%read_ionper()
  end subroutine dynamic_df

  !> @brief advance routine for dynamic loader
  !<
  subroutine dynamic_ad(this)
    class(Mf6FileDynamicPkgLoadType), intent(inout) :: this
    ! invoke loader advance
    call this%rp_loader%ad()
  end subroutine dynamic_ad

  !> @brief read and prepare routine for dynamic loader
  !<
  subroutine dynamic_rp(this)
    use TdisModule, only: kper, nper
    class(Mf6FileDynamicPkgLoadType), intent(inout) :: this

    ! check if ready to load
    if (this%ionper /= kper) return

    ! dynamic load
    call this%rp_loader%rp(this%parser)

    ! update loaded iper
    this%iper = kper

    ! read next iper
    if (kper < nper) then
      call this%read_ionper()
    else
      this%ionper = nper + 1
    end if
  end subroutine dynamic_rp

  !> @brief dynamic loader read ionper of next period block
  !<
  subroutine dynamic_read_ionper(this)
    use TdisModule, only: kper, nper
    class(Mf6FileDynamicPkgLoadType), intent(inout) :: this
    character(len=LINELENGTH) :: line
    logical(LGP) :: isblockfound
    integer(I4B) :: ierr
    character(len=*), parameter :: fmtblkerr = &
      &"('Looking for BEGIN PERIOD iper.  Found ', a, ' instead.')"

    call this%parser%GetBlock('PERIOD', isblockfound, ierr, &
                              supportOpenClose=.true., &
                              blockRequired=.false.)
    ! set first period block IPER
    if (isblockfound) then
      this%ionper = this%parser%GetInteger()
      if (this%ionper <= this%iper) then
        write (errmsg, '(a, i0, a, i0, a, i0, a)') &
          'Error in stress period ', kper, &
          '. Period numbers not increasing.  Found ', this%ionper, &
          ' but last period block was assigned ', this%iper, '.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      end if
    else
      ! PERIOD block not found
      if (ierr < 0) then
        ! End of file found; data applies for remainder of simulation.
        this%ionper = nper + 1
      else
        ! Found invalid block
        call this%parser%GetCurrentLine(line)
        write (errmsg, fmtblkerr) adjustl(trim(line))
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      end if
    end if
  end subroutine dynamic_read_ionper

  !> @brief allocate a dynamic loader based on load context
  !<
  subroutine dynamic_create_loader(this)
    use LayerArrayLoadModule, only: LayerArrayLoadType
    use GridArrayLoadModule, only: GridArrayLoadType
    use ListLoadModule, only: ListLoadType
    use Mf6FileSettingLoadModule, only: SettingLoadType
    use Mf6FileStoInputModule, only: StoInputType
    use FeatureFlagsModule, only: developmode
    class(Mf6FileDynamicPkgLoadType), intent(inout) :: this
    class(ListLoadType), pointer :: list_loader
    class(GridArrayLoadType), pointer :: arrgrid_loader
    class(LayerArrayLoadType), pointer :: arrlayer_loader
    class(SettingLoadType), pointer :: setting_loader
    class(StoInputType), pointer :: sto_loader

    ! allocate and set loader
    if (this%mf6_input%subcomponent_type == 'STO') then
      allocate (sto_loader)
      this%rp_loader => sto_loader
    else if (this%has_setting) then
      allocate (setting_loader)
      this%rp_loader => setting_loader
    else if (this%readasarrays) then
      allocate (arrlayer_loader)
      this%rp_loader => arrlayer_loader
    else if (this%readarraygrid) then
      call developmode('Input file "'//trim(this%input_name)// &
        '" READARRAYGRID option is still under development, install the &
        &nightly build or compile from source with IDEVELOPMODE = 1.', &
        this%iout)
      allocate (arrgrid_loader)
      this%rp_loader => arrgrid_loader
    else
      allocate (list_loader)
      this%rp_loader => list_loader
    end if

    ! set nc_vars pointer
    this%rp_loader%nc_vars => this%nc_vars

    ! initialize loader
    call this%rp_loader%ainit(this%mf6_input, &
                              this%component_name, &
                              this%component_input_name, &
                              this%input_name, &
                              this%iperblock, &
                              this%parser, &
                              this%iout)
  end subroutine dynamic_create_loader

  !> @brief dynamic loader destroy
  !<
  subroutine dynamic_destroy(this)
    use MemoryManagerModule, only: mem_deallocate
    class(Mf6FileDynamicPkgLoadType), intent(inout) :: this

    ! deallocate scalars
    call mem_deallocate(this%iper)
    call mem_deallocate(this%ionper)

    ! deallocate loader
    nullify (this%rp_loader%nc_vars)
    call this%rp_loader%destroy()
    deallocate (this%rp_loader)

    ! deallocate parser
    call this%parser%clear()
    deallocate (this%parser)

    ! deallocate input context
    call this%DynamicPkgLoadType%destroy()
  end subroutine dynamic_destroy

  !> @brief open a model package files
  !<
  function open_mf6file(filetype, filename, component_fname, iout) result(inunit)
    use InputOutputModule, only: openfile, getunit
    character(len=*), intent(in) :: filetype
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: component_fname
    integer(I4B), intent(in) :: iout
    integer(I4B) :: inunit

    ! initialize
    inunit = 0

    if (filename /= '') then
      ! get unit number and open file
      inunit = getunit()
      call openfile(inunit, iout, trim(adjustl(filename)), filetype, &
                    'FORMATTED', 'SEQUENTIAL', 'OLD')
    else
      write (errmsg, '(a,a,a)') &
        'File unspecified, cannot load model or package &
        &type "', trim(filetype), '".'
      call store_error(errmsg)
      call store_error_filename(component_fname)
    end if
  end function open_mf6file

end module IdmMf6FileModule
