!> @brief This module contains the atmospheric boundary condition
!!
!! This module contains the methods used to calculate the heat fluxes
!! for surface-water boundaries, like streams and lakes.  In its current form,
!! this class acts like a package to a package, similar to the TVK package that
!! can be invoked from the NPF package.  Once this package is completed in its
!! prototyped form, it will likely be moved around.
!<

! SFR flows (sfrbudptr)     index var     SFE term      Equation
! ---------------------------------------------------------------------------------
! -- PBST terms
! SHORTWAVE RADIATION       idxbudswr     SHORTWAVE     (1 - shd) * (1 - swrefl) * solr
! LONGWAVE RADIATION        idxbudlwr     LONGWAVE      longwv_in * (1 - lwrefl) + longwv_out
! LATENT HEAT FLUX          idxbudlhf     LATENT HEAT   evap_rate * latent_heat_vaporization * rho_w
! SENSIBLE HEAT FLUX        idxbudshf     SENS HEAT     cd * rho_a * C_p_a * wspd * (t_air - t_feat)

module AbcModule
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, DZERO, LENVARNAME, &
                             LENPACKAGENAME, TABLEFT, TABCENTER, LENMEMTYPE, &
                             DHUNDRED, DCTOK
  use KindModule, only: I4B, DP
  use MemoryManagerModule, only: mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use SimModule, only: store_error, count_errors
  use SimVariablesModule, only: errmsg
  use PbstBaseModule, only: PbstBaseType, pbstbase_da
  use SensHeatModule, only: ShfType, shf_cr
  use ShortwaveModule, only: SwrType, swr_cr
  use LatHeatModule, only: LhfType, lhf_cr
  use LongwaveModule, only: LwrType, lwr_cr
  use ObserveModule
  use BudgetObjectModule, only: BudgetObjectType, budgetobject_cr
  use NumericalPackageModule, only: NumericalPackageType
  use TimeSeriesManagerModule, only: TimeSeriesManagerType, tsmanager_cr
  use TableModule, only: TableType, table_cr
  use BndModule, only: BndType
  use GweInputDataModule, only: GweInputDataType

  implicit none

  private

  public :: AbcType
  public :: abc_cr

  character(len=LENVARNAME) :: text = '          ABC'

  type, extends(BndType) :: AbcType

    type(GweInputDataType), pointer :: gwecommon => null() !< pointer to shared gwe data used by multiple packages but set in est

    character(len=8), dimension(:), pointer, contiguous :: status => null() !< active, inactive, constant
    integer(I4B), pointer :: ncv => null() !< number of control volumes
    integer(I4B), dimension(:), pointer, contiguous :: iboundpbst => null() !< package ibound
    character(len=LINELENGTH), pointer, public :: inputFilename => null() !< a particular abc input file name, could be for sensible heat flux or latent heat flux subpackages, for example
    logical, pointer, public :: active => null() !< logical indicating if a atmospheric boundary condition object is active
    ! -- table objects
    !type(TableType), pointer :: inputtab => null() !< input table object
    !
    logical, pointer, public :: swr_active => null() !< logical indicating if a shortwave radiation heat flux object is active
    logical, pointer, public :: lwr_active => null() !< logical indicating if a longwave radiation heat flux object is active
    logical, pointer, public :: lhf_active => null() !< logical indicating if a latent heat flux object is active
    logical, pointer, public :: shf_active => null() !< logical indicating if a sensible heat flux object is active
    !
    type(ShfType), pointer :: shf => null() ! sensible heat flux (shf) object
    type(SwrType), pointer :: swr => null() ! shortwave radiation heat flux (swr) object
    type(LhfType), pointer :: lhf => null() ! latent heat flux (lhf) object
    type(LwrType), pointer :: lwr => null() ! longwave radiation heat flux (lwr) object
    !
    ! -- abc budget object
    type(BudgetObjectType), pointer :: budobj => null() !< ABC budget object
    !
    real(DP), pointer :: rhoa => null() !< desity of air
    real(DP), pointer :: cpa => null() !< heat capacity of air
    real(DP), pointer :: cd => null() !< drag coefficient
    real(DP), pointer :: wfslope => null() !< wind function slope
    real(DP), pointer :: wfint => null() !< wind function intercept
    real(DP), pointer :: lwrefl => null() !< reflectance of longwave radiation by the water surface
    real(DP), pointer :: emissw => null() !< emissivity of water
    real(DP), pointer :: emissr => null() !< emissivity of the riparian canopy
    !
    real(DP), dimension(:), pointer, contiguous :: wspd => null() !< wind speed
    real(DP), dimension(:), pointer, contiguous :: tatm => null() !< temperature of the atmosphere
    real(DP), dimension(:), pointer, contiguous :: solr => null() !< solar radiation
    real(DP), dimension(:), pointer, contiguous :: shd => null() !< shade fraction
    real(DP), dimension(:), pointer, contiguous :: swrefl => null() !< shortwave reflectance of water surface
    real(DP), dimension(:), pointer, contiguous :: rh => null() !< relative humidity
    real(DP), dimension(:), pointer, contiguous :: atmc => null() !< atmospheric composition adjustment
    real(DP), dimension(:), pointer, contiguous :: patm => null() !< atmospheric pressure (mbar)
    !
    real(DP), dimension(:), pointer, contiguous :: ea => null() !< ambient vapor pressure of the atmosphere, internally calculated and used by multiple heat calculations
    real(DP), dimension(:), pointer, contiguous :: es => null() !< saturation vapor pressure at air temperature, make available for multiple heat flux calculations
    real(DP), dimension(:), pointer, contiguous :: ew => null() !< saturation vapor pressure at water temperature, make available for multiple heat flux calculations

  contains

    procedure :: da => abc_da
    procedure :: ar
    procedure, public :: abc_rp
    procedure :: abc_check_valid
    procedure :: bnd_options => abc_read_options
    procedure :: abc_set_stressperiod
    procedure :: abc_allocate_arrays
    procedure, private :: abc_allocate_scalars
    procedure, public :: abc_cq
    procedure, private :: recalc_shared_vars
    procedure, private :: calc_eatm !< function for calculating ambient vapor pressure of the atmosphere

  end type AbcType

contains

  !> @brief Create a new AbcType object
  !!
  !! Create a new atmospheric boundary condition (AbcType) object. Initially for use with
  !! the SFE package.
  !<
  subroutine abc_cr(abc, name_model, inunit, iout, fname, ncv, gwecommon)
    ! -- dummy
    type(AbcType), pointer, intent(out) :: abc
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=LINELENGTH), intent(in) :: fname
    integer(I4B), target, intent(in) :: ncv
    type(GweInputDataType), intent(in), target :: gwecommon !< shared data container for use by multiple GWE packages
    !
    ! -- Create the object
    allocate (abc)
    !
    call abc%set_names(1, name_model, 'ABC', 'ABC')
    !
    !abc%text = text
    !
    ! -- allocate scalars
    call abc%abc_allocate_scalars()
    !
    abc%inunit = inunit
    abc%iout = iout
    abc%inputFilename = fname
    abc%ncv => ncv
    call abc%parser%Initialize(abc%inunit, abc%iout)
    !
    ! -- initialize associated abc utilities
    call shf_cr(abc%shf, name_model, inunit, iout, ncv)
    call swr_cr(abc%swr, name_model, inunit, iout, ncv)
    call lhf_cr(abc%lhf, name_model, inunit, iout, ncv)
    call lwr_cr(abc%lwr, name_model, inunit, iout, ncv)
    !
    ! -- Create time series manager
    call tsmanager_cr(abc%tsmanager, abc%iout, &
                      removeTsLinksOnCompletion=.true., &
                      extendTsToEndOfSimulation=.true.)
    !
    ! -- Store pointer to shared data module for accessing cpw, rhow
    !    for the heat flux calculations
    abc%gwecommon => gwecommon
  end subroutine abc_cr

  !> @brief Allocate and read
  !!
  !!  Method to allocate and read static data for the SHF, SWR, and LHF sub-utilities
  !<
  subroutine ar(this)
    ! -- dummy
    class(AbcType) :: this !< AbcType object
    ! -- formats
    character(len=*), parameter :: fmtapt = &
      "(1x,/1x,'ABC -- ATMOSPHERIC BOUNDARY CONDITION PACKAGE, VERSION 1, 7/20/2025', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! -- print a message identifying the apt package.
    write (this%iout, fmtapt) this%inunit
    !
  end subroutine ar

  !> @brief ABC read and prepare for setting stress period information
  !<
  subroutine abc_rp(this)
    ! -- module
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    use TdisModule, only: kper, nper
    ! -- dummy
    class(AbcType) :: this !< AbcType object
    ! -- local
    integer(I4B) :: ierr
    integer(I4B) :: n
    logical :: isfound, endOfBlock
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: line
    integer(I4B) :: itemno
    ! -- formats
    character(len=*), parameter :: fmtblkerr = &
      &"('Error.  Looking for BEGIN PERIOD iper.  Found ', a, ' instead.')"
    character(len=*), parameter :: fmtlsp = &
      &"(1X,/1X,'REUSING ',A,'S FROM LAST STRESS PERIOD')"
    !
    ! -- Set ionper to the stress period number for which a new block of data
    !    will be read.
    if (this%inunit == 0) return
    !
    ! -- get stress period data
    if (this%ionper < kper) then
      !
      ! -- get period block
      call this%parser%GetBlock('PERIOD', isfound, ierr, &
                                supportOpenClose=.true., &
                                blockRequired=.false.)
      if (isfound) then
        !
        ! -- read ionper and check for increasing period numbers
        call this%read_check_ionper()
      else
        !
        ! -- PERIOD block not found
        if (ierr < 0) then
          ! -- End of file found; data applies for remainder of simulation.
          this%ionper = nper + 1
        else
          ! -- Found invalid block
          call this%parser%GetCurrentLine(line)
          write (errmsg, fmtblkerr) adjustl(trim(line))
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
      end if
    end if
    !
    ! -- Read data if ionper == kper
    if (this%ionper == kper) then
      !
      ! -- setup table for period data
      if (this%iprpak /= 0) then
        !
        ! -- reset the input table object
        title = trim(adjustl(this%text))//' PACKAGE ('// &
                trim(adjustl(this%packName))//') DATA FOR PERIOD'
        write (title, '(a,1x,i6)') trim(adjustl(title)), kper
        call table_cr(this%inputtab, this%packName, title)
        call this%inputtab%table_df(1, 4, this%iout, finalize=.FALSE.)
        text = 'NUMBER'
        call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
        text = 'KEYWORD'
        call this%inputtab%initialize_column(text, 20, alignment=TABLEFT)
        do n = 1, 2
          write (text, '(a,1x,i6)') 'VALUE', n
          call this%inputtab%initialize_column(text, 15, alignment=TABCENTER)
        end do
      end if
      !
      ! -- read data
      stressperiod: do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        !
        ! -- get feature number
        itemno = this%parser%GetInteger()
        !
        ! -- read data from the rest of the line
        call this%abc_set_stressperiod(itemno)
        !
        ! -- write line to table
        if (this%iprpak /= 0) then
          call this%parser%GetCurrentLine(line)
          call this%inputtab%line_to_columns(line)
        end if
      end do stressperiod
      !
      if (this%iprpak /= 0) then
        call this%inputtab%finalize_table()
      end if
      !
      ! -- using stress period data from the previous stress period
    else
      write (this%iout, fmtlsp) trim(this%filtyp)
    end if
    !
    ! -- write summary of stress period error messages
    ierr = count_errors()
    if (ierr > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine abc_rp

  !> @brief Set options specific to the AbcType
  !!
  !! This routine overrides TspAptType%gc_options
  !<
  subroutine abc_read_options(this, option, found)
    ! -- modules
    use ConstantsModule, only: MAXCHARLEN, LGP
    use InputOutputModule, only: urword, getunit, assign_iounit, openfile
    ! -- dummy
    class(AbcType), intent(inout) :: this
    character(len=*), intent(inout) :: option
    logical, intent(inout) :: found
    ! -- formats
    character(len=*), parameter :: fmtaptbin = &
      "(4x, a, 1x, a, 1x, ' WILL BE SAVED TO FILE: ', a, &
      &/4x, 'OPENED ON UNIT: ', I0)"
    !
    found = .true.
    select case (option)
    case ('DENSITY_AIR')
      this%rhoa = this%parser%GetDouble()
      if (this%rhoa <= 0.0) then
        write (errmsg, '(a)') 'Specified value for the density of &
          &the atmosphere must be greater than 0.0.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      else
        write (this%iout, '(4x,a,1pg15.6)') &
          "The density of the atmosphere has been set to: ", this%rhoa
      end if
    case ('HEAT_CAPACITY_AIR')
      this%cpa = this%parser%GetDouble()
      if (this%cpa <= 0.0) then
        write (errmsg, '(a)') 'Specified value for the heat capacity of &
          &the atmosphere must be greater than 0.0.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      else
        write (this%iout, '(4x,a,1pg15.6)') &
          "The heat capacity of the atmosphere has been set to: ", this%cpa
      end if
    case ('DRAG_COEFFICIENT')
      this%cd = this%parser%GetDouble()
      write (this%iout, '(4x,a,1pg15.6)') &
        "The surface-atmosphere drag coefficient has been set to: ", this%cd
      !end if
    case ('WIND_FUNC_SLOPE')
      this%wfslope = this%parser%GetDouble()
      if (this%wfslope <= 0.0) then
        write (errmsg, '(a)') 'Specified value for the wind function slope &
          &must be greater than 0.0.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      else
        write (this%iout, '(4x,a,1pg15.6)') &
          "The evaporation wind function slope has been set to: ", this%wfslope
      end if
    case ('WIND_FUNC_INT')
      this%wfint = this%parser%GetDouble()
      if (this%wfint <= 0.0) then
        write (errmsg, '(a)') 'Specified value for the wind function intercept &
          &must be greater than 0.0.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      else
        write (this%iout, '(4x,a,1pg15.6)') &
          "The evaporation wind function intercept has been set to: ", this%wfint
      end if
    case ('LONGWAVE_REFLECTANCE')
      this%lwrefl = this%parser%GetDouble()
      if (this%lwrefl <= 0.0) then
        write (errmsg, '(a)') 'Specified value for the reflectance of longwave radiation &
          &must be greater than 0.0.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      else
        write (this%iout, '(4x,a,1pg15.6)') &
          "The reflectance of longwave radiation has been set to: ", this%lwrefl
      end if
    case ('EMISSIVITY_WATER')
      this%emissw = this%parser%GetDouble()
      if (this%emissw <= 0.0) then
        write (errmsg, '(a)') 'Specified value for the emissivity of water &
          &must be greater than 0.0.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      else
        write (this%iout, '(4x,a,1pg15.6)') &
          "The emissivity of water has been set to: ", this%emissw
      end if
    case ('EMISSIVITY_CANOPY')
      this%emissr = this%parser%GetDouble()
      if (this%emissr <= 0.0) then
        write (errmsg, '(a)') 'Specified value for the emissivity of the riparian canopy &
          &must be greater than 0.0.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      else
        write (this%iout, '(4x,a,1pg15.6)') &
          "The emissivity of the riparian canopy has been set to: ", this%emissw
      end if
    case ('SWR_OFF')
      this%swr_active = .false.
      write (this%iout, '(4x,a)') &
        "Shortwave thermal energy exchange between the stream reaches and the &
        &atmosphere has been turned off."
    case ('LWR_OFF')
      this%lwr_active = .false.
      write (this%iout, '(4x,a)') &
        "Longwave thermal energy exchange between stream reaches and the &
        &atmosphere has been turned off."
    case ('LHF_OFF')
      this%lhf_active = .false.
      write (this%iout, '(4x,a)') &
        "Latent heat exchange between stream reaches and the atmosphere &
        &has been turned off."
    case ('SHF_OFF')
      this%shf_active = .false.
      write (this%iout, '(4x,a)') &
        "Sensible heat exchange between the stream reaches and the &
        &atmosphere has been turned off."
    case default
      write (errmsg, '(a,a)') 'Unknown ABC option: ', trim(option)
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end select
  end subroutine abc_read_options

  !> @brief Allocate scalars specific to the streamflow energy transport (SFE)
  !! package.
  !<
  subroutine abc_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(AbcType) :: this
    !
    allocate (this%active)
    allocate (this%inputFilename)
    !
    ! -- initialize
    this%active = .false.
    this%inputFilename = ''
    !
    ! -- allocate
    call mem_allocate(this%shf_active, 'SHF_ACTIVE', this%memoryPath)
    call mem_allocate(this%swr_active, 'SWR_ACTIVE', this%memoryPath)
    call mem_allocate(this%lhf_active, 'LHF_ACTIVE', this%memoryPath)
    call mem_allocate(this%lwr_active, 'LWR_ACTIVE', this%memoryPath)
    !
    ! -- allocate SHF specific
    call mem_allocate(this%rhoa, 'RHOA', this%memoryPath)
    call mem_allocate(this%cpa, 'CPA', this%memoryPath)
    call mem_allocate(this%cd, 'CD', this%memoryPath)
    ! -- allocate LHF specific
    call mem_allocate(this%wfslope, 'WFSLOPE', this%memoryPath)
    call mem_allocate(this%wfint, 'WFINT', this%memoryPath)
    ! -- allocate LWR specific
    call mem_allocate(this%lwrefl, 'LWREFL', this%memoryPath)
    call mem_allocate(this%emissw, 'EMISSW', this%memoryPath)
    call mem_allocate(this%emissr, 'EMISSR', this%memoryPath)
    !
    ! -- initialize to default values
    this%shf_active = .true. ! Initialize to one for 'on'
    this%swr_active = .true.
    this%lhf_active = .true.
    this%lwr_active = .true.
    ! -- initialize to SHF specific default values
    this%rhoa = 1.225 ! kg/m3
    this%cpa = 717.0 ! J/kg/C
    this%cd = 0.002 ! unitless
    ! -- initialize to LHF specific default values
    this%wfslope = 1.383e-08 ! 1/mbar Fogg 2023 (change!)
    this%wfint = 3.445e-09 ! m/s Fogg 2023 (change!)
    ! -- initialize to LWR specific default values
    this%lwrefl = 0.03 ! unitless (Anderson, 1954)
    this%emissw = 0.95 ! unitless (Dingman, 2015)
    this%emissr = 0.97 ! unitless (Sobrino et al, 2005)
    !
    ! -- call standard NumericalPackageType allocate scalars
    call this%BndType%allocate_scalars()
    !
    ! -- allocate time series manager
    allocate (this%tsmanager)
  end subroutine abc_allocate_scalars

  !> @brief Allocate arrays specific to the atmspheric boundary package
  !<
  subroutine abc_allocate_arrays(this)
    ! -- modules
    !! use MemoryManagerModule, only: mem_allocate
    use MemoryManagerModule, only: mem_setptr, mem_checkin, mem_allocate, &
                                   get_mem_type, mem_reallocate
    ! -- dummy
    class(AbcType), intent(inout) :: this
    ! integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    !real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    ! -- local
    integer(I4B) :: n
    !
    ! -- allocate character array for status
    allocate (this%status(this%ncv))
    !
    ! -- initialize arrays
    do n = 1, this%ncv
      this%status(n) = 'ACTIVE'
    end do
    !
    ! -- allocate all atmospheric boundary condition vars, initialize to size 0
    call mem_allocate(this%wspd, 0, 'WSPD', this%memoryPath)
    call mem_allocate(this%tatm, 0, 'TATM', this%memoryPath)
    call mem_allocate(this%solr, 0, 'SOLR', this%memoryPath)
    call mem_allocate(this%shd, 0, 'SHD', this%memoryPath)
    call mem_allocate(this%swrefl, 0, 'SWREFL', this%memoryPath)
    call mem_allocate(this%rh, 0, 'RH', this%memoryPath)
    call mem_allocate(this%atmc, 0, 'ATMC', this%memoryPath)
    call mem_allocate(this%patm, 0, 'PATM', this%memoryPath)
    call mem_allocate(this%ea, 0, 'EA', this%memoryPath)
    call mem_allocate(this%es, 0, 'ES', this%memoryPath)
    call mem_allocate(this%ew, 0, 'EW', this%memoryPath)
    !
    ! -- reallocate abc variables based on which calculations are used
    if (this%shf_active .or. this%lhf_active) then
      call mem_reallocate(this%wspd, this%ncv, 'WSPD', this%memoryPath)
      call mem_reallocate(this%ew, this%ncv, 'ES', this%memoryPath)
      do n = 1, this%ncv
        this%wspd(n) = DZERO
        this%ew(n) = DZERO
      end do
    end if
    !
    if (this%shf_active .or. this%lhf_active .or. this%lwr_active) then
      call mem_reallocate(this%tatm, this%ncv, 'TATM', this%memoryPath)
      call mem_reallocate(this%ea, this%ncv, 'EA', this%memoryPath)
      call mem_reallocate(this%es, this%ncv, 'ES', this%memoryPath)
      call mem_reallocate(this%ew, this%ncv, 'EW', this%memoryPath)
      do n = 1, this%ncv
        this%tatm(n) = DZERO
        this%ea(n) = DZERO
        this%es(n) = DZERO
        this%ew(n) = DZERO
      end do
    end if
    if (this%swr_active) then
      call mem_reallocate(this%solr, this%ncv, 'SOLR', this%memoryPath)
      call mem_reallocate(this%swrefl, this%ncv, 'SWREFL', this%memoryPath)
      do n = 1, this%ncv
        this%solr(n) = DZERO
        this%swrefl(n) = DZERO
      end do
    end if
    if (this%swr_active .or. this%lwr_active) then
      call mem_reallocate(this%shd, this%ncv, 'SHD', this%memoryPath)
      do n = 1, this%ncv
        this%shd(n) = DZERO
      end do
    end if
    if (this%lhf_active .or. this%lwr_active) then
      call mem_reallocate(this%rh, this%ncv, 'RH', this%memoryPath)
      do n = 1, this%ncv
        this%rh(n) = DZERO
      end do
    end if
    if (this%lwr_active) then
      call mem_reallocate(this%atmc, this%ncv, 'ATMC', this%memoryPath)
      do n = 1, this%ncv
        this%atmc(n) = DZERO
      end do
    end if
    if (this%shf_active) then
      call mem_reallocate(this%patm, this%ncv, 'PATM', this%memoryPath)
      do n = 1, this%ncv
        this%patm(n) = DZERO
      end do
    end if
    !
    ! -- call utility ar routines if active
    if (this%swr_active) then
      call this%swr%pbst_ar()
    end if
    if (this%lwr_active) then
      call this%lwr%pbst_ar()
    end if
    if (this%lhf_active) then
      call this%lhf%pbst_ar()
    end if
    if (this%shf_active) then
      call this%shf%pbst_ar()
    end if
  end subroutine abc_allocate_arrays

  !> @brief Deallocate memory
  !<
  subroutine abc_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(AbcType) :: this
    !
    ! -- SHF (sensible heat flux)
    if (this%shf_active) then
      call this%shf%da()
      deallocate (this%shf)
    end if
    ! -- SWR (shortwave radiation heat flux)
    if (this%swr_active) then
      call this%swr%da()
      deallocate (this%swr)
    end if
    ! -- LHF (latent heat flux)
    if (this%lhf_active) then
      call this%lhf%da()
      deallocate (this%lhf)
    end if
    ! -- LWR (longwave radiation heat flux)
    if (this%lwr_active) then
      call this%lwr%da()
      deallocate (this%lwr)
    end if
    !
    ! -- Deallocate scalars
    call mem_deallocate(this%ncv)
    call mem_deallocate(this%shf_active)
    call mem_deallocate(this%swr_active)
    call mem_deallocate(this%lhf_active)
    call mem_deallocate(this%lwr_active)
    call mem_deallocate(this%rhoa)
    call mem_deallocate(this%cpa)
    call mem_deallocate(this%cd)
    call mem_deallocate(this%wfslope)
    call mem_deallocate(this%wfint)
    call mem_deallocate(this%lwrefl)
    call mem_deallocate(this%emissw)
    call mem_deallocate(this%emissr)
    !
    ! -- Deallocate time series manager
    deallocate (this%tsmanager)
    !
    ! -- Deallocate arrays
    call mem_deallocate(this%status)
    call mem_deallocate(this%iboundpbst)
    call mem_deallocate(this%wspd)
    call mem_deallocate(this%tatm)
    call mem_deallocate(this%shd)
    call mem_deallocate(this%swrefl)
    call mem_deallocate(this%solr)
    call mem_deallocate(this%rh)
    call mem_deallocate(this%atmc)
    call mem_deallocate(this%ea)
    call mem_deallocate(this%es)
    call mem_deallocate(this%ew)
    call mem_deallocate(this%patm)
    !
    ! -- Deallocate scalars in TspAptType
    call this%NumericalPackageType%da() ! this may not work -- revisit and cleanup !!!
  end subroutine abc_da

!  !> @brief Calculate observation value and pass it back to APT
!  !<
!  subroutine abc_bd_obs(this, obstypeid, jj, v, found)
!    ! -- dummy
!    class(AbcType), intent(inout) :: this
!    character(len=*), intent(in) :: obstypeid
!    real(DP), intent(inout) :: v
!    integer(I4B), intent(in) :: jj
!    logical, intent(inout) :: found
!    ! -- local
!    !integer(I4B) :: n1, n2
!    !
!    found = .true.
!    !select case (obstypeid)
!    !case ('SHF')
!    !  if (this%iboundpak(jj) /= 0) then
!    !    call this%sfe_shf_term(jj, n1, n2, v)
!    !  end if
!    !case ('SWR')
!    !  if (this%iboundpak(jj) /= 0) then
!    !    call this%sfe_swr_term(jj, n1, n2, v)
!    !  end if
!    !case default
!    !  found = .false.
!    !end select
!  end subroutine abc_bd_obs

  !> @brief Calculate Atmospheric-Stream Heat Flux
  !!
  !! Calculate and return the atmospheric heat flux for one reach
  !<
  subroutine abc_cq(this, ifno, tstrm, abcflx, obstype)
    ! -- dummy
    class(AbcType), intent(inout) :: this
    integer(I4B), intent(in) :: ifno !< stream reach integer id
    real(DP), intent(in) :: tstrm !< temperature of the stream reach
    real(DP), intent(inout) :: abcflx !< calculated atmospheric boundary flux amount
    character(len=*), optional, intent(in) :: obstype !< when present, subroutine will return a specific energy flux for an observation
    ! -- local
    real(DP) :: swrflx = DZERO
    real(DP) :: lwrflx = DZERO
    real(DP) :: lhfflx = DZERO
    real(DP) :: shfflx = DZERO
    !
    ! -- update shared variables
    call this%recalc_shared_vars(ifno, tstrm)
    !
    ! -- calculate shortwave radiation using HGS equation
    if (this%swr_active) then
      call this%swr%swr_cq(ifno, swrflx)
    end if
    !
    ! -- calculate longwave radiation
    if (this%lwr_active) then
      call this%lwr%lwr_cq(ifno, tstrm, lwrflx)
    end if
    !
    ! -- calculate latent heat flux using Dalton-like mass transfer equation
    if (this%lhf_active) then
      call this%lhf%lhf_cq(ifno, tstrm, this%gwecommon%gwerhow, lhfflx)
    end if
    !
    ! -- calculate sensible heat flux using HGS equation
    if (this%shf_active .and. .not. this%lhf_active) then
      call this%shf%shf_cq(ifno, tstrm, shfflx) ! default to Bowen ratio method ("1")
    else if (this%shf_active .and. this%lhf_active) then
      call this%shf%shf_cq(ifno, tstrm, shfflx, lhfflx) ! use Bowen ratio method ("2")
    end if
    !
    if (present(obstype)) then
      select case (obstype)
      case ('swr')
        abcflx = swrflx
      case ('lwr')
        abcflx = lwrflx
      case ('lhf')
        abcflx = lhfflx
      case ('shf')
        abcflx = shfflx
      case default
        errmsg = 'Unrecognized observation type "'// &
                 trim(obstype)//'" for '// &
                 trim(adjustl(this%text))//' utility.'
        call store_error(errmsg, terminate=.TRUE.)
      end select
    else
      abcflx = swrflx + lwrflx + shfflx - lhfflx
    end if
  end subroutine abc_cq

  !> @brief Recalculate variables that are used by various heat fluxes
  !!
  !! For parameters like ambient vapor pressure of the atmosphere that are
  !! used in the calculation of multiple heat fluxes, in this case longwave,
  !! latent, and sensible heat flux, need to update the value held in memory
  !<
  subroutine recalc_shared_vars(this, ifno, tstrm)
    ! -- dummy
    class(AbcType), intent(inout) :: this
    integer(I4B), intent(in) :: ifno !< reach id
    real(DP), intent(in) :: tstrm !< temperature of the stream reach
    !
    ! -- calculate saturation vapor pressure at air temperature
    this%es(ifno) = calc_sat_vap_pres(this%tatm(ifno))
    !
    ! -- calculate ambient vapor pressure at the atmospheric temperature
    this%ea(ifno) = this%calc_eatm(ifno)
    !
    ! -- calculate saturation vapor pressure at the water temperature
    this%ew(ifno) = calc_sat_vap_pres(tstrm)
  end subroutine recalc_shared_vars

  !> @brief Calculate saturated vapor pressure for a given temperature
  !!
  !! A function for calculating the saturated vapor pressure given a
  !! temperature, commonly either the stream temperature or atmospheric
  !! temperature.
  !<
  function calc_sat_vap_pres(temp) result(e)
    ! -- dummy
    real(DP), intent(in) :: temp
    ! -- return
    real(DP) :: e
    !
    e = 6.1275_DP * exp(17.2693882_DP * &
                        (temp / (temp + DCTOK - 35.86_DP)))
  end function calc_sat_vap_pres

  !> @brief Calculate ambient vapor pressure
  !!
  !! Calculate ambient vapor pressure of the atmosphere as a function of
  !! relative humidity and saturation vapor pressure of the atmosphere
  !<
  function calc_eatm(this, ifno) result(eatm)
    ! -- dummy
    class(AbcType) :: this
    integer(I4B) :: ifno
    ! -- return
    real(DP) :: eatm
    !
    eatm = this%rh(ifno) / DHUNDRED * this%es(ifno)
  end function calc_eatm

  !> @brief Set the stress period attributes based on the keyword
  !<
  subroutine abc_set_stressperiod(this, itemno) !, keyword, found) MAKING LIKE PBST SET STRESSPERIOD
    ! -- module
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(AbcType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    !logical, intent(inout) :: found
    ! -- local
    character(len=LINELENGTH) :: text
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    integer(I4B) :: jj
    real(DP), pointer :: bndElem => null()
    !
    ! <wspd> WIND SPEED
    ! <tatm> TEMPERATURE OF THE ATMOSPHERE
    ! <shd> SHADE
    ! <swrefl> REFLECTANCE OF SHORTWAVE RADIATION OFF WATER SURFACE
    ! <solr> SOLAR RADIATION
    ! <rh> RELATIVE HUMIDITY
    ! <atmc> ATMOSPHERIC COMPOSITION
    ! <patm> ATMOSPHERIC PRESSURE
    !
    ! -- read line
    call this%parser%GetStringCaps(keyword)
    select case (keyword)
    case ('STATUS')
      ierr = this%abc_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetStringCaps(text)
      this%status(itemno) = text(1:8)
      if (text == 'CONSTANT') then
        this%iboundpbst(itemno) = -1
      else if (text == 'INACTIVE') then
        this%iboundpbst(itemno) = 0
      else if (text == 'ACTIVE') then
        this%iboundpbst(itemno) = 1
      else
        write (errmsg, '(a,a)') &
          'Unknown '//trim(this%text)//' status keyword: ', text//'.'
        call store_error(errmsg)
      end if
    case ('WSPD')
      ierr = this%abc_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1
      bndElem => this%wspd(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'WSPD')
    case ('TATM')
      ierr = this%abc_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1
      bndElem => this%tatm(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'TATM')
    case ('SHD')
      ierr = this%abc_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1
      bndElem => this%shd(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'SHD')
    case ('SWREFL')
      ierr = this%abc_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1
      bndElem => this%swrefl(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'SWREFL')
    case ('SOLR')
      ierr = this%abc_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1
      bndElem => this%solr(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'SOLR')
    case ('RH')
      ierr = this%abc_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1
      bndElem => this%rh(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'RH')
    case ('ATMC')
      ierr = this%abc_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1
      bndElem => this%atmc(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'ATMC')
    case ('PATM')
      ierr = this%abc_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1
      bndElem => this%patm(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'PATM')
    case default
      !
      ! -- Keyword not recognized so return to caller with found = .false.
      !found = .false.
    end select
    !
999 continue
  end subroutine abc_set_stressperiod

  !> @brief Determine if a valid feature number has been specified.
  !<
  function abc_check_valid(this, itemno) result(ierr)
    ! -- return
    integer(I4B) :: ierr
    ! -- dummy
    class(AbcType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    !
    ierr = 0
    if (itemno < 1 .or. itemno > this%ncv) then
      write (errmsg, '(a,1x,i6,1x,a,1x,i6)') &
        'Featureno ', itemno, 'must be > 0 and <= ', this%ncv
      call store_error(errmsg)
      ierr = 1
    end if
  end function abc_check_valid

end module AbcModule
