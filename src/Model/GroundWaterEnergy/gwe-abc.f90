!> @brief This module contains the atmospheric boundary condition
!!
!! This module contains the methods used to calculate the heat fluxes
!! for surface-water boundaries, like streams and lakes.  In its current form,
!! this class acts like a package to a package, similar to the TVK package that
!! can be invoked from the NPF package.  Once this package is completed in its
!! prototyped form, it will likely be moved around.
!<
! SFR flows (sfrbudptr)     index var     SFE term              Transport Type
! ---------------------------------------------------------------------------------    
! -- PBST terms
! SENSIBLE HEAT FLUX        idxbudshf     SENS HEAT             cd * rho_a * C_p_a * wspd * (t_air - t_feat)
! SHORTWAVE RADIATION       idxbudsWR     SHORTWAVE             (1 - shd) * (1 - swrefl) * solr


module AbcModule
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, DZERO, LENVARNAME, &
                             LENPACKAGENAME, TABLEFT, TABCENTER, LENMEMTYPE
  use KindModule, only: I4B, DP
  use MemoryManagerModule, only: mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use SimModule, only: store_error, count_errors
  use SimVariablesModule, only: errmsg
  use PbstBaseModule, only: PbstBaseType, pbstbase_da
  use SensHeatModule, only: ShfType, shf_cr
  use ShortwaveModule, only: SwrType, swr_cr
  use LatHeatModule, only: LhfType, lhf_cr
  !use BndModule, only: BndType, AddBndToList, GetBndFromList
  !use TspAptModule, only: TspAptType
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

  !type, extends(NumericalPackageType) :: AbcType
  type, extends(BndType) :: AbcType
      
    type(GweInputDataType), pointer :: gwecommon => null() !< pointer to shared gwe data used by multiple packages but set in est    

    character(len=8), dimension(:), pointer, contiguous :: status => null() !< active, inactive, constant
    integer(I4B), pointer :: ncv => null() !< number of control volumes
    integer(I4B), dimension(:), pointer, contiguous :: iboundpbst => null() !< package ibound
    !character(len=LENPACKAGENAME) :: text = '' !< text string for package transport term
    character(len=LINELENGTH), pointer, public :: inputFilename => null() !< a particular abc input file name, could be for sensible heat flux or latent heat flux subpackages, for example
    logical, pointer, public :: active => null() !< logical indicating if a atmospheric boundary condition object is active
    ! -- table objects
    !type(TableType), pointer :: inputtab => null() !< input table object
  
    logical, pointer, public :: shf_active => null() !< logical indicating if a sensible heat flux object is active
    logical, pointer, public :: swr_active => null() !< logical indicating if a shortwave radition heat flux object is active
    logical, pointer, public :: lhf_active => null() !< logical indicating if a latent heat flux object is active

    type(ShfType), pointer :: shf => null() ! sensible heat flux (shf) object
    type(SwrType), pointer :: swr => null() ! shortwave radiation heat flux (swr) object
    type(LhfType), pointer :: lhf => null() ! latent heat flux (lhf) object
    ! -- abc budget object
    type(BudgetObjectType), pointer :: budobj => null() !< ABC budget object
 
    integer(I4B), pointer :: inshf => null() ! SHF (sensible heat flux utility) unit number (0 if unused)
    integer(I4B), pointer :: inswr => null() ! SWR (shortwave radiation heat flux utility) unit number (0 if unused)
    integer(I4B), pointer :: inlhf => null() ! LHF (latent heat flux utility) unit number (0 if unused)

    real(DP), pointer :: rhoa => null() !< desity of air
    real(DP), pointer :: cpa => null() !< heat capacity of air
    real(DP), pointer :: cd => null() !< drag coefficient
    real(DP), pointer :: wfslope => null() !< wind function slope
    real(DP), pointer :: wfint => null() !< wind function intercept
    
    real(DP), dimension(:), pointer, contiguous :: wspd => null() !< wind speed
    real(DP), dimension(:), pointer, contiguous :: tatm => null() !< temperature of the atmosphere
    real(DP), dimension(:), pointer, contiguous :: solr => null() !< solar radiation
    real(DP), dimension(:), pointer, contiguous :: shd => null() !< shade fraction
    real(DP), dimension(:), pointer, contiguous :: swrefl => null() !< shortwave reflectance of water surface
    real(DP), dimension(:), pointer, contiguous :: rh => null() !< relative humidity
    
  contains

    procedure :: da => abc_da
    !procedure :: init
    procedure :: ar
    procedure, public :: abc_rp
    procedure :: abc_check_valid
    procedure :: bnd_options => abc_read_options
    procedure :: read_option => abc_read_option ! reads stress period
    procedure :: abc_read_options ! read options block
    !procedure :: subpck_set_stressperiod => abc_set_stressperiod
    procedure :: abc_set_stressperiod
    procedure :: abc_allocate_arrays
    procedure, private :: abc_allocate_scalars
    procedure, public :: abc_cq
    procedure, private :: abc_shf_term
    procedure, private :: abc_swr_term
    procedure, private :: abc_lhf_term
    ! -- budget
    !procedure, private :: abc_setup_shfobj
    !procedure, private :: abc_setup_swrobj

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
    ! -- Set pointers to SHF package variables
    if (this%inshf) then
      call this%shf%ar_set_pointers()
    end if
    ! -- Set pointers to LHF package variables
    if (this%inlhf) then
      call this%lhf%ar_set_pointers()
    end if
    !
    ! -- Allocate arrays
    !call this%pbst_allocate_arrays()
    !
    ! -- Read options
    !call this%read_options()
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
    ! -- local
    character(len=LINELENGTH) :: fname
    character(len=MAXCHARLEN) :: keyword
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
      !if (this%cd <= 0.0) then
      !  write (errmsg, '(a)') 'Specified value for the drag coefficient &
      !    &must be greater than 0.0.'
      !  call store_error(errmsg)
      !  call this%parser%StoreErrorUnit()
      !else
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
    !use MemoryHelperModule, only: create_mem_path !! NEW KF 7/24 from npf pkg
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
    
    call mem_allocate(this%inshf, 'INSHF', this%memoryPath)
    call mem_allocate(this%inswr, 'INSWR', this%memoryPath)
    call mem_allocate(this%inlhf, 'INLHF', this%memoryPath)
    ! -- allocate SHF specific
    call mem_allocate(this%rhoa, 'RHOA', this%memoryPath)
    call mem_allocate(this%cpa, 'CPA', this%memoryPath)
    call mem_allocate(this%cd, 'CD', this%memoryPath)
    ! -- allocate LHF specific
    call mem_allocate(this%wfslope, 'WFSLOPE', this%memoryPath)
    call mem_allocate(this%wfint, 'WFINT', this%memoryPath)
    !
    ! -- initialize to default values
    this%shf_active = .false.
    this%swr_active = .false.
    this%lhf_active = .false.
    this%inshf = 1 ! Initialize to one for 'on'
    this%inswr = 1
    this%inlhf = 1
    ! -- initalize to SHF specific default values
    this%rhoa = 1.225 ! kg/m3
    this%cpa = 717.0 ! J/kg/C
    this%cd = 0.002 ! unitless
    ! -- initalize to LHF specific default values
    this%wfslope = 1.383e-01 ! 1/mbar Fogg 2023 (change!)
    this%wfint = 3.445e-09 ! m/s Fogg 2023 (change!)
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
                                   get_mem_type
    ! -- dummy
    class(AbcType), intent(inout) :: this
    ! integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    !real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    ! -- local
    integer(I4B) :: n
    character(len=LENMEMTYPE) :: var_type !< memory type
    !
    ! -- allocate character array for status
    allocate (this%status(this%ncv))
    !
    ! -- initialize arrays
    do n = 1, this%ncv
      this%status(n) = 'ACTIVE'
    end do
    !
    ! -- Call sub-package(s) allocate arrays
    if (this%inshf /= 0) then
       ! -- allocate base arrays
       !call this%BndExtType%allocate_arrays(nodelist, auxvar)
       !
       ! -- set WSPD and TATM context pointer
       call mem_allocate(this%wspd, this%ncv, 'WSPD', this%memoryPath)
       call mem_allocate(this%tatm, this%ncv, 'TATM', this%memoryPath)
       ! -- initialize
       do n = 1, this%ncv
         this%wspd(n) = DZERO
         this%tatm(n) = DZERO
       end do
       call this%shf%ar()
    end if
    if (this%inswr /= 0) then
       call mem_allocate(this%solr, this%ncv, 'SOLR', this%memoryPath)
       call mem_allocate(this%shd, this%ncv, 'SHD', this%memoryPath)
       call mem_allocate(this%swrefl, this%ncv, 'SWREFL', this%memoryPath)
       !
       ! -- initialize
       do n = 1, this%ncv
         this%solr(n) = DZERO
         this%shd(n) = DZERO
         this%swrefl(n) = DZERO
       end do
       call this%swr%ar()
    end if
    if (this%inlhf /= 0) then
       call get_mem_type('WSPD', this%memoryPath, var_type)
       if (var_type == 'UNKNOWN') then
         call mem_allocate(this%wspd, this%ncv, 'WSPD', this%memoryPath)
         do n = 1, this%ncv
           this%wspd(n) = DZERO
         end do
       end if
       call get_mem_type('TATM', this%memoryPath, var_type)
       if (var_type == 'UNKNOWN') then
         call mem_allocate(this%tatm, this%ncv, 'TATM', this%memoryPath)
         do n = 1, this%ncv
           this%tatm(n) = DZERO
         end do
       end if
       call mem_allocate(this%rh, this%ncv, 'RH', this%memoryPath)
       !
       ! -- initialize
       do n = 1, this%ncv
         this%wspd(n) = DZERO
         this%tatm(n) = DZERO
         this%rh(n) = DZERO
       end do
       call this%lhf%ar()
    end if
    
    !! -- allocate character array for status
    !allocate (this%status(this%ncv))
    !!
    !! -- initialize arrays
    !do n = 1, this%ncv
    !  this%status(n) = 'ACTIVE'
    !end do
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
    if (this%inshf /= 0) then
      call this%shf%da()
      deallocate (this%shf)
    end if
    ! -- SWR (shortwave radiation heat flux)
    if (this%inswr /= 0) then
      call this%swr%da()
      deallocate (this%swr)
    end if
    ! -- LHF (latent heat flux)
    if (this%inlhf /= 0) then
      call this%lhf%da()
      deallocate (this%lhf)
    end if
    !
    ! -- Deallocate scalars
    call mem_deallocate(this%ncv)
    call mem_deallocate(this%shf_active)
    call mem_deallocate(this%swr_active)
    call mem_deallocate(this%lhf_active)
    call mem_deallocate(this%inshf)
    call mem_deallocate(this%inswr)
    call mem_deallocate(this%inlhf)
    call mem_deallocate(this%rhoa)
    call mem_deallocate(this%cpa)
    call mem_deallocate(this%cd)
    call mem_deallocate(this%wfslope)
    call mem_deallocate(this%wfint)
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
    !
    ! -- Deallocate scalars in TspAptType
    call this%NumericalPackageType%da() ! this may not work -- revisit and cleanup !!!
    
    ! -- Deallocate parent
    !call pbstbase_da(this)
  end subroutine abc_da
  
   !> @brief Read a ABC-specific option from the OPTIONS block
  !!
  !! Process a single ABC-specific option. Used when reading the OPTIONS block
  !! of the ABC package input file.
  !<
  function abc_read_option(this, keyword) result(success)
    ! -- dummy
    class(AbcType) :: this
    character(len=*), intent(in) :: keyword
    ! -- return
    logical :: success
    !
    ! -- There are no ABC-specific options, so just return false
    success = .false.
  end function abc_read_option
  
  !> @brief Sensible Heat Flux (SHF) term
  !<
  subroutine abc_shf_term(this, ientry, n1, n2, rrate, rhsval, hcofval)
    ! -- dummy
    class(AbcType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout), optional :: rrate
    real(DP), intent(inout), optional :: rhsval
    real(DP), intent(inout), optional :: hcofval
    ! -- local
    real(DP) :: sensheat
    real(DP) :: strmtemp
    integer(I4B) :: auxpos
    real(DP) :: sa !< surface area of stream reach, different than wetted area
    !
    !n1 = this%flowbudptr%budterm(this%idxbudevap)%id1(ientry)
    !! -- For now, there is only 1 aux variable under 'EVAPORATION'
    !auxpos = this%flowbudptr%budterm(this%idxbudevap)%naux
    !sa = this%flowbudptr%budterm(this%idxbudevap)%auxvar(auxpos, ientry)
    !!
    !strmtemp = this%xnewpak(n1)
    !call this%shf%shf_cq(n1, strmtemp, sensheat)
    !!
    !if (present(rrate)) rrate = sensheat * sa
    !if (present(rhsval)) rhsval = -rrate
    !if (present(hcofval)) hcofval = DZERO
  end subroutine abc_shf_term

  !> @brief Shortwave Radiation (SWR) term
  !<
  subroutine abc_swr_term(this, ientry, n1, n2, rrate, rhsval, hcofval)
    ! -- dummy
    class(AbcType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout), optional :: rrate
    real(DP), intent(inout), optional :: rhsval
    real(DP), intent(inout), optional :: hcofval
    ! -- local
    real(DP) :: shrtwvheat
    real(DP) :: strmtemp
    integer(I4B) :: auxpos
    real(DP) :: sa !< surface area of stream reach, different than wetted area
    !
    !n1 = this%flowbudptr%budterm(this%idxbudevap)%id1(ientry)
    !! -- For now, there is only 1 aux variable under 'EVAPORATION'
    !auxpos = this%flowbudptr%budterm(this%idxbudevap)%naux
    !sa = this%flowbudptr%budterm(this%idxbudevap)%auxvar(auxpos, ientry)
    !!
    !strmtemp = this%xnewpak(n1)
    !call this%swr%swr_cq(n1, shrtwvheat)
    !!
    !if (present(rrate)) rrate = shrtwvheat * sa
    !if (present(rhsval)) rhsval = -rrate
    !if (present(hcofval)) hcofval = DZERO
  end subroutine abc_swr_term
  
  !> @brief Latent Heat Flux (LHF) term
  !<
  subroutine abc_lhf_term(this, ientry, n1, n2, rrate, rhsval, hcofval)
    ! -- dummy
    class(AbcType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout), optional :: rrate
    real(DP), intent(inout), optional :: rhsval
    real(DP), intent(inout), optional :: hcofval
    ! -- local
    real(DP) :: latheat
    real(DP) :: strmtemp
    integer(I4B) :: auxpos
    real(DP) :: sa !< surface area of stream reach, different than wetted area
    !
  end subroutine abc_lhf_term
  
   !> @brief Observations
  !!
  !! Store the observation type supported by the APT package and override
  !! BndType%bnd_df_obs
  !<
  subroutine abc_df_obs(this)
    ! -- modules
    ! -- dummy
    class(AbcType) :: this
    ! -- local
    integer(I4B) :: indx
    
    ! -- Store obs type and assign procedure pointer
    !    for sens-heat-flux observation type.
    !call this%obs%StoreObsType('shf', .true., indx)
    !this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !!
    !! -- Store obs type and assign procedure pointer
    !!    for shortwave-radiation-flux observation type.
    !call this%obs%StoreObsType('swr', .true., indx)
    !this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
  end subroutine abc_df_obs
  
  !> @brief Process package specific obs
  !!
  !! Method to process specific observations for this package.
  !<
  subroutine abc_rp_obs(this, obsrv, found)
    ! -- dummy
    class(AbcType), intent(inout) :: this !< package class
    type(ObserveType), intent(inout) :: obsrv !< observation object
    logical, intent(inout) :: found !< indicate whether observation was found
    ! -- local
    !
    found = .true.
    select case (obsrv%ObsTypeId)
    case ('SHF')
!      call this%rp_obs_byfeature(obsrv)
    case ('SWR')
!      call this%rp_obs_byfeature(obsrv)
    case ('LHF')
!      call this%rp_obs_byfeature(obsrv)
    case default
      found = .false.
    end select
  end subroutine abc_rp_obs
  
   !> @brief Calculate observation value and pass it back to APT
  !<
  subroutine abc_bd_obs(this, obstypeid, jj, v, found)
    ! -- dummy
    class(AbcType), intent(inout) :: this
    character(len=*), intent(in) :: obstypeid
    real(DP), intent(inout) :: v
    integer(I4B), intent(in) :: jj
    logical, intent(inout) :: found
    ! -- local
    integer(I4B) :: n1, n2
    !
    found = .true.
    !select case (obstypeid)
    !case ('SHF')
    !  if (this%iboundpak(jj) /= 0) then
    !    call this%sfe_shf_term(jj, n1, n2, v)
    !  end if
    !case ('SWR')
    !  if (this%iboundpak(jj) /= 0) then
    !    call this%sfe_swr_term(jj, n1, n2, v)
    !  end if
    !case default
    !  found = .false.
    !end select
  end subroutine abc_bd_obs

  !> @brief Calculate Atmospheric-Stream Heat Flux
  !!
  !! Calculate and return the atmospheric heat flux for one reach
  !<
  subroutine abc_cq(this, ifno, tstrm, abcflx)
    ! -- dummy
    class(AbcType), intent(inout) :: this
    integer(I4B), intent(in) :: ifno !< stream reach integer id
    real(DP), intent(in) :: tstrm !< temperature of the stream reach
    real(DP), intent(inout) :: abcflx !< calculated atmospheric boundary flux amount
    ! -- local
    real(DP) :: shflx
    real(DP) :: swrflx
    real(DP) :: lhflx
    !
    ! -- calculate sensible heat flux using HGS equation
    call this%shf%shf_cq(ifno, tstrm, shflx)
    !
    ! -- calculate shortwave radiation using HGS equation
    call this%swr%swr_cq(ifno, swrflx)
    !
    ! -- calculate latent heat flux using Dalton-like mass transfer equation
    call this%lhf%lhf_cq(ifno, tstrm, this%gwecommon%gwerhow, lhflx)
    
    abcflx = shflx + swrflx - lhflx
  end subroutine abc_cq

  !> @brief Set the stress period attributes based on the keyword
  !<
  subroutine abc_set_stressperiod(this, itemno) !, keyword, found) MAKING LIKE PBST SET STRESSPERIOD
    ! -- module
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(AbcType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    !character(len=*), intent(in) :: keyword
    !logical, intent(inout) :: found
    ! -- local
    character(len=LINELENGTH) :: text
    character(len=LINELENGTH) :: keyword
    !logical(LGP) :: found
    integer(I4B) :: ierr
    integer(I4B) :: jj
    real(DP), pointer :: bndElem => null()
    !
    ! <wspd> WIND SPEED
    ! <tatm> TEMPERATURE OF THE ATMOSPHERE
    ! <shd> SHADE
    ! <swrefl> REFLECTANCE OF SHORTWAVE RADIATION OFF WATER SURFACE
    ! <solr> SOLAR RADIATION
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
    ! -- formats
    ierr = 0
    if (itemno < 1 .or. itemno > this%ncv) then
      write (errmsg, '(a,1x,i6,1x,a,1x,i6)') &
        'Featureno ', itemno, 'must be > 0 and <= ', this%ncv
      call store_error(errmsg)
      ierr = 1
    end if
  end function abc_check_valid

end module AbcModule
