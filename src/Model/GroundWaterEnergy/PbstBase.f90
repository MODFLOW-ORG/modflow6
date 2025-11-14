!> @brief This module contains common process-based stream temperature functionality
!!
!! This module contains methods for implementing functionality associated with
!! heat fluxes to a stream reach.  Four sources of thermal energy commonly
!! accounted for in process-based stream temperature modeling include short-
!! wave radiation, long-wave radiation, sensible heat flux, and latent heat
!! flux.
!<
module PbstBaseModule
  use ConstantsModule, only: LINELENGTH, MAXCHARLEN, DZERO, LGP, &
                             LENPACKAGENAME, TABLEFT, TABCENTER, &
                             LENVARNAME, DCTOK, DHUNDRED
  use KindModule, only: I4B, DP
  use NumericalPackageModule, only: NumericalPackageType
  use SimModule, only: count_errors, store_error, ustop
  use SimVariablesModule, only: errmsg
  use TdisModule, only: kper, nper, kstp
  use TimeSeriesLinkModule, only: TimeSeriesLinkType
  use TimeSeriesManagerModule, only: TimeSeriesManagerType, tsmanager_cr
  use TableModule, only: TableType, table_cr

  implicit none

  private

  public :: PbstBaseType
  public :: pbstbase_da

  character(len=LENVARNAME) :: text = '         PBST'

  type, extends(NumericalPackageType) :: PbstBaseType

    integer(I4B), pointer :: ncv => null() !< number of control volumes
    integer(I4B), dimension(:), pointer, contiguous :: iboundpbst => null() !< package ibound
    logical, pointer, public :: active => null() !< logical indicating if a sensible heat flux object is active
    character(len=8), dimension(:), pointer, contiguous :: status => null() !< active, inactive, constant
    character(len=LENPACKAGENAME) :: text = '' !< text string for package transport term
    character(len=LINELENGTH), pointer, public :: inputFilename => null() !< a particular pbst input file name, could be for sensible heat flux or latent heat flux subpackages, for example
    type(TimeSeriesManagerType), pointer :: tsmanager => null()
    !
    ! -- table objects
    type(TableType), pointer :: inputtab => null() !< input table object

  contains

    procedure :: init
    procedure :: pbst_ar
    procedure :: rp
    procedure :: pbst_options
    procedure :: pbst_set_stressperiod
    procedure :: subpck_set_stressperiod
    procedure, private :: pbstbase_allocate_scalars
    procedure :: da => pbstbase_da
    procedure :: pbst_check_valid
    procedure, public :: eair !< function for calculating vapor pressure for a specified temperature
    procedure, public :: eatm !< funciton for calculating ambient vapor pressure of the atmosphere
    procedure, public :: epsa !< function for calculating atmospheric emissivity
    procedure, public :: epss !< function for calculating shade-weighted atmospheric emissivity

  end type PbstBaseType

contains

  !> @brief Initialize the PbstBaseType object
  !!
  !! Allocate and initialize data members of the object.
  !<
  subroutine init(this, name_model, pakname, ftype, inunit, iout, ncv)
    ! -- dummy
    class(PbstBaseType) :: this
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: pakname
    character(len=*), intent(in) :: ftype
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    integer(I4B), target, intent(in) :: ncv
    !
    call this%set_names(1, name_model, pakname, ftype)
    call this%pbstbase_allocate_scalars()
    this%inunit = inunit
    this%iout = iout
    this%ncv => ncv
    call this%parser%Initialize(this%inunit, this%iout)
  end subroutine init

  !> @brief Allocate and read
  !!
  !! Method to allocate and read static data for the PBST utility packages
  !<
  subroutine pbst_ar(this)
    ! -- dummy
    class(PbstBaseType) :: this !< ShfType object
    !
    ! -- will be overridden
  end subroutine pbst_ar

  !> @brief PaBST read and prepare for setting stress period information
  !<
  subroutine rp(this)
    ! -- module
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    use TdisModule, only: kper, nper
    ! -- dummy
    class(PbstBaseType) :: this !< PbstBaseType object
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
    ! -- set ionper to the stress period number for which a new block of data
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
          ! -- end of file found; data applies for remainder of simulation.
          this%ionper = nper + 1
        else
          ! -- found invalid block
          call this%parser%GetCurrentLine(line)
          write (errmsg, fmtblkerr) adjustl(trim(line))
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
      end if
    end if
    !
    ! --read data if ionper == kper
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
        call this%pbst_set_stressperiod(itemno)
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
  end subroutine rp

  !> @brief pbst_set_stressperiod()
  !!
  !! To be overridden by Pbst sub-packages (utilities)
  !<
  subroutine pbst_set_stressperiod(this, itemno)
    ! -- dummy
    class(PbstBaseType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    ! -- local
    integer(I4B) :: ierr
    character(len=LINELENGTH) :: text
    character(len=LINELENGTH) :: keyword
    logical(LGP) :: found
    !
    ! -- read line
    call this%parser%GetStringCaps(keyword)
    select case (keyword)
    case ('STATUS')
      ierr = this%pbst_check_valid(itemno)
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
    case default
      !
      ! -- call the specific package to deal with parameters specific to the
      !    process
      call this%subpck_set_stressperiod(itemno, keyword, found)
      ! -- terminate with error if data not valid
      if (.not. found) then
        write (errmsg, '(2a)') &
          'Unknown '//trim(adjustl(this%text))//' data keyword: ', &
          trim(keyword)//'.'
        call store_error(errmsg)
      end if
    end select
    !
    ! -- terminate if any errors were detected
999 if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine pbst_set_stressperiod

  !> @brief pbst_set_stressperiod()
  !!
  !! To be overridden by Pbst sub-packages
  !<
  subroutine subpck_set_stressperiod(this, itemno, keyword, found)
    ! -- dummy
    class(PbstBaseType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    character(len=*), intent(in) :: keyword
    logical, intent(inout) :: found
    ! -- to be overwritten by pbst subpackages (or "utilities")
  end subroutine subpck_set_stressperiod

  !> @brief Read additional options for sub-package
  !!
  !! Read additional options for the SFE boundary package. This method should
  !! be overridden by option-processing routine that is in addition to the
  !! base options available for all PbstBase packages.
  !<
  subroutine pbst_options(this, option, found)
    ! -- dummy
    class(PbstBaseType), intent(inout) :: this !< PbstBaseType object
    character(len=*), intent(inout) :: option !< option keyword string
    logical(LGP), intent(inout) :: found !< boolean indicating if the option was found
    !
    ! -- return with found = .false.
    found = .false.
  end subroutine pbst_options

  !> @brief Allocate scalar variables
  !!
  !! Allocate scalar data members of the object.
  !<
  subroutine pbstbase_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(PbstBaseType) :: this
    !
    allocate (this%active)
    allocate (this%inputFilename)
    !
    ! -- initialize
    this%active = .false.
    this%inputFilename = ''
    !
    ! -- call standard NumericalPackageType allocate scalars
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- allocate
    call mem_allocate(this%ncv, 'NCV', this%memoryPath)
    !
    ! -- initialize
    this%ncv = 0
    !
    ! -- allocate time series manager
    allocate (this%tsmanager)
  end subroutine pbstbase_allocate_scalars

  !> @brief Deallocate package memory
  !!
  !! Deallocate package scalars and arrays.
  !<
  subroutine pbstbase_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(PbstBaseType) :: this
    !
    deallocate (this%active)
    deallocate (this%inputFilename)
    !
    ! -- deallocate time series manager
    deallocate (this%tsmanager)
    !
    ! -- deallocate scalars
    call mem_deallocate(this%ncv)
    !
    ! -- deallocate arrays
    call mem_deallocate(this%iboundpbst)
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()
    !
    ! -- input table object
    if (associated(this%inputtab)) then
      call this%inputtab%table_da()
      deallocate (this%inputtab)
      nullify (this%inputtab)
    end if
  end subroutine pbstbase_da

  !> @brief Process-based stream temperature transport (or utility) routine
  !!
  !! Determine if a valid feature number has been specified.
  !<
  function pbst_check_valid(this, itemno) result(ierr)
    ! -- return
    integer(I4B) :: ierr
    ! -- dummy
    class(PbstBaseType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    !
    ! -- initialize
    ierr = 0
    !
    if (itemno < 1 .or. itemno > this%ncv) then
      write (errmsg, '(a,1x,i6,1x,a,1x,i6)') &
        'Featureno ', itemno, 'must be > 0 and <= ', this%ncv
      call store_error(errmsg)
      ierr = 1
    end if
  end function pbst_check_valid

  !> @brief Calculate saturation vapor pressure
  !!
  !! Calculate vapor pressure for a passed-in temperature.  Calculated using
  !! the Magnus-Tetens formula (Buck, 1981)
  !<
  function eair(this, tempc, tempk)
    ! -- dummy
    class(PbstBaseType) :: this
    real(DP) :: tempc !< temperature in Celcius
    real(DP) :: tempk !< temperature in Kelvin
    ! -- return
    real(DP) :: eair
    !
    eair = 6.1275_DP * exp(17.2693882_DP * (tempc / (tempk - 35.86_DP)))
  end function eair

  !> @brief Calculate ambient vapor pressure
  !!
  !! Calculate ambient vapor pressure of the atmosphere as a function of
  !! relative humidity and saturation vapor pressure of the atmosphere
  !<
  function eatm(this, rh, eair)
    ! -- dummy
    class(PbstBaseType) :: this
    real(DP) :: rh
    real(DP) :: eair
    ! -- return
    real(DP) :: eatm
    !
    eatm = rh / DHUNDRED * eair
  end function

  !> @brief Calculate atmospheric emissivity
  !!
  !! Atmospheric emissivity is highly dependent upon ambient vapor
  !! concentration of the air (Campbell and Norman 1989). Estimated
  !! clear-sky emissivity is calculated using (Brutsaert and Jirka 1984)
  !<
  function epsa(this, eatm, tatm, atmc)
    ! -- dummy
    class(PbstBaseType) :: this
    real(DP) :: atmc
    real(DP) :: eatm
    real(DP) :: tatm
    ! -- return
    real(DP) :: epsa !< atmospheric emissivity
    !
    epsa = 1.24_DP * (eatm / tatm)**(1 / 7) * atmc
  end function epsa

  !> @brief Calculate shade-altered emissivity
  !!
  !! Riparian shade can alter the above-channel emissivity which therefore
  !! affects the net longwave heat flux.  Emissivity above the stream channel
  !! is the combination of shade-weighted average of clear sky  atmospheric
  !! emissivity (epsa) and riparian canopy emissivity
  function epss(this, shd, epsa, epsr)
    ! -- dummy
    class(PbstBaseType) :: this
    real(DP) :: shd
    real(DP) :: epsa !< atmospheric emissivity
    real(DP) :: epsr !< riparian canopy emissivity
    ! -- return
    real(DP) :: epss !< shade weighted atmospheric emissivity
    !
    epss = (1 - shd) * epsa + shd * epsr
  end function epss

end module PbstBaseModule
