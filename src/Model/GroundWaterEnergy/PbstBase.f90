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
                             LENVARNAME, DCTOK, DHUNDRED, DONESEVENTH
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
    procedure :: pbst_options
    procedure :: subpck_set_stressperiod
    procedure, private :: pbstbase_allocate_scalars
    procedure :: da => pbstbase_da
    procedure :: pbst_check_valid
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
    epsa = 1.24_DP * (eatm / tatm)**(DONESEVENTH) * atmc
  end function epsa

  !> @brief Calculate shade-altered emissivity
  !!
  !! Riparian shade can alter the above-channel emissivity which therefore
  !! affects the net longwave heat flux.  Emissivity above the stream channel
  !! is the combination of shade-weighted average of clear sky  atmospheric
  !! emissivity (epsa) and riparian canopy emissivity
  !<
  function epss(this, shd, epsa, epsr)
    ! -- dummy
    class(PbstBaseType) :: this
    real(DP) :: shd !< percent shading, expressed as a fraction
    real(DP) :: epsa !< atmospheric emissivity
    real(DP) :: epsr !< riparian canopy emissivity
    ! -- return
    real(DP) :: epss !< shade weighted atmospheric emissivity
    !
    epss = (1 - shd) * epsa + shd * epsr
  end function epss

end module PbstBaseModule
