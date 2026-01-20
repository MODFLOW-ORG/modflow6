!> @brief This module contains the derived type ObsType
!!
!! This module defines type ObsType, which is the highest-level
!! derived type for implementing observations. All objects derived from
!! NumericalModelType or BndType already contain an ObsType member.
!!
!! Examples:
!!   NumericalModelType.obs
!!   BndType.obs
!!
!! Similarly, an ObsType member could be added to, say,
!! NumericalExchangeType or any other type that has DF, AR, RP, AD, BD, and OT
!! routines.
!!
!! IMPLEMENTATION OF OBSERVATIONS IN A MODEL OR PACKAGE
!!
!! For simple boundary packages like RIV and DRN, only steps 1-6 are
!! needed. For models and advanced packages like MAW and SFR, additional
!! steps are needed.
!!
!! 1. (package only) Override BndType.bnd_obs_supported to return true.
!!    bnd_obs_supported is called from various places in code.
!!
!! 2. (optional) Write a subroutine that implements abstract interface
!!    ObserveModule.ProcessIdSub. (Not needed if IDstring, which identifies
!!    location in model to be observed, is either a single node number or
!!    a single {lay, row, col} set of indices).
!!
!!    Examples:
!!    gwf_process_head_drawdown_obs_id, gwf_process_intercell_obs_id
!!
!!    A package can allow IDstring to be a boundary name.
!!    Example: ObsModule.DefaultObsIdProcessor
!!
!! 3. Override BndType.bnd_df_obs() to define string(s) to be
!!    recognized as observation type(s) and (optional) assign ProcessIdPtr
!!    (not needed  if IDstring is either a node number or a {lay, row, col}
!!    set of indices).
!!
!!    Examples: gwf_df_obs, drn_df_obs
!!
!!    When boundary names are allowed and developer wants simulated value
!!    to be cumulative (flow, for example) if user specifies multiple
!!    boundaries with the same BOUNDNAME, in bnd_df_obs call to
!!    ObsPackage.StoreObsType, provide cumulative argument as true.
!!    Otherwise, simulated values are not cumulative.
!!
!! 4. In DF routine: Call bnd_df_obs
!!
!! 5. In AR routine: Call ObsType.obs_ar. This reads the OBS input
!!    file.
!!    Example (gwf_ar): call this%obs%obs_ar()
!!    Example (lak_ar): call this%obs%obs_ar()
!!
!! 6. Override BndType.bnd_rp_obs for any package that needs to
!!    check user input or process observation input in any special way.
!!    If no special processing is needed, BndType.bnd_rp_obs can
!!    be used.  This routine also expands the ObserveType%indxbnds array for
!!    each observation in a package. ObserveType%indxbnds is used to sum
!!    simulated values from multiple boundaries when BOUNDNAMES is used.
!!    Equivalent routine may or may not be needed for model observations.
!!    If needed, call it from bottom of RP routine.
!!
!!    Examples:
!!        BndType.bnd_rp_obs, which is called from gwf_rp
!!
!! 7. In AD routine: Call ObsType.obs_ad
!!    Example: gwf_ad
!!
!! 8. Write a *_bd_obs routine. This is the routine that actually
!!    calculates the simulated value for each observation type supported
!!    by the model/package.  Call *_bd_obs from the bottom of the
!!    _bd routine.
!!     *_bd_obs needs to:
!!         Call ObsType.obs_bd_clear
!!         For each observation:
!!              Calculate the simulated value
!!              Call ObsType.SaveOneSimval
!!     Examples: gwf_bd_obs, maw_bd_obs, lak_bd_obs
!!
!! 9. In BD routine:
!!     Call BndType.bnd_bd_obs
!!     Examples: BndType.bnd_bd calls bnd_bd_obs
!!               GwfModelType.gwf_bd calls gwf_bd_obs
!!               MawType.maw_bd calls maw_bd_obs
!!               LakType.lak_bd calls lak_bd_obs
!!
!! 10. Ensure that ObsType.obs_ot is called. For packages, obs_ot is called
!!     from the model _ot procedure.  The model _ot procedure should also call
!!     obs_ot for its own observations.  Do not call obs_ot from a package _ot
!!     procedure because the package _ot procedure may not be called, depending
!!     on Output Control settings (ibudfl).
!!
!!     Note: BndType.bnd_ot_obs calls:
!!               ObsType.obs_ot
!!
!!     Note: ObsType.obs_ot calls:
!!               store_all_simvals
!!               write_continuous_simvals
!!               obsOutputList.WriteOutputLines
!!
!! BINARY OUTPUT:
!!
!! When observation-output files are written, the user has the option to have
!! output written to a binary file.  Binary obs output files start with a
!! 100-byte header structured as follows:
!!
!! bytes 1-4   (ascii): Observation type contained in file; options are:
!!                        "cont" -- Continuous observations
!! byte 5: blank
!! bytes 6-11  (ascii): Precision of all floating-point values; options are:
!!                        "single" -- Single precision
!!                        "double" -- Double precision
!! bytes 12-15 (ascii): LENOBSNAME (integer; length of observation names,
!!                                  in bytes)
!! bytes 16-100: blank
!!
!! IN A FILE OF CONTINUOUS OBSERVATIONS:
!!
!! The 100-byte header is followed by:
!! NOBS (4-byte integer) -- Number of observations.
!! NOBS repetitions of OBSNAME (ascii, LENOBSNAME bytes each).
!! Any number of repetitions of:
!! TIME SIMVAL-1 SIMVAL-2 ... SIMVAL-NOBS  (floating point)
!!
!<
module ObsModule

  use KindModule, only: DP, I4B
  use ArrayHandlersModule, only: ExpandArray
  use BaseDisModule, only: DisBaseType
  use BlockParserModule, only: BlockParserType
  use ConstantsModule, only: LENBIGLINE, LENFTYPE, LENOBSNAME, &
                             LENOBSTYPE, LENPACKAGENAME, LENBOUNDNAME, &
                             LINELENGTH, NAMEDBOUNDFLAG, MAXCHARLEN, &
                             MAXOBSTYPES, LENHUGELINE, DNODATA, &
                             TABLEFT
  use TableModule, only: TableType, table_cr
  use InputOutputModule, only: UPCASE, openfile, GetUnit, GetFileFromPath
  use ListModule, only: ListType
  use ObsContainerModule, only: ObsContainerType
  use ObserveModule, only: ConstructObservation, ObsDataType, &
                           ObserveType, GetObsFromList, &
                           AddObsToList
  use ObsOutputListModule, only: ObsOutputListType
  use ObsOutputModule, only: ObsOutputType
  use ObsUtilityModule, only: write_fmtd_obs, write_unfmtd_obs
  use OpenSpecModule, only: ACCESS, FORM
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, store_error_unit
  use TdisModule, only: totim

  implicit none

  private
  public :: ObsType, DefaultObsIdProcessor, obs_cr

  type :: ObsType
    ! -- Public members
    integer(I4B), public :: iout = 0 !< model list file unit
    integer(I4B), public :: npakobs = 0 !< number of observations
    integer(I4B), pointer, public :: inUnitObs => null() !< observation input file unit
    character(len=LINELENGTH), pointer, public :: inputFilename => null() !< observation input file name
    character(len=2*LENPACKAGENAME + 4), public :: pkgName = '' !< package name
    character(len=LENFTYPE), public :: filtyp = '' !< package file type
    logical, pointer, public :: active => null() !> logical indicating if a observation is active
    type(ObsContainerType), dimension(:), pointer, public :: pakobs => null() !< package observations
    type(ObsDataType), dimension(:), pointer, public :: obsData => null() !< observation data
    ! -- Private members
    integer(I4B), private :: iprecision = 2 ! 2=double; 1=single
    integer(I4B), private :: idigits = 0
    character(len=LINELENGTH), private :: outputFilename = ''
    character(len=LINELENGTH), private :: blockTypeFound = ''
    character(len=20), private :: obsfmtcont = ''
    logical, private :: echo = .false.
    logical, private :: more
    type(ListType), private :: obsList
    type(ObsOutputListType), pointer, private :: obsOutputList => null()
    class(DisBaseType), pointer, private :: dis => null()
    type(BlockParserType), private :: parser
    !
    ! -- table object
    type(TableType), pointer :: obstab => null()
  contains
    ! -- Public procedures
    procedure, public :: obs_df
    procedure, public :: obs_ar
    procedure, public :: obs_ad
    procedure, public :: obs_bd_clear
    procedure, public :: obs_ot
    procedure, public :: obs_da
    procedure, public :: SaveOneSimval
    procedure, public :: StoreObsType
    procedure, public :: allocate_scalars
    ! -- Private procedures
    procedure, private :: build_headers
    procedure, private :: define_fmts
    procedure, private :: get_num
    procedure, private :: get_obs
    procedure, private :: get_obs_array
    procedure, private :: get_obs_datum
    procedure, private :: obs_ar1
    procedure, private :: obs_ar2
    procedure, private :: set_obs_array
    procedure, private :: read_observations
    procedure, private :: read_obs_blocks
    procedure, private :: read_obs_options
    procedure, private :: write_obs_simvals
  end type ObsType

contains

  ! Non-type-bound procedures

  !> @ brief Create a new ObsType object
  !!
  !!  Subroutine to create a new ObsType object. Soubroutine
  !!
  !!  - creates object
  !!  - allocates pointer
  !!  - initializes values
  !!
  !<
  subroutine obs_cr(obs, inobs)
    ! -- dummy
    type(ObsType), pointer, intent(out) :: obs !< observation ObsType
    integer(I4B), pointer, intent(in) :: inobs !< observation input file unit
    !
    allocate (obs)
    call obs%allocate_scalars()
    obs%inUnitObs => inobs
  end subroutine obs_cr

  !> @ brief Process IDstring provided for each observation
  !!
  !!  Subroutine to process the IDstring provided for each observation. The
  !!  IDstring identifies the location in the model of the node(s) or feature(s)
  !!  where the simulated value is to be extracted and recorded. Subroutine
  !!
  !!  - interprets the IDstring
  !!  - stores the location of interest in the ObserveType object that
  !!    contains information about the observation
  !!
  !<
  subroutine DefaultObsIdProcessor(obsrv, dis, inunitobs, iout)
    ! -- dummy
    type(ObserveType), intent(inout) :: obsrv !< observation ObserveType
    class(DisBaseType), intent(in) :: dis !< discretization object
    integer(I4B), intent(in) :: inunitobs !< observation input file unit
    integer(I4B), intent(in) :: iout !< model list file
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: icol, istart, istop
    character(len=LINELENGTH) :: string
    logical :: flag_string
    !
    ! -- Initialize variables
    string = obsrv%IDstring
    icol = 1
    flag_string = .true. ! Allow string to contain a boundary name
    !
    n = dis%noder_from_string(icol, istart, istop, inunitobs, &
                              iout, string, flag_string)
    !
    if (n > 0) then
      obsrv%NodeNumber = n
    elseif (n == -2) then
      ! Integer can't be read from string; it's presumed to be a boundary
      ! name (already converted to uppercase)
      obsrv%FeatureName = string(istart:istop)
      ! -- Observation may require summing rates from multiple boundaries,
      !    so assign NodeNumber as a value that indicates observation
      !    is for a named boundary or group of boundaries.
      obsrv%NodeNumber = NAMEDBOUNDFLAG
    else
      errmsg = 'Error reading data from ID string'
      call store_error(errmsg)
      call store_error_unit(inunitobs)
    end if
  end subroutine DefaultObsIdProcessor

  ! Type-bound public procedures

  !> @ brief Define some members of an ObsType object
  !!
  !!  Subroutine to define some members of an ObsType object.
  !!
  !<
  subroutine obs_df(this, iout, pkgname, filtyp, dis)
    ! -- dummy
    class(ObsType), intent(inout) :: this
    integer(I4B), intent(in) :: iout !< model list file unit
    character(len=*), intent(in) :: pkgname !< package name
    character(len=*), intent(in) :: filtyp !< package file type
    class(DisBaseType), pointer :: dis !< discretization object
    !
    this%iout = iout
    this%pkgName = pkgname
    this%filtyp = filtyp
    this%dis => dis
    !
    ! -- Initialize block parser
    call this%parser%Initialize(this%inUnitObs, this%iout)
  end subroutine obs_df

  !> @ brief Allocate and read package observations
  !!
  !!  Subroutine to allocate and read observations for a package. Subroutine
  !!
  !!  - reads OPTIONS block of OBS input file
  !!  - reads CONTINUOUS blocks of OBS input file
  !!
  !<
  subroutine obs_ar(this)
    ! -- dummy
    class(ObsType) :: this
    !
    call this%obs_ar1(this%pkgName)
    if (this%active) then
      call this%obs_ar2(this%dis)
    end if
  end subroutine obs_ar

  !> @ brief Advance package observations
  !!
  !!  Subroutine to advance each package observations by resetting the
  !!  "current" value.
  !!
  !<
  subroutine obs_ad(this)
    ! -- dummy
    class(ObsType) :: this
    ! -- local
    integer(I4B) :: i, n
    class(ObserveType), pointer :: obsrv => null()
    !
    n = this%get_num()
    do i = 1, n
      obsrv => this%get_obs(i)
      call obsrv%ResetCurrentValue()
    end do
  end subroutine obs_ad

  !> @ brief Clear observation output lines
  !!
  !!  Subroutine to clear output lines in preparation for new rows of
  !!  continuous observations.
  !!
  !<
  subroutine obs_bd_clear(this)
    ! -- dummy
    class(ObsType), target :: this
    !
    call this%obsOutputList%ResetAllObsEmptyLines()
  end subroutine obs_bd_clear

  !> @ brief Output observation data
  !!
  !!  Subroutine to output observation data. Subroutine
  !!
  !!  - stores each simulated value into its ObserveType object
  !!  - writes each simulated value to it ObsOutputList object
  !!  _ writes contents of ObsOutputList to output file
  !!
  !! This procedure should NOT be called from a package's _ot procedure
  !! because the package _ot procedure may not be called every time step.
  !!
  !<
  subroutine obs_ot(this)
    ! -- dummy
    class(ObsType), intent(inout) :: this
    !
    if (this%npakobs > 0) then
      call this%write_obs_simvals()
      call this%obsOutputList%WriteAllObsLineReturns()
    end if
  end subroutine obs_ot

  !> @ brief Deallocate observation data
  !!
  !!  Subroutine to deallocate observation data.
  !!
  !<
  subroutine obs_da(this)
    ! -- dummy
    class(ObsType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i
    class(ObserveType), pointer :: obsrv => null()
    !
    deallocate (this%active)
    deallocate (this%inputFilename)
    deallocate (this%obsData)
    !
    ! -- observation table object
    if (associated(this%obstab)) then
      call this%obstab%table_da()
      deallocate (this%obstab)
      nullify (this%obstab)
    end if
    !
    ! -- deallocate pakobs components and pakobs
    if (associated(this%pakobs)) then
      do i = 1, this%npakobs
        obsrv => this%pakobs(i)%obsrv
        call obsrv%da()
        deallocate (obsrv)
        nullify (this%pakobs(i)%obsrv)
      end do
      deallocate (this%pakobs)
    end if
    !
    ! -- deallocate obsOutputList
    call this%obsOutputList%DeallocObsOutputList()
    deallocate (this%obsOutputList)
    !
    ! -- deallocate obslist
    call this%obslist%Clear()
    !
    ! -- nullify
    nullify (this%inUnitObs)
  end subroutine obs_da

  !> @ brief Save a simulated value
  !!
  !!  Subroutine to save or accumulate a simulated value to its ObserveType
  !!  object.
  !!
  !<
  subroutine SaveOneSimval(this, obsrv, simval)
    ! -- dummy
    class(ObsType) :: this
    class(ObserveType), intent(inout) :: obsrv !< observation ObserveType
    real(DP), intent(in) :: simval !< simulated value
    ! -- local
    character(len=LENOBSTYPE) :: obsTypeID
    type(ObsDataType), pointer :: obsDatum => null()
    !
    ! -- initialize variables
    obsTypeID = obsrv%ObsTypeId
    obsDatum => this%get_obs_datum(obsTypeID)
    !
    ! -- save current simulation time
    obsrv%CurrentTimeStepEndTime = totim
    !
    ! -- assign or accumulate simulated value
    if (obsDatum%Cumulative .and. simval /= DNODATA) then
      obsrv%CurrentTimeStepEndValue = obsrv%CurrentTimeStepEndValue + simval
    else
      obsrv%CurrentTimeStepEndValue = simval
    end if
  end subroutine SaveOneSimval

  !> @ brief Store observation type
  !!
  !!  Subroutine to store type name and related information for an
  !!  observation type that belongs to a package or model in the
  !!  obsData array.
  !!
  !<
  subroutine StoreObsType(this, obsrvType, cumulative, indx)
    ! -- dummy
    class(ObsType), intent(inout) :: this
    character(len=*), intent(in) :: obsrvType !< observation type
    ! cumulative:  Accumulate simulated values for multiple boundaries
    logical, intent(in) :: cumulative !< logical indicating if the observation should be accumulated
    integer(I4B), intent(out) :: indx !< observation index
    ! -- local
    integer(I4B) :: i
    character(len=LENOBSTYPE) :: obsTypeUpper
    character(len=100) :: msg
    !
    ! -- Ensure that obsrvType is not blank
    if (obsrvType == '') then
      msg = 'Programmer error: Invalid argument in store_obs_type.'
      call store_error(msg, terminate=.TRUE.)
    end if
    !
    ! -- Find first unused element
    indx = -1
    do i = 1, MAXOBSTYPES
      if (this%obsData(i)%ObsTypeID /= '') cycle
      indx = i
      exit
    end do
    !
    ! -- Ensure that array size is not exceeded
    if (indx == -1) then
      msg = 'Size of obsData array is insufficient; ' &
            //'need to increase MAXOBSTYPES.'
      call store_error(msg)
      call store_error_unit(this%inUnitObs)
    end if
    !
    ! -- Convert character argument to upper case
    obsTypeUpper = obsrvType
    call upcase(obsTypeUpper)
    !
    ! -- Assign members
    this%obsData(indx)%ObsTypeID = obsTypeUpper
    this%obsData(indx)%Cumulative = cumulative
  end subroutine StoreObsType

  ! Type-bound private procedures

  !> @ brief Allocate observation scalars
  !!
  !!  Subroutine to allocate and initialize memory for non-allocatable
  !   members (scalars).
  !!
  !<
  subroutine allocate_scalars(this)
    ! -- dummy
    class(ObsType) :: this
    !
    allocate (this%active)
    allocate (this%inputFilename)
    allocate (this%obsOutputList)
    allocate (this%obsData(MAXOBSTYPES))
    !
    ! -- Initialize
    this%active = .false.
    this%inputFilename = ''
  end subroutine allocate_scalars

  !> @ brief Read observation options and output formats
  !!
  !!  Subroutine to read the options block in the observation input file and
  !!  define output formats.
  !!
  !<
  subroutine obs_ar1(this, pkgname)
    ! -- dummy
    class(ObsType), intent(inout) :: this
    character(len=*), intent(in) :: pkgname !< package name
    ! -- formats
10  format(/, 'The observation utility is active for "', a, '"')
    !
    if (this%inUnitObs > 0) then
      this%active = .true.
      !
      ! -- Indicate that OBS is active
      write (this%iout, 10) trim(pkgname)
      !
      ! -- Read Options block
      call this%read_obs_options()
      !
      ! -- define output formats
      call this%define_fmts()
    end if
  end subroutine obs_ar1

  !> @ brief Call procedure provided by package
  !!
  !!  Subroutine to call procedure provided by package to interpret IDstring
  !!  and store required data.
  !!
  !<
  subroutine obs_ar2(this, dis)
    ! -- dummy
    class(ObsType), intent(inout) :: this
    class(DisBaseType) :: dis !< discretization object
    ! -- local
    integer(I4B) :: i
    type(ObsDataType), pointer :: obsDat => null()
    character(len=LENOBSTYPE) :: obsTypeID
    class(ObserveType), pointer :: obsrv => null()
    !
    call this%read_observations()
    ! -- allocate and set observation array
    call this%get_obs_array(this%npakobs, this%pakobs)
    !
    do i = 1, this%npakobs
      obsrv => this%pakobs(i)%obsrv
      ! -- Call IDstring processor procedure provided by package
      obsTypeID = obsrv%ObsTypeId
      obsDat => this%get_obs_datum(obsTypeID)
      if (associated(obsDat%ProcessIdPtr)) then
        call obsDat%ProcessIdPtr(obsrv, dis, &
                                 this%inUnitObs, this%iout)
      else
        call DefaultObsIdProcessor(obsrv, dis, &
                                   this%inUnitObs, this%iout)
      end if
    end do
    !
    if (count_errors() > 0) then
      call store_error_unit(this%inunitobs)
    end if
  end subroutine obs_ar2

  !> @ brief Read observation options block
  !!
  !!  Subroutine to read the options block in the observation input file.
  !!
  !<
  subroutine read_obs_options(this)
    ! -- dummy
    class(ObsType) :: this
    ! -- local
    integer(I4B) :: iin
    integer(I4B) :: ierr
    integer(I4B) :: localprecision
    integer(I4B) :: localdigits
    character(len=40) :: keyword
    character(len=LINELENGTH) :: fname
    type(ListType), pointer :: lineList => null()
    logical :: continueread, found, endOfBlock
    ! -- formats
10  format('No options block found in OBS input. Defaults will be used.')
40  format('Text output number of digits of precision set to: ', i2)
50  format('Text output number of digits set to internal representation (G0).')
60  format(/, 'Processing observation options:',/)
    !
    localprecision = 0
    localdigits = -1
    lineList => null()
    !
    ! -- Find and store file name
    iin = this%inUnitObs
    inquire (unit=iin, name=fname)
    this%inputFilename = fname
    !
    ! -- Read Options block
    continueread = .false.
    ierr = 0
    !
    ! -- get BEGIN line of OPTIONS block
    call this%parser%GetBlock('OPTIONS', found, ierr, &
                              supportOpenClose=.true., blockRequired=.false.)
    if (ierr /= 0) then
      ! end of file
      errmsg = 'End-of-file encountered while searching for'// &
               ' OPTIONS in OBS '// &
               'input file "'//trim(this%inputFilename)//'"'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    elseif (.not. found) then
      this%blockTypeFound = ''
      if (this%iout > 0) write (this%iout, 10)
    end if
    !
    ! -- parse OPTIONS entries
    if (found) then
      write (this%iout, 60)
      readblockoptions: do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('DIGITS')
          !
          ! -- error if digits already read
          if (localdigits /= -1) then
            errmsg = 'Error in OBS input: DIGITS has already been defined'
            call store_error(errmsg)
            exit readblockoptions
          end if
          !
          ! -- Specifies number of significant digits used writing simulated
          !    values to a text file. Default is stored digits.
          !
          ! -- Read integer value
          localdigits = this%parser%GetInteger()
          !
          ! -- Set localdigits to valid value: 0, or 2 to 16
          if (localdigits == 0) then
            write (this%iout, 50)
          else if (localdigits < 1) then
            errmsg = 'Error in OBS input: Invalid value for DIGITS option'
            call store_error(errmsg)
            exit readblockoptions
          else
            if (localdigits < 2) localdigits = 2
            if (localdigits > 16) localdigits = 16
            write (this%iout, 40) localdigits
          end if
        case ('PRINT_INPUT')
          this%echo = .true.
          write (this%iout, '(a)') 'The PRINT_INPUT option has been specified.'
        case default
          errmsg = 'Error in OBS input: Unrecognized option: '// &
                   trim(keyword)
          call store_error(errmsg)
          exit readblockoptions
        end select
      end do readblockoptions
    end if
    !
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    write (this%iout, '(1x)')
    !
    ! -- Assign type variables
    if (localprecision > 0) this%iprecision = localprecision
    if (localdigits >= 0) this%idigits = localdigits
  end subroutine read_obs_options

  !> @ brief Define observation output formats
  !!
  !!  Subroutine to define observation output formats.
  !!
  !<
  subroutine define_fmts(this)
    ! -- dummy
    class(ObsType) :: this
    ! formats
50  format('(g', i2.2, '.', i2.2, ')')
    !
    if (this%idigits == 0) then
      this%obsfmtcont = '(G0)'
    else
      write (this%obsfmtcont, 50) this%idigits + 7, this%idigits
    end if
  end subroutine define_fmts

  !> @ brief Read observations
  !!
  !!  Subroutine to read the observations from the observation input file
  !!  and build headers for the observation output files.
  !!
  !<
  subroutine read_observations(this)
    ! -- dummy
    class(ObsType) :: this
    ! -- local
    !
    ! -- Read CONTINUOUS blocks and store observations
    call this%read_obs_blocks(this%outputFilename)
    !
    ! -- build headers
    call this%build_headers()
  end subroutine read_observations

  !> @ brief Get the number of observations
  !!
  !!  Function to get the number of observationns in this ObsType object.
  !!
  !<
  function get_num(this)
    ! -- return
    integer(I4B) :: get_num !< number of observations
    ! -- dummy
    class(ObsType) :: this
    !
    get_num = this%obsList%Count()
  end function get_num

  !> @ brief Build observation headers
  !!
  !!  Subroutine to build headers for CSV-formatted and unformatted
  !!  continuous-observation output files and write them to those files.
  !!
  !!  Each formatted header will have the form: "time,obsname-1,obsname-2, ..."
  !!
  !<
  subroutine build_headers(this)
    ! -- module
    use iso_fortran_env, only: int32
    ! -- dummy
    class(ObsType), target :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: ii
    integer(I4B) :: idx
    integer(I4B) :: iu
    integer(I4B) :: num
    integer(int32) :: nobs
    character(len=4) :: clenobsname
    type(ObserveType), pointer :: obsrv => null()
    type(ObsOutputType), pointer :: obsOutput => null()
    !
    ! -- Cycle through ObsOutputList to write headers
    !    to formatted and unformatted file(s).
    idx = 1
    num = this%obsOutputList%Count()
    all_obsfiles: do i = 1, num
      obsOutput => this%obsOutputList%Get(i)
      nobs = obsOutput%nobs
      iu = obsOutput%nunit
      !
      ! -- write header information to the formatted file
      if (obsOutput%FormattedOutput) then
        write (iu, '(a)', advance='NO') 'time'
      else
        ! -- write header to unformatted file
        !    First 11 bytes are obs type and precision
        if (this%iprecision == 1) then
          ! -- single precision output
          write (iu) 'cont single'
        else if (this%iprecision == 2) then
          ! -- double precision output
          write (iu) 'cont double'
        end if
        ! -- write LENOBSNAME to bytes 12-15
        write (clenobsname, '(i4)') LENOBSNAME
        write (iu) clenobsname
        ! -- write blanks to complete 100-byte header
        do ii = 16, 100
          write (iu) ' '
        end do
        ! -- write NOBS
        write (iu) nobs
      end if
      !
      ! -- write observation name
      obsfile: do ii = 1, nobs
        obsrv => this%get_obs(idx)
        if (obsOutput%FormattedOutput) then
          write (iu, '(a,a)', advance='NO') ',', trim(obsrv%Name)
          !
          ! -- terminate the line on the last observation in file
          if (ii == nobs) then
            write (iu, '(a)', advance='YES') ''
          end if
        else
          write (iu) obsrv%Name
        end if
        idx = idx + 1
      end do obsfile
    end do all_obsfiles
  end subroutine build_headers

  !> @ brief Get an array of observations
  !!
  !!  Subroutine to get an array containing all observations in this
  !!  ObsType object.
  !!
  !<
  subroutine get_obs_array(this, nObs, obsArray)
    ! -- dummy
    class(ObsType), intent(inout) :: this
    integer(I4B), intent(out) :: nObs !< number of observations
    type(ObsContainerType), dimension(:), pointer, intent(inout) :: obsArray !< observation array
    !
    nObs = this%get_num()
    if (associated(obsArray)) deallocate (obsArray)
    allocate (obsArray(nObs))
    !
    ! set observations in obsArray
    if (nObs > 0) then
      call this%set_obs_array(nObs, obsArray)
    end if
  end subroutine get_obs_array

  !> @ brief Get an ObsDataType object
  !!
  !!  Function to get an ObsDataType object for the specified observation type.
  !!
  !<
  function get_obs_datum(this, obsTypeID) result(obsDatum)
    ! -- dummy
    class(ObsType) :: this
    character(len=*), intent(in) :: obsTypeID !< observation type
    ! -- return
    type(ObsDataType), pointer :: obsDatum !< observation ObsDataType
    ! -- local
    integer(I4B) :: i
    !
    obsDatum => null()
    do i = 1, MAXOBSTYPES
      if (this%obsData(i)%ObsTypeID == obsTypeID) then
        obsDatum => this%obsData(i)
        exit
      end if
    end do
    !
    if (.not. associated(obsDatum)) then
      errmsg = 'Observation type not found: '//trim(obsTypeID)
      call store_error(errmsg)
      call store_error_unit(this%inUnitObs)
    end if
  end function get_obs_datum

  !> @ brief Set observation array values
  !!
  !!  Subroutine to set values in an observation array.
  !!
  !<
  subroutine set_obs_array(this, nObs, obsArray)
    ! -- dummy
    class(ObsType), intent(inout) :: this
    integer(I4B), intent(in) :: nObs !< number of observations
    type(ObsContainerType), dimension(nObs), intent(inout) :: obsArray !< observation array
    !
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: n
    type(ObserveType), pointer :: obsrv => null()
    !
    n = this%get_num()
    do i = 1, n
      obsrv => this%get_obs(i)
      obsArray(i)%obsrv => obsrv
    end do
  end subroutine set_obs_array

  !> @ brief Get an ObserveType object
  !!
  !!  Subroutine to get an ObserveType object from the list of observations
  !!  using an list index.
  !!
  !<
  function get_obs(this, indx) result(obsrv)
    ! -- dummy
    class(ObsType) :: this
    integer(I4B), intent(in) :: indx !< observation list index
    class(ObserveType), pointer :: obsrv !< observation ObserveType
    !
    obsrv => GetObsFromList(this%obsList, indx)
  end function get_obs

  !> @ brief Read observation blocks
  !!
  !!  Subroutine to read CONTIGUOUS block from the observation input file.
  !!
  !<
  subroutine read_obs_blocks(this, fname)
    ! -- dummy
    class(ObsType), intent(inout) :: this
    character(len=*), intent(inout) :: fname
    ! -- local
    integer(I4B) :: ierr, indexobsout, numspec
    logical :: fmtd, found, endOfBlock
    character(len=LENBIGLINE) :: pnamein, fnamein
    character(len=LENHUGELINE) :: line
    character(len=LINELENGTH) :: btagfound, message, word
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: tag
    character(len=20) :: accarg, bin, fmtarg
    type(ObserveType), pointer :: obsrv => null()
    type(ObsOutputType), pointer :: obsOutput => null()
    integer(I4B) :: ntabrows
    integer(I4B) :: ntabcols
    !
    ! -- initialize local variables
    numspec = -1
    errmsg = ''
    !
    inquire (unit=this%parser%iuactive, name=pnamein)
    call GetFileFromPath(pnamein, fnamein)
    !
    if (this%echo) then
      !
      ! -- create the observation table
      ! -- table dimensions
      ntabrows = 1
      ntabcols = 5
      !
      ! -- initialize table and define columns
      title = 'OBSERVATIONS READ FROM FILE "'//trim(fnamein)//'"'
      call table_cr(this%obstab, fnamein, title)
      call this%obstab%table_df(ntabrows, ntabcols, this%iout, &
                                finalize=.FALSE.)
      tag = 'NAME'
      call this%obstab%initialize_column(tag, LENOBSNAME, alignment=TABLEFT)
      tag = 'TYPE'
      call this%obstab%initialize_column(tag, LENOBSTYPE + 12, alignment=TABLEFT)
      tag = 'TIME'
      call this%obstab%initialize_column(tag, 12, alignment=TABLEFT)
      tag = 'LOCATION DATA'
      call this%obstab%initialize_column(tag, LENBOUNDNAME + 2, alignment=TABLEFT)
      tag = 'OUTPUT FILENAME'
      call this%obstab%initialize_column(tag, 80, alignment=TABLEFT)
    end if
    !
    found = .true.
    readblocks: do
      if (.not. found) exit
      !
      call this%parser%GetBlock('*', found, ierr, .true., .false., btagfound)
      if (.not. found) then
        exit readblocks
      end if
      this%blockTypeFound = btagfound
      !
      ! Get keyword, which should be FILEOUT
      call this%parser%GetStringCaps(word)
      if (word /= 'FILEOUT') then
        call store_error('CONTINUOUS keyword must be followed by '// &
                         '"FILEOUT" then by filename.')
        cycle
      end if
      !
      ! -- get name of output file
      call this%parser%GetString(fname)
      ! Fname is the output file name defined in the BEGIN line of the block.
      if (fname == '') then
        message = 'Error reading OBS input file, likely due to bad'// &
                  ' block or missing file name.'
        call store_error(message)
        cycle
      else if (this%obsOutputList%ContainsFile(fname)) then
        errmsg = 'OBS outfile "'//trim(fname)// &
                 '" is provided more than once.'
        call store_error(errmsg)
        cycle
      end if
      !
      ! -- look for BINARY option
      call this%parser%GetStringCaps(bin)
      if (bin == 'BINARY') then
        fmtarg = FORM
        accarg = ACCESS
        fmtd = .false.
      else
        fmtarg = 'FORMATTED'
        accarg = 'SEQUENTIAL'
        fmtd = .true.
      end if
      !
      ! -- open the output file
      numspec = 0
      call openfile(numspec, 0, fname, 'OBS OUTPUT', fmtarg, &
                    accarg, 'REPLACE')
      !
      ! -- add output file to list of output files and assign its
      !    FormattedOutput member appropriately
      call this%obsOutputList%Add(fname, numspec)
      indexobsout = this%obsOutputList%Count()
      obsOutput => this%obsOutputList%Get(indexobsout)
      obsOutput%FormattedOutput = fmtd
      !
      ! -- process lines defining observations
      select case (btagfound)
      case ('CONTINUOUS')
        !
        ! -- construct a continuous observation from each line in the block
        readblockcontinuous: do
          call this%parser%GetNextLine(endOfBlock)
          if (endOfBlock) exit
          call this%parser%GetCurrentLine(line)
          call ConstructObservation(obsrv, line, numspec, fmtd, &
                                    indexobsout, this%obsData, &
                                    this%parser%iuactive)
          !
          ! -- increment number of observations
          !    to be written to this output file.
          obsOutput => this%obsOutputList%Get(indexobsout)
          obsOutput%nobs = obsOutput%nobs + 1
          call AddObsToList(this%obsList, obsrv)
          !
          ! -- write line to the observation table
          if (this%echo) then
            call obsrv%WriteTo(this%obstab, btagfound, fname)
          end if
        end do readblockcontinuous
      case default
        errmsg = 'Error: Observation block type not recognized: '// &
                 trim(btagfound)
        call store_error(errmsg)
      end select
    end do readblocks
    !
    ! -- finalize the observation table
    if (this%echo) then
      call this%obstab%finalize_table()
    end if
    !
    ! -- determine if error condition occurs
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine read_obs_blocks

  !> @ brief Write observation data
  !!
  !!  Subroutine to write observation data for a time step for each observation
  !!  to the observation output file.
  !!
  !<
  subroutine write_obs_simvals(this)
    ! -- dummy
    class(ObsType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: iprec
    integer(I4B) :: numobs
    character(len=20) :: fmtc
    real(DP) :: simval
    class(ObserveType), pointer :: obsrv => null()
    !
    ! Write simulated values for observations
    iprec = this%iprecision
    fmtc = this%obsfmtcont
    ! -- iterate through all observations
    numobs = this%obsList%Count()
    do i = 1, numobs
      obsrv => this%get_obs(i)
      ! -- continuous observation
      simval = obsrv%CurrentTimeStepEndValue
      if (obsrv%FormattedOutput) then
        call write_fmtd_obs(fmtc, obsrv, this%obsOutputList, simval)
      else
        call write_unfmtd_obs(obsrv, iprec, this%obsOutputList, simval)
      end if
    end do
  end subroutine write_obs_simvals

end module ObsModule
