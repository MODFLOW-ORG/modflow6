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
                             LENVARNAME
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

  type, abstract, extends(NumericalPackageType) :: PbstBaseType

    character(len=8), dimension(:), pointer, contiguous :: status => null() !< active, inactive, constant
    character(len=LENPACKAGENAME) :: text = '' !< text string for package transport term
    integer(I4B), pointer :: ncv => null() !< number of control volumes
    integer(I4B), dimension(:), pointer, contiguous :: iboundpbst => null() !< package ibound
    logical, pointer, public :: active => null() !< logical indicating if a sensible heat flux object is active
    character(len=LINELENGTH), pointer, public :: inputFilename => null() !< a particular pbst input file name, could be for sensible heat flux or latent heat flux subpackages, for example
    type(TimeSeriesManagerType), pointer :: tsmanager => null()
    !
    ! -- table objects
    type(TableType), pointer :: inputtab => null() !< input table object

  contains

    procedure :: init
    procedure :: ar
    procedure :: rp
    procedure, private :: read_options
    procedure :: pbst_options
    procedure :: pbst_set_stressperiod
    procedure :: subpck_set_stressperiod
    procedure(read_option), deferred :: read_option
    procedure, private :: pbstbase_allocate_scalars
    procedure :: pbst_allocate_arrays
    procedure :: da => pbstbase_da
    procedure :: pbst_check_valid

  end type PbstBaseType

  abstract interface

    !> @brief Announce package and set pointers to variables
    !!
    !! Deferred procedure called by the PbstBaseType code to process a single
    !! keyword from the OPTIONS block of the package input file.
    !<
    function read_option(this, keyword) result(success)
      ! -- modules
      import PbstBaseType
      ! -- dummy
      class(PbstBaseType) :: this
      character(len=*), intent(in) :: keyword
      ! -- return
      logical :: success
    end function

  end interface

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
  !!  Method to allocate and read static data for the SHF package
  !<
  subroutine ar(this)
    ! -- dummy
    class(PbstBaseType) :: this !< ShfType object
    ! -- formats
    character(len=*), parameter :: fmtapt = &
      "(1x,/1x,'SHF -- SENSIBLE HEAT FLUX TRANSPORT PACKAGE, VERSION 1, 3/12/2025', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! -- print a message identifying the apt package.
    write (this%iout, fmtapt) this%inunit
    !
    ! -- Allocate arrays
    !call this%pbst_allocate_arrays()
    !
    ! -- Create time series manager
    call tsmanager_cr(this%tsmanager, this%iout, &
                      removeTsLinksOnCompletion=.true., &
                      extendTsToEndOfSimulation=.true.)
    !
    ! -- Read options
    call this%read_options()
  end subroutine ar

  !> @brief PaBST read and prepare for setting stress period information
  !<
  subroutine rp(this)
    ! -- module
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    use TdisModule, only: kper, nper
    ! -- dummy
    class(PbstBaseType) :: this !< ShfType object
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
  !! To be overridden by Pbst sub-packages
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

  !> @brief Read the SHF-specific options from the OPTIONS block
  !<
  subroutine read_options(this)
    ! -- dummy
    class(PbstBaseType) :: this
    ! -- local
    character(len=LINELENGTH) :: keyword
    character(len=MAXCHARLEN) :: fname
    logical :: endOfBlock
    logical(LGP) :: found
    integer(I4B) :: ierr
    ! -- formats
    character(len=*), parameter :: fmtts = &
      &"(4x, 'TIME-SERIES DATA WILL BE READ FROM FILE: ', a)"
    !
    ! -- Get options block
    call this%parser%GetBlock('OPTIONS', found, ierr, &
                              blockRequired=.false., supportOpenClose=.true.)
    !
    ! -- Parse options block if detected
    if (found) then
      write (this%iout, '(1x,a)') &
        'PROCESSING '//trim(adjustl(this%packName))//' OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) then
          exit
        end if
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('PRINT_INPUT')
          this%iprpak = 1
          write (this%iout, '(4x,a)') 'TIME-VARYING INPUT WILL BE PRINTED.'
        case ('TS6')
          !
          ! -- Add a time series file
          call this%parser%GetStringCaps(keyword)
          if (trim(adjustl(keyword)) /= 'FILEIN') then
            errmsg = &
              'TS6 keyword must be followed by "FILEIN" then by filename.'
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
            call ustop()
          end if
          call this%parser%GetString(fname)
          write (this%iout, fmtts) trim(fname)
          call this%tsmanager%add_tsfile(fname, this%inunit)
        case default
          !
          ! -- Check for child class options
          call this%pbst_options(keyword, found)
          !
          ! -- Defer to subtype to read the option;
          ! -- if the subtype can't handle it, report an error
          if (.not. found) then
            write (errmsg, '(a,3(1x,a),a)') &
              'Unknown', trim(adjustl(this%packName)), "option '", &
              trim(keyword), "'."
            call store_error(errmsg)
          end if
        end select
      end do
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%packName))//' OPTIONS'
    end if
  end subroutine read_options

  !> @ brief Read additional options for sub-package
  !!
  !!  Read additional options for the SFE boundary package. This method should
  !!  be overridden by option-processing routine that is in addition to the
  !!  base options available for all PbstBase packages.
  !<
  subroutine pbst_options(this, option, found)
    ! -- dummy
    class(PbstBaseType), intent(inout) :: this !< PbstBaseType object
    character(len=*), intent(inout) :: option !< option keyword string
    logical(LGP), intent(inout) :: found !< boolean indicating if the option was found
    !
    ! Return with found = .false.
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

  !> @ brief Allocate arrays
  !!
  !! Allocate base process-based stream temperature package transport arrays
  !<
  subroutine pbst_allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(PbstBaseType), intent(inout) :: this
    ! -- local
    integer(I4B) :: n
    !
    ! -- Note: For the time-being, no call to parent class allocation of arrays
    !    as is done in tsp-apt.f90, for example.  Remember that this class
    !    extends NumericalPackage.f90 which doesn't have a standard set of
    !    arrays to be allocated, only scalars.
    !call this%BndType%allocate_arrays()
    !
    ! -- allocate character array for status
    allocate (this%status(this%ncv))
    !
    ! -- initialize arrays
    do n = 1, this%ncv
      this%status(n) = 'ACTIVE'
    end do
  end subroutine pbst_allocate_arrays

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
    ! -- Deallocate time series manager
    deallocate (this%tsmanager)
    !
    ! -- deallocate scalars
    call mem_deallocate(this%ncv)
    !
    ! -- deallocate arrays
    call mem_deallocate(this%iboundpbst)
    !
    ! -- Deallocate parent
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
    ! -- formats
    ierr = 0
    if (itemno < 1 .or. itemno > this%ncv) then
      write (errmsg, '(a,1x,i6,1x,a,1x,i6)') &
        'Featureno ', itemno, 'must be > 0 and <= ', this%ncv
      call store_error(errmsg)
      ierr = 1
    end if
  end function pbst_check_valid

end module PbstBaseModule
