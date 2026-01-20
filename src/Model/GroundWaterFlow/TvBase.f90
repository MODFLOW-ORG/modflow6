!> @brief This module contains common time-varying property functionality
!!
!! This module contains methods implementing functionality common to both
!! time-varying hydraulic conductivity (TVK) and time-varying storage (TVS)
!! packages.
!!
!<
module TvBaseModule
  use BaseDisModule, only: DisBaseType
  use ConstantsModule, only: LINELENGTH, MAXCHARLEN, DZERO
  use KindModule, only: I4B, DP
  use NumericalPackageModule, only: NumericalPackageType
  use SimModule, only: count_errors, store_error, ustop
  use SimVariablesModule, only: errmsg
  use TdisModule, only: kper, nper, kstp
  use TimeSeriesLinkModule, only: TimeSeriesLinkType
  use TimeSeriesManagerModule, only: TimeSeriesManagerType, tsmanager_cr, &
                                     read_value_or_time_series_adv

  implicit none

  private

  public :: TvBaseType
  public :: tvbase_da

  type, abstract, extends(NumericalPackageType) :: TvBaseType
    type(TimeSeriesManagerType), pointer, private :: tsmanager => null()
  contains
    procedure :: init
    procedure :: ar
    procedure :: rp
    procedure :: ad
    procedure :: da => tvbase_da
    procedure, private :: tvbase_allocate_scalars
    procedure, private :: read_options
    procedure(ar_set_pointers), deferred :: ar_set_pointers
    procedure(read_option), deferred :: read_option
    procedure(get_pointer_to_value), deferred :: get_pointer_to_value
    procedure(set_changed_at), deferred :: set_changed_at
    procedure(reset_change_flags), deferred :: reset_change_flags
    procedure(validate_change), deferred :: validate_change
  end type TvBaseType

  abstract interface

    !> @brief Announce package and set pointers to variables
    !!
    !! Deferred procedure called by the TvBaseType code to announce the
    !! specific package version and set any required array and variable
    !! pointers from other packages.
    !<
    subroutine ar_set_pointers(this)
      ! -- modules
      import TvBaseType
      ! -- dummy
      class(TvBaseType) :: this
    end subroutine

    !> @brief Announce package and set pointers to variables
    !!
    !! Deferred procedure called by the TvBaseType code to process a single
    !! keyword from the OPTIONS block of the package input file.
    !<
    function read_option(this, keyword) result(success)
      ! -- modules
      import TvBaseType
      ! -- dummy
      class(TvBaseType) :: this
      character(len=*), intent(in) :: keyword
      ! -- return
      logical :: success
    end function

    !> @brief Get an array value pointer given a variable name and node index
    !!
    !! Deferred procedure called by the TvBaseType code to retrieve a pointer
    !! to a given node's value for a given named variable.
    !<
    function get_pointer_to_value(this, n, varName) result(bndElem)
      ! -- modules
      use KindModule, only: I4B, DP
      import TvBaseType
      ! -- dummy
      class(TvBaseType) :: this
      integer(I4B), intent(in) :: n
      character(len=*), intent(in) :: varName
      ! -- return
      real(DP), pointer :: bndElem
    end function

    !> @brief Mark property changes as having occurred at (kper, kstp)
    !!
    !! Deferred procedure called by the TvBaseType code when a property value
    !! change occurs at (kper, kstp).
    !<
    subroutine set_changed_at(this, kper, kstp)
      ! -- modules
      use KindModule, only: I4B
      import TvBaseType
      ! -- dummy
      class(TvBaseType) :: this
      integer(I4B), intent(in) :: kper
      integer(I4B), intent(in) :: kstp
    end subroutine

    !> @brief Clear all per-node change flags
    !!
    !! Deferred procedure called by the TvBaseType code when a new time step
    !! commences, indicating that any previously set per-node property value
    !! change flags should be reset.
    !<
    subroutine reset_change_flags(this)
      ! -- modules
      import TvBaseType
      ! -- dummy
      class(TvBaseType) :: this
    end subroutine

    !> @brief Check that a given property value is valid
    !!
    !! Deferred procedure called by the TvBaseType code after a property value
    !! change occurs to perform any required validity checks on the value of
    !! the given variable at the given node. Perform any required updates to
    !! the property value if it is valid, or log an error if not.
    !<
    subroutine validate_change(this, n, varName)
      ! -- modules
      use KindModule, only: I4B
      import TvBaseType
      ! -- dummy
      class(TvBaseType) :: this
      integer(I4B), intent(in) :: n
      character(len=*), intent(in) :: varName
    end subroutine

  end interface

contains

  !> @brief Initialize the TvBaseType object
  !!
  !! Allocate and initialize data members of the object.
  !<
  subroutine init(this, name_model, pakname, ftype, inunit, iout)
    ! -- dummy
    class(TvBaseType) :: this
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: pakname
    character(len=*), intent(in) :: ftype
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    !
    call this%set_names(1, name_model, pakname, ftype)
    call this%tvbase_allocate_scalars()
    this%inunit = inunit
    this%iout = iout
    call this%parser%Initialize(this%inunit, this%iout)
  end subroutine init

  !> @brief Allocate scalar variables
  !!
  !! Allocate scalar data members of the object.
  !<
  subroutine tvbase_allocate_scalars(this)
    ! -- dummy
    class(TvBaseType) :: this
    !
    ! -- Call standard NumericalPackageType allocate scalars
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate time series manager
    allocate (this%tsmanager)
  end subroutine tvbase_allocate_scalars

  !> @brief Allocate and read method for package
  !!
  !! Allocate and read static data for the package.
  !<
  subroutine ar(this, dis)
    ! -- dummy
    class(TvBaseType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    !
    ! -- Set pointers to other package variables and announce package
    this%dis => dis
    call this%ar_set_pointers()
    !
    ! -- Create time series manager
    call tsmanager_cr(this%tsmanager, &
                      this%iout, &
                      removeTsLinksOnCompletion=.true., &
                      extendTsToEndOfSimulation=.true.)
    !
    ! -- Read options
    call this%read_options()
    !
    ! -- Now that time series will have been read, need to call the df routine
    ! -- to define the manager
    call this%tsmanager%tsmanager_df()
    !
    ! -- Terminate if any errors were encountered
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
      call ustop()
    end if
  end subroutine ar

  !> @brief Read OPTIONS block of package input file
  !!
  !! Reads the OPTIONS block of the package's input file, deferring to the
  !! derived type to process any package-specific keywords.
  !<
  subroutine read_options(this)
    ! -- dummy
    class(TvBaseType) :: this
    ! -- local
    character(len=LINELENGTH) :: keyword
    character(len=MAXCHARLEN) :: fname
    logical :: isfound
    logical :: endOfBlock
    integer(I4B) :: ierr
    ! -- formats
    character(len=*), parameter :: fmtts = &
      &"(4x, 'TIME-SERIES DATA WILL BE READ FROM FILE: ', a)"
    !
    ! -- Get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, &
                              blockRequired=.false., supportOpenClose=.true.)
    !
    ! -- Parse options block if detected
    if (isfound) then
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
          ! -- Defer to subtype to read the option;
          ! -- if the subtype can't handle it, report an error
          if (.not. this%read_option(keyword)) then
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

  !> @brief Read and prepare method for package
  !!
  !! Read and prepare stress period data for the package.
  !<
  subroutine rp(this)
    ! -- dummy
    class(TvBaseType) :: this
    ! -- local
    character(len=LINELENGTH) :: line, cellid, varName, text
    logical :: isfound, endOfBlock, haveChanges
    integer(I4B) :: ierr, node
    real(DP), pointer :: bndElem => null()
    ! -- formats
    character(len=*), parameter :: fmtblkerr = &
      &"('Looking for BEGIN PERIOD iper.  Found ', a, ' instead.')"
    character(len=*), parameter :: fmtvalchg = &
      "(a, ' package: Setting ', a, ' value for cell ', a, ' at start of &
      &stress period ', i0, ' = ', g12.5)"
    !
    if (this%inunit == 0) return
    !
    ! -- Get stress period data
    if (this%ionper < kper) then
      !
      ! -- Get PERIOD block
      call this%parser%GetBlock('PERIOD', isfound, ierr, &
                                supportOpenClose=.true., &
                                blockRequired=.false.)
      if (isfound) then
        !
        ! -- Read ionper and check for increasing period numbers
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
        end if
      end if
    end if
    !
    ! -- Read data if ionper == kper
    if (this%ionper == kper) then
      !
      ! -- Reset per-node property change flags
      call this%reset_change_flags()
      !
      haveChanges = .false.
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) then
          exit
        end if
        !
        ! -- Read cell ID
        call this%parser%GetCellid(this%dis%ndim, cellid)
        node = this%dis%noder_from_cellid(cellid, this%parser%iuactive, &
                                          this%iout)
        !
        ! -- Validate cell ID
        if (node < 1 .or. node > this%dis%nodes) then
          write (errmsg, '(a,2(1x,a))') &
            'CELLID', cellid, 'is not in the active model domain.'
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- Read variable name
        call this%parser%GetStringCaps(varName)
        !
        ! -- Get a pointer to the property value given by varName for the node
        ! -- with the specified cell ID
        bndElem => this%get_pointer_to_value(node, varName)
        if (.not. associated(bndElem)) then
          write (errmsg, '(a,3(1x,a),a)') &
            'Unknown', trim(adjustl(this%packName)), "variable '", &
            trim(varName), "'."
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- Read and apply value or time series link
        call this%parser%GetString(text)
        call read_value_or_time_series_adv(text, node, 0, bndElem, &
                                           this%packName, 'BND', &
                                           this%tsmanager, this%iprpak, &
                                           varName)
        !
        ! -- Report value change
        if (this%iprpak /= 0) then
          write (this%iout, fmtvalchg) &
            trim(adjustl(this%packName)), trim(varName), trim(cellid), &
            kper, bndElem
        end if
        !
        ! -- Validate the new property value
        call this%validate_change(node, varName)
        haveChanges = .true.
      end do
      !
      ! -- Record that any changes were made at the first time step of the
      ! -- stress period
      if (haveChanges) then
        call this%set_changed_at(kper, 1)
      end if
    end if
    !
    ! -- Terminate if errors were encountered in the PERIOD block
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
      call ustop()
    end if
  end subroutine rp

  !> @brief Advance the package
  !!
  !! Advance data for a new time step.
  !<
  subroutine ad(this)
    ! -- dummy
    class(TvBaseType) :: this
    ! -- local
    integer(I4B) :: i, n, numlinks
    type(TimeSeriesLinkType), pointer :: tsLink
    !
    ! -- Advance the time series manager;
    ! -- this will apply any time series changes to property values
    call this%tsmanager%ad()
    !
    ! -- If there are no time series property changes,
    ! -- there is nothing else to be done
    numlinks = this%tsmanager%CountLinks('BND')
    if (numlinks <= 0) then
      return
    end if
    !
    ! -- Record that changes were made at the current time step
    ! -- (as we have at least one active time series link)
    call this%set_changed_at(kper, kstp)
    !
    ! -- Reset node K change flags at all time steps except the first of each
    ! -- period (the first is done in rp(), to allow non-time series changes)
    if (kstp /= 1) then
      call this%reset_change_flags()
    end if
    !
    ! -- Validate all new property values
    do i = 1, numlinks
      tsLink => this%tsmanager%GetLink('BND', i)
      n = tsLink%iRow
      call this%validate_change(n, tsLink%Text)
    end do
    !
    ! -- Terminate if there were errors
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
      call ustop()
    end if
  end subroutine ad

  !> @brief Deallocate package memory
  !!
  !! Deallocate package scalars and arrays.
  !<
  subroutine tvbase_da(this)
    ! -- dummy
    class(TvBaseType) :: this
    !
    ! -- Deallocate time series manager
    deallocate (this%tsmanager)
    !
    ! -- Deallocate parent
    call this%NumericalPackageType%da()
  end subroutine tvbase_da

end module TvBaseModule
