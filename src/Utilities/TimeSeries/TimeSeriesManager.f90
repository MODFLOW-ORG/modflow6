module TimeSeriesManagerModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, LENPACKAGENAME, MAXCHARLEN, &
                             LINELENGTH, LENTIMESERIESNAME
  use HashTableModule, only: HashTableType, hash_table_cr, &
                             hash_table_da
  use InputOutputModule, only: same_word, UPCASE
  use ListModule, only: ListType
  use SimModule, only: store_error, store_error_unit
  use TdisModule, only: delt, kper, kstp, totim, totimc, &
                        totimsav
  use TimeSeriesFileListModule, only: TimeSeriesFileListType
  use TimeSeriesLinkModule, only: TimeSeriesLinkType, &
                                  ConstructTimeSeriesLink, &
                                  GetTimeSeriesLinkFromList, &
                                  AddTimeSeriesLinkToList
  use TimeSeriesModule, only: TimeSeriesContainerType, &
                              TimeSeriesFileType, &
                              TimeSeriesType

  implicit none

  private
  public :: TimeSeriesManagerType, read_value_or_time_series, &
            read_value_or_time_series_adv, &
            var_timeseries, tsmanager_cr

  type TimeSeriesManagerType
    integer(I4B), public :: iout = 0 ! output unit number
    type(TimeSeriesFileListType), pointer, public :: tsfileList => null() ! list of ts files objs
    type(ListType), pointer, public :: boundTsLinks => null() ! links to bound and aux
    integer(I4B) :: numtsfiles = 0 ! number of ts files
    character(len=MAXCHARLEN), allocatable, dimension(:) :: tsfiles ! list of ts files
    logical, private :: removeTsLinksOnCompletion = .false. ! flag indicating whether time series links should be removed in ad() once simulation time passes the end of the time series
    logical, private :: extendTsToEndOfSimulation = .false. ! flag indicating whether time series should be extended to provide their final value for all times after the series end time
    type(ListType), pointer, private :: auxvarTsLinks => null() ! list of aux links
    type(HashTableType), pointer, private :: BndTsHashTable => null() ! hash of ts to tsobj
    type(TimeSeriesContainerType), allocatable, dimension(:), &
      private :: TsContainers

  contains

    ! -- Public procedures
    procedure, public :: tsmanager_df
    procedure, public :: ad => tsmgr_ad
    procedure, public :: da => tsmgr_da
    procedure, public :: add_tsfile
    procedure, public :: CountLinks
    procedure, public :: GetLink
    procedure, public :: Reset
    procedure, public :: HashBndTimeSeries
    ! -- Private procedures
    procedure, private :: get_time_series
    procedure, private :: make_link
  end type TimeSeriesManagerType

contains

  !> @brief Create the tsmanager
  !<
  subroutine tsmanager_cr(this, iout, removeTsLinksOnCompletion, &
                          extendTsToEndOfSimulation)
    ! -- dummy
    type(TimeSeriesManagerType) :: this
    integer(I4B), intent(in) :: iout
    logical, intent(in), optional :: removeTsLinksOnCompletion
    logical, intent(in), optional :: extendTsToEndOfSimulation
    !
    this%iout = iout
    if (present(removeTsLinksOnCompletion)) then
      this%removeTsLinksOnCompletion = removeTsLinksOnCompletion
    end if
    if (present(extendTsToEndOfSimulation)) then
      this%extendTsToEndOfSimulation = extendTsToEndOfSimulation
    end if
    allocate (this%boundTsLinks)
    allocate (this%auxvarTsLinks)
    allocate (this%tsfileList)
    allocate (this%tsfiles(1000))
  end subroutine tsmanager_cr

  !> @brief Define time series manager object
  !<
  subroutine tsmanager_df(this)
    ! -- dummy
    class(TimeSeriesManagerType) :: this
    !
    if (this%numtsfiles > 0) then
      call this%HashBndTimeSeries()
    end if
  end subroutine tsmanager_df

  !> @brief Add a time series file to this manager
  !<
  subroutine add_tsfile(this, fname, inunit)
    ! -- modules
    use SimModule, only: store_error, store_error_unit
    use ArrayHandlersModule, only: ExpandArray
    ! -- dummy
    class(TimeSeriesManagerType) :: this
    character(len=*), intent(in) :: fname
    integer(I4B), intent(in) :: inunit
    ! -- local
    integer(I4B) :: isize
    integer(I4B) :: i
    class(TimeSeriesFileType), pointer :: tsfile => null()
    !
    ! -- Check for fname duplicates
    if (this%numtsfiles > 0) then
      do i = 1, this%numtsfiles
        if (this%tsfiles(i) == fname) then
          call store_error('Found duplicate time-series file name: '//trim(fname))
          call store_error_unit(inunit)
        end if
      end do
    end if
    !
    ! -- Save fname
    this%numtsfiles = this%numtsfiles + 1
    isize = size(this%tsfiles)
    if (this%numtsfiles > isize) then
      call ExpandArray(this%tsfiles, 1000)
    end if
    this%tsfiles(this%numtsfiles) = fname
    !
    ! --
    call this%tsfileList%Add(fname, this%iout, tsfile)
  end subroutine add_tsfile

  !> @brief Time step (or subtime step) advance. Call this each time step or
  !! subtime step
  !<
  subroutine tsmgr_ad(this)
    ! -- dummy
    class(TimeSeriesManagerType) :: this
    ! -- local
    type(TimeSeriesLinkType), pointer :: tsLink => null()
    type(TimeSeriesType), pointer :: timeseries => null()
    integer(I4B) :: i, nlinks, nauxlinks
    real(DP) :: begintime, endtime, tsendtime
    character(len=LENPACKAGENAME + 2) :: pkgID
    ! -- formats
    character(len=*), parameter :: fmt5 = &
      "(/,'Time-series controlled values in stress period: ', i0, &
      &', time step ', i0, ':')"
10  format(a, ' package: Boundary ', i0, ', entry ', i0, &
           ' value from time series "', a, '" = ', g12.5)
15  format(a, ' package: Boundary ', i0, ', entry ', i0, &
           ' value from time series "', a, '" = ', g12.5, ' (', a, ')')
20  format(a, ' package: Boundary ', i0, ', ', a, &
           ' value from time series "', a, '" = ', g12.5)
25  format(a, ' package: Boundary ', i0, ', ', a, &
           ' value from time series "', a, '" = ', g12.5, ' (', a, ')')
    !
    ! -- Initialize time variables
    begintime = totimc
    endtime = begintime + delt
    !
    ! -- Determine number of ts links
    nlinks = this%boundtslinks%Count()
    nauxlinks = this%auxvartslinks%Count()
    !
    ! -- Iterate through auxvartslinks and replace specified
    !    elements of auxvar with average value obtained from
    !    appropriate time series.  Need to do auxvartslinks
    !    first because they may be a multiplier column
    i = 1
    do while (i <= nauxlinks)
      tsLink => GetTimeSeriesLinkFromList(this%auxvarTsLinks, i)
      timeseries => tsLink%timeSeries
      !
      ! -- Remove time series link once its end time has passed, if requested
      if (this%removeTsLinksOnCompletion) then
        tsendtime = timeseries%FindLatestTime(.true.)
        if (tsendtime < begintime) then
          call this%auxvarTsLinks%RemoveNode(i, .TRUE.)
          nauxlinks = this%auxvartslinks%Count()
          cycle
        end if
      end if
      !
      if (i == 1) then
        if (tsLink%Iprpak == 1) then
          write (this%iout, fmt5) kper, kstp
        end if
      end if
      tsLink%BndElement = timeseries%GetValue(begintime, endtime, &
                                              this%extendTsToEndOfSimulation)
      !
      ! -- Write time series values to output file
      if (tsLink%Iprpak == 1) then
        pkgID = '"'//trim(tsLink%PackageName)//'"'
        if (tsLink%Text == '') then
          if (tsLink%BndName == '') then
            write (this%iout, 10) trim(pkgID), tsLink%IRow, tsLink%JCol, &
              trim(tsLink%timeSeries%Name), &
              tsLink%BndElement
          else
            write (this%iout, 15) trim(pkgID), tsLink%IRow, tsLink%JCol, &
              trim(tsLink%timeSeries%Name), &
              tsLink%BndElement, trim(tsLink%BndName)
          end if
        else
          if (tsLink%BndName == '') then
            write (this%iout, 20) trim(pkgID), tsLink%IRow, trim(tsLink%Text), &
              trim(tsLink%timeSeries%Name), &
              tsLink%BndElement
          else
            write (this%iout, 25) trim(pkgID), tsLink%IRow, trim(tsLink%Text), &
              trim(tsLink%timeSeries%Name), &
              tsLink%BndElement, trim(tsLink%BndName)
          end if
        end if
      end if
      !
      i = i + 1
    end do
    !
    ! -- Iterate through boundtslinks and replace specified
    !    elements of bound with average value obtained from
    !    appropriate time series. (For list-type packages)
    i = 1
    do while (i <= nlinks)
      tsLink => GetTimeSeriesLinkFromList(this%boundTsLinks, i)
      timeseries => tsLink%timeSeries
      !
      ! -- Remove time series link once its end time has passed, if requested
      if (this%removeTsLinksOnCompletion) then
        tsendtime = timeseries%FindLatestTime(.true.)
        if (tsendtime < begintime) then
          call this%boundTsLinks%RemoveNode(i, .TRUE.)
          nlinks = this%boundTsLinks%Count()
          cycle
        end if
      end if
      !
      if (i == 1 .and. nauxlinks == 0) then
        if (tsLink%Iprpak == 1) then
          write (this%iout, fmt5) kper, kstp
        end if
      end if
      ! this part needs to be different for MAW because MAW does not use
      ! bound array for well rate (although rate is stored in
      ! this%bound(4,ibnd)), it uses this%mawwells(n)%rate%value
      if (tsLink%UseDefaultProc) then
        timeseries => tsLink%timeSeries
        tsLink%BndElement = timeseries%GetValue(begintime, endtime, &
                                                this%extendTsToEndOfSimulation)
        !
        ! -- If multiplier is active and it applies to this element,
        !    do the multiplication.  This must be done after the auxlinks
        !    have been calculated in case iauxmultcol is being used.
        if (associated(tsLink%RMultiplier)) then
          tsLink%BndElement = tsLink%BndElement * tsLink%RMultiplier
        end if
        !
        ! -- Write time series values to output files
        if (tsLink%Iprpak == 1) then
          pkgID = '"'//trim(tsLink%PackageName)//'"'
          if (tsLink%Text == '') then
            if (tsLink%BndName == '') then
              write (this%iout, 10) trim(pkgID), tsLink%IRow, tsLink%JCol, &
                trim(tsLink%timeSeries%Name), &
                tsLink%BndElement
            else
              write (this%iout, 15) trim(pkgID), tsLink%IRow, tsLink%JCol, &
                trim(tsLink%timeSeries%Name), &
                tsLink%BndElement, trim(tsLink%BndName)
            end if
          else
            if (tsLink%BndName == '') then
              write (this%iout, 20) trim(pkgID), tsLink%IRow, trim(tsLink%Text), &
                trim(tsLink%timeSeries%Name), &
                tsLink%BndElement
            else
              write (this%iout, 25) trim(pkgID), tsLink%IRow, trim(tsLink%Text), &
                trim(tsLink%timeSeries%Name), &
                tsLink%BndElement, trim(tsLink%BndName)
            end if
          end if
        end if
        !
        ! -- If conversion from flux to flow is required, multiply by cell area
        if (tsLink%ConvertFlux) then
          tsLink%BndElement = tsLink%BndElement * tsLink%CellArea
        end if
      end if
      !
      i = i + 1
    end do
    !
    ! -- Finish with ending line
    if (nlinks + nauxlinks > 0) then
      if (tsLink%Iprpak == 1) then
        write (this%iout, '()')
      end if
    end if
  end subroutine tsmgr_ad

  !> @brief Deallocate memory
  !<
  subroutine tsmgr_da(this)
    ! -- dummy
    class(TimeSeriesManagerType) :: this
    ! -- local
    !
    ! -- Deallocate time-series links in boundTsLinks
    call this%boundTsLinks%Clear(.true.)
    deallocate (this%boundTsLinks)
    !
    ! -- Deallocate time-series links in auxvarTsLinks
    call this%auxvarTsLinks%Clear(.true.)
    deallocate (this%auxvarTsLinks)
    !
    ! -- Deallocate tsfileList
    call this%tsfileList%da()
    deallocate (this%tsfileList)
    !
    ! -- Deallocate the hash table
    if (associated(this%BndTsHashTable)) then
      call hash_table_da(this%BndTsHashTable)
    end if
    !
    deallocate (this%tsfiles)
  end subroutine tsmgr_da

  !> @brief Call this when a new BEGIN PERIOD block is read for a new stress
  !! period
  !<
  subroutine Reset(this, pkgName)
    ! -- dummy
    class(TimeSeriesManagerType) :: this
    character(len=*), intent(in) :: pkgName
    ! -- local
    integer(I4B) :: i, nlinks
    type(TimeSeriesLinkType), pointer :: tslink
    !
    ! Zero out values for time-series controlled stresses.
    ! Also deallocate all tslinks too.
    ! Then when time series are
    ! specified in this or another stress period,
    ! a new tslink would be set up.
    !
    ! Reassign all linked elements to zero
    nlinks = this%boundTsLinks%Count()
    do i = 1, nlinks
      tslink => GetTimeSeriesLinkFromList(this%boundTsLinks, i)
      if (associated(tslink)) then
        if (tslink%PackageName == pkgName) then
          tslink%BndElement = DZERO
        end if
      end if
    end do
    !
    ! Remove links belonging to calling package
    nlinks = this%boundTsLinks%Count()
    do i = nlinks, 1, -1
      tslink => GetTimeSeriesLinkFromList(this%boundTsLinks, i)
      if (associated(tslink)) then
        if (tslink%PackageName == pkgName) then
          call this%boundTsLinks%RemoveNode(i, .true.)
        end if
      end if
    end do
    nlinks = this%auxvarTsLinks%Count()
    do i = nlinks, 1, -1
      tslink => GetTimeSeriesLinkFromList(this%auxvarTsLinks, i)
      if (associated(tslink)) then
        if (tslink%PackageName == pkgName) then
          call this%auxvarTsLinks%RemoveNode(i, .true.)
        end if
      end if
    end do
  end subroutine Reset

  !> @brief Make link
  !<
  subroutine make_link(this, timeSeries, pkgName, auxOrBnd, bndElem, &
                       irow, jcol, iprpak, tsLink, text, bndName)
    ! -- dummy
    class(TimeSeriesManagerType), intent(inout) :: this
    type(TimeSeriesType), pointer, intent(inout) :: timeSeries
    character(len=*), intent(in) :: pkgName
    character(len=3), intent(in) :: auxOrBnd
    real(DP), pointer, intent(inout) :: bndElem
    integer(I4B), intent(in) :: irow, jcol
    integer(I4B), intent(in) :: iprpak
    type(TimeSeriesLinkType), pointer, intent(inout) :: tsLink
    character(len=*), intent(in) :: text
    character(len=*), intent(in) :: bndName
    !
    tsLink => null()
    call ConstructTimeSeriesLink(tsLink, timeSeries, pkgName, &
                                 auxOrBnd, bndElem, irow, jcol, iprpak)
    if (associated(tsLink)) then
      if (auxOrBnd == 'BND') then
        call AddTimeSeriesLinkToList(this%boundTsLinks, tsLink)
      elseif (auxOrBnd == 'AUX') then
        call AddTimeSeriesLinkToList(this%auxvarTsLinks, tsLink)
      else
        call store_error('programmer error in make_link', terminate=.TRUE.)
      end if
      tsLink%Text = text
      tsLink%BndName = bndName
    end if
  end subroutine make_link

  !> @brief Get link
  !<
  function GetLink(this, auxOrBnd, indx) result(tsLink)
    ! -- dummy
    class(TimeSeriesManagerType) :: this
    character(len=3), intent(in) :: auxOrBnd
    integer(I4B), intent(in) :: indx
    type(TimeSeriesLinkType), pointer :: tsLink
    ! -- local
    type(ListType), pointer :: list
    !
    list => null()
    tsLink => null()
    !
    select case (auxOrBnd)
    case ('AUX')
      list => this%auxvarTsLinks
    case ('BND')
      list => this%boundTsLinks
    end select
    !
    if (associated(list)) then
      tsLink => GetTimeSeriesLinkFromList(list, indx)
    end if
  end function GetLink

  !> @brief Count links
  !<
  function CountLinks(this, auxOrBnd)
    ! -- return
    integer(I4B) :: CountLinks
    ! -- dummy
    class(TimeSeriesManagerType) :: this
    character(len=3), intent(in) :: auxOrBnd
    !
    CountLinks = 0
    if (auxOrBnd == 'BND') then
      CountLinks = this%boundTsLinks%Count()
    elseif (auxOrBnd == 'AUX') then
      CountLinks = this%auxvarTsLinks%count()
    end if
  end function CountLinks

  !> @brief Get time series
  !<
  function get_time_series(this, name) result(res)
    ! -- dummy
    class(TimeSeriesManagerType) :: this
    character(len=*), intent(in) :: name
    ! -- result
    type(TimeSeriesType), pointer :: res
    ! -- local
    integer(I4B) :: indx
    !
    ! Get index from hash table, get time series from TsContainers,
    !     and assign result to time series contained in link.
    res => null()
    indx = this%BndTsHashTable%get(name)
    if (indx > 0) then
      res => this%TsContainers(indx)%timeSeries
    end if
  end function get_time_series

  !> @brief Store all boundary (stress) time series links in TsContainers
  !! and construct hash table BndTsHashTable
  !<
  subroutine HashBndTimeSeries(this)
    ! -- dummy
    class(TimeSeriesManagerType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i, j, k, numtsfiles, numts
    character(len=LENTIMESERIESNAME) :: name
    type(TimeSeriesFileType), pointer :: tsfile => null()
    !
    ! Initialize the hash table
    call hash_table_cr(this%BndTsHashTable)
    !
    ! Allocate the TsContainers array to accommodate all time-series links.
    numts = this%tsfileList%CountTimeSeries()
    allocate (this%TsContainers(numts))
    !
    ! Store a pointer to each time series in the TsContainers array
    ! and put its key (time-series name) and index in the hash table.
    numtsfiles = this%tsfileList%Counttsfiles()
    k = 0
    do i = 1, numtsfiles
      tsfile => this%tsfileList%Gettsfile(i)
      numts = tsfile%Count()
      do j = 1, numts
        k = k + 1
        this%TsContainers(k)%timeSeries => tsfile%GetTimeSeries(j)
        if (associated(this%TsContainers(k)%timeSeries)) then
          name = this%TsContainers(k)%timeSeries%Name
          call this%BndTsHashTable%add(name, k)
        end if
      end do
    end do
  end subroutine HashBndTimeSeries

  ! -- Non-type-bound procedures

  !> @brief Call this subroutine if the time-series link is available or needed
  !!
  !! This routine assumes that there is not an existing link for the
  !! specified package and array row and column.  For the standard
  !! boundary package, all links are removed by calling the tsmanager%reset()
  !! method.  For advanced packages, there is a separate routine called
  !! read_value_or_time_series_adv, which should be used instead of this one.
  !<
  subroutine read_value_or_time_series(textInput, ii, jj, bndElem, pkgName, &
                                       auxOrBnd, tsManager, iprpak, tsLink)
    ! dummy
    character(len=*), intent(in) :: textInput
    integer(I4B), intent(in) :: ii
    integer(I4B), intent(in) :: jj
    real(DP), pointer, intent(inout) :: bndElem
    character(len=*), intent(in) :: pkgName
    character(len=3), intent(in) :: auxOrBnd
    type(TimeSeriesManagerType), intent(inout) :: tsManager
    integer(I4B), intent(in) :: iprpak
    type(TimeSeriesLinkType), pointer, intent(inout) :: tsLink
    ! local
    type(TimeSeriesType), pointer :: timeseries => null()
    integer(I4B) :: istat
    real(DP) :: r
    character(len=LINELENGTH) :: errmsg
    character(len=LENTIMESERIESNAME) :: tsNameTemp

    read (textInput, *, iostat=istat) r
    if (istat == 0) then
      bndElem = r
    else

      ! Check to see if this is a time series name
      tsNameTemp = textInput
      call UPCASE(tsNameTemp)
      timeseries => tsManager%get_time_series(tsNameTemp)

      ! If this is a time series, then create a link between
      ! the time series and the row and column position of
      ! bndElem
      if (associated(timeseries)) then
        ! -- Assign value from time series to current
        !    array element
        r = timeseries%GetValue(totimsav, totim, &
                                tsManager%extendTsToEndOfSimulation)
        bndElem = r
        ! Make a new link between the time series and the row and
        ! column position in the array
        call tsManager%make_link(timeseries, pkgName, auxOrBnd, bndElem, &
                                 ii, jj, iprpak, tsLink, '', '')
      else
        errmsg = 'Error in list input. Expected numeric value or '// &
                 "time-series name, but found '"//trim(textInput)//"'."
        call store_error(errmsg)
      end if
    end if

  end subroutine read_value_or_time_series

  !> @brief Call this subroutine from advanced packages to define timeseries
  !! link for a variable (varName).
  !!
  !! Arguments are as follows:
  !!    textInput    : string that is either a float or a string name
  !!    ii           : column number
  !!    jj           : row number
  !!    bndElem      : pointer to a position in an array in package pkgName
  !!    pkgName      : package name
  !!    auxOrBnd     : 'AUX' or 'BND' keyword
  !!    tsManager    : timeseries manager object for package
  !!    iprpak       : integer flag indicating if interpolated timeseries values
  !!                   should be printed to package iout during TsManager%ad()
  !!    varName      : variable name
  !<
  subroutine read_value_or_time_series_adv(textInput, ii, jj, bndElem, pkgName, &
                                           auxOrBnd, tsManager, iprpak, varName)
    ! -- dummy
    character(len=*), intent(in) :: textInput
    integer(I4B), intent(in) :: ii
    integer(I4B), intent(in) :: jj
    real(DP), pointer, intent(inout) :: bndElem
    character(len=*), intent(in) :: pkgName
    character(len=3), intent(in) :: auxOrBnd
    type(TimeSeriesManagerType), intent(inout) :: tsManager
    integer(I4B), intent(in) :: iprpak
    character(len=*), intent(in) :: varName
    ! -- local
    integer(I4B) :: istat
    real(DP) :: v
    character(len=LINELENGTH) :: errmsg
    character(len=LENTIMESERIESNAME) :: tsNameTemp
    logical :: found
    type(TimeSeriesType), pointer :: timeseries => null()
    type(TimeSeriesLinkType), pointer :: tsLink => null()
    !
    ! -- attempt to read textInput as a real value
    read (textInput, *, iostat=istat) v
    !
    ! -- numeric value
    if (istat == 0) then
      !
      ! -- Numeric value was successfully read.
      bndElem = v
      !
      ! -- remove existing link if it exists for this boundary element
      found = remove_existing_link(tsManager, ii, jj, pkgName, &
                                   auxOrBnd, varName)
      !
      ! -- timeseries
    else
      !
      ! -- attempt to read numeric value from textInput failed.
      !    Text should be a time-series name.
      tsNameTemp = textInput
      call UPCASE(tsNameTemp)
      !
      ! -- if textInput is a time-series name, get average value
      !    from time series.
      timeseries => tsManager%get_time_series(tsNameTemp)
      !
      ! -- create a time series link and add it to the package
      !    list of time series links used by the array.
      if (associated(timeseries)) then
        !
        ! -- Assign average value from time series to current array element
        v = timeseries%GetValue(totimsav, totim, &
                                tsManager%extendTsToEndOfSimulation)
        bndElem = v
        !
        ! -- remove existing link if it exists for this boundary element
        found = remove_existing_link(tsManager, ii, jj, &
                                     pkgName, auxOrBnd, varName)
        !
        ! -- Add link to the list.
        call tsManager%make_link(timeseries, pkgName, auxOrBnd, bndElem, &
                                 ii, jj, iprpak, tsLink, varName, '')
        !
        ! -- not a valid timeseries name
      else
        errmsg = 'Error in list input. Expected numeric value or '// &
                 "time-series name, but found '"//trim(textInput)//"'."
        call store_error(errmsg)
      end if
    end if
  end subroutine read_value_or_time_series_adv

! -- private subroutines

  !> @brief Remove an existing timeseries link if it is defined.
  !!
  !! Arguments are as follows:
  !!    tsManager    : timeseries manager object for package
  !!    ii           : column number
  !!    jj           : row number
  !!    pkgName      : package name
  !!    auxOrBnd     : 'AUX' or 'BND' keyword
  !!    varName      : variable name
  !<
  function remove_existing_link(tsManager, ii, jj, &
                                pkgName, auxOrBnd, varName) result(found)
    ! -- return
    logical :: found
    ! -- dummy
    type(TimeSeriesManagerType), intent(inout) :: tsManager
    integer(I4B), intent(in) :: ii
    integer(I4B), intent(in) :: jj
    character(len=*), intent(in) :: pkgName
    character(len=3), intent(in) :: auxOrBnd
    character(len=*), intent(in) :: varName
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: nlinks
    integer(I4B) :: removeLink
    type(TimeSeriesLinkType), pointer :: tslTemp => null()
    !
    ! -- determine if link exists
    nlinks = tsManager%CountLinks(auxOrBnd)
    found = .FALSE.
    removeLink = -1
    csearchlinks: do i = 1, nlinks
      tslTemp => tsManager%GetLink(auxOrBnd, i)
      !
      ! -- Check ii against iRow, jj against jCol, and varName
      !    against Text member of link
      if (tslTemp%PackageName == pkgName) then
        !
        ! -- This array element is already linked to a time series.
        if (tslTemp%IRow == ii .and. tslTemp%JCol == jj .and. &
            same_word(tslTemp%Text, varName)) then
          found = .TRUE.
          removeLink = i
          exit csearchlinks
        end if
      end if
    end do csearchlinks
    !
    ! -- remove link if it was found
    if (removeLink > 0) then
      if (auxOrBnd == 'BND') then
        call tsManager%boundTsLinks%RemoveNode(removeLink, .TRUE.)
      else if (auxOrBnd == 'AUX') then
        call tsManager%auxvarTsLinks%RemoveNode(removeLink, .TRUE.)
      end if
    end if
  end function remove_existing_link

  !> @brief Determine if a timeseries link with varName is defined.
  !!
  !! Arguments are as follows:
  !!    tsManager    : timeseries manager object for package
  !!    pkgName      : package name
  !!    varName      : variable name
  !!    auxOrBnd     : optional 'AUX' or 'BND' keyword
  !<
  function var_timeseries(tsManager, pkgName, varName, auxOrBnd) result(tsexists)
    ! -- return
    logical :: tsexists
    ! -- dummy
    type(TimeSeriesManagerType), intent(inout) :: tsManager
    character(len=*), intent(in) :: pkgName
    character(len=*), intent(in) :: varName
    character(len=3), intent(in), optional :: auxOrBnd
    ! -- local
    character(len=3) :: ctstype
    integer(I4B) :: i
    integer(I4B) :: nlinks
    type(TimeSeriesLinkType), pointer :: tslTemp => null()
    !
    ! -- process optional variables
    if (present(auxOrBnd)) then
      ctstype = auxOrBnd
    else
      ctstype = 'BND'
    end if
    !
    ! -- initialize the return variable and the number of timeseries links
    tsexists = .FALSE.
    nlinks = tsManager%CountLinks(ctstype)
    !
    ! -- determine if link exists
    csearchlinks: do i = 1, nlinks
      tslTemp => tsManager%GetLink(ctstype, i)
      if (tslTemp%PackageName == pkgName) then
        !
        ! -- Check varName against Text member of link
        if (same_word(tslTemp%Text, varName)) then
          tsexists = .TRUE.
          exit csearchlinks
        end if
      end if
    end do csearchlinks
  end function var_timeseries

end module TimeSeriesManagerModule
