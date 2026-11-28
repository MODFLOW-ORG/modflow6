module TimeArraySeriesManagerModule

  use KindModule, only: DP, I4B
  use SimVariablesModule, only: errmsg
  use ConstantsModule, only: DZERO, LENTIMESERIESNAME, LINELENGTH, &
                             LENHUGELINE, LENMODELNAME, LENPACKAGENAME
  use ListModule, only: ListType
  use SimModule, only: store_error, store_error_unit
  use TdisModule, only: delt, totimc, kper, kstp
  use TimeArraySeriesLinkModule, only: TimeArraySeriesLinkType, &
                                       ConstructTimeArraySeriesLink, &
                                       GetTimeArraySeriesLinkFromList, &
                                       AddTimeArraySeriesLinkToList
  use TimeArraySeriesModule, only: TimeArraySeriesType
  use BaseDisModule, only: DisBaseType

  implicit none

  private
  public :: TimeArraySeriesManagerType, tasmanager_cr

  type TimeArraySeriesManagerType
    ! -- Public members
    integer(I4B), public :: iout = 0 ! output unit num
    class(DisBaseType), pointer :: dis => null() ! pointer to dis
    character(len=LENMODELNAME) :: modelname
    logical, public :: fullgrid = .false. ! series arrays for all user nodes
    ! -- Private members
    type(ListType), pointer, private :: boundTasLinks => null() ! list of TAS links
    character(len=LINELENGTH), allocatable, dimension(:) :: tasfiles ! list of TA file names
    type(TimeArraySeriesType), dimension(:), pointer, contiguous :: taslist ! array of TA pointers
    character(len=LENTIMESERIESNAME), allocatable, dimension(:) :: tasnames ! array of TA names

  contains

    ! -- Public procedures
    procedure, public :: tasmanager_df
    procedure, public :: ad => tasmgr_ad
    procedure, public :: da => tasmgr_da
    procedure, public :: add_tasfile
    procedure, public :: CountLinks
    procedure, public :: GetLink
    procedure, public :: MakeTasLink
    procedure, public :: Reset
    procedure, public :: remove_links
    procedure, public :: get_series
    ! -- Private procedures
    procedure, private :: tasmgr_add_link
    procedure, private :: tasmgr_convert_flux
  end type TimeArraySeriesManagerType

contains

! -- Type-bound procedures of TimeArraySeriesManagerType

  ! -- Public procedures

  !> @brief Create the time-array series manager
  !<
  subroutine tasmanager_cr(this, dis, modelname, iout, fullgrid)
    ! -- dummy
    type(TimeArraySeriesManagerType) :: this
    class(DisBaseType), pointer, optional :: dis
    character(len=*), intent(in) :: modelname
    integer(I4B), intent(in) :: iout
    logical, optional, intent(in) :: fullgrid !< series arrays for all user nodes
    !
    if (present(dis)) then
      this%dis => dis
    end if
    if (present(fullgrid)) this%fullgrid = fullgrid
    !
    this%modelname = modelname
    this%iout = iout
    allocate (this%boundTasLinks)
    allocate (this%tasfiles(0))
  end subroutine tasmanager_cr

  !> @brief Define the time-array series manager
  !<
  subroutine tasmanager_df(this)
    ! -- dummy
    class(TimeArraySeriesManagerType) :: this
    ! -- local
    type(TimeArraySeriesType), pointer :: tasptr => null()
    integer(I4B) :: nfiles
    integer(I4B) :: i
    !
    ! -- determine how many tasfiles.  This is the number of time array series
    !    so allocate arrays to store them
    nfiles = size(this%tasfiles)
    allocate (this%taslist(nfiles))
    allocate (this%tasnames(nfiles))
    !
    ! -- Setup a time array series for each file specified
    do i = 1, nfiles
      tasptr => this%taslist(i)
      call tasptr%tas_init(this%tasfiles(i), this%modelname, &
                           this%iout, this%tasnames(i), &
                           fullgrid=this%fullgrid)
    end do
  end subroutine tasmanager_df

  !> @brief Time step (or subtime step) advance.
  !!
  !! Call this each time step or subtime step.
  !<
  subroutine tasmgr_ad(this)
    ! -- dummy
    class(TimeArraySeriesManagerType) :: this
    ! -- local
    type(TimeArraySeriesLinkType), pointer :: tasLink => null()
    type(TimeArraySeriesType), pointer :: timearrayseries => null()
    integer(I4B) :: i, j, nlinks, nvals, isize1, isize2, inunit
    real(DP) :: begintime, endtime
    ! -- formats
    character(len=*), parameter :: fmt5 = &
      "(/,'Time-array-series controlled arrays in stress period ', &
      &i0, ', time step ', i0, ':')"
10  format('"', a, '" package: ', a, ' array obtained from time-array series "', &
           a, '"')
    !
    ! -- Initialize time variables
    begintime = totimc
    endtime = begintime + delt
    !
    ! -- Iterate through boundtaslinks and update specified
    !    array with array of average values obtained from
    !    appropriate time series.
    if (associated(this%boundTasLinks)) then
      nlinks = this%boundTasLinks%Count()
      do i = 1, nlinks
        tasLink => GetTimeArraySeriesLinkFromList(this%boundTasLinks, i)
        if (tasLink%Iprpak == 1 .and. i == 1) then
          write (this%iout, fmt5) kper, kstp
        end if
        if (tasLink%UseDefaultProc) then
          timearrayseries => tasLink%timeArraySeries
          !
          ! -- Fill the package array with integrated values
          nvals = link_nvals(tasLink)
          if (associated(tasLink%imap) .and. tasLink%keep_nodata) then
            ! -- record the elements set, keeping elements where the
            !    series is DNODATA
            if (allocated(tasLink%written)) then
              if (size(tasLink%written) /= nvals) deallocate (tasLink%written)
            end if
            if (.not. allocated(tasLink%written)) then
              allocate (tasLink%written(nvals))
            end if
            call timearrayseries%GetAverageValues( &
              nvals, tasLink%BndArray(1:nvals), begintime, endtime, &
              tasLink%imap(1:nvals), keep=.true., written=tasLink%written)
          else if (associated(tasLink%imap)) then
            ! -- package array elements map to series elements imap
            call timearrayseries%GetAverageValues( &
              nvals, tasLink%BndArray(1:nvals), begintime, endtime, &
              tasLink%imap(1:nvals))
          else
            call timearrayseries%GetAverageValues(nvals, tasLink%BndArray, &
                                                  begintime, endtime)
          end if
          !
          ! -- If conversion from flux to flow is required, multiply by cell area
          if (tasLink%ConvertFlux) then
            call this%tasmgr_convert_flux(tasLink)
          end if
          !
          ! -- If PRINT_INPUT is specified, write information
          !    regarding source of time-array series data
          if (tasLink%Iprpak == 1) then
            write (this%iout, 10) trim(tasLink%PackageName), &
              trim(tasLink%Text), &
              trim(tasLink%timeArraySeries%Name)
          end if
        end if
        if (i == nlinks) then
          write (this%iout, '()')
        end if
      end do
      !
      ! -- Now that all array values have been substituted, can now multiply
      !    an array by a multiplier array
      do i = 1, nlinks
        tasLink => GetTimeArraySeriesLinkFromList(this%boundTasLinks, i)
        if (tasLink%UseDefaultProc) then
          if (associated(tasLink%RMultArray)) then
            nvals = link_nvals(tasLink)
            isize1 = size(tasLink%BndArray)
            isize2 = size(tasLink%RMultArray)
            if (isize1 == isize2 .and. nvals <= isize1) then
              do j = 1, nvals
                tasLink%BndArray(j) = tasLink%BndArray(j) * tasLink%RMultArray(j)
              end do
            else
              errmsg = 'Size mismatch between boundary and multiplier arrays'// &
                       ' using time-array series: '// &
                       trim(tasLink%TimeArraySeries%Name)
              call store_error(errmsg)
              inunit = tasLink%TimeArraySeries%GetInunit()
              call store_error_unit(inunit)
            end if
          end if
        end if
      end do
    end if
  end subroutine tasmgr_ad

  !> @brief Deallocate
  !<
  subroutine tasmgr_da(this)
    ! -- dummy
    class(TimeArraySeriesManagerType) :: this
    ! -- local
    integer :: i, n
    type(TimeArraySeriesLinkType), pointer :: tasLink => null()
    !
    ! -- Deallocate contents of each TimeArraySeriesType object in list
    !    of time-array series links.
    n = this%boundTasLinks%Count()
    do i = 1, n
      tasLink => GetTimeArraySeriesLinkFromList(this%boundTasLinks, i)
      call tasLink%da()
    end do
    !
    ! -- Go through and deallocate individual time array series
    do i = 1, size(this%taslist)
      call this%taslist(i)%da()
    end do
    !
    ! -- Deallocate the list of time-array series links.
    call this%boundTasLinks%Clear(.true.)
    deallocate (this%boundTasLinks)
    deallocate (this%tasfiles)
    !
    ! -- Deallocate the time array series
    deallocate (this%taslist)
    deallocate (this%tasnames)
    !
    ! -- nullify pointers
    this%dis => null()
    this%boundTasLinks => null()
  end subroutine tasmgr_da

  !> @brief Add a time-array series file
  !<
  subroutine add_tasfile(this, fname)
    ! -- modules
    use ArrayHandlersModule, only: ExpandArray
    ! -- dummy
    class(TimeArraySeriesManagerType) :: this
    character(len=*), intent(in) :: fname
    ! -- local
    integer(I4B) :: indx
    !
    call ExpandArray(this%tasfiles, 1)
    indx = size(this%tasfiles)
    this%tasfiles(indx) = fname
  end subroutine add_tasfile

  !> @brief Zero out arrays that are represented with time series
  !!
  !! Delete all existing links from time array series to package arrays as they
  !! will need to be created with a new BEGIN PERIOD block
  !<
  subroutine Reset(this, pkgName)
    ! -- dummy
    class(TimeArraySeriesManagerType) :: this
    character(len=*), intent(in) :: pkgName
    ! -- local
    integer(I4B) :: i, j, nlinks
    type(TimeArraySeriesLinkType), pointer :: taslink
    !
    ! -- Reassign all linked elements to zero
    nlinks = this%boundTasLinks%Count()
    do i = 1, nlinks
      taslink => GetTimeArraySeriesLinkFromList(this%boundTasLinks, i)
      if (associated(taslink)) then
        do j = 1, size(taslink%BndArray)
          taslink%BndArray(j) = DZERO
        end do
      end if
    end do
    !
    ! -- Delete all existing time array links
    if (associated(this%boundTasLinks)) then
      ! Deallocate and remove all links belonging to package
      nlinks = this%boundTasLinks%Count()
      do i = nlinks, 1, -1
        taslink => GetTimeArraySeriesLinkFromList(this%boundTasLinks, i)
        if (associated(taslink)) then
          call taslink%da()
          call this%boundTasLinks%RemoveNode(i, .true.)
        end if
      end do
    end if
  end subroutine Reset

  !> @brief Remove time-array series links for a package
  !!
  !! Linked array values are not changed.
  !<
  subroutine remove_links(this, pkgName)
    ! -- modules
    use InputOutputModule, only: UPCASE
    ! -- dummy
    class(TimeArraySeriesManagerType) :: this
    character(len=*), intent(in) :: pkgName
    ! -- local
    integer(I4B) :: i
    character(len=LENPACKAGENAME) :: pkgNameTemp
    type(TimeArraySeriesLinkType), pointer :: taslink
    !
    if (.not. associated(this%boundTasLinks)) return
    !
    ! -- link package names are stored as all caps
    pkgNameTemp = pkgName
    call UPCASE(pkgNameTemp)
    !
    do i = this%boundTasLinks%Count(), 1, -1
      taslink => GetTimeArraySeriesLinkFromList(this%boundTasLinks, i)
      if (.not. associated(taslink)) cycle
      if (taslink%PackageName /= pkgNameTemp) cycle
      call taslink%da()
      call this%boundTasLinks%RemoveNode(i, .true.)
    end do
  end subroutine remove_links

  !> @brief Return a time-array series by name
  !<
  function get_series(this, tasName, inunit) result(tasptr)
    ! -- dummy
    class(TimeArraySeriesManagerType) :: this
    character(len=*), intent(in) :: tasName
    integer(I4B), intent(in) :: inunit !< unit for error reporting
    ! -- return
    type(TimeArraySeriesType), pointer :: tasptr
    ! -- local
    integer(I4B) :: i
    !
    tasptr => null()
    do i = 1, size(this%tasnames)
      if (this%tasnames(i) == tasName) then
        tasptr => this%taslist(i)
        exit
      end if
    end do
    if (.not. associated(tasptr)) then
      errmsg = 'Error: Time-array series "'//trim(tasName)//'" not found.'
      call store_error(errmsg)
      call store_error_unit(inunit)
    end if
  end function get_series

  !> @brief Make link from time-array series to package array
  !<
  subroutine MakeTasLink(this, pkgName, bndArray, iprpak, &
                         tasName, text, convertFlux, nodelist, inunit, &
                         imap, nmap, keep_nodata)
    ! -- dummy
    class(TimeArraySeriesManagerType) :: this
    character(len=*), intent(in) :: pkgName
    real(DP), dimension(:), pointer :: bndArray
    integer(I4B), intent(in) :: iprpak
    character(len=*), intent(in) :: tasName
    character(len=*), intent(in) :: text
    logical, intent(in) :: convertFlux
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: nodelist
    integer(I4B), intent(in) :: inunit
    integer(I4B), dimension(:), pointer, contiguous, optional, &
      intent(in) :: imap !< series element for each bndArray element
    integer(I4B), pointer, optional, intent(in) :: nmap !< elements in use, required with imap
    logical, optional, intent(in) :: keep_nodata !< with imap, keep values where series is DNODATA
    ! -- local
    type(TimeArraySeriesLinkType), pointer :: newTasLink
    type(TimeArraySeriesType), pointer :: tasptr => null()
    !
    ! -- Find the time array series
    tasptr => this%get_series(tasName, inunit)
    !
    ! -- Construct a time-array series link
    newTasLink => null()
    call ConstructTimeArraySeriesLink(newTasLink, tasptr, &
                                      pkgName, bndArray, iprpak, &
                                      text)
    newTasLink%ConvertFlux = convertFlux
    newTasLink%nodelist => nodelist
    if (present(imap)) newTasLink%imap => imap
    if (present(nmap)) newTasLink%nmap => nmap
    if (present(keep_nodata)) newTasLink%keep_nodata = keep_nodata
    !
    ! -- Add link to list of links
    call this%tasmgr_add_link(newTasLink)
  end subroutine MakeTasLink

  !> @brief Get link from the boundtaslinks list
  !<
  function GetLink(this, indx) result(tasLink)
    ! -- dummy
    class(TimeArraySeriesManagerType) :: this
    integer(I4B), intent(in) :: indx
    ! -- return
    type(TimeArraySeriesLinkType), pointer :: tasLink
    !
    tasLink => null()
    !
    if (associated(this%boundTasLinks)) then
      tasLink => GetTimeArraySeriesLinkFromList(this%boundTasLinks, indx)
    end if
  end function GetLink

  !> @brief Count number of links in the boundtaslinks list
  !<
  function CountLinks(this)
    ! -- return
    integer(I4B) :: CountLinks
    ! -- dummy
    class(TimeArraySeriesManagerType) :: this
    !
    if (associated(this%boundtaslinks)) then
      CountLinks = this%boundTasLinks%Count()
    else
      CountLinks = 0
    end if
  end function CountLinks

  ! -- Private procedures

  !> @brief Convert the array from a flux to a flow rate by multiplying by the
  !! cell area
  !<
  subroutine tasmgr_convert_flux(this, tasLink)
    ! -- dummy
    class(TimeArraySeriesManagerType) :: this
    type(TimeArraySeriesLinkType), pointer, intent(inout) :: tasLink
    ! -- local
    integer(I4B) :: i, n, noder
    real(DP) :: area
    !
    if (.not. (associated(this%dis) .and. &
               associated(tasLink%nodelist))) then
      errmsg = 'Programming error. Cannot convert flux. Verify that '&
               &'a valid DIS instance and nodelist were provided.'
      call store_error(errmsg)
      call store_error_unit(tasLink%TimeArraySeries%GetInunit())
    end if
    !
    n = size(tasLink%BndArray)
    do i = 1, n
      noder = tasLink%nodelist(i)
      if (noder > 0) then
        area = this%dis%get_area(noder)
        tasLink%BndArray(i) = tasLink%BndArray(i) * area
      end if
    end do
  end subroutine tasmgr_convert_flux

  !> @brief Return the number of package array elements a link fills
  !<
  function link_nvals(tasLink) result(nvals)
    ! -- dummy
    type(TimeArraySeriesLinkType), pointer, intent(in) :: tasLink
    ! -- return
    integer(I4B) :: nvals
    !
    if (associated(tasLink%imap)) then
      nvals = tasLink%nmap
    else
      nvals = size(tasLink%BndArray)
    end if
  end function link_nvals

  !> @brief Add a time arrays series link
  !<
  subroutine tasmgr_add_link(this, tasLink)
    ! -- dummy
    class(TimeArraySeriesManagerType) :: this
    type(TimeArraySeriesLinkType), pointer :: tasLink
    ! -- local
    !
    call AddTimeArraySeriesLinkToList(this%boundTasLinks, tasLink)
  end subroutine tasmgr_add_link

end module TimeArraySeriesManagerModule

