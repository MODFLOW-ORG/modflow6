!> @brief This module contains the TspSpc Module
!!
!! This module contains the code for reading and storing a
!! generic input file of source and sink concentrations or
!! temperatures.
!<
module TspSpcModule

  use KindModule, only: DP, LGP, I4B
  use ConstantsModule, only: LENPACKAGENAME, LENMODELNAME, &
                             LENMEMPATH, DZERO, LENFTYPE, &
                             LINELENGTH, TABLEFT, TABCENTER, &
                             LENVARNAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, count_errors
  use MemoryHelperModule, only: create_mem_path
  use BlockParserModule, only: BlockParserType
  use BaseDisModule, only: DisBaseType
  use TimeSeriesManagerModule, only: TimeSeriesManagerType, tsmanager_cr
  use TimeArraySeriesManagerModule, only: TimeArraySeriesManagerType, &
                                          tasmanager_cr
  use InputOutputModule, only: str_pad_left
  use TableModule, only: TableType, table_cr

  implicit none
  private
  public :: TspSpcType

  character(len=LENFTYPE) :: ftype = 'SPC'
  character(len=LENPACKAGENAME) :: text = 'STRESS PACK COMP'

  !> @brief Derived type for managing SPC input
  !!
  !! This derived type will read and process an SPC input file,
  !! make time series interpolations, and provide concentrations or
  !! temperatures to the SSM package that correspond to an individual
  !! GWF stress package.
  !<
  type :: TspSpcType

    character(len=LENMODELNAME) :: name_model = '' !< the name of the model that contains this package
    character(len=LENPACKAGENAME) :: packName = '' !< name of the package
    character(len=LENPACKAGENAME) :: packNameFlow = '' !< name of the corresponding flow package
    character(len=LENVARNAME) :: depvarname = '' !< name of the dependent variable (CONCENTRATION or TEMPERATURE)
    character(len=LENMEMPATH) :: memoryPath = '' !< the location in the memory manager where the variables are stored
    integer(I4B), pointer :: id => null() !< id number for this spc package
    integer(I4B), pointer :: inunit => null() !< unit number for input
    integer(I4B), pointer :: iout => null() !< unit number for output
    integer(I4B), pointer :: maxbound => null() !< length of dblvec
    integer(I4B), pointer :: ionper => null() !< stress period for next data
    integer(I4B), pointer :: lastonper => null() !< last value of ionper (for checking)
    integer(I4B), pointer :: iprpak => null() !< flag for printing input
    logical(LGP), pointer :: readasarrays => null() !< flag for reading concentrations as an array
    real(DP), dimension(:), pointer, contiguous :: dblvec => null() !< vector of floats read from file
    class(DisBaseType), pointer :: dis => null() !< model discretization object
    type(BlockParserType) :: parser !< parser object for reading blocks of information
    type(TimeSeriesManagerType), pointer :: TsManager => null() !< time series manager
    type(TimeArraySeriesManagerType), pointer :: TasManager => null() !< time array series manager
    type(TableType), pointer :: inputtab => null() !< input table object

  contains

    procedure :: initialize
    procedure :: allocate_scalars
    procedure :: read_options
    procedure :: read_dimensions
    procedure :: allocate_arrays
    procedure :: get_value
    procedure :: set_value
    procedure :: spc_rp
    procedure :: spc_rp_list
    procedure :: spc_rp_array
    procedure :: spc_ad
    procedure :: spc_da
    procedure :: read_check_ionper
    procedure :: check_flow_package

  end type TspSpcType

contains

  !> @ brief Initialize the SPC type
  !!
  !! Initialize the SPC object by setting up the parser,
  !! and time series manager, reading options and dimensions,
  !! and allocating memory.
  !!
  !<
  subroutine initialize(this, dis, id, inunit, iout, name_model, packNameFlow, &
                        dvn)
    ! -- dummy variables
    class(TspSpcType) :: this !<  TspSpcType
    class(DisBaseType), pointer, intent(in) :: dis !<  discretization package
    integer(I4B), intent(in) :: id !<  id number for this spc package
    integer(I4B), intent(in) :: inunit !<  unit number for input
    integer(I4B), intent(in) :: iout !<  unit number for output
    character(len=*), intent(in) :: name_model !<  character string containing model name
    character(len=*), intent(in) :: packNameflow !<  character string containing name of corresponding flow package
    character(len=*), intent(in) :: dvn !<  dependent variable name (CONCENTRATION or TEMPERATURE)
    ! -- local
    !
    ! -- construct the memory path
    write (this%packName, '(a, i0)') 'SPC'//'-', id
    this%name_model = name_model
    this%memoryPath = create_mem_path(this%name_model, this%packName)
    !
    ! -- allocate scalar variables
    call this%allocate_scalars()
    !
    ! -- assign member values
    this%id = id
    this%inunit = inunit
    this%iout = iout
    this%packNameFlow = packNameFlow
    this%depvarname = dvn
    !
    ! -- set pointers
    this%dis => dis
    !
    ! -- Setup the parser
    call this%parser%Initialize(this%inunit, this%iout)
    !
    ! -- Setup the time series manager
    call tsmanager_cr(this%TsManager, this%iout)
    call tasmanager_cr(this%TasManager, dis, name_model, this%iout)
    !
    ! -- read options
    call this%read_options()
    !
    ! -- read dimensions
    if (this%readasarrays) then
      this%maxbound = this%dis%get_ncpl()
    else
      call this%read_dimensions()
    end if
    !
    ! -- allocate arrays
    call this%allocate_arrays()
    !
    ! -- Now that time series are read, call define
    call this%tsmanager%tsmanager_df()
    call this%tasmanager%tasmanager_df()
  end subroutine initialize

  !> @ brief Allocate package scalars
  !!
  !!  Allocate and initialize package scalars.
  !!
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy variables
    class(TspSpcType) :: this !< TspSpcType object
    !
    ! -- allocate scalars in memory manager
    call mem_allocate(this%id, 'ID', this%memoryPath)
    call mem_allocate(this%inunit, 'INUNIT', this%memoryPath)
    call mem_allocate(this%iout, 'IOUT', this%memoryPath)
    call mem_allocate(this%maxbound, 'MAXBOUND', this%memoryPath)
    call mem_allocate(this%ionper, 'IONPER', this%memoryPath)
    call mem_allocate(this%lastonper, 'LASTONPER', this%memoryPath)
    call mem_allocate(this%iprpak, 'IPRPAK', this%memoryPath)
    call mem_allocate(this%readasarrays, 'READASARRAYS', this%memoryPath)
    !
    ! -- allocate special derived types
    allocate (this%TsManager)
    allocate (this%TasManager)
    !
    ! -- initialize
    this%id = 0
    this%inunit = 0
    this%iout = 0
    this%maxbound = 0
    this%ionper = 0
    this%lastonper = 0
    this%iprpak = 0
    this%readasarrays = .false.
  end subroutine allocate_scalars

  !> @ brief Read options for package
  !!
  !!  Read options for this package.
  !!
  !<
  subroutine read_options(this)
    ! -- modules
    ! -- dummy
    class(TspSpcType) :: this
    ! -- local
    character(len=LINELENGTH) :: keyword, fname
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    ! -- formats
    character(len=*), parameter :: fmtiprpak = &
      &"(4x,'SPC INFORMATION WILL BE PRINTED TO LISTING FILE.')"
    character(len=*), parameter :: fmtreadasarrays = &
      "(4x,'SPC INFORMATION WILL BE READ AS ARRAYS RATHER THAN IN LIST &
      &FORMAT.')"
    character(len=*), parameter :: fmtts = &
      &"(4x, 'TIME-SERIES DATA WILL BE READ FROM FILE: ', a)"
    character(len=*), parameter :: fmttas = &
      &"(4x, 'TIME-ARRAY SERIES DATA WILL BE READ FROM FILE: ', a)"
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, blockRequired=.false., &
                              supportOpenClose=.true.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING SPC OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('PRINT_INPUT')
          this%iprpak = 1
          write (this%iout, fmtiprpak)
        case ('READASARRAYS')
          this%readasarrays = .true.
          write (this%iout, fmtreadasarrays)
        case ('TS6')
          call this%parser%GetStringCaps(keyword)
          if (trim(adjustl(keyword)) /= 'FILEIN') then
            errmsg = 'TS6 keyword must be followed by "FILEIN" '// &
                     'then by filename.'
            call store_error(errmsg)
          end if
          call this%parser%GetString(fname)
          write (this%iout, fmtts) trim(fname)
          call this%TsManager%add_tsfile(fname, this%inunit)
        case ('TAS6')
          call this%parser%GetStringCaps(keyword)
          if (trim(adjustl(keyword)) /= 'FILEIN') then
            errmsg = 'TAS6 keyword must be followed by "FILEIN" '// &
                     'then by filename.'
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetString(fname)
          write (this%iout, fmttas) trim(fname)
          call this%TasManager%add_tasfile(fname)
        case default
          write (errmsg, '(a,a)') 'Unknown SPC option: ', trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (this%iout, '(1x,a)') 'END OF SPC OPTIONS'
    end if
  end subroutine read_options

  !> @ brief Read dimensions for package
  !!
  !!  Read dimensions for this package.
  !!
  !<
  subroutine read_dimensions(this)
    ! -- dummy variables
    class(TspSpcType), intent(inout) :: this !< TspSpcType object
    ! -- local variables
    character(len=LINELENGTH) :: keyword
    logical(LGP) :: isfound
    logical(LGP) :: endOfBlock
    integer(I4B) :: ierr
    !
    ! -- get dimensions block
    call this%parser%GetBlock('DIMENSIONS', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse dimensions block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(text))// &
        ' DIMENSIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('MAXBOUND')
          this%maxbound = this%parser%GetInteger()
          write (this%iout, '(4x,a,i7)') 'MAXBOUND = ', this%maxbound
        case default
          write (errmsg, '(a,3(1x,a))') &
            'Unknown', trim(text), 'dimension:', trim(keyword)
          call store_error(errmsg)
        end select
      end do
      !
      write (this%iout, '(1x,a)') 'END OF '//trim(adjustl(text))//' DIMENSIONS'
    else
      call store_error('Required DIMENSIONS block not found.')
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- verify dimensions were set
    if (this%maxbound <= 0) then
      write (errmsg, '(a)') 'MAXBOUND must be an integer greater than zero.'
      call store_error(errmsg)
    end if
    !
    ! -- terminate if there are errors
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine read_dimensions

  !> @ brief Allocate package arrays
  !!
  !!  Allocate and initialize package arrays.
  !!
  !<
  subroutine allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy variables
    class(TspSpcType) :: this !< TspSpcType object
    ! -- local
    integer(I4B) :: i
    !
    ! -- allocate array
    call mem_allocate(this%dblvec, this%maxbound, 'DBLVEC', this%memoryPath)
    !
    ! -- initialize dblvec to zero
    do i = 1, this%maxbound
      this%dblvec(i) = DZERO
    end do
  end subroutine allocate_arrays

  !> @ brief Get the data value from this package
  !!
  !!  Get the floating point value from the dblvec array.
  !!
  !<
  function get_value(this, ientry, nbound_flow) result(value)
    class(TspSpcType) :: this !< TspSpcType object
    integer(I4B), intent(in) :: ientry !< index of the data to return
    integer(I4B), intent(in) :: nbound_flow !< size of bound list in flow package
    real(DP) :: value
    integer(I4B) :: nu
    if (this%readasarrays) then
      ! Special handling for reduced grids and readasarrays
      ! if flow and transport are in the same simulation, then
      ! ientry is a user node number and it corresponds to the
      ! correct position in the dblvec array.  But if flow and
      ! transport are not in the same simulation, then ientry is
      ! a reduced node number, because the list of flows in the
      ! budget file do not include idomain < 1 entries. In this
      ! case, ientry must be converted to a user node number so
      ! that it corresponds to a user array, which includes
      ! idomain < 1 values.
      if (nbound_flow == this%maxbound) then
        ! flow and transport are in the same simulation or there
        ! are no idomain < 1 cells.
        value = this%dblvec(ientry)
      else
        ! This identifies case where flow and transport must be
        ! in a separate simulation, because nbound_flow is not
        ! the same as this%maxbound.  Under these conditions, we
        ! must assume that ientry corresponds to a flow list that
        ! would be of size ncpl if flow and transport were in the
        ! same simulation, but because boundary cells with
        ! idomain < 1 are not written to binary budget file, the
        ! list size is smaller.
        nu = this%dis%get_nodeuser(ientry)
        value = this%dblvec(nu)
      end if
    else
      value = this%dblvec(ientry)
    end if
  end function get_value

  !> @ brief Read and prepare
  !!
  !!  Read and prepare the period data block and fill dblvec
  !!  if the next period block corresponds to this time step.
  !!
  !<
  subroutine spc_rp(this)
    ! -- modules
    use TdisModule, only: kper, nper
    ! -- dummy
    class(TspSpcType), intent(inout) :: this !< TspSpcType object
    ! -- local
    character(len=LINELENGTH) :: line
    logical :: isfound
    integer(I4B) :: ierr
    ! -- formats
    character(len=*), parameter :: fmtblkerr = &
      &"('Looking for BEGIN PERIOD iper.  Found ', a, ' instead.')"
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
          call store_error(errmsg, terminate=.TRUE.)
        end if
      end if
    end if
    !
    ! -- Read data if ionper == kper
    if (this%ionper == kper) then
      !
      ! -- Remove all time-series and time-array-series links associated with
      !    this package.
      !    Do not reset as we are using a "settings" approach here in which the
      !    settings remain the same until the user changes them.
      ! call this%TsManager%Reset(this%packName)
      call this%TasManager%Reset(this%packName)
      if (this%readasarrays) then
        call this%spc_rp_array(line)
      else
        call this%spc_rp_list()
      end if
      !
      ! -- using data from the last stress period
    else
      write (this%iout, fmtlsp) trim(ftype)
    end if
    !
    ! -- write summary of maw well stress period error messages
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine spc_rp

  !> @ brief spc_rp_list
  !!
  !!  Read the stress period data in list format
  !!
  !<
  subroutine spc_rp_list(this)
    ! -- modules
    use TdisModule, only: kper
    ! -- dummy
    class(TspSpcType), intent(inout) :: this !< TspSpcType object
    ! -- local
    character(len=LINELENGTH) :: line
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: tabletext
    logical :: endOfBlock
    integer(I4B) :: ival
    !
    !
    ! -- setup table for period data
    if (this%iprpak /= 0) then
      !
      ! -- reset the input table object
      title = trim(adjustl(text))//' PACKAGE ('// &
              'SPC'//') DATA FOR PERIOD'
      write (title, '(a,1x,i6)') trim(adjustl(title)), kper
      call table_cr(this%inputtab, ftype, title)
      call this%inputtab%table_df(1, 3, this%iout, finalize=.FALSE.)
      tabletext = 'NUMBER'
      call this%inputtab%initialize_column(tabletext, 10, alignment=TABCENTER)
      tabletext = 'DATA TYPE'
      call this%inputtab%initialize_column(tabletext, 20, alignment=TABLEFT)
      write (tabletext, '(a,1x,i6)') 'VALUE'
      call this%inputtab%initialize_column(tabletext, 15, alignment=TABCENTER)
    end if
    !
    ! -- read data
    do
      call this%parser%GetNextLine(endOfBlock)
      if (endOfBlock) exit

      ival = this%parser%GetInteger()
      if (ival < 1 .or. ival > this%maxbound) then
        write (errmsg, '(2(a,1x),i0,a)') &
          'IVAL must be greater than 0 and', &
          'less than or equal to ', this%maxbound, '.'
        call store_error(errmsg)
        cycle
      end if
      !
      ! -- set stress period data
      call this%set_value(ival)
      !
      ! -- write line to table
      if (this%iprpak /= 0) then
        call this%parser%GetCurrentLine(line)
        call this%inputtab%line_to_columns(line)
      end if
    end do
    !
    ! -- finalize the table
    if (this%iprpak /= 0) then
      call this%inputtab%finalize_table()
    end if
  end subroutine spc_rp_list

  !> @ brief spc_rp_array
  !!
  !!  Read the stress period data in array format
  !!
  !<
  subroutine spc_rp_array(this, line)
    use ConstantsModule, only: LENTIMESERIESNAME, LENANAME
    use SimModule, only: store_error
    use ArrayHandlersModule, only: ifind
    ! -- dummy
    class(TspSpcType), intent(inout) :: this !< TspSpcType object
    character(len=LINELENGTH), intent(inout) :: line
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: ncolbnd
    integer(I4B) :: jauxcol, ivarsread
    integer(I4B), dimension(:), allocatable, target :: nodelist
    character(len=LENTIMESERIESNAME) :: tasName
    character(len=LENANAME) :: aname
    character(len=LINELENGTH) :: keyword
    logical :: endOfBlock
    logical :: convertFlux
    !
    ! -- these time array series pointers need to be non-contiguous
    !    because a slice of bound is passed
    real(DP), dimension(:), pointer :: bndArrayPtr => null()
    !
    write (aname, '(a)') str_pad_left(this%depvarname, LENANAME)
    !
    ! -- Initialize
    jauxcol = 0
    ivarsread = 0
    ncolbnd = 1
    allocate (nodelist(this%maxbound))
    do n = 1, size(nodelist)
      nodelist(n) = n
    end do
    !
    ! -- Read CONCENTRATION variables as arrays
    do
      call this%parser%GetNextLine(endOfBlock)
      if (endOfBlock) exit
      call this%parser%GetStringCaps(keyword)
      !
      ! -- Parse the keywords
      select case (keyword)
      case ('CONCENTRATION', 'TEMPERATURE')
        !
        ! -- Look for keyword TIMEARRAYSERIES and time-array series
        !    name on line, following RECHARGE
        call this%parser%GetStringCaps(keyword)
        if (keyword == 'TIMEARRAYSERIES') then
          ! -- Get time-array series name
          call this%parser%GetStringCaps(tasName)
          bndArrayPtr => this%dblvec(:)
          ! Make a time-array-series link and add it to the list of links
          ! contained in the TimeArraySeriesManagerType object.
          convertflux = .false.
          call this%TasManager%MakeTasLink(this%packName, bndArrayPtr, &
                                           this%iprpak, tasName, &
                                           this%depvarname, &
                                           convertFlux, nodelist, &
                                           this%parser%iuactive)
        else
          !
          ! -- Read the concentration array
          call this%dis%read_layer_array(nodelist, this%dblvec, ncolbnd, &
                                         this%maxbound, 1, aname, &
                                         this%parser%iuactive, this%iout)
        end if
        !
      case default
        call store_error('Looking for component name, either CONCENTRATION &
                         &or TEMPERATURE.  Found: '//trim(line))
        call this%parser%StoreErrorUnit()
      end select

    end do
  end subroutine spc_rp_array

  !> @ brief Advance
  !!
  !!  Call the advance method on the time series so that new values
  !!  are interpolated and entered into dblvec
  !!
  !<
  subroutine spc_ad(this, nbound_flowpack, budtxt)
    ! -- modules
    ! -- dummy
    class(TspSpcType), intent(inout) :: this !< TspSpcType object
    integer(I4B), intent(in) :: nbound_flowpack
    character(len=*), intent(in) :: budtxt
    ! -- local
    !
    ! -- Advance the time series
    call this%TsManager%ad()
    call this%TasManager%ad()
    !
    ! -- Check flow package consistency
    call this%check_flow_package(nbound_flowpack, budtxt)
  end subroutine spc_ad

  !> @ brief Deallocate variables
  !!
  !!  Deallocate and nullify package variables.
  !!
  !<
  subroutine spc_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy variables
    class(TspSpcType) :: this !< TspSpcType object
    !
    ! -- deallocate arrays in memory manager
    call mem_deallocate(this%dblvec)
    !
    ! -- deallocate scalars in memory manager
    call mem_deallocate(this%id)
    call mem_deallocate(this%inunit)
    call mem_deallocate(this%iout)
    call mem_deallocate(this%maxbound)
    call mem_deallocate(this%ionper)
    call mem_deallocate(this%lastonper)
    call mem_deallocate(this%iprpak)
    call mem_deallocate(this%readasarrays)
    !
    ! -- deallocate derived types
    call this%TsManager%da()
    deallocate (this%TsManager)
    nullify (this%TsManager)
  end subroutine spc_da

  !> @ brief Check ionper
  !!
  !!  Generic method to read and check ionperiod, which is used to determine
  !!  if new period data should be read from the input file. The check of
  !!  ionperiod also makes sure periods are increasing in subsequent period
  !!  data blocks.  Copied from NumericalPackage
  !!
  !<
  subroutine read_check_ionper(this)
    ! -- modules
    use TdisModule, only: kper
    ! -- dummy variables
    class(TspSpcType), intent(inout) :: this !< TspSpcType object
    !
    ! -- save last value and read period number
    this%lastonper = this%ionper
    this%ionper = this%parser%GetInteger()
    !
    ! -- make check
    if (this%ionper <= this%lastonper) then
      write (errmsg, '(a, i0, a, i0, a, i0, a)') &
        'Error in stress period ', kper, &
        '. Period numbers not increasing.  Found ', this%ionper, &
        ' but last period block was assigned ', this%lastonper, '.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
  end subroutine read_check_ionper

  !> @ brief Set the data value from the input file
  !!
  !!  Set the floating point value in the dblvec array using strings
  !!  returned from the parser.  Allow for time series names.
  !!
  !<
  subroutine set_value(this, ival)
    ! -- modules
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(TspSpcType), intent(inout) :: this !< TspSpcType object
    integer(I4B), intent(in) :: ival
    ! -- local
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: jj
    real(DP), pointer :: bndElem => null()
    !
    ! -- read remainder of variables on the line
    call this%parser%GetStringCaps(keyword)
    select case (keyword)
    case ('CONCENTRATION', 'TEMPERATURE')
      call this%parser%GetString(text)
      jj = 1 ! For CONCENTRATION
      bndElem => this%dblvec(ival)
      call read_value_or_time_series_adv(text, ival, jj, bndElem, this%packName, &
                                         'BND', this%tsManager, this%iprpak, &
                                         this%depvarname)

    end select
  end subroutine set_value

  !> @ brief check_flow_package
  !!
  !!  Check to make sure that flow package information is consistent
  !!  with this SPC information.
  !!
  !<
  subroutine check_flow_package(this, nbound_flowpack, budtxt)
    ! -- modules
    ! -- dummy
    class(TspSpcType), intent(inout) :: this !< TspSpcType object
    integer(I4B), intent(in) :: nbound_flowpack
    character(len=*), intent(in) :: budtxt
    ! -- local
    !
    ! -- Check and make sure MAXBOUND is not less than nbound_flowpack
    if (this%maxbound < nbound_flowpack) then
      write (errmsg, '(a, a, a, i0, a, i0, a)') &
            'The SPC Package corresponding to flow package ', &
            trim(this%packNameFlow), &
            ' has MAXBOUND set less than the number of boundaries &
            &active in this package.  Found MAXBOUND equal ', &
            this%maxbound, &
            ' and number of flow boundaries (NBOUND) equal ', &
            nbound_flowpack, &
            '. Increase MAXBOUND in the SPC input file for this package.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- If budtxt is RCHA or EVTA, then readasarrays must be used, otherwise
    !    readasarrays cannot be used
    select case (trim(adjustl(budtxt)))
    case ('RCHA')
      if (.not. this%readasarrays) then
        write (errmsg, '(a, a, a)') &
          'Array-based recharge must be used with array-based stress package &
          &concentrations.  GWF Package ', trim(this%packNameFlow), ' is being &
          &used with list-based SPC6 input.  Use array-based SPC6 input instead.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      end if
    case ('EVTA')
      if (.not. this%readasarrays) then
        write (errmsg, '(a, a, a)') &
          'Array-based evapotranspiration must be used with array-based stress &
          &package concentrations.  GWF Package ', trim(this%packNameFlow), &
          &' is being used with list-based SPC6 input.  Use array-based SPC6 &
          &input instead.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      end if
    case default
      if (this%readasarrays) then
        write (errmsg, '(a, a, a)') &
          'List-based packages must be used with list-based stress &
          &package concentrations.  GWF Package ', trim(this%packNameFlow), &
          &' is being used with array-based SPC6 input.  Use list-based SPC6 &
          &input instead.'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      end if
    end select
  end subroutine check_flow_package

end module TspSpcModule
