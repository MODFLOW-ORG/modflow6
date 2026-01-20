! This module implements all logic of the PreHeadsMF preprocessor.
module PreprocModule
  use ArrayHandlersModule,       only: ExpandArray
  use BlockParserModule,         only: BlockParserType
  use LineListModule,            only: LineListType
  use ConstantsModule,           only: DONE, DZERO, LENBIGLINE, &
                                       LINELENGTH, MAXCHARLEN,  &
                                       DHALF
  use ConstantsPHMFModule,       only: HUGEDBL, SINGLE, CONTINUOUS, &
                                       DEGREES_TO_RADIANS, HDRYDEFAULT, &
                                       LENOBSNAMENEW
  use DnmDis3dModule,            only: dis3d_cr, Dis3dType, CastAsDis3dType
  use DnmDisBaseModule,          only: DisBaseType
  use GLOBAL,                    only: NCOL, NROW, DELC, DELR
  use globalPHMF,                only: ioutPHMF, outfile
  use GlobalVariablesPHMFModule, only: prognamPHMF, verbose, vnam
  use InputOutputModule,         only: GetUnit, urword, GetUnit, openfile
  use ListModule,                only: ListType
  use ObsBlockModule,            only: ObsBlockType, ConstructObsBlockType, &
                                       AddObsBlockToList, GetObsBlockFromList
  use OpenSpecModule,            only: ACCESS, ACTION, FORM
  use SimPHMFModule,             only: count_errors, print_notes, store_error, &
                                       store_error_unit, ustop
  use UtilitiesModule,           only: get_extension, CalcContribFactors

  implicit none

  type :: PreprocType
    integer :: Idigits = 5
    integer :: Indis = 0
    integer :: Inunit = 0
    integer :: Iout = 0
    integer :: IoutMFobs = 0
    integer :: IoutPostObs = 0
    integer :: IPrecision = 2
    double precision :: Theta = DZERO
    double precision :: Hdry = hdrydefault
    double precision :: Xorigin = DZERO     ! X coordinate of upper left corner of grid
    double precision :: Yorigin = DZERO     ! Y coordinate of upper left corner of grid
    double precision :: XoriginLL = DZERO   ! X coordinate of lower left corner of grid
    double precision :: YoriginLL = DZERO   ! Y coordinate of lower left corner of grid
    character(len=MAXCHARLEN) :: BlockTypeFound = ''
    character(len=MAXCHARLEN) :: DisFilename = ''
    character(len=MAXCHARLEN) :: InputFilename = ''
    character(len=MAXCHARLEN) :: MFObsFilename = ''
    character(len=MAXCHARLEN) :: PostObsFilename = ''
    character(len=MAXCHARLEN), allocatable, dimension(:) :: PostObsOutputCsvFiles
    character(len=9) :: LengthUnit = 'UNDEFINED'
    logical :: Active = .false.
    logical :: Echo = .true.
    logical :: WriteMfOptions = .true.
    logical :: WriteBeginEnd = .true.
    class(DisBaseType), pointer :: dis => null()
    type(ListType), pointer     :: ObsBlockList => null()
    type(LineListType), private :: externalFiles
    type(BlockParserType) :: parser
  contains
    ! Public procedures
    procedure, public :: Run
    procedure, public :: find_lower_left
    ! Private procedures
    procedure, private :: add_obs_block
    procedure, private :: close_files
    procedure, private :: find_upper_left
    procedure, private :: initialize_preproc
    procedure, private :: open_files
    procedure, private :: process_options
    procedure, private :: read_any_block
    procedure, private :: read_blocks
    procedure, private :: read_options
    procedure, private :: write_mfobs_options
    procedure, private :: write_postobs_file
    procedure, private :: write_postobs_input_files
    procedure, private :: write_postobs_options
  end type PreprocType

  ! data
  character(len=10), dimension(2) :: btype
  data btype /'SINGLE    ','CONTINUOUS'/

contains

  subroutine Run(this, fname, WriteBeginEndArg)
    ! dummy
    class(PreprocType), intent(inout) :: this
    character(len=*), intent(in) :: fname
    logical, intent(in), optional :: WriteBeginEndArg
    ! local
    integer :: inpp, iu
    logical :: lop, WriteBeginEnd
    ! formats
    5 format(a)
    10 format(/,15x,a,/)
    20 format(/,'End of PreHeadsMF',/)
    !
    if (present(WriteBeginEndArg)) then
      WriteBeginEnd = WriteBeginEndArg
    else
      WriteBeginEnd = .true.
    endif
    !
    ! Open the main input file
    write(*,10)trim(prognamPHMF) // ' ' // trim(vnam)
    write(*,*)'Using input file: ',trim(fname)
    !
    ! Open an output file for messages, etc.
    inquire(file=outfile, opened=lop, number=iu)
    if (lop) then
      ioutPHMF = iu
    else
      ioutPHMF = GetUnit()
      open(unit=ioutPHMF,file=outfile,status='REPLACE')
    endif
    write(ioutPHMF,5)'Output listing from program ' // trim(prognamPHMF)
    !
    inpp = GetUnit()
    call openfile(inpp,0,fname,'PREPROC')
    !
    call this%initialize_preproc(inpp, ioutPHMF)
    !
    ! Read Options block
    call this%read_options()
    !
    ! Process options
    call this%process_options()
    !
    ! Open the MF obs input file and the PostObs input file
    call this%open_files()
    !
    call print_notes()
    !
    ! Read and process SINGLE and CONTINUOUS blocks,
    ! and write OBS8 input for MODFLOW 6.
    call this%read_blocks(WriteBeginEnd)
    !
    ! Write input for PostObsMF
    call this%write_postobs_file()
    !
    ! Close all files
    call this%close_files()
    !
    write(*,20)
    !
    return
  end subroutine Run

  subroutine initialize_preproc(this, inpp, iout)
    ! dummy
    class(PreprocType) :: this
    integer, intent(in) :: inpp, iout
    !
    allocate(this%ObsBlockList)
    allocate(this%PostObsOutputCsvFiles(0))
    call this%externalFiles%InitializeLineList()
    this%Inunit = inpp
    this%Iout = iout
    call this%parser%Initialize(inpp, iout)
    !
    return
  end subroutine initialize_preproc

  subroutine read_options(this)
    ! dummy
    class(PreprocType) :: this
    ! local variables
    integer :: iin
    integer :: ierr
    integer :: localprecision
    integer :: localdigits
    character(len=20) :: blockTypeWanted
    character(len=40) :: keyword
    character(len=LINELENGTH) :: ermsg
    character(len=LINELENGTH) :: errormessage, fname
    type(ListType), pointer :: lineList => null()
    double precision :: angle
    logical :: continueread, found, endOfBlock
    ! -- formats
 1  format()
10  format('No options block found in ',a,' input.')
20  format('Error reading begin/end block: ',a)
30  format('Binary output precision set to: ',a)
40  format('Text output number of digits of precision set to: ',i2)
55  format(a,' set to: ',g14.7)
60  format(/,'Processing options:',/)
70  format('DIS input file: ',a)
75  format(a,': ',a)
80  format(a)
    !
    localprecision = 2
    localdigits = 0
    lineList => null()
    !
    ! -- Find and store file name
    iin = this%inUnit
    inquire(unit=iin, name=fname)
    this%inputFilename = fname
    !
    ! -- Read Options block
    blockTypeWanted = 'OPTIONS'
    continueread = .false.
    ierr = 0
    !
    ! -- get BEGIN line of OPTIONS block
    call this%parser%GetBlock('OPTIONS', found, ierr, supportOpenClose=.true.)
    if (ierr /= 0) then
      ! end of file
      ermsg = 'End-of-file encountered while searching for' // &
              ' OPTIONS in ' // &
              'input file "' // trim(this%inputFilename) // '"'
      call store_error(ermsg)
      call this%parser%StoreErrorUnit()
      call ustop()
    elseif (.not. found) then
      this%blockTypeFound = ''
      write(ermsg,10)trim(prognamPHMF)
      call store_error(ermsg)
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    !
    ! -- parse OPTIONS entries
    if (found) then
      write(this%iout,60)
      readblockoptions: do
!        lloc = 1
        ! -- read a line from input
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
!        if (eof) then
!          ! End-of-file encountered. Encountered here, this is an error.
!          ermsg = 'End-of-file encountered in OPTIONS block.'
!          call store_error(ermsg)
!          call store_error_unit(iin)
!          call ustop()
!        endif
        !
        ! -- get the keyword
!        call urword(line, lloc, istart, istop, 1, idum, rdum, this%iout, &
!                    iin)
!        keyword = line(istart:istop)
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
! Following section commented out unless/until mf6 obs supports specifying
! single or double precision for binary output.
!        case ('PRECISION')
!          ! -- Specifies SINGLE or DOUBLE precision for writing simulated values
!          !    to a binary file. Default is DOUBLE.
!          ! -- get the word following the keyword (the key value)
!          call urword(line, lloc, istart, istop, 1, idum, rdum, this%iout, &
!                      iin)
!          keyvalue = line(istart:istop)
!          if (localprecision==0) then
!            if (keyvalue=='SINGLE') then
!              localprecision = 1
!              write(this%iout,30)'SINGLE'
!            elseif (keyvalue=='DOUBLE') then
!              localprecision = 2
!              write(this%iout,30)'DOUBLE'
!            else
!              errormessage = 'Error in input: "'//trim(keyvalue)// &
!                  '" is not a valid option for PRECISION'
!!              write(*,20)trim(errormessage)
!!              if (this%iout>0) write(this%iout,20)trim(errormessage)
!              call store_error(errormessage)
!              goto 900
!            endif
!          else
!            errormessage = 'Error in input: PRECISION has already been defined'
!!            write(*,20)trim(errormessage)
!!            if (this%iout>0) write(this%iout,20)trim(errormessage)
!            call store_error(errormessage)
!            goto 900
!          endif
        case ('DIGITS')
          ! -- Specifies number of significant digits used writing simulated
          !    values to a text file. Default is 5 digits.
          if (localdigits==0) then
!            call urword(line,lloc,istart,istop,2,localdigits,rdum,this%iout,iin)
            localdigits = this%parser%GetInteger()
            if (localdigits < 1) then
              errormessage = 'Error in ' // trim(prognamPHMF) // &
                             ' input: Invalid value for DIGITS option'
              call store_error(errormessage)
              goto 900
            endif
            if (localdigits < 2) localdigits = 2
            if (localdigits > 17) localdigits = 17
            write(this%iout,40)localdigits
          else
            errormessage = 'Error in ' // trim(prognamPHMF) // &
                           ' input: DIGITS has already been defined'
            call store_error(errormessage)
            goto 900
          endif
!        case ('EXTERNAL')
!          ! Get the file type keyword and verify that it is FILEIN
!!          call urword(line, lloc, istart, istop, 1, idum, rdum, this%iout, &
!!                      iin)
!          call this%parser%GetStringCaps(keyword)
!          if (keyword /= 'FILEIN') then
!            ermsg = 'Expected "FILEIN" but found "' // line(istart:istop) &
!                    // '"'
!            call store_error(ermsg)
!            call this%parser%StoreErrorUnit()
!            call ustop()
!          endif
!          ! -- get the word following the keyword (the key value)
!!          call urword(line, lloc, istart, istop, 0, idum, rdum, this%iout, &
!!                      iin)
!          call this%parser%GetString(extfilename)
!!          extfilename = line(istart:istop)
!          inquire(file=extfilename,exist=lex)
!          if (.not. lex) then
!            errormessage = 'Error in OBS input. External file does not' // &
!                           ' exist: ' // trim(extfilename)
!            call store_error(errormessage)
!            goto 900
!          endif
!          call this%externalFiles%AddLine(extfilename)
        case ('VERBOSE')
          verbose = .true.
          write(this%iout,'(a)')'The VERBOSE option has been specified.'
        case ('OMITOPTIONS')
          this%WriteMfOptions = .false.
          write(this%iout,'(a)')'The OMITOPTIONS option has been specified.'
        case ('DIS')
          if (this%DisFilename /= '') then
            ermsg = 'DIS has already been specified.'
            call store_error(ermsg)
            call store_error_unit(iin)
            call ustop()
          else
            ! Get the file type keyword and verify that it is FILEIN
!            call urword(line, lloc, istart, istop, 1, idum, rdum, this%iout, &
!                        iin)
            call this%parser%GetStringCaps(keyword)
            if (keyword /= 'FILEIN') then
              ermsg = 'Expected "FILEIN" but found "' // trim(keyword) &
                      // '"'
              call store_error(ermsg)
              call this%parser%StoreErrorUnit()
              call ustop()
            endif
            ! Get DIS file name
!            call urword(line, lloc, istart, istop, 0, idum, rdum, this%iout, &
!                        iin)
!            this%DisFilename = line(istart:istop)
            call this%parser%GetString(this%DisFilename)
            write(this%iout,70)trim(this%DisFilename)
          endif
        case ('XORIGIN')
          ! Get xorigin
!          call URWORD(line, lloc, istart, istop, 3, idum, this%Xorigin, &
!                      this%Iout, iin)
          this%XoriginLL = this%parser%GetDouble()
          write(this%iout,55)'Xorigin',this%XoriginLL
        case ('YORIGIN')
          ! Get yorigin
!          call URWORD(line, lloc, istart, istop, 3, idum, this%Yorigin, &
!                      this%Iout, iin)
          this%YoriginLL = this%parser%GetDouble()
          write(this%iout,55)'Yorigin',this%YoriginLL
        case ('ANGLE')
          ! Get angle
!          call URWORD(line, lloc, istart, istop, 3, idum, angle, &
!                      this%Iout, iin)
          angle = this%parser%GetDouble()
          write(this%iout,55)'Angle',angle          ! Rotation angle of grid, in degrees
          this%Theta = angle * DEGREES_TO_RADIANS   ! Rotation angle of grid, in radians
        case ('MFOBSFILE')
          ! Get the file type keyword and verify that it is FILEOUT
!          call urword(line, lloc, istart, istop, 1, idum, rdum, this%iout, &
!                      iin)
          call this%parser%GetStringCaps(keyword)
          if (keyword /= 'FILEOUT') then
            ermsg = 'Expected "FILEOUT" but found "' // trim(keyword) &
                    // '"'
            call store_error(ermsg)
            call this%parser%StoreErrorUnit()
            call ustop()
          endif
          ! Get name of MF6 OBS8 file to be created
!          call URWORD(line, lloc, istart, istop, 0, idum, rdum, &
!                      this%Iout, iin)
!          this%MFObsFilename = line(istart:istop)
          call this%parser%GetString(this%MFObsFilename)
          write(this%iout,75)'MODFLOW 6 head observations input file',trim(this%MFObsFilename)
        case ('POSTOBSFILE')
          ! Get the file type keyword and verify that it is FILEOUT
!          call urword(line, lloc, istart, istop, 1, idum, rdum, this%iout, &
!                      iin)
          call this%parser%GetStringCaps(keyword)
          if (keyword /= 'FILEOUT') then
            ermsg = 'Expected "FILEOUT" but found "' // trim(keyword) &
                    // '"'
            call store_error(ermsg)
            call this%parser%StoreErrorUnit()
            call ustop()
          endif
          ! Get name of PostObsMF input file to be created
!          call URWORD(line, lloc, istart, istop, 0, idum, rdum, &
!                      this%Iout, iin)
!          this%PostObsFilename = line(istart:istop)
          call this%parser%GetString(this%PostObsFilename)
          write(this%Iout,75)'PostObsMF input file',trim(this%PostObsFilename)
        case ('LENGTH_UNITS')
          ! Get geospatial length unit
!          call URWORD(line, lloc, istart, istop, 1, idum, rdum, &
!                      this%Iout, iin)
!          this%LengthUnit = line(istart:istop)
          call this%parser%GetStringCaps(this%LengthUnit)
          write(this%iout,75)'Geospatial length unit: ',trim(this%LengthUnit)
!        case ('END')
!          call uterminate_block(iin, this%iout, keyword, blockTypeWanted, lloc, &
!                                line, ierr, iuext)
!          write(this%iout,1)
!          write(this%iout,80)'End of options'
!          exit readblockoptions
        case default
          errormessage = 'Error in ' // trim(prognamPHMF) // &
                         ' input: Unrecognized option: ' // &
                         trim(keyword)
          call store_error(errormessage)
          goto 900
        end select
      enddo readblockoptions
    endif
    !
    write(this%iout,'(1x)')
    !
    ! -- Assign type variables
    if (localprecision>0) this%iprecision = localprecision
    if (localdigits>0) this%idigits = localdigits
    !
    ! -- normal return
    return
    !
    ! -- Errors transfer here
900 continue
    if (count_errors()>0) then
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    return
  end subroutine read_options

  subroutine process_options(this)
    ! Check for required options and define discretization.
    ! dummy
    class(PreprocType) :: this
    ! local
    character(len=MAXCHARLEN) :: ermsg
    character(len=11) :: gslenunit, gridlenunit
    character(len=11), dimension(0:3) :: lenunit
    data lenunit(0:3)/'UNDEFINED  ','FEET       ','METERS     ','CENTIMETERS'/
    ! format
    10 format('Geospatial length unit given as ',a, &
           ' but grid length unit is undefined in the DIS input file. The', &
           ' DIS file needs to be edited to set LENGTH_UNITS.')
    !
    if (this%MFObsFilename == '') then
      ermsg = 'MFOBSFILE not specified in OPTIONS block.'
      call store_error(ermsg)
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    !
    if (this%PostObsFilename == '') then
      ermsg = 'POSTOBSFILE not specified in OPTIONS block.'
      call store_error(ermsg)
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    !
    ! -- Read the DIS file
    if (this%DisFilename /= '') then
      ! Open the DIS file
      call openfile(this%Indis, this%Iout, this%DisFilename, 'DIS')
      ! Create discretization object
      call dis3d_cr(this%dis, 'unknown_model', this%indis, this%iout)
      ! Read discretization information
      call this%dis%dis_df()
    else
      ermsg = 'DIS not specified in OPTIONS block.'
      call store_error(ermsg)
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    !
    ! Find x,y coordinates of upper left corner of grid
    call this%find_upper_left()
    !
    this%dis%Xorigin = this%Xorigin
    this%dis%Yorigin = this%Yorigin
    this%dis%Theta = this%Theta
    !
    ! Check length units. Assign conversion factor
    ! to convert geospatial length unit to grid
    ! length unit if conversion is needed.
    gridlenunit = lenunit(this%dis%lenuni)
    gslenunit = this%LengthUnit
    ermsg = ''
    select case (gslenunit)
    case ('UNDEFINED')
      ! OK -- leave conversion factor = 1.0
    case ('FEET')
      select case (gridlenunit)
      case ('UNDEFINED')
        ! Not OK -- flag as error
        write(ermsg,10)trim(gslenunit)
      case ('FEET')
        ! OK -- leave conversion factor = 1.0
      case ('METERS')
        ! Need to convert FEET to METERS
        this%dis%ConvertFactor  = 0.3048d0
      case ('CENTIMETERS')
        ! Need to convert FEET to CENTIMETERS
        this%dis%ConvertFactor  = 30.48d0
      end select
    case ('METERS')
      select case (gridlenunit)
      case ('UNDEFINED')
        ! Not OK -- flag as error
        write(ermsg,10)trim(gslenunit)
      case ('FEET')
        ! Need to convert METERS to FEET
        this%dis%ConvertFactor  = DONE / 0.3048d0
      case ('METERS')
        ! OK -- leave conversion factor = 1.0
      case ('CENTIMETERS')
        ! Need to convert METERS to CENTIMETERS
        this%dis%ConvertFactor  = 100.0d0
      end select
    end select
    !
    if (ermsg /= '') then
      call store_error(ermsg)
      call this%parser%StoreErrorUnit()
      call ustop()
    endif
    !
    return
  end subroutine process_options

  subroutine find_upper_left(this)
    ! Find x,y coordinates of upper left corner of grid.
    ! dummy
    class(PreprocType), intent(inout) :: this
    ! local
    double precision :: a, b, sumdelc
    integer :: i
    character(len=MAXCHARLEN) :: ermsg
    type(Dis3dType), pointer :: dis3d
    class(*), pointer :: obj => null()
    !
    obj => this%dis
    dis3d => CastAsDis3dType(obj)
    if (associated(dis3d)) then
      ! Find height of grid (sum of DELC)
      sumdelc = DZERO
      do i=1,dis3d%nrow
        sumdelc = sumdelc + dis3d%delc(i)
      enddo
      ! Find upper left
      a = sumdelc * sin(this%Theta)
      b = sumdelc * cos((this%Theta))
      this%Xorigin = this%XoriginLL - a
      this%Yorigin = this%YoriginLL + b
    else
      ermsg = 'Error: Discretization is not 3D.'
      call store_error(ermsg)
      call ustop()
    endif
    !
    return
  end subroutine find_upper_left

  subroutine find_lower_left(this)
    ! Find x,y coordinates of lower left corner of grid, based on delr.
    ! dummy
    class(PreprocType), intent(inout) :: this
    ! local
    integer :: j
    double precision :: yo
    !
    this%XoriginLL = DZERO
    this%YoriginLL = DZERO
    !
    yo = DZERO
    do j=1,NROW
      yo = yo - DELR(j)
    enddo
    this%YoriginLL = yo
    !
    return
  end subroutine find_lower_left

  subroutine open_files(this)
    ! Open the MF obs input file and the PostObs input file
    ! dummy
    class(PreprocType) :: this
    ! local
    character(len=MAXCHARLEN) :: accarg, fmtarg
    !
    fmtarg = 'FORMATTED'
    accarg = 'SEQUENTIAL'
    !
    ! Open the MF obs input file and write an Options block
    call openfile(this%IoutMFobs, this%iout, this%MFObsFilename, &
                  'MF6 OBS INPUT', fmtarg, accarg, 'REPLACE')
    call this%write_mfobs_options()
    !
    ! Open the PostObs input file
    call openfile(this%IoutPostObs, this%iout, this%PostObsFilename, &
                  'PostObs INPUT', fmtarg, accarg, 'REPLACE')
    !
    return
  end subroutine open_files

  subroutine close_files(this)
    class(PreprocType) :: this
    !
    close(this%IoutMFobs)
    close(this%IoutPostObs)
!    close(this%Iout)
    !
    return
  end subroutine close_files

  subroutine write_mfobs_options(this)
    ! dummy
    class(PreprocType) :: this
    ! local
    integer :: ioutmf
    ! formats
    1 format()
    5 format('# File prepared by PreHeadsMF. Copy SINGLE and CONTINUOUS',/, &
             '# blocks from this file into the OBS8 input file listed ',/, &
             '# in your MODFLOW 6 model name file, or use this file as',/, &
             '# the OBS6 input file.',/)
    10 format(a)
    20 format(2x,a)
    30 format(2x,a,2x,i0)
    40 format(2x,a,2x,a)
    !
    ioutmf = this%IoutMFobs
    !
    if (.not. this%WriteMfOptions) return
    write(ioutmf,5)
    write(ioutmf,10)'BEGIN OPTIONS'
!    write(ioutmf,40)'PRECISION', prec(this%IPrecision)
    write(ioutmf,30)'DIGITS', this%Idigits
!    write(ioutmf,20)'STRUCTURED'
    write(ioutmf,20)'PRINT_INPUT'
    write(ioutmf,10)'END OPTIONS'
    !
    return
  end subroutine write_mfobs_options

  subroutine read_blocks(this, WriteBeginEnd)
    ! dummy
    class(PreprocType), intent(inout) :: this
    logical, intent(in) :: WriteBeginEnd
    ! local
    integer :: i, iu, k, numextfiles
    logical :: eof
    character(len=MAXCHARLEN) :: fname, msg
    type(Dis3dType), pointer :: dis3d => null()
    ! format
    20 format(/,a)
    !
    select type (d => this%dis)
    type is (Dis3dType)
      dis3d => d
    end select
    !
    msg = 'Processing SINGLE and/or CONTINUOUS blocks:'
    write(*,20)trim(msg)
    write(this%Iout,20)trim(msg)
    !
    k = 0
    eof = .false.
    ! Read and process SINGLE and CONTINUOUS blocks in all external files.
    numextfiles = this%externalFiles%CountLines()
    do i=1,numextfiles
      call this%externalFiles%GetLine(i, fname)
      iu = GetUnit()
      call openfile(iu, 0, fname, 'EXTERNAL', filstat_opt='OLD')
      ! Read all blocks in this external file
      loopext: do
        ! -- Read block from external file
        call this%read_any_block(iu, k, eof, dis3d, .true.)
        if (eof) exit loopext
      enddo loopext
      close(iu)
    enddo
    !
    ! Read and process all SINGLE and CONTINUOUS blocks in main input file.
    iu = this%Inunit
    loop: do
      call this%read_any_block(iu, k, eof, dis3d, WriteBeginEnd)
      if (eof) exit
    enddo loop
    !
    return
  end subroutine read_blocks

  subroutine read_any_block(this, iu ,k, eof, dis3d, WriteBeginEnd)
    ! dummy
    class(PreprocType), intent(inout) :: this
    integer, intent(in)    :: iu
    integer, intent(inout) :: k
    logical, intent(inout) :: eof
    type(Dis3dType), pointer, intent(inout) :: dis3d
    logical, intent(in) :: WriteBeginEnd
    ! local
    integer :: ierr, lloc
    logical :: binary, insertLine, isfound
    character(len=MAXCHARLEN) :: accarg, ctagfound, ermsg, fmtarg, &
                                 bname, fname, keyword
    type(ObsBlockType), pointer :: obsblock => null()
    integer :: iuext
    ! format
    10 format('Expecting block name SINGLE or CONTINUOUS, but found: ',a)
    20 format(/,a)
    30 format(/,'Observation type: ',a, &
              /,'Output base name: ',a, &
              /,'Binary or text: ',a)
    !
    ! -- Read any block as long as it's SINGLE or CONTINUOUS.
    lloc = 1
    call this%parser%GetBlock('*', isfound, ierr, .true., &
                              .false., ctagfound)
    if (.not. isfound) then
      eof = .true.
      return
    endif !exit loop
    if (ctagfound /= 'SINGLE' .and. ctagfound /= 'CONTINUOUS') then
      write(ermsg,10)trim(ctagfound)
      call store_error(ermsg)
      call store_error_unit(iu)
      call ustop()
    endif
    !
!    ! Check for FILEOUT keyword
!    call URWORD(line, lloc, istart, istop, 1, ndum, rdum, this%iout, iu)
!    if (line(istart:istop) /= 'FILEOUT') then
!      ermsg = 'Expected FILEOUT but found "' // line(istart:istop) //'"'
!      call store_error(ermsg)
!      call store_error_unit(iu)
!      call ustop()
!    endif
    !
    ! Get file basename and BINARY option.
!    call URWORD(line, lloc, istart, istop, 0, ndum, rdum, this%iout, iu)
!    bname = line(istart:istop)
    call this%parser%GetString(bname)
    binary = .false.
!    call URWORD(line, lloc, istart, istop, 1, ndum, rdum, this%iout, iu)
    call this%parser%GetStringCaps(keyword)
    if (keyword== 'BINARY') binary = .true.
    if (binary) then
      fname = trim(bname) // '.bsv'
      fmtarg = FORM
      accarg = ACCESS
      write(*,30)trim(ctagfound),trim(bname),'BINARY'
      write(this%Iout,30)trim(ctagfound),trim(bname),'BINARY'
    else
      fname = trim(bname) // '.csv'
      fmtarg = 'FORMATTED'
      accarg = 'SEQUENTIAL'
      write(*,30)trim(ctagfound),trim(bname),'TEXT'
      write(this%Iout,30)trim(ctagfound),trim(bname),'TEXT'
    endif
    !
    ! Allocate a new ObsBlock
    call ConstructObsBlockType(obsblock, ctagfound, bname, binary)
    obsblock%dis3d => dis3d
    obsblock%iout = this%Iout
    !obsblock%inunit = this%Inunit
    obsblock%inunit = this%parser%iuactive
    iuext = this%parser%iuactive
    obsblock%ioutMFobs = this%IoutMFobs
    obsblock%IoutPostObs = this%IoutPostObs
    !
    ! Process block and add it to list
    k = k + 1
    if (k == 1) then
      insertLine = this%WriteMfOptions
    else
      insertLine = .true.
    endif
    call obsblock%process_block(insertLine, WriteBeginEnd, this%parser)
    call this%add_obs_block(obsblock)
    !
    return
  end subroutine read_any_block

  subroutine add_obs_block(this, obsblock)
    ! dummy
    class(PreprocType), intent(inout) :: this
    type(ObsBlockType), pointer, intent(inout) :: obsblock
    !
    call AddObsBlockToList(this%ObsBlockList, obsblock)
    !
    return
  end subroutine add_obs_block

  subroutine write_postobs_file(this)
    ! Write all blocks of a PostObsMF input file
    ! dummy
    class(PreprocType) :: this
    ! local
    integer :: i, iout, n, nobsblocks
    type(ObsBlockType), pointer :: obsblock => null()
    character(len=MAXCHARLEN) :: outputfilename
    character(len=3) :: ext
    ! formats
    1 format()
    10 format(a,1x,a,1x,a)
    15 format(a,1x,a)
    20 format(2x,a)
    !
    iout = this%IoutPostObs
    !
    ! Write Options block
    call this%write_postobs_options()
    !
    ! Write Input_Files block
    call this%write_postobs_input_files()
    !
    ! Write SINGLE and CONTINUOUS blocks
    nobsblocks = this%ObsBlockList%Count()
    do i=1,nobsblocks
      obsblock => GetObsBlockFromList(this%ObsBlockList, i)
      call obsblock%write_postobs_input(outputfilename)
      call get_extension(outputfilename, ext)
      if (ext == 'CSV') then
        call ExpandArray(this%PostObsOutputCsvFiles)
        n = size(this%PostObsOutputCsvFiles)
        this%PostObsOutputCsvFiles(n) = outputfilename
      endif
    enddo
    !
    return
  end subroutine write_postobs_file

  subroutine write_postobs_options(this)
    ! dummy
    class(PreprocType), intent(inout) :: this
    ! local
    integer :: iout
    ! formats
    5 format('# PostObsMF input file prepared by PreHeadsMF.',/)
    10 format(a,1x,a)
    20 format(2x,a)
    30 format(2x,a,2x,a)
    40 format(2x,a,2x,i0)
    !
    iout = this%IoutPostObs
    write(iout,5)
    !
    ! Write BEGIN line
    write(iout,10)'BEGIN', 'OPTIONS'
    !
!    write(iout,30)'PRECISION', prec(this%IPrecision)
    write(iout,40)'DIGITS', this%Idigits
    !
    ! Write END line
    write(iout,10)'END', 'OPTIONS'
    !
    return
  end subroutine write_postobs_options

  subroutine write_postobs_input_files(this)
    ! Write the INPUT_FILES block to be read by PostObsMF
    ! dummy
    class(PreprocType), intent(inout) :: this
    ! local
    integer :: i, iout, nblocks
    character(len=MAXCHARLEN) :: fname
    type(ObsBlockType), pointer :: obsblock => null()
    ! formats
    1 format()
    10 format(a,1x,a)
    20 format(2x,a,2x,'FILEIN',2x,a)
    30 format(2x,a,2x,'FILEIN',2x,a,2x,'BINARY')
    !
    iout = this%IoutPostObs
    write(iout,1)
    !
    ! Write BEGIN line
    write(iout,10)'BEGIN', 'INPUT_FILES'
    !
    nblocks = this%ObsBlockList%Count()
    do i=1,nblocks
      obsblock => GetObsBlockFromList(this%ObsBlockList, i)
      if (obsblock%Binary) then
        fname = trim(obsblock%OutputBaseName) // '.bsv'
        write(iout,30)trim(obsblock%SorC), trim(fname)
      else
        fname = trim(obsblock%OutputBaseName) // '.csv'
        write(iout,20)trim(obsblock%SorC), trim(fname)
!        call ExpandArray(this%PostObsCsvFiles)
!        n = size(this%PostObsCsvFiles)
!        this%PostObsCsvFiles(n) = fname
      endif
    enddo
    !
    ! Write END line
    write(iout,10)'END', 'INPUT_FILES'
    !
    return
  end subroutine write_postobs_input_files

end module PreprocModule








