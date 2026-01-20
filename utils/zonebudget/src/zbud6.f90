program zonbudmf6
  use KindModule
  use ConstantsModule, only: LINELENGTH, LENHUGELINE
  use VersionModule, only: VERSION
  use SimVariablesModule, only: iout, errmsg
  use SimModule, only: store_error
  use MessageModule, only: write_message, write_message_centered
  use InputOutputModule, only: openfile

  implicit none

  character(len=10), parameter :: mfvnam = ' Version 6'
  character(len=LINELENGTH) :: line
  character(len=LENHUGELINE) :: fnam, flst, fcsv
  integer(I4B) :: iunit_lst = 20
  integer(I4B) :: iunit_csv = 21
  integer(I4B) :: iunit_nam = 22
  integer(I4B) :: iunit_bud = 23
  integer(I4B) :: iunit_zon = 24
  integer(I4B) :: iunit_grb = 25
  logical :: exists

  ! -- Write title to screen
  call write_message_centered('ZONEBUDGET'//mfvnam, 80)
  call write_message_centered('U.S. GEOLOGICAL SURVEY', 80)
  call write_message_centered('VERSION '//VERSION, 80)
  !
  ! -- Find name of zone budget name file and lst file
  fnam = 'zbud.nam'
  call parse_command_line(fnam, flst, fcsv)
  inquire (file=fnam, exist=exists)
  if (.not. exists) then
    write (errmsg, '(a)') 'Name file not found. Looking for: '//trim(fnam)
    call store_error(errmsg, terminate=.TRUE.)
  end if
  !
  ! -- Open list file and write title
  iout = iunit_lst
  call openfile(iunit_lst, 0, flst, 'LIST', filstat_opt='REPLACE')
  call write_message_centered('ZONEBUDGET'//mfvnam, 80, iunit=iout)
  call write_message_centered('U.S. GEOLOGICAL SURVEY', 80, iunit=iout)
  call write_message_centered('VERSION '//VERSION, 80, iunit=iout)
  !
  ! -- Open name file, read name file, and open csv file
  call openfile(iunit_nam, iout, fnam, 'NAM')
  call read_namefile(iunit_nam, iunit_bud, iunit_zon, iunit_grb)
  call openfile(iunit_csv, iout, fcsv, 'CSV', filstat_opt='REPLACE')
  !
  ! -- Process the budget file and write output
  call process_budget(iunit_csv, iunit_bud, iunit_zon, iunit_grb)
  !
  ! -- close output files
  write (iunit_lst, '(/, a)') 'Normal Termination'
  close (iunit_lst)
  close (iunit_csv)
  write (line, '(a)') 'Normal Termination'
  call write_message(line, skipbefore=1)
  !
  ! -- end of program
end program zonbudmf6

subroutine read_namefile(iunit_nam, iunit_bud, iunit_zon, iunit_grb)
! ******************************************************************************
! read_namefile
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- modules
  use KindModule
  use SimVariablesModule, only: iout, errmsg
  use SimModule, only: store_error
  use ConstantsModule, only: LENHUGELINE, LINELENGTH
  use InputOutputModule, only: openfile
  use OpenSpecModule, only: form, access
  use BlockParserModule, only: BlockParserType
  implicit none
  ! -- dummy
  integer, intent(in) :: iunit_nam
  integer, intent(in) :: iunit_bud
  integer, intent(in) :: iunit_zon
  integer, intent(inout) :: iunit_grb
  ! -- local
  type(BlockParserType) :: parser
  integer(I4B) :: ierr, iu
  logical :: isfound, endOfBlock
  character(len=LINELENGTH) :: keyword
  character(len=LENHUGELINE) :: filename
  character(len=20) :: fm, acc
! ------------------------------------------------------------------------------
  !
  call parser%Initialize(iunit_nam, iout)
  call parser%GetBlock('ZONEBUDGET', isfound, ierr)
  if (isfound) then
    do
      call parser%GetNextLine(endOfBlock)
      if (endOfBlock) exit
      call parser%GetStringCaps(keyword)
      fm = 'FORMATTED'
      acc = 'SEQUENTIAL'
      select case (keyword)
      case ('BUD')
        iu = iunit_bud
        fm = form
        acc = access
        call parser%GetString(filename)
      case ('ZON')
        iu = iunit_zon
        call parser%GetString(filename)
      case ('GRB')
        iu = iunit_grb
        fm = form
        acc = access
        call parser%GetString(filename)
      case default
        write (errmsg, '(4x,a,a)') 'ERROR. UNKNOWN ZONEBUDGET ENTRY: ', &
          trim(keyword)
        call store_error(errmsg)
        call parser%StoreErrorUnit()
      end select
      call openfile(iu, iout, trim(filename), trim(keyword), fm, acc)
    end do
  else
    write (errmsg, '(1x,a)') 'ERROR.  REQUIRED ZONEBUDGET BLOCK NOT FOUND.'
    call store_error(errmsg)
    call parser%StoreErrorUnit()
  end if
  !
  ! -- close name file
  close (iunit_nam)
  !
  ! -- return
  return
end subroutine read_namefile

subroutine process_budget(iunit_csv, iunit_bud, iunit_zon, iunit_grb)
! ******************************************************************************
! process_budget
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- modules
  use KindModule
  use ConstantsModule, only: LINELENGTH
  use SimVariablesModule, only: iout, errmsg
  use SimModule, only: store_error
  use BudgetDataModule, only: budgetdata_init, budgetdata_read, &
                              budgetdata_finalize, &
                              ia, ja, budtxt, nbudterms, &
                              nodesrc, nodedst, flowdata, flowja, kper, kstp, &
                              delt, totim, dstpackagename, hasimeth1flowja, &
                              srcmodelname, dstmodelname
  use ZoneModule, only: zone_init, clear_accumulators, &
                        flowja_accumulate, flowiaja_accumulate, &
                        flow_accumulate, &
                        flowch_setich, flowch_accumulate, &
                        zone_finalize, nmznfl, vbvl, vbznfl, maxzone
  use ZoneOutputModule, only: zoneoutput_init, zoneoutput_write, &
                              zoneoutput_finalize
  use GridFileReaderModule, only: GridFileReaderType
  use MessageModule, only: write_message, write_message_centered
  implicit none
  ! -- dummy
  integer, intent(in) :: iunit_csv
  integer, intent(in) :: iunit_bud
  integer, intent(in) :: iunit_zon
  integer, intent(in) :: iunit_grb
  ! -- local
  character(len=1) :: cdot
  character(len=16), dimension(:), allocatable :: budtxtarray
  character(len=16), dimension(:), allocatable :: packagenamearray
  integer, dimension(:), allocatable :: internalflow
  integer, allocatable, dimension(:) :: mshape
  integer(I4B) :: ibudterm
  integer(I4B) :: itime = 1
  integer(I4B) :: ncrgrb
  integer(I4B) :: ncrbud = 0
  integer(I4B) :: ncr
  logical :: opengrb
  logical :: success
  logical :: hasiaja = .false.
  logical :: foundchd = .false.
  type(GridFileReaderType) :: gfr
! ------------------------------------------------------------------------------
  !
  ! -- initialize local variables
  cdot = '.'
  !
  ! -- Initialize budget data
  call budgetdata_init(iunit_bud, iout, ncrbud)
  !
  ! -- Check to see if GRB is required, and read it if necessary
  ncrgrb = 0
  if (hasimeth1flowja) then
    inquire (unit=iunit_grb, opened=opengrb)
    if (opengrb) then
      hasiaja = .true.
      call gfr%initialize(iunit_grb)
      mshape = gfr%read_grid_shape()
      ia = gfr%read_int_1d("IA")
      ja = gfr%read_int_1d("JA")
      call gfr%finalize()
      ncrgrb = size(ia) - 1
    else
      errmsg = 'BUDGET FILE HAS "FLOW-JA-FACE" RECORD BUT NO GRB FILE SPECIFIED.'
      call store_error(errmsg)
      errmsg = 'ADD GRB ENTRY TO ZONE BUDGET NAME FILE.'
      call store_error(errmsg, terminate=.TRUE.)
    end if
  else
    inquire (unit=iunit_grb, opened=opengrb)
    if (opengrb) then
      errmsg = 'BINARY GRID FILE IS PRESENT, BUT BUDGET FILE DOES NOT HAVE &
               &"FLOW-JA-FACE" RECORD IN THE IMETH=1 FORMAT. CHECK TO MAKE SURE &
               &FLOWS ARE SAVED TO THE BUDGET FILE'
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- At this point, must be a budget file from an advanced package without
    !    the IMETH=1 flow-ja-face record.
    allocate (mshape(1))
    mshape(1) = ncrgrb
  end if
  !
  ! -- Read the zone file to get number of cells/reaches
  ncr = ncrgrb
  call zone_init(iunit_zon, nbudterms, ncr, mshape)
  !
  ! -- Initialize zone and zoneoutput modules
  !call zone_init(iunit_zon, nbudterms, ncr)
  call zoneoutput_init(iout, iunit_csv, maxzone, nbudterms)
  allocate (budtxtarray(nbudterms))
  allocate (packagenamearray(nbudterms))
  allocate (internalflow(nbudterms))
  !
  ! -- time loop
  timeloop: do
    !
    ! -- Clear budget accumulators and loop through budget terms
    call clear_accumulators()
    write (iout, '(/, a)') 'Reading records from budget file'
    do ibudterm = 1, nbudterms
      !
      ! -- read data
      call budgetdata_read(success, iout)
      if (.not. success) then
        write (iout, '(a)') 'Done reading records.  Exiting time loop.'
        exit timeloop
      end if
      !
      ! -- write message and check
      call write_message(text=cdot, advance=.FALSE.)
      if (itime == 1) then
        budtxtarray(ibudterm) = budtxt
        packagenamearray(ibudterm) = dstpackagename
        if (trim(adjustl(budtxt)) == 'FLOW-JA-FACE' .and. &
            srcmodelname == dstmodelname) then
          internalflow(ibudterm) = 1
        else
          internalflow(ibudterm) = 0
        end if
      else
        if (budtxt /= budtxtarray(ibudterm) .or. &
            dstpackagename /= packagenamearray(ibudterm)) then
          errmsg = 'Expecting '//trim(packagenamearray(itime))//'-'// &
                   trim(budtxtarray(itime))//' but found '//trim(dstpackagename) &
                   //'-'//trim(budtxt)
          call store_error(errmsg, terminate=.TRUE.)
        end if
      end if
      !
      ! -- Accumulate flow terms (or set ich for constant heads)
      if (internalflow(ibudterm) == 1) then
        if (hasiaja) then
          call flowiaja_accumulate(ia, ja, flowja)
        else
          call flowja_accumulate(nodesrc, nodedst, flowdata)
        end if
      else
        if (trim(adjustl(budtxt)) == 'CONSTANT HEAD') then
          call flowch_setich(ibudterm, nodesrc)
          foundchd = .true.
        else
          call flow_accumulate(ibudterm, nodesrc, flowdata)
        end if
      end if
      !
    end do
    write (iout, '(a)') 'Done reading records from budget file'
    !
    ! -- Now that all constant heads read, can process budgets for them
    if (hasiaja .and. foundchd) then
      call flowch_accumulate(ia, ja, flowja)
    end if
    !
    ! -- Write information for this time
    call zoneoutput_write(itime, kstp, kper, delt, totim, nbudterms, nmznfl, &
                          vbvl, vbznfl, packagenamearray, budtxtarray, &
                          internalflow)
    itime = itime + 1

  end do timeloop
  !
  ! -- Finalize
  call write_message(text=cdot)
  call budgetdata_finalize()
  call zoneoutput_finalize()
  call zone_finalize()
  !
  ! -- return
  return
end subroutine process_budget

subroutine parse_command_line(fnam, flst, fcsv)
! ******************************************************************************
! Parse command line arguments
!   Assign zone budget name file as first argument.
!   Assign file names for list and csv files based on root name of name file
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- modules
  use KindModule
  implicit none
  ! -- dummy
  character(len=*), intent(inout) :: fnam
  character(len=*), intent(inout) :: flst
  character(len=*), intent(inout) :: fcsv
  ! -- local
  integer(I4B) :: icountcmd
  integer(I4B) :: istart
  integer(I4B) :: istop
  integer(I4B) :: i
! ------------------------------------------------------------------------------
  !
  ! -- assign fnam to first command line argument
  icountcmd = command_argument_count()
  if (icountcmd > 0) then
    call get_command_argument(1, fnam)
  end if
  !
  ! -- Set lst and csv file names by replacing fnam suffix with .lst
  istart = 0
  istop = len_trim(fnam)
  do i = istop, 1, -1
    if (fnam(i:i) == '.') then
      istart = i
      exit
    end if
  end do
  if (istart == 0) istart = istop + 1
  !
  ! -- Create flst name
  flst = fnam(1:istart)
  istop = istart + 3
  flst(istart:istop) = '.lst'
  !
  ! -- Create fcsv name
  fcsv = fnam(1:istart)
  istop = istart + 3
  fcsv(istart:istop) = '.csv'
  !
  ! -- Return
  return
end subroutine parse_command_line

