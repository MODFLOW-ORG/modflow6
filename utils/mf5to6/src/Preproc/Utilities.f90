module UtilitiesModule

  use ConstantsModule, only: MAXCHARLEN, DZERO, MAXCHARLEN
  use GlobalVariablesModule, only: optfile, PathToPostObsMf, ScriptType, &
                                   verbose, echo
  use InputOutputModule, only: GetUnit, openfile, UPCASE, URWORD
  use SimPHMFModule, only: store_error, store_note, store_warning, ustop

  private
  public :: GetArgs, Write1Drel, Write2Drel, Write3Drel, &
            Write2Dint, Write3Dint, &
            ConstantInt1D, ConstantInt2D, ConstantInt3D, &
            ConstantReal1D, ConstantReal2D, ConstantReal3D, &
            BuildArrayFormat, Write1dValues, &
            Write2dValues, Write3dValues, findcell, &
            close_file, GreaterOf, GreatestOf, RemoveElement, &
            get_extension, count_file_records, &
            CalcContribFactors, PhmfOption

  interface RemoveElement
    module procedure :: remove_element_int
  end interface RemoveElement

contains

  subroutine GetArgs(fname, basename)
    !*******************************************************************
    ! Get the name of the MF2005 name file and base-name
    ! for MF6 input.
    ! Also, open an output file for conversion messages on iout.
    !*******************************************************************
    !    SPECIFICATIONS:
    !-------------------------------------------------------------------
    implicit none
    ! dummy arguments
    character(len=*),           intent(inout) :: fname
    character(len=*), optional, intent(inout) :: basename
    ! local variables
    integer :: istat, narg
    character(len=100)  :: text
    logical :: lex
    ! format
    10 format('Usage:  MF5to6  mf2005-name-file', &
              '  mf6-base-name',//, &
              'Or enter name-file and base-name in response to ', &
              'prompts:',/)
    20 format('Usage:  MF5to6  mf2005-name-file', &
              //, &
              'Or enter name-file in response to ', &
              'prompt:',/)
    !-------------------------------------------------------------------
    !
    ! Get name file and base-name from command line or user input.
    fname = ''
    if (present(basename)) then
      basename = ''
    endif
    optfile = ''
    narg = command_argument_count()
    if (narg > 0) then
      ! First argument is MF2005 name file
      call get_command_argument(1, fname, status=istat)
      if (istat/=0) then
        call store_error('Error getting command argument')
        goto 999
      endif
      if (present(basename)) then
        if (narg > 1) then
          ! Second argument is MF6 base-name (or option, like -phmf)
          call get_command_argument(2, basename, status=istat)
          if (basename(1:1) == '-') basename = ''
          if (istat/=0) then
            call store_error('Error getting command argument')
            goto 999
          endif
        else
          ! Get MF6 base-name interactively
          text = ' Enter the base-name for MODFLOW 6 input: '
          call getfilename(text, basename)
        endif
      endif
    else
      ! Explain usage
      if (present(basename)) then
        write(*,10)
      else
        write(*,20)
      endif
      ! Get MF2005 name file interactively
      text = ' Enter the name of a MODFLOW-2005, MODFLOW-LGR, or' // &
             ' MODFLOW-NWT NAME file: '
      call getfilename(text, fname)
      if (present(basename)) then
        ! Get MF6 base-name interactively
        text = ' Enter the base-name for MODFLOW 6 input: '
        call getfilename(text, basename)
        ! Get mf5to6 options file name interactively
!        text = ' Enter the name of an MF5to6 Options file: '
!        call getfilename(text, optfile)
      endif
    endif
    !
    ! Look for name of options file
    !if (narg > 2) then
    !  ! Third argument is mf5to6 options file
    !  call get_command_argument(3, optfile, status=istat)
    !  if (istat/=0) then
    !    call store_error('Error getting command argument')
    !    goto 999
    !  endif
    !endif
    !
!    ! look for NOECHO
!    if (narg > 2) then
!      call get_command_argument(3, arg3, status=istat)
!      if (istat/=0) then
!        call store_error('Error getting command argument')
!        goto 999
!      endif
!      call UPCASE(arg3)
!      select case (arg3)
!      case ('-NOECHO')
!        echo = .false.
!        msg = 'NOECHO option found on command line. ' // &
!              'Echoing of package input is turned off.'
!        call store_note(msg)
!      case ('-ECHO')
!        echo = .true.
!        msg = 'ECHO option found on command line. ' // &
!              'Echoing of package input is turned on.'
!        call store_note(msg)
!      case ('-VERBOSE')
!        verbose = .true.
!        msg = '-VERBOSE option found on command line. '
!        call store_note(msg)
!      case default
!        msg = 'Command-line option not recognized: ' // trim(arg3)
!        call store_warning(msg)
!      end select
!    endif
    !
    ! Check for existence of MF2005 name file
    inquire(file=fname, exist=lex)
    if (.not. lex) then
      text = 'Name file "' // trim(fname) // '" does not exist.'
      call store_error(text)
      call ustop()
    endif
    !
    ! Ensure that base-name is not blank
    if (present(basename)) then
      if (basename=='') then
        text = 'Base-name for MODFLOW 6 input is blank.'
        call store_error(text)
        call ustop()
      endif
    endif
    !
    return
    !
    ! Handle errors
999 continue
    call ustop()
  end subroutine GetArgs

  subroutine PhmfOption(phmfFound)
    ! Returns phmmfFound as true if "-phmf" is found on command line.
    ! If "-phmf" is found, the following  word is assigned to the
    ! optfile variable, which stores the name of a PreHeadsMF options
    ! file.
    !
    ! dummy
    logical, intent(inout) :: phmfFound
    ! local
    integer :: i, narg
    character(len=MAXCHARLEN) :: arg
    !
    phmfFound = .false.
    narg = command_argument_count()
    if (narg > 0) then
      do i=1,narg
        call get_command_argument(i, arg)
        call upcase(arg)
        if (arg == '-PHMF') then
          phmfFound = .true.
          if (i < narg) then
            call get_command_argument(i+1, optfile)
          endif
          exit
        endif
      enddo
    endif
    !
    return
  end subroutine PhmfOption

  subroutine getfilename(text, filname)
    ! Get a filename in response to a text prompt
    use InputOutputModule, only: URWORD
    implicit none
    ! dummy arguments
    character(len=*), intent(in) :: text
    character(len=*), intent(inout) :: filname
    ! local variables
    integer :: icol, istart, istop, n
    double precision :: r
    !
    write (*,'(a)') text
    read (*,'(a)') filname
    icol = 1
    call urword(filname,icol,istart,istop,0,n,r,0,0)
    filname=filname(istart:istop)
    return
  end subroutine getfilename

  subroutine Write1Drel(mout, n, array, constant, constantval, name, &
                        writenameline, iprn)
    ! Write input to define a 1D real array
    implicit none
    ! dummy arguments
    integer, intent(in)                        :: mout, n
    double precision, intent(in), dimension(n) :: array
    logical, intent(in)                        :: constant
    double precision, intent(in)               :: constantval
    character(len=*), intent(in)               :: name
    logical, intent(in)                        :: writenameline
    ! optional arguments
    integer, intent(in), optional              :: iprn
    ! local variables
    ! formats
    10 format(2x,a)
    20 format(4x,'CONSTANT',2x,g16.9,4x,a)
    30 format(4x,'INTERNAL  FACTOR',2x,g16.9,2x,'IPRN',2x,i0,4x,a)
    !
    ! write name if requested, indented 2 spaces
    if (writenameline) then
      write(mout,10) trim(name)
    endif
    !
    ! write array control record
    if (constant) then
      write(mout,20) constantval, trim(name)
    else
      write(mout,30) 1.0, iprn, trim(name)
    endif
    !
    ! write array values
    if (.not. constant) then
      call Write1dValues(mout, n, array)
    endif
    !
    return
  end subroutine Write1Drel

  subroutine Write2Drel(mout, nrow, ncol, array, constant, constantval, name, &
                        writenameline, iprn)
    ! Write input to define a 2D real array
    implicit none
    ! dummy arguments
    integer, intent(in)                    :: mout, nrow, ncol
    double precision, intent(in), dimension(ncol,nrow) :: array
    logical, intent(in)                    :: constant
    double precision, intent(in)           :: constantval
    character(len=*), intent(in)           :: name
    logical, intent(in)                    :: writenameline
    ! optional arguments
    integer, intent(in), optional          :: iprn
    ! local variables
    ! formats
    10 format(2x,a)
    20 format(4x,'CONSTANT',2x,g14.7)
    30 format(4x,'INTERNAL  FACTOR',2x,g14.7,2x,'IPRN',2x,i0)
    !
    ! write name if requested, indented 2 spaces
    if (writenameline) then
      write(mout,10) trim(name)
    endif
    !
    ! write array control record
    if (constant) then
      write(mout,20) constantval !, trim(name)
    else
      write(mout,30) 1.0, iprn   !, trim(name)
    endif
    !
    ! write array values
    if (.not. constant) then
      call Write2dValues(mout,nrow,ncol,array)
    endif
    !
    return
  end subroutine Write2Drel

  subroutine Write3Drel(mout, nrow, ncol, nlay, array, constant, &
                        constantval, name, writenameline, iprn)
    ! Write input to define a 3D real array
    implicit none
    ! dummy arguments
    integer, intent(in)                         :: mout, nrow, ncol, nlay
    double precision, intent(in), dimension(ncol,nrow,nlay) :: array
    logical, intent(in)                         :: constant
    double precision, intent(in)                            :: constantval
    character(len=*), intent(in)                :: name
    logical, intent(in)                         :: writenameline
    ! optional arguments
    integer, intent(in), optional               :: iprn
    ! local variables
    integer :: k
    ! formats
    10 format(2x,a)
    20 format(4x,'CONSTANT',2x,g16.9,4x,a)
    30 format(4x,'INTERNAL  FACTOR',2x,g16.9,2x,'IPRN',2x,i0,4x,a)
    !
    ! write name if requested, indented 2 spaces
    if (writenameline) then
      write(mout,10) trim(name)
    endif
    !
    ! write array control record
    if (constant) then
      write(mout,20) constantval, trim(name)
    else
      write(mout,30) 1.0, iprn, trim(name)
    endif
    !
    ! write array values
    if (.not. constant) then
      do k=1,nlay
        call Write2dValues(mout,nrow,ncol,array)
      enddo
    endif
    !
    return
  end subroutine Write3Drel

  subroutine ConstantReal1D(idim1, array, constant, constantVal)
    implicit none
    ! dummy
    integer, intent(in) :: idim1
    double precision, dimension(idim1), intent(in) :: array
    logical :: constant
    double precision :: constantVal
    ! local
    integer :: j
    !
    constant = .true.
    constantVal = array(1)
    do j=1,idim1
      if (array(j) /= constantVal) then
        constant = .false.
        constantVal = -999.99
        exit
      endif
    enddo
    !
    return
  end subroutine ConstantReal1D

  subroutine ConstantReal2D(idim1, idim2, array, constant, constantVal)
    implicit none
    ! dummy
    integer, intent(in) :: idim1, idim2
    double precision, dimension(idim1, idim2), intent(in) :: array
    logical :: constant
    double precision :: constantVal
    ! local
    integer :: i, j
    !
    constant = .true.
    constantVal = array(1,1)
    do i=1,idim2
      do j=1,idim1
        if (array(j,i) /= constantVal) then
          constant = .false.
          constantVal = -999.99
          exit
        endif
      enddo
    enddo
    !
    return
  end subroutine ConstantReal2D

  subroutine ConstantReal3D(idim1, idim2, idim3, array, constant, constantVal)
    implicit none
    ! dummy
    integer, intent(in) :: idim1, idim2, idim3
    double precision, dimension(idim1, idim2, idim3), intent(in) :: array
    logical :: constant
    double precision :: constantVal
    ! local
    integer :: i, j, k
    !
    constant = .true.
    constantVal = array(1,1,1)
    do k=1,idim3
      do i=1,idim2
        do j=1,idim1
          if (array(j,i,k) /= constantVal) then
            constant = .false.
            constantVal = -999.99
            exit
          endif
        enddo
      enddo
    enddo
    !
    return
  end subroutine ConstantReal3D

  subroutine ConstantInt1D(idim1, array, constant, constantVal)
    implicit none
    ! dummy
    integer, intent(in) :: idim1
    integer, dimension(idim1), intent(in) :: array
    logical, intent(inout) :: constant
    integer, intent(inout) :: constantVal
    ! local
    integer :: j
    !
    constant = .true.
    constantVal = array(1)
    do j=1,idim1
      if (array(j) /= constantVal) then
        constant = .false.
        constantVal = -999
        exit
      endif
    enddo
    !
    return
  end subroutine ConstantInt1D

  subroutine ConstantInt2D(idim1, idim2, array, constant, constantVal)
    implicit none
    ! dummy
    integer, intent(in) :: idim1, idim2
    integer, dimension(idim1, idim2), intent(in) :: array
    logical, intent(inout) :: constant
    integer, intent(inout) :: constantVal
    ! local
    integer :: i, j
    !
    constant = .true.
    constantVal = array(1,1)
    do i=1,idim2
      do j=1,idim1
        if (array(j,i) /= constantVal) then
          constant = .false.
          constantVal = -999
          exit
        endif
      enddo
    enddo
    !
    return
  end subroutine ConstantInt2D

  subroutine ConstantInt3D(idim1, idim2, idim3, array, constant, constantVal)
    implicit none
    ! dummy
    integer, intent(in) :: idim1, idim2, idim3
    integer, dimension(idim1, idim2, idim3), intent(in) :: array
    logical, intent(inout) :: constant
    integer, intent(inout) :: constantVal
    ! local
    integer :: i, j, k
    !
    constant = .true.
    constantVal = array(1,1,1)
    do k=1,idim3
      do i=1,idim2
        do j=1,idim1
          if (array(j,i, k) /= constantVal) then
            constant = .false.
            constantVal = -999
            exit
          endif
        enddo
      enddo
    enddo
    !
    return
  end subroutine ConstantInt3D

  subroutine Write2Dint(mout, nrow, ncol, array, constant, constantval, name, &
                        writenameline, iprn)
    ! Write input to define a 2D real array
    implicit none
    ! dummy arguments
    integer, intent(in)                       :: mout, nrow, ncol
    integer, intent(in), dimension(ncol,nrow) :: array
    logical, intent(in)                       :: constant
    integer, intent(in)                       :: constantval
    character(len=*), intent(in)              :: name
    logical, intent(in)                       :: writenameline
    ! optional arguments
    integer, intent(in), optional             :: iprn
    ! local variables
    integer :: i, j
    ! formats
    10 format(2x,a)
    20 format(4x,'CONSTANT',2x,i0,2x,a)
    30 format(4x,'INTERNAL  FACTOR',2x,i0,2x,'IPRN',2x,i0,4x,a)
    40 format(20(1x,i4))
    !
    ! write name if requested, indented 2 spaces
    if (writenameline) then
      write(mout,10) trim(name)
    endif
    !
    ! write array control record
    if (constant) then
      write(mout,20) constantval, trim(name)
    else
      write(mout,30) 1, iprn, trim(name)
    endif
    !
    ! write array values
    if (.not. constant) then
      do i=1,nrow
        write(mout,40) (array(j,i),j=1,ncol)
      enddo
    endif
    !
    return
  end subroutine Write2Dint

  subroutine Write3Dint(mout, nrow, ncol, nlay, array, constant, constantval, name, &
                        writenameline, iprn)
    ! Write input to define a 2D real array
    implicit none
    ! dummy arguments
    integer, intent(in)                             :: mout, nrow, ncol, nlay
    integer, intent(in), dimension(ncol,nrow, nlay) :: array
    logical, intent(in)                             :: constant
    integer, intent(in)                             :: constantval
    character(len=*), intent(in)                    :: name
    logical, intent(in)                             :: writenameline
    ! optional arguments
    integer, intent(in), optional                   :: iprn
    ! local variables
    integer :: i, j, k
    ! formats
    10 format(2x,a)
    20 format(4x,'CONSTANT',2x,i0,2x,a)
    30 format(4x,'INTERNAL  FACTOR',2x,i0,2x,'IPRN',2x,i0,4x,a)
    40 format(20(1x,i4))
    !
    ! write name if requested, indented 2 spaces
    if (writenameline) then
      write(mout,10) trim(name)
    endif
    !
    ! write array control record
    if (constant) then
      write(mout,20) constantval, trim(name)
    else
      write(mout,30) 1, iprn, trim(name)
    endif
    !
    ! write array values
    if (.not. constant) then
      do k=1,nlay
        do i=1,nrow
          write(mout,40) (array(j,i,k),j=1,ncol)
        enddo
      enddo
    endif
    !
    return
  end subroutine Write3Dint

  subroutine Write1dValues(mout, ndim, array)
    implicit none
    ! dummy arguments
    integer,          intent(in) :: mout, ndim
    double precision, intent(in), dimension(ndim) :: array
    ! local
    integer :: i
    ! format
    40 format(10(1x,g16.9))
    !
    ! write array values
    write(mout,40) (array(i),i=1,ndim)
    !
    return
  end subroutine Write1dValues

  subroutine Write2dValues(mout, nrow, ncol, array)
    implicit none
    ! dummy arguments
    integer, intent(in)                    :: mout, nrow, ncol
    double precision, intent(in), dimension(ncol,nrow) :: array
    ! local
    integer :: i, j
    ! format
    40 format(10(1x,g14.7))
    !
    ! write array values
    do i=1,nrow
      write(mout,40) (array(j,i),j=1,ncol)
    enddo
    !
    return
  end subroutine Write2dValues

  subroutine Write3dValues(mout, nlay, nrow, ncol, array)
    implicit none
    ! dummy arguments
    integer, intent(in) :: mout, nlay, nrow, ncol
    double precision, intent(in), dimension(ncol,nrow,nlay) :: array
    ! local
    integer :: i, j, k
    ! format
    40 format(10(1x,g16.9))
    !
    ! write array values
    do k=1,nlay
      do i=1,nrow
        write(mout,40) (array(j,i,k),j=1,ncol)
      enddo
    enddo
    !
    return
  end subroutine Write3dValues

! Ned todo: add subroutines Write1Dint?

  subroutine BuildArrayFormat(DataWidth, NDigits, LineLen, Fmat, FieldsPerLine)
    implicit none
    ! dummy
    integer,          intent(in)  :: DataWidth, NDigits, LineLen
    character(len=*), intent(out) :: Fmat
    integer,          intent(out) :: FieldsPerLine
    ! local
    integer :: FieldWidth
    character(len=10) :: cnum
    ! formats
    10  format(i0)
    !
    if (DataWidth < NDigits) then
      call store_error('In BuildArrayFormat, DataWidth < Digits')
      call ustop()
    endif
    FieldWidth = DataWidth + 1
    FieldsPerLine = LineLen / FieldWidth
    Fmat = ''
    cnum = ''
    write(cnum,10) FieldsPerLine
    Fmat = '(' // trim(cnum) // '(g'
    cnum = ''
    write(cnum,10) DataWidth
    fmat = trim(fmat) // trim(cnum) // '.'
    cnum = ''
    write(cnum,10) NDigits
    fmat = trim(fmat) // trim(cnum) // ',1x))'
    !
    return
  end subroutine BuildArrayFormat

  integer function GreaterOf(ival1, ival2)
    implicit none
    ! dummy
    integer, intent(in) :: ival1, ival2
    !
    if (ival1 >= ival2) then
      GreaterOf = ival1
    else
      GreaterOf = ival2
    endif
    return
  end function GreaterOf
  
  integer function GreatestOf(ival1, ival2, ival3, ival4, ival5)
    implicit none
    ! dummy
    integer, intent(in) :: ival1, ival2, ival3
    integer, intent(in), optional :: ival4, ival5
    ! local
    integer :: ivalmax
    !
    ivalmax = ival1
    if (ival2 > ivalmax) then
      ivalmax = ival2
    endif
    if (ival3 > ivalmax) then
      ivalmax = ival3
    endif
    if (present(ival4)) then
      if (ival4 > ivalmax) then
        ivalmax = ival4
      endif
    endif
    if (present(ival5)) then
      if (ival5 > ivalmax) then
        ivalmax = ival5
      endif
    endif
    !
    GreatestOf = ivalmax
    return
  end function GreatestOf

  subroutine findcell(node,nrow,ncol,nlay,i,j,k)
    ! Return row (I), column (J), and layer (K) that correspond to a
    ! node number
    implicit none
    integer, intent(in)  :: node, nrow, ncol, nlay
    integer, intent(out) :: i, j, k
    integer :: n, nodes, nrc
    !
    nrc=nrow*ncol
    nodes=nrc*nlay
    !
    if (node.le.0 .or. node.gt.nodes) then
      k = 0
      i = 0
      j = 0
    else
      k=(node-1)/nrc+1  ! Layer index
      n=node-(k-1)*nrc
      i=(n-1)/ncol+1    ! Row index
      j=n-(i-1)*ncol    ! Column index
    endif
    !
    return
  end subroutine findcell
  
  logical function close_file(NUNIT,MAXTRIES,WAITSECS)
    ! Ensure that unit NUNIT is closed.
    ! If closure of unit NUNIT is confirmed,
    ! return value is true.  If not confirmed, return value is false.
    !
    ! If unit NUNIT is open, the file connected to unit NUNIT is closed,
    ! regardless whether the connected file is FNAME.  
    !
!    use global_data, only: iverb, max_string_len
!    use utilities, only: utl_case
!    use ifport
    implicit none
    !
    !   Argument-list variables
    integer,                    intent(in) :: nunit
    integer, optional,          intent(in) :: maxtries
    double precision, optional, intent(in) :: waitsecs
    !
    !   Local variables
    integer :: istat, ktry, maxtrieslocal
    double precision :: waitsecslocal
    logical :: closed, lex, lop
    character(len=MAXCHARLEN) :: ermsg, locname, msg
    !
    ! Format statements
    120 format('Warning: Unable to close file "',A,'"')
    !
    ! Initialize variables
    if (present(maxtries)) then
      maxtrieslocal = maxtries
    else
      maxtrieslocal = 10
    endif
    if (present(waitsecs)) then
      waitsecslocal = waitsecs
    else
      waitsecslocal = 0.01d0
    endif
    !
    ! First, determine if file is open on unit NUNIT
    inquire(unit=nunit,name=locname,opened=lop)
    if (lop) then
      closed = .false.
    else
      closed = .true.
      close_file = .true.
      return
    endif
    !
    if (.not. closed) then
      ! Close file
      ktry = 0
      20 continue
      inquire(file=locname,iostat=istat,exist=lex)
      if (lex) then
        inquire(nunit,opened=lop)
        if (lop) then
          close(nunit,status='KEEP',iostat=istat)
          !if (istat==0) closed = .true.
          inquire(nunit,opened=lop)
          if (lop) then
            closed = .false.
          endif
        else
          closed = .true.
        endif
        if (.not. closed .and. ktry<maxtrieslocal) then
          ktry = ktry+1
!          call utl_sleep(waitsecslocal)
          goto 20
        endif
      endif
    endif
    !
    if (closed) then
      close_file = .true.
      msg = 'Closed file: ' // trim(locname)
      if (verbose) write(*,*)trim(msg)
    else
      msg = 'Failed to close file: ' // trim(locname)
      if (verbose) write(*,*)trim(msg)
      write(ermsg,120)trim(locname)
      call store_warning(ermsg)
    endif
    !
    return
  end function close_file

  subroutine remove_element_int(iarray, iremove)
    implicit none
    ! dummy
    integer, allocatable, dimension(:), intent(inout) :: iarray
    integer, intent(in) :: iremove
    ! local
    integer :: i, idimorig, idimnew, k
    integer, allocatable, dimension(:) :: itemp
    !
    idimorig = size(iarray)
    idimnew = idimorig - 1
    allocate(itemp(idimnew))
    ! Copy array elements to be kept in temporary array
    k = 0
    do i=1,idimorig
      if (i /= iremove) then
        k = k + 1
        itemp(k) = iarray(i)
      endif
    enddo
    ! Deallocate the original array
    deallocate(iarray)
    ! Reallocate and repopulate the original array
    allocate(iarray(idimnew))
    do i=1,idimnew
      iarray(i) = itemp(i)
    enddo
    !
    return
  end subroutine remove_element_int

  subroutine get_extension(name, ext)
    ! Return extension part of name, without "." and converted to uppercase.
    implicit none
    ! dummy
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: ext
    ! local
    integer :: i, locdot, namelen
    character(len=1), parameter :: dot='.'
    !
    namelen = len_trim(name)
    locdot = 0
    do i=namelen,1,-1
      if (name(i:i) == dot) then
        locdot = i
        exit
      endif
    enddo
    !
    if (locdot > 0) then
      ext = name(locdot+1:namelen)
      call UPCASE(ext)
    else
      ext = ''
    endif
    !
    return
  end subroutine get_extension

  function count_file_records(filename) result(nrecs)
    ! Open a text file, count the number of records in it, and close the file.
    ! dummy
    character(len=*), intent(in) :: filename
    integer :: nrecs
    ! local
    integer :: iu
    character(len=1) :: ch
    logical :: ok
    !
    nrecs = 0
    iu = 0
    call openfile(iu, 0, filename, 'CSV')
    do
      read(iu,'(a)',end=20)ch
      nrecs = nrecs + 1
    enddo
    20 continue
    ok = close_file(iu)
    !
    return
  end function count_file_records

  subroutine CalcContribFactors(bilinear, d1c, d1p, d2c, d2p, ibound12, &
                                cf1, cf2, cf12)
    ! Calculate contributing factors using either triangular or bilinear
    ! interpolation. 
    ! If this%bilinear is true and off-corner cell is active, bilinear is used.
    ! If this%bilinear is false, triangular is used.
    ! If off-corner cell is inactive, triangular is used.
    implicit none
    ! dummy
    logical, intent(in) :: bilinear
    double precision, intent(in) :: d1c, d1p, d2c, d2p
    ! Where: d1c = Distance from parent node (in cell where ghost node is
    !              needed) to child node in first dimension.
    !        d1p = Distance from parent node to contributing node in first
    !              dimension.
    !        d2c, d2p = as above, but for second dimension.
    !
    integer, intent(in) :: ibound12   ! Ibound for off-corner cell
    ! CF1 is contributing factor for contributing cell
    !     adjacent in dimension 1.
    double precision, intent(inout) :: cf1   
    ! CF2 is contributing factor for contributing cell
    !     adjacent in dimension 2.
    double precision, intent(inout) :: cf2
    ! CF12 is contributing factor for contributing cell
    !      off corner in dimensions 1 and 2.
    double precision, intent(inout) :: cf12
    !
    if (bilinear .and. ibound12 /= 0) then
      call CalcContribFactorsBilinear(d1c, d1p, d2c, d2p, cf1, cf2, cf12)
    else
      cf12 = DZERO
      call CalcContribFactorsTriangular(d1c, d1p, d2c, d2p, cf1, cf2)
    endif
    !
    return
  end subroutine CalcContribFactors
  
  subroutine CalcContribFactorsTriangular(d1c, d1p, d2c, d2p, cf1, cf2)
    ! Calculate contributing factors using barycentric interpolation, 2D case.
    ! Method:
    ! https://classes.soe.ucsc.edu/cmps160/Fall10/resources/barycentricInterpolation.pdf
    !
    ! This routine can be used for 1D case by providing:
    !       d2c = 0; results in: a2 = 0, contributing factor returned in cf2 = 0 
    !                (d2p value is arbitrary, non-zero). 
    ! OR:
    !       d1c = 0; results in: a1 = 0, contributing factor returned in cf1 = 0
    !                (d1p value is arbitrary, non-zero).
    implicit none
    ! dummy
    double precision, intent(in) :: d1c, d1p, d2c, d2p
    ! Where: d1c = Distance from parent node (in cell where ghost node is
    !              needed) to child node in first dimension.
    !        d1p = Distance from parent node to contributing node in first
    !              dimension.
    !        d2c, d2p = as above, but for second dimension.
    double precision, intent(inout) :: cf1   ! contributing factor for contributing
    !                                          cell adjacent in first dimension.
    double precision, intent(inout) :: cf2   ! as above, but for second dimension
    ! local
    !
    ! Barycentric calculations are:
    !    a = DHALF * d1p * d2p   ! area of enclosing triangle
    !    a1 = DHALF * d2p * d1c  ! area of subtriangle 1
    !    a2 = DHALF * d1p * d2c  ! area of subtriangle 2
    !    cf1 = a1 / a     ! contributing factor for cell 1
    !    cf2 = a2 / a     ! contributing factor for cell 2
    !
    ! These calculations are simplified to:
    if (d1p == DZERO) then
      ! linear interpolation
      cf2 = d2c / d2p
      cf1 = DZERO
    elseif (d2p == DZERO) then
      ! linear interpolation
      cf1 = d1c / d1p
      cf2 = DZERO
    else
      ! triangular interpolation
      cf1 = d1c / d1p
      cf2 = d2c / d2p
    endif
    ! 
    return
  end subroutine CalcContribFactorsTriangular

  subroutine CalcContribFactorsBilinear(d1c, d1p, d2c, d2p, cf1, cf2, cf12)
    ! Calculate contributing factors using bilinear interpolation.
    ! Method: 
    ! http://bmia.bmt.tue.nl/people/BRomeny/Courses/8C080/Interpolation.pdf
    !
    implicit none
    ! dummy
    double precision, intent(in) :: d1c, d1p, d2c, d2p
    ! Where: d1c = Distance from parent node (in cell where ghost node is
    !              needed) to child node in first dimension.
    !        d1p = Distance from parent node to contributing node in first
    !              dimension.
    !        d2c, d2p = as above, but for second dimension.
    !
    ! CF1 is contributing factor for contributing cell
    !     adjacent in dimension 1.
    double precision, intent(inout) :: cf1   
    ! CF2 is contributing factor for contributing cell
    !     adjacent in dimension 2.
    double precision, intent(inout) :: cf2
    ! CF12 is contributing factor for contributing cell
    !      off corner in dimensions 1 and 2.
    double precision, intent(inout) :: cf12
    ! local
    double precision :: ap  ! Area of rectangle formed by parent nodes
    !
    cf1 = DZERO
    cf2 = DZERO
    cf12 = DZERO
    ap = d1p * d2p
    !
    ! Allow d1p or d2p to equal zero if not needed.
    if (d1c > DZERO) then
      if (ap > DZERO) then
        cf1 = d1c * (d2p - d2c) / ap
      else
        ! Use linear interpolation in dimension 1.
        cf1 = d1c / d1p
      endif
    endif
    if (d2c > DZERO) then
      if (ap > DZERO) then
        cf2 = d2c * (d1p - d1c) / ap
      else
        ! Use linear interpolation in dimension 2.
        cf2 = d2c / d2p
      endif
    endif
    if (d1c > DZERO .and. d2c > DZERO) then
      cf12 = d1c * d2c / ap
    endif
    ! 
    return
  end subroutine CalcContribFactorsBilinear

end module UtilitiesModule
