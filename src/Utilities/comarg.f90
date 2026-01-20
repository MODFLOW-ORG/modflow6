module CommandArguments
  use KindModule
  use ConstantsModule, only: LINELENGTH, LENBIGLINE, LENHUGELINE, &
                             VSUMMARY, VALL, VDEBUG, &
                             MVALIDATE
  use VersionModule, only: VERSION, MFVNAM, IDEVELOPMODE, &
                           FMTDISCLAIMER, write_license
  use CompilerVersion
  use SimVariablesModule, only: istdout, isim_level, &
                                simfile, simlstfile, simstdout, &
                                isim_mode, simulation_mode
  use SimModule, only: store_error, ustop
  use MessageModule, only: write_message, write_message_counter
  use InputOutputModule, only: upcase, getunit
  !
  implicit none
  !
  private
  public :: GetCommandLineArguments
  !
contains

  !> @brief Get command line arguments
  !!
  !! Subroutine to get and write information on command line arguments.
  !!
  !<
  subroutine GetCommandLineArguments()
    ! -- dummy variables
    ! -- local variables
    character(len=LINELENGTH) :: tag
    character(len=LINELENGTH) :: uctag
    character(len=LENHUGELINE) :: line
    character(len=LINELENGTH) :: clevel
    character(len=LINELENGTH) :: cmode
    character(len=LINELENGTH) :: header
    character(len=LINELENGTH) :: errmsg
    character(len=LINELENGTH) :: cexe
    character(len=LENBIGLINE) :: compiler
    character(len=20) :: cdate
    character(len=LENBIGLINE) :: coptions
    logical :: ltyp
    logical :: lexist
    logical :: lstop
    integer(I4B) :: icountcmd
    integer(I4B) :: ipos
    integer(I4B) :: ilen
    integer(I4B) :: iarg
    !
    ! -- initialize local variables
    lstop = .FALSE.
    !
    ! -- set mf6 executable name
    icountcmd = command_argument_count()
    call get_command_argument(0, cexe)
    cexe = adjustl(cexe)
    !
    ! -- find the program basename, not including the path (this should be
    !    mf6.exe, mf6d.exe, etc.)
    ipos = index(cexe, '/', back=.TRUE.)
    if (ipos == 0) then
      ipos = index(cexe, '\', back=.TRUE.)
    end if
    if (ipos == 0) then
      ipos = 1
    else
      ipos = ipos + 1
    end if
    cexe = cexe(ipos:)
    !
    ! -- write header
    call get_compile_date(cdate)
    write (header, '(a,4(1x,a),a)') &
      trim(adjustl(cexe)), '- MODFLOW', &
      trim(adjustl(VERSION)), '(compiled', trim(adjustl(cdate)), ')'
    !
    ! -- check for silent option
    do iarg = 1, icountcmd
      call get_command_argument(iarg, uctag)
      call upcase(uctag)
      if (trim(adjustl(uctag)) == '-S' .or. &
          trim(adjustl(uctag)) == '--SILENT') then
        !
        ! -- get file unit and open mfsim.stdout
        istdout = getunit()
        open (unit=istdout, file=trim(adjustl(simstdout)))
        !
        ! -- exit loop
        exit
      end if
    end do
    !
    ! -- Read remaining arguments
    iarg = 0
    do
      !
      ! -- increment iarg and determine if loop should be terminated
      iarg = iarg + 1
      if (iarg > icountcmd) then
        exit
      end if
      !
      ! -- get command line argument
      call get_command_argument(iarg, tag)
      uctag = tag
      call upcase(uctag)
      !
      ! -- skip commands without - or --
      ipos = index(uctag, '-')
      if (ipos < 1) then
        cycle
      end if
      !
      ! -- parse level string, if necessary
      clevel = ' '
      ipos = index(uctag, '--LEVEL=')
      if (ipos > 0) then
        ipos = index(tag, '=')
        ilen = len_trim(tag)
        clevel = tag(ipos + 1:ilen)
        call upcase(clevel)
        uctag = tag(1:ipos - 1)
        call upcase(uctag)
      end if
      !
      ! -- parse mode string, if necessary
      cmode = ' '
      ipos = index(uctag, '--MODE=')
      if (ipos > 0) then
        ipos = index(tag, '=')
        ilen = len_trim(tag)
        cmode = tag(ipos + 1:ilen)
        call upcase(cmode)
        uctag = tag(1:ipos - 1)
        call upcase(uctag)
      end if
      !
      ! -- evaluate the command line argument (uctag)
      select case (trim(adjustl(uctag)))
      case ('-H', '-?', '--HELP')
        lstop = .TRUE.
        call write_usage(trim(adjustl(header)), trim(adjustl(cexe)))
      case ('-V', '--VERSION')
        lstop = .TRUE.
        write (line, '(2a,2(1x,a))') &
          trim(adjustl(cexe)), ':', trim(adjustl(VERSION))
        call write_message_counter(line, skipbefore=1, skipafter=1)
      case ('-DEV', '--DEVELOP')
        lstop = .TRUE.
        write (line, '(2a,g0)') &
          trim(adjustl(cexe)), ': develop version ', ltyp
        call write_message_counter(line, skipbefore=1, skipafter=1)
      case ('-C', '--COMPILER')
        lstop = .TRUE.
        call get_compiler(compiler)
        write (line, '(2a,1x,a)') &
          trim(adjustl(cexe)), ':', trim(adjustl(compiler))
        call write_message_counter(line, skipbefore=1, skipafter=1)
      case ('-S', '--SILENT')
        write (line, '(2a,1x,a)') &
          trim(adjustl(cexe)), ':', 'all screen output sent to mfsim.stdout'
        call write_message_counter(line, skipbefore=1, skipafter=1)
      case ('-D', '--DISCLAIMER')
        lstop = .TRUE.
        call write_message('', fmt=FMTDISCLAIMER)
      case ('-P', '--PARALLEL')
        simulation_mode = 'PARALLEL'
      case ('-LIC', '--LICENSE')
        lstop = .TRUE.
        call write_license()
      case ('-CO', '--COMPILER-OPT')
        lstop = .TRUE.
        call get_compile_options(coptions)
        call write_message_counter(coptions, skipbefore=1, skipafter=1)
      case ('-L', '--LEVEL')
        if (len_trim(clevel) < 1) then
          iarg = iarg + 1
          call get_command_argument(iarg, clevel)
          call upcase(clevel)
        end if
        select case (trim(adjustl(clevel)))
        case ('SUMMARY')
          isim_level = VSUMMARY
        case ('DEBUG')
          isim_level = VDEBUG
        case default
          call write_usage(trim(adjustl(header)), trim(adjustl(cexe)))
          write (errmsg, '(2a,1x,a)') &
            trim(adjustl(cexe)), ': illegal STDOUT level option -', &
            trim(adjustl(clevel))
          call store_error(errmsg)
        end select
        !
        ! -- write message to stdout
        write (line, '(2a,2(1x,a))') &
          trim(adjustl(cexe)), ':', 'stdout output level', &
          trim(adjustl(clevel))
        call write_message_counter(line, skipbefore=1, skipafter=1)
      case ('-M', '--MODE')
        if (len_trim(cmode) < 1) then
          iarg = iarg + 1
          call get_command_argument(iarg, cmode)
          call upcase(cmode)
        end if
        select case (trim(adjustl(cmode)))
        case ('VALIDATE')
          isim_mode = MVALIDATE
        case default
          call write_usage(trim(adjustl(header)), trim(adjustl(cexe)))
          errmsg = trim(adjustl(cexe))//': illegal MODFLOW 6 '// &
                   'simulation mode option - '//trim(adjustl(cmode))
          call store_error(errmsg, terminate=.TRUE.)
        end select
        !
        ! -- write message to stdout
        line = trim(adjustl(cexe))//': MODFLOW 6 simulation mode '// &
               trim(adjustl(cmode))//'. Model input will be checked for all '// &
               'stress periods but the matrix equations will not be '// &
               'assembled or solved.'
        call write_message_counter(line, skipbefore=1, skipafter=1)
      case default
        lstop = .TRUE.
        call write_usage(trim(adjustl(header)), trim(adjustl(cexe)))
        write (errmsg, '(2a,1x,a)') &
          trim(adjustl(cexe)), ': illegal option -', trim(adjustl(tag))
        call store_error(errmsg, terminate=.TRUE.)
      end select
    end do
    !
    ! -- check if simfile exists, only if the model should be run
    if (.not. lstop) then
      inquire (file=trim(adjustl(simfile)), exist=lexist)
      if (.NOT. lexist) then
        lstop = .TRUE.
        write (errmsg, '(2a,2(1x,a))') &
          trim(adjustl(cexe)), ':', trim(adjustl(simfile)), &
          'is not present in working directory.'
        call store_error(errmsg, terminate=.TRUE.)
      end if
    end if
    !
    ! -- terminate program if lstop
    if (lstop) then
      call ustop()
    end if
    !
    ! -- write blank line to stdout
    if (icountcmd > 0) then
      call write_message('')
    end if
  end subroutine GetCommandLineArguments

  !> @brief Write command line argument usage
  !!
  !! Subroutine to write usage information for command line arguments.
  !!
  !<
  subroutine write_usage(header, cexe)
    ! -- dummy variables
    character(len=*), intent(in) :: header !< header for usage
    character(len=*), intent(in) :: cexe !< executable name
    ! -- local variables
    character(len=LINELENGTH) :: line
    ! -- format
    character(len=*), parameter :: OPTIONSFMT = &
      "(/,&
      &'Options   GNU long option   Meaning ',/,&
      &' -h, -?    --help           Show this message',/,&
      &' -v        --version        Display program version information.',/,&
      &' -dev      --develop        Display program develop option mode.',/,&
      &' -d        --disclaimer     Display program disclaimer.',/,&
      &' -p        --parallel       Run program in parallel mode.',/,&
      &' -lic      --license        Display program license information.',/,&
      &' -c        --compiler       Display compiler information.',/,&
      &' -co       --compiler-opt   Display compiler options.',/,&
      &' -s        --silent         All STDOUT to mfsim.stdout.',/,"// &
      "' -l <str>  --level <str>    STDOUT output to screen based on <str>.',/,&
      &'                            <str>=summary Limited output to STDOUT.',/,&
      &'                            <str>=debug   Enhanced output to STDOUT.',/,&
      &' -m <str>  --mode <str>     MODFLOW 6 simulation mode based on <str>.',/,&
      &'                            <str>=validate Check model input for',/,&
      &'                                           errors but do not ',/,&
      &'                                           assemble or solve matrix ',/,&
      &'                                           equations or write ',/,&
      &'                                           solution output.',/,"// &
      "'                                                                    ',/,&
      &'Bug reporting and contributions are welcome from the community. ',/,&
      &'Questions can be asked on the issues page[1]. Before creating a new',/,&
      &'issue, please take a moment to search and make sure a similar issue',/,&
      &'does not already exist. If one does exist, you can comment (most',/,&
      &'simply even with just :+1:) to show your support for that issue.',/,&
      &'                                                                    ',/,&
      &'[1] https://github.com/MODFLOW-ORG/modflow6/issues',/)"
    !
    ! -- write command line usage information to the screen
    call write_message(header)
    write (line, '(a,1x,a,15x,a,2(1x,a),2a)') &
      'usage:', cexe, 'run MODFLOW', trim(adjustl(MFVNAM)), &
      'using "', trim(adjustl(simfile)), '"'
    call write_message(line)
    write (line, '(a,1x,a,1x,a,5x,a)') &
      '   or:', cexe, '[options]', &
      'retrieve program information'
    call write_message(line)
    call write_message('', fmt=OPTIONSFMT)
  end subroutine write_usage

end module CommandArguments
