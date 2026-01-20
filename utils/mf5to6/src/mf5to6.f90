program mf5to6
  ! Mf5to6, a MODFLOW-2005 to MODFLOW 6 converter
  use ConstantsModule, only: LINELENGTH, MAXCHARLEN
  use ConverterCommonModule, only: SupportPreproc
  use ExchangeModule, only: ExchangeType
  use ExchangeWriterModule, only: ExchangeWriterType
  use GLOBAL, only: IOUT
  use GlobalVariablesModule, only: prognamconv, prognamlong, mfvnam, &
                                   niunit, ilgr, ilunit, ngrids, &
                                   masteridomain, verbose
  use GwfLgrSubsModule, only: GETNAMFIL, GETNAMFILLGR
  use InputOutputModule, only: GetUnit
  use ListModule, only: ListType
  use ModelConverterModule, only: ModelConverterType
  use MoverModule, only: MoverType, GetMoverFromList
  use SfrPackageWriterModule, only: SfrPackageWriterType, &
                                    GetSfrPackageWriter, &
                                    AllSfrPkgWriters
  use SimFileWriterModule, only: SimFileWriterType
  use SimPHMFModule, only: ustop
  use SimListVariablesModule, only: SimMovers
  use UtilitiesModule, only: GetArgs, PhmfOption
  !
  implicit none
  integer :: iexg, igrid, ispw, iu
  integer :: nc, nr, nl, nexg, nsfrpw
  type(ModelConverterType), pointer :: parentConverter => null()
  type(ModelConverterType), pointer :: childConverter => null()
  type(ListType),           pointer :: modelConverters => null()
  type(ExchangeWriterType), pointer :: exchangeWriter => null()
  type(ExchangeType), pointer :: exchange => null()
  type(SimFileWriterType) :: SimFileWriter
  type(SfrPackageWriterType), pointer :: SfrWriter => null()
  class(*),                 pointer :: obj => null()
  character(len=MAXCHARLEN) :: fname, namfil, basnam, basnamtemp, msg
  logical :: WriteDisFile
  ! formats
10 format(a,a,i0)
20 format(a,i0)
30 format(/,1x,a)
490 format(a,a)
  !
  write(*,*)
  write(*,*)trim(prognamlong)
  write(*,*)trim(MFVNAM)
  !
  fname = ''
  namfil = ''
  basnam = ''
  basnamtemp = ''
  allocate(modelConverters)
  allocate(exchangeWriter)
  allocate(SimMovers)
  SimMovers%name = 'SimMovers'
  call SimFileWriter%InitializeSimWriter()
  ! The ExchangeWriter needs access to SimFileWriter%MvrWriter
  call exchangeWriter%InitializeExchangeWriter(SimFileWriter%MvrWriter)
  !
  ! Open name file and read first entry to determine if it's an LGR model.
  call GETNAMFIL(ilgr,ngrids,namfil,ilunit,basnam)
  !
  ! Check command line for "-phmf" option, which turns on support for the
  ! PreHeadsMF head-observations preprocessor.
  call PhmfOption(SupportPreproc)
  !
  ! Ned todo:
  ! Get name of options file from command line, open it, and assign options.
  ! One option will be path to PostObsMF executable, so converter can
  ! provide a command-prompt instruction that will run PostObsMF, or maybe
  ! generate a batch or python file (could also be an option) that would
  ! run PostObsMF (twice, if there are multilayer head observations).
  SimFileWriter%BaseName = basnam
  if (ilgr > 0) then
    ! LGR is active; read and initialize parent and all children.
    DO igrid = 1, ngrids
      ! Read name file from LGR name file
      call getnamfillgr(ilunit,fname,igrid)
      ! Initialize and convert the model associated with the current grid.
      if (igrid == 1) then
        ! parent
        basnamtemp = basnam
      else
        ! child
        basnamtemp = ''
        write(basnamtemp,10) trim(basnam), '_child_',igrid-1
      endif
      !
      ! Set up parent/child relations
      if (igrid == 1) then
        ! parent
        allocate(parentConverter)
        allocate(parentConverter%ChildConverters)
        ! Initialize and convert the parent model.
        call parentConverter%InitializeModel(fname, basnamtemp, igrid)
        WriteDisFile = .false.
        call parentConverter%ConvertModel(WriteDisFile)
        call parentConverter%model%SavePointers()
        ! Add parent converter to list of converters.
        obj => parentConverter
        call modelConverters%Add(obj)
      else
        ! child
        allocate(childConverter)
        ! Assign pointer to parent in this child.
        childConverter%ParentConverter => parentConverter
        ! Initialize and convert the child model.
        call childConverter%InitializeModel(fname, basnamtemp, igrid)
        WriteDisFile = .true.
        call childConverter%ConvertModel(WriteDisFile)
        call childConverter%model%SavePointers()
        ! Add child to list of children in parent.
        obj => childConverter
        call parentConverter%ChildConverters%Add(obj)
        ! Add child converter to list of converters.
        call modelConverters%Add(obj)
      endif
    enddo
    !
    call exchangeWriter%DefineExchanges(modelConverters)
!****************************************************************
!  Build Idomain for parent model with child areas deactivated
!****************************************************************
    call parentConverter%model%InitializeIdomain()
    nc = parentConverter%model%Ncol
    nr = parentConverter%model%Nrow
    nl = parentConverter%model%Nlaynew
    nexg = exchangeWriter%Exchanges%Count()
    do iexg=1,nexg
      obj => exchangeWriter%Exchanges%GetItem(iexg)
      select type (obj)
      type is (ExchangeType)
        exchange => obj
        call exchange%ModifyIdomainForChild(nc, nr, nl, &
                          parentConverter%model%Idomain)
      end select
    enddo
    !
    ! Write the rest of the DIS file, including Idomain
    fname = parentConverter%model%DisWriter%fileobj%FName
    iu = GetUnit()
    open(unit=iu, file=fname, Access = 'APPEND',Status='OLD')
    parentConverter%model%DisWriter%fileobj%IUnit = iu
    call parentConverter%model%PointToGrid()
    call parentConverter%model%DisWriter%WriteDisFile(nc, nr, nl, &
                                   parentConverter%model%Idomain)
    !
    ! Now that idomain for parent has been constructed, repeat
    ! conversion of parent model, eliminating stress boundaries
    ! within child domains.
    call parentConverter%model%Mf6Files%CloseAll()
    call parentConverter%model%Mf6Files%files%Clear(.true.)
    call parentConverter%model%Mf2005Files%CloseAll()
    call parentConverter%model%Mf2005Files%files%Clear(.true.)
    call parentConverter%model%PackageWriters%Clear(.true.)
    call parentConverter%model%IbChdWriter%IbChdList%Clear(.true.)
    parentConverter%model%IbChdWriter%Active = .false.
    call parentConverter%model%FhbWriter%ChdWriter%IbChdList%Clear(.true.)
    parentConverter%model%FhbWriter%ChdWriter%Active = .false.
    close(ilunit)
    ilgr = 0
    masteridomain => parentConverter%model%Idomain
    ! Do not need to rewrite the DIS file
    WriteDisFile = .false.
    close(iout)
    iout = 0
!****************************************************************
!  Reconvert parent model
!****************************************************************
    msg = 'Repeating conversion of parent model after modifications.'
    write(*,30)trim(msg)
    call parentConverter%ConvertModel(WriteDisFile)
    !
    ! Remove the first sfr writer if the parent SFR package was 
    ! recreated.  This was added 3/16/2021 to fix intermittent
    ! memory issues with the lgrex1.lgr test problem.
    nsfrpw = AllSfrPkgWriters%Count()
    if (nsfrpw >= 3) then
      call AllSfrPkgWriters%RemoveNode(1, .false.)
    end if
    !
    ! Build SFR-SFR movers
    nsfrpw = AllSfrPkgWriters%Count()
    do ispw=1,nsfrpw
      SfrWriter => GetSfrPackageWriter(ispw)
      if (associated(SfrWriter)) then
        if (verbose) then
          write(*,*)'Building mover(s) for ',trim(SfrWriter%Segments%name)
        endif
        call SfrWriter%BuildSimMovers()
      endif
    enddo
    !
    ! Assemble movers defined in parent
    call parentConverter%ValidateMovers('SIMULATION')
    !
    ! Write mover file of simulation scope
    call SimFileWriter%WriteMoverFile()
    !
    ! Write GWF-GWF Exchange files
    call exchangeWriter%WriteExchangeFiles()
    !
    ! Write simulation name file.
    call SimFileWriter%WriteSimFile(modelConverters, exchangeWriter)
  else
    ! LGR is not active; just initialize and convert the model.
    allocate(parentConverter)
    call parentConverter%InitializeModel(namfil, basnam, 1)
    WriteDisFile = .true.
    call parentConverter%ConvertModel(WriteDisFile)
    !
    ! Write simulation name file.
    call SimFileWriter%WriteSimFile(parentConverter)
  endif
  !
  call ustop()
  !
end program mf5to6
