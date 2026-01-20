module ModelConverterModule
  ! A class for converting one model (MF-2005, MODFLOW-NWT,
  ! MODFLOW-USG) to MODFLOW 6.
  !
  ! For converting a system of models (e.g. for MODFLOW-LGR),
  ! each model (parent, children) would be converted by its
  ! own ModelConverterType object.

  use ConstantsModule, only: LINELENGTH, MAXCHARLEN, LENBIGLINE, &
                             LENMODELNAME, LENPACKAGENAME, LENFTYPE
  use DrnObsWriterModule, only: createDrnObsWriter, DrnObsWriterType
  use FileTypeModule, only: FileType
  use GLOBAL, only: IFREFM, IUNIT, MXITER, iout, cbcfilename
  use GlobalVariablesModule, only: prognamconv, mfvnam, niunit, cunit, &
                                   ilgr, ilunit, ngrids, PathToPostObsMf, &
                                   ScriptType
  use GlobalVariablesPHMFModule, only: prognamPHMF
  use GwfBasSubs, only: GWF2BAS7AR
  use GwfBcfSubs, only: GWF2BCF7AR
  use GWFFHBMODULE, only: SGWF2FHB7PNT, SGWF2FHB7PSV
  use GwfFhbSubs, only: GWF2FHB7AR
  use GwfLgrSubsModule, only: GWF2LGR2AR
  use GwfLpfSubs, only: GWF2LPF7AR
  use HfBSubsNwt, only: GWF2HFB7AR, GWF2HFB7UPW
  use InputOutputModule, only: GetUnit, openfile
  use ListModule, only: ListType
  use ModelModule, only: ModelType
  use ModelPackageModule, only: ModelPackageType, ConstructModelPackageType
  use MoverModule, only: MoverType, AddMoverToList, GetMoverFromList
  use ObsWriterModule, only: ObsWriterType
  use OpenSpecModule, only: ACCESS, ACTION, FORM
  use PackageWriterModule, only: PackageWriterType
  use PCGN, only: PCGN2AR
  use RivObsWriterModule, only: createRivObsWriter, RivObsWriterType
  use GhbObsWriterModule, only: createGhbObsWriter, GhbObsWriterType
  use GlobalVariablesModule, only: echo
  use SimPHMFModule, only: store_error, store_note, store_warning, ustop
  use SimListVariablesModule, only: SimMovers
  use UpwSubsModule, only: GWF2UPW1AR
  use UtilitiesModule, only: GetArgs

  implicit none

  type :: ModelConverterType
    logical :: ConversionDone = .false.
    type(ModelType),          pointer :: model => null()
    type(ModelConverterType), pointer :: ParentConverter => null()
    type(ListType),           pointer :: ChildConverters => null()
    type(ListType),           pointer :: ModelMovers
    logical :: NeedWaterMover = .false.
    character(len=MAXCHARLEN) :: PostObsScriptName = ''
  contains
    procedure :: AddMover
    procedure :: ValidateMovers
    procedure :: ConvertModel
    procedure :: GetMover
    procedure :: InitializeModel
    procedure :: LookUpModelName
    procedure :: LookUpPkgName
    procedure, private :: WritePostObsScript
  end type ModelConverterType

contains

  subroutine InitializeModel(this, namefile, basename, igrid)
    implicit none
    ! dummy
    class(ModelConverterType), intent(inout) :: this
    character(len=*),          intent(in)    :: namefile
    character(len=*),          intent(in)    :: basename
    integer, intent(in) :: igrid
    !
    allocate(this%ModelMovers)
    allocate(this%model)
    this%model%ModelMovers => this%ModelMovers
    this%model%IGrid = igrid
    call this%model%Initialize(namefile, basename)
    !
    return
  end subroutine InitializeModel

  subroutine ConvertModel(this, WriteDisFile)
    implicit none
    ! dummy
    class(ModelConverterType) :: this
    logical, intent(in) :: WriteDisFile
    ! local
    integer :: fcode, ii, inunit, iu, iu15, numobsml
    integer :: maxunit, nsol, npkg, kpkg, lenbase
    integer :: nfiles, icbc, nmov
    integer :: ibotavg
    character(len=80)         :: headng(2)
    character(len=maxcharlen) :: msg
    character(len=MAXCHARLEN) :: fname, fname15
    character(len=12)         :: filtyp15
    character(len=LENFTYPE)   :: filtype
    logical                   :: lex, lop
    class(PackageWriterType), pointer :: pkgWriPtr => null()
    type(FileType),           pointer :: cbcfil => null(), filptr => null()
    class(ObsWriterType),     pointer :: obsWriter => null()
    type(RivObsWriterType),   pointer :: rivObsWriter => null()
    type(GhbObsWriterType),   pointer :: ghbObsWriter => null()
    type(DrnObsWriterType),   pointer :: drnObsWriter => null()
    type(ModelPackageType),   pointer :: newModelPack => null()
    !
    ! formats
    20 format('New name file is "',a,'"')
    30 format('Basename exceeds length limit, which is ',i0,'.')
    40 format(a,' has generated input for PostObsMF.', &
              ' After the MODFLOW 6 model is run, PostObsMF can', &
              ' be run by executing "[PostObsMF] ',a,'" at the', &
              ' command prompt, where [PostObsMF] is a pathname ', &
              'pointing to the PostObsMF.exe executable file.')
    50 format(a,' has generated input for PostObsMF for processing ', &
              'multilayer head observations.', &
              ' After PostObsMF has been run as described in the ', &
              'preceding note, PostObsMF can', &
              ' be run again to process those multilayer ', &
              'observations by executing "[PostObsMF] ',a,'" at the', &
              ' command prompt, where [PostObsMF] is a pathname ', &
              'pointing to the PostObsMF.exe executable file.')
    !
    iout = 0
    icbc = 0
    write(*,*)
    write(*,*)
    write(*,*)'Converting model based on name file: ', &
              trim(this%model%NameFile2005)
    write(*,*)
    echo = .true.
    ! Get MF-2005 name file and MF-2015 base-name, and open conversion-report file
    lenbase = len_trim(this%model%BaseName)
    if (lenbase > MAXCHARLEN) then
      write(msg,30)MAXCHARLEN
      call store_error(msg)
      call ustop()
    endif
    !
    INUNIT = 99
    NSOL = 1
    MAXUNIT = INUNIT
    !
    ! OpenMF2005 name file.
    OPEN (UNIT=INUNIT,FILE=this%model%NameFile2005,STATUS='OLD', &
          ACTION=ACTION(1))
    !
    ! Allocate BAS; open files listed in Name file; read DIS data;
    ! read BAS file; read Zone and Multiplier arrays; read PVAL file.
    CALL GWF2BAS7AR(INUNIT,CUNIT,prognamconv,24,31,32,MAXUNIT,12, &
                    HEADNG,26,MFVNAM, this%model)
    call this%model%OrderPackageWriters()
    IF (ILGR /= 0) CALL GWF2LGR2AR(ILUNIT, this%model%NameFile2005, NGRIDS, &
                                   this%model%IGrid)
    !
    if (.not. this%ConversionDone) then
      ! Write MF6 simulation tdis file
      call this%model%TdisWriter%WriteFile()
      !
      ! Convert MF2005 flow package to NPF
      IF(IUNIT(1).GT.0) then
        CALL GWF2BCF7AR(IUNIT(1), this%model)
        if (this%model%NpfWriter%Inpfcb > 0) then
          icbc = this%model%NpfWriter%Inpfcb
        endif
      endif
      IF(IUNIT(23).GT.0) then
        CALL GWF2LPF7AR(IUNIT(23), this%model)
      endif
      !
      ! NWT input has to be read before UPW input
      ibotavg = 0
      IF(IUNIT(63).GT.0) then
        CALL GWF2NWT1AR(IUNIT(63),MXITER, IUNIT(22),this%model%IGrid,ibotavg)
        this%model%NpfWriter%Newton = .true.
        this%model%StoWriter%Newton = .true.
        if (ibotavg /= 0) this%model%NewtonUnderRelaxation = .true.
      endif
      IF(IUNIT(62).GT.0) then
        CALL GWF2UPW1AR(IUNIT(62), this%model)
        this%model%Newton = .true.
      endif


    !skip  IF(IUNIT(37).GT.0) CALL GWF2HUF7AR(IUNIT(37),IUNIT(47),
    !     1                                     IUNIT(53),IGRID)
    !      IF(IUNIT(3).GT.0) CALL GWF2DRN7AR(IUNIT(3),IGRID)
     !     IF(IUNIT(4).GT.0) CALL GWF2RIV7AR(IUNIT(4),IGRID)
    !      IF(IUNIT(5).GT.0) CALL GWF2EVT7AR(IUNIT(5),IGRID)
    !      IF(IUNIT(7).GT.0) CALL GWF2GHB7AR(IUNIT(7),IGRID)
    !      IF(IUNIT(8).GT.0) CALL GWF2RCH7AR(IUNIT(8),IGRID)
    !  IF(IUNIT(16).GT.0) then
    !    write(*,*)'Processing FHB package input...'
    !    CALL GWF2FHB7AR(IUNIT(16),IGRID)
    !  endif
    !skip  IF(IUNIT(17).GT.0) CALL GWF2RES7AR(IUNIT(17),IGRID)
    !skip  IF(IUNIT(18).GT.0) CALL GWF2STR7AR(IUNIT(18),IGRID)
    !skip  IF(IUNIT(19).GT.0) CALL GWF2IBS7AR(IUNIT(19),IUNIT(54),IGRID)
    !      IF(IUNIT(20).GT.0) CALL GWF2CHD7AR(IUNIT(20),IGRID)
      !
      ! Call GWF2HFB7AR here so that GWF2HFB7UPW can be called after it.
      IF(IUNIT(21).GT.0) then
        write(*,*)'Processing HFB package input...'
        CALL GWF2HFB7AR(IUNIT(21),this%model%IGrid)
      endif
      ! Modify conductance for HFB when using UPW.
      IF (IUNIT(62).GT.0 ) THEN
        IF(IUNIT(21).GT.0) CALL GWF2HFB7UPW(this%model%IGrid)
      END IF
      !
    !skip  IF(IUNIT(22).GT.0 .OR. IUNIT(44).GT.0) CALL GWF2LAK7AR(
    !     1             IUNIT(22),IUNIT(44),IUNIT(15),IUNIT(55),NSOL,IGRID)
    !skip  IF(IUNIT(44).GT.0) CALL GWF2SFR7AR(IUNIT(44),IUNIT(1),IUNIT(23),
    !     1                           IUNIT(37),IUNIT(15),NSOL,IOUTS,
    !     2                           IUNIT(55),IGRID)
    !skip      IF(IUNIT(44).GT.0) CALL GWF2SFR7AR(IUNIT(44),IUNIT(1),IUNIT(23),
    !skip     1                           IUNIT(37),IUNIT(15),NSOL,IOUTS,
    !skip     2                           IUNIT(62),IUNIT(55),IGRID)
    !skip  IF(IUNIT(55).GT.0) CALL GWF2UZF1AR(IUNIT(55),IUNIT(1),
    !     1                                   IUNIT(23),IUNIT(37),IGRID)
    !skip  IF(IUNIT(46).GT.0) CALL GWF2GAG7AR(IUNIT(46),IUNIT(44),
    !     1                                     IUNIT(22),IGRID)
    !      IF(IUNIT(39).GT.0) CALL GWF2ETS7AR(IUNIT(39),IGRID)
    !skip  IF(IUNIT(40).GT.0) CALL GWF2DRT7AR(IUNIT(40),IGRID)
    !skip  IF(IUNIT(54).GT.0) CALL GWF2SUB7AR(IUNIT(54),IGRID)
    !skip      IF(IUNIT(64).GT.0) CALL GWF2SWR7AR(IUNIT(64),
    !skip     2                        IUNIT(1),IUNIT(23),IUNIT(37),
    !skip     3                        IUNIT(62),IUNIT(44),IUNIT(63),IGRID)  !SWR  - JDH
    !skip      IF(IUNIT(65).GT.0) CALL GWF2SWI2AR(IUNIT(65),
    !skip     2                        IUNIT(1),IUNIT(23),IUNIT(37),IUNIT(62),
    !skip     3                        IGRID)   !SWI2  - JDH
    !
    ! SOLVERS
      IF(IUNIT(9).GT.0) CALL SIP7AR(IUNIT(9),MXITER,this%model%IGrid)
      IF(IUNIT(10).GT.0) CALL DE47AR(IUNIT(10),MXITER,this%model%IGrid)
      IF(IUNIT(13).GT.0) CALL PCG7AR(IUNIT(13),MXITER,this%model%IGrid)
    !c      IF(IUNIT(14).GT.0) CALL LMG7AR(IUNIT(14),MXITER,IGRID)
      IF(IUNIT(42).GT.0) CALL GMG7AR(IUNIT(42),MXITER,this%model%IGrid)
      IF(IUNIT(59).GT.0) &
          CALL PCGN2AR(IUNIT(59),IFREFM,MXITER,this%model%IGrid)
    !
    !      IF(IUNIT(50).GT.0) CALL GWF2MNW27AR(IUNIT(50),IGRID)
    !      IF(IUNIT(51).GT.0) CALL GWF2MNW2I7AR(IUNIT(51),IUNIT(50),IGRID)
    !      IF(IUNIT(52).GT.0) CALL GWF2MNW17AR(IUNIT(52),IUNIT(9),
    !     1                     IUNIT(10),0,IUNIT(13),
    !     2                     IUNIT(42),IUNIT(59),FNAME,IGRID)
    !      IF(IUNIT(57).GT.0) CALL GWF2SWT7AR(IUNIT(57),IGRID)
    !      IF(IUNIT(65).GT.0) CALL GWF2SWI2AR(IUNIT(65),IUNIT(1),
    !     2                     IUNIT(23),IUNIT(37),IUNIT(62),IGRID)
    !skip  IF(IUNIT(43).GT.0) CALL GWF2HYD7BAS7AR(IUNIT(43),IGRID)
    !skip  IF(IUNIT(43).GT.0 .AND. IUNIT(19).GT.0)
    !     1                   CALL GWF2HYD7IBS7AR(IUNIT(43),IGRID)
    !skip  IF(IUNIT(43).GT.0 .AND. IUNIT(54).GT.0)
    !     1                   CALL GWF2HYD7SUB7AR(IUNIT(43),IGRID)
    !skip  IF(IUNIT(43).GT.0 .AND. IUNIT(18).GT.0)
    !     1                   CALL GWF2HYD7STR7AR(IUNIT(43),IGRID)
    !skip  IF(IUNIT(43).GT.0 .AND. IUNIT(44).GT.0)
    !     1                   CALL GWF2HYD7SFR7AR(IUNIT(43),IGRID)
    !skip  IF(IUNIT(49).GT.0) CALL LMT7BAS7AR(INUNIT,CUNIT,IGRID)
    endif
    !
    npkg = this%model%PackageWriters%Count()
    !
    ! Iterate through package writers and have each writer process
    ! the allocate-related input.
    call this%model%OcWriter%ProcessAllocate(this%model%IGrid)
    do kpkg=1,npkg
      pkgWriPtr => this%model%GetPackageWriter(kpkg)
      call pkgWriPtr%ProcessAllocate(this%model%IGrid)
      newModelPack => ConstructModelPackageType(this%model%IGrid, &
                      this%model%ModelName, pkgWriPtr%PackageName)
      pkgWriPtr%ModelPack => newModelPack
      fname15 = pkgWriPtr%fileobj%FName
      filtyp15 = pkgWriPtr%fileobj%FType
      iu15 = pkgWriPtr%fileobj%IUnit
      fcode = pkgWriPtr%fileobj%FCode
      if (pkgWriPtr%Active) then
        call this%model%Mf6Files%AddFile(fname15, filtyp15, iu15, fcode, &
                                            pkgWriPtr%fileobj%PkgName)
      endif
      if (icbc==0) then
        if (associated(pkgWriPtr%ICbc)) then
          if (pkgWriPtr%ICbc > 0) then
            icbc = pkgWriPtr%ICbc
          endif
        endif
      endif
    enddo
    !
    if (WriteDisFile) then
      call this%model%InitializeIdomain()
    endif
    !
    IF(IUNIT(16).GT.0) then
      write(*,*)'Processing FHB package input...'
      call this%model%FhbWriter%ProcessAllocate(this%model%IGrid)
      call this%model%FhbWriter%PrepareFiles(this%model%Mf6Files)
      ! Remove CHD duplicates
      call this%model%RemoveFhbChdDuplicates()
      call this%model%FhbWriter%WriteFiles(this%model%IGrid)
    endif
    !
    if (this%model%IbChdWriter%Active) then
      if (.not. this%model%ChdFileAdded) then
        call this%model%Mf6Files%AddFile( &
            this%model%IbChdWriter%fileobj%FName, &
            'CHD6', this%model%IbChdWriter%fileobj%IUnit, &
            this%model%IbChdWriter%fileobj%FCode, &
            this%model%IbChdWriter%fileobj%PkgName)
           this%model%ChdFileAdded = .true.
      endif
    endif
    !
    if (icbc > 0 .and. cbcfilename == '') then
      cbcfil => this%model%Mf6Files%GetFileByUnit(icbc)
      if (associated(cbcfil)) then
        cbcfilename = cbcfil%FName
      endif
    endif
    !
    if (WriteDisFile) then
      ! Write DIS file
      call this%model%DisWriter%WriteDisFile(this%model%Ncol, this%model%Nrow, &
                          this%model%Nlaynew, this%model%Idomain)
    endif
    !
    if (.not. this%ConversionDone) then
      ! Write SMS file
      call this%model%ImsWriter%WriteFile(this%model%IGrid)
      IF(IUNIT(1).GT.0) then
        ! Write NPF and STO input
        call this%model%NpfWriter%WriteFile()
        if (this%model%NpfWriter%Inpfcb > 0) then
          icbc = this%model%NpfWriter%Inpfcb
        endif
        call this%model%StoWriter%WriteFile()
      endif
      IF(IUNIT(23).GT.0) then
        ! Write NPF and STO input
        call this%model%NpfWriter%WriteFile()
        call this%model%StoWriter%WriteFile()
      endif
      !
      ! NWT input has to be read before UPW input
      IF(IUNIT(63).GT.0) then
        this%model%NpfWriter%Newton = .true.
      endif
      IF(IUNIT(62).GT.0) then
        call this%model%NpfWriter%WriteFile()
        call this%model%StoWriter%WriteFile()
      endif
      !
      ! Write Output Control file
      call this%model%OcWriter%WriteFile(this%model%Mf6Files, &
                                         this%model%IGrid)
    endif
    !
    ! Write Head Observations file
    if (this%model%HedObsWriter%Active) then
      call this%model%HedObsWriter%InitializeObs(this%model%BaseName, 'heads')
      call this%model%Mf6Files%AddFile( &
                      this%model%HedObsWriter%fileobj%FName, &
                      this%model%HedObsWriter%fileobj%FType, &
                      this%model%HedObsWriter%fileobj%IUnit, &
                      this%model%HedObsWriter%fileobj%FCode)
      this%model%HedObsWriter%Mf6Files => this%model%Mf6Files
      call this%model%HedObsWriter%SetDisFileName          &
                      (this%model%DisWriter%fileobj%FName)
      call this%model%HedObsWriter%SetHdry(this%model%NpfWriter%Hdry)
      call this%model%HedObsWriter%SetStressPeriods &
                      (this%model%TdisWriter%StressPeriods)
      call this%model%HedObsWriter%WriteObsFile(this%model%IGrid)
      if (this%model%HedObsWriter%Preproc%Active) then
        ! Run PreHeadsMF
        write(*,'(/,a)')trim(prognamconv) // ' is invoking ' // trim(prognamPHMF) &
            // ' to preprocess head observations...'
        call this%model%HedObsWriter%Preproc%Run( &
                 this%model%HedObsWriter%PreHeadsMfFile, .false.)
        if (PathToPostObsMf == '') then
          write(msg,40)trim(prognamPHMF), &
              trim(this%model%HedObsWriter%Preproc%PostObsFilename)
          call store_note(msg)
        endif
        numobsml = this%model%HedObsWriter%MLObsList%Count()
        if (numobsml > 0) then
          call this%model%HedObsWriter%write_ml_postobs_file()
          if (PathToPostObsMf == '') then
            write(msg,50)trim(prognamconv), &
                         trim(this%model%HedObsWriter%MlPostObsFileName)
            call store_note(msg)
          endif
        endif
        ! Write script to enable user to run PostObsMF for this converted model.
        call this%WritePostObsScript()
      endif
    endif
    !
    ! Write CHD Observations file
    if (associated(this%model%IbChdWriter)) then
      if (this%model%IbChdWriter%Active) then
        if (associated(this%model%IbChdWriter%PkgObsWriter)) then
          this%model%IbChdWriter%PkgObsWriter%Mf6Files => this%model%Mf6Files
        endif
      endif
    endif
    !
    ! Look for OBS input files in list of MF2005 input files
    nfiles = this%model%Mf2005Files%NCount()
    do ii=1,nfiles
      filptr => this%model%Mf2005Files%GetFile(ii)
      filtype = filptr%FType
      select case (filtype)
      case ('RVOB')
        do kpkg=1,npkg
          pkgWriPtr => this%model%GetPackageWriter(kpkg)
          if (pkgWriPtr%PkgType == 'RIV') then
            call createRivObsWriter(rivObsWriter, this%model%BaseName, &
                                    this%model%IuRivObs)
            pkgWriPtr%PkgObsWriter => rivObsWriter
            pkgWriPtr%ObsActive = .true.
            pkgWriPtr%PkgObsWriter%Mf6Files => this%model%Mf6Files
            exit
          endif
        enddo
      case ('GBOB')
        do kpkg=1,npkg
          pkgWriPtr => this%model%GetPackageWriter(kpkg)
          if (pkgWriPtr%PkgType == 'GHB') then
            call createGhbObsWriter(ghbObsWriter, this%model%BaseName,  &
                                    this%model%IuGhbObs)
            pkgWriPtr%PkgObsWriter => ghbObsWriter
            pkgWriPtr%ObsActive = .true.
            pkgWriPtr%PkgObsWriter%Mf6Files => this%model%Mf6Files
            exit
          endif
        enddo
      case ('DROB')
        do kpkg=1,npkg
          pkgWriPtr => this%model%GetPackageWriter(kpkg)
          if (pkgWriPtr%PkgType == 'DRN') then
            call createDrnObsWriter(drnObsWriter, this%model%BaseName, &
                                    this%model%IuDrnObs)
            pkgWriPtr%PkgObsWriter => drnObsWriter
            pkgWriPtr%ObsActive = .true.
            pkgWriPtr%PkgObsWriter%Mf6Files => this%model%Mf6Files
            exit
          endif
        enddo
      case ('CHOB')
        continue
      end select
    enddo
    !
    ! Iterate through package writers and have each writer process
    ! input for one package for all stress periods.
    do kpkg=1,npkg
      pkgWriPtr => this%model%GetPackageWriter(kpkg)
      call pkgWriPtr%ProcessStressLoop(this%model%IGrid)
      call pkgWriPtr%CloseFile()
    enddo
    !
    ! Build movers (between packages of this model)
    call this%model%BuildAllModelMovers()
    !
    nmov = this%model%ModelMovers%Count()
    if (nmov > 0) then
      !
      ! Valicate model movers
      call this%ValidateMovers('MODEL')
      !
      ! Write model mover file
      this%model%MvrWriter%ModelBasename = this%model%BaseName
      call this%model%MvrWriter%ProcessAllocate(this%model%IGrid)
      call this%model%WriteMvrFile()
      !
    endif
    !
    ! Iterate through package writers and have each writer process
    ! OBS input for one package. This will write the OBS6 input
    ! file for each package that has observations.
    do kpkg=1,npkg
      pkgWriPtr => this%model%GetPackageWriter(kpkg)
      obsWriter => pkgWriPtr%PkgObsWriter
      if (associated(obsWriter)) then
        if (obsWriter%Active) then
          call obsWriter%WriteObsFile(this%model%IGrid)
        endif
      endif
    enddo
    !
    ! Write CHD observations file if needed.
    if (associated(this%model%IbChdWriter%PkgObsWriter)) then
      if (this%model%IbChdWriter%PkgObsWriter%Active) then
        call this%model%IbChdWriter%PkgObsWriter%WriteObsFile(this%model%IGrid)
      endif
    endif
    !
    ! Write CHD file if needed.
    if (this%model%IbChdWriter%Active) then
      call this%model%IbChdWriter%WriteFile(this%model%IGrid)
    endif
    !
    ! If IBOUND CHD OBS is inactive, delete the input file that was started.
    if (.not. this%model%IbChdWriter%Active) then
      fname = this%model%IbChdWriter%PkgObsWriter%fileobj%FName
      if (fname /= ' ') then
        inquire(file=fname,exist=lex,opened=lop,number=iu)
        if (lex) then
          if (lop) then
            close(iu,status='DELETE')
          else
            iu = getunit()
            open(iu,file=fname)
            close(iu,status='DELETE')
          endif
        endif
      endif
    endif
    !
    if (.not. this%ConversionDone) then
      ! Write name file for this model
      call this%model%WriteNameFile()
    endif
    !
    ! Write package conversion table, first to IOUT then to the screen.
    call this%model%PackageConversionTable(iout)
    call this%model%PackageConversionTable(0)
    !
    ! Close all MF6 files.
    call this%model%Mf6Files%CloseAll()
    !
    write(msg,20)trim(this%model%BaseName) // '.nam'
    !
    this%ConversionDone = .true.
    this%model%ConversionDone = .true.
    !
    return
  end subroutine ConvertModel

  subroutine WritePostObsScript(this)
    implicit none
    ! dummy
    class(ModelConverterType) :: this
    ! local
    integer :: iu, numobsml
    character(len=MAXCHARLEN) :: ermsg, msg
    ! formats
    1 format()
    10 format(a)
    20 format(a,1x,a)
    !
    if (PathToPostObsMf /= '') then
      select case (ScriptType)
      case ('BATCH')
        this%PostObsScriptName = trim(this%model%ModelName) // '_PostObs.bat'
        iu = GetUnit()
        call openfile(iu,0,this%PostObsScriptName,'SCRIPT',filstat_opt='REPLACE')
        write(iu,10)'@echo off'
        write(iu,1)
        write(iu,10)'@echo.'
        write(iu,10)'REM Invoke PostObsMF to post-process single-layer head ' // &
                    'observations'
        msg = 'Batch file ' // trim(this%PostObsScriptName) // &
              ' is running PostObsMF for ' // &
              'single-layer head observations of model "' // &
              trim(this%model%ModelName) // '".'
        write(iu,20)'echo', trim(msg)
        write(iu,20)trim(PathToPostObsMf), trim(this%model%HedObsWriter%Preproc%PostObsFilename)
        numobsml = this%model%HedObsWriter%MLObsList%Count()
        if (numobsml > 0) then
          write(iu,1)
          write(iu,10)'@echo.'
          write(iu,10)'REM Invoke PostObsMF to post-process multi-layer head observations'
          msg = 'Batch file ' // trim(this%PostObsScriptName) // &
                ' is running PostObsMF for ' // &
                'multi-layer head observations of model "' // &
                trim(this%model%ModelName) // '".'
          write(iu,20)'echo', trim(msg)
          write(iu,20)trim(PathToPostObsMf), trim(this%model%HedObsWriter%MlPostObsFileName)
        endif
        msg = 'Batch file ' // trim(this%PostObsScriptName) // &
              ' can be used to post-process head' // &
              ' observations for model "' // trim(this%model%ModelName) // &
              '" after the model is run.'
        call store_note(msg)
      case ('PYTHON')
        ermsg = 'SCRIPT type PYTHON is not yet supported.'
        call store_error(ermsg)
        call ustop()
  !      this%PostObsScriptName = trim(this%model%ModelName) // '_PostObs.py'
      end select
    else
      msg = 'Mf5to6 can be setup to generate a BATCH (or PYTHON, [not' // &
            ' yet  supported]) script to automate post-processing of' // &
            ' head observations by ' // &
            'providing an MF5to6 options file -- see documentation of MF5to6.'
      call store_note(msg)
    endif
    !
    !
    return
  end subroutine WritePostObsScript

  subroutine ValidateMovers(this, scope)
    implicit none
    ! dummy
    class(ModelConverterType) :: this
    character(len=*), intent(in) :: scope
    ! local
    integer :: igridProv, igridRec, imov, nmov
    character(len=LENMODELNAME) :: modelname1, modelname2
    character(len=LENPACKAGENAME) :: pkgname1, pkgname2
    character(len=MAXCHARLEN) :: ermsg
    type(MoverType), pointer :: mover
    type(ListType), pointer :: moverlist => null()
    !
    select case (scope)
    case ('MODEL')
      moverlist => this%model%ModelMovers
    case ('SIMULATION')
      moverlist => SimMovers
    case default
      ermsg = 'Programmer error in ValidateMovers. Invalid scope: ' // &
               trim(scope)
      call store_error(ermsg)
      call ustop()
    end select
    !
    nmov = moverlist%Count()
    do imov=1,nmov
      mover => GetMoverFromList(moverlist, imov)
      igridProv = mover%IgridProvider
      if (mover%ProvPkgName == '') then
        call this%LookUpPkgName(igridProv, mover%ProvPkgType, pkgname1)
        mover%ProvPkgName = pkgname1
      endif
      if (mover%ProvModelName == '') then
        call this%LookUpModelName(igridProv, modelname1)
        mover%ProvModelName = modelname1
      endif
      igridRec = mover%IgridReceiver
      if (mover%RecPkgName == '') then
        call this%LookUpPkgName(igridRec, mover%RecPkgType, pkgname2)
        mover%RecPkgName = pkgname2
      endif
      if (mover%RecModelName == '') then
        call this%LookUpModelName(igridRec, modelname2)
        mover%RecModelName = modelname2
      endif
    enddo
    !
    return
  end subroutine ValidateMovers

  subroutine LookUpModelName(this, igrid, modelName)
    ! Look up model name, given grid number.
    implicit none
    ! dummy
    class(ModelConverterType) :: this
    integer, intent(in) :: igrid
    character(len=*), intent(out) :: modelName
    ! local
    integer :: i, nchild
    type(ModelType), pointer :: model
    type(ModelConverterType), pointer :: modelConverter => null()
    class(*), pointer :: obj => null()
    !
    model => null()
    modelName = ''
    !
    ! First, find model that has Igrid of interest
    if (this%model%IGrid == igrid) then
      model => this%model
    endif
    !
    if (.not. associated(model)) then
      if (associated(this%ParentConverter)) then
        if (this%ParentConverter%model%IGrid == igrid) then
          model => this%ParentConverter%model
        endif
      endif
    endif
    !
    if (.not. associated(model)) then
      nchild = this%ChildConverters%Count()
      childloop: do i=1,nchild
        obj => this%ChildConverters%GetItem(i)
        select type (obj)
        type is (ModelConverterType)
          modelConverter => obj
          if (modelConverter%model%IGrid == igrid) then
            model => modelConverter%model
            exit childloop
          endif
        end select
      enddo childloop
    endif
    !
    ! Get the model name.
    if (associated(model)) then
      modelName = model%ModelName
    endif
    !
    return
  end subroutine LookUpModelName

  subroutine LookUpPkgName(this, igrid, pkgType, pkgName)
    ! Look up package name, given package type and grid number.
    implicit none
    ! dummy
    class(ModelConverterType) :: this
    integer, intent(in) :: igrid
    character(len=*), intent(in)  :: pkgType
    character(len=*), intent(out) :: pkgName
    ! local
    integer :: i, nchild, npkg
    class(PackageWriterType), pointer :: pkgWriPtr => null()
    type(ModelType), pointer :: model
    type(ModelConverterType), pointer :: modelConverter => null()
    class(*), pointer :: obj => null()
    !
    pkgName = ''
    model => null()
    !
    ! First, find model that has Igrid of interest
    if (this%model%IGrid == igrid) then
      model => this%model
    endif
    !
    if (.not. associated(model)) then
      if (associated(this%ParentConverter)) then
        if (this%ParentConverter%model%IGrid == igrid) then
          model => this%ParentConverter%model
        endif
      endif
    endif
    !
    if (.not. associated(model)) then
      nchild = this%ChildConverters%Count()
      childloop: do i=1,nchild
        obj => this%ChildConverters%GetItem(i)
        select type (obj)
        type is (ModelConverterType)
          modelConverter => obj
          if (modelConverter%model%IGrid == igrid) then
            model => modelConverter%model
            exit childloop
          endif
        end select
      enddo childloop
    endif
    !
    ! Iterate through package writers in model and find
    ! matching package type, then return the package name.
    if (associated(model)) then
      npkg = model%PackageWriters%Count()
      pkgloop: do i=1,npkg
        pkgWriPtr => model%GetPackageWriter(i)
        if (pkgWriPtr%PkgType == pkgType) then
          pkgName = pkgWriPtr%PackageName
          exit pkgloop
        endif
      enddo pkgloop
    endif
    !
    return
  end subroutine LookUpPkgName

  subroutine AddMover(this, mover)
    implicit none
    ! dummy
    class(ModelConverterType) :: this
    type(MoverType), pointer :: mover
    !
    if (mover%ProvModelName == mover%RecModelName) then
      call AddMoverToList(this%ModelMovers, mover)
    else
      call AddMoverToList(SimMovers, mover)
    endif
    !
    return
  end subroutine AddMover

  function GetMover(this, idx) result (res)
    implicit none
    ! dummy
    class(ModelConverterType) :: this
    integer, intent(in) :: idx
    type(MoverType), pointer :: res
    !
    res => GetMoverFromList(this%ModelMovers, idx)
    !
    return
  end function GetMover

end module ModelConverterModule
