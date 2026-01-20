module ModelModule

  use ChdModule, only: ChdType
  use ChdPackageWriterModule, only: ChdPackageWriterType
  use ConstantsModule, only: LENMODELNAME, LENPACKAGENAME, MAXCHARLEN, DZERO
  use ConstantsPHMFModule, only: FCDATABIN, FCDATABOUT, FCDATAIN, FCDATAOUT, &
                             FCINPUT, FCOUTPUT, FCUNKNOWN, LENCTYPE
  use DisWriterModule, only: DisWriterType
  use DrnPackageWriterModule, only: DrnPackageWriterType
  use EvtPackageWriterModule, only: EvtPackageWriterType
  use FhbPackageWriterModule, only: FhbPackageWriterType
  use FileListModule, only: FileListType
  use FileTypeModule, only: FileType
  use FileWriterModule, only: FileWriterType
  use GhbPackageWriterModule, only: GhbPackageWriterType
  use GLOBAL, only: GLOBALTYPE, LAYCBD, NLAY, AllocateGlobalScalars, &
                    NCOL, NROW, IBOUND, BOTM, LBOTM, NBOTM, DELC, DELR
  use GlobalVariablesModule, only: GetNextIgrid
  use GwfBasModule, only: SGWF2BAS7PNT, SGWF2BAS7PSV
  use GWFBCFMODULE, only: SGWF2BCF7PNT, SGWF2BCF7PSV, LAYCON, WetdryBcf => WETDRY
  use GWFCHDMODULE, only: SGWF2CHD7PNT, SGWF2CHD7PSV
  use GwfDrnModule, only: SGWF2DRN7PNT, SGWF2DRN7PSV
  use GwfEtsModule, only: SGWF2ETS7PNT, SGWF2ETS7PSV
  use GwfEvtModule, only: SGWF2EVT7PNT, SGWF2EVT7PSV
  use GWFFHBMODULE, only: SGWF2FHB7PNT, SGWF2FHB7PSV
  use GwfGhbModule, only: SGWF2GHB7PNT, SGWF2GHB7PSV
  use GWFLAKMODULE, only: SGWF2LAK7PNT, SGWF2LAK7PSV
  use GWFLPFMODULE, only: SGWF2LPF7PNT, GWF2LPF7PSV, VKCB, LAYTYP, LAYWET, WetdryLpf => WETDRY
  use GWFMNW2MODULE, only: SGWF2MNW2PNT, SGWF2MNW2PSV
  use GwfRchModule, only: SGWF2RCH7PNT, SGWF2RCH7PSV
  use GwfResModule, only: SGWF2RES7PNT, SGWF2RES7PSV
  use GwfRivModule, only: SGWF2RIV7PNT, SGWF2RIV7PSV
  use GWFSFRMODULE, only: SGWF2SFR7PNT, SGWF2SFR7PSV
  use GWFUPWMODULE, only: LaytypUpw, LaywetUpw => LAYWET, VkcbUpw => VKCB, &
                          WetdryUpw => WETDRY, HaniUpw => HANI
  use GwfWelNwt, only: SGWF2WEL7PNT, SGWF2WEL7PSV
  use HfbPackageWriterModule, only: HfbPackageWriterType
  use HfBSubsNwt, only: SGWF2HFB7PNT, SGWF2HFB7PSV
  use IcWriterModule, only: IcWriterType
  use InputOutputModule, only: GetUnit, openfile
  use LakPackageWriterModule, only: LakPackageWriterType
  use LGRMODULE, only: IBFLG, SGWF2LGR2PNT
  use ListModule, only: ListType
  use MawPackageWriterModule, only: MawPackageWriterType
  use MvrPackageWriterModule, only: MvrPackageWriterType
  use NpfWriterModule, only: NpfWriterType
  use ObsWriterModule, only: ObsWriterType
  use OutputControlWriterModule, only: OutputControlWriterType
  use PackageWriterModule, only: PackageWriterType, CastAsPackageWriterClass
  use RchPackageWriterModule, only: RchPackageWriterType
  use RivPackageWriterModule, only: RivPackageWriterType
  use SfrPackageWriterModule, only: SfrPackageWriterType, &
                                    AddSfrPackageWriter
  use SimPHMFModule, only: store_error, store_note, ustop
  use ImsPackageWriterModule, only: ImsPackageWriterType
  use StoWriterModule, only: StoWriterType
  use TdisVariablesModule, only: GlobalTdisWriter
  use TdisWriterModule, only: TdisWriterType
  use UzfPackageWriterModule, only: UzfPackageWriterType
  use WelPackageWriterModule, only: WelPackageWriterType

  implicit none

  private
  public :: ModelType

  type :: ModelType
    character(len=MAXCHARLEN)   :: NameFile2005
    character(len=MAXCHARLEN)   :: NameFile2015
    character(len=MAXCHARLEN)   :: BaseName
    character(len=LENMODELNAME) :: ModelName
    character(len=MAXCHARLEN)   :: cbcfilename
    character(len=12),    private  :: LengthUnit = 'UNDEFINED'
    character(len=12),    private  :: TimeUnit = 'UNDEFINED'
    integer, pointer,     public   :: Nlaynew => null() ! = nlay + num quasi-3d units
    integer, pointer,     public   :: Ncol => null(), Nrow => null()
    integer, pointer, dimension(:) :: Layptr => null()  ! dim = nlay
    integer, pointer, dimension(:,:,:) :: Idomain => null()
    ! Grid number for this model
    integer :: IGrid = 0
    ! LGR interface flag for this model
    integer :: Ibflg = 0
    ! Unit numbers associated with observation input for individual packages
    integer :: IuChObs = 0
    integer :: IuDrnObs = 0
    integer :: IuDrtObs = 0
    integer :: IuGhbObs = 0
    integer :: IuHeadObs = 0
    integer :: IuRivObs = 0
    integer :: IuSfrObs = 0
    integer :: IuStrObs = 0
    integer :: iulist = 0
    logical :: ChdFileAdded = .false.
    logical :: ConversionDone = .false.
    logical :: Newton = .false.
    logical :: NewtonUnderRelaxation = .false.
    integer, pointer :: nbotm => null()
    double precision, pointer, dimension(:) :: delr => null()
    double precision, pointer, dimension(:) :: delc => null()
    double precision, pointer, dimension(:,:,:) :: botm => null()
    integer, pointer, dimension(:) :: lbotm => null()
    type(FileListType),            pointer :: Mf2005Files => null()
    type(FileListType),            pointer :: Mf6Files => null()
    type(DisWriterType),           pointer :: DisWriter => null()
    type(TdisWriterType),          pointer :: TdisWriter => null()
    type(IcWriterType),            pointer :: IcWriter => null()
    type(ChdPackageWriterType),    pointer :: IbChdWriter => null()
    type(NpfWriterType),           pointer :: NpfWriter => null()
    type(ObsWriterType),           pointer :: HedObsWriter => null()
    type(StoWriterType),           pointer :: StoWriter => null()
    type(FhbPackageWriterType),    pointer :: FhbWriter => null()
    type(ImsPackageWriterType),    pointer :: ImsWriter => null()
    type(OutputControlWriterType), pointer :: OcWriter => null()
    type(ListType),                pointer :: PackageWriters => null()
    type(MvrPackageWriterType),    pointer :: MvrWriter => null()
    type(ListType),                pointer :: ModelMovers => null()
  contains
    procedure :: AddPackage
    procedure :: AddPackageWriter
    procedure :: AssignIgrid
    procedure :: AssignLengthUnit
    procedure :: AssignTimeUnit
    procedure :: BuildAllModelMovers
    procedure :: ContainsPackageWriterOfType
    procedure :: GetLengthUnit
    procedure :: GetPackageWriter
    procedure :: GetTimeUnit
    procedure :: Initialize
    procedure :: InitializeMF6File
    procedure :: InitializeMF6Files
    procedure :: InitializeIdomain
    procedure :: InitializeLayptr
    procedure :: OrderPackageWriters
    procedure :: PackageConversionTable
    procedure :: PointToGrid
    procedure :: RemoveFhbChdDuplicates
    procedure :: SavePointers
    procedure :: WriteNameFile
    procedure :: WriteMvrFile
  end type ModelType

contains

  subroutine Initialize(this, namefile, basename)
    implicit none
    ! dummy
    class(ModelType), intent(inout) :: this
    character(len=*), intent(in) :: namefile, basename
    ! local
    character(len=MAXCHARLEN) :: fname
    !
    if (.not. this%ConversionDone) then
      allocate(this%ncol)
      allocate(this%nrow)
      allocate(this%nbotm)
      call AllocateGlobalScalars()
      call this%AssignIgrid()
      !
      ! Allocate and initialize Mf2005Files
      allocate(this%Mf2005Files)
      call this%Mf2005Files%Initialize()
      this%Mf2005Files%files%name = 'Mf2005Files'
      !
      ! Allocate and initialize Mf6Files
      allocate(this%Mf6Files)
      call this%Mf6Files%Initialize()
      this%Mf6Files%files%name = 'Mf6Files'
      !
      ! Allocate pointers
      allocate(this%DisWriter)
      allocate(this%TdisWriter)
      if (.not. associated(GlobalTdisWriter)) then
        GlobalTdisWriter => this%TdisWriter
      endif
      allocate(this%IcWriter)
      !
      allocate(this%IbChdWriter)
      this%IbChdWriter%ModelBasename = trim(basename) // '.ib'
      call this%IbChdWriter%AllocatePointers()
!      this%IbChdWriter%PackageName = 'CHD'
      this%IbChdWriter%source = 'IBOUND'
      this%IbChdWriter%PkgObsWriter%Source = 'IBOUND'
      !allocate(this%IbChdWriter%fileobj)
      call this%IbChdWriter%fileobj%Initialize()
      this%IbChdWriter%fileobj%PkgName = this%IbChdWriter%PackageName
      !
      allocate(this%FhbWriter)
      !this%FhbWriter%ModelBasename = basename
!      this%FhbWriter%PackageName = 'FHB'
      allocate(this%FhbWriter%fileobj)
      this%FhbWriter%ModelBaseName = basename
      call this%FhbWriter%fileobj%Initialize()
      this%FhbWriter%fileobj%PkgName = this%FhbWriter%PackageName
      allocate(this%FhbWriter%ChdWriter)
      allocate(this%FhbWriter%ChdWriter%IbChdList)
      allocate(this%FhbWriter%ChdWriter%TvChdList)
      this%FhbWriter%ChdWriter%ModelBasename = basename
      this%FhbWriter%ChdWriter%PackageName = 'CHD-FHB'
      allocate(this%FhbWriter%ChdWriter%fileobj)
      call this%FhbWriter%ChdWriter%fileobj%Initialize()
      this%FhbWriter%ChdWriter%fileobj%PkgName = this%FhbWriter%ChdWriter%PackageName
      allocate(this%FhbWriter%WellWriter)
      this%FhbWriter%WellWriter%ModelBasename = basename
      this%FhbWriter%WellWriter%PackageName = 'WEL-FHB'
      allocate(this%FhbWriter%WellWriter%fileobj)
      call this%FhbWriter%WellWriter%fileobj%Initialize()
      this%FhbWriter%WellWriter%fileobj%PkgName = this%FhbWriter%WellWriter%PackageName
      !
      allocate(this%OcWriter)
      allocate(this%NpfWriter)
      allocate(this%HedObsWriter)
      allocate(this%StoWriter)
      allocate(this%ImsWriter)
      this%ImsWriter%ModelBasename = basename
      call this%ImsWriter%AllocatePointers()
      allocate(this%PackageWriters)
      !
      ! Water mover stuff
      allocate(this%MvrWriter)
      this%MvrWriter%WriteModelNames = .false.
      this%MvrWriter%ModelBasename = this%BaseName
      allocate(this%ModelMovers)
      this%ModelMovers%name = 'ModelMovers'
      call this%MvrWriter%InitializeMvrWriter(this%ModelMovers)
      fname = trim(this%BaseName) // '.mvr6'
      call this%MvrWriter%InitializeFile(fname, 'MVR6', 'MVR')
      !
      ! Assign member pointers
      this%NpfWriter%Mf2005Files => this%Mf2005Files
      this%NpfWriter%Mf6Files => this%Mf6Files
      !
      ! Assign members
      this%cbcfilename = ''
      this%NameFile2005 = namefile
      this%BaseName = basename
      this%NameFile2015 = trim(basename) // '.nam'
      this%ModelName = basename
    endif
    !
    return
  end subroutine Initialize

  subroutine InitializeMF6Files(this)
    implicit none
    ! dummy
    class(ModelType), intent(inout) :: this
    !
    ! Initialize TdisWriter
    call this%InitializeMF6File('TDIS6')
    !
    ! Initialize file writers
    call this%InitializeMF6File('DIS6')
    call this%InitializeMF6File('IC6')
    call this%InitializeMF6File('LIST')
    call this%InitializeMF6File('NPF6')
    call this%InitializeMF6File('STO6')
    !
    ! Write Options blocks for file types that have no options
    call this%IcWriter%WriteOptions()
    !
    return
  end subroutine InitializeMF6Files

  subroutine InitializeMF6File(this, ftype, FCode)
    implicit none
    ! dummy
    class(ModelType), intent(inout) :: this
    character(len=*), intent(in) :: ftype
    integer, intent(in), optional :: FCode
    ! local
    class(FileWriterType), pointer :: fwt
    character(len=MAXCHARLEN) :: fname, msg
    character(len=LENPACKAGENAME) :: pkgname
    character(len=12) :: ftypelocal
    integer :: iu, FCodeLocal
    !
    FCodeLocal = -1
    if (present(FCode)) FCodeLocal = FCode
    !
    fwt => null()
    pkgname = ''
    !
    ftypelocal = ftype
    select case(ftype)
    !
    ! INPUT files
    case ('CHD6')
      ftypelocal = 'CHD6'
      fwt => this%IbChdWriter
      fwt%PkgType = 'CHD'
      fname = trim(this%BaseName) // '.chd'
      FCodeLocal = FCUNKNOWN
      pkgname = this%IbChdWriter%PackageName
      ! CHD (constant) will be activated only if needed
      if (fwt%Active) FCodeLocal = FCINPUT
    case ('DIS6')
      if (this%ConversionDone) return
      fwt => this%DisWriter
      fname = trim(this%BaseName) // '.dis'
      FCodeLocal = FCINPUT
      if (associated(this%DisWriter%fileobj)) then
        pkgname = this%DisWriter%fileobj%PkgName
      endif
      ! DIS is always active
      fwt%Active = .true.
    case ('IC6')
      if (this%ConversionDone) return
      fwt => this%IcWriter
      fname = trim(this%BaseName) // '.ic'
      FCodeLocal = FCINPUT
      if (associated(this%IcWriter%fileobj)) then
        pkgname = this%IcWriter%fileobj%PkgName
      endif
      ! IC is always active
      fwt%Active = .true.
    case ('NPF6')
      if (this%ConversionDone) return
      fwt => this%NpfWriter
      fname = trim(this%BaseName) // '.npf'
      FCodeLocal = FCINPUT
      if (associated(this%NpfWriter%fileobj)) then
        pkgname = this%NpfWriter%fileobj%PkgName
      endif
      ! NPF is always active
      fwt%Active = .true.
    case ('STO6')
      if (this%ConversionDone) return
      fwt => this%StoWriter
      fname = trim(this%BaseName) // '.sto'
      FCodeLocal = FCUNKNOWN
      if (associated(this%StoWriter%fileobj)) then
        pkgname = this%StoWriter%fileobj%PkgName
      endif
      ! STO will be activated only if needed
      if (fwt%Active) FCodeLocal = FCINPUT
    case ('TDIS6')
      if (this%ConversionDone) return
      fwt => this%TdisWriter
      fname = trim(this%BaseName) // '.tdis'
      FCodeLocal = FCINPUT
      if (associated(this%TdisWriter%fileobj)) then
        pkgname = this%TdisWriter%fileobj%PkgName
      endif
      ! TDIS is always active
      fwt%Active = .true.
    case ('OC6')
      if (this%ConversionDone) return
      fwt => this%OcWriter
      fname = trim(this%BaseName) // '.oc'
      FCodeLocal = FCUNKNOWN
      ! OC will be activated only if needed
      if (fwt%Active) FCodeLocal = FCINPUT
    !
    ! OUTPUT files
    case ('LIST')
      if (this%ConversionDone) return
      fname = trim(this%BaseName) // '.lst'
      iu = 0
      FCodeLocal = FCOUTPUT
    !
    case ('DATA')
      ! Ned todo: add support for DATA file type; need to determine if input or output
      msg = 'need to add support for DATA file type'
      call store_error(msg)
      call ustop()
    case ('DATA(BINARY)')
      ! Ned todo: add support for DATA(BINARY) file type; need to determine if input or output
      msg = 'need to add support for DATA(BINARY) file type'
      call store_error(msg)
      call ustop()
    case default
      msg = 'Unrecognized file type in ModelType.InitializeFile: ' // trim(ftype)
      call store_error(msg)
      call ustop()
    end select
    !
    ! Initialize MF6 input files being prepared by MF5to6
    if (associated(fwt)) then
      call fwt%InitializeFile(fname, ftypelocal)
      iu = fwt%fileobj%IUnit
    endif
    !
    ! Add input and output files for active file writers to list
    ! of files that will be written to name file.
    if (iu>=0 .and. FCodeLocal/=FCUNKNOWN) then
      call this%Mf6Files%AddFile(fname, ftypelocal, iu, FCodeLocal, &
                                    pkgname)
    endif
    !
    return
  end subroutine InitializeMF6File

  subroutine InitializeLayptr(this)
    ! Allocate and populate Layptr, which, for each layer of the
    ! original model, will contain the layer index in the new model.
    implicit none
    ! dummy
    class(ModelType) :: this
    ! local
    integer :: i, ilay, npak
    class(PackageWriterType), pointer :: pkgwriter
    !
    allocate(this%Nlaynew)
    allocate(this%Layptr(NLAY))
    ilay = 0
    do i = 1,nlay
      ilay = ilay + 1
      this%Layptr(i) = ilay
      if (LAYCBD(i) /= 0) then
        ! layer has an underlying quasi-3d confining unit
        ilay = ilay + 1
      endif
    enddo
    this%Nlaynew = ilay
    !
    ! -- Assign pointers needed for conversion of quasi-3d units
    this%DisWriter%Nlaynew => this%Nlaynew
    this%DisWriter%Layptr => this%Layptr
    !
    this%NpfWriter%Nlaynew => this%Nlaynew
    this%NpfWriter%Layptr => this%Layptr
    !
    this%StoWriter%Nlaynew => this%Nlaynew
    this%StoWriter%Layptr => this%Layptr
    !
    this%IbChdWriter%Nlaynew => this%Nlaynew
    this%IbChdWriter%Layptr => this%Layptr
    !
    npak = this%PackageWriters%Count()
    do i=1,npak
      pkgwriter => this%GetPackageWriter(i)
      pkgwriter%Nlaynew => this%Nlaynew
      pkgwriter%Layptr => this%Layptr
    enddo
    !
    return
  end subroutine InitializeLayptr

  subroutine WriteNameFile(this)
    implicit none
    ! dummy
    class(ModelType), intent(inout) :: this
    ! local
    character(len=MAXCHARLEN) :: fname, listfname
    integer :: iu, i, n
    type(FileType), pointer :: fil => null()
    character(len=LENPACKAGENAME+2) :: pname
    ! formats
    10 format('BEGIN ',a)
    15 format(2x,a,t14,a)
    20 format(2x,a,t14,a,2x,a)
    50 format('END ',a)
    60 format('')
    !
    if (this%ConversionDone) return
    !
    ! Define name for name file and open it
    fname = trim(this%BaseName) // '.nam'
    iu = GetUnit()
    open(iu, file=fname, status='REPLACE')
    !
    ! Write Options block containing name of LIST file
    fil => this%Mf6Files%GetFile('LIST')
    listfname = fil%FName
    write(iu,10)'Options'
    write(iu,15)'LIST',trim(listfname)
    if (this%Newton) then
      if (this%NewtonUnderRelaxation) then
        write(iu,15)'NEWTON UNDER_RELAXATION'
      else
        write(iu,15)'NEWTON'
      endif
    endif
    write(iu,50)'Options'
    write(iu,60)
    !
    ! Write Packages block
    write(iu,10)'Packages'
    n = this%Mf6Files%NCount()
    do i=1,n
      fil => this%Mf6Files%GetFile(i)
      if (fil%FCode==FCINPUT .or. fil%FCode==FCDATAIN .or. &
          fil%FCode==FCDATABIN) then
        if (fil%ftype/='TDIS6') then
          ! TDIS6 file needs to be listed in simulation name file
          if (fil%PkgName /= '') then
            pname = '''' // trim(fil%PkgName) // ''''
            write(iu,20)trim(fil%ftype), trim(fil%fname), pname
          else
            write(iu,15)trim(fil%ftype), trim(fil%fname)
          endif
        endif
      endif
    enddo
    write(iu,50)'Packages'
    write(iu,60)
    !
    return
  end subroutine WriteNameFile

  subroutine AssignTimeUnit(this, TimeUnit)
    implicit none
    class(ModelType), intent(inout) :: this
    character(len=*), intent(in) :: TimeUnit
    !
    this%TimeUnit = TimeUnit
    call this%TdisWriter%AssignTimeUnit(TimeUnit)
    !
    return
  end subroutine AssignTimeUnit

  subroutine AssignLengthUnit(this, LengthUnit)
    implicit none
    class(ModelType), intent(inout) :: this
    character(len=*), intent(in) :: LengthUnit
    !
    this%LengthUnit = LengthUnit
    call this%DisWriter%AssignLengthUnit(LengthUnit)
    !
    return
  end subroutine AssignLengthUnit

  function GetTimeUnit(this) result(tu)
    implicit none
    class(ModelType), intent(in) :: this
    character(len=12) :: tu
    tu = this%TimeUnit
    return
  end function GetTimeUnit

  function GetLengthUnit(this) result(lu)
    implicit none
    class(ModelType), intent(in) :: this
    character(len=12) :: lu
    lu = this%LengthUnit
    return
  end function GetLengthUnit

  function AddPackage(this, PType) result(PkgWriter)
    implicit none
    ! dummy
    class(ModelType) :: this
    character(len=*), intent(in) :: PType
    class(PackageWriterType), pointer :: PkgWriter, pw
    ! local
    type(ChdPackageWriterType), pointer :: ChdWriter => null()
    type(DrnPackageWriterType), pointer :: DrnWriter => null()
    type(EvtPackageWriterType), pointer :: EvtWriter => null()
    type(GhbPackageWriterType), pointer :: GhbWriter => null()
    type(HfbPackageWriterType), pointer :: HfbWriter => null()
    type(LakPackageWriterType), pointer :: LakWriter => null()
    type(MawPackageWriterType), pointer :: MawWriter => null()
    type(RchPackageWriterType), pointer :: RchWriter => null()
    type(RivPackageWriterType), pointer :: RivWriter => null()
    type(SfrPackageWriterType), pointer :: SfrWriter => null()
    type(UzfPackageWriterType), pointer :: UzfWriter => null()
    type(WelPackageWriterType), pointer :: WelWriter => null()
    character(len=100) :: msg
    ! format
    10 format(a,i0)
    !
    PkgWriter => null()
    select case (PType)
    case ('CHD')
      pw => ChdWriter
      if (.not. this%ContainsPackageWriterOfType(pw)) then
        allocate(ChdWriter)
        ChdWriter%ModelBasename = this%BaseName
        call ChdWriter%AllocatePointers()
        PkgWriter => ChdWriter
        PkgWriter%ModelMovers => this%ModelMovers
        PkgWriter%source = 'CHD'
        !PkgWriter%ModelBasename = this%BaseName
        !allocate(PkgWriter%fileobj)
        !call PkgWriter%fileobj%Initialize()
        PkgWriter%PkgType = 'CHD'
        call this%AddPackageWriter(PkgWriter)
      endif
    case ('DRN')
      pw => DrnWriter
      if (.not. this%ContainsPackageWriterOfType(pw)) then
        allocate(DrnWriter)
        PkgWriter => DrnWriter
        PkgWriter%ModelMovers => this%ModelMovers
        PkgWriter%source = 'DRN'
        PkgWriter%ModelBasename = this%BaseName
        allocate(PkgWriter%fileobj)
        call PkgWriter%fileobj%Initialize()
        PkgWriter%PkgType = 'DRN'
        call this%AddPackageWriter(PkgWriter)
      endif
    case ('EVT')
      pw => EvtWriter
      if (.not. this%ContainsPackageWriterOfType(pw)) then
        allocate(EvtWriter)
        PkgWriter => EvtWriter
        PkgWriter%ModelMovers => this%ModelMovers
        PkgWriter%source = 'EVT'
        PkgWriter%ModelBasename = this%BaseName
        allocate(PkgWriter%fileobj)
        call PkgWriter%fileobj%Initialize()
        PkgWriter%PkgType = 'EVT'
        call this%AddPackageWriter(PkgWriter)
      endif
    case ('ETS')
      ! Convert ETS to EVT8
      pw => EvtWriter
      if (.not. this%ContainsPackageWriterOfType(pw)) then
        allocate(EvtWriter)
        PkgWriter => EvtWriter
        PkgWriter%ModelMovers => this%ModelMovers
        PkgWriter%source = 'ETS'
        PkgWriter%ModelBasename = this%BaseName
        allocate(PkgWriter%fileobj)
        call PkgWriter%fileobj%Initialize()
        PkgWriter%PkgType = 'ETS'
        call this%AddPackageWriter(PkgWriter)
      endif
    case ('GHB')
      pw => GhbWriter
      if (.not. this%ContainsPackageWriterOfType(pw)) then
        allocate(GhbWriter)
        PkgWriter => GhbWriter
        PkgWriter%ModelMovers => this%ModelMovers
        PkgWriter%source = 'GHB'
        PkgWriter%ModelBasename = this%BaseName
        allocate(PkgWriter%fileobj)
        call PkgWriter%fileobj%Initialize()
        PkgWriter%PkgType = 'GHB'
        call this%AddPackageWriter(PkgWriter)
      endif
    case ('LAK')
      pw => LakWriter
      if (.not. this%ContainsPackageWriterOfType(pw)) then
        allocate(LakWriter)
        LakWriter%Igrid = this%IGrid
        PkgWriter => LakWriter
        PkgWriter%ModelMovers => this%ModelMovers
        PkgWriter%source = 'LAK'
        PkgWriter%ModelBasename = this%BaseName
        allocate(PkgWriter%fileobj)
        call PkgWriter%fileobj%Initialize()
        PkgWriter%PkgType = 'LAK'
        call this%AddPackageWriter(PkgWriter)
      endif
    case ('MAW')
      pw => MawWriter
      if (.not. this%ContainsPackageWriterOfType(pw)) then
        allocate(MawWriter)
        MawWriter%Mf6Files => this%Mf6Files
        MawWriter%source = 'MAW'
        MawWriter%ModelBasename = this%BaseName
        allocate(MawWriter%fileobj)
        call MawWriter%fileobj%Initialize()
        MawWriter%PkgType = 'MAW'
        MawWriter%Igrid = this%IGrid
        write(MawWriter%PackageName,10)'MAW_',this%IGrid
        PkgWriter => MawWriter
        PkgWriter%ModelMovers => this%ModelMovers
        call this%AddPackageWriter(PkgWriter)
      endif
    case ('RCH')
      pw => RchWriter
      if (.not. this%ContainsPackageWriterOfType(pw)) then
        allocate(RchWriter)
        PkgWriter => RchWriter
        PkgWriter%ModelMovers => this%ModelMovers
        PkgWriter%source = 'RCH'
        PkgWriter%ModelBasename = this%BaseName
        allocate(PkgWriter%fileobj)
        call PkgWriter%fileobj%Initialize()
        PkgWriter%PkgType = 'RCH'
        call this%AddPackageWriter(PkgWriter)
      endif
    case ('RIV')
      pw => RivWriter
      if (.not. this%ContainsPackageWriterOfType(pw)) then
        allocate(RivWriter)
        PkgWriter => RivWriter
        PkgWriter%ModelMovers => this%ModelMovers
        PkgWriter%source = 'RIV'
        PkgWriter%ModelBasename = this%BaseName
        allocate(PkgWriter%fileobj)
        call PkgWriter%fileobj%Initialize()
        PkgWriter%PkgType = 'RIV'
        call this%AddPackageWriter(PkgWriter)
      endif
    case ('SFR')
      pw => SfrWriter
      if (.not. this%ContainsPackageWriterOfType(pw)) then
        allocate(SfrWriter)
        SfrWriter%Mf6Files => this%Mf6Files
        SfrWriter%source = 'SFR'
        SfrWriter%ModelBasename = this%BaseName
        allocate(SfrWriter%fileobj)
        call SfrWriter%fileobj%Initialize()
        SfrWriter%PkgType = 'SFR'
        SfrWriter%Igrid = this%IGrid
        write(SfrWriter%PackageName,10)'SFR_',this%IGrid
        PkgWriter => SfrWriter
        PkgWriter%ModelMovers => this%ModelMovers
        call this%AddPackageWriter(PkgWriter)
        ! Add the SfrWriter to global list of all SfrPackageWriters
        call AddSfrPackageWriter(SfrWriter)
      endif
    case ('UZF')
      pw => UzfWriter
      if (.not. this%ContainsPackageWriterOfType(pw)) then
        allocate(UzfWriter)
        UzfWriter%Mf6Files => this%Mf6Files
        UzfWriter%source = 'UZF'
        UzfWriter%ModelBasename = this%BaseName
        allocate(UzfWriter%fileobj)
        call UzfWriter%fileobj%Initialize()
        UzfWriter%PkgType = 'UZF'
        UzfWriter%Igrid = this%IGrid
        write(UzfWriter%PackageName,10)'UZF_',this%IGrid
        PkgWriter => UzfWriter
        PkgWriter%ModelMovers => this%ModelMovers
        call this%AddPackageWriter(PkgWriter)
      endif
    case ('WEL')
      pw => WelWriter
      if (.not. this%ContainsPackageWriterOfType(pw)) then
        allocate(WelWriter)
        PkgWriter => WelWriter
        PkgWriter%ModelMovers => this%ModelMovers
        PkgWriter%source = 'WEL'
        PkgWriter%ModelBasename = this%BaseName
        allocate(PkgWriter%fileobj)
        call PkgWriter%fileobj%Initialize()
        PkgWriter%PkgType = 'WEL'
        call this%AddPackageWriter(PkgWriter)
      endif
    case ('HFB')
      pw => HfbWriter
      if (.not. this%ContainsPackageWriterOfType(pw)) then
        allocate(HfbWriter)
        PkgWriter => HfbWriter
        PkgWriter%ModelMovers => this%ModelMovers
        PkgWriter%source = 'HFB'
        PkgWriter%ModelBasename = this%BaseName
        allocate(PkgWriter%fileobj)
        call PkgWriter%fileobj%Initialize()
        PkgWriter%PkgType = 'HFB'
        call this%AddPackageWriter(PkgWriter)
      endif
    case default
      msg = 'In ModelType%AddPackage, PType not recognized: ' // trim(PType)
      call store_error(msg)
      call ustop()
    end select
    !
    return
  end function AddPackage

  subroutine AddPackageWriter(this, PkgWriter)
    implicit none
    ! dummy
    class(ModelType) :: this
    class(PackageWriterType), pointer :: PkgWriter
    ! local
    class(*), pointer :: obj => null()
    !
    PkgWriter%iulist = this%iulist
    obj => PkgWriter
    call this%PackageWriters%Add(obj)
    return
  end subroutine AddPackageWriter

  function GetPackageWriter(this, indx) result(pakwriter)
    implicit none
    ! dummy
    class(ModelType) :: this
    integer, intent(in) :: indx
    ! local
    class(*), pointer :: obj => null()
    class(PackageWriterType), pointer :: pakwriter
    !
    pakwriter => null()
    obj => this%PackageWriters%GetItem(indx)
    if (associated(obj)) then
      pakwriter => CastAsPackageWriterClass(obj)
    endif
    !
    return
  end function GetPackageWriter

  subroutine PackageConversionTable(this, iu)
    ! Write a table showing source and budget text for each package.
    implicit none
    ! dummy
    class(ModelType) :: this
    integer, intent(in) :: iu
    ! local
    integer :: i, npkg
    class(PackageWriterType), pointer :: pkg => null()
    character(len=LENPACKAGENAME) :: textb, texts
    ! formats
    1 format()
    5 format(1x,37('-'))
    10 format( &
    1x,'Data',    t13,'New'    ,t23,'Volume Budget' ,/, &
    1x,'Source',  t13,'Package',t23,'Table Text'  ,/, &
    1x,'--------',t13,'-------',t23,'----------------')
    15 format(1x,a)
    20 format(1x,a,t13,a,t23,a)
    30 format(1x,'Original name file: ',a)
    40 format(1x,'New name file: ',a)
    !
    if (iu>0) then
      write(iu,1)
      write(iu,15)'CONVERSION REPORT'
      write(iu,30)trim(this%NameFile2005)
      write(iu,40)trim(this%NameFile2015)
      write(iu,5)
      write(iu,10)
    else
      write(*,1)
      write(*,15)'CONVERSION REPORT'
      write(*,30)trim(this%NameFile2005)
      write(*,40)trim(this%NameFile2015)
      write(*,5)
      write(*,10)
    endif
    !
    if (this%NpfWriter%Active) then
      texts = this%NpfWriter%FlowPackage
      if (iu>0) then
        write(iu,20)trim(texts),'NPF','[no entry]'
      else
        write(*,20)trim(texts),'NPF','[no entry]'
      endif
    endif
    ! Special cases
    if (this%StoWriter%Active) then
      texts = this%NpfWriter%FlowPackage
      if (iu>0) then
        write(iu,20)trim(texts),'STO','STORAGE'
      else
        write(*,20)trim(texts),'STO','STORAGE'
      endif
    endif
    if (this%FhbWriter%Active) then
      if (associated(this%FhbWriter%WellWriter)) then
        call this%FhbWriter%WellWriter%GetBudgetText(textb)
        if (iu>0) then
          write(iu,20)'FHB','WEL',trim(textb)
        else
          write(*,20)'FHB','WEL',trim(textb)
        endif
      endif
      if (associated(this%FhbWriter%ChdWriter)) then
        call this%FhbWriter%ChdWriter%GetBudgetText(textb)
        if (iu>0) then
          write(iu,20)'FHB','CHD',trim(textb)
        else
          write(*,20)'FHB','CHD',trim(textb)
        endif
      endif
    endif
    call WritePackageEntry(this%IbChdWriter, iu)
    !
    ! Ordinary packages
    npkg = this%PackageWriters%Count()
    do i=1,npkg
      pkg => this%GetPackageWriter(i)
      if (associated(pkg)) then
        call WritePackageEntry(pkg, iu)
      endif
    enddo
    !
    if (iu>0) then
      write(iu,5)
!      write(iu,1)
    else
      write(*,5)
!      write(*,1)
    endif
    !
    return
  end subroutine PackageConversionTable

  subroutine WritePackageEntry(pkg, iu)
    implicit none
    ! dummy
    class(PackageWriterType) :: pkg
    integer, intent(in) :: iu
    ! local
    character(len=LENPACKAGENAME) :: textb, textn, texts
    ! formats
    20 format(1x,a,t13,a,t23,a)
    !
    if (pkg%Active) then
      texts = pkg%source
      textn = pkg%PkgType
      call pkg%GetBudgetText(textb)
      if (iu>0) then
        write(iu,20)trim(texts), textn, trim(textb)
      else
        write(*,20)trim(texts), textn, trim(textb)
      endif
    endif
    !
    return
  end subroutine WritePackageEntry

  subroutine AssignIgrid(this)
    implicit none
    ! dummy
    class(ModelType) :: this
    !
    if (this%IGrid < 1) then
      this%IGrid = GetNextIgrid()
    endif
    !
    return
  end subroutine AssignIgrid

  subroutine SavePointers(this)
    implicit none
    ! dummy
    class(ModelType) :: this
    ! local
    character(len=100) :: msg
    !
    if (this%IGrid <= 0) then
      msg = 'Error: IGrid not assigned for model from file: ' // trim(this%NameFile2005)
      call store_error(msg)
      call ustop()
    endif
    !
    call SGWF2BAS7PSV(this%IGrid)
    call SGWF2BCF7PSV(this%IGrid)
    call SGWF2CHD7PSV(this%IGrid)
    call SGWF2DRN7PSV(this%IGrid)
    call SGWF2ETS7PSV(this%IGrid)
    call SGWF2EVT7PSV(this%IGrid)
    call SGWF2FHB7PSV(this%IGrid)
    call SGWF2GHB7PSV(this%IGrid)
    call SGWF2HFB7PSV(this%IGrid)
    call SGWF2LAK7PSV(this%IGrid)
    call GWF2LPF7PSV(this%IGrid)
    call SGWF2RCH7PSV(this%IGrid)
    call SGWF2RES7PSV(this%IGrid)
    call SGWF2RIV7PSV(this%IGrid)
    call SGWF2SFR7PSV(this%IGrid)
    call SGWF2WEL7PSV(this%IGrid)
    !
    this%Nrow => NROW
    this%Ncol => NCOL
    this%nbotm => NBOTM
    this%delr => DELR
    this%delc => DELC
    this%lbotm => LBOTM
    this%botm => BOTM
    !
    return
  end subroutine SavePointers

  subroutine PointToGrid(this)
    implicit none
    ! dummy
    class(ModelType) :: this
    ! local
    character(len=100) :: msg
    !
    if (this%IGrid <= 0) then
      msg = 'Error: IGrid not assigned for model from file: ' // trim(this%NameFile2005)
      call store_error(msg)
      call ustop()
    endif
    !
    call SGWF2BAS7PNT(this%IGrid)
    call SGWF2LGR2PNT(this%IGrid)
    ! Ned todo: add more calls to *7PNT here
    !
    return
  end subroutine PointToGrid

  subroutine InitializeIdomain(this)
    implicit none
    ! dummy
    class(ModelType), intent(inout) :: this
    ! local
    integer :: i, ilay, j, k, idim1, idim2, idim3, kbot
    integer :: kk
    double precision :: b
    character(len=MAXCHARLEN) :: msg
    logical :: WetdryActive
    ! formats
    20 format(2x,a)
    30 format(a,1x,i0)
    40 format('Idomain for MODFLOW 6 model layer ',i0, &
        ', which replaces a quasi-3D unit, has been assigned equal' &
        ' to Ibound for layer ',i0,' of original model.')
    !
    ! Ensure IBOUND points to IBOUND of this model
    call this%PointToGrid()
    !
    WetdryActive = .false.
    ! Assign NpfWriter%WetDry pointer to appropriate WetDry array
    if (associated(WetdryBcf)) then
      this%NpfWriter%WetDry => WetdryBcf
      WetdryActive = .true.
    elseif (associated(WetdryLpf)) then
      this%NpfWriter%WetDry => WetdryLpf
      WetdryActive = .true.
    elseif (associated(WetdryUpw)) then
      this%NpfWriter%WetDry => WetdryUpw
      WetdryActive = .true.
    endif
    if (WetdryActive) then
      idim1 = size(this%NpfWriter%WetDry,1)
      idim2 = size(this%NpfWriter%WetDry,2)
      idim3 = size(this%NpfWriter%WetDry,3)
      if (idim1==1 .and. idim2==1 .and. idim3==1) then
        WetdryActive = .false.
      endif
    endif
    !
    ! Allocate and populate Idomain
    allocate(this%Idomain(ncol, nrow, this%Nlaynew))
    this%FhbWriter%Idomain => this%Idomain
    this%FhbWriter%ChdWriter%Idomain => this%Idomain
    this%IbChdWriter%Idomain => this%Idomain
    ilay = 0
    kbot = 0
    do k=1,nlay
      ilay = ilay + 1
      if (WetdryActive) then
        if (this%NpfWriter%FlowPackage == 'BCF') then
          if (laycon(k) == 1 .or. laycon(k) == 3) then
            kbot = kbot + 1
          endif
        elseif (this%NpfWriter%FlowPackage == 'LPF') then
          if (LAYTYP(k) /= 0 .and. LAYWET(k) /= 0) then
            kbot = kbot + 1
          endif
        endif
      endif
      do i=1,nrow
        do j=1,ncol
          ! Idomain < 0 does not have same meaning as IBOUND < 0
          if (IBOUND(j,i,k) /= 0) then
            this%Idomain(j,i,ilay) = 1
          else
            this%Idomain(j,i,ilay) = 0
            ! If rewetting is active and ibound=0, cell can be activated if wetdry /= 0.0
            if (WetdryActive) then
              if (kbot > 0 .and. kbot <= idim3) then
                if (this%NpfWriter%WetDry(j,i,kbot) /= DZERO) then
                  this%Idomain(j,i,ilay) = 1
                endif
              endif
            endif
          endif
        enddo
      enddo
      if (LAYCBD(k) /= 0) then
        ! Assign Idomain of nodes corresponding to underlying
        ! quasi-3d unit equal to Idomain of current layer.
        do i=1,nrow
          do j=1,ncol
            kk = LBOTM(k)
            b = BOTM(j,i,kk) - BOTM(j,i,kk+1)
            if (b > DZERO) then
              this%Idomain(j,i,kk+1) = this%Idomain(j,i,kk)
            else
              this%Idomain(j,i,kk+1) = -1
            end if
          enddo
        enddo
        ilay = ilay + 1
        write(msg,40)ilay,k
        call store_note(msg)
      endif
    enddo
    !
    !
    return
  end subroutine InitializeIdomain

  function ContainsPackageWriterOfType(this, pkgwriter) result (res)
    ! dummy
    class(ModelType), intent(inout) :: this
    class(PackageWriterType), pointer, intent(in) :: pkgwriter
    logical :: res
    ! local
    integer :: i, npkg
    class(PackageWriterType), pointer :: pw
    !
    res = .false.
    if (associated(pkgwriter)) then
      npkg = this%PackageWriters%Count()
      do i=1,npkg
        pw => this%GetPackageWriter(i)
        if (pw%SameType(pkgwriter)) then
          res = .true.
          exit
        endif
      enddo
    endif
    !
    return
  end function ContainsPackageWriterOfType

  subroutine BuildAllModelMovers(this)
    ! dummy
    class(ModelType) :: this
    ! local
    integer :: i, npkg
    class(PackageWriterType), pointer :: pkgwriter => null()
    !
    npkg = this%PackageWriters%Count()
    do i=1,npkg
      pkgwriter => this%GetPackageWriter(i)
      call pkgwriter%BuildModelMovers()
    enddo
    !
    return
  end subroutine BuildAllModelMovers

  subroutine RemoveFhbChdDuplicates(this)
    ! dummy
    class(ModelType) :: this
    ! local
    integer :: nchd
    !
    call this%FhbWriter%RemoveChdDuplicates(this%IbChdWriter%IbChdList)
    nchd = this%IbChdWriter%IbChdList%Count()
    if (nchd == 0) then
      this%IbChdWriter%Active = .false.
    endif
    !
    return
  end subroutine RemoveFhbChdDuplicates

  subroutine WriteMvrFile(this)
    ! dummy
    class(ModelType) :: this
    character(len=MAXCHARLEN) :: fname
    integer :: iu
    !
    call this%MvrWriter%WriteMvrFile()
    fname = this%MvrWriter%fileobj%FName
    iu = this%MvrWriter%fileobj%IUnit
    call this%Mf6Files%AddFile(fname, 'MVR6', iu, FCINPUT)
    !
    return
  end subroutine WriteMvrFile

  subroutine OrderPackageWriters(this)
    ! Ensure that LAK does not precede SFR in list of package writers
    ! dummy
    class(ModelType) :: this
    ! local
    class(PackageWriterType), pointer :: pkgWriter => null()
    integer :: i, indxLak, indxSfr, npkg
    class(*), pointer :: obj => null()
    logical :: destroyValue
    !
    indxLak = 0
    indxSfr = 0
    npkg = this%PackageWriters%Count()
    destroyValue = .false.
    !
    if (npkg > 1) then
      ! Determine positions of LAK and SFR, if present
      do i=1,npkg
        pkgWriter => this%GetPackageWriter(i)
        if (pkgWriter%PkgType == 'LAK') then
          indxLak = i
        elseif (pkgWriter%PkgType == 'SFR') then
          indxSfr = i
        endif
      enddo
      !
      ! Determine if reordering is needed
      if (indxLak > 0 .and. indxSfr > 0) then
        if (indxLak < indxSfr) then
          ! Reorder package writers
          pkgWriter => this%GetPackageWriter(indxLak)
          obj => pkgWriter
          call this%PackageWriters%InsertAfter(obj, indxSfr)
          call this%PackageWriters%RemoveNode(indxLak, destroyValue)
        endif
      endif
    endif
    !
    return
  end subroutine OrderPackageWriters

end module ModelModule
